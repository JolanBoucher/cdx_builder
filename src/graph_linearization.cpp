//
// Created by Jolan on 2026-07-24.
//

#include "graph_linearization.h"
#include "constant.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <boost/fusion/container/list/cons.hpp>
#include <boost/graph/filtered_graph.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif
static_assert(true, "_OPENMP test");

// this private function build bidirectional key of the CSRMatrix
// used for the topological relaxation
std::pair<std::vector<CSRMatrix::Offset>, std::vector<CSRMatrix::Offset>>
    build_csr_key(const std::unordered_map<uint64_t, uint32_t> &weights) {
    // initializing the key for both side
    std::vector<CSRMatrix::Offset> forward_key(cfg::ARRAY_SIZE + 1, 0);
    std::vector<CSRMatrix::Offset> backward_key(cfg::ARRAY_SIZE + 1, 0);

    // count the number of edge for each node
    for (const auto& [key, weight] : weights) {
        const auto u = static_cast<uint32_t>(key>>32);
        const auto v = static_cast<uint32_t>(key & 0xFFFFFFFF);

        if (v < cfg::ARRAY_SIZE) {
            forward_key[v]++;
        }
        if (u < cfg::ARRAY_SIZE) {
            backward_key[u]++;
        }
    }

    // transform key vectors into cumulative offset (prefix sum)
    CSRMatrix::Offset running_forward = 0;
    CSRMatrix::Offset running_backward = 0;

    for (size_t i = 0; i < forward_key.size(); ++i) {
        const CSRMatrix::Offset degree_forward = forward_key[i];
        const CSRMatrix::Offset degree_backward = backward_key[i];

        forward_key[i] = running_forward;
        backward_key[i] = running_backward;

        running_forward += degree_forward;
        running_backward += degree_backward;
    }
    return {std::move(forward_key), std::move(backward_key)};
}

// small structure to facilitate filling the CSRMatrix
struct CsrMatrix {
    std::vector<CSRMatrix::NodeId> forward_nodes;
    std::vector<CSRMatrix::Weight> forward_weights;
    std::vector<CSRMatrix::NodeId> backward_nodes;
    std::vector<CSRMatrix::Weight> backward_weights;
};

// private function to fill the csr matrix
// this matrix will be used for the topological relaxation
CsrMatrix fill_csr_matrix(
    const std::unordered_map<uint64_t, uint32_t> &weights,
    const std::vector<CSRMatrix::Offset> &forward_key,
    const std::vector<CSRMatrix::Offset> &backward_key)
{
    const size_t n_edges = weights.size();

    // initializing the data structure necessary to filling the matrix
    CsrMatrix matrix;
    matrix.forward_nodes.resize(n_edges, 0);
    matrix.forward_weights.resize(n_edges, 0);
    matrix.backward_nodes.resize(n_edges, 0);
    matrix.backward_weights.resize(n_edges, 0);

    // creating the copy that of the key that we will use as writing pointers
    std::vector<CSRMatrix::Offset> forward_ptr = forward_key;
    std::vector<CSRMatrix::Offset> backward_ptr = backward_key;

    // filling the matrix
    for (const auto& [key, weight] : weights) {

        // accessing the nid of the edges in the unordered_map
        const auto u = static_cast<uint32_t>(key>>32);
        const auto v = static_cast<uint32_t>(key & 0xFFFFFFFF);

        // forward view filling
        const CSRMatrix::Offset pos_forward = forward_ptr[v]++;
        matrix.forward_nodes[pos_forward] = u;
        matrix.forward_weights[pos_forward] = static_cast<CSRMatrix::Weight>(weight);

        // backward view filling
        const CSRMatrix::Offset pos_backward = backward_ptr[u]++;
        matrix.backward_nodes[pos_backward] = v;
        matrix.backward_weights[pos_backward] = static_cast<CSRMatrix::Weight>(weight);

    }
    return matrix;
}

// compute the alpha-factor (how much we should trust its position based on it's direct neighbor)
// of each individual node in both direction and fill the bidirectional alpha vectors
std::pair<std::vector<CSRMatrix::Alpha>, std::vector<CSRMatrix::Alpha>>
    compute_alpha_factors(
        const std::vector<CSRMatrix::Offset> &forward_key,
        const std::vector<CSRMatrix::Offset> &backward_key,
        const std::vector<CSRMatrix::Weight> &forward_weight,
        const std::vector<CSRMatrix::Weight> &backward_weight)
{
    // initialize the alpha's bidirectional vectors
    std::vector<CSRMatrix::Alpha> forward_alpha(cfg::ARRAY_SIZE, 0.95f);
    std::vector<CSRMatrix::Alpha> backward_alpha(cfg::ARRAY_SIZE, 0.95f);

    //local binding of the number of haplotype in the graph... used to fine-tune the alpha calculation
    const size_t n_haplo = cfg::N_HAPLO;

    // filling the alpha vectors
    for (size_t v = 1; v < cfg::ARRAY_SIZE; ++v) {
        // forward alpha filling
        // fetching the slice of the csr matrix corresponding to a specific nid
        const CSRMatrix::Offset low_slice_forward = forward_key[v];
        const CSRMatrix::Offset high_slice_forward = forward_key[v+1];
        const auto degrees_forward = static_cast<size_t> (high_slice_forward - low_slice_forward);

        // for slice that are not empty (no source nodes)
        if (degrees_forward > 0) {
            CSRMatrix::Weight sum_weight_forward = 0;
            CSRMatrix::Weight max_weight_forward = 0;

            // fetch the sum of weight entering this node and the max weight
            for (CSRMatrix::Offset i = low_slice_forward; i < high_slice_forward; ++i) {
                const CSRMatrix::Weight weight = forward_weight[i];
                sum_weight_forward += weight;
                if (weight > max_weight_forward) {
                    max_weight_forward = weight;
                }

                // get the ratio of the dominant incoming edge
                const float dominance_forward = static_cast<float>(max_weight_forward) /
                    static_cast<float>(sum_weight_forward);

                // if the node have extra connectivity it signals a loop in the graph
                // so we penalize it due to extra positional uncertainty created by the loop
                float complexity_forward = 0.0f;
                if (degrees_forward > n_haplo) {
                    complexity_forward =
                        static_cast<float>(degrees_forward - n_haplo) /
                        static_cast<float>(degrees_forward);
                }

                // Alpha remains within [0.50, 0.95].
                // Strong dominance increases alpha, while excess connectivity reduces it.
                forward_alpha[v] = std::clamp(
                    0.50f + 0.45f * dominance_forward * (1.0f - complexity_forward),
                    0.50f, 0.95f
                    );
            }
        }
        else {forward_alpha[v] = 0.95;}  // source node: maximal freedom

        // filling the backward alpha vector
        const CSRMatrix::Offset low_slice_backward = backward_key[v];
        const CSRMatrix::Offset high_slice_backward = backward_key[v+1];
        const auto degrees_backward = static_cast<float>(high_slice_backward);

        if (degrees_backward > 0.0f)
        {
            CSRMatrix::Weight sum_weight_backward = 0;
            CSRMatrix::Weight max_weight_backward = 0;

            // fetch the sum of weight entering this node and the max weight
            for (CSRMatrix::Offset i = low_slice_backward; i < high_slice_backward; ++i) {
                const CSRMatrix::Weight weight = backward_weight[i];
                sum_weight_backward += weight;
                if (weight > max_weight_backward) {
                    max_weight_backward = weight;
                }

                // get the ratio of the dominant incoming edge
                const float dominance_backward = static_cast<float>(max_weight_backward) /
                    static_cast<float>(sum_weight_backward);

                // if the node have extra connectivity it signals a loop in the graph
                // so we penalize it due to extra positional uncertainty created by the loop
                float complexity_backward = 0.0f;
                if (degrees_forward > n_haplo) {
                    complexity_backward =
                        (degrees_backward - static_cast<float>(n_haplo)) /
                        static_cast<float>(degrees_backward);
                }
                // Alpha remains within [0.50, 0.95].
                // Strong dominance increases alpha, while excess connectivity reduces it.
                backward_alpha[v] = std::clamp(
                    0.50f + 0.45f * dominance_backward * (1.0f - complexity_backward),
                    0.50f, 0.95f
                    );
            }
        }
        else {backward_alpha[v] = 0.95;}  // sink node: maximal freedom
    }
    return {std::move(forward_alpha), std::move(backward_alpha)};
}

// the public orchestator of construction of the linearization CSR matrix
CSRMatrix build_csr_matrix(const std::unordered_map<uint64_t, uint32_t> &weights) {

    // first pass building of the keys
    auto [forward_key, backward_key] =build_csr_key(weights);

    // second pas filling of the matrix
    auto matrix = fill_csr_matrix(weights, forward_key, backward_key);

    // lastly computing alpha of individual nodes in both direction
    auto [forward_alpha, backward_alpha] = compute_alpha_factors(
        forward_key,
        backward_key,
        matrix.forward_weights,
        matrix.backward_weights);

    // automatic construction and validation of the invariant
    return CSRMatrix(
        std::move(forward_key),
        std::move(matrix.forward_nodes),
        std::move(matrix.forward_weights),
        std::move(forward_alpha),
        std::move(backward_key),
        std::move(matrix.backward_nodes),
        std::move(matrix.backward_weights),
        std::move(backward_alpha)
    );
}

// this function round safely a vector of double into a vector of uint32_t
std::vector<uint32_t> double2uint32(
    const std::vector<uint32_t>& n_length,
    const std::vector<double> &input) {

    // initialize the result vector
    std::vector<uint32_t> result(cfg::ARRAY_SIZE, cfg::NODE_UNSEEN_32);
    for (size_t nid = 1; nid < cfg::ARRAY_SIZE; ++nid) {
        // Nonexistent nodes remain NODE_UNSEEN_32.
        if (n_length[nid] == cfg::NODE_UNSEEN_32) continue;
        const auto value = static_cast<double>(input[nid]);

        // some security for the rounding
        if (!std::isfinite(value)) {
            throw std::domain_error(
                "float2uint32: non-finite linearized coordinate for node " +
                std::to_string(nid) + ": " + std::to_string(input[nid]));
        }

        const double rounded_value = std::round(value);
        if (rounded_value < 0.0) {
            throw std::domain_error(
                "float2uint32: negative linearized coordinate for node " +
                std::to_string(nid) + ": " + std::to_string(input[nid]));
        }
        if (rounded_value >= static_cast<double>(cfg::NODE_UNSEEN_32)) {
            throw std::overflow_error(
                "float2uint32: linearized coordinate exceeds uint32_t capacity "
                "for node " + std::to_string(nid) + ": " + std::to_string(input[nid]));
        }
        result[nid] = static_cast<uint32_t>(rounded_value);
    }
    return result;
}


    std::pair<std::vector<uint32_t>, size_t> relax_topology(
    const CSRMatrix& matrix,
    const std::vector<uint32_t>& n_length,
    const std::vector<uint32_t>& median_pos,
    const float convergence_threshold,
    const size_t max_iterations,
    const float lambda_factor) {

    // copy the anchor vector. it will serve as the starting position of the node for the relaxation
    std::vector<double> nodes_pos(cfg::ARRAY_SIZE);
    for (size_t nid = 0; nid < cfg::ARRAY_SIZE; ++nid) {
        nodes_pos[nid] = static_cast<double>(median_pos[nid]);
    }
    const std::vector<double> anchor_pos = nodes_pos;

    // reusable buffer for the double-buffering in the hot-loop
    std::vector<double> next_node_pos = nodes_pos;

    // binding getters of csrMatrix for the hot-loop
    const auto& forward_key = matrix.forward_key();
    const auto& backward_key = matrix.backward_key();
    const auto& forward_nodes = matrix.forward_nodes();
    const auto& backward_nodes = matrix.backward_nodes();
    const auto& forward_weight = matrix.forward_weights();
    const auto& backward_weight = matrix.backward_weights();
    const auto& forward_alpha = matrix.forward_alphas();
    const auto& backward_alpha = matrix.backward_alphas();

    // early filtering of the node to calculate only on the active node in the vector
    std::vector<ActiveNode> active_nodes;
    active_nodes.reserve(cfg::ARRAY_SIZE);

    for (uint32_t nid = 1; nid < cfg::ARRAY_SIZE; ++nid) {
        if (n_length[nid] == cfg::NODE_UNSEEN_32) continue; // skip the nid associated to no nodes

        const CSRMatrix::Offset forward_slice_low = forward_key[nid];
        const CSRMatrix::Offset forward_slice_high = forward_key[nid + 1];
        const auto degrees_forward = forward_slice_high - forward_slice_low;

        const CSRMatrix::Offset backward_slice_low = backward_key[nid];
        const CSRMatrix::Offset backward_slice_high = backward_key[nid + 1];
        const auto degrees_backward = backward_slice_high - backward_slice_low;

        if (degrees_forward == 0 && degrees_backward == 0) {continue;} // skip node without edges

        //caching the metadata (precomputed degrees and offsets)
        active_nodes.push_back(ActiveNode{
            nid,
            degrees_forward, degrees_backward,
            forward_alpha[nid], backward_alpha[nid],
            forward_slice_low, backward_slice_low
        });
    }

        // Some more local binding
        const auto lambda = static_cast<double>(lambda_factor);
        const double one_minus_lambda = 1.0 - lambda;
        const size_t n_active = active_nodes.size();
        int executed_iteration = max_iterations;

    // the program function properly and the hot-loop
    // it will consume the cache ActiveNode N time
    for (int iter = 1; iter <= executed_iteration; ++iter) {
        double total_movement = 0.0;
        double total_squared_mouvement = 0.0; // mesure the system residual energy after each iteration

        #pragma omp parallel for reduction(+:total_movement, total_squared_mouvement) schedule(static) if(n_active > 1000)
        // calculated the new position for each active node
        for (size_t n = 0; n < n_active; ++n) {
            const auto& node = active_nodes[n];
            const uint32_t nid = node.nid;

            double forward_pred = anchor_pos[nid];
            double backward_pred = anchor_pos[nid];

            const uint32_t deg_forward = node.f_degree;
            const uint32_t deg_backward = node.b_degree;

            // 1. compute forward mouvement
            if (deg_forward > 0) {
                double sum_val = 0.0;
                uint32_t sum_weight = 0;

                const CSRMatrix::Offset forward_end = node.f_lo + deg_forward;
                for (CSRMatrix::Offset id = node.f_lo; id < forward_end; ++id) {
                    const uint32_t u = forward_nodes[id];
                    const uint32_t w = forward_weight[id];
                    sum_val += (nodes_pos[u] + static_cast<double>(n_length[u])) * static_cast<double>(w);
                    sum_weight += w;
                }

                // computing the mouvement in the forward direction for node 'nid'
                if (sum_weight > 0) {
                    const double barycenter = sum_val / static_cast<double>(sum_weight);
                    const auto alpha_forward = static_cast<double>(node.f_alpha);
                    const double local_pred = alpha_forward * barycenter + (1.0 - alpha_forward) * nodes_pos[nid];
                    forward_pred = one_minus_lambda * local_pred + lambda * anchor_pos[nid];
                }
            }

            // 2. compute backward mouvement
            if (deg_backward > 0) {
                double sum_val = 0.0;
                uint32_t sum_weight = 0;

                const CSRMatrix::Offset backward_end = node.b_lo + deg_backward;
                for (CSRMatrix::Offset id = node.b_lo; id < backward_end; ++id) {
                    const uint32_t s = backward_nodes[id];
                    const uint32_t w = backward_weight[id];
                    sum_val += (nodes_pos[s] - static_cast<double>(n_length[nid])) * static_cast<double>(w);
                    sum_weight += w;
                }
                if (sum_weight > 0) {
                    const double barycenter = sum_val / static_cast<double>(sum_weight);
                    const auto alpha_backward = static_cast<double>(node.b_alpha);
                    const double local_pred = alpha_backward * barycenter + (1.0 - alpha_backward) * nodes_pos[nid];
                    backward_pred = one_minus_lambda * local_pred + lambda * anchor_pos[nid];
                }
            }

            // 3. Asymmetric fusion of the bidirectional relaxation for this iteration
            double new_val;
            if (deg_forward > 0 && deg_backward > 0) {new_val = (forward_pred + backward_pred) * 0.5;}
            else if (deg_forward > 0)                {new_val = forward_pred;}
            else                                     {new_val = backward_pred;}

            // 4. Compute system energy and write in the reusable buffer
            const double old_val = nodes_pos[nid];
            const double diff = new_val - old_val;

            total_movement += std::abs(diff);
            total_squared_mouvement += diff * diff;
            next_node_pos[nid] = new_val;
        }

        // state synchronization
        nodes_pos.swap(next_node_pos);
        const double active_count = static_cast<double>(std::max<std::size_t>(1, n_active));
        const double mean_displacement = total_movement / active_count;
        const double mean_squared_displacement = total_squared_mouvement / active_count;

        const double rms_displacement = std::sqrt(mean_squared_displacement);

        std::cout << "Iteration " << std::setw(3) << iter
                  << " | mean: "  << std::fixed << std::setprecision(8) << std::setw(14) << mean_displacement << " bp"
                  << " | RMS: "   << std::setw(14) << rms_displacement << " bp" << '\n';

        if (rms_displacement < convergence_threshold) {
            executed_iteration = iter;
            std::cout << "Convergence reached at iteration " << iter
                      << ": RMS displacement " << std::fixed << std::setprecision(8) << rms_displacement
                      << " bp < threshold "<< std::setprecision(2)<< convergence_threshold << " bp." << '\n';
            break;
        }
    }

    std::vector<uint32_t> result = double2uint32(n_length, nodes_pos);
    return std::pair<std::vector<uint32_t>, size_t>(result, executed_iteration);
}