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

/**
 * @brief Builds bidirectional CSR (Compressed Sparse Row) key offset vectors from edge weights.
 *
 * Computes prefix sums for both forward and backward adjacency lookups based on the
 * source and destination node IDs packed into 64-bit edge keys.
 *
 * @param weights A reference to the unordered map storing packed node pairs as keys and weights as values.
 * @return std::pair<std::vector<CSRMatrix::Offset>, std::vector<CSRMatrix::Offset>>
 *         A pair containing the forward and backward cumulative offset key vectors.
 */
std::pair<std::vector<CSRMatrix::Offset>, std::vector<CSRMatrix::Offset>>
    build_csr_key(const std::unordered_map<uint64_t, uint32_t> &weights) {

    // Initialize forward and backward key vectors with zeros (+1 size for standard CSR trailing offset)
    std::vector<CSRMatrix::Offset> forward_key(cfg::ARRAY_SIZE + 1, 0);
    std::vector<CSRMatrix::Offset> backward_key(cfg::ARRAY_SIZE + 1, 0);

    // Count incoming edges (destination node counts) and outgoing edges (source node counts)
    for (const auto& [key, weight] : weights) {
        // Unpack source (u) and destination (v) node IDs from the 64-bit edge key
        const auto u = static_cast<uint32_t>(key >> 32);
        const auto v = static_cast<uint32_t>(key & 0xFFFFFFFF);

        // Increment forward degree counter for destination node
        if (v < cfg::ARRAY_SIZE) {
            forward_key[v]++;
        }
        // Increment backward degree counter for source node
        if (u < cfg::ARRAY_SIZE) {
            backward_key[u]++;
        }
    }

    // Transform raw degree counts into cumulative offset tables (prefix sums)
    CSRMatrix::Offset running_forward = 0;
    CSRMatrix::Offset running_backward = 0;

    for (size_t i = 0; i < forward_key.size(); ++i) {
        const CSRMatrix::Offset degree_forward = forward_key[i];
        const CSRMatrix::Offset degree_backward = backward_key[i];

        // Assign current running offset
        forward_key[i] = running_forward;
        backward_key[i] = running_backward;

        // Accumulate running totals for the next index
        running_forward += degree_forward;
        running_backward += degree_backward;
    }

    return {std::move(forward_key), std::move(backward_key)};
}


/**
 * @brief Helper structure containing temporary buffers to facilitate filling a CSRMatrix.
 */
struct CsrMatrix {
    std::vector<CSRMatrix::NodeId> forward_nodes;
    std::vector<CSRMatrix::Weight> forward_weights;
    std::vector<CSRMatrix::NodeId> backward_nodes;
    std::vector<CSRMatrix::Weight> backward_weights;
};


/**
 * @brief Fills the bidirectional CSR matrix buffers for topological relaxation using edge weights.
 *
 * Iterates through the edge weight map, unpacking source and destination node IDs,
 * and populates the forward and backward adjacency nodes and weights buffers using
 * pre-calculated offset pointers.
 *
 * @param weights Unmapped collection of packed 64-bit node pairs and their corresponding weights.
 * @param forward_key Cumulative prefix sum offset vector for incoming edges.
 * @param backward_key Cumulative prefix sum offset vector for outgoing edges.
 * @return CsrMatrix A structure containing the allocated and populated forward and backward CSR vectors.
 */
CsrMatrix fill_csr_matrix(
    const std::unordered_map<uint64_t, uint32_t> &weights,
    const std::vector<CSRMatrix::Offset> &forward_key,
    const std::vector<CSRMatrix::Offset> &backward_key)
{
    const size_t n_edges = weights.size();

    // Initialize temporary CSR matrix buffers to store exact edge capacities
    CsrMatrix matrix;
    matrix.forward_nodes.resize(n_edges, 0);
    matrix.forward_weights.resize(n_edges, 0);
    matrix.backward_nodes.resize(n_edges, 0);
    matrix.backward_weights.resize(n_edges, 0);

    // Copy offset key tables to use as active writing cursors during population
    std::vector<CSRMatrix::Offset> forward_ptr = forward_key;
    std::vector<CSRMatrix::Offset> backward_ptr = backward_key;

    // Populate the adjacency matrices using the edge weights and offset cursors
    for (const auto& [key, weight] : weights) {

        // Unpack source (u) and destination (v) node IDs from the 64-bit edge key
        const auto u = static_cast<uint32_t>(key >> 32);
        const auto v = static_cast<uint32_t>(key & 0xFFFFFFFF);

        // Fill forward view (indexed by destination node v, storing source u)
        const CSRMatrix::Offset pos_forward = forward_ptr[v]++;
        matrix.forward_nodes[pos_forward] = u;
        matrix.forward_weights[pos_forward] = static_cast<CSRMatrix::Weight>(weight);

        // Fill backward view (indexed by source node u, storing destination v)
        const CSRMatrix::Offset pos_backward = backward_ptr[u]++;
        matrix.backward_nodes[pos_backward] = v;
        matrix.backward_weights[pos_backward] = static_cast<CSRMatrix::Weight>(weight);

    }

    return matrix;
}

/**
 * @brief Computes bidirectional alpha regularization factors for each node based on local graph topology and edge weight dominance.
 *
 * Mathematically, the alpha factor modulates trust between path-derived anchor coordinates and local
 * neighborhood smoothing. For each node $v$, the algorithm evaluates incoming/outgoing edge weight distributions
 * to compute a confidence score bounded within $[0.50, 0.95]$.
 *
 * The alpha calculation combines two primary metrics:
 * 1. **Edge Dominance (D):** The ratio of the maximum weight edge to the total sum of edge weights:
 *    D = [max(w_i)]/[sum(w_i)]
 *    Stronger dominant edges increase confidence (alpha values approach 0.95).
 *
 * 2. **Structural Complexity Penalty (C):** When node degree exceeds the expected number of unique haplotypes (N_haplo),
 *    extra connectivity indicates loops or structural variations, introducing positional uncertainty. The complexity penalty is defined as:
 *    C = max(0, (degree - N_haplo) / degree)
 *
 * The final clamped alpha formula balances dominance and complexity:
 *    alpha_v = clamp(0.50 + [0.45 * D * (1.0 - C)], ; 0.50, ; 0.95)
 *
 * Source and sink nodes lacking directional connectivity default to maximal freedom (alpha = 0.95).
 *
 * @param forward_key Cumulative prefix sum offset vector for incoming edges.
 * @param backward_key Cumulative prefix sum offset vector for outgoing edges.
 * @param forward_weight Buffer containing weights for incoming edges.
 * @param backward_weight Buffer containing weights for outgoing edges.
 * @return std::pair<std::vector<CSRMatrix::Alpha>, std::vector<CSRMatrix::Alpha>>
 *         A pair containing the forward and backward alpha factor vectors for all graph nodes.
 */
std::pair<std::vector<CSRMatrix::Alpha>, std::vector<CSRMatrix::Alpha>>
    compute_alpha_factors(
        const std::vector<CSRMatrix::Offset> &forward_key,
        const std::vector<CSRMatrix::Offset> &backward_key,
        const std::vector<CSRMatrix::Weight> &forward_weight,
        const std::vector<CSRMatrix::Weight> &backward_weight)
{
    // Initialize bidirectional alpha vectors with default maximum trust (0.95)
    std::vector forward_alpha(cfg::ARRAY_SIZE, 0.95f);
    std::vector backward_alpha(cfg::ARRAY_SIZE, 0.95f);

    // Baseline haplotype count used to scale structural complexity penalties
    const size_t n_haplo = cfg::N_HAPLO;

    // Iterate through all nodes to compute topological alpha trust factors
    for (size_t v = 1; v < cfg::ARRAY_SIZE; ++v) {

        // --- 1. Forward Alpha Calculation (Incoming Edges) ---
        const CSRMatrix::Offset low_slice_forward = forward_key[v];
        const CSRMatrix::Offset high_slice_forward = forward_key[v+1];
        const auto degrees_forward = static_cast<size_t> (high_slice_forward - low_slice_forward);

        // Process nodes with incoming connections
        if (degrees_forward > 0) {
            CSRMatrix::Weight sum_weight_forward = 0;
            CSRMatrix::Weight max_weight_forward = 0;

            // Aggregate edge weights to find total volume and maximum dominant edge weight
            for (CSRMatrix::Offset i = low_slice_forward; i < high_slice_forward; ++i) {
                const CSRMatrix::Weight weight = forward_weight[i];
                sum_weight_forward += weight;
                if (weight > max_weight_forward) {
                    max_weight_forward = weight;
                }

                // Compute edge dominance ratio (D)
                const float dominance_forward = static_cast<float>(max_weight_forward) /
                    static_cast<float>(sum_weight_forward);

                // Compute complexity penalty (C) if node degree exceeds expected haplotype count
                float complexity_forward = 0.0f;
                if (degrees_forward > n_haplo) {
                    complexity_forward =
                        static_cast<float>(degrees_forward - n_haplo) /
                        static_cast<float>(degrees_forward);
                }

                // Calculate final alpha bounded strictly within [0.50, 0.95]
                forward_alpha[v] = std::clamp(
                    0.50f + 0.45f * dominance_forward * (1.0f - complexity_forward),
                    0.50f, 0.95f
                );
            }
        }
        else {
            forward_alpha[v] = 0.95f;  // Source node: assign maximal positional freedom
        }

        // --- 2. Backward Alpha Calculation (Outgoing Edges) ---
        const CSRMatrix::Offset low_slice_backward = backward_key[v];
        const CSRMatrix::Offset high_slice_backward = backward_key[v+1];
        const auto degrees_backward = static_cast<float>(high_slice_backward - low_slice_backward);

        // Process nodes with outgoing connections
        if (degrees_backward > 0.0f)
        {
            CSRMatrix::Weight sum_weight_backward = 0;
            CSRMatrix::Weight max_weight_backward = 0;

            // Aggregate edge weights for outgoing directions
            for (CSRMatrix::Offset i = low_slice_backward; i < high_slice_backward; ++i) {
                const CSRMatrix::Weight weight = backward_weight[i];
                sum_weight_backward += weight;
                if (weight > max_weight_backward) {
                    max_weight_backward = weight;
                }

                // Compute edge dominance ratio (D)
                const float dominance_backward = static_cast<float>(max_weight_backward) /
                    static_cast<float>(sum_weight_backward);

                // Compute complexity penalty (C) for backward orientation
                float complexity_backward = 0.0f;
                if (degrees_backward > static_cast<float>(n_haplo)) {
                    complexity_backward =
                        (degrees_backward - static_cast<float>(n_haplo)) /
                        degrees_backward;
                }

                // Calculate final alpha bounded strictly within [0.50, 0.95]
                backward_alpha[v] = std::clamp(
                    0.50f + 0.45f * dominance_backward * (1.0f - complexity_backward),
                    0.50f, 0.95f
                );
            }
        }
        else {
            backward_alpha[v] = 0.95f;  // Sink node: assign maximal positional freedom
        }
    }

    return {std::move(forward_alpha), std::move(backward_alpha)};
}

/**
 * @brief Safely rounds and converts a vector of double-precision coordinates into a vector of 32-bit unsigned integers.
 *
 * Iterates through all node indices, validating that each coordinate is finite, non-negative,
 * and within the valid range of `uint32_t` capacity (excluding sentinel values like `NODE_UNSEEN_32`).
 * Non-existent nodes identified via `n_length` retain their `NODE_UNSEEN_32` sentinel.
 *
 * @param n_length Reference vector indicating active/valid nodes in the graph.
 * @param input Vector of double precision values (e.g., continuous linearized coordinates) to convert.
 * @return std::vector<uint32_t> A vector containing safely rounded 32-bit unsigned integer values.
 * @throws std::domain_error If any coordinate is non-finite (NaN/inf) or negative.
 * @throws std::overflow_error If any rounded coordinate exceeds the `uint32_t` capacity.
 */
std::vector<uint32_t> double2uint32(
    const std::vector<uint32_t>& n_length,
    const std::vector<double> &input) {

    // Initialize the result vector with the global unseen sentinel value
    std::vector result(cfg::ARRAY_SIZE, cfg::NODE_UNSEEN_32);

    for (size_t nid = 1; nid < cfg::ARRAY_SIZE; ++nid) {
        // Non-existent nodes remain NODE_UNSEEN_32 and are skipped
        if (n_length[nid] == cfg::NODE_UNSEEN_32) continue;

        const auto value = static_cast<double>(input[nid]);

        // Validate finiteness to prevent undefined behavior during rounding (NaN or infinity)
        if (!std::isfinite(value)) {
            throw std::domain_error(
                "float2uint32: non-finite linearized coordinate for node " +
                std::to_string(nid) + ": " + std::to_string(input[nid]));
        }

        // Round the continuous coordinate to the nearest integer
        const double rounded_value = std::round(value);

        // Ensure coordinates are non-negative
        if (rounded_value < 0.0) {
            throw std::domain_error(
                "float2uint32: negative linearized coordinate for node " +
                std::to_string(nid) + ": " + std::to_string(input[nid]));
        }

        // Ensure coordinates do not overflow uint32_t or clash with sentinel values
        if (rounded_value >= static_cast<double>(cfg::NODE_UNSEEN_32)) {
            throw std::overflow_error(
                "float2uint32: linearized coordinate exceeds uint32_t capacity "
                "for node " + std::to_string(nid) + ": " + std::to_string(input[nid]));
        }

        // Store safely casted integer coordinate
        result[nid] = static_cast<uint32_t>(rounded_value);
    }

    return result;
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


// Main function of the module
std::pair<std::vector<uint32_t>, size_t> relax_topology(
    const CSRMatrix& matrix,
    const std::vector<uint32_t>& n_length,
    const std::vector<uint32_t>& median_pos,
    const float convergence_threshold,
    const size_t max_iterations,
    const float lambda_factor) {

    // Initialize node positions from median anchor positions; maintain a copy for anchor regularization
    std::vector<double> nodes_pos(cfg::ARRAY_SIZE);
    for (size_t nid = 0; nid < cfg::ARRAY_SIZE; ++nid) {
        nodes_pos[nid] = static_cast<double>(median_pos[nid]);
    }
    const std::vector<double> anchor_pos = nodes_pos;

    // Reusable double-buffering vector for position updates during the hot-loop
    std::vector<double> next_node_pos = nodes_pos;

    // Bind CSR matrix view getters for high-performance access inside the hot-loop
    const auto& forward_key = matrix.forward_key();
    const auto& backward_key = matrix.backward_key();
    const auto& forward_nodes = matrix.forward_nodes();
    const auto& backward_nodes = matrix.backward_nodes();
    const auto& forward_weight = matrix.forward_weights();
    const auto& backward_weight = matrix.backward_weights();
    const auto& forward_alpha = matrix.forward_alphas();
    const auto& backward_alpha = matrix.backward_alphas();

    // Filter and cache metadata for active graph nodes to optimize iteration performance
    std::vector<ActiveNode> active_nodes;
    active_nodes.reserve(cfg::ARRAY_SIZE);

    for (uint32_t nid = 1; nid < cfg::ARRAY_SIZE; ++nid) {
        if (n_length[nid] == cfg::NODE_UNSEEN_32) continue; // Skip nodes that do not exist

        const CSRMatrix::Offset forward_slice_low = forward_key[nid];
        const CSRMatrix::Offset forward_slice_high = forward_key[nid + 1];
        const auto degrees_forward = forward_slice_high - forward_slice_low;

        const CSRMatrix::Offset backward_slice_low = backward_key[nid];
        const CSRMatrix::Offset backward_slice_high = backward_key[nid + 1];
        const auto degrees_backward = backward_slice_high - backward_slice_low;

        if (degrees_forward == 0 && degrees_backward == 0) { continue; } // Skip isolated nodes with no edges

        // Cache precomputed degrees, offsets, and alpha trust factors
        active_nodes.push_back(ActiveNode{
            nid,
            degrees_forward, degrees_backward,
            forward_alpha[nid], backward_alpha[nid],
            forward_slice_low, backward_slice_low
        });
    }

    // Local constant bindings for optimization
    const auto lambda = static_cast<double>(lambda_factor);
    const double one_minus_lambda = 1.0 - lambda;
    const size_t n_active = active_nodes.size();
    int executed_iteration = max_iterations;

    // Main relaxation hot-loop: processes active node buffers across iterations
    for (int iter = 1; iter <= executed_iteration; ++iter) {
        double total_movement = 0.0;
        double total_squared_mouvement = 0.0; // Tracks system residual energy per iteration

        #pragma omp parallel for reduction(+:total_movement, total_squared_mouvement) schedule(static) if(n_active > 1000)
        // Compute updated continuous coordinates for each active node in parallel
        for (size_t n = 0; n < n_active; ++n) {
            const auto& node = active_nodes[n];
            const uint32_t nid = node.nid;

            double forward_pred = anchor_pos[nid];
            double backward_pred = anchor_pos[nid];

            const uint32_t deg_forward = node.f_degree;
            const uint32_t deg_backward = node.b_degree;

            // 1. Compute forward coordinate prediction based on incoming neighbor barycenter
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

                // Apply forward barycenter and regularize with node alpha and anchor position
                if (sum_weight > 0) {
                    const double barycenter = sum_val / static_cast<double>(sum_weight);
                    const auto alpha_forward = static_cast<double>(node.f_alpha);
                    const double local_pred = alpha_forward * barycenter + (1.0 - alpha_forward) * nodes_pos[nid];
                    forward_pred = one_minus_lambda * local_pred + lambda * anchor_pos[nid];
                }
            }

            // 2. Compute backward coordinate prediction based on outgoing neighbor barycenter
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

                // Apply backward barycenter and regularize with node alpha and anchor position
                if (sum_weight > 0) {
                    const double barycenter = sum_val / static_cast<double>(sum_weight);
                    const auto alpha_backward = static_cast<double>(node.b_alpha);
                    const double local_pred = alpha_backward * barycenter + (1.0 - alpha_backward) * nodes_pos[nid];
                    backward_pred = one_minus_lambda * local_pred + lambda * anchor_pos[nid];
                }
            }

            // 3. Asymmetric fusion of bidirectional predictions for the current iteration
            double new_val;
            if (deg_forward > 0 && deg_backward > 0) { new_val = (forward_pred + backward_pred) * 0.5; }
            else if (deg_forward > 0)                { new_val = forward_pred; }
            else                                     { new_val = backward_pred; }

            // 4. Measure displacement energy and store updated position in back-buffer
            const double old_val = nodes_pos[nid];
            const double diff = new_val - old_val;

            total_movement += std::abs(diff);
            total_squared_mouvement += diff * diff;
            next_node_pos[nid] = new_val;
        }

        // Synchronize state buffers for the next iteration
        nodes_pos.swap(next_node_pos);
        const double active_count = static_cast<double>(std::max<std::size_t>(1, n_active));
        const double mean_displacement = total_movement / active_count;
        const double mean_squared_displacement = total_squared_mouvement / active_count;
        const double rms_displacement = std::sqrt(mean_squared_displacement);

        // Report iteration metrics
        std::cerr << "  - Iteration " << std::setw(3) << iter
                  << " | mean: "  << std::fixed << std::setprecision(8) << std::setw(14) << mean_displacement << " bp"
                  << " | RMS: "   << std::setw(14) << rms_displacement << " bp" << '\n';

        // Check convergence criterion against Root Mean Square displacement
        if (rms_displacement < convergence_threshold) {
            executed_iteration = iter;
            std::cerr << "  - Convergence reached at iteration " << iter
                      << ": RMS displacement " << std::fixed << std::setprecision(8) << rms_displacement
                      << " bp < threshold " << std::setprecision(6) << convergence_threshold << " bp." << '\n';
            break;
        }
    }

    // Safely convert final continuous coordinates into 32-bit unsigned integers
    std::vector<uint32_t> result = double2uint32(n_length, nodes_pos);
    return std::pair<std::vector<uint32_t>, size_t>(result, executed_iteration);
}