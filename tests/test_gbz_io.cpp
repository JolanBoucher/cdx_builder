#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <algorithm>
#include "gbz_IO.h"
#include <gbwt/gbwt.h>
#include <gbwt/dynamic_gbwt.h>
#include <gbwt/support.h>
#include <gbwtgraph/gbwtgraph.h>

// ============================================================================
// 1. EXTRACTION ROBUSTE DES NOMS DE CONTIGS (extract_compo_robust)
// ============================================================================

TEST(ExtractCompoTest, PlainContigAndEmpty) {
    EXPECT_EQ(detail::extract_compo_robust("chr12"), "chr12");
    EXPECT_EQ(detail::extract_compo_robust(""), "");
}

TEST(ExtractCompoTest, IncompleteAndFourFieldPanSN) {
    // PanSN à 2 champs (Sample#Contig)
    EXPECT_EQ(detail::extract_compo_robust("HG001#chr1"), "chr1");
    // PanSN à 4 champs (Sample#Haplotype#Contig#Frag)
    EXPECT_EQ(detail::extract_compo_robust("HG001#0#chr1#0"), "chr1");
}

TEST(ExtractCompoTest, SubrangeCleaning) {
    EXPECT_EQ(detail::extract_compo_robust("chr1:1000-2000"), "chr1");
    EXPECT_EQ(detail::extract_compo_robust("chr2[500-1500]"), "chr2");
}

// ============================================================================
// 2. TRI ET COMPARAISON BIOLOGIQUE (compare_compo_names)
// ============================================================================

TEST(CompareCompoNamesTest, NaturalNumericAndRefSeq) {
    // Tri numérique naturel
    EXPECT_TRUE(detail::compare_compo_names("contig2", "contig10"));
    EXPECT_FALSE(detail::compare_compo_names("contig10", "contig2"));

    // Identifiants RefSeq
    EXPECT_TRUE(detail::compare_compo_names("NC_000002.12", "NC_000010.11"));
}

TEST(CompareCompoNamesTest, CaseInsensitivityHandling) {
    // S'assurer que le tri ne crash pas sur la casse
    // Soit la fonction ignore la casse, soit elle impose un ordre déterministe
    bool a_before_b = detail::compare_compo_names("Chr1", "chr1");
    bool b_before_a = detail::compare_compo_names("chr1", "Chr1");

    // L'un doit être vrai et l'autre faux (ordre strict)
    EXPECT_NE(a_before_b, b_before_a);
}

TEST(CompareCompoNamesTest, MitochondrialPositioning) {
    // Les deux formes doivent être placées après les chromosomes sexuels
    EXPECT_TRUE(detail::compare_compo_names("chrY", "chrM"));
    EXPECT_TRUE(detail::compare_compo_names("chrY", "chrMT"));

    // Ordre déterministe entre chrM et chrMT pour préserver la stabilité de std::sort
    EXPECT_TRUE(detail::compare_compo_names("chrM", "chrMT"));
}

TEST(CompareCompoNamesTest, SpecialLociPositioning) {
    // Alternate loci juste après le chromosome principal
    EXPECT_TRUE(detail::compare_compo_names("chr1", "chr1_GL383518v1_alt"));
    EXPECT_TRUE(detail::compare_compo_names("chr1_GL383518v1_alt", "chr2"));

    // Unplaced (chrUn) tout à la fin
    EXPECT_TRUE(detail::compare_compo_names("chrMT", "chrUn"));
    EXPECT_TRUE(detail::compare_compo_names("chrY", "chrUn"));
}

// ============================================================================
// 3. TRAITEMENT MULTI-ÉCHANTILLONS ET COHÉRENCE DU PANGÉNOME
// ============================================================================

TEST(GetSortedComponentsTest, DuplicateContigsFromMultipleSamples) {
    std::vector<std::string> paths = {
        "HG001#1#chr1",
        "HG002#1#chr1",
        "HG003#1#chr1",
        "HG001#1#chr2"
    };

    std::vector<std::string> cleaned;
    for (const auto& p : paths) {
        cleaned.push_back(detail::extract_compo_robust(p));
    }

    // Dédoublonnage
    std::sort(cleaned.begin(), cleaned.end());
    cleaned.erase(std::unique(cleaned.begin(), cleaned.end()), cleaned.end());

    // Tri biologique
    std::sort(cleaned.begin(), cleaned.end(), detail::compare_compo_names);

    std::vector<std::string> expected = {"chr1", "chr2"};
    EXPECT_EQ(cleaned, expected);
}

TEST(ComponentValidationTest, DetectMixedContigsInComponent) {
    std::vector<std::string> contigs = {"chr1", "chr1", "chr2"};

    auto first = contigs.front();
    bool is_coherent = std::all_of(contigs.begin(), contigs.end(),
        [&](const std::string& c) { return c == first; });

    // Doit être faux car la composante contient un mélange de chr1 et chr2
    EXPECT_FALSE(is_coherent);
}

// -----------------------------------------------------------------------------
// Tests unitaires pour detail::pack_node_pair
// -----------------------------------------------------------------------------
TEST(PackNodePairTest, BitShiftCorrectness) {
    // Cas de base
    EXPECT_EQ(detail::pack_node_pair(1, 2), 0x0000000100000002ULL);

    // Cas limites (32-bit uint)
    handlegraph::nid_t max_val = 0xFFFFFFFF; // UINT32_MAX

    EXPECT_EQ(detail::pack_node_pair(0, 0), 0ULL);
    EXPECT_EQ(detail::pack_node_pair(0, max_val), 0x00000000FFFFFFFFULL);
    EXPECT_EQ(detail::pack_node_pair(max_val, 0), 0xFFFFFFFF00000000ULL);
    EXPECT_EQ(detail::pack_node_pair(max_val, max_val), 0xFFFFFFFFFFFFFFFFULL);
}

// =============================================================================
//      Test pour compute_edge_weights
// =============================================================================

// Helper Fixture
class ComputeEdgeWeightsTest : public ::testing::Test {
protected:
    gbwtgraph::GBWTGraph build_test_graph(
        const std::vector<std::vector<handlegraph::nid_t>>& paths,
        const std::vector<handlegraph::nid_t>& node_ids)
    {
        gbwtgraph::NaiveGraph naive_graph;

        for (handlegraph::nid_t nid : node_ids) {
            naive_graph.create_node(nid, "A");
        }

        // Calcul de la largeur exacte nécessaire pour encoder le plus grand ID orienté
        gbwt::size_type node_width = 1;

        if (!node_ids.empty()) {
            const handlegraph::nid_t max_nid = *std::max_element(node_ids.begin(), node_ids.end());
            const gbwt::node_type max_encoded_node = gbwt::Node::encode(max_nid, true);
            node_width = gbwt::bit_length(max_encoded_node);
        }

        gbwt::GBWTBuilder builder(node_width, 1024);

        for (const auto& path : paths) {
            gbwt::vector_type gbwt_path;
            gbwt_path.reserve(path.size());

            for (handlegraph::nid_t nid : path) {
                gbwt_path.push_back(gbwt::Node::encode(nid, false));
            }

            builder.insert(gbwt_path);

            gbwt::vector_type reverse_path;
            gbwt::reversePath(gbwt_path, reverse_path);
            builder.insert(reverse_path);
        }

        builder.finish();
        builder.index.header.set(gbwt::GBWTHeader::FLAG_BIDIRECTIONAL);

        gbwt_index = builder.index;

        return gbwtgraph::GBWTGraph(gbwt_index, naive_graph);
    }

    gbwtgraph::GBWTGraph build_test_graph_with_reverse(
    const std::vector<
        std::pair<handlegraph::nid_t, bool>
    >& path,
    const std::vector<handlegraph::nid_t>& node_ids)
    {
        gbwtgraph::NaiveGraph naive_graph;

        for (handlegraph::nid_t nid : node_ids) {
            naive_graph.create_node(nid, "A");
        }

        // Calcul de la largeur requise pour le plus grand nœud orienté.
        gbwt::size_type node_width = 1;

        if (!node_ids.empty()) {
            const handlegraph::nid_t max_nid =
                *std::max_element(
                    node_ids.begin(),
                    node_ids.end()
                );

            const gbwt::node_type max_encoded_node =
                gbwt::Node::encode(max_nid, true);

            node_width =
                gbwt::bit_length(max_encoded_node);
        }

        gbwt::GBWTBuilder builder(node_width, 1024);

        // Construction du chemin orienté fourni par le test.
        gbwt::vector_type gbwt_path;
        gbwt_path.reserve(path.size());

        for (const auto& [node_id, is_reverse] : path) {
            gbwt_path.push_back(
                gbwt::Node::encode(
                    node_id,
                    is_reverse
                )
            );
        }

        if (!gbwt_path.empty()) {
            // Chemin original, ici 2- -> 1-.
            builder.insert(gbwt_path);

            // Reverse-complément, ici 1+ -> 2+.
            gbwt::vector_type reverse_path;
            gbwt::reversePath(gbwt_path, reverse_path);

            builder.insert(reverse_path);
        }

        builder.finish();

        builder.index.header.set(
            gbwt::GBWTHeader::FLAG_BIDIRECTIONAL
        );

        gbwt_index = builder.index;

        return gbwtgraph::GBWTGraph(
            gbwt_index,
            naive_graph
        );
    }

    gbwt::GBWT gbwt_index;
};

// Tests unitaires
// CAS LIMIT 1 : Graphe vide (aucun nœud, aucun chemin)
TEST_F(ComputeEdgeWeightsTest, EmptyGraph) {
    // On passe une liste de chemins vide ET une liste de nœuds vide
    gbwtgraph::GBWTGraph empty_graph = build_test_graph({}, {});

    auto weights = compute_edge_weights(empty_graph);
    EXPECT_TRUE(weights.empty());
}

// CAS LIMIT 2 : Nœud isolé sans arêtes ni chemins GBWT
TEST_F(ComputeEdgeWeightsTest, IsolatedNodes) {
    gbwtgraph::GBWTGraph graph = build_test_graph({}, {1});
    auto weights = compute_edge_weights(graph);
    EXPECT_TRUE(weights.empty());
}

// CAS LIMIT 3 : Discontinuité dans les IDs des nœuds (Trou d'IDs : 1 -> 100)
TEST_F(ComputeEdgeWeightsTest, DiscontinuousNodeIDs) {
    gbwtgraph::GBWTGraph graph = build_test_graph({{1, 100}}, {1, 100});

    auto weights = compute_edge_weights(graph);
    const uint64_t key = detail::pack_node_pair(1, 100);

    EXPECT_EQ(weights.size(), 1);
    EXPECT_EQ(weights.count(key), 1);
    EXPECT_EQ(weights.at(key), 1);
}

// CAS LIMIT 4 : Poids multiples sur une même arête (Multiples haplotypes traversants)
TEST_F(ComputeEdgeWeightsTest, MultipleHaplotypesSameEdge) {
    gbwtgraph::GBWTGraph graph = build_test_graph({{1, 2}, {1, 2}, {1, 2}}, {1, 2});

    auto weights = compute_edge_weights(graph);

    uint64_t key = detail::pack_node_pair(1, 2);
    EXPECT_EQ(weights.count(key), 1);
    EXPECT_EQ(weights[key], 3);
}

TEST_F(ComputeEdgeWeightsTest,ReverseRepresentationIsNotStoredSeparately){
    // Chemin original : 2- -> 1-
    // Reverse-complément : 1+ -> 2+
    gbwtgraph::GBWTGraph graph =
        build_test_graph_with_reverse(
            {
                {2, true},
                {1, true}
            },
            {1, 2}
        );

    const auto weights = compute_edge_weights(graph);

    const std::uint64_t forward_key =
        detail::pack_node_pair(1, 2);

    const std::uint64_t reverse_key =
        detail::pack_node_pair(2, 1);

    ASSERT_EQ(weights.size(), 1);

    const auto forward_edge =
        weights.find(forward_key);

    ASSERT_NE(forward_edge, weights.end());
    EXPECT_EQ(forward_edge->second, 1);

    EXPECT_EQ(weights.count(reverse_key), 0);
}

// CAS LIMIT 6 : Auto-boucle (Self-loop : 1+ -> 1+)
TEST_F(ComputeEdgeWeightsTest, SelfLoopEdge) {
    gbwtgraph::GBWTGraph graph = build_test_graph({{1, 1}}, {1});

    auto weights = compute_edge_weights(graph);

    uint64_t self_key = detail::pack_node_pair(1, 1);
    EXPECT_EQ(weights.count(self_key), 1);
    EXPECT_EQ(weights[self_key], 1);
}

// CAS LIMIT 7 : Arêtes en bulle (Divergence / Convergence)
TEST_F(ComputeEdgeWeightsTest, BubbleDivergenceAndConvergence) {
    gbwtgraph::GBWTGraph graph = build_test_graph({{1, 2, 4}, {1, 3, 4}}, {1, 2, 3, 4});

    auto weights = compute_edge_weights(graph);

    EXPECT_EQ(weights.size(), 4);
    EXPECT_EQ(weights[detail::pack_node_pair(1, 2)], 1);
    EXPECT_EQ(weights[detail::pack_node_pair(1, 3)], 1);
    EXPECT_EQ(weights[detail::pack_node_pair(2, 4)], 1);
    EXPECT_EQ(weights[detail::pack_node_pair(3, 4)], 1);
}