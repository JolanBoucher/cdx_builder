/**
 * @file test_gbz_io.cpp
 * @brief Unit test suite for the gbz_IO module (src/gbz_IO.h / gbz_IO.cpp).
 *
 * Covers, in order:
 *   - detail::extract_compo_robust   (PanSN / heuristic contig name extraction)
 *   - detail::compare_compo_names    (biological chromosome ordering comparator)
 *   - detail::pack_node_pair         (64-bit edge key packing)
 *   - compute_edge_weights           (path co-occurrence edge weighting)
 *   - assign_connected_components / get_graph_components (Union-Find over the graph)
 *   - compute_nodes_median_offset    (per-node median path offset)
 *   - bind_component_names           (component <-> contig name binding + biological reorder)
 *   - count_haplotypes               (unique (sample, haplotype) counting)
 *
 * Shared test infrastructure:
 *   - CfgFixture: saves/restores every mutable `cfg::` global around each
 *     test. Several functions under test read/write `cfg::ARRAY_SIZE`,
 *     `cfg::N_COMPO`, etc., and all tests share one gtest binary/process, so
 *     without this fixture a test could leak global state into the next one.
 *   - GbwtGraphTestFixture (extends CfgFixture): builds small, deterministic
 *     `gbwtgraph::GBZ` graphs from an explicit node list and a list of
 *     `PathSpec` path descriptions. Every path is inserted in both
 *     orientations (forward + reverse complement), as required for a valid
 *     bidirectional GBWT index, and is backed by a real `gbwt::Metadata`
 *     entry (sample/contig/haplotype/phase_block). This matters because
 *     `GBWTGraph::for_each_path_handle()`'s unfiltered fallback iterates
 *     `gbwt::Metadata::paths()` directly -- a GBWT with inserted sequences
 *     but no metadata entries yields *zero* iterable paths. Any test that
 *     walks paths (not just ones asserting on path names) therefore needs
 *     real metadata, not just inserted sequences.
 */

#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <unordered_map>
#include "gbz_IO.h"
#include "constant.h"
#include <gbwt/gbwt.h>
#include <gbwt/dynamic_gbwt.h>
#include <gbwt/support.h>
#include <gbwt/metadata.h>
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/gbz.h>

namespace {

// ============================================================================
// Shared fixtures
// ============================================================================

/**
 * @brief Saves and restores every mutable cfg:: global around each test.
 *
 * `constant.h` exposes process-wide mutable state (cfg::ARRAY_SIZE,
 * cfg::N_COMPO, cfg::N_HAPLO, cfg::NB_NODES) that several functions under
 * test depend on as an implicit precondition/postcondition. Since all
 * GoogleTest cases run in a single process, a test that sets cfg::ARRAY_SIZE
 * and forgets to reset it would silently corrupt every test that runs after
 * it. This fixture makes that impossible: tests are free to set whatever
 * cfg:: values they need in the test body, without caring what ran before or
 * will run after.
 */
class CfgFixture : public ::testing::Test {
protected:
    void SetUp() override {
        saved_nb_nodes_ = cfg::NB_NODES;
        saved_array_size_ = cfg::ARRAY_SIZE;
        saved_n_compo_ = cfg::N_COMPO;
        saved_n_haplo_ = cfg::N_HAPLO;
    }

    void TearDown() override {
        cfg::NB_NODES = saved_nb_nodes_;
        cfg::ARRAY_SIZE = saved_array_size_;
        cfg::N_COMPO = saved_n_compo_;
        cfg::N_HAPLO = saved_n_haplo_;
    }

private:
    size_t saved_nb_nodes_{};
    size_t saved_array_size_{};
    size_t saved_n_compo_{};
    size_t saved_n_haplo_{};
};

/// A single path to insert into a test graph, with the PanSN-relevant metadata attached to it.
struct PathSpec {
    std::string sample = "S";                                    // Sample name (metadata).
    std::string contig = "C";                                    // Contig/locus name (metadata).
    std::size_t haplotype = 0;                                    // Haplotype/phase number (metadata).
    std::size_t phase_block = 0;                                  // Phase block / fragment count (metadata).
    std::vector<std::pair<handlegraph::nid_t, bool>> steps;       // (node id, is_reverse) per step.
};

/// Convenience builder for the common case of an all-forward path.
std::vector<std::pair<handlegraph::nid_t, bool>> fwd(const std::vector<handlegraph::nid_t>& nodes) {
    std::vector<std::pair<handlegraph::nid_t, bool>> steps;
    steps.reserve(nodes.size());
    for (const handlegraph::nid_t nid : nodes) steps.emplace_back(nid, false);
    return steps;
}

/**
 * @brief Builds small, deterministic GBZ graphs for path-related and topology-related tests.
 *
 * See the file-level docstring for why every path is backed by real
 * `gbwt::Metadata`, even when the test does not care about path names.
 */
class GbwtGraphTestFixture : public CfgFixture {
protected:
    gbwtgraph::GBZ build_gbz(
        const std::vector<handlegraph::nid_t>& node_ids,
        const std::vector<PathSpec>& path_specs) {
        gbwtgraph::NaiveGraph naive_graph;
        for (const handlegraph::nid_t nid : node_ids) {
            naive_graph.create_node(nid, "A");
        }

        // Compute the exact bit width needed to encode the largest oriented node.
        gbwt::size_type node_width = 1;
        if (!node_ids.empty()) {
            const handlegraph::nid_t max_nid = *std::max_element(node_ids.begin(), node_ids.end());
            const gbwt::node_type max_encoded = gbwt::Node::encode(max_nid, true);
            node_width = gbwt::bit_length(max_encoded);
        }

        gbwt::GBWTBuilder builder(node_width, 1024);
        gbwt::Metadata metadata;

        // Collect unique sample/contig names in first-seen order to derive metadata indices.
        std::vector<std::string> sample_names;
        std::vector<std::string> contig_names;
        auto index_of = [](std::vector<std::string>& names, const std::string& name) -> std::size_t {
            const auto it = std::find(names.begin(), names.end(), name);
            if (it != names.end()) return static_cast<std::size_t>(it - names.begin());
            names.push_back(name);
            return names.size() - 1;
        };

        for (const auto& spec : path_specs) {
            gbwt::vector_type gbwt_path;
            gbwt_path.reserve(spec.steps.size());
            for (const auto& [nid, is_reverse] : spec.steps) {
                gbwt_path.push_back(gbwt::Node::encode(nid, is_reverse));
            }
            builder.insert(gbwt_path);

            // Bidirectional GBWT: also insert the reverse-complement orientation.
            gbwt::vector_type reverse_path;
            gbwt::reversePath(gbwt_path, reverse_path);
            builder.insert(reverse_path);

            const auto sample_idx = index_of(sample_names, spec.sample);
            const auto contig_idx = index_of(contig_names, spec.contig);
            metadata.addPath(
                static_cast<gbwt::size_type>(sample_idx),
                static_cast<gbwt::size_type>(contig_idx),
                static_cast<gbwt::size_type>(spec.haplotype),
                static_cast<gbwt::size_type>(spec.phase_block)
            );
        }

        metadata.setSamples(sample_names);
        metadata.setContigs(contig_names);
        metadata.setHaplotypes(sample_names.size());

        builder.finish();
        builder.index.header.set(gbwt::GBWTHeader::FLAG_BIDIRECTIONAL);
        builder.index.metadata = metadata;
        builder.index.addMetadata();

        return gbwtgraph::GBZ(builder.index, naive_graph);
    }

    /**
     * @brief Builds a Union-Find `parent`/`children` pair pre-initialized per the
     *        documented precondition of assign_connected_components()/get_graph_components():
     *        existing nodes must start as their own root, non-existent nodes must be
     *        marked with cfg::NODE_UNSEEN_32. Neither function initializes these arrays
     *        itself -- it is entirely the caller's responsibility.
     */
    static std::pair<std::vector<uint32_t>, std::vector<uint16_t>> make_identity_union_find(
        const gbwtgraph::GBWTGraph& graph, const std::size_t array_size) {
        std::vector<uint32_t> parent(array_size, cfg::NODE_UNSEEN_32);
        std::vector<uint16_t> children(array_size, 0);
        for (std::size_t nid = 0; nid < array_size; ++nid) {
            if (graph.has_node(static_cast<handlegraph::nid_t>(nid))) {
                parent[nid] = static_cast<uint32_t>(nid);
            }
        }
        return {std::move(parent), std::move(children)};
    }
};

} // namespace

// ============================================================================
// 1. extract_compo_robust
// ============================================================================
//
// Invariants / contract:
//   - PanSN naming (Sample#Haplotype#Contig[#Fragment] or Sample#Contig) is
//     parsed positionally by counting '#' separators; it is NOT validated
//     against real sample/contig dictionaries.
//   - Subrange suffixes ("[start-end]" or ":start-end") are stripped from
//     whatever the PanSN/heuristic step extracted as the contig name.
//   - When there is no '#' at all, a permissive, unanchored regex heuristic
//     is used as a fallback ("chrN", "chrM", "chrMT", or a bare 1-2 digit
//     number found *anywhere* in the string). Being unanchored, this
//     heuristic can match digits that are not really a chromosome number
//     (e.g. inside a scaffold/contig accession) -- this is documented
//     behavior, not a bug fix target, and is covered explicitly below so a
//     future change to the regex doesn't silently alter behavior.
//   - Empty input returns an empty string.

TEST(ExtractCompoTest, PlainContigAndEmpty) {
    EXPECT_EQ(detail::extract_compo_robust("chr12"), "chr12");
    EXPECT_EQ(detail::extract_compo_robust(""), "");
}

TEST(ExtractCompoTest, TwoFieldPanSN) {
    // PanSN à 2 champs (Sample#Contig): pas de 2e '#', branche "incomplete".
    EXPECT_EQ(detail::extract_compo_robust("HG001#chr1"), "chr1");
}

TEST(ExtractCompoTest, ThreeFieldPanSN) {
    // PanSN à 3 champs (Sample#Haplotype#Contig), sans 3e '#': third==npos,
    // donc compo = substr(second+1) jusqu'à la fin. Cas intermédiaire entre
    // le 2-champs et le 4-champs, absent de la suite d'origine.
    EXPECT_EQ(detail::extract_compo_robust("HG001#1#chr1"), "chr1");
}

TEST(ExtractCompoTest, FourFieldPanSN) {
    // PanSN à 4 champs (Sample#Haplotype#Contig#Fragment).
    EXPECT_EQ(detail::extract_compo_robust("HG001#0#chr1#0"), "chr1");
}

TEST(ExtractCompoTest, SubrangeCleaning) {
    EXPECT_EQ(detail::extract_compo_robust("chr1:1000-2000"), "chr1");
    EXPECT_EQ(detail::extract_compo_robust("chr2[500-1500]"), "chr2");
}

TEST(ExtractCompoTest, PanSNCombinedWithSubrange) {
    // La sous-plage doit être nettoyée APRÈS extraction du champ PanSN.
    EXPECT_EQ(detail::extract_compo_robust("HG001#1#chr1:100-200"), "chr1");
}

TEST(ExtractCompoTest, FallbackRegexCaseInsensitivePreservesCase) {
    // Le regex heuristique est insensible à la casse pour la détection, mais
    // `match.str(1)` retourne la sous-chaîne d'origine (casse non normalisée).
    EXPECT_EQ(detail::extract_compo_robust("CHR5"), "CHR5");
}

TEST(ExtractCompoTest, FallbackRegexNoDigitsNoChrPrefixReturnsWholeString) {
    // Aucun motif "chr..." ni chiffre isolé de 1-2 caractères -> la chaîne
    // nettoyée est retournée telle quelle.
    EXPECT_EQ(detail::extract_compo_robust("unplaced_contig"), "unplaced_contig");
}

TEST(ExtractCompoTest, FallbackRegexIsUnanchoredAndCanMatchStrayDigits) {
    // Documente une limite connue de l'heuristique: en l'absence de '#' et de
    // préfixe "chr", un nombre de 1-2 chiffres n'importe où dans la chaîne
    // est capturé, même s'il ne s'agit pas d'un numéro de chromosome.
    EXPECT_EQ(detail::extract_compo_robust("scaffold_5"), "5");
}

// ============================================================================
// 2. compare_compo_names
// ============================================================================
//
// Invariants / contract:
//   - Must be usable as a strict-weak-ordering comparator for std::sort:
//     irreflexive (compare(a,a) == false) and consistent enough not to
//     trigger UB on realistic chromosome-name sets.
//   - Ordering priority: standard/sex chromosomes < mitochondrial < unplaced.
//   - Among standard chromosomes, ordering is numeric (chr2 < chr10), not
//     lexicographic.
//   - The numeric extraction regex `(?:chr)?([0-9]{1,2})` matches digits with
//     or without a "chr" prefix, and only kicks in when *both* sides yield a
//     parsable number; otherwise names fall back to plain lexicographic
//     comparison, which is why e.g. "chrX" doesn't reliably fall in numeric
//     chromosome order relative to "chr2"/"chr10".
//   - When both names have the same chromosome number, the non-"_alt" one
//     sorts first.

TEST(CompareCompoNamesTest, NaturalNumericAndRefSeq) {
    // Tri numérique naturel
    EXPECT_TRUE(detail::compare_compo_names("contig2", "contig10"));
    EXPECT_FALSE(detail::compare_compo_names("contig10", "contig2"));

    // Identifiants RefSeq
    EXPECT_TRUE(detail::compare_compo_names("NC_000002.12", "NC_000010.11"));
}

TEST(CompareCompoNamesTest, NumericHeuristicWorksWithoutChrPrefix) {
    // Le regex d'extraction numérique fonctionne même sans préfixe "chr".
    EXPECT_TRUE(detail::compare_compo_names("scaffold_1", "scaffold_2"));
    EXPECT_FALSE(detail::compare_compo_names("scaffold_2", "scaffold_1"));
}

TEST(CompareCompoNamesTest, CaseInsensitivityHandling) {
    // S'assurer que le tri ne crash pas sur la casse
    // Soit la fonction ignore la casse, soit elle impose un ordre déterministe
    bool a_before_b = detail::compare_compo_names("Chr1", "chr1");
    bool b_before_a = detail::compare_compo_names("chr1", "Chr1");

    // L'un doit être vrai et l'autre faux (ordre strict)
    EXPECT_NE(a_before_b, b_before_a);
}

TEST(CompareCompoNamesTest, IsIrreflexive) {
    // Précondition d'un comparateur std::sort valide: compare(a, a) == false.
    for (const std::string& name : {std::string("chr1"), std::string("chrM"),
                                     std::string("chrUn_1"), std::string("chr1_alt")}) {
        EXPECT_FALSE(detail::compare_compo_names(name, name)) << "name = " << name;
    }
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

TEST(CompareCompoNamesTest, SortsMixedRealisticSetIntoExpectedOrder) {
    // Intégration: compose les invariants ci-dessus via un vrai std::sort et
    // vérifie qu'il n'y a pas de crash (UB de comparateur mal formé) en plus
    // de l'ordre final. Ensemble volontairement restreint à des paires déjà
    // validées individuellement plus haut, pour ne pas deviner un ordre total
    // non vérifié (chrX/chrY notamment retombent sur le fallback
    // lexicographique et ne sont pas inclus ici).
    std::vector<std::string> names = {"chrUn_1", "chrM", "chr10", "chr1", "chr2"};
    std::sort(names.begin(), names.end(), detail::compare_compo_names);

    const std::vector<std::string> expected = {"chr1", "chr2", "chr10", "chrM", "chrUn_1"};
    EXPECT_EQ(names, expected);
    EXPECT_TRUE(std::is_sorted(names.begin(), names.end(), detail::compare_compo_names));
}

// ============================================================================
// 3. EMPAQUETAGE DE PAIRES DE NOEUDS (detail::pack_node_pair)
// ============================================================================
//
// Invariant: combine deux IDs 32 bits (src, dst) en une clé 64 bits unique et
// réversible (src dans les 32 bits hauts, dst dans les 32 bits bas), sans
// validation de plage (suppose des entrées 32 bits valides).

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

// ============================================================================
// 4. POIDS DES ARETES PAR CO-OCCURRENCE (compute_edge_weights)
// ============================================================================
//
// Invariants / contract:
//   - Weights count consecutive-step co-occurrences per path, in the
//     orientation the path is stored in (forward view only: a path visiting
//     u then v increments weight[(u,v)], never weight[(v,u)] for the same
//     traversal).
//   - Because every test path is inserted together with its reverse
//     complement (see GbwtGraphTestFixture::build_gbz), the reverse-
//     complement traversal must NOT produce a separate reverse edge key --
//     it is expected to fold back onto the same node pair in the opposite
//     direction of travel, not double the forward key's weight.
//   - A path of a single node contributes no edges (no consecutive pair to
//     observe).
//   - Self-loops and bubble topologies (divergence/convergence) must be
//     tallied correctly, independent of node ID contiguity.

TEST_F(GbwtGraphTestFixture, EmptyGraph) {
    // CAS LIMITE: graphe vide (aucun nœud, aucun chemin)
    const auto gbz = build_gbz({}, {});
    const auto weights = compute_edge_weights(gbz.graph);
    EXPECT_TRUE(weights.empty());
}

TEST_F(GbwtGraphTestFixture, IsolatedNodeNoPaths) {
    // CAS LIMITE: nœud isolé, aucun chemin GBWT du tout.
    const auto gbz = build_gbz({1}, {});
    const auto weights = compute_edge_weights(gbz.graph);
    EXPECT_TRUE(weights.empty());
}

TEST_F(GbwtGraphTestFixture, SingleNodePathContributesNoEdges) {
    // CAS LIMITE: un chemin d'une seule étape n'a pas de paire consécutive à observer.
    const auto gbz = build_gbz({1}, {PathSpec{"S", "C", 0, 0, fwd({1})}});
    const auto weights = compute_edge_weights(gbz.graph);
    EXPECT_TRUE(weights.empty());
}

TEST_F(GbwtGraphTestFixture, DiscontinuousNodeIDs) {
    // CAS LIMITE: discontinuité dans les IDs des nœuds (trou d'IDs: 1 -> 100)
    const auto gbz = build_gbz({1, 100}, {PathSpec{"S", "C", 0, 0, fwd({1, 100})}});
    const auto weights = compute_edge_weights(gbz.graph);
    const uint64_t key = detail::pack_node_pair(1, 100);

    EXPECT_EQ(weights.size(), 1u);
    EXPECT_EQ(weights.count(key), 1u);
    EXPECT_EQ(weights.at(key), 1u);
}

TEST_F(GbwtGraphTestFixture, MultipleHaplotypesSameEdge) {
    // CAS LIMITE: poids multiples sur une même arête (plusieurs haplotypes traversants)
    const auto gbz = build_gbz({1, 2}, {
        PathSpec{"S1", "C", 0, 0, fwd({1, 2})},
        PathSpec{"S2", "C", 0, 0, fwd({1, 2})},
        PathSpec{"S3", "C", 0, 0, fwd({1, 2})},
    });
    const auto weights = compute_edge_weights(gbz.graph);

    const uint64_t key = detail::pack_node_pair(1, 2);
    EXPECT_EQ(weights.count(key), 1u);
    EXPECT_EQ(weights.at(key), 3u);
}

TEST_F(GbwtGraphTestFixture, ReverseRepresentationIsNotStoredSeparately) {
    // Original path (authored, inserted first -> this is sequence id 2p,
    // i.e. what compute_edge_weights actually extracts via
    // gbwt::Path::encode(p, false)) : 2- -> 1-.
    // Reverse complement (auto-inserted by the fixture right after,
    // sequence id 2p+1, never extracted directly by compute_edge_weights) :
    // 1+ -> 2+.
    const auto gbz = build_gbz({1, 2}, {
        PathSpec{"S", "C", 0, 0, {{2, true}, {1, true}}}
    });

    const auto weights = compute_edge_weights(gbz.graph);

    // "Forward" here means: the direction of the authored traversal (2 -> 1),
    // which is what gets counted. "Reverse" is the node pair that would only
    // appear if the auto-generated reverse-complement sequence were
    // (incorrectly) walked and counted separately.
    const std::uint64_t forward_key = detail::pack_node_pair(2, 1);
    const std::uint64_t reverse_key = detail::pack_node_pair(1, 2);

    ASSERT_EQ(weights.size(), 1u);

    const auto forward_edge = weights.find(forward_key);
    ASSERT_NE(forward_edge, weights.end());
    EXPECT_EQ(forward_edge->second, 1u);

    EXPECT_EQ(weights.count(reverse_key), 0u);
}

TEST_F(GbwtGraphTestFixture, SelfLoopEdge) {
    // CAS LIMITE: auto-boucle (Self-loop: 1+ -> 1+)
    const auto gbz = build_gbz({1}, {PathSpec{"S", "C", 0, 0, fwd({1, 1})}});
    const auto weights = compute_edge_weights(gbz.graph);

    const uint64_t self_key = detail::pack_node_pair(1, 1);
    EXPECT_EQ(weights.count(self_key), 1u);
    EXPECT_EQ(weights.at(self_key), 1u);
}

TEST_F(GbwtGraphTestFixture, BubbleDivergenceAndConvergence) {
    // CAS LIMITE: arêtes en bulle (divergence / convergence)
    const auto gbz = build_gbz({1, 2, 3, 4}, {
        PathSpec{"S1", "C", 0, 0, fwd({1, 2, 4})},
        PathSpec{"S2", "C", 0, 0, fwd({1, 3, 4})},
    });
    const auto weights = compute_edge_weights(gbz.graph);

    EXPECT_EQ(weights.size(), 4u);
    EXPECT_EQ(weights.at(detail::pack_node_pair(1, 2)), 1u);
    EXPECT_EQ(weights.at(detail::pack_node_pair(1, 3)), 1u);
    EXPECT_EQ(weights.at(detail::pack_node_pair(2, 4)), 1u);
    EXPECT_EQ(weights.at(detail::pack_node_pair(3, 4)), 1u);
}

// ============================================================================
// 5. UNION-FIND / COMPOSANTES CONNEXES (assign_connected_components, get_graph_components)
// ============================================================================
//
// Invariants / contract:
//   - CRITICAL PRECONDITION (undocumented by the type system): neither
//     function initializes `parent`/`children`. The caller must pre-fill
//     `parent[nid] = nid` for every node that exists in the graph and
//     `cfg::NODE_UNSEEN_32` for every index that does not. Violating this
//     precondition is undefined/nonsensical (path-halving would walk into
//     garbage indices). See GbwtGraphTestFixture::make_identity_union_find().
//   - Union-Find topology is derived purely from graph edges (both
//     `follow_edges` directions), independent of any GBWT metadata/paths --
//     an isolated node with no edges remains its own singleton component.
//   - `get_graph_components` compacts sparse roots into dense IDs in
//     [0, N_COMPO) and writes the discovered component count into
//     `cfg::N_COMPO` as a side effect.
//   - Component *identity* (which dense ID a given node ends up with) is not
//     guaranteed to follow any particular order, since it depends on
//     `for_each_handle`'s (implementation-defined) traversal order. Tests
//     therefore assert structural invariants (same/different component,
//     dense ID range, count) rather than hard-coded IDs.

TEST_F(GbwtGraphTestFixture, SingleIsolatedNodeIsItsOwnComponent) {
    cfg::ARRAY_SIZE = 2; // indices [0,1], node 1 exists
    // GBWTGraph only considers a node "real" (has_node() == true) once it is
    // referenced by at least one indexed sequence, so an isolated node still
    // needs a (trivial, edge-less) single-node path to register at all.
    const auto gbz = build_gbz({1}, {PathSpec{"S", "C", 0, 0, fwd({1})}});
    auto [parent, children] = make_identity_union_find(gbz.graph, cfg::ARRAY_SIZE);

    const auto compo = get_graph_components(gbz.graph, parent, children);

    EXPECT_EQ(cfg::N_COMPO, 1u);
    EXPECT_NE(compo[1], cfg::NODE_UNSEEN_16);
}

TEST_F(GbwtGraphTestFixture, TwoDisjointChainsFormTwoComponents) {
    cfg::ARRAY_SIZE = 6; // indices [0,5]
    const auto gbz = build_gbz({1, 2, 3, 4, 5}, {
        PathSpec{"S1", "C", 0, 0, fwd({1, 2, 3})},
        PathSpec{"S2", "C", 0, 0, fwd({4, 5})},
    });
    auto [parent, children] = make_identity_union_find(gbz.graph, cfg::ARRAY_SIZE);

    const auto compo = get_graph_components(gbz.graph, parent, children);

    EXPECT_EQ(cfg::N_COMPO, 2u);
    EXPECT_EQ(compo[1], compo[2]);
    EXPECT_EQ(compo[2], compo[3]);
    EXPECT_EQ(compo[4], compo[5]);
    EXPECT_NE(compo[1], compo[4]);
}

TEST_F(GbwtGraphTestFixture, CycleCollapsesIntoOneComponent) {
    cfg::ARRAY_SIZE = 4; // indices [0,3]
    // 1 -> 2 -> 3 -> 1: a cycle should not confuse Union-Find into multiple roots.
    const auto gbz = build_gbz({1, 2, 3}, {
        PathSpec{"S", "C", 0, 0, fwd({1, 2, 3, 1})},
    });
    auto [parent, children] = make_identity_union_find(gbz.graph, cfg::ARRAY_SIZE);

    const auto compo = get_graph_components(gbz.graph, parent, children);

    EXPECT_EQ(cfg::N_COMPO, 1u);
    EXPECT_EQ(compo[1], compo[2]);
    EXPECT_EQ(compo[2], compo[3]);
}

TEST_F(GbwtGraphTestFixture, SelfLoopStaysSingleton) {
    cfg::ARRAY_SIZE = 2;
    const auto gbz = build_gbz({1}, {PathSpec{"S", "C", 0, 0, fwd({1, 1})}});
    auto [parent, children] = make_identity_union_find(gbz.graph, cfg::ARRAY_SIZE);

    const auto compo = get_graph_components(gbz.graph, parent, children);

    EXPECT_EQ(cfg::N_COMPO, 1u);
    EXPECT_NE(compo[1], cfg::NODE_UNSEEN_16);
}

TEST_F(GbwtGraphTestFixture, ComponentIdsAreDenseAndContiguous) {
    cfg::ARRAY_SIZE = 6;
    const auto gbz = build_gbz({1, 2, 3, 4, 5}, {
        PathSpec{"S1", "C", 0, 0, fwd({1, 2})},
        PathSpec{"S2", "C", 0, 0, fwd({3, 4})},
        PathSpec{"S3", "C", 0, 0, fwd({5})},
    });
    auto [parent, children] = make_identity_union_find(gbz.graph, cfg::ARRAY_SIZE);

    const auto compo = get_graph_components(gbz.graph, parent, children);

    ASSERT_EQ(cfg::N_COMPO, 3u);
    std::vector<bool> seen(cfg::N_COMPO, false);
    for (const std::uint16_t cid : {compo[1], compo[3], compo[5]}) {
        ASSERT_LT(cid, cfg::N_COMPO);
        EXPECT_FALSE(seen[cid]) << "duplicate dense component id " << cid;
        seen[cid] = true;
    }
    EXPECT_TRUE(std::all_of(seen.begin(), seen.end(), [](bool b) { return b; }));
}

// ============================================================================
// 6. DECALAGE MEDIAN PAR NOEUD (compute_nodes_median_offset)
// ============================================================================
//
// Invariants / contract:
//   - For each node, collects every path-offset at which it is visited
//     (across ALL paths, not per-path), then picks the "right" median: for
//     `n` observations, `sorted_offsets[n / 2]` (i.e. for an even count this
//     is the upper of the two middle values, per `std::nth_element`'s exact
//     partition point used in the implementation).
//   - A node that exists in the graph but is never visited by any path
//     keeps the `cfg::NODE_UNSEEN_32` sentinel.
//   - A node visited multiple times within a single path (e.g. a loop)
//     contributes one offset entry per visit.

TEST_F(GbwtGraphTestFixture, SingleOccurrenceOffsetIsExact) {
    cfg::ARRAY_SIZE = 3;
    const auto gbz = build_gbz({1, 2}, {PathSpec{"S", "C", 0, 0, fwd({1, 2})}});

    const auto medians = compute_nodes_median_offset(gbz.graph);

    EXPECT_EQ(medians[1], 0u); // node 1 is the first step: offset 0
    // node 2's offset equals node 1's length; sequences are single-char "A" here.
    EXPECT_NE(medians[2], cfg::NODE_UNSEEN_32);
}

TEST_F(GbwtGraphTestFixture, NodeAbsentFromAnyPathKeepsSentinel) {
    cfg::ARRAY_SIZE = 3;
    // Node 2 exists in the graph but is never traversed by a path.
    const auto gbz = build_gbz({1, 2}, {PathSpec{"S", "C", 0, 0, fwd({1})}});

    const auto medians = compute_nodes_median_offset(gbz.graph);

    EXPECT_EQ(medians[2], cfg::NODE_UNSEEN_32);
}

TEST_F(GbwtGraphTestFixture, MedianAggregatesOffsetsAcrossMultiplePaths) {
    cfg::ARRAY_SIZE = 4; // indices [0,3]: nodes go up to id 3 here, unlike the two tests above
    // Node 2 is visited by three paths with different prefix lengths, so it
    // is seen at three different offsets.
    const auto gbz = build_gbz({1, 2, 3}, {
        PathSpec{"S1", "C", 0, 0, fwd({2})},          // node 2 at offset 0
        PathSpec{"S2", "C", 0, 0, fwd({1, 2})},        // node 2 at offset len(1)
        PathSpec{"S3", "C", 0, 0, fwd({1, 3, 2})},     // node 2 at offset len(1)+len(3)
    });

    const auto medians = compute_nodes_median_offset(gbz.graph);

    // Right median of 3 sorted offsets {0, L, 2L} (L=1 since sequences are "A")
    // is sorted_offsets[3/2] = sorted_offsets[1] = L = 1.
    EXPECT_EQ(medians[2], 1u);
}

// ============================================================================
// 7. LIAISON NOM DE COMPOSANTE <-> ID DE COMPOSANTE (bind_component_names)
// ============================================================================
//
// Invariants / contract:
//   - `nid2compo` is an INPUT here (pre-assigned component ids, e.g. from
//     get_graph_components), not something bind_component_names computes.
//   - Every step of every path must resolve to a single, consistent
//     component id; a path whose steps span more than one component is
//     considered corrupt data and throws std::runtime_error.
//   - Every connected component must be "covered" by at least one path
//     (i.e. get a name); an uncovered component throws.
//   - Two different paths landing in the same component must agree on the
//     extracted contig name; disagreement throws.
//   - `max_steps_to_check` caps how many steps of a path are inspected for
//     the single-component invariant. This is a real footgun: with a low
//     cap, a path that actually crosses components can silently pass
//     validation as long as its *checked* prefix stays within one
//     component and its extracted name doesn't conflict with what that
//     component was already named.
//   - On success, `nid2compo` is rewritten in-place so its values match the
//     *biologically sorted* order of the returned component name vector
//     (component 0 is whichever contig sorts first per compare_compo_names).

TEST_F(GbwtGraphTestFixture, TwoComponentsBoundAndBiologicallyReordered) {
    cfg::ARRAY_SIZE = 5; // indices [0,4]
    cfg::N_COMPO = 2;
    // Component 0 (raw) = nodes {1,2}, named "chr2".
    // Component 1 (raw) = nodes {3,4}, named "chr1".
    const auto gbz = build_gbz({1, 2, 3, 4}, {
        PathSpec{"S1", "chr2", 1, 0, fwd({1, 2})},
        PathSpec{"S2", "chr1", 1, 0, fwd({3, 4})},
    });

    std::vector<uint16_t> nid2compo(cfg::ARRAY_SIZE, cfg::NODE_UNSEEN_16);
    nid2compo[1] = 0; nid2compo[2] = 0;
    nid2compo[3] = 1; nid2compo[4] = 1;

    const auto names = bind_component_names(gbz, nid2compo, 0);

    ASSERT_EQ(names.size(), 2u);
    EXPECT_EQ(names[0], "chr1");
    EXPECT_EQ(names[1], "chr2");

    // nid2compo must be remapped to match the new (biological) ordering.
    EXPECT_EQ(nid2compo[1], 1);
    EXPECT_EQ(nid2compo[2], 1);
    EXPECT_EQ(nid2compo[3], 0);
    EXPECT_EQ(nid2compo[4], 0);
}

TEST_F(GbwtGraphTestFixture, PathCrossingComponentsThrows) {
    cfg::ARRAY_SIZE = 3;
    cfg::N_COMPO = 2;
    const auto gbz = build_gbz({1, 2}, {
        PathSpec{"S", "chr1", 1, 0, fwd({1, 2})}, // single path touching both components
    });

    std::vector<uint16_t> nid2compo(cfg::ARRAY_SIZE, cfg::NODE_UNSEEN_16);
    nid2compo[1] = 0;
    nid2compo[2] = 1; // different component than node 1, same path

    EXPECT_THROW(bind_component_names(gbz, nid2compo, 0), std::runtime_error);
}

TEST_F(GbwtGraphTestFixture, ComponentWithoutAnyPathThrows) {
    cfg::ARRAY_SIZE = 3;
    cfg::N_COMPO = 2; // component 1 will never be covered by a path
    const auto gbz = build_gbz({1}, {
        PathSpec{"S", "chr1", 1, 0, fwd({1})},
    });

    std::vector<uint16_t> nid2compo(cfg::ARRAY_SIZE, cfg::NODE_UNSEEN_16);
    nid2compo[1] = 0;
    // No node is ever assigned component 1.

    EXPECT_THROW(bind_component_names(gbz, nid2compo, 0), std::runtime_error);
}

TEST_F(GbwtGraphTestFixture, ConflictingNamesForSameComponentThrows) {
    cfg::ARRAY_SIZE = 3;
    cfg::N_COMPO = 1;
    // Two different paths, both fully contained in component 0, disagree on the contig name.
    const auto gbz = build_gbz({1, 2}, {
        PathSpec{"S1", "chr1", 1, 0, fwd({1, 2})},
        PathSpec{"S2", "chr5", 1, 0, fwd({1, 2})},
    });

    std::vector<uint16_t> nid2compo(cfg::ARRAY_SIZE, cfg::NODE_UNSEEN_16);
    nid2compo[1] = 0;
    nid2compo[2] = 0;

    EXPECT_THROW(bind_component_names(gbz, nid2compo, 0), std::runtime_error);
}

TEST_F(GbwtGraphTestFixture, MaxStepsToCheckCanSilentlyMissAComponentCrossing) {
    // Documents a real footgun rather than a "should" behavior: capping
    // validation at 1 step means a path that crosses from component 0 into
    // component 1 is never actually inspected past its first step, so the
    // crossing goes undetected as long as the (mis-)attributed name doesn't
    // conflict with what component 0 was already named.
    cfg::ARRAY_SIZE = 4;
    cfg::N_COMPO = 2;
    const auto gbz = build_gbz({1, 2, 3}, {
        PathSpec{"S1", "chr1", 1, 0, fwd({1, 2})},    // properly names component 0
        PathSpec{"S2", "chr2", 1, 0, fwd({3})},       // properly names component 1
        PathSpec{"Sbad", "chr1", 1, 0, fwd({1, 3})},  // truly crosses 0 -> 1, but agrees with component 0's name
    });

    std::vector<uint16_t> nid2compo(cfg::ARRAY_SIZE, cfg::NODE_UNSEEN_16);
    nid2compo[1] = 0; nid2compo[2] = 0;
    nid2compo[3] = 1;

    EXPECT_NO_THROW({
        const auto names = bind_component_names(gbz, nid2compo, /*max_steps_to_check=*/1);
        ASSERT_EQ(names.size(), 2u);
    });
}

// ============================================================================
// 8. COMPTAGE D'HAPLOTYPES UNIQUES (count_haplotypes)
// ============================================================================
//
// Invariants / contract:
//   - Uniqueness key is (sample, haplotype) ONLY -- the contig/locus is
//     irrelevant, so the same haplotype traversing several contigs is still
//     counted once.
//   - Paths whose name does not contain '#' fall back to being counted as
//     their own single haplotype (with a warning printed to stderr). This
//     fallback branch is not exercised here: constructing a graph where
//     `get_path_name()` yields a '#'-free string requires driving
//     GBWTGraph's internal path "sense" detection to GENERIC via an
//     internal sentinel sample name, which could not be verified against a
//     real build in this environment. Worth a follow-up test once the suite
//     can be compiled and iterated on locally.
//   - No paths at all returns 0.

TEST_F(GbwtGraphTestFixture, NoPathsReturnsZero) {
    const auto gbz = build_gbz({1}, {});
    EXPECT_EQ(count_haplotypes(gbz.graph), 0u);
}

TEST_F(GbwtGraphTestFixture, DistinctSamplesAndHaplotypesAreCountedSeparately) {
    const auto gbz = build_gbz({1, 2}, {
        PathSpec{"HG001", "chr1", 1, 0, fwd({1, 2})},
        PathSpec{"HG001", "chr1", 2, 0, fwd({1, 2})},
        PathSpec{"HG002", "chr1", 1, 0, fwd({1, 2})},
    });

    EXPECT_EQ(count_haplotypes(gbz.graph), 3u);
}

TEST_F(GbwtGraphTestFixture, SameHaplotypeAcrossDifferentContigsCountedOnce) {
    // Same (sample, haplotype) pair traversing two different contigs must
    // still count as a single haplotype: the contig is not part of the key.
    const auto gbz = build_gbz({1, 2, 3, 4}, {
        PathSpec{"HG001", "chr1", 1, 0, fwd({1, 2})},
        PathSpec{"HG001", "chr2", 1, 0, fwd({3, 4})},
    });

    EXPECT_EQ(count_haplotypes(gbz.graph), 1u);
}

TEST_F(GbwtGraphTestFixture, DuplicatePathsCollapseToOneHaplotype) {
    const auto gbz = build_gbz({1, 2}, {
        PathSpec{"HG001", "chr1", 1, 0, fwd({1, 2})},
        PathSpec{"HG001", "chr1", 1, 0, fwd({1, 2})}, // exact duplicate (sample, haplotype)
    });

    EXPECT_EQ(count_haplotypes(gbz.graph), 1u);
}
