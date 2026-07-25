import unittest
import warnings

warnings.simplefilter("ignore", SyntaxWarning)

from ete3 import Tree

from sgtree.marker_selection import loo_profile


TAXA = ("G", "A", "B", "C", "D", "E", "F", "H")


def _leaf(genome: str, marker: str, *, suffix: str = "") -> str:
    return f"{genome}|ctg|{marker}_{genome}{suffix}:1"


def _base_tree(
    marker: str,
    *,
    target_support: str | None = "0.95",
    far_stem: float = 1.0,
) -> Tree:
    support = "" if target_support is None else target_support
    newick = (
        f"((({_leaf('G', marker)},{_leaf('A', marker)}){support}:1,"
        f"({_leaf('B', marker)},{_leaf('C', marker)})0.95:1)0.95:1,"
        f"(({_leaf('D', marker)},{_leaf('E', marker)})0.95:1,"
        f"({_leaf('F', marker)},{_leaf('H', marker)})0.95:1)0.95:{far_stem});"
    )
    return Tree(newick, format=1)


def _localized_regraft_tree(marker: str) -> Tree:
    def name(genome: str) -> str:
        return f"{genome}|ctg|{marker}_{genome}"

    return Tree(
        (
            f"(({name('A')}:2,"
            f"({name('B')}:1,{name('C')}:1)0.95:1)0.95:1,"
            f"(({name('D')}:1,{name('E')}:1)0.95:1,"
            f"({name('H')}:1,"
            f"({name('F')}:0.5,{name('G')}:1)0.95:0.5)0.95:1)0.95:1);"
        ),
        format=1,
    )


def _scaled(tree: Tree, factor: float) -> Tree:
    result = tree.copy()
    for node in result.traverse():
        if not node.is_root():
            node.dist *= factor
    return result


def _without_genome(tree: Tree, genome: str) -> Tree:
    result = tree.copy()
    keep = [
        leaf.name
        for leaf in result.iter_leaves()
        if leaf.name.split("|", 1)[0] != genome
    ]
    result.prune(keep)
    return result


def _row(rows: list[dict], marker: str, genome: str = "G") -> dict:
    return next(
        row
        for row in rows
        if row["marker_name"] == marker and row["genome"] == genome
    )


def _nearest_genomes(tree: Tree, genome: str, count: int) -> list[str]:
    target = next(leaf for leaf in tree.iter_leaves() if leaf.name.split("|", 1)[0] == genome)
    neighbors = sorted(
        (
            float(target.get_distance(other)),
            other.name.split("|", 1)[0],
        )
        for other in tree.iter_leaves()
        if other is not target
    )
    return [name for _distance, name in neighbors[:count]]


class LeaveOneMarkerOutProfileTests(unittest.TestCase):
    def test_concordant_scaled_profiles_are_clean(self):
        trees = {
            f"M{index}": _scaled(_base_tree(f"M{index}"), index + 1)
            for index in range(6)
        }

        result = _row(loo_profile.score_loo_profiles(trees), "M0")

        self.assertEqual(result["loo_class"], "clean")
        self.assertIsNone(result["loo_abstention_reason"])
        self.assertEqual(result["loo_target_support"], 0.95)
        self.assertEqual(result["loo_attachment_taxa"], ["A"])
        self.assertEqual(result["loo_attachment_clade"], "A")
        self.assertEqual(result["loo_voter_count"], 5)
        self.assertEqual(result["loo_coordinate_count"], 7)
        self.assertEqual(result["loo_decision"], "kept_report_only")

    def test_strong_placement_conflict_is_discordant(self):
        trees = {"M0": _localized_regraft_tree("M0")}
        trees.update({f"M{index}": _base_tree(f"M{index}") for index in range(1, 6)})

        target_rows = [
            row
            for row in loo_profile.score_loo_profiles(trees)
            if row["marker_name"] == "M0"
        ]
        result = _row(target_rows, "M0")

        self.assertEqual(result["loo_class"], "discordant_marker")
        self.assertTrue(result["loo_conflict_beyond_dispersion"])
        self.assertGreater(result["loo_conflict_margin"], 0.0)
        self.assertEqual(result["loo_marker_rank"], 1)
        self.assertEqual(result["loo_attachment_taxa"], ["F"])
        self.assertEqual(result["loo_attachment_clade"], "F")
        self.assertGreaterEqual(result["loo_marker_margin"], loo_profile.MIN_MARKER_MARGIN)
        self.assertEqual(
            {row["genome"] for row in target_rows if row["loo_class"] == "discordant_marker"},
            {"G"},
        )

    def test_full_profile_detects_conflict_with_same_nearest_neighbors(self):
        target = _base_tree("M0", far_stem=8.0)
        trees = {"M0": target}
        trees.update({f"M{index}": _base_tree(f"M{index}") for index in range(1, 6)})

        self.assertEqual(
            _nearest_genomes(target, "G", 3),
            _nearest_genomes(trees["M1"], "G", 3),
        )
        target_rows = [
            row
            for row in loo_profile.score_loo_profiles(trees)
            if row["marker_name"] == "M0"
        ]
        result = _row(target_rows, "M0")

        self.assertGreater(result["loo_target_discordance"], 0.0)
        self.assertFalse(any(row["loo_class"] == "discordant_marker" for row in target_rows))
        self.assertIn(
            "marker_rank_not_unique",
            {row["loo_abstention_reason"] for row in target_rows},
        )

    def test_tiny_clean_branch_noise_does_not_create_discordant_marker(self):
        trees = {"M0": _base_tree("M0", far_stem=1.000000001)}
        trees.update({f"M{index}": _base_tree(f"M{index}") for index in range(1, 6)})

        target_rows = [
            row
            for row in loo_profile.score_loo_profiles(trees)
            if row["marker_name"] == "M0"
        ]

        self.assertFalse(any(row["loo_class"] == "discordant_marker" for row in target_rows))
        self.assertTrue(
            all(
                row["loo_abstention_reason"] == "discordance_below_effect_floor"
                for row in target_rows
            )
        )

    def test_ineligible_higher_score_does_not_mask_supported_conflict(self):
        for reason in (
            "missing_target_support",
            "target_support_below_threshold",
            "within_voter_dispersion",
        ):
            with self.subTest(reason=reason):
                valid = {
                    "marker_name": "M0",
                    "leaf_name": "G|c|g",
                    "loo_class": "discordant_marker",
                    "loo_abstention_reason": None,
                    "loo_target_discordance": 0.20,
                    "loo_marker_rank": None,
                    "loo_marker_margin": None,
                }
                ineligible = {
                    "marker_name": "M0",
                    "leaf_name": "H|c|h",
                    "loo_class": "ambiguous",
                    "loo_abstention_reason": reason,
                    "loo_target_discordance": 0.30,
                    "loo_marker_rank": None,
                    "loo_marker_margin": None,
                }

                loo_profile._apply_within_marker_gate([ineligible, valid])

                self.assertEqual(valid["loo_class"], "discordant_marker")
                self.assertEqual(valid["loo_marker_rank"], 1)
                self.assertEqual(valid["loo_marker_margin"], 0.20)
                self.assertIsNone(ineligible["loo_marker_rank"])

    def test_four_voters_abstains(self):
        trees = {f"M{index}": _base_tree(f"M{index}") for index in range(5)}

        result = _row(loo_profile.score_loo_profiles(trees), "M0")

        self.assertEqual(result["loo_abstention_reason"], "insufficient_voters")
        self.assertEqual(result["loo_voter_count"], 4)

    def test_five_coordinate_taxa_abstains(self):
        trees = {f"M{index}": _base_tree(f"M{index}") for index in range(6)}
        trees["M1"] = _without_genome(trees["M1"], "F")
        trees["M2"] = _without_genome(trees["M2"], "H")

        result = _row(loo_profile.score_loo_profiles(trees), "M0")

        self.assertEqual(result["loo_abstention_reason"], "insufficient_coordinates")
        self.assertEqual(result["loo_voter_count"], 5)
        self.assertEqual(result["loo_coordinate_count"], 5)

    def test_ragged_extra_voters_do_not_erase_valid_consensus(self):
        trees = {f"M{index}": _base_tree(f"M{index}") for index in range(8)}
        trees["M6"] = _without_genome(trees["M6"], "F")
        trees["M7"] = _without_genome(trees["M7"], "H")

        result = _row(loo_profile.score_loo_profiles(trees), "M0")

        self.assertNotEqual(result["loo_abstention_reason"], "insufficient_coordinates")
        self.assertEqual(result["loo_voter_count"], 6)
        self.assertEqual(result["loo_coordinate_count"], 6)

    def test_small_panel_subset_search_finds_non_greedy_consensus(self):
        coordinate_sets = (
            "0245678",
            "0245678",
            "0234678",
            "0123678",
            "01235678",
            "025678",
            "0245678",
        )
        target = {"taxa": {"G", *"012345678"}}
        voters = [
            (f"M{index}", {"taxa": {"G", *coordinates}})
            for index, coordinates in enumerate(coordinate_sets)
        ]

        selected, coordinates = loo_profile._select_voters(target, voters, "G")

        self.assertGreaterEqual(len(selected), loo_profile.MIN_VOTERS)
        self.assertGreaterEqual(len(coordinates), loo_profile.MIN_COORDINATES)
        self.assertEqual(
            set(coordinates),
            set.intersection(
                target["taxa"],
                *(voter["taxa"] for _name, voter in selected),
            )
            - {"G"},
        )

    def test_missing_and_low_target_support_abstain(self):
        voters = {f"M{index}": _base_tree(f"M{index}") for index in range(1, 6)}
        cases = (
            (None, "missing_target_support"),
            ("0.69", "target_support_below_threshold"),
            ("0.70", None),
            ("95", "missing_target_support"),
        )
        for support, reason in cases:
            with self.subTest(support=support):
                trees = {"M0": _base_tree("M0", target_support=support), **voters}
                result = _row(loo_profile.score_loo_profiles(trees), "M0")
                self.assertEqual(result["loo_abstention_reason"], reason)
                if reason is None:
                    self.assertEqual(result["loo_class"], "clean")

    def test_target_support_is_stable_across_equivalent_rootings(self):
        original = Tree(
            "(((G|c|g:1,A|c|a:1)0.8:1,(B|c|b:1,C|c|c:1)0.9:1)0.85:1,"
            "((D|c|d:1,E|c|e:1)0.9:1,(F|c|f:1,H|c|h:1)0.9:1)0.85:1);",
            format=1,
        )
        rerooted = Tree(
            "((D|c|d:1,E|c|e:1)0.9:1,"
            "(((G|c|g:1,A|c|a:1)0.8:1,(B|c|b:1,C|c|c:1)0.9:1)0.85:1,"
            "(F|c|f:1,H|c|h:1)0.9:1)0.9:1);",
            format=1,
        )

        left = loo_profile._target_attachment_support(
            loo_profile._prepare_tree(original),
            "G|c|g",
        )
        right = loo_profile._target_attachment_support(
            loo_profile._prepare_tree(rerooted),
            "G|c|g",
        )

        self.assertEqual(left, 0.8)
        self.assertEqual(right, 0.8)

    def test_non_cherry_target_keeps_all_minimal_support_sides(self):
        tree = Tree(
            "((A|c|a:1,B|c|b:1)0.81:1,G|c|g:1,"
            "((C|c|c:1,D|c|d:1)0.82:1,(E|c|e:1,F|c|f:1)0.83:1)0.84:1);",
            format=1,
        )

        support, sides, attachment_taxa = loo_profile._target_attachment_evidence(
            loo_profile._prepare_tree(tree),
            "G|c|g",
        )

        self.assertEqual(support, 0.81)
        self.assertEqual([side["support"] for side in sides], [0.84, 0.81])
        self.assertEqual(
            [side["taxa"] for side in sides],
            [["A", "B", "G"], ["C", "D", "E", "F", "G"]],
        )
        self.assertEqual(attachment_taxa, [])

    def test_complementary_support_merge_is_conservative_and_order_independent(self):
        left_first = Tree(
            "((A|c|a:1,B|c|b:1)0.6999999999998:1,"
            "(G|c|g:1,C|c|c:1)0.7000000000002:1);",
            format=1,
        )
        right_first = Tree(
            "((G|c|g:1,C|c|c:1)0.7000000000002:1,"
            "(A|c|a:1,B|c|b:1)0.6999999999998:1);",
            format=1,
        )

        supports = [
            loo_profile._target_attachment_support(
                loo_profile._prepare_tree(tree),
                "G|c|g",
            )
            for tree in (left_first, right_first)
        ]

        self.assertEqual(supports, [0.6999999999998, 0.6999999999998])

    def test_voter_dispersion_can_force_abstention(self):
        trees = {
            "M0": _base_tree("M0", far_stem=3.0),
            "M1": _base_tree("M1", far_stem=1.0),
            "M2": _base_tree("M2", far_stem=1.0),
            "M3": _base_tree("M3", far_stem=1.0),
            "M4": _base_tree("M4", far_stem=5.0),
            "M5": _base_tree("M5", far_stem=5.0),
        }

        result = _row(loo_profile.score_loo_profiles(trees), "M0")

        self.assertEqual(result["loo_abstention_reason"], "within_voter_dispersion")
        self.assertFalse(result["loo_conflict_beyond_dispersion"])

    def test_duplicate_target_genome_abstains(self):
        duplicate = Tree(
            "(((G|c|g1:1,G|c|g2:1)0.95:1,(B|c|b:1,C|c|c:1)0.95:1)0.95:1,"
            "((D|c|d:1,E|c|e:1)0.95:1,(F|c|f:1,H|c|h:1)0.95:1)0.95:1);",
            format=1,
        )
        trees = {"M0": duplicate}
        trees.update({f"M{index}": _base_tree(f"M{index}") for index in range(1, 6)})

        results = [
            row
            for row in loo_profile.score_loo_profiles(trees)
            if row["marker_name"] == "M0" and row["genome"] == "G"
        ]

        self.assertEqual(len(results), 2)
        self.assertEqual(
            {row["loo_abstention_reason"] for row in results},
            {"target_not_single_copy"},
        )

    def test_excluded_reference_targets_are_reported_but_not_scored(self):
        trees = {f"M{index}": _base_tree(f"M{index}") for index in range(6)}

        results = loo_profile.score_loo_profiles(
            trees,
            excluded_target_genomes={"G"},
        )
        reference_rows = [row for row in results if row["genome"] == "G"]

        self.assertEqual(len(reference_rows), 6)
        self.assertEqual(
            {row["loo_abstention_reason"] for row in reference_rows},
            {"reference_target"},
        )
        self.assertTrue(all(row["loo_decision"] == "kept_report_only" for row in reference_rows))

    def test_output_order_and_fields_are_deterministic(self):
        ascending = {f"M{index}": _base_tree(f"M{index}") for index in range(6)}
        descending = dict(reversed(list(ascending.items())))

        first = loo_profile.score_loo_profiles(ascending)
        second = loo_profile.score_loo_profiles(descending)

        self.assertEqual(first, second)
        self.assertEqual(
            [(row["marker_name"], row["leaf_name"]) for row in first],
            sorted((row["marker_name"], row["leaf_name"]) for row in first),
        )


if __name__ == "__main__":
    unittest.main()
