from __future__ import annotations

import math
import statistics
from collections.abc import Mapping
from itertools import combinations

from ete3 import Tree

from sgtree.id_schema import parse_sequence_id


MIN_VOTERS = 5
MIN_COORDINATES = 6
MIN_TARGET_SUPPORT = 0.70
MIN_TARGET_DISCORDANCE = 0.02
MIN_MARKER_MARGIN = 0.02
MAX_EXACT_VOTERS = 12
REVIEW_ROBUST_Z = 3.0
TOLERANCE = 1e-12


def _genome_id(leaf_name: str) -> str:
    return leaf_name.split("|", 1)[0]


def _canonical_split(
    left: set[str],
    right: set[str],
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    sides = [tuple(sorted(left)), tuple(sorted(right))]
    sides.sort(key=lambda side: (len(side), side))
    return sides[0], sides[1]


def _explicit_split_supports(tree: Tree) -> dict[
    tuple[tuple[str, ...], tuple[str, ...]],
    float | None,
]:
    all_leaves = {str(leaf.name) for leaf in tree.iter_leaves()}
    observations: dict[
        tuple[tuple[str, ...], tuple[str, ...]],
        list[float | None],
    ] = {}
    invalid: set[tuple[tuple[str, ...], tuple[str, ...]]] = set()

    for node in tree.traverse():
        if node.is_root() or node.is_leaf():
            continue
        left = {str(leaf.name) for leaf in node.iter_leaves()}
        right = all_leaves - left
        if len(left) < 2 or len(right) < 2:
            continue
        key = _canonical_split(left, right)
        label = str(node.name).strip()
        if not label:
            observations.setdefault(key, []).append(None)
            continue
        try:
            value = float(label)
        except ValueError:
            invalid.add(key)
            continue
        if not math.isfinite(value) or not 0.0 <= value <= 1.0:
            invalid.add(key)
            continue
        observations.setdefault(key, []).append(value)

    merged: dict[tuple[tuple[str, ...], tuple[str, ...]], float | None] = {}
    for key, values in observations.items():
        numeric = [value for value in values if value is not None]
        if key in invalid or not numeric:
            merged[key] = None
        elif max(numeric) - min(numeric) > TOLERANCE:
            merged[key] = None
        else:
            merged[key] = min(numeric)
    for key in invalid:
        merged.setdefault(key, None)
    return merged


def _prepare_tree(tree: Tree) -> dict:
    leaves = list(tree.iter_leaves())
    leaf_names = [str(leaf.name) for leaf in leaves]
    genome_to_leaves: dict[str, list[str]] = {}
    for leaf_name in leaf_names:
        genome_to_leaves.setdefault(_genome_id(leaf_name), []).append(leaf_name)
    for names in genome_to_leaves.values():
        names.sort()

    valid = bool(leaf_names) and len(leaf_names) == len(set(leaf_names))
    for node in tree.traverse():
        if node.is_root():
            continue
        distance = float(node.dist)
        if not math.isfinite(distance) or distance < 0.0:
            valid = False
            break

    duplicate_genomes = {
        genome for genome, names in genome_to_leaves.items() if len(names) != 1
    }
    leaf_by_genome = {
        genome: names[0]
        for genome, names in genome_to_leaves.items()
        if len(names) == 1
    }
    distances: dict[tuple[str, str], float] = {}
    if valid and not duplicate_genomes:
        node_by_genome = {
            _genome_id(str(leaf.name)): leaf
            for leaf in leaves
        }
        for left, right in combinations(sorted(node_by_genome), 2):
            distance = float(node_by_genome[left].get_distance(node_by_genome[right]))
            if not math.isfinite(distance) or distance < 0.0:
                valid = False
                distances = {}
                break
            distances[(left, right)] = distance

    return {
        "leaf_names": sorted(leaf_names),
        "genome_to_leaves": genome_to_leaves,
        "leaf_by_genome": leaf_by_genome,
        "taxa": set(genome_to_leaves),
        "duplicate_genomes": duplicate_genomes,
        "valid": valid,
        "distances": distances,
        "split_supports": _explicit_split_supports(tree) if valid else {},
    }


def _pair_distance(prepared: dict, left: str, right: str) -> float | None:
    if left == right:
        return 0.0
    key = tuple(sorted((left, right)))
    return prepared["distances"].get(key)


def _target_attachment_evidence(
    prepared: dict,
    leaf_name: str,
) -> tuple[float | None, list[dict], list[str]]:
    supports = prepared["split_supports"]
    target_sides: list[
        tuple[frozenset[str], tuple[tuple[str, ...], tuple[str, ...]]]
    ] = []
    for key in supports:
        left = frozenset(key[0])
        right = frozenset(key[1])
        if leaf_name in left:
            target_sides.append((left, key))
        elif leaf_name in right:
            target_sides.append((right, key))

    minimal = [
        (side, key)
        for side, key in target_sides
        if not any(other < side for other, _other_key in target_sides)
    ]
    if not minimal:
        return None, [], []

    distinct_sides = {tuple(sorted(side)) for side, _key in minimal}
    support_sides = []
    for side in sorted(distinct_sides, key=lambda item: (len(item), item)):
        key = next(key for candidate, key in minimal if tuple(sorted(candidate)) == side)
        support_sides.append(
            {
                "taxa": sorted(_genome_id(name) for name in side),
                "support": supports[key],
            }
        )
    attachment_taxa = (
        sorted(
            _genome_id(name)
            for name in next(iter(distinct_sides))
            if name != leaf_name
        )
        if len(distinct_sides) == 1
        else []
    )
    values = [supports[key] for _side, key in minimal]
    if any(value is None for value in values):
        return None, support_sides, attachment_taxa
    return (
        min(float(value) for value in values if value is not None),
        support_sides,
        attachment_taxa,
    )


def _target_attachment_support(prepared: dict, leaf_name: str) -> float | None:
    support, _support_sides, _attachment_taxa = _target_attachment_evidence(
        prepared,
        leaf_name,
    )
    return support


def _normalized_profile(
    prepared: dict,
    target_genome: str,
    coordinates: list[str],
) -> dict[str, float] | None:
    distances: dict[str, float] = {}
    for coordinate in coordinates:
        distance = _pair_distance(prepared, target_genome, coordinate)
        if distance is None or not math.isfinite(distance) or distance < 0.0:
            return None
        distances[coordinate] = distance
    total = sum(distances.values())
    if not math.isfinite(total) or total <= 0.0:
        return None
    return {
        coordinate: distances[coordinate] / total
        for coordinate in coordinates
    }


def _consensus(profiles: list[dict[str, float]], coordinates: list[str]) -> dict[str, float] | None:
    values = {
        coordinate: statistics.median(profile[coordinate] for profile in profiles)
        for coordinate in coordinates
    }
    total = sum(values.values())
    if not math.isfinite(total) or total <= 0.0:
        return None
    return {
        coordinate: values[coordinate] / total
        for coordinate in coordinates
    }


def _total_variation(
    left: dict[str, float],
    right: dict[str, float],
    coordinates: list[str],
) -> float:
    return 0.5 * sum(abs(left[item] - right[item]) for item in coordinates)


def _base_row(marker_name: str, leaf_name: str) -> dict:
    parsed_genome, contig_id, gene_id = parse_sequence_id(leaf_name)
    genome = _genome_id(leaf_name) or parsed_genome
    return {
        "marker_name": marker_name,
        "leaf_name": leaf_name,
        "genome": genome,
        "contig_id": contig_id,
        "gene_id": gene_id,
        "loo_class": "ambiguous",
        "loo_abstention_reason": None,
        "loo_target_support": None,
        "loo_target_support_sides": [],
        "loo_attachment_taxa": [],
        "loo_attachment_clade": None,
        "loo_voter_count": 0,
        "loo_voter_markers": [],
        "loo_voter_search_mode": None,
        "loo_coordinate_count": 0,
        "loo_coordinate_taxa": [],
        "loo_target_discordance": None,
        "loo_voter_center": None,
        "loo_voter_mad": None,
        "loo_voter_upper": None,
        "loo_conflict_margin": None,
        "loo_robust_z": None,
        "loo_score": None,
        "loo_conflict_beyond_dispersion": None,
        "loo_review_candidate": False,
        "loo_marker_rank": None,
        "loo_marker_margin": None,
        "loo_decision": "kept_report_only",
    }


def _shared_coordinates(target: dict, voters: list[tuple[str, dict]], genome: str) -> list[str]:
    coordinates = set(target["taxa"])
    for _voter_name, voter in voters:
        coordinates &= voter["taxa"]
    coordinates.discard(genome)
    return sorted(coordinates)


def _select_voters(
    target: dict,
    voters: list[tuple[str, dict]],
    genome: str,
) -> tuple[list[tuple[str, dict]], list[str], str]:
    """Find a shared-coordinate voter set without unbounded subset search."""
    if len(voters) <= MAX_EXACT_VOTERS:
        for size in range(len(voters), MIN_VOTERS - 1, -1):
            candidates = []
            for subset in combinations(voters, size):
                selected = list(subset)
                shared = _shared_coordinates(target, selected, genome)
                if len(shared) >= MIN_COORDINATES:
                    candidates.append(
                        (
                            -len(shared),
                            tuple(name for name, _voter in selected),
                            selected,
                            shared,
                        )
                    )
            if candidates:
                _negative_count, _names, selected, shared = min(candidates)
                return selected, shared, "exact"

    # Large production panels stay bounded. Greedy failure only causes a
    # report-only abstention; the fixed eight-marker benchmark uses exact search.
    selected = list(voters)
    coordinates = _shared_coordinates(target, selected, genome)
    while len(coordinates) < MIN_COORDINATES and len(selected) > MIN_VOTERS:
        choices = []
        for index, (voter_name, _voter) in enumerate(selected):
            remaining = selected[:index] + selected[index + 1 :]
            shared = _shared_coordinates(target, remaining, genome)
            choices.append((-len(shared), voter_name, index, shared))
        _negative_count, _voter_name, remove_index, coordinates = min(choices)
        selected.pop(remove_index)
    return selected, coordinates, "greedy"


def _apply_within_marker_gate(rows: list[dict]) -> None:
    by_marker: dict[str, list[dict]] = {}
    for row in rows:
        if (
            row["loo_class"] == "discordant_marker"
            and row["loo_abstention_reason"] is None
            and row["loo_target_discordance"] is not None
        ):
            by_marker.setdefault(str(row["marker_name"]), []).append(row)

    for marker_rows in by_marker.values():
        ranked = sorted(
            marker_rows,
            key=lambda row: (
                -float(row["loo_target_discordance"]),
                str(row["leaf_name"]),
            ),
        )
        for rank, row in enumerate(ranked, start=1):
            row["loo_marker_rank"] = rank
        top = ranked[0]
        runner_up_score = (
            float(ranked[1]["loo_target_discordance"])
            if len(ranked) > 1
            else 0.0
        )
        top["loo_marker_margin"] = (
            float(top["loo_target_discordance"]) - runner_up_score
        )

        for row in ranked:
            score = float(row["loo_target_discordance"])
            if score < MIN_TARGET_DISCORDANCE:
                row["loo_class"] = "ambiguous"
                row["loo_abstention_reason"] = "discordance_below_effect_floor"
            elif int(row["loo_marker_rank"]) != 1:
                row["loo_class"] = "ambiguous"
                row["loo_abstention_reason"] = "not_top_marker_outlier"
            elif float(row["loo_marker_margin"]) < MIN_MARKER_MARGIN:
                row["loo_class"] = "ambiguous"
                row["loo_abstention_reason"] = "marker_rank_not_unique"


def score_loo_profiles(
    marker_trees: Mapping[str, Tree],
    excluded_target_genomes: set[str] | None = None,
) -> list[dict]:
    """Score every marker leaf against the target genome's other markers.

    Trees must be parsed with ETE ``format=1`` so explicit internal support
    labels remain distinguishable from missing labels.
    """
    excluded_target_genomes = excluded_target_genomes or set()
    prepared = {
        str(marker_name): _prepare_tree(tree)
        for marker_name, tree in sorted(marker_trees.items(), key=lambda item: str(item[0]))
    }
    rows: list[dict] = []

    for marker_name in sorted(prepared):
        target = prepared[marker_name]
        for leaf_name in target["leaf_names"]:
            row = _base_row(marker_name, leaf_name)
            genome = row["genome"]
            if genome in excluded_target_genomes:
                row["loo_abstention_reason"] = "reference_target"
                rows.append(row)
                continue
            if genome in target["duplicate_genomes"]:
                row["loo_abstention_reason"] = "target_not_single_copy"
                rows.append(row)
                continue
            if not target["valid"] or target["duplicate_genomes"]:
                row["loo_abstention_reason"] = "invalid_target_tree"
                rows.append(row)
                continue

            support, support_sides, attachment_taxa = _target_attachment_evidence(
                target,
                leaf_name,
            )
            row["loo_target_support"] = support
            row["loo_target_support_sides"] = support_sides
            row["loo_attachment_taxa"] = attachment_taxa
            row["loo_attachment_clade"] = ",".join(attachment_taxa) or None
            voters: list[tuple[str, dict]] = []
            for voter_name in sorted(prepared):
                if voter_name == marker_name:
                    continue
                voter = prepared[voter_name]
                if not voter["valid"] or voter["duplicate_genomes"]:
                    continue
                if genome not in voter["leaf_by_genome"]:
                    continue
                shared = (target["taxa"] & voter["taxa"]) - {genome}
                if len(shared) >= MIN_COORDINATES:
                    voters.append((voter_name, voter))

            if len(voters) < MIN_VOTERS:
                row["loo_voter_count"] = len(voters)
                row["loo_voter_markers"] = [name for name, _voter in voters]
                row["loo_abstention_reason"] = "insufficient_voters"
                rows.append(row)
                continue

            voters, coordinates, search_mode = _select_voters(target, voters, genome)
            row["loo_voter_count"] = len(voters)
            row["loo_voter_markers"] = [name for name, _voter in voters]
            row["loo_voter_search_mode"] = search_mode
            row["loo_coordinate_count"] = len(coordinates)
            row["loo_coordinate_taxa"] = coordinates
            if len(coordinates) < MIN_COORDINATES:
                row["loo_abstention_reason"] = "insufficient_coordinates"
                rows.append(row)
                continue

            target_profile = _normalized_profile(target, genome, coordinates)
            voter_profiles = [
                _normalized_profile(voter, genome, coordinates)
                for _voter_name, voter in voters
            ]
            if target_profile is None or any(profile is None for profile in voter_profiles):
                row["loo_abstention_reason"] = "invalid_patristic_profile"
                rows.append(row)
                continue
            complete_voter_profiles = [
                profile for profile in voter_profiles if profile is not None
            ]
            voter_consensus = _consensus(complete_voter_profiles, coordinates)
            if voter_consensus is None:
                row["loo_abstention_reason"] = "invalid_patristic_profile"
                rows.append(row)
                continue

            target_discordance = _total_variation(
                target_profile,
                voter_consensus,
                coordinates,
            )
            voter_discordances: list[float] = []
            for index, voter_profile in enumerate(complete_voter_profiles):
                leave_one_out = complete_voter_profiles[:index] + complete_voter_profiles[index + 1 :]
                consensus = _consensus(leave_one_out, coordinates)
                if consensus is None:
                    voter_discordances = []
                    break
                voter_discordances.append(
                    _total_variation(voter_profile, consensus, coordinates)
                )
            if not voter_discordances:
                row["loo_abstention_reason"] = "invalid_patristic_profile"
                rows.append(row)
                continue

            center = statistics.median(voter_discordances)
            mad = statistics.median(abs(value - center) for value in voter_discordances)
            upper = max(max(voter_discordances), center + (3.0 * 1.4826 * mad))
            margin = target_discordance - upper
            beyond_dispersion = target_discordance > upper + TOLERANCE
            robust_z = (
                (target_discordance - center) / (1.4826 * mad)
                if mad > 0
                else None
            )
            row.update(
                {
                    "loo_target_discordance": target_discordance,
                    "loo_voter_center": center,
                    "loo_voter_mad": mad,
                    "loo_voter_upper": upper,
                    "loo_conflict_margin": margin,
                    "loo_robust_z": robust_z,
                    "loo_score": target_discordance,
                    "loo_conflict_beyond_dispersion": beyond_dispersion,
                }
            )

            support = row["loo_target_support"]
            if support is None:
                row["loo_abstention_reason"] = "missing_target_support"
            elif support < MIN_TARGET_SUPPORT:
                row["loo_abstention_reason"] = "target_support_below_threshold"
            elif target_discordance <= center + TOLERANCE:
                row["loo_class"] = "clean"
            elif beyond_dispersion:
                row["loo_class"] = "discordant_marker"
            else:
                row["loo_abstention_reason"] = "within_voter_dispersion"
            rows.append(row)

    _apply_within_marker_gate(rows)
    for row in rows:
        # Report-only review tier: dispersion-ceiling abstentions whose conflict
        # is still robustly outside the voter MAD band. Never a call by itself;
        # confirmation requires downstream contig evidence.
        row["loo_review_candidate"] = (
            row["loo_abstention_reason"] == "within_voter_dispersion"
            and row["loo_robust_z"] is not None
            and float(row["loo_robust_z"]) >= REVIEW_ROBUST_Z
        )
    return rows
