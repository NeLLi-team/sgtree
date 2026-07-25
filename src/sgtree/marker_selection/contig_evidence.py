from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Mapping


MIN_INFORMATIVE_GENES = 3
MIN_AGREEMENT_NUMERATOR = 2
MIN_AGREEMENT_DENOMINATOR = 3
MIN_CONFLICTING_GENES = 2


def _clean_id(value: object) -> str | None:
    if not isinstance(value, str):
        return None
    cleaned = value.strip()
    return cleaned or None


def contig_gene_vote_gate(
    votes: Iterable[Mapping[str, object]],
    proposed_clade: str | None,
) -> dict:
    """Test whether non-marker genes support the marker's attachment clade.

    Each informative vote must provide a distinct ``gene_id`` and a non-empty
    ``assigned_clade``. Two genes assigned to the same competing clade veto the
    proposed clade. The function only scores evidence; it does not remove a marker.
    """
    proposed = _clean_id(proposed_clade)
    informative: list[tuple[str, str]] = []
    seen_gene_ids: set[str] = set()
    invalid_gene_vote_count = 0
    duplicate_gene_vote_count = 0
    for vote in votes:
        if vote.get("informative") is not True:
            continue
        gene_id = _clean_id(vote.get("gene_id"))
        assigned_clade = _clean_id(vote.get("assigned_clade"))
        if gene_id is None or assigned_clade is None:
            invalid_gene_vote_count += 1
            continue
        if gene_id in seen_gene_ids:
            duplicate_gene_vote_count += 1
            continue
        seen_gene_ids.add(gene_id)
        informative.append((gene_id, assigned_clade))

    clades = [assigned_clade for _gene_id, assigned_clade in informative]
    counts = Counter(clades)
    top_clade = (
        min(counts, key=lambda clade: (-counts[clade], clade))
        if counts
        else None
    )
    agreement_count = counts.get(proposed, 0) if proposed is not None else 0
    informative_count = len(informative)
    agreement_fraction = (
        agreement_count / informative_count
        if informative_count
        else 0.0
    )
    conflicting_counts = {
        clade: count
        for clade, count in counts.items()
        if clade != proposed
    }
    strongest_conflict_clade = (
        min(
            conflicting_counts,
            key=lambda clade: (-conflicting_counts[clade], clade),
        )
        if conflicting_counts
        else None
    )
    strongest_conflict_count = (
        conflicting_counts[strongest_conflict_clade]
        if strongest_conflict_clade is not None
        else 0
    )
    strongest_conflict_fraction = (
        strongest_conflict_count / informative_count
        if informative_count
        else 0.0
    )
    strong_conflict_count = (
        strongest_conflict_count
        if strongest_conflict_count >= MIN_CONFLICTING_GENES
        else 0
    )

    reason = None
    passed = False
    if not proposed:
        reason = "missing_proposed_clade"
    elif invalid_gene_vote_count:
        reason = "invalid_gene_votes"
    elif duplicate_gene_vote_count:
        reason = "duplicate_gene_votes"
    elif informative_count < MIN_INFORMATIVE_GENES:
        reason = "insufficient_informative_genes"
    elif strong_conflict_count:
        reason = "strong_conflicting_clade"
    elif (
        agreement_count * MIN_AGREEMENT_DENOMINATOR
        < informative_count * MIN_AGREEMENT_NUMERATOR
    ):
        reason = "insufficient_donor_agreement"
    else:
        passed = True

    return {
        "proposed_clade": proposed,
        "top_clade": top_clade,
        "informative_gene_count": informative_count,
        "invalid_gene_vote_count": invalid_gene_vote_count,
        "duplicate_gene_vote_count": duplicate_gene_vote_count,
        "agreement_count": agreement_count,
        "agreement_fraction": agreement_fraction,
        "strongest_conflict_clade": strongest_conflict_clade,
        "strongest_conflict_count": strongest_conflict_count,
        "strongest_conflict_fraction": strongest_conflict_fraction,
        "strong_conflict_count": strong_conflict_count,
        "contig_gate_pass": passed,
        "contig_abstention_reason": reason,
    }
