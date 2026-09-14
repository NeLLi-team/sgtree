"""Benchmark-only contig evidence from protein-sequence matches."""

from __future__ import annotations

import hashlib
from collections.abc import Collection, Mapping
from dataclasses import dataclass
from typing import TypeAlias, TypedDict

from pyhmmer import easel, hmmer, plan7

from sgtree.id_schema import parse_sequence_id

AMINO_ACID_CHARACTERS = frozenset("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
PHMMER_MIN_QUERY_COUNT = 3
PHMMER_MAX_QUERY_COUNT = 12
PHMMER_MAX_EVALUE = 1e-5
PHMMER_MIN_QUERY_COVERAGE = 0.5
PHMMER_MIN_SCORE_MARGIN = 5.0
PHMMER_REPORT_EVALUE = 10.0
PHMMER_BUILDER_SEED = 42

_QueryItem: TypeAlias = tuple[str, object]


class _HitAudit(TypedDict):
    raw_hit_count: int
    evalue_filtered_hit_count: int
    coverage_filtered_hit_count: int
    qualifying_hit_record_count: int
    minimum_observed_evalue: float | None
    maximum_observed_query_coverage: float | None
    attachment_qualifying_genome_count: int
    background_qualifying_genome_count: int
    attachment_max_bitscore: float | None
    background_max_bitscore: float | None
    attachment_best_genome: str | None
    background_best_genome: str | None
    score_margin: float | None


class _SplitVote(_HitAudit):
    gene_id: str
    assigned_clade: str | None
    informative: bool
    match_status: str


class _SplitReport(TypedDict):
    input_status: str
    attachment_clade: str | None
    background_clade: str | None
    votes: list[_SplitVote]
    candidate_genome_id: str | None
    candidate_contig_id: str | None
    query_count: int
    eligible_non_marker_query_count: int
    query_selection_truncated_count: int
    candidate_record_count: int
    candidate_marker_record_count: int
    informative_vote_count: int
    attachment_vote_count: int
    background_vote_count: int
    reference_input_record_count: int
    reference_index_record_count: int
    attachment_reference_record_count: int
    background_reference_record_count: int
    excluded_recipient_reference_count: int
    excluded_marker_reference_count: int
    excluded_candidate_reference_count: int
    outside_split_reference_count: int
    invalid_reference_record_count: int


@dataclass(frozen=True, slots=True)
class _SplitInputs:
    recipient: str | None
    candidate_contig: str | None
    markers: set[str] | None
    attachment: set[str] | None
    background: set[str] | None

    @property
    def marker_ids(self) -> set[str]:
        return self.markers or set()

    @property
    def attachment_taxa(self) -> set[str]:
        return self.attachment or set()

    @property
    def background_taxa(self) -> set[str]:
        return self.background or set()

    @property
    def attachment_clade(self) -> str | None:
        return ",".join(sorted(self.attachment_taxa)) or None

    @property
    def background_clade(self) -> str | None:
        taxa = ",".join(sorted(self.background_taxa))
        return f"complement:{taxa}" if taxa else None


@dataclass(frozen=True, slots=True)
class _CandidateQueries:
    status: str
    items: list[_QueryItem]
    record_ids: set[str]


@dataclass(frozen=True, slots=True)
class _ReferenceIndex:
    alphabet: easel.AA
    sequences: easel.DigitalSequenceBlock[easel.AA]
    records_by_internal_id: dict[str, tuple[str, str]]


@dataclass(slots=True)
class _ReferenceRecords:
    sequences: list[tuple[str, str]]
    records_by_internal_id: dict[str, tuple[str, str]]


def _clean_id(value: object) -> str | None:
    if not isinstance(value, str):
        return None
    cleaned = value.strip()
    return cleaned or None


def _amino_acid_sequence(value: object) -> str | None:
    if not isinstance(value, str):
        return None
    sequence = value.strip().upper()
    if not sequence or any(
        character not in AMINO_ACID_CHARACTERS for character in sequence
    ):
        return None
    return sequence


def _clean_ids(values: object) -> set[str] | None:
    if isinstance(values, (str, bytes)) or not isinstance(values, Collection):
        return None
    cleaned = {_clean_id(value) for value in values}
    if None in cleaned:
        return None
    return {value for value in cleaned if value is not None}


def _opaque_gene_id(record_id: str) -> str:
    payload = record_id.encode()
    digest = hashlib.blake2s(payload, digest_size=6).hexdigest()
    return f"q{digest}"


def _split_report(status: str, inputs: _SplitInputs) -> _SplitReport:
    return {
        "input_status": status,
        "attachment_clade": inputs.attachment_clade,
        "background_clade": inputs.background_clade,
        "votes": [],
        "candidate_genome_id": inputs.recipient,
        "candidate_contig_id": inputs.candidate_contig,
        "query_count": 0,
        "eligible_non_marker_query_count": 0,
        "query_selection_truncated_count": 0,
        "candidate_record_count": 0,
        "candidate_marker_record_count": 0,
        "informative_vote_count": 0,
        "attachment_vote_count": 0,
        "background_vote_count": 0,
        "reference_input_record_count": 0,
        "reference_index_record_count": 0,
        "attachment_reference_record_count": 0,
        "background_reference_record_count": 0,
        "excluded_recipient_reference_count": 0,
        "excluded_marker_reference_count": 0,
        "excluded_candidate_reference_count": 0,
        "outside_split_reference_count": 0,
        "invalid_reference_record_count": 0,
    }


def _query_coverage(hit: plan7.Hit, query_length: int) -> float:
    intervals = sorted(
        (
            int(domain.alignment.hmm_from),
            int(domain.alignment.hmm_to),
        )
        for domain in hit.domains
        if domain.reported
    )
    covered = 0
    current_start: int | None = None
    current_end: int | None = None
    for start, end in intervals:
        if current_start is None or current_end is None:
            current_start, current_end = start, end
        elif start <= current_end + 1:
            current_end = max(current_end, end)
        else:
            covered += current_end - current_start + 1
            current_start, current_end = start, end
    if current_start is not None and current_end is not None:
        covered += current_end - current_start + 1
    return covered / query_length if query_length else 0.0


def _hit_name(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="strict")
    return str(value)


def _input_status(
    candidate_genes: object,
    panel_proteomes: object,
    inputs: _SplitInputs,
) -> str:
    if not isinstance(candidate_genes, Mapping):
        status = "invalid_candidate_genes"
    elif not isinstance(panel_proteomes, Mapping):
        status = "invalid_panel_proteomes"
    elif inputs.recipient is None:
        status = "invalid_recipient_genome"
    elif inputs.candidate_contig is None:
        status = "invalid_candidate_contig_id"
    elif inputs.markers is None:
        status = "invalid_marker_record_ids"
    elif not inputs.attachment:
        status = "invalid_attachment_taxa"
    elif not inputs.background:
        status = "invalid_background_taxa"
    elif inputs.attachment & inputs.background:
        status = "overlapping_split_taxa"
    elif inputs.recipient in inputs.attachment | inputs.background:
        status = "recipient_in_split_taxa"
    else:
        status = "ok"
    return status


def _candidate_record_status(
    record_id: str,
    record_ids: set[str],
    gene_ids: set[str],
    inputs: _SplitInputs,
) -> str:
    is_duplicate_record = record_id in record_ids
    record_ids.add(record_id)
    genome, contig, gene_id = parse_sequence_id(record_id)
    is_duplicate_gene = gene_id in gene_ids
    gene_ids.add(gene_id)
    if is_duplicate_record:
        return "duplicate_candidate_record_id"
    if genome != inputs.recipient:
        return "candidate_recipient_mismatch"
    if contig != inputs.candidate_contig:
        return "candidate_contig_mismatch"
    if not gene_id or gene_id == "unknown_gene":
        return "invalid_candidate_gene_id"
    if is_duplicate_gene:
        return "duplicate_candidate_gene_id"
    return "ok"


def _collect_candidate_queries(
    candidate_genes: Mapping[str, object],
    inputs: _SplitInputs,
    status: str,
    report: _SplitReport,
) -> _CandidateQueries:
    record_ids: set[str] = set()
    gene_ids: set[str] = set()
    query_items: list[_QueryItem] = []
    report["candidate_record_count"] = len(candidate_genes)
    for record_value, sequence_value in sorted(
        candidate_genes.items(),
        key=lambda item: str(item[0]),
    ):
        record_id = _clean_id(record_value)
        if record_id is None:
            if status == "ok":
                status = "invalid_candidate_record_id"
                report["input_status"] = status
            continue
        if status == "ok":
            status = _candidate_record_status(record_id, record_ids, gene_ids, inputs)
            report["input_status"] = status
        else:
            record_ids.add(record_id)
        if record_id in inputs.marker_ids:
            report["candidate_marker_record_count"] += 1
        else:
            query_items.append((record_id, sequence_value))

    report["eligible_non_marker_query_count"] = len(query_items)
    if status == "ok" and len(query_items) < PHMMER_MIN_QUERY_COUNT:
        status = "requires_at_least_three_non_marker_queries"
        report["input_status"] = status
    query_items.sort(
        key=lambda item: (
            hashlib.blake2s(item[0].encode(), digest_size=8).hexdigest(),
            item[0],
        )
    )
    report["query_selection_truncated_count"] = max(
        0,
        len(query_items) - PHMMER_MAX_QUERY_COUNT,
    )
    query_items = query_items[:PHMMER_MAX_QUERY_COUNT]
    report["query_count"] = len(query_items)
    return _CandidateQueries(status=status, items=query_items, record_ids=record_ids)


def _append_reference_records(
    genome: str,
    records: Mapping[str, object],
    inputs: _SplitInputs,
    candidate_record_ids: set[str],
    reference_records: _ReferenceRecords,
    report: _SplitReport,
) -> None:
    for record_value, sequence_value in sorted(
        records.items(),
        key=lambda item: str(item[0]),
    ):
        record_id = _clean_id(record_value)
        if record_id in inputs.marker_ids:
            report["excluded_marker_reference_count"] += 1
            continue
        if record_id in candidate_record_ids:
            report["excluded_candidate_reference_count"] += 1
            continue
        sequence = _amino_acid_sequence(sequence_value)
        if record_id is None or sequence is None:
            report["invalid_reference_record_count"] += 1
            continue
        internal_id = f"r{len(reference_records.sequences):08d}"
        reference_records.records_by_internal_id[internal_id] = (genome, record_id)
        reference_records.sequences.append((internal_id, sequence))
        report["reference_index_record_count"] += 1
        if genome in inputs.attachment_taxa:
            report["attachment_reference_record_count"] += 1
        else:
            report["background_reference_record_count"] += 1


def _collect_reference_sequences(
    panel_proteomes: Mapping[str, Mapping[str, object]],
    inputs: _SplitInputs,
    candidate_record_ids: set[str],
    report: _SplitReport,
) -> tuple[_ReferenceRecords, str]:
    reference_records = _ReferenceRecords(sequences=[], records_by_internal_id={})
    split_taxa = inputs.attachment_taxa | inputs.background_taxa
    for genome_value, records in sorted(
        panel_proteomes.items(),
        key=lambda item: str(item[0]),
    ):
        genome = _clean_id(genome_value)
        if not isinstance(records, Mapping):
            report["invalid_reference_record_count"] += 1
            continue
        report["reference_input_record_count"] += len(records)
        if genome is None:
            report["invalid_reference_record_count"] += len(records)
            continue
        if genome == inputs.recipient:
            report["excluded_recipient_reference_count"] += len(records)
            continue
        if genome not in split_taxa:
            report["outside_split_reference_count"] += len(records)
            continue
        _append_reference_records(
            genome,
            records,
            inputs,
            candidate_record_ids,
            reference_records,
            report,
        )

    if (
        report["attachment_reference_record_count"] == 0
        or report["background_reference_record_count"] == 0
    ):
        return reference_records, "missing_split_reference_side"
    return reference_records, "ok"


def _reference_index(
    reference_records: _ReferenceRecords,
) -> _ReferenceIndex:
    alphabet = easel.Alphabet.amino()
    digital_sequences = easel.DigitalSequenceBlock(
        alphabet,
        [
            easel.TextSequence(
                name=internal_id.encode(),
                sequence=sequence,
            ).digitize(alphabet)
            for internal_id, sequence in reference_records.sequences
        ],
    )
    return _ReferenceIndex(
        alphabet=alphabet,
        sequences=digital_sequences,
        records_by_internal_id=reference_records.records_by_internal_id,
    )


def _hit_audit() -> _HitAudit:
    return {
        "raw_hit_count": 0,
        "evalue_filtered_hit_count": 0,
        "coverage_filtered_hit_count": 0,
        "qualifying_hit_record_count": 0,
        "minimum_observed_evalue": None,
        "maximum_observed_query_coverage": None,
        "attachment_qualifying_genome_count": 0,
        "background_qualifying_genome_count": 0,
        "attachment_max_bitscore": None,
        "background_max_bitscore": None,
        "attachment_best_genome": None,
        "background_best_genome": None,
        "score_margin": None,
    }


def _score_hits(
    sequence: str,
    gene_id: str,
    reference_index: _ReferenceIndex,
) -> tuple[_HitAudit, dict[str, float]]:
    query = easel.TextSequence(
        name=gene_id.encode(),
        sequence=sequence,
    ).digitize(reference_index.alphabet)
    hits = next(
        hmmer.phmmer(
            [query],  # ty: ignore[invalid-argument-type]  # pyhmmer 0.12 types amino-acid sequences invariantly.
            reference_index.sequences,  # ty: ignore[invalid-argument-type]  # pyhmmer 0.12 accepts this block at runtime but types it invariantly.
            cpus=1,
            builder=plan7.Builder(
                reference_index.alphabet,
                seed=PHMMER_BUILDER_SEED,
            ),
            E=PHMMER_REPORT_EVALUE,
            domE=PHMMER_REPORT_EVALUE,
        )
    )
    audit = _hit_audit()
    scores_by_genome: dict[str, float] = {}
    for hit in hits:
        audit["raw_hit_count"] += 1
        evalue = float(hit.evalue)
        coverage = _query_coverage(hit, len(sequence))
        current_evalue = audit["minimum_observed_evalue"]
        current_coverage = audit["maximum_observed_query_coverage"]
        audit["minimum_observed_evalue"] = (
            evalue if current_evalue is None else min(current_evalue, evalue)
        )
        audit["maximum_observed_query_coverage"] = (
            coverage if current_coverage is None else max(current_coverage, coverage)
        )
        if evalue > PHMMER_MAX_EVALUE:
            audit["evalue_filtered_hit_count"] += 1
            continue
        if coverage < PHMMER_MIN_QUERY_COVERAGE:
            audit["coverage_filtered_hit_count"] += 1
            continue
        internal_id = _hit_name(hit.name)
        genome, _reference_record_id = reference_index.records_by_internal_id[
            internal_id
        ]
        audit["qualifying_hit_record_count"] += 1
        scores_by_genome[genome] = max(
            scores_by_genome.get(genome, float("-inf")),
            float(hit.score),
        )
    return audit, scores_by_genome


def _best_score(scores: Mapping[str, float]) -> tuple[str | None, float | None]:
    if not scores:
        return None, None
    best_genome = min(scores, key=lambda genome: (-scores[genome], genome))
    return best_genome, scores[best_genome]


def _summarize_side_scores(
    audit: _HitAudit,
    scores_by_genome: Mapping[str, float],
    inputs: _SplitInputs,
) -> tuple[float | None, float | None]:
    attachment_scores = {
        genome: score
        for genome, score in scores_by_genome.items()
        if genome in inputs.attachment_taxa
    }
    background_scores = {
        genome: score
        for genome, score in scores_by_genome.items()
        if genome in inputs.background_taxa
    }
    best_attachment, max_attachment = _best_score(attachment_scores)
    best_background, max_background = _best_score(background_scores)
    audit["attachment_qualifying_genome_count"] = len(attachment_scores)
    audit["background_qualifying_genome_count"] = len(background_scores)
    audit["attachment_best_genome"] = best_attachment
    audit["background_best_genome"] = best_background
    audit["attachment_max_bitscore"] = max_attachment
    audit["background_max_bitscore"] = max_background
    return max_attachment, max_background


def _assign_by_margin(
    audit: _HitAudit,
    max_attachment: float | None,
    max_background: float | None,
    inputs: _SplitInputs,
) -> tuple[str | None, bool, str]:
    if max_attachment is None or max_background is None:
        return None, False, "missing_split_side_hit"
    margin = max_attachment - max_background
    audit["score_margin"] = margin
    if margin >= PHMMER_MIN_SCORE_MARGIN:
        return inputs.attachment_clade, True, "attachment_margin"
    if -margin >= PHMMER_MIN_SCORE_MARGIN:
        return inputs.background_clade, True, "background_margin"
    return None, False, "score_margin_below_threshold"


def _split_vote(
    record_id: str,
    sequence_value: object,
    reference_index: _ReferenceIndex,
    inputs: _SplitInputs,
    status: str,
) -> _SplitVote:
    gene_id = _opaque_gene_id(record_id)
    assigned_clade = None
    informative = False
    match_status = status
    sequence = _amino_acid_sequence(sequence_value)
    if status == "ok" and sequence is not None:
        audit, scores_by_genome = _score_hits(sequence, gene_id, reference_index)
        max_attachment, max_background = _summarize_side_scores(
            audit,
            scores_by_genome,
            inputs,
        )
        assigned_clade, informative, match_status = _assign_by_margin(
            audit,
            max_attachment,
            max_background,
            inputs,
        )
    else:
        audit = _hit_audit()
        if status == "ok":
            match_status = "invalid_candidate_gene"
    return {
        "gene_id": gene_id,
        "assigned_clade": assigned_clade,
        "informative": informative,
        "match_status": match_status,
        **audit,
    }


def _record_vote(report: _SplitReport, vote: _SplitVote) -> None:
    if vote["informative"]:
        report["informative_vote_count"] += 1
    if vote["match_status"] == "attachment_margin":
        report["attachment_vote_count"] += 1
    elif vote["match_status"] == "background_margin":
        report["background_vote_count"] += 1
    report["votes"].append(vote)


def assign_contig_gene_split_votes(  # noqa: PLR0913  # Preserve the frozen benchmark API.
    candidate_genes: Mapping[str, object],
    panel_proteomes: Mapping[str, Mapping[str, object]],
    *,
    recipient_genome: str,
    candidate_contig_id: str,
    marker_record_ids: Collection[str],
    attachment_taxa: Collection[str],
    background_taxa: Collection[str],
) -> dict:
    """Score three non-marker proteins against a fixed attachment/background split.

    The function runs ``phmmer`` in one thread. It accepts no event, donor, or
    truth input. Panel genome IDs provide the only side assignments.
    """
    inputs = _SplitInputs(
        recipient=_clean_id(recipient_genome),
        candidate_contig=_clean_id(candidate_contig_id),
        markers=_clean_ids(marker_record_ids),
        attachment=_clean_ids(attachment_taxa),
        background=_clean_ids(background_taxa),
    )
    status = _input_status(candidate_genes, panel_proteomes, inputs)
    report = _split_report(status, inputs)
    if not isinstance(candidate_genes, Mapping):
        return dict(report)

    candidates = _collect_candidate_queries(candidate_genes, inputs, status, report)
    reference_records = _ReferenceRecords(sequences=[], records_by_internal_id={})
    status = candidates.status
    if status == "ok":
        reference_records, status = _collect_reference_sequences(
            panel_proteomes,
            inputs,
            candidates.record_ids,
            report,
        )
        report["input_status"] = status
    references = _reference_index(reference_records)
    for record_id, sequence_value in candidates.items:
        vote = _split_vote(record_id, sequence_value, references, inputs, status)
        _record_vote(report, vote)
    return dict(report)
