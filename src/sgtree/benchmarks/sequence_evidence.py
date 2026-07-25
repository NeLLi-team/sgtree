"""Benchmark-only contig evidence from protein-sequence matches."""

from __future__ import annotations

import hashlib
from collections.abc import Collection, Mapping

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
        character not in AMINO_ACID_CHARACTERS
        for character in sequence
    ):
        return None
    return sequence


def _marker_ids(values: object) -> set[str] | None:
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


def _taxa(values: object) -> set[str] | None:
    if isinstance(values, (str, bytes)) or not isinstance(values, Collection):
        return None
    cleaned = {_clean_id(value) for value in values}
    if None in cleaned:
        return None
    return {value for value in cleaned if value is not None}


def _split_report(status: str, attachment: set[str], background: set[str]) -> dict:
    return {
        "input_status": status,
        "attachment_clade": ",".join(sorted(attachment)) or None,
        "background_clade": (
            "complement:" + ",".join(sorted(background))
            if background
            else None
        ),
        "votes": [],
        "candidate_genome_id": None,
        "candidate_contig_id": None,
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


def _query_coverage(hit: object, query_length: int) -> float:
    intervals = sorted(
        (
            int(domain.alignment.hmm_from),
            int(domain.alignment.hmm_to),
        )
        for domain in hit.domains
        if domain.reported
    )
    covered = 0
    current_start = None
    current_end = None
    for start, end in intervals:
        if current_start is None:
            current_start, current_end = start, end
        elif start <= current_end + 1:
            current_end = max(current_end, end)
        else:
            covered += current_end - current_start + 1
            current_start, current_end = start, end
    if current_start is not None:
        covered += current_end - current_start + 1
    return covered / query_length if query_length else 0.0


def _hit_name(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="strict")
    return str(value)


def assign_contig_gene_split_votes(
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
    recipient = _clean_id(recipient_genome)
    candidate_contig = _clean_id(candidate_contig_id)
    markers = _marker_ids(marker_record_ids)
    attachment = _taxa(attachment_taxa)
    background = _taxa(background_taxa)
    attachment_for_report = attachment or set()
    background_for_report = background or set()

    if not isinstance(candidate_genes, Mapping):
        status = "invalid_candidate_genes"
    elif not isinstance(panel_proteomes, Mapping):
        status = "invalid_panel_proteomes"
    elif recipient is None:
        status = "invalid_recipient_genome"
    elif candidate_contig is None:
        status = "invalid_candidate_contig_id"
    elif markers is None:
        status = "invalid_marker_record_ids"
    elif attachment is None or not attachment:
        status = "invalid_attachment_taxa"
    elif background is None or not background:
        status = "invalid_background_taxa"
    elif attachment & background:
        status = "overlapping_split_taxa"
    elif recipient in attachment | background:
        status = "recipient_in_split_taxa"
    else:
        status = "ok"

    report = _split_report(status, attachment_for_report, background_for_report)
    report["candidate_genome_id"] = recipient
    report["candidate_contig_id"] = candidate_contig
    if not isinstance(candidate_genes, Mapping):
        return report
    report["candidate_record_count"] = len(candidate_genes)
    markers = markers or set()
    candidate_ids: set[str] = set()
    candidate_gene_ids: set[str] = set()
    query_items = []
    for record_value, sequence_value in sorted(
        candidate_genes.items(),
        key=lambda item: str(item[0]),
    ):
        record_id = _clean_id(record_value)
        if status == "ok" and record_id is None:
            status = "invalid_candidate_record_id"
            report["input_status"] = status
            continue
        if record_id is None:
            continue
        if status == "ok" and record_id in candidate_ids:
            status = "duplicate_candidate_record_id"
            report["input_status"] = status
        candidate_ids.add(record_id)
        genome, contig, gene_id = parse_sequence_id(record_id)
        if status == "ok" and genome != recipient:
            status = "candidate_recipient_mismatch"
            report["input_status"] = status
        elif status == "ok" and contig != candidate_contig:
            status = "candidate_contig_mismatch"
            report["input_status"] = status
        elif status == "ok" and (not gene_id or gene_id == "unknown_gene"):
            status = "invalid_candidate_gene_id"
            report["input_status"] = status
        elif status == "ok" and gene_id in candidate_gene_ids:
            status = "duplicate_candidate_gene_id"
            report["input_status"] = status
        candidate_gene_ids.add(gene_id)
        if record_id in markers:
            report["candidate_marker_record_count"] += 1
            continue
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

    reference_sequences = []
    reference_lookup: dict[str, tuple[str, str]] = {}
    if status == "ok":
        split_taxa = attachment | background
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
            if genome == recipient:
                report["excluded_recipient_reference_count"] += len(records)
                continue
            if genome not in split_taxa:
                report["outside_split_reference_count"] += len(records)
                continue
            for record_value, sequence_value in sorted(
                records.items(),
                key=lambda item: str(item[0]),
            ):
                record_id = _clean_id(record_value)
                if record_id in markers:
                    report["excluded_marker_reference_count"] += 1
                    continue
                if record_id in candidate_ids:
                    report["excluded_candidate_reference_count"] += 1
                    continue
                sequence = _amino_acid_sequence(sequence_value)
                if record_id is None or sequence is None:
                    report["invalid_reference_record_count"] += 1
                    continue
                internal_id = f"r{len(reference_sequences):08d}"
                reference_lookup[internal_id] = (genome, record_id)
                reference_sequences.append((internal_id, sequence))
                report["reference_index_record_count"] += 1
                side_count = (
                    "attachment_reference_record_count"
                    if genome in attachment
                    else "background_reference_record_count"
                )
                report[side_count] += 1
        if (
            report["attachment_reference_record_count"] == 0
            or report["background_reference_record_count"] == 0
        ):
            status = "missing_split_reference_side"
            report["input_status"] = status

    alphabet = easel.Alphabet.amino()
    digital_references = easel.DigitalSequenceBlock(
        alphabet,
        [
            easel.TextSequence(
                name=internal_id.encode(),
                sequence=sequence,
            ).digitize(alphabet)
            for internal_id, sequence in reference_sequences
        ],
    )
    votes = []
    for record_id, sequence_value in query_items:
        gene_id = _opaque_gene_id(record_id)
        sequence = _amino_acid_sequence(sequence_value)
        audit = {
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
        assigned_clade = None
        informative = False
        match_status = status

        if status == "ok" and (record_id is None or sequence is None):
            match_status = "invalid_candidate_gene"
        elif status == "ok":
            query = easel.TextSequence(
                name=gene_id.encode(),
                sequence=sequence,
            ).digitize(alphabet)
            hits = next(
                hmmer.phmmer(
                    [query],
                    digital_references,
                    cpus=1,
                    builder=plan7.Builder(alphabet, seed=PHMMER_BUILDER_SEED),
                    E=PHMMER_REPORT_EVALUE,
                    domE=PHMMER_REPORT_EVALUE,
                )
            )
            scores_by_genome: dict[str, float] = {}
            for hit in hits:
                audit["raw_hit_count"] += 1
                evalue = float(hit.evalue)
                coverage = _query_coverage(hit, len(sequence))
                current_evalue = audit["minimum_observed_evalue"]
                current_coverage = audit["maximum_observed_query_coverage"]
                audit["minimum_observed_evalue"] = (
                    evalue
                    if current_evalue is None
                    else min(float(current_evalue), evalue)
                )
                audit["maximum_observed_query_coverage"] = (
                    coverage
                    if current_coverage is None
                    else max(float(current_coverage), coverage)
                )
                if evalue > PHMMER_MAX_EVALUE:
                    audit["evalue_filtered_hit_count"] += 1
                    continue
                if coverage < PHMMER_MIN_QUERY_COVERAGE:
                    audit["coverage_filtered_hit_count"] += 1
                    continue
                internal_id = _hit_name(hit.name)
                genome, _reference_record_id = reference_lookup[internal_id]
                audit["qualifying_hit_record_count"] += 1
                scores_by_genome[genome] = max(
                    scores_by_genome.get(genome, float("-inf")),
                    float(hit.score),
                )

            attachment_scores = {
                genome: score
                for genome, score in scores_by_genome.items()
                if genome in attachment
            }
            background_scores = {
                genome: score
                for genome, score in scores_by_genome.items()
                if genome in background
            }
            audit["attachment_qualifying_genome_count"] = len(attachment_scores)
            audit["background_qualifying_genome_count"] = len(background_scores)
            if attachment_scores:
                best_a = min(
                    attachment_scores,
                    key=lambda genome: (-attachment_scores[genome], genome),
                )
                audit["attachment_best_genome"] = best_a
                audit["attachment_max_bitscore"] = attachment_scores[best_a]
            if background_scores:
                best_b = min(
                    background_scores,
                    key=lambda genome: (-background_scores[genome], genome),
                )
                audit["background_best_genome"] = best_b
                audit["background_max_bitscore"] = background_scores[best_b]

            max_a = audit["attachment_max_bitscore"]
            max_b = audit["background_max_bitscore"]
            if max_a is None or max_b is None:
                match_status = "missing_split_side_hit"
            else:
                margin = float(max_a) - float(max_b)
                audit["score_margin"] = margin
                if margin >= PHMMER_MIN_SCORE_MARGIN:
                    assigned_clade = report["attachment_clade"]
                    informative = True
                    match_status = "attachment_margin"
                    report["attachment_vote_count"] += 1
                elif -margin >= PHMMER_MIN_SCORE_MARGIN:
                    assigned_clade = report["background_clade"]
                    informative = True
                    match_status = "background_margin"
                    report["background_vote_count"] += 1
                else:
                    match_status = "score_margin_below_threshold"
        if informative:
            report["informative_vote_count"] += 1
        votes.append(
            {
                "gene_id": gene_id,
                "assigned_clade": assigned_clade,
                "informative": informative,
                "match_status": match_status,
                **audit,
            }
        )

    report["votes"] = votes
    return report
