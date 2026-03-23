from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path("/home/fschulz/dev/software/sgtree")
BENCH_ROOT = PROJECT_ROOT / "runs" / "contig_variant_proof_of_concept" / "benchmarks"
EVAL_ROOT = PROJECT_ROOT / "eval"
FULL_100_ROOT = EVAL_ROOT / "full" / "100gen"
PILOT_100_ROOT = EVAL_ROOT / "pilot" / "100gen"
FULL_MANIFEST = EVAL_ROOT / "full" / "manifest.json"
FULL_README = EVAL_ROOT / "full" / "README.md"
TOP_README = EVAL_ROOT / "README.md"
PILOT_MANIFEST = EVAL_ROOT / "pilot" / "pilot_manifest.json"

PILOT_100 = {
    "alpha": ("genus__seed_101", "combined", "solo_contig"),
    "bacteroidota": ("genus__seed_101", "combined", "solo_contig"),
}
VARIANTS = ("solo_contig", "merged_contam", "native_contig")


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def matrix_completeness(matrix_path: Path, genomes: list[str], extra_cols: dict[str, dict] | None = None) -> tuple[list[dict], list[str], dict[str, set[str]]]:
    df = pd.read_csv(matrix_path, index_col=0)
    truth_markers = [str(item) for item in df.index.tolist()]
    marker_presence: dict[str, set[str]] = {
        marker: {
            str(genome)
            for genome in df.columns
            if float(df.loc[marker, genome]) > 0
        }
        for marker in truth_markers
    }
    rows: list[dict] = []
    for genome in genomes:
        observed = 0
        in_matrix = genome in df.columns
        if in_matrix:
            observed = int((df[genome].fillna(0).astype(float) > 0).sum())
        row = {
            "genome_id": genome,
            "truth_marker_count": len(truth_markers),
            "observed_marker_count": observed,
            "truth_marker_completeness_fraction": (observed / len(truth_markers)) if truth_markers else 0.0,
            "present_in_marker_matrix": "yes" if in_matrix else "no",
        }
        if extra_cols and genome in extra_cols:
            row.update(extra_cols[genome])
        rows.append(row)
    return rows, truth_markers, marker_presence


def variant_contig_id(event: dict, variant: str) -> str:
    if variant == "solo_contig":
        return str(event.get("contaminant_contig_id", ""))
    if variant == "merged_contam":
        return f"contig__contam_merged__{event['recipient_genome']}"
    if variant == "native_contig":
        return str(event.get("native_contig_id", ""))
    raise ValueError(variant)


def rebuild_full_100gen() -> None:
    if FULL_100_ROOT.exists():
        shutil.rmtree(FULL_100_ROOT)
    FULL_100_ROOT.mkdir(parents=True, exist_ok=True)

    comparison = pd.read_csv(EVAL_ROOT / "full" / "summaries" / "100gen_comparison.tsv", sep="\t")
    run_index_rows: list[dict] = []
    panel_index_rows: list[dict] = []

    for benchmark_manifest in sorted(BENCH_ROOT.glob("*/*/benchmark_manifest.json")):
        bench_dir = benchmark_manifest.parent
        lineage = bench_dir.parent.name
        panel = bench_dir.name
        manifest = json.loads(benchmark_manifest.read_text())
        panel_target = FULL_100_ROOT / lineage / panel
        panel_target.mkdir(parents=True, exist_ok=True)

        selected_markers = [str(item) for item in manifest.get("selected_markers", [])]
        selected_genomes = [str(item) for item in manifest.get("selected_genomes", [])]
        selected_taxonomy = manifest.get("selected_genome_taxonomy", [])
        taxonomy_map = {row.get("genome_id"): row for row in selected_taxonomy if row.get("genome_id")}

        panel_manifest = {
            "source_type": "exact_benchmark_manifest",
            "lineage": lineage,
            "panel": panel,
            "scope": manifest.get("taxonomic_scope", ""),
            "seed": manifest.get("seed"),
            "selected_marker_count": len(selected_markers),
            "selected_markers_are_uni56_subset": True,
            "selected_genome_count": len(selected_genomes),
            "scenarios": [scenario["name"] for scenario in manifest.get("scenarios", [])],
            "variants": list(VARIANTS),
            "metadata_exactness": "exact",
            "reference_marker_presence_exact": True,
        }
        (panel_target / "panel_manifest.json").write_text(json.dumps(panel_manifest, indent=2) + "\n")
        (panel_target / "source_benchmark_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        write_tsv(panel_target / "selected_markers.tsv", [{"marker": marker} for marker in selected_markers])
        write_tsv(panel_target / "selected_genomes.tsv", selected_taxonomy)

        panel_summary = comparison[(comparison["lineage"] == lineage) & (comparison["panel"] == panel)].copy()
        panel_summary.to_csv(panel_target / "panel_run_summary.tsv", sep="\t", index=False)
        panel_index_rows.append(
            {
                "lineage": lineage,
                "panel": panel,
                "scope": manifest.get("taxonomic_scope", ""),
                "seed": manifest.get("seed"),
                "selected_marker_count": len(selected_markers),
                "selected_genome_count": len(selected_genomes),
                "scenario_count": len(manifest.get("scenarios", [])),
            }
        )

        for scenario_meta in manifest.get("scenarios", []):
            scenario = str(scenario_meta["name"])
            events_path = bench_dir / "scenarios" / scenario / "events.tsv"
            events_df = pd.read_csv(events_path, sep="\t")
            reference_matrix = bench_dir / "scenarios" / scenario / "reference_run" / "marker_count_matrix.csv"
            completeness_rows, truth_markers, marker_presence = matrix_completeness(reference_matrix, selected_genomes, taxonomy_map)

            for variant in VARIANTS:
                run_target = panel_target / "runs" / f"{scenario}__{variant}"
                run_target.mkdir(parents=True, exist_ok=True)
                write_tsv(run_target / "genome_truth_marker_completeness.tsv", completeness_rows)
                write_tsv(run_target / "genome_manifest.tsv", [{"genome_id": genome, **taxonomy_map.get(genome, {})} for genome in selected_genomes])

                introduced_rows: list[dict] = []
                for row in events_df.to_dict("records"):
                    marker = str(row.get("marker", ""))
                    present = sorted(marker_presence.get(marker, set()))
                    introduced_rows.append(
                        {
                            "event_index": row.get("event_index", ""),
                            "scenario": row.get("scenario", scenario),
                            "event_type": row.get("event_type", ""),
                            "taxonomic_scope": row.get("taxonomic_scope", ""),
                            "taxonomic_scope_label": row.get("taxonomic_scope_label", ""),
                            "recipient_genome": row.get("recipient_genome", ""),
                            "recipient_group": row.get("recipient_group", ""),
                            "marker": marker,
                            "donor_genome": row.get("donor_genome", ""),
                            "donor_group": row.get("donor_group", ""),
                            "source_relation": row.get("source_relation", ""),
                            "native_record_id": row.get("native_record_id", ""),
                            "native_contig_id": row.get("native_contig_id", ""),
                            "donor_record_id": row.get("donor_record_id", ""),
                            "donor_contig_id": row.get("donor_contig_id", ""),
                            "contaminant_record_id": row.get("contaminant_record_id", ""),
                            "contaminant_contig_id": row.get("contaminant_contig_id", ""),
                            "variant_contig_id": variant_contig_id(row, variant),
                            "expected_replacement_outcome": row.get("expected_replacement_outcome", ""),
                            "native_degrade_fraction": row.get("native_degrade_fraction", ""),
                            "reference_marker_present_genome_count": len(present),
                            "reference_marker_present_genomes": ",".join(present),
                        }
                    )
                write_tsv(run_target / "introduced_markers.tsv", introduced_rows)

                comparison_rows = panel_summary[
                    (panel_summary["scenario"] == scenario) & (panel_summary["variant"] == variant)
                ]
                comparison_row = comparison_rows.iloc[0] if not comparison_rows.empty else None
                run_manifest = {
                    "lineage": lineage,
                    "panel": panel,
                    "scenario": scenario,
                    "variant": variant,
                    "source_type": "exact_source_benchmark_plus_latest_full_run_summary",
                    "metadata_exactness": "exact",
                    "reference_marker_presence_exact": True,
                    "replacement_events": int(comparison_row.get("replacement_events", 0)) if comparison_row is not None else 0,
                    "duplicate_events": int(comparison_row.get("duplicate_events", 0)) if comparison_row is not None else 0,
                    "replacement_contaminant_removed": int(comparison_row.get("replacement_contaminant_removed", 0)) if comparison_row is not None else 0,
                    "duplicate_contaminant_removed": int(comparison_row.get("duplicate_contaminant_removed", 0)) if comparison_row is not None else 0,
                    "singleton_collateral_removed_count": int(comparison_row.get("singleton_collateral_removed_count", 0)) if comparison_row is not None else 0,
                    "singleton_collateral_genome_count": int(comparison_row.get("singleton_collateral_genome_count", 0)) if comparison_row is not None else 0,
                    "tree_rf_delta": float(comparison_row.get("tree_rf_delta", 0.0)) if comparison_row is not None else 0.0,
                    "result_status": "ok" if comparison_row is not None else "not_run_in_latest_full_recipient_consensus",
                    "event_source": "exact_source_scenario_events_tsv",
                    "completeness_source": "exact_source_reference_run_marker_count_matrix",
                }
                (run_target / "run_manifest.json").write_text(json.dumps(run_manifest, indent=2) + "\n")

                run_index_rows.append(
                    {
                        "lineage": lineage,
                        "panel": panel,
                        "scenario": scenario,
                        "variant": variant,
                        "replacement_events": int(comparison_row.get("replacement_events", 0)) if comparison_row is not None else 0,
                        "duplicate_events": int(comparison_row.get("duplicate_events", 0)) if comparison_row is not None else 0,
                        "replacement_contaminant_removed": int(comparison_row.get("replacement_contaminant_removed", 0)) if comparison_row is not None else 0,
                        "duplicate_contaminant_removed": int(comparison_row.get("duplicate_contaminant_removed", 0)) if comparison_row is not None else 0,
                        "singleton_collateral_removed_count": int(comparison_row.get("singleton_collateral_removed_count", 0)) if comparison_row is not None else 0,
                        "tree_rf_delta": float(comparison_row.get("tree_rf_delta", 0.0)) if comparison_row is not None else 0.0,
                        "result_status": "ok" if comparison_row is not None else "not_run_in_latest_full_recipient_consensus",
                        "metadata_exactness": "exact",
                    }
                )

        write_tsv(panel_target / "truth_markers.tsv", [{"marker": marker} for marker in truth_markers])

    write_tsv(FULL_100_ROOT / "panel_index.tsv", panel_index_rows)
    write_tsv(FULL_100_ROOT / "run_index_compact.tsv", run_index_rows)
    comparison.to_csv(FULL_100_ROOT / "run_index.tsv", sep="\t", index=False)


def rebuild_pilot_100gen() -> None:
    if PILOT_100_ROOT.exists():
        shutil.rmtree(PILOT_100_ROOT)
    PILOT_100_ROOT.mkdir(parents=True, exist_ok=True)
    for lineage, (panel, scenario, variant) in PILOT_100.items():
        src = FULL_100_ROOT / lineage / panel / "runs" / f"{scenario}__{variant}"
        dst = PILOT_100_ROOT / lineage / panel / f"{scenario}__{variant}"
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(src, dst)


def update_docs() -> None:
    manifest = json.loads(FULL_MANIFEST.read_text())
    manifest["datasets"]["100gen"]["metadata_exactness"] = "exact"
    manifest["datasets"]["100gen"]["reference_marker_presence_exact"] = True
    manifest["datasets"]["100gen"]["event_tables"] = "exact scenario events copied from regenerated source benchmark manifests"
    manifest["datasets"]["100gen"]["completeness_metric"] = "truth-marker panel completeness derived from exact source scenario reference_run marker_count_matrix.csv"
    manifest["datasets"]["100gen"]["source_root"] = str(BENCH_ROOT)
    manifest["datasets"]["100gen"]["limitations"] = []
    FULL_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")

    full_readme = FULL_README.read_text()
    full_readme = full_readme.replace(
        "100-gen event tables are reconstructed from contaminant headers preserved in the latest full run outputs because the original source benchmark directories were not available locally during this packaging pass.",
        "100-gen event tables are exact copies derived from regenerated source benchmark manifests and scenario event tables.",
    )
    full_readme = full_readme.replace(
        "In 100-gen `combined` runs, inserted contaminant headers do not always preserve enough information to label each row back to `duplicate` vs `replacement` with certainty. Those rows are kept as contaminant insertions with this limitation documented in the run manifests.",
        "In 100-gen panels, variant-specific contaminant contig placement is computed exactly from the source scenario events plus the documented rewrite rules for `solo_contig`, `merged_contam`, and `native_contig`.",
    )
    full_readme = full_readme.replace(
        "Exact per-marker reference-genome presence is available for 50-gen panels but not fully recoverable for 100-gen panels from the remaining local files.",
        "Exact per-marker reference-genome presence is available for both 50-gen and regenerated 100-gen panels.",
    )
    FULL_README.write_text(full_readme)

    top_readme = TOP_README.read_text()
    top_readme = top_readme.replace(
        "For `100gen/`, genome membership, truth-marker sets, run summaries, and introduced-marker headers are retained, but exact per-marker reference-genome presence is not fully recoverable from the surviving files.",
        "For `100gen/`, exact source manifests, exact scenario event tables, exact variant-specific contaminant placement, genome membership, and exact per-marker reference-genome presence are retained after regeneration.",
    )
    TOP_README.write_text(top_readme)

    pilot_manifest = json.loads(PILOT_MANIFEST.read_text())
    pilot_manifest["limitations"] = [
        "50-gen pilot ground truth is exact and includes per-marker reference-genome presence.",
        "100-gen pilot ground truth is exact and includes per-marker reference-genome presence after source benchmark regeneration.",
    ]
    PILOT_MANIFEST.write_text(json.dumps(pilot_manifest, indent=2) + "\n")


def main() -> None:
    rebuild_full_100gen()
    rebuild_pilot_100gen()
    update_docs()


if __name__ == "__main__":
    main()
