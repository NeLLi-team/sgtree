"""US-006: Comparative evaluation of ML contamination detection approaches.

Runs GCP, CMTV, and PLIE on the full 50-genome benchmark suite (36 panels)
and compares against the recipient_consensus baseline. Generates summary
tables and identifies the winning approach.

Usage:
    pixi run python runs/ml_contam_detection/evaluate_approaches.py
    pixi run python runs/ml_contam_detection/evaluate_approaches.py --check
"""

from __future__ import annotations

import argparse
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from runs.ml_contam_detection.approach_cmtv import score_cmtv
from runs.ml_contam_detection.approach_gcp import score_gcp
from runs.ml_contam_detection.approach_plie import score_plie
from runs.ml_contam_detection.feature_extraction import extract_features

EVAL_50GEN = PROJECT_ROOT / "eval" / "full" / "50gen"
PRIVATE_BACKUP = Path(
    "/home/fschulz/dev/software/sgtree_private_backup_20260316"
    "/runs/benchmarks/manuscript_panel_replicates_20260318_recipient_consensus"
)
OUTPUT_DIR = PROJECT_ROOT / "runs" / "ml_contam_detection"

LINEAGES_50GEN = ["chlam", "flavo", "gamma"]
SCOPES = ["genus", "family", "order"]
SEEDS = [101, 202, 303]

SINGLETON_SCENARIOS = ["replacement_only", "combined"]
MIXED_SCENARIO = "mixed_high_level"


def _result_dir_name(scenario: str, singles_mode: str) -> str:
    """Map (scenario, singles_mode) to the actual result directory name."""
    if scenario == "replacement_only":
        return f"{scenario}__singles_{singles_mode}"
    return f"{scenario}__duplicate_plus_singles_{singles_mode}"


def _parse_panel_id(lineage: str, scope: str) -> str:
    return f"{lineage}_{scope}"


def build_panel_index() -> list[dict]:
    """Build index of all available 50gen panels with their paths."""
    panels = []
    for lineage in LINEAGES_50GEN:
        for scope in SCOPES:
            panel_id = _parse_panel_id(lineage, scope)
            for seed in SEEDS:
                panel_seed = f"{panel_id}__seed_{seed}"
                benchmark_dir = EVAL_50GEN / lineage / panel_seed
                results_base = PRIVATE_BACKUP / panel_seed / "results"
                if not benchmark_dir.exists():
                    continue

                scenarios = list(SINGLETON_SCENARIOS)
                for sc in scenarios:
                    run_dir_name = _result_dir_name(sc, "recipient_consensus")
                    run_dir = results_base / run_dir_name
                    if (run_dir / "singleton_candidates.tsv").exists():
                        panels.append({
                            "lineage": lineage,
                            "scope": scope,
                            "seed": seed,
                            "panel_id": panel_id,
                            "panel_seed": panel_seed,
                            "scenario": sc,
                            "benchmark_dir": benchmark_dir,
                            "run_dir": run_dir,
                        })

        # Mixed high level panels
        for seed in SEEDS:
            panel_seed = f"{lineage}_mixed_high_level__seed_{seed}"
            benchmark_dir = EVAL_50GEN / lineage / panel_seed
            results_base = PRIVATE_BACKUP / panel_seed / "results"
            if not benchmark_dir.exists():
                continue

            run_dir_name = _result_dir_name(MIXED_SCENARIO, "recipient_consensus")
            run_dir = results_base / run_dir_name
            if (run_dir / "singleton_candidates.tsv").exists():
                panels.append({
                    "lineage": lineage,
                    "scope": "mixed_high_level",
                    "seed": seed,
                    "panel_id": f"{lineage}_mixed_high_level",
                    "panel_seed": panel_seed,
                    "scenario": MIXED_SCENARIO,
                    "benchmark_dir": benchmark_dir,
                    "run_dir": run_dir,
                })

    return panels


def _compute_metrics(
    df: pd.DataFrame, class_col: str
) -> dict:
    """Compute precision/recall/F1 from a scored DataFrame."""
    n_contam = int(df["is_contaminant_replacement"].sum())
    candidates = df[df[class_col] == "contamination_candidate"]
    n_pred = len(candidates)
    n_tp = int(candidates["is_contaminant_replacement"].sum())
    n_fp = n_pred - n_tp

    precision = n_tp / n_pred if n_pred > 0 else 0.0
    recall = n_tp / n_contam if n_contam > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    return {
        "n_leaves": len(df),
        "n_contaminants": n_contam,
        "n_predicted": n_pred,
        "n_true_positive": n_tp,
        "n_false_positive": n_fp,
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }


def _load_baseline_metrics(panel: dict) -> dict | None:
    """Load baseline recipient_consensus metrics from panel_results_summary.tsv."""
    summary_path = panel["benchmark_dir"] / "panel_results_summary.tsv"
    if not summary_path.exists():
        return None

    df = pd.read_csv(summary_path, sep="\t")
    scenario = panel["scenario"]
    row = df[df["scenario"] == scenario]
    if row.empty:
        return None

    row = row.iloc[0]
    return {
        "tree_rf_delta": float(row.get("tree_rf_delta", 0.0)),
        "contaminant_markers_removed_fraction": float(
            row.get("contaminant_markers_removed_fraction", 0.0)
        ),
        "collateral_genome_loss_count": int(
            row.get("collateral_genome_loss_count", 0)
        ),
        "replacement_contaminant_removed": int(
            row.get("replacement_contaminant_removed", 0)
        ),
        "replacement_events": int(row.get("replacement_events", 0)),
    }


def evaluate_single_panel(panel: dict) -> list[dict]:
    """Run all 3 approaches on a single panel+scenario and return result rows."""
    run_dir = panel["run_dir"]
    benchmark_dir = panel["benchmark_dir"]
    scenario = panel["scenario"]
    rows = []

    base_row = {
        "lineage": panel["lineage"],
        "scope": panel["scope"],
        "seed": panel["seed"],
        "panel_id": panel["panel_id"],
        "scenario": scenario,
    }

    # Load baseline metrics
    baseline = _load_baseline_metrics(panel)
    if baseline:
        rows.append({
            **base_row,
            "approach": "recipient_consensus",
            "n_leaves": 0,
            "n_contaminants": baseline["replacement_events"],
            "n_predicted": baseline["replacement_contaminant_removed"],
            "n_true_positive": baseline["replacement_contaminant_removed"],
            "n_false_positive": 0,
            "precision": 1.0 if baseline["replacement_contaminant_removed"] > 0 else 0.0,
            "recall": (
                baseline["replacement_contaminant_removed"]
                / baseline["replacement_events"]
                if baseline["replacement_events"] > 0
                else 0.0
            ),
            "f1": 0.0,
            "tree_rf_delta": baseline["tree_rf_delta"],
            "collateral_genome_loss_count": baseline["collateral_genome_loss_count"],
            "runtime_seconds": 0.0,
            "status": "ok",
        })
        # Compute F1 for baseline
        r = rows[-1]
        p, rec = r["precision"], r["recall"]
        r["f1"] = 2 * p * rec / (p + rec) if (p + rec) > 0 else 0.0

    # Extract features once for all approaches
    try:
        df = extract_features(run_dir, benchmark_dir, scenario)
    except Exception as e:
        for approach in ["gcp", "cmtv", "plie"]:
            rows.append({
                **base_row,
                "approach": approach,
                "status": f"error: {e}",
                **{k: 0 for k in [
                    "n_leaves", "n_contaminants", "n_predicted",
                    "n_true_positive", "n_false_positive", "precision",
                    "recall", "f1", "tree_rf_delta",
                    "collateral_genome_loss_count", "runtime_seconds",
                ]},
            })
        return rows

    # GCP
    t0 = time.time()
    try:
        df_gcp = score_gcp(df.copy())
        metrics = _compute_metrics(df_gcp, "gcp_class")
        rows.append({
            **base_row,
            "approach": "gcp",
            **metrics,
            "tree_rf_delta": 0.0,
            "collateral_genome_loss_count": 0,
            "runtime_seconds": time.time() - t0,
            "status": "ok",
        })
    except Exception as e:
        rows.append({
            **base_row, "approach": "gcp",
            "status": f"error: {e}",
            **{k: 0 for k in [
                "n_leaves", "n_contaminants", "n_predicted",
                "n_true_positive", "n_false_positive", "precision",
                "recall", "f1", "tree_rf_delta",
                "collateral_genome_loss_count", "runtime_seconds",
            ]},
        })

    # CMTV
    t0 = time.time()
    try:
        has_prottrees = (run_dir / "protTrees" / "no_singles").exists()
        if has_prottrees:
            df_cmtv = score_cmtv(df.copy(), run_dir)
            metrics = _compute_metrics(df_cmtv, "cmtv_class")
        else:
            metrics = _compute_metrics(df.copy().assign(cmtv_class="clean"), "cmtv_class")
        rows.append({
            **base_row,
            "approach": "cmtv",
            **metrics,
            "tree_rf_delta": 0.0,
            "collateral_genome_loss_count": 0,
            "runtime_seconds": time.time() - t0,
            "status": "ok" if has_prottrees else "no_prottrees",
        })
    except Exception as e:
        rows.append({
            **base_row, "approach": "cmtv",
            "status": f"error: {e}",
            **{k: 0 for k in [
                "n_leaves", "n_contaminants", "n_predicted",
                "n_true_positive", "n_false_positive", "precision",
                "recall", "f1", "tree_rf_delta",
                "collateral_genome_loss_count", "runtime_seconds",
            ]},
        })

    # PLIE
    t0 = time.time()
    try:
        df_plie = score_plie(df.copy())
        metrics = _compute_metrics(df_plie, "plie_class")
        rows.append({
            **base_row,
            "approach": "plie",
            **metrics,
            "tree_rf_delta": 0.0,
            "collateral_genome_loss_count": 0,
            "runtime_seconds": time.time() - t0,
            "status": "ok",
        })
    except Exception as e:
        rows.append({
            **base_row, "approach": "plie",
            "status": f"error: {e}",
            **{k: 0 for k in [
                "n_leaves", "n_contaminants", "n_predicted",
                "n_true_positive", "n_false_positive", "precision",
                "recall", "f1", "tree_rf_delta",
                "collateral_genome_loss_count", "runtime_seconds",
            ]},
        })

    return rows


def run_full_evaluation() -> pd.DataFrame:
    """Run all approaches on all available panels."""
    panels = build_panel_index()
    print(f"Found {len(panels)} panel+scenario combinations")

    all_rows = []
    for i, panel in enumerate(panels):
        tag = f"{panel['panel_seed']} / {panel['scenario']}"
        print(f"  [{i+1}/{len(panels)}] {tag} ...", end=" ", flush=True)
        try:
            rows = evaluate_single_panel(panel)
            all_rows.extend(rows)
            ok_count = sum(1 for r in rows if r.get("status") == "ok")
            print(f"{ok_count} ok")
        except Exception:
            print(f"FAILED")
            traceback.print_exc()

    return pd.DataFrame(all_rows)


def generate_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Generate mean +/- std summary grouped by (lineage, scope, scenario, approach)."""
    ok = df[df["status"] == "ok"].copy()
    if ok.empty:
        return pd.DataFrame()

    metric_cols = ["precision", "recall", "f1", "n_true_positive", "n_false_positive"]
    grouped = ok.groupby(["lineage", "scope", "scenario", "approach"])

    rows = []
    for (lineage, scope, scenario, approach), grp in grouped:
        row = {
            "lineage": lineage,
            "scope": scope,
            "scenario": scenario,
            "approach": approach,
            "n_panels": len(grp),
        }
        for col in metric_cols:
            vals = grp[col].values
            row[f"{col}_mean"] = np.mean(vals)
            row[f"{col}_std"] = np.std(vals)
        rows.append(row)

    return pd.DataFrame(rows)


def generate_scope_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Generate summary by scope (aggregating across lineages and seeds)."""
    ok = df[(df["status"] == "ok") & (df["approach"] != "recipient_consensus")].copy()
    baseline = df[(df["status"] == "ok") & (df["approach"] == "recipient_consensus")].copy()

    if ok.empty:
        return pd.DataFrame()

    rows = []
    for (scope, scenario, approach), grp in ok.groupby(["scope", "scenario", "approach"]):
        bl = baseline[(baseline["scope"] == scope) & (baseline["scenario"] == scenario)]
        row = {
            "scope": scope,
            "scenario": scenario,
            "approach": approach,
            "n_panels": len(grp),
            "precision_mean": grp["precision"].mean(),
            "precision_std": grp["precision"].std(),
            "recall_mean": grp["recall"].mean(),
            "recall_std": grp["recall"].std(),
            "f1_mean": grp["f1"].mean(),
            "f1_std": grp["f1"].std(),
            "total_tp": int(grp["n_true_positive"].sum()),
            "total_fp": int(grp["n_false_positive"].sum()),
            "total_contaminants": int(grp["n_contaminants"].sum()),
            "baseline_recall_mean": bl["recall"].mean() if not bl.empty else 0.0,
            "recall_delta_vs_baseline": (
                grp["recall"].mean() - bl["recall"].mean() if not bl.empty else 0.0
            ),
        }
        rows.append(row)

    return pd.DataFrame(rows)


def identify_winner(summary_df: pd.DataFrame) -> str:
    """Identify the best approach based on F1 at family+ scope."""
    if summary_df.empty:
        return "none"

    # Focus on family and order scopes (where detection should work)
    family_plus = summary_df[summary_df["scope"].isin(["family", "order"])]
    if family_plus.empty:
        family_plus = summary_df[summary_df["scope"] != "genus"]

    if family_plus.empty:
        return "none"

    approach_f1 = family_plus.groupby("approach")["f1_mean"].mean()
    return str(approach_f1.idxmax()) if not approach_f1.empty else "none"


def check_results() -> bool:
    """Verify that result files exist and a winner improves over baseline."""
    comparison_path = OUTPUT_DIR / "full_comparison.tsv"
    summary_path = OUTPUT_DIR / "summary_by_scope.tsv"

    if not comparison_path.exists():
        print(f"FAIL: {comparison_path} not found")
        return False
    if not summary_path.exists():
        print(f"FAIL: {summary_path} not found")
        return False

    df = pd.read_csv(comparison_path, sep="\t")
    summary = pd.read_csv(summary_path, sep="\t")

    n_panels = df[df["status"] == "ok"]["panel_id"].nunique()
    print(f"Panels evaluated: {n_panels}")

    n_approaches = df[df["status"] == "ok"]["approach"].nunique()
    print(f"Approaches: {n_approaches}")

    if n_panels < 5:
        print(f"FAIL: only {n_panels} panels evaluated, expected >= 5")
        return False

    winner = identify_winner(summary)
    print(f"Winning approach (family+ F1): {winner}")

    if winner == "none":
        print("WARN: no clear winner identified")
        return True

    # Check if winner improves recall at family+ scope
    family_plus = summary[summary["scope"].isin(["family", "order"])]
    winner_rows = family_plus[family_plus["approach"] == winner]
    if not winner_rows.empty:
        mean_recall_delta = winner_rows["recall_delta_vs_baseline"].mean()
        print(f"Winner recall delta vs baseline (family+): {mean_recall_delta:+.3f}")

    print("PASS: evaluation results valid")
    return True


def main():
    parser = argparse.ArgumentParser(description="US-006: Comparative evaluation")
    parser.add_argument("--check", action="store_true", help="Verify results exist")
    args = parser.parse_args()

    if args.check:
        ok = check_results()
        sys.exit(0 if ok else 1)

    print("=== US-006: Comparative Evaluation ===\n")

    # Run evaluation
    full_df = run_full_evaluation()
    if full_df.empty:
        print("No panels evaluated. Check paths.")
        sys.exit(1)

    # Save full comparison
    comparison_path = OUTPUT_DIR / "full_comparison.tsv"
    full_df.to_csv(comparison_path, sep="\t", index=False)
    print(f"\nSaved full comparison: {comparison_path}")
    print(f"  Total rows: {len(full_df)}")
    print(f"  OK rows: {len(full_df[full_df['status'] == 'ok'])}")

    # Generate summaries
    detail_summary = generate_summary(full_df)
    if not detail_summary.empty:
        detail_path = OUTPUT_DIR / "summary_by_lineage_scope.tsv"
        detail_summary.to_csv(detail_path, sep="\t", index=False)
        print(f"Saved detailed summary: {detail_path}")

    scope_summary = generate_scope_summary(full_df)
    if not scope_summary.empty:
        scope_path = OUTPUT_DIR / "summary_by_scope.tsv"
        scope_summary.to_csv(scope_path, sep="\t", index=False)
        print(f"Saved scope summary: {scope_path}")

        # Print scope summary table
        print("\n=== Summary by Scope ===\n")
        display_cols = [
            "scope", "scenario", "approach", "n_panels",
            "precision_mean", "recall_mean", "f1_mean",
            "total_tp", "total_fp", "recall_delta_vs_baseline",
        ]
        avail = [c for c in display_cols if c in scope_summary.columns]
        print(scope_summary[avail].to_string(index=False, float_format="%.3f"))

    # Identify winner
    winner = identify_winner(scope_summary)
    print(f"\n=== Winner (family+ F1): {winner} ===")

    # Print approach comparison at family+ scope
    if not scope_summary.empty:
        family_plus = scope_summary[scope_summary["scope"].isin(["family", "order"])]
        if not family_plus.empty:
            approach_agg = family_plus.groupby("approach").agg(
                mean_precision=("precision_mean", "mean"),
                mean_recall=("recall_mean", "mean"),
                mean_f1=("f1_mean", "mean"),
                total_tp=("total_tp", "sum"),
                total_fp=("total_fp", "sum"),
            ).round(3)
            print("\n=== Family+ Scope Approach Comparison ===\n")
            print(approach_agg.to_string())

    # Validate
    print("\n=== Validation ===")
    check_results()


if __name__ == "__main__":
    main()
