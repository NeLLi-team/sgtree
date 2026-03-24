{
  "version": 2,
  "project": "ML Contamination Detection",
  "overview": "Improve singleton/replacement contamination detection in sgtree by developing and testing 3+ creative ML/statistical approaches that leverage unsupervised anomaly detection on per-leaf marker features. Current approaches (delta_rf, topoknn, hybrid, composite, recipient_consensus, neighbor_clade, neighbor_ml) either lack sensitivity or cause excessive collateral damage. New approaches learn the 'normal' marker distribution from the majority of clean markers without requiring explicit reference/query labels, then flag outliers as contamination candidates. All testing uses the existing 50-genome (chlam/flavo/gamma) and 100-genome (alpha/bacteroidota) benchmark panels. The task is done when at least one approach strongly improves both contaminant detection rate and collateral damage compared to the current best (recipient_consensus).",
  "goals": [
    "Develop 3+ ML-based contamination detection approaches leveraging unsupervised anomaly detection",
    "Achieve higher contaminant_markers_removed_fraction with same or lower collateral_genome_loss_count vs recipient_consensus baseline",
    "Improve tree_rf_delta (final tree closer to truth) across all scenarios: duplicate_only, replacement_only, combined, mixed_high_level",
    "Evaluate on both 50-genome and 100-genome benchmark panels for generalizability",
    "Integrate the best-performing approach into the CLI as a new --singles-mode"
  ],
  "nonGoals": [
    "Changing the duplicate resolution pipeline (RF-distance coordinate descent) -- only singleton detection is in scope",
    "Requiring explicit --ref flags or labeled reference/query genome labels at runtime",
    "Deep learning or GPU-dependent approaches (keep to sklearn + lightweight libs)",
    "Rewriting the benchmark generation or evaluation infrastructure",
    "Achieving perfect detection at genus-level scope (known to be fundamentally hard due to sequence similarity)"
  ],
  "successMetrics": [
    "At least one approach achieves >80% replacement_contaminant_removed on family+ scope (vs ~70-90% for recipient_consensus)",
    "Collateral genome loss stays at 0 for >90% of panels (fewer false positive removals)",
    "Mean tree_rf_delta improvement of >10% over recipient_consensus across combined/replacement_only scenarios",
    "Approach generalizes across lineages (chlam/flavo/gamma/alpha/bacteroidota) without per-lineage tuning"
  ],
  "openQuestions": [
    "Should cross-marker topology voting use full RF distance between pruned subtrees, or a lighter proxy like shared-neighbor count?",
    "What is the minimum number of markers per genome needed for genome-consistency profiling to work reliably?",
    "Can pseudo-labeling iteration converge in a reasonable number of rounds (2-3) without overfitting?"
  ],
  "stack": {
    "language": "Python 3.10+",
    "environment": "pixi (conda-based, for trimal/VeryFastTree/hmmer)",
    "ml_libraries": "scikit-learn (existing), xgboost, umap-learn, hdbscan",
    "testing": "unittest (existing pattern)",
    "data": "pandas, numpy (existing)"
  },
  "routes": [],
  "uiNotes": [],
  "dataModel": [
    {
      "entity": "LeafFeatureMatrix",
      "fields": ["panel_id", "scenario", "genome", "marker_name", "leaf_name", "contig_id", "delta_rf", "neighbor_overlap", "topoknn_score", "branch_outlier", "bitscore_outlier", "recipient_consensus_score", "species_anchor_score", "purity_drop", "anchor_knn_agreement", "attachment_gap", "is_contaminant (truth label)"]
    },
    {
      "entity": "ApproachResult",
      "fields": ["panel_id", "scenario", "approach_name", "contaminant_markers_removed_fraction", "collateral_genome_loss_count", "tree_rf_delta", "singleton_precision", "singleton_recall", "runtime_seconds"]
    }
  ],
  "importFormat": {
    "description": "Existing benchmark panel structure with events.tsv ground truth and SGTree output directories",
    "example": {
      "benchmark_dir": "eval/full/50gen/chlam/chlam_family__seed_101",
      "events_tsv": "scenarios/replacement_only/events.tsv",
      "run_dir": "runs/replacement_only__singles_recipient_consensus"
    }
  },
  "rules": [
    "All approaches must be unsupervised at runtime -- no explicit reference/query genome labels required",
    "The majority of markers in any real dataset are clean; approaches exploit this distributional assumption",
    "Per-genome budget constraints apply: max 1 prune per genome, must retain min_markers_per_genome",
    "RF guard gate remains: proposed removals must not worsen marker-tree RF distance",
    "Evaluation uses the existing evaluate_benchmark_run() infrastructure for consistent comparison",
    "Feature extraction must handle variable panel sizes (20-100 genomes, 10-56 markers)",
    "New pixi dependencies (xgboost, umap-learn, hdbscan) added to pixi.toml"
  ],
  "qualityGates": [
    "pixi run python -m unittest tests.test_benchmark -v",
    "pixi run python runs/ml_contam_detection/evaluate_approaches.py --check"
  ],
  "stories": [
    {
      "id": "US-001",
      "title": "Establish baseline results for all existing singleton modes",
      "status": "done",
      "dependsOn": [],
      "description": "As a researcher, I want baseline results for all 7 existing singleton modes on both 50-genome and 100-genome benchmark panels so that I have a clear comparison target for new approaches.",
      "acceptanceCriteria": [
        "Run all 7 singleton modes (delta_rf, topoknn, hybrid, composite, recipient_consensus, neighbor_clade, neighbor_ml) on at least 6 representative panels: 3 from 50gen (one each of chlam/flavo/gamma, family scope, seed 101) and 3 from 100gen (alpha + bacteroidota, family scope, seed 101)",
        "For each run, collect: contaminant_markers_removed_fraction, collateral_genome_loss_count, tree_rf_delta, singleton_precision, singleton_recall",
        "Run across all applicable scenarios: duplicate_only, replacement_only, combined (and mixed_high_level for 50gen panels that have it)",
        "Save results to runs/ml_contam_detection/baseline_results.tsv with columns: panel_id, scenario, singles_mode, and all metrics",
        "Example: chlam_family_seed101 + replacement_only + recipient_consensus should show the current detection rate",
        "Negative case: duplicate_only scenario should show 0 singleton activity (no replacements to detect)",
        "Reuse existing evaluated results from eval/full/50gen/ and eval/full/100gen/ where available instead of re-running"
      ]
    },
    {
      "id": "US-002",
      "title": "Build unified feature extraction framework",
      "status": "done",
      "dependsOn": ["US-001"],
      "description": "As a developer, I want a feature extraction module that extracts per-leaf candidate metrics from SGTree marker selection outputs into a unified DataFrame, so that all ML approaches can share the same feature matrix.",
      "acceptanceCriteria": [
        "Create runs/ml_contam_detection/feature_extraction.py with a function extract_features(run_dir, benchmark_dir, scenario_name) -> pd.DataFrame",
        "Output DataFrame has one row per (genome, marker, leaf) with columns: all _score_singleton_candidates() metrics (delta_rf, neighbor_overlap, topoknn_score, branch_outlier, bitscore_outlier, recipient_consensus_score, species_anchor_score, purity_drop, anchor_knn_agreement, attachment_gap, etc.)",
        "Add derived cross-marker features: genome_median_overlap, genome_marker_count, marker_genome_count, genome_consistency_score (MAD of neighbor_overlap across markers for same genome)",
        "Add ground-truth labels from events.tsv: is_contaminant_replacement, is_contaminant_duplicate, is_native",
        "Feature matrix can be extracted from existing benchmark run directories (reads singleton_proposals or marker_selection intermediates)",
        "Example: extract_features() on chlam_family_seed101 replacement_only produces a DataFrame with ~500 rows (50 genomes x 10 markers) with 4 rows marked as replacement contaminants",
        "Negative case: if run_dir has no singleton proposal data, function raises a clear error with instructions",
        "Add pixi dependencies: xgboost, umap-learn, hdbscan to pixi.toml"
      ]
    },
    {
      "id": "US-003",
      "title": "Approach 1: Genome Consistency Profiling (GCP)",
      "status": "done",
      "dependsOn": ["US-002"],
      "description": "As a researcher, I want a genome-consistency-based anomaly detector that scores each marker leaf by how inconsistent it is with the genome's other marker placements, so that replacement contaminants (which disrupt a genome's topological profile) are detected without explicit labels.",
      "acceptanceCriteria": [
        "Create runs/ml_contam_detection/approach_gcp.py implementing Genome Consistency Profiling",
        "Core idea: for each genome, collect feature vectors across all its markers. Clean markers form a consistent cluster; contaminant markers deviate. Score = Mahalanobis distance of each marker from its genome's robust center",
        "Use MinCovDet (robust covariance) per genome when marker count >= 5, fall back to median absolute deviation for smaller panels",
        "Features used: neighbor_overlap, topoknn_score, branch_outlier, bitscore_outlier, recipient_consensus_score, species_anchor_score, purity_drop",
        "Apply HDBSCAN on per-genome feature matrices as a secondary signal: markers not in the main cluster get a penalty",
        "Output: scored proposals with gcp_score (higher = more likely contamination) and gcp_class (contamination_candidate / clean / ambiguous)",
        "Run on 6 baseline panels from US-001, compute precision/recall/F1 for contamination detection",
        "Example: a replacement contaminant in chlam_family should have high gcp_score because its topology features differ from the genome's other 9 clean markers",
        "Negative case: genomes with only 1-2 markers cannot compute consistency -- these should be classified as 'ambiguous', not false-positived"
      ]
    },
    {
      "id": "US-004",
      "title": "Approach 2: Cross-Marker Topology Voting (CMTV)",
      "status": "done",
      "dependsOn": ["US-002"],
      "description": "As a researcher, I want a cross-marker voting approach where each marker tree independently 'votes' on whether a genome's placement is consistent, so that a contaminated marker is outvoted by the genome's clean markers across other trees.",
      "acceptanceCriteria": [
        "Create runs/ml_contam_detection/approach_cmtv.py implementing Cross-Marker Topology Voting",
        "Core idea: for each genome G and candidate marker M, check G's neighborhood in marker M's tree vs G's neighborhoods in all other marker trees. If M places G near different taxa than the consensus of other markers, M is suspicious",
        "For each (genome, marker) pair, compute a 'topology agreement score': fraction of other marker trees where the genome's k-nearest neighbors overlap with this marker tree's neighbors (reuse _leaf_neighbor_overlap logic adapted for cross-tree comparison)",
        "Aggregate cross-marker votes: if a marker's placement disagrees with >60% of other markers for the same genome, flag it",
        "Weight votes by the quality of each voter tree (trees with low mean RF get higher weight)",
        "Output: scored proposals with cmtv_score and cmtv_class",
        "Run on 6 baseline panels, compute precision/recall/F1",
        "Example: in combined scenario, a replacement marker should show low topology agreement because it was placed by a donor genome's sequence, not the recipient's",
        "Negative case: HGT events (real biological transfers) might also show disagreement -- the approach should use a stricter threshold to avoid flagging these"
      ]
    },
    {
      "id": "US-005",
      "title": "Approach 3: Pseudo-Label Iterative Ensemble (PLIE)",
      "status": "done",
      "dependsOn": ["US-002"],
      "description": "As a researcher, I want an iterative pseudo-labeling approach that bootstraps contamination labels from high-confidence outlier detection, then trains a supervised XGBoost classifier on these soft labels, so that we get the power of supervised learning without requiring actual labels at runtime.",
      "acceptanceCriteria": [
        "Create runs/ml_contam_detection/approach_plie.py implementing Pseudo-Label Iterative Ensemble",
        "Round 0: Run robust outlier detection (IsolationForest + LOF + robust Mahalanobis) on all leaf features. High-confidence inliers (bottom 80% of anomaly scores) become pseudo-clean. High-confidence outliers (top 5%) become pseudo-contaminated",
        "Round 1+: Train XGBoost binary classifier on pseudo-labels. Use genome-level and marker-level aggregated features in addition to raw features. Predict contamination probability for all leaves. Update pseudo-labels from high-confidence predictions",
        "Iterate for 2-3 rounds or until labels stabilize (< 1% change)",
        "Features: all raw features from US-002 + genome_consistency_score + per-genome z-scores of each feature + per-marker z-scores",
        "Final output: scored proposals with plie_score (contamination probability) and plie_class",
        "Run on 6 baseline panels, compute precision/recall/F1",
        "Example: the pseudo-labeling should correctly identify replacement contaminants at family+ scope even without true labels, because they are genuine statistical outliers",
        "Negative case: at genus scope, contaminants are similar to clean markers -- approach should assign low confidence (score < 0.5) rather than confidently misclassifying"
      ]
    },
    {
      "id": "US-006",
      "title": "Comparative evaluation on full benchmark suite",
      "status": "done",
      "dependsOn": ["US-003", "US-004", "US-005"],
      "description": "As a researcher, I want a comprehensive comparative evaluation of all 3 new approaches against all existing modes across the full 50-genome and 100-genome benchmark suites, so that I can identify the best-performing approach with statistical confidence.",
      "acceptanceCriteria": [
        "Create runs/ml_contam_detection/evaluate_approaches.py that runs all 3 new approaches + existing baseline modes across the full benchmark matrix",
        "50-genome panels: 3 lineages x 4 scopes x 3 seeds = 36 panels, each with 3-4 scenarios",
        "100-genome panels: 2 lineages x 4 scopes x 3 seeds = 24 panels, each with 3 scenarios x 3 contig variants",
        "For each (panel, scenario, approach), collect: contaminant_markers_removed_fraction, collateral_genome_loss_count, tree_rf_delta, singleton_precision, singleton_recall, F1, runtime_seconds",
        "Generate summary tables: mean +/- std across seeds for each (lineage, scope, scenario, approach)",
        "Identify the winning approach: highest mean tree_rf_delta with collateral_genome_loss_count <= baseline",
        "Save results to runs/ml_contam_detection/full_comparison.tsv and runs/ml_contam_detection/summary_by_scope.tsv",
        "Add --check flag to evaluate_approaches.py that verifies the result files exist and the winner improves over baseline",
        "Example: at family scope combined scenario, the winning approach should show measurably better RF delta than recipient_consensus",
        "Negative case: if no approach beats baseline on a specific scope (e.g., genus), document this clearly in the summary rather than cherry-picking results"
      ]
    },
    {
      "id": "US-007",
      "title": "Integrate winning approach into CLI as new singles_mode",
      "status": "open",
      "dependsOn": ["US-006"],
      "description": "As a user, I want to run the best-performing contamination detection approach via sgtree --singles-mode <new_mode> so that I can use it in my phylogenomic workflows.",
      "acceptanceCriteria": [
        "Add the winning approach as a new singleton mode in sgtree/marker_selection.py",
        "The mode is selectable via --singles-mode in the CLI (sgtree/cli.py)",
        "Implementation reuses the existing _score_singleton_candidates() feature extraction pipeline",
        "New ML scoring is applied as a post-processing step on the standard candidate features",
        "Any new pixi dependencies are already in pixi.toml from US-002",
        "The mode works on arbitrary panel sizes (not just benchmark panels)",
        "Add the new mode to SINGLETON_MODES list and DEFAULT_CLEANUP_PROFILES",
        "Example: pixi run python -m sgtree --genomedir proteomes/ --modeldir models/ --singles --singles-mode <new_mode> produces tree_final.nwk with contamination removed",
        "Negative case: if the panel has < 3 genomes or < 3 markers, fall back to recipient_consensus mode with a warning"
      ]
    },
    {
      "id": "US-008",
      "title": "Final validation and documentation",
      "status": "open",
      "dependsOn": ["US-007"],
      "description": "As a researcher, I want the integrated mode validated on the full benchmark suite through the standard SGTree pipeline (not standalone scripts), confirming that end-to-end performance matches the standalone evaluation.",
      "acceptanceCriteria": [
        "Run the integrated mode through run_baseline.py (modified to include the new mode) on at least 6 representative panels",
        "Confirm tree_rf_delta matches standalone evaluation results within 0.01 tolerance",
        "Confirm collateral_genome_loss_count matches standalone results exactly",
        "Update eval/full/50gen/ and eval/full/100gen/ panel_results_summary.tsv files with the new mode's results alongside existing recipient_consensus results",
        "All quality gates pass: unittest suite and evaluate_approaches.py --check",
        "Example: chlam_family_seed101 combined scenario should produce identical results whether run via standalone script or via CLI",
        "Negative case: if any panel shows regression vs standalone evaluation, investigate and fix the integration before marking complete"
      ]
    }
  ]
}
