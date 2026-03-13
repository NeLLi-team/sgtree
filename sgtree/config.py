import os
import glob
import datetime
from dataclasses import dataclass, field


@dataclass
class Config:
    genomedir: str
    modeldir: str
    outdir: str
    num_cpus: int
    percent_models: int
    lflt_fraction: float
    aln_method: str
    tree_method: str
    iqtree_fast: bool
    iqtree_model: str
    hmmsearch_cutoff: str
    hmmsearch_evalue: float
    selection_mode: str
    selection_max_rounds: int
    selection_global_rounds: int
    lock_references: bool
    max_sdup: int
    max_dupl: float
    ref: str | None
    ref_concat: str
    marker_selection: bool
    singles: bool
    singles_mode: str
    num_nei: int
    singles_min_rfdist: float
    keep_intermediates: bool
    is_ref: bool
    start_time: str
    input_format: str = "auto"
    ani_cluster: bool = False
    snp: bool = False
    ani_threshold: float = 95.0
    ani_backend: str = "auto"
    ani_mcl_inflation: float = 2.0
    snp_tree_min_cluster_size: int = 3
    original_genomedir: str | None = None
    original_ref: str | None = None
    model_count: int = 0

    # derived paths (set in __post_init__)
    models_path: str = field(init=False)
    proteomes_path: str = field(init=False)
    staged_proteomes_dir: str = field(init=False)
    gene_call_map_path: str = field(init=False)
    genome_manifest_path: str = field(init=False)
    tables_dir: str = field(init=False)
    extracted_dir: str = field(init=False)
    extracted_seqs_dir: str = field(init=False)
    aligned_dir: str = field(init=False)
    aln_spectree_dir: str = field(init=False)
    trimmed_dir: str = field(init=False)
    concat_dir: str = field(init=False)
    ref_proteomes_path: str = field(init=False)
    ani_dir: str = field(init=False)
    ani_inputs_path: str = field(init=False)
    ani_pairs_path: str = field(init=False)
    ani_cluster_members_path: str = field(init=False)
    ani_representatives_path: str = field(init=False)
    ani_keep_list_path: str = field(init=False)
    ani_selected_query_dir: str = field(init=False)
    ani_selected_ref_dir: str = field(init=False)
    snp_trees_dir: str = field(init=False)
    snp_tree_summary_path: str = field(init=False)

    def __post_init__(self):
        self.models_path = os.path.join(self.outdir, "models")
        self.proteomes_path = os.path.join(self.outdir, "proteomes")
        self.staged_proteomes_dir = os.path.join(self.outdir, "staged_proteomes")
        self.gene_call_map_path = os.path.join(self.outdir, "gene_calls.tsv")
        self.genome_manifest_path = os.path.join(self.outdir, "genome_manifest.tsv")
        self.tables_dir = os.path.join(self.outdir, "tables")
        self.extracted_dir = os.path.join(self.outdir, "extracted")
        self.extracted_seqs_dir = os.path.join(self.outdir, "extracted_seqs")
        self.aligned_dir = os.path.join(self.outdir, "aligned")
        self.aln_spectree_dir = os.path.join(self.outdir, "aln_SpecTree")
        self.trimmed_dir = os.path.join(self.outdir, "trimmed_SpeciesTree")
        self.concat_dir = os.path.join(self.outdir, "concat")
        self.ref_proteomes_path = os.path.join(self.outdir, "ref_and_query_proteomes")
        self.ani_dir = os.path.join(self.outdir, "ani")
        self.ani_inputs_path = os.path.join(self.ani_dir, "inputs.tsv")
        self.ani_pairs_path = os.path.join(self.ani_dir, "ani_pairwise.tsv")
        self.ani_cluster_members_path = os.path.join(self.ani_dir, "ani_clusters.tsv")
        self.ani_representatives_path = os.path.join(self.ani_dir, "ani_representatives.tsv")
        self.ani_keep_list_path = os.path.join(self.ani_dir, "ani_kept_genomes.txt")
        self.ani_selected_query_dir = os.path.join(self.ani_dir, "query_representatives")
        self.ani_selected_ref_dir = os.path.join(self.ani_dir, "ref_representatives")
        self.snp_trees_dir = os.path.join(self.outdir, "snp_trees")
        self.snp_tree_summary_path = os.path.join(self.snp_trees_dir, "snp_tree_summary.tsv")
        if self.original_genomedir is None:
            self.original_genomedir = self.genomedir
        if self.original_ref is None:
            self.original_ref = self.ref

    @property
    def hitsoutdir(self):
        return os.path.join(self.outdir, "hits.hmmout")

    @property
    def min_models_fraction(self):
        return self.percent_models / 100

    @property
    def genome_count(self):
        if os.path.isdir(self.genomedir):
            return len(glob.glob(os.path.join(self.genomedir, "*")))
        if os.path.isfile(self.genomedir):
            return 1
        return 0

    @property
    def model_file_count(self):
        if os.path.isfile(self.modeldir):
            count = 0
            with open(self.modeldir, "rb") as handle:
                for line in handle:
                    if line.startswith(b"NAME"):
                        count += 1
            return count
        return len(glob.glob(os.path.join(self.modeldir, "*.hmm")))

    def ref_dir_path(self):
        """Path to the reference concat directory for this ref+model combination."""
        if self.ref is None:
            return None
        ref_name = self.ref.rstrip("/").split("/")[-1]
        model_name = os.path.basename(self.modeldir.rstrip("/"))
        if model_name.endswith(".hmm"):
            model_name = model_name[:-4]
        return os.path.join(self.ref_concat, f"{ref_name}_{model_name}")

    def print_banner(self):
        sep = "=" * 80
        print(f"{self.outdir}\n{sep}")
        print(f"Sg_Tree v.2\nstart time: {self.start_time}\n{sep}")
        print(f"Genomes database {self.genomedir} contains {self.genome_count} genomes\n{sep}")
        print(f"Marker database {self.modeldir} contains {self.model_file_count} models\n{sep}\n")
