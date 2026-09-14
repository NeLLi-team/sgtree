"""Store run settings and derive the paths used by pipeline stages."""

from dataclasses import dataclass, field
from pathlib import Path

from sgtree._version import DISPLAY_VERSION


@dataclass
class Config:
    """Hold validated CLI options and the derived paths for one run."""

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
    ani_cluster_members_path: str = field(init=False)
    ani_keep_list_path: str = field(init=False)
    ani_selected_query_dir: str = field(init=False)
    ani_selected_ref_dir: str = field(init=False)
    snp_trees_dir: str = field(init=False)

    def __post_init__(self) -> None:
        outdir = Path(self.outdir)
        self.models_path = str(outdir / "models")
        self.proteomes_path = str(outdir / "proteomes")
        self.staged_proteomes_dir = str(outdir / "staged_proteomes")
        self.gene_call_map_path = str(outdir / "gene_calls.tsv")
        self.genome_manifest_path = str(outdir / "genome_manifest.tsv")
        self.tables_dir = str(outdir / "tables")
        self.extracted_dir = str(outdir / "extracted")
        self.extracted_seqs_dir = str(outdir / "extracted_seqs")
        self.aligned_dir = str(outdir / "aligned")
        self.aln_spectree_dir = str(outdir / "aln_SpecTree")
        self.trimmed_dir = str(outdir / "trimmed_SpeciesTree")
        self.concat_dir = str(outdir / "concat")
        self.ref_proteomes_path = str(outdir / "ref_and_query_proteomes")
        ani_dir = outdir / "ani"
        self.ani_dir = str(ani_dir)
        self.ani_cluster_members_path = str(ani_dir / "ani_clusters.tsv")
        self.ani_keep_list_path = str(ani_dir / "ani_kept_genomes.txt")
        self.ani_selected_query_dir = str(ani_dir / "query_representatives")
        self.ani_selected_ref_dir = str(ani_dir / "ref_representatives")
        self.snp_trees_dir = str(outdir / "snp_trees")
        if self.original_genomedir is None:
            self.original_genomedir = self.genomedir
        if self.original_ref is None:
            self.original_ref = self.ref

    @property
    def hitsoutdir(self) -> str:
        """Return the unfiltered HMM hit table path."""
        return str(Path(self.outdir) / "hits.hmmout")

    @property
    def min_models_fraction(self) -> float:
        """Return the minimum marker completeness as a fraction."""
        return self.percent_models / 100

    @property
    def genome_count(self) -> int:
        """Count visible files in the genome input or staged proteome directory."""
        directory = Path(self.genomedir)
        if directory.is_dir():
            return sum(
                path.is_file() and not path.name.startswith(".")
                for path in directory.iterdir()
            )
        if directory.is_file():
            return 1
        return 0

    @property
    def model_file_count(self) -> int:
        """Count models in a concatenated HMM file or model directory."""
        models = Path(self.modeldir)
        if models.is_file():
            with models.open("rb") as handle:
                return sum(line.startswith(b"NAME") for line in handle)
        return sum(path.is_file() for path in models.glob("*.hmm"))

    def ref_dir_path(self) -> str | None:
        """Return the reference cache path for this reference and model pair."""
        if self.ref is None:
            return None
        ref_name = Path(self.ref).name or "root"
        model_name = Path(self.modeldir).name or "root"
        model_name = model_name.removesuffix(".hmm")
        return str(Path(self.ref_concat) / f"{ref_name}_{model_name}")

    def print_banner(self) -> None:
        """Print the input summary used by the command-line runner."""
        sep = "=" * 80
        print(f"{self.outdir}\n{sep}")
        print(f"{DISPLAY_VERSION}\nstart time: {self.start_time}\n{sep}")
        print(
            f"Genomes database {self.genomedir} contains "
            f"{self.genome_count} genomes\n{sep}"
        )
        print(
            f"Marker database {self.modeldir} contains "
            f"{self.model_file_count} models\n{sep}\n"
        )
