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
    ref: str | None
    ref_concat: str
    marker_selection: bool
    singles: bool
    is_ref: bool
    start_time: str
    model_count: int = 0

    # derived paths (set in __post_init__)
    models_path: str = field(init=False)
    proteomes_path: str = field(init=False)
    tables_dir: str = field(init=False)
    extracted_dir: str = field(init=False)
    extracted_seqs_dir: str = field(init=False)
    aligned_dir: str = field(init=False)
    aln_spectree_dir: str = field(init=False)
    trimmed_dir: str = field(init=False)
    concat_dir: str = field(init=False)
    ref_proteomes_path: str = field(init=False)

    def __post_init__(self):
        self.models_path = os.path.join(self.outdir, "models")
        self.proteomes_path = os.path.join(self.outdir, "proteomes")
        self.tables_dir = os.path.join(self.outdir, "tables")
        self.extracted_dir = os.path.join(self.outdir, "extracted")
        self.extracted_seqs_dir = os.path.join(self.outdir, "extracted_seqs")
        self.aligned_dir = os.path.join(self.outdir, "aligned")
        self.aln_spectree_dir = os.path.join(self.outdir, "aln_SpecTree")
        self.trimmed_dir = os.path.join(self.outdir, "trimmed_SpeciesTree")
        self.concat_dir = os.path.join(self.outdir, "concat")
        self.ref_proteomes_path = os.path.join(self.outdir, "ref_and_query_proteomes")

    @property
    def hitsoutdir(self):
        return os.path.join(self.outdir, "hits.hmmout")

    @property
    def min_models_fraction(self):
        return self.percent_models / 100

    @property
    def genome_count(self):
        return len(glob.glob(os.path.join(self.genomedir, "*")))

    @property
    def model_file_count(self):
        return len(glob.glob(os.path.join(self.modeldir, "*")))

    def ref_dir_path(self):
        """Path to the reference concat directory for this ref+model combination."""
        if self.ref is None:
            return None
        ref_name = self.ref.rstrip("/").split("/")[-1]
        model_name = self.modeldir.rstrip("/").split("/")[-1]
        return os.path.join(self.ref_concat, f"{ref_name}_{model_name}")

    def print_banner(self):
        sep = "=" * 80
        print(f"{self.outdir}\n{sep}")
        print(f"Sg_Tree v.2\nstart time: {self.start_time}\n{sep}")
        print(f"Genomes database {self.genomedir} contains {self.genome_count} genomes\n{sep}")
        print(f"Marker database {self.modeldir} contains {self.model_file_count} models\n{sep}\n")
