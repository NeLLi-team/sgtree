import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from sgtree import ani


class AniTests(unittest.TestCase):
    def test_python_mcl_clusters_thresholded_pairs(self):
        labels = ["A", "B", "C"]
        rows = [
            {"genome_a": "A", "genome_b": "B", "ani": 0.99},
            {"genome_a": "A", "genome_b": "C", "ani": 0.80},
            {"genome_a": "B", "genome_b": "C", "ani": 0.79},
        ]

        clusters = ani._run_python_mcl(labels, rows, ani_threshold=0.95, inflation=2.0)

        self.assertIn(["A", "B"], clusters)
        self.assertIn(["C"], clusters)

    def test_choose_cluster_representative_prefers_reference_then_contiguity(self):
        cluster = [
            ani.GenomeRecord("QueryA", "query", "fna", "QueryA.fna", "QueryA.fna", None, 1000, 10),
            ani.GenomeRecord("RefA", "ref", "fna", "RefA.fna", "RefA.fna", None, 900, 20),
            ani.GenomeRecord("RefB", "ref", "fna", "RefB.fna", "RefB.fna", None, 800, 5),
        ]

        representative = ani.choose_cluster_representative(cluster)

        self.assertEqual(representative.genome_id, "RefB")

    def test_build_snp_trees_writes_alignment_and_summary(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            query_manifest = tmp / "query_manifest.tsv"
            clusters = tmp / "clusters.tsv"
            outdir = tmp / "snp_trees"

            query_manifest.write_text(
                "\n".join(
                    [
                        "genome_id\tinput_format\tsource_file\tassembly_path\tstaged_proteome_path\tcontigs\ttotal_bases",
                        f"Rep\tfna\t{tmp / 'Rep.fna'}\t{tmp / 'Rep.fna'}\t\t1\t4",
                        f"Q1\tfna\t{tmp / 'Q1.fna'}\t{tmp / 'Q1.fna'}\t\t1\t4",
                        f"Q2\tfna\t{tmp / 'Q2.fna'}\t{tmp / 'Q2.fna'}\t\t1\t4",
                    ]
                )
                + "\n"
            )
            clusters.write_text(
                "\n".join(
                    [
                        "cluster_id\trepresentative_genome\tgenome_id\tsource_role\tcluster_size\ttotal_bases\tcontigs\tkept_for_species_tree\tbackend\tani_threshold\tmcl_inflation",
                        "ani_cluster_001\tRep\tRep\tquery\t3\t4\t1\tyes\tminimap2\t0.9500\t2.00",
                        "ani_cluster_001\tRep\tQ1\tquery\t3\t4\t1\tno\tminimap2\t0.9500\t2.00",
                        "ani_cluster_001\tRep\tQ2\tquery\t3\t4\t1\tno\tminimap2\t0.9500\t2.00",
                    ]
                )
                + "\n"
            )
            for genome_id, sequence in [("Rep", "AAAA"), ("Q1", "AATA"), ("Q2", "AACA")]:
                (tmp / f"{genome_id}.fna").write_text(f">contig1\n{sequence}\n")

            projections = {
                "Q1.fna": ani.SamProjection(coverage={"contig1": [(0, 4)]}, mismatches={("contig1", 2): "T"}),
                "Q2.fna": ani.SamProjection(coverage={"contig1": [(0, 4)]}, mismatches={("contig1", 2): "C"}),
            }

            def fake_projection(reference_path: str, query_path: str):
                return projections[Path(query_path).name]

            def fake_fasttree(alignment_path: Path, tree_path: Path):
                tree_path.write_text("(Rep,Q1,Q2);\n")

            with patch("sgtree.ani._project_alignment_to_reference", side_effect=fake_projection), patch(
                "sgtree.ani._run_fasttree_nt",
                side_effect=fake_fasttree,
            ):
                summary_path = ani.build_snp_trees(
                    clusters_path=clusters,
                    query_manifest=query_manifest,
                    ref_manifest=None,
                    outdir=outdir,
                    min_cluster_size=3,
                    cpus=1,
                )

            summary = pd.read_csv(summary_path, sep="\t")
            self.assertEqual(summary.loc[0, "status"], "built")
            self.assertEqual(int(summary.loc[0, "snp_sites"]), 1)
            alignment = (outdir / "ani_cluster_001" / "core_snps.fna").read_text()
            self.assertIn(">Rep\nA\n", alignment)
            self.assertIn(">Q1\nT\n", alignment)
            self.assertIn(">Q2\nC\n", alignment)


if __name__ == "__main__":
    unittest.main()
