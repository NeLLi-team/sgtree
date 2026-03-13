import unittest

import pandas as pd

from sgtree.benchmark_dataset import (
    DEFAULT_BURKHOLDERIACEAE_SPECIES_PLAN,
    select_species_rows,
)


class BenchmarkDatasetSelectionTests(unittest.TestCase):
    def test_default_plan_has_expected_shape(self):
        self.assertEqual(len(DEFAULT_BURKHOLDERIACEAE_SPECIES_PLAN), 20)
        self.assertEqual(
            sum(int(row["target_count"]) for row in DEFAULT_BURKHOLDERIACEAE_SPECIES_PLAN),
            50,
        )
        self.assertEqual(
            len({str(row["genus"]) for row in DEFAULT_BURKHOLDERIACEAE_SPECIES_PLAN}),
            20,
        )

    def test_select_species_rows_prefers_quality_and_respects_counts(self):
        candidates = pd.DataFrame(
            [
                {
                    "assembly_accession": "GCA_000001111.1",
                    "organism_name": "Achromobacter xylosoxidans strain A",
                    "ftp_path": "https://example.org/A",
                    "genome_size_bp": 6500000,
                    "assembly_level": "Scaffold",
                    "refseq_category": "na",
                    "phylum": "Pseudomonadota",
                    "class": "Gammaproteobacteria",
                    "order_name": "Burkholderiales",
                    "family": "Burkholderiaceae",
                    "genus": "Achromobacter",
                    "species": "Achromobacter xylosoxidans",
                    "taxonomy_source": "gtdb",
                },
                {
                    "assembly_accession": "GCA_000001112.1",
                    "organism_name": "Achromobacter xylosoxidans strain B",
                    "ftp_path": "https://example.org/B",
                    "genome_size_bp": 6499000,
                    "assembly_level": "Complete Genome",
                    "refseq_category": "reference genome",
                    "phylum": "Pseudomonadota",
                    "class": "Gammaproteobacteria",
                    "order_name": "Burkholderiales",
                    "family": "Burkholderiaceae",
                    "genus": "Achromobacter",
                    "species": "Achromobacter xylosoxidans",
                    "taxonomy_source": "gtdb",
                },
                {
                    "assembly_accession": "GCA_000001211.1",
                    "organism_name": "Bordetella pertussis strain A",
                    "ftp_path": "https://example.org/C",
                    "genome_size_bp": 4106000,
                    "assembly_level": "Complete Genome",
                    "refseq_category": "na",
                    "phylum": "Pseudomonadota",
                    "class": "Gammaproteobacteria",
                    "order_name": "Burkholderiales",
                    "family": "Burkholderiaceae",
                    "genus": "Bordetella",
                    "species": "Bordetella pertussis",
                    "taxonomy_source": "gtdb",
                },
                {
                    "assembly_accession": "GCA_000001212.1",
                    "organism_name": "Bordetella pertussis strain B",
                    "ftp_path": "https://example.org/D",
                    "genome_size_bp": 4106100,
                    "assembly_level": "Complete Genome",
                    "refseq_category": "reference genome",
                    "phylum": "Pseudomonadota",
                    "class": "Gammaproteobacteria",
                    "order_name": "Burkholderiales",
                    "family": "Burkholderiaceae",
                    "genus": "Bordetella",
                    "species": "Bordetella pertussis",
                    "taxonomy_source": "gtdb",
                },
                {
                    "assembly_accession": "GCA_000001213.1",
                    "organism_name": "Bordetella pertussis strain C",
                    "ftp_path": "https://example.org/E",
                    "genome_size_bp": 4106200,
                    "assembly_level": "Scaffold",
                    "refseq_category": "na",
                    "phylum": "Pseudomonadota",
                    "class": "Gammaproteobacteria",
                    "order_name": "Burkholderiales",
                    "family": "Burkholderiaceae",
                    "genus": "Bordetella",
                    "species": "Bordetella pertussis",
                    "taxonomy_source": "gtdb",
                },
            ]
        )
        species_plan = [
            {
                "genus": "Achromobacter",
                "species": "Achromobacter xylosoxidans",
                "target_count": 1,
                "role": "singleton",
            },
            {
                "genus": "Bordetella",
                "species": "Bordetella pertussis",
                "target_count": 2,
                "role": "strain_set_2",
            },
        ]

        selection = select_species_rows(candidates, species_plan=species_plan, prefix="TEST__")

        self.assertEqual(len(selection), 3)
        self.assertEqual(
            selection.groupby(["genus", "species"]).size().to_dict(),
            {
                ("Achromobacter", "Achromobacter xylosoxidans"): 1,
                ("Bordetella", "Bordetella pertussis"): 2,
            },
        )
        achromobacter_row = selection[selection["genus"] == "Achromobacter"].iloc[0]
        self.assertEqual(achromobacter_row["assembly_accession"], "GCA_000001112.1")
        self.assertEqual(achromobacter_row["genome_id"], "TEST__GCA-000001112-1")
        self.assertIn("Burkholderiaceae", achromobacter_row["taxonomy_lookup"])


if __name__ == "__main__":
    unittest.main()
