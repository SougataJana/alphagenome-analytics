import unittest

from pydantic import ValidationError

from app.models.variants import VariantRequest
from app.services.evidence import EvidenceService
from app.services.statistics import compare_groups, gene_enrichment


class VariantValidationTests(unittest.TestCase):
    def test_normalizes_chromosome_and_alleles(self):
        variant = VariantRequest(
            chromosome=" 7 ", position=140_453_136, reference="a", alternate="g"
        )
        self.assertEqual(variant.chromosome, "chr7")
        self.assertEqual(variant.reference, "A")
        self.assertEqual(variant.alternate, "G")

    def test_rejects_same_allele_and_out_of_bounds_position(self):
        with self.assertRaises(ValidationError):
            VariantRequest(chromosome="chr7", position=1, reference="A", alternate="a")
        with self.assertRaises(ValidationError):
            VariantRequest(
                chromosome="chrM", position=16_570, reference="A", alternate="C"
            )

    def test_rejects_unknown_chromosome_and_non_snv(self):
        with self.assertRaises(ValidationError):
            VariantRequest(chromosome="chrUn", position=1, reference="A", alternate="C")
        with self.assertRaises(ValidationError):
            VariantRequest(chromosome="chr1", position=1, reference="AA", alternate="C")


class StatisticsTests(unittest.TestCase):
    def test_seed_reproduces_resampling_statistics(self):
        first = compare_groups([1, 2, 3], [4, 5, 6], 500, 7)
        second = compare_groups([1, 2, 3], [4, 5, 6], 500, 7)
        for key in ("case_mean", "control_mean", "mean_difference", "empirical_p_value", "confidence_interval_95"):
            self.assertEqual(first[key], second[key])
        self.assertEqual(first["mean_difference"], -3.0)
        self.assertTrue(0 < first["empirical_p_value"] <= 1)

    def test_compare_requires_independent_samples_and_finite_values(self):
        with self.assertRaises(ValueError):
            compare_groups([1], [2, 3], 100, 0)
        with self.assertRaises(ValueError):
            compare_groups([1, float("nan")], [2, 3], 100, 0)

    def test_enrichment_rejects_selected_genes_outside_background(self):
        with self.assertRaisesRegex(ValueError, "occur in the stated background"):
            gene_enrichment(["A", "outside"], ["A", "B"], {"set": ["A"]})

    def test_enrichment_reports_background_and_adjusted_probability(self):
        result = gene_enrichment(
            ["a", "b"], ["a", "b", "c", "d"],
            {"contains_a": ["a", "c"], "contains_none": ["d"]},
        )
        by_name = {row["gene_set"]: row for row in result["results"]}
        self.assertEqual(result["selected_genes_in_background"], 2)
        self.assertEqual(by_name["contains_a"]["overlap_n"], 1)
        self.assertAlmostEqual(by_name["contains_a"]["p_value"], 5 / 6)
        self.assertGreaterEqual(by_name["contains_a"]["q_value_bh"], by_name["contains_a"]["p_value"])


class RegionValidationTests(unittest.TestCase):
    def setUp(self):
        self.service = EvidenceService()

    def test_rejects_invalid_region_chromosome_and_grch38_overflow(self):
        with self.assertRaisesRegex(ValueError, "Chromosome must"):
            self.service.region("chrUn", 1, 100)
        with self.assertRaisesRegex(ValueError, "chromosome length"):
            self.service.encode_region("chrM", 16_000, 16_570)

    def test_rejects_zero_based_or_reversed_interval(self):
        with self.assertRaisesRegex(ValueError, "positive"):
            self.service.encode_region("chr8", 0, 100)
        with self.assertRaisesRegex(ValueError, "positive"):
            self.service.encode_region("chr8", 200, 100)


if __name__ == "__main__":
    unittest.main()
