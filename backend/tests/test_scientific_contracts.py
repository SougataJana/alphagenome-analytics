import unittest
import tempfile
import sys
import types
import os
from pathlib import Path
from unittest.mock import Mock, patch

import httpx
import numpy as np
import pandas as pd

from pydantic import ValidationError

from app.models.variants import ModalityEffect, VariantRequest
from app.services.evidence import EvidenceProviderError, EvidenceService
from app.services.alphagenome import AlphaGenomeService, _serialize_atlas_matrix
from app.services.reference import ReferenceGenomeService, ReferenceValidationError
from app.services.statistics import _bh, compare_groups, gene_enrichment
from app.services import store


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

    def test_accepts_first_and_last_grch38_positions(self):
        first = VariantRequest(chromosome="chr1", position=1, reference="A", alternate="C")
        last = VariantRequest(chromosome="chr1", position=248_956_422, reference="A", alternate="C")
        self.assertEqual(first.position, 1)
        self.assertEqual(last.position, 248_956_422)


class ReferenceAlleleTests(unittest.TestCase):
    def variant(self, chromosome="chr7", position=140_453_136, reference="A", alternate="G"):
        return VariantRequest(chromosome=chromosome, position=position, reference=reference, alternate=alternate)

    def service(self, status_code=200, text="A"):
        response = httpx.Response(status_code, text=text, request=httpx.Request("GET", "https://rest.ensembl.org"))
        client = Mock()
        client.get.return_value = response
        return ReferenceGenomeService(client=client), client

    def test_verifies_matching_ref_and_records_provider_provenance(self):
        service, client = self.service(text="a\n")
        result = service.verify(self.variant())
        self.assertEqual(result["status"], "verified")
        self.assertEqual(result["assembly"], "GRCh38")
        self.assertEqual(result["observed_reference"], "A")
        self.assertEqual(result["source"], "Ensembl REST sequence/region")
        self.assertEqual(result["source_url"], "https://rest.ensembl.org/sequence/region/human/7:140453136..140453136:1?coord_system_version=GRCh38")
        self.assertEqual(client.get.call_args.kwargs["params"], {"coord_system_version": "GRCh38"})
        client.get.assert_called_once()

    def test_rejects_mismatch_without_correcting_input(self):
        service, _ = self.service(text="C")
        with self.assertRaisesRegex(ReferenceValidationError, "submitted REF=A, Ensembl reports C") as error:
            service.verify(self.variant())
        self.assertEqual(error.exception.status_code, 422)

    def test_invalid_provider_response_fails_closed(self):
        service, _ = self.service(text="AC")
        with self.assertRaisesRegex(ReferenceValidationError, "invalid single-base"):
            service.verify(self.variant())

    def test_ensembl_failure_does_not_silently_skip_verification(self):
        service, _ = self.service(status_code=503, text="unavailable")
        with self.assertRaises(ReferenceValidationError) as error:
            service.verify(self.variant())
        self.assertEqual(error.exception.status_code, 502)

    def test_mitochondrial_alias_maps_to_ensembl_mt_sequence(self):
        service, client = self.service(text="A")
        service.verify(self.variant(chromosome="chrM", position=1))
        self.assertIn("/MT:1..1:1", client.get.call_args.args[0])

    def test_alpha_query_is_not_called_when_reference_check_fails(self):
        reference = Mock()
        reference.verify.side_effect = ReferenceValidationError("REF mismatch", 422)
        service = AlphaGenomeService(reference_service=reference)
        with patch("app.services.alphagenome._atlas_client") as atlas_client:
            from app.services.alphagenome import AlphaGenomeRequestError
            with self.assertRaises(AlphaGenomeRequestError) as error:
                service.analyze(self.variant())
        self.assertEqual(error.exception.status_code, 422)
        atlas_client.assert_not_called()

    def test_effect_quantiles_preserve_signed_direction(self):
        effect = ModalityEffect(
            modality="RNA_SEQ", raw_score=-1.36, quantile_score=-0.9999998
        )
        self.assertAlmostEqual(effect.quantile_score, -0.9999998)
        with self.assertRaises(ValidationError):
            ModalityEffect(modality="RNA_SEQ", raw_score=1.0, quantile_score=-1.01)


class AtlasResponsePreservationTests(unittest.TestCase):
    def test_matrix_serializer_keeps_metadata_quantiles_and_json_safe_values(self):
        class Frame:
            def to_json(self, **kwargs):
                return '[{"label":"biosample"}]'

        class Matrix:
            shape = (1, 2)
            X = [[1.25, float("nan")]]
            obs = Frame()
            var = Frame()
            layers = {"quantiles": [[-0.98, 0.4]], None: [[123]]}

        result = _serialize_atlas_matrix(Matrix())
        self.assertEqual(result["shape"], [1, 2])
        self.assertEqual(result["X"], [[1.25, None]])
        self.assertEqual(result["obs"], [{"label": "biosample"}])
        self.assertEqual(result["layers"], {"quantiles": [[-0.98, 0.4]]})

    def test_compressed_provider_response_round_trip(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(store, "DB_PATH", Path(directory) / "test.sqlite3"):
                payload = {"RNA_SEQ": {"X": [[-1.36]], "layers": {"quantiles": [[-0.99]]}}}
                self.assertTrue(store.save_provider_response("analysis-1", "AlphaGenome Atlas", payload))
                restored = store.get_provider_response("analysis-1")
        self.assertEqual(restored["response"], payload)
        self.assertIn("not raw HTTP bytes", restored["representation"])

    def test_atlas_adapter_preserves_signed_quantile_and_full_decoded_matrices(self):
        class Matrix:
            def __init__(self, x, obs, var, layers=None):
                self.X = np.asarray(x, dtype=float)
                self.obs = pd.DataFrame(obs)
                self.var = pd.DataFrame(var)
                self.layers = {key: np.asarray(value, dtype=float) for key, value in (layers or {}).items()}
                self.shape = self.X.shape

        response = {
            "AVI_SCORE": Matrix([[1.5548]], {"variant": ["v1"]}, {"name": ["AVI"]}, {"quantiles": [[0.99894]]}),
            "AVI_SCORE_FEATURE_IMPORTANCE": Matrix([[0.7, -0.2]], {"variant": ["v1"]}, {"name": ["RNA_SEQ", "ATAC"]}),
            "RNA_SEQ": Matrix(
                [[-1.36, 0.2]],
                {"gene_name": ["TP53"], "gene_id": ["ENSG00000141510"], "strand": [1]},
                {"name": ["track-a", "track-b"], "biosample_name": ["lung", "liver"], "gtex_tissue": [None, None]},
                {"quantiles": [[-0.99998, 0.15]]},
            ),
        }
        client = Mock()
        client.query_variant.return_value = response
        parent = types.ModuleType("alphagenome")
        data = types.ModuleType("alphagenome.data")
        genome = types.ModuleType("alphagenome.data.genome")
        genome.Variant = lambda **kwargs: kwargs
        parent.data = data
        data.genome = genome
        variant = VariantRequest(chromosome="chr22", position=36_201_698, reference="A", alternate="C")
        service = AlphaGenomeService(reference_service=Mock(verify=Mock(return_value={
            "status": "verified", "assembly": "GRCh38", "source": "Ensembl REST sequence/region",
            "source_url": "https://rest.ensembl.org/sequence/region/human/22:36201698..36201698:1",
            "requested_reference": "A", "observed_reference": "A", "checked_at": "2026-10-09T00:00:00+00:00",
        })))
        with patch("app.services.alphagenome._atlas_client", return_value=client), patch.dict(
            sys.modules,
            {"alphagenome": parent, "alphagenome.data": data, "alphagenome.data.genome": genome},
        ):
            result = service.analyze(variant)

        self.assertEqual(result.avi_score, 1.5548)
        self.assertAlmostEqual(result.avi_quantile, 0.99894)
        self.assertEqual(result.reference_validation.observed_reference, "A")
        self.assertEqual(result.effects[0].quantile_score, -0.99998)
        self.assertEqual(result.effects[0].gene, "TP53")
        self.assertIsNone(result.provider_model_version)
        self.assertEqual(result.atlas_sdk_response["RNA_SEQ"]["shape"], [1, 2])
        self.assertEqual(result.atlas_sdk_response["RNA_SEQ"]["layers"]["quantiles"], [[-0.99998, 0.15]])


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

    def test_permutation_p_value_agrees_with_exact_small_sample_reference(self):
        result = compare_groups([1, 2], [3, 4], 10_000, 4)
        # Six possible allocations of two observations to the case group; two
        # are at least as extreme as the observed absolute mean difference.
        self.assertAlmostEqual(result["empirical_p_value"], 2 / 6, delta=0.025)

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

    def test_bh_matches_known_reference_vector(self):
        observed = _bh([0.01, 0.04, 0.03, 0.002])
        expected = [0.02, 0.04, 0.04, 0.008]
        np.testing.assert_allclose(observed, expected, rtol=0, atol=1e-12)

    def test_enrichment_zero_overlap_has_unit_p_value(self):
        result = gene_enrichment(["A"], ["A", "B", "C"], {"other": ["B"]})
        row = result["results"][0]
        self.assertEqual(row["overlap_n"], 0)
        self.assertEqual(row["p_value"], 1.0)


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


class ProviderAdapterTests(unittest.TestCase):
    def setUp(self):
        self.service = EvidenceService()

    def test_ensembl_region_adapter_normalizes_annotation_and_provenance(self):
        response = httpx.Response(
            200,
            json=[{"id": "ENSG00000141510", "feature_type": "gene", "external_name": "TP53", "start": 100, "end": 200}],
            request=httpx.Request("GET", "https://rest.ensembl.org/overlap/region"),
        )
        with patch.object(self.service, "_get", return_value=response):
            result = self.service.region("chr1", 100, 200)
        self.assertEqual(result["region"]["assembly"], "GRCh38")
        self.assertEqual(result["features"][0]["symbol"], "TP53")
        self.assertEqual(result["features"][0]["source"], "Ensembl")
        self.assertTrue(result["source_url"].startswith("https://rest.ensembl.org/"))

    def test_gnomad_adapter_keeps_cohorts_separate_and_computes_frequency(self):
        response = httpx.Response(
            200,
            json={"data": {"variant": {"variantId": "22-36201698-A-C", "genome": {"ac": 10, "an": 1000}, "exome": {"ac": 2, "an": 200}}}},
            request=httpx.Request("POST", "https://gnomad.broadinstitute.org/api"),
        )
        with patch("app.services.evidence.httpx.post", return_value=response):
            result = self.service.gnomad_variant("chr22", 36_201_698, "A", "C")
        self.assertEqual(result["dataset"], "gnomad_r4 (GRCh38)")
        self.assertEqual(result["result"]["genome"]["af_calculated"], 0.01)
        self.assertEqual(result["result"]["exome"]["af_calculated"], 0.01)
        self.assertIn("kept separate", result["note"])

    def test_gwas_adapter_preserves_records_and_source(self):
        response = httpx.Response(
            200,
            json={"_embedded": {"associations": [{"id": "association-1", "pvalue": 1e-8}]}},
            request=httpx.Request("GET", "https://www.ebi.ac.uk/gwas/rest/api/v2/associations"),
        )
        with patch.object(self.service, "_get", return_value=response):
            result = self.service.gwas_associations("rs123")
        self.assertEqual(result["associations"][0]["id"], "association-1")
        self.assertEqual(result["source"], "GWAS Catalog API v2")

    def test_clinvar_adapter_returns_empty_records_without_inventing_assertions(self):
        search = httpx.Response(
            200, json={"esearchresult": {"idlist": []}},
            request=httpx.Request("GET", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"),
        )
        with patch.object(self.service, "_get", return_value=search):
            result = self.service.clinvar("rs123")
        self.assertEqual(result["records"], [])
        self.assertIsNone(result["summary_url"])
        self.assertEqual(result["source"], "NCBI ClinVar E-utilities")

    def test_gtex_adapter_preserves_provider_payload_and_dataset(self):
        response = httpx.Response(
            200, json={"associations": [{"tissueSiteDetailId": "Lung", "pValue": 1e-5}]},
            request=httpx.Request("GET", "https://gtexportal.org/api/v2/association/singleTissueEqtl"),
        )
        with patch.object(self.service, "_get", return_value=response):
            result = self.service.gtex_eqtls("rs123")
        self.assertEqual(result["dataset"], "gtex_v10")
        self.assertEqual(result["associations"]["associations"][0]["tissueSiteDetailId"], "Lung")
        self.assertEqual(result["source"], "GTEx Portal API v2")

    def test_screen_adapter_separates_cres_from_prediction_scores(self):
        response = httpx.Response(
            200,
            json={"data": {"cCRESCREENSearch": [{"chrom": "chr8", "start": 100, "len": 200, "info": {"accession": "EH38E1"}}]}},
            request=httpx.Request("POST", "https://screen.api.wenglab.org/graphql"),
        )
        with patch.dict(os.environ, {"SCREEN_API_KEY": "test-key"}), patch(
            "app.services.evidence.httpx.post", return_value=response
        ):
            result = self.service.encode_region("chr8", 100, 200)
        self.assertEqual(result["ccres"][0]["info"]["accession"], "EH38E1")
        self.assertEqual(result["source"], "SCREEN GraphQL API (ENCODE Registry of cCREs)")
        self.assertNotIn("avi_score", result)

    def test_provider_failure_is_not_returned_as_an_empty_biological_result(self):
        with patch.object(self.service, "region", side_effect=EvidenceProviderError("Ensembl timed out")):
            result = self.service.synthesize_variant("chr1", 100, "A", "C")
        provider = result["providers"]["Ensembl region"]
        self.assertEqual(provider["status"], "unavailable")
        self.assertIn("timed out", provider["error"])
        self.assertNotIn("result", provider)


class VariantToGeneScoringTests(unittest.TestCase):
    def setUp(self):
        self.service = EvidenceService()

    def test_effect_strength_can_outrank_row_count(self):
        features = [
            {"id": "ENSG00000000001", "feature_type": "gene", "symbol": "GENEA", "start": 1500, "end": 1600},
            {"id": "ENSG00000000002", "feature_type": "gene", "symbol": "GENEB", "start": 1900, "end": 1950},
        ]
        effects = [
            {"gene": "GENEA", "modality": "RNA_SEQ", "quantile_score": 0.1}
            for _ in range(12)
        ] + [{"gene": "GENEB", "modality": "RNA_SEQ", "quantile_score": -0.9}]
        with patch.object(self.service, "region", return_value={"features": features, "source_url": "source"}):
            result = self.service.prioritize_variant_genes("chr1", 1000, effects, window=1000)
        self.assertEqual(result["candidates"][0]["gene"], "GENEB")
        by_gene = {row["gene"]: row for row in result["candidates"]}
        self.assertGreater(by_gene["GENEA"]["effect_count"], by_gene["GENEB"]["effect_count"])
        self.assertLess(by_gene["GENEA"]["research_priority_score"], by_gene["GENEB"]["research_priority_score"])
        self.assertEqual(by_gene["GENEB"]["best_effect_quantile"], -0.9)
        self.assertEqual(by_gene["GENEB"]["alphagenome_support_component"], 0.9)

    def test_unmapped_gene_with_missing_distance_is_unscored(self):
        features = [
            {"id": "ENSG00000000003", "feature_type": "gene", "symbol": "NEARBY", "start": 1000, "end": 1010},
        ]
        effects = [{"gene": "UNMAPPED", "modality": "RNA_SEQ", "quantile_score": 1.0}]
        with patch.object(self.service, "region", return_value={"features": features, "source_url": "source"}):
            result = self.service.prioritize_variant_genes("chr1", 1005, effects, window=1000)
        by_gene = {row["gene"]: row for row in result["candidates"]}
        self.assertEqual(by_gene["NEARBY"]["research_priority_score"], 0.5)
        self.assertIsNone(by_gene["UNMAPPED"]["research_priority_score"])
        self.assertIsNone(by_gene["UNMAPPED"]["rank"])

    def test_symbol_and_versioned_ensembl_id_merge_to_one_candidate(self):
        features = [
            {"id": "ENSG00000141510", "feature_type": "gene", "symbol": "TP53", "start": 900, "end": 1100},
        ]
        effects = [
            {"gene": "tp53", "modality": "RNA_SEQ", "quantile_score": 0.2},
            {"gene": "ENSG00000141510.4", "modality": "RNA_SEQ", "quantile_score": 0.7},
        ]
        with patch.object(self.service, "region", return_value={"features": features, "source_url": "source"}):
            result = self.service.prioritize_variant_genes("chr1", 1000, effects, window=1000)
        self.assertEqual(len(result["candidates"]), 1)
        candidate = result["candidates"][0]
        self.assertEqual(candidate["gene"], "TP53")
        self.assertEqual(candidate["effect_count"], 2)
        self.assertEqual(candidate["best_effect_quantile"], 0.7)

    def test_ambiguous_symbol_does_not_merge_effects_to_arbitrary_gene(self):
        features = [
            {"id": "ENSG00000000011", "feature_type": "gene", "symbol": "DUP", "start": 900, "end": 950},
            {"id": "ENSG00000000012", "feature_type": "gene", "symbol": "DUP", "start": 1050, "end": 1100},
        ]
        effects = [{"gene": "DUP", "modality": "RNA_SEQ", "quantile_score": 0.9}]
        with patch.object(self.service, "region", return_value={"features": features, "source_url": "source"}):
            result = self.service.prioritize_variant_genes("chr1", 1000, effects, window=1000)
        genes = [row for row in result["candidates"] if row["gene"] == "DUP"]
        unmapped = [row for row in result["candidates"] if row["gene"] == "DUP" and row["ensembl_id"] is None]
        self.assertEqual(len(genes), 3)
        self.assertEqual(len(unmapped), 1)
        self.assertIsNone(unmapped[0]["research_priority_score"])
        self.assertIsNone(unmapped[0]["rank"])


if __name__ == "__main__":
    unittest.main()
