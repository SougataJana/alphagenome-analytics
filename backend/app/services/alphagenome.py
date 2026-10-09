import logging
import json
import math
import os
from functools import lru_cache
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from dotenv import load_dotenv

from app.models.variants import (
    FeatureAttribution,
    ModalityEffect,
    VariantAnalysisResult,
    VariantRequest,
)
from app.services.reference import ReferenceGenomeService, ReferenceValidationError

logger = logging.getLogger(__name__)
BACKEND_ROOT = Path(__file__).resolve().parents[2]
load_dotenv(BACKEND_ROOT / ".env")
REQUESTED_SCORERS = [
    "AVI_SCORE",
    "AVI_SCORE_FEATURE_IMPORTANCE",
    "RNA_SEQ",
    "SPLICE_SITES",
    "ATAC",
]
MAX_RETURNED_EFFECTS = 200


def _serialize_atlas_matrix(matrix) -> dict:
    def frame_records(frame):
        # Pandas performs JSON-safe conversion for NumPy scalars and missing values.
        return json.loads(frame.to_json(orient="records", double_precision=15))

    def dense_values(values):
        if hasattr(values, "toarray"):
            values = values.toarray()
        if hasattr(values, "tolist"):
            values = values.tolist()

        def json_safe(value):
            if isinstance(value, list):
                return [json_safe(item) for item in value]
            if isinstance(value, float) and not math.isfinite(value):
                return None
            return value

        return json_safe(values)

    layers = {}
    for name in matrix.layers.keys():
        if isinstance(name, str):
            layers[name] = dense_values(matrix.layers[name])
    return {
        "shape": list(matrix.shape),
        "X": dense_values(matrix.X),
        "obs": frame_records(matrix.obs),
        "var": frame_records(matrix.var),
        "layers": layers,
    }


def serialize_atlas_sdk_response(results) -> dict:
    """Preserve the SDK-decoded scorer matrices and metadata without truncation."""
    return {
        scorer: _serialize_atlas_matrix(matrix)
        for scorer, matrix in results.items()
        if matrix is not None
    }


class AlphaGenomeNotConfigured(Exception):
    pass


class AlphaGenomeRequestError(Exception):
    def __init__(self, message: str, status_code: int):
        super().__init__(message)
        self.status_code = status_code


@lru_cache(maxsize=1)
def _atlas_client():
    api_key = os.getenv("ALPHAGENOME_API_KEY")
    if not api_key:
        raise AlphaGenomeNotConfigured(
            "Set ALPHAGENOME_API_KEY in backend/.env to enable AlphaGenome Atlas."
        )

    from alphagenome.atlas import atlas

    return atlas.create(api_key, timeout=30)


class AlphaGenomeService:
    """Fetch normalized single-SNV scores from the official AlphaGenome Atlas API."""

    def __init__(self, reference_service: ReferenceGenomeService | None = None):
        self.reference_service = reference_service or ReferenceGenomeService()

    def analyze(self, variant: VariantRequest) -> VariantAnalysisResult:
        try:
            reference_validation = self.reference_service.verify(variant)
        except ReferenceValidationError as exc:
            raise AlphaGenomeRequestError(str(exc), exc.status_code) from exc

        from alphagenome.data import genome

        client = _atlas_client()
        try:
            results = client.query_variant(
                genome.Variant(
                    chromosome=variant.chromosome,
                    position=variant.position,
                    reference_bases=variant.reference,
                    alternate_bases=variant.alternate,
                ),
                requested_scorers=REQUESTED_SCORERS,
            )
        except PermissionError as exc:
            logger.warning("AlphaGenome Atlas rejected the configured credentials or access.")
            raise AlphaGenomeRequestError(
                "AlphaGenome Atlas rejected the configured key or account access.", 502
            ) from exc
        except TimeoutError as exc:
            raise AlphaGenomeRequestError("AlphaGenome Atlas request timed out.", 504) from exc
        except (ValueError, IndexError) as exc:
            raise AlphaGenomeRequestError(
                "AlphaGenome Atlas could not find scores for this variant. Confirm the "
                "GRCh38 coordinate and reference/alternate alleles.",
                404,
            ) from exc
        except Exception as exc:
            logger.exception("AlphaGenome Atlas request failed.")
            raise AlphaGenomeRequestError("AlphaGenome Atlas request failed.", 502) from exc

        atlas_sdk_response = serialize_atlas_sdk_response(results)

        avi = results.get("AVI_SCORE")
        if avi is None or avi.X.shape[0] == 0 or avi.X.shape[1] == 0:
            raise AlphaGenomeRequestError(
                "AlphaGenome Atlas returned no AVI score for this variant.", 404
            )

        avi_score = float(avi.X[0, 0])
        if not math.isfinite(avi_score):
            raise AlphaGenomeRequestError("AlphaGenome Atlas returned a non-finite AVI score.", 502)
        quantiles = avi.layers.get("quantiles")
        avi_quantile = float(quantiles[0, 0]) if quantiles is not None else None
        if avi_quantile is not None and (
            not math.isfinite(avi_quantile) or not 0 <= avi_quantile <= 1
        ):
            logger.warning("AlphaGenome Atlas returned an invalid AVI quantile; omitting it.")
            avi_quantile = None

        feature_attributions: list[FeatureAttribution] = []
        feature_scores = results.get("AVI_SCORE_FEATURE_IMPORTANCE")
        if feature_scores is not None and feature_scores.X.shape[0] > 0:
            feature_names = (
                feature_scores.var["name"].tolist()
                if "name" in feature_scores.var
                else [str(index) for index in range(feature_scores.X.shape[1])]
            )
            for index, name in enumerate(feature_names):
                value = float(feature_scores.X[0, index])
                if math.isfinite(value):
                    feature_attributions.append(FeatureAttribution(feature=str(name), value=value))

        effects: list[ModalityEffect] = []
        for scorer_name in ("RNA_SEQ", "SPLICE_SITES", "ATAC"):
            matrix = results.get(scorer_name)
            if matrix is None or matrix.X.shape[0] == 0 or matrix.X.shape[1] == 0:
                continue
            quantile_matrix = matrix.layers.get("quantiles")
            for row_index in range(matrix.X.shape[0]):
                obs = matrix.obs.iloc[row_index]
                for column_index in range(matrix.X.shape[1]):
                    var = matrix.var.iloc[column_index]
                    value = matrix.X[row_index, column_index]
                    if hasattr(value, "toarray"):
                        value = value.toarray()[0, 0]
                    raw_score = float(value)
                    if not (raw_score == raw_score and abs(raw_score) != float("inf")):
                        continue
                    feature = var.get("name")
                    def clean_label(*values):
                        for candidate in values:
                            if candidate is None:
                                continue
                            try:
                                if candidate != candidate:
                                    continue
                            except (TypeError, ValueError):
                                continue
                            text = str(candidate).strip()
                            if text and text.lower() not in {"nan", "none"}:
                                return text
                        return None

                    gene = clean_label(obs.get("gene_name"), obs.get("gene_id"))
                    tissue = clean_label(var.get("biosample_name"), var.get("gtex_tissue"), var.get("ontology_curie"))
                    quantile = (
                        float(quantile_matrix[row_index, column_index])
                        if quantile_matrix is not None
                        else None
                    )
                    if quantile is not None and (
                        not math.isfinite(quantile) or not -1 <= quantile <= 1
                    ):
                        logger.warning("AlphaGenome Atlas returned an invalid effect quantile; omitting it.")
                        quantile = None
                    effects.append(
                        ModalityEffect(
                            modality=scorer_name,
                            feature=str(feature) if feature is not None else None,
                            gene=str(gene) if gene is not None else None,
                            tissue=str(tissue) if tissue is not None else None,
                            raw_score=raw_score,
                            quantile_score=quantile,
                        )
                    )
        total_effects = len(effects)
        effects_by_scorer = {
            scorer: sorted(
                (effect for effect in effects if effect.modality == scorer),
                key=lambda effect: abs(effect.raw_score),
                reverse=True,
            )
            for scorer in ("RNA_SEQ", "SPLICE_SITES", "ATAC")
        }
        effects_total_by_scorer = {scorer: len(items) for scorer, items in effects_by_scorer.items()}
        # Cap each scorer independently: absolute raw values have different scales
        # across scorers and must not compete for a shared top-N cutoff.
        effects = [
            effect
            for scorer in ("RNA_SEQ", "SPLICE_SITES", "ATAC")
            for effect in effects_by_scorer[scorer][:MAX_RETURNED_EFFECTS]
        ]

        return VariantAnalysisResult(
            status="complete",
            provider="AlphaGenome Atlas",
            provider_sdk_version=self._sdk_version(),
            provider_model_version=None,
            variant=variant,
            reference_validation=reference_validation,
            avi_score=avi_score,
            avi_quantile=avi_quantile,
            requested_scorers=REQUESTED_SCORERS.copy(),
            feature_importance=feature_attributions,
            effects=effects,
            effects_total=total_effects,
            effects_total_by_scorer=effects_total_by_scorer,
            effects_truncated=any(
                count > MAX_RETURNED_EFFECTS for count in effects_total_by_scorer.values()
            ),
            provider_response_available=True,
            atlas_sdk_response=atlas_sdk_response,
        )

    @staticmethod
    def _sdk_version() -> str | None:
        try:
            return version("alphagenome")
        except PackageNotFoundError:
            return None
