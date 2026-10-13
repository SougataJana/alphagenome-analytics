"""Optional, on-demand AlphaGenome REF/ALT signal tracks for visualization.

These predictions are deliberately separate from precomputed Atlas lookups.
"""

import logging
import math
import os
from functools import lru_cache
from importlib.metadata import PackageNotFoundError, version

from app.models.variants import VariantRequest

logger = logging.getLogger(__name__)

SUPPORTED_TRACK_MODALITIES = (
    "RNA_SEQ",
    "CAGE",
    "ATAC",
    "DNASE",
    "CHIP_HISTONE",
    "CHIP_TF",
    "SPLICE_SITES",
    "SPLICE_SITE_USAGE",
)
MAX_TRACKS = 8
MAX_PLOT_POINTS = 1200
PLOT_WINDOW_BP = 50_000


class TrackPredictionError(Exception):
    """A useful, non-secret-bearing error from the optional track request."""


@lru_cache(maxsize=1)
def _model_client():
    api_key = os.getenv("ALPHAGENOME_API_KEY")
    if not api_key:
        raise TrackPredictionError(
            "Set ALPHAGENOME_API_KEY on the backend to enable detailed predictions."
        )
    from alphagenome.models import dna_client

    return dna_client.create(api_key)


def _get_sdk_version() -> str | None:
    try:
        return version("alphagenome")
    except PackageNotFoundError:
        return None


def _plain(value):
    if value is None:
        return None
    try:
        if value != value:
            return None
    except (TypeError, ValueError):
        pass
    return str(value)


@lru_cache(maxsize=1)
def get_track_options() -> dict:
    """Return tissue/assay options from the SDK's current output metadata."""
    try:
        from alphagenome.models import dna_client

        metadata = _model_client().output_metadata(
            organism=dna_client.Organism.HOMO_SAPIENS
        )
        options: dict[str, list[dict[str, str]]] = {}
        for modality in SUPPORTED_TRACK_MODALITIES:
            frame = getattr(metadata, modality.lower(), None)
            if frame is None or frame.empty:
                options[modality] = []
                continue
            columns = set(frame.columns)
            rows = []
            seen = set()
            for _, row in frame.iterrows():
                ontology_term = _plain(row.get("ontology_curie")) if "ontology_curie" in columns else None
                biosample = _plain(row.get("biosample_name")) if "biosample_name" in columns else None
                if not ontology_term:
                    continue
                label = biosample or ontology_term
                key = (ontology_term, label)
                if key in seen:
                    continue
                seen.add(key)
                rows.append({"ontology_term": ontology_term, "label": label})
            options[modality] = sorted(rows, key=lambda item: item["label"].casefold())
        return {"provider": "AlphaGenome SDK output metadata", "options": options}
    except TrackPredictionError:
        raise
    except Exception as exc:
        logger.exception("Could not retrieve AlphaGenome output metadata.")
        raise TrackPredictionError(
            "AlphaGenome could not provide track metadata right now. Retry shortly."
        ) from exc


def _downsample_pair(reference, alternate, track_data, variant_position: int) -> dict:
    """Crop a TrackData pair to 50 kb and mean-bin to a bounded SVG payload."""
    import numpy as np

    ref_values = np.asarray(reference.values)
    alt_values = np.asarray(alternate.values)
    if ref_values.ndim != 2 or alt_values.ndim != 2 or ref_values.shape != alt_values.shape:
        raise TrackPredictionError("AlphaGenome returned an unsupported track matrix shape.")

    interval = reference.interval
    resolution = max(1, int(reference.resolution))
    left = max(int(interval.start), variant_position - 1 - PLOT_WINDOW_BP // 2)
    right = min(int(interval.end), variant_position + PLOT_WINDOW_BP // 2)
    index_start = max(0, (left - int(interval.start)) // resolution)
    index_end = min(ref_values.shape[0], math.ceil((right - int(interval.start)) / resolution))
    ref_values = ref_values[index_start:index_end]
    alt_values = alt_values[index_start:index_end]
    if ref_values.shape[0] == 0:
        raise TrackPredictionError("AlphaGenome returned no signal around the selected variant.")

    bins = min(MAX_PLOT_POINTS, ref_values.shape[0])
    edges = np.linspace(0, ref_values.shape[0], bins + 1, dtype=int)

    def reduce(values, column):
        series = []
        for start, end in zip(edges[:-1], edges[1:]):
            segment = values[start:end, column]
            segment = segment[np.isfinite(segment)]
            series.append(float(segment.mean()) if segment.size else None)
        return series

    metadata = track_data.metadata.reset_index(drop=True)
    result_tracks = []
    count = min(ref_values.shape[1], MAX_TRACKS)
    for column in range(count):
        row = metadata.iloc[column] if column < len(metadata) else {}
        row_dict = row.to_dict() if hasattr(row, "to_dict") else {}
        result_tracks.append(
            {
                "name": _plain(row_dict.get("name")) or f"Track {column + 1}",
                "biosample": _plain(row_dict.get("biosample_name")),
                "strand": _plain(row_dict.get("strand")),
                "resolution_bp": resolution,
                "reference": reduce(ref_values, column),
                "alternate": reduce(alt_values, column),
            }
        )
    return {
        "start_1based": left + 1,
        "end_1based_inclusive": min(right, left + ref_values.shape[0] * resolution),
        "variant_position": variant_position,
        "downsampling": "mean signal per bin",
        "tracks_truncated": ref_values.shape[1] > MAX_TRACKS,
        "tracks": result_tracks,
    }


def predict_variant_tracks(
    variant: VariantRequest, ontology_term: str, modality: str
) -> dict:
    """Compute one selected on-demand modality for one biosample ontology term."""
    if modality not in SUPPORTED_TRACK_MODALITIES:
        raise TrackPredictionError("Choose a supported AlphaGenome track modality.")
    if not ontology_term or len(ontology_term) > 80:
        raise TrackPredictionError("Choose a tissue or biosample ontology term.")
    try:
        from alphagenome.data import genome
        from alphagenome.models import dna_client

        model = _model_client()
        output_type = getattr(dna_client.OutputType, modality)
        model_interval = genome.Variant(
            chromosome=variant.chromosome,
            position=variant.position,
            reference_bases=variant.reference,
            alternate_bases=variant.alternate,
        ).reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
        output = model.predict_variant(
            interval=model_interval,
            variant=genome.Variant(
                chromosome=variant.chromosome,
                position=variant.position,
                reference_bases=variant.reference,
                alternate_bases=variant.alternate,
            ),
            requested_outputs=[output_type],
            ontology_terms=[ontology_term],
        )
        ref = getattr(output.reference, modality.lower(), None)
        alt = getattr(output.alternate, modality.lower(), None)
        if ref is None or alt is None or ref.values.size == 0 or alt.values.size == 0:
            raise TrackPredictionError(
                "No tracks were returned for this modality and tissue. Try another available biosample."
            )
        plot = _downsample_pair(ref, alt, ref, variant.position)
        if not plot["tracks"]:
            raise TrackPredictionError(
                "No tracks were returned for this modality and tissue. Try another available biosample."
            )
        return {
            "provider": "AlphaGenome API · on-demand predict_variant",
            "provider_sdk_version": _get_sdk_version(),
            "model_version": "Not exposed by the API response",
            "variant": variant.model_dump(mode="json"),
            "ontology_term": ontology_term,
            "modality": modality,
            "assembly": "GRCh38",
            "output_kind": "on-demand REF/ALT model prediction; separate from precomputed Atlas scores",
            **plot,
        }
    except TrackPredictionError:
        raise
    except TimeoutError as exc:
        raise TrackPredictionError(
            "Detailed AlphaGenome prediction timed out. Atlas summary results remain available; try again later."
        ) from exc
    except Exception as exc:
        logger.exception("Detailed AlphaGenome track prediction failed.")
        raise TrackPredictionError(
            "AlphaGenome could not generate detailed tracks for this request. Atlas summary results remain available."
        ) from exc
