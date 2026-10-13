"""Optional, on-demand AlphaGenome REF/ALT signal tracks for visualization.

These predictions are deliberately separate from precomputed Atlas lookups.
"""

import logging
import math
import os
import base64
from io import BytesIO
from functools import lru_cache
from importlib.metadata import PackageNotFoundError, version
from datetime import datetime, timezone
import httpx

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
    "SPLICE_JUNCTIONS",
    "CONTACT_MAPS",
)
MAX_TRACKS = 8
MAX_PLOT_POINTS = 1200
PLOT_WINDOW_BP = 50_000
CONTACT_WINDOW_BP = 250_000
VISUALIZATION_PROFILES = {
    "expression_chromatin": ("RNA_SEQ", "CAGE", "ATAC", "DNASE"),
    "splicing": ("RNA_SEQ", "SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"),
    "contact_maps": ("CONTACT_MAPS",),
}
GENCODE_FEATHER_URL = (
    "https://storage.googleapis.com/alphagenome/reference/gencode/"
    "hg38/gencode.v46.annotation.gtf.gz.feather"
)
ENABLE_GENCODE_TRANSCRIPTS_ENV = "AGA_ENABLE_GENCODE_TRANSCRIPTS"


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
            "generated_at": datetime.now(timezone.utc).isoformat(),
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


def _selected_tracks(track_data, ontology_term: str):
    metadata = track_data.metadata
    if "ontology_curie" not in metadata.columns:
        return track_data
    mask = metadata["ontology_curie"].astype(str) == ontology_term
    return track_data.filter_tracks(mask.to_numpy())


@lru_cache(maxsize=1)
def _transcript_extractor():
    """Load the official GENCODE v46 GRCh38 annotation once per process."""
    import pandas as pd
    import pyarrow as pa
    from alphagenome.data import gene_annotation, transcript

    response = httpx.get(GENCODE_FEATHER_URL, timeout=httpx.Timeout(30.0, connect=8.0))
    response.raise_for_status()
    gtf = pd.read_feather(pa.BufferReader(response.content))
    supported = gene_annotation.filter_transcript_support_level(
        gene_annotation.filter_protein_coding(gtf), ["1"]
    )
    return transcript.TranscriptExtractor(supported)


def _render_advanced_plot(output, variant, profile: str, ontology_term: str) -> tuple[str, str | None]:
    """Render official AlphaGenome plotting components from returned model data."""
    import matplotlib

    matplotlib.use("Agg")
    from alphagenome.data import genome
    from alphagenome.visualization import plot_components
    import matplotlib.pyplot as plt

    components = []
    if profile == "expression_chromatin":
        for modality in VISUALIZATION_PROFILES[profile]:
            ref = getattr(output.reference, modality.lower(), None)
            alt = getattr(output.alternate, modality.lower(), None)
            if ref is None or alt is None or not ref.values.size:
                continue
            ref = _selected_tracks(ref, ontology_term)
            alt = _selected_tracks(alt, ontology_term)
            if not ref.values.size or not alt.values.size:
                continue
            components.append(
                plot_components.OverlaidTracks(
                    tdata={"REF": ref, "ALT": alt},
                    colors={"REF": "#63766b", "ALT": "#c76d55"},
                    ylabel_template=f"{modality}: {{biosample_name}} ({{strand}})\\n{{name}}",
                    max_num_tracks=6,
                )
            )
        window = PLOT_WINDOW_BP
        title = "AlphaGenome REF vs ALT · expression and chromatin predictions"
    elif profile == "splicing":
        for modality in ("RNA_SEQ", "SPLICE_SITES", "SPLICE_SITE_USAGE"):
            ref = getattr(output.reference, modality.lower(), None)
            alt = getattr(output.alternate, modality.lower(), None)
            if ref is None or alt is None or not ref.values.size:
                continue
            ref = _selected_tracks(ref, ontology_term)
            alt = _selected_tracks(alt, ontology_term)
            if not ref.values.size or not alt.values.size:
                continue
            components.append(
                plot_components.OverlaidTracks(
                    tdata={"REF": ref, "ALT": alt},
                    colors={"REF": "#63766b", "ALT": "#c76d55"},
                    ylabel_template=f"{modality}: {{biosample_name}} ({{strand}})\\n{{name}}",
                    max_num_tracks=6,
                )
            )
        ref_junctions = getattr(output.reference, "splice_junctions", None)
        alt_junctions = getattr(output.alternate, "splice_junctions", None)
        if ref_junctions is not None and ref_junctions.values.size:
            ref_junctions = _selected_tracks(ref_junctions, ontology_term)
        if ref_junctions is not None and ref_junctions.values.size:
            components.append(
                plot_components.Sashimi(
                    _selected_tracks(ref_junctions, ontology_term),
                    ylabel_template="REF splice junctions · {biosample_name} ({strand})\\n{name}",
                )
            )
        if alt_junctions is not None and alt_junctions.values.size:
            alt_junctions = _selected_tracks(alt_junctions, ontology_term)
        if alt_junctions is not None and alt_junctions.values.size:
            components.append(
                plot_components.Sashimi(
                    _selected_tracks(alt_junctions, ontology_term),
                    ylabel_template="ALT splice junctions · {biosample_name} ({strand})\\n{name}",
                )
            )
        window = PLOT_WINDOW_BP
        title = "AlphaGenome REF vs ALT · splicing predictions"
    elif profile == "contact_maps":
        for label, side in (("REF", output.reference), ("ALT", output.alternate)):
            contact = getattr(side, "contact_maps", None)
            if contact is None or not contact.values.size:
                continue
            contact = _selected_tracks(contact, ontology_term)
            if not contact.values.size:
                continue
            components.append(
                plot_components.ContactMaps(
                    tdata=contact,
                    ylabel_template=f"{label} · {{biosample_name}}\\n{{name}}",
                    cmap="viridis",
                )
            )
        window = CONTACT_WINDOW_BP
        title = "AlphaGenome REF and ALT · predicted contact maps"
    else:
        raise TrackPredictionError("Choose a supported visualization profile.")

    if not components:
        raise TrackPredictionError(
            "AlphaGenome returned no plottable tracks for this profile and biosample. Try another profile or biosample."
        )
    genome_variant = genome.Variant(
        chromosome=variant.chromosome,
        position=variant.position,
        reference_bases=variant.reference,
        alternate_bases=variant.alternate,
    )
    display_interval = genome_variant.reference_interval.resize(window)
    # Loading the full GENCODE Feather reference can use several hundred MB.
    # Keep it opt-in because the public Render free instance has a tight memory limit.
    include_transcripts = os.getenv(ENABLE_GENCODE_TRANSCRIPTS_ENV, "false").lower() in {
        "1", "true", "yes"
    }
    if include_transcripts:
        transcript_status = "included"
        try:
            transcripts = _transcript_extractor().extract(display_interval)
            if transcripts:
                components.insert(0, plot_components.TranscriptAnnotation(transcripts))
            else:
                transcript_status = "no protein-coding, support-level-1 transcripts overlap this display window"
        except Exception:
            logger.exception("Could not load the GENCODE transcript annotations for plotting.")
            transcript_status = "GENCODE annotation unavailable; prediction tracks are still shown"
    else:
        transcript_status = "omitted to keep the public backend within its memory limit"
    figure = plot_components.plot(
        components,
        interval=display_interval,
        annotations=[plot_components.VariantAnnotation([genome_variant])],
        title=title,
    )
    buffer = BytesIO()
    try:
        figure.savefig(buffer, format="svg", bbox_inches="tight")
        return base64.b64encode(buffer.getvalue()).decode("ascii"), transcript_status
    finally:
        plt.close(figure)


def predict_advanced_visualization(variant: VariantRequest, ontology_term: str, profile: str) -> dict:
    """Return plotted SVG from real REF/ALT predictions for one supported profile."""
    if profile not in VISUALIZATION_PROFILES:
        raise TrackPredictionError("Choose a supported AlphaGenome visualization profile.")
    if not ontology_term or len(ontology_term) > 80:
        raise TrackPredictionError("Choose a tissue or biosample ontology term.")
    try:
        from alphagenome.data import genome
        from alphagenome.models import dna_client

        model = _model_client()
        outputs = [getattr(dna_client.OutputType, name) for name in VISUALIZATION_PROFILES[profile]]
        variant_data = genome.Variant(
            chromosome=variant.chromosome,
            position=variant.position,
            reference_bases=variant.reference,
            alternate_bases=variant.alternate,
        )
        interval = variant_data.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
        prediction = model.predict_variant(
            interval=interval,
            variant=variant_data,
            requested_outputs=outputs,
            ontology_terms=[ontology_term],
        )
        plot_svg, transcript_status = _render_advanced_plot(
            prediction, variant, profile, ontology_term
        )
        return {
            "provider": "AlphaGenome API · on-demand predict_variant",
            "provider_sdk_version": _get_sdk_version(),
            "model_version": "Not exposed by the API response",
            "variant": variant.model_dump(mode="json"),
            "ontology_term": ontology_term,
            "profile": profile,
            "modalities_requested": list(VISUALIZATION_PROFILES[profile]),
            "assembly": "GRCh38",
            "output_kind": "Official AlphaGenome visualization components rendered from on-demand REF/ALT model predictions; separate from precomputed Atlas scores",
            "window_bp": CONTACT_WINDOW_BP if profile == "contact_maps" else PLOT_WINDOW_BP,
            "plot_svg_base64": plot_svg,
            "transcript_annotation": transcript_status,
            "transcript_source": GENCODE_FEATHER_URL if transcript_status != "omitted to keep the public backend within its memory limit" else None,
        }
    except TrackPredictionError:
        raise
    except TimeoutError as exc:
        raise TrackPredictionError(
            "Detailed AlphaGenome prediction timed out. Atlas summary results remain available; try again later."
        ) from exc
    except Exception as exc:
        logger.exception("Advanced AlphaGenome visualization failed.")
        raise TrackPredictionError(
            "AlphaGenome could not generate this detailed visualization. Atlas summary results remain available."
        ) from exc
