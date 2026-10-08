import os

from fastapi import APIRouter, BackgroundTasks, Header, HTTPException, Query
from pydantic import BaseModel, Field

from app.models.variants import ModalityEffect, VariantAnalysisResult, VariantRequest
from app.services.alphagenome import (
    AlphaGenomeNotConfigured,
    AlphaGenomeRequestError,
    AlphaGenomeService,
)
from app.services import jobs
from app.services.jobs import JobCapacityError
from app.services.evidence import EvidenceProviderError, EvidenceService
from app.services.store import get_analysis, list_analyses, save_analysis
from app.services.statistics import compare_groups, gene_enrichment
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version

router = APIRouter(tags=["variants"])
alphagenome = AlphaGenomeService()
evidence = EvidenceService()


def _sdk_version() -> str | None:
    try:
        return version("alphagenome")
    except PackageNotFoundError:
        return None


class BatchRequest(BaseModel):
    variants: list[VariantRequest] = Field(min_length=1, max_length=5000)


class CompareRequest(BaseModel):
    case_scores: list[float] = Field(min_length=2, max_length=50000)
    control_scores: list[float] = Field(min_length=2, max_length=50000)
    iterations: int = Field(default=10000, ge=100, le=50000)
    seed: int = Field(default=42, ge=0, le=2147483647)


class EnrichmentRequest(BaseModel):
    genes: list[str] = Field(min_length=1, max_length=5000)
    background: list[str] = Field(min_length=2, max_length=50000)
    gene_sets: dict[str, list[str]]


class EvidenceSynthesisRequest(VariantRequest):
    rsid: str | None = Field(default=None, pattern=r"^rs[0-9]+$")


class VariantGeneRequest(VariantRequest):
    effects: list[ModalityEffect] = Field(default_factory=list, max_length=5000)
    window_bp: int = Field(default=100_000, ge=1, le=1_000_000)


@router.post("/variants/analyze", response_model=VariantAnalysisResult)
def analyze_variant(request: VariantRequest) -> VariantAnalysisResult:
    try:
        result = alphagenome.analyze(request)
    except AlphaGenomeNotConfigured as exc:
        raise HTTPException(status_code=503, detail=str(exc)) from exc
    except AlphaGenomeRequestError as exc:
        raise HTTPException(status_code=exc.status_code, detail=str(exc)) from exc
    save_analysis(
        result.analysis_id,
        result.created_at.isoformat(),
        "single_variant",
        "complete",
        {"variant": request.model_dump(mode="json"), "genome_assembly": "GRCh38", "provider": "AlphaGenome Atlas", "sdk_version": _sdk_version(), "provider_model_version": "not exposed by the Atlas response"},
        result.model_dump(mode="json"),
    )
    return result


@router.post("/variants/analyze-batch", status_code=202)
def analyze_batch(request: BatchRequest, background_tasks: BackgroundTasks) -> dict:
    if os.getenv("PUBLIC_DEPLOYMENT", "false").lower() in {"1", "true", "yes"} and len(request.variants) > 100:
        raise HTTPException(status_code=422, detail="Public batch submissions are limited to 100 SNVs.")
    try:
        job = jobs.submit(request.variants)
    except JobCapacityError as exc:
        raise HTTPException(status_code=429, detail=str(exc)) from exc
    background_tasks.add_task(jobs.run, job["job_id"], request.variants)
    return job


@router.get("/jobs/{job_id}")
def get_job(job_id: str, authorization: str | None = Header(default=None)) -> dict:
    access_token = authorization.removeprefix("Bearer ") if authorization else None
    job = jobs.get(job_id, access_token)
    if job is None:
        raise HTTPException(status_code=404, detail="Analysis job not found or access token is invalid.")
    return job


@router.get("/analyses")
def analyses(limit: int = Query(default=50, ge=1, le=200)) -> list[dict]:
    return list_analyses(limit)


@router.get("/analyses/{analysis_id}")
def analysis(analysis_id: str) -> dict:
    record = get_analysis(analysis_id)
    if record is None:
        raise HTTPException(status_code=404, detail="Analysis not found.")
    return record


def _evidence_call(callable_, *args):
    try:
        return callable_(*args)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except EvidenceProviderError as exc:
        raise HTTPException(status_code=502, detail=str(exc)) from exc


@router.get("/evidence/gene/{symbol}")
def gene_evidence(symbol: str) -> dict:
    return _evidence_call(evidence.search_gene, symbol)


@router.get("/evidence/region")
def region_evidence(chromosome: str, start: int, end: int) -> dict:
    return _evidence_call(evidence.region, chromosome, start, end)


@router.get("/evidence/variant/{rsid}")
def variant_evidence(rsid: str) -> dict:
    return _evidence_call(evidence.variant_evidence, rsid)


@router.get("/evidence/gwas/{rsid}")
def gwas_evidence(rsid: str) -> dict:
    return _evidence_call(evidence.gwas_associations, rsid)


@router.get("/evidence/clinvar/{rsid}")
def clinvar_evidence(rsid: str) -> dict:
    return _evidence_call(evidence.clinvar, rsid)


@router.get("/evidence/gtex/{rsid}")
def gtex_evidence(rsid: str) -> dict:
    return _evidence_call(evidence.gtex_eqtls, rsid)


@router.get("/evidence/gnomad")
def gnomad_evidence(chromosome: str, position: int, reference: str, alternate: str) -> dict:
    return _evidence_call(evidence.gnomad_variant, chromosome, position, reference, alternate)


@router.get("/evidence/encode/region")
def encode_region_evidence(chromosome: str, start: int, end: int) -> dict:
    return _evidence_call(evidence.encode_region, chromosome, start, end)


@router.post("/evidence/synthesize")
def synthesize_evidence(request: EvidenceSynthesisRequest) -> dict:
    return _evidence_call(
        evidence.synthesize_variant,
        request.chromosome,
        request.position,
        request.reference,
        request.alternate,
        request.rsid,
    )


@router.post("/analysis/variant-to-gene")
def variant_to_gene(request: VariantGeneRequest) -> dict:
    return _evidence_call(
        evidence.prioritize_variant_genes,
        request.chromosome,
        request.position,
        [effect.model_dump() for effect in request.effects],
        request.window_bp,
    )


@router.post("/statistics/compare")
def compare(request: CompareRequest) -> dict:
    try:
        result = compare_groups(request.case_scores, request.control_scores, request.iterations, request.seed)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    analysis_id = result["analysis_id"]
    save_analysis(analysis_id, result["created_at"], "statistical_comparison", "complete", request.model_dump(), result)
    return result


@router.post("/statistics/enrichment")
def enrichment(request: EnrichmentRequest) -> dict:
    try:
        result = gene_enrichment(request.genes, request.background, request.gene_sets)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    created = datetime.now(timezone.utc).isoformat()
    analysis_id = result["analysis_id"]
    result["created_at"] = created
    save_analysis(analysis_id, created, "gene_enrichment", "complete", request.model_dump(), result)
    return result
