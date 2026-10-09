"""Background VCF jobs with local persistence for review and export."""

from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
import hmac
import secrets
from threading import Lock
from uuid import uuid4
from datetime import timedelta

from app.models.variants import VariantAnalysisResult, VariantFailure, VariantRequest
from app.services.alphagenome import AlphaGenomeNotConfigured, AlphaGenomeRequestError, AlphaGenomeService
from app.services.store import save_analysis

_jobs: dict[str, dict] = {}
_lock = Lock()
_service = AlphaGenomeService()
_MAX_RETAINED_JOBS = 100
_MAX_ACTIVE_JOBS = 3
_JOB_RETENTION = timedelta(hours=2)


class JobCapacityError(Exception):
    pass


def _sdk_version() -> str | None:
    try:
        return version("alphagenome")
    except PackageNotFoundError:
        return None


def submit(variants: list[VariantRequest]) -> dict:
    job_id = str(uuid4())
    created_at = datetime.now(timezone.utc).isoformat()
    job = {
        "analysis_id": job_id,
        "access_token": secrets.token_urlsafe(32),
        "created_at": created_at,
        "analysis_type": "vcf_batch",
        "status": "queued",
        "requested": len(variants),
        "completed": 0,
        "results": [],
        "failures": [],
        "manifest": {"genome_assembly": "GRCh38", "provider": "AlphaGenome Atlas", "sdk_version": _sdk_version(), "provider_model_version": None, "requested_scorers": ["AVI_SCORE", "AVI_SCORE_FEATURE_IMPORTANCE", "RNA_SEQ", "SPLICE_SITES", "ATAC"], "input_order_preserved": True, "input_variants": [v.model_dump(mode="json") for v in variants]},
    }
    with _lock:
        now = datetime.now(timezone.utc)
        for existing_id, existing in list(_jobs.items()):
            created = datetime.fromisoformat(existing["created_at"])
            if existing["status"] not in {"queued", "running"} and now - created > _JOB_RETENTION:
                _jobs.pop(existing_id, None)
        active_count = sum(item["status"] in {"queued", "running"} for item in _jobs.values())
        if active_count >= _MAX_ACTIVE_JOBS:
            raise JobCapacityError("The public analysis queue is busy. Please try again shortly.")
        while len(_jobs) >= _MAX_RETAINED_JOBS:
            removable = next((key for key, item in _jobs.items() if item["status"] not in {"queued", "running"}), None)
            if removable is None:
                raise JobCapacityError("The public analysis queue is busy. Please try again shortly.")
            _jobs.pop(removable, None)
        _jobs[job_id] = job
    save_analysis(job_id, created_at, "vcf_batch", "queued", job["manifest"], {})
    return {
        "job_id": job_id,
        "access_token": job["access_token"],
        "status": "queued",
        "requested": len(variants),
    }


def get(job_id: str, access_token: str | None) -> dict | None:
    with _lock:
        job = _jobs.get(job_id)
        if not job or not access_token or not hmac.compare_digest(job["access_token"], access_token):
            return None
        return {key: value for key, value in job.items() if key != "access_token"}


def run(job_id: str, variants: list[VariantRequest]) -> None:
    with _lock:
        job = _jobs.get(job_id)
        if not job:
            return
        job["status"] = "running"

    results_by_index: dict[int, VariantAnalysisResult] = {}
    failures_by_index: dict[int, VariantFailure] = {}
    try:
        with ThreadPoolExecutor(max_workers=5) as executor:
            futures = {
                executor.submit(_service.analyze, variant): (index, variant)
                for index, variant in enumerate(variants)
            }
            for future in as_completed(futures):
                index, variant = futures[future]
                try:
                    results_by_index[index] = future.result()
                except AlphaGenomeNotConfigured:
                    with _lock:
                        job["status"] = "failed"
                        job["error"] = "AlphaGenome API key is not configured."
                    break
                except AlphaGenomeRequestError as exc:
                    failures_by_index[index] = VariantFailure(variant=variant, error=str(exc))
                except Exception:
                    failures_by_index[index] = VariantFailure(variant=variant, error="Unexpected provider error.")
                with _lock:
                    job["completed"] = len(results_by_index) + len(failures_by_index)
                    job["results"] = [results_by_index[i].model_dump(mode="json") for i in sorted(results_by_index)]
                    job["failures"] = [failures_by_index[i].model_dump(mode="json") for i in sorted(failures_by_index)]
                    if job["completed"] % 25 == 0:
                        save_analysis(
                            job_id, job["created_at"], "vcf_batch", "running", job["manifest"],
                            {"results": job["results"], "failures": job["failures"], "completed": job["completed"]},
                        )
        with _lock:
            if job["status"] != "failed":
                job["status"] = "complete"
            job["completed"] = len(results_by_index) + len(failures_by_index)
            job["results"] = [results_by_index[i].model_dump(mode="json") for i in sorted(results_by_index)]
            job["failures"] = [failures_by_index[i].model_dump(mode="json") for i in sorted(failures_by_index)]
            save_analysis(
                job_id, job["created_at"], "vcf_batch", job["status"], job["manifest"],
                {"results": job["results"], "failures": job["failures"], "completed": job["completed"], "input_order_preserved": True},
            )
    except Exception:
        with _lock:
            job["status"] = "failed"
            job["error"] = "Batch analysis failed."
            save_analysis(job_id, job["created_at"], "vcf_batch", "failed", job["manifest"], {"error": job["error"]})
