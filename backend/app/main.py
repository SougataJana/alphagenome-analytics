import os
import time
from collections import defaultdict, deque
from threading import Lock

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from app.routers.variants import router as variants_router

app = FastAPI(
    title="AlphaGenome Analytics API",
    version="0.1.0",
    description="A research interface for genomic prediction and evidence review.",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        origin.strip()
        for origin in os.getenv(
            "CORS_ORIGINS", "http://localhost:3000,http://127.0.0.1:3000"
        ).split(",")
        if origin.strip()
    ],
    allow_credentials=False,
    allow_methods=["GET", "POST"],
    allow_headers=["Content-Type", "Authorization"],
)

_rate_lock = Lock()
_rate_windows: dict[tuple[str, str], deque[float]] = defaultdict(deque)


def _rate_rule(method: str, path: str) -> tuple[str, int, int] | None:
    if method == "OPTIONS":
        return None
    if method == "POST" and path == "/api/v1/variants/analyze":
        return "single-variant", 12, 60
    if method == "POST" and path == "/api/v1/variants/analyze-batch":
        return "batch", 2, 60
    if method == "POST" and path.startswith("/api/v1/statistics/"):
        return "statistics", 10, 60
    if path.startswith("/api/v1/"):
        return "api", 60, 60
    return None


@app.middleware("http")
async def limit_public_requests(request, call_next):
    """Best-effort per-process rate limits for the single-instance public service."""
    rule = _rate_rule(request.method, request.url.path)
    if rule is None:
        return await call_next(request)

    bucket, limit, window_seconds = rule
    client_ip = request.client.host if request.client else "unknown"
    now = time.monotonic()
    key = (client_ip, bucket)
    with _rate_lock:
        hits = _rate_windows[key]
        while hits and now - hits[0] >= window_seconds:
            hits.popleft()
        if len(hits) >= limit:
            retry_after = max(1, int(window_seconds - (now - hits[0])))
            return JSONResponse(
                status_code=429,
                content={"detail": "Request rate limit reached. Wait briefly and retry."},
                headers={"Retry-After": str(retry_after), "X-RateLimit-Limit": str(limit)},
            )
        hits.append(now)
        if len(_rate_windows) > 10_000:
            expired = [
                old_key for old_key, old_hits in _rate_windows.items()
                if not old_hits or now - old_hits[-1] >= window_seconds
            ]
            for old_key in expired:
                _rate_windows.pop(old_key, None)

    response = await call_next(request)
    response.headers["X-RateLimit-Limit"] = str(limit)
    with _rate_lock:
        response.headers["X-RateLimit-Remaining"] = str(max(0, limit - len(_rate_windows[key])))
    return response


@app.get("/health", tags=["system"])
def health() -> dict[str, str]:
    return {"status": "ok", "service": "alphagenome-analytics-api"}


app.include_router(variants_router, prefix="/api/v1")
