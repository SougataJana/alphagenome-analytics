import os

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

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
    allow_headers=["Content-Type"],
)


@app.get("/health", tags=["system"])
def health() -> dict[str, str]:
    return {"status": "ok", "service": "alphagenome-analytics-api"}


app.include_router(variants_router, prefix="/api/v1")
