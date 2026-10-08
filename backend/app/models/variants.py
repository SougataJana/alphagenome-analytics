from typing import Literal
from datetime import datetime, timezone
from uuid import uuid4

from pydantic import BaseModel, Field, model_validator


class VariantRequest(BaseModel):
    """A single reference-genome variant, using a 1-based position."""

    chromosome: str = Field(pattern=r"^(chr)?([0-9]+|X|Y|M|MT)$", examples=["chr7"])
    position: int = Field(gt=0, examples=[140453136])
    reference: str = Field(min_length=1, max_length=1, pattern=r"^[ACGT]$", examples=["A"])
    alternate: str = Field(min_length=1, max_length=1, pattern=r"^[ACGT]$", examples=["G"])
    genome_assembly: Literal["GRCh38"] = Field(default="GRCh38", examples=["GRCh38"])

    @model_validator(mode="after")
    def alleles_must_differ(self) -> "VariantRequest":
        if self.reference == self.alternate:
            raise ValueError("Reference and alternate alleles must differ.")
        return self


class FeatureAttribution(BaseModel):
    feature: str
    value: float


class ModalityEffect(BaseModel):
    modality: str
    feature: str | None = None
    gene: str | None = None
    tissue: str | None = None
    raw_score: float
    quantile_score: float | None = None


class VariantAnalysisResult(BaseModel):
    analysis_id: str = Field(default_factory=lambda: str(uuid4()))
    created_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))
    status: Literal["complete"]
    provider: Literal["AlphaGenome Atlas"]
    provider_sdk_version: str | None = None
    provider_model_version: str | None = None
    variant: VariantRequest
    avi_score: float
    avi_quantile: float | None = None
    feature_importance: list[FeatureAttribution]
    effects: list[ModalityEffect] = Field(default_factory=list)


class VariantFailure(BaseModel):
    variant: VariantRequest
    error: str
