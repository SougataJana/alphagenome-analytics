from typing import Any, Literal
from datetime import datetime, timezone
from uuid import uuid4

from pydantic import BaseModel, Field, field_validator, model_validator


GRCH38_CHROMOSOME_LENGTHS = {
    "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555,
    "5": 181538259, "6": 170805979, "7": 159345973, "8": 145138636,
    "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
    "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345,
    "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167,
    "21": 46709983, "22": 50818468, "X": 156040895, "Y": 57227415,
    "M": 16569, "MT": 16569,
}


def normalize_grch38_chromosome(chromosome: str) -> str:
    chrom = chromosome.strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    chrom = chrom.upper()
    if chrom.isdigit():
        chrom = str(int(chrom))
    if chrom not in GRCH38_CHROMOSOME_LENGTHS:
        raise ValueError("Chromosome must be chr1–chr22, chrX, chrY, or chrM.")
    return chrom


class VariantRequest(BaseModel):
    """A single reference-genome variant, using a 1-based position."""

    chromosome: str = Field(examples=["chr7"])
    position: int = Field(gt=0, examples=[140453136])
    reference: str = Field(min_length=1, max_length=1, pattern=r"^[ACGT]$", examples=["A"])
    alternate: str = Field(min_length=1, max_length=1, pattern=r"^[ACGT]$", examples=["G"])
    genome_assembly: Literal["GRCh38"] = Field(default="GRCh38", examples=["GRCh38"])

    @field_validator("chromosome", mode="before")
    @classmethod
    def normalize_chromosome(cls, value: str) -> str:
        if not isinstance(value, str):
            raise ValueError("Chromosome must be a GRCh38 chromosome name.")
        chromosome = normalize_grch38_chromosome(value)
        return f"chr{chromosome}"

    @field_validator("reference", "alternate", mode="before")
    @classmethod
    def normalize_allele(cls, value: str) -> str:
        if isinstance(value, str):
            return value.upper()
        return value

    @model_validator(mode="after")
    def alleles_must_differ(self) -> "VariantRequest":
        if self.reference == self.alternate:
            raise ValueError("Reference and alternate alleles must differ.")
        chromosome = self.chromosome.removeprefix("chr")
        if self.position > GRCH38_CHROMOSOME_LENGTHS[chromosome]:
            raise ValueError(f"Position exceeds the GRCh38 length of chromosome {chromosome}.")
        return self


class FeatureAttribution(BaseModel):
    feature: str
    value: float = Field(allow_inf_nan=False)


class ModalityEffect(BaseModel):
    modality: str
    feature: str | None = None
    gene: str | None = None
    tissue: str | None = None
    raw_score: float = Field(allow_inf_nan=False)
    # Signed scorer quantiles preserve direction and range from -1 to 1;
    # unsigned scorer quantiles use the non-negative portion of that range.
    quantile_score: float | None = Field(default=None, ge=-1, le=1, allow_inf_nan=False)


class ReferenceValidation(BaseModel):
    status: Literal["verified"]
    assembly: Literal["GRCh38"]
    source: str
    source_url: str
    requested_reference: str = Field(pattern=r"^[ACGT]$")
    observed_reference: str = Field(pattern=r"^[ACGT]$")
    checked_at: datetime


class VariantAnalysisResult(BaseModel):
    analysis_id: str = Field(default_factory=lambda: str(uuid4()))
    created_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))
    status: Literal["complete"]
    provider: Literal["AlphaGenome Atlas"]
    provider_sdk_version: str | None = None
    provider_model_version: str | None = None
    variant: VariantRequest
    reference_validation: ReferenceValidation
    avi_score: float = Field(allow_inf_nan=False)
    avi_quantile: float | None = Field(default=None, ge=0, le=1, allow_inf_nan=False)
    requested_scorers: list[str] = Field(default_factory=list)
    feature_importance: list[FeatureAttribution]
    effects: list[ModalityEffect] = Field(default_factory=list)
    effects_total: int = Field(default=0, ge=0)
    effects_total_by_scorer: dict[str, int] = Field(default_factory=dict)
    effects_truncated: bool = False
    provider_response_available: bool = False
    provider_response_endpoint: str | None = None
    atlas_sdk_response: dict[str, Any] | None = Field(default=None, exclude=True)


class VariantFailure(BaseModel):
    variant: VariantRequest
    error: str
