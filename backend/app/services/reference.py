"""GRCh38 single-base validation using Ensembl's reference sequence service."""

from datetime import datetime, timezone
from functools import lru_cache

import httpx

from app.models.variants import VariantRequest

ENSEMBL_SEQUENCE_URL = "https://rest.ensembl.org/sequence/region/human"


class ReferenceValidationError(Exception):
    """Reference validation failed or could not be completed safely."""

    def __init__(self, message: str, status_code: int = 502):
        super().__init__(message)
        self.status_code = status_code


class ReferenceGenomeService:
    """Verify the requested REF allele against the Ensembl GRCh38 sequence."""

    def __init__(self, client: httpx.Client | None = None):
        self._client = client

    def verify(self, variant: VariantRequest) -> dict:
        chromosome = variant.chromosome.removeprefix("chr")
        # Ensembl calls the mitochondrial reference sequence MT (the UI also
        # accepts the common chrM alias).
        chromosome = "MT" if chromosome == "M" else chromosome
        coordinate = (
            f"{chromosome}:{variant.position}"
            f"..{variant.position}:1"
        )
        url = f"{ENSEMBL_SEQUENCE_URL}/{coordinate}?coord_system_version=GRCh38"
        try:
            if self._client is None:
                observed, retrieved_at = self._fetch_reference_base(coordinate)
            else:
                observed = self._fetch_with_client(coordinate)
                retrieved_at = datetime.now(timezone.utc).isoformat()
        except ReferenceValidationError:
            raise
        except (httpx.TimeoutException, TimeoutError) as exc:
            raise ReferenceValidationError(
                "GRCh38 reference verification timed out; AlphaGenome was not queried.", 504
            ) from exc
        except httpx.HTTPError as exc:
            raise ReferenceValidationError(
                "GRCh38 reference verification is unavailable; AlphaGenome was not queried.", 502
            ) from exc

        if observed != variant.reference:
            raise ReferenceValidationError(
                f"Reference allele mismatch at {variant.chromosome}:{variant.position} "
                f"(GRCh38): submitted REF={variant.reference}, Ensembl reports {observed}. "
                "The input was not corrected and AlphaGenome was not queried.",
                422,
            )
        return {
            "status": "verified",
            "assembly": "GRCh38",
            "source": "Ensembl REST sequence/region",
            "source_url": url,
            "requested_reference": variant.reference,
            "observed_reference": observed,
            "checked_at": retrieved_at,
        }

    def _fetch_with_client(self, coordinate: str) -> str:
        response = self._client.get(
            f"{ENSEMBL_SEQUENCE_URL}/{coordinate}",
            params={"coord_system_version": "GRCh38"},
            headers={"Content-Type": "text/plain", "Accept": "text/plain"},
            timeout=8.0,
        )
        if response.status_code == 400 or response.status_code == 404:
            raise ReferenceValidationError(
                "Ensembl could not resolve this GRCh38 coordinate; AlphaGenome was not queried.",
                422,
            )
        response.raise_for_status()
        sequence = response.text.strip().upper()
        if len(sequence) != 1 or sequence not in {"A", "C", "G", "T", "N"}:
            raise ReferenceValidationError(
                "Ensembl returned an invalid single-base GRCh38 sequence response; "
                "AlphaGenome was not queried.",
                502,
            )
        return sequence

    @staticmethod
    @lru_cache(maxsize=8192)
    def _fetch_reference_base(coordinate: str) -> tuple[str, str]:
        try:
            response = httpx.get(
                f"{ENSEMBL_SEQUENCE_URL}/{coordinate}",
                params={"coord_system_version": "GRCh38"},
                headers={"Content-Type": "text/plain", "Accept": "text/plain"},
                timeout=httpx.Timeout(8.0, connect=4.0),
            )
            if response.status_code in {400, 404}:
                raise ReferenceValidationError(
                    "Ensembl could not resolve this GRCh38 coordinate; AlphaGenome was not queried.",
                    422,
                )
            response.raise_for_status()
        except ReferenceValidationError:
            raise
        except httpx.TimeoutException:
            raise
        except httpx.HTTPError:
            raise
        sequence = response.text.strip().upper()
        if len(sequence) != 1 or sequence not in {"A", "C", "G", "T", "N"}:
            raise ReferenceValidationError(
                "Ensembl returned an invalid single-base GRCh38 sequence response; "
                "AlphaGenome was not queried.",
                502,
            )
        return sequence, datetime.now(timezone.utc).isoformat()
