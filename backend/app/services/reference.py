"""GRCh38 single-base validation using independent public reference services."""

from datetime import datetime, timezone
from functools import lru_cache
import time
from urllib.parse import urlencode

import httpx

from app.models.variants import VariantRequest

ENSEMBL_SEQUENCE_URL = "https://rest.ensembl.org/sequence/region/human"
UCSC_SEQUENCE_URL = "https://api.genome.ucsc.edu/getData/sequence"
_DNA_BASES = {"A", "C", "G", "T"}


class ReferenceValidationError(Exception):
    """Reference validation failed or could not be completed safely."""

    def __init__(self, message: str, status_code: int = 502):
        super().__init__(message)
        self.status_code = status_code


class ReferenceGenomeService:
    """Verify the requested REF allele against GRCh38 before Atlas analysis.

    Ensembl is queried first. Transient Ensembl errors are retried once and then
    fall back to UCSC's hg38 sequence endpoint. A mismatch from a successful
    provider remains a hard error; the service never changes the user's allele.
    """

    def __init__(self, client: httpx.Client | None = None):
        self._client = client

    def verify(self, variant: VariantRequest) -> dict:
        chromosome = variant.chromosome.removeprefix("chr")
        # Ensembl and UCSC use MT and chrM, respectively.
        ensembl_chromosome = "MT" if chromosome in {"M", "MT"} else chromosome
        ucsc_chromosome = "chrM" if chromosome in {"M", "MT"} else f"chr{chromosome}"
        coordinate = f"{ensembl_chromosome}:{variant.position}..{variant.position}:1"
        ensembl_url = (
            f"{ENSEMBL_SEQUENCE_URL}/{coordinate}?coord_system_version=GRCh38"
        )

        if self._client is None:
            observed, retrieved_at, source, source_url = self._fetch_reference_base(
                ensembl_chromosome, ucsc_chromosome, variant.position
            )
        else:
            observed, source, source_url = self._fetch_with_client(
                coordinate,
                ensembl_chromosome,
                ucsc_chromosome,
                variant.position,
            )
            retrieved_at = datetime.now(timezone.utc).isoformat()

        if observed != variant.reference:
            provider_label = "Ensembl" if source.startswith("Ensembl") else "UCSC"
            raise ReferenceValidationError(
                f"Reference allele mismatch at {variant.chromosome}:{variant.position} "
                f"(GRCh38): submitted REF={variant.reference}, {provider_label} reports {observed}. "
                "The input was not corrected and AlphaGenome was not queried.",
                422,
            )
        return {
            "status": "verified",
            "assembly": "GRCh38",
            "source": source,
            "source_url": source_url or ensembl_url,
            "requested_reference": variant.reference,
            "observed_reference": observed,
            "checked_at": retrieved_at,
        }

    def _fetch_with_client(
        self,
        coordinate: str,
        ensembl_chromosome: str,
        ucsc_chromosome: str,
        position: int,
    ) -> tuple[str, str, str]:
        getter = self._client.get
        ensembl_error: Exception | None = None
        for attempt in range(2):
            try:
                return self._fetch_ensembl(getter, coordinate)
            except ReferenceValidationError:
                raise
            except (httpx.TimeoutException, httpx.NetworkError, httpx.HTTPStatusError) as exc:
                ensembl_error = exc
            if attempt == 0:
                time.sleep(0.25)

        try:
            return self._fetch_ucsc(
                getter, ucsc_chromosome, position, ensembl_error
            )
        except ReferenceValidationError:
            raise

    @staticmethod
    @lru_cache(maxsize=8192)
    def _fetch_reference_base(
        ensembl_chromosome: str,
        ucsc_chromosome: str,
        position: int,
    ) -> tuple[str, str, str, str]:
        coordinate = f"{ensembl_chromosome}:{position}..{position}:1"
        ensembl_error: Exception | None = None
        for attempt in range(2):
            try:
                result = ReferenceGenomeService._fetch_ensembl(httpx.get, coordinate)
                return (*result, datetime.now(timezone.utc).isoformat())
            except ReferenceValidationError:
                raise
            except (httpx.TimeoutException, httpx.NetworkError, httpx.HTTPStatusError) as exc:
                ensembl_error = exc
            if attempt == 0:
                time.sleep(0.25)

        try:
            result = ReferenceGenomeService._fetch_ucsc(
                httpx.get, ucsc_chromosome, position, ensembl_error
            )
            return (*result, datetime.now(timezone.utc).isoformat())
        except ReferenceValidationError:
            raise

    @staticmethod
    def _fetch_ensembl(getter, coordinate: str) -> tuple[str, str, str]:
        url = f"{ENSEMBL_SEQUENCE_URL}/{coordinate}"
        response = getter(
            url,
            params={"coord_system_version": "GRCh38"},
            headers={"Content-Type": "text/plain", "Accept": "text/plain"},
            timeout=httpx.Timeout(5.0, connect=3.0),
        )
        if response.status_code in {400, 404}:
            raise ReferenceValidationError(
                "Ensembl could not resolve this GRCh38 coordinate; AlphaGenome was not queried.",
                422,
            )
        response.raise_for_status()
        sequence = response.text.strip().upper()
        if len(sequence) != 1 or sequence not in _DNA_BASES:
            raise ReferenceValidationError(
                "Ensembl returned an invalid single-base GRCh38 sequence response; "
                "AlphaGenome was not queried.",
                502,
            )
        return sequence, "Ensembl REST sequence/region", url + "?coord_system_version=GRCh38"

    @staticmethod
    def _fetch_ucsc(
        getter,
        chromosome: str,
        position: int,
        ensembl_error: Exception | None,
    ) -> tuple[str, str, str]:
        # UCSC start is zero-relative and end is one-relative, so [p-1, p)
        # retrieves exactly the GRCh38 base at the app's 1-based position p.
        params = {
            "genome": "hg38",
            "chrom": chromosome,
            "start": position - 1,
            "end": position,
        }
        url = f"{UCSC_SEQUENCE_URL}?{urlencode(params)}"
        try:
            response = getter(
                UCSC_SEQUENCE_URL,
                params=params,
                headers={"Accept": "application/json"},
                timeout=httpx.Timeout(6.0, connect=3.0),
            )
            response.raise_for_status()
            sequence = str(response.json().get("dna", "")).strip().upper()
        except (httpx.HTTPError, ValueError, AttributeError) as exc:
            raise ReferenceValidationError(
                "GRCh38 reference verification is unavailable from Ensembl and UCSC; "
                "AlphaGenome was not queried. Please retry later.",
                502,
            ) from (exc if exc is not None else ensembl_error)

        if len(sequence) != 1 or sequence not in _DNA_BASES:
            raise ReferenceValidationError(
                "UCSC did not return one unambiguous GRCh38 base after Ensembl "
                "verification failed; AlphaGenome was not queried.",
                502,
            )
        return sequence, "UCSC Genome Browser API (hg38 sequence)", url
