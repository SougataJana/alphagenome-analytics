"""Read-only public evidence adapters. Every result carries its source."""

from urllib.parse import quote

import httpx
import os


class EvidenceProviderError(Exception):
    pass


class EvidenceProviderNotConfigured(EvidenceProviderError):
    pass


class EvidenceService:
    def __init__(self) -> None:
        self.timeout = httpx.Timeout(20.0, connect=8.0)
        self.headers = {
            "Accept": "application/json",
            "User-Agent": "AlphaGenomeAnalytics/0.1 (research prototype)",
        }

    def _get(
        self,
        url: str,
        params: dict | list[tuple[str, str]] | None = None,
        *,
        provider: str,
    ) -> httpx.Response:
        try:
            response = httpx.get(url, params=params, headers=self.headers, timeout=self.timeout)
            response.raise_for_status()
            return response
        except httpx.TimeoutException as exc:
            raise EvidenceProviderError(
                f"{provider} did not respond before the timeout. Please retry in a moment."
            ) from exc
        except httpx.HTTPStatusError as exc:
            raise EvidenceProviderError(
                f"{provider} returned HTTP {exc.response.status_code}. Please retry later."
            ) from exc
        except httpx.RequestError as exc:
            raise EvidenceProviderError(
                f"Could not connect to {provider}. Please retry in a moment."
            ) from exc

    @staticmethod
    def _json(response: httpx.Response, provider: str) -> dict | list:
        try:
            return response.json()
        except ValueError as exc:
            raise EvidenceProviderError(
                f"{provider} returned a response that could not be read. Please retry later."
            ) from exc

    def search_gene(self, symbol: str) -> dict:
        url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{quote(symbol, safe='')}"
        response = self._get(url, {"expand": "1"}, provider="Ensembl")
        gene = self._json(response, "Ensembl")
        if not isinstance(gene, dict) or not gene.get("id"):
            raise EvidenceProviderError("Ensembl did not resolve this gene symbol.")
        start = int(gene["start"])
        end = int(gene["end"])
        region = self.region(gene["seq_region_name"], max(1, start - 100_000), end + 100_000)
        return {
            "query": symbol,
            "gene": {
                "id": gene["id"],
                "symbol": gene.get("display_name", symbol),
                "biotype": gene.get("biotype"),
                "chromosome": gene.get("seq_region_name"),
                "start": start,
                "end": end,
                "strand": gene.get("strand"),
                "assembly": gene.get("assembly_name", "GRCh38"),
                "source": "Ensembl",
            },
            "nearby_features": region["features"],
            "source_url": response.url.__str__(),
        }

    def region(self, chromosome: str, start: int, end: int) -> dict:
        if start < 1 or end < start or end - start > 5_000_000:
            raise ValueError("Region must be positive, ordered, and no wider than 5 Mb.")
        chrom = chromosome.removeprefix("chr")
        region = f"{chrom}:{start}-{end}"
        url = f"https://rest.ensembl.org/overlap/region/homo_sapiens/{quote(region, safe=':-')}"
        response = self._get(
            url,
            [("feature", "gene"), ("feature", "regulatory")],
            provider="Ensembl",
        )
        features = []
        for item in self._json(response, "Ensembl"):
            features.append(
                {
                    "id": item.get("id"),
                    "feature_type": item.get("feature_type") or item.get("biotype") or "feature",
                    "symbol": item.get("external_name") or item.get("gene_name"),
                    "biotype": item.get("biotype"),
                    "start": item.get("start"),
                    "end": item.get("end"),
                    "strand": item.get("strand"),
                    "gene_id": item.get("Parent") or item.get("gene_id"),
                    "source": "Ensembl",
                }
            )
        return {
            "region": {"chromosome": chromosome, "start": start, "end": end, "assembly": "GRCh38"},
            "features": features,
            "source_url": response.url.__str__(),
        }

    def variant_evidence(self, rsid: str) -> dict:
        if not rsid.lower().startswith("rs") or not rsid[2:].isdigit():
            raise ValueError("Enter an rsID such as rs12345.")
        url = f"https://rest.ensembl.org/variation/homo_sapiens/{quote(rsid.lower(), safe='')}"
        response = self._get(url, {"phenotypes": "1", "pops": "1"}, provider="Ensembl Variation")
        data = self._json(response, "Ensembl Variation")
        return {
            "rsid": rsid.lower(),
            "mappings": data.get("mappings", []),
            "phenotypes": data.get("phenotypes", []),
            "populations": data.get("populations", []),
            "source": "Ensembl Variation",
            "source_url": response.url.__str__(),
        }

    def gwas_associations(self, rsid: str) -> dict:
        if not rsid.lower().startswith("rs") or not rsid[2:].isdigit():
            raise ValueError("Enter an rsID such as rs12345.")
        response = self._get(
            "https://www.ebi.ac.uk/gwas/rest/api/v2/associations",
            {"rs_id": rsid.lower(), "size": "100"},
            provider="GWAS Catalog",
        )
        payload = self._json(response, "GWAS Catalog")
        records = payload.get("_embedded", {}).get("associations", payload.get("associations", []))
        return {
            "rsid": rsid.lower(),
            "associations": records,
            "source": "GWAS Catalog API v2",
            "source_url": response.url.__str__(),
        }

    def clinvar(self, rsid: str) -> dict:
        if not rsid.lower().startswith("rs") or not rsid[2:].isdigit():
            raise ValueError("Enter an rsID such as rs12345.")
        term = f'{rsid.lower()}[External Allele ID]'
        search = self._get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            {"db": "clinvar", "term": term, "retmode": "json", "retmax": "20"},
            provider="NCBI ClinVar",
        )
        search = self._json(search, "NCBI ClinVar")
        ids = search.get("esearchresult", {}).get("idlist", [])
        summaries = {}
        summary_url = None
        if ids:
            response = self._get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                {"db": "clinvar", "id": ",".join(ids), "retmode": "json"},
                provider="NCBI ClinVar",
            )
            summaries = self._json(response, "NCBI ClinVar").get("result", {})
            summary_url = str(response.url)
        return {
            "rsid": rsid.lower(),
            "records": [summaries.get(identifier) for identifier in ids if identifier in summaries],
            "source": "NCBI ClinVar E-utilities",
            "search_url": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            "summary_url": summary_url,
            "note": "ClinVar contains submitted assertions and is not a clinical interpretation.",
        }

    def gtex_eqtls(self, rsid: str) -> dict:
        if not rsid.lower().startswith("rs") or not rsid[2:].isdigit():
            raise ValueError("Enter an rsID such as rs12345.")
        response = self._get(
            "https://gtexportal.org/api/v2/association/singleTissueEqtl",
            {"snpId": rsid.lower(), "datasetId": "gtex_v10", "page": "0", "itemsPerPage": "100"},
            provider="GTEx",
        )
        return {
            "rsid": rsid.lower(),
            "dataset": "gtex_v10",
            "associations": self._json(response, "GTEx"),
            "source": "GTEx Portal API v2",
            "source_url": str(response.url),
        }

    def gnomad_variant(self, chromosome: str, position: int, reference: str, alternate: str) -> dict:
        if position < 1 or len(reference) != 1 or len(alternate) != 1:
            raise ValueError("gnomAD lookup currently accepts positive-position GRCh38 SNVs only.")
        if reference.upper() not in "ACGT" or alternate.upper() not in "ACGT" or reference.upper() == alternate.upper():
            raise ValueError("Provide different A/C/G/T reference and alternate alleles.")
        chrom = chromosome.removeprefix("chr").upper()
        if chrom == "MT":
            chrom = "M"
        if chrom not in {str(value) for value in range(1, 23)} | {"X", "Y", "M"}:
            raise ValueError("Enter a valid human chromosome.")
        variant_id = f"{chrom}-{position}-{reference.upper()}-{alternate.upper()}"
        query = """
        query VariantFrequency($variantId: String!) {
          variant(variantId: $variantId, dataset: gnomad_r4) {
            variantId
            genome { ac an }
            exome { ac an }
          }
        }
        """
        try:
            response = httpx.post(
                "https://gnomad.broadinstitute.org/api",
                json={"query": query, "variables": {"variantId": variant_id}},
                headers={**self.headers, "Content-Type": "application/json"},
                timeout=self.timeout,
            )
            response.raise_for_status()
            payload = response.json()
        except httpx.TimeoutException as exc:
            raise EvidenceProviderError(
                "gnomAD did not respond before the timeout. Please retry in a moment."
            ) from exc
        except httpx.HTTPStatusError as exc:
            raise EvidenceProviderError(
                f"gnomAD returned HTTP {exc.response.status_code}. Please retry later."
            ) from exc
        except httpx.RequestError as exc:
            raise EvidenceProviderError(
                "Could not connect to gnomAD. Please retry in a moment."
            ) from exc
        except ValueError as exc:
            raise EvidenceProviderError(
                "gnomAD returned a response that could not be read. Please retry later."
            ) from exc
        if payload.get("errors"):
            raise EvidenceProviderError("gnomAD returned a query error for this variant.")
        result = payload.get("data", {}).get("variant")
        if isinstance(result, dict):
            for cohort in ("genome", "exome"):
                counts = result.get(cohort)
                if isinstance(counts, dict):
                    ac, an = counts.get("ac"), counts.get("an")
                    counts["af_calculated"] = ac / an if isinstance(ac, (int, float)) and isinstance(an, (int, float)) and an > 0 else None
        return {
            "variant": variant_id,
            "dataset": "gnomad_r4 (GRCh38)",
            "result": result,
            "source": "gnomAD GraphQL API",
            "source_url": "https://gnomad.broadinstitute.org/api",
            "note": "Exome and genome allele counts are kept separate; missing data is not evidence of absence.",
        }

    def encode_region(self, chromosome: str, start: int, end: int) -> dict:
        if start < 1 or end <= start or end - start > 1_000_000:
            raise ValueError("ENCODE regions must be ordered, positive, and no wider than 1 Mb.")
        base_url = os.getenv("ENCODE_GDS_BASE_URL", "").rstrip("/")
        if not base_url:
            raise EvidenceProviderNotConfigured(
                "ENCODE region lookup needs an ENCODE genomic-data-service base URL. "
                "Set ENCODE_GDS_BASE_URL in backend/.env to a configured ENCODE-DCC genomic-data-service instance."
            )
        chrom = chromosome.removeprefix("chr")
        params = {
            "query": f"chr{chrom}:{start}-{end}",
            "assembly": "GRCh38",
            "format": "json",
            "page": "1",
            "limit": "100",
            "expand": "0",
            "interval": "intersects",
        }
        response = self._get(f"{base_url}/region-search/", params, provider="ENCODE")
        return {
            "region": {"chromosome": f"chr{chrom}", "start": start, "end": end, "assembly": "GRCh38"},
            "encode": self._json(response, "ENCODE"),
            "source": "ENCODE DCC genomic-data-service region search",
            "source_url": str(response.url),
        }
