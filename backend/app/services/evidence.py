"""Read-only public evidence adapters. Every result carries its source."""

from urllib.parse import quote
from pathlib import Path
import os

import httpx
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from dotenv import load_dotenv


BACKEND_ROOT = Path(__file__).resolve().parents[2]
load_dotenv(BACKEND_ROOT / ".env")


class EvidenceProviderError(Exception):
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
            raise ValueError("SCREEN regions must be ordered, positive, and no wider than 1 Mb.")
        chrom = chromosome.removeprefix("chr")
        if chrom.upper() not in {*(str(n) for n in range(1, 23)), "X", "Y", "M", "MT"}:
            raise ValueError("Enter a human chromosome from chr1–chr22, chrX, chrY, or chrM.")
        api_key = os.getenv("SCREEN_API_KEY", "").strip()
        if not api_key:
            raise EvidenceProviderError(
                "SCREEN region search is not configured. Add SCREEN_API_KEY to backend/.env locally or to the backend service environment in Render."
            )
        query = f'''query RegionCcreSearch {{
          cCRESCREENSearch(
            assembly: "grch38"
            coordinates: [{{ chromosome: "chr{chrom}", start: {start}, end: {end} }}]
          ) {{
            chrom start len pct ctcf_zscore dnase_zscore atac_zscore enhancer_zscore promoter_zscore
            info {{ accession }}
          }}
        }}'''
        url = "https://screen.api.wenglab.org/graphql"
        try:
            response = httpx.post(
                url,
                json={"query": query},
                headers={**self.headers, "Authorization": f"Bearer {api_key}"},
                timeout=self.timeout,
            )
            response.raise_for_status()
            payload = response.json()
        except httpx.TimeoutException as exc:
            raise EvidenceProviderError("SCREEN did not respond before the timeout. Please retry in a moment.") from exc
        except httpx.HTTPStatusError as exc:
            if exc.response.status_code in {401, 403}:
                detail = "SCREEN rejected the API key. Check that SCREEN_API_KEY is current and active."
            else:
                detail = f"SCREEN returned HTTP {exc.response.status_code}. Please retry later."
            raise EvidenceProviderError(detail) from exc
        except httpx.RequestError as exc:
            raise EvidenceProviderError("Could not connect to SCREEN. Please retry in a moment.") from exc
        except ValueError as exc:
            raise EvidenceProviderError("SCREEN returned a response that could not be read as JSON.") from exc
        if not isinstance(payload, dict):
            raise EvidenceProviderError("SCREEN returned an unexpected response format.")
        if payload.get("errors"):
            raise EvidenceProviderError("SCREEN could not complete this cCRE query. Check the interval and retry.")
        data = payload.get("data")
        ccres = data.get("cCRESCREENSearch", []) if isinstance(data, dict) else []
        if not isinstance(ccres, list):
            raise EvidenceProviderError("SCREEN returned an unexpected cCRE result format.")
        return {
            "region": {"chromosome": f"chr{chrom}", "start": start, "end": end, "assembly": "GRCh38"},
            "ccres": ccres,
            "screen": payload,
            "source": "SCREEN GraphQL API (ENCODE Registry of cCREs)",
            "source_url": url,
        }

    def synthesize_variant(self, chromosome: str, position: int, reference: str, alternate: str, rsid: str | None = None) -> dict:
        """Gather independent provider records without inferring a consensus classification."""
        calls = {
            "Ensembl region": lambda: self.region(chromosome, max(1, position - 50_000), position + 50_000),
            "gnomAD": lambda: self.gnomad_variant(chromosome, position, reference, alternate),
        }
        if rsid:
            if not rsid.lower().startswith("rs") or not rsid[2:].isdigit():
                raise ValueError("Enter an rsID such as rs12345, or leave it blank.")
            calls.update({
                "Ensembl Variation": lambda: self.variant_evidence(rsid),
                "GWAS Catalog": lambda: self.gwas_associations(rsid),
                "ClinVar": lambda: self.clinvar(rsid),
                "GTEx": lambda: self.gtex_eqtls(rsid),
            })

        providers: dict[str, dict] = {}
        with ThreadPoolExecutor(max_workers=min(6, len(calls))) as pool:
            futures = {pool.submit(call): name for name, call in calls.items()}
            for future in as_completed(futures):
                name = futures[future]
                try:
                    providers[name] = {"status": "complete", "retrieved_at": datetime.now(timezone.utc).isoformat(), "result": future.result()}
                except (EvidenceProviderError, ValueError) as exc:
                    providers[name] = {"status": "unavailable", "retrieved_at": datetime.now(timezone.utc).isoformat(), "error": str(exc)}

        return {
            "query": {"chromosome": chromosome, "position": position, "reference": reference, "alternate": alternate, "assembly": "GRCh38", "rsid": rsid.lower() if rsid else None},
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "providers": providers,
            "interpretation": "Provider records are shown independently. An rsID query is not assumed to identify the supplied coordinates; verify variant mapping before combining them. No consensus or clinical classification is inferred.",
        }

    def prioritize_variant_genes(self, chromosome: str, position: int, effects: list[dict], window: int = 100_000) -> dict:
        if window < 1 or window > 1_000_000:
            raise ValueError("Gene-prioritization window must be between 1 bp and 1 Mb.")
        region = self.region(chromosome, max(1, position - window), position + window)
        candidates: dict[str, dict] = {}
        for feature in region["features"]:
            if feature.get("feature_type") != "gene" and not str(feature.get("id", "")).startswith("ENSG"):
                continue
            symbol = feature.get("symbol") or feature.get("gene_id") or feature.get("id")
            if not symbol:
                continue
            start, end = feature.get("start"), feature.get("end")
            distance = 0 if isinstance(start, int) and isinstance(end, int) and start <= position <= end else (
                min(abs(position - start), abs(position - end)) if isinstance(start, int) and isinstance(end, int) else None
            )
            candidates[str(symbol).casefold()] = {
                "gene": str(symbol), "ensembl_id": feature.get("id") or feature.get("gene_id"),
                "start": start, "end": end, "distance_bp": distance,
                "effect_count": 0, "modalities": [], "best_effect_quantile": None,
            }

        for effect in effects:
            gene = effect.get("gene")
            if not gene:
                continue
            key = str(gene).casefold()
            candidate = candidates.setdefault(key, {
                "gene": str(gene), "ensembl_id": None, "start": None, "end": None,
                "distance_bp": None, "effect_count": 0, "modalities": [], "best_effect_quantile": None,
            })
            candidate["effect_count"] += 1
            modality = effect.get("modality")
            if modality and modality not in candidate["modalities"]:
                candidate["modalities"].append(modality)
            quantile = effect.get("quantile_score")
            if isinstance(quantile, (int, float)):
                current = candidate["best_effect_quantile"]
                candidate["best_effect_quantile"] = max(current, quantile) if current is not None else quantile

        max_effect_count = max((item["effect_count"] for item in candidates.values()), default=0)
        rows = list(candidates.values())
        for row in rows:
            distance = row["distance_bp"]
            row["proximity_component"] = max(0.0, 1.0 - distance / window) if distance is not None else None
            row["alphagenome_support_component"] = row["effect_count"] / max_effect_count if max_effect_count else None
            components = [value for value in (row["proximity_component"], row["alphagenome_support_component"]) if value is not None]
            row["research_priority_score"] = sum(components) / len(components) if components else None
            row["modalities"].sort()
        rows.sort(key=lambda item: (item["research_priority_score"] is None, -(item["research_priority_score"] or 0), item["gene"].casefold()))
        return {
            "variant": {"chromosome": chromosome, "position": position, "assembly": "GRCh38"},
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "window_bp": window,
            "method": "Equal-weight mean of available components: linear proximity (1 - distance/window) and AlphaGenome effect-record count normalized by the maximum count among returned candidates. Missing components are omitted.",
            "candidates": rows,
            "source": "Ensembl REST overlap/region plus submitted AlphaGenome Atlas effect records",
            "source_url": region.get("source_url"),
            "note": "Research prioritization heuristic, not a probability, causal-gene determination, or clinical classification. Candidates and scores depend on the returned Ensembl interval and AlphaGenome effect labels.",
        }
