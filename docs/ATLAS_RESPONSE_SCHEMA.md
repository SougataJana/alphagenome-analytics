# Observed AlphaGenome Atlas response

This note records one real GRCh38 Atlas query inspected on 2026-10-09 through the official Python SDK (`alphagenome` 0.10.0). It documents the SDK-decoded response for schema discovery; it is not a formal provider guarantee and does not expose or store credentials.

## Query used

- Variant: `chr22:36201698 A>C` (GRCh38)
- Requested scorers: `AVI_SCORE`, `AVI_SCORE_FEATURE_IMPORTANCE`, `RNA_SEQ`, `SPLICE_SITES`, `ATAC`
- Provider model version: not returned by the response inspected

## SDK response structure observed

The SDK returned scorer-keyed AnnData matrices. Each matrix exposes `X`, row metadata (`obs`), column metadata (`var`), and, for score matrices, a `quantiles` layer.

| Scorer | Observed matrix shape | Row metadata examples | Column metadata examples | Quantiles layer |
| --- | ---: | --- | --- | --- |
| `AVI_SCORE` | 1 × 1 | `variant`, `interval` | `name`, `strand` | yes |
| `AVI_SCORE_FEATURE_IMPORTANCE` | 1 × 18 | `variant`, `interval` | `name`, `strand` | no |
| `RNA_SEQ` | 37 × 371 | `gene_id`, `gene_name`, `strand`, `variant`, `interval` | `name`, `strand`, `Assay title`, `ontology_curie`, `biosample_name`, `biosample_type`, `biosample_life_stage`, `gtex_tissue`, `data_source`, `endedness`, `genetically_modified`, `nonzero_mean` | yes |
| `SPLICE_SITES` | 1 × 2 | `gene_id`, `gene_name`, `strand`, `variant`, `interval` | `name`, `strand` | yes |
| `ATAC` | 1 × 167 | `variant`, `interval` | `name`, `strand`, `Assay title`, `ontology_curie`, `biosample_name`, `biosample_type`, `biosample_life_stage`, `data_source`, `endedness`, `genetically_modified`, `nonzero_mean` | yes |

Shapes and available metadata can vary by variant and scorer. The application records the requested scorer names and SDK version; Atlas did not return a model version in this response.

## Quantile interpretation

In this response, `AVI_SCORE` returned raw score `1.5548` and quantile `0.99894`. Signed scorer quantiles included values near `−1` and positive values near `+1` (for example, RNA_SEQ values near `−0.99998`). Therefore, effect quantiles cannot be constrained to `[0, 1]`: for signed scorers, the sign retains direction, and absolute value represents extremeness. The application must preserve both. AVI quantile remains a non-negative tail quantile.

## Preservation status and limitations

The standard result keeps AVI, feature-attribution values, selected gene/tissue labels, raw effect scores, and signed effect quantiles. Its displayed effects are sorted and capped at 200 per scorer while retaining total counts. The complete SDK-decoded scorer matrices (`X`, `obs`, `var`, and all named layers) are now stored separately as a compressed provider-response artifact and can be retrieved from `/api/v1/analyses/{analysis_id}/atlas-response`. This artifact is not the original HTTP wire payload: the SDK does not expose that through this call. Analyses saved with persistence disabled will not have a retrievable artifact.

This single query is schema reconnaissance, not a benchmark or biological validation. The example variant is known from the AlphaGenome publication; values and annotations should still be verified against the current Atlas website before being used as a parity test. The automated adapter test uses synthetic matrices to exercise conversion and preservation; it does not count as Atlas website parity. A parity fixture must be captured independently from the Atlas UI and include the exact variant, assembly, scorer names, raw values, quantiles, dimensions, metadata, and retrieval date before a comparison test can make a parity claim.
