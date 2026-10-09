# Methods and scientific scope

This document describes what the software computes and the limits on scientific interpretation. AlphaGenome Analytics is a research interface, not a validated diagnostic or clinical decision tool.

## Inputs and coordinate rules

- Variant requests are GRCh38 single nucleotide substitutions using 1-based chromosome positions. Chromosome aliases are normalized to `chr1`–`chr22`, `chrX`, `chrY`, or `chrM`/`chrMT`.
- The service rejects non-ACGT alleles, identical reference/alternate alleles, and positions outside the corresponding GRCh38 chromosome length. The reference allele is not independently checked against a reference FASTA; users must verify that the supplied allele matches GRCh38.
- VCF parsing is limited to supported SNVs. Indels, structural variants, genotype-based analyses, and liftover are not supported.

## AlphaGenome prediction results

- The backend queries the configured AlphaGenome Atlas service with `AVI_SCORE`, `AVI_SCORE_FEATURE_IMPORTANCE`, `RNA_SEQ`, `SPLICE_SITES`, and `ATAC` scorers.
- The application displays provider-returned raw values and quantiles. It does not rescale an AVI raw score into a probability. A quantile is a reference-distribution rank, not a probability of pathogenicity or causality.
- Effects are sorted by absolute raw score within each scorer and capped at 200 rows per scorer. The response reports the original total and scorer-specific counts, plus whether any scorer was truncated. This ordering is a display choice, not a cross-scorer ranking; raw scores from different scorers may not be comparable.
- Atlas does not expose a model version in the response used by this application. The manifest records that as unavailable instead of inventing a value. The installed SDK version and requested scorer names are recorded where available.
- Prediction labels, feature attributions, and proximity do not establish mechanism or causal target genes.

## Evidence sources

External lookups are run only when requested. Each response should retain its provider name, query, retrieval time, and source URL when the provider supplies them. The evidence pages must distinguish a provider failure from an empty result; an empty response is not evidence that an association or feature is absent.

- Ensembl supplies gene, variant, and interval annotations.
- GWAS Catalog, ClinVar, GTEx, and gnomAD provide distinct evidence types and are not interchangeable. rsID-based lookups may refer to multiple alleles or studies; verify the record's assembly, allele, phenotype, tissue, and study context before interpretation.
- SCREEN supplies ENCODE Registry candidate cis-regulatory elements (cCREs). cCRE overlap is an annotation, not a functional assay for the user's variant. The configured SCREEN endpoint requires a backend API key.
- Evidence is not automatically joined across batch variants. A user-entered combined ranking is keyed by exact variant text and uses only the data supplied by the user.

## Ranking and visualization

- AVI sorting orders the returned variants by the provider's AVI raw score. It is not a calibrated pathogenicity ranking.
- Combined evidence ranking transforms each available input component into a percentile within the submitted dataset, reverses direction for lower-is-stronger components, and computes a user-weighted mean over available components. The resulting percentile is dataset-relative and changes with the submitted rows and weights; it is not an absolute evidence score.
- Variant-to-gene scores are fixed equal-weight heuristics combining linear distance decay and the maximum RNA_SEQ effect quantile per Ensembl-mapped gene. Missing RNA_SEQ support contributes zero under fixed weights; effect-row counts are context only. Unmapped effect labels are shown without rank or score, so missing distance cannot increase their score. The linear distance window and weighting are not empirically calibrated. Scores are dimensionless ranking values, not probabilities, and do not demonstrate causality.
- Tissue summaries and tissue-to-tissue comparisons describe returned model-effect rows. They do not run a hypothesis test, calculate a confidence interval, or treat those rows as independent biological replicates.
- Region displays locate the selected coordinate and show annotations returned by the selected provider. A coordinate line by itself is not an annotation track.

## Statistics

- User-supplied score comparison reports a two-sided permutation test for the difference in means with plus-one correction and a percentile bootstrap confidence interval for the mean difference. The permutation test assumes exchangeability under the null. The bootstrap assumes independent sampling units. These assumptions can fail for repeated measures, related variants, or model outputs; do not use model-effect rows as biological replicates.
- Gene-set over-representation uses a one-sided hypergeometric survival probability and Benjamini–Hochberg adjustment across the gene sets submitted in that run. The selected gene list must be a subset of the explicitly supplied background; otherwise the request is rejected. Gene identifiers are uppercased for exact matching, not alias-mapped; the output reports gene-set members excluded from the background. The result depends on the chosen background and identifier mapping.
- Statistical outputs are not proof of biological mechanism. Analyses should be designed before looking at results, include appropriate controls, and report the tested universe and preprocessing.

## Reproducibility and validation status

Inputs, timestamps, analysis IDs, requested scorers, SDK version when available, evidence-provider provenance, random seed, iteration count, and methods are included in analysis records or exports. The application does not claim that external data snapshots are immutable; provider databases may change over time.

Validation is incomplete. The current repository does not contain a curated benchmark showing agreement with the Atlas website, a ClinVar pathogenic-versus-benign evaluation, or independent user evaluation. The software's backend checks validate input handling and deterministic calculations only; they cannot establish biological validity. Before making performance claims, compare exact known variants against Atlas with assembly/alleles/scorers recorded, evaluate a pre-specified and appropriately controlled variant set, and report coverage, failures, and discrepancies.

No repository code license has been selected yet. Until a license is added, users should not assume that code reuse or redistribution is permitted.

## Terms, privacy, and intended use

AlphaGenome is provided for non-commercial research and theoretical modelling; its outputs must not be used for clinical decision-making or to train machine-learning models. The public app sends submitted variants to its configured backend and then to the provider. Before operating a shared-key public service, the owner must confirm with the provider that this hosted proxy and its intended audience comply with current account/API terms. Do not submit identifiable or sensitive genomic data without an appropriate legal and ethics review. See the [official AlphaGenome repository and terms](https://github.com/google-deepmind/alphagenome).

The GRCh38 coordinate limits in the API are based on the [NCBI Genome Reference Consortium assembly data](https://www.ncbi.nlm.nih.gov/grc/human/data). SCREEN source documentation: [SCREEN API query guide](https://weng-lab.github.io/SCREEN2.0/queries/ccres/).
