# AlphaGenome Analytics roadmap

The product makes AlphaGenome predictions searchable, interpretable, comparable, and biologically useful. It should expose model outputs and evidence transparently rather than create another opaque score.

## Product feature groups

1. **Explore** — variant, gene, and region entry points; genomic context and candidate regulatory links.
2. **Variant Effect Analyzer** — expression, splicing, chromatin/regulatory effects, tissue context, and an overall effect profile.
3. **Prioritization and Dataset Analysis** — transparent evidence-backed rankings; VCF processing, distributions, heatmaps, and mechanism clustering.
4. **Variant-to-Gene and Tissue Discovery** — candidate target genes and cell/tissue context from model predictions and external regulatory evidence.
5. **Research Workflow** — statistical comparisons, evidence integration, regulatory networks, deterministic Ask AGA orchestration, reports, and reproducibility.

## Delivery stages

### Implemented in this build

- Explore a GRCh38 single-nucleotide variant by `chr:position REF>ALT`.
- Retrieve AlphaGenome Atlas AVI and feature-attribution results.
- Visualize the returned effect profile.
- Navigate a basic coordinate window around the variant.

- Single-SNV AlphaGenome Atlas query and effect profile, including returned AVI, feature attribution, scorer rows, genes, tissues, and quantiles.
- Local VCF parsing and background batch analysis (up to 5,000 GRCh38 SNVs per job), with explicit consent before coordinates are sent, job progress, raw AVI ranking, and a relative batch magnitude display.
- Ensembl gene, region, and rsID lookups; GWAS Catalog association, ClinVar, and GTEx single-tissue eQTL lookups by rsID; gnomAD GRCh38 exome/genome allele-count lookup for a loaded SNV; ENCODE public portal region search. Each response includes source information.
- Tissue summaries, annotated prediction gene labels, a feature-to-gene co-occurrence view, score-group permutation/bootstrap analysis, and supplied gene-set enrichment with Benjamini-Hochberg correction.
- Combined ranking can include pasted GWAS, GTEx, and gnomAD values with user-selected weights, within-dataset percentile normalization, and downloaded raw inputs/method. It does not auto-join those sources across batch variants.
- Deterministic Ask AGA summaries over loaded results, JSON/HTML exports, and a local SQLite history of manifests/results.

### Data integrations and limitations

- ENCODE region search calls the public ENCODE portal and is subject to its availability and request limits. It returns overlapping records, not a per-feature genomic annotation track. The current implementation has not been verified against a live ENCODE response.
- GTEx eQTL and ClinVar are individual rsID lookups; they are not automatically merged into batch rankings. Combined ranking uses a user-provided keyed score table; ClinVar assertions are displayed separately.
- The candidate-gene view reports AlphaGenome-returned labels and Ensembl locus annotations. It does not infer a causal target gene.
- The network shows feature/gene co-occurrence from returned predictions, not experimentally established regulatory edges.
- The region view has a coordinate navigator and Ensembl annotations; richer prediction/evidence multi-track visualization remains planned. Batch plots show relative AVI magnitude and descriptive dominant-scorer groups, not statistical clustering or inferred mechanisms.
- VCF parsing supports uncompressed SNVs only. Genotype columns are discarded by the browser parser. Variant coordinates are sent to the configured AlphaGenome Atlas service after consent and saved in the local SQLite database.
- Live provider requests were not made as part of this implementation. Availability depends on the saved AlphaGenome key, account permissions, network access, and installed dependencies.

## Scientific and privacy rules

- Do not show fabricated predictions, gene links, pathway enrichments, or statistical significance.
- Keep calculations in deterministic analytical code; any language model may organize calls and explain returned outputs, but must not invent results.
- Make variant transmission to external services explicit, especially for VCFs or other potentially identifiable genetic data.
- Label AlphaGenome outputs as research predictions, not clinical evidence.
