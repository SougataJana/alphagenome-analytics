# Atlas website parity validation

## Status

**Pending.** No Atlas website values have been entered as expected results. Do not treat the synthetic adapter tests or one SDK schema inspection as Atlas UI parity.

## Capture protocol

For each case, record the exact GRCh38 chromosome, 1-based position, REF and ALT, query date/time (UTC), Atlas UI URL or page identity, requested/displayed scorers, and any UI filters. Save a screenshot for context and transcribe the underlying values shown by Atlas. Do not estimate values from graph pixels. Keep provider display rounding separate from the underlying value when the UI exposes both.

Use multiple cases selected before comparing results:

1. A published example with a non-empty RNA_SEQ matrix.
2. A variant with at least one negative and one positive signed quantile.
3. A variant with non-empty SPLICE_SITES and ATAC output, if the Atlas UI exposes those scorers.
4. A low-signal or sparse-output case to check missing scorers and empty results.

For every case compare, before UI formatting:

- AVI raw value and quantile, including the displayed scale.
- Scorer names and matrix dimensions.
- Raw scorer matrix values and signed quantiles (sign and magnitude).
- Feature-attribution labels and values.
- Row/column metadata used for gene, tissue, and feature labels.
- Missing-value handling, ordering, and any truncation or filters.

Do not compare values from different model or data releases as if they were the same run. Record the SDK version and Atlas model/version metadata when available. If the UI does not reveal a field, mark it `not exposed` rather than inferring it.

## Case log

| Case | Exact input | Atlas retrieval time (UTC) | Atlas fields captured | AGA analysis ID | Result | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | — | — | — | — | Pending | — |
| 2 | — | — | — | — | Pending | — |
| 3 | — | — | — | — | Pending | — |
| 4 | — | — | — | — | Pending | — |

## Interpretation

Parity means agreement in faithfully representing the same provider output under matched inputs and versions. It does not validate AlphaGenome's biological accuracy, gene causality, clinical value, or the tool's added statistical methods. Report any disagreement, exclusions, missing data, and rounding tolerance. Do not deploy or make an accuracy claim while parity remains pending.
