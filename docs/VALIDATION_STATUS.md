# Validation status

Last updated: 2026-10-09

This is an implementation-validation record, not a biological performance claim.

| Check | Status | Evidence and scope |
| --- | --- | --- |
| GRCh38 input syntax, chromosome aliases, SNV-only restriction, and coordinate bounds | Passed | Backend unit tests cover accepted aliases, invalid chromosomes/alleles, and first/last chromosome positions. |
| REF allele verification behavior | Passed in unit tests; live provider check pending | Mocked Ensembl sequence responses cover a matching REF, mismatch rejection, malformed response, provider failure, and explicit GRCh38 coordinate-system parameter. Live reference access has not been confirmed from this environment. |
| Atlas request blocked on failed REF validation | Passed | Unit test confirms the Atlas client is not constructed when reference verification fails. |
| AlphaGenome SDK adapter | Passed for synthetic fixture; schema reconnaissance only for real response | Adapter test exercises AVI, signed quantile, feature attribution, labels, and preservation. The one real SDK response summarized in `ATLAS_RESPONSE_SCHEMA.md` is not an Atlas UI parity test. |
| Evidence provider adapters | Passed for mocked provider responses | Ensembl, gnomAD, GWAS Catalog, ClinVar, GTEx, and SCREEN tests check response mapping, source labels, and independent records. Live provider availability and current response contracts remain unverified. |
| Statistical calculations | Passed for reference examples | Tests check a known BH vector, exact small-sample hypergeometric examples, and permutation output against an exhaustive allocation reference. This does not validate biological study design or the choice of input units/background. |
| Atlas website parity | Pending | No Atlas UI values have been captured or compared. See `ATLAS_PARITY_CHECKLIST.md`. |
| Biological benchmark or predictive performance | Pending | No curated variant benchmark, ClinVar evaluation, or external user evaluation is included. |
| Frontend production build | Passed | `npm run build` completed successfully. |

## Current automated checks

- `python -m unittest discover -s tests -v` in `backend/`: 34 tests passed.
- `npm run build` in `frontend/`: passed.
- `git diff --check`: passed.

The live Ensembl-to-Atlas probe could not be completed from the restricted command environment because outbound DNS/network access was unavailable before the Atlas call. The Atlas website also requires the account holder to complete a Google verification step in the browser. Do not describe UI parity, live provider connectivity, or biological validity as passed until those checks are performed and the case log is filled with observed values.
