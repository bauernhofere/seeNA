# ichorViz contributor and coding-agent guide

## Scope and architecture

Read existing ichorCNA output; do not call CNAs, choose a fitted solution, infer
patient pairing, or import unrelated clinical spreadsheets. User-supplied manifest
annotations are supported but are not automatically safe to share.

`files -> parsers -> validated sample -> cohort -> coverage-aware matrix -> plots`

- `R/read.R`: exact, single-sample schemas; normalize without guessing.
- `R/sample.R`: source identity checks, interval/value validation, explicit bounds.
- `R/provenance.R`: immutable import-time fingerprints; paths opt-in.
- `R/genome.R`: coordinate reference and region parsing.
- `R/cohort.R`: atomic manifest import, sample/metadata order validation.
- `R/matrix.R`: continuous means, categorical modes, coverage and mixed flags.
- `R/plot.R`, `R/heatmap.R`: standard plot objects, no file writes.
- `docs/data-contracts.md`: scientific semantics and compatibility limits.
- `tests/testthat/test-adversarial.R`: regressions for independent-review findings.

## Non-negotiable scientific rules

- Coordinates: 1-based closed; normalized chromosomes 1-22, X, Y; explicit build.
- Missing calls/flags are not neutral/false. Unknown tokens fail explicitly.
- Internal TF is a fraction in [0,1]; percent conversion is presentation only.
- Raw event and corrected_call are distinct. Never silently substitute them.
- Call codes are categories: no arithmetic averaging or silent tie-breaking.
- Coverage and heterogeneous-call flags travel with the matrix.
- Match source component IDs before applying aliases. Identical IDs do not prove
  the files came from the same fitted run; that remains an input-selection duty.
- Inputs are read-only; bounds trimming is opt-in, auditable and not liftover.
- No patient-specific code patches. Fix source data or selection manifests.

## Privacy and testing

Do not commit clinical files, sample identifiers, local paths, credentials, or
rendered patient figures. Real data stay outside the repository. Existing tiny
format fixtures are test-only examples, not biological validation. New malformed
or boundary-case test inputs must be clearly test-only, never manuscript data.
No test or example may require private data, credentials, or a network connection.
Path redaction and identifier-pattern tests do not de-identify genomic data.

## Change checklist

1. Write a failing regression tied to an invariant or verified upstream format.
2. Make a focused change; avoid compatibility shims during this pre-release stage.
3. Update contracts, examples and NEWS for scientific/API changes.
4. Run `make test`, then `make check`. Generated `man/` and `NAMESPACE` come from
   roxygen: do not edit them manually.
5. Test actual rendering for plot changes; plot-object class tests are insufficient.
6. Keep real-data checks local and report aggregates, never identifying data.

`make bootstrap` installs development dependencies (network required).
`make check` generates docs, tests, and performs R CMD check. Passing it is
necessary, not evidence of scientific correctness. Manuscript release also needs
real-input reconciliation, visual review, a pinned commit, and author approval.
Do not publish the private repository, tag a release, or mint citation metadata
without that approval. Review reports are hypotheses to verify, not authority.
