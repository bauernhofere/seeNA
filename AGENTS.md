# seeNA contributor and coding-agent guide

## Scope and architecture

Read existing ichorCNA output; do not call CNAs, automatically choose a fitted solution, infer
patient pairing, or import unrelated clinical spreadsheets. User-supplied manifest
annotations are supported but are not automatically safe to share.

`files -> parsers -> validated sample -> cohort -> coverage-aware matrix -> plots`

- `R/read.R`: exact, single-sample schemas; normalize without guessing.
- `R/sample.R`, `R/bounds.R`: source identity, interval validation, audited window padding.
- `R/accessors.R`: public parameters, upstream display transform, NEUT-bin evidence.
- `R/provenance.R`: immutable import-time fingerprints; paths opt-in.
- `R/genome.R`: coordinate reference and region parsing.
- `R/cohort.R`: atomic manifest import, sample/metadata order validation.
- `R/matrix.R`: bin/segment continuous means, categorical modes, coverage and mixed flags.
- `R/plot.R`, `R/heatmap.R`: standard plot objects, no file writes.
- `R/concordance.R`, `R/plot-concordance.R`: explicit-pair directional comparison,
  separate call/height evidence and sex-reference flags; no inferred pairing.
- `docs/data-contracts.md`: scientific semantics and compatibility limits.
- `docs/decision-register.md`: rationale, sources, alternatives and open study gates.
- `README.Rmd`: executable source for README.md and the four fixture-only pictures.
- `tests/testthat/test-adversarial.R`: regressions for independent-review findings.

## Non-negotiable scientific rules

- Coordinates: 1-based closed; normalized chromosomes 1-22, X, Y; explicit build.
- Missing calls/flags are not neutral/false. Unknown tokens fail explicitly.
- Internal TF is a fraction in [0,1]; percent conversion is presentation only.
- Raw event and corrected_call are distinct. Never silently substitute them.
- Call codes are categories: no arithmetic averaging or silent tie-breaking.
- Coverage and heterogeneous-call flags travel with the matrix.
- `segment_median` coverage is finite segment-span coverage, not bin/read coverage.
  Never overwrite source bin logR with segment medians.
- Match source component IDs before applying aliases. Identical IDs do not prove
  the files came from the same fitted run; that remains an input-selection duty.
- Inputs are read-only. Default window-padding handling requires regular source
  grid evidence; broader trimming is opt-in. Neither policy proves assembly identity.
- Adjusted logR is an explicit display transform, not a new call or guarantee
  that every loss/gain has negative/positive y. Raw output stays the default.
- NEUT-bin CN evidence never falls back across chromosome classes or to CN 2.
- No patient-specific code patches. Fix source data or selection manifests.

## Privacy and testing

Do not commit clinical files, sample identifiers, local paths, credentials, or
rendered patient figures. Real data stay outside the repository. Existing tiny
format fixtures are test-only examples, not biological validation. New malformed
or boundary-case test inputs must be clearly test-only, never manuscript data.
No test or example may require private data, credentials, or a network connection.
`make readme` regenerates only the four allowlisted fixture PNGs in docs/figures;
never replace them with patient-derived images. `make audit-reference` is an
explicit optional network audit of public UCSC chromosome tables, not a test.
Path redaction and identifier-pattern tests do not de-identify genomic data.
`.gitignore` allows only the named public fixtures; do not widen its exceptions
for real inputs. Ignoring a file does not remove it from Git history. Audit
tracked files and history before publishing new data or changing visibility.
Keep this guide tracked; `.Rbuildignore` excludes it from the installed package.

## Change checklist

1. Write a failing regression tied to an invariant or verified upstream format.
2. Make a focused change; avoid compatibility shims during this pre-release stage.
3. Update contracts, examples and NEWS for scientific/API changes.
4. Run `make test`, then `make check`. Generated `man/` and `NAMESPACE` come from
   roxygen: do not edit them manually.
5. Test actual rendering for plot changes; plot-object class tests are insufficient.
6. Keep real-data checks local and report aggregates, never identifying data.
7. Change README.Rmd, then run `make readme`; inspect the generated figures.
   Keep source/markdown/images together. For scientific/default changes, update
   the decision ID and methods, distinguishing policy from external evidence.

`make bootstrap` installs development dependencies (network required).
`make check` generates docs, tests, and performs R CMD check. Passing it is
necessary, not evidence of scientific correctness. Manuscript release also needs
real-input reconciliation, visual review, a pinned commit, and author approval.
Public repository visibility does not approve manuscript conclusions or a
software release. Do not tag a release or create an archive/DOI without author
approval. Review reports are hypotheses to verify, not authority.
