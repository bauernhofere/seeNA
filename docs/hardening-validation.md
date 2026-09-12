# Hardening validation — 2026-09-12

Scope: private development version 0.0.0.9001; not a public release or claim of
clinical/scientific validation of ichorCNA itself. See
[downstream-validation.md](downstream-validation.md) for subsequent API changes,
new evidence, and the diagnosis of CI failures not seen in this local run.

## Automated evidence

- 35 test cases, 145 passing assertions, no failures/errors in the local run.
- `make check`: 0 errors, 0 warnings; one environment-only NOTE:
  `unable to verify current time`.
- Examples execute in R CMD check. Methods and upstream notice ship in `inst/`.
- Actual PDF heatmap and headless ragg PNG rendering are exercised, including
  rasterized genome-scale heatmaps and ggplot profiles.
- Regression cases cover malformed/missing fields, sample mismatches, duplicate
  semantics, both segment schemas, appended parameter diagnostics, terminal
  clipping, identity aliases, file drift, worker failure, allocation limits,
  categorical mode/ties, partial coverage, matrix metadata drift and ordering.

## Real-input compatibility and smoke check (local only)

- 1,272 `.cna.seg`, 1,272 `.seg`, and 1,272 `.params.txt` files parsed with the
  hardened readers, with no parse failures. These include repeated runs and are
  **not** 1,272 distinct patients.
- All 71 samples in the supplied canonical directory imported with explicit
  `bounds='trim'`; 213 terminal-interval changes were recorded.
- A 71 x 3,102 call matrix was constructed with full-coverage requirement.
  Values were checked to be only the documented codes or NA.
- Reference local run: reading ~2.0 seconds, matrix construction ~0.75 seconds;
  returned matrix object ~5.8 MB including coverage/mixed layers and provenance.
  These are warm local smoke timings, not a general performance guarantee.
- Real-input profile, paired overlay, chromosome-8 region, and cohort heatmap
  were rendered and inspected. These are illustrative smoke plots, not the final
  manuscript figure set. Clinical inputs, selection records and output images
  remain outside this repository.

Local environment: R 4.4.2, Apple Silicon macOS. CI has been expanded to Windows,
macOS, Linux and multiple R versions; remote results must be checked separately.
The minimum declared R 4.1 is not independently verified by this local run.

## Still required before a manuscript release

1. Explicitly freeze the manuscript's source-run and sample-pair manifest.
2. Reconcile expected/plotted samples; do not silently drop failed pairs.
3. Reproduce the existing CN-minus-ploidy heatmap with explicit display thresholds,
   ordering and saturation (see review-triage.md); a categorical call plot is not
   an equivalent substitute.
4. Compare package-driven final figures/matrices with source results and obtain
   author approval of the measurement and presentation.
5. Confirm study-code copyright/contributor attribution before public distribution.
6. Obtain successful cross-platform CI, then pin a reviewed version/commit and
   publish/archive only with author approval. No DOI or release date is claimed.
