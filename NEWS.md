# ichorViz 0.0.0.9003 (development)

* Added explicit-pair directional call comparison and a stacked agreement plot:
  shared grid, independent call/logR coverage, unknown mixed/missing calls,
  baseline-evidence flags and optional conservative sex-chromosome gating.
* Paired plots retain both altered heights by default; the manuscript's
  representative-height policy is an explicit option, not a difference metric.
  Added arithmetic, rendering and vdiffr regressions plus a fixture-only example.
* Sample identifiers are now shown by default in every plot, with an explicit
  opt-out: `show_sample_id = FALSE` for profile, comparison and region plots and
  `show_row_names = FALSE` for heatmaps (whose default changed from FALSE to TRUE).
* Replaced the custom `%||%` operator, whose NA/empty semantics differed from
  base R, with an internal NULL-only helper; empty and NA arguments are now
  handled explicitly (an NA `sample_id` is an error).
* Added a runnable example to every exported function and rewrote roxygen
  descriptions to lead with what each function does; limits moved to `@details`
  and the installed methods contract.
* Restructured the README (install, quick start, comparison, cohort) and trimmed
  repeated caveats in the README and vignette.
* Removed hard-coded commit pins from the README and docs; the install snippet
  now tracks the branch and records `RemoteSha`.
* Added vdiffr snapshot tests for the ggplot outputs (`vdiffr` in Suggests).

# ichorViz 0.0.0.9002

* Added a source/evidence-linked decision register, distinguishing upstream
  contracts, package policies and unresolved live-manuscript assumptions.
* Added an executable illustrated README with three fixture-only figures, a
  CI execution step, and an optional UCSC reference-length audit.
* Documented the ordinal-distance assumption of optional call-code clustering.
* Added default bounds='window' for source-grid-supported terminal padding;
  arbitrary overhangs still fail. Window acceptance does not verify genome build.
* Aggregated coordinate notices at cohort level, without blanket warning suppression.
* Exported TF/ploidy accessors, call-state mapping, reference genome layout,
  upstream-formula adjusted logR, and diagnostic neutral-CN evidence.
* Added explicit ploidy_adjust to profile, comparison and region plots. Raw
  output remains the default; transforms never modify source values or calls.
* Added an executable workflow vignette and an auditable run-selection recipe.
* Fixed cross-platform rendering tests to close PDF devices before checking bytes.
* Moved development tooling out of Suggests; version remains pre-release.

# ichorViz 0.0.0.9001

* Strict numeric/flag/call parsing preserves missingness and rejects ambiguity.
* Reconciled source identities before aliasing, validated segments and bounds,
  and added explicit audited terminal-interval trimming.
* Captured immutable import fingerprints with opt-in paths.
* Made cohort imports atomic, including worker failures.
* Required explicit matrix values; categorical rebinning now uses base-pair mode
  with NA ties, coverage thresholds, and heterogeneous-call flags.
* Added cohort/matrix validators and allocation limits.
* Removed silent plot-call fallback and region midpoint shifts; comparisons
  tolerate different optional schemas, and raw segments do not share a corrected
  call legend.
* Made heatmap ordering deterministic, clustering explicit, missingness safe,
  and continuous color ranges inclusive of observed extremes.
* Added adversarial scientific and rendering regression tests, methods contracts
  and a manuscript release gate. This remains a private development version.

# ichorViz 0.0.0.9000

* Added parsers for `.cna.seg`, `.seg`, and `.params.txt` output.
* Added validated sample and cohort objects.
* Added genome-wide, regional, and multi-sample comparison plots.
* Added cohort bin matrices and optional ComplexHeatmap rendering.
