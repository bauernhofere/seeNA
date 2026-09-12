# ichorViz 0.0.0.9002

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
