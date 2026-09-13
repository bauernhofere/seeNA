# Downstream workflow follow-up — development version 0.0.0.9002

This follows `hardening-validation.md`. It is not a release or approval of final
manuscript figures. The repository remains private and the hardening PR is draft.
A subsequent [decision audit](decision-register.md#live-manuscript-port-audit)
records changes in the live manuscript wrapper since the earlier review; the
historical CN-minus-ploidy description below is not its current measurement.

## Decisions from the downstream review

- Added default `bounds='window'` using observed fixed-bin grid evidence, not a
  generic sub-megabase tolerance. All changes remain recorded. Strict `error`
  and explicitly broader `trim` remain available. Some hg19/hg38 chromosome
  lengths round to the same terminal window: passing this check does not verify
  assembly identity. Regression tests demonstrate that limitation.
- Coordinate notices are aggregated once per cohort. Only the package's specific
  padding conditions are intercepted; an unrelated-warning regression verifies
  other warnings remain visible.
- Exported `ichor_tf`, `ichor_ploidy`, `ichor_call_state`, `ichor_genome_layout`,
  `ichor_adjusted_logr`, and `ichor_neutral_cn` for custom manuscript wrappers.
- Verified the positive additive logR shift against upstream v0.4.0
  `plotGWSolution` and `plotCNlogRByChr` in `R/plotting.R`. It uses mixture ploidy
  and applies identically to bins and segment medians. Exported rounded fit
  parameters can differ from upstream's full-precision in-memory values.
- Added `ploidy_adjust=TRUE` to all three profile/overlay/region APIs. Kept raw
  default FALSE deliberately: samples without parameters still plot, and a new
  version does not silently alter y values. Labels and transform metadata make
  the selected scale explicit; adjustment is not a new CN call.
- Neutral-CN helper reports observed evidence, not a presumed baseline. Missing
  or conflicting NEUT-bin CN returns NA with status/counts/candidates. It never
  borrows an autosomal reference for X or assumes male X must be CN 1.
- Added an installed, executable workflow vignette. Run-selection support is a
  documented, auditable recipe using externally specified TF/precision, not an
  automatic scientific optimum or first matching directory. Ambiguity aborts.
- Documented expected NA coverage due to upstream filtering. The missing fraction
  is dataset-specific, not a package-wide expected percentage.
- Moved devtools/roxygen2 out of Suggests into developer tooling configuration.

Optional colored raw segment states remain deferred: grey currently avoids
merging raw segment events into the corrected bin-call legend. A future option
needs a clearly separate legend, not silent reuse of the bin interpretation.
The manuscript's CN-minus-ploidy landscape still needs reconciliation with its
explicit dead band, saturation and ordering; call heatmaps are not a substitute.

## CI failures diagnosed

The prior remote run was not green despite local success:

1. Four jobs failed a heatmap PDF byte-size assertion while its device was still
   open. Tests now close the device in `finally` before checking output size.
2. Linux release failed generated-documentation verification: CI used roxygen2
   8.1.0 while local/generated docs used 7.3.3. The generator is now pinned and
   checked against RoxygenNote; documentation drift checks remain enabled.

Subsequently confirmed: both [PR](https://github.com/bauernhofere/ichorViz/actions/runs/34724852864)
and [push](https://github.com/bauernhofere/ichorViz/actions/runs/34724851027) runs
for `1dd2d5d` succeeded on all five configurations. This does not approve the
manuscript workflow.

## Local evidence

- 44 test cases / 203 passing assertions; no failures or errors.
- R CMD check (including installed vignette execution and rebuilding): zero
  errors/warnings; environment NOTE `unable to verify current time`.
- All 71 canonical samples imported under the new default, with one summary for
  213 recorded terminal changes. No blanket warning suppression was needed.
- A 71 x 3,102 full-coverage call matrix retained 39,382 missing cells (~17.9%).
  Its values remained documented categorical codes or NA.
- Warm smoke timings: ~2.17 seconds for import and ~0.67 seconds for matrix
  construction; matrix object ~5.8 MB including diagnostic layers/provenance.
- Adjusted bin vectors and segment medians for every canonical sample agreed
  with independently written arithmetic using the exported TF/ploidy.
- Raw/adjusted real-input profile, paired overlay and regional figures rendered;
  adjusted profile and region were visually inspected. The cohort heatmap also
  rendered. This is a smoke check, not author approval or final run selection.
- Clinical files, images, manifests and detailed validation records remain
  outside the package repository.

No release tag, public visibility change, DOI, or final citation is authorized
by these checks. See the release gates in `hardening-validation.md`.
