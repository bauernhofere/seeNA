# Decision register: rationale, evidence and limits

**Scope:** the current seeNA development branch (see `git log` for the
exact commit).
This register makes decisions defensible and reviewable; it does **not** certify
clinical validity or retroactively turn engineering choices into published
methods. No literature reference below endorses seeNA or all its defaults.

## How to read the evidence

- **Contract:** behavior traceable to a source format, reference assembly or
  upstream implementation. Compatibility is limited to the inspected versions.
- **Policy:** an explicit package design choice. Tests establish that we implement
  it consistently, not that it is uniquely correct for every scientific question.
- **Study choice:** a decision the study investigators must justify and approve.
  The package must not make it silently.

Status "retained" means retained in the development implementation, **not**
author approval for a manuscript. Tests, source review, real-file compatibility,
visual smoke checks and scientific study validation are different evidence types.
The [methods contract](../inst/methods.md) specifies the exact behavior; this
register explains *why* and what remains uncertain.

## Primary sources

**U1 — ichorCNA formats.** [output.R at the inspected commit](https://github.com/GavinHaLab/ichorCNA/blob/d31ed52e9e9aa225084a6d380ba0a81a21e1943a/R/output.R):
`outputHMM()` distinguishes `.seg`, `.seg.txt` and `.cna.seg`;
`outputParametersToFile()` reports TF as `1 - n` and tumor ploidy `phi`, both to
four significant digits. The upstream tag is v0.4.0; its DESCRIPTION says 0.3.4.
Use the commit/configuration rather than assuming the two version labels agree.

**U2 — ichorCNA plotting.** [plotting.R at the same commit](https://github.com/GavinHaLab/ichorCNA/blob/d31ed52e9e9aa225084a6d380ba0a81a21e1943a/R/plotting.R):
`plotGWSolution()` computes mixture ploidy; `plotCNlogRByChr()` adds its log2
relative-to-two shift to bins **and** segment medians. These are display operations.

**U3 — ichorCNA CN correction.** [utils.R at the same commit](https://github.com/GavinHaLab/ichorCNA/blob/d31ed52e9e9aa225084a6d380ba0a81a21e1943a/R/utils.R):
`correctIntegerCN()` keeps raw and corrected fields, with purity/configuration-
dependent and chromosome-X-specific handling. `logRbasedCN()` divides by
purity and cellular prevalence; its diagnostic output differs from plotted logR.

**G — Reference lengths.** UCSC primary chromosome lengths from
[hg19.chrom.sizes](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes)
and [hg38.chrom.sizes](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes).
All 24 supported chromosome lengths per build were compared numerically against
these downloads in this audit. Full-download MD5s were
`b3b0fcf79b5477ab0b3af02e81eac8dc` (hg19) and
`c42e9f75fff906c6a7143cb2fac86602` (hg38). These URLs are not immutable archives;
the implemented integer lengths and commit are the reproducible reference.
`make audit-reference` repeats the comparison (network required; not a package test).
This verifies the lookup table, **not the assembly of a user's samples**.

**P — Original method publication.** Adalsteinsson, Ha, Freeman et al. (2017),
*Scalable whole-exome sequencing of cell-free DNA reveals high concordance with
metastatic tumors*, Nature Communications 8:1324,
[doi:10.1038/s41467-017-00965-y](https://doi.org/10.1038/s41467-017-00965-y).
Cite this for ichorCNA's original method, not as evidence for our rebinning,
terminal-padding heuristic, plotting palettes or manuscript-specific thresholds.

## Input and coordinate decisions

### D01 — R, file-based scope and small public objects — retained policy

**Decision:** stay in R, reuse ggplot2 and optional ComplexHeatmap, return ordinary
validated lists/data frames/matrices, and expose small helpers. No Python bridge,
implicit cache, upstream model fitting or clinical-spreadsheet discovery.

**Why:** this matches the downstream R workflow, avoids cross-language state and
lets users style/export standard objects. A wrapper or shared mutable analysis
state would add installation and provenance boundaries without solving the core
file-to-figure task. This is an engineering trade-off, not evidence R is globally
faster. Optional heatmap dependencies keep the basic reader/profile path smaller;
development-only tooling is not a runtime dependency.

**Evidence/limit:** [DESCRIPTION](../DESCRIPTION), [public API](../NAMESPACE),
local smoke timings in [validation](downstream-validation.md). Dense matrices
still have memory limits; Windows loading is sequential, not parallel.

### D02 — Strict single-sample formats and identity before aliases — retained contract/policy

**Decision:** explicit known aliases; reject ambiguous prefixes/duplicate semantic
columns and malformed values. Support both segment formats in U1. Preserve optional
fields; never select the first sample from a multi-sample export. Reconcile source
IDs before applying user aliases.

**Why:** U1 permits multiple samples and different segment layouts; choosing by
column order or a permissive substring risks plausible wrong-sample plots. Strict
rejection is preferred to guessing. This does not mean upstream multi-sample
exports are invalid—only outside the supported sample API.

**Evidence:** U1; [read.R](../R/read.R), [sample.R](../R/sample.R);
[adversarial regressions](../tests/testthat/test-adversarial.R) for component
ambiguity, identity, both segment formats and conflicting parameters. Appended
candidate-diagnostic tables are not another selected sample.
**Limit:** matching IDs cannot prove all files came from the same fitted run.

### D03 — Missingness and numerical domains — retained contract/policy

**Decision:** preserve missing calls/logR/flags; reject unknown nonmissing call or
flag tokens. TF is a fraction in [0,1], ploidy positive, coordinates finite
integers. Diagnostic inverse-CN infinity is allowed only in that diagnostic
field—not logR, plotted medians or supported matrix measurements.

**Why:** unknown is not neutral and missing is not false. U3's inverse diagnostic
can be singular at zero TF; rejecting that field would reject otherwise usable
output, while plotting it as logR would be wrong. No implicit zero imputation.

**Evidence:** U1/U3; [read.R](../R/read.R), [sample.R](../R/sample.R);
[tests](../tests/testthat/test-adversarial.R), including zero-TF diagnostics.
**Limit:** accepted missing fields can leave a particular plot/analysis unavailable.

### D04 — Explicit assemblies, closed intervals and shared axes — retained contract/policy

**Decision:** hg19/hg38 primary chromosomes 1–22/X/Y, 1-based inclusive integer
intervals; normalize `chr` and numeric sex-chromosome aliases. Validate sorted,
nonoverlapping intervals. No liftover, assembly guessing or silent contig removal.
Use reference-length cumulative axes, not sample-observed extents. Layout/profile
defaults include 1–22/X (profiles add observed Y); matrix defaults include Y too.
Choose/reuse the same chromosome set for custom layouts.

**Why:** inclusive widths match U1 (`end - start + 1`); shared references prevent
sparse inputs from shifting chromosomes between overlays. Whole-reference axes
can look empty for sparse data: that is preferable to implying coverage.

**Evidence:** U1/G; [genome.R](../R/genome.R), [sample.R](../R/sample.R);
[region tests](../tests/testthat/test-region.R) and
[contract regressions](../tests/testthat/test-contract-regressions.R).
**Limit:** unsupported assemblies/contigs need a deliberate future API, not aliases.

### D05 — Default fixed-window terminal padding — retained, qualified policy

**Decision:** `bounds='window'` clips only supported grid padding, records each
original/new interval, and emits a notice. Require at least two equal-width
interior bins, consistent 1-based grid geometry, less-than-one-window terminal
overhangs, and matching observed bin boundaries for overhanging segments. Filtering
may leave gaps. Keep strict `error` and broader explicit `trim` available.

**Why:** the earlier strict default rejected routine fixed-width terminal windows
in the real-input workflow. Arbitrary "less than 1 Mb" trimming would be too
permissive and tied to one bin width. The observed-grid rule is narrower and
works at other regular widths. Two interior bins is a minimal engineering
support rule, **not** a calibrated confidence threshold or upstream requirement.

**Evidence:** [bounds.R](../R/bounds.R), [downstream tests](../tests/testthat/test-downstream-workflow.R),
and 71 canonical imports/213 recorded changes in [validation](downstream-validation.md).
Tests show that two builds can share a rounded terminal window.
**Limit:** this heuristic cannot establish assembly correctness, recover missing
bases, or handle every variable-width input. Strict rejection remains appropriate
when original source geometry or assembly is not trusted. Final reference choice
requires run-configuration evidence.

### D06 — Atomic cohort loading and specific notices — retained policy

**Decision:** any failed sample aborts the cohort; preserve manifest/metadata order.
Aggregate only classed coordinate notices. Unix can use explicit fork workers;
Windows warns and loads sequentially.

**Why:** silent row dropping changes denominators and pair sets. Blanket warning
suppression hides unrelated issues. Ordering and failures must be independent of
worker completion order. A sequential Windows fallback avoids another worker
framework dependency, but is a disclosed performance limitation.

**Evidence:** [cohort.R](../R/cohort.R);
[atomic-failure tests](../tests/testthat/test-adversarial.R) and
[unrelated-warning test](../tests/testthat/test-downstream-workflow.R).
**Limit:** importing successfully does not adjudicate sample eligibility.

## Measurement and aggregation decisions

### D07 — Raw versus corrected calls and display categories — retained contract/policy

**Decision:** require explicit matrix measurement; never substitute `event` for
`corrected_call`. Map HOMD/HETD/NEUT to Deep loss/Loss/Neutral; collapse
GAIN/AMP/HLAMP variants to Gain for the four-state display.

**Why:** U3 distinguishes raw and corrected fields. The four-state vocabulary
is a visualization simplification; it is **our policy**, not an upstream claim
that all amplification levels are equivalent. Source fields remain available.
Collapsed categories are formed **before** categorical rebinning; their combined
support can therefore win over another category.

**Evidence:** U1/U3; [utils.R](../R/utils.R), [matrix.R](../R/matrix.R),
[call tests](../tests/testthat/test-cohort.R).
**Limit:** to study amplification magnitude, use an appropriate CN measurement
or the original fields, not four-state codes.

### D08 — Optional upstream logR shift, raw default — retained contract/policy

For fitted TF `f` and tumor ploidy `phi`:

```
mixture_ploidy = f * phi + (1 - f) * 2
adjusted_logR  = raw_logR + log2(mixture_ploidy / 2)
```

**Decision:** apply exactly this positive additive display shift to bins and
segment medians when requested; never mutate source values or reassign calls.
Keep `ploidy_adjust=FALSE` by default, label the scale, store shift metadata,
and error on missing parameters when TRUE.

**Why:** U2 supplies the formula. Raw default is a deliberate deviation from the
upstream plotting default: it preserves exported values and supports input without
parameter files. A silent default change would alter existing plots. Requiring
opt-in is not a claim the raw scale is biologically superior.

**Evidence:** U1/U2; [accessors.R](../R/accessors.R), [plot.R](../R/plot.R),
[vector/plot regressions](../tests/testthat/test-downstream-workflow.R),
real-input arithmetic reconciliation in [validation](downstream-validation.md).
**Limit:** U1 rounds fitted parameters, so bitwise agreement with full-precision
upstream plots is not promised. This is not purity deconvolution; gains/losses,
particularly on X, are not guaranteed to have a particular y sign.

### D09 — Neutral-CN evidence, not baseline inference — retained cautious policy

**Decision:** report CN among NEUT bins by autosome/X/Y, pairing raw or corrected
call and CN fields. One distinct observed value is reported; none or several
yield NA with status/counts/candidates. No modal tie-break, gender-based guess,
automatic CN 2, or borrowing from autosomes for X.

**Why:** U3 makes X handling configuration/purity dependent. A hardcoded reference
can silently change downstream CN-minus-reference measurements. Reporting evidence
is more honest than claiming to recover an unobserved run-level baseline.

**Evidence:** U3; [accessors.R](../R/accessors.R),
[missing/ambiguous/X tests](../tests/testthat/test-downstream-workflow.R).
**Limit:** `supported` means *one distinct observed CN*, not statistical confidence.
One observed bin can suffice even if others have missing CN; inspect counts.
Autosomal pooling can be ambiguous. Even unanimous observed values do not prove
biological neutrality or all run-level thresholds.

### D10 — Base-pair means for continuous rebinning — retained mathematical policy

**Decision:** weight each finite source value by inclusive interval overlap:
`w = max(0, min(end_i,end_j) - max(start_i,start_j) + 1)`;
return `sum(w * value) / sum(w)` over observed support.

**Why:** unweighted bin means distort unequal-width/split boundary contributions.
This summarizes the selected numerical field over genomic span; it is not a
refit of ichorCNA. Mean logR is a mean **on the log scale**, not the log of mean
linear ratio. Averaged CN can be fractional and is not a new integer CN call.

**Evidence:** [matrix.R](../R/matrix.R), hand-calculated split/weight tests in
[adversarial regressions](../tests/testthat/test-adversarial.R).
**Limit:** assumes a source-bin value represents its interval uniformly. This is
not read-depth-, uncertainty- or purity-weighted inference; narrow events can be
diluted. Bin width must suit the biological question.

### D11 — Base-pair mode for categorical rebinning — retained display policy

**Decision:** greatest observed BP support wins; ties yield NA, heterogeneous
support sets `mixed=TRUE`. Never average category codes to produce a state.

**Why:** equal gain/loss support must not cancel into a fabricated neutral call.
Alternatives include "any alteration", "most severe state" or a separate mixed
category; those answer different questions and may exaggerate small events.
Mode is a conservative summary of dominant observed span, **not** a unique
biological definition or an upstream calling algorithm.

**Evidence:** [matrix.R](../R/matrix.R), [mode/tie tests](../tests/testthat/test-adversarial.R).
**Limit:** minority events can disappear from the displayed value; inspect `mixed`
and use smaller bins/region plots. Invariant under splitting an interval into
same-state subintervals, but not under changing target width or category mapping.

### D12 — Coverage and no imputation — retained policy

**Decision:** `coverage = observed selected-value overlap / actual target width`.
Default `min_coverage=1`; zero support always yields NA. Shorten terminal target
bins to reference length. Keep values/coverage/mixed aligned. The implementation
allows a `1e-12` floating-point comparison tolerance, not a biological threshold.

**Why:** filling filtered windows with neutral states invents observations. Full
coverage prevents partial intervals from silently representing entire target
bins. A lower threshold can be useful, but must be explicit and recorded.

**Evidence:** [matrix.R](../R/matrix.R), [coverage tests](../tests/testthat/test-adversarial.R).
**Limit:** coverage concerns the chosen field, not sequencing depth, certainty or
fraction of tumor cells. Full coverage can discard useful partial evidence;
filtering may be systematic, not missing at random. NA can represent no coverage,
insufficient coverage or a categorical tie; inspect diagnostic layers. The ~17.9%
missingness in one local matrix is not a universal expected rate.

### D13 — Allocation guards and 1 Mb default grid — retained engineering policy

**Decision:** a configurable 1 Mb target grid; default 50 million cell guard;
values/coverage/mixed need at least 20 bytes per cell before other overhead.

**Why:** an unbounded dense matrix can terminate a session. A familiar coarse grid
is convenient, not a claim of validated biological resolution. Sparse/chunked
storage is a future option if actual workloads justify its complexity.

**Evidence:** [matrix.R](../R/matrix.R), [allocation tests](../tests/testthat/test-adversarial.R).
**Limit:** the guard is not a peak-RAM guarantee. Choose bin width/chromosomes and
resources explicitly; neither the cell cap nor 1 Mb is a clinical cutoff.

## Display, ordering and provenance decisions

### D14 — Original midpoints and distinct segment semantics — retained policy

**Decision:** plot original source-bin midpoints; clip only segment endpoints in
regional views (after any audited import trimming). Do not move points to the
middle of their visible overlap. Grey profile segments avoid merging raw segment
events into the corrected-bin legend. Comparison colors identify samples.

**Why:** moving/rebinning points just to align overlays changes the displayed
observations. Segment states and bin calls can differ (U3). Colored raw segments
could be useful, but require an explicit separate interpretation/legend.

**Evidence:** U3; [plot.R](../R/plot.R),
[midpoint/schema tests](../tests/testthat/test-adversarial.R).
**Limit:** a bin overlapping a region can have its original midpoint outside the
view. Absence of a visible point does not establish absence of overlapping input.
Raw-segment state coloring remains deferred rather than silently conflated.

### D15 — Palettes, scales and rasterization — retained display policies

**Decision:** named palettes, missing-state grey, discrete call colors; continuous
heatmap defaults cover observed extremes. CN 2 is a color reference only; callers
can set explicit saturation. Annotation colors are deterministic for the supplied
metadata levels. Large heatmaps use ragg if available, without numerical matrix
resizing/averaging of categorical codes.

**Why:** clipping outliers or guessing a neutral baseline in the default color
scale hides assumptions. Rasterization is a rendering optimization, not a new
matrix transformation. Independent panels need fixed palettes/scales for comparison.

**Evidence:** [palettes.R](../R/palettes.R), [heatmap.R](../R/heatmap.R),
[headless rendering tests](../tests/testthat/test-contract-regressions.R).
**Limit:** image resolution can conceal narrow events even without matrix resizing.
Color interpolation, hue choices and perceptual accessibility are presentation
choices, not biological evidence. Defaults are not certified color-vision-safe.
Neutral profile points are blue while neutral heatmap cells are near-white;
legends are authoritative. Specify fixed annotation palettes across subsets.

### D16 — Manifest order; clustering only by explicit request — retained, qualified policy

**Decision:** no default clustering or genomic column reorder. Optional row
clustering uses Euclidean distance/complete linkage on bins observed in **every**
sample, with no imputation; requires at least two such bins and 2–2000 samples.

**Why:** stable manifest order preserves known pairing/timepoint structure.
Common observed bins provide a common distance feature set; pairwise missingness
handling would compare different regions for different pairs. Two bins is a
minimal implementation guard, not sufficient evidence for meaningful clustering;
2000 rows limits quadratic work, not a statistical boundary.

**Evidence:** [heatmap.R](../R/heatmap.R), [ordering/clustering tests](../tests/testthat/test-adversarial.R).
**Important limit:** for calls, Euclidean distance treats codes -2/-1/0/1 as
numerically spaced. That is an **additional ordinal-distance assumption**, not
implied by using categorical modes. It is exploratory and not a validated
biological similarity metric. Prefer a justified continuous measurement for
Euclidean clustering, or supply an externally justified row order. A dedicated
categorical distance would need a separate explicit API/validation. Missingness
can also remove biologically important regions from the common-bin set.

### D17 — Immutable import fingerprints and privacy boundaries — retained policy

**Decision:** record MD5 and bytes before/after import; abort if they disagree;
retrieve stored fingerprints rather than rehashing current files. Keep build,
version, schema and transformation metadata. Paths require explicit opt-in;
metadata allowlists are available but non-file annotations are retained by default.

**Why:** late hashing can describe files different from those that produced a plot.
Input-time checks identify common accidental drift. Minimal path retention avoids
unnecessary path disclosure without disabling intentional annotations.

**Evidence:** [provenance.R](../R/provenance.R), [sample.R](../R/sample.R),
[drift/privacy tests](../tests/testthat/test-adversarial.R).
**Limit:** MD5 is not a security signature; before/after checks are not filesystem
locking or protection against adversarial edits. Aliases, hashes, coordinates,
metadata and figures can remain identifying. No automatic de-identification or
permission to publish clinical data is implied. Save the commit separately;
package version alone cannot distinguish all development checkouts.

### D18 — Evidence gates, tooling and public example images — retained policies

**Decision:** unit/adversarial tests, actual device rendering, package checks,
pinned documentation generation and cross-platform CI are required. Keep public
example images reproducible from the existing tiny format fixtures only; exclude
clinical files and figures from version control. Public repository visibility
has author approval following the [publication preflight](publication-audit.md);
a tagged/archived release requires separate approval. GPL-3-or-later and upstream
attribution are recorded; final study-code contributor/copyright review is still
required.

**Why:** object-class tests miss rendering failures; the previous PDF test checked
bytes before device closure, and a roxygen version change caused documentation
drift. Fixtures prove API behavior, not ichorCNA accuracy. Clinical realism cannot
justify republishing patient profiles without permission.

**Evidence:** [CI](../.github/workflows/R-CMD-check.yaml),
[validation history](downstream-validation.md), [NOTICE](../NOTICE.md),
[fixture provenance](../inst/extdata/README.md), [README source](../README.Rmd).
An earlier development commit passed all five CI configurations:
[PR run](https://github.com/bauernhofere/seeNA/actions/runs/34724852864) and
[push run](https://github.com/bauernhofere/seeNA/actions/runs/34724851027).
**Limit:** minimum declared R 4.1 has not been separately verified by these CI
configurations. Rendering/tests/source compatibility are not independent clinical
validation or author endorsement. No release DOI/date is invented.

### D19 — Paired directional call agreement on a shared grid — retained policy

**Decision:** `ichor_pair_concordance(a,b)` takes an explicit, same-build pair with
unique IDs. Reuse fixed-grid BP aggregation with independent call/logR coverage.
Also require overlapping observed call support (`joint_coverage`) to meet the
same threshold; disjoint partially observed source spans are not co-observations.
Partial logR means can still represent different observed portions, so they are
descriptive heights, not a paired difference.
Unknown, insufficiently covered, tied **or heterogeneous** call support makes the
comparison unknown, even if a categorical mode has a unique winner. Retain every
target bin and its coverage, mixed flags, reason, source-scale and selected-scale
logR, settings and provenance. No inner join that silently drops unmatched bins.

Classify both-neutral, a-only, b-only, concordant (both altered in the same
direction), discordant (gain versus loss), or unknown. Deep loss and loss agree
in direction, not magnitude. Classification uses the selected source-call layer,
never the sign of logR; `event` is an explicit alternative to `corrected_call`.

**Why:** the manuscript's exact-coordinate inner join is adequate only for matching
source grids. It discards unmatched/missing calls without an availability track.
Coverage-aware alignment makes that loss of evidence visible. Independent dominant
modes could agree while hiding opposing minority events, hence mixed bins are
unknown in this comparison. Same-direction agreement is not statistical agreement,
truth, clinical sensitivity or proof of fluid-specific absence in the other sample.

**Sex-chromosome policy:** no NEUT X/Y bins does **not** prove a gain call is false:
a whole chromosome can genuinely be altered. Default `sex_chromosomes='flag'`
retains comparable source calls and reports unresolved reference evidence.
The explicitly conservative `require_neutral` option instead makes X/Y comparisons
unknown when either sample lacks one distinct observed NEUT-bin CN. This records
an evidence requirement, not biological reclassification. No source calls are
erased and no autosomal/diploid baseline is substituted. Even a supported NEUT
value is only observed evidence (D09), not verification of calling calibration.

**Evidence:** [concordance.R](../R/concordance.R),
[classification/alignment/sex-policy tests](../tests/testthat/test-concordance.R).
**Limit:** common bins do not eliminate differing TF, calling thresholds, noise or
sex-chromosome conventions. Aggregation assumes uniform source-bin values; select
a grid appropriate to the question. This is a new display measurement, not exact
reproduction of every manuscript join or a validated diagnostic comparison.

### D20 — Separate agreement colors from bar heights — retained display policy

**Decision:** `plot_ichor_concordance()` stacks original profiles above a target-grid
track on shared fixed y limits. Default `height='both'` retains both altered
heights in side-by-side half-bin bars; one-sided categories show the altered sample.
Half-bin positioning is visual packing, not a sub-bin measurement or breakpoint.

The opt-in `representative` rule reproduces the manuscript's *height policy*:

- both neutral: zero (hidden);
- one altered sample: that sample's logR;
- same-direction alterations: mean of the two logR values;
- opposite-direction alterations: value with larger **absolute** logR, preserving
  its sign; exact absolute ties choose a.

This is not the larger signed value, not a subtraction, and not a measure of
"amount of discordance". It can suppress one opposite-direction signal; swapping
samples reverses the selected sign in equal-absolute discordant ties. Both heights
must be available for a two-sided representative; NA is never silently ignored.
The default retains both to avoid this information loss. Averaging continuous
logR is permitted here, unlike averaging categorical codes; concordant source
calls do not guarantee matching logR signs or a nonzero mean.

Grey background marks unavailable comparisons/heights, gold strips flag unresolved
sex-reference evidence, and both-neutral bins stay blank. Exact-zero altered
heights have a marker without y nudging. Sample colors and agreement colors have
separate legends. Fitted TF/ploidy remain visible; identifiers can be hidden, but
plot data/attributes are not anonymized. Raw logR remains default; optional
ploidy adjustment applies the same run-specific shift to bins, segments and track
heights. `ylim=c(-2,2)` clips the view only and warns for out-of-range heights.
Relative `panel_heights`, `segment_linewidth` and `point_stroke` are presentation
controls recorded in `ichor_view`; they do not change measurements or limits.
Defaults retain equal panels, segment width 0.55 and inherited point stroke.
The package supports ggplot2 >= 3.5.0. Unequal proportions use native ggplot2
>= 4.0.0 panel sizing; on older versions the request errors rather than silently
losing the requested layout. Equal panels and all other styling remain supported.
This avoids a package-wide 4.0 requirement and preserves a themeable ggplot without
a composition dependency. A custom-style snapshot and numeric layout/data
invariance tests cover these controls. One-sided default track colors follow
`sample_colors`; an explicit agreement palette takes precedence.

The existing linear joint-coverage sweep is retained. Fine-bin profiling should
precede an optimization; any replacement must retain interval/NA semantics and
avoid unexpectedly large overlap joins. No fine-grid speed improvement is claimed.

**Evidence:** [plot-concordance.R](../R/plot-concordance.R), numeric/device tests in
[test-concordance.R](../tests/testthat/test-concordance.R), and new vdiffr cases in
[test-concordance-snapshots.R](../tests/testthat/test-concordance-snapshots.R).
**Limit:** fixed limits/finite pixels can hide magnitude or narrow events. Region
cropping follows aggregation, so a displayed partial target bin still summarizes
its complete target interval. Upper points retain source midpoints, lower bars
use the target grid. No automatic equivalence with the manuscript is claimed.

### D21 — Grouped heatmap rows and presentation arguments — retained display policy

**Decision:** ComplexHeatmap stays the only heatmap renderer. `plot_ichor_heatmap()`
maps `group` to `row_split` with slice order fixed (factor levels, otherwise first
appearance in `row_order`), no slice clustering, and rejects `group` with
`cluster_rows`. `row_labels` replace displayed text while rows stay keyed by sample
ID. `legend_direction` and `text_style` (`fontfamily`, `fontsize`, `col` only) are
presentation controls; call legend keys are outlined so a white state stays visible.

**Why:** paired and longitudinal designs need blocks of rows in a caller-supplied
order, and figures need their text to match surrounding panels. A separate ggplot
landscape renderer was prototyped and rejected: it duplicated this function's
semantics for the sake of theming. Callers pass their own style values instead of
the package translating ggplot themes.

**Evidence:** [heatmap.R](../R/heatmap.R), block-order, label, style, legend and
validation tests in [test-heatmap-grouping.R](../tests/testthat/test-heatmap-grouping.R).
**Limit:** block membership and order are study choices, not results. `text_style`
does not cover gaps, borders or colours. Fonts must exist on the device that
builds and draws the heatmap, because label widths are measured when the object is
created. Placing the drawn heatmap in a ggplot layout is left to the caller.

### D22 — Exported segment-median matrices — retained contract/policy

**Decision:** `ichor_matrix(value = "segment_median")` aggregates the exported
segment `median` directly onto the requested grid (1 Mb by default), using D10's
bp-weighted continuous mean and coverage threshold. The source bin logR remains
unchanged. Every sample needs an explicitly supplied selected-run segment table;
there is no fallback to bins, refitting, purity correction or ploidy shift.
Missing medians and segment gaps contribute no support. Coverage refers to finite
segment spans, which may bridge absent/filtered source bins; it is **not**
observed bin/read coverage and no bin mask is applied. At boundaries this is a
weighted mean of exported medians, not a newly computed median. Continuous
`mixed` flags remain FALSE. `value`, `aggregation` and source-file fingerprints
record the operation and inputs without a new renderer or data overwrite.

**Why:** continuous segment-level signal retains amplitude while reducing
bin-level speckling. Discrete calls discard amplitude; tumor-CN estimates are a
different quantity. Overwriting bin logR with segment values would falsely label
the source layer and make coverage semantics opaque. Segment-span coverage is
explicit rather than silently borrowing a mask from a different measurement.

**Evidence:** U1/U2 export and plot segment medians; [matrix.R](../R/matrix.R)
and [test-segment-matrix.R](../tests/testthat/test-segment-matrix.R) cover exported
values, weighted boundaries, NA/coverage, absent segment files/chromosomes,
short terminal bins, independent segment/bin support, input immutability,
validation/allocation guards and actual ComplexHeatmap rendering.
**Limit:** segmentation is fitted smoothing, not independent evidence for an
alteration or validation of a low-TF fit. LogR amplitude depends on TF as well as
CNA magnitude; weaker urine signal does not establish smaller tumor alterations
or their biological absence. Shared display limits, fitted TF labels and final
manuscript interpretation remain study choices. This option is unshifted logR,
not the optional upstream plotting shift in D08.

## Study decisions that this audit cannot approve

| Decision | Current boundary | Required evidence / owner |
|---|---|---|
| Source assembly and chosen fit | Explicit caller inputs; no automatic ranking | Investigators freeze run configuration and component manifest |
| External reported TF to resolve candidates | Documented recipe, predeclared tolerance, fail on ambiguity | Document external provenance/reporting precision and any adjudication; do not widen tolerance to obtain a preferred result |
| Pairing, timepoints, eligibility and exclusions | Supplied by study manifest | Reconcile expected/loaded/plotted sets and explain every exclusion |
| Measurement for manuscript heatmap | Earlier snapshot: CN minus fitted ploidy. Live port: CN minus a NEUT-bin reference with wrapper-specific fallbacks | Investigators approve the change and reference evidence; neither is interchangeable with calls or absolute CN |
| Dead band 0.35, clipping at ±2, palette saturation at ±1 | Existing study display choices, **not** package defaults or CNA thresholds | State rationale and sensitivity; no biological optimum has been demonstrated here |
| Plasma-burden row ordering | Existing study presentation | Define burden/denominator/missingness and freeze tie/pair ordering policy |
| Sex-chromosome reference and interpretation | Evidence helper does not guess | Inspect run configuration and available observations; report unsupported/ambiguous cases |
| Final figures, attribution and citation | Draft manuscript figures; public development code, no archive or release tag | Author visual/scientific approval, contributor review, pinned artifact and approved archive/publication |

See [review triage](review-triage.md) for the existing study workflow and
[vignette](../vignettes/manuscript-workflow.Rmd) for the explicit selection recipe.
Thresholds, fitting, pairing and manuscript statistical models are not validated
by this package's tests or by a citation to the original ichorCNA paper.

## Live manuscript port audit

A fresh source inspection found that the manuscript wrapper has evolved since
the earlier hardening review. This is a **code audit, not a rerun or numerical
validation**. No patient-level content is reproduced here.

1. **Measurement changed.** The live `R/ichor_io.R` calculates corrected CN minus
   a reference from NEUT-bin CN; `R/cn_landscape.R` consumes that quantity. Earlier
   review notes described CN minus fitted ploidy. Some current comments still
   mention the old measurement. Freeze the actual code/input version, align
   captions/comments and explicitly approve the change before comparing figures.
2. **Wrapper baseline is less conservative than the package helper.** The local
   function takes a mode, falls back to autosomal CN 2 if absent, and borrows that
   value for X if X evidence is absent. Those are assumptions, not observed
   evidence. Its name also masks the package export; merely updating seeNA
   does not adopt the package's diagnostic behavior. Migrate deliberately using
   `seeNA::ichor_neutral_cn()` and an explicit missing/ambiguous-evidence policy;
   its data-frame result is not a drop-in replacement for the wrapper's vector.
3. **Run ambiguity remains.** The resolver chooses the first candidate within
   `tol=0.01` (one percentage point of TF), not a unique eligible solution. This
   threshold is not automatically justified by U1's four-significant-digit
   output. A directory-priority rule may be legitimate external adjudication,
   but requires evidence and approval; it is not proof of fit equivalence.
4. **Silent exclusions/warnings remain outside the package.** The wrapper uses
   broad warning suppression, and the landscape loader catches failed pair loads
   and returns NULL. Package-level atomic imports do not protect a downstream
   caller that explicitly catches errors. Reconcile expected versus plotted pairs
   and record every failed/excluded sample privately; remove blanket suppression
   when migrating to classed package coordinate notices.
5. **The shift does not prove every call's y sign.** Wrapper prose asserting that
   all corrected gains/losses must fall above/below zero overstates U2. Update that
   interpretation when reconciling the manuscript methods; see D08/D09.

These are unresolved investigator-owned choices, not package behavior endorsed
by this register. This documentation pass does not change the manuscript's
calculations, source data or sample selection.

## What would change a retained decision?

A reproducible counterexample, a newly supported upstream format, measured
workload constraints or a justified study requirement can motivate a change.
Record the affected decision ID, alternative and trade-off; add a regression and
update the methods, user-facing defaults and validation evidence. Prefer an
explicit new measurement/option over changing the meaning of an existing one.
The current register is a reviewable rationale, not a claim that every choice is
optimal, universal, or final.
