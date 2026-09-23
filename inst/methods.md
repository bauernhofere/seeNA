# Scientific contracts (pre-release schema 1)

This file defines behavior, not clinical validity. The repository's
[decision register](../docs/decision-register.md) documents primary sources,
alternatives, policy limitations and unresolved manuscript choices. For an
installed copy of these methods, consult that register in the matching source
checkout; `docs/` is repository-only.

## Supported formats and boundaries

Reference: ichorCNA v0.4.0 `R/output.R`, commit
`d31ed52e9e9aa225084a6d380ba0a81a21e1943a`. Its package DESCRIPTION says 0.3.4;
record upstream commit/run configuration rather than relying on that number.

- `.cna.seg`: coordinates plus plain or single-sample-prefixed measurement columns.
  Multi-sample exports are rejected, not guessed or silently reduced to one sample.
- `.seg`: sample, chr, start, end, event, copy.number, bins, median.
- `.seg.txt`: the distinct ID/chrom/start/end/num.mark/seg.median.logR export.
- `.params.txt`: one selected sample in a tabular header and/or key-value section;
  repeated fields must agree. The appended init/n_est/phi_est candidate table is
  diagnostic, not a second selected solution.

Public parser aliases are explicit in `R/read.R`. Unknown extra columns are not
used; ambiguous semantic matches and malformed numeric values raise classed
errors. Valid missing values stay NA. The diagnostic `logR_Copy_Number` inverse
transform may be infinite at zero TF; it is retained but never used as a plotted
log ratio or as a supported matrix value. This exception does not apply to logR,
segment medians or corrected copy numbers.

The test fixtures exercise schemas, not ichorCNA inference. Reading many files
successfully demonstrates compatibility only with those files, not every version.

## Coordinates and identity

1-based, inclusive integer intervals. Builds hg19/hg38 use primary chromosome
lengths (GRCh37/GRCh38; UCSC chr sizes), chromosomes 1-22, X, Y; chr prefixes and
23/24 aliases normalize to X/Y. No liftover, build guessing, or contig dropping.
Bins/segments must be sorted, nonoverlapping and assembly-bounded. NA logR is
allowed in sample objects; a plotting request without finite logR fails.

Fixed-width terminal bins can exceed reference lengths. Default `bounds='window'`
requires at least two equal-width observed interior bins, consistent 1-based grid
starts/endpoints, and overhangs ending at the rounded terminal window. A segment
overhang must end at an observed overhanging bin boundary. Matched padding is
clipped, messaged, and stored in `coordinate_changes`; unsupported cases fail.
`bounds='error'` rejects all overruns. Explicit `bounds='trim'` remains available
for broader, manually inspected straddling intervals. Its warnings and window
messages are classed and aggregated once per cohort, leaving unrelated warnings
untouched. Wholly out-of-bounds intervals still fail. Several assemblies can
share a rounded terminal window; this is not a genome-build verification test.

A sample object reconciles source IDs across all supplied files before applying
an alias. Aliases replace source IDs in returned tables. With unprefixed bin
columns, source identity falls back to the .cna.seg basename. Renamed standalone
files therefore need matching source naming. Aliases are not multi-sample selectors.
Identity agreement cannot distinguish multiple runs for the same sample: the
caller must explicitly choose one run for all components. The workflow vignette
shows an auditable comparison against externally adjudicated TF, with explicit
tolerance and failure on non-unique matches. The package does not select runs by
file order or optimize fits to a desired manuscript result.

## Measurements and categorical calls

- `logR`: unmodified source signal. Optional upstream-formula plotting adjustment
  is a separate display transform; it does not mutate this column.
- `copy_number`, `event`: raw ichorCNA output fields.
- `corrected_copy_number`, `corrected_call`: caller-corrected output fields, not
  recomputed by seeNA. Corrected calls, absolute CN, logR sign and neutrality
  relative to fitted ploidy must not be conflated.
- `logR_copy_number`: diagnostic only.
- TF: fraction [0,1]; ploidy positive. Missing values remain missing.
- Logical subclone_status recognizes explicit true/false tokens; unknowns error.

Call vocabulary: HOMD -> Deep loss (-2), HETD -> Loss (-1), NEUT -> Neutral (0),
GAIN/AMP/HLAMP and numbered HLAMP variants -> Gain (1). NA/empty -> missing;
other tokens error. There is no fallback from corrected_call to event.

Public helpers: `ichor_tf()` and `ichor_ploidy()` return scalars (NA if absent);
`ichor_call_state()` implements the mapping above; `ichor_genome_layout()` returns
selected reference chromosomes with lengths/offsets/midpoints. Reuse the same
chromosome selection for all samples to preserve coordinate alignment.

`ichor_neutral_cn()` reports observed NEUT-bin CN evidence by autosome/X/Y, paired
with the corresponding raw or corrected CN field. Exactly one distinct observed
CN yields a supported value. Missing or conflicting evidence yields NA plus
status/counts/candidates, never a diploid default or cross-chromosome fallback.
This evidence does not identify every run-level threshold or biological baseline.

## Fixed-bin matrices

`value` is required. `call_column` explicitly chooses corrected_call (default)
or event when value='call'. Metadata/sample order is validated at every boundary.
The continuous bin measurements are `logR`, `copy_number` and
`corrected_copy_number`. `value='segment_median'` instead uses the normalized
segment table's exported `median` (raw log2 ratio). A selected-run segment table
is required for every sample; missing files error, with no bin-logR fallback.
No source values are overwritten and no purity/ploidy adjustment is applied.

For source interval i and target bin j, overlap in bp is:

`w_ij = max(0, min(end_i,end_j) - max(start_i,start_j) + 1)`.

Source intervals do not overlap; only non-missing selected values contribute.
Continuous result: `sum(w_ij * value_i) / sum(w_ij)`.
Categorical result: category with greatest total bp support; equal support for
multiple winners -> NA. Mixed categories do not cancel into neutral. `mixed`
records any heterogeneous observed call support, including unequal modes. The
mode is a display aggregation, not a new biological CNA call.

`coverage = observed overlap / target width`. The final target bin is shortened
to the actual chromosome end. Default min_coverage=1 requires full observed
coverage; zero coverage always gives NA even if min_coverage=0. Lowering the
threshold is explicit and recorded. The numeric values, coverage and logical
mixed matrices share dimensions/dimnames. No implicit NA imputation. Dense matrix
allocation is guarded by max_cells (50 million default); memory use also includes
inputs, temporary calculations, bin labels and metadata. Source bin filtering
can produce substantial NA fractions even for valid input (the exact proportion
is dataset-specific); inspect coverage rather than assuming the import failed.

For `segment_median`, the source intervals in these formulas are **segments**.
A target bin crossing a segment boundary receives a bp-weighted mean of the
exported medians, not a newly calculated median. NA medians, gaps between
segments and absent chromosomes contribute no support; they are not zero.
Coverage measures finite **segment-span coverage**, not bin/read coverage:
a fitted segment can bridge filtered/absent bins or bins with missing logR.
No intersection with the original bin grid is silently imposed. `mixed` is
FALSE for this continuous measurement, even at segment boundaries. The matrix
`value` and `aggregation` identify this operation, and import provenance retains
the selected segment-file fingerprints. Smoothed fitted summaries do not rescue
unreliable low-TF fits. LogR amplitude reflects TF as well as CNA magnitude, so
weaker signal cannot establish a smaller tumor alteration or biological absence.

## Plotting

Profiles color bins by the specified call field; segment medians are grey to
avoid merging raw segment events into a corrected-call legend. Missing calls use
grey, never neutral blue. Comparison overlays use original genomic midpoints,
not aligned/rebinned observations; only segment endpoints are clipped in region
views. Common genome axes use reference lengths even for sparse inputs. No
vertical nudging, hidden downsampling or biological reclassification.

`ploidy_adjust=TRUE` applies `log2((TF*ploidy + (1-TF)*2)/2)` to both bins and
segment medians in profile/comparison/region plots. `ichor_adjusted_logr()` exposes
the same vectors. This reproduces the v0.4 plotting formula using exported,
possibly rounded parameters. Raw output remains the default; missing parameters
when adjustment is requested are an error. Plots label the scale and store shift
metadata in `ichor_transform`. This is not purity deconvolution, integer CN, or a
promise that every gain/loss will lie above/below zero. Chromosome-X conventions
and noise still matter; neither y=0 nor gender alone identifies neutral CN.

Heatmaps preserve manifest order by default. `group` splits rows into blocks in factor-level
or first-appearance order and never reorders rows within a block; it cannot be
combined with clustering. `row_labels` and `text_style` change displayed text only. Every plot shows user-chosen sample
identifiers by default; `show_sample_id = FALSE` (profiles, comparisons, regions)
or `show_row_names = FALSE` (heatmaps) omits them. Categorical colors are discrete. Continuous default scales include
observed extremes (CN 2 is a reference color, not a neutrality decision). Custom
colors may intentionally saturate: authors must report their scale. Unknown,
tied and insufficient-coverage cells use the NA color. Inspect coverage/mixed
layers to distinguish these causes. Rasterized figures have finite display
resolution; use a regional view or more pixels for bin-level inspection.

Optional clustering: 2-2000 samples, at least two bins observed in every sample,
Euclidean distance and complete linkage on these common bins, no imputation.
Sample ordering is exploratory, not an inferential result. Euclidean distance on
call codes additionally assumes ordinal spacing of -2/-1/0/1; categorical modes
do not justify that distance. Prefer a scientifically justified continuous
measurement or an externally supplied row order when this assumption is unwanted.
The common-bin subset can be biased by systematic filtering. No global clustering
by default. Annotation colors are deterministic; supply fixed palettes when
comparing figures from different subsets.

## Paired directional agreement

`ichor_pair_concordance(a,b)` requires an explicit pair with unique IDs and the
same build. Calls and logR are aligned independently with the matrix BP rules
and separate coverage. Joint observed call support must also meet `min_coverage`;
disjoint partial source spans cannot become apparent co-observations. With relaxed
coverage, logR means can still summarize different observed portions of a bin.
All reference-grid bins remain. Missing, under-covered,
tied or mixed calls make a comparison unknown; a known call with missing logR
remains classified but has an unavailable height. No calls are inferred from
logR sign, fitted ploidy or filenames.

Both-neutral, a-only, b-only, same-direction altered (concordant), opposite-
direction altered (discordant), and unknown are distinct classes. Loss and deep
loss agree in direction, not amplitude. Baseline status is reported per sample
and chromosome class. Default sex policy `flag` retains source-call comparisons
but flags missing/ambiguous X/Y NEUT-bin evidence. Opt-in `require_neutral` makes
those X/Y comparisons unknown instead. No NEUT bins is not proof of an invalid
call: genuinely altered whole chromosomes can lack them. No baseline is guessed.

The paired plot shows original profiles above target-bin bars, not a numerical
subtraction. Default `height='both'` packs a/b into left/right half-bin bars when
both are altered, retaining both values; for one-sided calls only the altered
height contributes. `representative` instead uses the one-sided height, mean
for same-direction alterations, or larger absolute height for opposite directions
(ties choose a). Both values are required for a two-sided representative; its
sign can change when swapping tied discordant samples. Calls are never averaged.

Both-neutral bins are blank; grey marks unavailable calls/selected heights;
gold marks unresolved sex-reference evidence. Zero altered heights receive a
marker without nudging. Separate sample/agreement legends prevent color-semantic
mixing. Both panels share fixed y limits (default -2,2) and the selected raw or
upstream-adjusted scale. Out-of-view values warn; source values are retained.
Regions crop the display after aggregation, with original upper-bin midpoints.
Data/settings/provenance and view settings are attached to the ggplot object.
Relative `panel_heights` (profiles, agreement), `segment_linewidth` and
`point_stroke` only control layout/marks and are recorded in `ichor_view`.
Defaults preserve equal panels, width 0.55 and inherited ggplot point stroke.
Unequal panel proportions require ggplot2 >= 4.0.0; older supported versions error
on those requests, rather than ignore them. Equal panels work on ggplot2 3.5.
Default one-sided track colors follow the sample palette; an explicit agreement
palette overrides that derivation.
These are display policies, not validated biological agreement statistics; see
D19–D20 in the repository decision register.

## Reproducibility and privacy

Schema-versioned sample/cohort/matrix objects are ordinary R lists, frames and
matrices. Validation is provided for all three. Sample fingerprints (MD5 and byte
counts) are captured before/after reading; input drift aborts import. Provenance
retrieval never re-hashes changed files. Absolute paths are omitted by default;
retain_paths=TRUE opts in. Schema/build/package version and matrix transformation
settings travel with the objects. Original coordinates for terminal clipping are
retained. Saving objects is the caller's choice; there is no implicit cache.

MD5 is a reproducibility checksum, not a security guarantee. Genomic signals,
aliases, metadata, hashes and paths can still be sensitive. Metadata allowlisting
helps minimize accidental inclusion but does not anonymize a cohort. Do not
publish manuscript inputs, objects or session logs without institutional approval.
