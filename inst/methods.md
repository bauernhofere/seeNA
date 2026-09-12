# Scientific contracts (pre-release schema 1)

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

Some real ichorCNA fixed-width terminal bins exceed reference lengths. Default
`bounds = 'error'` rejects them. Explicit `bounds = 'trim'` clips only intervals
that straddle a reference chromosome end, warns, and stores role/row/original/new
coordinates in `coordinate_changes`. Intervals wholly outside the reference still
fail. Trimming is not evidence that the declared build is correct.

A sample object reconciles source IDs across all supplied files before applying
an alias. Aliases replace source IDs in returned tables. With unprefixed bin
columns, source identity falls back to the .cna.seg basename. Renamed standalone
files therefore need matching source naming. Aliases are not multi-sample selectors.
Identity agreement cannot distinguish multiple runs for the same sample: the
caller must explicitly choose one run for all components. Never select a run by
matching a desired manuscript TF inside the package.

## Measurements and categorical calls

- `logR`: source log2-ratio signal; never re-centered or purity-adjusted here.
- `copy_number`, `event`: raw ichorCNA output fields.
- `corrected_copy_number`, `corrected_call`: caller-corrected output fields, not
  recomputed by ichorViz. Corrected calls, absolute CN, logR sign and neutrality
  relative to fitted ploidy must not be conflated.
- `logR_copy_number`: diagnostic only.
- TF: fraction [0,1]; ploidy positive. Missing values remain missing.
- Logical subclone_status recognizes explicit true/false tokens; unknowns error.

Call vocabulary: HOMD -> Deep loss (-2), HETD -> Loss (-1), NEUT -> Neutral (0),
GAIN/AMP/HLAMP and numbered HLAMP variants -> Gain (1). NA/empty -> missing;
other tokens error. There is no fallback from corrected_call to event.

## Fixed-bin matrices

`value` is required. `call_column` explicitly chooses corrected_call (default)
or event when value='call'. Metadata/sample order is validated at every boundary.

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
inputs, temporary calculations, bin labels and metadata.

## Plotting

Profiles color bins by the specified call field; segment medians are grey to
avoid merging raw segment events into a corrected-call legend. Missing calls use
grey, never neutral blue. Comparison overlays use original genomic midpoints,
not aligned/rebinned observations; only segment endpoints are clipped in region
views. Common genome axes use reference lengths even for sparse inputs. No
vertical nudging, hidden downsampling or biological reclassification.

Heatmaps preserve manifest order by default; identifiers are hidden unless
requested. Categorical colors are discrete. Continuous default scales include
observed extremes (CN 2 is a reference color, not a neutrality decision). Custom
colors may intentionally saturate: authors must report their scale. Unknown,
tied and insufficient-coverage cells use the NA color. Inspect coverage/mixed
layers to distinguish these causes. Rasterized figures have finite display
resolution; use a regional view or more pixels for bin-level inspection.

Optional clustering: 2-2000 samples, at least two bins observed in every sample,
Euclidean distance and complete linkage on these common bins, no imputation.
Sample ordering is exploratory, not an inferential result. No global clustering
by default. Annotation colors are deterministic; supply fixed palettes when
comparing figures from different subsets.

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
