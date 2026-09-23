# Synthetic gallery data

`make-demo.R` defines `write_demo_cohort(directory = tempfile("seena-demo-"))`.
It writes a manifest and bin, segment and parameter files for 24 synthetic
samples, then returns the manifest path. Use a new or empty directory; existing
files are never overwritten. The script is distributed under the package's
GPL-3-or-later license.

## Recipe

- hg38 reference lengths define a 500 kb grid on chromosomes 1–22 and X. Terminal
  intervals end at the chromosome boundary.
- Three explicitly authored copy-state templates, with pair-specific breakpoints
  and focal changes, give each example group a recognizable pattern. The sample
  IDs, groups, pair IDs, timepoints and TF values are all invented.
- Bin logR is `log2((2 * (1 - TF) + CN * TF) / 2)` plus Gaussian variation with
  SD 0.055. The generator fixes seed 41027 and all RNG kinds, and restores the
  caller's RNG state on exit.
- Calls are assigned from the designed integer CN states, not inferred from
  noisy logR. Raw and corrected calls are intentionally identical.
- Segment medians summarize the generated bins within contiguous state runs.
  Deliberately omitted intervals split segments as well as bins, preserving gaps
  in both types of matrix. The first pair includes shared gain, one-sided loss,
  one-sided gain, opposite directions and missing support on chromosome 1.
- Parameter files contain specified TF and nominal ploidy 2, not fitted estimates.
  This is an illustrative mixture formula, not a model of tumor evolution,
  realistic sequencing noise, sex-chromosome biology or ichorCNA estimation.

The files follow the supported ichorCNA input schemas but **were not produced by
an ichorCNA fit**. They contain no patient data and are not calibrated to the
manuscript cohort. The generator uses only the package's public chromosome-length
helper and base/recommended R packages; it reads no clinical inputs or network
resources. Generated raw files stay in the caller's temporary directory and are
not committed. The original small `inst/extdata` parser fixtures remain separate.
