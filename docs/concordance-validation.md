# Paired concordance feature — local validation

Scope: development working tree, version 0.0.0.9003. This extends the existing
README/help/identifier/vdiffr work without changing manuscript code or patient
inputs. No public release, manuscript approval or new CI result is implied.

## Source review and intentional differences

The manuscript overlay joins exact coordinates, drops missing/unmatched calls,
uses corrected-call direction, then draws a one-sided/mean/largest-absolute logR
height (absolute ties select the first fluid). The useful feature is retained,
but not copied as an implicit numerical difference or diagnostic truth:

- Explicit same-build pair; no sample/fit/fluid inference.
- Reference-grid alignment with separate call/logR coverage and joint observed
  call support. Mixed/tied/unsupported comparisons are unknown, not dropped.
- Both-neutral blanks are distinct from grey unavailable evidence.
- Default two-height display preserves both altered values. Representative
  heights are an explicit option, including the documented tie asymmetry.
- X/Y reference evidence is flagged independently of source-call direction.
  Opt-in conservative gating is available; absence of NEUT bins alone does not
  disprove an altered whole chromosome. No named-sample corrections or fallbacks.
- Original upper-bin midpoints; separate sample/agreement legends; shared limits
  with clipping warnings, explicit scale, fitted parameters, and privacy opt-out.

D19 and D20 in [the decision register](decision-register.md) give the policies,
alternatives and limitations. The [methods](../inst/methods.md) define the contract.

## Automated evidence

- The suite covers the pre-existing profile, comparison and region APIs alongside
  the paired track. Visual baselines cover ordinary, styled and unresolved
  sex-reference views. Test totals are reported by the run, not maintained here.
- Numeric and device tests run on ggplot2 3.5 as well as current versions; visual
  snapshots target ggplot2 4. Unequal panel sizing is tested on supporting versions
  and explicitly rejected on older versions. CI pins a ggplot2 3.5.2 job.
- Added numeric regressions cover all agreement classes, logR sign independence,
  deep-loss versus loss, unequal/shifted grids, missing calls versus missing
  heights, mixed winning modes, disjoint partial support, representative ties,
  missing X reference evidence, conservative sex policy, zero-height markers,
  parameter/build/ID/allocation errors, clipping, original midpoints and legends.
- `make readme` executes the expanded gallery, including its fourth fixture-only
  PNG. `make check` passes with zero errors/warnings and the existing environment
  NOTE `unable to verify current time`. A PDF-device smoke test also renders the
  new track, closing the device before the byte-size assertion.

## Real-input reconciliation (private, local only)

One explicitly selected pair from the earlier canonical smoke-test manifest was
loaded. An independent exact-coordinate source join and directly written
call-direction/height arithmetic matched the new output on all **2,548 shared
source bins matching target coordinates**. This validates the representative
rule on comparable bins, not equivalence of unmatched/mixed support policies.

The pair table contained 3,044 target bins (default autosomes/X, excluding Y),
including unknown bins. Warm pair construction was about 0.20 seconds locally;
this is a smoke timing, not a general benchmark. Both two-height and representative
regional real-input plots rendered and were visually inspected. Full clinical
records and rendered images remain outside the package repository.

This is not an adjudication of alternative fitted runs or the user's individual
sex-chromosome cases. Biological validity, final sample selection and investigator
approval remain separate from API/mathematical/rendering validation.
