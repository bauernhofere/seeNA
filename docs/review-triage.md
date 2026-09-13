# Independent-review triage and manuscript release gate

The initial commit is a prototype, not a scientifically validated release.
Passing R CMD check or parsing files without errors does not establish that
calls, intervals, pairing, annotations, or provenance are correct.

## Accepted priorities

1. Preserve missing calls/flags; reject malformed fields and ambiguous schemas.
2. Match source identities before applying user aliases; validate segments too.
3. Never average categorical call codes. Return coverage and mixed/tied states.
4. Freeze input fingerprints at ingestion and omit source paths by default.
5. Fail cohort loading atomically; validate metadata/matrix identity and order.
6. Preserve manifest order in heatmaps; explicitly define clustering missingness.
7. Exercise rendering and hand-calculated scientific regressions, not just types.
8. Document methods, evidence, limitations, and a pinned manuscript workflow.

## Reviewer claims that need qualification

- ichorCNA v0.4 output.R writes our existing `sample/chr/start/end/.../median`
  format to `.seg`. Its separate `.seg.txt` export uses `ID/chrom/.../
  seg.median.logR`. Do not claim the original parser rejected canonical `.seg`.
- Legitimate logR values may be missing. Never require every bin to be finite.
- Clinical plot annotations supplied intentionally in a manifest are in scope.
  Metadata privacy is a user responsibility; a regex is not a de-identification
  guarantee. The package must not import unrelated clinical spreadsheets.
- Chromosome bounds require care: the real canonical input files contain
  terminal fixed-width bins extending beyond the reference chromosome length.
  The initial hardening used strict rejection with opt-in recorded trimming.
  This was superseded by source-grid-supported `bounds='window'` after the
  downstream port; see D05 in [decision-register.md](decision-register.md).
  Neither policy verifies assembly identity.
- An agent-friendly package is an explicit API, small helpers, documented
  scientific invariants, and executable tests—not a dependence on a particular
  model, harness, or guessed API examples.
- The earlier generated fixtures are test-only examples, not biological data or
  evidence of ichorCNA accuracy. Manuscript verification must use real inputs
  outside this repository, without patient-specific code corrections.

## Manuscript heatmap is a different measurement

**Historical snapshot:** the initial review found corrected copy number minus
fitted ploidy, not categorical Corrected_Call. A fresh audit found the live port
now subtracts a NEUT-bin reference, with wrapper-specific fallback assumptions.
See [the live-port audit](decision-register.md#live-manuscript-port-audit); the
measurement change is not automatically approved by passing package checks.

The earlier manuscript display path was recorded as follows. Its display path clips values at
+/-2, blanks deviations below 0.35, and saturates the tighter palette at +/-1.
It also orders paired rows by plasma alteration burden. These choices must not
be silently replaced by the package's categorical call heatmap or an absolute-CN
heatmap. This remains true for the newer NEUT-reference measurement. The display
dead band is not an ichorCNA neutrality call.

`load_cn_landscape()` also catches failed pair loads and returns NULL, allowing
pairs to disappear. The next manuscript integration pass must reconcile the
expected versus plotted pair set explicitly and document the final input-run
manifest. Keep study-specific pair selection outside the package; implement any
CN-minus-ploidy transformation as an explicitly named, tested measurement rather
than changing what existing matrix values mean. Final selection, display
thresholds, ordering and scale limits need author approval.

## Manuscript release gate

Keep the repository private until the author approves publication. Before a
citable release: complete regression tests and checks, rerun the actual figure
workflow against source data, visually inspect plots, reconcile matrix/figure
semantics with the manuscript, document source ichorCNA run selection, and freeze
the package commit and analysis inputs. A private GitHub URL is not an accessible
software release; public release/archive and citation metadata require author
approval. No release tag, DOI, or numerical results should be invented.
