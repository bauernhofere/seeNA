# Public development repository preflight

The author approved public repository visibility. This does not create a tagged
release, archive or DOI, validate clinical conclusions, or approve manuscript
figures. Existing manuscript integration and scientific review gates remain.

## Checks performed before changing visibility

- Inspected tracked and untracked source paths and paths across reachable Git
  history, including remote branches. Historical binary paths were limited to
  the four fixture-only README images.
- Scanned reachable Git history (`--all`) and the current source tree with
  Gitleaks 8.30.1; no secrets were detected. The scanner archive checksum was
  verified against its published release checksums. Raw reports remain local.
- Scanned historical text blobs and current sources for the study's identifier
  pattern and user-specific local paths; no matches were found. These are
  heuristic checks, not de-identification or a guarantee against all disclosures.
- Reviewed GitHub issue/PR descriptions and comments. No clinical inputs or
  credentials were found. Author identity/contact in Git and package citation
  metadata is intentional public attribution.
- Retained only hand-authored format fixtures, their explicit manifest, and
  fixture-derived images/snapshots. Real-input validation data remain external.
- Ran the test suite, regenerated the README/help, and completed R CMD check;
  no errors or warnings, with the existing clock-verification NOTE only.

## Ignore and contributor-guide policy

`.gitignore` now covers local credentials/environment files, agent/editor state,
clinical tables and common genomic exports, fitted segment files, serialized
objects, renderings, build/check output and archives. Exceptions name each public
input fixture explicitly; there is no blanket exemption for segment files under
`inst/extdata`. The repository-only regression checks these rules without relying
on the maintainer's personal Git exclusions.

`AGENTS.md` stays tracked: it contains project architecture, scientific invariants
and privacy/review instructions, not private session notes. `.Rbuildignore`
excludes it, local configuration and local data/output directories from the R
package build. Public visibility is not authorization for future data uploads.

Ignore rules do not untrack existing files, erase history, inspect file contents,
or prevent `git add -f`. Review future staged content, fixture exceptions and
history rather than treating these patterns or a clean secret scan as a security
guarantee.
