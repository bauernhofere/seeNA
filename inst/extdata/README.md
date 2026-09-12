# Test-only format fixtures

These two tiny illustrative profiles and their manifest were authored in the
initial package commit c7849b6. They are not outputs of an ichorCNA fit, patient
data, biological simulations, or evidence about copy-number calling accuracy.
The values are deliberately small test cases for file parsing and visualization.

Tests modify copies in temporary directories to exercise malformed inputs,
missing calls, overlap weighting and identity mismatches. Real manuscript
compatibility/render checks use source files outside the package repository.
No test fixture should ever be substituted for manuscript analysis data.

Header reference: ichorCNA v0.4.0 R/output.R, commit
`d31ed52e9e9aa225084a6d380ba0a81a21e1943a`.
