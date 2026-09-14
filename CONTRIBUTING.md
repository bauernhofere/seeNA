# Contributing

1. Create a focused branch.
2. Add or update tests for behavior changes.
3. Run `make check` before opening a pull request.
4. Do not commit patient data or real clinical sample identifiers.
5. Keep genome build, coordinate, call-column, and tumor-fraction semantics
   explicit. Update the affected entry in [the decision register](docs/decision-register.md)
   with rationale, evidence, alternatives and limitations; a passing test does
   not establish biological validity.
6. Edit `README.Rmd`, not generated `README.md`, and run `make readme`. Inspect
   all four regenerated fixture-only PNGs; never substitute clinical profiles.

`make audit-reference` optionally compares the implemented primary chromosome
lengths against public UCSC tables. It requires network access and does not
verify the assembly of user inputs.

Bug reports should include the ichorCNA version, file headers, genome build, and
a minimal de-identified reproduction. Do not attach protected clinical data.
