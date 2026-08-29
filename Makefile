.PHONY: document test check install

document:
	Rscript -e 'roxygen2::roxygenise()'

test:
	Rscript -e 'devtools::test()'

check: document
	Rscript -e 'devtools::check(document = FALSE, error_on = "warning")'

install: document
	R CMD INSTALL .
