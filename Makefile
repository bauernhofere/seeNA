.PHONY: bootstrap document test check install

bootstrap:
	Rscript -e 'if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak", repos = "https://cloud.r-project.org"); pak::local_install_dev_deps()'

document:
	Rscript -e 'roxygen2::roxygenise()'

test:
	Rscript -e 'devtools::test()'

check: document
	Rscript -e 'devtools::check(document = FALSE, error_on = "warning")'

install: document
	R CMD INSTALL .
