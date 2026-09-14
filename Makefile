.PHONY: bootstrap document readme audit-reference test check install

bootstrap:
	Rscript -e 'if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak", repos = "https://cloud.r-project.org"); pak::pkg_install(c("devtools", paste0("roxygen2@", read.dcf("DESCRIPTION")[1, "RoxygenNote"]))); pak::local_install_dev_deps()'

document:
	Rscript -e 'roxygen2::roxygenise()'

readme:
	Rscript -e 'knitr::knit("README.Rmd", output = "README.md", quiet = TRUE)'

audit-reference:
	Rscript docs/check-reference-lengths.R

test:
	Rscript -e 'devtools::test()'

check: document
	Rscript -e 'devtools::check(document = FALSE, error_on = "warning")'

install: document
	R CMD INSTALL .
