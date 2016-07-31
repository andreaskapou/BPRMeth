## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"---------------------------------
BiocStyle::latex()

## ----eval=TRUE, echo=TRUE, message=FALSE, results="asis"--------------------------------
library(BPRMeth)

## ----eval=TRUE, echo=TRUE, results="asis"-----------------------------------------------
rrbs_file <- system.file("extdata", "rrbs.bed", package = "BPRMeth")
rnaseq_file <- system.file("extdata", "rnaseq.bed", package = "BPRMeth")

## ----eval=TRUE, echo=TRUE, message=FALSE, results="asis"--------------------------------
HTS_data <- process_haib_caltech_wrap(rrbs_file, rnaseq_file)

