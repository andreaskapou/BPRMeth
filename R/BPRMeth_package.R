#' @title \code{BPRMeth}: Extracting higher order methylation features
#'
#' @description Higher order methylation features for clustering and prediction
#'   in epigenomic studies
#' @docType package
#' @name BPRMeth
#'
#' @return BPRMeth main package documentation.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @rawNamespace importFrom(magrittr,"%>%")
#' @rawNamespace importFrom(data.table,":=")
#' @import GenomicRanges ggplot2
#' @importFrom stats pnorm dbinom dnorm
#' @importFrom Rcpp evalCpp
#' @importFrom cowplot plot_grid
#' @useDynLib BPRMeth
#'
.datatable.aware <- TRUE
NULL
#> NULL


.onLoad <- function(libname = find.package("BPRMeth"), pkgname = "BPRMeth"){
    # CRAN Note avoidance
    if (getRversion() >= "2.15.1")
        utils::globalVariables(
            # sample file names from taxstats
            c(# we use the magrittr pipe
                "."
            )
        )
    invisible()
}
