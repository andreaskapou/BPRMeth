% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_data.R
\name{read_expr}
\alias{read_expr}
\title{Read expression data file}
\usage{
read_expr(file, log2_transf = FALSE, delimiter = "\\t")
}
\arguments{
\item{file}{File name}

\item{log2_transf}{Logical, whether to log2 transform the expression data.}

\item{delimiter}{Delimiter format the columns are splitted. Default is tab.}
}
\value{
A \code{\link[data.table]{data.table}} object.
}
\description{
\code{read_expr} reads a file containing expression data. It
  should contain two columns: "id", "expr" and by default columns are
  considered in tab delimited format.
}
\examples{
# Obtain the path to files
file <- system.file("extdata", "dummy_expr.bed", package = "BPRMeth")
expr_dt <- read_expr(file)

# Extract feature id
expr_ids <- expr_dt$id

}
\seealso{
\code{\link{read_met}}, \code{\link{read_anno}}
}
\author{
C.A.Kapourani \email{C.A.Kapourani@ed.ac.uk}
}
