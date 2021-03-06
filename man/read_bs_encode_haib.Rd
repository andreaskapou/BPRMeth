% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated_process_functions.R
\name{read_bs_encode_haib}
\alias{read_bs_encode_haib}
\title{(DEPRECATED) Read ENCODE HAIB bed formatted BS-Seq file}
\usage{
read_bs_encode_haib(file, chr_discarded = NULL, is_GRanges = TRUE)
}
\arguments{
\item{file}{The name of the file to read data values from.}

\item{chr_discarded}{A vector with chromosome names to be discarded.}

\item{is_GRanges}{Logical: if TRUE a GRanges object is returned, otherwise a
data.frame object is returned.}
}
\value{
A \code{GRanges} object if \code{is_GRanges} is
  TRUE, otherwise a \code{data.table} object.

  The GRanges object contains two additional metadata columns: \itemize{
  \item \code{total_reads}: total reads mapped to each genomic location.
  \item \code{meth_reads}: methylated reads mapped to each genomic location.
  } These columns can be accessed as follows:
  \code{granges_object$total_reads}
}
\description{
(DEPRECATED) \code{read_bs_encode_haib} reads a file containing
  methylation data from BS-Seq experiments using the \code{\link{scan}}
  function. The BS-Seq file should be in ENCODE HAIB \code{bed} format. Read
  the Important section below on when to use this function.
}
\section{Important}{
 Unless you want to create a different workflow when
  processing the BS-Seq data, you should NOT call this function, since this
  is a helper function. Instead you should call the
  \code{\link{preprocess_bs_seq}} function.
}

\examples{
# Obtain the path to the file and then read it
bs_file <- system.file("extdata", "rrbs.bed", package = "BPRMeth")
bs_data <- read_bs_encode_haib(bs_file)

}
\seealso{
\code{\link{pool_bs_seq_rep}}, \code{\link{preprocess_bs_seq}}
}
\author{
C.A.Kapourani \email{C.A.Kapourani@ed.ac.uk}
}
