#' Read file containing Bismark Cov formatted BS-Seq data
#'
#' \code{read_bs_bismark_cov} reads a file containing methylation data from
#'  BS-Seq experiments using the \code{\link[data.table]{fread}} function.
#'  The BS-Seq file should be in Bismark Cov format.
#'
#' @inheritParams read_bs_encode_haib
#'
#' @return a \code{\link[GenomicRanges]{GRanges}} object if \code{is_GRanges}
#'  is TRUE, otherwise a \code{\link[data.table]{data.table}} object.
#'
#' @seealso \code{\link{pool_bs_bismark_cov_rep}},
#'  \code{\link{preprocess_bs_bismark_cov}}
#'
#' @references
#'   \url{http://rnbeads.mpi-inf.mpg.de/data/RnBeads.pdf}
#'
#' @export
read_bs_bismark_cov <- function(file, chr_discarded = NULL, is_GRanges = TRUE){
  message("Reading file ", file, " ...")
  bs_data <- data.table::fread(input = file,
                               sep = "\t",
                               header = FALSE,
                               col.names = c("chr", "start", "meth_reads",
                                             "unmeth_reads"))


  # Remove selected chromosomes  -------------------------------
  bs_data <- discard_chr(x = bs_data, chr_discarded = chr_discarded)


  # Sorting data -----------------------------------------------
  # With order priority: 1. chr, 2. start
  message("Sorting BS-Seq data ...")
  bs_data <- bs_data[order(bs_data$chr, bs_data$start)]


  if (is_GRanges){
    # Create a GRanges object ---------------------------------
    message("Creating GRanges object ...")
    bs_data <- GenomicRanges::GRanges(seqnames = bs_data$chr,
                      ranges = IRanges::IRanges(start=bs_data$start, width=1),
                      meth_reads   = bs_data$meth_reads,
                      unmeth_reads = bs_data$unmeth_reads)
  }
  message("Done!\n")
  return(bs_data)
}
