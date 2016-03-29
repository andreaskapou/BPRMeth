#' Process BS-Seq data in Bismark cov format
#'
#' \code{preprocess_bs_bismark_cov} pre-processes BS-Seq files in Bismark Cov
#'  format. If a vector of files is given, the replicates are pooled together,
#'  then the total reads are extracted and finally noisy reads are discarded.
#'
#' @param files A vector of filenames containing replicate experiments.
#' @param chr_discarded A vector with chromosome names to be discarded.
#' @param min_bs_cov Minimum number of reads mapping to each CpG site.
#' @param max_bs_cov Maximum number of reads mapping to each CpG site.
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object.
#'
#' @export
preprocess_bs_bismark_cov <- function(files, chr_discarded = NULL,
                                      min_bs_cov = 2, max_bs_cov = 1000){

  # Read BS-Seq data in Bismark Cov format ---------------
  if (length(files) > 1){
    bs_data <- pool_bs_bismark_cov_rep(files         = files,
                                       chr_discarded = chr_discarded)
  }else{
    bs_data <- read_bs_bismark_cov(file          = files,
                                   chr_discarded = chr_discarded,
                                   is_GRanges    = TRUE)
  }

  # Get total reads from meth_reads + unmeth_reads -------
  total_reads <- bs_data$meth_reads + bs_data$unmeth_reads

  # Create final GRanges object --------------------------
  message("Creating final GRanges object for BS-Seq data ...")
  bs_data <- GenomicRanges::GRanges(seqnames    = bs_data@seqnames,
                                    ranges      = bs_data@ranges,
                                    total_reads = total_reads,
                                    meth_reads  = bs_data$meth_reads)

  bs_data <- discard_bs_noise_reads(bs_data     = bs_data,
                                    min_bs_cov  = min_bs_cov,
                                    max_bs_cov  = max_bs_cov)
  message("Done!\n")
  return(bs_data)
}
