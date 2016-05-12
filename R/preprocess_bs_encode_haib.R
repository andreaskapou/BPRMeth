#' Process BS-Seq data in ENCODE HAIB format
#'
#' \code{preprocess_bs_encode_haib} pre-processes BS-Seq files in ENCODE HAIB
#'  format. If a vector of files is given, the replicates are pooled together,
#'  and finally noisy reads are discarded.
#'
#' @inheritParams preprocess_bs_bismark_cov
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object.
#'
#' @export
preprocess_bs_encode_haib <- function(files, chr_discarded = NULL,
                                      min_bs_cov = 2, max_bs_cov = 1000){

  # Read BS-Seq data in ENCODE HAIB format ---------------
  if (length(files) > 1){
    bs_data <- pool_bs_encode_haib_rep(files         = files,
                                       chr_discarded = chr_discarded)
  }else{
    bs_data <- read_bs_encode_haib(file          = files,
                                   chr_discarded = chr_discarded,
                                   is_GRanges    = TRUE)
  }

  bs_data <- .discard_bs_noise_reads(bs_data     = bs_data,
                                     min_bs_cov  = min_bs_cov,
                                     max_bs_cov  = max_bs_cov)
  message("Done!\n")
  return(bs_data)
}
