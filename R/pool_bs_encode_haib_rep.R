#' Read and pool replicates from Bismark Cov formatted BS-Seq data
#'
#' \code{pool_bs_bismark_cov_rep} reads and pools replicate methylation data
#'  from  BS-Seq experiments that are in Bismark cov format.
#'
#' @inheritParams preprocess_bs_bismark_cov
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object.
#'
#' @export
pool_bs_encode_haib_rep <- function(files, chr_discarded = NULL){
  assertthat::assert_that(length(files) > 1)
  message("Pooling BS-Seq replicates ...")

  # Read first file
  pooled_bs <- read_bs_encode_haib(file          = files[1],
                                   chr_discarded = chr_discarded,
                                   is_GRanges    = TRUE)
  for (i in 2:length(files)){
    # Read replicate file i
    bs_data_rep <- read_bs_encode_haib(file          = files[i],
                                       chr_discarded = chr_discarded,
                                       is_GRanges    = TRUE)

    # Find overlaps between BS-Seq replicates.
    # A Hits object containing in the 1st column the query indices and in
    # the 2nd column the corresponding subject indices that overlap.
    overlaps <- GenomicRanges::findOverlaps(query   = pooled_bs,
                                            subject = bs_data_rep)

    # Get only the subset of overlapping CpG sites
    tmp_bs <- pooled_bs[S4Vectors::queryHits(overlaps)]

    # Add the total reads from the two distinct replicates
    tmp_bs$total_reads <- tmp_bs$total_reads +
                  bs_data_rep[S4Vectors::subjectHits(overlaps)]$total_reads

    # Add the methylated reads from the two distinct replicates
    tmp_bs$meth_reads <- tmp_bs$meth_reads +
                  bs_data_rep[S4Vectors::subjectHits(overlaps)]$meth_reads

    # Add CpG sites that were not overlapping between replicates
    tmp_bs <- c(tmp_bs,
                pooled_bs[-S4Vectors::queryHits(overlaps)],
                bs_data_rep[-S4Vectors::subjectHits(overlaps)])

    # Sort the pooled GRanges object
    pooled_bs <- sort(tmp_bs, ignore.strand=TRUE)
  }
  message("Finished pooling BS-Seq replicates!\n")
  return(pooled_bs)
}
