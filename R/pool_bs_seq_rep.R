#' Read and pool replicates from BS-Seq data
#'
#' \code{pool_bs_seq_rep} reads and pools replicate methylation data from BS-Seq
#' experiments that are either in Encode RRBS or Bismark Cov format. Read the
#' Important section below on when to use this function.
#'
#' @inheritParams preprocess_bs_seq
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object. The GRanges object
#'   contains two additional metadata columns: \itemize{ \item
#'   \code{total_reads}: total reads mapped to each genomic location. \item
#'   \code{meth_reads}: methylated reads mapped to each genomic location. }
#'   These columns can be accessed as follows: granges_object$total_reads
#'
#' @section Important: Unless you want to create a different workflow when
#'   processing the BS-Seq data, you should NOT call this function, since this
#'   is a helper function. Instead you should call the
#'   \code{\link{preprocess_bs_seq}} function.
#'
#'   Information about the file formats can be found in the following links:
#'
#'   Encode RRBS format:
#'   \url{http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgTables?db=hg19&hgta_group=regulation&hgta_track=wgEncodeHaibMethylRrbs&hgta_table=wgEncodeHaibMethylRrbsBcbreast0203015BiochainSitesRep2&hgta_doSchema=describe+table+schema}
#'
#'   Bismark Cov format: \url{http://rnbeads.mpi-inf.mpg.de/data/RnBeads.pdf}
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{read_bs_bismark_cov}},
#'   \code{\link{read_bs_encode_haib}}, \code{\link{preprocess_bs_seq}}
#'
#' @examples
#' \dontrun{
#' # Download the files and change the working directory to that location
#' files <- c("name_of_bs_encode_file1", "name_of_bs_encode_file2")
#' bs_data <- pool_bs_seq_rep(files, file_format = "encode_rrbs")
#'
#' # Extract the total reads and methylated reads
#' total_reads <- bs_data$total_reads
#' meth_reads <- bs_data$meth_reads
#' }
#'
#' @export
pool_bs_seq_rep <- function(files, file_format = "encode_rrbs",
                            chr_discarded = NULL){
  assertthat::assert_that(length(files) > 1)
  message("Pooling BS-Seq replicates ...")

  # Read first file
  if (file_format == "encode_rrbs"){
    pooled_bs <- read_bs_encode_haib(file          = files[1],
                                     chr_discarded = chr_discarded,
                                     is_GRanges    = TRUE)
  }else if (file_format == "bismark_cov"){
    pooled_bs <- read_bs_bismark_cov(file          = files[1],
                                     chr_discarded = chr_discarded,
                                     is_GRanges    = TRUE)
  }else{
    stop("Wrong file format. Please check the available file formats!")
  }
  for (i in 2:length(files)){
    # Read replicate file i
    if (file_format == "encode_rrbs"){
      bs_data_rep <- read_bs_encode_haib(file          = files[i],
                                         chr_discarded = chr_discarded,
                                         is_GRanges    = TRUE)
    }else if (file_format == "bismark_cov"){
      bs_data_rep <- read_bs_bismark_cov(file          = files[i],
                                         chr_discarded = chr_discarded,
                                         is_GRanges    = TRUE)
    }

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
    pooled_bs <- sort(tmp_bs, ignore.strand = TRUE)
  }
  message("Finished pooling BS-Seq replicates!\n")
  return(pooled_bs)
}
