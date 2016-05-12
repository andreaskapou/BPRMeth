#' Read file containing CpG island locations
#'
#' \code{read_encode_cgi} reads a file containing CpG island (CGI) locations
#'  in the human genome using the \code{\link[data.table]{fread}} function.
#'
#' @inheritParams read_bs_encode_haib
#'
#' @return A \code{\link[data.table]{data.table}} object.
#'
#' @seealso \code{\link{read_rna_encode_caltech}},
#'  \code{\link{read_bs_encode_haib}}
#'
#' @export
read_encode_cgi <- function(file, is_GRanges = TRUE){
  message("Reading file ", file, " ...")
  cgi_data <- data.table::fread(input = file,
                                sep = "\t",
                                header = FALSE,
                                col.names = c("chr", "start", "end", "cgi_id"))

  N <- NROW(cgi_data)
  cgi_id <- paste0("cgi_", seq(1, N))
  cgi_data$cgi_id <- cgi_id

  if (is_GRanges){
    # Create a GRanges object -----------------------------------
    message("Creating GRanges object ...")
    cgi_data <- GenomicRanges::GRanges(seqnames = cgi_data$chr,
                   ranges = IRanges::IRanges(start = cgi_data$start,
                                             end = cgi_data$end),
                   cgi_id = cgi_data$cgi_id)
  }
  message("Done!\n")
  return(cgi_data)
}
