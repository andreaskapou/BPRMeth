#' Read file containing CpG island locations
#'
#' \code{read_encode_cgi} reads a file containing CpG island (CGI) locations in
#' the human genome using the \code{\link[data.table]{fread}} function.
#'
#' @inheritParams read_bs_encode_haib
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object if \code{is_GRanges} is
#'   TRUE, otherwise a \code{\link[data.table]{data.table}} object.
#'
#'   The GRanges object contains one additional metadata column: \itemize{ \item
#'   \code{cgi_id}: Unique ID of the CpG Island. } This column can be accessed
#'   as follows: \code{granges_object$cgi_id}
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{read_rna_encode_caltech}},
#'   \code{\link{read_bs_encode_haib}}
#'
#' @examples
#' \dontrun{
#' # Download the file and change the working directory to that location
#' file <- "name_of_CGI_file")
#' cgi_data <- read_encode_cgi(file)
#'
#' # Extract the CGI ID
#' cgi_id <- cgi_data$cgi_id
#' }
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
                                             end   = cgi_data$end),
                   cgi_id = cgi_data$cgi_id)
  }
  message("Finished reading CGI file!\n")
  return(cgi_data)
}
