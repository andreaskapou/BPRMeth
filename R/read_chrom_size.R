#' Read genome chromosome sizes file.
#'
#' \code{read_chrom_size} reads a file containing genome chromosome sizes using
#'  the \code{\link[data.table]{fread}} function.
#'
#' @param file The name of the file to read data values from.
#'
#' @return A \code{\link[data.table]{data.table}} object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{read_rna_encode_caltech}},
#'  \code{\link{read_bs_encode_haib}}
#'
#' @examples
#' chr_file <- system.file("extdata", "hg19.chr.sizes", package = "BPRMeth")
#' chr_data <- read_chrom_size(chr_file)
#'
#' # Extract the size of the chr1
#' chr_data[1]
#'
#' @export
read_chrom_size <- function(file){
  message("Reading file ", file, " ...")
  data_raw <- data.table::fread(input = file,
                                sep = "\t",
                                header = FALSE,
                                col.names = c("chr", "size"))

  message("Done!\n")
  return(data_raw)
}
