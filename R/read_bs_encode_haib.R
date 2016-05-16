#' Read ENCODE HAIB bed formatted BS-Seq file
#'
#' \code{read_bs_encode_haib} reads a file containing methylation data from
#' BS-Seq experiments using the \code{\link{scan}} function. The BS-Seq file
#' should be in ENCODE HAIB \code{bed} format. Read the Important section below
#' on when to use this function.
#'
#' @param file The name of the file to read data values from.
#' @param chr_discarded A vector with chromosome names to be discarded.
#' @param is_GRanges Logical: if TRUE a GRanges object is returned, otherwise a
#'   data.frame object is returned.
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object if \code{is_GRanges} is
#'   TRUE, otherwise a \code{\link[data.table]{data.table}} object.
#'
#'   The GRanges object contains two additional metadata columns: \itemize{
#'   \item \code{total_reads}: total reads mapped to each genomic location.
#'   \item \code{meth_reads}: methylated reads mapped to each genomic location.
#'   } These columns can be accessed as follows:
#'   \code{granges_object$total_reads}
#'
#' @section Important: Unless you want to create a different workflow when
#'   processing the BS-Seq data, you should NOT call this function, since this
#'   is a helper function. Instead you should call the
#'   \code{\link{preprocess_bs_seq}} function.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @references
#' \url{http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgTables?db=hg19&hgta_group=regulation&hgta_track=wgEncodeHaibMethylRrbs&hgta_table=wgEncodeHaibMethylRrbsBcbreast0203015BiochainSitesRep2&hgta_doSchema=describe+table+schema}
#'
#' @seealso \code{\link{pool_bs_seq_rep}}, \code{\link{preprocess_bs_seq}}
#'
#' @examples
#' \dontrun{
#' # Download the files and change the working directory to that location
#' file <- "name_of_bs_encode_file"
#' bs_data <- read_bs_encode_haib(file)
#'
#' # Extract the total reads and methylated reads
#' total_reads <- bs_data$total_reads
#' meth_reads <- bs_data$meth_reads
#' }
#' @export
read_bs_encode_haib <- function(file, chr_discarded = NULL, is_GRanges = TRUE){
  message("Reading file ", file, " ...")
  data_raw <- scan(file = file,
                   skip = 1,
                   sep = "\t",
                   what = list("character",  # Reference chromosome or scaffold
                               integer(),    # Start position in chromosome
                               NULL,         # End position in chromosome
                               NULL,         # Name of item
                               integer(),    # Score from 0-1000. Capped number
                               "character",  # Strand : + or - or . for unknown
                               NULL,         # Start position
                               NULL,         # End position
                               NULL,         # Color value R,G,B
                               NULL,         # Number of reads or coverage
                               integer()     # Methylation percentage
                   ))


  # Convert to actual methylated reads -------------------------
  data_raw[[11]] <- as.integer(round(0.01 * data_raw[[5]] * data_raw[[11]]))


  # Store only required fields
  bs_data <- data.table::data.table(chr = data_raw[[1]],
                                    start = data_raw[[2]],
                                    strand = data_raw[[6]],
                                    total_reads = data_raw[[5]],
                                    meth_reads = data_raw[[11]])
  rm(data_raw)


  # Remove selected chromosomes  -------------------------------
  bs_data <- .discard_chr(x = bs_data, chr_discarded = chr_discarded)


  # Sorting data -----------------------------------------------
  # With order priority: 1. chr, 2. start, 3. strand
  message("Sorting BS-Seq data ...")
  bs_data <- bs_data[order(bs_data$chr, bs_data$start, bs_data$strand)]


  if (is_GRanges){
    # Create a GRanges object -----------------------------------
    message("Creating GRanges object ...")
    bs_data <- GenomicRanges::GRanges(seqnames = bs_data$chr,
                    strand = bs_data$strand,
                    ranges = IRanges::IRanges(start = bs_data$start, width = 1),
                    total_reads = bs_data$total_reads,
                    meth_reads  = bs_data$meth_reads)
  }
  message("Finished reading BS-Seq file!\n")
  return(bs_data)
}
