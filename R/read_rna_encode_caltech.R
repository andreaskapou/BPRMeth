#' Read ENCODE Caltech bed formatted RNA-Seq file
#'
#' \code{read_rna_encode_caltech} reads a file containing promoter annotation
#' data together with gene expression levels from RNA-Seq experiments using the
#' \code{\link{scan}} function. The RNA-Seq file should be in ENCODE Caltech
#' \code{bed} format, e.g. use \code{gtf2bed} tool if your initial file is in
#' \code{gtf} format.
#'
#' @inheritParams read_bs_encode_haib
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object if \code{is_GRanges} is
#'   TRUE, otherwise a \code{\link[data.table]{data.table}} object.
#'
#'   The GRanges object contains three additional metadata columns: \itemize{
#'   \item \code{ensembl_id}: Ensembl IDs of each gene promoter. \item
#'   \code{gene_name}: Gene name. \item \code{gene_fpkm}: Expression level in
#'   FPKM. } These columns can be accessed as follows:
#'   \code{granges_object$ensembl_id}
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{read_chrom_size}}, \code{\link{read_bs_encode_haib}}
#'
#' @examples
#' # Obtain the path to the file and then read it
#' rnaseq_file <- system.file("extdata", "rnaseq.bed", package = "BPRMeth")
#' rna_data <- read_rna_encode_caltech(rnaseq_file)
#'
#' # Extract the gene name and gene expression in fpkm
#' gene_name <- rna_data$gene_name
#' gene_fpkm <- rna_data$gene_fpkm
#'
#' @export
read_rna_encode_caltech <- function(file, chr_discarded = NULL,
                                    is_GRanges = TRUE){
    message("Reading file ", file, " ...")
    data_raw <- scan(file = file,
             sep = "\t",
             what = list("character",  # Reference chromosome
                         integer(),    # Start position in chromosome
                         integer(),    # End position in chromosome
                         "character",  # Gene ENSEMBL id
                         numeric(),    # Expression level
                         "character",  # Strand : + or - or . for unknown
                         NULL,         # Source, e.g. HAVANA
                         NULL,         # Type of feature, e.g. gene
                         NULL,         # No information
                         "character"   # Metadata
             ))


    # Store only required fields
    rna_data <- data.table::data.table(chr = data_raw[[1]],
                                       start = data_raw[[2]],
                                       end = data_raw[[3]],
                                       strand = data_raw[[6]],
                                       ensembl_id = data_raw[[4]],
                                       gene_expr = data_raw[[5]])


    # Extract FPKM and gene name from each gene ------------------
    gene_info <- data_raw[[10]]
    fpkm      <- vector(mode = "numeric")
    gene_name <- vector(mode = "character")

    message("Extracting FPKM and gene names ...")
    for (i in 1:length(gene_info)){
        fpkm[i]      <- .extract_fpkm(gene_info[i])
        gene_name[i] <- .extract_gene_name(gene_info[i])
    }

    # Store extracted data to data.frame
    rna_data <- data.table::data.table(rna_data,
                                       fpkm = fpkm,
                                       gene_name = gene_name)
    rm(data_raw)


    # Remove selected chromosomes  -------------------------------
    rna_data <- .discard_chr(x = rna_data, chr_discarded = chr_discarded)


    # Sorting data -----------------------------------------------
    # With order priority: 1. chr, 2. start, 3. strand
    message("Sorting RNA-Seq data ...")
    rna_data <- rna_data[order(rna_data$chr, rna_data$start, rna_data$strand)]


    if (is_GRanges){
        # Create a GRanges object -----------------------------------
        message("Creating GRanges object ...")
        rna_data <- GenomicRanges::GRanges(seqnames = rna_data$chr,
                           strand = rna_data$strand,
                           ranges = IRanges::IRanges(start = rna_data$start,
                                                     end = rna_data$end),
                           ensembl_id = rna_data$ensembl_id,
                           gene_name  = rna_data$gene_name,
                           gene_fpkm  = rna_data$fpkm)
    }
    message("Finished reading RNA-Seq file!\n")
    return(rna_data)
}
