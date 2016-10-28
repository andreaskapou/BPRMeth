#' Create promoter regions from gene annotation data.
#'
#' \code{create_prom_region} creates promoter region from gene annotation data.
#' Using the TSS of gene annotation data as ground truth labels we create
#' promoter regions \code{N} bp upstream and \code{M} bp downstream of TSS.
#'
#' @param annot_data A \code{\link[GenomicRanges]{GRanges}} object containing
#'   the gene annotation data. This for example can be RNA-Seq data output from
#'   \code{\link{read_rna_encode_caltech}}.
#' @param chrom_size Optional \code{\link[data.table]{data.table}} containing
#'   chromosome sizes, e.g. using the \code{\link{read_chrom_size}} function.
#' @param upstream Integer defining the length of bp upstream of TSS.
#' @param downstream Integer defining the length of bp downstream of TSS.
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object containing the promoter
#'   regions data.
#'
#'   The GRanges object contains one additional metadata column: \itemize{ \item
#'   \code{tss}: TSS of each gene promoter.} This column can be accessed as
#'   follows: \code{granges_object$tss}
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{create_methyl_region}}, \code{\link{read_chrom_size}},
#'   \code{\link{read_rna_encode_caltech}}
#'
#' @examples
#' # Obtain the path to the file and then read it
#' rnaseq_file <- system.file("extdata", "rnaseq.bed", package = "BPRMeth")
#' annot_data <- read_rna_encode_caltech(rnaseq_file)
#' prom_region <- create_prom_region(annot_data)
#'
#' # Extract the TSS
#' tss <- prom_region$tss
#'
#' @importFrom methods is
#' @export
create_prom_region <- function(annot_data, chrom_size = NULL, upstream = -7000,
                               downstream = 7000){

    message("Creating promoter regions ...")
    assertthat::assert_that(methods::is(annot_data, "GRanges"))
    N <- NROW(annot_data)  # Number of genes
    if (upstream > 0 ){
        upstream <- -upstream
    }

    # Create empty vectors
    tss       <- vector(mode = "integer", N)  # TSS location
    up_prom   <- vector(mode = "integer", N)  # Start location in chromosome
    down_prom <- vector(mode = "integer", N)  # End location in chromosome

    # Extract chromosome information
    annot_chr    <- as.character(annot_data@seqnames)
    # Extract strand information
    annot_strand <- as.character(GenomicRanges::strand(annot_data))
    # Extract start information
    annot_start  <- GenomicRanges::ranges(annot_data)@start
    # Extract end information
    annot_end    <- annot_start + GenomicRanges::ranges(annot_data)@width - 1

    for (i in 1:N){
        # Depending on the strand we change regions up or downstream of TSS
        if (identical(annot_strand[i], "+")){
            # Set TSS location
            tss[i] <- annot_start[i]
            # Set upstream bp promoter region
            up_prom[i] <- max(0, annot_start[i] + upstream)
            # Set downstream bp promoter region
            if (is.null(chrom_size)){
                down_prom[i] <- annot_start[i] + downstream
            }else{
                down_prom[i] <- min(chrom_size[chrom_size$chr ==
                                                   annot_chr[i]]$size,
                                    annot_start[i] + downstream)
            }
        }else if (identical(annot_strand[i], "-")){
            # Set TSS location
            tss[i] <- annot_end[i]
            # Set downstream bp promoter region
            up_prom[i] <- max(0, annot_end[i] - downstream)
            # Set upstream bp promoter region
            if (is.null(chrom_size)){
                down_prom[i] <- annot_end[i] - upstream
            }else{
                down_prom[i] <- min(chrom_size[chrom_size$chr ==
                                                   annot_chr[i]]$size,
                                    annot_end[i] - upstream)
            }
        }
    }

    # Create GRanges object
    message("Creating GRanges object for promoter regions ...")
    prom_region <- GenomicRanges::GRanges(seqnames = annot_chr,
                              ranges   = IRanges::IRanges(start = up_prom,
                                                          end = down_prom),
                              strand   = annot_strand,
                              tss      = tss)
    message("Finished creating promoter regions!\n")
    return(prom_region)
}
