#' Keep CGIs that intersect with promoter regions
#'
#' \code{intersect_cgi_tss} keeps only CpG islands that overlap with promoter
#'  regions \code{N} bp upstream and \code{M} bp downstream of TSS.
#'
#' @param cgi_file The name of the CGI '.bed' formatted file to read
#'  data values from.
#' @param rna_files The name of the RNA-Seq '.bed' formatted file to
#'  read data values from.
#' @param chrom_size_file Optional name of the file containing genome
#'  chromosome sizes
#' @param chr_discarded A vector with chromosome names to be discarded.
#' @param upstream Integer defining the length of bp upstream of TSS.
#' @param downstream Integer defining the length of bp downstream of TSS.
#' @param unique_tss Logical, indicating if CGIs should match only to unique
#'  TSS, default false, so as to consider strand direction.
#'
#' @return The remaining CGIs which intersect with promoter regions stored in a
#'  \code{\link[GenomicRanges]{GRanges}} object.
#'
#' @seealso \code{\link{create_prom_region}}
#'
#' @export
intersect_cgi_tss <- function(cgi_file, rna_files, chrom_size_file = NULL,
                              chr_discarded = NULL, upstream = -100,
                              downstream = 100, unique_tss = FALSE){

  cgi_data <- read_encode_cgi(file = cgi_file,
                              is_GRanges = TRUE)

  # Read RNA-Seq BED file
  rna_data <- read_rna_encode_caltech(file          = rna_files,
                                      chr_discarded = chr_discarded,
                                      is_GRanges    = TRUE)

  # Read the chromosome size file, if it is supplied
  if (!is.null(chrom_size_file)){
    chrom_size <- read_chrom_size(file = chrom_size_file)
  }else{
    chrom_size <- NULL
  }

  # Create promoter regions
  prom_region <- create_prom_region(annot_data = rna_data,
                                    chrom_size = chrom_size,
                                    upstream   = upstream,
                                    downstream = downstream)


  # Find overlaps between promoter regions and CGI data
  overlaps <- GenomicRanges::findOverlaps(query   = cgi_data,
                                          subject = prom_region,
                                          ignore.strand = TRUE)

  if (! unique_tss){
  # Get only the subset of overlapping CpG sites
  cgi_data <- cgi_data[S4Vectors::queryHits(overlaps)]
  # Add the corresponding ensembl ids
  cgi_data$ensembl_id <- rna_data[S4Vectors::subjectHits(overlaps)]$ensembl_id
  }else{
    # Find overlaps between promoter regions and CGI data
    overlaps <- GenomicRanges::findOverlaps(query   = cgi_data,
                                            subject = prom_region,
                                            select = "first",
                                            ignore.strand = TRUE)

    # Get only the subset of overlapping CpG islands
    keep_non_na <- which(!is.na(overlaps))
    cgi_data <- cgi_data[keep_non_na]
    # Add the corresponding ensembl ids
    cgi_data$ensembl_id <- rna_data[overlaps[keep_non_na]]$ensembl_id
  }

  return(cgi_data)
}
