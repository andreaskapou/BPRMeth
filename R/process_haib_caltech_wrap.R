#' Wrapper method for processing ENCODE HAIB and Caltech HTS data
#'
#' \code{process_haib_caltech_wrap} is a wrapper method for processing HTS data
#' and returning the methylation promoter regions and the corresponding gene
#' expression data for those promoter regions. Note that the format of BS-Seq
#' data should be in the Encode Haib bed format and for the RNA-Seq data in
#' Encode Caltech bed format.
#'
#' @param bs_files Filename (or vector of filenames if there are replicates) of
#'   the BS-Seq '.bed' formatted data to read values from.
#' @param rna_files Filename of the RNA-Seq '.bed' formatted data to read values
#'   from. Currently, this version does not support pooling RNA-Seq replicates.
#' @param chrom_size_file Optional filename containing genome chromosome sizes.
#' @param chr_discarded A vector with chromosome names to be discarded.
#' @param upstream Integer defining the length of bp upstream of TSS for
#'   creating the promoter region.
#' @param downstream Integer defining the length of bp downstream of TSS for
#'   creating the promoter region.
#' @param min_bs_cov The minimum number of reads mapping to each CpG site. CpGs
#'   with less reads will be considered as noise and will be discarded.
#' @param max_bs_cov The maximum number of reads mapping to each CpG site. CpGs
#'   with more reads will be considered as noise and will be discarded.
#' @inheritParams create_methyl_region
#'
#' @return A \code{processHTS} object which contains following information:
#'   \itemize{ \item{ \code{methyl_region}: A list containing methylation data,
#'   where each entry in the list is an \eqn{L_{i} X 3} dimensional matrix,
#'   where \eqn{L_{i}} denotes the number of CpGs found in region \code{i}. The
#'   columns contain the following information: \enumerate{ \item{ 1st column:
#'   Contains the locations of CpGs relative to TSS. Note that the actual
#'   locations are scaled to the (fmin, fmax) region. } \item{ 2nd column:
#'   Contains the total reads of each CpG in the corresponding location.} \item{
#'   3rd column: Contains the methylated reads each CpG in the corresponding
#'   location.} } } \item{ \code{prom_region}: A
#'   \code{\link[GenomicRanges]{GRanges}} object containing corresponding
#'   annotated promoter regions for each entry of the \code{methyl_region} list.
#'   The GRanges object has one additional metadata column named \code{tss},
#'   which stores the TSS of each promoter. } \item{ \code{rna_data}: A
#'   \code{\link[GenomicRanges]{GRanges}} object containing the corresponding
#'   RNA-Seq data for each entry of the \code{methyl_region} list. The GRanges
#'   object has three additional metadata columns which are explained in
#'   \code{\link{read_rna_encode_caltech}}} \item{ \code{upstream}: Integer
#'   defining the length of bp upstream of TSS.} \item{ \code{downstream}:
#'   Integer defining the length of bp downstream of TSS.} \item{
#'   \code{cpg_density}: Integer defining the minimum number of CpGs that have
#'   to be in a methylated region. Regions with less than \code{n} CpGs are
#'   discarded.} \item{ \code{sd_thresh}: Numeric defining the minimum standard
#'   deviation of the methylation change in a region. This is used to filter
#'   regions with no methylation change.} \item{ \code{fmin}: Minimum range
#'   value for region location scaling.} \item{ \code{fmax}: Maximum range value
#'   for region location scaling.} }
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @examples
#' # Get the location of the files
#' rrbs_file <- system.file("extdata", "rrbsH1hESC.bed", package = "processHTS")
#' rnaseq_file <- system.file("extdata", "rnaseqH1hESC.bed", package = "processHTS")
#' data <- process_haib_caltech_wrap(rrbs_file, rnaseq_file)
#'
#' @export
process_haib_caltech_wrap <- function(bs_files, rna_files,
                                      chrom_size_file = NULL,
                                      chr_discarded = NULL, upstream = -100,
                                      downstream = 100, min_bs_cov = 0,
                                      max_bs_cov = 1000, cpg_density = 1,
                                      sd_thresh = 0, ignore_strand = FALSE,
                                      fmin = -1, fmax = 1){


  # Process BS-Seq file and return data in the required format
  bs_data <- preprocess_bs_seq(files         = bs_files,
                               file_format   = "encode_rrbs",
                               chr_discarded = chr_discarded,
                               min_bs_cov    = min_bs_cov,
                               max_bs_cov    = max_bs_cov)

  # Read the chromosome size file, if it is supplied
  if (!is.null(chrom_size_file)){
    chrom_size <- read_chrom_size(file = chrom_size_file)
  }else{
    chrom_size <- NULL
  }

  # Read RNA-Seq BED file
  rna_data <- read_rna_encode_caltech(file          = rna_files,
                                      chr_discarded = chr_discarded,
                                      is_GRanges    = TRUE)

  # Create promoter regions
  prom_reg <- create_prom_region(annot_data = rna_data,
                                 chrom_size = chrom_size,
                                 upstream   = upstream,
                                 downstream = downstream)

  # Create methylation regions data
  methyl_reg <- create_methyl_region(bs_data       = bs_data,
                                     prom_region   = prom_reg,
                                     cpg_density   = cpg_density,
                                     sd_thresh     = sd_thresh,
                                     ignore_strand = ignore_strand,
                                     fmin          = fmin,
                                     fmax          = fmax)

  # Keep only the corresponding gene expression data
  rna_data <- rna_data[methyl_reg$prom_ind]
  # Keep only the corresponding gene annotation data
  prom_reg <- prom_reg[methyl_reg$prom_ind]

  # Create object
  obj <- structure(list(methyl_region = methyl_reg$meth_data,
                        prom_region   = prom_reg,
                        rna_data      = rna_data,
                        upstream      = upstream,
                        downstream    = downstream,
                        cpg_density   = cpg_density,
                        sd_thresh     = sd_thresh,
                        fmin          = fmin,
                        fmax          = fmax),
                   class = "processHTS")
  return(obj)
}
