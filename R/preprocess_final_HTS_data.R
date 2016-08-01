#' Pre-process final HTS data for downstream analysis
#'
#' \code{preprocess_final_HTS_data} performs a final filtering and preprocessing
#' on the data for use in downstream analysis. These include, removing noisy
#' gene expression data, removing or not un-expressed genes and log2-transorming
#' of the FPKM values.
#'
#' @param methyl_region Methylation region data, which are the output of the
#'   "\code{create_methyl_region}" function.
#' @param rna_data A \code{\link[GenomicRanges]{GRanges}} object containing
#'   corresponding RNA-Seq data for each entry of the \code{methyl_region} list.
#'   This is the output of the "\code{read_rna_encode_caltech} function.
#' @param prom_reg A \code{\link[GenomicRanges]{GRanges}} object containing
#'   corresponding annotated promoter regions for each entry of the
#'   \code{methyl_region} list.
#' @inheritParams process_haib_caltech_wrap
#'
#' @return  An object which contains following information: \itemize{ \item
#'   \code{methyl_region}: The subset of promoter methylation region data after
#'   the filtering process. \item \code{gex}: A vectoring storing only the
#'   corresponding gene expression values for each promoter region. \item
#'   \code{rna_data}: The corresponding gene expression data stored as a GRanges
#'   object.}
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{read_rna_encode_caltech}}
#'   \code{\link{process_haib_caltech_wrap}}
#'
#' @export
preprocess_final_HTS_data <- function(methyl_region, prom_reg, rna_data,
                                      gene_log2_transf = TRUE,
                                      gene_outl_thresh = TRUE,
                                      gex_outlier = 300){

  # Extract only the gene expression data in FPKM
  gex <- as.numeric(rna_data$gene_fpkm)

  # Option to discard possible outliers / noisy data
  if (gene_outl_thresh){
    ind <- which(gex > gex_outlier)
    if (length(ind) > 0){
      gex <- gex[-ind]
      methyl_region <- methyl_region[-ind]
      prom_reg <- prom_reg[-ind]
      rna_data <- rna_data[-ind]
    }
  }

  # Option to log-transform gene expression data
  if (gene_log2_transf){
    gex <- gex + 0.1
    gex <- log2(gex)
  }

  return(list(methyl_region = methyl_region,
              gex = gex,
              prom_reg = prom_reg,
              rna_data = rna_data))
}
