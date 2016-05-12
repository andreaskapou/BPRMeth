#' Create methylation regions for each gene promoter.
#'
#' \code{create_methyl_region} creates methylation regions using BS-Seq and
#'  pre-processed annotated gene promoter regions. BS-Seq data give information
#'  for the methylated CpGs individually, and promoter data are used to locate
#'  the TSS of each gene and its promoter region.
#'
#' @param bs_data \code{\link[GenomicRanges]{GRanges}} object containing the
#'  BS-Seq data.
#' @param prom_region \code{\link[GenomicRanges]{GRanges}} object containing
#'  promoter regions, i.e. N bp upstream and M bp downstream of TSS location.
#' @param cpg_density Optional integer defining the minimum number of CpGs that
#'  have to be in a methylated region. Regions with less than \code{n} CpGs
#'  are discarded.
#' @param sd_thresh Optional numeric defining the minimum standard deviation
#'  of the methylation change in a region.
#' @param ignore_strand Logical, whether or not to ignore strand information.
#' @param fmin Optional minimum range value for region location scaling.
#' @param fmax Optional maximum range value for region location scaling.
#'
#' @return A \code{meth_regions} object containing the following information:
#'  \itemize{
#'    \item{ \code{meth_data}: A list containing methylation data, where each
#'      each entry in the list consists of an L X 3 dimensional matrix, where:
#'      \enumerate{
#'        \item{ 1st column: Contains the locations of CpGs relative to TSS,
#'          where the range (min, max) of possible values is given, by the
#'          inputs fmin and fmax.
#'        }
#'        \item{ 2nd column: Contains the total reads of the CpG in the
#'          corresponding location.}
#'        \item{ 3rd column: Contains the methylated reads of the CpG in the
#'          corresponding location.}
#'      }
#'    }
#'    \item{ \code{prom_ind}: A vector storing the corresponding promoter
#'      indices.}
#'  }
#'
#' @seealso \code{\link{read_rna_encode_caltech}},
#'  \code{\link{read_bs_encode_haib}},
#'  \code{\link{create_prom_region}}
#'
#' @importFrom stats sd
#' @importFrom methods is
#' @export
create_methyl_region <- function(bs_data, prom_region, cpg_density = 1,
                                 sd_thresh = 0, ignore_strand = FALSE,
                                                  fmin = -1, fmax = 1){

  message("Creating methylation regions ...")
  assertthat::assert_that(methods::is(bs_data, "GRanges"))
  assertthat::assert_that(methods::is(prom_region, "GRanges"))

  # Find overlaps between promoter regions and BS-Seq data -------
  overlaps <- GenomicRanges::findOverlaps(query   = prom_region,
                                          subject = bs_data,
                                          ignore.strand = ignore_strand)
  if (length(overlaps) < 2){
    stop("Not enough matches between the BS-Seq and RNA-Seq data.")
  }


  # Convert data in vector format for faster lookup --------------
  query_hits <- S4Vectors::queryHits(overlaps)
  subj_hits  <- S4Vectors::subjectHits(overlaps)

  # Indices of promoter locations
  prom_loc   <- unique(query_hits)
  tss_loc    <- prom_region$tss        # TSS locations
  tss_strand <- as.character(GenomicRanges::strand(prom_region))
  cpg_loc    <- GenomicRanges::ranges(bs_data)@start  # CpG locations
  tot_reads  <- bs_data$total_reads    # Total reads
  meth_reads <- bs_data$meth_reads     # Methylated reads

  # Extract upstream and downstream lengths in bps
  width        <- GenomicRanges::ranges(prom_region)@width[1]
  if (identical(tss_strand[1], "+")){
    upstream   <- GenomicRanges::ranges(prom_region)@start[1] - tss_loc[1]
    downstream <- width + upstream - 1
  }else if (identical(tss_strand[1], "-")){
    downstream <- abs(GenomicRanges::ranges(prom_region)@start[1] - tss_loc[1])
    upstream   <- downstream - width + 1
  }


  # Initialize variables -----------------------------------------
  n            <- 1                         # Data points counter
  LABEL        <- FALSE                     # Flag variable
  meth_data    <- list()                    # List where data will be stored
  prom_counter <- 0                         # Promoter counter
  prom_ind     <- vector(mode = "integer")  # Vector of promoter indices
  cpg_ind      <- vector(mode = "integer")  # Vector of CpG indices
  cpg_ind      <- c(cpg_ind, subj_hits[1])  # Add the first subject hit

  for (i in 2:NROW(query_hits)){
    # If query hits is the same as the previous one
    if (query_hits[i] == query_hits[i - 1]){
      cpg_ind <- c(cpg_ind, subj_hits[i])  # Add subject hit
    }else{
      prom_counter <- prom_counter + 1  # Increase promoter counter
      LABEL <- TRUE
    }

    if(LABEL){
      # Only keep regions that have more than 'n' CpGs
      if (length(cpg_ind) > cpg_density){

        # If standard deviation of the methylation level is above threshold
        if (sd(meth_reads[cpg_ind] / tot_reads[cpg_ind]) > sd_thresh){

          # Promoter indices
          prom_ind <- c(prom_ind, prom_loc[prom_counter])

          # Locations of CpGs in the genome
          region <- cpg_loc[cpg_ind]

          # TSS location for promoter 'promCount'
          tss <- tss_loc[prom_loc[prom_counter]]

          # Extract strand information, i.e. direction
          strand_direction <- tss_strand[prom_loc[prom_counter]]

          # Shift CpG locations relative to TSS
          center_data  <- center_loc(region = region,
                                     tss = tss,
                                     strand_direction = strand_direction)

          # In the "-" strand the order of the locations should change
          Order <- base::order(center_data)

          meth_data[[n]] <- matrix(0, nrow = length(cpg_ind), ncol = 3)

          # Store normalized locations of methylated CpGs in (fmin,fmax) region
          meth_data[[n]][ ,1] <- .minmax_scaling(data = center_data[Order],
                                                 xmin = upstream,
                                                 xmax = downstream,
                                                 fmin = fmin,
                                                 fmax = fmax)

          # Store total reads in the corresponding locations
          meth_data[[n]][ ,2] <- tot_reads[cpg_ind][Order]

          # Store methylated reads in the corresponding locations
          meth_data[[n]][ ,3] <- meth_reads[cpg_ind][Order]

          # Increase data points counter
          n <- n + 1
        }
      }
      LABEL   <- FALSE
      cpg_ind <- vector(mode="integer")
      cpg_ind <- c(cpg_ind, subj_hits[i])
    }
  }
  methyl_region <- structure(list(meth_data = meth_data,
                                  prom_ind = prom_ind),
                             class = "methyl_region")
  message("Done!\n")
  return(methyl_region)
}


#' Center CpG locations relative to TSS
#'
#' \code{center_loc} centera CpG locations relative to TSS
#'
#' @param region CpG locations
#' @param tss TSS location
#' @param strand_direction Strand direction
#'
#' @return Centered location data relative to TSS
#'
center_loc <- function(region, tss, strand_direction){
  assertthat::assert_that(is.character(strand_direction))

  center <- region - tss
  if (identical(strand_direction, "-")){
    center  <- (-center)  # If '-' strand, swap CpG locations
  }
  return(center)
}
