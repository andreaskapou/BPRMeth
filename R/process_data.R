#' @title Read methylation file
#'
#' @description \code{read_met} reads a file containing methylation data using
#'   the \code{\link[data.table]{fread}} function. Since there are different
#'   technologies, e.g. array, bulk BS-Seq, scBS-Seq, and still there is no
#'   standard file format, different options are available, check the Important
#'   section below on the file format for each type you choose. If a file format
#'   is not availabe, you need to read the file and create a
#'   \code{\link[GenomicRanges]{GRanges}} object, where the data are ordered by
#'   chromosome and genomic location.
#'
#' @param file File name.
#' @param type Type of technology as character. Either "bulk_seq", "sc_seq" or
#'   "array". Check the Important section below for more details.
#' @param strand_info Logical, whether or not the file contains strand
#'   information.
#' @param chr_discarded Optional vector with chromosomes to be discarded.
#' @param min_bulk_cov Minimum number of reads mapping to each CpG site. Used
#'   only for "bulk_seq" and CpGs with less reads will be discarded as noise.
#' @param max_bulk_cov Maximum number of reads mapping to each CpG site. Used
#'   only for "bulk_seq" and CpGs with less reads will be discarded as noise.
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object.
#'
#'   The GRanges object contains one or two additional metadata columns:
#'   \itemize{ \item{ \code{met}: Methylation level. \itemize{ \item{For "array"
#'   this is the Beta or M-values} \item{For "sc_seq" this is either 0 or 1
#'   (unmethylated or methylated)} \item{For "bulk_seq" this contains the number
#'   of methylated reads for each CpG.} } }  \item{ \code{total}: Total number
#'   of reads for each CpG. Present only for "bulk_seq" type. } } These columns
#'   can be accessed as follows: \code{granges_obj$met}
#'
#' @section Important: Depending on technology \code{type} we assume different
#'   file formats. \itemize{ \item{ \code{"array"} File format: "chromosome",
#'   "start", "strand" (optional), "met" .} \item{ \code{"sc_seq"} File format:
#'   "chromosome", "start", "strand" (optional), "met" . Where "met" should
#'   contain only 1s or 0s. } \item{ \code{"bulk_seq"} File format: "chromosome",
#'   "start", "strand" (optional), "met", "total". } } Columns should be in tab
#'   delimited format.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{preprocess_bs_seq}}
#'
#' @examples
#' \dontrun{
#' # Download the files and change the working directory to that location
#' file <- "name_of_file"
#' met_dt <- read_met(file)
#'
#' # Extract methylation level
#' met <- met_dt$met
#' }
#'
#' @export
read_met <- function(file, type = "sc_seq", strand_info = FALSE,
                     chr_discarded = NULL, min_bulk_cov = 4,
                     max_bulk_cov = 1000){

    chr <- NULL # So RMD CHECK does not complain
    message("Reading file ", file, " ...")
    met_dt <- data.table::fread(input = file, sep = "\t",
                                stringsAsFactors = TRUE, showProgress = FALSE)
    if (identical(type, "sc_seq") || identical(type, "array")) {
        if (strand_info) {
            met_dt <- met_dt %>%
                data.table::setnames(c("chr", "start", "strand", "met")) %>%
                subset(!(chr %in% chr_discarded)) %>%
                .[, c("end") := list(start)] %>%
                .[,c("chr", "start", "end", "strand", "met")] %>%
                data.table::setkey(chr, start, end) %>% GRanges()
        }else {
            met_dt <- met_dt %>%
                data.table::setnames(c("chr", "start", "met")) %>%
                subset(!(chr %in% chr_discarded)) %>%
                .[, c("end") := list(start)] %>%
                .[,c("chr", "start", "end", "met")] %>%
                data.table::setkey(chr, start, end) %>% GRanges()
        }
    }else if (identical(type, "bulk_seq")) {
        if (strand_info) {
            met_dt <- met_dt %>%
                data.table::setnames(c("chr", "start", "strand", "met",
                                       "total")) %>%
                subset(!(chr %in% chr_discarded)) %>%
                .[, c("end") := list(start)] %>%
                .[, c("chr", "start", "end", "strand", "met", "total")] %>%
                data.table::setkey(chr, start, end) %>% GRanges()
        }else {
            met_dt <- met_dt %>%
                data.table::setnames(c("chr", "start", "met", "total")) %>%
                subset(!(chr %in% chr_discarded)) %>%
                .[, c("end") := list(start)] %>%
                .[, c("chr", "start", "end", "met", "total")] %>%
                data.table::setkey(chr, start, end) %>% GRanges()
        }
        met_dt <- subset(met_dt, met_dt$total >= min_bulk_cov)
        met_dt <- subset(met_dt, met_dt$total <= max_bulk_cov)
    }
    return(met_dt)
}


#' @title Read annotation file
#'
#' @description \code{read_anno} reads a file containing annotation data, which
#'   will be used to select genomic regions of interest for downstream analysis.
#'   The annotation file format is the following: "chromosome", "start", "end",
#'   "strand", "id", "name" (optional). Read the Important section below for
#'   more details.
#'
#' @param file File name.
#' @param chrom_size_file Optional file name to read genome chromosome sizes.
#' @param chr_discarded Optional vector with chromosomes to be discarded.
#' @param is_centre Logical, whether 'start' and 'end' locations are
#'   pre-centred. If TRUE, the mean of the locations will be chosen as centre.
#'   If FALSE, the 'start' will be chosen as the center; e.g. for genes the
#'   'start' denotes the TSS and we use this as centre to obtain K-bp upstream
#'   and downstream of TSS.
#' @param upstream Integer defining the length of bp upstream of 'centre' for
#'   creating the genomic region.
#' @param downstream Integer defining the length of bp downstream of 'centre'
#'   for creating the genomic region.
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object.
#'
#'   The GRanges object contains two or three additional metadata columns:
#'   \itemize{ \item{ \code{id}: Feature ID, e.g. Ensembl IDs for each gene.}
#'   \item{ \code{centre}: The central (middle) location of the genomic region.
#'   This is required when transforming 'methylation regions' in the (-1, 1)
#'   interval, the 'centre' would be at 0.} \item{\code{name} (Optional) the
#'   feature name.} } These columns can be accessed as follows:
#'   \code{granges_obj$id}
#'
#' @section Important: \itemize{ \item{When having strand information, and we
#'   are in the antisense strand ("-"), the 'start' is automatically switched
#'   with the 'end' location.} \item{Columns should be in tab delimited format.}
#'   \item{ The "name" feature is optional.} \item{When there is no "strand"
#'   info, this can be set to "*".}}
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{read_met}}
#'
#' @examples
#' \dontrun{
#' # Download the files and change the working directory to that location
#' file <- "name_of_anno_file"
#' anno_dt <- read_anno(file)
#'
#' # Extract feature id
#' anno_ids <- anno_dt$id
#' }
#'
#' @export
read_anno <- function(file, chrom_size_file = NULL, chr_discarded = NULL,
                      is_centre = FALSE, upstream = -5000, downstream = 5000){

    chr <- NULL # So RMD CHECK does not complain
    message("Reading file ", file, " ...")
    anno_dt <- data.table::fread(input = file, sep = "\t",
                                 stringsAsFactors = TRUE, showProgress = FALSE)
    # Read the chromosome size file, if it is supplied
    if (!is.null(chrom_size_file)) {
        chrom_size <- read_chrom_size(file = chrom_size_file)
    }else {chrom_size <- NULL }
    if (NCOL(anno_dt) == 5) {
        anno_dt <- anno_dt %>%
            data.table::setnames(c("chr", "start", "end", "strand", "id"))
    }else{
        anno_dt <- anno_dt %>%
            data.table::setnames(c("chr","start","end","strand","id","name"))
    }
    anno_dt <- anno_dt %>% subset(!(chr %in% chr_discarded)) %>%
        data.table::setkey(chr, start, end) %>% GRanges()
    # Create genomic region
    anno_dt <- create_genomic_region(anno = anno_dt, chrom_size = chrom_size,
          is_centre = is_centre, upstream = upstream, downstream = downstream)
    return(anno_dt)
}


#' @name create_met_region
#' @rdname create_met_region
#' @aliases met_region, methylation_region, create_region
#'
#' @title Create methylation regions
#'
#' @description \code{create_met_region} creates methylation regions using as
#'   input methylation and annotation data with genomic regions of interest.
#'
#' @param met_dt A \code{\link[GenomicRanges]{GRanges}} object with methylation
#'   data, whose format should be similar to \code{\link{read_met}} function.
#' @param anno_dt A \code{\link[GenomicRanges]{GRanges}} object with annotation
#'   data, whose format should be similar to  \code{\link{read_anno}}.
#' @param cov Integer defining the minimum coverage of CpGs that each region
#'   must contain.
#' @param sd_thresh Optional numeric defining the minimum standard deviation of
#'   the methylation change in a region. This is used to filter regions with no
#'   methylation variability.
#' @param ignore_strand Logical, whether or not to ignore strand information.
#' @param filter_empty_region Logical, whether to discard genomic regions that
#'   have no CpG coverage or do not pass filtering options.
#' @param fmin Minimum range value for location scaling. Under this version, it
#'   should be left to its default value -1.
#' @param fmax Maximum range value for location scaling. Under this version, it
#'   should be left to its default value 1.
#'
#' @return A \code{list} object containing the two elements: \itemize{ \item{
#'   \code{met}: A list containing methylation regions, where each entry in the
#'   list is an \eqn{L_{i} X D} dimensional matrix, where \eqn{L_{i}} denotes
#'   the number of CpGs found in region \code{i}. The columns contain the
#'   following information: \enumerate{ \item{ 1st column: Contains the
#'   locations of CpGs relative to centre. Note that the actual locations are
#'   scaled to the (fmin, fmax) region. } \item{ 2nd column: If "bulk" data
#'   (i.e. binomial) it contains the total number of reads at each CpG location,
#'   otherwise the methylation level.} \item{ 3rd column: If "bulk" data, the
#'   methylated reads at each CpG location, otherwise this D = 2 and this column
#'   is absent.} } } \item{ \code{anno}: The annotation object.} } Note: The
#'   lengths of \code{met} and \code{anno} should match.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @examples
#' \dontrun{
#' # Download the files and change the working directory to that location
#' met_dt <- read_met("name_of_met_file")
#' anno_dt <- read_met("name_of_anno_file")
#'
#' obj <- create_met_region(met_dt, anno_dt)
#'
#' # Extract methylation regions
#' met <- obj$met
#' }
#'
#' @seealso \code{\link{read_met}}, \code{\link{read_anno}}
#'
#' @importFrom methods is
#' @export
create_met_region <- function(met_dt, anno_dt, cov = 5, sd_thresh = 1e-1,
                              ignore_strand = TRUE, filter_empty_region = TRUE,
                              fmin = -1, fmax = 1){
    message("Creating methylation regions ...")
    assertthat::assert_that(methods::is(met_dt, "GRanges"))
    assertthat::assert_that(methods::is(anno_dt, "GRanges"))
    # Find overlaps between met and anno
    overlaps <- GenomicRanges::findOverlaps(query = anno_dt, subject = met_dt,
                                            ignore.strand = ignore_strand)
    if (length(overlaps) < 2) { stop("Low overlap between anno and met data!")}

    # Convert data in vector format for efficiency
    query_hits <- S4Vectors::queryHits(overlaps)
    subj_hits  <- S4Vectors::subjectHits(overlaps)
    gen_ind    <- unique(query_hits) # Indices of genomic locations
    centre     <- anno_dt$centre     # Central locations
    id         <- anno_dt$id         # (Ensembl) IDs
    strand     <- as.character(GenomicRanges::strand(anno_dt))
    cpg_loc    <- GenomicRanges::ranges(met_dt)@start  # CpG locations
    met        <- met_dt$met         # Methylated read
    # Number of columns for each matrix
    if (is.null(met_dt$total)) { D <- 2 } else {total <- met_dt$total; D <- 3 }
    # Extract upstream and downstream lengths in bps
    width <- GenomicRanges::ranges(anno_dt)@width[1]
    if (identical(strand[1], "-")) {
        downstream <- abs(GenomicRanges::ranges(anno_dt)@start[1] - centre[1])
        upstream   <- downstream - width + 1
    }else{
        upstream   <- GenomicRanges::ranges(anno_dt)@start[1] - centre[1]
        downstream <- width + upstream - 1
    }

    # Initialize variables -----------------------------------------
    met_region <- list()
    # Initilize list to NA
    for (j in 1:length(id)) { met_region[[id[j]]] <- NA }
    LABEL     <- FALSE                     # Flag variable
    reg_count <- 0                         # Promoter counter
    cpg_ind   <- vector(mode = "integer")  # Vector of CpG indices
    cpg_ind   <- c(cpg_ind, subj_hits[1])  # Add the first subject hit

    total_iter <- NROW(query_hits)
    for (i in 2:total_iter) {
        # If query hits is the same as the previous one
        if (query_hits[i] == query_hits[i - 1]) {
            cpg_ind <- c(cpg_ind, subj_hits[i])  # Add subject hit
            # In case we have the last region
            if (i == total_iter) { reg_count <- reg_count + 1; LABEL <- TRUE }
        }else {reg_count <- reg_count + 1; LABEL <- TRUE } # Increase counter

        if (LABEL) {
            # Central location for region 'reg_count'
            idx <- id[gen_ind[reg_count]]
            # Only keep regions that have more than 'n' CpGs
            if (length(cpg_ind) > cov) {
                # If sd of the methylation level is above threshold
                if (D == 2) { obs_var <- stats::sd(met[cpg_ind])
                }else {obs_var <- stats::sd(met[cpg_ind]/total[cpg_ind]) }
                if (obs_var > sd_thresh) {
                    # Locations of CpGs in the genome
                    region <- cpg_loc[cpg_ind]
                    # Middle location for region 'reg_count'
                    middle <- centre[gen_ind[reg_count]]
                    # Extract strand information, i.e. direction
                    strand_direction <- strand[gen_ind[reg_count]]
                    # Shift CpG locations relative to TSS
                    center_data  <- .do_centre_loc(region = region,
                        centre = middle, strand_direction = strand_direction)
                    # In the "-" strand the order of the locations should change
                    Order <- base::order(center_data)
                    met_region[[idx]] <- matrix(0, length(cpg_ind), D)
                    # Store normalized locations of methylated CpGs
                    met_region[[idx]][, 1] <- round(.minmax_scaling(
                        data = center_data[Order], xmin = upstream,
                        xmax = downstream, fmin = fmin, fmax = fmax), 4)
                    # Store methylated reads in the corresponding locations
                    if (D == 2) { met_region[[idx]][, 2] <- met[cpg_ind][Order]
                    }else{
                        # Store total reads in the corresponding locations
                        met_region[[idx]][, 2] <- total[cpg_ind][Order]
                        # Store methylated reads in the corresponding locations
                        met_region[[idx]][, 3] <- met[cpg_ind][Order]
                    }
                }
            }
            LABEL   <- FALSE
            cpg_ind <- vector(mode = "integer")
            cpg_ind <- c(cpg_ind, subj_hits[i])
        }
    }
    # Keep only covered genomic regions
    if (filter_empty_region) {
        cov_ind <- which(!is.na(met_region))
        met_region <- met_region[cov_ind]
        anno_dt <- anno_dt[cov_ind, ]
    }
    return(list(met = met_region, anno = anno_dt))
}


#' @title Read expression data file
#'
#' @description \code{read_expr} reads a file containing expression data. It
#'   should contain two columns: "id", "expr" and columns should be in tab
#'   delimited format.
#'
#' @param file File name
#' @param log2_transf Logical, whether to log2 transform the expression data.
#'
#' @return A \code{\link[data.table]{data.table}} object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{read_met}}, \code{\link{read_anno}}
#'
#' @examples
#' \dontrun{
#' # Download the files and change the working directory to that location
#' file <- "name_of_expr_file"
#' anno_dt <- read_expr(file)
#'
#' # Extract feature id
#' expr_ids <- anno_dt$id
#' }
#'
#' @export
read_expr <- function(file, log2_transf = FALSE){
    expr_dt <- data.table::fread(input = file, sep = "\t",
                                 col.names = c("id", "expr"))
    if (log2_transf) { expr_dt$expr <- log2(expr_dt$expr + 0.1) }
    return(expr_dt)
}


#' @title Read genome chromosome sizes file
#'
#' @description \code{read_chrom_size} reads a file containing genome chromosome
#'   sizes using the \code{\link[data.table]{fread}} function.
#'
#' @param file File name
#'
#' @return A \code{\link[data.table]{data.table}} object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{read_met}}, \code{\link{read_anno}}
#'
#' @examples
#' chr_file <- system.file("extdata", "hg19.chr.sizes", package = "BPRMeth")
#' chr_data <- read_chrom_size(chr_file)
#'
#' # Extract the size of chr1
#' chr_data[1]
#'
#' @export
read_chrom_size <- function(file){
    chr_size <- data.table::fread(input = file, sep = "\t",
                                  col.names = c("chr", "size"))
    return(chr_size)
}


#' @title Create genomic regions from annotation data.
#'
#' @description \code{create_genomic_region} creates genomic region from
#'   annotation data, using the central point of the annotation features as
#'   ground truth labels we create genomic regions \code{N} bp upstream and
#'   \code{M} bp downstream of central location.
#'
#' @param anno A \code{\link[GenomicRanges]{GRanges}} object containing the
#'   annotation data, this normally would be the output from
#'   \code{\link{read_anno}} function.
#' @param chrom_size Object containing genome chromosome sizes, normally would
#'   be the output of \code{\link{read_chrom_size}} function.
#' @inheritParams read_anno
#' @return A \code{\link[GenomicRanges]{GRanges}} object containing the genomic
#'   regions.
#'
#'   The GRanges object contains two or three additional metadata column:
#'   \itemize{\item{\code{id}: Genomic region id.} \item{\code{centre}: Central
#'   location of each genomic region.} \item{\code{name}: (Optional) Genomic
#'   region name.} } This column can be accessed as follows:
#'   \code{granges_object$tss}
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{create_met_region}}, \code{\link{read_anno}}
#'
#' @examples
#' \dontrun{
#' # Download the files and change the working directory to that location
#' anno_file <- "name_of_anno_file"
#' annot_data <- read_anno(anno_file)
#' gen_region <- create_genomic_region(annot_data)
#'
#' # Extract ID
#' id <- gen_region$id
#' }
#'
#' @importFrom methods is
#' @export
create_genomic_region <- function(anno, chrom_size = NULL, is_centre = FALSE,
                                  upstream = -5000, downstream = 5000){
    assertthat::assert_that(methods::is(anno, "GRanges"))
    N <- NROW(anno)  # Number of features
    if (upstream > 0 ) { upstream <- -upstream }
    # Create empty vectors
    centre <- vector(mode = "integer", N)  # Centre location
    up     <- vector(mode = "integer", N)  # Start location in chromosome
    down   <- vector(mode = "integer", N)  # End location in chromosome
    chr    <- as.character(anno@seqnames)                   # Chrom info
    strand <- as.character(GenomicRanges::strand(anno))     # Strand info
    start  <- GenomicRanges::ranges(anno)@start             # Start info
    end    <- start + GenomicRanges::ranges(anno)@width - 1 # End info

    for (i in 1:N) {
        # Depending on the strand we change regions up or downstream of centre
        if (identical(strand[i], "-")) {
            centre[i] <- end[i]  # Set centre location
            # Set downstream bp promoter region
            up[i] <- max(0, end[i] - downstream)
            # Set upstream bp promoter region
            if (is.null(chrom_size)) {down[i] <- end[i] - upstream
            }else {down[i] <- min(chrom_size[chrom_size$chr == chr[i]]$size,
                                  end[i] - upstream) }
        }else {
            centre[i] <- start[i]  # Set centre location
            # Set upstream bp promoter region
            up[i] <- max(0, start[i] + upstream)
            # Set downstream bp promoter region
            if (is.null(chrom_size)) { down[i] <- start[i] + downstream
            }else {down[i] <- min(chrom_size[chrom_size$chr == chr[i]]$size,
                                  start[i] + downstream) }
        }
    }

    # Create GRanges object
    if (is.null(anno$name)) {
        genomic_region <- GenomicRanges::GRanges(seqnames = chr,
            ranges = IRanges::IRanges(start = up, end = down),
            strand = strand, id = anno$id, centre = centre)
    }else {
        genomic_region <- GenomicRanges::GRanges(seqnames = chr,
            ranges = IRanges::IRanges(start = up, end = down),
            strand = strand, id = anno$id, centre = centre, name = anno$name)
    }
    return(genomic_region)
}
