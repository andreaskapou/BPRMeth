#' @title (DEPRECATED) Wrapper for processing ENCODE HAIB and Caltech HTS
#'
#' @description (DEPRECATED) \code{process_haib_caltech_wrap} is a wrapper
#'   method for processing HTS data and returning the methylation promoter
#'   regions and the corresponding gene expression data for those promoter
#'   regions. Note that the format of BS-Seq data should be in the Encode Haib
#'   bed format and for the RNA-Seq data in Encode Caltech bed format.
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
#' @param cpg_density Optional integer defining the minimum number of CpGs that
#'   have to be in a methylated region. Regions with less than \code{n} CpGs are
#'   discarded.
#' @param sd_thresh Optional numeric defining the minimum standard deviation of
#'   the methylation change in a region. This is used to filter regions with no
#'   methylation change.
#' @param ignore_strand Logical, whether or not to ignore strand information.
#' @param gene_log2_transf Logical, whether or not to log2 transform the gene
#'   expression data.
#' @param gene_outl_thresh Logical, whehter or not to remove outlier gene
#'   expression data.
#' @param gex_outlier Numeric, denoting the threshold above of which the gene
#'   expression data (before the log2 transformation) are considered as noise.
#' @param fmin Optional minimum range value for region location scaling. Under
#'   this version, this parameter should be left to its default value.
#' @param fmax Optional maximum range value for region location scaling. Under
#'   this version, this parameter should be left to its default value.
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
#'   location.} } } \item{\code{gex}: A vector containing the corresponding gene
#'   expression levels for each entry of the \code{methyl_region} list.} \item{
#'   \code{prom_region}: A \code{GRanges} object
#'   containing corresponding annotated promoter regions for each entry of the
#'   \code{methyl_region} list. The GRanges object has one additional metadata
#'   column named \code{tss}, which stores the TSS of each promoter. } \item{
#'   \code{rna_data}: A \code{GRanges} object containing
#'   the corresponding RNA-Seq data for each entry of the \code{methyl_region}
#'   list. The GRanges object has three additional metadata columns which are
#'   explained in \code{\link{read_rna_encode_caltech}}} \item{ \code{upstream}:
#'   Integer defining the length of bp upstream of TSS.} \item{
#'   \code{downstream}: Integer defining the length of bp downstream of TSS.}
#'   \item{ \code{cpg_density}: Integer defining the minimum number of CpGs that
#'   have to be in a methylated region. Regions with less than \code{n} CpGs are
#'   discarded.} \item{ \code{sd_thresh}: Numeric defining the minimum standard
#'   deviation of the methylation change in a region. This is used to filter
#'   regions with no methylation change.} \item{ \code{fmin}: Minimum range
#'   value for region location scaling.} \item{ \code{fmax}: Maximum range value
#'   for region location scaling.} }
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @examples
#' # Obtain the path to the files
#' rrbs_file <- system.file("extdata", "rrbs.bed", package = "BPRMeth")
#' rnaseq_file <- system.file("extdata", "rnaseq.bed", package = "BPRMeth")
#' proc_data <- process_haib_caltech_wrap(rrbs_file, rnaseq_file)
#'
#' @export
process_haib_caltech_wrap <- function(bs_files, rna_files,
                                      chrom_size_file = NULL,
                                      chr_discarded = NULL, upstream = -7000,
                                      downstream = 7000, min_bs_cov = 4,
                                      max_bs_cov = 1000, cpg_density = 10,
                                      sd_thresh = 10e-02, ignore_strand = TRUE,
                                      gene_log2_transf = TRUE,
                                      gene_outl_thresh = TRUE,
                                      gex_outlier = 300,
                                      fmin = -1, fmax = 1){
    # Process BS-Seq file and return data in the required format
    bs_data <- preprocess_bs_seq(files = bs_files, file_format = "encode_rrbs",
                        chr_discarded = chr_discarded, min_bs_cov = min_bs_cov,
                        max_bs_cov = max_bs_cov)
    # Read the chromosome size file, if it is supplied
    if (!is.null(chrom_size_file)) {
        chrom_size <- read_chrom_size(file = chrom_size_file)
    }else {chrom_size <- NULL }
    # Read RNA-Seq BED file
    rna_data <- read_rna_encode_caltech(file = rna_files,
                chr_discarded = chr_discarded, is_GRanges = TRUE)
    # Create promoter regions
    prom_reg <- create_anno_region(anno = rna_data, chrom_size = chrom_size,
                                   upstream = upstream, downstream = downstream)
    # Create methylation regions data
    methyl_reg <- create_region_object(met_dt = bs_data, anno_dt = prom_reg,
        cov = cpg_density, sd_thresh = sd_thresh, ignore_strand = ignore_strand,
        filter_empty_region = FALSE, fmin = fmin, fmax = fmax)$met

    # Keep only covered genomic regions
    cov_ind <- which(!is.na(methyl_reg))
    methyl_reg <- methyl_reg[cov_ind]
    prom_reg <- prom_reg[cov_ind, ]
    rna_data <- rna_data[cov_ind, ]
    proc_data <- preprocess_final_HTS_data(methyl_region = methyl_reg,
                                           prom_reg = prom_reg,
                                           rna_data = rna_data,
                                           gene_log2_transf = gene_log2_transf,
                                           gene_outl_thresh = gene_outl_thresh,
                                           gex_outlier = gex_outlier)
    # Create object
    obj <- structure(list(methyl_region = proc_data$methyl_region,
                          gex           = proc_data$gex,
                          prom_region   = proc_data$prom_reg,
                          rna_data      = proc_data$rna_data,
                          upstream      = upstream,
                          downstream    = downstream,
                          cpg_density   = cpg_density,
                          sd_thresh     = sd_thresh,
                          fmin          = fmin,
                          fmax          = fmax),
                     class = "processHTS")
    return(obj)
}


#' @title (DEPRECATED) Pre-process BS-Seq data in any given format
#'
#' @description (DEPRECATED) \code{preprocess_bs_seq} is a general function for
#'   reading and preprocessing BS-Seq data. If a vector of files is given, these
#'   are considered as replicates and are pooled together. Finally, noisy reads
#'   are discarded.
#'
#' @param files A vector of filenames containing replicate experiments. This can
#'   also be just a single replicate.
#' @param file_format A string denoting the file format that the BS-Seq data are
#'   stored. Current version allows "\code{encode_rrbs}" or "\code{bismark_cov}"
#'   formats.
#' @param chr_discarded A vector with chromosome names to be discarded.
#' @param min_bs_cov The minimum number of reads mapping to each CpG site. CpGs
#'   with less reads will be considered as noise and will be discarded.
#' @param max_bs_cov The maximum number of reads mapping to each CpG site. CpGs
#'   with more reads will be considered as noise and will be discarded.
#'
#' @return A \code{GRanges} object. The GRanges object contains two additional
#'   metadata columns: \itemize{ \item \code{total_reads}: total reads mapped to
#'   each genomic location. \item \code{meth_reads}: methylated reads mapped to
#'   each genomic location. } These columns can be accessed as follows:
#'   \code{granges_object$total_reads}
#'
#' @section Additional Info: Information about the file formats can be found in
#'   the following links:
#'
#'   Encode RRBS format:
#'   \url{http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgTables?db=
#'   hg19&hgta_group=regulation&hgta_track=wgEncodeHaibMethylRrbs&hgta_table=
#'   wgEncodeHaibMethylRrbsBcbreast0203015BiochainSitesRep2&hgta_doSchema=
#'   describe+table+schema}
#'
#'   Bismark Cov format: \url{http://rnbeads.mpi-inf.mpg.de/data/RnBeads.pdf}
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{read_bs_encode_haib}} \code{\link{pool_bs_seq_rep}}
#'
#' @examples
#' # Obtain the path to the files
#' bs_file <- system.file("extdata", "rrbs.bed", package = "BPRMeth")
#' bs_data <- preprocess_bs_seq(bs_file, file_format = "encode_rrbs")
#'
#' @export
preprocess_bs_seq <- function(files, file_format = "encode_rrbs",
                              chr_discarded = NULL, min_bs_cov = 4,
                              max_bs_cov = 1000){
    # If we have more than one replicates
    if (length(files) > 1) {
        bs_data <- pool_bs_seq_rep(files = files, file_format = file_format,
                                   chr_discarded = chr_discarded)
    }else{
        if (file_format == "encode_rrbs") {
            bs_data <- read_bs_encode_haib(file = files,
                            chr_discarded = chr_discarded, is_GRanges = TRUE)
        }else if (file_format == "bismark_cov") {
            bs_data <- .read_bs_bismark_cov(file = files,
                            chr_discarded = chr_discarded, is_GRanges = TRUE)
        }
        else {stop("Wrong file format.") }
    }
    bs_data <- .discard_bs_noise_reads(bs_data = bs_data,
                            min_bs_cov = min_bs_cov, max_bs_cov = max_bs_cov)
    return(bs_data)
}


#' @title (DEPRECATED) Pre-process final HTS data for downstream analysis
#'
#' @description (DEPRECATED) \code{preprocess_final_HTS_data} performs a final
#'   filtering and preprocessing on the data for use in downstream analysis.
#'   These include, removing noisy gene expression data, removing or not
#'   un-expressed genes and log2-transorming of the FPKM values.
#'
#' @param methyl_region Methylation region data, which are the output of the
#'   "\code{create_region_object}" function.
#' @param rna_data A \code{GRanges} object containing corresponding RNA-Seq data
#'   for each entry of the \code{methyl_region} list. This is the output of the
#'   "\code{read_rna_encode_caltech} function.
#' @param prom_reg A \code{GRanges} object containing corresponding annotated
#'   promoter regions for each entry of the \code{methyl_region} list.
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
#' @examples
#' # Obtain the path to the BS file and then read it
#' bs_file <- system.file("extdata", "rrbs.bed", package = "BPRMeth")
#' bs_data <- read_bs_encode_haib(bs_file)
#'
#' # Create promoter regions
#' rnaseq_file <- system.file("extdata", "rnaseq.bed", package = "BPRMeth")
#' annot_data <- read_rna_encode_caltech(rnaseq_file)
#' prom_region <- create_anno_region(annot_data)
#'
#' # Create methylation regions
#' methyl_reg <- create_region_object(bs_data, prom_region,
#'   filter_empty_region = FALSE)
#'
#' # Keep only covered genomic regions
#' cov_ind <- which(!is.na(methyl_reg))
#' methyl_reg <- methyl_reg[cov_ind]
#' prom_region <- prom_region[cov_ind, ]
#' annot_data <- annot_data[cov_ind, ]
#'
#' # Finally preprocess the HTS data
#' res <- preprocess_final_HTS_data(methyl_reg, prom_region, annot_data)
#'
#' @export
preprocess_final_HTS_data <- function(methyl_region, prom_reg, rna_data,
                                      gene_log2_transf = TRUE,
                                      gene_outl_thresh = TRUE,
                                      gex_outlier = 300){
    # Extract only the gene expression data in FPKM
    gex <- as.numeric(rna_data$gene_fpkm)
    # Option to discard possible outliers / noisy data
    if (gene_outl_thresh) {
        ind <- which(gex > gex_outlier)
        if (length(ind) > 0) {
            gex <- gex[-ind]
            methyl_region <- methyl_region[-ind]
            prom_reg <- prom_reg[-ind]
            rna_data <- rna_data[-ind]
        }
    }
    # Option to log-transform gene expression data
    if (gene_log2_transf) { gex <- log2(gex + 0.1) }
    return(list(methyl_region = methyl_region, gex = gex,
                prom_reg = prom_reg, rna_data = rna_data))
}




#' @title (DEPRECATED) Read and pool replicates from BS-Seq data
#'
#' @description (DEPRECATED) \code{pool_bs_seq_rep} reads and pools replicate
#'   methylation data from BS-Seq experiments that are either in Encode RRBS or
#'   Bismark Cov format. Read the Important section below on when to use this
#'   function.
#'
#' @inheritParams preprocess_bs_seq
#'
#' @return A \code{GRanges} object. The GRanges object contains two additional
#'   metadata columns: \itemize{ \item \code{total_reads}: total reads mapped to
#'   each genomic location. \item \code{meth_reads}: methylated reads mapped to
#'   each genomic location. } These columns can be accessed as follows:
#'   granges_object$total_reads
#'
#' @section Important: Unless you want to create a different workflow when
#'   processing the BS-Seq data, you should NOT call this function, since this
#'   is a helper function. Instead you should call the
#'   \code{\link{preprocess_bs_seq}} function.
#'
#'   Information about the file formats can be found in the following links:
#'
#'   Encode RRBS format:
#'   \url{http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgTables?db=hg19&hgta_group=
#'   regulation&hgta_track=wgEncodeHaibMethylRrbs&hgta_table=
#'   wgEncodeHaibMethylRrbsBcbreast0203015BiochainSitesRep2&hgta_doSchema=
#'   describe+table+schema}
#'
#'   Bismark Cov format: \url{http://rnbeads.mpi-inf.mpg.de/data/RnBeads.pdf}
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{read_bs_encode_haib}}, \code{\link{preprocess_bs_seq}}
#'
#' @examples
#' # Obtain the path to the file
#' bs_file1 <- system.file("extdata", "rrbs.bed", package = "BPRMeth")
#' bs_file2 <- system.file("extdata", "rrbs.bed", package = "BPRMeth")
#'
#' # Concatenate the files
#' bs_files <-  c(bs_file1, bs_file2)
#' # Pool the replicates
#' pooled_data <- pool_bs_seq_rep(bs_files)
#'
#' @export
pool_bs_seq_rep <- function(files, file_format = "encode_rrbs",
                            chr_discarded = NULL){
    assertthat::assert_that(length(files) > 1)
    message("Pooling BS-Seq replicates ...")

    # Read first file
    if (file_format == "encode_rrbs") {
        pooled_bs <- read_bs_encode_haib(file = files[1],
                    chr_discarded = chr_discarded, is_GRanges = TRUE)
    }else if (file_format == "bismark_cov") {
        pooled_bs <- .read_bs_bismark_cov(file = files[1],
                    chr_discarded = chr_discarded, is_GRanges = TRUE)
    }else{
        stop("Wrong file format. Please check the available file formats!")
    }
    for (i in 2:length(files)) {
        # Read replicate file i
        if (file_format == "encode_rrbs") {
            bs_data_rep <- read_bs_encode_haib(file = files[i],
                     chr_discarded = chr_discarded, is_GRanges = TRUE)
        }else if (file_format == "bismark_cov") {
            bs_data_rep <- .read_bs_bismark_cov(file = files[i],
                     chr_discarded = chr_discarded, is_GRanges = TRUE)
        }
        # Find overlaps between BS-Seq replicates.
        # A Hits object containing in the 1st column the query indices and in
        # the 2nd column the corresponding subject indices that overlap.
        overlaps <- GenomicRanges::findOverlaps(query   = pooled_bs,
                                                subject = bs_data_rep)
        # Get only the subset of overlapping CpG sites
        tmp_bs <- pooled_bs[S4Vectors::queryHits(overlaps)]
        # Add the total reads from the two distinct replicates
        tmp_bs$total <- tmp_bs$total +
            bs_data_rep[S4Vectors::subjectHits(overlaps)]$total
        # Add the methylated reads from the two distinct replicates
        tmp_bs$met <- tmp_bs$met +
            bs_data_rep[S4Vectors::subjectHits(overlaps)]$met
        # Add CpG sites that were not overlapping between replicates
        tmp_bs <- c(tmp_bs,
                    pooled_bs[-S4Vectors::queryHits(overlaps)],
                    bs_data_rep[-S4Vectors::subjectHits(overlaps)])
        # Sort the pooled GRanges object
        pooled_bs <- sort(tmp_bs, ignore.strand = TRUE)
    }
    return(pooled_bs)
}



#' @title (DEPRECATED) Read ENCODE HAIB bed formatted BS-Seq file
#'
#' @description (DEPRECATED) \code{read_bs_encode_haib} reads a file containing
#'   methylation data from BS-Seq experiments using the \code{\link{scan}}
#'   function. The BS-Seq file should be in ENCODE HAIB \code{bed} format. Read
#'   the Important section below on when to use this function.
#'
#' @param file The name of the file to read data values from.
#' @param chr_discarded A vector with chromosome names to be discarded.
#' @param is_GRanges Logical: if TRUE a GRanges object is returned, otherwise a
#'   data.frame object is returned.
#'
#' @return A \code{GRanges} object if \code{is_GRanges} is
#'   TRUE, otherwise a \code{data.table} object.
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
#' @seealso \code{\link{pool_bs_seq_rep}}, \code{\link{preprocess_bs_seq}}
#'
#' @examples
#' # Obtain the path to the file and then read it
#' bs_file <- system.file("extdata", "rrbs.bed", package = "BPRMeth")
#' bs_data <- read_bs_encode_haib(bs_file)
#'
#' @export
read_bs_encode_haib <- function(file, chr_discarded = NULL, is_GRanges = TRUE){
    message("Reading file ", file, " ...")
    data_raw <- scan(file = file, skip = 1, sep = "\t",
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
    bs_data <- data.table::data.table(chr = data_raw[[1]],start = data_raw[[2]],
                          strand = data_raw[[6]], total_reads = data_raw[[5]],
                          meth_reads = data_raw[[11]])
    rm(data_raw)
    # Remove selected chromosomes  -------------------------------
    bs_data <- .discard_chr(x = bs_data, chr_discarded = chr_discarded)
    # Sorting data -----------------------------------------------
    # With order priority: 1. chr, 2. start, 3. strand
    message("Sorting BS-Seq data ...")
    bs_data <- bs_data[order(bs_data$chr, bs_data$start, bs_data$strand)]
    if (is_GRanges) {
        # Create a GRanges object -----------------------------------
        message("Creating GRanges object ...")
        bs_data <- GenomicRanges::GRanges(seqnames = bs_data$chr,
                  strand = bs_data$strand,
                  ranges = IRanges::IRanges(start = bs_data$start, width = 1),
                  met = bs_data$meth_reads,
                  total = bs_data$total_reads)
    }
    return(bs_data)
}


# (DEPRECATED) Read Bismark Cov formatted BS-Seq file
#
# (DEPRECATED) \code{read_bs_bismark_cov} reads a file containing
#   methylation data from BS-Seq experiments using the
#   \code{\link[data.table]{fread}} function. The BS-Seq file should be in
#   Bismark Cov format. Read the Important section below on when to use this
#   function.
#
# @inheritParams read_bs_encode_haib
#
# @return A \code{GRanges} object if \code{is_GRanges} is
#   TRUE, otherwise a \code{data.table} object.
#
#   The GRanges object contains two additional metadata columns: \itemize{
#   \item \code{total_reads}: total reads mapped to each genomic location.
#   \item \code{meth_reads}: methylated reads mapped to each genomic location.
#   } These columns can be accessed as follows:
#   \code{granges_object$total_reads}
#
# @section Important: Unless you want to create a different workflow when
#   processing the BS-Seq data, you should NOT call this function, since this
#   is a helper function. Instead you should call the
#   \code{\link{preprocess_bs_seq}} function.
#
# @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#
# @references \url{http://rnbeads.mpi-inf.mpg.de/data/RnBeads.pdf}
#
# @seealso \code{\link{pool_bs_seq_rep}}, \code{\link{preprocess_bs_seq}}
#
# @examples
# \dontrun{
# # Download the files and change the working directory to that location
# file <- "name_of_bismark_file"
# bs_data <- read_bs_bismark_cov(file)
#
# # Extract the total reads and methylated reads
# total_reads <- bs_data$total_reads
# meth_reads <- bs_data$meth_reads
# }
.read_bs_bismark_cov <- function(file, chr_discarded = NULL, is_GRanges = TRUE){
    message("Reading file ", file, " ...")
    bs_data <- data.table::fread(input = file, sep = "\t", header = FALSE,
                col.names = c("chr", "start", "meth_reads", "unmeth_reads"))
    # Remove selected chromosomes  -------------------------------
    bs_data <- .discard_chr(x = bs_data, chr_discarded = chr_discarded)
    # Sorting data -----------------------------------------------
    # With order priority: 1. chr, 2. start
    message("Sorting BS-Seq data ...")
    bs_data <- bs_data[order(bs_data$chr, bs_data$start)]
    if (is_GRanges) {
        # Create a GRanges object ---------------------------------
        message("Creating GRanges object ...")
        bs_data <- GenomicRanges::GRanges(seqnames = bs_data$chr,
                  ranges = IRanges::IRanges(start = bs_data$start, width = 1),
                  met = bs_data$meth_reads,
                  total = bs_data$meth_reads + bs_data$unmeth_reads)
    }
    return(bs_data)
}


#' @title (DEPRECATED) Read ENCODE Caltech bed formatted RNA-Seq file
#'
#' @description (DEPRECATED) \code{read_rna_encode_caltech} reads a file
#'   containing promoter annotation data together with gene expression levels
#'   from RNA-Seq experiments using the \code{\link{scan}} function. The RNA-Seq
#'   file should be in ENCODE Caltech \code{bed} format, e.g. use \code{gtf2bed}
#'   tool if your initial file is in \code{gtf} format.
#'
#' @inheritParams read_bs_encode_haib
#'
#' @return A \code{GRanges} object if \code{is_GRanges} is TRUE, otherwise a
#'   \code{data.table} object.
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
#' @export
read_rna_encode_caltech <- function(file, chr_discarded = NULL,
                                    is_GRanges = TRUE){
    message("Reading file ", file, " ...")
    data_raw <- scan(file = file, sep = "\t",
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
        start = data_raw[[2]], end = data_raw[[3]], strand = data_raw[[6]],
        ensembl_id = data_raw[[4]],gene_expr = data_raw[[5]])
    # Extract FPKM and gene name from each gene ------------------
    gene_info <- data_raw[[10]]
    fpkm      <- vector(mode = "numeric")
    gene_name <- vector(mode = "character")

    message("Extracting FPKM and gene names ...")
    for (i in 1:length(gene_info)) {
        fpkm[i]      <- .extract_fpkm(gene_info[i])
        gene_name[i] <- .extract_gene_name(gene_info[i])
    }
    # Store extracted data to data.frame
    rna_data <- data.table::data.table(rna_data, fpkm = fpkm,
                                       gene_name = gene_name)
    rm(data_raw)
    # Remove selected chromosomes  -------------------------------
    rna_data <- .discard_chr(x = rna_data, chr_discarded = chr_discarded)
    # Sorting data -----------------------------------------------
    # With order priority: 1. chr, 2. start, 3. strand
    message("Sorting RNA-Seq data ...")
    rna_data <- rna_data[order(rna_data$chr, rna_data$start, rna_data$strand)]

    if (is_GRanges) {
        # Create a GRanges object -----------------------------------
        message("Creating GRanges object ...")
        rna_data <- GenomicRanges::GRanges(seqnames = rna_data$chr,
                       strand = rna_data$strand,
                       ranges = IRanges::IRanges(start = rna_data$start,
                                                 end = rna_data$end),
                       id = rna_data$ensembl_id,
                       gene_name  = rna_data$gene_name,
                       gene_fpkm  = rna_data$fpkm)
    }
    return(rna_data)
}
