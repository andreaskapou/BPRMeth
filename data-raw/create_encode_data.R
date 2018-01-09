library(BPRMeth)
rrbs_file <- system.file("extdata", "rrbs.bed", package = "BPRMeth")
rnaseq_file <- system.file("extdata", "rnaseq.bed", package = "BPRMeth")
chr_discarded <- NULL
min_bs_cov <- 4
max_bs_cov <- 1000
chrom_size <- NULL

# Process BS-Seq file and return data in the required format
met_data  <- preprocess_bs_seq(files = rrbs_file, file_format = "encode_rrbs",
                             chr_discarded = chr_discarded,
                             min_bs_cov = min_bs_cov, max_bs_cov = max_bs_cov)
# Read RNA-Seq BED file
rna_data <- read_rna_encode_caltech(file = rnaseq_file,
                                    chr_discarded = chr_discarded,
                                    is_GRanges = TRUE)
# Create promoter regions
prom_reg <- BPRMeth:::.genomic_region(anno = rna_data, chrom_size = chrom_size,
                                      upstream = -7000, downstream = 7000)
encode_met <- create_met_region(met_data, prom_reg, cov = 10, sd_thresh = 10e-2)
encode_expr <- data.table(id = rna_data$id, expr = log2(rna_data$gene_fpkm + 0.1))
# encode_data <- list(met_region = met_region, expr = expr)
set.seed(1)
devtools::use_data(encode_met, encode_expr, overwrite = TRUE)
