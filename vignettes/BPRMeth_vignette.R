## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"---------------------------------
BiocStyle::latex()

## ----eval=TRUE, echo=TRUE, message=FALSE, results="asis"--------------------------------
library(BPRMeth)
rrbs_file <- system.file("extdata", "rrbs.bed", package = "BPRMeth")
rnaseq_file <- system.file("extdata", "rnaseq.bed", package = "BPRMeth")

## ----eval=TRUE, echo=TRUE, message=FALSE, results="asis"--------------------------------
# Preprocess both RRBS and RNA-Seq files
HTS_data <- process_haib_caltech_wrap(rrbs_file, rnaseq_file)

## ----eval=TRUE, echo=TRUE---------------------------------------------------------------
HTS_data$methyl_region[[16]]

## ----eval=TRUE, echo=TRUE---------------------------------------------------------------
head(HTS_data$gex, 10)

## ----eval=TRUE, echo=TRUE---------------------------------------------------------------
HTS_data$rna_data

## ----eval=TRUE, echo=TRUE---------------------------------------------------------------
# Obtain the number of gene promoters
length(HTS_data$gex)

## ----eval=TRUE, echo=TRUE---------------------------------------------------------------
# Create basis object with 9 RBFs
basis_profile <- create_rbf_object(M = 9)

# Create basis object with 0 RBFs, i.e. constant function
basis_mean <- create_rbf_object(M = 0)

## ----eval=TRUE, echo=TRUE---------------------------------------------------------------
# Show the slots of the 'rbf' object
basis_profile

## ----eval=TRUE, message=FALSE, echo=TRUE------------------------------------------------

# Set seed for reproducible results
set.seed(1234)


# Perform predictions using methylation profiles
res_profile <- bpr_predict_wrap(x = HTS_data$methyl_region, y = HTS_data$gex,
                                basis = basis_profile, fit_feature = "RMSE",
                                cpg_dens_feat = TRUE, is_parallel = FALSE,
                                is_summary = FALSE)

# Perform predictions using mean methylation level
res_mean <- bpr_predict_wrap(x = HTS_data$methyl_region, y = HTS_data$gex,
                             basis = basis_mean, is_parallel = FALSE,
                             is_summary = FALSE)

## ----eval=TRUE, echo=TRUE---------------------------------------------------------------
# Test errors for methylation profiles, PCC = Pearson's r
res_profile$test_errors$pcc

# Test errors for mean methylation levels
res_mean$test_errors$pcc

## ----figureexample, fig.show='hide', fig.width=7.8, fig.height=5------------------------
# Choose promoter region 21 -> i.e. LEPREL2 gene
gene_name <- HTS_data$rna_data$gene_name[21]
plot_fitted_profiles(region = 21, X = HTS_data$methyl_region, fit_prof = res_profile,
                     fit_mean = res_mean, title = paste0("Gene ", gene_name))

## ----figureexample2, fig.show='hide', fig.width=12.8, fig.height=5.8--------------------
par(mfrow=c(1,2))
plot_scatter_gex(bpr_predict_obj = res_profile)
plot_scatter_gex(bpr_predict_obj = res_mean, main_lab = "Mean Methylation")

## ----eval=TRUE, message=FALSE, warning=FALSE, echo=TRUE---------------------------------
# Set seed for reproducible results
set.seed(1234)
# Create basis object with 4 RBFs
basis_obj <- create_rbf_object(M = 4)
# Set number of clusters K = 5
K <- 5
# Perform clustering
res <- bpr_cluster_wrap(x = HTS_data$methyl_region, K = K, basis = basis_obj,
                        em_max_iter = 15, opt_itnmax = 30, is_parallel = TRUE)

## ----figureexample3, warning=FALSE, fig.show='hide', fig.width=13.5, fig.height=5-------
par(mfrow=c(1,2))
plot_cluster_prof(bpr_cluster_obj = res)
boxplot_cluster_gex(bpr_cluster_obj = res, gex = HTS_data$gex)

## ---------------------------------------------------------------------------------------
sessionInfo()

