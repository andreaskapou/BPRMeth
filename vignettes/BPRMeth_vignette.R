## ----installation, echo=TRUE, eval=FALSE-----------------------------------
#  ## try http:// if https:// URLs are not supported
#  if (!requireNamespace("BiocManager", quietly=TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("BPRMeth")
#  
#  ## Or download from Github repository
#  # install.packages("devtools")
#  devtools::install_github("andreaskapou/BPRMeth", build_vignettes = TRUE)

## ---- fig.retina = NULL, fig.align='center', fig.cap="Illustration of the process for infering methylation profiles using the `BPRMeth` package.", echo=FALSE----
knitr::include_graphics("../inst/figures/model.png")

## ---- fig.retina = NULL, fig.align='center', fig.cap="Workflow diagram and functionalities of the BPRMeth package.", echo=FALSE----
knitr::include_graphics("../inst/figures/workflow.png")

## ----load_met, echo=TRUE, message=FALSE, warning=FALSE---------------------
library(BPRMeth)  # Load package
met_region <- encode_met

## ---- echo=TRUE, message=FALSE, warning=FALSE------------------------------
head(met_region$met[[20]])

## ---- echo=TRUE, message=FALSE, warning=FALSE------------------------------
head(met_region$anno)

## ---- echo=TRUE, message=FALSE, warning=FALSE------------------------------
NROW(met_region$anno)

## ----read_met_dt_pipeline, echo=TRUE, message=FALSE, eval=FALSE------------
#  # Read bulk methylation data
#  met_dt <- read_met(file = "met_file", type = "bulk_seq")
#  # Read annotation file, which will create annotation regions as well
#  anno_dt <- read_anno(file = "anno_file")
#  # Create methylation regions that have CpG coverage
#  met_region <- create_region_object(met_dt = met_dt, anno_dt = anno_dt)

## ---- echo=TRUE, message=FALSE, warning=FALSE------------------------------
expr <- encode_expr

## ---- echo=TRUE, message=FALSE, warning=FALSE------------------------------
head(expr)

## ---- echo=TRUE, message=FALSE, eval=FALSE, warning=FALSE------------------
#  # Read expression data and log2 transform them
#  expr <- read_expr(file = "expr_file", log2_transf = TRUE)

## ----create_basis, echo=TRUE, message=FALSE, warning=FALSE-----------------
# Create RBF basis object with 5 RBFs
basis_profile <- create_rbf_object(M = 5)
# Create RBF basis object with 0 RBFs, i.e. constant function
basis_mean <- create_rbf_object(M = 0)

## ----show_basis, echo=TRUE, message=FALSE, warning=FALSE-------------------
# Show the slots of the 'rbf' object
basis_profile

## ----infer_profiles, echo=TRUE, message=FALSE, warning=FALSE---------------
fit_profiles <- infer_profiles_vb(X = met_region$met, model = "binomial",
                                basis = basis_profile, is_parallel = FALSE)

fit_mean <- infer_profiles_vb(X = met_region$met, model = "binomial", 
                              basis = basis_mean, is_parallel = FALSE)

## ----show_infer_object, echo=TRUE, message=FALSE, warning=FALSE------------
# Show structure of fit_profiles object
str(fit_profiles, max.level = 1)

## ----plotprofile, fig.wide=TRUE, echo=TRUE, fig.cap="Inferring methylation profiles. Methylation pattern for specific genomic region over +/-7kb promoter region.", message=FALSE, eval=TRUE, warning=FALSE----
# Choose promoter region 21 -> i.e. LEPREL2 gene
gene_id <- met_region$anno$id[21]
p <- plot_infer_profiles(region = 21, obj_prof = fit_profiles,
                         obj_mean = fit_mean, obs = met_region$met, 
                         title = paste0("Gene ID ", gene_id))
print(p)

## ----predictexpr, echo=TRUE, message=FALSE, warning=FALSE, warning=FALSE----
# Set seed for reproducible results
set.seed(22)
# Perform predictions using methylation profiles
pred_profile <- predict_expr(prof_obj = fit_profiles, expr = expr,
                             anno = met_region$anno, model_name = "svm",
                             fit_feature = "RMSE", cov_feature = TRUE, 
                             is_summary = FALSE)
# Perform predictions using mean methylation rate
pred_mean <- predict_expr(prof_obj = fit_mean, expr = expr, 
                          anno = met_region$anno, model_name = "svm",
                          is_summary = FALSE)

## ----l9, echo=TRUE, message=FALSE, warning=FALSE---------------------------
print(paste("Profile r: ", round(pred_profile$test_errors$pcc, 2)))
print(paste("Mean r:", round(pred_mean$test_errors$pcc, 2)))

## ----plotpred, echo=TRUE, fig.cap="Relationship between DNA methylation patterns and gene expression. Scatter plots of predicted versus measured (log2-transformed) gene expression values using profiles (left) and rates (right); each shaded blue dot represents a different gene, The red line indicates the linear fit between the predicted and measured expression values.", fig.wide=TRUE, message=FALSE, eval=TRUE, warning=FALSE----
g1 <- plot_predicted_expr(pred_obj = pred_profile, 
                          title = "Methylation profile")
g2 <- plot_predicted_expr(pred_obj = pred_mean, 
                          title = "Methylation rate")
g <- cowplot::plot_grid(g1, g2, ncol = 2, nrow = 1)
print(g)

## ----clustering, echo=TRUE, message=FALSE, eval=TRUE, warning=FALSE--------
# Set seed for reproducible results
set.seed(12) 
# Create basis object
basis_obj <- create_rbf_object(M = 3)
# Set number of clusters K
K <- 4
# Perform clustering
cl_obj <- cluster_profiles_vb(X = met_region$met, K = K, model = "binomial", 
                              alpha_0 = .5, beta_0 = .1,
                              basis = basis_obj, vb_max_iter = 20)

## ----clusterplot, fig.wide=TRUE, echo=TRUE, fig.cap="Clustering methylation profiles across promoter-proximal regions.", message=FALSE, eval=TRUE, warning=FALSE----
cluster_plot <- plot_cluster_profiles(cluster_obj = cl_obj)
print(cluster_plot)

## ----infer_profiles_beta, fig.wide = TRUE, echo=TRUE, message=FALSE, warning=FALSE----
basis_obj <- create_rbf_object(M = 5)
# Infer beta profiles
beta_prof <- infer_profiles_mle(X = beta_data, model = "beta",
                                basis = basis_obj, is_parallel = FALSE)
p <- plot_infer_profiles(region = 15, obj_prof = beta_prof,
                         obs = beta_data, title = "Beta profile")
print(p)

## ----infer_profiles_bernoulli, fig.wide=TRUE, echo=TRUE, message=FALSE, warning=FALSE----
basis_prof <- create_rbf_object(M = 5)
basis_mean <- create_rbf_object(M = 0)
bern_prof <- infer_profiles_vb(X = bernoulli_data, model = "bernoulli",
                               basis = basis_prof, is_parallel = FALSE)
bern_mean <- infer_profiles_vb(X = bernoulli_data, model = "bernoulli",
                               basis = basis_mean, is_parallel = FALSE)
p <- plot_infer_profiles(region = 50, obj_prof = bern_prof, 
                         obj_mean = bern_mean, obs = bernoulli_data, 
                         title = "Bernoulli profile")
print(p)

## ----cluster_bernoulli, fig.wide=TRUE, echo=TRUE, fig.cap="Clustering Bernoulli methylation profiles.", message=FALSE, eval=TRUE, warning=FALSE----
# Perform clustering
cl_obj <- cluster_profiles_vb(X = bernoulli_data, K = 3, model = "bernoulli", 
                              basis = basis_obj, vb_max_iter = 40)
                              
cluster_plot <- plot_cluster_profiles(cluster_obj = cl_obj)
print(cluster_plot)

## ----session_info, echo=TRUE, message=FALSE--------------------------------
sessionInfo()

