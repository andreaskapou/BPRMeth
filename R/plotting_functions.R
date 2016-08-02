#' Plot the fit of methylation profiles across a region
#'
#' \code{plot_fitted_profiles} is a simple function for plotting the methylation
#' data across a give region, together with the fit of the methylation profiles.
#'
#' @param region Promoter region number
#' @param X Methylation data observations
#' @param fit_prof Fitted profile
#' @param fit_mean Fitted mean function
#' @param title Title of the plot
#' @param ... Additional parameters
#'
#' @return The figure to be plotted in the device.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{plot_cluster_prof}}, \code{\link{plot_scatter_gex}},
#'   \code{\link{boxplot_cluster_gex}}
#'
#' @examples
#' # Fit methylation profiles using 8 RBFs
#' obs <- meth_data
#' y   <- gex_data
#' basis <- create_rbf_object(M = 8)
#' out   <- bpr_predict_wrap(x = obs, y = y, basis = basis,
#'                           is_parallel = FALSE, opt_itnmax = 10)
#'
#' # Create the plot
#' plot_fitted_profiles(region = 16, X = meth_data, fit_prof = out)
#'
#' @export
plot_fitted_profiles <- function(region, X, fit_prof, fit_mean = NULL,
                                 title = "Gene promoter", ...){

    graphics::par(cex=1.05, mai=c(1.37,1.37,.7,.3) )
    x <- X[[region]][,1]
    y <- X[[region]][,3]/X[[region]][,2]
    xs <- seq(from = -1, to = 1, by = 0.01)
    graphics::plot(x, y, col = "blue2", pch = 21, ylim = c(0,1),
                   xlim = c(-1,1), lwd = 0.8, xlab = NA, ylab = NA,
                   cex.axis = 1.1, xaxt = "n")
    graphics::mtext(side = 1, "genomic region", line = 3, cex = 1.2)
    graphics::mtext(side = 2, "methylation level", line = 3, cex = 1.2)
    graphics::axis(side = 1, at = c(-1, 0, 1), labels=c("-7kb", "TSS", "+7kb"))
    graphics::title(main=title, line = 1, cex.main=1.4)
    if(!is.null(fit_mean)){
        graphics::lines(x = xs,
                        y = eval_probit_function(fit_mean$basis, xs,
                                                 fit_mean$W_opt[region, ]),
                        col = 'coral', lwd = 2, lty = 2)
    }
    graphics::lines(x = xs,
                    y = eval_probit_function(fit_prof$basis, xs,
                                             fit_prof$W_opt[region, ]),
                    col = 'red2', lwd = 2)
}


#' Scatter plot of predicted vs measured gene expression levels
#'
#' \code{plot_scatter_gex} creates a scatter plot of predicted gene expression
#' values on the x-axis versus the measured gene expression values on the
#' y-axis.
#'
#' @param bpr_predict_obj The output of the \code{bpr_predict_wrap} function.
#' @param main_lab The title of the plot
#' @param is_margins Use specific margins or not.
#'
#' @return The figure to be plotted in the device.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{plot_cluster_prof}}, \code{\link{plot_fitted_profiles}},
#'   \code{\link{boxplot_cluster_gex}}
#'
#' @examples
#' # Fit methylation profiles using 8 RBFs
#' obs <- meth_data
#' y   <- gex_data
#' basis <- create_rbf_object(M = 8)
#' res   <- bpr_predict_wrap(x = obs, y = y, basis = basis,
#'                           is_parallel = FALSE, opt_itnmax = 10)
#'
#' # Create the scatter plot
#' plot_scatter_gex(bpr_predict_obj = res)
#'
#' @export
plot_scatter_gex <- function(bpr_predict_obj,
                             main_lab = "Methylation Profile",
                             is_margins = TRUE){
    output <- bpr_predict_obj
    if (is_margins){
        ylim=c(-3.55, 8.02)
        xlim=c(-3.90, 7.02)
    }else{
        max_data <- max(output$test_pred, output$test$y)
        min_data <- min(output$test_pred, output$test$y)
        ylim=c(min_data, max_data)
        xlim=c(min_data, max_data)
    }

    # Compute correlation
    r <- round(stats::cor(output$test_pred, output$test$y), 2)
    # Compute RMSE
    rmse   <- round(output$test_errors$rmse, 2)
    # Perform simple linear regression using lm
    my_lm = stats::lm(output$test$y ~ output$test_pred)

    graphics::plot(output$test_pred, output$test$y, ylab = NA, xlab = NA,
                   ylim = ylim, xlim = xlim, cex.axis = 1, col = "#0000ff56",
                   pch = 16, cex = 1.3)
    graphics::mtext(side = 1, "predicted expression (log2)",
                    line = 2.2, cex = 1.1)
    graphics::mtext(side = 2, "measured expression (log2)",
                    line = 2.2, cex = 1.1)
    graphics::title(main=main_lab, line = 1.1, cex.main=1.45)
    graphics::abline(stats::coef(my_lm)[1], stats::coef(my_lm)[2],
                     col = "red", lty = 4, lwd = 3)
    # Create legend inside the figure
    graphics::legend(3.7, -1.85,
                     legend = paste0("r = ", r, "\n", "RMSE = ", rmse),
                     bty = "n", cex = 1.2)
}


#' Plot of clustered methylation profiles
#'
#' \code{plot_cluster_prof} creates a plot of cluster methylation profiles,
#' where each colour denotes a different cluster.
#'
#' @param bpr_cluster_obj The output of the \code{bpr_cluster_wrap} function.
#' @param main_lab The title of the plot
#'
#' @return The figure to be plotted in the device.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{plot_scatter_gex}},
#'   \code{\link{plot_fitted_profiles}}, \code{\link{boxplot_cluster_gex}}
#'
#' @examples
#' # Cluster methylation profiles using 4 RBFs
#' obs <- meth_data
#' basis <- create_rbf_object(M = 4)
#' res   <- bpr_cluster_wrap(x = obs, K = 3, em_max_iter = 5, opt_itnmax = 4,
#'                           init_opt_itnmax = 5, is_parallel = FALSE)
#'
#' # Create the plot
#' plot_cluster_prof(bpr_cluster_obj = res)
#'
#' @export
plot_cluster_prof <- function(bpr_cluster_obj,
                              main_lab = "Clustered methylation profiles"){
    graphics::par(mar=c(4.2, 4.1, 3.1, 2), xpd=TRUE)
    cols <- c("darkolivegreen4", "cornflowerblue",
              "coral", "firebrick","#E69F00")
    xs <- seq(-1,1,len=2000) # create some values
    graphics::plot(x = xs,
                   y = eval_probit_function(bpr_cluster_obj$basis, xs,
                                            bpr_cluster_obj$w[, 1]),
                   xlim = c(-1, 1), ylim = c(0, 1),
                   type = "l", col = cols[1], lwd = 4,
                   xlab = "promoter region",
                   ylab = "methylation level",
                   main = main_lab)
    K <- 5
    if (bpr_cluster_obj$K < 5){
        K <- bpr_cluster_obj$K
    }
    for (k in 2:K){
        graphics::lines(x = xs,
                        y = eval_probit_function(bpr_cluster_obj$basis, xs,
                                                 bpr_cluster_obj$w[, k]),
                        col = cols[k], lwd = 4)
    }
    #   graphics::legend("right", inset=c(-0.18,0), legend=seq(1,K),
    #                    lty = 1, lwd=4, col=cols[1:K], title="Cluster")
}


#' Boxplot of clustered expression levels
#'
#' \code{boxplot_cluster_gex} creates a boxplot of clustered gene expression
#' levels which depend on the clustered methylation profiles. Each colour
#' denotes a different cluster.
#'
#' @param bpr_cluster_obj The output of the \code{bpr_cluster_wrap} function.
#' @param gex The vector of gene expression data for each promoter region.
#' @param main_lab The title of the plot
#'
#' @return The figure to be plotted in the device.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{plot_cluster_prof}}, \code{\link{plot_scatter_gex}},
#'   \code{\link{plot_fitted_profiles}}
#'
#' @examples
#' # Cluster methylation profiles using 4 RBFs
#' obs <- meth_data
#' basis <- create_rbf_object(M = 4)
#' res   <- bpr_cluster_wrap(x = obs, K = 3, em_max_iter = 5, opt_itnmax = 4,
#'                           init_opt_itnmax = 5, is_parallel = FALSE)
#'
#' # Create the plot
#' boxplot_cluster_gex(bpr_cluster_obj = res, gex = gex_data)
#'
#' @export
boxplot_cluster_gex <- function(bpr_cluster_obj, gex,
                                main_lab = "Gene expression levels"){
    graphics::par(mar=c(4.2, 4.1, 3.1, 5.5), xpd=TRUE)

    cols <- c("darkolivegreen4", "cornflowerblue",
              "coral", "firebrick","#E69F00")
    gex_list <- list()
    for (k in 1:bpr_cluster_obj$K){
        gex_list[[k]] <- gex[which(bpr_cluster_obj$labels == k)]
    }

    graphics::boxplot(gex_list, col = cols[1:bpr_cluster_obj$K], notch = TRUE,
                      xlab = "Cluster K", ylab = "expression level",
                      main = main_lab)
    graphics::legend("right", inset=c(-0.18,0),
                     legend = seq(1,bpr_cluster_obj$K),
                     lty = 1, lwd = 3, col = cols[1:bpr_cluster_obj$K],
                     title = "Cluster")
}
