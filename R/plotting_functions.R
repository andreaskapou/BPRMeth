# Define ggplot2 theme
.gg_theme <- function(){
    p <- theme(
        plot.title = element_text(size = 20,face = 'bold',
                                  margin = margin(0,0,3,0), hjust = 0.5),
        axis.text = element_text(size = rel(1.05), color = 'black'),
        axis.title = element_text(size = rel(1.45), color = 'black'),
        axis.title.y = element_text(margin = margin(0,10,0,0)),
        axis.title.x = element_text(margin = margin(10,0,0,0)),
        axis.ticks.x = element_line(colour = "black", size = rel(0.8)),
        axis.ticks.y = element_blank(),
        legend.position = "right",
        legend.key.size = unit(1.4, 'lines'),
        legend.title = element_text(size = 15, face = 'bold'),
        legend.text = element_text(size = 13),
        panel.border = element_blank(),
        panel.grid.major = element_line(colour = "gainsboro"),
        panel.background = element_blank()
    )
    return(p)
}

#' @title Plot inferred methylation profiles across a region
#'
#' @description Function for plotting the inferred methylation profiles across a
#'   given region, and optionally the mean methylation rate together with the
#'   observed methylation data, using \code{\link{ggplot2}}.
#'
#' @param region Genomic region number
#' @param obj_prof Inferred profile, i.e. output from
#'   \code{\link{infer_profiles_vb}} or \code{\link{infer_profiles_mle}}
#' @param obj_mean Inferred mean function, i.e. output from
#'   \code{\link{infer_profiles_vb}} or \code{\link{infer_profiles_mle}}
#' @param obs Methylation data observations, if a list, will extract the
#'   specific \code{region}, if a matrix will plot directly its observations.
#' @param title Plot title
#' @param x_axis x axis label
#' @param y_axis x axis label
#' @param x_labels x axis ticks labels
#' @param ... Additional parameters
#'
#' @return A ggplot2 object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @examples
#' # Fit methylation profiles using 5 RBFs
#' basis <- create_rbf_object(M = 5)
#' prof <- infer_profiles_vb(X = encode_met$met, model = "binomial",
#'     basis = basis, is_parallel = FALSE, vb_max_iter = 5)
#' # Create the plot
#' g <- plot_infer_profiles(region = 16, obj_prof = prof, obs = encode_met$met)
#'
#' @seealso \code{\link{plot_predicted_expr}},
#' \code{\link{plot_cluster_profiles}},
#'   \code{\link{boxplot_cluster_expr}}
#'
#' @export
plot_infer_profiles <- function(region = 1, obj_prof, obj_mean = NULL,
                                obs = NULL, title = "Inferred profiles",
                                x_axis = "genomic region", y_axis = "met level",
                                x_labels = c("Upstream", "", "Centre", "", "Downstream"),
                                ...) {
    aes_xs <- seq(from = -1, to = 1, by = 0.01)
    # For RMD CHECK to pass without NOTEs
    aes_ys = x = y = Model = ys_prof = ys_mean <- NULL
    ys_low = ys_high <- 0
    ys_low_prof = ys_high_prof = ys_low_mean = ys_high_mean <- 0
    if (methods::is(obj_prof, "infer_profiles_mle")) {
        if (methods::is(obj_prof, "infer_profiles_mle_gaussian")) {
            ys_prof <- eval_function(obj_prof$basis,aes_xs,obj_prof$W[region,])
        }else{
            ys_prof <- eval_probit_function(obj_prof$basis, aes_xs,
                                            obj_prof$W[region, ])
        }
        if (!is.null(obj_mean)) {
            if (methods::is(obj_mean, "infer_profiles_mle_gaussian")) {
                ys_mean <- eval_function(obj_mean$basis, aes_xs,
                                         obj_mean$W[region, ])
            }else{
                ys_mean <- eval_probit_function(obj_mean$basis, aes_xs,
                                                obj_mean$W[region, ])
            }
        }
    }else if (methods::is(obj_prof, "infer_profiles_vb")) {
        tmp <- .predictive_infer_profile(obj_prof, aes_xs, region)
        ys_prof <- tmp$W_pred
        if (methods::is(obj_prof, "infer_profiles_vb_binomial") ||
            methods::is(obj_prof, "infer_profiles_vb_bernoulli")) {
            ys_low_prof <- ys_prof - ys_prof*(1 - ys_prof);
            ys_high_prof <- ys_prof + ys_prof*(1 - ys_prof)
        }else if (methods::is(obj_prof, "infer_profiles_vb_gaussian")) {
            ys_low_prof <- ys_prof - 2 * tmp$W_sd_pred;
            ys_high_prof <- ys_prof + 2 * tmp$W_sd_pred
        }

        if (!is.null(obj_mean)) {
            tmp <- .predictive_infer_profile(obj_mean, aes_xs, region)
            ys_mean <- tmp$W_pred
            if (methods::is(obj_mean, "infer_profiles_vb_binomial") ||
                methods::is(obj_mean, "infer_profiles_vb_bernoulli")) {
                ys_low_mean <- ys_mean - ys_mean*(1 - ys_mean);
                ys_high_mean <- ys_mean + ys_mean*(1 - ys_mean)
            }else if (methods::is(obj_mean, "infer_profiles_vb_gaussian")) {
                ys_low_mean <- ys_mean - 2 * tmp$W_sd_pred;
                ys_high_mean <- ys_mean + 2 * tmp$W_sd_pred
            }
        }
    }else{
        stop("No plotting function for this model!")
    }

    dt <- data.table::data.table(aes_xs = aes_xs, aes_ys = ys_prof,
                                 ys_low = ys_low_prof, ys_high = ys_high_prof)
    if (is.null(obj_mean)) {
        p <- ggplot(dt, aes(x = aes_xs, y = aes_ys)) +
            geom_line(aes(x = aes_xs, y = aes_ys), size = 1.5, col = "darkblue")
        if (methods::is(obj_prof, "infer_profiles_vb") ||
            methods::is(obj_prof, "infer_profiles_gibbs") ) {
            p <- p + geom_ribbon(dt, mapping = aes(ymin = ys_low,
             ymax = ys_high), alpha = 0.2, size = 0.1, fill = "cornflowerblue")
        }
    }else{
        dt <- dt %>% .[, c("Model") := list("Profile")]
        dt_mean <- data.table::data.table(aes_xs = aes_xs, aes_ys = ys_mean,
                                          ys_low = ys_low_mean,
                                          ys_high = ys_high_mean,
                                          Model = "Mean")
        dt <- rbind(dt, dt_mean)
        if (methods::is(obj_prof, "infer_profiles_vb") ||
            methods::is(obj_prof, "infer_profiles_gibbs") ) {
            p <- ggplot(dt, aes(x = aes_xs, y = aes_ys, color = Model)) +
                geom_line(size = 1.5) +
                geom_ribbon(dt, mapping = aes(ymin = ys_low, ymax = ys_high,
                                              fill = Model),
                            alpha = 0.2, size = 0.1) +
                scale_color_brewer(palette = "Set1") +
                scale_fill_brewer(palette = "Set1")
        }else{
            p <- ggplot(dt, aes(x = aes_xs, y = aes_ys, color = Model)) +
                geom_line(size = 1.5) +
                scale_color_brewer(palette = "Set1") +
                scale_fill_brewer(palette = "Set1")
        }
    }
    if (!is.null(obs)) {
        # If we are given a list object
        if (is.list(obs)) { obs <- obs[[region]] }
        # If we have binomial observations
        if (NCOL(obs) == 3) {
            dt_obs <- data.table::data.table(x = obs[,1], y = obs[, 3]/obs[, 2])
        }else{
            dt_obs <- data.table::data.table(x = obs[,1], y = obs[, 2])
        }
        p <- p + geom_point(data = dt_obs, mapping = aes(x = x, y = y),
                            shape = 1, color = "red", size = 3)
    }
    # If nto gaussian data set y-lim to (0, 1)
    if (!methods::is(obj_prof, "infer_profiles_mle_gaussian") &
        !methods::is(obj_prof, "infer_profiles_vb_gaussian")) {
        p <- p + scale_y_continuous(limits = c(0, 1))
    }
    p <- p + scale_x_continuous(limits = c(-1, 1), labels = x_labels) +
        labs(title = title, x = x_axis, y = y_axis) +
        .gg_theme()
    return(p)
}


#' @title Plot clustered methylation profiles across a region
#'
#' @description Function for plotting the clusterd methylation profiles across a
#'   given region where each colour denotes a different cluster.
#'
#' @param cluster_obj Clustered profiles object, i.e. output from
#'   \code{\link{cluster_profiles_vb}} or \code{\link{cluster_profiles_mle}}
#'   functions.
#' @inheritParams plot_infer_profiles
#'
#' @return A ggplot2 object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{plot_infer_profiles}},
#'   \code{\link{plot_predicted_expr}}, \code{\link{boxplot_cluster_expr}}
#'
#' @examples
#' # Cluster methylation profiles using 3 RBFs
#' basis <- create_rbf_object(M = 3)
#' # Perform clustering
#' cl_obj <- cluster_profiles_vb(X = encode_met$met, K = 3, model = "binomial",
#'            basis = basis, vb_max_iter = 5)
#' # Create plot
#' g <- plot_cluster_profiles(cluster_obj = cl_obj)
#'
#' @export
plot_cluster_profiles <- function(cluster_obj, title = "Clustered profiles",
                                x_axis = "genomic region", y_axis = "met level",
                                x_labels = c("Upstream", "", "Centre", "", "Downstream"),
                                ...) {
    # Test data
    aes_xs <- seq(from = -1, to = 1, by = 0.01)
    # Number of clusters
    K <- NCOL(cluster_obj$W)
    ys <- matrix(0, ncol = K, nrow = length(aes_xs))
    # For RMD CHECK to pass without NOTEs
    aes_ys = x = y = Cluster <- NULL
    ys_low = ys_high <- NULL
    if (methods::is(cluster_obj, "cluster_profiles_mle")) {
        if (methods::is(cluster_obj, "cluster_profiles_mle_gaussian")) {
            for (k in 1:K) {
                ys[, k] <- eval_function(cluster_obj$basis, aes_xs,
                                         cluster_obj$W[,k])
            }
        }else{
            for (k in 1:K) {
                ys[, k] <- eval_probit_function(cluster_obj$basis, aes_xs,
                                         cluster_obj$W[,k])
            }
        }
    }else if (methods::is(cluster_obj, "cluster_profiles_vb")) {
        tmp <- .predictive_cluster_profile(cluster_obj, aes_xs)
        ys <- tmp$W_pred
        if (methods::is(cluster_obj, "cluster_profiles_vb_binomial") ||
            methods::is(cluster_obj, "cluster_profiles_vb_bernoulli")) {
            ys_low <- ys - ys*(1 - ys);
            ys_high <- ys + ys*(1 - ys)
        }else if (methods::is(cluster_obj, "cluster_profiles_vb_gaussian")) {
            ys_low <- ys - 2 * tmp$W_sd_pred;
            ys_high <- ys + 2 * tmp$W_sd_pred
        }
    }else{
        stop("No plotting function for this model!")
    }

    dt <- data.table::data.table(aes_xs = numeric(), aes_ys = numeric(),
            ys_low = numeric(), ys_high = numeric(), Cluster = numeric())
    if (methods::is(cluster_obj, "cluster_profiles_vb") ||
        methods::is(cluster_obj, "cluster_profiles_gibbs") ) {
        for (k in 1:K) {
            dt <- rbind(dt, data.table::data.table(aes_xs = aes_xs,
                aes_ys = ys[,k], ys_low = ys_low[,k], ys_high = ys_high[,k],
                Cluster = as.factor(k)))
        }
    }else{
        for (k in 1:K) {
            dt <- rbind(dt, data.table::data.table(aes_xs = aes_xs,
                                                   aes_ys = ys[,k],
                                                   ys_low = 0, ys_high = 0,
                                                   Cluster = as.factor(k)))
        }
    }

    p <- ggplot(dt, aes(x = aes_xs, y = aes_ys, color = Cluster)) +
        geom_line(size = 2)
    if (methods::is(cluster_obj, "cluster_profiles_vb") ||
        methods::is(cluster_obj, "cluster_profiles_gibbs") ) {
        p <- p + geom_ribbon(dt, mapping = aes(ymin = ys_low, ymax = ys_high,
                 fill = Cluster), alpha = 0.2, size = 0.1)
    }
    # If nto gaussian data set y_lim to (0, 1)
    if (!methods::is(cluster_obj, "cluster_profiles_mle_gaussian") &
        !methods::is(cluster_obj, "cluster_profiles_vb_gaussian")) {
        p <- p + scale_y_continuous(limits = c(0, 1))
    }

    p <- p + scale_x_continuous(limits = c(-1, 1), labels = x_labels) +
        scale_color_brewer(palette = "Dark2") +
        scale_fill_brewer(palette = "Dark2") +
        labs(title = title, x = x_axis, y = y_axis) + .gg_theme()
    return(p)
}




#' @title Scatter plot of predicted vs measured gene expression levels
#'
#' @description \code{plot_predicted_expr} creates a scatter plot of predicted
#'   gene expression values on the x-axis versus the measured gene expression
#'   values on the y-axis.
#'
#' @param pred_obj The output of the \code{predict_expr} function.
#' @param title The title of the plot.
#' @param is_margins Use specified margins or not.
#'
#' @return A ggplot2 object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{plot_infer_profiles}},
#'   \code{\link{plot_cluster_profiles}}, \code{\link{boxplot_cluster_expr}}
#'
#' @examples
#' # Fit methylation profiles using 5 RBFs
#' basis <- create_rbf_object(M = 5)
#' prof <- infer_profiles_vb(X = encode_met$met, model = "binomial",
#'     basis = basis, is_parallel = FALSE, vb_max_iter = 5)
#' # Predict expression
#' pred_obj <- predict_expr(prof_obj = prof, expr = encode_expr,
#'   anno = encode_met$anno, model_name = "lm", is_summary = FALSE)
#' # Create plot
#' g <- plot_predicted_expr(pred_obj = pred_obj)
#'
#' @export
plot_predicted_expr <- function(pred_obj, title = "Predicted expression",
                                is_margins = FALSE){
    pred = meas <- NULL
    if (is_margins) {
        ylim = c(-3.35, 8.02)
        xlim = c(-3.85, 7.02)
    }else{
        max_data <- max(pred_obj$test_pred, pred_obj$test$y)
        min_data <- min(pred_obj$test_pred, pred_obj$test$y)
        ylim = c(min_data, max_data)
        xlim = c(min_data, max_data)
    }
    # Compute correlation
    r <- round(stats::cor(pred_obj$test_pred, pred_obj$test$y), 2)
    # Compute RMSE
    rmse <- round(pred_obj$test_errors$rmse, 2)
    # CReate data.frame
    out_plot <- data.frame(pred = pred_obj$test_pred, meas = pred_obj$test$y)
    # Create gpplot
    g <- ggplot(out_plot, aes(x = pred, y = meas)) +
        geom_point(pch = 16, col = "#0000ff56", cex = 3) +
        labs(x = "predicted expression (log2)",y = "measured expression (log2)",
             title = title) +
        geom_smooth(method = lm, se = FALSE, col = "red") +
        scale_y_continuous(limits = ylim) +
        scale_x_continuous(limits = xlim) +
        geom_text(data = data.frame(),
                  aes(label = paste0("  r = ",r,"\n","RMSE = ",rmse)),
                  x = -Inf, y = Inf, hjust = -0.2, vjust = 1.2) +
        .gg_theme()
    return(g)
}

#' @title Boxplot of clustered expression levels
#'
#' @description Create a boxplot of clustered gene expression levels which
#'   depend on the clustered methylation profiles. Each colour denotes a
#'   different cluster.
#'
#' @param cluster_obj The output from \code{\link{cluster_profiles_vb}} or
#'   \code{\link{cluster_profiles_mle}} functions.
#' @param expr The expression data object.
#' @param anno The annotation data object.
#' @inheritParams plot_infer_profiles
#'
#' @return A ggplot2 object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{plot_cluster_profiles}},
#' \code{\link{plot_infer_profiles}}, \code{\link{plot_predicted_expr}}
#'
#' @examples
#' # Cluster methylation profiles using 3 RBFs
#' basis <- create_rbf_object(M = 3)
#' # Perform clustering
#' cl_obj <- cluster_profiles_vb(X = encode_met$met, K = 3, model = "binomial",
#'            basis = basis, vb_max_iter = 5)
#' # Create plot
#' g <- boxplot_cluster_expr(cluster_obj = cl_obj, expr = encode_expr,
#'        anno = encode_met$anno)
#'
#' @export
boxplot_cluster_expr <- function(cluster_obj, expr, anno,
                                 title = "Expression levels"){
    my_min <- function(x){ xx <- min(x); return(xx - 1.3) }
    # Init so RMD CHECK does not complain
    id = N = Cluster <- NULL
    # Append feature ID so we can merge with expression data later
    dt <- data.table::data.table(cbind(anno$id, cluster_obj$labels))
    colnames(dt) <- c("id", "Cluster")
    # Make sure that the first column name of expression matrix is called "id"
    if (!identical(colnames(expr)[1], "id")) {
        stop("First column of expression matrix must have name 'id'.")
    }
    # Merge met profiles with expression data by id
    dt <- merge(dt, expr, by = "id") %>% .[, N := length(id), by = "Cluster"]

    ggplot(dt, aes(x = Cluster, y = expr)) +
        geom_boxplot(aes(fill = Cluster), outlier.shape = 1, notch = FALSE) +
        scale_color_brewer(palette = "Dark2") +
        scale_fill_brewer(palette = "Dark2") +
        stat_summary(fun.y = my_min,
                     aes(label = paste0("", N)),
                     geom = 'text', lwd = 5, col = 'black', cex = 5) +
        .gg_theme() + theme(axis.text.x = element_blank()) +
        labs(list(title = title, y = "expression level"))
}
