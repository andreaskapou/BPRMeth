% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{plot_cluster_profiles}
\alias{plot_cluster_profiles}
\title{Plot clustered methylation profiles across a region}
\usage{
plot_cluster_profiles(
  cluster_obj,
  title = "Clustered profiles",
  x_axis = "genomic region",
  y_axis = "met level",
  x_labels = c("Upstream", "", "Centre", "", "Downstream"),
  ...
)
}
\arguments{
\item{cluster_obj}{Clustered profiles object, i.e. output from
\code{\link{cluster_profiles_vb}} or \code{\link{cluster_profiles_mle}}
functions.}

\item{title}{Plot title}

\item{x_axis}{x axis label}

\item{y_axis}{x axis label}

\item{x_labels}{x axis ticks labels}

\item{...}{Additional parameters}
}
\value{
A ggplot2 object.
}
\description{
Function for plotting the clusterd methylation profiles across a
  given region where each colour denotes a different cluster.
}
\examples{
# Cluster methylation profiles using 3 RBFs
basis <- create_rbf_object(M = 3)
# Perform clustering
cl_obj <- cluster_profiles_vb(X = encode_met$met, K = 3, model = "binomial",
           basis = basis, vb_max_iter = 5)
# Create plot
g <- plot_cluster_profiles(cluster_obj = cl_obj)

}
\seealso{
\code{\link{plot_infer_profiles}},
  \code{\link{plot_predicted_expr}}, \code{\link{boxplot_cluster_expr}}
}
\author{
C.A.Kapourani \email{C.A.Kapourani@ed.ac.uk}
}
