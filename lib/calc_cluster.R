#' Calculate the cluster
#'
#' https://stat.ethz.ch/R-manual/R-devel/library/cluster/html/pam.html
#' https://stat.ethz.ch/R-manual/R-devel/library/cluster/html/silhouette.html
#' http://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
calc_cluster <- function(dist) {
  cluster::pam(dist)
  fpc::pamk(dist)

}
