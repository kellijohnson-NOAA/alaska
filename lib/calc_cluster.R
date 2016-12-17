#' Calculate the cluster
#'
#' https://stat.ethz.ch/R-manual/R-devel/library/cluster/html/pam.html
#' https://stat.ethz.ch/R-manual/R-devel/library/cluster/html/silhouette.html
#' http://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
calc_cluster <- function(data, ...) {
  get <- c("Longitude", "Latitude", "Omega_em")
  if (any(!get %in% colnames(data))) stop(get, " need to be columns",
    " of data.")
  temp <- data[data$Year == 1 & data$local, get]
  if ("weight" %in% names(list(...))) {
    dist <- calc_dist(temp, ...)
  } else {
    dist <- calc_dist(temp, weight = 1.0)
  }
  # Could also use cluster::pam(), but have to specify k
  pamk <- fpc::pamk(dist)
  temp$cluster <- pamk$pamobject$clustering
  temp <- temp[, -which(colnames(temp) == "Omega_em")]
  data <- merge(data, temp, by = c("Longitude", "Latitude"))

  return(data)
}
