#' Calculate the cluster
#'
#' https://stat.ethz.ch/R-manual/R-devel/library/cluster/html/pam.html
#' https://stat.ethz.ch/R-manual/R-devel/library/cluster/html/silhouette.html
#' http://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
calc_cluster <- function(data, return = c("basic", "detailed"), ...) {

  return <- match.arg(return, choices = c("basic", "detailed"))

  get <- c("Longitude", "Latitude", "Omega_em")
  if (any(!c(get, "Year", "local") %in% colnames(data))) {
    stop(c(get, "Year", "local"), " need to be columns of data.")
  }
  temp <- data[data$Year == 1 & data$local, get]
  if ("weight" %in% names(list(...))) {
    dist <- calc_dist(temp, ...)
  } else {
    dist <- calc_dist(temp, weight = 1.0)
  }

  # Could also use cluster::pam(), but have to specify k
  pamk <- fpc::pamk(dist)
  temp$cluster <- pamk$pamobject$clustering
  # If the true Omega exists calculate best clusters
  if ("Omega_om" %in% colnames(data)) {
    dist <- data[data$Year == 1 & data$local,
      c("Longitude", "Latitude", "Omega_om")]
    if ("weight" %in% names(list(...))) {
      dist <- calc_dist(dist, ...)
    } else {
      dist <- calc_dist(dist, weight = 1.0)
    }
    object <- fpc::pamk(dist)
    temp$groupcluster <- object$pamobject$clustering
  }
  temp <- temp[, -which(colnames(temp) == "Omega_em")]
  data <- merge(data, temp, by = c("Longitude", "Latitude"))

  # If the correct columns exist return a better indexing
  if ("group" %in% colnames(data)) {
    temp <- data[, c("group", "cluster")]
    colnames(temp) <- c("true", "est")
    data$cluster <- calc_matches(temp)$estgroup
    if ("groupcluster" %in% colnames(data)) {
      temp <- data[, c("group", "groupcluster")]
      colnames(temp) <- c("true", "est")
      data$groupcluster <- calc_matches(temp)$estgroup
    }
  }

  if (return == "basic") {
    return(data)
  }
  if (return == "detailed") {
    return(list("em" = pamk, "om" = ifelse(exists("object"), object, NULL)))
  }
}
