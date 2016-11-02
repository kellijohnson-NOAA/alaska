#' Read a saved R data set into a new environment and return the environment.
#'
#' The new environment can be subsetted using the dollar sign and it will
#' not replace any names in the global environment.
#'
#' @param RData A saved \code{.RData} file using the full path, or the
#' file must be in the current working directory.
#' @param env A specified environment for which to save the read in file.

read_env <- function(RData, env = new.env()){
  load(RData, env)
  return(env)
}
