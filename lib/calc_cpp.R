#' Load a \code{\link{TMB}} library.
#'
#' The function first checks to see if the library is loaded and then
#' if so, it removes the library prior to attempting to compile and load it.
#'
#' @param cpp A character entry giving the name of the \code{.cpp} file
#' @param loc The directory location where \code{cpp} is located
calc_cpp <- function(cpp, loc) {
  ###############################################################################
  #### Unload and load cpp
  ###############################################################################
  # Make sure that the cpp character value does not have the extension
  if(grepl("\\.[[:alpha:]]{3}", cpp)) {
    cpp <- substr(cpp, 1, nchar(basename(cpp)) - 4)
  }
  # First attempt
  firsttry <- try(dyn.unload(TMB::dynlib(file.path(loc, cpp))), silent = TRUE)
  # Garbage Collection
  if(is(firsttry, "try-error")) gc(verbose = FALSE)
  # Second attempt
  secondtry <- try(dyn.unload(TMB::dynlib(file.path(loc, cpp))), silent = TRUE)

  # Compile
  TMB::compile(file.path(loc, paste0(cpp, ".cpp")))
  dyn.load(TMB::dynlib(file.path(loc, cpp)))
}
