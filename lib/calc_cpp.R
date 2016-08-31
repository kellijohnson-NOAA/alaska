#' @param cpp A character entry giving the name of the \code{.cpp} file
#' @param loc The directory location where \code{cpp} is located
calc_cpp <- function(cpp, loc) {
  ###############################################################################
  #### Unload and load cpp
  ###############################################################################
  # First attempt
  firsttry <- try(dyn.unload(dynlib(file.path(cpp))), silent = TRUE)
  # Garbage Collection
  if(is(firsttry, "try-error")) gc(verbose = FALSE)
  # Second attempt
  secondtry <- try(dyn.unload(dynlib(file.path(cpp))), silent = TRUE)

  # Check if current directory is same as file location
  if (getwd() != loc) {
    firsttry <- try(dyn.unload(dynlib(file.path(loc, cpp))), silent = TRUE)
    if(is(firsttry, "try-error")) gc(verbose = FALSE)
    secondtry <- try(dyn.unload(dynlib(file.path(loc, cpp))), silent = TRUE)
    # Copy the file over to the current directory
    ignore <- file.copy(file.path(loc, paste0(cpp, ".cpp")),
      paste0(cpp, ".cpp"))
  }

  # Compile
  compile(paste0(cpp, ".cpp"))
  dyn.load(dynlib(cpp))
}
