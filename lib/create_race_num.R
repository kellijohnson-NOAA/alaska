#' Create a table with RACE numbers for a given species.
#'
#' Creates a single column table with RACE numbers for a vector
#' of common names.
#' @author Kelli Faye Johnson
#' @param dir The filepath to a folder that contains the \textitt{filename}.
#' @param filename The \textitt{.csv} file that holds the lookup information.
#' @param speciesname The common names for the species that you want
#' the RACE numbers for.
#' @export

create_race_num <- function(dir = dir.data,
                            filename = "EBS_RACE_Look_2012.csv",
                            speciesnames){
    my.file <- file.path(dir, filename)
    if (!file.exists(my.file)) {
        stop(paste("The file either does not exist \n",
                   "or the wrong PATH was specified",
                   "PATH = ", my.file))
    }
    lookup.table <- read.csv(my.file)
    race.table <- data.frame(
      lookup.table[lookup.table$COMMON %in% speciesnames, "RACE"])
    colnames(race.table) <- "RACE"
    if (dim(race.table)[1] != length(speciesnames)) {
        stop(paste("One or more of the specified species names",
                   "were not in the specified lookup table."))
    }
    row.names(race.table) <- speciesnames
    return(race.table)
}
