#' Subset dataframe for certain species by common name
#'
#' Methods to subset a dataframe by the "SID" column.
#' Users can supply a single column dataframe or a vector of racenums,
#' to the function which will then be used to subset a dataframe.
#' @author Kelli Faye Johnson
#' @param racenums A single column dataframe or a vector of numeric values.
#' @param data A dataframe containing a column called "SID" containing
#' @param racenums.all A vector of all RACE numbers that you care about for this study
#' race numbers.
#' @export

    create_data_spp <- function(racenums, data, racenums.all = race.num$RACE,
                                zeroobs = "NA") {
    	if (is.data.frame(racenums)) {
            if (dim(racenums)[2] > 1) {
            	stop("racenums must be a vector or a dataframe of 1 column.")
            }
    		racenums <- as.vector(racenums[, 1])
    	}
        ## Subset the data for the current racenum and values inside the 
        ## study area
        data.new <- droplevels(subset(data,
                                      SID %in% racenums & !is.na(inside)))
        ## Get all station x year combinations in the dataframe
        ## that were sampled
        data.long <- subset(melt(reshape::cast(data, 
                                   formula = STATION ~ YEAR,
                                   fun.aggregate = length, 
                                   value = "WTCPUE", fill = NA), 
                                 id = "STATION", variable_name = "Year"),
                                 !is.na(value) & 
                                 STATION %in% data.new$STATION)

        notthere <- list(); counter <- 0
        for(q in seq(dim(data.long)[1])) {
            if(!(data.long[q, "STATION"] %in% 
                 subset(data.new, YEAR == data.long[q, "YEAR"])$STATION)) {
                counter <- counter + 1
                notthere[[counter]] <- data.long[q, ]
            }
        }
        if(length(notthere) > 0 & zeroobs == "0.001") {
          data.needed <- data[match(with(do.call("rbind", notthere), 
                                         paste(STATION, YEAR, sep = "_")),
                                    with(data.inside, 
                                         paste(STATION, YEAR, sep = "_"))), ]
          data.needed$SID <- racenums
          data.needed$id <- with(data.needed, paste(STATION, SID, sep = "_"))
          data.needed$WTCPUE <- as.numeric(zeroobs)
          data.needed$NUMCPUE <- 0
          data.needed$COMMON <- data.new$COMMON[1]
          data.needed$SCIENTIFIC <- data.new$SCIENTIFIC[1]
        
          data.return <- rbind(data.new, data.needed)
          return(data.return)
        } else {
          return(data.new)
        }

    }