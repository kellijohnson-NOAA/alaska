###############################################################################
###############################################################################
## Purpose:    Stock analysis of cod in Alaska
##             Initial inputs, no analysis
## Author:     Kelli Faye Johnson
## Contact:    kellifayejohnson@gmail.com
## Date:       2016-04-13
## Comments:
###############################################################################
###############################################################################

###############################################################################
# Set initial inputs
###############################################################################
# Must set my.base or have loaded in R
my.tmb       <- "spatial_gompertz"

my.filetype <- "png"
my.width <- c(3.34646, 6.69291) # 85 mm, 170mm
my.height <- my.width[1]
my.height.map <- 4.5
my.resolution <- 500

desired.areas <- c("ai", "goa")
desired.spp <- "Pacific cod"
desired.years <- 1990:2015 #2015 = last year of data in the AI

###############################################################################
#### Initial inputs for plots
###############################################################################
my.textcolour <- "black"
my.font <- 2
my.fontsize <- c(0.4, 0.7)

my.spp.line <- -5
col.gridlines <- "dark gray"

lty.eez <- 2
lty.currents <- c(1, 3, 4)
plot.colours <- c("grey90", "grey30")
grayscale <- colorRampPalette(plot.colours)

###############################################################################
#### Load packages
###############################################################################
dir.results <- file.path(my.base, "results")
dir.data    <- file.path(my.base, "data")
dir.create(dir.results, showWarnings = FALSE)
dir.create(dir.data, showWarnings = FALSE)

if (getweb) {
  # Install personal R package: "r4kfj"
  devtools::install_github("kellijohnson/r4kfj@master")
  # Install the stable version of INLA
  install.packages("INLA", repos = "http://www.math.ntnu.no/inla/R/stable")
  # Install a needed package for rgeos
  install.packages("gpclib", repos = "http://cran.fhcrc.org/", type = "source")
}

# Install TMB
if (file.exists(file.path("c:", "adcomp"))) {
  old_wd <- getwd()
  setwd(file.path("c:", "adcomp"))
  if (getweb) {
    system("git fetch")
    system("git rebase origin/master")
  }
  remove.packages("TMB")
  source("install_windows.R", verbose = FALSE, echo = FALSE, print.eval = FALSE)
  setwd(old_wd)
  rm(old_wd)
} else {
  stop("TMB is not installed on this computer.")
}

# Load .R files specific for the alaska analysis, located in "lib" folder
ignore <- sapply(dir(file.path(my.base, "lib"), full.names = TRUE), source)
library("r4kfj", quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE)
library(INLA, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE)
load_packages(c(
  "doParallel",
  "fields",
  "foreach",
  "fpc", "ggplot2", "gstat", "igraph",
  "maps", "maptools", "mapproj", "Matrix",
  "plyr", "RandomFields",
  "raster", "reshape", "reshape2", "rgdal",
  "rgeos", "sp", "SPODT", "splancs",
  "spdep", "stats4", "TMB", "xtable"))

###############################################################################
#### Define custom inputs, such as projections and plotting
###############################################################################
define_projections()

# Save custom ggplot theme
my.theme <- theme(plot.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(colour = "white", fill = NA, size = 1),
  legend.key = element_rect(colour = "white"),
  legend.title = element_text(size = 0, face = "bold"),
  legend.text = element_text(size = 7, face = "bold")
)

#EndOfFile
