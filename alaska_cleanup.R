###############################################################################
###############################################################################
## Purpose:    Clean up some objects and do garbage collection
## Author:     Kelli Faye Johnson
## Contact:    kellifayejohnson@gmail.com
## Date:       2016-09-08
## Comments:
###############################################################################
###############################################################################

rm(list = ls(pattern = "maps\\."))
gcinfo(verbose = FALSE)
ignore <- gc(verbose = FALSE)
