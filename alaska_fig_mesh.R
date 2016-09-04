###############################################################################
###############################################################################
## Purpose:  Stock analysis of cod in Alaska
##           Base map: AlaskaData.png
## Author:   Kelli Faye Johnson
## Contact:  kellifayejohnson@gmail.com
## Date:     2016-09-03
## Comments:
###############################################################################
###############################################################################

###############################################################################
## Alaska and survey data: AlaskaData.png
###############################################################################
par(oma = c(0,0,0,0), mar = c(4,5.5,0,0), xpd = TRUE)
plot(spTransform(maps.ak, akCRS), col = rgb(0,0,0,0))
plot(mesh, add = TRUE)
par(new = TRUE)
plot(maps.ak, col = rgb(0,0,0,0))
lines(maps.ak, col = plot.colours[2])
r4kfj::llgridlines(maps.eez, recenter = TRUE, lty = 1, col = col.gridlines)
label_major(font = my.font, colour = my.textcolour, size = my.fontsize[2])

dev.copy(png, file.path(dir.results, "AlaskaMesh.png"), units = "in",
  width = my.width[2], height = my.height.map, res = my.resolution)
dev.off(); dev.off()
