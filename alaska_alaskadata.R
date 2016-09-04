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
    text(x = 200000, y = 1800000, "Alaska", font = my.font, col = my.textcolour, cex = my.fontsize[2])
    text(x = -1000000, y = 1300000, "Bering \n Sea", font = my.font, col = my.textcolour, cex = my.fontsize[2])
    text(x = -1200000, y = 250000, "Aleutian Islands", font = my.font, col = my.textcolour, cex = my.fontsize[2])
    text(x = 600000, y = 500000, "Gulf of \n Alaska", font = my.font, col = my.textcolour, cex = my.fontsize[2])
    par(new = TRUE)
    plot.data <- subset(data.all, !is.na(inside))
    points(spTransform(plot.data, efhCRS), pch = 1, col = rgb(0,0,0, alpha = 0.25),
           cex = (plot.data$WTCPUE / (max(plot.data$WTCPUE) / 4)))

  dev.copy(png, file.path(dir.results, "AlaskaData.png"), units = "in",
    width = my.width[2], height = my.height.map, res = my.resolution)
  dev.off(); dev.off()
