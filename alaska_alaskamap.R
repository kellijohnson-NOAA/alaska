###############################################################################
###############################################################################
## Purpose:  Stock analysis of cod in Alaska
##           Base map: AlaskaMap.png
## Author:   Kelli Faye Johnson
## Contact:  kellifayejohnson@gmail.com
## Date:     2016-09-03
## Comments:
###############################################################################
###############################################################################

png(filename = file.path(dir.results, "AlaskaMap.png"),
  width = my.width[2], height = my.height.map, res = my.resolution, units = "in")

par(xpd = TRUE, oma = c(1, 0, 0, 0), mar = c(3, 2, 0, 0))
plot(maps.eez, lty = lty.eez)
r4kfj::llgridlines(maps.eez, recenter = TRUE, lty = 1, col = col.gridlines)
for(q in seq_along(names.nt)){
  plot(eval(parse(text = paste0("names.nt_", q))),
       col = plot.colours[1], border = plot.colours[1], add = TRUE)
}
for(q in seq_along(names.nf)){
  plot(eval(parse(text = paste0("names.nf_", q))),
       col = plot.colours[2], border = plot.colours[2], add = TRUE)
}
lines(maps.ak, col = plot.colours[2])
plot(maps.eez, lty = lty.eez, add = TRUE)
legend(x = -2800000, y = 2800000,
  c("EEZ", "No trawling", "No fishing", "Alaska Coastal Current", "Alaskan Stream", "North Slope Current"),
  lty = c(lty.eez, 0, 0, lty.currents),
  pch = c(NA, 15, 15, rep(NA, length(lty.currents))),
  col = c("black", plot.colours[1], plot.colours[2], rep("black", length(lty.currents))),
  bty = "n",
  cex = my.fontsize[1], pt.cex = 1.2)
label_major(font = my.font, colour = my.textcolour, size = my.fontsize[2])
text(x = 1300000, y = 500000, "Dixon \nEntrance", font = my.font, col = my.textcolour, cex = my.fontsize[1])
## TODO:
## 2. Add lines for the currents, where the current is specified w/ a line type
text(x = -600000, y = 490000, "Unimak Pass", font = my.font, col = my.textcolour, cex = my.fontsize[1])
text(x = -1100000, y = 380000, "Samalga Pass", font = my.font, col = my.textcolour, cex = my.fontsize[1])
text(x = -1800000, y = 350000, "Amchitka\nPass", font = my.font, col = my.textcolour, cex = my.fontsize[1])
text(x = -2100000, y = 650000, "Buldir Strait", font = my.font, col = my.textcolour, cex = my.fontsize[1])

map.arrows <- data.frame(x = c(1351213, 75106.41, -444787.1, -876262.4, 361313.4,  -1636016,  -1890990, -1976426,
                               115538.2, -296584, -752666.3, -1037039, -1103498.1, -1806016,  -2024768, -944654.6),
                         y = c(715993.2, 967678.3, 599844.7, 422556.2, 815718.1, 421632.8, 505505.1, 792254.6,
                               968683, 620295, 576070.5, 499378.6, 320066.8, 542503.4, 630913.3 , 546415),
  curve = rep(c(-0.47, 0.00, 0.50, 0.76, 0.20, 0.78, 0.10, -0.20), 2),
  lty = rep(c(rep(lty.currents[1], 4), rep(lty.currents[2], 3), rep(lty.currents[3], 1)), 2))
for(q in 1:(dim(map.arrows)[1] / 2)){
  z <- (dim(map.arrows)[1] / 2) + q
  igraph:::igraph.Arrows(map.arrows[q, 1], map.arrows[q, 2], map.arrows[z, 1], map.arrows[z, 2],
    h.lwd = 1.5, sh.lwd = 1.5, curve = map.arrows[q, 3], width = 1, size = 0.4,
    h.lty = 1, sh.lty = map.arrows[q, "lty"])
}
dev.off()
