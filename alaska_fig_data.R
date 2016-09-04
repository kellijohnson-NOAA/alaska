###############################################################################
###############################################################################
## Purpose:  Stock analysis of cod in Alaska
##           Base map: fig_data.png
## Author:   Kelli Faye Johnson
## Contact:  kellifayejohnson@gmail.com
## Date:     2016-09-03
## Comments:
###############################################################################
###############################################################################

###############################################################################
#### figure_data
###############################################################################
png(filename = file.path(dir.results, "fig_data.png"),
  width = my.width[2], height = my.height, res = my.resolution, units = "in")
data.plot <- subset(data.all, SID == race.num[1,1] & YEAR %in% desired.years)
data.plot$year.f <- factor(data.plot$YEAR, levels = desired.years)

col.year <- rep(NA, length(desired.years))
for(y in seq_along(desired.years)){
  if(!desired.years[y] %in% data.plot$YEAR) {next}
  temp <- subset(data.plot, YEAR == desired.years[y])@data[1, "survey"]
  col.year[y] <- ifelse(length(grep("ai", temp)) > 0, 1, 2)
}
par(mgp = c(1, 0.5, 0), oma = c(0, 0, 0, 0), mar = c(3, 3, 0, 0))
boxplot(log(WTCPUE) ~ year.f,
  data = data.plot, lty = col.year, las = 1, ann = FALSE)
mtext("year", side = 1, line = 1.5)
mtext("ln(CPUE)", side = 2, line = 1.5)
text(seq_along(desired.years), rep(par("yaxp")[1] - 1, length(desired.years)),
       labels = summary(data.plot$year.f), cex = 0.65)
text(0.15, par("yaxp")[1] - 1, "n =", cex = 0.65)
legend("topleft", lty = 1:2, legend = c("AIs", "GOA"),
       bty = "n", horiz = TRUE)
dev.off()
