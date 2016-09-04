###############################################################################
###############################################################################
## Purpose:  Stock analysis of cod in Alaska
##           Create spde for spatial analysis
## Author:   Kelli Faye Johnson
## Contact:  kellifayejohnson@gmail.com
## Date:     2014-09-18
## Comments: use \code{calc_table} to create word tables.
##           Latex and pandoc notes:
##           1. parameter names == "$\\sigma_{\\epsilon}$"
###############################################################################
###############################################################################
n.areas <- length(unique(data.spp$inside))



logfile.spp <- sapply(desired.spp, function(x){
                      paste(tolower(substring(unlist(strsplit(x, " ")),
                                              1, 1)),
                            collapse = "")
                      })
logfile.spp <- logfile.spp[order(logfile.spp)]
###############################################################################
#### Read in results data
###############################################################################
rdata_files <- dir(dir.results, pattern = "log", full.names = TRUE)
rdata_species <- do.call("rbind", strsplit(rdata_files, "_"))[, 2]
rdata_files <- rdata_files[rdata_species %in% logfile.spp]
for(q in seq_along(logfile.spp)){
  all <- grep(logfile.spp[q], rdata_files)
  remove <- rev(order(file.info(rdata_files[all])$mtime))[-1]
  rdata_files <- rdata_files[-all[remove]]
}

my.res <- list()
for(q in seq_along(desired.spp)){
  load(rdata_files[q])
  assign(logfile.spp[q], saved)
  my.res[[q]] <- eval(parse(text = logfile.spp[q]))
}
names(my.res) <- logfile.spp

###############################################################################
#### Results
###############################################################################
sink(file.path(dir.results, "spde_results.txt"), append = FALSE)
cat("number of mesh nodes", "\n")
cat(my.res[[1]]$mesh$n, "\n\n\n")
cat("Decimal degrees to minutes seconds for stock splits", "\n",
    "of my.res[[1]], which is currently Pacific cod")
apply(cbind(my.res[[1]]$stock_all$splits[, "index"], 60),
          1, function(x) dd2dms(as.vector(x)))
sink()

###############################################################################
###############################################################################
#### fig_stock.png
###############################################################################
# make_file(my.filetype, file.path(dir.results, "fig_stock.png"),
#           width = my.width[1], res = my.resolution)
par(mfrow = c(length(my.res), 1))
par(xpd = TRUE, fig = c(0,1,0,1))
for(q in seq_along(my.res)) {
  plotcp(my.res[[q]]$stock_all, las = 1, upper = "splits")
  points(tail(my.res[[q]]$stock$cptable, 1)[, c("nsplit", "xerror")] + c(1, 0),
         col = "black", pch = 19, cex = 3)

  par(new = TRUE, fig = c(0.2,0.9,0.45,0.95))
  plot(my.res[[q]]$stock_all, uniform = TRUE)
  text(my.res[[q]]$stock_all, use.n = TRUE, all = FALSE,
       digits = 3, cex = 0.8, splits = TRUE,
       fancy = FALSE, fwidth = 0.6, fheight = 0.75)
}
dev.copy(png, file.path(dir.results, "fig_stock.png"),
         res = my.resolution)
dev.off();dev.off()
#make_file_off()

###############################################################################
#### Mean abundance
###############################################################################
make_file(my.filetype, file.path(dir.results, "hat_abundance.png"), width = 4,
    height=2.5*length(desired.spp), res = my.resolution)
  par(mfrow = c(length(desired.spp), 1),
      mar = c(0, 0, 0, 0), oma = c(3, 3, 0, 0),
      mgp = c(1.5,0.5,0), tck = -0.02, xaxs = "i")
  for(q in seq_along(desired.spp)){
    # Load old results
    ylim <- c(0, range(my.res[[q]]$B_conf_spatial)[2])
    xlim <- c(desired.years, rev(desired.years)[1] + 1)
    plot(1, type = "n", xlab = "", ylab = "", lwd = 3, ylim = ylim,
         xlim = range(xlim), log = "", las = 1, xaxt = "n")
    Which = match(c("rho", "phi"), rownames(my.res[[q]]$opt$summary) )
    val <- formatC(round(my.res[[q]]$opt$summary[Which[1], "Estimate"], 3), 3)
    se <- formatC(round(my.res[[q]]$opt$summary[Which[1], "Std. Error"], 3), 3)
    mtext(side = 1, adj = 0, line = -1.5, text = bquote(rho == .(val) (.(se))))
    val <- formatC(round(my.res[[q]]$opt$summary[Which[2],
                               "Estimate"], 3), 3)
    se <- formatC(round(my.res[[q]]$opt$summary[Which[2],
                               "Std. Error"], 3), 3)
    mtext(side = 1, adj = 0, line = -2.5, text = bquote(phi == .(val) (.(se))))
    if(length(desired.spp) > 1){
      plot.text <- paste0("(", letters[q], ")", " ", desired.spp[q])
      mtext(plot.text, side = 3, adj = 1, line = -2, font = my.font)
    }
    polygon(x = c(xlim, rev(xlim)),
            y = c(my.res[[q]]$B_conf_spatial[,1],
                  rev(my.res[[q]]$B_conf_spatial[,2])), lty = "solid",
            col=rgb(0, 0, 0, 0.2), border = NA)
    lines(x = xlim, y = my.res[[q]]$B_mean_spatial, lwd = 3)

  Dji <- my.res[[q]]$opt$summary[rownames(opt$summary) == "Dji", ]
  stockset <- my.res[[q]]$stock$where
  Dji <- cbind(rep(xlim, each = length(stockset)), Dji)
  colnames(Dji)[1] <- "year"

  stock.1 <- Dji[my.res[[q]]$stock$where == unique(stockset)[1], ]
  stock.1[, "Estimate"] <- exp(stock.1[, "Estimate"])
  stock.1.u <- tapply(stock.1[, "Estimate"], stock.1[, "year"], mean)
  stock.lines <- list(stock.1.u)
  if(length(stockset) > 1){
    stock.2 <- Dji[my.res[[q]]$stock$where == unique(stockset)[2], ]
    stock.2[, "Estimate"] <- exp(stock.2[, "Estimate"])
    stock.2.u <- tapply(stock.2[, "Estimate"], stock.2[, "year"], mean)
    stock.lines[[2]] <- stock.2.u
    lines(x = xlim,
          y = tapply(stock.1[, "Estimate"], stock.1[, "year"], mean),
          lty = 2)
    lines(x = xlim,
          y = tapply(stock.2[, "Estimate"], stock.2[, "year"], mean),
          lty = 3)
  }
 }
 axis(1)
 mtext(side = 1, "years", outer = TRUE, line = 2)
 mtext(side = 2, expression(bar(lnBiomass)), outer = TRUE, line = 1.5)
 legend("bottomright", c(paste("longitude >=",
                               round(my.res[[q]]$stock$splits[1, "index"], 2)),
                         paste("longitude <",
                               round(my.res[[q]]$stock$splits[1, "index"], 2))),
        lty = c(2, 3), bty = "n")
make_file_off()


###############################################################################
#### fig spatial variation in productivity
###############################################################################

plot.coords <- data.frame("x" = my.res[[q]]$x_stations,
                          "y" = my.res[[q]]$y_stations,
                          "omega" = my.res[[q]]$Report[["Omega"]])
coordinates(plot.coords) <- ~ x + y
proj4string(plot.coords) <- akCRS
plot.grid  <- SpatialGrid(GridTopology(cellcentre.offset = bbox(plot.coords)[, "min"],
                                       cellsize = c(50, 50),
                                       cells.dim = c(62, 17)),
                          proj4string = akCRS)

us <- getData("GADM", country="USA", level=1)
# extract states (need to uppercase everything)
ak <- us[match(toupper("Alaska"), toupper(us$NAME_1)),]
ak <- spTransform(ak, akCRS)

# bubble(plot.coords, "omega")
idw <- krige(omega ~ 1, plot.coords, plot.grid)

make_file(my.filetype, file.path(dir.results, "fig_omega.png"),
          width = my.width, height = my.height, res = my.resolution)
par(oma = c(0,0,0,0), mar = c(0,0,0,0))
image(idw["var1.pred"], col= gray.colors(101),
      ylim = c(0, 1500))
plot(ak, add = TRUE, col = "white")
abline(v = my.res[[1]]$stock_orig$splits[, "index"],
       lty = c(1, rep(2, dim(my.res[[1]]$stock_orig$splits)[1])))
# r4kfj::llgridlines(ak, recenter = TRUE)
par(new = TRUE)
plot.new()
par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
plot.window(xlim = round(range(idw["var1.pred"]@data), 1), ylim = c(0, 1))
points(x = seq(round(range(idw["var1.pred"]@data)[1], 1),
               round(range(idw["var1.pred"]@data)[2], 1),
               length.out = 101),
       y = rep(0.33, 101), col = gray.colors(101), pch = 15, cex = 1.2)
axis(1, line = -5.4, cex.axis = 0.5, mgp = c(0.3, 0.25, 0))
make_file_off()

plot(crop.image(idw["var1.pred"], rbind( c(-2000, 300),
             c(500, 1200))))
axis(2)
#### TABLES ####
################
###############################################################################
###############################################################################



###############################################################################
## Code that is not currently used
## But too good to throw away
###############################################################################
NPFMCstrata.afsc <- spTransform(readOGR(file.path(dir.data, "stratum"),
  "strata_geographic"), akCRS)


###############################################################################
## population trend for each STRATUM
###############################################################################
make_file(my.filetype, file.path(dir.results, "trend.png"),
          width = my.width[2], height = my.height.map * length(desired.spp),
          res = my.resolution)
par(mfrow = c(length(desired.spp), 1), oma = c(0,5,0,0), mar = c(4,5.5,0,0),
    xpd = TRUE)
  for(q in seq(desired.spp)) {
    growth <- list()
    my.coords <- list()
    use.data <- subset(data.all, SID == race.num[q, 1] & !is.na(inside))
    use.data <- spTransform(use.data, efhCRS)
    stations.unique.RACE <- unique(use.data@data$STRATUM)
 for(s in seq_along(stations.unique.RACE)) {
   md.data <- use.data[which(use.data@data$STRATUM == stations.unique.RACE[s]), ]
   my.coords[[s]] <- coordinates(md.data[1, ])
   md.data <- md.data@data
   if(dim(md.data)[1] < 2) next
   if(length(unique(md.data$YEAR)) == 1) next
     mylm <- lm(log(WTCPUE) ~ YEAR, data = md.data)
   growth[[s]] <- summary(mylm)$coefficients["YEAR",]
   names(growth)[s] <- stations.unique.RACE[s]
 }
  growth <- unlist(lapply(growth, "[[", 1))
  points.shape <- ifelse(growth > 0, 0, 1) # 0 = square 1 = circle
   plot(maps.ak, col = plot.colours[2])
    r4kfj::llgridlines(maps.eez, recenter = TRUE, lty = 1, col = col.gridlines)
    text(x = 200000, y = 1800000, "Alaska", font = my.font, col = my.textcolour, cex = my.fontsize[2])
    text(x = -1000000, y = 1300000, "Bering \n Sea", font = my.font, col = my.textcolour, cex = my.fontsize[2])
    text(x = -1200000, y = 250000, "Aleutian Islands", font = my.font, col = my.textcolour, cex = my.fontsize[2])
    text(x = 600000, y = 500000, "Gulf of \n Alaska", font = my.font, col = my.textcolour, cex = my.fontsize[2])
    points(do.call("rbind", my.coords),
           pch = points.shape, cex = abs(growth*15))
      plot.text <- paste0("(", letters[q], ")", " ", desired.spp[q])
      if(length(desired.spp) > 1){
            mtext(plot.text, side = 3, adj = 1, line = my.spp.line, font = my.font,
                  col = my.textcolour, cex = 1.5)
          }
      if(q == 1){
        legend("topleft", legend = c("positive", "negative"),
               pch = c(0, 1), bty = "n")
      }
}
make_file_off()

###############################################################################
#### table_data
###############################################################################
# split data between GOA and AI, and place in a list with length 2

# create the table_data dataframe
  f <- function(x) c(mean(x), length(x))
table_data <- apply(race.num, 1, function(d) {
    data <- subset(data.spp, SID == d & YEAR %in% desired.years)
    area <- do.call("rbind", strsplit(rownames(data@data), ".", fixed = TRUE))[, 1]
    temp <- list("ai" = data[grep("i", area), ],
                 "goa" = data[grep("o", area), ])

    data.mid <- lapply(temp, function(x) {
      temp.2 <- do.call(rbind, tapply(x@data[, c("WTCPUE")], x@data$YEAR, f))
      data.frame("area" = eval.parent(quote(names(X)))[substitute(x)[[3]]],
                 "year" = rownames(temp.2),
                 "n" = temp.2[, 2],
                 "meanCPUE" = temp.2[, 1])
    })
    do.call("rbind", data.mid)
  })
table_data <- do.call("rbind", lapply(table_data, function(d) {
              species <- rep(eval.parent(quote(names(X)))[substitute(d)[[3]]],
                             dim(d)[1])
              cbind(species, d)
              }))
colnames(table_data)[match("meanCPUE", colnames(table_data))] <- "$\\bar{CPUE}$"

table_data$year <- factor(table_data$year, min(as.numeric(as.character(table_data$year))) :
                                           max(as.numeric(as.character(table_data$year))))
table_data <- table_data[order(table_data$species, table_data$year), ]
xtable_data <- xtable(table_data, digits = c(0,0,0,0,0,2))
print(xtable_data, file = file.path(dir.results, "table_data.tex"),
      include.rownames = FALSE, sanitize.text.function = function(x){x},
      )
currwd <- getwd()
setwd(dir.results)
system("pandoc table_data.tex -o table_data.docx")
setwd(currwd)

###############################################################################
#### Stock partitioning table
###############################################################################
table_stock <- lapply(my.res, function(x) x$stock$cptable)
xtable_stock <- xtable(table_stock)
digits(xtable_stock) <- 4
print(xtable_stock, file = file.path(dir.results, "table_stock.tex"),
      include.rownames = FALSE, sanitize.text.function = function(x){x})
currwd <- getwd()
setwd(dir.results)
system("pandoc table_stock.tex -o table_stock.docx")
setwd(currwd)

df.stock <- list()
for(r in seq_along(desired.spp)){
df.stock[[r]] <- my.res[[r]]$stock$frame[, c("n", "dev", "yval", "complexity")]
if(dim(my.res[[r]]$stock$frame)[1] > 1){
  pos <- data.frame(x = my.res[[r]]$stock$splits[1, "index"],
                    y = mean(my.res[[r]]$y_stations))
  pos <- SpatialPoints(pos, proj4string = akCRS)
  posll <- spTransform(pos, llCRS)

  df.stock[[r]]$name <- c(desired.spp[r],
                          paste0(c(">=", "<"),
                                 round(coordinates(posll)[1, 1],
                                              digits = 2)))
} else {df.stock[[r]]$name <- desired.spp[r]}
df.stock[[r]] <- df.stock[[r]][, c("name", "n", "yval", "dev", "complexity")]
}
(df.stock <- do.call("rbind", df.stock))
colnames(df.stock) <- c("name", "n", "mean sigmaO", "dev", "complexity")

df.stock.xtable <- xtable(df.stock)
 print(df.stock.xtable, type = "html",
       file = file.path(dir.results, "stock.html"),
       include.colnames = TRUE, include.rownames = FALSE,
       html.table.attributes = getOption("xtable.html.table.attributes",
                                    "border = 0"))



 # Calculate esimates of random effects (RE)
  if(TRUE){

    y.lim <- range(mesh$loc[, 2])
    x.lim <- range(mesh$loc[, 1])
    col.use <- rev(heat.colors(100)[26:75])

    # Omega
    Omega_est <- Report[["Omega"]]
      Rel <- ((Omega_est - min(Omega_est)) /
              diff(range(Omega_est)))
    png(file = file.path(dir.results, paste0("Omega_est_", logfile.spp,".png")),
        width = 4, height = 4, res = 200, units = "in")

      plot(eval(parse(text = my.shape)), ylim = y.lim, xlim = x.lim,
           col = "grey90", border = "grey90", main = "", mar = c(0, 0, 2.5, 0))

      points(x = x_stations, y = y_stations, pch = 20, col = col.use[round(length(col.use)*Rel)])
        box(bty = "o", lwd = 2)
        legend("topleft", legend = c("high", "low"), bty = "n",
               pch = 20, col = col.use[c(length(col.use), 1)])
        r4kfj::llgridlines(eval(parse(text = my.shape)), recenter = TRUE)
    dev.off()

    # Equilibrium field
    Equil_est = Report[["Equil"]]
    # Epsilon
    Epsilon_est <- Report[["Epsilon"]]
    Nrow <- ceiling(sqrt(n_years))
    Ncol <- ceiling(n_years/Nrow)
    png(file = file.path(dir.results, paste0("Epsilon_est_", logfile.spp, ".png")),
        width = 2 * Ncol, height = 2 * Nrow, res = 200, units = "in")
      par(mfrow = c(Nrow, Ncol), mar = c(0, 0, 0, 0), mgp = c(1.5, 0.5, 0), tck = -0.02)
      for(i in 1:n_years){
        Rel = ((Epsilon_est[, i] - min(Epsilon_est[, i])) / diff(range(Epsilon_est[, i])))
      plot(eval(parse(text = my.shape)), ylim = y.lim, xlim = x.lim,
           col = "grey90", border = "grey90", main = "", mar = c(0, 0, 2.5, 0))
        points(x = x_stations, y = y_stations, pch = 20,
               col = col.use[round(length(col.use)*Rel)])
        box(bty = "o", lwd = 2)
        if(i == 1){
          legend("topleft", legend = c("high", "low"), bty = "n",
                 pch = 20, col = col.use[c(length(col.use), 1)])
        }
      }
    dev.off()

    D_est <- Report[["Dji"]]
    Nrow <- ceiling(sqrt(n_years))
    Ncol <- ceiling(n_years / Nrow)
    png(file = file.path(dir.results, paste0("D_est_", logfile.spp, ".png")),
        width = 2 * Ncol, height = 2 * Nrow, res = 200, units = "in")
      par(mfrow = c(Nrow, Ncol), mar = c(0, 0, 0, 0), mgp = c(1.5, 0.5, 0), tck = -0.02)
      for(i in 1:n_years){
        Rel <- ((D_est[, i]- min(D_est[, i])) / diff(range(D_est[, i])))
      plot(eval(parse(text = my.shape)), ylim = y.lim, xlim = x.lim,
           col = "grey90", border = "grey90", main = "", mar = c(0, 0, 2.5, 0))
        points(x = x_stations, y = y_stations, pch = 20, col = col.use[round(length(col.use) * Rel)])
         box(bty = "o", lwd = 2)
        if(i == 1){
          legend("topleft", legend = c("high", "low"), bty = "n",
                 pch = 20, col = col.use[c(length(col.use), 1)])
        }
      }
    dev.off()

    # Site specific abundance
    png(file = file.path(dir.results, paste0("Abundance_bysite_", logfile.spp, ".png")),
        width = 4, height = 4, res = 200, units = "in")
      par(mar = c(4, 4, 2, 0), mgp = c(2.25, 0.5, 0), tck = -0.01)
      matplot(t(D_est), type = "l", col = "black", lty = "solid", las = 1,
              xlab = "year", ylab = "ln(abundance)")
    dev.off()
    # Compare my calculations and TMB
    plot(y = colMeans(exp(D_est)), x = B_mean_spatial)
    #rowMeans( exp(Report[["Dji"]]) )
  }



###############################################################################
#### JTT code
###############################################################################

png(file.path(DateFile,"Abundance_all.png"), width=1*3, height=2.5*length(SpeciesSet), res=200, units="in")
  par(mfrow=c(length(SpeciesSet),1), mar=c(2,3,0,0), mgp=c(1.5,0.5,0), tck=-0.02, oma=c(2,2,1.5,0), xaxs="i")
  for(SpeciesI in 1:length(SpeciesSet)){
    # Load old results
    spp = SpeciesSet[SpeciesI]   # Rho = 0.98, 0.94, 0.53/0.29
    load(file=paste(DateFile,"/Save_spatial_",logfile.spp,".RData",sep=""))
    Ylim = range(saved$B_conf_spatial)
    plot(1, type="n", xlab="", ylab="", lwd=3, ylim=Ylim, xlim=2003+c(1,10), log="")
    text(x=2003+2, y=Ylim[2], pos=1, labels="rho = \nphi = \nAIC =")
    mtext(side=2, line=2, outer=FALSE, text=switch(SpeciesSet[SpeciesI], "NumVermC"="Vermilion", "NumGpot"="Greenspotted", "NumBoc"="Boccacio", "NumVerm"="Vermilion", "NumSunset"="Sunset"))
    # Spatial
    polygon( x=2003+c(1:10,10:1), y=c(saved$B_conf_spatial[,1],rev(saved$B_conf_spatial[,2])), lty="solid", col=rgb(1,0,0,0.2), border=NA)
    lines( x=2003+c(1:10), y=saved$B_mean_spatial, col="red", lwd=3)
    Which = match( c("rho","phi"), rownames(saved$opt$summary) )
    text(x=2003+4, y=Ylim[2], pos=1, col="red", labels=paste(formatC(round(saved$opt$summary[Which[1],"Estimate"],3),3)," (",formatC(round(saved$opt$summary[Which[1],"Std. Error"],3),3),")\n",formatC(round(saved$opt$summary[Which[2],"Estimate"],3),3)," (",formatC(round(saved$opt$summary[Which[2],"Std. Error"],3),3),")\n",formatC(round(2*saved$opt$objective+2*length(saved$opt),1),-1),sep=""))
    n_years = length(unique(Data$Year))
    Ymat = matrix(NA, nrow=length(unique(Data$Site)), ncol=n_years)
    for(YearI in 1:n_years){
    for(SiteI in 1:length(unique(Data$Site))){
      Which = which(Data$Year==unique(Data$Year)[YearI]&Data$Site==unique(Data$Site)[SiteI])
      if(length(Which)>=1) Ymat[SiteI,YearI] = Data[Which[1],logfile.spp]
    }}
    Y = as.vector(Ymat)    # Use Gaussian measurement error
    Site = as.vector(row(Ymat))
    Year = as.vector(col(Ymat))
    ( Glm = glm(Y ~ 0 + factor(Year) + factor(Site), family="poisson") )
    B_pred = predict(Glm, newdata=data.frame( "Year"=Year, "Site"=Site ), type="response")
    B_hat = tapply( B_pred, INDEX=Year, FUN=mean)
    lines( x=2003+c(1:10), y=B_hat, lwd=2, col="black", lty="dotted")
    if(SpeciesI==1) mtext(side=3, expression("Average density"), line=0, outer=FALSE)
  }
  mtext(side=1, line=0, outer=TRUE, text="Year")
dev.off()

# Plot together -- Omega
y.lim = c( 32, 34.8)     # range(Data[,"Lat..DD.DDDDD."])
x.lim = c( -121, -117 )  # range(Data[,"Lon..DDD.DDDDD."])
Omega_range = c(-1.2,0.9)
redblue2 <- colorRampPalette(colors=(c("red","purple","blue","green","yellow")))
Seq = seq(-1.2,0.9,length=5)
png(file.path(DateFile, "Omega_all.png"), width = 1 * 3,
    height = 2.5 * length(SpeciesSet), res = 200, units = "in")
  par(mfrow=c(length(SpeciesSet), 1), mar = c(0,0,0,0), mgp = c(1.5,0.5,0),
      tck = -0.02, oma = c(3,3,1.5,0), xaxs = "i")
  for(SpeciesI in 1:length(SpeciesSet)){
    spp = SpeciesSet[SpeciesI]   # Rho = 0.98, 0.94, 0.53/0.29
    load(file = paste0(DateFile, "/Save_spatial_", spp, ".RData"))
    if(TRUE){
      saved[["x_stations"]] = x_stations
      saved[["y_stations"]] = y_stations
      saved[["mesh"]] = mesh
    }
    Omega_est = saved$Report[["Omega"]][saved$mesh$idx$loc]
    print(range(Omega_est))
    Rel = ((Omega_est[saved$mesh$idx$loc]-min(Omega_est[saved$mesh$idx$loc]))/diff(range(Omega_est[saved$mesh$idx$loc])))
    #Rel = (Omega_est - min(Omega_range)) / diff(Omega_range)
    map("worldHires", ylim=y.lim, xlim=x.lim, col="grey90", fill=T, mar=c(0,0,0,0), main="", myborder=0.01) # myborder is buffer around limits
    points(x=saved$x_stations, y=saved$y_stations, pch=20,
           col = redblue2(21)[round(21*Rel)])
    print( table( round(21*Rel) ))
    if(SpeciesI==length(SpeciesSet)) axis(1)                                                                    # mar=c(0,0,2,3),
    axis(2)
    box(bty="o",lwd=1)
  }
  mtext(side=1, "Longitude", line=1.5, outer=TRUE)
  mtext(side=2, "Latitude", line=1.5, outer=TRUE)
  mtext(side=3, expression(Omega), line=0, outer=TRUE)
dev.off()
# Now with GLM
png(file.path(DateFile,"Omega_all_GLM.png"), width=1*3, height=2.5*length(SpeciesSet), res=200, units="in")
  par(mfrow=c(length(SpeciesSet),1), mar=c(0,0,0,0), mgp=c(1.5,0.5,0), tck=-0.02, oma=c(3,3,1.5,0), xaxs="i")
  for(SpeciesI in 1:length(SpeciesSet)){
    spp = SpeciesSet[SpeciesI]   # Rho = 0.98, 0.94, 0.53/0.29
    map("worldHires", ylim=y.lim, xlim=x.lim, col="grey90", fill=T, mar=c(0,0,0,0), main="", myborder=0.01) # myborder is buffer around limits
    #points(x=saved$x_stations, y=saved$y_stations, pch=20, col=redblue2(21)[round(21*Rel)])
    print( table( round(21*Rel) ))
    if(SpeciesI==length(SpeciesSet)) axis(1)                                                                    # mar=c(0,0,2,3),
    axis(2)
    box(bty="o",lwd=1)
    # GLM
    n_years = length(unique(Data$Year))
    Ymat = matrix(NA, nrow=length(unique(Data$Site)), ncol=n_years)
    for(YearI in 1:n_years){
    for(SiteI in 1:length(unique(Data$Site))){
      Which = which(Data$Year==unique(Data$Year)[YearI]&Data$Site==unique(Data$Site)[SiteI])
      if(length(Which)>=1) Ymat[SiteI,YearI] = Data[Which[1],spp]
    }}
    Y = as.vector(Ymat)    # Use Gaussian measurement error
    Site = as.vector(row(Ymat))
    Year = as.vector(col(Ymat))
    ( Glm = glm(Y ~ 0 + factor(Site) + factor(Year), family="poisson") )
    Match = match(unique(Data$Site),Data$Site)
    Rel = (Glm$coef[1:length(unique(Data$Site))] - mean(Glm$coef[1:length(unique(Data$Site))])) / diff(Omega_range)
    Rel = round(21*Rel)
    Rel = ifelse(Rel>21,21,Rel)
    Rel = ifelse(Rel<1,1,Rel)
    points( y=Data[Match,"Lat..DD.DDDDD."], x=Data[Match,"Lon..DDD.DDDDD."], pch=20, col=redblue2(21)[Rel] )
  }
  mtext(side=1, "Longitude", line=1.5, outer=TRUE)
  mtext(side=2, "Latitude", line=1.5, outer=TRUE)
  mtext(side=3, expression(Omega), line=0, outer=TRUE)
dev.off()


# Plot individually
for(SpeciesI in 1:length(SpeciesSet)){
  png(paste0(DateFile,"/Abundance_pretty_",SpeciesSet[SpeciesI],".png"), width=1*4, height=3*1, res=200, units="in")
    par(mfrow=c(1,1), mar=c(2,2,0,0), mgp=c(1.5,0.5,0), tck=-0.02, oma=c(1,0,1,0), xaxs="i")
    # Load old results
    spp = SpeciesSet[SpeciesI]   # Rho = 0.98, 0.94, 0.53/0.29
    load(file=paste0(DateFile,"/Save_spatial_",spp,".RData"))
    Ylim = range(saved$B_conf_spatial)
    plot( 1, type="n", xlab="", ylab="", lwd=3, ylim=Ylim, xlim=2003+c(1,10), log="")
    text(x=2003+2, y=Ylim[2], pos=1, labels="rho = \nphi = \nAIC =")
    mtext(side=3, line=0, outer=FALSE, text=c("Simulated")[SpeciesI])
    # Spatial
    polygon( x=2003+c(1:10,10:1), y=c(saved$B_conf_spatial[,1],rev(saved$B_conf_spatial[,2])), lty="solid", col=rgb(1,0,0,0.2), border=NA)
    lines( x=2003+c(1:10), y=saved$B_mean_spatial, col="red", lwd=3)
    Which = match( c("rho","phi"), rownames(saved$opt$summary) )
    text(x=2003+4, y=Ylim[2], pos=1, col="red", labels=paste(formatC(round(saved$opt$summary[Which[1],"Estimate"],3),3)," (",formatC(round(saved$opt$summary[Which[1],"Std. Error"],3),3),")\n",formatC(round(saved$opt$summary[Which[2],"Estimate"],3),3)," (",formatC(round(saved$opt$summary[Which[2],"Std. Error"],3),3),")\n",formatC(round(2*saved$opt$objective+2*length(saved$opt),1),-1),sep=""))
   # GLM
    n_years = length(unique(Data$Year))
    Ymat = matrix(NA, nrow = length(unique(Data$Site)), ncol=n_years)
    for(YearI in 1:n_years){
    for(SiteI in 1:length(unique(Data$Site))){
      Which = which(Data$Year==unique(Data$Year)[YearI]&Data$Site==unique(Data$Site)[SiteI])
      if(length(Which)>=1) Ymat[SiteI,YearI] = Data[Which[1],spp]
    }}
    Y = as.vector(Ymat)    # Use Gaussian measurement error
    Site = as.vector(row(Ymat))
    Year = as.vector(col(Ymat))
    ( Glm = glm(Y ~ 0 + factor(Year) + factor(Site), family="poisson") )
    B_pred = predict(Glm, newdata=data.frame( "Year"=Year, "Site"=Site ), type="response")
    B_hat = tapply( B_pred, INDEX=Year, FUN=mean)
    lines( x=2003+c(1:10), y=B_hat, lwd=2, col="black", lty="dotted")
    mtext(side=1, line=1.3, outer=FALSE, text="Year")
  dev.off()
}
