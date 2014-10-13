###############################################################################
###############################################################################
## Purpose:    Stock analysis using gompertz model that includes space
## Author:     James T. Thorson
## Contact:    kellifayejohnson@gmail.com
## Date:       2014-05-10
###############################################################################
###############################################################################
jtt.base <- "c:/alaska/tmb"
setwd(jtt.base)
dataFldr <- file.path(jtt.base, "Data")
  DateFile = file.path(jtt.base, paste("sptlGmprtzSim", Sys.Date(), sep = "_"))
  dir.create(DateFile)

pkgs.need <- c("INLA", "TMB", "lme4", "maps", "mapdata", "mapproj")
pkgs.toinstl <- pkgs.need %in% installed.packages()
  sapply(pkgs.need[!pkgs.toinstl], install.packages, quiet = TRUE)
  sapply(pkgs.need, require, character.only = TRUE)
## mapdata == high-resolution world data
## mapproj == extra map data (e.g. gridlines)

SpeciesSet = c('NumSim')
  MeshType = c("Minimal","Recommended")[2] # Minimal: Refine=FALSE; Recommended: much slower!
  RunNew = TRUE

# Save record
RecordList = list("MeshType" = MeshType, "SpeciesSet" = SpeciesSet)
  capture.output(RecordList, 
                 file = file.path(DateFile, "RecordList.txt", sep=""))

# Loop across species
SpeciesI = 1
Species = SpeciesSet[SpeciesI]   # Rho = 0.98, 0.94, 0.53/0.29
Data <- read.csv(file = file.path(dataFldr, "DataSim.csv"), header = TRUE)

  tapply(Data[,Species], Data[,'Year'], function(Vec){sum(!is.na(Vec))})

  n_years    <- length(unique(Data$Year))
  n_stations <- length(unique(Data$Site))
  n_data     <- n_stations * n_years
  x_stations <- Data[match(unique(Data$Site),Data$Site),'Lon..DDD.DDDDD.']
  y_stations <- Data[match(unique(Data$Site),Data$Site),'Lat..DD.DDDDD.']
  area <- Data[match(unique(Data$Site),Data$Site),'Area']

  # Visualize distances
  plot(x_stations, y_stations, col = rainbow(length(unique(area)))[area],
       main = "Station locations w/ colors by area", las = 1)
    legend("topright", legend = unique(area), cex = 0.5, 
           pch = 1, col = rainbow(length(unique(area))))
    data.stations <- Data[!duplicated(Data$Area), ]
    text(data.stations$'Lon..DDD.DDDDD.', data.stations$'Lat..DD.DDDDD.',
         labels = data.stations$Area, font = 2, cex = 0.5, adj = 1,
         col = rainbow(length(unique(area)))[data.stations$Area])
  
# Compute distances between stations for every station
# Summary of distances between stations, 
# number of NAs should be equal to the number of stations
  Dist <- as.matrix(dist(cbind(x_stations, y_stations), upper = TRUE))
  summary(as.vector(ifelse(Dist==0, NA, Dist)) )
# Determine which stations are within a given area
# Summarize the distances of stations within an area
  Within <- outer(area, area, FUN=function(Num1,Num2){ifelse(Num1==Num2,1,NA)})
  summary( as.vector(ifelse((Dist*Within)==0,NA,Dist*Within)) )
  hist( log10(Dist) )
  abline(v = log10(0.1), lwd = 3)

# Determine the number of station observations per year for each station
# Retrieves the number of simulated fish at a given station in a given year
# and places them in a matrix with rows == station and column == year
  Ymat <- matrix(NA, nrow = n_stations, ncol = n_years)
  for(YearI in 1:n_years){
  for(SiteI in 1:n_stations){
    Which <- which(Data$Year == unique(Data$Year)[YearI] & 
                   Data$Site == unique(Data$Site)[SiteI])
    if(length(Which) >= 1) Ymat[SiteI, YearI] = Data[Which[1], Species]
  }}
  Y    <- as.vector(Ymat)		# Use Gaussian measurement error
  X    <- cbind(rep(1, n_data))
  Site <- as.vector(row(Ymat))
  Year <- as.vector(col(Ymat))
  #An indicator whethere the site was sampled or not.
  #1 if not sampled and 0 otherwise
  NAind = as.integer(ifelse(is.na(Y), 1, 0))
  # Inspect data
  mean(Y == 0, na.rm = TRUE)
  hist(log10(Y))         #
  (Glm = glm(Y ~ 0 + factor(Year) + factor(Site), family = "poisson"))

  ## Build SPDE object using INLA
  if(MeshType == "Recommended") {
    # loc_samp;  min.angle = 26,max.edge.data=0.08, max.edge.extra=0.2
    # defaults: min.angle = 21, max.edge = Inf
    my.refine <- list(min.angle = 26)
  } else{if(MeshType == "Minimal") {
    # loc_samp
    my.refine <- FALSE
  }}

  mesh <- inla.mesh.create(cbind(x_stations, y_stations), plot.delay = NULL,  
                           extend = list(n = 8, offset = -0.15),
                           refine = my.refine)  

  ## Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
  spde <- inla.spde2.matern(mesh, alpha = 2)
  plot(mesh)

  # Compile spatial
  if(FALSE){
    if(FALSE){
      dyn.unload(dynlib("spacetime_v7c"))
      file.remove(c("spacetime_v7c.o","spacetime_v7c.dll"))
    }
    compile("spacetime_v7c.cpp")
    # Compile nonspatial
    if(FALSE){
      dyn.unload(dynlib("time_v2b"))
      file.remove(c("time_v2b.o","time_v2b.dll"))
    }
    compile("time_v2c.cpp")
  }

  # Settings
  newtonOption(smartsearch = TRUE)
  ## NegBin specifies whether or not to use the negative
  ## binomial distribution in the cpp file
  NegBin = 1 

  # Run spatial version
  if(RunNew == FALSE & 
     (paste("Save_spatial_", Species, ".RData", sep = "") %in% list.files(DateFile))){
    if("Save_spatial" %in% search()) detach(Save_spatial)
    load(file = paste(DateFile, "/Save_spatial_", Species, ".RData", sep = ""))
    Match = which(names(Save_spatial) %in% ls())
    if(length(Match) >= 1) remove(list = names(Save_spatial)[Match])
    attach(Save_spatial)
  }else{
    # Run spatial model
    dyn.load(dynlib("TMB code/2014-02-09/spacetime_v7c"))
    Data_spatial = list(n_data = n_stations * n_years, Y = Y, NAind = NAind, 
                        NegBin = as.integer(NegBin), n_stations = n_stations, 
                        meshidxloc = mesh$idx$loc - 1, n_years = n_years, 
                        n_p = ncol(X), X = X, G0 = spde$param.inla$M0, 
                        G1 = spde$param.inla$M1, G2 = spde$param.inla$M2)
    if(FALSE) Parameters_spatial = list(alpha=c(0.0), phi=0.0, log_tau_E=0.0, log_tau_O=0.0, log_kappa=0.0,	rho=0.5, ln_VarInfl=c(0.0,0.0), Epsilon_input=matrix(0.0,nrow=spde$n.spde,ncol=n_years), Omega_input=rep(0.0,spde$n.spde))
    Parameters_spatial = list(alpha=c(0.0), phi=0.0, 
                              log_tau_E=0.0, log_tau_O=0.0, 
                              log_kappa=0.0,	rho=0.5, 
                              ln_VarInfl=c(0.0,0.0), 
                              Epsilon_input=matrix(rnorm(spde$n.spde*n_years),nrow=spde$n.spde,ncol=n_years), 
                              Omega_input=rnorm(spde$n.spde))
    if(FALSE) Parameters_spatial = list(alpha=c(0.0), phi=0.0, log_tau_E=0.0, log_tau_O=0.0, log_kappa=0.0,	rho=0.5, ln_VarInfl=c(0.0,0.0), Epsilon_input=matrix(rnorm(spde$n.spde*n_years),nrow=spde$n.spde,ncol=n_years), Omega_input=c(rowMeans(Epsilon_est),rnorm(spde$n.spde-nrow(Epsilon_est))))
    obj_spatial <- MakeADFun(data=Data_spatial, parameters=Parameters_spatial, random=c("Epsilon_input","Omega_input"))
    obj_spatial$fn(obj_spatial$par)
    # Run optimizer
    obj_spatial$control <- list(trace=1,parscale=rep(1,13),REPORT=1,reltol=1e-12,maxit=100)
    obj_spatial$hessian <- F
    opt_spatial = nlminb(obj_spatial$par, obj_spatial$fn, obj_spatial$gr, lower=c(rep(-20,2),rep(-10,3),-0.999,rep(-10,2)), upper=c(rep(20,2),rep(10,3),0.999,rep(10,2)), control=list(eval.max=1e4, iter.max=1e4))
    sdreport_spatial = try( sdreport(obj_spatial) )
    Report_spatial = obj_spatial$report()
    if( !("condition" %in% names(attributes(sdreport_spatial))) ) opt_spatial[["summary"]] = summary(sdreport_spatial)
    # spatial indices
    B_mean_spatial = opt_spatial$summary[which(rownames(opt_spatial$summary)=="mean_abundance"),'Estimate']
    B_conf_spatial = exp(opt_spatial$summary[which(rownames(opt_spatial$summary)=="log(mean_abundance)"),'Estimate']%o%rep(1,2) + opt_spatial$summary[which(rownames(opt_spatial$summary)=="log(mean_abundance)"),'Std. Error']%o%qnorm(c(0.1,0.9)))
    # Save for later
    capture.output( opt_spatial, file=paste(DateFile,"opt_spatial_",Species,".txt",sep=""))
    Save_spatial = list( "opt_spatial"=opt_spatial, "obj_spatial"=obj_spatial, "B_mean_spatial"=B_mean_spatial, "B_conf_spatial"=B_conf_spatial, "sdreport_spatial"=sdreport_spatial, "Report_spatial"=Report_spatial, "mesh"=mesh, "x_stations"=x_stations, "y_stations"=y_stations)
    save( Save_spatial, file=paste(DateFile,"Save_spatial_",Species,".RData",sep=""))
  }
  
  # Run nonspatial version
  if( RunNew==FALSE & ( paste("Save_nonspatial",Species,".RData",sep="") %in% list.files(DateFile)) ){
    if( "Save_nonspatial" %in% search()) detach( Save_nonspatial )
    load( file=paste(DateFile,"Save_nonspatial",Species,".RData",sep=""))
    Match = which( names(Save_nonspatial) %in% ls() )
    if( length(Match) >= 1) remove( list=names(Save_nonspatial)[Match] )
    attach( Save_nonspatial )
  }else{
    # Run non-spatial model
    dyn.load(dynlib("time_v2c"))
    Data_nonspatial = list(n_data=n_stations*n_years, Y=Y, NAind=NAind, NegBin=as.integer(NegBin), n_stations=n_stations, n_years=n_years, n_p=ncol(X), X=X)
    Parameters_nonspatial = list(alpha=c(0.0), phi=0.0, log_tau=0.0, rho=0.5, ln_VarInfl=c(0.0,0.0), Epsilon_input=rep(0,n_years))
    obj_nonspatial <- MakeADFun(data=Data_nonspatial, parameters=Parameters_nonspatial, random=c("Epsilon_input"))
    obj_nonspatial$fn(obj_nonspatial$par)
    # Run optimizer
    opt_nonspatial = nlminb(obj_nonspatial$par, obj_nonspatial$fn, obj_nonspatial$gr, lower=c(rep(-20,2),rep(-10,1),-0.999,rep(-10,2)),upper=c(rep(20,2),rep(10,1),0.999,rep(10,2)), control=list(eval.max=1e4, iter.max=1e4))
    sdreport_nonspatial = try( sdreport(obj_nonspatial) )
    Report_nonspatial = obj_nonspatial$report()
    if( !("condition" %in% names(attributes(opt_nonspatial))) ) opt_nonspatial[["summary"]] = summary(sdreport_nonspatial)
    # non-spatial indices
    B_mean_nonspatial = opt_nonspatial$summary[which(rownames(opt_nonspatial$summary)=="mean_abundance"),'Estimate']
    B_conf_nonspatial = exp(opt_nonspatial$summary[which(rownames(opt_nonspatial$summary)=="log(mean_abundance)"),'Estimate']%o%rep(1,2) + opt_nonspatial$summary[which(rownames(opt_nonspatial$summary)=="log(mean_abundance)"),'Std. Error']%o%qnorm(c(0.1,0.9)))
    # Save estimates
    capture.output( opt_nonspatial, file=paste(DateFile,"opt_nonspatial_",Species,".txt",sep=""))
    Save_nonspatial = list( "opt_nonspatial"=opt_nonspatial, "obj_nonspatial"=obj_nonspatial, "B_mean_nonspatial"=B_mean_nonspatial, "B_conf_nonspatial"=B_conf_nonspatial, "sdreport_nonspatial"=sdreport_nonspatial, "Report_nonspatial"=Report_nonspatial)
    save( Save_nonspatial, file=paste(DateFile,"Save_nonspatial",Species,".RData",sep=""))
  }
  
  # Range of correlation (Lindgren and Rue 2013, immediately before Eq. 4)
    Nu = 1
    sqrt(8*Nu)/exp(opt_spatial$par['log_kappa'])
  # Marginal Standard Deviation  (Lindgren and Rue 2013 Eq. 2)
    1 / sqrt(4*pi*exp(2*opt_spatial$par['log_tau_E'])*exp(2*opt_spatial$par['log_kappa']))
  # Marginal Standard Deviation  (Lindgren and Rue 2013 Eq. 2)
    1 / sqrt(4*pi*exp(2*opt_spatial$par['log_tau_O'])*exp(2*opt_spatial$par['log_kappa']))

  # Calculate esimates of random effects (RE)
  if(TRUE){
    y.lim = c( 32, 34.8)     # range(Data[,'Lat..DD.DDDDD.'])
    x.lim = c( -121, -117 )  # range(Data[,'Lon..DDD.DDDDD.'])
    redblue2 <- colorRampPalette(colors=c("blue","purple","red"))
    # Omega
    Omega_est = Report_spatial[['Omega']][mesh$idx$loc]
    png(file=paste(DateFile,"Omega_est_",Species,".png",sep=""), width=4, height=4, res=200, units="in")
      Rel = ((Omega_est[mesh$idx$loc]-min(Omega_est[mesh$idx$loc]))/diff(range(Omega_est[mesh$idx$loc])))
      map("worldHires", ylim=y.lim,xlim=x.lim,col="grey90",fill=T, main="", mar=c(0,0,2.5,0))
      #plot(x=x_stations, y=y_stations, cex=3*Rel)
      points(x=x_stations, y=y_stations, pch=20, col=redblue2(21)[round(21*Rel)])
     	box(bty="o",lwd=2)
    dev.off()
    # Equilibrium field
    Equil_est = Report_spatial[['Equil']][mesh$idx$loc]
    # Epsilon
    Epsilon_est = matrix( Report_spatial[['Epsilon']], ncol=n_years)[mesh$idx$loc,]
    Nrow = ceiling(sqrt(n_years)); Ncol=ceiling(n_years/Nrow)
    png(file=paste(DateFile,"Epsilon_est_",Species,".png",sep=""), width=2*Ncol, height=2*Nrow, res=200, units="in")
      par(mfrow=c(Nrow,Ncol), mar=c(2,2,1,0), mgp=c(1.5,0.5,0), tck=-0.02)
      for(i in 1:n_years){
        Rel = ((Epsilon_est[mesh$idx$loc,i]-min(Epsilon_est[mesh$idx$loc,]))/diff(range(Epsilon_est[mesh$idx$loc,])))
        map("worldHires", ylim=y.lim,xlim=x.lim,col="grey90",fill=T, main="", mar=c(0,0,2.5,0))
        points(x=x_stations, y=y_stations, pch=20, col=redblue2(21)[round(21*Rel)])
     	  box(bty="o",lwd=2)
      }
    dev.off()
    # lnB
    pow = function(a,b) a^b
    D_est = Report_spatial[['Dji']]
    Nrow = ceiling(sqrt(n_years)); Ncol=ceiling(n_years/Nrow)
    png(file=paste(DateFile,"D_est_",Species,".png",sep=""), width=2*Ncol, height=2*Nrow, res=200, units="in")
      par(mfrow=c(Nrow,Ncol), mar=c(2,2,1,0), mgp=c(1.5,0.5,0), tck=-0.02)
      for(i in 1:n_years){
        Rel = ((D_est[mesh$idx$loc,i]-min(D_est[mesh$idx$loc,]))/diff(range(D_est[mesh$idx$loc,])))
        #plot(x=x_stations, y=y_stations, cex=3*sqrt(Rel))
        map("worldHires", ylim=y.lim,xlim=x.lim,col="grey90",fill=T, main="", mar=c(0,0,2.5,0))
        points(x=x_stations, y=y_stations, pch=20, col=redblue2(21)[round(21*Rel)])
     	  box(bty="o",lwd=2)
      }
    dev.off()
    # Site specific abundance
    png(file=paste(DateFile,"Abundance_bysite_",Species,".png",sep=""), width=4, height=4, res=200, units="in")
      par(mar=c(3,3,2,0), mgp=c(1.5,0.5,0), tck=-0.02)
      matplot( t(D_est), type="l", col="black", lty="solid" )
    dev.off()
    # Compare my calculations and TMB
    plot( y=colMeans(exp(D_est)), x=B_mean_spatial)
    #rowMeans( exp(Report_spatial[['Dji']]) )
  }

  # Plotting abundance
  png( paste(DateFile,"Abundance_",Species,".png",sep=""), width=4, height=4, res=200, units="in")
    par(mfrow=c(1,1), mar=c(3,3,2,0), mgp=c(1.5,0.5,0), tck=-0.02)
    # Plotting
    plot( 1, type="n", xlab="Years", ylab="Abundance", lwd=3, ylim=range(range(B_conf_nonspatial),range(B_conf_spatial)), xlim=c(1,10), log="")
    polygon( x=c(1:nrow(B_conf_nonspatial),nrow(B_conf_nonspatial):1), y=c(B_conf_nonspatial[,1],rev(B_conf_nonspatial[,2])), lty="solid", col=rgb(1,0,0,0.2), border=NA)
    lines( B_mean_nonspatial, col="red", lwd=3)
    polygon( x=c(1:nrow(B_conf_nonspatial),nrow(B_conf_nonspatial):1), y=c(B_conf_spatial[,1],rev(B_conf_spatial[,2])), lty="solid", col=rgb(0,0,1,0.2), border=NA)
    lines( B_mean_spatial, col="blue", lwd=3)
  dev.off()
}

#################################
# Plot all together
#################################

#Data = read.csv( file=paste(dataFldr,"2014 H&L data request to Jim Thorson.csv",sep=""), header=TRUE)
Data <- read.csv(file = file.path(dataFldr, "DataSim.csv"), header = TRUE)
SpeciesSet <- c("NumSim")

# Plot together
png(file.path(DateFile,"Abundance_all.png"), width=1*3, height=2.5*length(SpeciesSet), res=200, units="in")
  par(mfrow=c(length(SpeciesSet),1), mar=c(2,3,0,0), mgp=c(1.5,0.5,0), tck=-0.02, oma=c(2,2,1.5,0), xaxs="i")
  for(SpeciesI in 1:length(SpeciesSet)){
    # Load old results
    Species = SpeciesSet[SpeciesI]   # Rho = 0.98, 0.94, 0.53/0.29
    load(file=paste(DateFile,"/Save_spatial_",Species,".RData",sep=""))
    #load(file=paste(DateFile,"/Save_nonspatial",Species,".RData",sep=""))
    # Plotting
    # Ylim = range(range(Save_nonspatial$B_conf_nonspatial),range(Save_spatial$B_conf_spatial))
    Ylim = range(Save_spatial$B_conf_spatial)
    plot(1, type="n", xlab="", ylab="", lwd=3, ylim=Ylim, xlim=2003+c(1,10), log="")
    text(x=2003+2, y=Ylim[2], pos=1, labels="rho = \nphi = \nAIC =")
    mtext(side=2, line=2, outer=FALSE, text=switch(SpeciesSet[SpeciesI], "NumVermC"="Vermilion", "NumGpot"="Greenspotted", "NumBoc"="Boccacio", "NumVerm"="Vermilion", "NumSunset"="Sunset"))
    # Spatial
    polygon( x=2003+c(1:10,10:1), y=c(Save_spatial$B_conf_spatial[,1],rev(Save_spatial$B_conf_spatial[,2])), lty="solid", col=rgb(1,0,0,0.2), border=NA)
    lines( x=2003+c(1:10), y=Save_spatial$B_mean_spatial, col="red", lwd=3)
    Which = match( c("rho","phi"), rownames(Save_spatial$opt_spatial$summary) )
    text(x=2003+4, y=Ylim[2], pos=1, col="red", labels=paste(formatC(round(Save_spatial$opt_spatial$summary[Which[1],'Estimate'],3),3)," (",formatC(round(Save_spatial$opt_spatial$summary[Which[1],'Std. Error'],3),3),")\n",formatC(round(Save_spatial$opt_spatial$summary[Which[2],'Estimate'],3),3)," (",formatC(round(Save_spatial$opt_spatial$summary[Which[2],'Std. Error'],3),3),")\n",formatC(round(2*Save_spatial$opt_spatial$objective+2*length(Save_spatial$opt_spatial),1),-1),sep=""))
    # non-spatial
    # polygon( x=2003+c(1:10,10:1), y=c(Save_nonspatial$B_conf_nonspatial[,1],rev(Save_nonspatial$B_conf_nonspatial[,2])), lty="solid", col=rgb(0,0,1,0.2), border=NA)
    # lines( x=2003+c(1:10), y=Save_nonspatial$B_mean_nonspatial, col="blue", lwd=3)
    # Which = match( c("rho","phi"), rownames(Save_nonspatial$opt_nonspatial$summary) )
    # text(x=2003+8, y=Ylim[2], pos=1, col="blue", labels=paste(formatC(round(Save_nonspatial$opt_nonspatial$summary[Which[1],'Estimate'],3),3)," (",formatC(round(Save_nonspatial$opt_nonspatial$summary[Which[1],'Std. Error'],3),3),")\n",formatC(round(Save_nonspatial$opt_nonspatial$summary[Which[2],'Estimate'],3),3)," (",formatC(round(Save_nonspatial$opt_nonspatial$summary[Which[2],'Std. Error'],3),3),")\n",formatC(round(2*Save_nonspatial$opt_nonspatial$objective+2*length(Save_spatial$opt_spatial),1),-1),sep=""))
    # GLM
    n_years = length(unique(Data$Year))
    n_stations = length(unique(Data$Site))
    n_data = n_stations*n_years
    Ymat = matrix(NA, nrow=n_stations, ncol=n_years)
    for(YearI in 1:n_years){
    for(SiteI in 1:n_stations){
      Which = which(Data$Year==unique(Data$Year)[YearI]&Data$Site==unique(Data$Site)[SiteI])
      if(length(Which)>=1) Ymat[SiteI,YearI] = Data[Which[1],Species]
    }}
    Y = as.vector(Ymat)		# Use Gaussian measurement error
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
y.lim = c( 32, 34.8)     # range(Data[,'Lat..DD.DDDDD.'])
x.lim = c( -121, -117 )  # range(Data[,'Lon..DDD.DDDDD.'])
Omega_range = c(-1.2,0.9)
redblue2 <- colorRampPalette(colors=(c("red","purple","blue","green","yellow")))
Seq = seq(-1.2,0.9,length=5)
png(file.path(DateFile, "Omega_all.png"), width = 1 * 3, 
    height = 2.5 * length(SpeciesSet), res = 200, units = "in")
  par(mfrow=c(length(SpeciesSet), 1), mar = c(0,0,0,0), mgp = c(1.5,0.5,0), 
      tck = -0.02, oma = c(3,3,1.5,0), xaxs = "i")
  for(SpeciesI in 1:length(SpeciesSet)){
    Species = SpeciesSet[SpeciesI]   # Rho = 0.98, 0.94, 0.53/0.29
    load(file = paste0(DateFile, "/Save_spatial_", Species, ".RData"))
    if(TRUE){ 
      Save_spatial[["x_stations"]] = x_stations
      Save_spatial[["y_stations"]] = y_stations
      Save_spatial[["mesh"]] = mesh
    } 
    Omega_est = Save_spatial$Report_spatial[['Omega']][Save_spatial$mesh$idx$loc]
    print(range(Omega_est))
    Rel = ((Omega_est[Save_spatial$mesh$idx$loc]-min(Omega_est[Save_spatial$mesh$idx$loc]))/diff(range(Omega_est[Save_spatial$mesh$idx$loc])))
    #Rel = (Omega_est - min(Omega_range)) / diff(Omega_range)
    map("worldHires", ylim=y.lim, xlim=x.lim, col="grey90", fill=T, mar=c(0,0,0,0), main="", myborder=0.01) # myborder is buffer around limits
    points(x=Save_spatial$x_stations, y=Save_spatial$y_stations, pch=20, 
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
    Species = SpeciesSet[SpeciesI]   # Rho = 0.98, 0.94, 0.53/0.29
    map("worldHires", ylim=y.lim, xlim=x.lim, col="grey90", fill=T, mar=c(0,0,0,0), main="", myborder=0.01) # myborder is buffer around limits
    #points(x=Save_spatial$x_stations, y=Save_spatial$y_stations, pch=20, col=redblue2(21)[round(21*Rel)])
    print( table( round(21*Rel) ))
	  if(SpeciesI==length(SpeciesSet)) axis(1)                                                                    # mar=c(0,0,2,3), 
    axis(2)
    box(bty="o",lwd=1)
    # GLM
    n_years = length(unique(Data$Year))
    n_stations = length(unique(Data$Site))
    n_data = n_stations*n_years
    Ymat = matrix(NA, nrow=n_stations, ncol=n_years)
    for(YearI in 1:n_years){
    for(SiteI in 1:n_stations){
      Which = which(Data$Year==unique(Data$Year)[YearI]&Data$Site==unique(Data$Site)[SiteI])
      if(length(Which)>=1) Ymat[SiteI,YearI] = Data[Which[1],Species]
    }}
    Y = as.vector(Ymat)		# Use Gaussian measurement error
    Site = as.vector(row(Ymat))
    Year = as.vector(col(Ymat))
    ( Glm = glm(Y ~ 0 + factor(Site) + factor(Year), family="poisson") )
    Match = match(unique(Data$Site),Data$Site)
    Rel = (Glm$coef[1:n_stations] - mean(Glm$coef[1:n_stations])) / diff(Omega_range)
    Rel = round(21*Rel)
    Rel = ifelse(Rel>21,21,Rel)
    Rel = ifelse(Rel<1,1,Rel)
    points( y=Data[Match,'Lat..DD.DDDDD.'], x=Data[Match,'Lon..DDD.DDDDD.'], pch=20, col=redblue2(21)[Rel] )  
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
    Species = SpeciesSet[SpeciesI]   # Rho = 0.98, 0.94, 0.53/0.29
    load(file=paste0(DateFile,"/Save_spatial_",Species,".RData"))
    #load( file=paste(DateFile,"Save_nonspatial",Species,".RData",sep=""))
    # Plotting
    # Ylim = range(range(Save_nonspatial$B_conf_nonspatial),range(Save_spatial$B_conf_spatial))
    Ylim = range(Save_spatial$B_conf_spatial)
    plot( 1, type="n", xlab="", ylab="", lwd=3, ylim=Ylim, xlim=2003+c(1,10), log="")
    text(x=2003+2, y=Ylim[2], pos=1, labels="rho = \nphi = \nAIC =")
    mtext(side=3, line=0, outer=FALSE, text=c("Simulated")[SpeciesI])
    # Spatial
    polygon( x=2003+c(1:10,10:1), y=c(Save_spatial$B_conf_spatial[,1],rev(Save_spatial$B_conf_spatial[,2])), lty="solid", col=rgb(1,0,0,0.2), border=NA)
    lines( x=2003+c(1:10), y=Save_spatial$B_mean_spatial, col="red", lwd=3)
    Which = match( c("rho","phi"), rownames(Save_spatial$opt_spatial$summary) )
    text(x=2003+4, y=Ylim[2], pos=1, col="red", labels=paste(formatC(round(Save_spatial$opt_spatial$summary[Which[1],'Estimate'],3),3)," (",formatC(round(Save_spatial$opt_spatial$summary[Which[1],'Std. Error'],3),3),")\n",formatC(round(Save_spatial$opt_spatial$summary[Which[2],'Estimate'],3),3)," (",formatC(round(Save_spatial$opt_spatial$summary[Which[2],'Std. Error'],3),3),")\n",formatC(round(2*Save_spatial$opt_spatial$objective+2*length(Save_spatial$opt_spatial),1),-1),sep=""))
    # non-spatial
    # polygon( x=2003+c(1:10,10:1), y=c(Save_nonspatial$B_conf_nonspatial[,1],rev(Save_nonspatial$B_conf_nonspatial[,2])), lty="solid", col=rgb(0,0,1,0.2), border=NA)
    # lines( x=2003+c(1:10), y=Save_nonspatial$B_mean_nonspatial, col="blue", lwd=3)
    # Which = match( c("rho","phi"), rownames(Save_nonspatial$opt_nonspatial$summary) )
    # text(x=2003+8, y=Ylim[2], pos=1, col="blue", labels=paste(formatC(round(Save_nonspatial$opt_nonspatial$summary[Which[1],'Estimate'],3),3)," (",formatC(round(Save_nonspatial$opt_nonspatial$summary[Which[1],'Std. Error'],3),3),")\n",formatC(round(Save_nonspatial$opt_nonspatial$summary[Which[2],'Estimate'],3),3)," (",formatC(round(Save_nonspatial$opt_nonspatial$summary[Which[2],'Std. Error'],3),3),")\n",formatC(round(2*Save_nonspatial$opt_nonspatial$objective+2*length(Save_spatial$opt_spatial),1),-1),sep=""))
    # GLM
    n_years = length(unique(Data$Year))
    n_stations = length(unique(Data$Site))
    n_data = n_stations*n_years
    Ymat = matrix(NA, nrow=n_stations, ncol=n_years)
    for(YearI in 1:n_years){
    for(SiteI in 1:n_stations){
      Which = which(Data$Year==unique(Data$Year)[YearI]&Data$Site==unique(Data$Site)[SiteI])
      if(length(Which)>=1) Ymat[SiteI,YearI] = Data[Which[1],Species]
    }}
    Y = as.vector(Ymat)		# Use Gaussian measurement error
    Site = as.vector(row(Ymat))
    Year = as.vector(col(Ymat))
    ( Glm = glm(Y ~ 0 + factor(Year) + factor(Site), family="poisson") )
    B_pred = predict(Glm, newdata=data.frame( "Year"=Year, "Site"=Site ), type="response")
    B_hat = tapply( B_pred, INDEX=Year, FUN=mean)
    lines( x=2003+c(1:10), y=B_hat, lwd=2, col="black", lty="dotted")
    mtext(side=1, line=1.3, outer=FALSE, text="Year")
  dev.off()
}
