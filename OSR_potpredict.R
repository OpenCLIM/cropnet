library(raster)
library(rgdal)
library(rgeos)
library(geosphere)
library(data.table)
library(ncdf4)
library(abind)

### UKCP18 NetCDF processing

# List NetCDF files
nlist <- list.files("\\\\nercwlsmb01\\data\\UKCP18\\RCM_12km\\daily_timeseries\\downloaded_data", ".nc", recursive = TRUE, full.names = T)

# List date range
dranges <- unique(substr(basename(nlist), nchar(basename(nlist))-19, nchar(basename(nlist))-3))

# List simulation runs
simr <- unique(substr(basename(nlist), nchar(basename(nlist))-26, nchar(basename(nlist))-25))

for(simx in simr[-1:-2]){
#simx <- "01"

# Set date range and subset list
for(dmet in dranges){
#dmet <- "19901201-20001130"
snlist <- nlist[grep(paste(simx,dmet, sep = "_day_"),  nlist)]

# Subset to required variables
vars <-  c("pr","tasmax","tasmin","rls" ,"rss","tas")
snlist <- snlist[sapply(strsplit(basename(snlist), "_"), function(x)x[1])  %in% vars]

# Read NetCDF and rearrange to required format
for(V in 1:length(vars)){
	# Get file index (helps with checking date ranges)
	vlist <- nlist[grep(paste(vars[V], "_", sep = ""),  nlist)]
	vlist <- vlist[grep(paste("_",simx,"_",sep = ""), vlist)]
	vind <- which(vlist == snlist[V])

	# Open netCDF
	ncin <- nc_open(vlist[vind])

	# Get coordiantes
	X <- ncvar_get(ncin, "projection_x_coordinate")
	Y <- ncvar_get(ncin, "projection_y_coordinate")
	
	# Get data
	tmp.array <- ncvar_get(ncin, vars[V])
	dlname <- ncatt_get(ncin, vars[V], "long_name")
	dunits <- ncatt_get(ncin, vars[V], "units")
	#fillvalue <- ncatt_get(ncin, vars[V], "_FillValue")
	#tmp.array[tmp.array == fillvalue$value] <- NA
	dim(tmp.array)

	# Get dates, years and growing years (based on benchmark sowing date of 01 September)
	dates <- gsub("  ", "",ncvar_get(ncin, "yyyymmdd"))
	groyears <- ifelse(as.numeric(substr(dates , 5,8)) >= as.numeric("0901"),   as.numeric(substr(dates ,1,4)) +1, as.numeric(substr(dates ,1,4)))

	if(vind-1 != 0){
		# Get previous date range to add start of first growing year
		ncprev <- nc_open(vlist[vind-1])
		prevdates <- gsub("  ", "",ncvar_get(ncprev, "yyyymmdd"))
		prevyears <- as.numeric(substr(prevdates,1,4))
		prevgroyears <- ifelse(as.numeric(substr(prevdates, 5,8)) >= as.numeric("0901"),   prevyears +1, prevyears)
	
		pvind <- which(prevgroyears == min(groyears))
		prev.array <- ncvar_get(ncprev, vars[V], start = c(1,1, min(pvind),1), count = c(-1,-1, length(pvind),-1))
		tmp.array <- abind(prev.array, tmp.array)
		dates <- c(prevdates[pvind], dates)
	}

	if(vind +1 <= length(vlist)){
		# Get next date range to add end of last growing year
		ncnext <- nc_open(vlist[vind+1])
		nxtdates <- gsub("  ", "",ncvar_get(ncnext, "yyyymmdd"))
		nxtyears <- as.numeric(substr(nxtdates,1,4))
		nxtgroyears <- ifelse(as.numeric(substr(nxtdates, 5,8)) >= as.numeric("0901"),   nxtyears +1, nxtyears)
	
		nvind <- which(nxtgroyears == max(groyears))
		nxt.array <- ncvar_get(ncnext, vars[V], start = c(1,1, min(nvind), 1), count = c(-1,-1, length(nvind),-1))
		tmp.array <- abind(tmp.array, nxt.array)
		dates <- c(dates, nxtdates[nvind])
	}

	# Calculate date and growing year for each time step
	allyears <-  as.numeric(substr(dates,1,4))
	allgroyears <- ifelse(as.numeric(substr(dates, 5,8)) >= as.numeric("0901"),   allyears +1, allyears)

	# Drop years without full year's data 
	byears <- names(table(allgroyears)[table(allgroyears) < 360])
	tmp.array <- tmp.array[,,!allgroyears %in% byears]
	dates <- dates[!allgroyears %in% byears]	
	allgroyears <- allgroyears[!allgroyears %in% byears]	

	# Check lengths
	dim(tmp.array)[3]
	length(dates)

	dimnames(tmp.array)[[1]] <- X
	dimnames(tmp.array)[[2]] <- Y
	dimnames(tmp.array)[[3]] <- dates
	#tmp.array[tmp.array == fillvalue$value] <- NA
	assign(paste(vars[V], "UKCP18", sep = "_"), tmp.array)
	save(list = paste(vars[V], "UKCP18", sep = "_"), file = paste("P:\\NEC07148_Crop-Net\\Workfiles\\Redhead crop modelling work\\UKCP18 data\\OSR\\",vars[V],"_UKCP18_",min(allgroyears),"_",max(allgroyears),"_",simx,".RData",sep = ""))
}
}	
}



### Can run from here

# Read soils data, reproject to grid and aggregate
	ukgrid <- "+init=epsg:27700"
	latlong <- "+init=epsg:4326"
AWCrast <- raster("N:\\ASSIST\\WP1 Yield data\\DEFRA Yield\\Landscape data\\WaterStorage\\MaxWet1.tif")
projection(AWCrast) <- projection(ukgrid)
AWCrastag <- aggregate(AWCrast, 10, fun=mean)

# Read CO2 data
pCO2 <- read.csv("P:\\NEC07148_Crop-Net\\Workfiles\\Redhead crop modelling work\\UKCP18_CO2_RCP85.csv")

# Choose growing year (demo to test growth model)
rsim <- "01"
gyear <- 2015

# Find and load required files and subset to year
rlist <- list.files("P:\\NEC07148_Crop-Net\\Workfiles\\Redhead crop modelling work\\UKCP18 data\\OSR", full.names = T)
dateranges <- t(sapply(unique(sapply(strsplit(rlist, "_"), function(x) paste(x[4], gsub(".RData","",x[5]), sep = "_"))), function(z){
	c(as.numeric(strsplit(z,"_")[[1]][1]),as.numeric(strsplit(z,"_")[[1]][2]))
}))
grange <- rownames(dateranges)[dateranges[,1] <= gyear & dateranges[,2] >= gyear][1]
gfiles <- rlist[grep(grange,rlist)]
gfiles <- gfiles[grep(paste("_",rsim,".RData",sep = ""), gfiles)]

findgyear <- function(xn, ggyear) {
	ifelse(as.numeric(substr(dimnames(xn)[[3]], 5,8)) >= as.numeric("0901"),   as.numeric(substr(dimnames(xn)[[3]],1,4)) +1, as.numeric(substr(dimnames(xn)[[3]],1,4))) == ggyear
}

lapply(gfiles, load, environment())
temp_tas <- tas_UKCP18[,,findgyear(tas_UKCP18, gyear)]
temp_rls <- rls_UKCP18[,,findgyear(rls_UKCP18, gyear)]
temp_rss <- rss_UKCP18[,,findgyear(rss_UKCP18, gyear)]
temp_pr <- pr_UKCP18[,,findgyear(pr_UKCP18, gyear)]
temp_tasmax <- tasmax_UKCP18[,,findgyear(tasmax_UKCP18, gyear)]
temp_tasmin <- tasmin_UKCP18[,,findgyear(tasmin_UKCP18, gyear)]
temp_cconc <- pCO2[pCO2$YEAR == gyear,2]

# Set growth model function
osr_py <- function(tmean, tmax, tmin, prec, solarrad, cconc = NULL,
	AWCraster = AWCrastag, Tbase = 4.5, GAItab = NULL, 
	HarvestJday = 212, k = 0.75, FCO2 = TRUE, savepath = NULL){

	# Defaults for testing
	# GAItab = NULL
	# tmean = temp_tas
	# tmax = temp_tasmax
	# tmin = temp_tasmin
	# prec = temp_pr 
	# solarrad = temp_rss	
	# AWCraster = AWCrastag
	# HarvestJday <- 212 # julian day of harvest
	# k <- 0.75 # Light extinciton coefficent
	# cconc = temp_cconc

	X <- as.numeric(dimnames(tmean)[[1]] )
      Y <- as.numeric(dimnames(tmean)[[2]] )

	### Degree of vernalisiton
	# Vernalisation rate (unit increase per day)
	TVmin <- -3.7
	TVmax <- 17.20
	TVop1	<- 0.72
	TVop2 <- 5.37
	RVmax <- 0.0145

	Vrate <- ifelse(tmean > TVmin, 
		ifelse(tmean <= TVop1, ((tmean - TVmin) / (TVop1 - TVmin)) * RVmax,
			ifelse(tmean <= TVop2, RVmax,
				ifelse(tmean < TVmax, (1-((tmean - TVop2) / (TVmax - TVop2))) * RVmax,
					0))), 0)
	Vf <- aperm(apply(Vrate, 1:2, cumsum), c(2,3,1)) 		
	Vf[Vf > 1] <- 1	

	### Photoperiod 
	# Convert OScoordinates to latitudes
	coords <- expand.grid(as.numeric(dimnames(tmean)[[1]]), as.numeric(dimnames(tmean)[[2]]))
	dat_SP <- SpatialPoints(coords, proj4string = CRS("+init=epsg:27700"))
	dat_SP_LL <- spTransform(dat_SP, CRS(latlong))
	lats <- matrix(coordinates(dat_SP_LL)[,2], nrow = dim(tmean)[1], ncol= dim(tmean)[2])

	# Julian days
	Jday <- as.POSIXlt(strptime(dimnames(tmean)[[3]], format = "%Y%m%d"))$yday + 1

	# Day length
	DL <- sapply(Jday, function(d) daylength(lats, d), simplify = "array")

	# Photoperiod factor
	Pbase <- 5.7
	Pmax <- 14.8
	Pf <- ifelse(DL > Pbase, ifelse(DL < Pmax, (DL-Pbase) / (Pmax - Pbase),1), 0) 
	Pf <- aperm(apply(Pf, 1:2, function(x){
		x[which(is.na(x))] <- (x[(min(which(is.na(x)))-1)] + x[(max(which(is.na(x)))+1)])/2
		return(x)
	}), c(2,3,1))

	### Effective daily temperature
	Tbase <- 4.5
	Teff <- ifelse(tmean > Tbase, tmean - Tbase, 0)

	# Growth response to daily temperature
	at <- c(7.6, 2.0, 5.1, 1.5)*(10^-3)

	### Daily development rate

	# Stage 0-1 (sowing to emrgence)
	devRate <- Teff * at[1]
	devStage <- aperm(apply(devRate, 1:2, cumsum), c(2,3,1)) 

	# Stage 1-2 (emergence to onset of flowering)
	devRate <- ifelse(devStage >= 1,  Teff * at[2] * Vf * Pf  ,devRate)
	devStage <- aperm(apply(devRate, 1:2, cumsum), c(2,3,1)) 

	# Stage 2-3 (flowering)
	devRate <- ifelse(devStage >= 2,  Teff * at[3]  ,devRate)
	devStage <- aperm(apply(devRate, 1:2, cumsum), c(2,3,1)) 

	# Stage 3-4 (end of flowering to maturity)
	devRate <- ifelse(devStage >= 3,  Teff * at[3]  ,devRate)
	devStage <- aperm(apply(devRate, 1:2, cumsum), c(2,3,1)) 
	devStage <- ifelse(devStage >= 4,  4 ,devStage)

	#plot(Pf[60,30,], type = "l", col = "orange", ylim = c(0,1))
	#points(Vf[60,30,], type = "l", col = "blue")
	#points(devRate[60,30,], type = "l", col = "green")
	#points(devStage[60,30,]/4, type = "l", col = "darkgreen")

	# Inital LAI
	LAIs <- 0.53
	RGRlai <- 0.00732
	LAI <- LAIs *  exp(RGRlai * aperm(apply(Teff, 1:2, cumsum), c(2,3,1)))
	LAI <- ifelse(LAI > 0.75, 2, LAI)

	# Residual light 
	Klai <- 0.75
	PAR <- temp_rss*(0.0036*24)/0.5
	REFg <- ifelse(devStage < 2, 0.046, ifelse(devStage > 3, 0.072, 0)) # Green canopy reflectance
	TR   <- ifelse(devStage > 2 & devStage <3,  1 - ((devStage - 2) * 0.5)   ,1)   # Flowering canopy transmission

	lt <- PAR * TR *(1-REFg) * exp(-Klai * LAI)
	lt <-ifelse(devStage > 3.5, PAR * (exp(-Klai * LAI) + 0.184*(devStage - 3.5)), lt)

	# Light absorbtion	
	la <- (PAR * TR *(1-REFg) ) - lt # MJ per m2

	# Crop growth rate
	RUE <- (ifelse(devStage < 2, 0.6, ifelse(devStage < 3, 1.4, 0.7)))

	# CO2 fertilisation effect
   	if(FCO2 == TRUE){
		RUE <-  ifelse(cconc < 350, RUE, ifelse(cconc > 750, RUE + RUE * (750/350 - 1) * 0.333, RUE + RUE * (cconc/350 - 1) * 0.333)) 
	}
	Wd <- la*RUE 

	#  Adjust LAI
	SLA <- 0.0177/10      # cm/g -> m/g
	LAI2 <- ifelse(LAI > 0.75,  aperm(apply(Wd * SLA, 1:2, cumsum), c(2,3,1)), LAI)
	LAI2[LAI > 0.75 & LAI2 < 0.75] <- 0.75
	LAI2[LAI2 > 4]  <- 4 
	maxLAI <- array(apply(LAI2, 1:2, max, na.rm = T), dim(LAI))
	LAI2[devStage > 3] <-  (maxLAI -((devStage - 3)/1)*(maxLAI-2.5))[devStage > 3]

	# Adjust light absorbtion
	lt2 <- PAR * TR *(1-REFg) * exp(-Klai * LAI)
	lt2 <- ifelse(devStage > 3.5, PAR * (exp(-Klai * LAI2) + 0.184*(devStage - 3.5)), lt2)
	la2 <- (PAR * TR *(1-REFg) ) - lt2 # MJ per m2

	# plot(PAR[60,30,], type = "l", col = "orange")
	# points(lt2[60,30,], type = "l", col = "red")
	# points(la2[60,30,], type = "l", col = "green")
	# abline(v = which(diff(cut(devStage[60,30,],c(0,1,2,3,4,5)-0.0001, labels = F), type= "l") == 1), lty = 2, col = "grey")
	# plot(LAI2[60,30,], type = "l", col = "darkgreen")

	# Crop dry matter biomass (g/m2) (water unlimited)
	WU_Biomass <- la2 * RUE

	###  Calculate availible water 

	# Extract soil AWC per grid cell 
	AWC <- array(extract(AWCraster,  dat_SP, fun = mean, na.rm = TRUE), dim(tmean))  

	# Assume saturation prior to Julian day 77
	Jarray <- array(rep(Jday,each = dim(tmean)[1]*dim(tmean)[2]), dim(tmean))
	Cday <- array(rep(1:dim(tmean)[3],each = dim(tmean)[1]*dim(tmean)[2]), dim(tmean))
	
	# Calcuate daily PAWC
	# Method 1 (apply/reduce)
	#PAWC <- array(0, dim(AWC))
	#PAWC[,,1] <- AWC[,,1]
	#PAWC[which((Jarray <= 77 | Cday < 150) & Cday != 1)] <- 0
	#Biouse <- WU_Biomass/ifelse(devStage < 3, 5, 3.8)
	#PAWC[which(Jarray > 77 & Cday > 150)] <- (prec-Biouse)[which(Jarray > 77 & Cday > 150)] 
	
	#PAWC <- aperm(apply(PAWC, 1:2, function(x){
	#	maxX <- x[1]
	#	Reduce(function(u, v){
	#		pmax(pmin(u + v, maxX), 0, na.rm = T)
	#	}, x, accumulate = TRUE)
	#}), c(2,3,1), init = 0)
	#plot(PAWC[60,30,], type = "l")
	
	# Method 2 (for loop)
	PAWC <- array(0, dim(AWC))
	PAWC[which((Jarray <= 77 | Cday < 150))] <- AWC[which((Jarray <= 77 | Cday < 150))]
	PAWC[which(is.na(Jarray) == T)] <- AWC[which(is.na(Jarray) == T)]

	for(i in 1:dim(PAWC)[3]){
		ji <- which(Jarray[,,i] > 77 & Cday[,,i] > 150)
		PAWC[,,i][ji] <- pmin(AWC, pmax((PAWC[,,max(1,i-1)] + prec[,,i]) - (WU_Biomass[,,max(1,i-1)]/ifelse(devStage < 3, 5, 3.8)[,,max(1,i-1)]), 0))[ji]
	}
	#points(PAWC[60,30,], type = "l", col = "red")
	
	#Calculate water-limited biomass
	WL_Biomass <- ifelse(PAWC>0, WU_Biomass, 0)

	# Consecutive days of waterlogging
	# For each waterlogging event beyond 14 days, lose 50% biomass if pre-stem elongation... 
	wlog1x <- array(0,dim(PAWC))
	wlog1x[which(Jarray <= HarvestJday & (PAWC+prec) > (AWC + 10) & devStage < 2)] <- 1
	wlog1x <- aperm(apply(wlog1x, 1:2, function(x){
		ifelse(x == 1, unlist(sapply(rle(x)$lengths, seq)), 0)
	}), c(2,3,1))

	pwloss1 <- ifelse(wlog1x > 5 & wlog1x <= 14, (wlog1x-5)/(14-5)*0.5, ifelse(wlog1x > 14,0.5,0))
	WLWL_Biomass <- WL_Biomass*(1-pwloss1)

	# ...15% if flowering (Wollmer et al 2017)
	wlog2x <- array(0,dim(PAWC))
	wlog2x[which(Jarray <= HarvestJday & (PAWC+prec) > (AWC + 10) & devStage >= 2)] <- 1
	wlog2x <- aperm(apply(wlog2x, 1:2, function(x){
		ifelse(x == 1, unlist(sapply(rle(x)$lengths, seq)), 0)
	}), c(2,3,1))

	pwloss2 <- ifelse(wlog1x > 5 & wlog1x <= 14, (wlog1x-5)/(14-5)*0.15, ifelse(wlog1x > 14,0.5,0))
	WLWL_Biomass <- WLWL_Biomass*(1-pwloss2)

	# plot(PAWC[60,30,], type ="l", col = "blue", ylim = c(0,450))
	# points(prec[60,30,], type ="l", col = "lightblue")
	# plot(wlog1x[60,30,], type ="l", col = "lightblue")
	# plot(wlog2x[60,30,], type ="l", col = "lightblue")
	# plot(WU_Biomass[60,30,], type ="l", col ="red")
	# points(WL_Biomass[60,30,], type ="l", col = "blue")
	# points(WLWL_Biomass[60,30,], type ="l", col = "purple")

	### Partitioning water limited biomass to seed via harvest index 
	WLWL_Biomass[Jarray > HarvestJday & Cday > 150] <- 0
	WLyield <- apply(WLWL_Biomass,1:2, sum) * 0.35  #(g/m2)
	WLYield <- (WLyield/1000000)*10000  #(t/Ha2)
	#hist(WLYield)

	### Alternative apporach assuming only biomass post-flowering availible for yield
	AWLWL_Biomass <- ifelse(devStage > 3,  WLWL_Biomass, 0.1 * WLWL_Biomass)
	# points(AWLWL_Biomass[60,30,], type ="l", col = "green")
	WLYield2 <- (apply(AWLWL_Biomass, 1:2, sum)/1000000)*10000
	#hist(WLYield2)

	# par(mfrow = c(2,2))
	# brk <- seq(1000, 8000, 500)
	# clr <- colorRampPalette(c("pink", "yellow", "green", "darkgreen"))(14)
	#image(X,Y, apply(WU_Biomass,1:2,sum), asp = 1, breaks = brk , col = clr) 
	#image(X,Y, apply(WL_Biomass,1:2,sum), asp = 1, breaks = brk , col = clr) 
	#image(X,Y, apply(WLWL_Biomass,1:2,sum), asp = 1, breaks = brk , col = clr) 
	#image(X,Y, WLYield, asp = 1, breaks = 1:15, col = clr) 
	#image(X,Y, WLYield2, asp = 1, breaks = 1:15, col = clr) 

	# Convert final yield to raster
	WLPY=list()
	WLPY$x= X
 	WLPY$y= Y
 	WLPY$z= WLYield2
	WLPYrast <- raster(WLPY, crs = ukgrid)

	if(is.null(savepath) == FALSE){
		writeRaster(WLPYrast, savepath,"GTiff", overwrite = TRUE)
	} 	
	return(WLPYrast)
}

# Run for each 
yrs <- 1982:2080
sms <- c("01","04","05","06","07","08","09","10","11","12","13","15")

allc <- expand.grid(yrs, sms)

findgyear <- function(xn, ggyear) {
	ifelse(as.numeric(substr(dimnames(xn)[[3]], 5,8)) >= as.numeric("0901"),   as.numeric(substr(dimnames(xn)[[3]],1,4)) +1, as.numeric(substr(dimnames(xn)[[3]],1,4))) == ggyear
}

for(XC in 1077:nrow(allc)){
	# Find and load required files and subset to year
	grange <- rownames(dateranges)[dateranges[,1] <= allc[XC,1] & dateranges[,2] >= allc[XC,1]][1]
	gfiles <- rlist[grep(grange,rlist)]
	gfiles <- gfiles[grep(paste("_", allc[XC,2], sep = ""), gfiles)]

	lapply(gfiles, load, environment())
	temp_tas <- tas_UKCP18[,,findgyear(tas_UKCP18, allc[XC,1])]
	temp_rls <- rls_UKCP18[,,findgyear(rls_UKCP18, allc[XC,1])]
	temp_rss <- rss_UKCP18[,,findgyear(rss_UKCP18, allc[XC,1])]
	temp_pr <- pr_UKCP18[,,findgyear(pr_UKCP18, allc[XC,1])]
	temp_tasmax <- tasmax_UKCP18[,,findgyear(tasmax_UKCP18, allc[XC,1])]
	temp_tasmin <- tasmin_UKCP18[,,findgyear(tasmin_UKCP18, allc[XC,1])]
	temp_cconc <- pCO2[pCO2$YEAR == allc[XC,1],2]

	osr_py(temp_tas,temp_tasmax, temp_tasmin, temp_pr, temp_rss, cconc = temp_cconc, FCO2 = TRUE, savepath = paste("P:\\NEC07148_Crop-Net\\Workfiles\\Redhead crop modelling work\\UKCP18 Yield\\COSR2", allc[XC,1],allc[XC,2], sep = "_"))
}





#### Compare with suveyed yield data

# Individual hectads
yw <- readOGR("N:\\ASSIST\\WP1 Yield data\\DEFRA Yield\\Analysis outputs", "YieldStats_Update_winter oilseed")
yw$Mean_Yield

# Load yield trends from survey data
load("N:\\ASSIST\\WP1 Yield data\\DEFRA Yield\\Analysis outputs\\CollatedYieldData_withPre2010.R")
sag <- aggregate(coll_ydat$Yield[coll_ydat$Crop == "winter oilseed"], list(coll_ydat$Year[coll_ydat$Crop == "winter oilseed"]), mean)

# Modelled potential yield for surveyd years
ppat <- "P:\\NEC07148_Crop-Net\\Workfiles\\Redhead crop modelling work\\UKCP18 Yield\\OSR2"
ystack <- stack(paste(ppat,  "\\OSR_2010_01.tif", sep = ""),
	paste(ppat,  "\\OSR_2011_01.tif", sep = ""),
	paste(ppat,  "\\OSR_2012_01.tif", sep = ""),
	paste(ppat,  "\\OSR_2013_01.tif", sep = ""),
	paste(ppat,  "\\OSR_2014_01.tif", sep = ""),
	paste(ppat,  "\\OSR_2015_01.tif", sep = ""),
	paste(ppat,  "\\OSR_2016_01.tif", sep = ""),
	paste(ppat,  "\\OSR_2017_01.tif", sep = ""))

# Averages per year (whole dataset)
natav <- cellStats(ystack, mean)

# Averages across surveyed cells and per cell
survav <- extract(ystack, yw, fun = mean, na.rm= TRUE)
spy <- apply(survav,2,mean)
spc <- apply(survav,1,mean)

# Compare with plots
par(mfrow = c(2,2), las = 1, mar = c(4,4,2,1), cex = 1.1)

plot(yw$Mean_Yield, spc, xlab = "Mean yield over time per hectad \n(surveyed)", ylab = "Mean yield over time \n(UKCP18 modelled)")
plot(sag$x, spy, xlab = "Mean yield per year \n(surveyed)", ylab = "Mean yield per year \n(UKCP18 modelled)", xlim = c(2,5), ylim = c(12, 14))
text(sag$x, spy, 2008:2017, col = "red", pos= 3)
cor.test(sag$x, spy)





























