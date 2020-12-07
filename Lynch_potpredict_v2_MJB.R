library(raster)
library(rgdal)
library(rgeos)
library(geosphere)
library(data.table)
library(ncdf4)
library(abind)

formatdrivingdatawheat <- function(datadir, outdir, xname, yname){
### UKCP18 NetCDF processing

# List NetCDF files
nlist <- list.files(datadir, ".nc", recursive = TRUE, full.names = T)

# List date range
dranges <- unique(substr(basename(nlist), nchar(basename(nlist))-19, nchar(basename(nlist))-3))

# List simulation runs
simr <- unique(substr(basename(nlist), nchar(basename(nlist))-26, nchar(basename(nlist))-25))

for(simx in simr){
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
	X <- ncvar_get(ncin, xname)
	Y <- ncvar_get(ncin, yname)	
	# Get data
	tmp.array <- ncvar_get(ncin, vars[V])
	dlname <- ncatt_get(ncin, vars[V], "long_name")
	dunits <- ncatt_get(ncin, vars[V], "units")
	#fillvalue <- ncatt_get(ncin, vars[V], "_FillValue")
	#tmp.array[tmp.array == fillvalue$value] <- NA
	dim(tmp.array)

	# Get dates, years and growing years (based on AHDB benchmark sowing date of 01 October)
	dates <- gsub("  ", "",ncvar_get(ncin, "yyyymmdd"))
	groyears <- ifelse(as.numeric(substr(dates , 5,8)) >= as.numeric("1001"),   as.numeric(substr(dates ,1,4)) +1, as.numeric(substr(dates ,1,4)))

	if(vind-1 != 0){
		# Get previous date range to add start of first growing year
		ncprev <- nc_open(vlist[vind-1])
		prevdates <- gsub("  ", "",ncvar_get(ncprev, "yyyymmdd"))
		prevyears <- as.numeric(substr(prevdates,1,4))
		prevgroyears <- ifelse(as.numeric(substr(prevdates, 5,8)) >= as.numeric("1001"),   prevyears +1, prevyears)
	
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
		nxtgroyears <- ifelse(as.numeric(substr(nxtdates, 5,8)) >= as.numeric("1001"),   nxtyears +1, nxtyears)
	
		nvind <- which(nxtgroyears == max(groyears))
		nxt.array <- ncvar_get(ncnext, vars[V], start = c(1,1, min(nvind), 1), count = c(-1,-1, length(nvind),-1))
		tmp.array <- abind(tmp.array, nxt.array)
		dates <- c(dates, nxtdates[nvind])
	}

	# Calculate date and growing year for each time step
	allyears <-  as.numeric(substr(dates,1,4))
	allgroyears <- ifelse(as.numeric(substr(dates, 5,8)) >= as.numeric("1001"),   allyears +1, allyears)

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
	save(list = paste(vars[V], "UKCP18", sep = "_"), file = paste(outdir,vars[V],"_UKCP18_",min(allgroyears),"_",max(allgroyears),"_",simx,".RData",sep = ""))
}
}	
}
}


loaddata <- function(rsim, gyear, Rdatpath, CO2file){

# Read CO2 data
pCO2 <- read.csv(CO2file)

# Choose growing year (demo to test growth model)
#rsim <- "11"
#gyear <- 2015

# Find and load required files and subset to year
rlist <- list.files(Rdatpath, full.names = T)
dateranges <- t(sapply(unique(sapply(strsplit(rlist, "_"), function(x) paste(x[4], gsub(".RData","",x[5]), sep = "_"))), function(z){
	c(as.numeric(strsplit(z,"_")[[1]][1]),as.numeric(strsplit(z,"_")[[1]][2]))
}))
grange <- rownames(dateranges)[dateranges[,1] <= gyear & dateranges[,2] >= gyear][1]
gfiles <- rlist[grep(grange,rlist)]
gfiles <- gfiles[grep(paste("_",rsim,".RData",sep = ""), gfiles)]

findgyear <- function(xn, ggyear) {ifelse(as.numeric(substr(dimnames(xn)[[3]], 5,8)) >= as.numeric("1001"),
                                          as.numeric(substr(dimnames(xn)[[3]],1,4)) +1, as.numeric(substr(dimnames(xn)[[3]],1,4))) == ggyear}

lapply(gfiles, load, environment())
temp_tas <- tas_UKCP18[,,findgyear(tas_UKCP18, gyear)]
#temp_rls <- rls_UKCP18[,,findgyear(rls_UKCP18, gyear)]
temp_rss <- rss_UKCP18[,,findgyear(rss_UKCP18, gyear)]
temp_pr <- pr_UKCP18[,,findgyear(pr_UKCP18, gyear)]
temp_tasmax <- tasmax_UKCP18[,,findgyear(tasmax_UKCP18, gyear)]
temp_tasmin <- tasmin_UKCP18[,,findgyear(tasmin_UKCP18, gyear)]
temp_cconc <- pCO2[pCO2$YEAR == gyear,2]

coords <- expand.grid(as.numeric(dimnames(temp_tas)[[1]]), as.numeric(dimnames(temp_tas)[[2]]))
dat_SP <- SpatialPoints(coords, proj4string = CRS("+init=epsg:27700"))
# Read soils data, reproject and aggregate to coarser resolution
ukgrid <- "+init=epsg:27700"
latlong <- "+init=epsg:4326"
AWCrast <- raster("/users/sgsys/matbro/cropNET/data/MaxWet1.tif")
#print(dim(AWCrast))
projection(AWCrast) <- projection(ukgrid)
AWCraster <- aggregate(AWCrast, 10, fun=mean)
#print(AWCraster)
#print(dim(AWCraster))
# Extract soil AWC per grid cell - constant over time. 
# Below function regrids the AWC data on to the model grid
AWC <- array(extract(AWCraster,  dat_SP, fun = mean, na.rm = TRUE), dim(temp_tas))
#print(AWC)
#print(dim(dat_SP))
#print(dim(AWC))

X <- as.numeric(dimnames(temp_tas)[[1]] )
Y <- as.numeric(dimnames(temp_tas)[[2]] )
    
data <- list('temp_tas'=temp_tas, 'temp_rss'=temp_rss, 'temp_pr'=temp_pr, 'temp_tasmax'=temp_tasmax,'temp_tasmin'= temp_tasmin, 'temp_cconc'=temp_cconc, 'AWC'=AWC, 'X'=X, 'Y'=Y)
return(data)
}

# Set growth model function
GAI <- function(tmean, tmax, tmin, prec, solarrad, Dt = 14, Tbase = 0, GAItab = NULL, HarvestJday = 243){

	# Defaults for testing for this and wheat_yield function
	#GAItab = NULL
	#tmean = temp_tas
	#tmax = temp_tasmax
	#tmin = temp_tasmin
	#prec = temp_pr 
	#solarrad = temp_rss	
	#AWCraster = AWCrastag
	#Dt <- 14 # Day length threshold for onset of GS31 (construction phase)
	#Tbase <- 0 # Threshold for degree days
	#HarvestJday <- 243 # julian day of harvest
	#RUE <- 3.1 # Radiation use efficiency
	#WLP <- 117.7 # Waterlogging loss penalty (kg/ha yield loss per day of waterlogging) (Olgun et al)
	#cconc = temp_cconc
        #FCO2 = TRUE

        latlong <- "+init=epsg:4326"

	X <- as.numeric(dimnames(tmean)[[1]] )
        Y <- as.numeric(dimnames(tmean)[[2]] )
        #print(X)
        #print(Y)

	# Set up thermal time to growth stages and associated GAI (if not specified)
	if(is.null(GAItab)){
		
		GS <- c("Sowing-GS30","GS30-GS31","GS31-GS61","GS61-GS69","GS69-GS87","GS87-senescense")
		TT <- c(1100,100,900,50,750,200) # thresholds of cumulative degree days
		GAIstart <- c(0, 1.6, 2.0, 6.3, 6.3, 1.3)
		GAIend <- c(1.6, 2.0, 6.3, 6.3, 0.7, 0)	
	}

        # Calculate day length
        dates <- paste(substr(dimnames(tmean)[[3]],1,4),substr(dimnames(tmean)[[3]],5,6), substr(dimnames(tmean)[[3]],7,8),sep = "-")
        #print(dimnames(tmean)[[3]])

	# Convert OScoordinates to latitudes
	coords <- expand.grid(as.numeric(dimnames(tmean)[[1]]), as.numeric(dimnames(tmean)[[2]]))
	dat_SP <- SpatialPoints(coords, proj4string = CRS("+init=epsg:27700"))
	dat_SP_LL <- spTransform(dat_SP, CRS(latlong))
        lats <- matrix(coordinates(dat_SP_LL)[,2], nrow = dim(tmean)[1], ncol= dim(tmean)[2])
        #print(dim(lats))

	# Julian days
        Jday <- as.POSIXlt(strptime(dimnames(tmean)[[3]], format = "%Y%m%d"))$yday + 1
        #print(Jday)
	# Days since sowing
        Cday <- array(rep(1:dim(tmean)[3],each = dim(tmean)[1]*dim(tmean)[2]), dim(tmean))
        #print(Cday[1,1,1:360])

	# Day length
	DL <- sapply(Jday, function(d) daylength(lats, d), simplify = "array")

 	# Day on which day length > threshold
	DL14 <- apply(DL, 1:2, function(DX) min(which(Jday > 1 & DX >= Dt) ))
        #print(dim(DL14))
	# Degree days
        Tbase <- 0
        #print(tmin[50,50,1:360])
        #print(tmax[50,50,1:360])
	DD <-  ifelse(tmin > Tbase, ((tmax + tmin)/2)-Tbase,
				ifelse(tmax < Tbase, 0,
					ifelse(tmin <= Tbase, ((tmax-Tbase)/2) - ((Tbase - tmin)/4), 999)))
        #print(DD[50,50,1:360])
	# Vernalisation days
	VD <- ifelse(tmean > 3 & tmean  <= 10, 1,
			ifelse(tmean  > -4 & tmean  <= 3, (tmean +4)/7,
				ifelse(tmean  > 10 & tmean  <= 17, (17-tmean )/7,
					ifelse(tmean  < -4, 0, 
						ifelse(tmean > 17, 0, 999)))))
        #print(dim(VD))
	# Set starting values for GAI and growth stage arrays
	GAI <- array(0, dim(tmean))
	GSS <- array("X", dim(tmean))

        ### GAI for Growth stages 0 - GS30
        # for each latlon point, do a cumulative sum over the time dimension
        # for some reason this rearranges the dimensions, so the aperm(..., c(2,3,1)) bit rearranges them back
	CDD <- aperm(apply(DD, 1:2, cumsum), c(2,3,1)) 
	p1i <- which(CDD <= TT[1])
	GAI[p1i] <- ((CDD/TT[1])* GAIend[1])[p1i]
        GSS[p1i] <- "Sowing-GS30"
        #print(CDD[50,50,1:360])
        #print(GAI[50,50,1:360])
        #print(GSS[50,50,1:360])

        ### GAI for Growth stages GS30 - GS31
	p2i <- which(CDD > TT[1] & CDD <= (TT[1] + TT[2]))
	DD2 <- DD
	DD2[-p2i] <- 0
	CDD2 <- aperm(apply(DD2, 1:2, cumsum), c(2,3,1))

	GAI[p2i] <- ((CDD2/TT[2]) * (GAIend[2] - GAIstart[2]) + GAIstart[2])[p2i]
	GSS[p2i] <- "GS30-GS31"

	### GAI for Growth stages GS31-GS61
	p3i <- which(CDD > (TT[1] + TT[2]) & CDD <= (TT[1] + TT[2] + TT[3]))
	
	# Check vernalisaiton
	GAI[p3i][aperm(apply(VD, 1:2, cumsum), c(2,3,1))[p3i] < 50] <- GAIstart[3]

	# Check day length
	GAI[p3i][Cday[p3i] < DL14[p3i]] <- GAIstart[3]

	# Reduce max GAI for GS61 if failed to accumulate sufficent degree days by DL14
	mGAI61 <- CDD
	mGAI61[Cday < array(DL14,dim(CDD))] <- NA
	mGAI61 <- apply(mGAI61, 1:2, min, na.rm = T)
	mGAI61 <- apply((1- ((TT[1] + TT[2]) - mGAI61)/(TT[1] + TT[2])) * GAIend[3], 1:2, function(x) min(x,GAIend[3]))

	# Set degree days in GS31-61 before adequate vernalisation or DL14 to zero
	DD[CDD > (TT[1] + TT[2]) & (aperm(apply(VD, 1:2, cumsum), c(2,3,1)) < 50)] <- 0
	DD[CDD > (TT[1] + TT[2]) & (Cday < array(DL14,dim(CDD)))] <- 0

 	# Calculate GAI for following days
	CDD <- aperm(apply(DD, 1:2, cumsum), c(2,3,1)) 
	p4i <- which(CDD > (TT[1] + TT[2]) & CDD <= (TT[1] + TT[2] + TT[3])) 
	DD4 <- DD
	DD4[-p4i] <- 0
	CDD4 <- aperm(apply(DD4, 1:2, cumsum), c(2,3,1))
	GAI[setdiff(p3i, p4i)] <- GAIstart[3]
	GAI[p4i] <- (((CDD4/TT[3]) * (array(mGAI61, dim(CDD)) - GAIstart[3])) + GAIstart[3])[p4i]
	GSS[setdiff(p3i, p4i)] <- "GS31-GS61"
	GSS[p4i] <- "GS31-GS61"

	# Freeze GAI if temperature drops below -5 after GS31
	if(length(which(tmin <= -5 & GSS == "GS31-GS61")) >0){
		b5i <- which(tmin <= -5 & GSS == "GS31-GS61") 
		B5 <- GAI
		B5[-b5i] <- NA
		B5GAI <- array(apply(B5,1:2, min, na.rm = T), dim(GAI))
		GAI[b5i] <- B5GAI[b5i]
	}
	
	### GAI for Growth stages GS61-GS69
	p5i <- which(CDD > (TT[1] + TT[2]+ TT[3]) & CDD <= (TT[1] + TT[2] + TT[3] + TT[4]))
	GAI[p5i] <- array(apply(GAI,1:2, max), dim(GAI))[p5i]
	GSS[p5i] <- "GS61-GS69"

	### GAI for Growth stages GS69 - 87
	p6i <- which(CDD > (TT[1] + TT[2]+ TT[3] +TT[4]) & CDD <= (TT[1] + TT[2] + TT[3] + TT[4] + TT[5]))
	DD6 <- DD
	DD6[-p6i] <- 0
	CDD6 <- aperm(apply(DD6, 1:2, cumsum), c(2,3,1))
	GAI[p6i] <- ((CDD6/TT[5]) * (GAIend[5] - array(apply(GAI,1:2, max, na.rm = T), dim(GAI))) + array(apply(GAI,1:2, max, na.rm = T), dim(GAI)))[p6i] 
	GSS[p6i] <- "GS69-GS87"

	### GAI for Growth stages GS87 - senescense
	p7i <- which(CDD > (TT[1] + TT[2]+ TT[3] +TT[4] + TT[5]) & CDD <= (TT[1] + TT[2] + TT[3] + TT[4] +TT[5]  + TT[6]))
	DD7 <- DD
	DD7[-p7i] <- 0
	CDD7 <- aperm(apply(DD7, 1:2, cumsum), c(2,3,1))
	GAI6 <- GAI
	GAI6[-p6i] <- NA
	GAI6 <- array(apply(GAI6, 1:2, min, na.rm = T), dim(GAI))
	GAI[p7i] <- ((CDD7/TT[6]) * (GAIend[6] - GAI6 ) + GAI6)[p7i] 
	GSS[p7i] <- "GS87-senescense"


#plot(Cday[50,40,], GAI[50,40,], type= "o", col = as.factor(GSS[50,40,]))
#abline(v = 330)

	### Account for failure to senesce pre-harvest (ANNOTATE THIS SECTION)
	Jarray <- array(rep(Jday,each = dim(tmean)[1]*dim(tmean)[2]), dim(tmean))

	# Find day since sowing for harvest
	harvestCDD <- CDD          
	harvestCDD[Jarray > HarvestJday & Cday > 100] <- NA
	harvestCDD <- array(apply(harvestCDD, 1:2, max, na.rm = T),dim(CDD))

	# Find julian day of maximum GAI
	maxGAI <- array(apply(GAI, 1:2, max),dim(GAI))	
	MJarray <- Jarray
	MJarray[GAI != maxGAI] <- NA
	MJarray <- array(apply(MJarray, 1:2, max, na.rm = T),dim(Jarray))

	# Porptional progress from julian day of maximum GAI to harvest
	d <- Jarray - MJarray
	d[Cday < 100 | d < 0] <- 0
	d[Jarray > HarvestJday] <- 0
	d <- d/array(apply(d,1:2,max, na.rm = T),dim(d))

	p8i <- which(harvestCDD < (TT[1] + TT[2] + TT[3] + TT[4] + TT[5] + TT[6]) & Jarray <= HarvestJday & Jarray > MJarray)
	GAI[p8i] <- ((d * (GAIend[6] - maxGAI)) + maxGAI)[p8i]

	GAI[Jarray > HarvestJday & Cday > 100] <- 0

	### Account for failure to flower pre-harvest!
	GAI[array(apply(CDD,1:2, max),dim(CDD)) < (TT[1] + TT[2] + TT[3] + TT[4])] <- 0

        
        # split function here
                                        #points(Cday[50,40,], GAI[50,40,], type= "l", col = "green")
        FCO2 <- TRUE
        #print(X)
        #print(Y)

        print(GAI[50,50,1:360])
        data <- list('GAI'=GAI, 'tmean'=tmean, 'prec'=prec, 'solarrad'=solarrad, 'dat_SP'=dat_SP, 'Jarray'=Jarray, 'Cday'=Cday, 'GSS'=GSS, 'HarvestJday'=HarvestJday, 'x'=X, 'y'=Y, 't'=dates)
      	return(data)
      }
        # AFAIK GAI=LAI
	### Converting GAI to biomass
                                        # % PAR intercepted by canopy


GAI_point<- function(tmean, tmax, tmin, prec, solarrad, X, Y, lat, T, Dt = 14, Tbase = 0, GAItab = NULL, HarvestJday = 243){

	# Defaults for testing for this and wheat_yield function
	#GAItab = NULL
	#tmean = temp_tas
	#tmax = temp_tasmax
	#tmin = temp_tasmin
	#prec = temp_pr 
	#solarrad = temp_rss	
	#AWCraster = AWCrastag
	#Dt <- 14 # Day length threshold for onset of GS31 (construction phase)
	#Tbase <- 0 # Threshold for degree days
	#HarvestJday <- 243 # julian day of harvest
	#RUE <- 3.1 # Radiation use efficiency
	#WLP <- 117.7 # Waterlogging loss penalty (kg/ha yield loss per day of waterlogging) (Olgun et al)
	#cconc = temp_cconc
        #FCO2 = TRUE

        latlong <- "+init=epsg:4326"

        #print(X)
        #print(Y)

	# Set up thermal time to growth stages and associated GAI (if not specified)
	if(is.null(GAItab)){
		
		GS <- c("Sowing-GS30","GS30-GS31","GS31-GS61","GS61-GS69","GS69-GS87","GS87-senescense")
		TT <- c(1100,100,900,50,750,200) # thresholds of cumulative degree days
		GAIstart <- c(0, 1.6, 2.0, 6.3, 6.3, 1.3)
		GAIend <- c(1.6, 2.0, 6.3, 6.3, 0.7, 0)	
	}

	# Calculate day length
	dates <- paste(substr(T,1,4), substr(T,5,6), substr(T,7,8), sep = "-")
        #print(dates)

        #print(lat)

	# Julian days
	Jday <- as.POSIXlt(strptime(T, format = "%Y%m%d"))$yday + 1
	# Days since sowing
        Cday <- array(1:dim(T)[1])
                                        #print(Cday)
    

	# Day length
	DL <- sapply(Jday, function(d) daylength(lat, d), simplify = "array")
        #print(DL)
 	# Day on which day length > threshold
	DL14 <- min(which(Jday > 1 & DL >= Dt))
        #print(DL14)
	# Degree days
        Tbase <- 0
        #print(tmin>Tbase)
        tmin <- array(tmin)
        tmax <- array(tmax)
	DD <-  ifelse(tmin > Tbase, ((tmax + tmin)/2)-Tbase,
				ifelse(tmax < Tbase, 0,
					ifelse(tmin <= Tbase, ((tmax-Tbase)/2) - ((Tbase - tmin)/4), 999)))
        #print(DD)
	# Vernalisation days
        tmean <- array(tmean)
        VD <- ifelse(tmean > 3 & tmean  <= 10, 1,
			ifelse(tmean  > -4 & tmean  <= 3, (tmean +4)/7,
				ifelse(tmean  > 10 & tmean  <= 17, (17-tmean )/7,
					ifelse(tmean  < -4, 0, 
						ifelse(tmean > 17, 0, 999)))))
        #print(dim(VD))
	# Set starting values for GAI and growth stage arrays
	GAI <- array(0, dim(tmean))
	GSS <- array("X", dim(tmean))
        #print(GAI)
        #print(GSS)
    
        ### GAI for Growth stages 0 - GS30
	CDD <- cumsum(DD)
	p1i <- which(CDD <= TT[1])
	GAI[p1i] <- ((CDD/TT[1])* GAIend[1])[p1i]
	GSS[p1i] <- "Sowing-GS30"
        #print(CDD)
        #print(GAI)
        #print(GSS)
    
	### GAI for Growth stages GS30 - GS31
	p2i <- which(CDD > TT[1] & CDD <= (TT[1] + TT[2]))
	DD2 <- DD
	DD2[-p2i] <- 0
	CDD2 <- cumsum(DD2)

	GAI[p2i] <- ((CDD2/TT[2]) * (GAIend[2] - GAIstart[2]) + GAIstart[2])[p2i]
	GSS[p2i] <- "GS30-GS31"

	### GAI for Growth stages GS31-GS61
	p3i <- which(CDD > (TT[1] + TT[2]) & CDD <= (TT[1] + TT[2] + TT[3]))
	
	# Check vernalisaiton
	GAI[p3i][cumsum(VD)[p3i] < 50] <- GAIstart[3]

	# Check day length
	GAI[p3i][Cday[p3i] < DL14[p3i]] <- GAIstart[3]

	# Reduce max GAI for GS61 if failed to accumulate sufficent degree days by DL14
	mGAI61 <- CDD
	mGAI61[Cday < DL14] <- NA
	#mGAI61 <- apply(mGAI61, 1, min, na.rm = T)
        mGAI61 <- min(mGAI61, na.rm=T)
        mGAI61 <- array(mGAI61)
        #print(dim(mGAI61))
        mGAI61 <- apply((1- ((TT[1] + TT[2]) - mGAI61)/(TT[1] + TT[2])) * GAIend[3], 1, function(x) min(x,GAIend[3]))

	# Set degree days in GS31-61 before adequate vernalisation or DL14 to zero
	DD[CDD > (TT[1] + TT[2]) & cumsum(VD) < 50] <- 0
	DD[CDD > (TT[1] + TT[2]) & Cday < DL14] <- 0

 	# Calculate GAI for following days
	CDD <- cumsum(DD)
	p4i <- which(CDD > (TT[1] + TT[2]) & CDD <= (TT[1] + TT[2] + TT[3])) 
	DD4 <- DD
	DD4[-p4i] <- 0
	CDD4 <- cumsum(DD4)
	GAI[setdiff(p3i, p4i)] <- GAIstart[3]
	GAI[p4i] <- (((CDD4/TT[3]) * (mGAI61 - GAIstart[3])) + GAIstart[3])[p4i]
	GSS[setdiff(p3i, p4i)] <- "GS31-GS61"
	GSS[p4i] <- "GS31-GS61"

	# Freeze GAI if temperature drops below -5 after GS31
	if(length(which(tmin <= -5 & GSS == "GS31-GS61")) >0){
		b5i <- which(tmin <= -5 & GSS == "GS31-GS61") 
		B5 <- GAI
		B5[-b5i] <- NA
		B5GAI <- array(min(B5, na.rm = T), dim(GAI))
		GAI[b5i] <- B5GAI[b5i]
	}
	
	### GAI for Growth stages GS61-GS69
	p5i <- which(CDD > (TT[1] + TT[2]+ TT[3]) & CDD <= (TT[1] + TT[2] + TT[3] + TT[4]))
	GAI[p5i] <- array(max(GAI), dim(GAI))[p5i]
	GSS[p5i] <- "GS61-GS69"

	### GAI for Growth stages GS69 - 87
	p6i <- which(CDD > (TT[1] + TT[2]+ TT[3] +TT[4]) & CDD <= (TT[1] + TT[2] + TT[3] + TT[4] + TT[5]))
	DD6 <- DD
	DD6[-p6i] <- 0
	CDD6 <- cumsum(DD6)
	GAI[p6i] <- ((CDD6/TT[5]) * (GAIend[5] - array(max(GAI, na.rm = T), dim(GAI))) + array(max(GAI, na.rm = T), dim(GAI)))[p6i] 
	GSS[p6i] <- "GS69-GS87"

	### GAI for Growth stages GS87 - senescense
	p7i <- which(CDD > (TT[1] + TT[2]+ TT[3] +TT[4] + TT[5]) & CDD <= (TT[1] + TT[2] + TT[3] + TT[4] +TT[5]  + TT[6]))
	DD7 <- DD
	DD7[-p7i] <- 0
	CDD7 <- cumsum(DD7)
	GAI6 <- GAI
	GAI6[-p6i] <- NA
	GAI6 <- array(min(GAI6, na.rm = T), dim(GAI))
	GAI[p7i] <- ((CDD7/TT[6]) * (GAIend[6] - GAI6 ) + GAI6)[p7i] 
	GSS[p7i] <- "GS87-senescense"


#plot(Cday[50,40,], GAI[50,40,], type= "o", col = as.factor(GSS[50,40,]))
#abline(v = 330)

	### Account for failure to senesce pre-harvest (ANNOTATE THIS SECTION)
	Jarray <- array(Jday)

	# Find day since sowing for harvest
	harvestCDD <- CDD          
        harvestCDD[Jarray > HarvestJday & Cday > 100] <- NA
        harvestCDD <- array(harvestCDD)
        CDD <- array(CDD)
	harvestCDD <- array(max(harvestCDD, na.rm = T),dim(CDD))

	# Find julian day of maximum GAI
	maxGAI <- array(max(GAI),dim(GAI))	
	MJarray <- Jarray
	MJarray[GAI != maxGAI] <- NA
	MJarray <- array(max(MJarray, na.rm = T),dim(Jarray))

	# Porptional progress from julian day of maximum GAI to harvest
	d <- Jarray - MJarray
	d[Cday < 100 | d < 0] <- 0
	d[Jarray > HarvestJday] <- 0
	d <- d/array(max(d, na.rm = T),dim(d))

	p8i <- which(harvestCDD < (TT[1] + TT[2] + TT[3] + TT[4] + TT[5] + TT[6]) & Jarray <= HarvestJday & Jarray > MJarray)
	GAI[p8i] <- ((d * (GAIend[6] - maxGAI)) + maxGAI)[p8i]

	GAI[Jarray > HarvestJday & Cday > 100] <- 0

	### Account for failure to flower pre-harvest!
	GAI[array(max(CDD),dim(CDD)) < (TT[1] + TT[2] + TT[3] + TT[4])] <- 0

        
        # split function here
                                        #points(Cday[50,40,], GAI[50,40,], type= "l", col = "green")
        FCO2 <- TRUE
        #print(X)
        #print(Y)
        #print(GAI)
        data <- list('GAI'=GAI, 'tmean'=tmean, 'prec'=prec, 'solarrad'=solarrad, 'Jarray'=Jarray, 'Cday'=Cday, 'GSS'=GSS, 'HarvestJday'=HarvestJday, 'x'=X, 'y'=Y, 't'=dates)
      	return(data)
      }
        # AFAIK GAI=LAI
	### Converting GAI to biomass
                                        # % PAR intercepted by canopy



wheat_yield <- function(GAI, tmean, prec, solarrad, AWC, Jarray, Cday, GSS, HarvestJday, X, Y, cconc=NULL, FCO2=TRUE, RUE=3.1, WLP=117.7, savepath=NULL) {

    
        LI <- 1-exp(-0.5*GAI)

        # CO2 fertilisation effect
   	if(FCO2 == TRUE){
		RUE <-  ifelse(cconc < 350, RUE, ifelse(cconc > 750, RUE + RUE * (750/350 - 1) * 0.333, RUE + RUE * (cconc/350 - 1) * 0.333)) 
	}
	# Convert intercepted PAR to water-unlimited biomass (g/m2).  Note conversion of Watts to MJ
        PAR <- (solarrad*0.0036*24) * 0.5
	WU_Biomass <- (PAR * LI) * RUE

#par(mfrow = c(1,2))
#image(X,Y, AWC[,,1], asp =1)
#image(AWCraster, asp =1)

	# Set up PAWC array
	PAWC <- array(0, dim(AWC))
    
	# Assume saturation prior to Julian day 77
	PAWC[which(Jarray <= 77 | Cday < 100)] <- AWC[which(Jarray <= 77 | Cday < 100)]
    
	# Calcuate daily PAWC
	for(i in 1:dim(PAWC)[3]){
		ji <- which(Jarray[,,i] >= 77)
		PAWC[,,i][ji] <- pmin(AWC, pmax((PAWC[,,max(1,i-1)] + prec[,,i]) - (WU_Biomass[,,max(1,i-1)]/5), 0))[ji]
	}
    
	# Calculate water-limited biomass
	WL_Biomass <- ifelse(PAWC>0, WU_Biomass, 0)
    
	### Partitioning water-unlimited biomass to grain (g/m2)
        WUbmy <- apply(array(ifelse(GSS %in% c("GS69-GS87","GS87-senescense"), WU_Biomass, 0),dim(GSS)), 1:2, sum)
	WUwsc <-  apply(array(ifelse(GSS %in% c("GS31-GS61"), WU_Biomass, 0),dim(GSS)), 1:2, sum) * 0.265
	# Calculate water-unlimited yield (adjusted for water content)
	WUyield <- ((WUbmy+WUwsc)/1000000 * 10000) * 1.15 
	### Partitioning water limited biomass to grain (g/m2)
        bmy <- apply(array(ifelse(GSS %in% c("GS69-GS87","GS87-senescense"), WL_Biomass, 0),dim(GSS)), 1:2, sum)
	wsc <- apply(array(ifelse(GSS %in% c("GS31-GS61"), WL_Biomass, 0),dim(GSS)), 1:2, sum) * 0.265
    
	# Calculate water limited yield (adjusted for water content)
	WLyield <- ((bmy+wsc)/1000000 * 10000) * 1.15
    
	### Calculate consecutive days of waterlogging after onset of flowering 
	# Start of flowering
	fstart <- Cday
        fstart[GSS != "GS61-GS69"] <- NA
        #print(dim(fstart))
        fstart <- array(apply(fstart, 1:2, min, na.rm = T), dim(GAI))
    
	# Consecutive days of waterlogging
	wlog <- array(0,dim(fstart))
	wlog[which(Cday >= fstart & Jarray <= HarvestJday & (PAWC+prec) > AWC)] <- 1
	wlog <- aperm(apply(wlog, 1:2, function(x){
		ifelse(x == 1, unlist(sapply(rle(x)$lengths, seq)), 0)
	}), c(2,3,1))
    
	# For each consecutive day of waterlogging beyond 5 days, lose 117.7 kg/ha
	wloss <- apply(ifelse(wlog > 5, WLP/1000, 0), 1:2, sum)
    
	WLPyield <- (WLyield/1.15 - wloss) * 1.15	
#image(X,Y, WLPyield, asp =1)

        # Convert final yield to raster
        ukgrid <- "+init=epsg:27700"
 	WLPY=list()
 	WLPY$x= X
 	WLPY$y= Y
 	WLPY$z= WLPyield
	WLPYrast <- raster(WLPY, crs = ukgrid)

	if(is.null(savepath) == FALSE){
		writeRaster(WLPYrast, savepath,"GTiff", overwrite = TRUE)
	} 	
	return(WLPY)
}


# as wheat_yield but for point driving data
wheat_yield_point <- function(GAI, tmean, prec, solarrad, AWC, Jarray, Cday, GSS, HarvestJday, cconc=NULL, FCO2=TRUE, RUE=3.1, WLP=117.7, savepath=NULL) {

    
        LI <- 1-exp(-0.5*GAI)
    
	# CO2 fertilisation effect
   	if(FCO2 == TRUE){
		RUE <-  ifelse(cconc < 350, RUE, ifelse(cconc > 750, RUE + RUE * (750/350 - 1) * 0.333, RUE + RUE * (cconc/350 - 1) * 0.333)) 
	}
    
	# Convert intercepted PAR to water-unlimited biomass (g/m2).  Note conversion of Watts to MJ
        PAR <- (solarrad*0.0036*24) * 0.5
	WU_Biomass <- (PAR * LI) * RUE
          

#par(mfrow = c(1,2))
#image(X,Y, AWC[,,1], asp =1)
#image(AWCraster, asp =1)

	# Set up PAWC array
	PAWC <- array(0, dim(AWC))
        # Assume saturation prior to Julian day 77
        PAWC[which(Jarray <= 77 | Cday < 100)] <- AWC[which(Jarray <= 77 | Cday < 100)]
        
	# Calcuate daily PAWC
        for(i in 1:dim(PAWC)){
            if(!is.na(Jarray[i])) {
		if(Jarray[i] >= 77) { #which(Jarray[i] >= 77)
                    PAWC[i] <- pmin(AWC, pmax((PAWC[max(1,i-1)] + prec[i]) - (WU_Biomass[max(1,i-1)]/5), 0))
                    }}
	}
	# Calculate water-limited biomass
	WL_Biomass <- ifelse(PAWC>0, WU_Biomass, 0)
	### Partitioning water-unlimited biomass to grain (g/m2)

        GSS <- array(GSS)
	WUbmy <- sum(array(ifelse(GSS %in% c("GS69-GS87","GS87-senescense"), WU_Biomass, 0),dim(GSS)))
        WUwsc <- sum(array(ifelse(GSS %in% c("GS31-GS61"), WU_Biomass, 0),dim(GSS))) * 0.265
        WUbmy <- array(WUbmy)
        WUwsc <- array(WUwsc)
	# Calculate water-unlimited yield (adjusted for water content)
        WUyield <- ((WUbmy+WUwsc)/1000000 * 10000) * 1.15
        WUyield <- array(WUyield)
	### Partitioning water limited biomass to grain (g/m2)
	bmy <- sum(array(ifelse(GSS %in% c("GS69-GS87","GS87-senescense"), WL_Biomass, 0),dim(GSS)))
        wsc <- sum(array(ifelse(GSS %in% c("GS31-GS61"), WL_Biomass, 0),dim(GSS)) * 0.265)
        bmy <- array(bmy)
        wsc <- array(wsc)
    
	# Calculate water limited yield (adjusted for water content)
        WLyield <- ((bmy+wsc)/1000000 * 10000) * 1.15
        WLyield <- array(WLyield)
        ### Calculate consecutive days of waterlogging after onset of flowering 
	# Start of flowering
	fstart <- Cday
	fstart[GSS != "GS61-GS69"] <- NA
	fstart <- array(apply(fstart, 1, min, na.rm = T), dim(GAI))	
	# Consecutive days of waterlogging
	wlog <- array(0,dim(fstart))
	wlog[which(Cday >= fstart & Jarray <= HarvestJday & (PAWC+prec) > AWC)] <- 1
	wlog <- apply(wlog, 1, function(x){
		ifelse(x == 1, unlist(sapply(rle(x)$lengths, seq)), 0)
	})
    
        # For each consecutive day of waterlogging beyond 5 days, lose 117.7 kg/ha
        wlog <- array(wlog)
        #print(wlog)
        wloss <- sum(array(ifelse(wlog > 5, WLP/1000, 0)))
        wloss <- array(wloss)
        WLPyield <- (WLyield/1.15 - wloss) * 1.15
        WLPyield <- array(WLPyield)
        
#image(X,Y, WLPyield, asp =1)

        # Convert final yield to raster
        #ukgrid <- "+init=epsg:27700"
 	WLPY=list()
 	#WLPY$x= X
        #WLPY$y= Y
        #print(WLPyield)
 	WLPY$z= WLPyield
	#WLPYrast <- raster(WLPY, crs = ukgrid)#

	#if(is.null(savepath) == FALSE){
	#	writeRaster(WLPYrast, savepath,"GTiff", overwrite = TRUE)
	#} 	
	return(WLPY)
}

# Run for each # MODIFIED TO LOOP THROUGH SCENARIO RUNS
run_all <- function(){
rlist <- list.files("/users/sgsys/matbro/cropNET/data/wheat_driving", full.names = T)
dateranges <- t(sapply(unique(sapply(strsplit(rlist, "_"), function(x) paste(x[4], gsub(".RData","",x[5]), sep = "_"))), function(z){
	c(as.numeric(strsplit(z,"_")[[1]][1]),as.numeric(strsplit(z,"_")[[1]][2]))
}))

yrs <- 1982:2079
sms <- c("01","04","05","06","07","08","09","10","11","12","13","15")

# list of years (yrs) to run for each scenario (sms)
allc <- expand.grid(yrs, sms)

findgyear <- function(xn, ggyear) {
	ifelse(as.numeric(substr(dimnames(xn)[[3]], 5,8)) >= as.numeric("0901"),   as.numeric(substr(dimnames(xn)[[3]],1,4)) +1, as.numeric(substr(dimnames(xn)[[3]],1,4))) == ggyear
}

for(XC in 1:nrow(allc)){
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
 
	wheat_py(temp_tas,temp_tasmax, temp_tasmin, temp_pr, temp_rss, cconc = temp_cconc, FCO2 = TRUE, savepath = paste("/users/sgsys/matbro/cropNET/data/output/wheat", allc[XC,1],allc[XC,2], sep = "_"))
}
}








#### Compare with suveyed yield data
compare <- function(){
# Individual hectads
yw <- readOGR("N:\\ASSIST\\WP1 Yield data\\DEFRA Yield\\Analysis outputs", "YieldStats_Update_wheat")
yw$Mean_Yield

# Load yield trends from survey data
load("N:\\ASSIST\\WP1 Yield data\\DEFRA Yield\\Analysis outputs\\CollatedYieldData_withPre2010.R")
sag <- aggregate(coll_ydat$Yield[coll_ydat$Crop == "wheat"], list(coll_ydat$Year[coll_ydat$Crop == "wheat"]), mean)

# Modelled potential yield for surveyd years
ystack <- stack("N:\\CropNet\\UKCP18 Yield\\wheat_2008.tif",
	"N:\\CropNet\\UKCP18 Yield\\wheat_2009.tif",
	"N:\\CropNet\\UKCP18 Yield\\wheat_2010.tif",
	"N:\\CropNet\\UKCP18 Yield\\wheat_2011.tif",
	"N:\\CropNet\\UKCP18 Yield\\wheat_2012.tif",
	"N:\\CropNet\\UKCP18 Yield\\wheat_2013.tif",
	"N:\\CropNet\\UKCP18 Yield\\wheat_2014.tif",
	"N:\\CropNet\\UKCP18 Yield\\wheat_2015.tif",
	"N:\\CropNet\\UKCP18 Yield\\wheat_2016.tif",
	"N:\\CropNet\\UKCP18 Yield\\wheat_2017.tif")

# Averages per year (whole dataset)
natav <- cellStats(ystack, mean)

# Averages across surveyed cells and per cell
survav <- extract(ystack, yw, fun = mean, na.rm= TRUE)
spy <- apply(survav,2,mean)
spc <- apply(survav,1,mean)

# Compare with plots
par(mfrow = c(2,2), las = 1, mar = c(4,4,2,1), cex = 1.1)

plot(yw$Mean_Yield, spc, xlab = "Mean yield over time per hectad \n(surveyed)", ylab = "Mean yield over time \n(UKCP18 modelled)")
plot(sag$x, spy, xlab = "Mean yield per year \n(surveyed)", ylab = "Mean yield per year \n(UKCP18 modelled)", xlim = c(6,10), ylim = c(9, 16))
text(sag$x, spy, 2008:2017, col = "red", pos= 3)
cor.test(sag$x, spy)

### Look at projected yield over time
allrasts <- list.files("N:\\CropNet\\UKCP18 Yield", "tif", full.names =T)
allrasts <- stack(allrasts)

meany <- cellStats(allrasts, mean)
plot(1982:2080, meany, type = "l", col = "blue", xlab = "Year", ylab = "Mean potential yield (GB)")
abline(v = 2019, col ="red", lty = 2)

meany1 <- extract(allrasts, yw[yw$Hectad == "TL49",], fun = mean, na.rm= TRUE)
plot(1982:2080, meany1, type = "l", col ="blue", xlab = "Year", ylab = "Mean potential yield (TL49)")
abline(v = 2019, col ="red", lty = 2)
}


testfunc <- function(a){
  b=a+1
  return(b)
}




















