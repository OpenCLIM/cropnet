library(raster)
library(rgdal)
library(rgeos)
library(geosphere)
library(data.table)
library(ncdf4)
library(abind)

### UKCP18 NetCDF processing

# List NetCDF files
nlist <- list.files("/data/UKCP18/RCM_12km/daily_timeseries/downloaded_data", ".nc", recursive = TRUE, full.names = T)

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
	X <- ncvar_get(ncin, "projection_x_coordinate")
	Y <- ncvar_get(ncin, "projection_y_coordinate")
	
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
	save(list = paste(vars[V], "UKCP18", sep = "_"), file = paste("/users/sgsys/matbro/cropNET/data/wheat_driving/",vars[V],"_UKCP18_",min(allgroyears),"_",max(allgroyears),"_",simx,".RData",sep = ""))
}
}	
}
