library(raster)
library(rgdal)
library(rgeos)
library(geosphere)
library(data.table)
library(ncdf4)
library(abind)

# Read soils data, reproject to grid and aggregate
load_AWC <- function(X, Y, loc, datares, dims){
coords <- expand.grid(as.numeric(X), as.numeric(Y))
dat_SP <- SpatialPoints(coords, proj4string = CRS("+init=epsg:27700"))
## Read soils data, reproject and aggregate to coarser resolution
ukgrid <- "+init=epsg:27700"
latlong <- "+init=epsg:4326"
AWCrast <- raster(loc)
projection(AWCrast) <- projection(ukgrid)
AWCres <- res(AWCrast)[0]
resrat <- datares/AWCres
if (resrat >= 2) {
   AWCraster <- aggregate(AWCrast, floor(resrat), fun=mean)
}
##print(AWCraster)
##print(dim(AWCraster))
## Extract soil AWC per grid cell - constant over time. 
## Below function regrids the AWC data on to the model grid
AWC <- array(extract(AWCraster,  dat_SP, fun = mean, na.rm = TRUE), dims)
return(AWC)
}

load_AWC_no_agg <- function(X, Y, loc, dims){
coords <- expand.grid(as.numeric(X), as.numeric(Y))
dat_SP <- SpatialPoints(coords, proj4string = CRS("+init=epsg:27700"))
## Read soils data, reproject and aggregate to coarser resolution
ukgrid <- "+init=epsg:27700"
latlong <- "+init=epsg:4326"
AWCrast <- raster(loc)
projection(AWCrast) <- projection(ukgrid)
#AWCraster <- aggregate(AWCrast, 10, fun=mean)
##print(AWCraster)
##print(dim(AWCraster))
## Extract soil AWC per grid cell - constant over time. 
## Below function regrids the AWC data on to the model grid
AWC <- array(extract(AWCrast,  dat_SP, fun = mean, na.rm = TRUE), dims)
return(AWC)
}


# Set growth model function
osr_py <-function(tmean, tmax, tmin, prec, solarrad, AWC, X, Y, T, lats,
                   cconc = NULL, Tbase = 4.5, GAItab = NULL, HarvestJday = 212,
                   k = 0.75, datasetname = 'ukcp18', precipname = 'None',
                   radname = 'None', FCO2 = TRUE){
  
  if (datasetname %in% "ukcp18") {
    print('Using UKCP18 data, units:')
    print('Temperature: Celsius')
    if (precipname %in% "aphrodite") {
      print('Using aphrodite for precip, units mm')
    } else {
      print('Precip: mm')
    }
    if (grepl('ceres', radname, fixed=TRUE)) {
      print('Using ceres for solar rad, units W/m^2')
    } else {
      print('Solar rad: W/m^2')
    }
  }
  if (grepl("chess-scape", datasetname, fixed=TRUE)) {
    print('Using UKCP18 1km (CHESS-SCAPE) data, units:')
    print('Temperature: Converting from Kelvin to Celsius')
    tmean <- tmean - 273.15                                             
    tmax <- tmax - 273.15
    tmin <- tmin - 273.15 ## all Kelvin --> Celsius
    if (precipname %in% "aphrodite") {
      print('Using aphrodite for precip, units mm')
    } else {
      print('Precip: Converting from kg/m^2/s to mm/day')
      prec <- prec*86400
    }
    if (grepl('ceres', radname, fixed=TRUE)) {
      print('Using ceres for solar rad, units W/m^2')
    } else {
      print('Solar rad: W/m^2')
    }
  }


  if (datasetname %in% "chess_and_haduk") {
    print('Using chess-met and HadUK data, units:')
    print('Temperature: Celsius')
    if (precipname %in% "aphrodite") {
      print('Using aphrodite for precip, units mm')
    } else {
      print('Precip: Converting from kg/m^2/s to mm/day')
      prec <- prec*86400
    }
    if (grepl('ceres', radname, fixed=TRUE)) {
      print('Using ceres for solar rad, units W/m^2')
    } else {
      print('Solar rad: W/m^2')
    }
  }  

  
  ## These conversions only necessary for era5 data
  if (datasetname %in% "era5") {
    print('Using era5 data, units:')
    print('Temperature: Converting from Kelvin to Celsius')
    tmean <- tmean - 273.15
    tmax <- tmax - 273.15
    tmin <- tmin - 273.15 ## all Kelvin --> Celsius
    if (!(precipname %in% "aphrodite")) { ## only if we're not using aphrodite for precip
      print('Precip: Converting from m to mm')
      prec <- prec*1000. ## m --> mm
    } else {
      print('Using aphrodite for precip, units mm')
    }
    if (!(grepl('ceres', radname, fixed=TRUE))) { ## only if we're not using ceres for radiation
      print('Solar rad: converting from J/m^2/day to W/m^2')
      solarrad <- solarrad*(1/86400) ## J/m^2/day --> W/m^2
    } else {
      print('Using ceres for solar rad, units W/m^2')
    }
  }
  
  ukgrid <- "+init=epsg:27700"
  latlong <- "+init=epsg:4326"
    
  ## Degree of vernalisiton
  ## Vernalisation rate (unit increase per day)
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

  ## Photoperiod ------------- 

  ## Julian days
  dates <- paste(substr(T,1,4),substr(T,5,6), substr(T,7,8),sep = "-")
  Jday <- as.POSIXlt(strptime(T, format = "%Y%m%d"))$yday + 1
  
  ## Day length
  DL <- sapply(Jday, function(d) daylength(lats, d), simplify = "array")
  print(dim(DL))
    
  ## Photoperiod factor
  Pbase <- 5.7
  Pmax <- 14.8
  Pf <- ifelse(DL > Pbase, ifelse(DL < Pmax, (DL-Pbase) / (Pmax - Pbase),1), 0)
  print(dim(Pf))
  Pf <- aperm(apply(Pf, 1:2, function(x){
        x[which(is.na(x))] <- (x[(min(which(is.na(x)))-1)] + x[(max(which(is.na(x)))+1)])/2
        return(x)
	}), c(2,3,1))

  ## Effective daily temperature
  Tbase <- 4.5
  Teff <- ifelse(tmean > Tbase, tmean - Tbase, 0)
  
  ## Growth response to daily temperature
  at <- c(7.6, 2.0, 5.1, 1.5)*(10^-3)
    
  ### Daily development rate
  
  ## Stage 0-1 (sowing to emrgence)
  devRate <- Teff * at[1]
  devStage <- aperm(apply(devRate, 1:2, cumsum), c(2,3,1)) 
  
  ## Stage 1-2 (emergence to onset of flowering)
  devRate <- ifelse(devStage >= 1,  Teff * at[2] * Vf * Pf  ,devRate)
  devStage <- aperm(apply(devRate, 1:2, cumsum), c(2,3,1)) 
  
  ## Stage 2-3 (flowering)
  devRate <- ifelse(devStage >= 2,  Teff * at[3]  ,devRate)
  devStage <- aperm(apply(devRate, 1:2, cumsum), c(2,3,1)) 
  
  ## Stage 3-4 (end of flowering to maturity)
  devRate <- ifelse(devStage >= 3,  Teff * at[3]  ,devRate)
  devStage <- aperm(apply(devRate, 1:2, cumsum), c(2,3,1)) 
  devStage <- ifelse(devStage >= 4,  4 ,devStage)

  ## Inital LAI
  LAIs <- 0.53
  RGRlai <- 0.00732
  LAI <- LAIs *  exp(RGRlai * aperm(apply(Teff, 1:2, cumsum), c(2,3,1)))
  LAI <- ifelse(LAI > 0.75, 2, LAI)
  
  ## Residual light 
  Klai <- 0.75
  PAR <- temp_rss*(0.0036*24)/0.5
  # Green canopy reflectance
  REFg <- ifelse(devStage < 2, 0.046, ifelse(devStage > 3, 0.072, 0))
  # Flowering canopy transmission
  TR   <- ifelse(devStage > 2 & devStage <3,  1 - ((devStage - 2) * 0.5)   ,1)

  lt <- PAR * TR *(1-REFg) * exp(-Klai * LAI)
  lt <-ifelse(devStage > 3.5, PAR * (exp(-Klai * LAI) + 0.184*(devStage - 3.5)), lt)

  ## Light absorbtion	
  la <- (PAR * TR *(1-REFg) ) - lt # MJ per m2
  
  ## Crop growth rate
  RUE <- (ifelse(devStage < 2, 0.6, ifelse(devStage < 3, 1.4, 0.7)))
  
  ## CO2 fertilisation effect
  if(FCO2 == TRUE){
     RUE <-  ifelse(cconc < 350, RUE,
             ifelse(cconc > 750, RUE + RUE * (750/350 - 1) * 0.333, RUE + RUE * (cconc/350 - 1) * 0.333)) 
  }
  Wd <- la*RUE 
  
  ##  Adjust LAI
  SLA <- 0.0177/10      # cm/g -> m/g
  LAI2 <- ifelse(LAI > 0.75,  aperm(apply(Wd * SLA, 1:2, cumsum), c(2,3,1)), LAI)
  LAI2[LAI > 0.75 & LAI2 < 0.75] <- 0.75
  LAI2[LAI2 > 4]  <- 4 
  maxLAI <- array(apply(LAI2, 1:2, max, na.rm = T), dim(LAI))
  LAI2[devStage > 3] <-  (maxLAI -((devStage - 3)/1)*(maxLAI-2.5))[devStage > 3]
  
  ## Adjust light absorbtion
  lt2 <- PAR * TR *(1-REFg) * exp(-Klai * LAI)
  lt2 <- ifelse(devStage > 3.5, PAR * (exp(-Klai * LAI2) + 0.184*(devStage - 3.5)), lt2)
  la2 <- (PAR * TR *(1-REFg) ) - lt2 # MJ per m2

  ## Crop dry matter biomass (g/m2) (water unlimited)
  WU_Biomass <- la2 * RUE

  ### Partitioning water unlimited biomass to seed via harvest index 
  WU_Biomass[Jarray > HarvestJday & Cday > 150] <- 0
  WUyield <- apply(WU_Biomass,1:2, sum) * 0.35  #(g/m2)
  WUyield <- (WUyield/1000000)*10000  #(t/Ha2)
  
  ###  Calculate availible water 

  ## Assume saturation prior to Julian day 77
  Jarray <- array(rep(Jday,each = dim(tmean)[1]*dim(tmean)[2]), dim(tmean))
  Cday <- array(rep(1:dim(tmean)[3],each = dim(tmean)[1]*dim(tmean)[2]), dim(tmean))
  
  ## Method 2 (for loop)
  PAWC <- array(0, dim(AWC))
  PAWC[which((Jarray <= 77 | Cday < 150))] <- AWC[which((Jarray <= 77 | Cday < 150))]
  PAWC[which(is.na(Jarray) == T)] <- AWC[which(is.na(Jarray) == T)]
  
  for(i in 1:dim(PAWC)[3]){
      ji <- which(Jarray[,,i] > 77 & Cday[,,i] > 150)
      PAWC[,,i][ji] <- pmin(AWC, pmax((PAWC[,,max(1,i-1)] + prec[,,i]) - (WU_Biomass[,,max(1,i-1)]/ifelse(devStage < 3, 5, 3.8)[,,max(1,i-1)]), 0))[ji]
	}
	
  ## Calculate water-limited biomass
  WL_Biomass <- ifelse(PAWC>0, WU_Biomass, 0)
  
  ## Consecutive days of waterlogging
  ## For each waterlogging event beyond 14 days, lose 50% biomass if pre-stem elongation... 
  wlog1x <- array(0,dim(PAWC))
  wlog1x[which(Jarray <= HarvestJday & (PAWC+prec) > (AWC + 10) & devStage < 2)] <- 1
  wlog1x <- aperm(apply(wlog1x, 1:2, function(x){
                  ifelse(x == 1, unlist(sapply(rle(x)$lengths, seq)), 0)
                  }), c(2,3,1))
  
  pwloss1 <- ifelse(wlog1x > 5 & wlog1x <= 14, (wlog1x-5)/(14-5)*0.5, ifelse(wlog1x > 14,0.5,0))
  WLWL_Biomass <- WL_Biomass*(1-pwloss1)
  
  ## ...15% if flowering (Wollmer et al 2017)
  wlog2x <- array(0,dim(PAWC))
  wlog2x[which(Jarray <= HarvestJday & (PAWC+prec) > (AWC + 10) & devStage >= 2)] <- 1
  wlog2x <- aperm(apply(wlog2x, 1:2, function(x){
		  ifelse(x == 1, unlist(sapply(rle(x)$lengths, seq)), 0)
                  }), c(2,3,1))
  
  pwloss2 <- ifelse(wlog1x > 5 & wlog1x <= 14, (wlog1x-5)/(14-5)*0.15, ifelse(wlog1x > 14,0.5,0))
  WLWL_Biomass <- WLWL_Biomass*(1-pwloss2)
  
  ### Partitioning water limited biomass to seed via harvest index 
  WLWL_Biomass[Jarray > HarvestJday & Cday > 150] <- 0
  WLyield <- apply(WLWL_Biomass,1:2, sum) * 0.35  #(g/m2)
  WLYield <- (WLyield/1000000)*10000  #(t/Ha2)

  ### Alternative apporach assuming only biomass post-flowering availible for yield
  AWLWL_Biomass <- ifelse(devStage > 3,  WLWL_Biomass, 0.1 * WLWL_Biomass)
  WLYield2 <- (apply(AWLWL_Biomass, 1:2, sum)/1000000)*10000

  ### Calculate HDDs during Stage 2-3
    
  ## Set all non-stage2-3 temps to 0
  tmin_s2  <- ifelse(devStage > 2 & devStage <3, tmin, tmin-tmin)
  tmean_s2 <- ifelse(devStage > 2 & devStage <3, tmean, tmean-tmean)
  tmax_s2  <- ifelse(devStage > 2 & devStage <3, tmax, tmax-tmax)
  ## Calculate HDDs with baseline temp of 29.5C
  HDD <- ifelse(tmin_s2 > 29.5, ((tmin_s2 + tmax_s2)/2) -29.5,
         ifelse(tmax_s2 < 29.5, 0,
         ifelse(tmax_s2 >= 29.5, (((29.5 + tmax_s2)/2) -29.5) * ((tmax_s2-29.5)/(tmax_s2-tmin_s2)), 999)))

  print('Heat degree days dimensions:')
  print(dim(HDD))
  # Calculate total HDDs across stage2
  CHDD <- apply(HDD, 1:2, sum)
  print('Summed heat degree days dimensions:')
  print(dim(CHDD))
  print('Yield dimensions:')
  print(dim(WLyield))
  # Apply penalty of 0.2T/Ha per HDD to the yields
  WLHLyield <- WLyield - 0.2*CHDD
  WUHLyield <- WUyield - 0.2*CHDD
  
  data <- list('WUyield'=WUyield, 'WLyield'=WLyield, 'WLHLyield'=WLHLyield, 'WUHLyield'=WUHLyield, 'WU_Biomass'=WU_Biomass, 'WL_Biomass'=WL_Biomass, 'PAWC'=PAWC, 'AWC'=AWC, 'tmean'=tmean, 'prec'=prec)
  return(data)
}
