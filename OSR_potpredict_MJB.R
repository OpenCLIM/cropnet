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
osr_py <-function(tmean, tmax, tmin, prec, solarrad, AWC, X, Y, times, lats,
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
  print('CO2 concentration: ')
  print(cconc)
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
  print('Any Vf NAs:')
  print(any(is.na(Vf)))
  print('All NAs?:')
  print(all(is.na(Vf)))
  print('Vf max:')
  print(max(Vf, na.rm=T))

  ## Photoperiod ------------- 

  ## Julian days
  ## dates <- paste(substr(T,1,4),substr(T,5,6), substr(T,7,8),sep = "-")
  Jday <- as.POSIXlt(strptime(times, format = "%Y%m%d"))$yday + 1
    
  ## Day length
  DL <- sapply(Jday, function(d) daylength(lats, d), simplify = "array")
  print('Any DL NAs:')
  print(any(is.na(DL)))
  print('All NAs?:')
  print(all(is.na(DL)))
    
  ## Photoperiod factor
  Pbase <- 5.7
  Pmax <- 14.8
  Pf <- ifelse(DL > Pbase, ifelse(DL < Pmax, (DL-Pbase) / (Pmax - Pbase),1), 0)
  Pf <- aperm(apply(Pf, 1:2, function(x){
        x[which(is.na(x))] <- (x[(min(which(is.na(x)))-1)] + x[(max(which(is.na(x)))+1)])/2
        return(x)
	}), c(2,3,1))
  print('Any Pf NAs:')
  print(any(is.na(Pf)))
  print('All NAs?:')
  print(all(is.na(Pf)))
  print('Pf max:')
  print(max(Pf, na.rm=T))
  
  ## Effective daily temperature
  Tbase <- 4.5
  Teff <- ifelse(tmean > Tbase, tmean - Tbase, 0)
  print('Any Teff NAs:')
  print(any(is.na(Teff)))
  print('All NAs?:')
  print(all(is.na(Teff)))
  print('Teff max:')
  print(max(Teff, na.rm=T))
    
  ## Growth response to daily temperature
  at <- c(7.6, 2.0, 5.1, 1.5)*(10^-3)
    
  ### Daily development rate
  
  ## Stage 0-1 (sowing to emrgence)
  devRate <- Teff * at[1]
  devStage <- aperm(apply(devRate, 1:2, cumsum), c(2,3,1)) 

  print('Any devStage1 NAs:')
  print(any(is.na(devStage)))
  print('All NAs?:')
  print(all(is.na(devStage)))
  print('max:')
  print(max(devStage, na.rm=T))
    
  ## Stage 1-2 (emergence to onset of flowering)
  devRate <- ifelse(devStage >= 1,  Teff * at[2] * Vf * Pf  ,devRate)
  devStage <- aperm(apply(devRate, 1:2, cumsum), c(2,3,1)) 

  print('Any devStage2 NAs:')
  print(any(is.na(devStage)))
  print('All NAs?:')
  print(all(is.na(devStage)))
  print('max:')
  print(max(devStage, na.rm=T))
    
  ## Stage 2-3 (flowering)
  devRate <- ifelse(devStage >= 2,  Teff * at[3]  ,devRate)
  devStage <- aperm(apply(devRate, 1:2, cumsum), c(2,3,1))

  print('Any devStage3 NAs:')
  print(any(is.na(devStage)))
  print('All NAs?:')
  print(all(is.na(devStage)))
  print('max:')
  print(max(devStage, na.rm=T))    
    
  ## Stage 3-4 (end of flowering to maturity)
  devRate <- ifelse(devStage >= 3,  Teff * at[3]  ,devRate)
  devStage <- aperm(apply(devRate, 1:2, cumsum), c(2,3,1)) 
  devStage <- ifelse(devStage >= 4,  4 ,devStage)

  print('Any devStage4 NAs:')
  print(any(is.na(devStage)))
  print('All NAs?:')
  print(all(is.na(devStage)))
  print('max:')
  print(max(devStage, na.rm=T))
    
  ## Inital LAI
  LAIs <- 0.53
  RGRlai <- 0.00732
  LAI <- LAIs *  exp(RGRlai * aperm(apply(Teff, 1:2, cumsum), c(2,3,1)))
  LAI <- ifelse(LAI > 0.75, 2, LAI)

  print('Any Initial LAI NAs:')
  print(any(is.na(LAI)))
  print('All NAs?:')
  print(all(is.na(LAI)))
  print('max:')
  print(max(LAI, na.rm=T))
  
  ## Residual light 
  Klai <- 0.75
  PAR <- solarrad*(0.0036*24)/0.5

  print('Any PAR NAs')
  print(any(is.na(PAR)))
  print('All NAs?:')
  print(all(is.na(PAR)))
  print('max:')
  print(max(PAR, na.rm=T))
    
  # Green canopy reflectance
  REFg <- ifelse(devStage < 2, 0.046, ifelse(devStage > 3, 0.072, 0))

  print('Any REFg NAs:')
  print(any(is.na(REFg)))
  print('All NAs?:')
  print(all(is.na(REFg)))
  print('max:')
  print(max(REFg, na.rm=T))
    
  # Flowering canopy transmission
  TR   <- ifelse(devStage > 2 & devStage <3,  1 - ((devStage - 2) * 0.5)   ,1)
  print('Any TR NAs:')
  print(any(is.na(TR)))
  print('All NAs?:')
  print(all(is.na(TR)))
  print('max:')
  print(max(TR, na.rm=T))
    
  lt <- PAR * TR *(1-REFg) * exp(-Klai * LAI)

  print('Any lt1 NAs:')
  print(any(is.na(lt)))
  print('All NAs?:')
  print(all(is.na(lt)))
  print('max:')
  print(max(lt, na.rm=T))
    
  lt <-ifelse(devStage > 3.5, PAR * (exp(-Klai * LAI) + 0.184*(devStage - 3.5)), lt)

  print('Any lt2 NAs:')
  print(any(is.na(lt)))
  print('All NAs?:')
  print(all(is.na(lt)))
  print('max:')
  print(max(lt, na.rm=T))

  ## Light absorbtion	
  la <- (PAR * TR *(1-REFg) ) - lt # MJ per m2

  print('Any la NAs:')
  print(any(is.na(la)))
  print('All NAs?:')
  print(all(is.na(la)))
  print('max:')
  print(max(la, na.rm=T))
  
  ## Crop growth rate
  RUE <- (ifelse(devStage < 2, 0.6, ifelse(devStage < 3, 1.4, 0.7)))

  print('Any RUE NAs:')
  print(any(is.na(RUE)))
  print('All NAs?:')
  print(all(is.na(RUE)))
  print('max:')
  print(max(RUE, na.rm=T))
    
  ## CO2 fertilisation effect
  # make cconc 3D so ifelse RUE statement works as expected
  cconc <- array(rep(cconc, each = dim(tmean)[1]*dim(tmean)[2]), dim(tmean))
  print('cconc dims:')
  print(dim(cconc))
  if(FCO2 == TRUE){
     RUE <-  ifelse(cconc < 350, RUE,
             ifelse(cconc > 750, RUE + RUE * (750/350 - 1) * 0.333, RUE + RUE * (cconc/350 - 1) * 0.333))
     print('Any CO2 RUE NAs:')
     print(any(is.na(RUE)))
     print('All NAs?:')
     print(all(is.na(RUE)))
     print('max:')
     print(max(RUE, na.rm=T))
  }
  Wd <- la*RUE

  print('Any Wd NAs:')
  print(any(is.na(Wd)))
  print('All NAs?:')
  print(all(is.na(Wd)))
  print('max:')
  print(max(Wd, na.rm=T))
  
  ##  Adjust LAI
  SLA <- 0.0177/10      # cm/g -> m/g
  LAI2 <- ifelse(LAI > 0.75,  aperm(apply(Wd * SLA, 1:2, cumsum), c(2,3,1)), LAI)

  print('Any LAI2 stage 1 NAs:')
  print(any(is.na(LAI2)))
  print('All NAs?:')
  print(all(is.na(LAI2)))
  print('max:')
  print(max(LAI2, na.rm=T))
    
  LAI2[LAI > 0.75 & LAI2 < 0.75] <- 0.75

  print('Any LAI2 stage 2 NAs:')
  print(any(is.na(LAI2)))
  print('All NAs?:')
  print(all(is.na(LAI2)))
  print('max:')
  print(max(LAI2, na.rm=T))

  LAI2[LAI2 > 4]  <- 4
  maxLAI <- array(apply(LAI2, 1:2, max, na.rm = T), dim(LAI))

  print('Any maxLAI NAs:')
  print(any(is.na(maxLAI)))
  print('All NAs?:')
  print(all(is.na(maxLAI)))
  print('max:')
  print(max(maxLAI, na.rm=T))
    
  LAI2[!is.na(devStage) & devStage > 3] <-  (maxLAI -((devStage - 3)/1)*(maxLAI-2.5))[!is.na(devStage) & devStage > 3]

  print('Any LAI2 stage 3 NAs:')
  print(any(is.na(LAI2)))
  print('All NAs?:')
  print(all(is.na(LAI2)))
  print('max:')
  print(max(LAI2, na.rm=T))
  
  ## Adjust light absorbtion
  lt2 <- PAR * TR *(1-REFg) * exp(-Klai * LAI)

  print('Any lt2 stage 1 NAs:')
  print(any(is.na(lt2)))
  print('All NAs?:')
  print(all(is.na(lt2)))
  print('max:')
  print(max(lt2, na.rm=T))
   
  lt2 <- ifelse(devStage > 3.5, PAR * (exp(-Klai * LAI2) + 0.184*(devStage - 3.5)), lt2)

  print('Any lt2 stage 2 NAs:')
  print(any(is.na(lt2)))
  print('All NAs?:')
  print(all(is.na(lt2)))
  print('max:')
  print(max(lt2, na.rm=T))
    
  la2 <- (PAR * TR *(1-REFg) ) - lt2 # MJ per m2

  print('Any la2 NAs:')
  print(any(is.na(la2)))
  print('All NAs?:')
  print(all(is.na(la2)))
  print('max:')
  print(max(la2, na.rm=T))
    
  ## Crop dry matter biomass (g/m2) (water unlimited)
  WU_Biomass <- la2 * RUE
  print('Any WU_Biomass NAs:')
  print(any(is.na(WU_Biomass)))
  print('All NAs?:')
  print(all(is.na(WU_Biomass)))
  print('max:')
  print(max(WU_Biomass, na.rm=T))    

  Jarray <- array(rep(Jday,each = dim(tmean)[1]*dim(tmean)[2]), dim(tmean))
  Cday <- array(rep(1:dim(tmean)[3],each = dim(tmean)[1]*dim(tmean)[2]), dim(tmean))
    
  ### Partitioning water unlimited biomass to seed via harvest index 
  WU_Biomass[Jarray > HarvestJday & Cday > 150] <- 0
  print('Any WU_Biomass NAs:')
  print(any(is.na(WU_Biomass)))
  print('All NAs?:')
  print(all(is.na(WU_Biomass)))
  print('max:')
  print(max(WU_Biomass, na.rm=T))    

  WUyield <- apply(WU_Biomass,1:2, sum) * 0.35  #(g/m2)
  WUyield <- (WUyield/1000000)*10000  #(t/Ha2)
  print('Any WUyield NAs:')
  print(any(is.na(WUyield)))
  print('All NAs?:')
  print(all(is.na(WUyield)))
  print('max:')
  print(max(WUyield, na.rm=T))    
    
  ###  Calculate availible water 
  
  ## Method 2 (for loop)
  # Assume saturation prior to Julian day 77
  print('AWC max:')
  print(max(AWC, na.rm=T))
  PAWC <- array(0, dim(AWC))
  PAWC[which((Jarray <= 77 | Cday < 150))] <- AWC[which((Jarray <= 77 | Cday < 150))]
  PAWC[which(is.na(Jarray) == TRUE)] <- AWC[which(is.na(Jarray) == TRUE)]
    
  for(i in 1:dim(PAWC)[3]){
      ji <- which(Jarray[,,i] > 77 & Cday[,,i] > 150)
      PAWC[,,i][ji] <- pmin(AWC, pmax((PAWC[,,max(1,i-1)] + prec[,,i]) - (WU_Biomass[,,max(1,i-1)]/ifelse(devStage < 3, 5, 3.8)[,,max(1,i-1)]), 0))[ji]
	}
  print('PAWC max:')
  print(max(PAWC, na.rm=T))
    
  ## Calculate water-limited biomass
  WL_Biomass <- ifelse(PAWC>0, WU_Biomass, 0)
  print('WL_Biomass max:')
  print(max(WL_Biomass, na.rm=T))
    
  ## Consecutive days of waterlogging
  ## For each waterlogging event beyond 14 days, lose 50% biomass if pre-stem elongation... 
  wlog1x <- array(0,dim(PAWC))
  wlog1x[which(Jarray <= HarvestJday & (PAWC+prec) > (AWC + 10) & devStage < 2)] <- 1
  wlog1x <- aperm(apply(wlog1x, 1:2, function(x){
                  ifelse(x == 1, unlist(sapply(rle(x)$lengths, seq)), 0)
                  }), c(2,3,1))
  print('Any wlog1x NAs:')
  print(any(is.na(wlog1x)))
  print('All NAs?:')
  print(all(is.na(wlog1x)))
  print('max:')
  print(max(wlog1x, na.rm=T))    
  
  pwloss1 <- ifelse(wlog1x > 5 & wlog1x <= 14, (wlog1x-5)/(14-5)*0.5, ifelse(wlog1x > 14,0.5,0))
  WLWL_Biomass <- WL_Biomass*(1-pwloss1)
    
  ## ...15% if flowering (Wollmer et al 2017)
  wlog2x <- array(0,dim(PAWC))
  wlog2x[which(Jarray <= HarvestJday & (PAWC+prec) > (AWC + 10) & devStage >= 2)] <- 1
  wlog2x <- aperm(apply(wlog2x, 1:2, function(x){
		  ifelse(x == 1, unlist(sapply(rle(x)$lengths, seq)), 0)
                  }), c(2,3,1))
  print('Any wlog1x NAs:')
  print(any(is.na(wlog2x)))
  print('All NAs?:')
  print(all(is.na(wlog2x)))
  print('max:')
  print(max(wlog2x, na.rm=T))    
  
  pwloss2 <- ifelse(wlog1x > 5 & wlog1x <= 14, (wlog1x-5)/(14-5)*0.15, ifelse(wlog1x > 14,0.5,0))
    
  WLWL_Biomass <- WLWL_Biomass*(1-pwloss2)
  print('WLWL_Biomass max:')
  print(max(WLWL_Biomass, na.rm=T))


    
  ### Partitioning water limited biomass to seed via harvest index 
  WLWL_Biomass[Jarray > HarvestJday & Cday > 150] <- 0
  WLyield <- apply(WLWL_Biomass,1:2, sum) * 0.35  #(g/m2)
  WLYield <- (WLyield/1000000)*10000  #(t/Ha2)
  print('WLYield max:')
  print(max(WLYield, na.rm=T))
    
  ### Alternative apporach assuming only biomass post-flowering availible for yield
  AWLWL_Biomass <- ifelse(devStage > 3,  WLWL_Biomass, 0.1 * WLWL_Biomass)
  WLYield2 <- (apply(AWLWL_Biomass, 1:2, sum)/1000000)*10000

  ### Calculate HDDs during Stage 2-3
    
  ## Set all non-stage2-3 temps to 0
  tmin_s2  <- ifelse(devStage >= 2 & devStage <=3, tmin, 0)
  tmean_s2 <- ifelse(devStage >= 2 & devStage <=3, tmean, 0)
  tmax_s2  <- ifelse(devStage >= 2 & devStage <=3, tmax, 0)
  print('tmin_s2 max:')
  print(max(tmin_s2, na.rm=T))
  print('tmean_s2 max:')
  print(max(tmean_s2, na.rm=T))
  print('tmax_s2 max:')
  print(max(tmax_s2, na.rm=T))    
    
  ## Calculate HDDs with baseline temp of 29.5C
  HDD <- ifelse(tmin_s2 > 29.5, ((tmin_s2 + tmax_s2)/2) -29.5,
         ifelse(tmax_s2 < 29.5, 0,
         ifelse(tmax_s2 >= 29.5, (((29.5 + tmax_s2)/2) -29.5) * ((tmax_s2-29.5)/(tmax_s2-tmin_s2)), 999)))
  print('HDD max:')
  print(max(HDD, na.rm=T))    

  # Calculate total HDDs across stage2
  CHDD <- apply(HDD, 1:2, sum)
  print('CHDD max:')
  print(max(CHDD, na.rm=T))    
    
  # Apply penalty of 0.2T/Ha per HDD to the yields
  WLHLyield <- WLYield - 0.2*CHDD
  WUHLyield <- WUyield - 0.2*CHDD
  print('Max yield penalty:')
  print(max(0.2*CHDD, na.rm=T))
  print('Max HL yield difference:')
  print(max(WLYield-WLHLyield, na.rm=T))
  print('WLHLyield max:')
  print(max(WLHLyield, na.rm=T))
  print('WUHLyield max:')
  print(max(WUHLyield, na.rm=T))    
  
  data <- list('WUyield'=WUyield, 'WLyield'=WLYield, 'WLHLyield'=WLHLyield, 'WUHLyield'=WUHLyield, 'WU_Biomass'=WU_Biomass, 'WL_Biomass'=WL_Biomass, 'PAWC'=PAWC, 'AWC'=AWC, 'tmean'=tmean, 'prec'=prec, 'HDD'=HDD, 'CHDD'=CHDD, 'devStage'=devStage)
  return(data)
}
