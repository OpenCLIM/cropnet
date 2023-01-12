library(ncdf4)
#library(Evapotranspiration)
library(rgdal)
library(rgeos)
library(raster)
library(abind)

load_elev <- function (X, Y, loc)
{
  ## Get elevation data
  ##print(X)
  ##print(Y)
  coords <- expand.grid(as.numeric(X), as.numeric(Y))
  ##print(coords)
  dat_SP <- SpatialPoints(coords, proj4string = CRS("+init=epsg:27700"))
  elev <- raster(loc)
  elevs <- matrix(extract(elev,  dat_SP, fun = mean, na.rm = TRUE), c(dim(X), dim(Y)))
  return(elevs)
}

grass_py <- function(Tt, Ttmx, Ttmn, prec, Rr,  rh, wind, X, Y, T, datasetname='ukcp18', precipname = 'None', radname = 'None', cconc = NULL, sfcP = NULL, elevfile = NULL, FCO2 = TRUE, Q = 18.81, SMDmax = 110, SMDc = 10, savepath = NULL) {

  if (datasetname %in% "era5") {
    print('using')
    print(datasetname)
    Tdp <- rh
  }
  
  if (datasetname %in% "ukcp18") {
    print('Using UKCP18 data, units:')
    print('Temperature: Celsius')
    if (precipname %in% "aphrodite") {
      print('Using aphrodite for precip, units mm')
    } else {
      print('Precip: mm')
    }
    if (grepl('ceres', radname, fixed=TRUE)) {
      print('Using ceres for solar rad, units W/m^2, converting to MJ/m^2/day')
      Rr <- Rr*0.0036*24 ## W/m^2 --> MJ/m^2/day		
    } else {
      Rr <- Rr*0.0036*24 ## W/m^2 --> MJ/m^2/day	
      print('Solar rad: W/m^2, converting to MJ/m^2/day')
    }
  }

  if (datasetname %in% "ukcp18bc") {
    print('Using UKCP18 Bias Corrected (CHESS-SCAPE) data, units:')
    print('Temperature: Converting from Kelvin to Celsius')
    Tt <- Tt - 273.15                                             
    Ttmx <- Ttmx - 273.15
    Ttmn <- Ttmn - 273.15 ## all Kelvin --> Celsius
    if (precipname %in% "aphrodite") {
      print('Using aphrodite for precip, units mm')
    } else {
      print('Precip: Converting from kg/m^2/s to mm/day')
      prec <- prec*86400
    }
    if (grepl('ceres', radname, fixed=TRUE)) {
      print('Using ceres for solar rad')
      print('Converting rad units from W/m^2 --> MJ/m^2/day')
      Rr <- Rr*0.0036*24 ## W/m^2 --> MJ/m^2/day
    } else {
      print('Converting rad units from W/m^2 --> MJ/m^2/day')
      Rr <- Rr*0.0036*24 ## W/m^2 --> MJ/m^2/day
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
      print('Using ceres for solar rad')
      print('Converting rad units from W/m^2 --> MJ/m^2/day')
      Rr <- Rr*0.0036*24 ## W/m^2 --> MJ/m^2/day
    } else {
      print('Converting rad units from W/m^2 --> MJ/m^2/day')
      Rr <- Rr*0.0036*24 ## W/m^2 --> MJ/m^2/day
    }
  }  

  ## These conversions only necessary for era5 data
  if (datasetname %in% "era5") {
    Tt <- Tt - 273.15
    Ttmx <- Ttmx - 273.15
    Ttmn <- Ttmn - 273.15
    Tdp <- Tdp - 273.15 ## all Kelvin --> Celsius
    if (!(precipname %in% 'aphrodite')) { ## only if we're not using aphrodite for precip
      prec <- prec*1000. ## m --> mm
    }
    if (!(grepl('ceres', radname, fixed=TRUE))) { ## only if we're not using ceres for solar radiation
      print('Converting rad units from J/m^2/day --> MJ/m^2/day')
      Rr <- Rr/1000000
    }
  }


################################################
  
  
  ## Get elevation data
  ##print(dim(X))
  ##print(dim(Y))
  coords <- expand.grid(as.numeric(X), as.numeric(Y))
  ##print(coords)
  if (!is.null(elevfile)){
    dat_SP <- SpatialPoints(coords, proj4string = CRS("+init=epsg:27700"))
    elev <- raster(elevfile)
    elevs <- matrix(extract(elev,  dat_SP, fun = mean, na.rm = TRUE), c(dim(X), dim(Y)))
    elevs <- t(elevs)
    print(elevs)
  }
  ##print(dim(elevs))
  ## Get julian days
  ##print(T)
  Jday <- as.POSIXlt(strptime(T, format="%Y%m%d"))$yday + 1  ## Julian day
  Jarray <- array(rep(Jday,each = dim(Tt)[1]*dim(Tt)[2]), dim(Tt))	

  ## efficiency of radiation use at different temperature ranges
  etab <- data.frame(Tmn = c(-20,4.5,9.5,19,25,40), 
                     Tmx = c(4.5,9.5,19,25,40,100),
                     aT = c(0,-0.00893,0.00349,0.01750,0.04680,0),
                     bT = c(0,0.00204,0.00070,0,-0.00117,0) )

  ## Efficiency of conversion of radiation energy to herbage energy
  TtX <- cut(Tt, breaks = unique(c(etab$Tmn, etab$Tmx)), labels = F)
  ##print(length(TtX))
  e <- array(etab$aT[TtX] + (etab$bT[TtX] * Tt) , dim(Tt))
  ##print(dim(e))
  ## Seasonality in growth potential
  x <- array(ifelse(Jarray < 166, 1, 0.4), dim(Jarray))

  ## Radiation and temperature dependent daily growth rate (kg Ha-1 day-1)
  RUE <- (e * 1/Q) *1000
  preRUE <- RUE
  ##print(dim(RUE))
  ## CO2 fertilisation effect
  if(FCO2 == TRUE & cconc >350){
    print('Doing CO2 fertilisation')
    if(cconc > 750){
      RUE <- RUE + RUE * (750/350 - 1) * 0.333 }
    else if(cconc <750){
      RUE <- RUE + RUE * (cconc/350 - 1) * 0.333 }
  }
  ##print(dim(RUE))
  ##Yp <- (e *  * x)/Q ## Note conversion from m-2 to Ha
  # 'Potential' yield, if PET == AET, I think?
  Yp <-(RUE/1000 * ((Rr*10000)*0.5) * x) ## Note conversion from m-2 to Ha
  
  ## PET via FAO Penman-Monteith
  G <- 0 ## soil heat flux
  u <- wind*0.75 ## wind speed at 2m (convert from 10m)
  lambda <- 2.45 ## latent heat of vaporisation

  Eo <- (0.6108 * exp((17.27*Tt) / (Tt + 237.3))) ## Saturation vapour pressure
  if(datasetname %in% 'era5'){
    Ea <- (0.6108 * exp((17.27*Tdp) / (Tdp + 237.3))) } ## Actual vapour pressure = Saturation vapour pressure at dew point
   else {
    Ea <- (rh/100) * Eo }   
  VPD <- Eo - Ea ## Vapour pressure deficit
  delta <- (4098 * Eo) / (Tt + 237.3)^2 ## slope of saturation vapour pressure-temperature curve
  if (is.null(sfcP)){
    print('Using elevfile')
    P <- 101.3*((293-0.0065*array(elevs, dim(Tt))) / 293)^5.26 ## atmospheric pressure
  }
  if (is.null(elevfile)){
    print('Using sfcP')
    P <- sfcP/1000.
  }
  pgamma <- ((1.013 *10^-3) * P)/(0.622*lambda) ## Psychrometric constant

  PET <- ((0.408*delta*(Rr-G)) + (pgamma*(900/(Tt + 273))* u *VPD))/(delta + pgamma*(1 + 0.34 * u))

  ## Soil moisture deficit
  AET <- PET
  SMD <- array(-10, dim(PET))
  for(i in 2:dim(PET)[3]){
    AET[,,i] <- PET[,,i] * ifelse(SMD[,,i-1] > SMDc, ((SMDmax - SMD[,,i-1]) / (SMDmax - SMDc)), 1)
    SMD[,,i] <- ifelse(SMD[,,i-1] + AET[,,i-1] - prec[,,i-1] > -10, SMD[,,i-1] + AET[,,i-1] - prec[,,i-1], -10)
  }

  ##plot(1:dim(Tt)[3], prec[60, 30,], type = "l", col = "blue", ylim = c(-10,100))
  ##points(1:dim(Tt)[3], PET[60, 30,], type ="l", col = "orange")
  ##points(1:dim(Tt)[3], AET[60, 30,], type ="l", col = "red")
  ##points(1:dim(Tt)[3], SMD[60, 30,], type ="l", col = "brown")

  ## Translate to effect on dry matter growth per day (kg ha-1 day-1)
  EaEp <- AET/PET
  YaYp <- 0.2 + (0.8*EaEp)
  Ya <- Yp * YaYp 
  
  ##plot(Yp[60,30,], type  ="l")
  ##points(Ya[60,30,], type  ="l", col = "green")

  ## Total yield per year (convert to T per Ha) 
  YaSum <- apply(Ya, 1:2, sum, na.rm = T)/1000

  ## Convert final yield to raster
  #WLPY = list()
  #WLPY$x= X
  #WLPY$y= Y
  #WLPY$z= YaSum
  #WLPYrast <- raster(WLPY, crs = "+init=epsg:27700")
  ##plot(WLPYrast)

  ##if(is.null(savepath) == FALSE){
  ##  writeRaster(WLPYrast, savepath,"GTiff", overwrite = TRUE)
  ##}
  datalist <- list('PET'=PET, 'AET'=AET, 'SMD'=SMD, 'Yp'=Yp, 'Ya'=Ya, 'YaSum'=YaSum)
  return(datalist)
}

grass_py_point <- function(Tt, Ttmx, Ttmn, prec, Rr,  rh, wind, X, Y, T, datasetname, precipname = 'None', radname = 'None', cconc = NULL, sfcP = NULL, elevs = NULL, FCO2 = TRUE, Q = 18.81, SMDmax = 110, SMDc = 10, savepath = NULL) {

  if (datasetname %in% "era5"){
    Tdp <- rh
  }
  
  ## Set weather datasets for testing
  ##Tt <- temp_tas ## Mean daily temperature data (degrees C)
  ##Ttmn <- temp_tasmin
  ##Ttmx <-  temp_tasmax
  ##Rr <- temp_rss*0.0036*24 ## shortwave radiation data (MJ m-2 day-1) 
  ##prec <- temp_pr ## daily rainfall (mm)
  ##rh <- temp_hurs ## Relative humidity (%)
  ##wind <- temp_wind ## Wind speed 
  ## Q <- 18.81 ## energy content of herbage dry matter (18.81 MJ kg-1) 
  ## SMDmax <- 110 
  ## SMDc <- 10
  ## cconc = temp_cconc

  ## units conversion
  if (datasetname %in% "ukcp18" | grepl('ceres', radname, fixed=TRUE)) {
    print('Converting rad units from W/m^2 --> MJ/m^2/day')
    Rr <- Rr*0.0036*24 ## W/m^2 --> MJ/m^2/day
  }

  ## These conversions only necessary for era5 data
  if (datasetname %in% "era5") {
    Tt <- Tt - 273.15
    Ttmx <- Ttmx - 273.15
    Ttmn <- Ttmn - 273.15
    Tdp <- Tdp - 273.15 ## all Kelvin --> Celsius
    if (!(precipname %in% 'aphrodite')) { ## only if we're not using aphrodite for precip
      prec <- prec*1000. ## m --> mm
    }
    if (!(grepl('ceres', radname, fixed=TRUE))){ ## only if we're not using ceres for solar radiation
      Rr <- Rr/1000000
    }
  }
  
  ## Get julian days
  ##print(T)
  Jarray <- array(as.POSIXlt(strptime(T, format="%Y%m%d"))$yday + 1)  ## Julian day
  ##print('Jarray')
  ##print(dim(Jarray))

  ## efficiency of radiation use at different temperature ranges
  etab <- data.frame(Tmn = c(-20,4.5,9.5,19,25,40), 
                     Tmx = c(4.5,9.5,19,25,40,100),
                     aT = c(0,0.00893,0.00349,0.01750,0.04680,0),
                     bT = c(0,0.00204,0.00070,0,-0.00117,0) )

  ## Efficiency of conversion of radiation energy to herbage energy
  TtX <- array(cut(Tt, breaks = unique(c(etab$Tmn, etab$Tmx)), labels = F))
  ##print('TtX')
  ##print(dim(TtX))
  e <- array(etab$aT[TtX] + (etab$bT[TtX] * Tt) , dim(Tt))
  ##print('e')
  ##print(dim(e))
  
  ## Seasonality in growth potential
  x <- array(ifelse(Jarray < 166, 1, 0.4), dim(Jarray))
  ##print('x')
  ##print(dim(x))
  ## Radiation and temperature dependent daily growth rate (kg Ha-1 day-1)
  RUE <- (e * 1/Q) *1000
  ##print('RUE')
  ##print(dim(RUE))
  ##print(RUE)
  ## CO2 fertilisation effect
  if(FCO2 == TRUE & cconc >350){
    if(cconc > 750){
      RUE <- RUE + RUE * (750/350 - 1) * 0.333 }
    else if(cconc <750){
      RUE <- RUE + RUE * (cconc/350 - 1) * 0.333 }
  }
  ##print(dim(RUE))
  ##print(RUE)

  ##Yp <- (e *  * x)/Q ## Note conversion from m-2 to Ha	
  Yp <-(RUE/1000 * ((Rr*10000)*0.5) * x) ## Note conversion from m-2 to Ha
  ##print('Yp')
  ##print(dim(Yp))
  
  ## PET via FAO Penman-Monteith
  G <- 0 ## soil heat flux
  u <- wind*0.75 ## wind speed at 2m (convert from 10m)
  lambda <- 2.45 ## latent heat of vaporisation

  Eo <- (0.6108 * exp((17.27*Tt) / (Tt + 273.3))) ## Saturation vapour pressure
  ##print('Eo')
  ##print(dim(Eo))
  if(datasetname=='ukcp18'){
    Ea <- (rh/100) * Eo ## Actual vapour pressure
  }
  if(datasetname %in% 'era5'){
    Ea <- (0.6108 * exp((17.27*Tdp) / (Tdp + 237.3))) ## Actual vapour pressure = Saturation vapour pressure at dew point
  }
  ##print('Ea')
  ##print(dim(Ea))
  VPD <- Eo - Ea ## Vapout pressure deficit
  delta <- (4098 * Eo) / (Tt + 273.3)^2 ## slope of saturation vapour pressure-temperature curve
  ##print('delta')
  ##print(dim(delta))
  if (!is.null(elevs)){
    print('Using elevs')
    P <- 101.3*((293-0.0065*array(elevs, dim(Tt))) / 293)^5.26 ## atmospheric pressure
  }
  if (!is.null(sfcP)){
    print('Using sfcP')
    P <- sfcP/1000.
  }
  ##print('P')
  ##print(dim(P))
  pgamma <- ((1.013 *10^-3) * P)/(0.622*lambda) ## Psychrometric constant
  ##print('pgamma')
  ##print(dim(pgamma))

  PET <- ((0.408*delta*(Rr-G)) + (pgamma*(900/(Tt + 273))* u *VPD))/(delta + pgamma*(1 + 0.34 * u))
  ##print('PET')
  ##print(dim(PET))
  
  ## Soil moisture deficit
  AET <- PET
  SMD <- array(-10, dim(PET))
  ##print('SMD')
  ##print(dim(SMD))
  for(i in 2:dim(PET)[1]){
    AET[i] <- PET[i] * ifelse(SMD[i-1] > SMDc, ((SMDmax - SMD[i-1]) / (SMDmax - SMDc)), 1)
    SMD[i] <- ifelse(SMD[i-1] + AET[i-1] - prec[i-1] > -10, SMD[i-1] + AET[i-1] - prec[i-1], -10)
  }
  ##print('AET')
  ##print(dim(AET))
  ##print('SMD')
  ##print(dim(SMD))

  ##plot(1:dim(Tt)[3], prec[60, 30,], type = "l", col = "blue", ylim = c(-10,100))
  ##points(1:dim(Tt)[3], PET[60, 30,], type ="l", col = "orange")
  ##points(1:dim(Tt)[3], AET[60, 30,], type ="l", col = "red")
  ##points(1:dim(Tt)[3], SMD[60, 30,], type ="l", col = "brown")

  ## Translate to effect on dry matter growth per day (kg ha-1 day-1)
  EaEp <- AET/PET
  YaYp <- 0.2 + (0.8*EaEp)
  ##print('YaYp')
  ##print(dim(YaYp))
  Ya <- Yp * YaYp 
  ##print('Ya')
  ##print(dim(Ya))
  
  ##plot(Yp[60,30,], type  ="l")
  ##points(Ya[60,30,], type  ="l", col = "green")

  ## Total yield per year (convert to T per Ha) 
  YaSum <- sum(Ya, na.rm=TRUE)/1000
  ##print('YaSum')
  ##print(dim(YaSum))
  
  ##if(is.null(savepath) == FALSE){
  ##  writeRaster(WLPYrast, savepath,"GTiff", overwrite = TRUE)
  ##}
  datalist <- list('PET'=PET, 'AET'=AET, 'SMD'=SMD, 'Yp'=Yp, 'Ya'=Ya, 'YaSum'=YaSum)
  return(datalist)
}


## Run for each 
##yrs <- 1982:2079
##sms <- c("01","04","05","06","07","08","09","10","11","12","13","15")##
##
##allc <- expand.grid(yrs, sms)
##
##findgyear <- function(xn, ggyear) {
##	as.numeric(substr(dimnames(xn)[[3]], 1,4)) == ggyear
##}
##
##for(XC in 536:nrow(allc)){
##	## Find and load required files and subset to year
##	grange <- rownames(dateranges)[dateranges[,1] <= allc[XC,1] & dateranges[,2] >= allc[XC,1]][1]
##	gfiles <- rlist[grep(grange,rlist)]
##	gfiles <- gfiles[grep(paste("_", allc[XC,2], sep = ""), gfiles)]
##
##	lapply(gfiles, load, environment())
##	temp_tas <- tas_UKCP18[,,findgyear(tas_UKCP18, allc[XC,1])]
##	temp_rls <- rls_UKCP18[,,findgyear(rls_UKCP18, allc[XC,1])]
##	temp_rss <- rss_UKCP18[,,findgyear(rss_UKCP18, allc[XC,1])]
##	temp_pr <- pr_UKCP18[,,findgyear(pr_UKCP18, allc[XC,1])]
##	temp_tasmax <- tasmax_UKCP18[,,findgyear(tasmax_UKCP18, allc[XC,1])]
##	temp_tasmin <- tasmin_UKCP18[,,findgyear(tasmin_UKCP18, allc[XC,1])]
##	temp_hurs <- hurs_UKCP18[,,findgyear(hurs_UKCP18, allc[XC,1])]
##	temp_wind <- sfcWind_UKCP18[,,findgyear(sfcWind_UKCP18, allc[XC,1])]
##	temp_cconc <- pCO2[pCO2$YEAR == allc[XC,1],2]
##
##	grass_py(temp_tas, temp_tasmax, temp_tasmin, temp_pr, temp_rss, temp_hurs, temp_wind,cconc = temp_cconc, FCO2 = TRUE,
##		 savepath = paste("P:\\NEC07148_Crop-Net\\Workfiles\\Redhead crop modelling work\\UKCP18 Yield\\CGrass", allc##[XC,1],allc[XC,2], sep = "_"))
##}
##
##
##
##
##
##
##
##
##
##
##
########################################## OLD PET CALCULATIONS
####
##
##	## Set up soil moisture array
##	sarry <- array(NA, dim(Yp))
##	sgrid <- expand.grid(1:dim(sarry)[1], 1:dim(sarry)[2])
##
##	## Set up PET dataframe 
##	XX <- rep(0,length(dimnames(Tt)[[3]]))
##	pmdata <- data.frame(Year = XX, Month = XX, Day = XX, Tmax = XX, Tmin = XX, RHmax = XX, RHmin = XX, Rs = XX, uz = XX))##
##
##	## Loop over cells
##	for(test in 1:(dim(sarry)[1]*dim(sarry)[2])){
##		
##		testx <- sgrid[test,1]
##		testy <- sgrid[test,2]
##
##		if(sum(is.na(Ttmx[testx, testy,])) < 10){
##		## Set data for PET model
##		pmdata$Year = as.numeric(substr(dimnames(Tt)[[3]],1,4)) 
##		pmdata$Month = as.numeric(substr(dimnames(Tt)[[3]],5,6)) 
##		pmdata$Day = as.numeric(substr(dimnames(Tt)[[3]],7,8))
##		pmdata$Tmax = Ttmx[testx, testy,]
##		pmdata$Tmin = Ttmn[testx, testy,] 
##		pmdata$RHmax = rh[testx, testy,]
##		pmdata$RHmin = rh[testx, testy,]
##		pmdata$Rs = Rr[testx, testy,]
##		pmdata$uz = wind[testx, testy,]
##		##pmdata <- pmdata[dimnames(Tt)[[3]] != "20000230",] 
##
##		## Set constants for PET model
##		pmc <- pmc1
##		pmc$Elev <- elevs[testx, testy]
##		pmc$lat_rad <- latsrad[testx, testy]
##			
##
##		## Prep inputs for PET model
##		pin <- ReadInputs(c("Tmax","Tmin","RHmax","RHmin","Rs","uz"), pmdata, pmc, stopmissing=c(10,10,3), timestep = ##"daily", message = "no")
##
##		## Run PET model
##		PET <- as.numeric(ET.PenmanMonteith(pin, pmc, ts="daily", solar= "data", 
##			wind= "yes", crop="short", message="yes", AdditionalStats="no", save.csv="no")$ET.Daily)
##		PET <- Reduce(function(x,v) ifelse(is.na(v) == TRUE, x, v), PET, accumulate = TRUE)
##
##		## Soil moisture deficit
##		tws <- mapply(function(x,y,z,w) c(x,y,z,w), PET, prec[testx, testy,], 0, -10, SIMPLIFY = FALSE)
##		S <- Reduce(f = function(x, v){
##				
##			aet <- v[1] * ifelse(x[4] > SMDc, ((SMDmax - x[4]) / (SMDmax - SMDc)), 1)
##			Sdef <- ifelse(x[4] + aet - x[2] > -10, x[4] + aet - x[2], -10) 
##			c(v[1],v[2], aet, Sdef)
##		}, x = tws, accumulate = TRUE, init = c(0,0,0, -10))[-1]
##		AET <- data.frame(matrix(unlist(S), nrow=length(S), byrow=T))[,3]
##		S <- data.frame(matrix(unlist(S), nrow=length(S), byrow=T))[,4]
##		sarry[testx,testy,] <- S
##









