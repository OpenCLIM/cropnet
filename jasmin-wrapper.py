import sys
from utils import *
from MODIS import *
import warnings
warnings.filterwarnings('ignore')
from rpy2.robjects import numpy2ri
import numpy as np
import xarray as xr
from dateutil.relativedelta import relativedelta
import datetime as dt
import cftime as cft
import pyproj
import os
numpy2ri.activate()

'''
Wrapper to run the CropNET/EcoCrop crop models on the JASMIN LOTUS HPC cluster

Input options:
dataloc - Location of all the meteorological driving data files. Use wildcards
          to ensure all files are selected. '**' means 'all folders'
saveloc - Location to store the processed driving data files that are created
          by the data read-in routine
outloc  - Where to store the final yield outputs of the model
AWCrast - File path of the available water content file as a geotiff raster
CO2file - File path of the CO2 concentrations file as a csv file with two columns,
          the year and the CO2 conc in ppm. Set to 'None' (including quotes) to 
          disable the CO2 fertlisation option.
simx    - Only used if basedatasetname is 'ukcp18'. Ignored otherwise. Identifies the
          ensemble member to use by selecting out the driving data files with the 
          specified simx in the their name
crop    - Currently redundant as only wheat crop is supported
basedatasetname - Which driving dataset to use. This is used to determine the file and
                  variable names to look for in the data read-in routine. These can be
                  configured if necessary by editing the 'getnames' function directly.
                  Pre-configured options are 'ukcp18', 'ukcp18bc', 'era5', 'chess_and_haduk'.
'''

dataloc = "/gws/nopw/j04/ceh_generic/matbro/cropnet/drivingdata/chess-scape_RCP26_test/**/*.nc"
outloc  = "/gws/nopw/j04/ceh_generic/matbro/cropnet/outputs/no_assim/ukcp18bc_rcp26"
AWCrast = "/gws/nopw/j04/ceh_generic/matbro/cropnet/data/AWC/G2G_derived_AWC/G2G_AWC_UK.nc"
elevfile = '/gws/nopw/j04/ceh_generic/matbro/cropnet/data/NextMap_DTM_50m.tif'
CO2file = 'None'
simx = '01'
crop='wheat'
basedatasetname = 'ukcp18bc'

startyear=1981#int(sys.argv[1])
startmonth=10#int(sys.argv[2])
startday=1

if crop=='grass':
    AWCrast = None

if basedatasetname == 'chess_and_haduk':
    caltype = 'gregorian'
else:
    caltype = '360_day'

if caltype == '360_day':
    startdate = cft.datetime(startyear,startmonth,startday,calendar='360_day')
    enddate = startdate + dt.timedelta(359)
    fnenddate = startdate + dt.timedelta(360)
elif caltype == 'gregorian':
    startdate = dt.datetime(startyear,startmonth,startday)
    enddate = startdate + relativedelta(years=1) - relativedelta(days=1)
    fnenddate = startdate + relativedelta(years=1)

endyear=enddate.year
endmonth=enddate.month
endday=enddate.day

if not os.path.exists(outloc):
    os.makedirs(outloc)

if startyear == endyear:
    years = [startyear]
elif endyear > startyear+1:
    if endmonth < startmonth:
        years = list(np.arange(startyear+1, endyear+1))[:-1]
    elif endmonth==startmonth and endday < startday:
        years = list(np.arange(startyear+1, endyear+1))[:-1]
    else:
        years = list(np.arange(startyear+1, endyear+1))
else:
    years = list(np.arange(startyear+1, endyear+1))

startmonthstr = '{:02.0f}'.format(startmonth)
startdaystr = '{:02.0f}'.format(startday)
endmonthstr = '{:02.0f}'.format(endmonth)
enddaystr = '{:02.0f}'.format(endday)

for year in years:
    # extract out a year's worth of obs
    if startyear==endyear:
        gyeardates = xr.cftime_range(str(year) + '-' + startmonthstr + '-' + startdaystr,
                                     str(year) + '-' + endmonthstr + '-' + enddaystr,
                                     calendar=caltype, freq='D', name='date').values
    elif len(years)==1 and endmonth < startmonth:
        gyeardates = xr.cftime_range(str(year-1) + '-' + startmonthstr + '-' + startdaystr,
                                     str(year) + '-' + endmonthstr + '-' + enddaystr,
                                     calendar=caltype, freq='D', name='date').values
    elif len(years)==1 and endmonth > startmonth:
        gyeardates = xr.cftime_range(str(year-1) + '-' + startmonthstr + '-' + startdaystr,
                                     str(year) + '-' + startmonthstr + '-' + startdaystr,
                                     calendar=caltype, freq='D', name='date').values
    elif len(years)==1 and endday < startday:
        gyeardates = xr.cftime_range(str(year-1) + '-' + startmonthstr + '-' + startdaystr,
                                     str(year) + '-' + endmonthstr + '-' + enddaystr,
                                     calendar=caltype, freq='D', name='date').values
    elif len(years)==1 and endday >= startday:
        gyeardates = xr.cftime_range(str(year-1) + '-' + startmonthstr + '-' + startdaystr,
                                     str(year) + '-' + startmonthstr + '-' + startdaystr,
                                     calendar=caltype, freq='D', name='date').values
    else:
        gyeardates = xr.cftime_range(str(year-1) + '-' + startmonthstr + '-' + startdaystr,
                                     str(year) + '-' + startmonthstr + '-' + startdaystr,
                                     calendar=caltype, freq='D', name='date').values

    print('Running for ' + str(gyeardates[0]) + ' to ' + str(gyeardates[-1]))

    times = np.array([str(t.year)+'{:02.0f}'.format(t.month)+'{:02.0f}'.format(t.day) for t in list(gyeardates)])

    if crop == 'wheat':
        r['source']('Lynch_potpredict_v2_MJB.R')
    elif crop == 'grass':
        r['source']('Grass_potpredict_MJB.R')
    else:
        SyntaxError('crop must be either wheat or grass, currently it is ' + str(crop))
    dd = load_driving_data(basedatasetname, times,
                           None, dataloc, simx, crop, AWCrast, CO2file)
    
    tmean = dd['tmean']
    tmax  = dd['tmax']
    tmin  = dd['tmin']
    prec  = dd['prec']
    solarrad = dd['solarrad']
    x = dd['x']
    y = dd['y']
    t = dd['t']
    if crop == 'wheat':
        AWC = dd['AWC']
    elif crop == 'grass':
        if basedatasetname == 'era5':
            # era5 has dew point instead of rh
            tdp = dd['tdp']
            # we use sfcP to determine elevation when using era5, instead of a raster file
            sfcP = dd['sfcP']
        else:
            tdp = dd['rh']
            sfcP = 0
        wind  = dd['wind']
    
    FCO2 = True
    try:
        cconc = dd['cconc']
    except KeyError:
        print('CO2 fertilisation disabled')
        cconc = 0
        FCO2 = False
    
    # Convert these to lon,lat
    if crop == 'wheat':
        proj = pyproj.Transformer.from_crs(27700, 4326, always_xy=True)
        xx,yy = np.meshgrid(x,y)
        try:
            lons,lats = proj.transform(xx, yy, errcheck=True)
        except pyproj.exceptions.ProjError:
            lons,lats = proj.transform(xx, yy, errcheck=True)
    
    if crop == 'grass':
        grassfunc = r['grass_py']
        print('Running grass model')
        datalist = grassfunc(tmean, tmax, tmin, prec, solarrad, tdp, wind,
                             x, y, t, basedatasetname, cconc=cconc, FCO2=FCO2, 
                             elevfile=elevfile, sfcP=sfcP)
    
        PET = np.array(datalist.rx2('PET'))
        AET = np.array(datalist.rx2('AET'))
        SMD = np.array(datalist.rx2('SMD'))
        Yp  = np.array(datalist.rx2('Yp'))
        Ya  = np.array(datalist.rx2('Ya'))
        YaSum = np.array(datalist.rx2('YaSum'))
        
        outfile = os.path.join(outloc, 'yields_grass_' + basedatasetname + '_' + str(endyear) + '.nc')
        outputsave(YaSum, [y, x], ['y', 'x'], endyear, 'yield', 'tn/hc', outfile)
        outfile = os.path.join(outloc, 'Ya_grass_' + basedatasetname + '_' + str(endyear) + '.nc')
        outputsave(Ya.transpose(2,0,1), [t, y, x], ['t', 'y', 'x'], endyear, 'Actual_yield', 'tn/hc', outfile)
        outfile = os.path.join(outloc, 'Yp_grass_' + basedatasetname + '_' + str(endyear) + '.nc')
        outputsave(Yp.transpose(2,0,1), [t, y, x], ['t', 'y', 'x'], endyear, 'Potential_yield', 'tn/hc', outfile)
    
    elif crop == 'wheat':
        GAIfunc = r['GAI']
        print('Running model, calculating GAI')
        datalist = GAIfunc(tmean, tmax, tmin, prec, solarrad, x, y, t, lats, basedatasetname)
        GAI = np.array(datalist.rx2('GAI'))
        tmean       = np.array(datalist.rx2('tmean'))
        tmin        = np.array(datalist.rx2('tmin'))
        tmax        = np.array(datalist.rx2('tmax'))
        prec        = np.array(datalist.rx2('prec'))
        solarrad    = np.array(datalist.rx2('solarrad'))
        Jarray      = np.array(datalist.rx2('Jarray'))
        Cday        = np.array(datalist.rx2('Cday'))
        CDD         = np.array(datalist.rx2('CDD'))
        TT          = np.array(datalist.rx2('TT'))
        GSS_r       = datalist.rx2('GSS')
        GSS         = np.array(datalist.rx2('GSS')) # As it's a string array it comes out as a flattened array               
        GSS         = GSS.reshape((GSS_r.dim[2], GSS_r.dim[1], GSS_r.dim[0]))
        GSS         = GSS.transpose(2,1,0) # therefore we need to manually reshape it                                         
        HarvestJday = datalist.rx2('HarvestJday')
        
        yieldfunc = r['wheat_yield']
        print('Calculating yield')
        datalist2 = yieldfunc(GAI, tmean, tmin, tmax, prec, solarrad, AWC, Jarray, Cday, GSS, HarvestJday, CDD, TT, x, y, cconc, FCO2=FCO2)
        WUyield = np.array(datalist2.rx2('WUyield'))
        WLyield = np.array(datalist2.rx2('WLyield'))
        WLHLyield = np.array(datalist2.rx2('WLHLyield'))
        WUHLyield = np.array(datalist2.rx2('WUHLyield'))
        WUyieldxr = xr.DataArray(WUyield, [y, x], ['y', 'x'])
        WLyieldxr = xr.DataArray(WLyield, [y, x], ['y', 'x'])
        WUHLyieldxr = xr.DataArray(WUHLyield, [y, x], ['y', 'x'])
        WLHLyieldxr = xr.DataArray(WLHLyield, [y, x], ['y', 'x'])
        WUyieldxr.name = 'water_unlimited_potential_wheat_yield'
        WLyieldxr.name = 'water_limited_potential_wheat_yield'
        WUHLyieldxr.name = 'water_unlimited_heat_stressed_potential_wheat_yield'
        WLHLyieldxr.name = 'water_limited_heat_stressed_potential_wheat_yield'
        enddatestr = fnenddate.strftime('%b%d')
        WUyieldname   = 'UK_WUpotyield_' + basedatasetname + '_' + enddatestr + '_' + str(endyear) + '.nc'
        WLyieldname   = 'UK_WLpotyield_' + basedatasetname + '_' + enddatestr + '_' + str(endyear) + '.nc'
        WUHLyieldname = 'UK_WUHLpotyield_' + basedatasetname + '_' + enddatestr + '_' + str(endyear) + '.nc'
        WLHLyieldname = 'UK_WLHLpotyield_' + basedatasetname + '_' + enddatestr + '_' + str(endyear) + '.nc'
        WUyieldxr.to_netcdf(os.path.join(outloc, WUyieldname))
        WLyieldxr.to_netcdf(os.path.join(outloc, WLyieldname))
        WUHLyieldxr.to_netcdf(os.path.join(outloc, WUHLyieldname))
        WLHLyieldxr.to_netcdf(os.path.join(outloc, WLHLyieldname))
    
