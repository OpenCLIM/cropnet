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
Wrapper to run the CropNet crop models on the JASMIN LOTUS HPC cluster

Input options:
datalocbase ----- Location of all the meteorological driving data files. 
                  With subfolders within for the different rcp/ensmem permutations.
outlocbase ------ Where to store the final yield outputs of the model. Subfolders based
                  on the crop, ensmem and basedataset will be created under here.
AWCrast --------- File path of the available water content file as a geotiff raster or netcdf file.
CO2filesloc ----- Location of the folder containing the CO2 concentrations files as csv files with two columns,
                  the year and the CO2 conc in ppm. Set to 'None' (including quotes) to 
                  disable the CO2 fertlisation option.
simx ------------ Only used if basedatasetname is 'ukcp18'. Ignored otherwise. Identifies the
                  ensemble member to use by selecting out the driving data files with the 
                  specified simx in the their name
crop ------------ Crop to model. Can be 'wheat', 'grass' or 'OSR'
output_biomass -- Whether or not to output the crop biomass evolution in addition to the yields.
                  This significantly increases the disk size of the output as it is a 3D variable
                  (time, y, x). Can be True or False.
basedatasetname - Which driving dataset to use, along with possible rcp and ensmem identifiers.
                  This is used to determine the file and variable names to look for in the 
                  data read-in routine. These can be configured if necessary by editing the 'getnames' function directly.
                  Pre-configured options are 'ukcp18', 'era5', 'chess_and_haduk',
                  'chess-scape_8.5_01', 'chess-scape_8.5_04', 'chess-scape_8.5_06', 'chess-scape_8.5_15', 
                  'chess-scape_2.6_01', 'chess-scape_2.6_04', 'chess-scape_2.6_06', 'chess-scape_2.6_15'.
'''

datalocbase = "/gws/nopw/j04/ceh_generic/matbro/cropnet/drivingdata"
outlocbase  = "/gws/nopw/j04/ceh_generic/matbro/cropnet/outputs/no_assim/"
CO2filesloc = '/gws/nopw/j04/ceh_generic/matbro/cropnet/openclim-cropnet/CO2_files/'
AWCrast = "/gws/nopw/j04/ceh_generic/matbro/cropnet/data/AWC/G2G_derived_AWC/G2G_AWC_UK.nc"
elevfile = '/gws/nopw/j04/ceh_generic/matbro/cropnet/data/NextMap_DTM_50m.tif'
CO2 = sys.argv[3] # '8.5', '8.5_01', '8.5_04', '8.5_06', '8.5_15', '2.6' or 'None'
simx = '01'
crop= sys.argv[1] # 'wheat', 'grass' or 'OSR'
basedatasetname = sys.argv[2] 
# 'ukcp18', 'chess-scape_8.5_01', 'chess-scape_8.5_04', 'chess-scape_8.5_06', 'chess-scape_8.5_15', 
# 'chess-scape_2.6_01', 'chess-scape_2.6_04', 'chess-scape_2.6_06', 'chess-scape_2.6_15', 'era5', 'chess_and_haduk'
output_biomass = False # True or False

startyear=int(sys.argv[4])
startmonth=int(sys.argv[5])
startday=1

if CO2 == '8.5':
    CO2file = os.path.join(CO2filesloc, 'CHESS-SCAPE_RCP85_01.csv')
elif CO2 == '8.5_01':
    CO2file = os.path.join(CO2filesloc, 'CHESS-SCAPE_RCP85_01.csv')
elif CO2 == '8.5_04':
    CO2file = os.path.join(CO2filesloc, 'CHESS-SCAPE_RCP85_04.csv')
elif CO2 == '8.5_06':
    CO2file = os.path.join(CO2filesloc, 'CHESS-SCAPE_RCP85_06.csv')
elif CO2 == '8.5_15':
    CO2file = os.path.join(CO2filesloc, 'CHESS-SCAPE_RCP85_15.csv')
elif CO2 == '2.6':
    CO2file = os.path.join(CO2filesloc, 'CHESS-SCAPE_RCP26.csv')
else:
    print('Unrecognised or no CO2 option specified, running without CO2 fertilisation')
    CO2file = 'None'

if basedatasetname == 'ukcp18':
    insubfolder = 'ukcp18/**/*.nc'
    outsubfolder = crop + '/ukcp18'
    dataloc = os.path.join(datalocbase, insubfolder)
    outloc = os.path.join(outlocbase, outsubfolder)
elif basedatasetname == 'chess-scape_8.5_01':
    insubfolder = 'chess-scape/rcp85/01/**/*.nc'
    if CO2file == 'None':
        outsubfolder = crop + '/chess-scape/rcp85_noCO2/01'
    else:
        outsubfolder = crop + '/chess-scape/rcp85/01'
    dataloc = os.path.join(datalocbase, insubfolder)
    outloc = os.path.join(outlocbase, outsubfolder)
elif basedatasetname == 'chess-scape_8.5_04':
    insubfolder = 'chess-scape/rcp85/04/**/*.nc'
    if CO2file == 'None':
        outsubfolder = crop + '/chess-scape/rcp85_noCO2/04'
    else:
        outsubfolder = crop + '/chess-scape/rcp85/04'
    dataloc = os.path.join(datalocbase, insubfolder)
    outloc = os.path.join(outlocbase, outsubfolder)
elif basedatasetname == 'chess-scape_8.5_06':
    insubfolder = 'chess-scape/rcp85/06/**/*.nc'
    if CO2file == 'None':
        outsubfolder = crop + '/chess-scape/rcp85_noCO2/06'
    else:
        outsubfolder = crop + '/chess-scape/rcp85/06'
    dataloc = os.path.join(datalocbase, insubfolder)
    outloc = os.path.join(outlocbase, outsubfolder)
elif basedatasetname == 'chess-scape_8.5_15':
    insubfolder = 'chess-scape/rcp85/15/**/*.nc'
    if CO2file == 'None':
        outsubfolder = crop + '/chess-scape/rcp85_noCO2/15'
    else:
        outsubfolder = crop + '/chess-scape/rcp85/15'
    dataloc = os.path.join(datalocbase, insubfolder)
    outloc = os.path.join(outlocbase, outsubfolder)
elif basedatasetname == 'chess-scape_2.6_01':
    insubfolder = 'chess-scape/rcp26/01/**/*.nc'
    if CO2file == 'None':
        outsubfolder = crop + '/chess-scape/rcp26_noCO2/01'
    else:
        outsubfolder = crop + '/chess-scape/rcp26/01'
    dataloc = os.path.join(datalocbase, insubfolder)
    outloc = os.path.join(outlocbase, outsubfolder)
elif basedatasetname == 'chess-scape_2.6_04':
    insubfolder = 'chess-scape/rcp26/04/**/*.nc'
    if CO2file == 'None':
        outsubfolder = crop + '/chess-scape/rcp26_noCO2/04'
    else:
        outsubfolder = crop + '/chess-scape/rcp26/04'
    dataloc = os.path.join(datalocbase, insubfolder)
    outloc = os.path.join(outlocbase, outsubfolder)
elif basedatasetname == 'chess-scape_2.6_06':
    insubfolder = 'chess-scape/rcp26/06/**/*.nc'
    if CO2file == 'None':
        outsubfolder = crop + '/chess-scape/rcp26_noCO2/06'
    else:
        outsubfolder = crop + '/chess-scape/rcp26/06'
    dataloc = os.path.join(datalocbase, insubfolder)
    outloc = os.path.join(outlocbase, outsubfolder)
elif basedatasetname == 'chess-scape_2.6_15':
    insubfolder = 'chess-scape/rcp26/15/**/*.nc'
    if CO2file == 'None':
        outsubfolder = crop + '/chess-scape/rcp26_noCO2/15'
    else:
        outsubfolder = crop + '/chess-scape/rcp26/15'
    dataloc = os.path.join(datalocbase, insubfolder)
    outloc = os.path.join(outlocbase, outsubfolder)
elif basedatasetname == 'era5':
    insubfolder = 'era5/**/*.nc'
    outsubfolder = 'era5'
    dataloc = os.path.join(datalocbase, insubfolder)
    outloc = os.path.join(outlocbase, outsubfolder)
elif basedatasetname == 'chess_and_haduk':
    insubfolder = 'chess_and_haduk/**/*.nc'
    outsubfolder = 'chess_and_haduk'
    dataloc = os.path.join(datalocbase, insubfolder)
    outloc = os.path.join(outlocbase, outsubfolder)

print('################## RUN INFO ###################')
print('Crop: ' + crop)
print('Dataset: ' + basedatasetname)
print('CO2 Scenario: ' + CO2)
print('Starting year: ' + str(startyear))
print('Starting month: ' + str(startmonth))
print('Starting day: ' + str(startday))
print('Data input path: ' + dataloc)
print('Output path: ' + outloc)
print('###############################################')

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
    elif crop == 'OSR':
        r['source']('OSR_potpredict_MJB.R')
    elif crop == 'grass':
        r['source']('Grass_potpredict_MJB.R')
    else:
        SyntaxError('crop must be either wheat, grass or OSR. Currently it is ' + str(crop))
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
    if crop == 'wheat' or crop == 'OSR':
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
        print('CO2 concentration: ' + str(cconc))
    except KeyError:
        print('CO2 fertilisation disabled')
        cconc = 0
        FCO2 = False
    
    # Convert these to lon,lat
    if crop == 'wheat' or crop == 'OSR':
        proj = pyproj.Transformer.from_crs(27700, 4326, always_xy=True)
        xx,yy = np.meshgrid(x,y)
        try:
            lons,lats = proj.transform(xx, yy, errcheck=True)
        except pyproj.exceptions.ProjError:
            lons,lats = proj.transform(xx, yy, errcheck=True)


    # run crop models
    if crop == 'grass':
        grassfunc = r['grass_py']
        print('Running grass model')
        if basedatasetname == 'era5':
            datalist = grassfunc(tmean, tmax, tmin, prec, solarrad, tdp, wind,
                                 x, y, t, datasetname = basedatasetname, cconc=cconc, 
                                 FCO2=FCO2, elevfile=elevfile, sfcP=sfcP)
        else:
            datalist = grassfunc(tmean, tmax, tmin, prec, solarrad, tdp, wind,
                                 x, y, t, basedatasetname, cconc=cconc, FCO2=FCO2, 
                                 elevfile=elevfile)
    
        PET = np.array(datalist.rx2('PET'))
        AET = np.array(datalist.rx2('AET'))
        SMD = np.array(datalist.rx2('SMD'))
        Yp  = np.array(datalist.rx2('Yp'))
        Ya  = np.array(datalist.rx2('Ya'))
        YaSum = np.array(datalist.rx2('YaSum'))
        
        outfile = os.path.join(outloc, 'yields_grass_' + basedatasetname + '_' + str(endyear) + '.nc')
        outputsave(YaSum, [y, x], ['y', 'x'], endyear, 'yield', 'tn/hc', outfile)
        if output_biomass:
            outfile = os.path.join(outloc, 'Ya_grass_' + basedatasetname + '_' + str(endyear) + '.nc')
            outputsave(Ya.transpose(2,0,1), [t, y, x], ['t', 'y', 'x'], endyear, 'Actual_yield', 'tn/hc', outfile)
            outfile = os.path.join(outloc, 'Yp_grass_' + basedatasetname + '_' + str(endyear) + '.nc')
            outputsave(Yp.transpose(2,0,1), [t, y, x], ['t', 'y', 'x'], endyear, 'Potential_yield', 'tn/hc', outfile)
    
    elif crop == 'wheat':
        GAIfunc = r['GAI']
        print('Running model, calculating GAI')
        datalist = GAIfunc(tmean, tmax, tmin, prec, solarrad, x, y, t, lats, datasetname = basedatasetname)
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
        if output_biomass:
            WU_Biomass = np.array(datalist2.rx2('WU_Biomass'))
            WL_Biomass = np.array(datalist2.rx2('WL_Biomass'))
            WU_Biomassxr = xr.DataArray(WU_Biomass, [y, x, t], ['y', 'x', 't'])
            WL_Biomassxr = xr.DataArray(WL_Biomass, [y, x, t], ['y', 'x', 't'])
            WU_Biomassxr.name = 'water_unlimited_biomass'
            WL_Biomassxr.name = 'water_limited_biomass'
            WU_Biomassname = 'UK_WUbiomass_' + basedatasetname + '_' + enddatestr + '_' + str(endyear) + '.nc'
            WL_Biomassname = 'UK_WLbiomass_' + basedatasetname + '_' + enddatestr + '_' + str(endyear) + '.nc'
            WU_Biomassxr.to_netcdf(os.path.join(outloc, WU_Biomassname))
            WL_Biomassxr.to_netcdf(os.path.join(outloc, WL_Biomassname))
    elif crop == 'OSR':
        yieldfunc = r['osr_py']
        print('Any tmean NaNs?:')
        print(np.any(np.isnan(tmean)))
        print('All NaNs?:')
        print(np.all(np.isnan(tmean)))
        datalist2 = yieldfunc(tmean, tmax, tmin, prec, solarrad, 
                              AWC, x, y, t, lats, cconc = cconc, 
                              datasetname = basedatasetname, FCO2=FCO2)
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
        if output_biomass:
            WU_Biomass = np.array(datalist2.rx2('WU_Biomass'))
            WL_Biomass = np.array(datalist2.rx2('WL_Biomass'))
            WU_Biomassxr = xr.DataArray(WU_Biomass, [y, x, t], ['y', 'x', 't'])
            WL_Biomassxr = xr.DataArray(WL_Biomass, [y, x, t], ['y', 'x', 't'])
            WU_Biomassxr.name = 'water_unlimited_biomass'
            WL_Biomassxr.name = 'water_limited_biomass'
            WU_Biomassname = 'UK_WUbiomass_' + basedatasetname + '_' + enddatestr + '_' + str(endyear) + '.nc'
            WL_Biomassname = 'UK_WLbiomass_' + basedatasetname + '_' + enddatestr + '_' + str(endyear) + '.nc'
            WU_Biomassxr.to_netcdf(os.path.join(outloc, WU_Biomassname))
            WL_Biomassxr.to_netcdf(os.path.join(outloc, WL_Biomassname))
