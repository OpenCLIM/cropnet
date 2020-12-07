import sys
sys.path.insert(1, '/users/sgsys/matbro/cropNET/code')
from utils import *
#import sys
#import time
#import requests
#import urllib.request
import warnings
warnings.filterwarnings('ignore')
#import rpy2.robjects as robjects
#from rpy2.robjects.packages import importr
#import rpy2.robjects.packages as rpackages
#from rpy2.robjects.vectors import StrVector
from rpy2.robjects import numpy2ri
import numpy as np
import xarray as xr
#import matplotlib.pyplot as plt
#import scipy.optimize as spo
#import datetime as dt
import pandas as pd
#import geopandas as gpd
#import pyproj
#import netCDF4 as nc4
#import cftime as cft
#import os
#import glob
#import shutil
numpy2ri.activate()
#from memory_profiler import profile


##################################################################################
# Edit these:
#
##################################################################################

# UKCP18 ensemble members to use. Note that '02', '03' and '14' don't exist.
# Note that because the error in the modelled value is based on the ensemble
# spread, reducing the number of ensemble members used will result in the
# assimilated/optimum GAI being closer to the original modelled GAI.
# To counteract this, change the moderrinfl variable below. You may need to
# experiment a bit to get the results you want!
ensmems = ['01', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '15']

# years to run for, from the first specified year to the year before the second
years = np.arange(2016, 2020)

# location of the R-formatted driving data                                                                                    
# If this needs to be generated set 'genRdata' to 1. This will take a few hours,
# and will use the UKCP18 data in UKCP18dir. 
# The UKCP18 data must follow a certain format at the moment, see README for details.
Rdatpath = '/users/sgsys/matbro/cropNET/data/wheat_driving'
# location of UKCP18 netcdf files. Must encompass the range specified by years,
# plus a year in front if the growing season does not start on 1st Jan
UKCP18dir = '/data/UKCP18/RCM_12km/daily_timeseries/downloaded_data'
genRdata = 0

# location and filename of the csv file containing the CO2 concentration projections
CO2file = '/users/sgsys/matbro/cropNET/data/UKCP18_CO2_RCP85.csv'

# location to save outputted yields in
saveloc = '/users/sgsys/matbro/cropNET/outputs/'

# switch controlling plotting of outputs.
# may not work outside of a jupyter notebook example...
# 0: no plots, 1: 1 plot containing the timeseries at each iteration, 
# 2: 2 plots, one showing the value of the cost function at each iteration 
plot = 0
# output a plot every interval iterations, if saving the plots, the smaller
# this number is, the slower the code will run. The slowdown when not saving
# plots is smaller, but still significant.
interval = 50
# directory to store plots in. Set to 'None' to not save plots and just display them.
# Saving plots makes the code significantly slower.
# Ignored if plot==0.
plotdir = 'None'


# verify switch. If set to 1, the code will run for all the locations for which we
# have precision yield data, and produce plots comparing the original modelled
# yield to that after assimilation. The yield data is read in from the
# yieldfile.
# As there are ~1300 fields of precision yield data, the code will take a few days
# to run. To reduce this time, change the number of precision yield locations we
# run for in howmanyobs below.
verify = 0
howmanyobs = -1 # set to -1 for all
yieldfile='/users/sgsys/matbro/data/cropyield/Mean_wheat_yields_v2/mean_wheat_yields_OSGB.shp'

# obs switch. Set to 1 to download MODIS LAI, 2 to use already downloaded data,
# otherwise fake obs will be generated, options 0 and 1 use the obscoords specified below
# unless verify==1
obs = 1

# directory containing MODIS LAI data, or data to download it to
# Ignored if obs==0
# Note that all files in this directory will be used, regardless
# of whether they are newly downloaded or not.
# Create separate directories for different sets of locations
MODISdir = '/users/sgsys/matbro/cropNET/data/MODIS/test'

# Coordinates to run over. Alternatively, a csv file with one x,y coord per line can
# be specfied in coordsfile. Set coordsfile to None to use obscoords.
# Both ignored if verify==1
#             x     , y
obscoords = [[580000, 300000], # E Anglia
             [460000, 420000], # C. Eng
             [260000,  90000]] # SW. Eng
coordsfile = None

# Assimilation controls:
obserrtype = 0 # if 1, obsstd = obs*obserr for each tstep (relative)
               # if 0, obsstd = obserr for each tstep (actual)
obserr = 0.1 # the obs stdev to assume (either relative or actual)
moderrinfl = 1 # factor to multiply the model error by

# maximum value by which the optimum/assimilated timeseries should not change by more than
# for one timestep to the next, for the smoothing constraint.
# See README for more details
tsvar=0.25
# To smooth the timeseries based on 1st order or 2nd, 3rd... order differences.
# 1st order seems to work best
order=1
# The strength of the smoothing. The higher this is, the smoother the timeseries
# will be that the smoothing component of the cost function is trying to push the
# optimum/assimilated timeseries to. 
power=10

##########################################################################################
# Main code:
#
##########################################################################################

r['source']('Lynch_potpredict_v2_MJB.R')
fddw = r['formatdrivingdatawheat']
loaddata = r['loaddata']
#GAIfunc = r['GAI']
GAIfunc_point = r['GAI_point']

if genRdata == 1:
    fddw(UKCP18dir, Rdatpath, 'projection_x_coordinate', 'projection_y_coordinate')

# use yield shape file to generate locations to request MODIS data for if verify mode
if verify==1 and obs==1:
    datacodes = MODIS_request(yieldfile, obscoords, years[0]-1, years[-1], MODIScode='MCD15A2H')
    MODIS_download(MODISdir, datacodes, product_name='Lai_500m')

# otherwise use provided x,y coords
if verify==0 and obs==1:
    datacodes = MODIS_request(coordsfile, obscoords, years[0]-1, years[-1], MODIScode='MCD15A2H')
    MODIS_download(MODISdir, datacodes, product_name='Lai_500m')

if obs==0:
    obspdall = fakeobsgen(obscoords, moddata)
else:
    obspdall, obscoords = MODIS_process(os.path.join(MODISdir, '*.csv*'))

if verify==1:
    #obspdall = obspdall.iloc[:, :howmanyobs]
    #obscoords = obscoords[:howmanyobs]
    obspdall = obspdall.iloc[:, :howmanyobs]
    obscoords = obscoords[:howmanyobs]

cfyears_vari = []
# create xr datasets to store data
ensmemsint = [int(ensmem) for ensmem in ensmems]
oldyields_dict = {}
newyields_dict = {}
for tob in obscoords:
    oldyields_dict[str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((len(ensmems), len(years))), coords=[ensmemsint, years], dims=['ensmem', 'year'])
    newyields_dict[str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((len(ensmems), len(years))), coords=[ensmemsint, years], dims=['ensmem', 'year'])
oldyields_vari = xr.Dataset(oldyields_dict)
newyields_vari = xr.Dataset(newyields_dict)
for year in years:

    # extract out a year's worth of obs
    gyeardates = xr.cftime_range(str(year-1)+'-10-01', str(year)+'-09-30', calendar='360_day', freq='D', name='date').values
    obspd = obspdall.loc[gyeardates]
    
    if plotdir != 'None':
        if not os.path.exists(plotdir):
            os.makedirs(plotdir)

    GAI_p_all, GAI_p_ensmean, GAI_p_ensstd, \
    tmean_p_all, prec_p_all, solarrad_p_all, \
    Jarray_p_all, Cday_p_all, GSS_p_all, \
    HarvestJday, AWC_allp, temp_cconc = ensgen_point(ensmems, year, Rdatpath, CO2file, obscoords)

    # now we've done the ensgen, do the assimilation
    cfxrs_vari = []
    for ensmem in ensmems:
        print('Doing DA for GAI using variational method for ensmem ' + ensmem)
        mergedall_vari, cfall_vari = vari_method_point(obspd, GAI_p_all.sel(ensmem=int(ensmem)), GAI_p_ensstd, 
                                                       obscoords, year, ensmem, obserrtype, obserr, moderrinfl,
                                                       tsvar, order, power, plot, plotdir, interval)

        cfxr_vari = cfall_vari.to_xarray()
        del cfall_vari
        cfxr_vari = cfxr_vari.expand_dims({'ensmem': [int(ensmem)]})
        cfxr_vari = xr.where(np.isinf(cfxr_vari), np.nan, cfxr_vari) # set infs as nans
        cfxrs_vari.append(cfxr_vari)

        oldyields, newyields = update_yield_points_point(GAI_p_all.sel(ensmem=int(ensmem)), mergedall_vari, obscoords, 
                                         tmean_p_all.sel(ensmem=int(ensmem)), prec_p_all.sel(ensmem=int(ensmem)), 
                                         solarrad_p_all.sel(ensmem=int(ensmem)), Jarray_p_all.sel(ensmem=int(ensmem)), 
                                         Cday_p_all.sel(ensmem=int(ensmem)), GSS_p_all.sel(ensmem=int(ensmem)), 
                                         HarvestJday, AWC_allp, temp_cconc, ensmem)

        del mergedall_vari
        counter=0
        for tob in obscoords:
            oldyields_vari[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmem), year=year)] = oldyields[counter]
            newyields_vari[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmem), year=year)] = newyields[counter]
            counter+=1

    cfallensmems_vari = xr.concat(cfxrs_vari, dim='ensmem')
    cfyears_vari.append(cfallensmems_vari)

print('Calculating climatological conversion factor using years: ' + str(years) + '\n')
cfallyears_vari = xr.concat(cfyears_vari, dim='date')
del cfyears_vari
cfensmean_vari = cfallyears_vari.mean(dim='ensmem')
del cfallyears_vari
cfclim_vari = cfensmean_vari.groupby('date.dayofyear').mean().roll(dayofyear=-270, roll_coords=True)

print('Saving to netcdf')
if not os.path.exists(saveloc):
    os.makedirs(saveloc)
now = dt.datetime.now().strftime('%Y%m%d%H%M%s')[:14]
oldyieldfile = os.path.join(saveloc, now + '_oldyield.nc')
newyieldfile = os.path.join(saveloc, now + '_newyield.nc')
cfclimfile   = os.path.join(saveloc, now + '_cfclimo.nc')
oldyields_vari.to_netcdf(oldyieldfile)
newyields_vari.to_netcdf(newyieldfile)
cfclim_vari.to_netcdf(cfclimfile)

if verify==1:
    verifyyield(yieldfile, obscoords, oldyields_vari, newyields_vari, years, saveloc)
