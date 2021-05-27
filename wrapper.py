import sys
from utils import *
from MODIS import *
import warnings
warnings.filterwarnings('ignore')
from rpy2.robjects import numpy2ri
import numpy as np
import xarray as xr
import pandas as pd
numpy2ri.activate()
#from memory_profiler import profile


##################################################################################
# Edit these:
#
##################################################################################

# name of crop and dataset name to use
# names of alternative precipitation and radiation datasets to use
# set to 'None' to use the precip and radiation from the main dataset
# lastly, name of the variable to assimilate, either fPAR or LAI
crop = 'wheat'
datasetname = 'ukcp18'
precname = 'None'
radname = 'None'
assimvar = 'fPAR'

# location of the driving data netcdf files. Must encompass the time range specified 
# Use wildcards to match all the required data files. '**' means every directory 
dataloc = "/data/UKCP18/RCM_12km/daily_timeseries/downloaded_data/**/*.nc"

# names of the driving data variables to use
# CURRENTLY NOT IMPLEMENTED - TO DO
varnames = ["pr", "tasmax", "tasmin", "rss", "tas"]

# location and filename of elevation data
# only used for grass crop currently
elevfile = 'None'

# location and filename of the csv file containing the CO2 concentration projections
CO2file = '/users/sgsys/matbro/cropnet/data/UKCP18_CO2_RCP85.csv'

# location and filename of the AWC data
AWCrast = "/users/sgsys/matbro/cropnet/data/MaxWet1.tif"

# location to save outputted yields in                                                                                       
yieldsaveloc = '/users/sgsys/matbro/cropnet/outputs/'

# where to store the formatted driving data
saveloc = '/users/sgsys/matbro/cropnet/driving_datafiles'

# switch controlling plotting of outputs.
# may not work outside of a jupyter notebook...
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
verify = 1
howmanyobs = -1 # set to -1 for all
yieldshapefile='/users/sgsys/matbro/cropnet/data/cropyield/Mean_wheat_yields_v3/hayman_etal_val_data.shp'

# obs switch. Set to 1 to download MODIS LAI, 2 to use already downloaded data,                                               
# option 1 uses the coords specified below
obs = 2

# directory containing MODIS LAI data, or data to download it to
# Ignored if obs==0
# Note that all files in this directory will be used, regardless
# of whether they are newly downloaded or not.
# Create separate directories for different sets of locations
MODISdir = '/users/sgsys/matbro/cropnet/data/MODIS/fPAR/'

# Coordinates to run over. Alternatively, a csv file with one x,y coord per line can
# be specfied in coordsfile. Set coordsfile to None to use obscoords.
# Both ignored if verify==1
#          x     , y
coords = [[580000, 300000], # E Anglia
          [460000, 420000], # C. Eng
          [260000,  90000]] # SW. Eng
coordsfile = None

# UKCP18 ensemble members to use. Note that '02', '03' and '14' don't exist.
# Note that because the error in the modelled value is based on the ensemble
# spread, reducing the number of ensemble members used will result in the
# assimilated/optimum GAI being closer to the original modelled GAI.
# To counteract this, change the moderrinfl variable below. You may need to
# experiment a bit to get the results you want!
ensmems = ['01', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '15']
#ensmems = ['01']

# timespan to run for. 
# If the timespan spans more than 1 whole year,
# the end of the timespan will be chopped off to
# make a whole number of years.
startyear = 2015
startmonth = 10
startday = 1
endyear = 2016
endmonth = 10
endday = 1

# Assimilation controls:
obserrtype = 0 # if 1, obsstd = obs*obserr for each tstep (relative)
               # if 0, obsstd = obserr for each tstep (actual)
obserr = 0.025 # the obs stdev to assume (either relative or actual)

moderrtype = 1 # as obserrtype
moderr = 0.5 # as obserr. These two variables only used if len(ensmems)==1
moderrinfl = 3 # factor to multiply the model error by. 

# maximum value by which the optimum/assimilated timeseries should not change by more than
# for one timestep to the next, for the smoothing constraint.
# See README for more details
tsvar=0.1
# To smooth the timeseries based on 1st order or 2nd, 3rd... order differences.
# 1st order seems to work best
order=1
# The strength of the smoothing. The higher this is, the smoother the timeseries
# will be that the smoothing component of the cost function is trying to push the
# optimum/assimilated timeseries to. 
power=20

##########################################################################################
# Main code:
#
##########################################################################################

r['source']('Lynch_potpredict_v2_MJB.R')
#GAIfunc = r['GAI']
GAIfunc_point = r['GAI_point']

if assimvar == 'LAI':
    mname = 'Lai_500m'
elif assimvar == 'fPAR':
    mname = 'Fpar_500m'
# use yield shape file to generate locations to request MODIS data for if verify mode
if verify==1 and obs==1:
    datacodes = MODIS_request(yieldshapefile, coords, startyear, endyear, MODIScode='MCD15A2H')
    MODIS_download(MODISdir, datacodes, product_name=mname)

# otherwise use provided x,y coords
if verify==0 and obs==1:
    datacodes = MODIS_request(coordsfile, coords, startyear, endyear, MODIScode='MCD15A2H')
    MODIS_download(MODISdir, datacodes, product_name=mname)

print('Processing obs')
if assimvar == 'LAI':
    ft = 0.5
elif assimvar == 'fPAR':
    ft = 0
obspdall, coords = MODIS_process(os.path.join(MODISdir, '*.csv*'), filter_threshold=ft)

if verify==1:
    obspdall = obspdall.iloc[:, :howmanyobs]
    coords = coords[:howmanyobs]

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

cfyears_vari = []
post_errs_years = []
# create xr datasets to store data
ensmemsint = [int(ensmem) for ensmem in ensmems]
oldyields_dict = {}
newyields_dict = {}
for tob in coords:
    oldyields_dict[str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((len(ensmems), len(years))), coords=[ensmemsint, years], dims=['ensmem', 'year'])
    newyields_dict[str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((len(ensmems), len(years))), coords=[ensmemsint, years], dims=['ensmem', 'year'])
oldyields_vari = xr.Dataset(oldyields_dict)
newyields_vari = xr.Dataset(newyields_dict)

for year in years:
    # extract out a year's worth of obs
    # TODO adapt to handle non-360day calendars
    # can use days_in_month code from PET project
    if startyear==endyear:
        gyeardates = xr.cftime_range(str(year) + '-' + startmonthstr + '-' + startdaystr, 
                                     str(year) + '-' + endmonthstr + '-' + enddaystr, 
                                     calendar='360_day', freq='D', name='date').values
    elif len(years)==1 and endmonth < startmonth:
        gyeardates = xr.cftime_range(str(year-1) + '-' + startmonthstr + '-' + startdaystr, 
                                     str(year) + '-' + endmonthstr + '-' + enddaystr, 
                                     calendar='360_day', freq='D', name='date').values
    elif len(years)==1 and endmonth > startmonth:
        gyeardates = xr.cftime_range(str(year-1) + '-' + startmonthstr + '-' + startdaystr, 
                                     str(year) + '-' + startmonthstr + '-' + startdaystr, 
                                     calendar='360_day', freq='D', name='date').values
    elif len(years)==1 and endday < startday:
        gyeardates = xr.cftime_range(str(year-1) + '-' + startmonthstr + '-' + startdaystr, 
                                     str(year) + '-' + endmonthstr + '-' + enddaystr, 
                                     calendar='360_day', freq='D', name='date').values
    elif len(years)==1 and endday >= startday:
        gyeardates = xr.cftime_range(str(year-1) + '-' + startmonthstr + '-' + startdaystr, 
                                     str(year) + '-' + startmonthstr + '-' + startdaystr, 
                                     calendar='360_day', freq='D', name='date').values
    else:
        gyeardates = xr.cftime_range(str(year-1) + '-' + startmonthstr + '-' + startdaystr, 
                                     str(year) + '-' + startmonthstr + '-' + startdaystr, 
                                     calendar='360_day', freq='D', name='date').values

    print('Running for ' + str(gyeardates[0]) + ' to ' + str(gyeardates[-1]))
    try:
        obspd = obspdall.loc[gyeardates]
    except KeyError:
        print('Required times not present in downloaded MODIS data, re-run code with obs=1 to download')
        sys.exit()
    times = np.array([str(t.year)+'{:02.0f}'.format(t.month)+'{:02.0f}'.format(t.day) for t in list(gyeardates)])    

    if plotdir != 'None':
        if not os.path.exists(plotdir):
            os.makedirs(plotdir)

    prior_p_all, prior_p_ensmean, prior_p_ensstd, nonprior_p_all, \
    tmean_p_all, prec_p_all, solarrad_p_all, \
    Jarray_p_all, Cday_p_all, GSS_p_all, \
    HarvestJday, AWC_allp, CDD, TT, temp_cconc = \
    ensgen_point(ensmems, assimvar, times, moderrtype, moderr,
                 datasetname, precname, radname, crop, coords, 
                 dataloc, saveloc, AWCrast, elevfile, CO2file)

    # now we've done the ensgen, do the assimilation
    cfxrs_vari = []
    post_errs_all = []
    for ensmem in ensmems:
        print('Doing DA for ' + assimvar + ' using variational method for ensmem ' + ensmem)
        mergedall_vari, cfall_vari, post_errs = vari_method_point(obspd, prior_p_all.sel(ensmem=int(ensmem)), prior_p_ensstd, 
                                                       assimvar, coords, year, ensmem, obserrtype, obserr, moderrinfl,
                                                       tsvar, order, power, plot, plotdir, interval)

        cfxr_vari = cfall_vari.to_xarray()
        post_errs_xr = post_errs.to_xarray()
        del cfall_vari
        del post_errs
        cfxr_vari = cfxr_vari.expand_dims({'ensmem': [int(ensmem)]})
        post_errs_xr = post_errs_xr.expand_dims({'ensmem': [int(ensmem)]})
        cfxr_vari = xr.where(np.isinf(cfxr_vari), np.nan, cfxr_vari) # set infs as nans
        post_errs_xr = xr.where(np.isinf(post_errs_xr), np.nan, post_errs_xr) # set infs as nans
        cfxrs_vari.append(cfxr_vari)
        post_errs_all.append(post_errs_xr)

        oldyields, newyields = update_yield_points_point(prior_p_all.sel(ensmem=int(ensmem)), mergedall_vari, coords, assimvar, 
                                         nonprior_p_all.sel(ensmem=int(ensmem)), tmean_p_all.sel(ensmem=int(ensmem)), 
                                         prec_p_all.sel(ensmem=int(ensmem)), solarrad_p_all.sel(ensmem=int(ensmem)), 
                                         Jarray_p_all.sel(ensmem=int(ensmem)), Cday_p_all.sel(ensmem=int(ensmem)), 
                                         GSS_p_all.sel(ensmem=int(ensmem)), HarvestJday, AWC_allp, CDD, TT, temp_cconc, ensmem)

        del mergedall_vari
        counter=0
        for tob in coords:
            oldyields_vari[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmem), year=year)] = oldyields[counter]
            newyields_vari[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmem), year=year)] = newyields[counter]
            counter+=1

    cfallensmems_vari = xr.concat(cfxrs_vari, dim='ensmem')
    post_errs_allensmems = xr.concat(post_errs_all, dim='ensmem')
    cfyears_vari.append(cfallensmems_vari)
    post_errs_years.append(post_errs_allensmems)

if len(years) > 1:
    print('Calculating climatological conversion factor using years: ' + str(years) + '\n')
    cfallyears_vari = xr.concat(cfyears_vari, dim='date')
    del cfyears_vari
    cfensmean_vari = cfallyears_vari.mean(dim='ensmem')
    del cfallyears_vari
    dayroll = -1*(startmonth*30 + startday-1)
    cfclim_vari = cfensmean_vari.groupby('date.dayofyear').mean().roll(dayofyear=dayroll, roll_coords=True)

print('Saving to netcdf')
if not os.path.exists(yieldsaveloc):
    os.makedirs(yieldsaveloc)
now = dt.datetime.now().strftime('%Y%m%d%H%M%s')[:14]
oldyieldfile = os.path.join(yieldsaveloc, now + '_oldyield.nc')
newyieldfile = os.path.join(yieldsaveloc, now + '_newyield.nc')
cfclimfile   = os.path.join(yieldsaveloc, now + '_cfclimo.nc')
oldyields_vari.to_netcdf(oldyieldfile)
newyields_vari.to_netcdf(newyieldfile)
if len(years) > 1:
    cfclim_vari.to_netcdf(cfclimfile)

post_errs_years = xr.concat(post_errs_years, dim='date')

#if verify==1:
#    verifyyield(yieldfile, coords, oldyields_vari, newyields_vari, years, yieldsaveloc)
