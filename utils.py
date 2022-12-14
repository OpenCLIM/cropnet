import os
import sys
import glob
import pyproj
import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd
import subprocess as subp
import matplotlib.pyplot as plt
import scipy.optimize as spo

from cdo import *
from affine import Affine
from IPython import display
from rasterio import features
from dask.diagnostics import  ProgressBar


import rpy2
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages

from rpy2.robjects import numpy2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector
numpy2ri.activate()


# Set up R environment on import of this module
rgdal = importr('rgdal')
rgeos = importr('rgeos')
ncdf4 = importr('ncdf4')
abind = importr('abind')
raster = importr('raster')
#Evapotranspiration = importr('Evapotranspiration')
r = robjects.r
try:
    r['source']('Grass_potpredict_MJB.R') # define R functions
    r['source']('Lynch_potpredict_v2_MJB.R')
except rpy2.rinterface.embedded.RRuntimeError:
    print('Warning: R crop model files not available')


def load_driving_data(basedatasetname, times,
                      saveloc, dataloc, simx, crop, AWCrast='None', CO2file='None',
                      lonmin=0, lonmax=0, latmin=0, latmax=0, downloaddata=0,
                      mask=0, sfname='ne_10m_admin_0_countries.shp', countries=['China'],
                      cropmask=0, cmfile=None, cropthresh=0, irrmask=0, irrfile=None,
                      precipname='None', radname='None'):

    '''
    Read in the various data files and merge them into a single
    netcdf file on the same grid, apply masking and load 
    ancillary data files. 
    Inputs:
    basedatasetname: Which driving dataset we're using
    times: 1D np array that defines the time period that the data will be cut down to
    dataloc: Directory containing all the data files
    simx: Only used for ukcp18 data, to identify the ensemble member used
    crop: Which crop to load data for
    AWCrast: Location of the AWC dataset used for wheat crops. Defaults to 'None'.
    CO2file: Location of csv file containing yearly CO2 concentrations. Defaults to 'None'.
    lonmin, lonmax, latmin, latmax: Region to download if downloaddata==1. All default to 0.
    downloaddata: Switch to control whether the basedataset is downloaded
                  or not. Currently only works for era5. Defaults to 0 (off).
    mask: Switch to enable masking of data by country. Defaults to 0 (off).
    sfname: The location of the shapefile containing at least the countries 
            to mask to. Defaults to one in current directory. 
    countries: List of countries to mask for. Defaults to China. 
    cropmask: Switch controlling masking based on crop. Defaults to 0 (off).
    cmfile: Location of dataset defining location of crop for masking. Defaults to None.
    cropthresh: Threshold in cmfile to use for masking. Defaults to 0.
    irrmask: Switch controlling whether to account for irrigation or not.
             Defaults to 0 (off).
    irrfile: File containing proportion of each pixel that is irrigated.
             Should be on same grid as cmfile, ideally from same dataset.
             Defaults to None.
    precipname: The name of the precip dataset to use
    radname: The name of the solar radiation dataset to use   

    Outputs:
    alldata: Dictionary containing all the variables needed to run the crop 
             model as np arrays. Keys are from dnames (see getnames function). 
    '''
    
    fnames,vnames,dnames,xnames,ynames,tnames,interpvar = getnames(basedatasetname, precipname, radname, crop)
    
    if downloaddata == 1:
        download_data(os.path.dirname(dataloc), np.arange(startyear,endyear+1), vnames,
                      lonmin, lonmax, latmin, latmax)

    print('Producing subsetted driving data nc file(s):')
    data = process_driving_data(basedatasetname, fnames, vnames, dnames, xnames, ynames, tnames, interpvar,
                                saveloc, dataloc, simx, times, precipname, radname,
                                mask, sfname, countries, crop, cropmask, cmfile, cropthresh)


    interpind = dnames.index(interpvar)

    x = data['x'].values # these can still be lon/lat even though named x/y
    y = data['y'].values

    alldata = {}
    print('Processing data into format ready for model')
    for v in range(0,len(vnames)):
        print('Reading in files')
        if type(vnames[v]) == list:
            if radname=='ceres-syn':
                alldata[dnames[v]] = (data[vnames[v][0]] - data[vnames[v][1]]).values.squeeze().transpose(1,2,0)
            else:
                raise TypeError(vnames[v] + ' should not be a list')
            #print(alldata[dnames[v]])
        else:
            alldata[dnames[v]] = data[vnames[v]].values.squeeze().transpose(1,2,0)
            
    alldata['x'] = x
    alldata['y'] = y
    alldata['t'] = times
    
    if basedatasetname=='era5':
        if crop=='grass':
            alldata['wind'] = np.sqrt(alldata['uwind']**2 + alldata['vwind']**2)        

    if crop=='wheat':
        if irrmask == 1:
            print('Reading in irrmask')
            irrprop = xr.load_dataarray(irrfile)
            irrprop = irrprop.interp({'lon': data['x'], 'lat': data['y']}, method='nearest')
            alldata['irrprop'] = irrprop.values

    # process AWC
    if crop=='wheat':
        print('Reading in AWC')
        interpdata = data
        interpx = 'x' # this is always x
        interpy = 'y' # and this is always y even if lonlat data
        # they are just variable names
        
        if AWCrast[-2:] == 'nc':
            # netcdf data
            AWCdataset = xr.open_dataset(AWCrast)
            
            for key in AWCdataset.keys():
                if 'AWC' in key:
                    AWC = AWCdataset[key]
                    break
                elif 'awc' in key:
                    AWC = AWCdataset[key]
                    break
                elif 'Awc' in key:
                    AWC = AWCdataset[key]
                    break
            else:
                raise KeyError('AWC variable not found in ' + AWCrast)
                    

            # find out coordinate type (lon/lat or x/y)
            for key in AWCdataset.dims:
                if key == 'lon': 
                    ctype = 'll'
                    AWCx = 'lon'
                    AWCy = 'lat'
                    break
                elif key == 'longitude':
                    ctype = 'll'
                    AWCx = 'longitude'
                    AWCy = 'latitude'
                    break
                elif key == 'x': 
                    ctype = 'xy'
                    AWCx = 'x'
                    AWCy = 'y'
                    break
                elif key == 'projection_x_coordinate':
                    ctype = 'xy'
                    AWCx = 'projection_x_coordinate'
                    AWCy = 'projection_y_coordinate'
                    break
            else:
                print('WARNING: Could not determine coord type, assuming x/y')
                ctype = 'xy'
                AWCx = 'x'
                AWCy = 'y'

            # coarsen/aggregate if necessary
            # agg=True
            datares = abs(data['x'][1] - data['x'][0]).values
            AWCres  = abs(AWC[AWCx].values[1] - AWC[AWCx].values[0])
            resrat = datares/AWCres
            if resrat >= 2:
                agg = True
            else:
                agg = False

            if ctype == 'll':
                # lon/lat
                if agg:
                    # coarsen by ratio of resolutions?
                    resratagg = np.floor(resrat)
                    if AWCx == 'lon':
                        AWC = AWC.coarsen(lon=resratagg, lat=resratagg, boundary='pad').mean()
                    elif AWCx == 'longitude':
                        AWC = AWC.coarsen(longitude=resratagg, latitude=resratagg, boundary='pad').mean()
                    
            elif ctype == 'xy':
                # x/y
                if agg:
                    # coarsen by ratio of resolutions?
                    resratagg = np.floor(resrat)
                    if AWCx == 'x':
                        AWC = AWC.coarsen(x=resratagg, y=resratagg, boundary='pad').mean()
                    elif AWCx == 'projection_x_coordinate':
                        AWC = AWC.coarsen(projection_x_coordinate=resratagg, 
                                          projection_y_coordinate=resratagg, boundary='pad').mean()

            # Make the AWC three dimensional (add a time dimension)
            AWClist=[]
            for ts in range(0,len(times)):
                AWCtemp = AWC.expand_dims('t')
                AWCtemp['t'] = [ts]
                AWClist.append(AWCtemp)
            AWC = xr.concat(AWClist, dim='t')

        else:
            # raster data
            # raster data is unlikely to be lon/lat so we can ignore this case
            
            # aggregation logic built into R func
            datares = abs(alldata['x'].values[1] - alldata['x'].values[0])
            # use loadAWC R function
            loadAWC = r['load_AWC']
            AWC = loadAWC(x, y, AWCrast, datares, np.array(alldata[dnames[0]].transpose(1,0,2).shape))
            AWC = np.array(AWC).transpose(1,0,2)
            AWC = xr.DataArray(AWC, coords=[y, x, times], dims=['y', 'x', 't'])
        

        ## ADAPTING THIS BLOCK ABOVE
        #if basedatasetname=='era5':
        #    AWCx = 'lon'
        #    AWCy = 'lat'
        #else:
        #    AWCx = 'x'
        #    AWCy = 'y'
        #if basedatasetname=='ukcp18bc' or basedatasetname=='chess_and_haduk':
        #    loadAWCnoagg = r['load_AWC_no_agg']
        #    AWC = loadAWCnoagg(x, y, AWCrast, np.array(alldata['tmean'].transpose(1,0,2).shape))
        #    AWC = np.array(AWC).transpose(1,0,2)
        #    AWC = xr.DataArray(AWC, coords=[y, x, times], dims=['y', 'x', 't'])
        #elif basedatasetname == 'era5':
        #    AWC = xr.load_dataarray(AWCrast)
        #    AWC = AWC.coarsen(lon=30, lat=30, boundary='pad').mean()
        #    AWClist = []
        #    for ts in range(0,len(times)):
        #        AWCtemp = AWC.expand_dims('t')
        #        AWCtemp['t'] = [ts]
        #        AWClist.append(AWCtemp)
        #    AWC = xr.concat(AWClist, dim='t')
        #else:
        #    loadAWC = r['load_AWC']
        #    AWC = loadAWC(x, y, AWCrast, np.array(alldata['tmean'].transpose(1,0,2).shape))
        #    AWC = np.array(AWC).transpose(1,0,2)
        #    AWC = xr.DataArray(AWC, coords=[y, x, times], dims=['y', 'x', 't'])

        print('Interpolating AWC onto common grid')
        AWC = AWC.interp({AWCx: interpdata[interpx], AWCy: interpdata[interpy]}, kwargs={"fill_value": None})
        if not AWCx==interpx:
            AWC=AWC.drop([AWCx])
        if not AWCy==interpy:
            AWC=AWC.drop([AWCy])

        if mask==1:
            print('Masking AWC to ' + str(countries))
            AWC = country_subset_shapefile(data=AWC, sfname=sfname, IDname='ADMIN', IDs=countries,
                                           xname=interpx, yname=interpy, drop=False)

    if crop=='wheat':
        AWC = AWC.transpose(interpy, interpx, 't')
        alldata['AWC'] = AWC.values#.transpose(1,0,2)
        #print(alldata['AWC'].shape)
        
    if CO2file!='None':
        pCO2 = pd.read_csv(CO2file)
        pCO2.set_index('YEAR', inplace=True)
        cconc = pCO2.loc[int(times[0][:4])].values[0]
        alldata['cconc'] = cconc
        
    # mask data to where we have AWC values
    #print(data)
    #print(AWC)
    #if crop=='wheat':
    #    data = data.where(AWC > 0)


    print('Done')
    print('Passing data to model')     
        
    return alldata


# not used
def load_driving_data_point(basedatasetname, times, coords, dataloc, saveloc,
                            simx, crop, elevfile='None', AWCrast='None', CO2file='None',
                            lonmin=0, lonmax=0, latmin=0, latmax=0, downloaddata=0,
                            mask=0, sfname='ne_10m_admin_0_countries.shp', countries=['China'],
                            cropmask=0, cmfile='None', cropthresh=0, irrmask=0, irrfile='None',
                            precipname='None', radname='None'):

    '''
    Read in the various data files and merge them into a single
    netcdf file on the same grid, apply masking and load 
    ancillary data files. 
    Inputs:
    basedatasetname: Which driving dataset we're using
    times: 1D np array of yyyymmdd strings. This defines the time period that the data will be cut down to
    coords: List of (x,y) coordinate pairs to extract data for
    dataloc: Directory containing all the data files
    saveloc: Where to output the single netcdf file containing all the variables
    simx: Only used for ukcp18 data, to identify the ensemble member used
    crop: Which crop to load data for
    elevfile: Location of raster dataset containing elevation data for the grass crop
    AWCrast: Location of the AWC dataset used for wheat crops. Defaults to 'None'.
    CO2file: Location of csv file containing yearly CO2 concentrations. Defaults to 'None'.
    lonmin, lonmax, latmin, latmax: Region to download if downloaddata==1. All default to 0.
    downloaddata: Switch to control whether the basedataset is downloaded
                  or not. Currently only works for era5. Defaults to 0 (off).
    mask: Switch to enable masking of data by country. Defaults to 0 (off).
    sfname: The location of the shapefile containing at least the countries 
            to mask to. Defaults to one in current directory. 
    countries: List of countries to mask for. Defaults to China. 
    cropmask: Switch controlling masking based on crop. Defaults to 0 (off).
    cmfile: Location of dataset defining location of crop for masking. Defaults to None.
    cropthresh: Threshold in cmfile to use for masking. Defaults to 0.
    irrmask: Switch controlling whether to account for irrigation or not.
             Defaults to 0 (off).
    irrfile: File containing proportion of each pixel that is irrigated.
             Should be on same grid as cmfile, ideally from same dataset.
             Defaults to None.
    precipname: The name of the precip dataset to use
    radname: The name of the solar radiation dataset to use   

    Outputs:
    alldata: Dictionary containing all the variables needed to run the crop 
             model as xr datasets for each variable. The variables within each
             xr dataset is the value(s) at each coordinate. These variable names
             are the coordinates (x,y). 
    '''
        
    fnames,vnames,dnames,xnames,ynames,tnames,interpvar = getnames(basedatasetname, precipname, radname, crop)
           
    print('Producing subsetted driving data nc file(s):')
    datafile = process_driving_data(basedatasetname, fnames, vnames, dnames, xnames, ynames, tnames, interpvar,
                                    dataloc, saveloc, simx, times, precipname, radname,
                                    mask, sfname, countries, crop, cropmask, cmfile, cropthresh)

    # Read in the data from the file into a dictionary with keys specified by dnames
    # Each variable in the dictionary will be an xarray dataarray

    alldata = {}
    data = xr.load_dataset(datafile)
    for v in range(0, len(fnames)):
        print('Standardising names ' + dnames[v])
        if type(vnames[v]) == list:
            if radname=='ceres-syn':
                alldata[dnames[v]] = data[vnames[v][0]] - data[vnames[v][1]]
            else:
                raise TypeError(vnames[v] + ' should not be a list')
        else:
            alldata[dnames[v]] = data[vnames[v]]

    print('Getting x,y coords')
    x = alldata[dnames[0]]['x'].values
    y = alldata[dnames[0]]['y'].values

    print('Formatting data as point variables')
    alldata_p = format_point_data(coords, dnames, alldata, times)

    #for v in range(0, len(fnames)):
    #    savefile = os.path.join(saveloc, dnames[v] + simx + '_' + times[0] + \
    #                            '-' + times[-1] + '.nc')
    #    # save to netcdf
    #    print('Processing and saving ' + dnames[v] + ' to disk')
    #    alldata_p[dnames[v]].to_netcdf(savefile)

    #alldata_p2 = {}
    #for v in range(0, len(filenames)):
    #    savefile = os.path.join(saveloc, dnames[v] + simx + '_' + times[0] + \
    #                            '-' + times[-1] + '.nc')
    #    
    #    print('Reading in processed ' + dnames[v] + ' data')
    #    alldata_p2[dnames[v]] = xr.load_dataset(savefile)

    if crop=='grass':
        if basedatasetname == 'ukcp18' or basedatasetname=='ukcp18bc':
            loadelev = r['load_elev']
            elevs = loadelev(x, y, elevfile)
            elevs = np.array(elevs).squeeze()
            alldata['elevs'] = xr.DataArray(elevs, coords=[('x', x), ('y', y)])
            alldata_p = format_point_data(coords, ['elevs'], alldata, datadict_p=alldata_p)
        elif basedatasetname == 'era5':
            alldata_p['wind'] = np.sqrt(alldata_p['uwind']**2 + alldata_p['vwind']**2)
        
    elif crop=='wheat':
        if basedatasetname=='ukcp18bc' or basedatasetname=='chess_and_haduk':
            # The load_AWC function will pull out the x,y coords supplied, using aggregation and averaging
            # The AWC data doesn't have to be on the same grid as the other data
            loadAWC = r['load_AWC_no_agg']
            AWC = loadAWC(x, y, AWCrast, np.array(alldata['tmean'].values.squeeze().transpose(2,1,0).shape))
            AWC = np.array(AWC).squeeze()
            alldata['AWC'] = xr.DataArray(AWC, coords=[('x', x), ('y', y), ('time', np.arange(0, len(times)))])
            alldata_p = format_point_data(coords, ['AWC'], alldata, times, alldata_p)
        if basedatasetname=='era5':
            AWC = xr.load_dataarray(AWCrast)
            AWC = AWC.coarsen(lon=30, lat=30, boundary='pad').mean()
            AWClist = []
            for ts in range(0,len(times)):
                AWCtemp = AWC.expand_dims('t')
                AWCtemp['t'] = [ts]
                AWClist.append(AWCtemp)
            AWC = xr.concat(AWClist, dim='t')
            AWC = AWC.interp({'lon': data['x'], 'lat': data['y']})
            alldata['AWC'] = AWC
            alldata_p = format_point_data(coords, ['AWC'], alldata, times, alldata_p)
        else:
            loadAWC = r['load_AWC']
            AWC = loadAWC(x, y, AWCrast, np.array(alldata['tmean'].values.squeeze().transpose(2,1,0).shape))
            AWC = np.array(AWC).squeeze()
            alldata['AWC'] = xr.DataArray(AWC, coords=[('x', x), ('y', y), ('time', np.arange(0, len(times)))])
            alldata_p = format_point_data(coords, ['AWC'], alldata, times, alldata_p)
            
    if CO2file:
        pCO2 = pd.read_csv(CO2file)
        pCO2.set_index('YEAR', inplace=True)
        cconc = pCO2.loc[int(times[0][:4])].values[0]
        alldata_p['cconc'] = cconc

    alldata_p['t'] = times
    del alldata
    return alldata_p


def process_driving_data(basedatasetname, filenames, vnames, dnames, xnames, ynames, tnames, interpvar, saveloc, dataloc, simx,
                         times, precipname, radname, mask, sfname, countries, crop, cropmask, cmfile, cropthresh):
    '''
    Read in data from various datasets and merge them into one netcdf file. 
    Inputs:
    basedatasetname: Which driving dataset we're using
    filenames: Part of the filenames of each file that identifies the variable within it.
               Must correspond with vnames and dnames.
    vnames: The names of the variables in the netcdf files to be read in.
            Must correspond with filenames and dnames.
    dnames: The names of the variables to be used throughout the rest of the code
            Must correspond with filenames and vnames.
    xnames: The x dimnames of the variables, must correspond with the above
    ynames: The y dimnames of the variables, must correspond with the above
    tnames: The time dimnames of the variables, must correspond with the above
    dataloc: Directory containing all the data files
    simx: Only used for ukcp18 data, to identify the ensemble member used
    times: 1D np array of yyyymmdd strings defining the time period that the data will be cut down to
    precipname: The name of the precip dataset to use
    radname: The name of the solar radiation dataset to use
    Outputs: 
    savefile: The location of the created netcdf file containing all the variables required to run the cropmodel
    '''

    # get list of extracted nc files
    nlist = glob.glob(dataloc, recursive=True)
    #print(nlist)

    # pull out only files corresponding to this ensemble member
    if basedatasetname == 'ukcp18':
        searchstr = simx + '_day_'
        snlist = [name for name in nlist if searchstr in name]
    else:
        snlist = nlist

    # extract only those files corresponding to the variables in filenames
    snlist2 = []
    for fname in filenames:
        for filename in snlist:
            if fname in os.path.basename(filename):
                snlist2.append(filename)

    snlist = snlist2
    #print(snlist)
        
    # extract out coordinates to interpolate all variables on to 
    intind = dnames.index(interpvar)
    vlist = [name for name in snlist if filenames[intind] in os.path.basename(name)]
    dataset = xr.open_mfdataset(vlist, parallel=True, combine='by_coords')
    subsetds = dataset.sel(time=slice(times[0][:4] + '-' + times[0][4:6] + '-' + times[0][6:8], \
                                      times[-1][:4] + '-' + times[-1][4:6] + '-' + times[-1][6:8]))
    interpxs = subsetds[xnames[intind]]
    interpys = subsetds[ynames[intind]]
    interpts = subsetds[tnames[intind]]
    
    # for each variable in varnames...
    data = []
    for v in range(0, len(filenames)):
        #print('Processing ' + filenames[v])
        # pull out only the files corresponding to this particular variable
        vlist = [name for name in snlist if filenames[v] in os.path.basename(name)]
        
        # open all the files selected into one xarray dataset using dask
        #print('Opening ' + str(vlist))
        dataset = xr.open_mfdataset(vlist, parallel=True, combine='by_coords')
        #print(dataset)
        # subset to required time period
        subsetds = dataset.sel(time=slice(times[0][:4] + '-' + times[0][4:6] + '-' + times[0][6:8], \
                                          times[-1][:4] + '-' + times[-1][4:6] + '-' + times[-1][6:8]))
        
        # interpolate onto common grid (if already on the common grid, nothing will happen)
        subsetinterp = subsetds.interp({xnames[v]: interpxs, ynames[v]: interpys},
                                       kwargs={'fill_value': None})
        
        # ensure time coordinate metadata and labels are the same between datasets
        subsetinterp[tnames[v]] = interpts
        #print(subsetds)
        
        if mask==1:
            print('Masking data to ' + str(countries))
            subsetinterp = country_subset_shapefile(data=subsetinterp, sfname=sfname, IDname='ADMIN', IDs=countries,
                                                    xname=xnames[intind], yname=ynames[intind], drop=False)
            
        if crop=='wheat':
            # mask out non-crop areas
            # currently only used for wheat crop
            # and for lon/lat driving data - will need to adapt this to work with UKCP18 data
            # or create a new cropmap/irrmap in the UKCP18 x/y coords
            if cropmask==1:
                cropmap = xr.load_dataarray(cmfile)
                cropmap = cropmap.interp({'lon': subsetinterp[xnames[intind]], 'lat': subsetinterp[ynames[intind]]})
                subsetinterp = subsetinterp.where(cropmap > cropthresh)
                
        subsetinterp = subsetinterp.rename({xnames[intind]: 'x', ynames[intind]: 'y'})
        data.append(subsetinterp)
        
    print('Merging data')
    mergedata = xr.merge(data)
    
    print('Processing, extracting and loading data to memory')
    with ProgressBar():
        mergedata.load()

    return mergedata


def getnames(basedatasetname, precipname, radname, crop):
    '''
    Generate the file, variable and dictionary names associated
    with each dataset.
    Users can edit this function to set up the code to use new meteorological driving datasets
    - fnames are the variable names that appear in the filenames
    - vnames are the variable names in the netcdf files
    - dnames are the variable names that the code uses to refer to different variables.
      For the following variables, the corresponding dname should be set to:
      precipitation: prec
      Max temperature: tmax
      Min temperature: tmin
      Relative humidity: rh
      Solar radiation: solarrad
      x component of wind speed: uwind
      y component of wind speed: vwind
      Wind speed: wind
      Mean temperature: tmean
      Dew point temperature: tdp
      Surface air pressure: sfcP
    - xnames are the variable names of the x coordinate of each variable
    - ynames are the variable names of the y coordinate of each variable
    - tnames are the variable names of the t coordinate of each variable
    For xnames, ynames and tnames, if all variables have the same coordinate variables, the 
    coordinate variable name can just be listed once.
    For the other lists, the order of the variables must be the same, i.e. the first item 
    in the lists corresponds to the same variable, the second items to another etc. 
    '''
    
    if basedatasetname == 'ukcp18':
        if crop == 'grass':
            fnames = ["pr", "tasmax", "tasmin", "hurs", "rss", "sfcWind", "tas"]
            vnames = fnames
            dnames = ['prec', 'tmax', 'tmin', 'rh', 'solarrad', 'wind', 'tmean']
            xnames = ['projection_x_coordinate']
            ynames = ['projection_y_coordinate']
            tnames = ['time']
        elif crop == 'wheat':
            fnames = ["pr", "tasmax", "tasmin", "rss", "tas"]
            vnames = fnames
            dnames = ['prec', 'tmax', 'tmin', 'solarrad', 'tmean']
            xnames = ['projection_x_coordinate']
            ynames = ['projection_y_coordinate']
            tnames = ['time']
    elif basedatasetname == 'era5':
        if crop == 'grass':
            fnames = ['total_precipitation', 'maximum_2m_temperature_since_previous_post_processing',
                      'minimum_2m_temperature_since_previous_post_processing', '2m_dewpoint_temperature',
                      'surface_net_solar_radiation', '10m_u_component_of_wind', '10m_v_component_of_wind',
                      '2m_temperature','surface_pressure']
            vnames = ['tp', 'mx2t', 'mn2t', 'd2m', 'ssr', 'u10', 'v10', 't2m' ,'sp']
            dnames = ['prec', 'tmax', 'tmin', 'tdp', 'solarrad', 'uwind', 'vwind', 'tmean', 'sfcP']
            xnames = ['longitude']
            ynames = ['latitude']
            tnames = ['time']
        elif crop == 'wheat':
            fnames = ['total_precipitation', 'maximum_2m_temperature_since_previous_post_processing',
                      'minimum_2m_temperature_since_previous_post_processing', 'surface_net_solar_radiation',
                      '2m_temperature']
            vnames = ['tp' ,'mx2t', 'mn2t', 'ssr', 't2m']
            dnames = ['prec', 'tmax', 'tmin', 'solarrad', 'tmean']
            xnames = ['longitude']
            ynames = ['latitude']
            tnames = ['time']
    elif basedatasetname == 'ukcp18bc':
        if crop == 'grass':
            fnames = ["pr", "tasmax", "tasmin", "hurs", "rss", "sfcWind", "tas"]
            vnames = fnames
            dnames = ['prec', 'tmax', 'tmin', 'rh', 'solarrad', 'wind', 'tmean']
            xnames = ['x']
            ynames = ['y']
            tnames = ['time']
        elif crop == 'wheat':
            fnames = ["pr", "tasmax", "tasmin", "rss", "tas"]
            vnames = ['pr', 'tasmax', 'tasmin', 'rss', 'tas']
            dnames = ['prec', 'tmax', 'tmin', 'solarrad', 'tmean']
            xnames = ['x']
            ynames = ['y']
            tnames = ['time']
    elif basedatasetname == 'chess_and_haduk':
        if crop == 'wheat':
            fnames = ['precip', 'tasmax', 'tasmin', 'netsolar', 'tasmean']
            vnames = fnames
            dnames = ['prec', 'tmax', 'tmin', 'solarrad', 'tmean']
            xnames = ['x']
            ynames = ['y']
            tnames = ['time']

    if len(xnames)==1:
        xnames2 = [xnames[0] for n in fnames]
    else:
        xnames2 = xnames
    if len(ynames)==1:
        ynames2 = [ynames[0] for n in fnames]
    else:
        ynames2 = ynames
    if len(tnames)==1:
        tnames2 = [tnames[0] for n in fnames]
    else:
        tnames2 = tnames

    precindex = dnames.index('prec')
    radindex  = dnames.index('solarrad')

    # determine an index that isn't rad or precip,
    # to be used for interpolation later
    for i in range(0, len(dnames)):
        if i == precindex or i == radindex:
            pass
        else:
            interpvar = dnames[i]
            break
    if precipname == 'aphrodite':
        fnames[precindex] = 'APHRO_MA_025deg_V1101'
        vnames[precindex] = 'precip'
        xnames2[precindex] = 'lon'
        ynames2[precindex] = 'lat'
        tnames2[precindex] = 'time'
    if radname == 'ceres-syn':
        fnames[radindex] = 'CERES_SYN1deg-Day'
        vnames[radindex] = ['adj_atmos_sw_down_all_surface_daily', 'adj_atmos_sw_up_all_surface_daily']
        xnames2[radindex] = 'lon'
        ynames2[radindex] = 'lat'
        tnames2[radindex] = 'time'
    if radname == 'ceres-flash':
        fnames[radindex] = 'CERES_FLASH_TISA'
        vnames[radindex] = 'sfc_net_sw_all_daily'
        xnames2[radindex] = 'lon'
        ynames2[radindex] = 'lat'
        tnames2[radindex] = 'time'

        
    return fnames,vnames,dnames,xnames2,ynames2,tnames2,interpvar


def format_point_data(coords, varnames, data, t=[None], datadict_p=None):
    '''
    Take the gridded data and format it into an xarray dataset with
    one variable per point we want to extract.

    Inputs:
    coords: List of (x,y) or (lon,lat) pairs, the coords we want to extract
            data for
    varnames: The names of the variables we want to extract data for
              These must correspond to the names of the variables in "data"
    data: Dictionary with names (keys) the same as varnames containing the 
          xarrays with the gridded data variables
    t: Time coordinates as a list
    datadict_p: If provided, the new variables in varnames will be added
                to this preexisting dictionary, outputted from a previous
                call of this function

    Outputs:
    datadict_p: A dictionary with names (keys) the same as varnames, containing
                the point formatted data as above
    '''
    
    datadict = {}
    if not datadict_p:
        datadict_p = {}

    # create dictionary of empty xarrays to store point data
    print('Creating dictionary to of empty xarrays to store point data')
    for var in varnames:
        datadict[var] = {}
        for tob in coords:
            if t[0] != None: 
                datadict[var][str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((1, len(t))).squeeze(), coords=[t], dims=['time'])
            else: 
                datadict[var][str(tob[0])+','+str(tob[1])] = xr.DataArray(np.array(0.))

        datadict_p[var] = xr.Dataset(datadict[var])
        
        # fill empty xarray with data
        for tob in coords:
            #print('Processing ' + str(tob[0]) + ',' + str(tob[1]))
            datadict_p[var][str(tob[0])+','+str(tob[1])].values = data[var].sel(x=tob[0], y=tob[1], method='nearest').values.squeeze()
    return datadict_p


## NOT USED ##
def ensgen(ensmems, year, Rdatpath, CO2file):
    '''
    Generates a set of ensemble members of GAI. Returns this as an xarray dataset
    along with the ensemble mean and stdev. 

    Inputs: ensmems - List of 0-padded 2-digit strings corresponding to the ensemble 
                      member numbers
            year - Growing year to work with (int/float). Note this corresponds
                   to the year that harvest is made in.
            Rdatpath - Location of .Rdata driving data files (string)
            CO2file - Location of file containing CO2 conc. data (string)
    
    Outputs: GAI_all_merged - xarray dataset of GAI, dims (ensmem,x,y,t)
             GAI_ensmean - xarray dataset containing the ensmean of GAI (x,y,t)
             GAI_ensstd  - xarray dataset containing the ensstd of GAI (x,y,t)
    '''

    # Load R funcs
    loaddata = r['loaddata']
    GAIfunc = r['GAI']
    
    GAI_all=[]
    tmean_all=[]
    prec_all=[]
    solarrad_all=[]
    Jarray_all=[]
    Cday_all=[]
    GSS_all=[]
    for ensmem in ensmems:
        print('Loading data for ensmem ' + ensmem + ' for gyear ' + str(year))

        # read in driving data for this growing year from the .Rdata files
        datalist = loaddata(ensmem, year, Rdatpath, CO2file)

        temp_tas = datalist.rx2('temp_tas')
        temp_rls = datalist.rx2('temp_rls')
        temp_rss = datalist.rx2('temp_rss')
        temp_pr  = datalist.rx2('temp_pr')
        temp_tasmax = datalist.rx2('temp_tasmax')
        temp_tasmin = datalist.rx2('temp_tasmin')
        temp_cconc  = datalist.rx2('temp_cconc')
        if ensmem==ensmems[0]:
            AWC = datalist.rx2('AWC')
            AWC = np.array(AWC)
            
        # run crop growth model to calculate GAI
        datalist2 = GAIfunc(temp_tas, temp_tasmax, temp_tasmin, temp_pr, temp_rss)
        del temp_tas
        del temp_rls
        del temp_rss
        del temp_pr
        del temp_tasmax
        del temp_tasmin
        
        GAI         = np.array(datalist2.rx2('GAI')) # convert to python numpy
        tmean       = np.array(datalist2.rx2('tmean')) 
        prec        = np.array(datalist2.rx2('prec'))
        solarrad    = np.array(datalist2.rx2('solarrad'))
        Jarray      = np.array(datalist2.rx2('Jarray'))
        Cday        = np.array(datalist2.rx2('Cday'))
        GSS_r       = datalist2.rx2('GSS')
        GSS         = np.array(datalist2.rx2('GSS')) # As it's a string array it comes out as a flattened array
        GSS         = GSS.reshape((GSS_r.dim[2], GSS_r.dim[1], GSS_r.dim[0]))
        GSS         = GSS.transpose(2,1,0) # therefore we need to manually reshape it
        HarvestJday = datalist2.rx2('HarvestJday')

        # get the data coordinates
        x = np.array(datalist2.rx2('x'))
        y = np.array(datalist2.rx2('y'))
        t = np.array(datalist2.rx2('t'))
        attrs = {'units': 'days since ' + t[0], 'calendar': '360_day'}
        # Convert GAI and driving vars to xarray to make the space and time handling/wrangling easier
        # Could definitely put this in a separate function...
        GAI      = GAI.transpose(1,0,2)
        tmean    = tmean.transpose(1,0,2)
        prec     = prec.transpose(1,0,2)
        solarrad = solarrad.transpose(1,0,2)
        Jarray   = Jarray.transpose(1,0,2)
        Cday     = Cday.transpose(1,0,2)
        GSS      = GSS.transpose(1,0,2)

        if ensmem==ensmems[0]:
            AWC      = AWC.transpose(1,0,2)
            AWC      = xr.DataArray(AWC, coords=[('y', y), ('x', x), ('time', np.arange(0, len(t)))])
            AWC.name = 'AWC'
            AWC      = AWC.to_dataset()
            AWC.time.attrs = attrs
            AWC      = xr.decode_cf(AWC, decode_coords=False)
            AWC      = AWC.assign_coords(x=x)
            AWC      = AWC.assign_coords(y=y)        
        
        GAI      = GAI.reshape(1, GAI.shape[0], GAI.shape[1], GAI.shape[2])
        tmean    = tmean.reshape(1, tmean.shape[0], tmean.shape[1], tmean.shape[2])
        prec     = prec.reshape(1, prec.shape[0], prec.shape[1], prec.shape[2])
        solarrad = solarrad.reshape(1, solarrad.shape[0], solarrad.shape[1], solarrad.shape[2])
        Jarray   = Jarray.reshape(1, Jarray.shape[0], Jarray.shape[1], Jarray.shape[2])
        Cday     = Cday.reshape(1, Cday.shape[0], Cday.shape[1], Cday.shape[2])
        GSS      = GSS.reshape(1, GSS.shape[0], GSS.shape[1], GSS.shape[2])
        
        
        GAI      = xr.DataArray(GAI, coords=[('ensmem', [int(ensmem)]), ('y', y), ('x', x), ('time', np.arange(0, len(t)))])
        tmean    = xr.DataArray(tmean, coords=[('ensmem', [int(ensmem)]), ('y', y), ('x', x), ('time', np.arange(0, len(t)))])
        prec     = xr.DataArray(prec, coords=[('ensmem', [int(ensmem)]), ('y', y), ('x', x), ('time', np.arange(0, len(t)))])
        solarrad = xr.DataArray(solarrad, coords=[('ensmem', [int(ensmem)]), ('y', y), ('x', x), ('time', np.arange(0, len(t)))])
        Jarray   = xr.DataArray(Jarray, coords=[('ensmem', [int(ensmem)]), ('y', y), ('x', x), ('time', np.arange(0, len(t)))])
        Cday     = xr.DataArray(Cday, coords=[('ensmem', [int(ensmem)]), ('y', y), ('x', x), ('time', np.arange(0, len(t)))])
        GSS      = xr.DataArray(GSS, coords=[('ensmem', [int(ensmem)]), ('y', y), ('x', x), ('time', np.arange(0, len(t)))])
        GAI.name      = 'GAI'
        tmean.name    = 'tmean'
        prec.name     = 'prec'
        solarrad.name = 'solarrad'
        Jarray.name   = 'Jarray'
        Cday.name     = 'Cday'
        GSS.name      = 'GSS'
        GAI      = GAI.to_dataset()
        tmean    = tmean.to_dataset()
        prec     = prec.to_dataset()
        solarrad = solarrad.to_dataset()
        Jarray   = Jarray.to_dataset()
        Cday     = Cday.to_dataset()
        GSS      = GSS.to_dataset()
        GAI.time.attrs      = attrs
        tmean.time.attrs    = attrs
        prec.time.attrs     = attrs
        solarrad.time.attrs = attrs
        Jarray.time.attrs   = attrs
        Cday.time.attrs     = attrs
        GAI.time.attrs      = attrs
        GAI      = xr.decode_cf(GAI, decode_coords=False)
        tmean    = xr.decode_cf(tmean, decode_coords=False)
        prec     = xr.decode_cf(prec, decode_coords=False)
        solarrad = xr.decode_cf(solarrad, decode_coords=False)
        Jarray   = xr.decode_cf(Jarray, decode_coords=False)
        Cday     = xr.decode_cf(Cday, decode_coords=False)
        GSS      = xr.decode_cf(GSS, decode_coords=False)
        GAI      = GAI.assign_coords(x=x)
        tmean    = tmean.assign_coords(x=x)
        prec     = prec.assign_coords(x=x)
        solarrad = solarrad.assign_coords(x=x)
        Jarray   = Jarray.assign_coords(x=x)
        Cday     = Cday.assign_coords(x=x)
        GSS      = GSS.assign_coords(x=x)
        GAI      = GAI.assign_coords(y=y) # the cf decoding sometimes mangles the coords, so reset them
        tmean    = tmean.assign_coords(y=y)
        prec     = prec.assign_coords(y=y)
        solarrad = solarrad.assign_coords(y=y)
        Jarray   = Jarray.assign_coords(y=y)
        Cday     = Cday.assign_coords(y=y)
        GSS      = GSS.assign_coords(y=y)
        GAI_all.append(GAI)
        tmean_all.append(tmean)
        prec_all.append(prec)
        solarrad_all.append(solarrad)
        Jarray_all.append(Jarray)
        Cday_all.append(Cday)
        GSS_all.append(GSS)
        

    GAI_all_merged = xr.concat(GAI_all, dim='ensmem')
    tmean_all_merged = xr.concat(tmean_all, dim='ensmem')
    prec_all_merged = xr.concat(prec_all, dim='ensmem')
    solarrad_all_merged = xr.concat(solarrad_all, dim='ensmem')
    Jarray_all_merged = xr.concat(Jarray_all, dim='ensmem')
    Cday_all_merged = xr.concat(Cday_all, dim='ensmem')
    GSS_all_merged = xr.concat(GSS_all, dim='ensmem')
    del GAI_all
    del tmean_all
    del prec_all
    del solarrad_all
    del Jarray_all
    del Cday_all
    del GSS_all

    #print(type(tmean_all_merged))
    #print(tmean_all_merged)
    # Calculate ensemble mean and stdev at each timestep
    GAI_ensmean = GAI_all_merged.copy()
    GAI_ensmean['GAI'] = GAI_ensmean['GAI'].mean(axis=0)

    GAI_ensstd = GAI_all_merged.copy()
    GAI_ensstd['GAI'] = GAI_ensstd['GAI'].std(axis=0)

    return GAI_all_merged, GAI_ensmean, GAI_ensstd, tmean_all_merged, prec_all_merged, solarrad_all_merged, Jarray_all_merged, Cday_all_merged, GSS_all_merged, HarvestJday, AWC, temp_cconc, x, y, t


def ensgen_point(ensmems, assimvar, times, moderrtype, moderr, moderrswitch, year, datasetname, precname, radname, caltype, crop, coords, dataloc, saveloc, AWCrast, elevfile, CO2file):
    '''
    Generates a set of ensemble members of GAI at specified coordinates,
    or the closest model grid point to them.
    Returns this as an xarray dataset of dimensions time and ensmem, with
    each variable in the dataset corresponding to a different coordinate.
    Also returned are the ensemble mean and stdev, in a similar format.

    Inputs: ensmems - List of 0-padded 2-digit strings corresponding to the ensemble 
                      member numbers
            assimvar - The variable to perform the assimilation on. 'LAI' or 'fPAR'
            moderrtype - Defines the type of modstd used if only one ensemble member
                         if 1, modstd = mod*moderr for each tstep (relative)
                         if 0, modstd = moderr for each tstep (actual)
            moderrswitch - If 0, generate stdev from ensemble driving data
                           If 1, read this in from file
                           If 2, use a single number (fixed-in-time)
            moderr - The mod stdev to assume (either relative or actual) if moderrswitch==2
            year - The growing year (year of yield harvest)
            CO2file - Location of file containing CO2 conc. data (string)
            coords - List of 2-element lists containing the x,y coordinates
                        of interest on the OSGB eastings/northings grid. 
            times - 1D np array of yyyymmdd strings defining the time period to be run
    
    Outputs: GAI_p_all - xarray dataset of GAI, dims (ensmem,t), one variable
                         for each spatial point, named after the x,y coords
             As above for all the driving data variables 
             GAI_p_ensmean - xarray dataset containing the ensmean of GAI (t)
                             one variable per x,y coord 
             GAI_p_ensstd  - xarray dataset containing the ensstd of GAI (t)
                             one variable per x,y coord
             HarvestJday - Julian day of harvest, in R variable format, constant
                           for each ensmem and coordinate. Needed for later in the 
                           main code
             AWC_allp1 - Something to do with average water content of soil for each
                         model gridbox. Format as GAI_p_all but without the ensmem dim.
             temp_cconc - CO2 concentration timeseries, in R format. Needed for later.
    '''

    # Load R funcs
    GAIfunc_p = r['GAI_point']

    # Load the gridded data for each ensmem, then generate the ensemble,
    # mean and stdev for each obscoord
    print('Generating ensemble data')
    # gen xarrays to store data
    if moderrswitch==2:
        if assimvar == 'LAI':
            prior_p_all, fPAR_p_all, tmean_p_all, tmin_p_all, tmax_p_all, \
            prec_p_all, solarrad_p_all, Jarray_p_all, Cday_p_all, \
            GSS_p_all, AWC_allp1, prior_p_ensstd = \
            create_xrs(coords, ensmems, len(times), 1)
        elif assimvar == 'fPAR':
            GAI_p_all, prior_p_all, tmean_p_all, tmin_p_all, tmax_p_all, \
            prec_p_all, solarrad_p_all, Jarray_p_all, Cday_p_all, \
            GSS_p_all, AWC_allp1, prior_p_ensstd = \
            create_xrs(coords, ensmems, len(times), 1)
    else:
        if assimvar == 'LAI':
            prior_p_all, fPAR_p_all, tmean_p_all, tmin_p_all, tmax_p_all, \
            prec_p_all, solarrad_p_all, Jarray_p_all, Cday_p_all, GSS_p_all, AWC_allp1 = \
            create_xrs(coords, ensmems, len(times), 0)
        elif assimvar == 'fPAR':
            GAI_p_all, prior_p_all, tmean_p_all, tmin_p_all, tmax_p_all, \
            prec_p_all, solarrad_p_all, Jarray_p_all, Cday_p_all, GSS_p_all, AWC_allp1 = \
            create_xrs(coords, ensmems, len(times), 0)

    for ensmem in ensmems:
        print('Ensemble ' + str(ensmem))
        dd = load_driving_data_point(datasetname, times, coords, dataloc, saveloc, ensmem, crop, 
                                     elevfile=elevfile, AWCrast=AWCrast, CO2file=CO2file, precipname=precname, radname=radname)
        
        if ensmem == ensmems[0]:
            AWC_allp1 = dd['AWC']
            temp_cconc = dd['cconc']

        # Select out grid points nearest to each coord, calculate GAI at each point
        print('Calculating GAI at each coord')
        counter=1
        totalobs = len(coords)
        
        for tob in coords:
            print('Calculating GAI for ensmem ' + str(ensmem) + ' for obs ' + str(counter) + ' of ' + str(totalobs))
            counter+=1
            coordname = str(tob[0]) + ',' + str(tob[1])
            tmean_p    = dd['tmean'][coordname].values
            prec_p     = dd['prec'][coordname].values
            solarrad_p = dd['solarrad'][coordname].values
            tmax_p     = dd['tmax'][coordname].values
            tmin_p     = dd['tmin'][coordname].values

            # Convert to lon,lat
            proj = pyproj.Transformer.from_crs(27700, 4326, always_xy=True)
            try:
                lon,lat = proj.transform(tob[0],tob[1], errcheck=True)
            except pyproj.exceptions.ProjError:
                lon,lat = proj.transform(tob[0],tob[1], errcheck=True)

            datalist2 = GAIfunc_p(tmean_p, tmax_p, tmin_p, prec_p, solarrad_p, tob[0], tob[1], lat, times, datasetname, precname, radname)
            
            fPAR1        = np.array(datalist2.rx2('LI')) # convert to python numpy
            GAI1         = np.array(datalist2.rx2('GAI')) 
            tmean1       = np.array(datalist2.rx2('tmean')) 
            tmin1        = np.array(datalist2.rx2('tmin')) 
            tmax1        = np.array(datalist2.rx2('tmax')) 
            prec1        = np.array(datalist2.rx2('prec'))
            solarrad1    = np.array(datalist2.rx2('solarrad'))
            Jarray1      = np.array(datalist2.rx2('Jarray'))
            Cday1        = np.array(datalist2.rx2('Cday'))
            GSS1         = np.array(datalist2.rx2('GSS')) 
            CDD          = np.array(datalist2.rx2('CDD'))
            TT           = np.array(datalist2.rx2('TT'))
            HarvestJday = datalist2.rx2('HarvestJday')
            
            # get the data coordinates                                                                                       
            x = np.array(datalist2.rx2('x'))
            y = np.array(datalist2.rx2('y'))
            t = np.array(datalist2.rx2('t'))
            attrs = {'units': 'days since ' + t[0], 'calendar': caltype}

            if assimvar == 'LAI':
                prior_p_all[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmem))] = GAI1
                fPAR_p_all[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmem))] = fPAR1
            elif assimvar == 'fPAR':
                prior_p_all[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmem))] = fPAR1
                GAI_p_all[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmem))] = GAI1
            tmean_p_all[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmem))] = tmean1
            tmin_p_all[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmem))] = tmin1
            tmax_p_all[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmem))] = tmax1
            prec_p_all[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmem))] = prec1
            solarrad_p_all[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmem))] = solarrad1
            Jarray_p_all[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmem))] = Jarray1
            Cday_p_all[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmem))] = Cday1
            GSS_p_all[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmem))] = GSS1
            #if ensmem==ensmems[0]:
            #    AWC1      = AWC.sel(x=tob[0], y=tob[1], method='nearest').values
            #    AWC_allp1[str(tob[0])+','+str(tob[1])].values = AWC1

    if assimvar == 'LAI':
        fPAR_p_all.time.attrs = attrs
    elif assimvar == 'fPAR':
        GAI_p_all.time.attrs  = attrs
    prior_p_all.time.attrs    = attrs
    tmean_p_all.time.attrs    = attrs
    tmin_p_all.time.attrs     = attrs
    tmax_p_all.time.attrs     = attrs
    prec_p_all.time.attrs     = attrs
    solarrad_p_all.time.attrs = attrs
    Jarray_p_all.time.attrs   = attrs
    Cday_p_all.time.attrs     = attrs
    GSS_p_all.time.attrs      = attrs
    #AWC_allp1.time.attrs      = attrs
    if moderrswitch==2:
        prior_p_ensstd.time.attrs = attrs

    if assimvar == 'LAI':
        fPAR_p_all = xr.decode_cf(fPAR_p_all, decode_coords=False)
    elif assimvar == 'fPAR':
        GAI_p_all      = xr.decode_cf(GAI_p_all, decode_coords=False)
    prior_p_all    = xr.decode_cf(prior_p_all, decode_coords=False)
    tmean_p_all    = xr.decode_cf(tmean_p_all, decode_coords=False)
    tmin_p_all     = xr.decode_cf(tmin_p_all, decode_coords=False)
    tmax_p_all     = xr.decode_cf(tmax_p_all, decode_coords=False)
    prec_p_all     = xr.decode_cf(prec_p_all, decode_coords=False)
    solarrad_p_all = xr.decode_cf(solarrad_p_all, decode_coords=False)
    Jarray_p_all   = xr.decode_cf(Jarray_p_all, decode_coords=False)
    Cday_p_all     = xr.decode_cf(Cday_p_all, decode_coords=False)
    GSS_p_all      = xr.decode_cf(GSS_p_all, decode_coords=False)
    #AWC_allp1      = xr.decode_cf(AWC_allp1, decode_coords=False)
    if moderrswitch==2:
        prior_p_ensstd = xr.decode_cf(GAI_p_ensstd, decode_coords=False)

    prior_p_ensmean = prior_p_all.mean(axis=0)
    if moderrswitch==2:
        for tob in coords:
            if moderrtype == 0:
                modstd = moderr
            elif moderrtype == 1:
                modstd = moderr * prior_p_all[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmems[0]))]
            prior_p_ensstd[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmems[0]))] = modstd
        prior_p_ensstd = prior_p_ensstd.squeeze()
    elif moderrswitch==1:
        infile = str(year).join(moderr.split('????'))
        print('Reading in ' + infile)
        prior_p_ensstd = xr.load_dataset(infile)
    elif moderrswitch==0:
        prior_p_ensstd = prior_p_all.std(axis=0)
        prior_p_ensstd.to_netcdf(os.path.join(saveloc, datasetname + '_' + str(year) + 'ensstd.nc'))
    del tmean1
    del tmin1
    del tmax1
    del prec1
    del solarrad1
    del tmax_p
    del tmin_p
    del tmean_p
    del prec_p
    del solarrad_p

    if assimvar == 'LAI':
        return prior_p_all, prior_p_ensmean, prior_p_ensstd, fPAR_p_all, tmean_p_all, tmin_p_all, tmax_p_all, prec_p_all, solarrad_p_all, Jarray_p_all, Cday_p_all, GSS_p_all, HarvestJday, AWC_allp1, CDD, TT, temp_cconc
    elif assimvar == 'fPAR':
        return prior_p_all, prior_p_ensmean, prior_p_ensstd, GAI_p_all, tmean_p_all, tmin_p_all, tmax_p_all, prec_p_all, solarrad_p_all, Jarray_p_all, Cday_p_all, GSS_p_all, HarvestJday, AWC_allp1, CDD, TT, temp_cconc


def create_xrs(obscoords, ensmems, timlen, stdswitch):
    '''
    Create the xarray datasets that are filled in the ensgen routine

    Inputs: 
    obscoords - List of 2-element lists containing the x,y coordinates
                of interest on the OSGB eastings/northings grid. 
    ensmems - List of 0-padded 2-digit strings corresponding to the ensemble 
              member numbers
    timlen - Length of the time dimension
    stdswitch - 1 if we need to create an xarray to store the 'stdev' in.
                This is only needed if we only have 1 ensemble member
    
    Outputs:
    GAI_p_all, tmean_p_all, prec_p_all, solarrad_p_all, 
    Jarray_p_all, Cday_p_all, GSS_p_all, AWC_allp1:
    Xarray datasets with ensmem,t dimensions, one variable for each obscoord,
    full of zeros.
    '''
    
    times = np.arange(0, timlen)
    ensmemsint = [int(ensmem) for ensmem in ensmems]
    GAI_dict = {}
    extra_dict = {}
    tmean_dict = {}
    tmin_dict = {}
    tmax_dict = {}
    prec_dict = {}
    solarrad_dict = {}
    Jarray_dict = {}
    Cday_dict = {}
    GSS_dict = {}
    AWC_dict = {}
    if stdswitch==1:
        GAI_p_ensstd_dict = {}
    for tob in obscoords:
        GAI_dict[str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((len(ensmems), timlen)), coords=[ensmemsint, times], dims=['ensmem', 'time'])
        extra_dict[str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((len(ensmems), timlen)), coords=[ensmemsint, times], dims=['ensmem', 'time'])
        tmean_dict[str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((len(ensmems), timlen)), coords=[ensmemsint, times], dims=['ensmem', 'time'])
        tmin_dict[str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((len(ensmems), timlen)), coords=[ensmemsint, times], dims=['ensmem', 'time'])
        tmax_dict[str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((len(ensmems), timlen)), coords=[ensmemsint, times], dims=['ensmem', 'time'])
        prec_dict[str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((len(ensmems), timlen)), coords=[ensmemsint, times], dims=['ensmem', 'time'])
        solarrad_dict[str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((len(ensmems), timlen)), coords=[ensmemsint, times], dims=['ensmem', 'time'])
        Jarray_dict[str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((len(ensmems), timlen)), coords=[ensmemsint, times], dims=['ensmem', 'time'])
        Cday_dict[str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((len(ensmems), timlen)), coords=[ensmemsint, times], dims=['ensmem', 'time'])
        GSS_dict[str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((len(ensmems), timlen), dtype='object'), coords=[ensmemsint, times], dims=['ensmem', 'time'])
        if stdswitch==1:
            GAI_p_ensstd_dict[str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((len(ensmems), timlen)), coords=[ensmemsint, times], dims=['ensmem', 'time'])
        AWC_dict[str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((1, timlen)).squeeze(), coords=[times], dims=['time'])
    GAI_p_all = xr.Dataset(GAI_dict)
    extra_p_all = xr.Dataset(extra_dict)
    tmean_p_all = xr.Dataset(tmean_dict)
    tmin_p_all = xr.Dataset(tmean_dict)
    tmax_p_all = xr.Dataset(tmean_dict)
    prec_p_all = xr.Dataset(prec_dict)
    solarrad_p_all = xr.Dataset(solarrad_dict)
    Jarray_p_all = xr.Dataset(Jarray_dict)
    Cday_p_all = xr.Dataset(Cday_dict)
    GSS_p_all = xr.Dataset(GSS_dict)
    AWC_allp1 = xr.Dataset(AWC_dict)
    if stdswitch==1:
        GAI_p_ensstd = xr.Dataset(GAI_p_ensstd_dict)

    if stdswitch==1:
        return GAI_p_all, extra_p_all, tmean_p_all, tmin_p_all, tmax_p_all, prec_p_all, solarrad_p_all, Jarray_p_all, Cday_p_all, GSS_p_all, AWC_allp1, GAI_p_ensstd
    else:
        return GAI_p_all, extra_p_all, tmean_p_all, tmin_p_all, tmax_p_all, prec_p_all, solarrad_p_all, Jarray_p_all, Cday_p_all, GSS_p_all, AWC_allp1

def load_drivingdata(loaddata, ensmem, ensmems, year, Rdatpath, CO2file):
    '''
    Wrapper around the R function that reads in the R-formatted driving data 
    needed for the crop model. This function puts it into an xarray dataset
    format for use in the rest of the assimilation code. 

    Inputs:
    loaddata - The R function object of the R function that reads in the data (r['loaddata'])
    ensmem - The ensemble member we're reading in (string)
    ensmems - List of all the ensemble members (strings)
    year - Growing year to load data for (int)
    Rdatpath - Directory in which the R driving data files are stored (string)
    CO2file - Location of the file containing the CO2 concentration data (string)

    Outputs:
    tmean, prec, solarrad, tmax, tmin - Xarray data arrays of dimension y,x,time
    temp_cconc - R formatted vector containing the CO2 concentration data
    AWC1 - Xarray data array of dimension (y,x,time). Only outputted when ensmem==ensmems[0]   
    '''

    datalist = loaddata(ensmem, year, Rdatpath, CO2file)

    temp_cconc  = datalist.rx2('temp_cconc')
    xall = np.array(datalist.rx2('X'))
    yall = np.array(datalist.rx2('Y'))
    if ensmem==ensmems[0]:
        AWC = datalist.rx2('AWC')
        AWC = np.array(AWC)
        AWC1 = AWC.transpose(1,0,2)
        AWC1 = xr.DataArray(AWC1, coords=[('y', yall), ('x', xall), ('time', np.arange(0, 360))])
    
    
    tmean    = np.array(datalist.rx2('temp_tas')).transpose(1,0,2)
    prec     = np.array(datalist.rx2('temp_pr')).transpose(1,0,2)
    solarrad = np.array(datalist.rx2('temp_rss')).transpose(1,0,2)
    tmax     = np.array(datalist.rx2('temp_tasmax')).transpose(1,0,2)
    tmin     = np.array(datalist.rx2('temp_tasmin')).transpose(1,0,2)
    tmean    = xr.DataArray(tmean, coords=[('y', yall), ('x', xall), ('time', np.arange(0, 360))])
    prec     = xr.DataArray(prec, coords=[('y', yall), ('x', xall), ('time', np.arange(0, 360))])
    solarrad = xr.DataArray(solarrad, coords=[('y', yall), ('x', xall), ('time', np.arange(0, 360))])
    tmax     = xr.DataArray(tmax, coords=[('y', yall), ('x', xall), ('time', np.arange(0, 360))])
    tmin     = xr.DataArray(tmin, coords=[('y', yall), ('x', xall), ('time', np.arange(0, 360))])
    del datalist

    if ensmem==ensmems[0]:
        return tmean, prec, solarrad, tmax, tmin, temp_cconc, AWC1
    else:
        return tmean, prec, solarrad, tmax, tmin, temp_cconc

## NOT USED ##
def fakeobsgen(obscoords, moddata):
    '''
    Generate some fake observation-like data at specified coords from the model data
    Selects out the nearest model grid point, reduces it, adds some random noise
    and randomly removes days from the timeseries.
    Inputs: obscoords - List of 2-element lists containing the x,y coords of the 
                        desired points in eastings,northings metres. 
            moddata - xarray dataarray of model data on an x,y,t grid

    Outputs: testobspdall - Pandas dataframe containing timeseries at each obscoords
    '''
    
    for tob in obscoords:    
        # select out nearest model grid point to obscoords
        testobs = moddata.sel(x=tob[0], y=tob[1], method='nearest')

        # multiply each time point by a random number between 0.6 and 0.85
        multiplier = [np.random.uniform(low=0.55, high=0.9) for t in range(0, len(testobs['time']))]

        # knock out all but every eighth day of obs
        ynoi = np.zeros((len(testobs['time'])//8, 8))
        ynoi[:,0] = 1
        ynoi = list(ynoi.reshape(8*(len(testobs['time'])//8)))
        exti = len(testobs['time']) - len(ynoi)
        exto = [0. for a in range(0,exti)]
        ynoi.extend(exto)
        yesnoobs = ynoi
        #yesnoobs = [round(np.random.uniform()) for t in range(0, len(testobs['time']))]

        # generate the timeseries
        testobs2 = testobs * multiplier * yesnoobs

        # set zeros to missing (nans)
        testobs2 = testobs2.where(testobs2 != 0)

        # convert to pandas series
        testobspd = testobs2.to_pandas()

        # add coordinates as col names
        testobspd = testobspd.rename(str(tob[0]) + ',' + str(tob[1]))

        # merge into dataframe
        if tob == obscoords[0]:
            testobspdall = testobspd
        else:
            testobspdall = pd.merge(testobspdall, testobspd, on='time')
            
    return testobspdall


## NOT USED ##
def nudg_method(obsall, moddata, mod_ensstd, obscoords):
    '''
    A basic DA method that nudges the model values towards the observed 
    based on the relative sizes of the assumed errors/stdevs of each.
    A stdev at each time step for the modelled values should be provided.
    The stdev for the errors is assumed to be 10% of the obs values
    
    Inputs: obsall - Pandas dataframe containing the timeseries of the 
                     observations at various x,y points (output of fakeobsgen)
            moddata - xarray dataarray on an x,y,t grid, containing the model 
                      data
            mod_ensstd - xarray dataarray on an x,y,t grid, containing the model 
                         stdev
            obscoords - List of 2-element lists containing the x,y coords of the 
                        obs points in eastings,northings metres. 

    Outputs: mergedall - Pandas dataframe containing the merged timeseries at the obscoords
             cfall - Pandas dataframe containing the factor needed to generate the merged
                     timeseries from the model one, for each obscoord. 
    '''
    
    for tob in obscoords:    
        # obs point
        obs = obsall[str(tob[0]) + ',' + str(tob[1])]

        # model point
        # do we do the ensemble nudging with all the model ensemble members or just 1?
        # using just 1 for now, as a test case. 
        mod = moddata.sel(x=tob[0], y=tob[1], method='nearest')

        # Smooth the obs timeseries in time, running mean?, 
        # so that we have a timeseries with one observation per model timestep.

        #testobs.plot()
        obsfilled = obs.interpolate(method='linear', limit_direction='both')
        #testobsfilled.plot()

        obssmooth = obsfilled.rolling(10).mean().interpolate(method='linear', limit_direction='both')
        #testobssmooth.plot()

        # Calculate the errors - model stdev at that timestep, and ~10% stdev for obs,
        # or a better value found in literature

        # model stdev for each timestep
        modstd = mod_ensstd.sel(x=tob[0], y=tob[1], method='nearest')

        # calc 10% stdev for obs
        obsstd = obssmooth*0.05

        #uppermod = testmod + modstd
        #lowermod = testmod - modstd
        #upperobs = testobssmooth + obsstd
        #lowerobs = testobssmooth - obsstd
        #plt.figure(figsize=(20,6))
        #ax1 = plt.subplot(1,2,1)
        #testmod.plot(ax=ax1)
        #uppermod.plot(ax=ax1)
        #lowermod.plot(ax=ax1)
        #ax2 = plt.subplot(1,2,2)
        #testobssmooth.plot(ax=ax2)
        #upperobs.plot(ax=ax2)
        #lowerobs.plot(ax=ax2)

        # Calculate a weighted mean of the two values/a correction factor to multiply the model value by,
        # using the above calculated errors

        obsweights = obsstd/obsstd
        modweights = modstd/obsstd
        
        # the timeseries with the highest stdev should have the smallest weight
        # the '1 - ...' ensures this
        obsnormweights = 1 - (obsweights/(obsweights+modweights))
        modnormweights = 1 - (modweights/(obsweights+modweights))

        allweights = np.zeros((len(obsweights), 2))
        allweights[:,0] = obsnormweights
        allweights[:,1] = modnormweights

        alldata = np.zeros((len(obsweights), 2))
        alldata[:,0] = obssmooth
        alldata[:,1] = mod

        merged = np.average(alldata, axis=1, weights=allweights)

        if tob == obscoords[0]:
            mergedall = pd.DataFrame(np.zeros((len(merged), len(obscoords))), index=obsall.index, columns=list(obsall.columns.values))
            cfall = pd.DataFrame(np.zeros((len(merged), len(obscoords))), index=obsall.index, columns=list(obsall.columns.values))
        mergedall[str(tob[0]) + ',' + str(tob[1])] = merged
        #ax1 = plt.subplot(1,1,1)
        #plt.plot(merged)
        #plt.plot(testobssmooth.values)
        #plt.plot(testmod.values)

        # Store each multiplicative factor used for each timestep
        cf = merged/mod
        cfall[str(tob[0]) + ',' + str(tob[1])] = cf

    return mergedall, cfall


## NOT USED ##
allcosts=[]
fig1 = []
def vari_method(obsall, moddata, mod_ensstd, obscoords, year, ensmem, obserrtype=0, obserr=0.1, tsvar=0.25, order=1, power=10,plot=0, nudgedtsall=None):
    '''
    Merge each observation timeseries with it's nearest model grid point timeseries,
    and keep it smooth, using a variational, cost function approach. This method
    attempts to minimize the cost function which consists of terms that grow in magnitude 
    the further the assimilated/guessed timeseries is away from a) the modelled GAI,
    b) the observations on the days that they exist and c) a smoothed version of the 
    timeseries itself, to find the optimum timeseries given these constraints.

    Inputs:
    obsall: A pandas table of the observation data. Rows corresponding to time and 
            columns to different locations. Column headings in format x,y (eastings,
            northings). Size(nxl)
    moddata: xarray dataarray on an x,y,t grid, containing the model 
             data
    mod_ensstd: xarray dataarray on an x,y,t grid, containing the model 
                stdev
    obscoords: List of 2-element lists containing the x,y coords of the 
               obs points in eastings,northings metres. 
    obserrtype: If 0, obs stdev is actual, i.e. obsstd = obserr for each ts
                If 1, obs stdev is relative, i.e. obsstd = obserr*obs for each ts
                0 by default.
    obserr: The obs stdev. Either actual value for each ts, or a relative multiplier
            between 0 and 1.]
            0.1 by default.
    tsvar: this is approximately the max value you think GAI should
           be changing from one timestep (day) to the next when changing by 1stdev,
           e.g 0.25 (default). The smoothing term will try to constrain the assimilated
           series to not vary more than this day-to-day.
           This is the maximum this value can be, for each term in the smoothed
           timeseries it will vary from 0 to this value depending on value of 
           the model timeseries relative to its max. 
           It also needs to have dimensions of length 1 less than the 
           length of the timeseries if we are only considering first-order
           differences. 2 less if considering 2nd order diffs and so on. 
    power: Smoothing matrix power. The greater this number the smoother the smoothed version of the
           timeseries (which makes up one of the three constraints of the cost function) will be.
           This has the effect of increasing the relative proportion of this term compared to the 
           others in the cost function, meaning this term will start to become prioritised in the
           minimisation. 10 by default.
    order: Order of differences to consider when doing the smoothing. 1st order differences will 
           only consider the differences from one timestep either side of each point, 2nd order
           differences, two steps, etc. 1 by default.

    Outputs:
    mergedall: Pandas dataframe containing the merged timeseries at the obscoords
    cfall: Pandas dataframe containing the factor needed to generate the merged
           timeseries from the model one, for each obscoord. 
    '''
    
    for tob in obscoords:
        print('Processing data and generating covariance matrices for coords ' + str(tob) + ' for  year ' + str(year) + ' ensmem ' + str(ensmem))

        # extract out the nudged timeseries for comparison when plotting
        try:
            if nudgedtsall:
                nudgedts = None
        except ValueError:
            nudgedts = nudgedtsall[str(tob[0]) + ',' + str(tob[1])].values
        
        # remove the zeros from either end of the modvars time series. 
        # and produce the indices needed to chop the ends of all the
        # timeseries. The assimilation routine will fail if zero stdevs
        # are used for the model or obs timeseries.
        modstds_trim, frontchop, backchop = remove_zeros(mod_ensstd.sel(x=tob[0], y=tob[1], method='nearest').values)
        # fill in any zeros in the middle of the timeseries
        modstds_trim = pd.Series(modstds_trim).where(pd.Series(modstds_trim)!=0).interpolate(method='linear').values
        
        # select out the obs point from the fake obs table
        obs = obsall[str(tob[0]) + ',' + str(tob[1])]

        # chop off the front and end of the timeseries where the modvar is 0
        obsGAI = obs.values[frontchop:-backchop]

        # Generate the matrix that multiplies the GAI timeseries to extract out only the timesteps
        # that also have observations
        nobs = np.where(~np.isnan(obsGAI))[0].shape[0]
        H = np.zeros((nobs, len(obsGAI)))
        rowinds = np.arange(0, nobs)
        colinds = np.where(~np.isnan(obsGAI))[0]
        H[rowinds, colinds] = 1.

        # generate the vector that says where there's an obs value or not (0 for no, 1 for yes)
        yno = np.where(np.isnan(obsGAI), 0, 1.)
        # generate a timeseries of the obs with the timesteps where there are no obs removed
        obsGAI_trim = np.asarray([obsGAI[i] for i in range(0,len(obsGAI)) if yno[i]==1])

        # select out the model grid point nearest the obs point
        modGAI = moddata.sel(x=tob[0], y=tob[1], method='nearest').values

        # chop off the front and end of the timeseries where the modstd is 0
        modGAI = modGAI[frontchop:-backchop]

        # assume the 5% stdev for the obs for now, or a user supplied %
        if obserrtype == 0:
            obsstds_trim = np.array([obserr for a in obsGAI_trim])
        elif obserrtype == 1:
            obsstds_trim = obsGAI_trim*obserr

        # generate the covariance matrices for the obs and the mod
        modcov = np.eye(len(modstds_trim))*(modstds_trim**2)
        obscov = np.eye(len(obsstds_trim))*(obsstds_trim**2)

        # generate the covariance matrix for the smoothed term
        # this is approximately the max value you think GAI should
        # be changing from one timestep (day) to the next. 
        # e.g 0.25. The smoothing term will try to constrain the assimilated
        # series to not vary more than this day-to-day.
        # It also needs to have dimensions of length 1 less than the 
        # length of the timeseries if we are only considering first-order
        # differences. 2 less if considering 2nd order diffs and so on. 
        tsvar = tsvar
        order = order
        #smoothcov = np.abs((modGAI_trim - smoother(modGAI_trim, 1000., 1))*np.eye(len(modGAI_trim)))
        modGAI_norm = modGAI/modGAI.max()
        modGAI_norm = np.where(modGAI_norm < 0.01, 0.01, modGAI_norm)
        tsvari = tsvar*modGAI_norm
        smoothcov = tsvari[:-order]**2*np.eye(len(modGAI)-order)


        # calculate the inverses of the covariance matrices
        obscovinv = np.linalg.inv(obscov)
        modcovinv = np.linalg.inv(modcov)
        smoothcovinv = np.linalg.inv(smoothcov)

        # set the first guess for the minimisation routine
        firstguess = modGAI

        # set bounds
        #bounds = [(a*0.5, a*1.5) for a in list(nudgedts[frontchop:-backchop])]
        
        # run the minimisation routine and plot the evolution of the cost function
        print('Running minimisation routine')
        global allcosts
        global fig1
        allcosts=[]
        fig1 = []
        #x,f,d = spo.fmin_l_bfgs_b(cost_function, firstguess, args=(obsGAI_trim2, modGAI_trim, obscovinv, modcovinv, smoothcovinv, yno), approx_grad=True, disp=5, maxfun=100000, pgtol=2e-02, epsilon=1e-04, factr=1e12, bounds=bounds)
        merged = spo.fmin_bfgs(cost_function, firstguess, grad_cost_function, args=(obsGAI_trim, obsGAI, modGAI, obscovinv, modcovinv, smoothcovinv, H, power, order, plot, nudgedts))
        # append and prepend zeros that were cut off
        #print(merged)
        merged = list(merged)
        merged.extend([0 for a in range(0, backchop)])
        zeros = [0 for a in range(0, frontchop)]
        zeros.extend(merged)
        merged = np.asarray(zeros)

        modGAI = list(modGAI)
        modGAI.extend([0 for a in range(0, backchop)])
        zeros = [0 for a in range(0, frontchop)]
        zeros.extend(modGAI)
        modGAI = np.asarray(zeros)

        if tob == obscoords[0]:
            mergedall = pd.DataFrame(np.zeros((len(merged), len(obscoords))), index=obsall.index, columns=list(obsall.columns.values))
            cfall = pd.DataFrame(np.zeros((len(merged), len(obscoords))), index=obsall.index, columns=list(obsall.columns.values))
        mergedall[str(tob[0]) + ',' + str(tob[1])] = merged

        # Store each multiplicative factor used for each timestep
        cf = merged/modGAI
        cfall[str(tob[0]) + ',' + str(tob[1])] = cf

        # save figure
        #if not os.path.exists('plots'):
        #    os.makedirs('plots')
        #fig1.savefig('plots/vari_year_' + str(year) + '_ensmem_' + str(ensmem) + '_loc_' + str(tob[0]) + ',' + str(tob[1]) + '.png', dpi=300)
        
    return mergedall, cfall


allcosts=[]
fig1=[]
def vari_method_point(obsall, moddata, mod_ensstd, assimvar, obscoords, year, ensmem, obserrtype=0, obserr=0.1, moderrinfl=1, tsvar=0.25, order=1, power=10, plot=0, plotdir='None', interval=10, nudgedtsall=None):
    '''
    Merge each observation timeseries with it's nearest model grid point timeseries,
    and keep it smooth, using a variational, cost function approach. This method
    attempts to minimize the cost function which consists of terms that grow in magnitude 
    the further the assimilated/guessed timeseries is away from a) the modelled GAI,
    b) the observations on the days that they exist and c) a smoothed version of the 
    timeseries itself, to find the optimum timeseries given these constraints.

    Inputs:
    obsall: A pandas table of the observation data. Rows corresponding to time and 
            columns to different locations. Column headings in format x,y (eastings,
            northings). Size(nxl)
    moddata: xarray dataset containing the model (GAI) data, one variable per coord
             Each variable containing the timeseries. 
    mod_ensstd: As moddata but containing the ensemble stdev instead
    assimvar: Name of the variable to be assimilated
    obscoords: List of 2-element lists containing the x,y coords of the 
               obs points in OSGB eastings,northings metres. 
    obserrtype: If 0, obs stdev is actual, i.e. obsstd = obserr for each ts
                If 1, obs stdev is relative, i.e. obsstd = obserr*obs for each ts
                0 by default.
    obserr: The obs stdev. Either actual value for each ts, or a relative multiplier
            between 0 and 1.
            0.1 by default.
    moderrinfl: Value to multiply the model error by. Defaults to 1. 
    tsvar: this is approximately the max value you think GAI should
           be changing from one timestep (day) to the next when changing by 1stdev,
           e.g 0.25 (default). The smoothing term will try to constrain the assimilated
           series to not vary more than this day-to-day.
           This is the maximum this value can be, for each term in the smoothed
           timeseries it will vary from 0 to this value depending on value of 
           the model timeseries relative to its max. 
    power: Smoothing matrix power. The greater this number the smoother the smoothed version of the
           timeseries (which makes up one of the three constraints of the cost function) will be.
           This has the effect of increasing the relative proportion of this term compared to the 
           others in the cost function, meaning this term will start to become prioritised in the
           minimisation. 10 by default.
    order: Order of differences to consider when doing the smoothing. 1st order differences will 
           only consider the differences from one timestep either side of each point, 2nd order
           differences, two steps, etc. 1 by default.
    plot: Switch to control plotting of timeseries and cost function convergence. May only work
          when used in a jupyter notebook, hence zero by default.
    nudgedtsall: For use in the above plotting, shows another timeseries for comparison.
                 Only used if plot==1, defaults to None. 

    Outputs:
    mergedall: Pandas dataframe containing the merged timeseries at the obscoords
    cfall: Pandas dataframe containing the factor needed to generate the merged
           timeseries from the model one, for each obscoord. 
    '''
    
    for tob in obscoords:
        print('Processing data and generating covariance matrices for coords ' + str(tob) + ' for  year ' + str(year) + ' ensmem ' + str(ensmem))

        # extract out the nudged timeseries for comparison when plotting
        try:
            if not nudgedtsall:
                nudgedts = None
        except ValueError:
            nudgedts = nudgedtsall[str(tob[0]) + ',' + str(tob[1])].values

        # select out the model grid point nearest the obs point
        modGAI = moddata[str(tob[0]) + ',' + str(tob[1])].values
        
        # remove the zeros from either end of the mod time series. 
        # and produce the indices needed to chop the ends of all the
        # timeseries. The assimilation routine will fail if zero stdevs
        # are used for the model or obs timeseries.
        #modGAI_trim, frontchop, backchop = remove_zeros(modGAI)
        modGAI_trim = modGAI

        # select out the obs point from the obs table
        obs = obsall[str(tob[0]) + ',' + str(tob[1])]

        # chop off the front and end of the timeseries where the modvar is 0
        #obsGAI = obs.values[frontchop:-backchop]
        obsGAI = obs.values

        # implement the obs error on the full length obs list
        if obserrtype == 0:
            obsstds = np.array([obserr if a!=0 else 0 for a in obsGAI])
        elif obserrtype == 1:
            obsstds = obsGAI*obserr

        # Generate the matrix that multiplies the GAI timeseries to extract out only the timesteps
        # that also have observations
        nobs = np.where(~np.isnan(obsGAI))[0].shape[0]
        H = np.zeros((nobs, len(obsGAI)))
        rowinds = np.arange(0, nobs)
        colinds = np.where(~np.isnan(obsGAI))[0]
        H[rowinds, colinds] = 1.

        # generate the vector that says where there's an obs value or not (0 for no, 1 for yes)
        yno = np.where(np.isnan(obsGAI), 0, 1.)
        # generate a timeseries of the obs with the timesteps where there are no obs removed
        obsGAI_trim = np.asarray([obsGAI[i] for i in range(0,len(obsGAI)) if yno[i]==1])

        modstds_trim = mod_ensstd[str(tob[0]) + ',' + str(tob[1])].values#[frontchop:-backchop]
        # set modstds to 0.001 where is modGAI==0
        #modstds_trim = list(np.where(modGAI_trim<=0., 0.001, modstds_trim))
        # set any remaining stds of 0s to 0.001
        modstds_trim2 = np.array([0.001 if x==0. else x for x in modstds_trim])
        # multiply the mod error by the model error inflation value (default 1)
        modstds_trim = modstds_trim2 * moderrinfl

        # implement the obs error
        if obserrtype == 0:
            obsstds_trim = np.array([obserr for a in obsGAI_trim])
        elif obserrtype == 1:
            obsstds_trim = obsGAI_trim*obserr

        # generate the covariance matrices for the obs and the mod
        modcov = np.eye(len(modstds_trim))*(modstds_trim**2)
        obscov = np.eye(len(obsstds_trim))*(obsstds_trim**2)
        print(modcov.dtype)
        print(obscov.dtype)

        # generate the covariance matrix for the smoothed term
        # tsvar is approximately the max value you think GAI should
        # be changing from one timestep (day) to the next. 
        # e.g 0.25. The smoothing term will try to constrain the assimilated
        # series to not vary more than this day-to-day.
        # It also needs to have dimensions of length 1 less than the 
        # length of the timeseries if we are only considering first-order
        # differences. 2 less if considering 2nd order diffs and so on. 
        tsvar = tsvar
        order = order
        #smoothcov = np.abs((modGAI_trim - smoother(modGAI_trim, 1000., 1))*np.eye(len(modGAI_trim)))
        modGAI_norm = modGAI_trim/modGAI_trim.max()
        modGAI_norm = np.where(modGAI_norm < 0.01, 0.01, modGAI_norm)
        tsvari = tsvar*modGAI_norm
        #maxind = tsvari.argmax()
        #tsvari[maxind-20:] = tsvar
        #tsteps = []
        #for a in range(1, len(modGAI)+1):
        #    b = a/len(modGAI)
        #    tsteps.append(b)
        #tsvari = tsvar*np.asarray(tsteps)
        smoothcov = tsvari[:-order]**2*np.eye(len(modGAI_trim)-order)


        # calculate the inverses of the covariance matrices
        obscovinv = np.linalg.inv(obscov)
        modcovinv = np.linalg.inv(modcov)
        smoothcovinv = np.linalg.inv(smoothcov)

        # set the first guess for the minimisation routine
        firstguess = modGAI_trim

        # set bounds
        #bounds = [(a*0.5, a*1.5) for a in list(nudgedts[frontchop:-backchop])]
        
        # run the minimisation routine and plot the evolution of the cost function
        print('Running minimisation routine')
        global allcosts
        global fig1
        allcosts=[]
        fig1=[]
        coord = str(tob[0]) + ',' + str(tob[1])
        #x,f,d = spo.fmin_l_bfgs_b(cost_function, firstguess, args=(obsGAI_trim2, modGAI_trim, obscovinv, modcovinv, smoothcovinv, yno), approx_grad=True, disp=5, maxfun=100000, pgtol=2e-02, epsilon=1e-04, factr=1e12, bounds=bounds)
        merged = spo.fmin_bfgs(cost_function, firstguess, grad_cost_function, args=(obsGAI_trim, obsGAI, modGAI_trim, obscovinv, modcovinv, smoothcovinv, H, power, order, plot, plotdir, coord, interval, ensmem, nudgedts))

        # calculate error on the posterior
        post_err = posterior_error(merged, H, obscovinv, modcovinv, power, smoothcovinv, order)
        #print(post_err)

        # plot the final result of the assimilation
        plt.rcParams['font.size'] = '16'
        obserrs = []
        for a in range(0, len(obsGAI)):
            if obsGAI[a]-obsstds[a] > 0:
                obserr1 = [[a,a,a], [obsGAI[a]-obsstds[a], obsGAI[a], obsGAI[a]+obsstds[a]]]
            else:
                obserr1 = [[a,a,a], [0, obsGAI[a], obsGAI[a]+obsstds[a]]]
            obserrs.append(obserr1)
        if plot==1 or plot==2:
            fig, ax1 = plt.subplots(1,1)
            fig.set_size_inches(8,8)
            ax1.plot(obsGAI, 'bx', label='Observed LAI', zorder=2)
            ax1.plot(modGAI, 'r', label='Modelled GAI', zorder=4)
            ax1.plot(merged, 'k', label='Assimilated GAI', zorder=3)
            upperposterr = np.array(merged) + np.array(post_err[0])
            lowerposterr = np.array(merged) - np.array(post_err[0])
            uppermoderr = np.array(modGAI_trim) + np.array(modstds_trim)
            lowermoderr = np.array(modGAI_trim) - np.array(modstds_trim)
            lowermoderr = np.where(lowermoderr<=0, 0, lowermoderr)
            #ax1.plot(upperposterr[0], 'g-', zorder=4)
            #ax1.plot(lowerposterr[0], 'g-', zorder=4)
            #ax1.plot(uppermoderr[0], 'r-', zorder=4)
            #ax1.plot(lowermoderr[0], 'r-', zorder=4)
            x=np.arange(0, len(merged))
            ax1.fill_between(x, upperposterr[0], lowerposterr[0], 'k', zorder=1, alpha=0.35) 
            ax1.fill_between(x, uppermoderr, lowermoderr, 'r', zorder=0, alpha=0.35) 
            for a in range(0, len(obserrs)):
                ax1.plot(obserrs[a][0], obserrs[a][1], 'b-', linewidth=0.5)
            ax1.set_xlabel('Day of growing year')
            ax1.set_ylabel('LAI/GAI')
            plt.legend(loc='upper left')
            if plotdir!='None':
                plt.savefig(os.path.join(plotdir, str(tob[0])+'_'+str(tob[1])+'_'+str(ensmem)+'_final.png'), dpi=300)
            #display.clear_output(wait=True)
            #display.display(fig)

        # append and prepend zeros that were cut off
        #print(merged)
        #merged = list(merged)
        #merged.extend([0 for a in range(0, backchop)])
        #zeros = [0 for a in range(0, frontchop)]
        #zeros.extend(merged)
        #merged = np.asarray(zeros)

        #modGAI_trim = list(modGAI_trim)
        #modGAI_trim.extend([0 for a in range(0, backchop)])
        #zeros = [0 for a in range(0, frontchop)]
        #zeros.extend(modGAI_trim)
        #modGAI_trim = np.asarray(zeros)

        #post_err = list(np.array(post_err[0])[0])
        #post_err.extend([0 for a in range(0, backchop)])
        #zeros = [0 for a in range(0, frontchop)]
        #zeros.extend(post_err)
        #post_err = np.asarray(zeros)

        if tob == obscoords[0]:
            mergedall = pd.DataFrame(np.zeros((len(merged), len(obscoords))), index=obsall.index, columns=list(obsall.columns.values))
            cfall = pd.DataFrame(np.zeros((len(merged), len(obscoords))), index=obsall.index, columns=list(obsall.columns.values))
            post_err_all = pd.DataFrame(np.zeros((len(merged), len(obscoords))), index=obsall.index, columns=list(obsall.columns.values))
        mergedall[str(tob[0]) + ',' + str(tob[1])] = merged

        # Store each multiplicative factor used for each timestep
        cf = merged/modGAI_trim
        cfall[str(tob[0]) + ',' + str(tob[1])] = cf

        # Store each posterior error vector
        post_err_all[str(tob[0]) + ',' + str(tob[1])] = np.array(post_err[0])[0]
        
    return mergedall, cfall, post_err_all


## NOT USED ##
def calc_GAI(ensmem, year, Rdatpath, CO2file):
    '''
    Calculates the GAI for a given driving-data-ensemble-member and year

    Inputs: ensmem - 0-padded 2-digit string corresponding to the ensemble 
                     member number
            year -  Int/float of the year of the data
            Rdatpath - Location of .Rdata driving data files (string)
            CO2file - Location of file containing CO2 conc. data (string)

    Outputs: GAIds - Xarray dataset containing the GAI (x,y,t)
             tmean, prec, solarrad, temp_cconc, dat_sp, Jarray, Cday, GSS, HarvestJday, x, y, t
             driving data and derived variables from the R script, needed to
             calculate yield later.
    '''
    
    print('Loading data for ensmem ' + ensmem + ' for gyear ' + str(year))

    loaddata = r['loaddata']
    GAIfunc = r['GAI']
    
    # read in driving data for this growing year from the .Rdata files
    datalist = loaddata(ensmem, year, Rdatpath, CO2file)

    temp_tas = datalist.rx2('temp_tas')
    temp_rls = datalist.rx2('temp_rls')
    temp_rss = datalist.rx2('temp_rss')
    temp_pr  = datalist.rx2('temp_pr')
    temp_tasmax = datalist.rx2('temp_tasmax')
    temp_tasmin = datalist.rx2('temp_tasmin')
    temp_cconc  = datalist.rx2('temp_cconc')

    # run crop growth model to calculate GAI
    print('Calculating GAI over season')
    datalist2 = GAIfunc(temp_tas, temp_tasmax, temp_tasmin, temp_pr, temp_rss)

    GAI         = np.array(datalist2.rx2('GAI')) # convert to python numpy
    tmean       = datalist2.rx2('tmean') # all others can stay as R type as not editing these
    prec        = datalist2.rx2('prec')
    solarrad    = datalist2.rx2('solarrad')
    dat_sp      = datalist2.rx2('dat_SP')
    Jarray      = datalist2.rx2('Jarray')
    Cday        = datalist2.rx2('Cday')
    GSS         = datalist2.rx2('GSS')
    #GSS         = np.array(datalist2.rx2('GSS')) # As it's a string array it comes out as a flattened array
    #GSS         = GSS.reshape((GSS_r.dim[2], GSS_r.dim[1], GSS_r.dim[0]))
    #GSS         = GSS.transpose(2,1,0) # therefore we need to manually reshape it
    HarvestJday = datalist2.rx2('HarvestJday')

    # get the data coordinates
    x = np.array(datalist2.rx2('x'))
    y = np.array(datalist2.rx2('y'))
    t = np.array(datalist2.rx2('t'))


    # Convert GIA to xarray to make the space and time handling/wrangling easier
    GAI = GAI.transpose(1,0,2)
    #print('GAI data dims are ' + str(GAI.shape))
    attrs = {'units': 'days since ' + t[0], 'calendar': '360_day'}
    GAIxr = xr.DataArray(GAI, coords=[('y', y), ('x', x), ('time', np.arange(0, len(t)))])
    GAIxr.name = 'GAI'
    GAIds = GAIxr.to_dataset()
    GAIds.time.attrs = attrs
    GAIds = xr.decode_cf(GAIds, decode_coords=False)
    GAIds = GAIds.assign_coords(x=x)
    GAIds = GAIds.assign_coords(y=y) # the cf time decoding sometimes mangles the coords, so reset them
    
    return GAIds, tmean, prec, solarrad, temp_cconc, dat_sp, Jarray, Cday, GSS, HarvestJday, x, y, t

## NOT USED ##
def update_GAI_points(GAIold, obscoords, mergedall):
    '''
    Replace gridpoints of the modeled GAI data with the merged obs&model timeseries.
    The nearest model gridpoints are used. 

    Inputs: GAIold - Xarray dataset containing the modelled GAI (x,y,t)
            obscoords - List of 2-element lists containing the x,y coords of the 
                        obs points in eastings,northings metres. 
            mergedall - Pandas dataframe containg the timeseries of the merged timeseries

    Outputs: GAInew - Xarray dataset with the relevant points updated with the merged data
    '''
    
    GAInew = GAIold['GAI'].copy()
    for tob in obscoords:
        nearestx = GAInew.sel(x=tob[0], y=tob[1], method='nearest')['x'].values
        nearesty = GAInew.sel(x=tob[0], y=tob[1], method='nearest')['y'].values
        GAInew.loc[dict(x=nearestx, y=nearesty)] = mergedall[str(tob[0]) + ',' + str(tob[1])].values
    return GAInew

## NOT USED ##
def update_yield_points(GAIold, GAInew, obscoords, tmean, prec, solarrad, Jarray, Cday, GSS, HarvestJday, AWC, x, y, temp_cconc, ensmem):
    '''
    Run updated array of GAI through the model to produce updated array of yield

    Inputs: GAIold - Original modelled GAI. Xarray datarray of dim (x,y,t)
            GAInew - New merged/assimilated GAI  
                     Pandas dataframe containing the merged timeseries at each
                     obscoord. 
            The rest as before.

    Outputs: oldyield - Original modelled yield at the modelled grid points nearest obscoords
             newyield - Updated yield at the same points
                        Both as a list

    '''
    
    #wheat_yield = r['wheat_yield']
    wheat_yield_point = r['wheat_yield_point']

    newyields = []
    oldyields = []
    counter=1
    totalobs = len(obscoords)
    for tob in obscoords:
        print('Processing obs point ' + str(counter) + ' of ' + str(totalobs))
        counter+=1
        # Select out the timeseries at the obs points, or the nearest model grid point
        GAInew1 = GAInew[str(tob[0])+','+str(tob[1])].values
        GAIold1 = GAIold.sel(x=tob[0], y=tob[1], method='nearest').values
        tmean1 = tmean.sel(x=tob[0], y=tob[1], method='nearest').values
        prec1 = prec.sel(x=tob[0], y=tob[1], method='nearest').values
        solarrad1 = solarrad.sel(x=tob[0], y=tob[1], method='nearest').values
        Jarray1 = Jarray.sel(x=tob[0], y=tob[1], method='nearest').values
        Cday1 = Cday.sel(x=tob[0], y=tob[1], method='nearest').values
        GSS1 = GSS.sel(x=tob[0], y=tob[1], method='nearest').values
        AWC1 = AWC.sel(x=tob[0], y=tob[1], method='nearest').values
        
        #print(type(GAInew))
        #print(type(tmean))
        #print(type(prec))
        #print(type(solarrad))
        #print(type(AWC))
        #print(type(Jarray))
        #print(type(Cday))
        #print(type(GSS))
        #print(type(HarvestJday))
        #print(type(x))
        #print(type(y))
        #print(type(temp_cconc))
        yieldlist = wheat_yield_point(GAInew1, tmean1, prec1, solarrad1, AWC1, Jarray1, Cday1, GSS1, HarvestJday, x, y, cconc=temp_cconc)
        oldyieldlist = wheat_yield_point(GAIold1, tmean1, prec1, solarrad1, AWC1, Jarray1, Cday1, GSS1, HarvestJday, x, y, cconc=temp_cconc)

        # Extract out the grid points we replaced (using the same xr, nearest, method)
        
        yieldval = np.asarray(yieldlist.rx2('z')).copy()[0]
        oldyieldval = np.asarray(oldyieldlist.rx2('z')).copy()[0]
        #x = np.asarray(yieldlist.rx2('x')).copy()
        #y = np.asarray(yieldlist.rx2('y')).copy()
        #yieldvals = yieldvals.transpose(1,0)
        #oldyieldvals = oldyieldvals.transpose(1,0)
        #yieldxr = xr.DataArray(yieldvals, coords=[('y', y), ('x', x)])
        #oldyieldxr = xr.DataArray(oldyieldvals, coords=[('y', y), ('x', x)])

        newyields.append(yieldval)
        oldyields.append(oldyieldval)
        
    return oldyields, newyields

def update_yield_points_point(priorold, priornew, obscoords, assimvar, nonprior, tmean, tmin, tmax, prec, solarrad, Jarray, Cday, GSS, HarvestJday, AWC, CDD, TT, temp_cconc, ensmem):
    '''
    Run a series of GAI timeseries through the crop model to calculate the corresponding yields

    Inputs: GAIold - Original modelled GAI. Xarray dataset. One variable (timeseries) per obscoord
            GAInew - New merged/assimilated GAI. Pandas dataframe containing the merged timeseries 
                     at each obscoord. 
            tmean, prec, solarrad, Jarray, Cday, GSS - Driving data/ancilaries needed for yield
                                                       calculation. Same format as GAI_old. 
            HarvestJday, AWC - Further ancillaries, but these without an ensemble member dimension.
            The rest as in above functions. 

    Outputs: oldyield - Original modelled yield at the modelled grid points nearest obscoords
             newyield - Updated yield at the same points
                        Both as a list
    '''
    
    #wheat_yield = r['wheat_yield']
    wheat_yield_point = r['wheat_yield_point']

    newyields = []
    oldyields = []
    counter=1
    totalobs = len(obscoords)
    for tob in obscoords:
        print('Processing obs point ' + str(counter) + ' of ' + str(totalobs))
        counter+=1
        # Select out the timeseries at the obs points, or the nearest model grid point
        priornew1 = priornew[str(tob[0])+','+str(tob[1])].values
        priorold1 = priorold[str(tob[0])+','+str(tob[1])].values
        nonprior1 = nonprior[str(tob[0])+','+str(tob[1])].values
        tmean1 = tmean[str(tob[0])+','+str(tob[1])].values
        tmin1 = tmean[str(tob[0])+','+str(tob[1])].values
        tmax1 = tmean[str(tob[0])+','+str(tob[1])].values
        prec1 = prec[str(tob[0])+','+str(tob[1])].values
        solarrad1 = solarrad[str(tob[0])+','+str(tob[1])].values
        Jarray1 = Jarray[str(tob[0])+','+str(tob[1])].values
        Cday1 = Cday[str(tob[0])+','+str(tob[1])].values
        GSS1 = GSS[str(tob[0])+','+str(tob[1])].values
        AWC1 = AWC[str(tob[0])+','+str(tob[1])].values
        
        #print(type(GAInew))
        #print(type(tmean))
        #print(type(tmin))
        #print(type(tmax))
        #print(type(prec))
        #print(type(solarrad))
        #print(type(AWC))
        #print(type(Jarray))
        #print(type(Cday))
        #print(type(GSS))
        #print(type(HarvestJday))
        #print(type(x))
        #print(type(y))
        #print(type(temp_cconc))
        if assimvar == 'LAI':
            yieldlist = wheat_yield_point(priornew1, nonprior1, assimvar, tmean1, tmin1, tmax1, prec1, solarrad1, AWC1, Jarray1, Cday1, GSS1, HarvestJday, CDD, TT, cconc=temp_cconc)
            oldyieldlist = wheat_yield_point(priorold1, nonprior1, assimvar, tmean1, tmin1, tmax1, prec1, solarrad1, AWC1, Jarray1, Cday1, GSS1, HarvestJday, CDD, TT, cconc=temp_cconc)
        elif assimvar == 'fPAR':
            yieldlist = wheat_yield_point(nonprior1, priornew1, assimvar, tmean1, tmin1, tmax1, prec1, solarrad1, AWC1, Jarray1, Cday1, GSS1, HarvestJday, CDD, TT, cconc=temp_cconc)
            oldyieldlist = wheat_yield_point(nonprior1, priorold1, assimvar, tmean1, tmin1, tmax1, prec1, solarrad1, AWC1, Jarray1, Cday1, GSS1, HarvestJday, CDD, TT, cconc=temp_cconc)

        yieldval = np.asarray(yieldlist.rx2('finalyield')).copy()[0]
        oldyieldval = np.asarray(oldyieldlist.rx2('finalyield')).copy()[0]

        newyields.append(yieldval)
        oldyields.append(oldyieldval)

    del yieldlist
    del oldyieldlist
    
    return oldyields, newyields


# newGAI is the only value that changes with each iteration. The others are all fixed
def cost_function(newGAI, obsGAI, obsGAInontrim, modGAI, obscovinv, modcovinv, smoothcovinv, H, power, order, plot, plotdir, coord, interval, ensmem, nudgedts):
    '''
    The cost function to minimize, consists of observation, model and
    smoothing constraints/terms.
    Inputs: 
    newGAI: The vector whose value that minimizes this function we're iteratively 
            solving for. Size (1xn)
    obsGAI: A vector of observations, with missing timesteps removed. Size (1xp)
    obsGAInontrim: A vector of observations with the missing timesteps (nans)
                   kept in. Size (1xn)
    modGAI: A vector of the modelled values. Size (1xn)
    obscovinv: The inverted matrix of observation covariances. Before inversion, this 
               matrix will typically be a diagonal matrix, with each term the variance
               (uncertainty) of the observations at that timestep (may all be the same). 
               Size (pxp)
    modcovinv: The inverted matrix of model timeseries covariances. As obscovinv, the
               pre-inverted matrix will typically be diagonal, with each term the variance
               (uncertainty) of the modelled value at each timestep. Size (nxn)
    smoothcovinv: The inverted matrix of the smoothed model covariances. Pre-inversion,
                  this will again likely be a diagonal matrix, with each value the variance
                  (uncertainty) of the smoothed value at each timestep. An easier-to-understand
                  interpretation of these values is that they are how the squared value of
                  the stdev of the distribution within which you would like the values to change
                  in the smoothed timeseries from {order} timestep(s) to the next. Size (nxn)
    H: The matrix that removes the timesteps of newGAI that don't have observations when multiplied
       by it. Size (pxn): Each column should be all zeros, except for columns corresponding to
       timesteps where there are observations, in which case there should be a 1 in the row that
       corresponds to the position in the obsGAI vector that this observation is. 
    power: Smoothing matrix power. The greater this number the smoother the smoothed version of the
           timeseries (which makes up one of the three constraints of the cost function) will be.
           This has the effect of increasing the relative proportion of this term compared to the 
           others in the cost function, meaning this term will start to become prioritised in the
           minimisation. 
    order: Order of differences to consider when doing the smoothing. 1st order differences will 
           only consider the differences from one timestep either side of each point, 2nd order
           differences, two steps, etc. 
    plot, nudgedts: See vari_method_point(...)

    Outputs:
    Cost: Value of the cost function, given the above inputs. 
    '''
    
    #start = time.time()
    obs_cost = obs_cost_func(newGAI, obsGAI, obscovinv, H)
    mod_cost = mod_cost_func(newGAI, modGAI, modcovinv)
    smooth_cost = smooth_cost_func(newGAI, smoothcovinv, power, order)
    cost = mod_cost + obs_cost + smooth_cost
    cost = np.asarray(cost)[0][0]
    if plot > 0:
        global allcosts
        allcosts.append(cost)
        plot_costs(allcosts, newGAI, obsGAInontrim, modGAI, nudgedts, plot, plotdir, coord, ensmem, interval)
    #end = time.time()
    #print(cost)
    return cost

def grad_cost_function(newGAI, obsGAI, obsGAInontrim, modGAI, obscovinv, modcovinv, smoothcovinv, H, power, order, plot, plotdir, coord, interval, ensmem, nudgedts):
    '''
    The gradient of the cost function with respect to newGAI. Requires all the same inputs 
    as cost_function() (even though many are not used) for the purposes of the 
    minimisation routine. 
    '''
    
    grad_obs_cost = grad_obs_cost_func(newGAI, obsGAI, obscovinv, H)
    grad_mod_cost = grad_mod_cost_func(newGAI, modGAI, modcovinv)
    grad_smooth_cost = grad_smooth_cost_func(newGAI, smoothcovinv, power, order)
    grad_cost = grad_obs_cost + grad_mod_cost + grad_smooth_cost
    return grad_cost

def obs_cost_func(newGAI, obsGAI, obscovinv, H):
    '''
    The component of the cost function whose value is proportional to the
    difference between newGAI and obsGAI.
    Inputs:
    newGAI: The vector whose value that minimizes this function we're iteratively 
            solving for. Size (1xn)
    obsGAI: A vector of observations, with missing timesteps removed. Size (1xp)
    obscovinv: The inverted matrix of observation covariances. Before inversion, this 
               matrix will typically be a diagonal matrix, with each term the variance
               (uncertainty) of the observations at that timestep (may all be the same). 
               Size (pxp)
    H: The matrix that removes the timesteps of newGAI that don't have observations when multiplied
       by it. Size (pxn): Each column should be all zeros, except for columns corresponding to
       timesteps where there are observations, in which case there should be a 1 in the row that
       corresponds to the position in the obsGAI vector that this observation is. 
    '''
    
    Hx = np.matmul(H, newGAI)
    diff = obsGAI - Hx
    diffTC = np.matmul(diff[:,None].T, obscovinv)
    obs_cost = 0.5 * np.matmul(diffTC, diff)
    return obs_cost


def grad_obs_cost_func(newGAI, obsGAI, obscovinv, H):
    '''
    The gradient of the obs term in the cost function, with respect to newGAI
    Inputs:
    newGAI: The vector whose value that minimizes this function we're iteratively 
            solving for. Size (1xn)
    obsGAI: A vector of observations, with missing timesteps removed. Size (1xp)
    obscovinv: The inverted matrix of observation covariances. Before inversion, this 
               matrix will typically be a diagonal matrix, with each term the variance
               (uncertainty) of the observations at that timestep (may all be the same). 
               Size (pxp)
    H: The matrix that removes the timesteps of newGAI that don't have observations when multiplied
       by it. Size (pxn): Each column should be all zeros, except for columns corresponding to
       timesteps where there are observations, in which case there should be a 1 in the row that
       corresponds to the position in the obsGAI vector that this observation is. 
    '''
    
    Hx = np.matmul(H, newGAI)
    HTC = np.matmul(H.T, obscovinv)
    diff = obsGAI - Hx
    grad_obs_cost = -1.*np.matmul(HTC, diff)
    return np.asarray(grad_obs_cost)


def mod_cost_func(newGAI, modGAI, modcovinv):
    '''
    The component of the cost function that is proportional to the difference between 
    newGAI and modGAI.
    Inputs:
    newGAI: The vector whose value that minimizes this function we're iteratively 
            solving for. Size (1xn)
    modGAI: A vector of the modelled values. Size (1xn)
    modcovinv: The inverted matrix of model timeseries covariances. As obscovinv, the
               pre-inverted matrix will typically be diagonal, with each term the variance
               (uncertainty) of the modelled value at each timestep. Size (nxn)    
    '''
    
    mod_cost = 0.5 * np.matmul((newGAI-modGAI).T, np.matmul(modcovinv, newGAI-modGAI))
    return mod_cost


def grad_mod_cost_func(newGAI, modGAI, modcovinv):
    '''
    The gradient of the model term of the cost function with respect to newGAI.
    Inputs:
    newGAI: The vector whose value that minimizes this function we're iteratively 
            solving for. Size (1xn)
    modGAI: A vector of the modelled values. Size (1xn)
    modcovinv: The inverted matrix of model timeseries covariances. As obscovinv, the
               pre-inverted matrix will typically be diagonal, with each term the variance
               (uncertainty) of the modelled value at each timestep. Size (nxn)   
    '''
    
    grad_mod_cost = np.matmul(modcovinv, (newGAI - modGAI))
    return np.asarray(grad_mod_cost)


def smooth_cost_func(newGAI, smoothcovinv, power, order):
    '''
    The component of the cost function that is proportional to the difference
    between newGAI and a smoothed version of newGAI. 
    Inputs:
    newGAI: The vector whose value that minimizes this function we're iteratively 
            solving for. Size (1xn)
    smoothcovinv: The inverted matrix of the smoothed model covariances. Pre-inversion,
                  this will again likely be a diagonal matrix, with each value the variance
                  (uncertainty) of the smoothed value at each timestep. An easier-to-understand
                  interpretation of these values is that they are how the squared value of
                  the stdev of the distribution within which you would like the values to change
                  in the smoothed timeseries from {order} timestep(s) to the next. Size (nxn)
    power: Smoothing matrix power. The greater this number the smoother the smoothed version of the
           timeseries (which makes up one of the three constraints of the cost function) will be.
           This has the effect of increasing the relative proportion of this term compared to the 
           others in the cost function, meaning this term will start to become prioritised in the
           minimisation. 
    order: Order of differences to consider when doing the smoothing. 1st order differences will 
           only consider the differences from one timestep either side of each point, 2nd order
           differences, two steps, etc. 
    '''
    
    x=newGAI # timeseries to smooth
    N=len(x)
    
    #order = order of difference
    #(int; 1=1st differences, 2=2nd differences etc)
    D=np.matrix(np.diff(np.eye(N, dtype='int64'),n=order, axis=0)*-1)
    Dx = np.matmul(D, x[:, None])
    smoothedx = np.matmul(Dx.T, smoothcovinv)
    smoothcost = (power**2)*np.matmul(smoothedx, Dx)
    
    #plot
    #plt.plot(x, label='modGAI')
    #plt.plot(smoothedx.T, label='smooth modGAI')
    #plt.legend()
    #plt.title('g = ' + str(g) + ', power = ' + str(power))
    #print(str(power))
    #print(str(smoothcost))
    return smoothcost


def grad_smooth_cost_func(newGAI, smoothcovinv, power, order):
    '''
    The gradient of the smoothing component of the cost function with respect to newGAI.
    Inputs:
    newGAI: The vector whose value that minimizes this function we're iteratively 
            solving for. Size (1xn)
    smoothcovinv: The inverted matrix of the smoothed model covariances. Pre-inversion,
                  this will again likely be a diagonal matrix, with each value the variance
                  (uncertainty) of the smoothed value at each timestep. An easier-to-understand
                  interpretation of these values is that they are how the squared value of
                  the stdev of the distribution within which you would like the values to change
                  in the smoothed timeseries from {order} timestep(s) to the next. Size (nxn)
    power: Smoothing matrix power. The greater this number the smoother the smoothed version of the
           timeseries (which makes up one of the three constraints of the cost function) will be.
           This has the effect of increasing the relative proportion of this term compared to the 
           others in the cost function, meaning this term will start to become prioritised in the
           minimisation. 
    order: Order of differences to consider when doing the smoothing. 1st order differences will 
           only consider the differences from one timestep either side of each point, 2nd order
           differences, two steps, etc. 
    '''
    
    x=newGAI # timeseries to smooth
    N=len(x)
    
    #order = order of difference
    #(int; 1=1st differences, 2=2nd differences etc)
    D=np.matrix(np.diff(np.eye(N, dtype='int64'),n=order, axis=0)*-1)
    Dx = np.matmul(D, x[:, None])
    CDx = np.matmul(smoothcovinv, Dx)
    grad_smooth_cost = (power**2)*np.matmul(D.T, CDx)
    return np.asarray(grad_smooth_cost).squeeze()

def posterior_error(newGAI, H, obscovinv, modcovinv, power, smoothcovinv, order):
    
    N = len(newGAI)
    D=np.matrix(np.diff(np.eye(N, dtype='int64'),n=order, axis=0)*-1)

    obs_comp = np.matmul(H.T, np.matmul(obscovinv, H))
    mod_comp = modcovinv
    smooth_comp = power**2 * np.matmul(D.T, np.matmul(smoothcovinv, D))
    
    grad_grad_cost_func = obs_comp + mod_comp + smooth_comp
    posterior_error_matrix = np.linalg.inv(grad_grad_cost_func)
    posterior_error = posterior_error_matrix.diagonal()
    return np.sqrt(posterior_error)

def plot_costs(allcosts, newGAI, obsGAI, modGAI, nudgedts, plot, plotdir, coord, ensmem, interval):
    '''
    Plot the value of the cost function on one plot and the current value of newGAI,
    along with obsGAI and modGAI on another. Called by cost_function(...) every interval
    iterations. 
    Inputs:
    allcosts: A list of all the values of the cost function from each iteration so far
    newGAI: Current value
    obsGAI: Obs timeseries with missing values (nans) retained
    modGAI: Model timeseries
    nudgedts: Nudged (or any) timeseries for comparison
    '''

    global fig1
    if len(allcosts) >= 1:
        fig1 = plt.gcf()
        
    if len(allcosts) == interval:
        if plot == 2:
            ax1, ax2 = setup_plot(plot)
            if max(allcosts) > 100000:
                ylim2 = 100000
            else:
                ylim2 = max(allcosts)*1.1
            ax1.set_ylim([0, ylim2])
            ax1.plot(allcosts, 'b')
        elif plot == 1:
            ax2 = setup_plot(plot)

        try:
            if nudgedts:
                pass
        except ValueError:
            ax2.plot(nudgedts, 'k', label='nudged')    
        ax2.plot(obsGAI, 'bx', label='Observed LAI')
        ax2.plot(modGAI, 'r', label='Modelled GAI')
        ax2.plot(newGAI, 'g', label='Assimilated GAI')
        ax2.set_xlabel('Day of growing year')
        ax2.set_ylabel('GAI/LAI')
        #ax2.set_ylims([0, ax2.get_ylims()[1]])
        plt.legend()
    
    if len(allcosts) % interval == 0:
        fig1 = plt.gcf()
        if plot == 2:
            ax1 = fig1.get_axes()[0]
            ax2 = fig1.get_axes()[1]
            if len(allcosts) > 1000 and len(allcosts)%1000 == 100:
                ax1.set_xlim([0, len(allcosts)+900])
            if len(ax1.get_lines()) > 0:
                lines = ax1.get_lines()
                lines.pop(0).remove()
            ax1.plot(allcosts, 'b')
        elif plot == 1:
            ax2 = fig1.get_axes()[0]
        #ax.set_ylim([0, max(allcosts)*1.1])
        if len(ax2.get_lines()) > 0:
            lines = ax2.get_lines()
            lines.pop(-1).remove()
        ax2.plot(newGAI, 'g', label='Assimilated GAI')
        #ax2.set_ylims([0, ax2.get_ylims()[1]])
        display.clear_output(wait=True)
        display.display(fig1)
        if plotdir != 'None':
            itern = '{:05.0f}'.format(len(allcosts))
            plotname = coord + '_' + str(ensmem) + '_iter' + itern + '.png'
            plotpath = os.path.join(plotdir, plotname)
            plt.savefig(plotpath, dpi=300)
    return

        
def setup_plot(plot):
    '''
    Called by plot_costs(...) to set up the plot the first time it is called
    '''
    
    if plot == 2:
        fig, (ax1, ax2) = plt.subplots(1,2)
        fig.set_size_inches(15, 8, forward=True)
        ax1.set_xlim([0, 1000])
        return ax1, ax2
    elif plot == 1:
        fig, ax1 = plt.subplots(1,1)
        fig.set_size_inches(8,8)
        return ax1


def remove_zeros(ts):
    '''
    Remove any zeros from the beginning or end of the supplied timeseries,
    and return this along with the indices needed to chop the same timesteps
    off any other timeseries of the same length using ts[frontchop:-backchop]
    '''
    
    frontchop = 0
    backchop  = len(ts)
    if ts[0] == 0:
        frontchop = np.where(ts!=0)[0].min()
    if ts[-1] == 0:
        backchop = np.where(ts!=0)[0].max()+1
    newts = ts[frontchop:backchop]
    return newts, frontchop, len(ts)-backchop


def plot_1box_yield(dat,pos,col,it):
    # just draw box (no whiskers etc.), and make colours etc. right
    #bp = plt.boxplot(dat,whis=0,sym='',positions=[pos],widths=0.37,patch_artist=True)
    bp = plt.boxplot(dat,whis=0,sym='',positions=[pos],widths=0.5,patch_artist=True)   
    plt.setp(bp['medians'],color='red')
    for pt in ['whiskers','caps']:
        plt.setp(bp[pt],color=col,lw=0)
    #if(it < 2):
    plt.setp(bp['boxes'], alpha=0.75)
    #else:
    #    plt.setp(bp['boxes'],facecolor='None',color=col)
    # add my own whiskers, plus overall min and max
    minvals = np.min(dat[:])
    p10s = np.percentile(dat[:],10)
    p25s = np.percentile(dat[:],25)
    p75s = np.percentile(dat[:],75)
    p90s = np.percentile(dat[:],90)
    maxvals = np.max(dat[:])
    size=40
    plt.scatter(pos,minvals,marker='_',color=col,s=size) # draw minima
    plt.scatter(pos,maxvals,marker='_',color=col,s=size) # draw maxima
    plt.scatter(pos,p10s,marker='_',color=col,s=size)    # draw lower caps
    plt.scatter(pos,p90s,marker='_',color=col,s=size)    # draw upper caps
    plt.vlines(pos,p25s,p10s,color=col,lw=1.5)           # draw lower whiskers
    plt.vlines(pos,p75s,p90s,color=col,lw=1.5)           # draw upper whiskers
    return bp

def plot_boxes_yield(data,stats,posis,xlabs, ylabs, ylims, title, actyield=0, savefile='plot.png'):
    '''
    Adapted from g2g stats code, hence some funny variable names and weird formatting

    Create 2 box plots per plot.
    Used for comparing the original modelled yield to the assimilated.

    Inputs:
    data: Pandas dataframe containing columns with the names in stats
    stats: Names of 'stats' to plot. Should correspond to column names in data
           List of lists. Each element of the main list corresponds to a 
           particular plot, with the sublist a pair of names corresponding to
           the boxplots on the plot in question.
           E.g. [['NS daily', 'NS hrly'], ['Bias daily', 'Bias hrly']] to compare the
                stats of two datasets.
    posis: x axis positions of the boxplots on each plot. e.g. [0,1] for 2 boxplots
    xlabs: x axis labels for each plot e.g. ['oldyield', 'newyield']
    ylabs: ylabel for each plot e.g. ['tonnes/hectare']
    ylims: y axis limits for each plot
    title: Plot title
    actyield: Value to plot as horizontal line across plots
    savefile: Direcotry and name of file to save plot as. Set to None if using notebook. 

    Outputs:
    datameanold: The mean of the left hand box plot data (ensemble mean of original yield)
    datameannew: The mean of the right hand box plot data (ensemble mean of asssimilated yield)
    For use in the verify(...) function. 
    '''
    
    # loop over stats
    counter=0
    for stat in stats:
        #print(counter,stat)
        fig = plt.figure(figsize=(5,4))
        fig.add_subplot(1,1,1)
        ## add boxes
        # boxplot1
        bp0 = plot_1box_yield(data[stat[0]],posis[0],'0.5',0)
        data25 = np.percentile(data[stat[0]],25)
        datamed = np.percentile(data[stat[0]],50)
        datameanold = np.mean(data[stat[0]])
        data75 = np.percentile(data[stat[0]],75)
        # boxplot2
        bp1 = plot_1box_yield(data[stat[1]],posis[1],'red',1)
        datameannew = np.mean(data[stat[1]])
        #print('25thptile:',data25,np.percentile(data[stat[1]],25))
        #print('med:',datamed,np.percentile(data[stat[1]],50))
        #print('75thptile:',data75,np.percentile(data[stat[1]],75))
        
        ## add lines, labels etc
        if len(actyield) == 1:
            plt.axhline(y=actyield[0], ls='-.',c='0.4',lw=1.5)
        else:
            if len(posis) > 2:
                bp2 = plot_1box_yield(actyield, posis[2], 'black', 2)
            else:
                print('Not plotting 3rd box plot because position (posis) not supplied')
        
        plt.xlim(posis[0]-0.5,posis[-1]+0.5)
        plt.xticks(posis,xlabs,fontsize=12)
        for xl in np.arange(0.5,0,1):
            plt.axvline(x=xl,ls=':',c='0.85') # add vertical line at x=xl
        plt.ylim(ylims[counter])
        plt.yticks(fontsize=12)
        plt.ylabel(ylabs[counter],fontsize=12)
        counter+=1
    # finish fig
    plt.title(title, fontsize=14)
    if savefile:
        plt.savefig(savefile, dpi=300)
    #plt.legend([bp0["boxes"][0],bp1["boxes"][0]],['data','dscalepe'],fontsize=11,bbox_to_anchor=(0.2,-0.06),ncol=3,loc=1)
    #plt.subplots_adjust(left=0.1,right=0.98,top=0.98,bottom=0.08,wspace=0.36,hspace=0.2)
    return datameanold, datameannew


def verifyyield(yieldfile, obscoords, oldyields_vari, newyields_vari, years, ensmem, savedir, allplots=0, newyields_vari_alt=None):
    '''
    Produce plots comparing the original modelled yield against the assimilated and
    observed yield values. 

    yieldfile: File containing the yield data as a shape file
    obscoords: List of 2-element lists containing the x,y coordinates of the observed yield data
    oldyields_vari: Xarray dataset containing the original modelled yields. 
                    One variable per obscoord, dims: ensmems, years
    newyields_vari: As oldyields_vari but for the assimilated yield
    years: List of years of we ran the assimilation routine for.  
    savedir: Location to save the plots. Set to None if using notebook.
    '''
    if savedir:
        if not os.path.exists(savedir):
            os.makedirs(savedir)
    
    # read in the actual yield data
    yielddata = gpd.read_file(yieldfile)
    # sort by descending field area
    for areaname in ['area', 'area_ha']:
        try:
            asortyielddata = yielddata.sort_values([areaname])[::-1]
            aname = areaname
            break
        except KeyError:
            continue
    # create pandas geodataframe for the obscoords points (which are taken from the satellite data)
    obscoordspd = pd.DataFrame(obscoords, columns=['x','y'])
    obspoints = gpd.GeoDataFrame(obscoordspd, geometry=gpd.points_from_xy(obscoordspd.x, obscoordspd.y))
    obspoints.crs = yielddata.crs
    # for each point in obspoints, found out which field from the yielddata it is inside,
    # and attach the yields to this coord. If a point is not inside any field, it is dropped
    actyields = gpd.sjoin(obspoints, asortyielddata, op = 'within')
    actyields['x'] = [str(x) for x in actyields.x]
    actyields['y'] = [str(y) for y in actyields.y]
    actyields['xy'] = actyields['x'] + ',' + actyields['y']
    actyields = actyields.set_index('xy')
    #print(actyields)

    nmissing = len(obscoords) - len(actyields)
    print(str(nmissing) + ' obs points not inside a field/with no yield data')

    if newyields_vari_alt:
        newyields2 = newyields_vari_alt
    else:
        newyields2 = newyields_vari
    # for each obs point...
    meanolds = []
    meannews = []
    obsyields=[]
    fields = []
    nobs = 0 # used to count no. of data points used in verification
    for tob in list(actyields.index.values):
        # check if we have yield data for each year we ran for
        # and if so, plot 
        for year in years:
            if np.isnan(actyields.loc[tob][str(year)]) == False:
                # check if we have model data
                if np.any(np.isnan(newyields2[tob].sel(year=year))) == False:
                    fields.append(tob)
                    nobs+=1
                    actyield = actyields.loc[tob][str(year)]
                    newyield = newyields_vari[tob].sel(year=year).to_pandas()
                    oldyield = oldyields_vari[tob].sel(year=year).to_pandas()
                    newyield.name = 'newyield'
                    oldyield.name = 'oldyield'
                    oldandnew = pd.merge(oldyield, newyield, on='ensmem')
                    if allplots:
                        if not savedir:
                            meanold, meannew = plot_boxes_yield(oldandnew, [['oldyield','newyield']], [0,1], ['oldyield', 'newyield'], ['yield (tn/hec)'], [[7,17]], str(year) + ': ' + tob, actyield, None)
                        else:
                            savefile = os.path.join(savedir, str(year) + '_' + tob + '.png')
                            meanold, meannew = plot_boxes_yield(oldandnew, [['oldyield','newyield']], [0,1], ['oldyield', 'newyield'], ['yield (tn/hec)'], [[7,17]], str(year) + ': ' + tob, actyield, savefile)
                    else:
                        #print(oldandnew)
                        if ensmem:
                            meanold = oldandnew['oldyield'].loc[int(ensmem)]
                            meannew = oldandnew['newyield'].loc[int(ensmem)]
                        else:
                            meanold = np.mean(oldandnew['oldyield'])
                            meannew = np.mean(oldandnew['newyield'])
                    meanolds.append(meanold)
                    meannews.append(meannew)
                    obsyields.append(actyields.loc[tob][str(year)])

    nfields = len(np.unique(np.asarray(fields)))
    print('Total number of fields used in verification: ' + str(nfields))
    print('Total data points used in verification: ' + str(nobs))
    # plot mean over gridpoints and years
    #meanolds = np.asarray(meanolds)
    #meannews = np.asarray(meannews)
    #meanolds = np.where(meanolds==0., np.nan, meanolds)
    #meannews = np.where(meannews==0., np.nan, meannews)
    meanoldspd = pd.DataFrame(meanolds)
    meannewspd = pd.DataFrame(meannews)
    meanoldspd.index.name = 'point'
    meannewspd.index.name = 'point'
    meanoldspd.columns = ['mod_yield']
    meannewspd.columns = ['assim_yield']
    meanmeanold = meanoldspd.mean()[0]
    meanmeannew = meannewspd.mean()[0]
    print('Mean modelled yield: ' + str(meanmeanold))
    print('Mean assimilated yield: ' + str(meanmeannew))
    oldsandnews = pd.merge(meanoldspd, meannewspd, on='point')
    #yieldsf = gpd.read_file(yieldfile)
    #yieldsf.drop(['field_id', aname, 'geometry'], axis=1, inplace=True)
    #allyields = yieldsf.values.flatten()
    #yields = allyields[~np.isnan(allyields)]
    yields = obsyields
    #actyields.drop(['x', 'y', 'geometry', 'index_right', 'field_id', 'area'], axis=1, inplace=True)
    #actyields.index = oldsandnews.index
    #oldsnewsacts = pd.merge(oldsandnews, actyields, on='point')
    #actyield = actyields.mean(axis=0).mean()
    #print(actyield)

    savefile = os.path.join(savedir, 'mean_yield_BC.png')
    junk1, junk2 = plot_boxes_yield(oldsandnews, [['mod_yield','assim_yield']], [0,1,2], ['mod_yield', 'assim_yield', 'actual_yield'], ['yield (' + r'$t\: ha^{-1}$' + ')'], [[7,17]], 'Yields', yields, savefile)

    return fields


def download_data(outdir, years, variables, lonmin, lonmax, latmin, latmax):

    downloadera5(outdir, years, variables, lonmin, lonmax, latmin, latmax)


def outputsave(data, coords, dims, year, name, units, outfile):
    data = xr.DataArray(data, coords=coords, dims=dims)
    if len(coords)==2:
        data = data.expand_dims({'time': [str(year+1)]})
    data.name = name
    data.attrs = {'units': units}
    data.to_netcdf(outfile)

def country_subset_shapefile(data=None, datafile=None, multifile=0, sfname=None, xname='x', yname='y', IDname=None, IDs=None, drop=True):
    '''
    Function to subset an xarray dataarray or dataset, or netcdf dataset, to selected shapes from 
    a shapefile. Returns an xarray dataset with of the same shape as the input
    datafile but with the data outside the selected shapes
    set to nans. Also returns the shapes so these can be plotted.

    data:     An xarray DataArray or DataSet
    datafile: The filename of the netcdf file to subset. Multiple files can be selected with * etc.
              If this is the case multifile should be set to 1. Defaults to 0.
    sfname:   The filename of the shapefile
    IDname:   The name of the catgeory to search over for selecting shapes to subset to (e.g. 'RIVER')
    IDs:      The values of the category to select (e.g. ['Thames', 'Severn'])
    multifile:Are multiple files specfied in datafile? Set to 1 if so. In this case the files are read in
              using dask, which can process data that exceeds the memory capacity of the machine by
              processing in parallel.
    xname: Name of the x-coordinate in the netcdf file(s). 'x' by default.
    yname: Name of the y-coordinate in the netcdf file(s). 'y' by default
    '''

    # Read in data
    if datafile:
        print('Reading in ' + datafile)
        if multifile == 1:
            data = xr.open_mfdataset(datafile, parallel=True)
        else:
            data = xr.load_dataset(datafile)

    subset = add_shape_coord_from_data_array(data, sfname, IDname, IDs, yname, xname)
    if drop==True:
        subset = subset.where(subset[IDname]==1, drop=True)
    else:
        subset = subset.where(subset[IDname]==1)
        
        
    return subset



def add_shape_coord_from_data_array(xr_da, shp_path, IDname, IDs, latname, lonname):
    """ Create a new coord for the xr_da indicating whether or not it 
    is inside the shapefile
    
    Creates a new coord - "coord_name" which will have integer values
    used to subset xr_da for plotting / analysis/
    
    Usage:
    -----
    precip_da = add_shape_coord_from_data_array(precip_da, "awash.shp", "awash")
    awash_da = precip_da.where(precip_da.awash==0, other=np.nan) 
    """
    # 1. read in shapefile
    shp_gpd = gpd.read_file(shp_path)
    
    # 2. create a list of tuples (shapely.geometry, id)
    #    this allows for many different polygons within a .shp file (e.g. States of US)
    shapes = []
    counter = 0
    for ID in shp_gpd[IDname]:
        if ID in IDs:
            shapes.append((shp_gpd['geometry'][counter], 1))
            print('Found: ' + str(shp_gpd[IDname][counter]))
        counter+=1

    if len(shapes) == 0:
        raise AttributeError(IDname + ' ' + str(IDs) + ' not found in shapefile')

    xr_da[IDname] = rasterize(shapes, xr_da.coords,
                              longitude=lonname, latitude=latname)
    
        
    return xr_da


def transform_from_latlon(lat, lon):
    """ 
    input 1D array of lat / lon and output an Affine transformation
    """

    lat = np.asarray(lat)
    lon = np.asarray(lon)
    trans = Affine.translation(lon[0], lat[0])
    scale = Affine.scale(lon[1] - lon[0], lat[1] - lat[0])
    return trans * scale

def rasterize(shapes, coords, latitude='y', longitude='x',
              fill=np.nan, **kwargs):
    """
    Rasterize a list of (geometry, fill_value) tuples onto the given
    xarray coordinates. This only works for 1d latitude and longitude
    arrays.

    usage:
    -----
    1. read shapefile to geopandas.GeoDataFrame
          `states = gpd.read_file(shp_dir+shp_file)`
    2. encode the different shapefiles that capture those lat-lons as different
        numbers i.e. 0.0, 1.0 ... and otherwise np.nan
          `shapes = (zip(states.geometry, range(len(states))))`
    3. Assign this to a new coord in your original xarray.DataArray
          `ds['states'] = rasterize(shapes, ds.coords, longitude='X', latitude='Y')`

    arguments:
    ---------
    : **kwargs (dict): passed to `rasterio.rasterize` function

    attrs:
    -----
    :transform (affine.Affine): how to translate from latlon to ...?
    :raster (numpy.ndarray): use rasterio.features.rasterize fill the values
      outside the .shp file with np.nan
    :spatial_coords (dict): dictionary of {"X":xr.DataArray, "Y":xr.DataArray()}
      with "X", "Y" as keys, and xr.DataArray as values

    returns:
    -------
    :(xr.DataArray): DataArray with `values` of nan for points outside shapefile
      and coords `Y` = latitude, 'X' = longitude.
    """
                  
    print('Adding mask to xarray')
    transform = transform_from_latlon(coords[latitude], coords[longitude])
    out_shape = (len(coords[latitude]), len(coords[longitude]))
    raster = features.rasterize(shapes, out_shape=out_shape,
                                fill=fill, transform=transform,
                                dtype=float, **kwargs)
    spatial_coords = {latitude: coords[latitude], longitude: coords[longitude]}
    return xr.DataArray(raster, coords=spatial_coords, dims=(latitude, longitude))


def calc_daily_climo(inputs, outname, varname, tname, xname='x', yname='y', roll=None):
    '''
    Function to calculate a daily climatology from several years worth of daily data
    Only really designed for outputs from the wheat model assuming a growing season
    running from 1st Oct. May generalise this later. 
    
    Inputs:
    inputs - String that identifies all the files e.g. files_????.nc
    outname - Output name and folder of climo file, excluding the file extension.
              climo.nc and stdev.nc will be added to the respective output files
    varname - Name of variable in file to calculate climo of
    tname - Name of the time dimension in the input files
    roll - The number of days to calculate a rolling mean of the mean and stdev statistics.
           Default is None which means no rolling mean of the mean or stdev stat is produced.
    '''

    coutname = outname + '_climo.nc'
    soutname = outname + '_stdev.nc'
    if roll:
        rcoutname = outname + '_' + str(roll) + 'dayclimo.nc'
        rsoutname = outname + '_' + str(roll) + 'daystdev.nc'
    
    cdo = Cdo()
    cdo.ydaymean(input='-cat ' + inputs,
                 output=coutname)
    cdo.ydaystd(input='-cat ' + inputs,
                 output=soutname)
    
    tc = xr.load_dataset(coutname)
    tc = tc.rename({tname: 't'})
    tc = tc.roll(t=-274, roll_coords=True)
    newtc = pd.date_range('2003-10-01', freq='D', periods=366)
    tc[varname + '2'] = xr.DataArray(tc[varname].values, [newtc, tc[yname], tc[xname]], ['time', yname, xname])
    tc = tc.drop_vars([varname, 't'])
    tc = tc.rename({varname + '2': varname})
    tc.to_netcdf(coutname)
    if roll:
        tcr = tc.rolling(time=roll, min_periods=1).mean()
        tcr.to_netcdf(rcoutname)

    ts = xr.load_dataset(soutname)
    ts = ts.rename({tname: 't'})
    ts = ts.roll(t=-274, roll_coords=True)
    ts[varname + '2'] = xr.DataArray(ts[varname].values, [newtc, ts[yname], ts[xname]], ['time', yname, xname])
    ts = ts.drop_vars([varname, 't'])
    ts = ts.rename({varname + '2': varname})
    ts.to_netcdf(soutname)
    if roll:
        tsr = ts.rolling(time=roll, min_periods=1).mean()
        tsr.to_netcdf(rsoutname)
