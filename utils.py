import sys
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
from rpy2.robjects import numpy2ri
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 20})
import scipy.optimize as spo
#import cartopy as cp
import datetime as dt
import pandas as pd
import pyproj
import geopandas as gpd
#from rasterio import features
#from affine import Affine
#from dateutil.relativedelta import relativedelta
import netCDF4 as nc4
import cftime as cft
import os
import shutil
import glob
import time
import requests
import urllib.request
#import warnings
#warnings.filterwarnings.ignore()
from IPython import display
numpy2ri.activate()
#import csv
#import shapefile
#from shapely.geometry import Polygon, MultiPolygon, Point, MultiPoint, shape
#from dask.diagnostics import Profiler, ResourceProfiler, CacheProfiler, ProgressBar, visualize
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#import matplotlib.animation as animation

# Set up R environment on import of this module
raster = importr('raster')
rgdal = importr('rgdal')
rgeos = importr('rgeos')
geosphere = importr('geosphere')
abind = importr('abind')
datatable = importr('data.table')
r = robjects.r
r['source']('Lynch_potpredict_v2_MJB.R') # define R functions


def MODIS_request(infile, obscoordsin, startyear, endyear, MODIScode='MCD15A2H'):
    '''
    Submit requests to the MODIS data servers to download data from
    modis.ornl.gov for the centre point of each of the fields we have 
    yield data for, or provided coordinates.
    A maximum of 100 requests can be submitted per email address per day,
    but this can be ignored by using different, fake, email addresses.
    Several email addresses are written in the 'emails' variable list
    below. Just add some more nonsense ones if you need more, substituting
    %40 for @

    Inputs:
    infile: Location of file containing the shapes of the fields 
            we have yield data for, or a csv file containing x,y coords,
            or None
    obscoordsin: List of of 2-element lists containing x,y coords
    startyear: First year you want to download data for
    endyear: Last year you want to download data for
    MODIScode: Product code of MODIS product you want

    Outputs:
    datacodes: List of datacodes for each point that can be used to 
               download the data once ready
    '''

    # check if input file provided
    if infile:
        ocswitch=0
        # check if it's a shape file or not
        if infile.split('.')[-1]=='shp':
            shapeswitch=1
        else:
            shapeswitch=0
    # if no file at all, use passed obscoords
    else:
        oscswitch=1
        shapeswitch=0

    if shapeswitch==1:
        # Read in shapefile containing yield data into geopandas dataframe
        obsyields = gpd.read_file(yieldshapefile)
        
        # Extract out the geometries and find the 'center' coordinates of each field
        obscoordsxy = [point.coords[0] for point in list(obsyields.centroid.values)]

    if shapeswitch==0 and oscswitch==0:
        obscoordspd = pd.read_csv(infile)
        xs = list(obscoordspd.iloc[:, 0].values)
        ys = list(obscoordspd.iloc[:, 1].values)
        obscoordsxy = [[xs[i],ys[i]] for i in range(0, len(xs))]

    if oscswitch==1:
        obscoordsxy=obscoordsin
    
    # Convert these to lon,lat
    proj = pyproj.Transformer.from_crs(27700, 4326, always_xy=True)
    obscoordslonlat = [proj.transform(x,y) for x,y in obscoordsxy]

    # Submit request to MODIS servers for data for each location
    datacodes=[]
    emails=['a%40b.c', 'd%40e.f', 'g%40h.i', 'j%40k.l', 'm%40n.o', 'p%40q.r', 's%40t.u', 'v%40w.x', 'y%40z.z', 'test%40test.c', 'test%40test.b', 'test%40test.a', 'fun%40be.com', 'matt%40me.com', 'matt%40me.be', 'mat%40think.c', 'bla%40blabla.com', 'bladebla%40evenmorebla.com', 'testy%40me.com', 'cdb%40de.com', 'mccool%40as.ice', 'cold%40hot.water', '%check%40mate.mb', 'need%40idea.s']

    for counter in range(0, len(obscoordsxy)):
        lon=str(obscoordslonlat[counter][0])
        lat=str(obscoordslonlat[counter][1])
        productcode = MODIScode
        uid='cnet' + str(counter)
        startdate='A' + str(startyear) + '001'
        enddate='A' + str(endyear) + '365'
        yarea=str(0)
        xarea=str(0)

        headers = {
            'Accept': 'text/csv',
        }

        url = 'https://modis.ornl.gov/rst/api/v1/' + productcode + '/subsetOrder?' + \
              'latitude=' + lat + '&longitude=' + lon + \
              '&email=' + emails[0] + '&uid=' + uid + \
              '&startDate=' + startdate + '&endDate=' + enddate + \
              '&kmAboveBelow=' + yarea + '&kmLeftRight=' + xarea

        print('Requesting: ' + url)
        datacode = requests.get(url, headers=headers).text

        # If requests limit reached (100 per day) try a different email address
        if 'Cannot' in datacode:
            emails.pop(0)
            url = 'https://modis.ornl.gov/rst/api/v1/' + productcode + '/subsetOrder?' + \
                  'latitude=' + lat + '&longitude=' + lon + \
                  '&email=' + emails[0] + '&uid=' + uid + \
                  '&startDate=' + startdate + '&endDate=' + enddate + \
                  '&kmAboveBelow=' + yarea + '&kmLeftRight=' + xarea
            print('Request limit reached, trying new email address: ' + emails[0])
            print('Requesting: ' + url)
            datacode = requests.get(url, headers=headers).text

        datacodes.append(datacode)
        print(datacode)
        time.sleep(2)
    return datacodes

def MODIS_download(datasavedir, datacodes, product_name='Lai_500m'):
    '''
    Download the MODIS data once the requested data is ready.
    
    Inputs:
    datasavedir: Where to save the data files
    datacodes: The data codes obtained from the data requests
    product_name: The MODIS product name as it appears in the url
    '''
    
    if not os.path.exists(datasavedir):
        os.makedirs(datasavedir)

    for dc in datacodes:
        status = 404
        url = 'https://modis.ornl.gov/subsetdata/' + dc + '/filtered_scaled_' + product_name + '.csv'
        lat, lon = dc.split('L')[1], dc.split('L')[2].split('S')[0]
        outfile = url.split('/')[-1].split('.')[0] + '_lat_' + lat + '_lon_' + lon + '.csv'
        savefile = os.path.join(datasavedir, outfile)

        while status==404:
            #status = requests.get(url, headers=headers).status_code
            status = requests.get(url).status_code
            if status!=404:
                break
            print('Data ' + dc + ' not ready, trying again in 2mins')
            time.sleep(120)

        print('Saving ' + dc + ' to ' + savefile)
        urllib.request.urlretrieve(url, savefile)
        time.sleep(1)

def MODIS_process(datasavedir, filter_threshold=0.5):
    '''
    Process the MODIS data from the downloaded CSV files into
    a pandas dataframe usable by the code.
    The data are translated onto a 360day calendar by shifting
    anything on the 31st of a month onto the 30th so that no 
    data is lost. Values below a default threshold of 0.5 are
    removed between 1st March and 1st July as unrealistic. 

    Inputs:
    datasavedir: Where the data are saved and the wildcard filenames
                 that identifies all of the files you want to load.
    filter_threshold: Detailed above.

    Outputs:
    obspdall: Pandas dataframe containing the obs data. Cols are
              each downloaded pixel (named using the x,y coords)
              Rows are days of a 360day calendar (to match the model)
    obscoords: List of 2-element lists containing the [x,y] coords
               of the obs pixels. OSGB eastings,northings. 
    '''
    
    datafiles = glob.glob(datasavedir)

    obscoords=[]
    counter=1
    totalfiles = len(datafiles)
    for filein in datafiles:
        print('Processing ' + str(counter) + ' of ' + str(totalfiles))
        # Read in lon,lat of pixel
        latlon = pd.read_csv(filein, header=None, usecols=[3], nrows=1).values[0][0]
        lat = float(latlon.split('S')[0].split('L')[1][2:])
        lon = float(latlon.split('S')[0].split('L')[2][2:])

        # Convert to OSGB eastings,northings
        proj = pyproj.Transformer.from_crs(4326, 27700, always_xy=True)
        x,y = proj.transform(lon,lat)
        # Set name of column for data table
        colname = str(x).split('.')[0] + ',' + str(y).split('.')[0]
        coord = [int(str(x).split('.')[0]), int(str(y).split('.')[0])]
        obscoords.append(coord)

        # Read in data file and do a whole load of nasty date wrangling
        pixel = pd.read_csv(filein, header=None, index_col=0, usecols=[2,6], names=['date', colname], na_values='F')
        # Convert date index column to datetimes
        pixel.index = pd.to_datetime(pixel.index, format='%Y%j', exact=False)
        # Change any '31st's of the month days to 30ths for 360day calendar
        startdate = str(pixel.index[0])[:10]
        enddate   = str(pixel.index[-1])[:10]
        oldindex = pixel.index.values
        newindex1 = [dt.datetime(int(str(d)[:4]), int(str(d)[5:7]), int(str(d)[8:10])) for d in list(oldindex)]
        newindex2 = [d if d.day<=30 else dt.datetime(d.year, d.month, 30) for d in newindex1]
        pixel.index = pd.DatetimeIndex(newindex2)
        # Convert index to 360day calendar
        datetimes = [dt.datetime.strptime(str(d), '%Y-%m-%dT%H:%M:%S.000000000') for d in pixel.index.values]
        cfdatetimes = [cft.Datetime360Day(d.year, d.month, d.day) for d in datetimes]
        cfdatetimesidx = xr.coding.cftimeindex.CFTimeIndex(cfdatetimes)
        pixel.index = cfdatetimesidx
        # Add in all the days on which there are no obs as NaNs
        alldaysidx = xr.cftime_range(startdate, enddate, calendar='360_day', freq='D', name='date')
        pixel = pixel.reindex(alldaysidx)
        # merge each pixel into one pandas table
        if filein == datafiles[0]:
            obspdall = pixel
        else:
            obspdall = pd.merge(obspdall, pixel, on='date')
        counter+=1

    # filter out observations below a certain threshold between
    # 1st March and 1st July each year as unrealistic
    thresh = filter_threshold
    startyear = obspdall.index.values[0].year + 1
    if obspdall.index.values[-1].month < 7:
        endyear = obspdall.index.values[-1].year - 1
    else:
        endyear = obspdall.index.values[-1].year

    for year in range(startyear, endyear+1):
        seldates = xr.cftime_range(str(year)+'-03-01', str(year)+'-07-30', calendar='360_day', freq='D', name='date').values
        obspdall.loc[seldates] = obspdall.loc[seldates].where(obspdall.loc[seldates]>thresh)

    return obspdall, obscoords

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

def ensgen_point(ensmems, year, Rdatpath, CO2file, obscoords):
    '''
    Generates a set of ensemble members of GAI at specified coordinates,
    or the closest model grid point to them.
    Returns this as an xarray dataset of dimensions time and ensmem, with
    each variable in the dataset corresponding to a different coordinate.
    Also returned are the ensemble mean and stdev, in a similar format.

    Inputs: ensmems - List of 0-padded 2-digit strings corresponding to the ensemble 
                      member numbers
            year - Growing year to work with (int/float). Note this corresponds
                   to the year that harvest is made in.
            Rdatpath - Location of .Rdata driving data files (string)
            CO2file - Location of file containing CO2 conc. data (string)
            obscoords - List of 2-element lists containing the x,y coordinates
                        of interest on the OSGB eastings/northings grid. 
    
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
    loaddata = r['loaddata']
    GAIfunc_point = r['GAI_point']
    
    cfxrs_vari = []
    ensmemoldyields_vari = []
    ensmemnewyields_vari = []    
    #GAI_p_all = []
    #tmean_p_all = []
    #prec_p_all = []
    #solarrad_p_all = []
    #Jarray_p_all = []
    #Cday_p_all = []
    #GSS_p_all = []

    # Load the gridded data for each ensmem, then generate the ensemble,
    # mean and stdev for each obscoord
    print('Generating ensemble data')
    # gen xarrays to store data
    GAI_p_all, tmean_p_all, prec_p_all, solarrad_p_all, Jarray_p_all, Cday_p_all, GSS_p_all, AWC_allp1 = \
    create_xrs(obscoords, ensmems)
    for ensmem in ensmems:
        print('Ensemble ' + str(ensmem))
        if ensmem==ensmems[0]:
            tmean, prec, solarrad, tmax, tmin, temp_cconc, AWC = \
            load_drivingdata(loaddata, ensmem, ensmems, year, Rdatpath, CO2file)
        else:
            tmean, prec, solarrad, tmax, tmin, temp_cconc = \
            load_drivingdata(loaddata, ensmem, ensmems, year, Rdatpath, CO2file)
        
        times = []
        for month in ['10','11','12']:
            for day in np.arange(1,10):
                day = '0' + str(day)
                times.append(str(year-1)+month+day)
            for day in np.arange(10,31):
                times.append(str(year-1)+month+str(day))
        for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09']:
            for day in np.arange(1,10):
                day = '0' + str(day)
                times.append(str(year)+month+day)
            for day in np.arange(10,31):
                times.append(str(year)+month+str(day))
        times = np.array(times)
        
        # Select out grid points nearest to each obscoord, calculate GAI at each point
        print('Calculating GAI at each obspoint')
        counter=1
        totalobs = len(obscoords)
        
        for tob in obscoords:
            print('Calculating GAI for ensmem ' + str(ensmem) + ' for obs ' + str(counter) + ' of ' + str(totalobs))
            counter+=1
            tmean_p    = tmean.sel(x=tob[0], y=tob[1], method='nearest')
            x_p = float(tmean_p['x'].values)
            y_p = float(tmean_p['y'].values)
            tmean_p = tmean_p.values
            prec_p     = prec.sel(x=tob[0], y=tob[1], method='nearest').values
            solarrad_p = solarrad.sel(x=tob[0], y=tob[1], method='nearest').values
            tmax_p     = tmax.sel(x=tob[0], y=tob[1], method='nearest').values
            tmin_p     = tmin.sel(x=tob[0], y=tob[1], method='nearest').values

            proj = pyproj.Transformer.from_crs(27700, 4326, always_xy=True)
            lon,lat = proj.transform(x_p,y_p)

            datalist2 = GAIfunc_point(tmean_p, tmax_p, tmin_p, prec_p, solarrad_p, x_p, y_p, lat, times)
            
            GAI1         = np.array(datalist2.rx2('GAI')) # convert to python numpy                                          
            tmean1       = np.array(datalist2.rx2('tmean'))
            prec1        = np.array(datalist2.rx2('prec'))
            solarrad1    = np.array(datalist2.rx2('solarrad'))
            Jarray1      = np.array(datalist2.rx2('Jarray'))
            Cday1        = np.array(datalist2.rx2('Cday'))
            GSS1         = np.array(datalist2.rx2('GSS')) 
            HarvestJday = datalist2.rx2('HarvestJday')
            
            # get the data coordinates                                                                                       
            x = np.array(datalist2.rx2('x'))
            y = np.array(datalist2.rx2('y'))
            t = np.array(datalist2.rx2('t'))
            attrs = {'units': 'days since ' + t[0], 'calendar': '360_day'}

            GAI_p_all[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmem))] = GAI1
            tmean_p_all[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmem))] = tmean1
            prec_p_all[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmem))] = prec1
            solarrad_p_all[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmem))] = solarrad1
            Jarray_p_all[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmem))] = Jarray1
            Cday_p_all[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmem))] = Cday1
            GSS_p_all[str(tob[0])+','+str(tob[1])].loc[dict(ensmem=int(ensmem))] = GSS1
            if ensmem==ensmems[0]:
                AWC1      = AWC.sel(x=tob[0], y=tob[1], method='nearest').values
                AWC_allp1[str(tob[0])+','+str(tob[1])].values = AWC1

    GAI_p_all.time.attrs      = attrs
    tmean_p_all.time.attrs    = attrs
    prec_p_all.time.attrs     = attrs
    solarrad_p_all.time.attrs = attrs
    Jarray_p_all.time.attrs   = attrs
    Cday_p_all.time.attrs     = attrs
    GSS_p_all.time.attrs      = attrs
    AWC_allp1.time.attrs      = attrs
    GAI_p_all      = xr.decode_cf(GAI_p_all, decode_coords=False)
    tmean_p_all    = xr.decode_cf(tmean_p_all, decode_coords=False)
    prec_p_all     = xr.decode_cf(prec_p_all, decode_coords=False)
    solarrad_p_all = xr.decode_cf(solarrad_p_all, decode_coords=False)
    Jarray_p_all   = xr.decode_cf(Jarray_p_all, decode_coords=False)
    Cday_p_all     = xr.decode_cf(Cday_p_all, decode_coords=False)
    GSS_p_all      = xr.decode_cf(GSS_p_all, decode_coords=False)
    AWC_allp1      = xr.decode_cf(AWC_allp1, decode_coords=False)

    GAI_p_ensmean = GAI_p_all.mean(axis=0)
    GAI_p_ensstd = GAI_p_all.std(axis=0)
    
    del tmean
    del prec
    del solarrad
    del tmax
    del tmin

    return GAI_p_all, GAI_p_ensmean, GAI_p_ensstd, tmean_p_all, prec_p_all, solarrad_p_all, Jarray_p_all, Cday_p_all, GSS_p_all, HarvestJday, AWC_allp1, temp_cconc

def create_xrs(obscoords, ensmems):
    '''
    Create the xarray datasets that are filled in the ensgen routine

    Inputs: 
    obscoords - List of 2-element lists containing the x,y coordinates
                of interest on the OSGB eastings/northings grid. 
    ensmems - List of 0-padded 2-digit strings corresponding to the ensemble 
              member numbers
    
    Outputs:
    GAI_p_all, tmean_p_all, prec_p_all, solarrad_p_all, 
    Jarray_p_all, Cday_p_all, GSS_p_all, AWC_allp1:
    Xarray datasets with ensmem,t dimensions, one variable for each obscoord,
    full of zeros.
    '''
    
    times = np.arange(0, 360)
    ensmemsint = [int(ensmem) for ensmem in ensmems]
    GAI_dict = {}
    tmean_dict = {}
    prec_dict = {}
    solarrad_dict = {}
    Jarray_dict = {}
    Cday_dict = {}
    GSS_dict = {}
    AWC_dict = {}
    for tob in obscoords:
        GAI_dict[str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((len(ensmems), 360)), coords=[ensmemsint, times], dims=['ensmem', 'time'])
        tmean_dict[str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((len(ensmems), 360)), coords=[ensmemsint, times], dims=['ensmem', 'time'])
        prec_dict[str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((len(ensmems), 360)), coords=[ensmemsint, times], dims=['ensmem', 'time'])
        solarrad_dict[str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((len(ensmems), 360)), coords=[ensmemsint, times], dims=['ensmem', 'time'])
        Jarray_dict[str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((len(ensmems), 360)), coords=[ensmemsint, times], dims=['ensmem', 'time'])
        Cday_dict[str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((len(ensmems), 360)), coords=[ensmemsint, times], dims=['ensmem', 'time'])
        GSS_dict[str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((len(ensmems), 360), dtype='object'), coords=[ensmemsint, times], dims=['ensmem', 'time'])
        AWC_dict[str(tob[0])+','+str(tob[1])] = xr.DataArray(np.zeros((1, 360)).squeeze(), coords=[times], dims=['time'])
    GAI_p_all = xr.Dataset(GAI_dict)
    tmean_p_all = xr.Dataset(tmean_dict)
    prec_p_all = xr.Dataset(prec_dict)
    solarrad_p_all = xr.Dataset(solarrad_dict)
    Jarray_p_all = xr.Dataset(Jarray_dict)
    Cday_p_all = xr.Dataset(Cday_dict)
    GSS_p_all = xr.Dataset(GSS_dict)
    AWC_allp1 = xr.Dataset(AWC_dict)

    return GAI_p_all, tmean_p_all, prec_p_all, solarrad_p_all, Jarray_p_all, Cday_p_all, GSS_p_all, AWC_allp1

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
def vari_method_point(obsall, moddata, mod_ensstd, obscoords, year, ensmem, obserrtype=0, obserr=0.1, moderrinfl=1, tsvar=0.25, order=1, power=10, plot=0, plotdir='None', interval=10, nudgedtsall=None):
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
        
        # remove the zeros from either end of the modvars time series. 
        # and produce the indices needed to chop the ends of all the
        # timeseries. The assimilation routine will fail if zero stdevs
        # are used for the model or obs timeseries.
        modstds_trim, frontchop, backchop = remove_zeros(mod_ensstd[str(tob[0]) + ',' + str(tob[1])].values)
        # fill in any zeros in the middle of the timeseries
        modstds_trim = pd.Series(modstds_trim).where(pd.Series(modstds_trim)!=0).interpolate(method='linear').values
        
        # multiply the mod error by the model error inflation value (default 1)
        modstds_trim = modstds_trim * moderrinfl
        
        # select out the obs point from the obs table
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
        modGAI = moddata[str(tob[0]) + ',' + str(tob[1])].values

        # chop off the front and end of the timeseries where the modstd is 0
        modGAI = modGAI[frontchop:-backchop]

        # implement the obs error
        if obserrtype == 0:
            obsstds_trim = np.array([obserr for a in obsGAI_trim])
        elif obserrtype == 1:
            obsstds_trim = obsGAI_trim*obserr

        # generate the covariance matrices for the obs and the mod
        modcov = np.eye(len(modstds_trim))*(modstds_trim**2)
        obscov = np.eye(len(obsstds_trim))*(obsstds_trim**2)

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
        fig1=[]
        coord = str(tob[0]) + ',' + str(tob[1])
        #x,f,d = spo.fmin_l_bfgs_b(cost_function, firstguess, args=(obsGAI_trim2, modGAI_trim, obscovinv, modcovinv, smoothcovinv, yno), approx_grad=True, disp=5, maxfun=100000, pgtol=2e-02, epsilon=1e-04, factr=1e12, bounds=bounds)
        merged = spo.fmin_bfgs(cost_function, firstguess, grad_cost_function, args=(obsGAI_trim, obsGAI, modGAI, obscovinv, modcovinv, smoothcovinv, H, power, order, plot, plotdir, coord, interval, ensmem, nudgedts))
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

def update_yield_points_point(GAIold, GAInew, obscoords, tmean, prec, solarrad, Jarray, Cday, GSS, HarvestJday, AWC, temp_cconc, ensmem):
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
        GAInew1 = GAInew[str(tob[0])+','+str(tob[1])].values
        GAIold1 = GAIold[str(tob[0])+','+str(tob[1])].values
        tmean1 = tmean[str(tob[0])+','+str(tob[1])].values
        prec1 = prec[str(tob[0])+','+str(tob[1])].values
        solarrad1 = solarrad[str(tob[0])+','+str(tob[1])].values
        Jarray1 = Jarray[str(tob[0])+','+str(tob[1])].values
        Cday1 = Cday[str(tob[0])+','+str(tob[1])].values
        GSS1 = GSS[str(tob[0])+','+str(tob[1])].values
        AWC1 = AWC[str(tob[0])+','+str(tob[1])].values
        
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
        yieldlist = wheat_yield_point(GAInew1, tmean1, prec1, solarrad1, AWC1, Jarray1, Cday1, GSS1, HarvestJday, cconc=temp_cconc)
        oldyieldlist = wheat_yield_point(GAIold1, tmean1, prec1, solarrad1, AWC1, Jarray1, Cday1, GSS1, HarvestJday, cconc=temp_cconc)

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
    plt.title(title)
    if savefile:
        plt.savefig(savefile, dpi=300)
    #plt.legend([bp0["boxes"][0],bp1["boxes"][0]],['data','dscalepe'],fontsize=11,bbox_to_anchor=(0.2,-0.06),ncol=3,loc=1)
    #plt.subplots_adjust(left=0.1,right=0.98,top=0.98,bottom=0.08,wspace=0.36,hspace=0.2)
    return datameanold, datameannew


def verifyyield(yieldfile, obscoords, oldyields_vari, newyields_vari, years, savedir, allplots=0):
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
    asortyielddata = yielddata.sort_values(['area'])[::-1]
    # create pandas geodataframe for the obscoords points (which are taken from the modis data)
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

    # for each obs point...
    meanolds = []
    meannews = []
    for tob in list(actyields.index.values):
        # check if we have yield data for each year we ran for
        # and if so, plot 
        for year in years:
            if np.isnan(actyields.loc[tob][str(year)]) == False:
                actyield = actyields.loc[tob][str(year)]
                newyield = newyields_vari[tob].sel(year=year).to_pandas()
                oldyield = oldyields_vari[tob].sel(year=year).to_pandas()
                newyield.name = 'newyield'
                oldyield.name = 'oldyield'
                oldandnew = pd.merge(oldyield, newyield, on='ensmem')
                if allplots:
                    if not savedir:
                        meanold, meannew = plot_boxes_yield(oldandnew, [['oldyield','newyield']], [0,1], ['oldyield', 'newyield'], ['yield (tn/hec)'], [[7,18]], str(year) + ': ' + tob, actyield, None)
                    else:
                        savefile = os.path.join(savedir, str(year) + '_' + tob + '.png')
                        meanold, meannew = plot_boxes_yield(oldandnew, [['oldyield','newyield']], [0,1], ['oldyield', 'newyield'], ['yield (tn/hec)'], [[7,18]], str(year) + ': ' + tob, actyield, savefile)
                else:
                    meanold = np.mean(oldandnew['oldyield'])
                    meannew = np.mean(oldandnew['newyield'])
                meanolds.append(meanold)
                meannews.append(meannew)

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
    oldsandnews = pd.merge(meanoldspd, meannewspd, on='point')
    yieldsf = gpd.read_file(yieldfile)
    yieldsf.drop(['field_id', 'area', 'geometry'], axis=1, inplace=True)
    allyields = yieldsf.values.flatten()
    yields = allyields[~np.isnan(allyields)]
    #actyields.drop(['x', 'y', 'geometry', 'index_right', 'field_id', 'area'], axis=1, inplace=True)
    #actyields.index = oldsandnews.index
    #oldsnewsacts = pd.merge(oldsandnews, actyields, on='point')
    #actyield = actyields.mean(axis=0).mean()
    #print(actyield)

    savefile = os.path.join(savedir, 'mean_yield.png')
    junk1, junk2 = plot_boxes_yield(oldsandnews, [['mod_yield','assim_yield']], [0,1,2], ['mod_yield', 'assim_yield', 'actual_yield'], ['yield (tn/hec)'], [[7,18]], 'Ensemble-mean yields', yields, savefile)
