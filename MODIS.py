import os
import time
import glob
import pyproj
import requests
import numpy as np
import xarray as xr
import pandas as pd
import cftime as cft
import datetime as dt
import urllib.request
import geopandas as gpd

from utils import *

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
        oscswitch=0
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
        obsyields = gpd.read_file(infile)
        
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
    try:
        obscoordslonlat = [proj.transform(x,y, errcheck=True) for x,y in obscoordsxy]
    except pyproj.exceptions.ProjError:
        obscoordslonlat = [proj.transform(x,y, errcheck=True) for x,y in obscoordsxy]

    # Submit request to MODIS servers for data for each location
    datacodes=[]
    emails=['a%40b.c', 'd%40e.f', 'g%40h.i', 'j%40k.l', 'm%40n.o', 'p%40q.r', 's%40t.u', 'v%40w.x', 'y%40z.z', 'test%40test.c', 'test%40test.b', 'test%40test.a', 'fun%40be.com', 'matt%40me.com', 'matt%40me.be', 'mat%40think.c', 'bla%40blabla.com', 'bladebla%40evenmorebla.com', 'testy%40me.com', 'cdb%40de.com', 'mccool%40as.ice', 'cold%40hot.water', 'check%40mate.mb', 'need%40idea.s']

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
        if 'error' in datacode:
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

def MODIS_process(datasavedir, caltype, filter_threshold=0.5):
    '''
    Process the MODIS data from the downloaded CSV files into
    a pandas dataframe usable by the code.
    The data are translated onto a 360day calendar by shifting
    anything on the 31st of a month onto the 30th so that no 
    data is lost, if the driving data is 360day. Otherwise they
    are kept on a gregorian calendar. 
    Values below a default threshold of 0.5 are
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
        #print(lon,lat)

        # Convert to OSGB
        proj = pyproj.Transformer.from_crs(4326, 27700, always_xy=True)
        try:
            x,y = proj.transform(lon,lat, errcheck=True)
        except pyproj.exceptions.ProjError:
            x,y = proj.transform(lon,lat, errcheck=True)
        #print(x,y)

        # Set name of column for data table
        colname = str(x).split('.')[0] + ',' + str(y).split('.')[0]
        coord = [int(str(x).split('.')[0]), int(str(y).split('.')[0])]
        obscoords.append(coord)
        print(coord)

        # Read in data file 
        pixel = pd.read_csv(filein, header=None, index_col=0, usecols=[2,6], names=['date', colname], na_values='F')
        # Convert date index column to datetimes
        pixel.index = pd.to_datetime(pixel.index, format='%Y%j', exact=False)

        # if calendar is 360day, then we need to do some nasty date wrangling to get the 
        # MODIS data on this calendar
        if caltype == '360_day':
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
        else:
            startdate = str(pixel.index[0])[:10]
            enddate   = str(pixel.index[-1])[:10]
            #oldindex = pixel.index.values
            #newindex1 = [dt.datetime(int(str(d)[:4]), int(str(d)[5:7]), int(str(d)[8:10])) for d in list(oldindex)]
            #pixel.index = pd.DatetimeIndex(newindex1)
            alldaysidx = xr.cftime_range(startdate, enddate, calendar='standard', freq='D', name='date')
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
        seldates = xr.cftime_range(str(year)+'-03-01', str(year)+'-07-30', calendar=caltype, freq='D', name='date').values
        obspdall.loc[seldates] = obspdall.loc[seldates].where(obspdall.loc[seldates]>thresh)

    return obspdall, obscoords


def get_S2tilebounds(S2dir):
    tiles = ['30UYE', '30UYD', '30UYC', '30UYB',
             '30UXE', '30UXD', '30UXC', '30UXB',
             '30UWE', '30UWD', '30UWC', '30UWB']
    
    tilebounds = {}
    for tile in tiles:
        testtile = glob.glob(os.path.join(S2dir, '*' + tile + '*'))[0]
        testtilexr = xr.open_rasterio(testtile)
        
        xcoords = testtilexr['x'].values
        ycoords = testtilexr['y'].values
        xcoords = np.sort(xcoords)
        ycoords = np.sort(ycoords)
        xres = xcoords[1]-xcoords[0]
        yres = ycoords[1]-ycoords[0]
    
        xbounds = [xcoords[0]-xres/2, xcoords[-1]+xres/2]
        ybounds = [ycoords[0]-yres/2, ycoords[-1]+yres/2]
        tilebounds[tile] = [xbounds, ybounds]
    return tilebounds


def process_S2(S2dir, yieldfile=None, startdate="2015-01-01", enddate="2019-12-31", cloudthresh=40, fieldno=None, calendar='gregorian'):

    print('Getting S2 tile bounds')
    tilebounds = get_S2tilebounds(S2dir)
    tiles = ['30UYE', '30UYD', '30UYC', '30UYB',
             '30UXE', '30UXD', '30UXC', '30UXB',
             '30UWE', '30UWD', '30UWC', '30UWB']

    if yieldfile:
        print('Reading in yields')
        # Read in shapefile containing yield data into geopandas dataframe
        obsyields = gpd.read_file(yieldfile)

        # Extract out the geometries and find the 'center' coordinates of each field
        if fieldno:
            fieldIDs = [list(obsyields['field_id'])[fieldno-1]]
            obscoordsxy = [[point.coords[0] for point in list(obsyields.centroid.values)][fieldno-1]]
        else:
            fieldIDs = list(obsyields['field_id'])
            obscoordsxy = [point.coords[0] for point in list(obsyields.centroid.values)]
        print(fieldIDs)
        print(obscoordsxy)
        
        # Convert these to Sentinel2 coords
        UTM30U = pyproj.CRS("+proj=utm +zone=30U, +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
        OSGB = pyproj.CRS(init='epsg:27700')
        OSGBtoUTM30U = pyproj.Transformer.from_crs(OSGB, UTM30U)
        print('Set up coord transform')
        # this try/except block is to get around a weird pyproj error/bug where it fails
        # for the first transformation but not subsequent ones in a code block
        try:
            obscoordsUTM30 = [OSGBtoUTM30U.transform(x,y, errcheck=True) for x,y in obscoordsxy]
            print('Trying coordinate transform')
        except pyproj.exceptions.ProjError:
            obscoordsUTM30 = [OSGBtoUTM30U.transform(x,y, errcheck=True) for x,y in obscoordsxy]
            print('Trying coordinate transform again')
        print(obscoordsUTM30)

        # repeat the top row at the bottom to get around the same pyproj error
        # for some reason this means we also have to reset the crs info on the dataframe
        obsyieldstemp = obsyields.append(obsyields.iloc[0,:])
        temp = obsyields.crs
        obsyieldstemp.crs = temp

        obsyieldsUTM30 = obsyieldstemp.to_crs(epsg=32630).iloc[1:,:]
        obsyieldsUTM30 = obsyieldsUTM30.iloc[[-1],:].append(obsyieldsUTM30.iloc[:-1,:])

        sfUTM30name = yieldfile[:-4] + 'UTM30.shp'
        if not os.path.exists(sfUTM30name):
            obsyieldsUTM30.to_file(sfUTM30name)

        print('Creating pandas table to store LAI data in')
        # create the pandas table to store the data in
        # index (rows)
        # generate list of dates
        tempindex = list(pd.date_range(startdate, enddate))
        if calendar=='360day':
            # change any 31sts to 30ths, to match 360day calendar of UKCP18 data
            tempindex2 = [d if d.day<=30 else dt.datetime(d.year, d.month, 30) for d in tempindex]
            tempindex3 = list(pd.DatetimeIndex(tempindex2))
            # remove duplicate 30ths
            tempindex4 = pd.DatetimeIndex(list(np.unique(np.asarray(tempindex3))))
            # convert to 360day-calendar datetime index
            cfdatetimes = [cft.Datetime360Day(d.year, d.month, d.day) for d in tempindex4]
            cfdatetimesidx = xr.coding.cftimeindex.CFTimeIndex(cfdatetimes)
        else:
            tempindex3 = list(pd.DatetimeIndex(tempindex))
            cfdatetimes = [cft.DatetimeGregorian(d.year, d.month, d.day) for d in tempindex3]
            cfdatetimesidx = xr.coding.cftimeindex.CFTimeIndex(cfdatetimes)
        # column names (the OSGB coords, minus all the digits after the decimal points)
        colnames = [str(x).split('.')[0] + ',' + str(y).split('.')[0] for x,y in obscoordsxy]
        # dataframe
        LAIall = pd.DataFrame(index=cfdatetimesidx, columns=colnames)

        totalfields = len(fieldIDs)
        for fID in range(0, totalfields):
            print('Processing field ' + str(fID+1) + ' of ' + str(totalfields))
            tfccs2 = obscoordsUTM30[fID]
            print(tfccs2)
            
            # find what S2 tile this field is in
            for tile in tiles:
                if tilebounds[tile][0][1] > tfccs2[0] > tilebounds[tile][0][0]:
                    if tilebounds[tile][1][1] > tfccs2[1] > tilebounds[tile][1][0]:
                        fieldtile = tile
                        print('Field in tile ' + fieldtile)
                        break

            # read in the tile (for each timestep) and subset to field
            fieldtiles = glob.glob(os.path.join(S2dir, '*' + fieldtile + '*'))
            for ftile in fieldtiles:
                print('Processing ' + os.path.basename(ftile) + ' for field ' + str(fID+1) + ' of ' + str(totalfields))
                data = xr.open_rasterio(ftile)
                test = country_subset_shapefile(data=data, sfname=sfUTM30name, IDname='field_id', IDs=[fieldIDs[fID]])
                maskedtest = test/1000.
                maskedtest = xr.where(maskedtest < 10, maskedtest, np.nan)
                maskedtest = xr.where(maskedtest > 0, maskedtest, np.nan)

                # calculate % of field that is not covered by cloud
                # then only use value if % is > 40% or something
                nfieldpoints = xr.where(np.isnan(test), np.nan, 1).sum().values
                nunmasked = xr.where(np.isnan(maskedtest), np.nan, 1).sum().values
                percentunmasked = (nunmasked/nfieldpoints)*100
                print(str(percentunmasked) + '% of field is not covered by cloud')

                if percentunmasked > cloudthresh:
                    fieldavgLAI = maskedtest.mean().values
                else:
                    fieldavgLAI = np.nan
                
                obsdate = os.path.basename(ftile).split('_')[1][:8]
                year = obsdate[:4]
                month = obsdate[4:6]
                day = obsdate[6:8]
                if calendar=='360day':
                    if int(day) > 30:
                        obsdate = year+month+'30'
                
                LAIall.loc[obsdate,colnames[fID]] = fieldavgLAI
                LAIall.to_csv('field_' + str(fieldno) + '.csv')
    return LAIall


def S2_readin(filein, caltype):
    '''
    Read in csv file of all fields to run for into pandas table 
    like the one produced in the above function
    '''

    allfields = pd.read_csv(filein, index_col=0, parse_dates=True)
    coordstemp = list(allfields.columns.values)
    coords = [[c.split(',')[0], c.split(',')[1]] for c in coordstemp]

    # create index with all the days of the calendar used (including Feb 29th, 30th if 360day)
    d = allfields.index[0]
    e = allfields.index[-1]
    if caltype == '360day':
        startdate = cft.Datetime360Day(d.year, d.month, d.day)
        enddate = cft.Datetime360Day(e.year, e.month, e.day)
        alldaysidx = xr.cftime_range(startdate, enddate, calendar='360_day', freq='D', name='date')
    elif caltype == 'gregorian':
        startdate = cft.DatetimeGregorian(d.year, d.month, d.day)
        enddate = cft.DatetimeGregorian(e.year, e.month, e.day)
        alldaysidx = xr.cftime_range(startdate, enddate, calendar='gregorian', freq='D', name='date')

    if caltype == '360day':
        # convert original index of data to 360day calendar datetimes
        tempindex = list(pd.date_range(allfields.index[0], allfields.index[-1]))
        # change any 31sts to 30ths, to match 360day calendar of UKCP18 data
        tempindex2 = [d if d.day<=30 else dt.datetime(d.year, d.month, 30) for d in tempindex]
        tempindex3 = list(pd.DatetimeIndex(tempindex2))
        # remove duplicate 30ths
        tempindex4 = pd.DatetimeIndex(list(np.unique(np.asarray(tempindex3))))
        # convert to 360day-calendar datetime index
        cfdatetimes = [cft.Datetime360Day(d.year, d.month, d.day) for d in tempindex4]
        cfdatetimesidx = xr.coding.cftimeindex.CFTimeIndex(cfdatetimes)
        allfields.index = cfdatetimesidx
    elif caltype == 'gregorian':
        # convert original index of data to CFtime gregorian calendar datetimes
        tempindex = list(pd.date_range(allfields.index[0], allfields.index[-1]))
        tempindex3 = list(pd.DatetimeIndex(tempindex))
        cfdatetimes = [cft.DatetimeGregorian(d.year, d.month, d.day) for d in tempindex3]
        cfdatetimesidx = xr.coding.cftimeindex.CFTimeIndex(cfdatetimes)
        allfields.index = cfdatetimesidx

    allfields = allfields.reindex(alldaysidx)    
                                            
    return allfields, coords
