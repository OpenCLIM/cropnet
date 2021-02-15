import geopandas as gpd
import pandas as pd
import pyproj
import requests
import time
import os
import urllib.request
import glob
import datetime as dt
import xarray as xr
import cftime as cft

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

        # Convert to OSGB
        proj = pyproj.Transformer.from_crs(27700, 4326, always_xy=True)
        x,y = proj.transform(lon,lat, errcheck=True)

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
