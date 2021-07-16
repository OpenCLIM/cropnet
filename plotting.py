import numpy as np
import xarray as xr
import cartopy as cp
import geopandas as gpd
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.io.shapereader as shapereader
import matplotlib.animation as animation
import matplotlib.dates as mdates
import os
from dask.diagnostics import Profiler, ResourceProfiler, CacheProfiler, ProgressBar, visualize
import pyproj
import geopandas as gpd
from rasterio import features
from affine import Affine
from shapely import geometry

def plot_yearanim(data, lonname='longitude', latname='latitude', timname='time', startdate='1980-01-01', plotlab=None, title='', outname='test/test', sfile=None, vmin='None', vmax='None', colmap='Reds', ext='gif', climo=0):
    '''
    Plots a timeseries animation of daily data from an xarray
    Currently only plots the first year of the data

    Inputs:
    data - The xarray dataarray to plot
    lonname - Name of the x dimension
    latname - Name of the y dimension
    timname - Name of the t dimension
    startdate - Start date of data (and the animation)
    plotlab - Label for the colourbar
    title - Plot title
    outname - Save path of animation
    sfile - Shape file to plot in addition to the data (e.g. rivers)
    vmin - Lower bound of colourbar
    vmax - Upper bound of colourbar
    colmap - Colormap to use
    ext - gif or mp4
    climo - Controls how the title looks. 1 if plotting a climatology, 0 otherwise
    '''

    # Create output directory
    dirname = os.path.split(outname)[0]
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    
    # Transpose data so that tim dimension is last
    data = data.transpose(latname, lonname, timname)
    lons = data[lonname].values
    lats = data[latname].values

    if len(data[timname]) > 365:
        endstep = 366
    else:
        endstep = 365

    # read in shapefile data
    if sfile:
        rivers = cp.feature.ShapelyFeature(cp.io.shapereader.Reader(sfile).geometries(), cp.crs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.25)

    # set up first frame of animation
    fig = plt.figure()
    ax1 = plt.axes(projection = cp.crs.PlateCarree())
    ax1.coastlines(resolution='10m') # plot coastlines
    # add shapefile 
    if sfile:
        ax1.add_feature(rivers)
    # add country borders
    ax1.add_feature(cp.feature.BORDERS, linestyle='-', alpha=.5)
    # set extent of plot to match the data
    ax1.set_extent((lons[0], lons[-1], lats[0], lats[-1]))
    # set axes labels
    gl = ax1.gridlines(draw_labels=True, alpha=0)
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    # format date for title
    startdate = dt.datetime.strptime(startdate, '%Y-%m-%d')
    if climo == 1:
        ax1.set_title(title + ' ' + startdate.strftime('%b %d') + ' climo')
    else:
        ax1.set_title(title + ' ' + startdate.strftime('%b %d %Y'))
    # set colourbar limits
    if vmax is not 'None':
        maxval=vmax
    else:
        maxval = data.max()
    if vmin is not 'None':
        minval=vmin
    else:
        minval = data.min()
    # plot
    pcm = ax1.pcolormesh(lons, lats, data[:,:,0].values, vmin=minval, vmax=maxval, cmap=colmap)
    # manually add a colourbar
    cbaxes=fig.add_axes([0.05, 0.1, 0.05, 0.8])
    plt.colorbar(pcm, cax=cbaxes, extend='both')
    cbaxes.yaxis.set_label_position("left")
    cbaxes.yaxis.tick_left()
    # add label to colourbar
    if not plotlab:
        try:
            plotlab = data.name
        except AttributeError:
            plotlab = ''
    cbaxes.set_ylabel(plotlab)
    
    # Function to run to plot for all the other frames of the animation.
    # f is just the frame number
    def animate(f):
        print('Plotting frame ' + str(f+1))
        # clear the current frame
        ax1.clear()
        # then repeat the above plotting code for each frame
        ax1.coastlines(resolution='10m')
        if sfile:
            ax1.add_feature(rivers)
        ax1.add_feature(cp.feature.BORDERS, linestyle='-', alpha=.5)
        ax1.set_extent((lons[0], lons[-1], lats[0], lats[-1]))
        gl = ax1.gridlines(draw_labels=True, alpha=0)
        gl.xlabels_top = False
        gl.ylabels_left = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        if climo == 1:
            ax1.set_title(title + ' ' + (startdate + dt.timedelta(f)).strftime('%b %d') + ' climo')
        else:
            ax1.set_title(title + ' ' + (startdate + dt.timedelta(f)).strftime('%b %d %Y'))
        pcm = ax1.pcolormesh(lons, lats, data[:,:,f].values, vmin=minval, vmax=maxval, cmap=colmap)
        #pcm.cmap.set_under('white')
        return pcm,

    # call animation function
    anim = animation.FuncAnimation(fig, animate, range(0,endstep), blit=True)
    
    # save to disk
    if ext=='mp4':
        anim.save(outname + '.mp4', fps=10, dpi=300, bitrate=10000, writer='ffmpeg')
    elif ext=='gif':
        anim.save(outname + '.gif', writer='imagemagick', fps=10, dpi=300)
    else:
        raise Error('format ' + ext + 'not supported')

    print('Animation saved to ' + outname + '.' + ext)

    plt.close()
    return


def plot_points(datalist, var, coords, xname, yname, tname, outdir, title='None', labels='None', year='None', startdate=None):
    '''
    Plot timeseries of datasets at various coords, or at the same coordinate over various years
    For the latter, set the start date to '10-01' (assuming your time series start 1st Oct)
    
    Inputs:
    datalist - List of the filepaths of the netcdf datasets to plot
               Must all be on the same grid
    var      - Variable name from the netcdf datasets
               to plot. 
    coords - List of x,y pairs, coords to plot
    xname - Name of x dim
    yname - Name of y dim 
    tname - Name of t dim in the datatsets provided
    outdir - Location to save plots
    title - Title of plot
    labels - Labels for legend
    year - Only used for the output filename
    startdate - When to start the xaxis
              - Note the time axis of the data will be changed to fit this
              - This is to allow plotting data on different time axes on top of each other

    Outputs:
    mean and stdev of plotted coords/years are also returned
    '''
    
    # read in data
    datasets = [xr.load_dataset(ds) for ds in datalist]

    datasets = [dataset.rename({xname: 'x', yname: 'y'}) for dataset in datasets]

    # to plot datasets with different time axes on top of each other,
    # set the time variable of all the datasets to be the same,
    # starting at the specified month and day in startdate
    if startdate:
        newleaptc = pd.date_range('2003-' + startdate, freq='D', periods=366)
        newtc = pd.date_range('2003-' + startdate, freq='D', periods=366).drop(np.datetime64('2004-02-29'))

        for ds in range(0, len(datasets)):
            if len(datasets[ds][tname]) > 365:
                datasets[ds][tname] = newleaptc
            else:
                datasets[ds][tname] = newtc

    
    dataarrs = [dataset[var] for dataset in datasets]
    tlen = len(dataarrs[0].sel(x=coords[0][0], y=coords[0][1], method='nearest'))
    datastore = np.zeros((tlen, len(datasets)*len(coords)))
    counter = 0
    ax = plt.axes()
    for dataarr in dataarrs:
        for coord in coords:
            data = dataarr.sel(x=coord[0], y=coord[1], method='nearest')
            try:
                data.plot.line(x=tname, ax=ax)
            except AttributeError:
                data.plot(ax=ax)
            if startdate:
                if len(data) > 365:
                    data = data.drop_sel({tname: np.datetime64('2004-02-29')})
            datastore[:,counter] = data.values
            counter += 1 
    mean = datastore.mean(axis=1)
    stdev = datastore.std(axis=1)
            
    if tname is not 'dayofyear':
        months = mdates.MonthLocator()  # every month
        months_fmt = mdates.DateFormatter('%b')
        ax.xaxis.set_major_locator(months)
        ax.xaxis.set_major_formatter(months_fmt)
    if title is 'None':
        ax.set_title(dataarr.name)
    else:
        ax.set_title(title)
    if labels is 'None':
        ax.legend(coords)
    else:
        ax.legend(labels)
    if year is 'None':
        try:
            year = dataarr[tname].dt.year.values[-1]
        except TypeError:
            year = ''
    savename = dataarr.name + '_' + str(year) + '.png'

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outname = os.path.join(outdir, savename)
    plt.savefig(outname, dpi=300)
    return mean,stdev


def spatial_plot(data, title, savename='None', countries=['China'], maskval='None'):
    sfname = shapereader.natural_earth(resolution='10m', category='cultural', name='admin_0_countries')
    maskeddata = country_subset_shapefile(data=data, sfname=sfname, IDname='ADMIN', IDs=countries)
    ax = plt.axes(projection=cp.crs.PlateCarree())
    borders = cp.feature.NaturalEarthFeature('cultural', 'admin_0_countries', '10m', edgecolor='black', facecolor='none')
    ax.add_feature(borders)
    if maskval is not 'None':
        maskeddata = maskeddata.where(maskeddata > maskval)
    maskeddata.plot()
    gl = ax.gridlines(crs=cp.crs.PlateCarree(), draw_labels=True, alpha=0)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = cp.mpl.gridliner.LONGITUDE_FORMATTER
    gl.yformatter = cp.mpl.gridliner.LATITUDE_FORMATTER
    plt.title(title)
    if savename is not 'None':
        plt.savefig(savename, dpi=300)


def plot6_summary(files, year, titles, lonnames, latnames, tnames, startdate, minval, maxval, colmaps, outdir, boxes=[]):

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # read in data
    dataarrays = [xr.load_dataarray(file) for file in files]
    sfname = shapereader.natural_earth(resolution='10m',
                                       category='cultural',
                                       name='admin_0_countries')
    # subset temp and precip to China
    dataarrays[4] = country_subset_shapefile(data=dataarrays[4], sfname=sfname, xname=lonnames[4], yname=latnames[4],
                                                                             IDname='ADMIN', IDs=['China'], drop=1)
    dataarrays[5] = country_subset_shapefile(data=dataarrays[5], sfname=sfname, xname=lonnames[5], yname=latnames[5],
                                                                             IDname='ADMIN', IDs=['China'], drop=1)

    # Ensure all data has the same dimension order (shape)
    for d in range(0, len(dataarrays)):
        dataarrays[d] = dataarrays[d].transpose(latnames[d], lonnames[d], tnames[d])

    # Read in lons and lats
    lons = []
    lats = []
    for d in range(0, len(dataarrays)):
        lons.append(dataarrays[d][lonnames[d]])
        lats.append(dataarrays[d][latnames[d]])

    if len(dataarrays[0][tnames[0]]) > 365:
        endstep = 366
    else:
        endstep = 365
    
    # Set up plot
    plt.rcParams['figure.figsize'] = [14, 12]
    fig, axes = plt.subplots(3,2, subplot_kw=dict(projection = cp.crs.PlateCarree()))
    axes = axes.flatten()

    # Function for plotting one of the six plots
    def plot_single(a, t):
        # add coastlines and borders and set plot extent to match data
        axes[a].coastlines(resolution='10m')
        axes[a].add_feature(cp.feature.BORDERS, linestyle='-', alpha=.5)
        axes[a].set_extent((lons[a][0], lons[a][-1], lats[a][0], lats[a][-1]))
        # Set coord labels on the left or right depending on which plot we're doing
        gl = axes[a].gridlines(draw_labels=True, alpha=0)
        gl.xlabels_top = False
        if a%2 == 0:
            gl.ylabels_left = False
        else:
            gl.ylabels_right = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        # Format the date for the title (except for the yield plot)
        startdatet = dt.datetime.strptime(startdate, '%Y-%m-%d')+dt.timedelta(t)
        if a != 3:
            axes[a].set_title(titles[a] + ' ' + startdatet.strftime('%b %d %Y'))
        else:
            axes[a].set_title(titles[a])

        # plot data
        if a != 3:
            pcm = axes[a].pcolormesh(lons[a], lats[a], dataarrays[a][:,:,t].values, vmin=minval, vmax=maxval, cmap=colmaps[a])
        else:
            pcm = axes[a].pcolormesh(lons[a], lats[a], dataarrays[a][:,:,0].values, vmin=minval, vmax=maxval, cmap=colmaps[a])

        # add boxes
        if len(boxes)>0:
            geoms = [geometry.box(minx=box[0], maxx=box[1], miny=box[2], maxy=box[3]) for box in boxes]
            axes[a].add_geometries(geoms, crs=cp.crs.PlateCarree(), facecolor='None', edgecolor='k')
            
        return pcm
                   
    # Call the function for each plot, for the first frame of the animation, and add colourbars on alternating sides
    for a in range(0, len(axes)):
        pcm = plot_single(a, 0)
        
        if a%2 == 0:
            cb = plt.colorbar(pcm, ax=[axes[a]], location='left', pad=0.05, shrink=0.93)
        else:
            cb = plt.colorbar(pcm, ax=[axes[a]], location='right', pad=0.05, shrink=0.93)

    # Function to plot all the other frames
    def animate(t):
        print('Plotting frame ' + str(t+1))
        pcms = []
        for a in range(0, len(axes)):
            axes[a].clear()
            pcm = plot_single(a, t)
            pcms.append(pcm)
        return pcms

    # Call animation function
    anim = animation.FuncAnimation(fig, animate, range(0, endstep), blit=True)

    # Save animation
    outname = os.path.join(outdir, '6plot_summary_' + str(year) + '.mp4')
    anim.save(outname, fps=5, dpi=300, bitrate=10000, writer='ffmpeg')

    plt.close()
    return
                            

    
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
            data = xr.load_dataset(filein)

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
