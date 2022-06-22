# Simple script to test the imports of all the packages we need
# to run the cropnet model
# MJB 21/06/22

print('Importing packages')

# import local packages
from utils import *
from MODIS import *

# import environment (anaconda) packages
# R
import rpy2
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector
from rpy2.robjects import numpy2ri
numpy2ri.activate()

# python
import os
import sys
import time
import mgrs
import glob
import shutil
import pyproj
import requests
import rasterio
import numpy as np
import xarray as xr
import pandas as pd
import cartopy as cp
import cftime as cft
import urllib.request
import netCDF4 as nc4
import datetime as dt
import geopandas as gpd
import scipy.optimize as spo
import matplotlib.pyplot as plt

from cdo import *
from affine import Affine
from rasterio import features
from dateutil.relativedelta import relativedelta
from dask.diagnostics import Profiler, ResourceProfiler
from dask.diagnostics import CacheProfiler, ProgressBar, visualize

#import warnings
#warnings.filterwarnings('ignore')
