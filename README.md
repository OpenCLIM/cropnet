CropNET Data Assimilation tool for UKCP18 wheat yield projections
-----------------------------------------------------------------

Overview
--------
This tool produces estimates of crop yield using a wheat crop model and MODIS or Sentinel2 LAI (leaf-area-index) satellite data,
for each pixel satellite data is provided for. The final product is NOT gridded. This is under development. 
The crop model produces estimates of potential yield, which are usually higher than the actual, measured yields.
The crop model is therefore combined with satellite derived measures of crop growth to produce more accurate values.
For a year of the UKCP18 driving (meteorological) data that we also have satellite data for:
- The satellite data is requested and downloaded if necessary (MODIS)
- The raw satellite data is downloaded and processed to produce LAI estimates on a 20mx20m grid (Sentinel2)*
- The grids are processed to produce field-scale estimates of LAI where there are enough non-cloud pixels (Sentinel2)*
- The modelled GAI over the year is calculated for each ensemble member
- The ensemble variance is used to generate a baseline estimate of the uncertainty in the modelled GAI (optional)
- A variational data assimilation method is applied to the timeseries of the satellite LAI at each pixel (MODIS) or field (Sentinel2) downloaded
  and the model GAI for the corresponding grid point within which it resides. This produces an optimum, smoothed
  timeseries that is the best guess of the true GAI, given the errors associated with the model and the observed values.
  See below for more details on the method.
- The optimal GAI timeseries is fed back in to the crop model and used to calculate a new yield value

* These steps are not automatic and have to be run separately, ideally on a cluster such as JASMIN LOTUS. 

TO BE IMPLEMENTED:
Once this is complete, the average correction factor vector is calculated over each year and ensemble member, for each
observation pixel used.
Using the data assimilation for future years, using a simple linear regression approach, or machine learning to tailor the
assimilation to each field. 

Details of the winter wheat crop model used can be found in Redhead et al. (2020): The influence of landscape composition and configuration on crop yield resilience, Journal of Applied Ecology, 57, 2180-2190, https://doi.org/10.1111/1365-2664.13722

Data Assimilation details
-------------------------
A method a bit like 4Dvar that attempts to minimize a cost function which consists of terms that grow in magnitude the
further the assimilated/guessed timeseries is away from a) the modelled GAI, b) the observations on the days that they
exist and c) a smoothed version of the timeseries itself, to find the optimum timeseries given these constraints.

The strength of each constraint is dependent on the errors (variances) associated with each timeseries. The smaller the
error, the more strongly the optimal timeseries is pulled towards the timeseries. These are user configurable.
By default, the error in the modelled timeseries is assumed to be the variance of the ensemble. There is a multiplier (default 1)
that can be changed to account for errors in the model not accounted for by errors in the meteorological driving data. 
The error in the observed timeseries/points is assumed to be 0.4 regardless of the magnitude of the observation. This
absolute error can be switched to a percentage error, whereby the error is assumed to be, e.g. 10% of the obs value.
However, in testing it was found that the absolute error produces much smoother results, particularly early on in the
timeseries when the observation values are small.

The error in the smoothed timeseries is a little more complex, but can be understood as the value by which the timeseries
should not change by more than for one timestep to the next. This is set to vary depending on the magnitude of each value in the
modelled GAI. If the values at a part of the modelled GAI timeseries are small, then this error is smaller, and thus the 
smoothing constraint is stronger for this part of the timeseries, and vice versa.

Other controls on the strength of the smoothing constraint exist. These are power and order.
Power controls the overall strength of the constraint throughout the timeseries, so the higher this is, the smoother
the optimal timeseries will be.
Order controls how it does the smoothing. An order of 1 will only consider 1st order differences, 2, 2nd order differences
etc. For example, order 2 will try to make the derivative of the optimal timeseries smooth, which can result in the
actual timeseries being very variable!
The defaults of power 10 and order 1 have been set through trial and error to provide the best results.
They are user configurable, but change them at your peril!!

Data requirements
-----------------
The code has been designed to work with UKCP18 data, but in theory any ensemble meteorological driving data can be used, 
provided it is in a similar format to the UKCP18 data. Note though, that this has not been tested.
The UKCP18 data that the model has been tested with is available on the CEDA archive at 
https://catalogue.ceda.ac.uk/uuid/589211abeb844070a95d061c8cc7f604
The format requirements are:
 - The daily variables required are total precipitation in mm (pr), maximum 2m temperature in degC (tasmax),
minimum 2m temperature in degC (tasmin), average net surface shortwave radiation in W/m^2 (rss), average 2m temperature in degC 
(tas).
 - The variable names are user editable using the variable 'varnames', but the varnames specified must be both in the filenames 
and the name of the variables within the netcdf files. 
 - The end of each filename must also be written _nn_day_startdate_enddate.nc with startdate and enddate formatted yyyymmdd. 
nn is the ensemble number, but can be any number string, so long as these are input to the code using the 'ensmems' variable.
- The data must be on an x,y grid where x and y are OSGB eastings and northings. The variable names of the x and y
variables must be projection_x_coordinate and projection_y_coordinate. 
- The data must be on a 360day calendar.

Available Water Content (AWC) is also required by the model. This data is provided as a geotiff in the 'MaxWet1.tif' and 
'MaxWet1.tfw' files. The data was derived from the UKCEH Grid-to-Grid model, taking the difference between the max and min
soil moisture content at each (1km) grid point over the UK from a multi-year simulation. 

Annual CO2 concentration is an optional requirement. This is provided in the 'UKCP18_CO2_RCP85.csv' file.

If you wish to verify the yield produced by the assimilation, you will need to set the verify switch to 1 and 
supply yield data in the form of a shapefile containing the polygons of each field with coordinates in OSGB eastings 
and northings. We cannot supply this data due to confidentiality arrangements with the data suppliers. 


Running the code
----------------
The code is written in python, which calls the crop model written in R.
The code uses several python/R packages that are not installed by default on the UKCEH Linux systems, or JASMIN.
The easiest way to install these is to use anaconda - a package management system for python and R modules/packages.
Anaconda installs packages in your home directory, and can be managed separately to any other python/R installations
you have, i.e. it won't interfere with these, and can be enabled only when you want to run this programme.
See below for instructions on installing and using anaconda.

MODIS
-----
If using MODIS satellite data, then everything can be run from the wrapper script or notebook. There are extra steps if using S2.
The code can either be run from a python script or a jupyter notebook. If you have not heard of the latter, then it
is probably easier to use the former.
To run from the script:
- Open wrapper.py in a text editor and change any of the user-editable values at the top that you need to
- Load the anaconda environment as below
- Load ipython (or just python, but ipython looks nicer!) from the directory with the script in
- Then type 'run wrapper.py'

To run from the notebook:
- Load the anaconda environment as below
- Run 'jupyter notebook assimilation.ipynb' in the directory containing this ipynb file. This will open up a web browser
  with the notebook loaded.
- Run the cell with all the import statements in
- Change the user-editable variables in the next code cell as necessary and run it
- Run the main code in the cell below

The code is set up to run for the x,y coordinates specified by the obscoords variable, for a single year or multiple years,
for each UKCP18 ensemble member. To download the required MODIS LAI data for assimilation set obs=1.
After this, if the xy coords and years have not been changed, the obs switch can be set to 2, so that it uses what is already
downloaded. 

Sentinel2
---------
If using Sentinel2 satellite data, there are separate steps that need to be taken before running the main wrapper:
1 - Downloading of the raw satellite data and calculating the LAI from this using the Sensor Invariate Atmospheric Correct (SIAC)
    and S2_TOA_TO_LAI python packages which are included in this repo
2 - Subsetting the produced LAI tiles to obtain field-scale LAI values

If running for more than a handful of fields, step 1 will need to be run in parallel on a cluster, otherwise it will take 
many days. Example scripts have been provided for the JASMIN LOTUS cluster - S2LAIwrapper.py, the python script that does
the processing and sbatch_multiplot_template which submits the python script to the cluster and contains the job-control settings.
The S2LAIwrapper.py script is set up to take values for the S2 UTM grid code, year and month as inputs and download all matching 
tiles. The S2 UTM grid code is the code given to the specific projected-coordinate grid of each of the 100x100km tiles that all 
Sentinel2 data is split into. The easiest way to find the one you want is to go to:
https://s2.boku.eodc.eu/ and click on 'pick from the map' then click where you want to download, then 'Find'. In the list of files
you can see the tile name as one of the columns. It may be that there are several tiles that cover the point you have clicked on,
you can see the area spanned by each of them by clicking on the row containing the tile of interest. It is also better to narrow
the date range a little before clicking 'Find', otherwise it will return tiles from the entire catalogue, which can get quite slow
Sometimes the 100x100km tiles are truncated diagonally for some reason, so you may have to click through a few to find a square
one that is actually representative of the area usually encompassed by that tile/UTM grid code. 
Once you have the UTM grid code and times of the tiles you want to download and process, generate several
'sbatch_multiplot_template' files, each with inputs spanning one month and UTM grid code, and submit these all to a cluster.



The final output of the code will be the original, modelled yield and the updated, assimilated yield in xarray datasets,
which are similar in structure to netcdf files. The data will also be saved on disk as netcdf files. The datasets will
have one variable per obs pixel, with the name of the variable the x,y coordinates in OSGB eastings/northings. Each
variable will have ensemble member and year dimensions.

For information on the other user configurable variables, outputs and internal functions, see the comments written in the
notebook and scripts. 


Anaconda installation instructions
----------------------------------
To install anaconda and my python/R packages:

1. Download anaconda from https://www.anaconda.com
You want a python>=3.7, 64bit (x86) installer for linux.
As of 08/06/2020 this can be obtained by running: wget https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh
But note that the link to the most up to date version will likely change, so don't rely on this.

2. Install anaconda by running: bash the-script-you-just-downloaded.sh
By default it will install in your home directory, requiring no admin/root permissions.
At the end of the installation it will ask you whether you want to run conda-init.
Say no to this.

3. Add the following line to your .bashrc file in your home directory. If you don't have one, create one.
alias loadconda='export PATH="/path/to/homedir/anaconda3/bin:$PATH"'
Obtain the full path to your home directory by running pwd in your home directory.

4. Start a new bash shell by running bash, then load anaconda by running loadconda (the command we just created above).

5. Setup a 'DAenv' environment that contains all the modules you need to run my scripts,
using the conda_env_file_DAenv.txt. (You can change the name if you want). 
In your home directory run 'bash' if you are not already using the bash shell (by default UKCEH uses the csh shell)
Then conda create --name DAenv --file conda_env_file_DAenv.txt (after copying the txt file to your home dir).

6. Run 'source activate DAenv' to load this environment, and you should be good to go.

Now, instead of sourcing the python and R modules from the default setup, or your own, it will source them
from ~/anaconda3/envs/DAenv/.
To get back to your own python/R setup, simply open a new terminal/shell.

Anytime you want to run these scripts, or use my python setup, run:
bash
loadconda
source activate DAenv

