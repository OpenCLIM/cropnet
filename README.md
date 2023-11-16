CropNET Models with Data Assimilation 
-------------------------------------

Overview
--------
This tool produces estimates of potential yield using a grass, wheat or oilseed rape (OSR) crop model for a given year.

Details of the winter wheat crop model used can be found in Redhead et al. (2020): The influence of landscape composition and configuration on crop yield resilience, Journal of Applied Ecology, 57, 2180-2190, https://doi.org/10.1111/1365-2664.13722

Data requirements
-----------------
The code has been designed to work with CHESS-SCAPE, UKCP18 or ERA5 data, with the optional addition of APHRODITE precipitation data 
for running over China.
12km UKCP18 data is available on the CEDA archive at 
https://catalogue.ceda.ac.uk/uuid/589211abeb844070a95d061c8cc7f604
ERA5 data is available from the copernicus climate data store:
https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview
APHRODITE data is available from:
http://aphrodite.st.hirosaki-u.ac.jp/
There are switches in the code that allows you to specify whether CHESS-SCAPE, UKCP18 or ERA5 data is being used, and whether you want to
use APHRODITE precipitation data. The variable names, filenames and coordinate names are all preset for these datasets.
If you are using ERA5 you will need to apply a daily average or sum (as appropriate) to the downloaded data. There is an option
to do this automatically in the code, but this feature should be considered as under development.

If you wish to use a different meteorological driving dataset you can edit the 'getnames' function in utils.py, using the
docstring as guidance. The dataset name is provided to the code as an argument and the translation of this to the folder
containing the data will also need to be manually edited in the jasmin-wrapper.py script. See 'Running the code' section for more details on the arguments required to run the model. 
The data required by the models are:
Wheat and OSR:
- Daily precipitation
- Daily maximum temperature
- Daily minimum temperature
- Daily mean temperature
- Daily average net solar radiation
Grass:
- Daily precipitation
- Daily maximum temperature
- Daily minimum temperature
- Daily mean temperature
- Daily average surface air pressure
- Daily average net solar radiation
- One of: daily mean dewpoint temperature or relative humidity
- One of: daily mean windspeed, or daily mean x and y components of windspeed

Available Water Content (AWC) is also required by the model. This data is provided for the UK as in the data/AWC folder. 
The data was derived from the UKCEH Grid-to-Grid model, as the average annual maximum soil moisture in the top 1.6m of soil at each 1km grid point over the UK, multiplied by 0.65 to account for the proportion of water typically available to the wheat plant (Barraclough and Leigh, 1984).
For China, data can be obtained using the 'get_AWC_China.sh' bash script. 

Annual CO2 concentration is an optional requirement. Concentrations corresponding to the ensemble members availabe in the CHESS-SCAPE dataset are available in the CO2_files folder.
They are derived from the RCP database at https://tntcat.iiasa.ac.at/RcpDb/dsd?Action=htmlpage&page=download 

Crop and irrigation maps are also available for use in the code for China in the wheatmapchina.nc and 
wheat_irrigated_proportion.nc files. They were generated from MAPSPAM data, example code for generating
equivalent maps for other areas is in the verdant.ipynb notebook.


Running the code
----------------
The code is written in python, which calls the crop model written in R.
It is run from the jasmin-wrapper.py script which is designed to be called either from the shell or from a shell script using the following call:
jasmin-wrapper.py arg1 arg2 arg3 arg4 arg5
- arg1: crop to model. 'grass', 'wheat' or 'OSR'
- arg2: driving dataset name. 'ukcp18', 'chess-scape_8.5_01', 'chess-scape_8.5_04', 'chess-scape_8.5_06', 'chess-scape_8.5_15', 'chess-scape_2.6_01', 'chess-scape_2.6_04', 'chess-scape_2.6_06', 'chess-scape_2.6_15', 'era5', 'chess_and_haduk'
- arg3: The CO2 Scenario to use for CO2 fertilisation. '8.5', '8.5_01', '8.5_04', '8.5_06', '8.5_15', '2.6' or 'None'. 'None' disables CO2 fertilisation.
- arg4: Year in which to start the simulation
- arg5: Month in which to start the simulation
arg2 determines which folder of input data will be read in and defines the variable names found in the filename and netcdf variable names. The code in jasmin-wrapper.py which handles the translation of arg2 to the folder containing the data would need to be edited to handle any other input/driving dataset than the ones described here. 


The code uses several python/R packages that are not installed by default on most systems.
The easiest way to install these is to use anaconda - a package management system for python and R modules/packages.
Anaconda installs packages in your home directory, and can be managed separately to any other python/R installations
you have, i.e. it won't interfere with these, and can be enabled only when you want to run this programme.
See below for instructions on installing and using anaconda.


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

