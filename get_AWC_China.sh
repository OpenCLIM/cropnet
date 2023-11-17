#!/bin/bash

# Edit or comment out the below two lines to use your own environment
loadconda
source activate ceh

wget http://globalchange.bnu.edu.cn/download/data/ptf/TH33.zip
wget http://globalchange.bnu.edu.cn/download/data/ptf/TH1500.zip

unzip TH1500.zip
unzip TH33.zip
rm *.zip

cdo delname,TH33cv -chname,TH33,TH1500 TH33.nc TH_FC.nc
cdo delname,TH1500cv TH1500.nc TH_WP.nc

cdo sub TH_FC.nc TH_WP.nc AWC_temp.nc
cdo chname,TH1500,AWC AWC_temp.nc AWC_China.nc

rm AWC_temp.nc TH1500.nc TH33.nc
