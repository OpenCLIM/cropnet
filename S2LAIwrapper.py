import calendar
import sys
import os
print(sys.executable)
sys.stdout.flush()
from S2_TOA_TO_LAI import TOA2LAI_S2



def days_in_month(month, year):
    # simple function to find the number of days in a given month of a given year
    month = int(month)
    year = int(year)
    if month in [9, 4, 6, 11]:
        days = 30
    if month in [1, 3, 5, 7, 8, 10, 12]:
        days = 31
    if month == 2:
        if calendar.isleap(year):
            days = 29
        else:
            days = 28
    return str(days)

month = sys.argv[3]
year = sys.argv[2]
tile = sys.argv[1]
dim = days_in_month(month, year)
start = str(year) + '-' + str(month) + '-01'
end = str(year) + '-' + str(month) + '-' + dim
currdir = os.getcwd()
mcd43 = os.path.join(currdir, 'MCD43')
vrt_dir = os.path.join(currdir, 'MCD43_VRT')
if not os.path.exists(mcd43):
    os.makedirs(mcd43)
if not os.path.exists:
    os.makedirs(vrt_dir)

ret = TOA2LAI_S2(tiles = [tile], start=start, end=end, mcd43=mcd43, vrt_dir=vrt_dir)
