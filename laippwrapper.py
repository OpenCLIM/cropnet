from MODIS import *
import multiprocessing as mp

S2dir = '/users/sgsys/matbro/cropnet/data/S2/LAI'
yieldfile = '/users/sgsys/matbro/cropnet/data/cropyield/Mean_wheat_yields_v3/hayman_etal_val_data.shp'

pool = mp.Pool(10) # no. of concurrent processes
results = []
for fieldno in range(1100,1169):
    results.append(pool.apply_async(process_S2, [S2dir, yieldfile, "2015-01-01", "2019-12-31", 20, fieldno]))

print([res.get() for res in results])
pool.close()
pool.join()

    
