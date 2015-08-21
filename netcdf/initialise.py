from scipy import interpolate
import numpy as np
from netCDF4 import Dataset
species='N2O'
dataset=Dataset('initial.nc','r+',format='NETCDF4')
height=np.linspace(80,0,81)
#dataset.createDimension('height',len(height))
#heights=dataset.createVariable('height',np.float32,('height',))
#heights[:]=height
spec=dataset.variables[species]
marray=Dataset('M.nc','r')
m=marray.variables['M'][:]
m=m[::-1]
array_spec=np.ones(81)
for i in range(10,50):
	array_spec[i]=1E-9*(320-(i-10)*8)*m[i]
for i in range(0,10):
	array_spec[i]=320E-9*m[i]
array_spec=array_spec[::-1]
spec[:]=array_spec
