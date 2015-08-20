from scipy import interpolate
import numpy as np
from netCDF4 import Dataset
species='ClO'
dataset=Dataset(species+'.nc','w',format='NETCDF4')
height=np.linspace(80,0,81)
dataset.createDimension('height',len(height))
heights=dataset.createVariable('height',np.float32,('height',))
heights[:]=height
spec=dataset.createVariable(species,np.float32,('height',))
marray=Dataset('M.nc','r')
m=marray.variables['M'][:]
m=m[::-1]
array_spec=np.zeros(81)
for i in range(81):
	array_spec[i]=3E-9*m[i]
array_spec=array_spec[::-1]
spec[:]=array_spec
