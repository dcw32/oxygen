# These functions construct the D arrays
# D1 holds the concentrations of each chemical species in species.dat
# This is initialised from values calculated in the steady state,
# or from a given initial concentration, all in netCDF files
# D2 holds the array of rates calculated for each rxn
# D3 holds the array of chemical tendencies d[species]/dt, for each species.
import numpy as np
from netCDF4 import Dataset
def d1_init(d1defs,nlevs,ratio):
	d1=np.zeros([len(d1defs[:,0]),nlevs])
	for i in range(len(d1defs[:,0])):
		file=Dataset('netcdf/'+d1defs[i,0]+'.nc')
		species=file.variables[d1defs[i,0]][:]
		d1[i,:]=species
		file.close()
	return d1
