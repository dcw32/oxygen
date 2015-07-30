# Chemical Solver
# Using the Backward Euler Method
import numpy as np
from netCDF4 import Dataset
import os, sys
def d1_init(d1defs,nlevs,ratio):
	d1=np.zeros([len(d1defs[:,0]),nlevs])
	for i in range(len(d1defs[:,0])):
#Perhaps this could be substantially tidied up...
		if d1defs[i,1]=='y':
			print d1defs[i,0]
			file=Dataset('netcdf/'+d1defs[i,0]+'.nc')
			species=file.variables[d1defs[i,0]][:]
			d1[i,:]=species
			file.close()
                elif d1defs[i,0]=='O2':
                        file=Dataset('netcdf/O2.nc')
                        species=file.variables['O2'][:]
                        d1[i,:]=species
                        file.close()
		elif d1defs[i,0]=='M':
			file=Dataset('netcdf/O2.nc')
			species=file.variables['O2'][:]
			species=species/ratio
			d1[i,:]=species
			file.close()
		else:
			#Something's going catastrophically wrong
			print >> sys.stderr, "ERROR: STEADY STATE SPECIES NOT DEFINED IN D1_INIT ARRAY"
			sys.exit(1)
	print d1
	return d1
def chem_ten(d1,d1defs,nlevs):
	bimol=np.genfromtxt("bimol.dat",dtype='str',skiprows=2)
	photo=np.genfromtxt("photol.dat",dtype='str',skiprows=2)
	nrxns=len(bimol)+len(photo)
	#Initialise array of chemical tendencies
	#For each reaction in the scheme
	chem_tend=np.zeros([nrxns,nlevs])
	
