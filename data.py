from scipy import interpolate
import numpy as np
from netCDF4 import Dataset
def extract_csec():
	#Extraction of solar data
	solar_data=np.loadtxt("solar.dat")
	lambda_solr,solr=solar_data[:,0],solar_data[:,1]
	#Extraction of cross-sections
	o2_data=np.loadtxt("o2_csec.dat",skiprows=2)
	lambda_o2,o2_csec=o2_data[:,0],o2_data[:,1]
	o3_data=np.loadtxt("o3_csec.dat",skiprows=2)
	lambda_bin,o3_binned_csec=o3_data[:,0],o3_data[:,1]
	# Interpolate solar fluxes and o2 onto o3 lambda vals
	lambda_max_o2=np.max(lambda_o2)
	o2_fn=interpolate.UnivariateSpline(lambda_o2,o2_csec)
	o2_binned_csec=o2_fn(lambda_bin)
	for i in range(lambda_bin.shape[0]):
        	if lambda_bin[i]>lambda_max_o2:
                	o2_binned_csec[i]=0
	sol_fn=interpolate.interp1d(lambda_solr,solr)
	sol_binned=sol_fn(lambda_bin)
	return lambda_bin,o2_binned_csec,o3_binned_csec,sol_binned
def ncsave(height,species,species_arr):
	dataset=Dataset('netcdf/'+species+'.nc','w',format='NETCDF4')
	dataset.createDimension('height',len(height))
	heights=dataset.createVariable('height',np.float32,('height',))
	heights[:]=height
	spec=dataset.createVariable(species,np.float32,('height',))
	spec[:]=species_arr
