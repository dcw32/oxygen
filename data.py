from scipy import interpolate
import numpy as np
from netCDF4 import Dataset
def extract_csec(heights):
	#Extraction of solar data
	solar_data=np.loadtxt("solar.dat")
	lambda_solr,solr=solar_data[:,0],solar_data[:,1]
	#Extraction of cross-sections
	o2_data=np.loadtxt("o2_csec.dat",skiprows=2)
	lambda_o2,o2_csec=o2_data[:,0],o2_data[:,1]
	o3_data=np.loadtxt("o3_csec.dat",skiprows=2)
	lambda_bin,o3_binned_csec=o3_data[:,0],o3_data[:,1]
	o3phi_data=np.loadtxt("o3_qyield.dat",skiprows=2)
	o3_qyield=o3phi_data[:,1]
	h2o2_data=np.loadtxt("h2o2_csec.dat",skiprows=2)
	lambda_h2o2,h2o2_csec=h2o2_data[:,0],h2o2_data[:,1]
        clono2_data=np.loadtxt("clono2_csec.dat",skiprows=2)
        lambda_clono2,clono2_csec=clono2_data[:,0],clono2_data[:,1]
	# Interpolate solar fluxes and o2 onto o3 lambda vals
	lambda_max_o2=np.max(lambda_o2)
	o2_fn=interpolate.UnivariateSpline(lambda_o2,o2_csec)
	o2_binned_csec=o2_fn(lambda_bin)
	for i in range(lambda_bin.shape[0]):
        	if lambda_bin[i]>lambda_max_o2:
                	o2_binned_csec[i]=0
	lambda_max_h2o2=np.max(lambda_h2o2)
	h2o2_fn=interpolate.interp1d(lambda_h2o2,h2o2_csec,bounds_error=False)
	h2o2_binned_csec=h2o2_fn(lambda_bin)
	for i in range(lambda_bin.shape[0]):
		if lambda_bin[i]>lambda_max_h2o2:
			h2o2_binned_csec[i]=0
	h2o2_binned_csec[h2o2_binned_csec<0]=0
	lambda_max_clono2=np.max(lambda_clono2)
	clono2_fn=interpolate.interp1d(lambda_clono2,clono2_csec,bounds_error=False)
        clono2_binned_csec=clono2_fn(lambda_bin)
        for i in range(lambda_bin.shape[0]):
                if lambda_bin[i]>lambda_max_clono2:
                        clono2_binned_csec[i]=0
        clono2_binned_csec[clono2_binned_csec<0]=0
	sol_fn=interpolate.interp1d(lambda_solr,solr)
	sol_binned=sol_fn(lambda_bin)
	#convert sol to photons cm-2 s-1 nm-1 for binning
	sol_binned=sol_binned/0.05
	sol_bin_width=np.zeros(len(lambda_bin))
	sol_bin_width[0]=2
	sol_bin_width[len(lambda_bin)-1]=5
	for i in range(1,len(lambda_bin)-2):
        	sol_bin_width[i]=0.5*(lambda_bin[i+1]-lambda_bin[i-1])
	temp_data=np.loadtxt("temp.dat")
	ht,temp=temp_data[:,0],temp_data[:,1]
	temp_fn=interpolate.UnivariateSpline(ht,temp)
	T=temp_fn(heights)
	o2_binned_csec=np.array(o2_binned_csec)
	o3_binned_csec=np.array(o3_binned_csec)
	constants={}
	key=['sol_bin_width','o2_c','o3_c','sol','T','o3_q','h2o2_c','clono2_c','lambda']
	vals=[sol_bin_width,o2_binned_csec,o3_binned_csec,sol_binned,T,o3_qyield,h2o2_binned_csec,clono2_binned_csec,lambda_bin]
	for i in range(len(key)):
		constants[key[i]]=vals[i]
	return constants
def ncsave(height,species,species_arr):
	dataset=Dataset('netcdf/'+species+'.nc','w',format='NETCDF4')
	dataset.createDimension('height',len(height))
	heights=dataset.createVariable('height',np.float32,('height',))
	heights[:]=height
	spec=dataset.createVariable(species,np.float32,('height',))
	spec[:]=species_arr
def savetonetcdf(constants,d1defs,d1,d3):
	height=constants['heights']
	file1=Dataset('netcdf/results.nc','w',format='NETCDF4')
	file1.createDimension('height',len(height))
        heights=file1.createVariable('height',np.float32,('height',))
        heights[:]=height
	for i in range(len(d1defs)):
        	spec=file1.createVariable(d1defs[i,0],np.float32,('height',))
		spec[:]=d3[i,:]
