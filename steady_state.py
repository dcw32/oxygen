# Calculates ozone profile in steady-state
# 
from j_rates import j
from k_rates import k
import numpy as np
import matplotlib.pyplot as plt
import sys
from netCDF4 import Dataset

def steady(constants,init):
	o3=np.zeros(constants['nlevs'])
	o3_sum=np.zeros(constants['nlevs'])
	o2=np.zeros(constants['nlevs'])
	o2_sum=np.zeros(constants['nlevs'])
	o=np.zeros(constants['nlevs'])
	J_o2=np.zeros(constants['nlevs'])
	J_o3=np.zeros(constants['nlevs'])
	k3M=np.zeros(constants['nlevs'])
	k4=np.zeros(constants['nlevs'])
	q=np.zeros(constants['nlevs'])
	o2_running=0
	o3_running=0
	for i in range(constants['nlevs']):
        	o2[i]=constants['ratio']*constants['M'][i]
	        optical_depth=constants['o2_c']*o2_running+constants['o3_c']*o3_running
	        I_factor=np.exp(-optical_depth)
        	I=np.array(constants['sol']*I_factor)
        	J_o2[i]=j(I,'JO2',constants)
        	J_o3[i]=j(I,'JO3',constants)
        	k3M[i]=k('k3M',constants,i)
        	k4[i]=k('k4',constants,i)
        	o3[i],q[i]=ozone(J_o2[i],J_o3[i],k3M[i],k4[i],o2[i],constants['M'][i])
        	o[i]=otp(o3[i],J_o3[i],k3M[i],o2[i])
        	o2_running=o2_running+o2[i]*constants['box_h'][i]
        	o2_sum[i]=o2_running
        	o3_running=o3_running+o3[i]*constants['box_h'][i]
        	o3_sum[i]=o3_running
	if init==True:
		file=Dataset('netcdf/initial.nc','w')
		file.createDimension('height',constants['nlevs'])
		spec_in_init=['O3','O','O2','M','CH4','N2O','H2O']
		val_in_init=[o3,o,o2,constants['M'],1.7E-6*constants['M'],3.2E-7*constants['M'],4.E-6*constants['M']]
		for i in range(len(spec_in_init)):
			file.createVariable(spec_in_init[i],np.float64,'height')
			file.variables[spec_in_init[i]][:]=val_in_init[i]
		#file.variables['O3'][:]=o3
		#file.variables['O'][:]=o
		#file.variables['O2'][:]=o2
		#file.variables['M'][:]=constants['M']
		#file.variables['CH4'][:]=1.7E-6*constants['M']
		#file.variables['N2O'][:]=3.2E-7*constants['M']
		#file.variables['H2O'][:]=4.E-6*constants['M']
	#plt.plot(constants['heights'],o3)
	#plt.savefig('ozone.png')
	return o3,o2,o,J_o2,J_o3,o3_running
# This routine calculates the ozone at the model level, given an overhead ozone column, photolysis rates and rate constants
def ozone(J_o2,J_o3,k3M,k4,o2,M):
	q=(J_o2/J_o3)*(k3M/k4)*o2*M*0.42
	if q<0:
		print >> sys.stderr, "NEGATIVE Q VALUE"
		sys.exit(1)
	o3=np.sqrt(q)
	return o3,q
def otp(o3,J_o3,k3M,o2):
	#[O]ss=[O3]*J(O3)/k3*[M]*[O2]
	ot=o3*J_o3/(k3M*o2)
	return ot
