#This script calculates the J(O2) and J(O3) photolysis rates required in calculating the steady state ozone concentration
import numpy as np
import sys
def j(I,jrate,constants):
#to add a rxn:
#elif jrate=='JRATE':
#	j_spec=value
	if jrate=='JO2':
		j_spec=np.multiply(constants['o2_c'],I)
		j_spec=np.multiply(j_spec,constants['sol_bin_width'])
		j_spec=np.sum(j_spec)
	elif jrate=='JO3':
		j_spec=np.multiply(constants['o3_c'],1-constants['o3_q'])
                j_spec=np.multiply(j_spec,I)
                j_spec=np.multiply(j_spec,constants['sol_bin_width'])
                j_spec=np.sum(j_spec)
	elif jrate=='JO1D':
		j_spec=np.multiply(constants['o3_c'],constants['o3_q'])
		j_spec=np.multiply(j_spec,I)
		j_spec=np.multiply(j_spec,constants['sol_bin_width'])
		j_spec=np.sum(j_spec)
	elif jrate=='JH2O2':
		j_spec=np.multiply(constants['h2o2_c'],I)
		j_spec=np.multiply(j_spec,constants['sol_bin_width'])
		j_spec=np.sum(j_spec)
	elif jrate=='JClONO2':
		j_spec=np.multiply(constants['clono2_c'],I)
		j_spec=np.multiply(j_spec,constants['sol_bin_width'])
		j_spec=np.sum(j_spec)
	elif jrate=='JNO2':
		j_spec=1E-2
	elif jrate=='JNO3':
		j_spec=0.156
	else:
		j_spec=0
		sys.exit("ERROR: J RATE NOT DEFINED")
	return j_spec	
