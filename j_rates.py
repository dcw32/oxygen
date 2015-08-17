#This script calculates the J(O2) and J(O3) photolysis rates required in calculating the steady state ozone concentration
import numpy as np
def j(I,jrate,constants):
#to add a rxn:
#elif jrate=='JRATE':
#	j_spec=value
	if jrate=='JO2':
		j_spec=np.multiply(constants['o2_c'],I)
		j_spec=np.multiply(j_spec,constants['sol_bin_width'])
		j_spec=np.sum(j_spec)
		j_spec=j_spec
	elif jrate=='JO3':
                j_spec=np.multiply(constants['o3_c'],I)
                j_spec=np.multiply(j_spec,constants['sol_bin_width'])
                j_spec=np.sum(j_spec)
		j_spec=j_spec
	return j_spec	
