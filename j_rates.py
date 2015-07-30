#This script calculates the J(O2) and J(O3) photolysis rates required in calculating the steady state ozone concentration
import numpy as np
def j(I,o2_c,o3_c,jrate,bin_w):
	if jrate=='JO2':
		j_spec=np.multiply(o2_c,I)
		j_spec=np.multiply(j_spec,bin_w)
		j_spec=np.sum(j_spec)
	elif jrate=='JO3':
                j_spec=np.multiply(o3_c,I)
                j_spec=np.multiply(j_spec,bin_w)
                j_spec=np.sum(j_spec)
	return j_spec
