#This script calculates the J(O2) and J(O3) photolysis rates required in calculating the steady state ozone concentration
import numpy as np
def j(I,csec,bin_w):
	j_spec=np.multiply(csec,I)
	j_spec=np.multiply(j_spec,bin_w)
	j_spec=np.sum(j_spec)
	return j_spec
