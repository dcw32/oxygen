#This script calculates the J(O2) and J(O3) photolysis rates required in calculating the steady state ozone concentration
import numpy as np
def j(I,csec):
	j_spec=np.multiply(csec,I)
	return j_spec
