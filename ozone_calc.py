# This routine calculates the ozone at the model level, given an overhead ozone column, photolysis rates and rate constants
import numpy as np
def ozone(J_o2,J_o3,T,o2):
	k3=6E-34*(300/T)*(300/T)
	k4=1E-11*np.exp(-2100/T)
	q=(J_o2/J_o3)*(k3/k4)*np.power(o2,3)/0.21
	o3=np.sqrt(q)
	return o3
