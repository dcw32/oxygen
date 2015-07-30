# This routine calculates the ozone at the model level, given an overhead ozone column, photolysis rates and rate constants
import numpy as np
import sys
def ozone(J_o2,J_o3,T,o2,ratio,k3,k4):
	q=(J_o2/J_o3)*(k3/k4)*np.power(o2,3)/ratio
	if q<0:
		print >> sys.stderr, "NEGATIVE Q VALUE"
		sys.exit(1)
	o3=np.sqrt(q)
	return o3,q
def otp(o3,J_o3,k3,o2,ratio):
	ot=o3*J_o3*ratio/(k3*o2*o2)
	return ot
