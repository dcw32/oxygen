# This routine calculates the ozone at the model level, given an overhead ozone column, photolysis rates and rate constants
import numpy as np
import sys
def ozone(J_o2,J_o3,T,o2,ratio,k3M,k4):
	#[O3]ss=(J(O2)*k3*[M]*[O2]^2/J(O3)*k4)^0.5
	q=(J_o2/J_o3)*(k3M/k4)*np.power(o2,2)
	if q<0:
		print >> sys.stderr, "NEGATIVE Q VALUE"
		sys.exit(1)
	o3=np.sqrt(q)
	return o3,q
def otp(o3,J_o3,k3M,o2,ratio):
	#[O]ss=[O3]*J(O3)/k3*[M]*[O2]
	ot=o3*J_o3/(k3M*o2)
	return ot
