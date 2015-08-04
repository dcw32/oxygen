import numpy as np
import sys
def k(krt,T,M):
	#k3 is really k3*[M] in order to make the equation quasi bimolecular
	if krt=='k3M':
	        k3=6E-34*(300/T)*(300/T)*M
        	return k3
	elif krt=='k4':
	        k4=1E-11*np.exp(-2100/T)
	        return k4
	else:
		print "ERROR: K RATE NOT DEFINED"
