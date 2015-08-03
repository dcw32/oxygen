import numpy as np
import sys
def k(krt,T,M):
	if krt=='k3':
	        k3=6E-34*(300/T)*(300/T)*M
        	return k3
	elif krt=='k4':
	        k4=1E-11*np.exp(-2100/T)
	        return k4
	else:
		print "ERROR: K RATE NOT DEFINED"
