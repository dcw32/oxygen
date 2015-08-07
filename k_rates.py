import numpy as np
import sys
def k(krt,T,M):
#to add new rate:
#elif krt=='RATE':
#	k=...
	#k3 is really k3*[M] in order to make the equation quasi bimolecular
	if krt=='k3M':
	        k=6E-34*(300/T)*(300/T)*M
	elif krt=='k4':
	        k=1E-11*np.exp(-2100/T)
	else:
		print "ERROR: K RATE NOT DEFINED"
	return k
