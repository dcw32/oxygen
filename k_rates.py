import numpy as np
import sys
def k(krt,constants,i):
#to add new rate:
#elif krt=='RATE':
#	k=...
	#k3 is really k3*[M] in order to make the equation quasi bimolecular
	if krt=='k3M':
	        k=6E-34*(300/constants['T'][i])*(300/constants['T'][i])*constants['M'][i]
	elif krt=='k4':
	        k=1E-11*np.exp(-2100/constants['T'][i])
	else:
		print "ERROR: K RATE NOT DEFINED"
	return k
