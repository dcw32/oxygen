import numpy as np
import sys
def k3_i(T):
        k3=6E-34*(300/T)*(300/T)
	return k3
def k4_i(T):
        k4=1E-11*np.exp(-2100/T)
	return k4
