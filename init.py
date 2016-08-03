import numpy as np
from data import extract_csec
#from data import extract_csec
def init(ratio,new_heights):
	# Initialises the model
	#nlevs=41
	#Height, km
	h_min=0
	h_max=80
	#heights=np.linspace(h_max,h_min,nlevs)
	heights=np.array([0.5,2.,4.5,8.,12.5,18.,24.5,32.0,40.5,50.,60.5,72.,84.5])
	heights=heights[::-1]
	nlevs=heights.shape[0]
	box_h=1E5*np.array([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.])
	box_h=box_h[::-1]
	#scale height, km
	if new_heights==True:
		heights2=heights
		nlevs=71
		heights=np.linspace(h_max,h_min,nlevs)
	H=7
	#surface M, cm-3
	M_surf=2.5E19*(ratio+0.79)
	#M_surf=2.5E19*0.5
	M=M_surf*np.exp(-heights/H)
	#box_h=(1E5*(h_max-h_min)/(nlevs-1))
	if new_heights==True:
		M_fix=M_surf*np.exp(-heights2/H)
		box_hold=(1E5*(h_max-h_min)/(80))
	else:
		box_hold=0.
	#
	# wavelength bins, interpolating o2 & solar onto o3 csec wavelengths, extract temps
	constants=extract_csec(heights)
	#
	# Creates a dictionary with all the constants
	key=['H','nlevs','heights','M_surf','ratio','ratio_copy','box_h','M','soln','box_hold']
	vals=[H,nlevs,heights,M_surf,ratio,ratio,box_h,M,constants['sol'],box_hold]
	#
	for i in range(len(key)):
	        constants[key[i]]=vals[i]
	return constants
