# Calculated the optical depth
import numpy as np
def od_chems(d1,d1defs,o2_c,o3_c,sol,nlevs,wbins,h_max,h_min):
	o2col=np.zeros([nlevs,wbins])
	o3col=np.zeros([nlevs,wbins])
	for i in range(nlevs):
		for j in range(wbins):
			o2col[i,j]=(i+1)*(1E5*(h_max-h_min)/(nlevs-1))*np.sum(d1[np.where(d1defs=='O2')[0][0],:i+1])
			o3col[i,j]=(i+1)*(1E5*(h_max-h_min)/(nlevs-1))*np.sum(d1[np.where(d1defs=='O3')[0][0],:i+1])
	opd=o2col*o2_c+o3col*o3_c
	I_factor=np.exp(-opd)
	I=I_factor*sol
	return I_factor
