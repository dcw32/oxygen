# Regrids data from one grid spacing to another, correcting for total mass
# Units of molec are molecules km cm-2... caveat emptor
from scipy.interpolate import interp1d
import numpy as np
def regrid(species,constants):
	molec=np.sum(constants['box_hold']*species)
	f=interp1d(np.linspace(0,80,81),species,kind='cubic')
	spec_new=f(constants['heights'])
	spec_new[spec_new<0]=0.
	molec_n=np.sum(constants['box_h']*spec_new)
	if molec_n>0:
		spec_new=molec*spec_new/molec_n
	return spec_new
