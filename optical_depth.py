# Calculated the optical depth
import numpy as np
def od(o2_c,o3_c,o2_col,o3_col):
	optical_depth=o2_c*o2_col+o3_c*o3_col
	I_factor=np.exp(-optical_depth)
	return I_factor
