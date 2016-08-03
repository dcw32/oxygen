# This file is the main control routine for the stratospheric ozone code
# July 2015, dcw32@cam.ac.uk
#
# Change the running options as desired 
# Then run the code as python oxygen.py
#
import numpy as np
import matplotlib.pyplot as plt
import sys

from init import init
from data import savetonetcdf
from j_rates import j
from k_rates import k
from steady_state import steady,ozone,otp
from plotting import *
from solver import solve
from driver import *

print "STARTING PROGRAM OXYGEN"

################################################################################
# RUNNING OPTIONS
################################################################################

# Oxygen ratio
ratio = 0.21
print str(ratio)+" OXYGEN CONTENT"
# Steady state mode - can set up 'initial values' for chem scheme
steadystate = True
# Plotting scripts - plots oxygen vs ozone for a range of oxygen vals
plot_on = True
# Interactive plotting tool for steady state chapman chemistry
interactive = False
# Numerical ODE solver
new_numerics = True

################################################################################
# MAIN ROUTINES
################################################################################

def main():
	# Sets up the model grid and constants
	constants=init(ratio,False)
	# Steady State Chapman Chemistry
	if steadystate==True:
		#if constants['nlevs']==81:
		steady_state_chapman(constants)
		#else:
		#	print "NEEDS 81 EVEN LEVS FOR SS"
		#	sys.exit()
		#Interactive plot of ozone vs o2
		if interactive==True:
	        	print "INTERACTIVE PLOTTING"
			inter_plot2(o3,J_o2,J_o3,constants)

		if plot_on==True:
			print "RUNNING PLOTTING SCRIPTS"
			logoxyoz(constants)
			linoxyoz(constants)

	if new_numerics==True:
		numerical(constants)
	print "OXYGEN.PY COMPLETE"

################################################################################
# END
################################################################################

if __name__ ==  "__main__":
	main()
