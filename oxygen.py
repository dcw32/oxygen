# This file is the main control routine for the stratospheric ozone code
# July 2015, dcw32@cam.ac.uk
# run as python oxygen.py
import numpy as np
import matplotlib.pyplot as plt

from data import extract_csec,savetonetcdf
from j_rates import j
from k_rates import k
from steady_state import steady,ozone,otp
from plotting import altconc,linoxyoz,logoxyoz,inter_plot2
from solver import solve

print "STARTING PROGRAM OXYGEN"

################################################################################
# RUNNING OPTIONS
################################################################################

# Oxygen ratio
ratio = 0.21
print ratio
# Steady state mode - can set up 'initial values' for chem scheme
steadystate = True
# Plotting scripts - plots oxygen vs ozone for a range of oxygen vals
plot_on = True
# Interactive plotting tool for steady state chapman chemistry
interactive = False
# Numerical ODE solver
new_numerics = False

################################################################################
# CONSTANTS
################################################################################

# Change these at your peril!
#Number of atmospheric levels
nlevs=81
#Height, km
h_min=0
h_max=80
heights=np.linspace(h_max,h_min,nlevs)
#scale height, km
H=7
#surface M, cm-3
M_surf=2.5E19*(ratio+0.79)
#M_surf=2.5E19*0.5
M=M_surf*np.exp(-heights/H)

box_h=(1E5*(h_max-h_min)/(nlevs-1))

# wavelength bins, interpolating o2 & solar onto o3 csec wavelengths, extract temps
constants=extract_csec(heights)

# Creates a dictionary with all the constants
key=['H','nlevs','heights','M_surf','ratio','ratio_copy','box_h','M','soln']
vals=[H,nlevs,heights,M_surf,ratio,ratio,box_h,M,constants['sol']]

for i in range(len(key)):
	constants[key[i]]=vals[i]

################################################################################
# MAIN ROUTINES
################################################################################

def main():
	if steadystate==True:
		print "CALCULATING STEADY STATE OZONE COLUMN (CHAPMAN)"
		o3,o2,o,J_o2,J_o3,o3_running=steady(constants,True)
		print "J(O2)="+str(J_o2[40])
		print "J(O3)="+str(J_o3[40])
		altconc(constants,o3,o3_running)

	#Interactive plot of ozone vs o2
		if interactive==True:
	        	print "INTERACTIVE PLOTTING"
			inter_plot2(o3,J_o2,J_o3,constants)

		if plot_on==True:
			print "RUNNING PLOTTING SCRIPTS"
			logoxyoz(constants)
			linoxyoz(constants)

	if new_numerics==True:
		print "RUNNING NUMERICAL SOLVER"
		d1defs,d1,d3=solve(constants)
		savetonetcdf(constants,d1defs,d1,d3)
		
	        o3=d1[np.where(d1defs=='O3')[0][0],:]
	        plt.plot(heights,d3[np.where(d1defs=='O3')[0][0],:])
		plt.plot(heights,d1[np.where(d1defs=='O3')[0][0],:])
		#plt.savefig('ozone_out.png')
		plt.show()
		o3=o3*constants['box_h']
	        du=np.sum(o3)/2.69E16
	        print du
		o3=d3[np.where(d1defs=='O3')[0][0],:]
		o3=o3*constants['box_h']
		du=np.sum(o3)/2.69E16
		print du
		#print d3[np.where(d1defs=='H2O')[0][0],:]

	print "OXYGEN.PY COMPLETE"

################################################################################
# END
################################################################################

if __name__ ==  "__main__":
	main()
