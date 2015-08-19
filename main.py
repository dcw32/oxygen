# This file is the main control routine for the stratospheric ozone code
# July 2015, dcw32@cam.ac.uk
# run as ipython --pylab
# %run ./main.py
import numpy as np
import matplotlib.pyplot as plt

from data import extract_csec,ncsave
from j_rates import j
from k_rates import k
from steady_state import steady,ozone,otp
from plotting import altconc,linoxyoz,logoxyoz,inter_plot
from solver import solve

print "STARTING PROGRAM OXYGEN"

#################################################################################
# RUNNING OPTIONS
#################################################################################

# Oxygen ratio
ratio = 0.21
# Steady state mode - can set up 'initial values' for chem scheme
steadystate = False
# Plotting scripts - plots oxygen vs ozone for a range of oxygen vals
plot_on = False
# Interactive plotting tool for steady state chapman chemistry
interactive = False
# Numerical ODE solver
new_numerics = True

#Constants
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
M=M_surf*np.exp(-heights/H)

box_h=(1E5*(h_max-h_min)/(nlevs-1))

# wavelength bins, interpolating o2 & solar onto o3 csec wavelengths, extract temps
sol_bin_width,o2_c,o3_c,sol,T=extract_csec(heights)

# Creates a dictionary with all the constants
constants={}
key=['H','nlevs','heights','M_surf','ratio','o2_c','o3_c','T','sol','sol_bin_width','box_h','M']
vals=[H,nlevs,heights,M_surf,ratio,o2_c,o3_c,T,sol,sol_bin_width,box_h,M]

for i in range(len(key)):
	constants[key[i]]=vals[i]

###################################################################################
# MAIN ROUTINES
###################################################################################

if steadystate==True:
	print "CALCULATING STEADY STATE OZONE COLUMN (CHAPMAN)"
	o3,o2,o,J_o2,J_o3,o3_running=steady(constants,True)
	altconc(constants,o3,o3_running)

#Interactive plot of ozone vs o2
	if interactive==True:
	        print "INTERACTIVE PLOTTING"
		inter_plot(o3,J_o2,J_o3,constants)

	if plot_on==True:
		print "RUNNING PLOTTING SCRIPTS"
		logoxyoz(constants)
		linoxyoz(constants)

if new_numerics==True:
	constants['ratio']=0.21
	constants['M_surf']=2.5E19*(constants['ratio']+0.79)
	constants['M']=constants['M_surf']*np.exp(-constants['heights']/constants['H'])
	print "RUNNING NUMERICAL SOLVER"
	d1defs,d1,d3=solve(constants)
	#plt.plot(heights,(d3[np.where(d1defs=='O3')[0][0],:]\
        #          -d1[np.where(d1defs=='O3')[0][0],:])\
        #          /d1[np.where(d1defs=='O3')[0][0],:])
	#plt.savefig('new_numerics.png')
        o3=d1[np.where(d1defs=='O3')[0][0],:]
        plt.plot(heights,d3[np.where(d1defs=='O3')[0][0],:])
	plt.plot(heights,d1[np.where(d1defs=='O3')[0][0],:])
	plt.show()
	o3=o3*constants['box_h']
        du=np.sum(o3)/2.69E16
        print du
	o3=d3[np.where(d1defs=='O3')[0][0],:]
	o3=o3*constants['box_h']
	du=np.sum(o3)/2.69E16
	print du
	#print d3[np.where(d1defs=='H2O')[0][0],:]

#######################################################################################
# END
#######################################################################################
print "END OF PROGRAM"
