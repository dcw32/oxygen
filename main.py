# This file is the main control routine for the stratospheric ozone code
# July 2015, dcw32@cam.ac.uk
# run as ipython --pylab
# %run ./main.py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.text as txt
from matplotlib.widgets import Slider, Button, RadioButtons
import sys
from pylab import *
from netCDF4 import Dataset

from data import extract_csec,ncsave
from optical_depth import od_chems_new,od_chems
from j_rates import j
from k_rates import k
from steady_state import steady,ozone,otp
from plotting import altconc,linoxyoz,logoxyoz
from solver import solve

print "STARTING PROGRAM OXYGEN"

# Running options

# Oxygen ratio
ratio = 0.21
# Steady state mode - can set up 'initial values' for chem scheme
steadystate = True
# Plotting scripts - plots oxygen vs ozone for a range of oxygen vals
plot_on = True
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

box_h=(1E5*(h_max-h_min)/(nlevs-1))

# wavelength bins, interpolating o2 & solar onto o3 csec wavelengths, extract temps
sol_bin_width,o2_c,o3_c,sol,T=extract_csec(heights)

###################################################################################
# MAIN ROUTINES
###################################################################################

if steadystate==True:
	print "CALCULATING STEADY STATE OZONE COLUMN (CHAPMAN)"
	height,o3,o2,o,J_o2,J_o3,o3_running=steady(nlevs,h_max,h_min,H,M_surf,ratio,o2_c,o3_c,T,\
                                                   sol,sol_bin_width,True)
	altconc(height,o3,o3_running)

#Interactive plot of ozone vs o2
	if interactive==True:
		print "INTERACTIVE PLOTTING"
		fig, (ax2,ax1,ax0) = plt.subplots(3, sharex=True)
		#fig1=plt.subplots()
		plt.subplots_adjust(bottom=0.25)
		ax0.set_xlabel('Height / km')
		ax0.set_ylabel('Ozone / 10^-12 cm-3')
		ax0.set_ylim([0,12])
		ax1.set_ylabel('JO2 / 10^10 s-1')
		ax2.set_ylabel('JO3 / 10^3 s-1')
		l,=ax0.plot(height,o3/1E12)
		m,=ax1.plot(height,J_o2*1E10,color='red')
		n,=ax2.plot(height,J_o3*1E3,color='purple')
		#ann=txt.Annotation(str(int(round(du))), xy=(50,10),xycoords='data')		
		#ax0.add_artist(ann)
		#Sets up sliders
		axo2=plt.axes([0.25,0.1,0.65,0.03],axisbg='lightgoldenrodyellow')
		so2=Slider(axo2, 'O\mathrm{_{2}}', 0.001, 1.00, valinit=0.21)
		#Updates graphs based on slider value
		def update(val):
			ratio=so2.val
			M_surf=2.5E19*(ratio+0.79)
			height,o3,o2,o,J_o2,J_o3,o3_running=steady(nlevs,h_max,h_min,H,\
                                               M_surf,ratio,o2_c,o3_c,T,sol,sol_bin_width,False)
			m.set_ydata(J_o2*1E10)
			l.set_ydata(o3/1E12)
			n.set_ydata(J_o3*1E3)
			du=o3_running/2.69E16
	                #ax0.annotate(str(int(round(du))), xy=(0.99,0.01),xycoords='axes fraction',horizontalalignment='right', verticalalignment='bottom')
			fig.canvas.draw_idle()
		so2.on_changed(update)
		plt.show()

	if plot_on==True:
		print "RUNNING PLOTTING SCRIPTS"
		logoxyoz(nlevs,h_max,h_min,H,o2_c,o3_c,T,sol,sol_bin_width)
		linoxyoz(nlevs,h_max,h_min,H,o2_c,o3_c,T,sol,sol_bin_width)

if new_numerics==True:
	d1defs,d1,d3=solve(o2_c,o3_c,sol,sol_bin_width,box_h,nlevs,ratio,T)
	plt.plot(heights,(d1[np.where(d1defs=='O3')[0][0],:]\
                  -d3[np.where(d1defs=='O3')[0][0],:])\
                  /d1[np.where(d1defs=='O3')[0][0],:])
	plt.savefig('new_numerics.png')
	print d3[np.where(d1defs=='H2O')[0][0],:]

#######################################################################################
# END
#######################################################################################
print "END OF PROGRAM"
