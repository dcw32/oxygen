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
from scipy.integrate import odeint,ode

from data import extract_csec,ncsave
from optical_depth import od_chems_new,od_chems
from j_rates import j
from k_rates import k
from steady_state import steady,ozone,otp
from d_arrays import d2_calc,d3_calc,d1_init
from plotting import altconc,linoxyoz,logoxyoz

#Running options
# oxygen ratio
ratio = 0.21
# Steady state mode - can set up 'initial values' for chem scheme
steadystate = True
interactive = False
new_numerics = True
# Plotting scripts - plots oxygen vs ozone for a range of oxygen vals
plot_on = True
## Chemical scheme: ChapmanSS, Chapman
#chem_scheme = 'Chapman'

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

if steadystate==True:
	print "CALCULATING STEADY STATE OZONE COLUMN"
	height,o3,o2,o,J_o2,J_o3,o3_running=steady(nlevs,h_max,h_min,H,M_surf,ratio,o2_c,o3_c,T,sol,sol_bin_width,True)
	altconc(height,o3,o3_running)

#Interactive plot of ozone vs o2
	if interactive==True:
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
			height,o3,o2,o,J_o2,J_o3,o3_running=steady(nlevs,h_max,h_min,H,M_surf,ratio,o2_c,o3_c,T,sol,sol_bin_width,False)
			m.set_ydata(J_o2*1E10)
			l.set_ydata(o3/1E12)
			n.set_ydata(J_o3*1E3)
			du=o3_running/2.69E16
	                #ax0.annotate(str(int(round(du))), xy=(0.99,0.01),xycoords='axes fraction',horizontalalignment='right', verticalalignment='bottom')
			fig.canvas.draw_idle()
		so2.on_changed(update)
		plt.show()

if plot_on==True:
	logoxyoz(nlevs,h_max,h_min,H,o2_c,o3_c,T,sol,sol_bin_width)
	linoxyoz(nlevs,h_max,h_min,H,o2_c,o3_c,T,sol,sol_bin_width)

if new_numerics==True:
	o2_running=0
	o3_running=0
	o3=np.zeros(nlevs)
	o=np.zeros(nlevs)
        d1defs=np.genfromtxt("species.dat",dtype='str',skiprows=2)
        bimol=np.genfromtxt("bimol.dat",dtype='str',skiprows=2)
        photo=np.genfromtxt("photol.dat",dtype='str',skiprows=2)
        nrxns=len(bimol)+len(photo)
        nspec=len(d1defs[:,0])
	d1=d1_init(d1defs,nlevs,ratio)
	o2=d1[np.where(d1defs=='O2')[0][0],:]
	M=d1[np.where(d1defs=='M')[0][0],:]
        d3=np.zeros([len(d1defs[:,0]),nlevs])
	for i in range(nlevs):
		#Create dictionary of rates
                rates={}
                for a in range(len(bimol)):
                        rates[a]=k(bimol[a,2],T[i],M[i])
                I=od_chems(o2_running,o3_running,o2_c,o3_c,sol)
                for a in range(len(photo)):
                        rates[a+len(bimol)]=j(I,o2_c,o3_c,photo[a,1],sol_bin_width)
		#Define function for chemical tendencies
		def f(y, t):
			g=[rates[0]*y[1]*o2[i]-(rates[3]+rates[1]*y[1])*y[0],2*rates[2]*o2[i]+rates[3]*y[0]-(rates[0]*o2[i]+rates[1]*y[0])*y[1]]
			return g
		#initial conditions
		y0=[d1[np.where(d1defs=='O3')[0][0],i],d1[np.where(d1defs=='O')[0][0],i]]
		time_tot=8640
		t=np.linspace(0,time_tot,time_tot/3600)
		#solution provides a time series, take the final value
		soln=odeint(f,y0,t)
		soln=soln[len(t)-1,:]
		if i==0:
			print soln
		d3[np.where(d1defs=='O3')[0][0],i]=soln[0]
                d3[np.where(d1defs=='O')[0][0],i]=soln[1]
		o3_running=o3_running+d3[np.where(d1defs=='O3')[0][0],i]*box_h
		o2_running=o2_running+d1[np.where(d1defs=='O2')[0][0],i]*box_h
	#print d3[np.where(d1defs=='O3')[0][0],:]
	#print o3_init
	plt.plot(heights,(d1[np.where(d1defs=='O3')[0][0],:]-d3[np.where(d1defs=='O3')[0][0],:])/d1[np.where(d1defs=='O3')[0][0],:])
	plt.savefig('new_numerics.png')
