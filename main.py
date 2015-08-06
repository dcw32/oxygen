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
from scipy.integrate import odeint

from data import extract_csec,ncsave
from optical_depth import od_chems
from j_rates import j
from k_rates import k
from steady_state import steady,ozone,otp
from d_arrays import d1_init,d2_calc,d3_calc
from plotting import altconc,linoxyoz,logoxyoz

#Running options
# oxygen ratio
ratio = 0.21
# Steady state mode - can set up 'initial values' for chem scheme
steadystate = True
interactive = False
numerics = False
new_numerics = True
# Plotting scripts - plots oxygen vs ozone for a range of oxygen vals
plot_on = True
# Chemical scheme: ChapmanSS, Chapman
chem_scheme = 'Chapman'

#Constants
#Number of atmospheric levels
nlevs=41
#Height, km
h_min=0
h_max=50
heights=np.linspace(h_max,h_min,nlevs)
#scale height, km
H=7
#surface M, cm-3
M_surf=2.5E19*(ratio+0.79)

# wavelength bins, interpolating o2 & solar onto o3 csec wavelengths, extract temps
wlen,o2_c,o3_c,sol,T=extract_csec(heights)

#convert sol to photons cm-2 s-1 nm-1 for binning
sol=sol/0.05

sol_bin_width=np.zeros(len(wlen))
sol_bin_width[0]=2
sol_bin_width[len(wlen)-1]=5
for i in range(1,len(wlen)-2):
        sol_bin_width[i]=0.5*(wlen[i+1]-wlen[i-1])

o2_c=np.array(o2_c)
o3_c=np.array(o3_c)

if steadystate==True:
	print "CALCULATING STEADY STATE OZONE COLUMN"
	height,o3,o2,o,J_o2,J_o3,o3_running=steady(nlevs,h_max,h_min,H,M_surf,ratio,o2_c,o3_c,T,sol,sol_bin_width)
	ncsave(heights,'O3',o3)
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
			height,o3,o2,o,J_o2,J_o3,o3_running=steady(nlevs,h_max,h_min,H,M_surf,ratio,o2_c,o3_c,T,sol,sol_bin_width)
			m.set_ydata(J_o2*1E10)
			l.set_ydata(o3/1E12)
			n.set_ydata(J_o3*1E3)
			du=o3_running/2.69E16
			print du
	                #ax0.annotate(str(int(round(du))), xy=(0.99,0.01),xycoords='axes fraction',horizontalalignment='right', verticalalignment='bottom')
			fig.canvas.draw_idle()
		so2.on_changed(update)
		plt.show()
else:
	if interactive==True:
		print "Must run Steady State calc for Interactive Plotting routine"

if plot_on==True:
	logoxyoz(nlevs,h_max,h_min,H,o2_c,o3_c,T,sol,sol_bin_width)
	linoxyoz(nlevs,h_max,h_min,H,o2_c,o3_c,T,sol,sol_bin_width)

if new_numerics==True:
        file=Dataset('netcdf/O3.nc')
        o3_init=file.variables['O3'][:]
	print o3_init
        file=Dataset('netcdf/O2.nc')
        o2=file.variables['O2'][:]
	M=o2/ratio
        file=Dataset('netcdf/O.nc')
        o_init=file.variables['O'][:]
	o2_running=0
	o3_running=0
	o3=np.zeros(nlevs)
	o=np.zeros(nlevs)
	for i in range(nlevs):
		k3M=k('k3M',T[i],M[i])
		k4=k('k4',T[i],M[i])
		optical_depth=o2_c*o2_running+o3_c*o3_running
		I_factor=np.exp(-optical_depth)
		I=np.array(sol*I_factor)
		J2=j(I,o2_c,o3_c,'JO2',sol_bin_width)
		J3=j(I,o2_c,o3_c,'JO3',sol_bin_width)
		def f(y, t):
			o3i=y[0]
			oi=y[1]
			f0=k3M*oi*o2[i]-(J3+k4*oi)*o3i
			f1=2*J2*o2[i]+J3*o3i-(k3M*o2[i]+k4*o3i)*oi
			return [f0, f1]
		#initial conditions
		o0=o_init[i]
		o30=o3_init[i]
		y0=[o30,o0]
		t=np.linspace(0,3.1104E8,3600)
		soln=odeint(f,y0,t)
		soln=soln[len(t)-1,:]
		o3[i]=soln[0]
		o[i]=soln[1]
		o3_running=o3_running+o3[i]*(1E5*(h_max-h_min)/(nlevs-1))
		o2_running=o2_running+o2[i]*(1E5*(h_max-h_min)/(nlevs-1))
	print o3
	plt.plot(heights,o3_init)
	plt.plot(heights,o3)
	plt.savefig('new_numerics.png')
if numerics==True:
#The D1 array contains the chemical species concentrations
#D1 metadata is contained in D1DEFS

	print "CONSTRUCTING D1 ARRAY"
	d1defs=np.genfromtxt("species.dat",dtype='str',skiprows=2)
	d1=d1_init(d1defs,nlevs,ratio)
	bimol=np.genfromtxt("bimol.dat",dtype='str',skiprows=2)
	photo=np.genfromtxt("photol.dat",dtype='str',skiprows=2)
	nrxns=len(bimol)+len(photo)
#The D2 array contains the chemical tendencies for each reaction
#There is a chemical tendency arising from each chemical reaction
	for i in range(100):
		d2=d2_calc(nrxns,nlevs,d1,bimol,T,photo,o2_c,o3_c,sol,sol_bin_width,d1defs,d1[np.where(d1defs=='M')[0][0],:],h_max,h_min)
		d3=d3_calc(d1defs,nlevs,bimol,photo,d2)
		#print d3[1,:]
		for i in range(len(d1defs[:,0])):
			for j in range(nlevs):
				d1[i,j]=d1[i,j]+d3[i,j]
				if d1[i,j]<0:
					d1[i,j]=0
		print d1
		doz=d1[1,:]
		doz=doz*(1E5*(h_max-h_min)/(nlevs-1))
		doz=np.sum(doz)
		doz=doz/2.69E16
		print doz
du=o3_running/2.69E16
print str(du)+" DU Ozone Column"
#plt.plot(np.linspace(h_min,h_max,nlevs),d1[1,:])
#plt.show()
