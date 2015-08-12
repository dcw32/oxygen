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
from optical_depth import od_chems_new
from j_rates import j
from k_rates import k
from steady_state import steady,ozone,otp
from d_arrays import d2_calc_new,d3_calc_new,d2_calc,d3_calc
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
## Chemical scheme: ChapmanSS, Chapman
#chem_scheme = 'Chapman'

#Constants
#Number of atmospheric levels
nlevs=161
#Height, km
h_min=0
h_max=80
heights=np.linspace(h_max,h_min,nlevs)
#scale height, km
H=7
#surface M, cm-3
M_surf=2.5E19*(ratio+0.79)

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
        file=Dataset('netcdf/O3.nc')
        o3_init=file.variables['O3'][:]
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
			g=[k3M*y[1]*o2[i]-(J3+k4*y[1])*y[0],2*J2*o2[i]+J3*y[0]-(k3M*o2[i]+k4*y[0])*y[1]]
			return g
		#initial conditions
		o0=o_init[i]
		o30=o3_init[i]
		y0=[o30,o0]
		t=np.linspace(0,8640000,2400)
		#There's probably a more efficient way to do this but meh
		soln=odeint(f,y0,t)
		soln=soln[len(t)-1,:]
		o3[i]=soln[0]
		o[i]=soln[1]
		o3_running=o3_running+o3[i]*(1E5*(h_max-h_min)/(nlevs-1))
		o2_running=o2_running+o2[i]*(1E5*(h_max-h_min)/(nlevs-1))
	plt.plot(heights,(o3_init-o3)/o3_init)
	plt.savefig('new_numerics.png')
if numerics==True:
#The D1 array contains the chemical species concentrations
#D1 metadata is contained in D1DEFS

	print "CONSTRUCTING D1 ARRAY"
	d1defs=np.genfromtxt("species.dat",dtype='str',skiprows=2)
	#d1=d1_init(d1defs,nlevs,ratio)
	d1={}
	y=0
	while y<len(d1defs[:,0]):
	        D=Dataset('netcdf/'+d1defs[y,0]+'.nc')
	        value=D.variables[d1defs[y,0]][:]
	        key=d1defs[y,0]
	        d1[key]=value
	        y+=1
	print d1.__class__
	bimol=np.genfromtxt("bimol.dat",dtype='str',skiprows=2)
	photo=np.genfromtxt("photol.dat",dtype='str',skiprows=2)
	nrxns=len(bimol)+len(photo)
	nspec=len(d1defs[:,0])
#The D2 array contains the chemical tendencies for each reaction
#There is a chemical tendency arising from each chemical reaction
	for i in range(1):
		rates={}
		for a in range(len(bimol)):
			rates[a]=k(bimol[a,2],T[i],M[i])
		I=od_chems_new(d1,d1defs,o2_c,o3_c,sol,nlevs,len(sol),h_max,h_min)
		for a in range(len(photo)):
			rates[a+len(bimol)]=j(I,o2_c,o3_c,photo[a,1],sol_bin_width)
		print rates
		def f(y, t):
			for a in range(nspec):
				sp=d1defs[a,0]
				d3[sp]=y[a]
				d2[sp]=0
				for b in range(np.where(bimol==d1defs[a,0])[:][0].shape[0]):
					rxno=np.where(bimol==d1defs[i,0])[0][j]
					if np.where(bimol==d1defs[a,0])[1][b]<2:
						d2[sp]=d2[sp]-d3[bimol[rxno,0]]*d3[bimol[rxno,1]]*rates[bimol[rxno,2]]
					else:
						d2[sp]=d2[sp]+d3[bimol[rxno,0]]*d3[bimol[rxno,1]]*rates[bimol[rxno,2]]
			return d2
		t=[0,3600,7200]
		soln=ode(f).set_integrator('lsoda')
		print soln
		soln.integrate(soln,3600)
#		soln=odeint(f,d1,t)
		#d2=d2_calc_new(nrxns,nlevs,d1,bimol,T,photo,o2_c,o3_c,sol,sol_bin_width,d1defs,d1['M'],h_max,h_min)
		#d3=d3_calc_new(d1defs,nlevs,bimol,photo,d2)
		#print d3[1,:]
		#for i in range(len(d1defs[:,0])):
		#	d1[d1defs[i,0]]=d1[d1defs[i,0]]+d3[d1defs[i,0]]
		#doz=d1['O3']
		#doz=doz*(1E5*(h_max-h_min)/(nlevs-1))
		#doz=np.sum(doz)
		#doz=doz/2.69E16
du=o3_running/2.69E16
#plt.plot(np.linspace(h_min,h_max,nlevs),d1[1,:])
#plt.show()
