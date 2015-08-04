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

from data import extract_csec,ncsave
from optical_depth import od_chems
from j_rates import j
from k_rates import k
from ozone_calc import ozone,otp
from chemical_solver import d1_init
from steady_state import steady

#Running options
# oxygen ratio
ratio = 0.21
# debug mode
debug = True
# Steady state mode - can set up 'initial values' for chem scheme
steadystate = True
interactive = False
numerics = False
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

	fig0,axis=plt.subplots()
	plt.plot(height,o3)
	axis.set_xlabel('Altitude / km')
	axis.set_ylabel('O3 concentration / cm-3')
	plt.savefig('ozone.png')
	fig0.clear
	du=o3_running/2.69E16
	print du

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

fig, ax=plt.subplots()
for i in range(100):
	ratio=0.01+0.01*i
	M_surf=2.5E19*(ratio+0.79)
	height,o3,o2,o,J_o2,J_o3,o3_running=steady(nlevs,h_max,h_min,H,M_surf,ratio,o2_c,o3_c,T,sol,sol_bin_width)
	du=o3_running/2.69E16
	plt.scatter(ratio,du,marker=".")
ax.set_xlim([0,1])
plt.xlabel(r'O2 ratio')
plt.ylabel(r'Ozone Column / DU')
#plt.show()
plt.savefig('dobson.png')

#height=np.linspace(h_max,h_min,nlevs)
#o=np.zeros(nlevs)
#ncsave(height,'O',o)
#plt.semilogy(height,o3)
#print np.divide(J_o2,J_o3)
#print np.power(o2,3)/ratio
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
	from d_arrays import d2_calc,d3_calc
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
