# This file is the main control routine for the stratospheric ozone code
# July 2015, dcw32@cam.ac.uk
# run as ipython --pylab
# %run ./main.py
import numpy as np
import matplotlib.pyplot as plt
import sys

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
steadystate = False
numerics = True
# Chemical scheme: ChapmanSS, Chapman
chem_scheme = 'Chapman'

#Constants
#Number of atmospheric levels
nlevs=101
#Height, km
h_min=0
h_max=50
#scale height, km
H=7
#surface M, cm-3
M_surf=2.5E19*(ratio+0.79)


#Temperature
#T=np.zeros(nlevs)
#T=np.linspace(270,180,nlevs)
#T[-6:]=np.linspace(T[-6],290,6)
T=291*np.ones(nlevs)

# wavelength bins, interpolating o2 & solar onto o3 csec wavelengths
wlen,o2_c,o3_c,sol=extract_csec()

#convert sol to photons cm-2 s-1 nm-1 for binning
sol=sol/0.05

sol_bin_width=np.zeros(len(wlen))
sol_bin_width[0]=2
sol_bin_width[len(wlen)-1]=5
for i in range(1,len(wlen)-2):
	sol_bin_width[i]=0.5*(sol_bin_width[i+1]-sol_bin_width[i-1])

o2_c=np.array(o2_c)
o3_c=np.array(o3_c)

if steadystate==True:
	print "CALCULATING STEADY STATE OZONE COLUMN"
	o3_running=steady(nlevs,h_max,h_min,H,M_surf,ratio,o2_c,o3_c,T,sol,sol_bin_width)
	
#plt.semilogy(height,o3)
#print np.divide(J_o2,J_o3)
#print np.power(o2,3)/ratio
if numerics==True:
#The D1 array contains the chemical species concentrations
#D1 metadata is contained in D1DEFS

	print "CONSTRUCTING D1 ARRAY"
	d1defs=np.genfromtxt("species.dat",dtype='str',skiprows=2)
	d1=d1_init(d1defs,nlevs,ratio)
	print d1
	bimol=np.genfromtxt("bimol.dat",dtype='str',skiprows=2)
	photo=np.genfromtxt("photol.dat",dtype='str',skiprows=2)
	nrxns=len(bimol)+len(photo)
#The D2 array contains the chemical tendencies for each reaction
#There is a chemical tendency arising from each chemical reaction
	from d_arrays import d2_calc,d3_calc
	for i in range(1000):
		d2=d2_calc(nrxns,nlevs,d1,bimol,T,photo,o2_c,o3_c,sol,sol_bin_width,d1defs)
		d3=d3_calc(d1defs,nlevs,bimol,photo,d2)
		d1=d1+3600000*d3
		for i in range(len(d1defs[:,0])):
			for j in range(nlevs):
				if d1[i,j]<0:
					d1[i,j]=0
print d1
#du=o3_running/2.69E16
#print str(du)+" DU Ozone Column"
doz=d1[1,:]
doz=doz*(1E5*h_max/(nlevs-1))
doz=np.sum(doz)
doz=doz/2.69E16
print str(doz)
#plt.plot(height,o3)
#plt.show()
