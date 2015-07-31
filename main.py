# This file is the main control routine for the stratospheric ozone code
# July 2015, dcw32@cam.ac.uk
# run as ipython --pylab
# %run ./main.py
import numpy as np
import matplotlib.pyplot as plt

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
numerics = True
# Chemical scheme: ChapmanSS, Chapman
chem_scheme = 'Chapman'

#Constants
#Number of atmospheric levels
nlevs=1001
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
	o3_running=steady(nlevs,h_max,h_min,H,M_surf,ratio,o2_c,o3_c,T,sol,sol_bin_width)

#plt.semilogy(height,o3)
#print np.divide(J_o2,J_o3)
#print np.power(o2,3)/ratio

if numerics==True:
#The D1 array contains the chemical species concentrations
#D1 metadata is contained in D1DEFS
	d1defs=np.genfromtxt("species.dat",dtype='str',skiprows=2)
	d1=d1_init(d1defs,nlevs,ratio)
	bimol=np.genfromtxt("bimol.dat",dtype='str',skiprows=2)
	photo=np.genfromtxt("photol.dat",dtype='str',skiprows=2)
	nrxns=len(bimol)+len(photo)
#The D2 array contains the chemical tendencies for each reaction
#There is a chemical tendency arising from each chemical reaction
	d2=np.empty([nrxns,nlevs])
	for i in range(len(bimol)):
		spec_1=d1[np.where(d1defs==bimol[i,0])[0][0],:]
		spec_2=d1[np.where(d1defs==bimol[i,1])[0][0],:]
		k_rt=k(bimol[i,2],T)
		spec=np.multiply(spec_1,spec_2)
		d2[i]=np.multiply(spec,k_rt)
	for i in range(len(photo)):
		spec=d1[np.where(d1defs==photo[i,0])[0][0],:]
		j_rt=np.zeros(nlevs)
		I=od_chems(d1,d1defs,o2_c,o3_c,sol,nlevs,len(sol))
		j_rt=j(I,o2_c,o3_c,photo[i,1],sol_bin_width)
		d2[i+len(bimol)]=np.multiply(spec,j_rt)

du=o3_running/2.69E16
print str(du)+" DU Ozone Column"
#plt.plot(height,o3)
#plt.show()
