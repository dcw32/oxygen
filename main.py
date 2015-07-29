# This file is the main control routine for the stratospheric ozone code
# July 2015, dcw32@cam.ac.uk
# run as ipython --pylab
# %run ./main.py
import numpy as np
import matplotlib.pyplot as plt

from data import extract_csec
from optical_depth import od
from j_rates import j
from k_rates import k3_i,k4_i
from ozone_calc import ozone

#Running options
# oxygen ratio
ratio = 0.20
# debug mode
debug = True

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
T=np.zeros(nlevs)
T=np.linspace(270,180,nlevs)
T[-6:]=np.linspace(T[-6],290,6)

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

#Initialise arrays
#arrays run 0 to nlevs-1 where 0 is TOA and nlevs-1 is at surface
init=np.zeros(nlevs)
o3=np.zeros(nlevs)
o3_sum=np.zeros(nlevs)
o2=np.zeros(nlevs)
o2_sum=np.zeros(nlevs)
J_o2=np.zeros(nlevs)
J_o3=np.zeros(nlevs)
k3=np.zeros(nlevs)
k4=np.zeros(nlevs)
q=np.zeros(nlevs)
height=np.linspace(h_max,h_min,nlevs)

o2_running=0
o3_running=0

#Calculation of j rates, ozone)
for i in range(nlevs):
        o2[i]=ratio*M_surf*np.exp(-height[i]/H)
	I_factor=od(o2_c,o3_c,o2_running,o3_running)
	I=np.array(sol*I_factor)
	J_o2[i]=j(o2_c,I,sol_bin_width)
	J_o3[i]=j(o3_c,I,sol_bin_width)
	k3[i]=k3_i(T[i])
	k4[i]=k4_i(T[i])
	o3[i],q[i]=ozone(J_o2[i],J_o3[i],T[i],o2[i],ratio,k3[i],k4[i])
	o2_running=o2_running+o2[i]*(1E5*h_max/(nlevs-1))	
	o2_sum[i]=o2_running
	o3_running=o3_running+o3[i]*(1E5*h_max/(nlevs-1))
	o3_sum[i]=o3_running
#plt.semilogy(height,o3)
#print np.divide(J_o2,J_o3)
#print np.power(o2,3)/ratio
print o3_running
plt.plot(height,np.divide(J_o2,J_o3))
plt.show()
