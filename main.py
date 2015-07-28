# This file is the main control routine for the stratospheric ozone code
# July 2015, dcw32@cam.ac.uk
# run as ipython --pylab
# %run ./main.py
import numpy as np
import matplotlib.pyplot as plt

from data import extract_csec
from optical_depth import od
from j_rates import j
from ozone_calc import ozone

#Constants
#Number of atmospheric levels
nlevs=41
#Height, km
h_min=0
h_max=80
#scale height, km
H=7
#surface M, cm-3
M_surf=2.5E19

#Temperature
T=np.zeros(nlevs)
T=np.linspace(260,150,nlevs)
T[-6:]=np.linspace(290,220,6)

# wavelength bins, interpolating o2 & solar onto o3 csec wavelengths
wlen,o2_c,o3_c,sol=extract_csec()

o2_c=np.array(o2_c)
o3_c=np.array(o3_c)

#Initialise arrays
#arrays run 0 to nlevs-1 where 0 is TOA and nlevs-1 is at surface
o3=np.zeros(nlevs)
o2=np.zeros(nlevs)
J_o2=np.zeros(nlevs)
J_o3=np.zeros(nlevs)
o3_running=0
o2_running=0
height=np.linspace(h_max,h_min,nlevs)

#Calculation of j rates, ozone)
for i in range(nlevs):
        o2[i]=0.21*M_surf*np.exp(-height[i]/H)
	I_factor=od(o2_c,o3_c,o2_running,o3_running)
	I=sol*I_factor
	I=np.array(I)
	J_o2_s=np.multiply(o2_c,I)
	J_o2[i]=np.sum(J_o2_s)
	print J_o2[i]
	J_o3_s=np.multiply(o3_c,I)
	J_o3[i]=np.sum(J_o3_s)
	print J_o3[i]
	o3[i]=ozone(J_o2[i],J_o3[i],T[i],o2[i])
	o2_running=o2_running+o2[i]*(1E5*h_max/(nlevs-1))	
	o3_running=o3_running+o3[i]*(1E5*h_max/(nlevs-1))
plt.semilogy(height,o3)
#plt.plot(height,o2*1E-5)
plt.show()
