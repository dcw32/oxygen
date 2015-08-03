# Calculates ozone profile in steady-state
# 
from j_rates import j
from k_rates import k
from ozone_calc import ozone,otp
from data import ncsave
import numpy as np
import matplotlib.pyplot as plt

def steady(nlevs,h_max,h_min,H,M_surf,ratio,o2_c,o3_c,T,sol,sol_bin_width):
	o3=np.zeros(nlevs)
	o3_sum=np.zeros(nlevs)
	o2=np.zeros(nlevs)
	o2_sum=np.zeros(nlevs)
	o=np.zeros(nlevs)
	J_o2=np.zeros(nlevs)
	J_o3=np.zeros(nlevs)
	k3=np.zeros(nlevs)
	k4=np.zeros(nlevs)
	q=np.zeros(nlevs)
	height=np.linspace(h_max,h_min,nlevs)
	o2_running=0
	o3_running=0
	for i in range(nlevs):
        	o2[i]=ratio*M_surf*np.exp(-height[i]/H)
	        optical_depth=o2_c*o2_running+o3_c*o3_running
	        I_factor=np.exp(-optical_depth)
        	I=np.array(sol*I_factor)
        	J_o2[i]=j(I,o2_c,o3_c,'JO2',sol_bin_width)
        	J_o3[i]=j(I,o2_c,o3_c,'JO3',sol_bin_width)
        	k3[i]=k('k3',T[i],o2[i]/ratio)
        	k4[i]=k('k4',T[i],o2[i]/ratio)
        	o3[i],q[i]=ozone(J_o2[i],J_o3[i],T[i],o2[i],ratio,k3[i],k4[i])
        	o[i]=otp(o3[i],J_o3[i],k3[i],o2[i],ratio)
        	o2_running=o2_running+o2[i]*(1E5*h_max/(nlevs-1))
        	o2_sum[i]=o2_running
        	o3_running=o3_running+o3[i]*(1E5*h_max/(nlevs-1))
        	o3_sum[i]=o3_running
	plt.plot(height,J_o2)
	plt.show()
	ncsave(height,'O3',o3)
	ncsave(height,'O2',o2)
	ncsave(height,'O',o)
	return o3_running
