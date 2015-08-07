# Calculates ozone profile in steady-state
# 
from j_rates import j
from k_rates import k
from data import ncsave
import numpy as np
import matplotlib.pyplot as plt
import sys

def steady(nlevs,h_max,h_min,H,M_surf,ratio,o2_c,o3_c,T,sol,sol_bin_width):
	o3=np.zeros(nlevs)
	o3_sum=np.zeros(nlevs)
	o2=np.zeros(nlevs)
	o2_sum=np.zeros(nlevs)
	o=np.zeros(nlevs)
	J_o2=np.zeros(nlevs)
	J_o3=np.zeros(nlevs)
	k3M=np.zeros(nlevs)
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
        	k3M[i]=k('k3M',T[i],o2[i]/ratio)
        	k4[i]=k('k4',T[i],o2[i]/ratio)
        	o3[i],q[i]=ozone(J_o2[i],J_o3[i],T[i],o2[i],ratio,k3M[i],k4[i])
        	o[i]=otp(o3[i],J_o3[i],k3M[i],o2[i],ratio)
        	o2_running=o2_running+o2[i]*(1E5*(h_max-h_min)/(nlevs-1))
        	o2_sum[i]=o2_running
        	o3_running=o3_running+o3[i]*(1E5*(h_max-h_min)/(nlevs-1))
        	o3_sum[i]=o3_running
        ncsave(height,'O3',o3)
        ncsave(height,'O',o)
        ncsave(height,'O2',o2)
        M=o2/ratio
        ncsave(height,'M',M)
	#plt.plot(height,o3)
	#plt.savefig('ozone.png')
	return height,o3,o2,o,J_o2,J_o3,o3_running
# This routine calculates the ozone at the model level, given an overhead ozone column, photolysis rates and rate constants
def ozone(J_o2,J_o3,T,o2,ratio,k3M,k4):
	#[O3]ss=(J(O2)*k3*[M]*[O2]^2/J(O3)*k4)^0.5
	q=(J_o2/J_o3)*(k3M/k4)*np.power(o2,2)
	if q<0:
		print >> sys.stderr, "NEGATIVE Q VALUE"
		sys.exit(1)
	o3=np.sqrt(q)
	return o3,q
def otp(o3,J_o3,k3M,o2,ratio):
	#[O]ss=[O3]*J(O3)/k3*[M]*[O2]
	ot=o3*J_o3/(k3M*o2)
	return ot
