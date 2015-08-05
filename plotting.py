import matplotlib.pyplot as plt
import numpy as np
from steady_state import steady
def altconc(height,o3,o3_running):
        fig0,axis=plt.subplots()
        plt.plot(height,o3)
        axis.set_xlabel('Altitude / km')
        axis.set_ylabel('O3 concentration / cm-3')
        fig0.clear
        du=o3_running/2.69E16
        du=int(round(du))
        axis.text(0.95,0.95,str(du)+" DU",horizontalalignment='right',verticalalignment='top',transform=axis.transAxes)
        plt.title('Steady State Chapman - PAL Oxygen')
        plt.savefig('ozone.png')
	plt.close()
def linoxyoz(nlevs,h_max,h_min,H,o2_c,o3_c,T,sol,sol_bin_width):
	ratios=np.linspace(0.01,1.00,100)
	fig, ax=plt.subplots()
	for i in range(len(ratios)):
	        ratio=ratios[i]
	        M_surf=2.5E19*(ratio+0.79)
	        height,o3,o2,o,J_o2,J_o3,o3_running=steady(nlevs,h_max,h_min,H,M_surf,ratio,o2_c,o3_c,T,sol,sol_bin_width)
	        du=o3_running
	        du=o3_running/2.69E16
	        ratio=ratio/0.21
	        plt.scatter(ratio,du,marker="x",color='purple')
	x,y=(1,1),(0,1000)
	plt.plot(x,y,color='k')
	ax.set_xlim([0,5])
	ax.set_ylim([550,800])
	plt.xlabel(r'Oxygen Content (PAL)')
	plt.ylabel(r'Ozone Column / DU')
	plt.savefig('dobson.png')
	plt.close()
def logoxyoz(nlevs,h_max,h_min,H,o2_c,o3_c,T,sol,sol_bin_width):
	ratios=np.logspace(-6,0,50)
	fig, ax=plt.subplots()
	for i in range(len(ratios)):
	        ratio=ratios[i]
	        M_surf=2.5E19*(ratio+0.79)
	        height,o3,o2,o,J_o2,J_o3,o3_running=steady(nlevs,h_max,h_min,H,M_surf,ratio,o2_c,o3_c,T,sol,sol_bin_width)
	        du=o3_running
	#       du=o3_running/2.69E16
	        ratio=ratio/0.21
	        plt.scatter(ratio,du,marker="x",color='r')
	ax.set_xscale('log')
	#ax.set_xlim([0,1])
	ax.set_yscale('log')
	plt.xlabel(r'Oxygen Content (PAL)')
	plt.ylabel(r'Ozone Column / cm-2')
	#plt.show()
	plt.savefig('dobson_log.png')
	plt.close()

