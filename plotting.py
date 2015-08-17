import matplotlib.pyplot as plt
import numpy as np
from steady_state import steady
from matplotlib.widgets import Slider, Button, RadioButtons

def altconc(constants,o3,o3_running):
        fig0,axis=plt.subplots()
        plt.plot(constants['heights'],o3)
        axis.set_xlabel('Altitude / km')
        axis.set_ylabel('O3 concentration / cm-3')
        fig0.clear
        du=o3_running/2.69E16
        du=int(round(du))
        axis.text(0.95,0.95,str(du)+" DU",horizontalalignment='right',verticalalignment='top',transform=axis.transAxes)
        plt.title('Steady State Chapman - PAL Oxygen')
        plt.savefig('ozone.png')
	plt.close()
def linoxyoz(constants):
	ratios=np.linspace(0.01,1.00,100)
	fig, ax=plt.subplots()
	for i in range(len(ratios)):
	        constants['ratio']=ratios[i]
	        constants['M_surf']=2.5E19*(constants['ratio']+0.79)
                constants['M']=constants['M_surf']*np.exp(-constants['heights']/constants['H'])
	        o3,o2,o,J_o2,J_o3,o3_running=steady(constants,False)
	        du=o3_running
	        du=o3_running/2.69E16
	        ratio=constants['ratio']/0.21
	        plt.scatter(ratio,du,marker="x",color='purple')
	x,y=(1,1),(0,1000)
	plt.plot(x,y,color='k')
	ax.set_xlim([0,5])
	ax.set_ylim([550,800])
	plt.xlabel(r'Oxygen Content (PAL)')
	plt.ylabel(r'Ozone Column / DU')
	plt.savefig('dobson.png')
	plt.close()
def logoxyoz(constants):
	ratios=np.logspace(-6,0,50)
	fig, ax=plt.subplots()
	for i in range(len(ratios)):
                constants['ratio']=ratios[i]
                constants['M_surf']=2.5E19*(constants['ratio']+0.79)
                constants['M']=constants['M_surf']*np.exp(-constants['heights']/constants['H'])
	        o3,o2,o,J_o2,J_o3,o3_running=steady(constants,False)
	        du=o3_running
	#       du=o3_running/2.69E16
	        ratio=constants['ratio']/0.21
	        plt.scatter(ratio,du,marker="x",color='r')
	ax.set_xscale('log')
	#ax.set_xlim([0,1])
	ax.set_yscale('log')
	plt.xlabel(r'Oxygen Content (PAL)')
	plt.ylabel(r'Ozone Column / cm-2')
	#plt.show()
	plt.savefig('dobson_log.png')
	plt.close()
def inter_plot(o3,J_o2,J_o3,constants):
        fig, (ax2,ax1,ax0) = plt.subplots(3, sharex=True)
        #fig1=plt.subplots()
        plt.subplots_adjust(bottom=0.25)
        ax0.set_xlabel('Height / km')
        ax0.set_ylabel('Ozone / 10^-12 cm-3')
        ax0.set_ylim([0,12])
        ax1.set_ylabel('JO2 / 10^10 s-1')
        ax2.set_ylabel('JO3 / 10^3 s-1')
        l,=ax0.plot(constants['heights'],o3/1E12)
        m,=ax1.plot(constants['heights'],J_o2*1E10,color='red')
        n,=ax2.plot(constants['heights'],J_o3*1E3,color='purple')
        #ann=txt.Annotation(str(int(round(du))), xy=(50,10),xycoords='data')            
        #ax0.add_artist(ann)
        #Sets up sliders
        axo2=plt.axes([0.25,0.1,0.65,0.03],axisbg='lightgoldenrodyellow')
        so2=Slider(axo2, 'O\mathrm{_{2}}', 0.001, 1.00, valinit=0.21)
        #Updates graphs based on slider value
        def update(val):
		ratio=so2.val
                constants['M_surf']=2.5E19*(ratio+0.79)
		constants['M']=constants['M_surf']*np.exp(-constants['heights']/constants['H'])
                o3,o2,o,J_o2,J_o3,o3_running=steady(constants,False)
                m.set_ydata(J_o2*1E10)
                l.set_ydata(o3/1E12)
                n.set_ydata(J_o3*1E3)
                du=o3_running/2.69E16
                #ax0.annotate(str(int(round(du))), xy=(0.99,0.01),xycoords='axes fraction',horizontalalignment='right', verticalalignment='bottom')
                fig.canvas.draw_idle()
        so2.on_changed(update)
        plt.show()

