
from data import savetonetcdf
from j_rates import j
from k_rates import k
from steady_state import steady,ozone,otp
from plotting import *
from solver import solve
# Driver for the model, calls the required routines for each feature
def steady_state_chapman(constants):
        print "CALCULATING STEADY STATE OZONE COLUMN (CHAPMAN)"
        o3,o2,o,J_o2,J_o3,o3_running=steady(constants,True)
        altconc(constants,o3,o3_running)
def numerical(constants):
                print "RUNNING NUMERICAL SOLVER"
                d1defs,d1,d3=solve(constants)
                savetonetcdf(constants,d1defs,d1,d3)
                
                o3=d1[np.where(d1defs=='O3')[0][0],:]
                plt.plot(constants['heights'],d3[np.where(d1defs=='O3')[0][0],:])
                plt.plot(constants['heights'],d1[np.where(d1defs=='O3')[0][0],:])
                #plt.savefig('ozone_out.png')
                plt.show()
                o3=o3*constants['box_h']
                du=np.sum(o3)/2.69E16
                print du
                o3=d3[np.where(d1defs=='O3')[0][0],:]
                o3=o3*constants['box_h']
                du=np.sum(o3)/2.69E16
                print du
                #print d3[np.where(d1defs=='H2O')[0][0],:]

