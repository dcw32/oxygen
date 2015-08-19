import numpy as np
from scipy.integrate import odeint
from netCDF4 import Dataset

from k_rates import k
from j_rates import j

def solve(constants):
 o2_running=0
 o3_running=0
 d1defs=np.genfromtxt("species.dat",dtype='str',skiprows=2)
 bimol=np.genfromtxt("bimol.dat",dtype='str',skiprows=2)
 photo=np.genfromtxt("photol.dat",dtype='str',skiprows=2)
 nspec=len(d1defs[:,0])
 d1=np.zeros([len(d1defs[:,0]),constants['nlevs']])
 for i in range(len(d1defs[:,0])):
  file=Dataset('netcdf/'+d1defs[i,0]+'.nc')
  species=file.variables[d1defs[i,0]][:]
  d1[i,:]=species
  file.close()
 #M=d1[np.where(d1defs=='M')[0][0],:]
 d3=np.zeros([nspec,constants['nlevs']])
 for i in range(constants['nlevs']):
# Create dictionary of rates
  rates={}
  for a in range(len(bimol)):
   rates[bimol[a,2]]=k(bimol[a,2],constants,i)
  optical_depth=constants['o2_c']*o2_running+constants['o3_c']*o3_running
  I_factor=np.exp(-optical_depth)
  I=np.array(I_factor*constants['sol'])
  #I=od_chems(o2_running,o3_running,constants)
  for a in range(len(photo)):
   rates[photo[a,1]]=j(I,photo[a,1],constants)
# Define function for chemical tendencies
  def f(y, t):
   g=np.empty(nspec)
# Loop through chemical species
   for spec in range(nspec):
# Set to zero if species constant, e.g. M, O2
    if d1defs[spec,1]=='n':
     g[spec]=0.0
    else:
     g[spec]=0.0
# Add chemical tendency for a bimolecular rxn
# np.where(bimol==d1defs[spec,0])[:][0].shape[0] is the number of times spec appears in bimol   
     for bi in range(np.where(bimol==d1defs[spec,0])[:][0].shape[0]):
      rxno=np.where(bimol==d1defs[spec,0])[0][bi]
      if np.where(bimol==d1defs[spec,0])[1][bi]<2:
# If appears on LHS of eqn, results in loss
       g[spec]-=y[np.where(d1defs==bimol[rxno,0])[0][0]]\
                *y[np.where(d1defs==bimol[rxno,1])[0][0]]\
                *rates[bimol[rxno,2]]
      else:
# If appears on RHS of eqn, results in production
       g[spec]+=y[np.where(d1defs==bimol[rxno,0])[0][0]]\
                *y[np.where(d1defs==bimol[rxno,1])[0][0]]\
                *rates[bimol[rxno,2]]
     for pho in range(np.where(photo==d1defs[spec,0])[:][0].shape[0]):
      rxno=np.where(photo==d1defs[spec,0])[0][pho]
      if np.where(photo==d1defs[spec,0])[1][pho]<1:
       g[spec]-=y[np.where(d1defs==photo[rxno,0])[0][0]]\
                *rates[photo[rxno,1]]
      else:
       g[spec]+=y[np.where(d1defs==photo[rxno,0])[0][0]]\
                *rates[photo[rxno,1]]
   return g
# Set up times
  time_tot=3.1104E9
  t=np.linspace(0,time_tot,time_tot/3600)
#solution provides a time series, take the final value
  soln=odeint(f,d1[:,i],t)
  soln=soln[len(t)-1,:]
  d3[:,i]=soln
  o3_running=o3_running+d3[np.where(d1defs=='O3')[0][0],i]*constants['box_h']
  o2_running=o2_running+d1[np.where(d1defs=='O2')[0][0],i]*constants['box_h']
 return d1defs,d1,d3
