import numpy as np
from scipy.integrate import odeint
from netCDF4 import Dataset
import sys
import os

from k_rates import k
from j_rates import j

def solve(constants):
 o2_running=0
 o3_running=0
 constants['M_surf']=2.5E19*(constants['ratio']+0.79)
 constants['M']=constants['M_surf']*np.exp(-constants['heights']/constants['H'])
 d1defs=np.genfromtxt("species.dat",dtype='str',skiprows=2)
 bimol=np.genfromtxt("bimol.dat",dtype='str',skiprows=2)
 photo=np.genfromtxt("photol.dat",dtype='str',skiprows=2)
 nspec=len(d1defs[:,0])
 d1=np.zeros([len(d1defs[:,0]),constants['nlevs']])
 file=Dataset('netcdf/initial.nc')
 for i in range(len(d1defs[:,0])):
  #if os.path.isfile('netcdf/'+d1defs[i,0]+'.nc')==True:
  #file=Dataset('netcdf/initial.nc')
  species=file.variables[d1defs[i,0]][:]
  d1[i,:]=species
  #else:
  # sys.exit("SPECIES "+d1defs[i,0]+" IS NOT INITIALISED")
 for i in range(len(bimol[:,0])):
  for a in [0,1,3,4,5,6]:
   if bimol[i,a] in d1defs[:,0] or bimol[i,a]=='X':
    pass
   else:
    sys.exit("SPECIES "+bimol[i,a]+" IN BIMOL RXN BUT NOT IN SPECIES.DAT")
 for i in range(len(photo[:,0])):
  for a in [0,2,3,4,5]:
   if photo[i,a] in d1defs[:,0] or photo[i,a]=='X':
    pass
   else:
    sys.exit("SPECIES "+photo[i,a]+" IN PHOTOL RXN BUT NOT IN SPECIES.DAT")
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
  t=np.linspace(0,5.174E8/60.,5.174E8/3600.)
#solution provides a time series, take the final value
  soln,hi=odeint(f,d1[:,i],t,full_output=1)
  soln=soln[len(t)-1,:]
  d3[:,i]=soln
  o3_running=o3_running+d3[np.where(d1defs=='O3')[0][0],i]*constants['box_h']
  o2_running=o2_running+d1[np.where(d1defs=='O2')[0][0],i]*constants['box_h']
 return d1defs,d1,d3
