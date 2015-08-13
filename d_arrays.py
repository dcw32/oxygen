# These functions construct the D arrays
# D1 holds the concentrations of each chemical species in species.dat
# This is initialised from values calculated in the steady state,
# or from a given initial concentration, all in netCDF files
# D2 holds the array of rates calculated for each rxn
# D3 holds the array of chemical tendencies d[species]/dt, for each species.
import numpy as np
from netCDF4 import Dataset
import os, sys
from k_rates import k
from j_rates import j
from optical_depth import od_chems,od_chems_new
def d2_calc_new(nrxns,nlevs,d1,bimol,T,photo,o2_c,o3_c,sol,sol_bin_width,d1defs,M,h_max,h_min):
#        print "CONSTRUCTING D2 ARRAY OF CHEMICAL TENDENCIES"
        #d2=np.zeros([nrxns,nlevs])
	d2={}
        for i in range(len(bimol)):
 #               print "BIMOLECULAR RXN"
                spec_1=d1[bimol[i,0]]
                spec_2=d1[bimol[i,1]]
                k_rt=k(bimol[i,2],T,M)
                spec=np.multiply(spec_1,spec_2)
                d2[i]=np.multiply(spec,k_rt)
        for i in range(len(photo)):
  #              print "PHOTOLYTIC RXN"
                spec=d1[photo[i,0]]
                j_rt=np.zeros(nlevs)
                I=od_chems_new(d1,d1defs,o2_c,o3_c,sol,nlevs,len(sol),h_max,h_min)
                j_rt=j(I,o2_c,o3_c,photo[i,1],sol_bin_width)
                d2[i+len(bimol)]=np.multiply(spec,j_rt)
        return d2
def d3_calc_new(d1defs,nlevs,bimol,photo,d2):
#        d3=np.zeros([len(d1defs[:,0]),nlevs])
	d3={}
	for i in range(len(d1defs[:,0])):
		#i.e. which species are we dealing with?
		key=d1defs[i,0]
		d3[key]=np.zeros(nlevs)
		print key
		try:
			for j in range(np.where(bimol==d1defs[i,0])[:][0].shape[0]):
				if np.where(bimol==d1defs[i,0])[1][j]<2:
					d3[key]=d3[key]-d2[np.where(bimol==d1defs[i,0])[0][j],:]
				else:
					d3[key]=d3[key]+d2[np.where(bimol==d1defs[i,0])[0][j],:]
		except:
			print "WARNING: SPECIES "+d1defs[i,0]+" NOT IN BIMOL"
		try:
			for j in range(np.where(photo==d1defs[i,0])[:][0].shape[0]):
				if np.where(photo==d1defs[i,0])[1][j]<1:
					d3[key]=d3[key]-d2[len(bimol)+np.where(photo==d1defs[i,0])[0][j],:]
				else:
					d3[key]=d3[key]+d2[len(bimol)+np.where(photo==d1defs[i,0])[0][j],:]
		except:
			print "WARNING: SPECIES "+d1defs[i,0]+" NOT IN PHOTOL"
	return d3
# Old code
# might be new code..
def d2_calc(nrxns,nlevs,d1,bimol,T,photo,o2_c,o3_c,sol,sol_bin_width,d1defs,M,h_max,h_min):
#        print "CONSTRUCTING D2 ARRAY OF CHEMICAL TENDENCIES"
        d2=np.zeros([nrxns,nlevs])
        for i in range(len(bimol)):
 #               print "BIMOLECULAR RXN"
                spec_1=d1[np.where(d1defs==bimol[i,0])[0][0],:]
                spec_2=d1[np.where(d1defs==bimol[i,1])[0][0],:]
                k_rt=k(bimol[i,2],T,M)
                spec=np.multiply(spec_1,spec_2)
                d2[i]=np.multiply(spec,k_rt)
        for i in range(len(photo)):
  #              print "PHOTOLYTIC RXN"
                spec=d1[np.where(d1defs==photo[i,0])[0][0],:]
                j_rt=np.zeros(nlevs)
                I=od_chems(d1,d1defs,o2_c,o3_c,sol,nlevs,len(sol),h_max,h_min)
                j_rt=j(I,o2_c,o3_c,photo[i,1],sol_bin_width)
                d2[i+len(bimol)]=np.multiply(spec,j_rt)
	return d2
def d3_calc(d1defs,nlevs,bimol,photo,d2):
        d3=np.zeros([len(d1defs[:,0]),nlevs])
        for i in range(len(d1defs[:,0])):
                        #calc n times spec occurs
         try:
          for j in range(np.where(bimol==d1defs[i,0])[:][0].shape[0]):
           if np.where(bimol==d1defs[i,0])[1][j]<2:
            d3[i,:]=d3[i,:]-d2[np.where(bimol==d1defs[i,0])[0][j],:]
           else:
            d3[i,:]=d3[i,:]+d2[np.where(bimol==d1defs[i,0])[0][j],:]
         except:
          print "WARNING: SPECIES "+d1defs[i,0]+" NOT IN BIMOL"
         try:
          for j in range(np.where(photo==d1defs[i,0])[:][0].shape[0]):
           if np.where(photo==d1defs[i,0])[1][j]<1:
            d3[i,:]=d3[i,:]-d2[len(bimol)+np.where(photo==d1defs[i,0])[0][j],:]
           else:
            d3[i,:]=d3[i,:]+d2[len(bimol)+np.where(photo==d1defs[i,0])[0][j],:]
         except:
          print "WARNING: SPECIES "+d1defs[i,0]+" NOT IN PHOTOL"
	return d3
def d1_init(d1defs,nlevs,ratio):
	d1=np.zeros([len(d1defs[:,0]),nlevs])
	for i in range(len(d1defs[:,0])):
#Perhaps this could be substantially tidied up...
		if d1defs[i,1]=='y':
			file=Dataset('netcdf/'+d1defs[i,0]+'.nc')
			species=file.variables[d1defs[i,0]][:]
			d1[i,:]=species
			file.close()
                elif d1defs[i,0]=='O2':
                        file=Dataset('netcdf/O2.nc')
                        species=file.variables['O2'][:]
                        d1[i,:]=species
                        file.close()
		elif d1defs[i,0]=='M':
			file=Dataset('netcdf/O2.nc')
			species=file.variables['O2'][:]
			species=species/ratio
			d1[i,:]=species
			file.close()
		else:
			#Something's going catastrophically wrong
			print >> sys.stderr, "ERROR: STEADY STATE SPECIES NOT DEFINED IN D1_INIT ARRAY"
			sys.exit(1)
	return d1
def d2_calc(nrxns,nlevs,d1,bimol,T,photo,o2_c,o3_c,sol,sol_bin_width,d1defs,M,h_max,h_min):
#        print "CONSTRUCTING D2 ARRAY OF CHEMICAL TENDENCIES"
        d2=np.zeros([nrxns,nlevs])
        for i in range(len(bimol)):
 #               print "BIMOLECULAR RXN"
                spec_1=d1[np.where(d1defs==bimol[i,0])[0][0],:]
                spec_2=d1[np.where(d1defs==bimol[i,1])[0][0],:]
                k_rt=k(bimol[i,2],T,M)
                spec=np.multiply(spec_1,spec_2)
                d2[i]=np.multiply(spec,k_rt)
        for i in range(len(photo)):
  #              print "PHOTOLYTIC RXN"
                spec=d1[np.where(d1defs==photo[i,0])[0][0],:]
                j_rt=np.zeros(nlevs)
                I=od_chems(d1,d1defs,o2_c,o3_c,sol,nlevs,len(sol),h_max,h_min)
                j_rt=j(I,o2_c,o3_c,photo[i,1],sol_bin_width)
                d2[i+len(bimol)]=np.multiply(spec,j_rt)
	return d2
def d3_calc(d1defs,nlevs,bimol,photo,d2):
        d3=np.zeros([len(d1defs[:,0]),nlevs])
        for i in range(len(d1defs[:,0])):
                        #calc n times spec occurs
         try:
          for j in range(np.where(bimol==d1defs[i,0])[:][0].shape[0]):
           if np.where(bimol==d1defs[i,0])[1][j]<2:
            d3[i,:]=d3[i,:]-d2[np.where(bimol==d1defs[i,0])[0][j],:]
           else:
            d3[i,:]=d3[i,:]+d2[np.where(bimol==d1defs[i,0])[0][j],:]
         except:
          print "WARNING: SPECIES "+d1defs[i,0]+" NOT IN BIMOL"
         try:
          for j in range(np.where(photo==d1defs[i,0])[:][0].shape[0]):
           if np.where(photo==d1defs[i,0])[1][j]<1:
            d3[i,:]=d3[i,:]-d2[len(bimol)+np.where(photo==d1defs[i,0])[0][j],:]
           else:
            d3[i,:]=d3[i,:]+d2[len(bimol)+np.where(photo==d1defs[i,0])[0][j],:]
         except:
          print "WARNING: SPECIES "+d1defs[i,0]+" NOT IN PHOTOL"
	return d3
def d1_init(d1defs,nlevs,ratio):
	d1=np.zeros([len(d1defs[:,0]),nlevs])
	for i in range(len(d1defs[:,0])):
		file=Dataset('netcdf/'+d1defs[i,0]+'.nc')
		species=file.variables[d1defs[i,0]][:]
		d1[i,:]=species
		file.close()
	return d1
