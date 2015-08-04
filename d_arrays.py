#
#
import numpy as np
from k_rates import k
from j_rates import j
from optical_depth import od_chems
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
