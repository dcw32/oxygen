if numerics==True:
#The D1 array contains the chemical species concentrations
#D1 metadata is contained in D1DEFS

        print "CONSTRUCTING D1 ARRAY"
        d1defs=np.genfromtxt("species.dat",dtype='str',skiprows=2)
        #d1=d1_init(d1defs,nlevs,ratio)
        d1={}
        y=0
        while y<len(d1defs[:,0]):
                D=Dataset('netcdf/'+d1defs[y,0]+'.nc')
                value=D.variables[d1defs[y,0]][:]
                key=d1defs[y,0]
                d1[key]=value
                y+=1
        print d1.__class__
        bimol=np.genfromtxt("bimol.dat",dtype='str',skiprows=2)
        photo=np.genfromtxt("photol.dat",dtype='str',skiprows=2)
        nrxns=len(bimol)+len(photo)
        nspec=len(d1defs[:,0])
#The D2 array contains the chemical tendencies for each reaction
#There is a chemical tendency arising from each chemical reaction
        for i in range(1):
                rates={}
                for a in range(len(bimol)):
                        rates[a]=k(bimol[a,2],T[i],M[i])
                I=od_chems_new(d1,d1defs,o2_c,o3_c,sol,nlevs,len(sol),h_max,h_min)
                for a in range(len(photo)):
                        rates[a+len(bimol)]=j(I,o2_c,o3_c,photo[a,1],sol_bin_width)
                print rates
                def f(y, t):
                        for a in range(nspec):
                                sp=d1defs[a,0]
                                d3[sp]=y[a]
                                d2[sp]=0
                                for b in range(np.where(bimol==d1defs[a,0])[:][0].shape[0]):
                                        rxno=np.where(bimol==d1defs[i,0])[0][j]
                                        if np.where(bimol==d1defs[a,0])[1][b]<2:
                                                d2[sp]=d2[sp]-d3[bimol[rxno,0]]*d3[bimol[rxno,1]]*rates[bimol[rxno,2]]
                                        else:
                                                d2[sp]=d2[sp]+d3[bimol[rxno,0]]*d3[bimol[rxno,1]]*rates[bimol[rxno,2]]
                        return d2
                t=[0,3600,7200]
                soln=ode(f).set_integrator('lsoda')
                print soln
                soln.integrate(soln,3600)
#               soln=odeint(f,d1,t)
                #d2=d2_calc_new(nrxns,nlevs,d1,bimol,T,photo,o2_c,o3_c,sol,sol_bin_width,d1defs,d1['M'],h_max,h_min)
                #d3=d3_calc_new(d1defs,nlevs,bimol,photo,d2)
                #print d3[1,:]
                #for i in range(len(d1defs[:,0])):
                #       d1[d1defs[i,0]]=d1[d1defs[i,0]]+d3[d1defs[i,0]]
                #doz=d1['O3']
                #doz=doz*(1E5*(h_max-h_min)/(nlevs-1))
                #doz=np.sum(doz)
                #doz=doz/2.69E16
du=o3_running/2.69E16
#plt.plot(np.linspace(h_min,h_max,nlevs),d1[1,:])
#plt.show()

