import numpy as np
import sys
def termol(k0,n,kinf,m,constants,i):
	k0T=k0*np.power((constants['T'][i]/300),-n)
	kinfT=kinf*np.power((constants['T'][i]/300),-m)
	coeff=1/(1+np.power(np.log10(k0T*constants['M'][i]/kinfT),2))
	k_termol=k0*constants['M'][i]*np.power(0.6,coeff)/(1+(k0T*constants['M'][i])/kinfT)
	return k_termol
def bimol(A,E,constants,i):
	k_bimol=A*np.exp(-E/constants['T'][i])
	return k_bimol
def k(krt,constants,i):
#to add new rate:
#elif krt=='RATE':
#	k=...
	#k3 is really k3*[M] in order to make the equation quasi bimolecular
	if krt=='k2':
		k=2.15E-11*np.exp(110/constants['T'][i])
	elif krt=='k3M':
	        k=6E-34*(300/constants['T'][i])*(300/constants['T'][i])*constants['M'][i]
	elif krt=='k4':
	        k=1E-11*np.exp(-2100/constants['T'][i])
	elif krt=='k5':
		k=3E-12*np.exp(-1500/constants['T'][i])
	elif krt=='k6':
		k=5.6E-12*np.exp(180/constants['T'][i])
	elif krt=='k7':
		k=1.63E-10*np.exp(60/constants['T'][i])
	elif krt=='k8':
		k=1.7E-12*np.exp(-940/constants['T'][i])
	elif krt=='k9':
		k=3.0E-11*np.exp(200/constants['T'][i])
	elif krt=='k10':
		k=3.3E-12*np.exp(270/constants['T'][i])
	elif krt=='k11':
		k=1E-14*np.exp(-490/constants['T'][i])
	elif krt=='k12':
		k=3E-13*np.exp(460/constants['T'][i])
	elif krt=='k13':
		k=1.8E-12
	elif krt=='k14':
		k=2.3E-11*np.exp(-200/constants['T'][i])
	elif krt=='k15':
		k=2.8E-11*np.exp(85/constants['T'][i])
	elif krt=='k16':
		k=6.4E-12*np.exp(290/constants['T'][i])
	elif krt=='k17M':
		k=2.66E-33*constants['M'][i]*np.exp(1140/constants['T'][i])
	elif krt=='k18':
		k=7.25E-11*np.exp(20/constants['T'][i])
	elif krt=='k19':
		k=4.63E-11*np.exp(20/constants['T'][i])
	elif krt=='k20':
		k=termol(1.8E-30,3.0,2.8E-11,0,constants,i)
	elif krt=='k21':
		k=bimol(2.45E-12,1775,constants,i)
	elif krt=='k22':
		k0=2.4E-14*np.exp(460/constants['T'][i])
		k2=2.7E-17*np.exp(2199/constants['T'][i])
		k3=6.5E-34*np.exp(1335/constants['T'][i])
		k=k0+k3*constants['M'][i]/(1+(k3*constants['M'][i])/k2)
	else:
		print krt+" ERROR: K RATE NOT DEFINED"
	return k
