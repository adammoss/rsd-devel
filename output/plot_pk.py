import matplotlib.pyplot as plt
import numpy as np

fileroot = 'test'

plt.figure()

data = np.loadtxt(fileroot+'_l0.dat')
k0 = data[:,0]
dd0 = data[:,1]
dt0 = data[:,2]
tt0 = data[:,3]
Acorr0 = data[:,4]
Bcorr0 = data[:,5]

#plt.plot(k0,dd0)
#plt.plot(k0,dd0+dt0)
plt.plot(k0,dd0+dt0+tt0,linestyle='--',color='k')
#plt.plot(k0,dd0+dt0+tt0+Acorr0,color='k')
plt.plot(k0,dd0+dt0+tt0+Acorr0+Bcorr0,color='k')

data = np.loadtxt(fileroot+'_l2.dat')
k2 = data[:,0]
dd2 = data[:,1]
dt2 = data[:,2]
tt2 = data[:,3]
Acorr2 = data[:,4]
Bcorr2 = data[:,5]

#plt.plot(k2,dd2)
#plt.plot(k2,dd2+dt2)
plt.plot(k2,dd2+dt2+tt2,linestyle='--',color='r')
#plt.plot(k2,dd2+dt2+tt2+Acorr2,color='r')
plt.plot(k2,dd2+dt2+tt2+Acorr2+Bcorr2,color='r')

data = np.loadtxt(fileroot+'_l4.dat')
k4 = data[:,0]
dd4 = data[:,1]
dt4 = data[:,2]
tt4 = data[:,3]
Acorr4 = data[:,4]
Bcorr4 = data[:,5]

#plt.plot(k4,dd4)
#plt.plot(k4,dd4+dt4)
plt.plot(k4,dd4+dt4+tt4,linestyle='--',color='b')
#plt.plot(k4,dd4+dt4+tt4+Acorr4,color='b')
plt.plot(k4,dd4+dt4+tt4+Acorr4+Bcorr4,color='b')


plt.xlim(0,0.2)
plt.ylim(500,100000)
plt.yscale('log')
plt.xlabel(u'$k [h \, Mpc^{-1}]$')
plt.ylabel(u'$P_L (k) [Mpc\, h^{-1}]^3$')
plt.savefig('pk_rsd.pdf', bbox_inches='tight')
