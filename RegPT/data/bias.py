import numpy as np
import scipy.integrate as int
import scipy.interpolate as intp
import matplotlib.pyplot as plt
import sys

class Bias():

     def __init__(self,rootname):

          CAMB = np.loadtxt(rootname)
          k = CAMB[:,0]
          P_k = CAMB[:,1]
          self.kmin = np.min(k)
          self.kmax = np.max(k)
          self.Pk=intp.InterpolatedUnivariateSpline(k,P_k,k=3)

          mu_min = -1.0
          mu_max = 0.999
          n_mu = 100
          n_k = 300
          dmu = (mu_max-mu_min)/n_mu
          dk = (kmax-kmin)/n_k

     def q1dotq2(self,k1,q1,mu1):
          return q1*k1*mu1-q1**2

     def q2(self,k1,q1,mu1):
          return np.sqrt(q1**2+k1**2-2*k1*q1*mu1)

          def FS2(self,k1,q1,mu1):
               return (5.0/7.0)+(q1dotq2(k1,q1,mu1)/(2.0*q1*q2(k1,q1,mu1)))*((q1/q2(k1,q1,mu1))+(q2(k1,q1,mu1)/q1))+(2.0/7.0)*(q1dotq2(k1,q1,mu1)/(q1*q2(k1,q1,mu1)))**2

          def GS2(self,k1,q1,mu1):
               return 3.0/7.0+(q1dotq2(k1,q1,mu1)/(2*q1*q2(k1,q1,mu1)))*((q1/q2(k1,q1,mu1))+(q1/q2(k1,q1,mu1)))+(4.0/7.0)*(q1dotq2(k1,q1,mu1)/(q1*q2(k1,q1,mu1)))**2

          def S2(self,k1,q1,mu1):
               return (q1dotq2(k1,q1,mu1)/(q1*q2(k1,q1,mu1)))**2-1.0/3.0

          def Pb2d(self,mu1,q1,k1):
               return (1/(2*np.pi)**2)*q1**2*np.float(Pk(q1))*np.float(Pk(q2(k1,q1,mu1)))*FS2(k1,q1,mu1)

          def Pb2dint(self,ka):
               return int.dblquad(Pb2d,kmin,kmax,lambda mu: -1, lambda mu: 1, args=(ka,))

          def Pb2t(self,mu1,q1,k1):
               return (1/(2*np.pi)**2)*q1**2*np.float(Pk(q1))*np.float(Pk(q2(k1,q1,mu1)))*GS2(k1,q1,mu1)

          def Pb2tint(self,ka):
               return int.dblquad(Pb2t,kmin,kmax,lambda mu: -1, lambda mu: 1, args=(ka,))

          def Pbs2d(self,mu1,q1,k1):
               return (1/(2*np.pi)**2)*q1**2*np.float(Pk(q1))*np.float(Pk(q2(k1,q1,mu1)))*FS2(k1,q1,mu1)*S2(k1,q1,mu1)

          def Pbs2dint(self,ka):
               return int.dblquad(Pbs2d ,kmin,kmax,lambda mu: -1, lambda mu: 1,args=(ka,))

          def Pbs2t(self,mu1,q1,k1):
               return (1/(2*np.pi)**2)*q1**2*np.float(Pk(q1))*np.float(Pk(q2(k1,q1,mu1)))*GS2(k1,q1,mu1)*S2(k1,q1,mu1)

          def Pbs2tint(self,ka):
               return int.dblquad(Pbs2t ,kmin,kmax,lambda mu: -1, lambda mu: 1,args=(ka,))

          def Pb2s2(self,mu1,q1,k1):
               return (-1/2)*(1/(2*np.pi)**2)*q1**2*np.float(Pk(q1))*((2/3)*np.float(Pk(q1))-np.float(Pk(q2(k1,q1,mu1)))*S2(k1,q1,mu1))

          def Pb2s2int(self,ka):
               return int.dblquad(Pb2s2 ,kmin,kmax,lambda mu: -1, lambda mu: 1,args=(ka,))

          def Pbs22(self,mu1,q1,k1):
               return (-1/2)*(1/(2*np.pi)**2)*q1**2*np.float(Pk(q1))*((4/9)*np.float(Pk(q1))-np.float(Pk(q2(k1,q1,mu1)))*(S2(k1,q1,mu1))**2)

          def Pbs22int(self,ka):
               return int.dblquad(Pbs22 ,kmin,kmax,lambda mu: -1, lambda mu: 1,args=(ka,))

def int_trap(f,ka):
     s1, s2, s3, s4, s5 = 0.0, 0.0, 0.0, 0.0, 0.0
     for i in range(1,n_mu):
          mu_i = mu_min + i*dmu 
          s1 += f(mu_i,kmin,ka)
          s2 += f(mu_i,kmax,ka)
     for i in range(1,n_k):
          k_i = kmin + i*dk 
          s3 += f(mu_min,k_i,ka)
          s4 += f(mu_max,k_i,ka)
     for i in range(1,n_mu):
          for j in range(1,n_k):
               mu_i = mu_min + i*dmu 
               k_i = kmin + j*dk 
               s5 +=  f(mu_i,k_i,ka)
     return 0.25*dmu*dk*(f(mu_min,kmin,ka)+f(mu_max,kmin,ka)+f(mu_min,kmax,ka)+f(mu_max,kmax,ka)+2.0*(s1+s2+s3+s4)+4.0*s5)

n = 10
karr = np.zeros(n)
Pa = np.zeros(n)
Pb2da = np.zeros(n)
Pb2ta = np.zeros(n)
Pbs2da = np.zeros(n)
Pbs2ta = np.zeros(n)
Pb2s2a = np.zeros(n)
Pbs22a = np.zeros(n)

for i in range(0,n):
     print i
     kval = (i+1)*0.01
     karr[i] = kval
     Pa[i] = Pk(kval)
     Pb2da[i] =  int_trap(Pb2d,kval)
     Pb2ta[i] =  int_trap(Pb2t,kval)
     Pbs2da[i] =  int_trap(Pbs2d,kval)
     Pbs2ta[i] =  int_trap(Pbs2t,kval)
     Pb2s2a[i] =  int_trap(Pb2s2,kval)
     Pbs22a[i] =  int_trap(Pbs22,kval)

plt.plot(karr,Pb2da/Pa)
plt.plot(karr,Pb2ta/Pa)
plt.plot(karr,Pbs2da/Pa)
plt.plot(karr,Pbs2ta/Pa)
plt.plot(karr,Pb2s2a/Pa)
plt.plot(karr,Pbs22a/Pa)
plt.show()
