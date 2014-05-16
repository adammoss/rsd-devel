import numpy as np
import scipy.integrate as int
import scipy.interpolate as intp
import matplotlib.pyplot as plt
import sys

class Bias():

     def __init__(self,rootname,mu_min=-1.0,mu_max=0.999,n_mu=100,n_k=300,n_k_sample=10,kmax_sample=0.2,debug=False,hubble=0.7):

          CAMB = np.loadtxt(rootname)
          # CAMB outputs k/h
          k = CAMB[:,0]#*hubble
          P_k = CAMB[:,1]
          self.kmin = np.min(k)
          self.kmax = np.max(k)          
          self.mu_min = mu_min
          self.mu_max = mu_max
          self.n_mu = n_mu
          self.n_k = n_k
          self.dmu = (self.mu_max-self.mu_min)/self.n_mu
          self.dk = (self.kmax-self.kmin)/self.n_k
          self.Pk=intp.InterpolatedUnivariateSpline(k,P_k,k=3)
          
          print 'k_min: ',self.kmin
          print 'k_max: ',self.kmax

          self.kmin_sample = self.kmin
          self.kmax_sample = kmax_sample
          self.n_k_sample = n_k_sample
          self.dk_sample = (self.kmax_sample-self.kmin_sample)/(self.n_k_sample-1)
          karr = np.zeros(self.n_k_sample)
          Pa = np.zeros(self.n_k_sample)
          Pb2da = np.zeros(self.n_k_sample)
          Pb2ta = np.zeros(self.n_k_sample)
          Pbs2da = np.zeros(self.n_k_sample)
          Pbs2ta = np.zeros(self.n_k_sample)
          Pb2s2a = np.zeros(self.n_k_sample)
          Pbs22a = np.zeros(self.n_k_sample)
          sigma3a = np.zeros(self.n_k_sample)

          for i in range(0,self.n_k_sample):
               print 'k number: ',i
               kval = self.kmin+i*self.dk_sample
               karr[i] = kval
               Pa[i] = self.Pk(kval)
               Pb2da[i] =  self.int_trap(self.Pb2d,kval)
               Pb2ta[i] =  self.int_trap(self.Pb2t,kval)
               Pbs2da[i] =  self.int_trap(self.Pbs2d,kval)
               Pbs2ta[i] =  self.int_trap(self.Pbs2t,kval)
               Pb2s2a[i] =  self.int_trap(self.Pb2s2,kval)
               Pbs22a[i] =  self.int_trap(self.Pbs22,kval)
               sigma3a[i] = self.int_trap(self.sigma3,kval)
               if debug: print karr[i],Pa[i],Pb2da[i],Pb2ta[i],Pbs2da[i],Pbs2ta[i],Pb2s2a[i],Pbs22a[i],sigma3a[i]

          self.Pk_b2d=intp.InterpolatedUnivariateSpline(karr,Pb2da,k=3)
          self.Pk_b2t=intp.InterpolatedUnivariateSpline(karr,Pb2ta,k=3)
          self.Pk_bs2d=intp.InterpolatedUnivariateSpline(karr,Pbs2da,k=3)
          self.Pk_bs2t=intp.InterpolatedUnivariateSpline(karr,Pbs2ta,k=3)
          self.Pk_b2s2=intp.InterpolatedUnivariateSpline(karr,Pb2s2a,k=3)
          self.Pk_bs22=intp.InterpolatedUnivariateSpline(karr,Pbs22a,k=3)
          self.Pk_sigma3=intp.InterpolatedUnivariateSpline(karr,sigma3a,k=3)

     def Pk_dd(self,k,b1,b2,N,bs2,b3nl,determinate=False):
          if determinate:
               bs2 = -4.0/7.0*(b1-1.0)
               b3nl = 32.0/315.0*(b1-1.0)
          
          # Note ps22 -> In Florian's paper b22 isn't defined - should it be ps22?
          return b1**2*self.Pk(k)+2.0*b1*b2*self.Pk_b2d(k)+2.0*bs2*b1*self.Pk_bs2d(k)+2.0*b3nl*b1*self.Pk_sigma3(k)*self.Pk(k)+b2**2*self.Pk_bs22(k)+2*b2*bs2*self.Pk_b2s2(k)+bs2**2*self.Pk_bs22(k)+N
     

 #    def Pk_plot(self,b1,b2,N):
 #         karr = np.linspace(self.kmin_sample,self.kmax_sample,100)
 #         plt.plot(karr,self.Pk_dd(karr,b1,0,0)/self.Pk_dd(karr,b1,b2,N))
 #         plt.show()
          
     # q2 = k - q
     def q1dotq2(self,k1,q1,mu1):
          return q1*k1*mu1-q1**2

     def q2(self,k1,q1,mu1):
          return np.sqrt(q1**2+k1**2-2*k1*q1*mu1)

     # Eqn 48
     def FS2(self,k1,q1,mu1):
          return (5.0/7.0)+(self.q1dotq2(k1,q1,mu1)/(2.0*q1*self.q2(k1,q1,mu1)))*((q1/self.q2(k1,q1,mu1))+(self.q2(k1,q1,mu1)/q1))+(2.0/7.0)*(self.q1dotq2(k1,q1,mu1)/(q1*self.q2(k1,q1,mu1)))**2

     # Eqn 49
     def GS2(self,k1,q1,mu1):
          return 3.0/7.0+(self.q1dotq2(k1,q1,mu1)/(2*q1*self.q2(k1,q1,mu1)))*((q1/self.q2(k1,q1,mu1))+(q1/self.q2(k1,q1,mu1)))+(4.0/7.0)*(self.q1dotq2(k1,q1,mu1)/(q1*self.q2(k1,q1,mu1)))**2

     # Eqn 50
     def S2(self,k1,q1,mu1):
          return (self.q1dotq2(k1,q1,mu1)/(q1*self.q2(k1,q1,mu1)))**2-1.0/3.0

     # Eqn 42
     def Pb2d(self,mu1,q1,k1):
          return (1.0/(2*np.pi)**2)*q1**2*np.float(self.Pk(q1))*np.float(self.Pk(self.q2(k1,q1,mu1)))*self.FS2(k1,q1,mu1)

     # Eqn 43
     def Pb2t(self,mu1,q1,k1):
          return (1.0/(2*np.pi)**2)*q1**2*np.float(self.Pk(q1))*np.float(self.Pk(self.q2(k1,q1,mu1)))*self.GS2(k1,q1,mu1)
         
     # Eqn 44
     def Pbs2d(self,mu1,q1,k1):
          return (1.0/(2*np.pi)**2)*q1**2*np.float(self.Pk(q1))*np.float(self.Pk(self.q2(k1,q1,mu1)))*self.FS2(k1,q1,mu1)*self.S2(k1,q1,mu1)

     # Eqn 45
     def Pbs2t(self,mu1,q1,k1):
          return (1.0/(2*np.pi)**2)*q1**2*np.float(self.Pk(q1))*np.float(self.Pk(self.q2(k1,q1,mu1)))*self.GS2(k1,q1,mu1)*self.S2(k1,q1,mu1)

     # Eqn 46
     def Pb2s2(self,mu1,q1,k1):
          return -0.5*(1.0/(2*np.pi)**2)*q1**2*np.float(self.Pk(q1))*((2.0/3.0)*np.float(self.Pk(q1))-np.float(self.Pk(self.q2(k1,q1,mu1)))*self.S2(k1,q1,mu1))

     # Eqn 47
     def Pbs22(self,mu1,q1,k1):
          return -0.5*(1.0/(2*np.pi)**2)*q1**2*np.float(self.Pk(q1))*((4.0/9.0)*np.float(self.Pk(q1))-np.float(self.Pk(self.q2(k1,q1,mu1)))*(self.S2(k1,q1,mu1))**2)

     # Eqn 51
     def sigma3(self,mu1,q1,k1):
          return 105.0/16.0*(1.0/(2*np.pi)**2)*q1**2*np.float(self.Pk(q1))*((2.0/7.0*(mu1**2-1.0))*self.S2(k1,q1,mu1)+8.0/63.0)

     # 2D trap intergration
     def int_trap(self,f,ka):
          s1, s2, s3, s4, s5 = 0.0, 0.0, 0.0, 0.0, 0.0
          for i in range(1,self.n_mu):
               mu_i = self.mu_min + i*self.dmu 
               s1 += f(mu_i,self.kmin,ka)
               s2 += f(mu_i,self.kmax,ka)
          for i in range(1,self.n_k):
               k_i = self.kmin + i*self.dk 
               s3 += f(self.mu_min,k_i,ka)
               s4 += f(self.mu_max,k_i,ka)
          for i in range(1,self.n_mu):
               for j in range(1,self.n_k):
                    mu_i = self.mu_min + i*self.dmu 
                    k_i = self.kmin + j*self.dk 
                    s5 +=  f(mu_i,k_i,ka) 
          return 0.25*self.dmu*self.dk*(f(self.mu_min,self.kmin,ka)+f(self.mu_max,self.kmin,ka)+f(self.mu_min,self.kmax,ka)+f(self.mu_max,self.kmax,ka)+2.0*(s1+s2+s3+s4)+4.0*s5)
