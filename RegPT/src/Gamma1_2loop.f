c ******************************************************* c
c
      subroutine two_loop_Gamma1(a, k, ss)
c
c ******************************************************* c
c
c     label a = 1: density
c               2: velocity divergence
c
      implicit none
c
      integer a, ikmax, ik_max
      integer iq1, iq2, iqmax
      parameter(ikmax=2000)
      parameter(iqmax=500)
      real*8   ak(ikmax), pk(ikmax)
      real*8   kmin, kmax, k
      real*8   qq(iqmax), wq(iqmax), q1, q2
      real*8   integ_pkcorr_G1, integ_G1_2loop
      real*8   pi, ss, pklin_q1, pklin_q2
      common /pk_data/ ak, pk, ik_max
      pi = 4.d0 * datan(1.d0)
c     --------------------------------------------------------
c
      kmin = ak(1)
      kmax = ak(ik_max)
c
      ss = 0.d0
c     
c     ////// Gauss-Legendre integration over x (=q/k) //////  c
c
      call gauleg(dlog(kmin),dlog(kmax),qq,wq,iqmax)
c
      do iq1=1, iqmax
         q1 = dexp(qq(iq1))
         call find_pk(q1, pklin_q1)
c
         integ_pkcorr_G1 = 0.d0
c
         do iq2=1, iqmax
            q2 = dexp(qq(iq2))
            call find_pk(q2, pklin_q2)
            integ_pkcorr_G1 = integ_pkcorr_G1 + wq(iq2) * 
     *           integ_G1_2loop(a, k, q1, q2) * pklin_q2 * 
     *           q1**3 * q2**3
         enddo
c
         ss = ss + wq(iq1) * integ_pkcorr_G1 * pklin_q1
c
      enddo
c
      ss = ss / (2.d0*pi*pi)**2
c
      end
c
c ******************************************************* c
c
      function integ_G1_2loop(a, k, q1, q2)
c
c ******************************************************* c
c
c     Fitting form of the Integrand of G2 2-loop 
c                          (Angular integration has been done) 
c
c     Here, the second term of the integrand is considered
c
      implicit none
c
      integer a
      real*8  integ_G1_2loop
      real*8  k, q1, q2, x1, x2
      real*8  alpha
c     --------------------------------------------------------
c
      x1 = q1 / k
      x2 = q2 / k
c
      integ_G1_2loop = - alpha(a, x1, x2) / (x1**2 + x2**2)
c
      end
c
c ******************************************************* c
c
      function alpha(a, q1, q2)
c
c ******************************************************* c
c
      implicit none
c
      integer a
      real*8  alpha, q1, q2
      real*8  qr, qs
      real*8  beta, cc, dd, delta
c     --------------------------------------------------------
c
      qr = dsqrt( q1**2 + q2**2 )
      qs = q1*q2 / qr
c
      alpha = ( beta(a, q1, q2) - cc(a, qs) - dd(a, qr) ) * delta(a, qs)
c
c     ///// check /////
ccc      write(6,'(A,1p7e18.10)')  'check', q1, q2, beta(a, q1, q2), 
ccc     &     cc(a, qs),  dd(a, qr), delta(a, qs), alpha
c     ///// check /////
c
      end
c
c ******************************************************* c
c
      function beta(a, q1, q2) 
c
c ******************************************************* c
c
      implicit none
      integer a
      real*8  q1, q2, beta, y, beta_func5
      real*8  beta1, beta2
c     --------------------------------------------------------
c
      y = (q1 - q2) / (q1 + q2) 
c
      if(a.eq.1) then
         if(dabs(dabs(y)-1.d0).le.5.d-2) then
            beta = 120424.d0 / 3009825.d0 + 
     &           (2792.d0*(1.d0 - dabs(y))**2)/429975.d0 + 
     &           (2792.d0*(1.d0 - dabs(y))**3)/429975.d0 + 
     &           (392606.d0*(1.d0 - dabs(y))**4)/9.9324225d7
         elseif(dabs(y).le.1.d-2) then
            beta = 22382.d0/429975.d0 - 57052.d0*y**2/429975.d0
         else
            beta = (2.d0*(1.d0 + y**2)*(-11191.d0 + 118054.d0*y**2 - 
     &           18215.d0*y**4 + 18215.d0*y**8 - 118054.d0*y**10 + 
     &           11191.d0*y**12 + 60.d0*y**4*(3467.d0 - 790.d0*y**2 + 
     &           3467.d0*y**4)*dlog(y**2))) / 
     &           ( 429975.d0*(-1.d0 + y**2)**7 )
         endif
      elseif(a.eq.2) then
         if(dabs(dabs(y)-1.d0).le.5.d-2) then
            beta = 594232.d0 / 33108075.d0 + 
     &           (91912.d0*(1.d0 - dabs(y))**2)/33108075.d0 - 
     &           (91912.d0*(1.d0 - dabs(y))**3)/33108075.d0 + 
     &           (1818458.d0*(1.d0 - dabs(y))**4)/1092566475.d0
         elseif(dabs(y).le.1.d-2) then
            beta = 9886.d0/429975.d0 - 254356.d0*y**2/4729725.d0
         else
            beta = (2.d0*(1.d0 + y**2)*(-54373.d0 + 562162.d0*y**2 - 
     &           408245.d0*y**4 + 408245.d0*y**8 - 562162.d0*y**10 + 
     &           54373.d0*y**12 + 60.d0*y**4*(14561.d0 - 10690.d0*y**2 + 
     &           14561.d0*y**4)*dlog(y**2))) / 
     &           ( 4.729725d6*(-1.d0 + y**2)**7)
         endif
      endif
c
cc      write(6,'(A,1p3e18.10)') 'beta', q1, q2, beta
c
      end
c
c ******************************************************* c
c
      function cc(a, qs) 
c
c ******************************************************* c
c
      implicit none
      integer a
      real*8  qs, cc
c     --------------------------------------------------------
c
      if(a.eq.1) 
     &     cc = 0.02088557981734545 / (35.09866396385646*qs**4 + 
     &     4.133811416743832*qs**2 + 1.d0) - 
     &     0.076100391588544*qs**3 / (77.79670692480381*qs**6 + 1.d0)
c
      if(a.eq.2) 
     &     cc = - 0.008217140060512867 / (42.14072830553836*qs**4 + 
     &     1.367564560397748*qs**2 + 1.d0) + 
     &     0.01099093588476197*qs**3 / (28.490424851390667*qs**8 + 1.d0)
c
ccc      write(6,*) 'cc', qs, cc
c
      end
c
c ******************************************************* c
c
      function dd(a, qr)
c
c ******************************************************* c
c
      implicit none
c
      integer a
      real*8  qr, dd
c     --------------------------------------------------------
c
      if(a.eq.1) 
     &     dd = - 0.022168478217299517 / (7.030631093970638*qr**4 + 
     &     2.457866449142683*qr**2 + 1.d0) + 
     &     0.009267495321465601*qr**2 / (4.11633699497035*qr**10 + 1.d0)
c
      if(a.eq.2) 
     &     dd = 0.008023147297149955 / (2.238261369090066*qr**5 + 
     &     1.d0) - 0.006173880966928251*qr**2 / 
     &     (0.4711737436482179*qr**5 + 1.d0)
c
cc      write(6,*) 'dd', qr, dd
c
      end
c
c ******************************************************* c
c
      function delta(a, qs)
c
c ******************************************************* c
c
      implicit none
c
      integer a
      real*8  qs, delta
c     --------------------------------------------------------
c
      qs = dsqrt(2.d0) * qs
c
      if(a.eq.1) 
     &     delta = 0.3191221482038663*qs**4 / 
     &     (1.3549058352752525*qs**4 + 1.d0) 
     &     + 1.2805575495849764 / (18.192939946270577*qs**4 + 
     &     3.98817716852858*qs**2 + 1.d0) + 0.764469131436698
c
      if(a.eq.2) 
     &     delta = 1.528058751211026 * (2.4414566000839355*qs**4 + 
     &     1.8616263354608626*qs**2) / 
     &     (2.4414566000839355*qs**4 + 1.d0) + 
     &     2.5227965281961247 / (0.0028106312591877226*qs**4 + 
     &     1.0332351481570086*qs**2 + 1.d0) - 0.528058751211026
c
ccc      write(6,*) 'delta', qs, delta
c
      end
