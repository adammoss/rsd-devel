c ******************************************************* c
c
      subroutine calc_pkcorr_from_Gamma2(a, b, k, pkcorr_G2_tree_tree, 
     &     pkcorr_G2_tree_1loop, pkcorr_G2_1loop_1loop)
c
c ******************************************************* c
c
c     \int d^3q /(2*pi)^3 [Gamma^(2)(k-q,q,k)]^2 * P(k-q) * P(q)
c
c     without exp. cutoff
c
c     pkcorr_G2_tree_tree : Gamma_tree^(2) * Gamma_tree^(2)
c     pkcorr_G2_tree_1loop : Gamma_tree^(2) * Gamma_1loop^(2)
c     pkcorr_G2_1loop_1loop : Gamma_1loop^(2) * Gamma_1loop^(2)
c
      implicit none
c
      integer ik, ikmax, ik_max
      integer ix, ixmax
      integer a, b, aa, bb
      parameter(ikmax=2000)
ccc      parameter(ixmax=300)
      parameter(ixmax=200)
      real*8   ak(ikmax), pk(ikmax)
      real*8   kmin, kmax, k, kk
      real*8   xx(ixmax), wx(ixmax), xmin, xmax, x
      real*8   pkcorr_G2_tree_tree, pkcorr_G2_1loop_1loop
      real*8   pkcorr_G2_tree_1loop
      real*8   integ_pkcorr1, integ_pkcorr2, integ_pkcorr3
      real*8   pi, sigmav, exp_factor
      common /pk_data/ ak, pk, ik_max
      common /wave_number/  kk, xmin, xmax
      common /sigma_v/ sigmav
      common /label_a_b/ aa, bb 
      pi = 4.d0 * datan(1.d0)
c     --------------------------------------------------------
c
      kmin = ak(1)
      kmax = ak(ik_max)
      aa = a
      bb = b
c
      kk = k
      pkcorr_G2_tree_tree = 0.d0
      pkcorr_G2_tree_1loop = 0.d0
      pkcorr_G2_1loop_1loop = 0.d0
      exp_factor = 0.5d0 * (k*sigmav)**2
c
      if(k.lt.1.d-2 .or. exp_factor.gt.50.d0) then
         goto 5
      else
         xmin = kmin / k
         xmax = kmax / k
c
c     ////// Gauss-Legendre integration over x (=q/k) //////  c
c
         call gauleg(dlog(xmin),dlog(xmax),xx,wx,ixmax)
c
         do ix=1, ixmax
            x = dexp(xx(ix))
            call calc_integ_pkcorr_G2(x, integ_pkcorr1,
     &           integ_pkcorr2, integ_pkcorr3)
            pkcorr_G2_tree_tree = pkcorr_G2_tree_tree + 
     *           wx(ix) * integ_pkcorr1 * x**3
            pkcorr_G2_tree_1loop = pkcorr_G2_tree_1loop + 
     *           wx(ix) * integ_pkcorr2 * x**3
            pkcorr_G2_1loop_1loop = pkcorr_G2_1loop_1loop + 
     *           wx(ix) * integ_pkcorr3 * x**3
         enddo
c
         endif
c
 5       pkcorr_G2_tree_tree = 2.d0 * pkcorr_G2_tree_tree * 
     *        k**3 / (2.d0*pi)**2
         pkcorr_G2_tree_1loop = 2.d0 * pkcorr_G2_tree_1loop * 
     *        k**3 / (2.d0*pi)**2
         pkcorr_G2_1loop_1loop = 2.d0 * pkcorr_G2_1loop_1loop *
     *        k**3 / (2.d0*pi)**2
c     
      end
c
c ******************************************************* c
c
      subroutine calc_integ_pkcorr_G2(x, integ_pkcorr1,
     &     integ_pkcorr2, integ_pkcorr3)
c
c ******************************************************* c
c
c     Integration of pkcorr_G2 over angle, mu
c
      implicit none
c
      integer imu, imu_max, itype
ccc      parameter(imu_max=20)
      parameter(imu_max=10)
      real*8  integ_pkcorr1, integ_pkcorr2, integ_pkcorr3
      real*8  xmin, xmax, mumin, mumax
      real*8  kernel1, kernel2, kernel3
      real*8  k, x, mmu(imu_max), wmu(imu_max)
      common /wave_number/  k, xmin, xmax
c     --------------------------------------------------------
c
      integ_pkcorr1 = 0.d0
      integ_pkcorr2 = 0.d0
      integ_pkcorr3 = 0.d0
c
      mumin = max(-1.d0, (1.d0+x**2-xmax**2)/2.d0/x)
      mumax = min( 1.d0, (1.d0+x**2-xmin**2)/2.d0/x)
c
      if(x.ge.0.5d0) mumax= 0.5d0/x
c
      call gauleg(mumin, mumax, mmu, wmu, imu_max)
c
      do imu=1, imu_max
         call kernel_pkcorr_G2(x, mmu(imu), kernel1, kernel2, kernel3)
         integ_pkcorr1 = integ_pkcorr1 + wmu(imu) * kernel1
         integ_pkcorr2 = integ_pkcorr2 + wmu(imu) * kernel2
         integ_pkcorr3 = integ_pkcorr3 + wmu(imu) * kernel3
      enddo
c
      end
c
c ******************************************************* c
c
      subroutine kernel_pkcorr_G2(x, mu, kernel1, kernel2, kernel3)
c
c ******************************************************* c
c
c     Integration kernel for pkcorr_G2
c
      implicit none
c
      integer a, b
      real*8  kernel1, kernel2, kernel3, x, mu
      real*8  k, xmin, xmax, pk_q, pk_kq
      real*8  Gamma2_stdPT, q, kq
      real*8  G2a_tree, G2b_tree, G2a_1loop, G2b_1loop
      common /wave_number/  k, xmin, xmax
      common /label_a_b/ a, b 
c     --------------------------------------------------------
c
      kq = k*dsqrt(1.d0+x**2-2.d0*mu*x)
      q = k*x
c
      call find_pk(q, pk_q)
      call find_pk(kq, pk_kq)
c
      if(a.eq.b) then 
         G2a_tree = Gamma2_stdPT(a, 0, kq, q, k) 
         G2a_1loop = Gamma2_stdPT(a, 1, kq, q, k) 
         G2b_tree = G2a_tree
         G2b_1loop = G2a_1loop
      else
         G2a_tree = Gamma2_stdPT(a, 0, kq, q, k) 
         G2a_1loop = Gamma2_stdPT(a, 1, kq, q, k) 
         G2b_tree = Gamma2_stdPT(b, 0, kq, q, k) 
         G2b_1loop = Gamma2_stdPT(b, 1, kq, q, k) 
      endif
c
      kernel1 = 2.d0 * pk_q * pk_kq * G2a_tree * G2b_tree
      kernel2 = 2.d0 * pk_q * pk_kq * 
     &     ( G2a_tree * G2b_1loop + G2a_1loop * G2b_tree )
      kernel3 = 2.d0 * pk_q * pk_kq * G2a_1loop * G2b_1loop
c
      end
c
c ******************************************************* c
c
      function Gamma2_stdPT(a, n_loop, k1, k2, k3)
c
c ******************************************************* c
c
c     Gamma2 in standard PT 
c
c     Note--.  vec{k3} = vec{k1} + vec{k2} 
c
c     n_loop = 0 :  tree,    1 : one-loop
c
      implicit none
c
      integer a, n_loop
      real*8  Gamma2_stdPT, k1, k2, k3
      real*8  Gamma2_tree, Gamma2_1loop
c     --------------------------------------------------------
c
      if(n_loop.eq.0) then 
         Gamma2_stdPT = Gamma2_tree(a, k1, k2, k3)
      elseif(n_loop.eq.1) then
         if(a.eq.1) call one_loop_Gamma2d(k1, k2, k3, Gamma2_1loop)
         if(a.eq.2) call one_loop_Gamma2v(k1, k2, k3, Gamma2_1loop)
         Gamma2_stdPT = Gamma2_1loop
      else
         Gamma2_stdPT = 0.d0
      endif
c
      end
c
c ******************************************************* c
c
      function Gamma2_tree(a, k1, k2, k3)
c
c ******************************************************* c
c
      implicit none
c
      integer a
      real*8  Gamma2_tree, k1, k2, k3, k1dk2
c     --------------------------------------------------------
c
      k1dk2 = ( k3**2 - k1**2 - k2**2 ) / 2.d0
c
      if(a.eq.1) 
     &     Gamma2_tree = 5.d0/7.d0 + 0.5d0*k1dk2*
     &     (1.d0/k1**2 + 1.d0/k2**2) + 2.d0/7.d0*(k1dk2/(k1*k2))**2
      if(a.eq.2) 
     &     Gamma2_tree = 3.d0/7.d0 + 0.5d0*k1dk2*
     &     (1.d0/k1**2 + 1.d0/k2**2) + 4.d0/7.d0*(k1dk2/(k1*k2))**2
c
      end
c
c ******************************************************* c
c
      function LFunc(k,q)
c
c ******************************************************* c
c
      implicit none
      real*8 LFunc, k, q
c     -------------------------------------------------
c
      LFunc = dlog((k + q)**2/(k - q)**2)
c
      end
c
c ******************************************************* c
c
      function WFunc(k1, k2, k3, q)
c
c ******************************************************* c
c
      implicit none
      real*8 aa, bb
      real*8 WFunc, k1, k2, k3, q
c     -------------------------------------------------
c
      aa = -4.d0*k3**2*q**2 - 2.d0*(k1**2 - q**2)*(k2**2 - q**2) 
      bb = 4.d0*k3*q*dsqrt(k1**2*k2**2 + (-k1**2 - k2**2 + k3**2)*q**2 + 
     -     q**4)
c
      WFunc = dlog( (aa-bb)/(aa+bb) )
c
      end
c
c ******************************************************* c
c
      function betafunc(i,z)
c
c ******************************************************* c
c
c     Incomplete beta function of the type, B(z, i, 0)
c     (i=2, 4, 6)
c
c     this is indeed expressed in terms of elementary functions
c
      implicit none
      integer i
      real*8  z, a, betafunc
c     ----------------------------------
c
      betafunc = 0.d0
c
      if(z.lt.0.1d0) then
         a = dble(i)
         betafunc = z**i*(1.d0/a + z/(1.d0 + a) + z**2/(2.d0 + a) + 
     &        z**3/(3.d0 + a) + z**4/(4.d0 + a) + z**5/(5.d0 + a) + 
     &        z**6/(6.d0 + a))
      else
         if(i.eq.2) betafunc =
     &        (z**2*(-2.d0/z - (2.d0*dlog(1.d0 - z))/z**2))/2.d0
         if(i.eq.4) betafunc = 
     &        (z**4*((-2.d0*(6.d0 + 3.d0*z + 2.d0*z**2))/(3.d0*z**3) - 
     &        (4.d0*dlog(1.d0 - z))/z**4))/4.d0
         if(i.eq.6) betafunc =             
     &           (z**6*((-60.d0 - 30.d0*z - 20.d0*z**2 - 15.d0*z**3 - 
     &           12.d0*z**4)/(10.d0*z**5) - 
     &           (6.d0*dlog(1 - z))/z**6))/6.d0
      endif
c
      end
c
c ******************************************************* c
c
      function small_beta(k, q)
c
c ******************************************************* c
c
c     Regular function involving incomplete beta function: 
c
c     (k**2 - q**2) * Beta( 4*k*q/(k+q)**2, 4, 0) / (k + q)**2
c
      implicit none
      real*8 small_beta, k, q
      real*8 y, betafunc
c     ----------------------------------
c
      y = (k - q) / (k + q)
c
      if(dabs(y).lt.1.d-6) then
         small_beta = 0.d0
      else
         small_beta = y * betafunc(4, 1.d0 - y**2) 
      endif
c
cc      write(6,'(A,1p2e18.10)') 'small_beta=',small_beta, y
c
      end
c
c ******************************************************* c
c
      function big_beta(k1, k2, k3, q)
c
c ******************************************************* c
c
c     Regular function involving incomplete beta function:  
c
c     x * Beta(1-X**2, 2, 0) / sqrt( x + a**2 )
c
c     where a, x and X are defined as
c
c     a = k3 * q
c     x = (k1**2 - q**2) * (k2**2 - q**2)
c     X = {1 - sqrt(x + a**2)} / {1 + sqrt(x + a**2)} 
c

      implicit none
      real*8 big_beta, k1, k2, k3, q
      real*8 x, a, y, betafunc
c     ----------------------------------
c
      x = (k1**2-q**2) * (k2**2-q**2)
      a = k3 * q
      y = dsqrt(x+a*a)/a
c     
      if(dabs(y-1).lt.1.d-5) then
         big_beta = 0.d0
      elseif(dabs(y-1).lt.1.d-2) then
         big_beta = 2.d0*a*(-1.d0 + y)*(-1.d0 + dlog(4.d0) - 
     &        2.d0*dlog(dabs(-1.d0 + y))) + a*(-1.d0 + y)**2*
     &        (3.d0 - dlog(4.d0) + 2.d0*dlog(dabs(-1.d0 + y)))
cc         write(6,*) 'asymptotic expansion'
      else 
         big_beta = x * betafunc(2, 1-(1.d0-y)**2/(1.d0+y)**2) / (a*y)
cc         write(6,*) 'exact expression'
      endif
c
ccc      write(6,'(A,1p1e18.10)') 'big_beta=',big_beta
c
      end
