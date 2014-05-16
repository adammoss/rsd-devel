c ******************************************************* c
c
      subroutine calc_pkcorr_from_Gamma3(a, b, k, pkcorr_G3_tree)
c
c ******************************************************* c
c
c     \int d^3p d^3q Gamma_a^(3)(p,q,k-p-q) * 
c                    Gamma_b^(3)(p,q,k-p-q) * P(p)*P(q)*P(|k-p-q|)
c     without exp. cutoff
c
c
c     labels a,  b = 1 : density
c                  = 2 : velocity divergence
c
      implicit none
c
      integer ndim6, ndim5, ncomp, mineval, maxeval
      integer verbose, last
      double precision epsrel, epsabs
      parameter (ndim5 = 5)
      parameter (ncomp = 1)
      parameter (epsrel = 0.005)
      parameter (epsabs = 0)
      parameter (verbose = 0)
      parameter (mineval = 0)
      parameter (maxeval=110000000)
      integer nstart, nincrease
      parameter (nstart = 4000)
      parameter (nincrease = 700)
      double precision integral(ncomp), error(ncomp), prob(ncomp)
      integer neval, fail
c     ---------------------------------------------------
      integer a, b, aa, bb
      integer ijob, ik, ik_max, ikmax, ik_init
      parameter(ikmax=2000)
      real*8  ak(ikmax), pk(ikmax)
      real*8  pkcorr_G3_tree
      real*8  error_pkcorr
      real*8  k, kk, kmin, kmax, xmin, xmax, pi
      real*8  sigmav, exp_factor
      external integ_pkcorr_G3
      common /pk_data/ ak, pk, ik_max
      common /wave_number/  kk, xmin, xmax
      common /G3_label/  aa, bb
      common /sigma_v/ sigmav
      pi = 4.d0 * datan(1.d0)
c     ---------------------------------------------------
c
      kmin = ak(1)
      kmax = ak(ik_max)
      aa = a
      bb = b
      kk = k
c
      xmin = kmin / k
      xmax = kmax / k
c
      exp_factor = 0.5d0 * (k*sigmav)**2
      pkcorr_G3_tree= 0.0d0
      integral(1) = 0.d0
      error(1) = 0.d0
      fail = 0
c
      if(k.lt.1.d-2 .or. exp_factor.gt.50.d0) then
         goto 5
      else
c
c     ////// 5D Monte Carlo integration //////  c
c
         call vegas(ndim5, ncomp, integ_pkcorr_G3,
     &        epsrel, epsabs, verbose + last, mineval, maxeval,
     &        nstart, nincrease, 
     &        neval, fail, integral, error, prob)
c
      endif
c     
 5    pkcorr_G3_tree = 6.d0 * integral(1) / (2.d0*pi)**6
      error_pkcorr = 6.d0 * error(1) / (2.d0*pi)**6
c
      end
c
c ******************************************************* c
c
      subroutine integ_pkcorr_G3(ndim, xx, ncomp, ff)
c
c ******************************************************* c
c
c     integrand of (Gamma3_tree)^2 for 5-dim integration 
c
      implicit none
c
      integer ndim, ncomp
      double precision xx(*), ff(*)
      double precision pi
      parameter (pi = 3.14159265358979323846D0)
c --------------------------------------- c
      integer a, b, i
      real*8  theta1, theta2, phi1, phi2, jacobian
      real*8  k, kmin, kmax, xmin, xmax
      real*8  F3_sym, kp(3), k_p, prod_q_kp
      real*8  pp(3), qq(3), kpq(3), p, q, k_p_q
      real*8  pklin_p, pklin_q, pklin_kpq
      common /wave_number/  k, xmin, xmax
      common /G3_label/  a, b
c
      kmin = k * xmin
      kmax = k * xmax
c
      p = dexp(dlog(kmin) + (dlog(kmax) - dlog(kmin)) * xx(1))
      theta1 = xx(2) * pi
      phi1 = xx(3) * 2.d0 * pi
      q = dexp(dlog(kmin) + (dlog(kmax) - dlog(kmin)) * xx(4))
      theta2 = xx(5) * pi
      jacobian = ( (dlog(kmax)-dlog(kmin)) * pi * 2.d0 * pi )**2
c
      pp(1) = p * dsin(theta1) * dcos(phi1)
      pp(2) = p * dsin(theta1) * dsin(phi1)
      pp(3) = p * dcos(theta1) 
      qq(1) = q * dsin(theta2) 
      qq(2) = 0.d0
      qq(3) = q * dcos(theta2) 
c
      kpq(1) =  -pp(1)-qq(1)
      kpq(2) =  -pp(2)-qq(2)
      kpq(3) = k-pp(3)-qq(3)
c
      k_p_q = dsqrt( kpq(1)**2 + kpq(2)**2 + kpq(3)**2 )
c
      kp(1) = -pp(1)
      kp(2) = -pp(2)
      kp(3) = k-pp(3)
      k_p = dsqrt( kp(1)**2 + kp(2)**2 + kp(3)**2 )
      prod_q_kp = qq(1)*kp(1) + qq(2)*kp(2) + qq(3)*kp(3)
c
      if(k_p_q.ge.kmin .and. k_p_q.le.kmax) then 
c
         if(prod_q_kp .le. 0.5d0*k_p*k_p) then
c
            call find_pk(p, pklin_p)
            call find_pk(q, pklin_q)
            call find_pk(k_p_q, pklin_kpq)
c
            ff(1) = F3_sym(a, pp, qq, kpq) * F3_sym(b, pp, qq, kpq)
     &           * pklin_p * pklin_q * pklin_kpq
            ff(1) = ff(1) * p**3 * q**3 * dsin(theta1) * dsin(theta2)
c
         else
            ff(1) = 0.d0
         endif
      else
         ff(1) = 0.d0
      endif
c
      ff(1) = 2.d0 * ff(1) * jacobian
c
      end
c
