! Code to compute redshift space anisotropic P(k) 

module pk_rsd 
  use am_routines
  implicit none
  integer, parameter :: DP = kind(1.0D0)
  real(DP), parameter :: pi = 3.1415926535897932384626433832795d0, twopi=2*pi, fourpi=4*pi
     
  integer :: nk = 100
  real(DP), allocatable :: ak(:), pk(:)
  
  integer :: ik_max
  real(DP), allocatable :: ak_camb(:), pk_camb(:)

  real(DP) :: sigmav ! Sigma_v
  real(DP) :: growth ! Growth factor 
  real(DP) :: ff ! Growth rate 

  real(DP), allocatable :: pk0corr(:), pk2corr(:), pk4corr(:)

contains
   
  subroutine init_pk()
    implicit none

    allocate(pk0corr(nk),pk2corr(nk),pk4corr(nk))
    pk0corr(:) = 0.0d0
    pk2corr(:) = 0.0d0
    pk4corr(:) = 0.0d0

  end subroutine init_pk
 
! ******************************************************* 

      subroutine load_matterpower_data(filename)

!     input file is assumed to be the matter power spectrum data 
!     created by CAMB code. 

      implicit none
      character(len=100) filename

      integer ikmax
      parameter(ikmax=3000)
      integer ik, ikk
      real(DP) :: ak_temp(ikmax), pk_temp(ikmax), dlnk
      
      open(9, file=trim(filename), status='unknown')
      
      do ik=1, ikmax
         read(9,*,END=10) ak_temp(ik), pk_temp(ik)
      enddo
 10   continue
      close(9)

      ik_max = ik-1
      allocate(ak_camb(ik_max),pk_camb(ik_max))
      ak_camb(1:ik_max) = ak_temp(1:ik_max)
      pk_camb(1:ik_max) = pk_temp(1:ik_max)

      !******************************************
      ! REMOVE THIS

      pk_camb(:)=pk_camb(:)*3.83

      !**************************************

      dlnk=log(ak_temp(ik-1)/ak_temp(1))/dble(nk-1)

      allocate(ak(nk),pk(nk))
      do ik=1,nk
         ak(ik) = ak_temp(1)*exp(dble(ik-1)*dlnk)
         pk(ik) = find_pk(ak(ik))
         write(*,*) ak(ik),pk(ik)
      end do

      write(6,*) 'ak(1)=', ak(1)
      write(6,*) 'ak(ik_max)=', ak(nk)
      
      call init_pk()

    end subroutine load_matterpower_data

! ******************************************************* 

    function find_pk(kk)

      implicit none
      integer j, jmin, jmax
      real(DP) :: kk, s, ds, find_pk

      call hunt(ak_camb(1:ik_max), kk, j)
      jmin = j - 2
      jmax = j + 2
      if(jmin.lt.1) jmin = 1
      if(jmax.ge.ik_max) jmax = ik_max
      call polint(ak_camb(jmin:jmax),pk_camb(jmin:jmax),kk,s,ds)
      find_pk = s
      
    end function find_pk

! ******************************************************* 

    function fp(ip, x, mu,k,xmin,xmax)

!     ip=1 for kernel of pk_B111
!     ip=2 for kernel of pk_B112
!     ip=3 for kernel of pk_B121
!     ip=4 for kernel of pk_B122
!     ip=5 for kernel of pk_B211
!     ip=6 for kernel of pk_B212
!     ip=7 for kernel of pk_B221
!     ip=8 for kernel of pk_B222
!     ip=9 for kernel of pk_B312
!     ip=10 for kernel of pk_B321
!     ip=11 for kernel of pk_B322
!     ip=12 for kernel of pk_B422
!
      implicit none
      integer ip
      real(DP) fp
      real(DP) mu, k, x, xmin, xmax
      real(DP) mumin, mumax

      if(ip.eq.1) then
         fp = x**2 * (mu*mu-1.) / 2.
      elseif(ip.eq.2) then
         fp = 3.*x**2 * (mu*mu-1.)**2 / 8.
      elseif(ip.eq.3) then
         fp = 3.*x**4 * (mu*mu-1.)**2 / (1.+x*x-2.*mu*x) / 8.
      elseif(ip.eq.4) then
         fp = 5.*x**4 * (mu*mu-1.)**3 / (1.+x*x-2.*mu*x) / 16.
      elseif(ip.eq.5) then
         fp = x * (x+2.*mu-3.*x*mu*mu) / 2.
      elseif(ip.eq.6) then
         fp = - 3.*x * (mu*mu-1.) * (-x-2.*mu+5.*x*mu*mu) / 4.
      elseif(ip.eq.7) then
         fp = 3.*x**2 * (mu*mu-1.) * (-2.+x*x+6.*x*mu-5.*x*x*mu*mu) &
             / (1.+x*x-2.*mu*x) / 4.
      elseif(ip.eq.8) then
         fp = - 3.*x**2 * (mu*mu-1.)**2 &
             * (6.-5.*x*x-30.*x*mu+35.*x*x*mu*mu) &
             / (1.+x*x-2.*mu*x) / 16.
      elseif(ip.eq.9) then
         fp = x * (4.*mu*(3.-5.*mu*mu) + x*(3.-30.*mu*mu+35.*mu**4) )/ 8.
      elseif(ip.eq.10) then
         fp = x * (-8.*mu + x*(-12.+36.*mu*mu+12.*x*mu*(3.-5.*mu*mu)+ &
             x**2*(3.-30.*mu*mu+35.*mu**4) ) ) / (1.+x*x-2.*mu*x) / 8.
      elseif(ip.eq.11) then
         fp = 3.*x * (mu*mu-1.) * (-8.*mu + x*(-12.+60.*mu*mu+ &
             20.*x*mu*(3.-7.*mu*mu)+5.*x*x*(1.-14.*mu*mu+21.*mu**4)) ) &
             / (1.+x*x-2.*mu*x) / 16.
      elseif(ip.eq.12) then
         fp = x * (8.*mu*(-3.+5.*mu*mu) - 6.*x*(3.-30.*mu*mu+35.*mu**4) &
             + 6.*x*x*mu*(15.-70.*mu*mu+63*mu**4) + x**3*(5.-21.*mu*mu* &
             (5.-15.*mu*mu+11.*mu**4)) ) / (1.+x*x-2.*mu*x) / 16.
      endif

      fp = fp * x  / (1.+x*x-2.*mu*x)* find_pk(k*x)*find_pk(k*sqrt(1.+x*x-2.*mu*x))
      
    end function fp


! ******************************************************* 

    function integ_fp(ip, x,k,xmin,xmax)

      implicit none
      integer ip, imu, imu_max
      parameter(imu_max=10)
      real(DP)  integ_fp, xmin, xmax, mumin, mumax
      real(DP)  k, x, mu, wmu(imu_max), mmu(imu_max)
     
      integ_fp = 0.d0

      mumin = max(-1.0, (1.+x**2-xmax**2)/2./x)
      mumax = min( 1.0, (1.+x**2-xmin**2)/2./x)

      if(x.ge.0.5d0) mumax= 0.5d0/x

      call gaulegf(mumin, mumax, mmu, wmu, imu_max)

      do imu=1, imu_max
         integ_fp = integ_fp + wmu(imu) * fp(ip, x, mmu(imu),k,xmin,xmax)
      enddo

    end function integ_fp

! ******************************************************* 

    subroutine calc_correction
      
      implicit none
      integer ik, isub
      integer ix, ixmax
      parameter(ixmax=600)
      real(DP) :: kmin, kmax, xmin, xmax, mumin, mumax
      real(DP) :: k, ww(ixmax), xx(ixmax)
      real(DP) :: alpha
      real(DP) :: pk_B111, pk_B112, pk_B121
      real(DP) :: pk_B122, pk_B211, pk_B212
      real(DP) :: pk_B221, pk_B222, pk_B312
      real(DP) :: pk_B321, pk_B322, pk_B422
      real(DP) :: pk_B1, pk_B2, pk_B3, pk_B4

      kmin = ak(1)
      kmax = ak(nk) 

      do ik=1, nk

         pk_B111 = 0.0d0
         pk_B112 = 0.0d0
         pk_B121 = 0.0d0
         pk_B122 = 0.0d0
         pk_B211 = 0.0d0
         pk_B212 = 0.0d0
         pk_B221 = 0.0d0
         pk_B222 = 0.0d0
         pk_B312 = 0.0d0
         pk_B321 = 0.0d0
         pk_B322 = 0.0d0
         pk_B422 = 0.0d0
         pk_B1 = 0.0d0
         pk_B2 = 0.0d0
         pk_B3 = 0.0d0
         pk_B4 = 0.0d0

         k = ak(ik)

         xmin = kmin / k
         xmax = kmax / k

!     ////// Gauss-Legendre integration //////  c

         if(k.lt.0.2) isub =200 
         if(k.ge.0.2) isub =0 
         !       
         call gaulegf(log(xmin),log(xmax),xx,ww,ixmax-isub)

         do ix=1, ixmax-isub
            xx(ix)= dexp(xx(ix))
            pk_B111 = pk_B111+ww(ix)*integ_fp(1,xx(ix),k,xmin,xmax)
            pk_B112 = pk_B112+ww(ix)*integ_fp(2,xx(ix),k,xmin,xmax)
            pk_B121 = pk_B121+ww(ix)*integ_fp(3,xx(ix),k,xmin,xmax)
            pk_B122 = pk_B122+ww(ix)*integ_fp(4,xx(ix),k,xmin,xmax)
            pk_B211 = pk_B211+ww(ix)*integ_fp(5,xx(ix),k,xmin,xmax)
            pk_B212 = pk_B212+ww(ix)*integ_fp(6,xx(ix),k,xmin,xmax)
            pk_B221 = pk_B221+ww(ix)*integ_fp(7,xx(ix),k,xmin,xmax)
            pk_B222 = pk_B222+ww(ix)*integ_fp(8,xx(ix),k,xmin,xmax)
            pk_B312 = pk_B312+ww(ix)*integ_fp(9,xx(ix),k,xmin,xmax)
            pk_B321 = pk_B321+ww(ix)*integ_fp(10,xx(ix),k,xmin,xmax)
            pk_B322 = pk_B322+ww(ix)*integ_fp(11,xx(ix),k,xmin,xmax)
            pk_B422 = pk_B422+ww(ix)*integ_fp(12,xx(ix),k,xmin,xmax)
         enddo

         pk_B111 = 2.d0 * pk_B111 * k**3 / (2.*pi)**2
         pk_B112 = - 2.d0 * pk_B112 * k**3 / (2.*pi)**2
         pk_B121 = - 2.d0 * pk_B121 * k**3 / (2.*pi)**2
         pk_B122 = 2.d0 * pk_B122 * k**3 / (2.*pi)**2
         pk_B211 = 2.d0 * pk_B211 * k**3 / (2.*pi)**2
         pk_B212 = - 2.d0 * pk_B212 * k**3 / (2.*pi)**2
         pk_B221 = - 2.d0 * pk_B221 * k**3 / (2.*pi)**2
         pk_B222 = 2.d0 * pk_B222 * k**3 / (2.*pi)**2
         pk_B312 = - 2.d0 * pk_B312 * k**3 / (2.*pi)**2
         pk_B321 = - 2.d0 * pk_B321 * k**3 / (2.*pi)**2
         pk_B322 = 2.d0 * pk_B322 * k**3 / (2.*pi)**2
         pk_B422 = 2.d0 * pk_B422 * k**3 / (2.*pi)**2

         !write(6,'(i4,1p4e18.10)') ik,k,pk(ik),pk_B111(ik),pk_B112(ik)
         
         pk_B1 = ff**2*pk_B111 + ff**3*pk_B112 + ff**3*pk_B121 + ff**4*pk_B122
         pk_B2 = ff**2*pk_B211 + ff**3*pk_B212 + ff**3*pk_B221 + ff**4*pk_B222
         pk_B3 = ff**3*pk_B312 + ff**3*pk_B321 + ff**4*pk_B322
         pk_B4 = ff**4*pk_B422

         alpha = (k*ff*sigmav)**2.0

         pk0corr(ik) = fact(1,0,alpha) * pk_B1 &
              + fact(2,0,alpha) * pk_B2 + fact(3,0,alpha) * pk_B3 &
              + fact(4,0,alpha) * pk_B4
         
         pk2corr(ik) = fact(1,2,alpha) * pk_B1 &
              + fact(2,2,alpha) * pk_B2 + fact(3,2,alpha) * pk_B3 &
              + fact(4,2,alpha) * pk_B4

         pk4corr(ik) = fact(1,4,alpha) * pk_B1 &
              + fact(2,4,alpha) * pk_B2 + fact(3,4,alpha) * pk_B3 &
              + fact(4,4,alpha) * pk_B4
         
         write(6,'(i4,1p4e18.10)') ik,k,pk0corr(ik),pk2corr(ik),pk4corr(ik)
         
      end do

    end subroutine calc_correction

! ******************************************************* 

      subroutine calc_pkred
        
!     Summing up all contributions to redshift P(k) in PT
!     and calculating monopole, quadrupole and hexadecapole
!     spectra


      end subroutine calc_pkred

! ************************************************ 

      function fact(n, l, alpha)

!     (2l+1)/2 * integ dmu  mu^(2n) * exp(-alpha*mu^2) * P_l(mu)

      implicit none
      integer n, l
      real(DP) fact, nn, alpha
      nn = dble(n)

      if(alpha.gt.0.05) then

         if(l.eq.0) then
            fact = gamhalf(n) * gammp(0.5+nn,alpha)
            fact = fact / alpha**(nn+0.5) / 4.d0
         elseif(l.eq.2) then
            fact = alpha * gamhalf(n) * gammp(0.5+nn,alpha) &
                - 3.d0 * gamhalf(n+1) * gammp(1.5+nn,alpha) 
            fact = fact / alpha**(nn+1.5) * (-5.d0/8.d0)
         elseif(l.eq.4) then
            fact = 12.*gamhalf(n)*gammp(0.5+nn,alpha)/alpha**(n+0.5)  &
                -120.*gamhalf(n+1)*gammp(1.5+nn,alpha)/alpha**(n+1.5) &
                +140.*gamhalf(n+2)*gammp(2.5+nn,alpha)/alpha**(n+2.5)
            fact = fact * 9./128.
         endif

      else

         if(l.eq.0) then
            fact = 1./(2.+4.*nn) - alpha/(6.+4.*nn) + alpha**2/(20.+8.*nn) 
         elseif(l.eq.2) then
            fact = nn/(3.+8.*nn+4.*nn**2) &
                - (nn+1.)*alpha/(15.+16.*nn+4.*nn**2) &
                + (nn+2.)*alpha**2/(70.+48.*nn+8.*nn**2)
            fact = fact * 5.d0
         elseif(l.eq.4) then
            fact = dble(n*(n-1))/dble(15+46*n+36*n**2+8*n**3) &
                - dble(n*(n+1))/dble(105+142*n+60*n**2+8*n**3)*alpha &
                + dble((n+1)*(n+2))/dble(315+286*n+84*n**2+8*n**3) &
                *alpha**2/2.d0
            fact = fact * 18.d0
         endif
         
      endif

      fact = fact * (1.d0 + (-1.d0)**(2.*nn)) 
      
    end function fact

  end module pk_rsd


