module pk_rsd 
  use am_routines
  implicit none
  integer, parameter :: DP = kind(1.0D0)
  real(DP), parameter :: pi = 3.1415926535897932384626433832795d0, twopi=2*pi, fourpi=4*pi
     
  integer :: nk = 100
  real(DP), allocatable :: ak(:), pk(:)

  integer imu_max
  parameter(imu_max=10)
  integer ixmax
  parameter(ixmax=2000)
  
  integer :: ik_max
  real(DP), allocatable :: ak_camb(:), pk_camb(:)

  real(DP) :: sigmav ! Sigma_v
  real(DP) :: growth ! Growth factor 
  real(DP) :: ff ! Growth rate 
  real(DP) :: sigma_8 = 0.0d0

  real(DP) :: b1

  real(DP), allocatable :: pk0dd(:),pk2dd(:),pk4dd(:)
  real(DP), allocatable :: pk0dt(:),pk2dt(:),pk4dt(:)
  real(DP), allocatable :: pk0tt(:),pk2tt(:),pk4tt(:)
  real(DP), allocatable :: pk0corr_A(:), pk2corr_A(:), pk4corr_A(:)
  real(DP), allocatable :: pk0corr_B(:), pk2corr_B(:), pk4corr_B(:)

contains
   
  subroutine init_pk()
    implicit none

    allocate(pk0dd(nk),pk2dd(nk),pk4dd(nk))
    allocate(pk0dt(nk),pk2dt(nk),pk4dt(nk))
    allocate(pk0tt(nk),pk2tt(nk),pk4tt(nk))
    allocate(pk0corr_A(nk),pk2corr_A(nk),pk4corr_A(nk))
    allocate(pk0corr_B(nk),pk2corr_B(nk),pk4corr_B(nk))
    pk0dd(:) = 0.0d0
    pk2dd(:) = 0.0d0
    pk4dd(:) = 0.0d0
    pk0dt(:) = 0.0d0
    pk2dt(:) = 0.0d0
    pk4dt(:) = 0.0d0
    pk0tt(:) = 0.0d0
    pk2tt(:) = 0.0d0
    pk4tt(:) = 0.0d0
    pk0corr_A(:) = 0.0d0
    pk2corr_A(:) = 0.0d0
    pk4corr_A(:) = 0.0d0
    pk0corr_B(:) = 0.0d0
    pk2corr_B(:) = 0.0d0
    pk4corr_B(:) = 0.0d0

  end subroutine init_pk
 
! ******************************************************* 

      subroutine load_matterpower_data(filename)

!     input file is assumed to be the matter power spectrum data 
!     created by CAMB code. 

      implicit none
      character(len=200) filename

      integer ikmax
      parameter(ikmax=10000)
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

      if (sigma_8 .ne. 0.0d0) call normalization_trapez(8.0d0)

      dlnk=log(ak_temp(ik-1)/ak_temp(1))/dble(nk-1)

      allocate(ak(nk),pk(nk))
      do ik=1,nk
         ak(ik) = ak_temp(1)*exp(dble(ik-1)*dlnk)
         pk(ik) = find_pk(ak(ik))
      end do

      write(6,*) 'k_min =', ak(1)
      write(6,*) 'k_max =', ak(nk)
      
      call init_pk()

    end subroutine load_matterpower_data

! ******************************************************* 

    subroutine normalization_trapez(r_th)
      implicit none
      real(DP) r_th
      integer ik
      real(DP) :: W_TH,sigma_a,sigma_b,x,const

      x = ak_camb(1) * r_th
      if(x.lt.1.d-3) then
         W_TH = 1.d0 - x*x / 10.d0 + x**4 / 280.d0 
      else
         W_TH = 3.d0 * (sin(x) - x * cos(x))/x/x/x
      endif
      sigma_a = W_TH * W_TH * pk_camb(1) * ak_camb(1) * ak_camb(1)
      sigma_a = sigma_a / (2.d0 * pi * pi)

      const = 0.d0 
      do ik=2, ik_max
         x = ak_camb(ik) * r_th
         if(x.lt.1.d-3) then
            W_TH = 1.d0 - x*x / 10.d0 + x**4 / 280.d0 
         else
            W_TH = 3.d0 * (sin(x) - x * cos(x))/x/x/x
         endif
         sigma_b = W_TH * W_TH * pk_camb(ik) * ak_camb(ik) * ak_camb(ik) 
         sigma_b = sigma_b / (2.d0 * pi * pi)
         const = const + (sigma_a + sigma_b) * ( ak_camb(ik) - ak_camb(ik-1) )/ 2.d0
         sigma_a = sigma_b
      enddo

      write(*,*) 'Sigma_8 = ',sqrt(const)
      pk_camb(:) = pk_camb(:) * (sigma_8**2.0)/const
      write(*,*) 'Rescaling P(k) to sigma_8 = ',sigma_8

    end subroutine normalization_trapez

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

      function  fp_A(ip, x, mu,k,xmin,xmax)

!     ip=1 for kernel of pk_A11
!     ip=2 for kernel of pk_A12
!     ip=3 for kernel of pk_A22
!     ip=4 for kernel of pk_A23
!     ip=5 for kernel of pk_A33

      implicit none
      integer ip
      real*8 fp_A
      real*8 mu, k, x, xmax, xmin

      if(ip.eq.1) then
         fp_A = - x**3 * ( mu+6.*mu**3+x*x*mu*(-3.+10.*mu*mu)+ &
             x*(-3.+mu*mu-12.*mu**4) ) / 7.
      elseif(ip.eq.2 .or. ip.eq.4) then
         fp_A = x**4 * (mu*mu-1.) * (-1.+7.*x*mu-6.*mu*mu) / 14.
      elseif(ip.eq.3) then
         fp_A = x**3 * ( x*x*mu*(13.-41.*mu*mu) -4.*(mu+6.*mu**3) &
             + x*(5.+9.*mu*mu+42.*mu**4) ) / 14.
      elseif(ip.eq.5) then
         fp_A = x**3 * (1.-7.*x*mu+6.*mu*mu) &
             * (-2.*mu+x*(-1.+3.*mu*mu)) / 14.
      endif

      fp_A = fp_A * x * find_pk(k) * find_pk(k*sqrt(1.+x*x-2.*mu*x)) / (1.+x*x-2.*mu*x)**2

      end function fp_A

! ******************************************************* 

      function  fp_tA(ip, x, mu,k,xmin,xmax)

!     ip=6 for kernel of pk_tA11
!     ip=7 for kernel of pk_tA12
!     ip=8 for kernel of pk_tA22
!     ip=9 for kernel of pk_tA23
!     ip=10 for kernel of pk_tA33

      implicit none
      integer ip
      real*8 fp_tA
      real*8 mu, k, x, xmax, xmin
      
      if(ip.eq.6) then
         fp_tA = (-mu+x*(2.*mu*mu-1.)) * (-3.*x-7.*mu+10.*x*mu*mu) / 7.
      elseif(ip.eq.7) then
         fp_tA = x * (mu*mu-1.) * (3.*x+7.*mu-10.*x*mu*mu) / 14.
      elseif(ip.eq.8) then
         fp_tA = ( 28.*mu*mu + x*mu*(25.-81.*mu*mu) &
             + x**2*(1.-27.*mu*mu+54.*mu**4) ) / 14.
      elseif(ip.eq.9) then
         fp_tA = - x * (mu*mu-1.) * (x-7.*mu+6.*x*mu*mu) / 14. 
      elseif(ip.eq.10) then
         fp_tA = ( x-7.*mu+6.*x*mu*mu ) * ( -2.*mu + &
             x*(3.*mu*mu-1.) ) / 14.
      endif

      fp_tA = fp_tA * x * find_pk(k*x) * find_pk(k*sqrt(1.+x*x-2.*mu*x)) / (1.+x*x-2.*mu*x)**2

      end function fp_tA

! ******************************************************* 

    function  fp_aa(ip, x, mu,k,xmin,xmax)

!     ip=11 for kernel of pk_aa11
!     ip=12 for kernel of pk_aa12
!     ip=13 for kernel of pk_aa22
!     ip=14 for kernel of pk_aa23
!     ip=15 for kernel of pk_aa33

      implicit none
      integer ip
      real*8 fp_aa
      real*8 mu, k, x, xmax, xmin

      if(ip.eq.11) then
         fp_aa = ( -7.*mu*mu + x**3*mu * (-3.+10.*mu*mu) + 3.*x &
             * (mu+6.*mu**3) + x*x * (6.-19.*mu*mu-8.*mu**4) ) / 7.
      elseif(ip.eq.12 .or. ip.eq.14) then
         fp_aa = x * (-1.+mu*mu) * (6.*x - 7.*(1.+x*x)*mu + 8.*x*mu*mu) &
             / 14.
      elseif(ip.eq.13) then
         fp_aa = ( -28.*mu*mu + x**3*mu* (-13.+41.*mu*mu) + &
             x*mu*(11.+73.*mu*mu) - 2.*x*x*(-9.+31.*mu*mu+20.*mu**4) ) &
             / 14.
      elseif(ip.eq.15) then
         fp_aa = ( 7.*mu + x * (-6.+7.*x*mu-8.*mu*mu) ) * &
             ( -2.*mu + x*(-1.+3.*mu*mu) ) / 14.
      endif

      fp_aa = fp_aa * x * find_pk(k) * find_pk(k*x) / (1.+x*x-2.*mu*x)

      end function fp_aa

! ******************************************************* 

    function integ_fp_B(ip, x,k,xmin,xmax)

      implicit none
      integer ip, imu
      real(DP)  integ_fp_B, xmin, xmax, mumin, mumax
      real(DP)  k, x, mu, wmu(imu_max), mmu(imu_max)
     
      integ_fp_B = 0.d0

      mumin = max(-1.0, (1.+x**2-xmax**2)/2./x)
      mumax = min( 1.0, (1.+x**2-xmin**2)/2./x)

      if(x.ge.0.5d0) mumax= 0.5d0/x

      call gaulegf(mumin, mumax, mmu, wmu, imu_max)

      do imu=1, imu_max
         integ_fp_B = integ_fp_B + wmu(imu) * fp(ip, x, mmu(imu),k,xmin,xmax)
      enddo

    end function integ_fp_B


    function integ_fp_A(ip, x,k,xmin,xmax)

      implicit none
      integer ip, imu
      real*8  integ_fp_A, xmin, xmax, mumin, mumax
      real*8  k, x, mu, wmu(imu_max), mmu(imu_max)

      integ_fp_A = 0.d0

      mumin = max(-1.0, (1.+x**2-xmax**2)/2./x)
      mumax = min( 1.0, (1.+x**2-xmin**2)/2./x)

      if(x.ge.0.5d0) mumax= 0.5d0/x

      call gaulegf(mumin, mumax, mmu, wmu, imu_max)

         do imu=1, imu_max
            if(ip.le.5)  integ_fp_A = integ_fp_A + wmu(imu) * fp_A(ip, x, mmu(imu),k,xmin,xmax)
            if(ip.ge.6 .and. ip.le.10) integ_fp_A = integ_fp_A + wmu(imu) * fp_tA(ip, x, mmu(imu),k,xmin,xmax)
            if(ip.ge.11) integ_fp_A = integ_fp_A + wmu(imu) * fp_aa(ip, x, mmu(imu),k,xmin,xmax)
         enddo
        
     end function integ_fp_A

! ******************************************************* 

    subroutine calc_correction(ik)
      
      implicit none
      integer ik, isub
      integer ix
      real(DP) :: kmin, kmax, xmin, xmax, mumin, mumax
      real(DP) :: k, ww(ixmax), xx(ixmax),kfact, pk_val
      real(DP) :: alpha
      real(DP) :: pk_B111, pk_B112, pk_B121
      real(DP) :: pk_B122, pk_B211, pk_B212
      real(DP) :: pk_B221, pk_B222, pk_B312
      real(DP) :: pk_B321, pk_B322, pk_B422
      real(DP) :: pk_B1, pk_B2, pk_B3, pk_B4
      real(DP) ::  pk_A11, pk_A12, pk_A22
      real(DP) ::  pk_A23, pk_A33
      real(DP) ::  pk_tA11, pk_tA12, pk_tA22
      real(DP) ::  pk_tA23, pk_tA33
      real(DP) ::  pk_aa11, pk_aa12, pk_aa22
      real(DP) ::  pk_aa23, pk_aa33
      real(DP) :: pk_A1, pk_A2, pk_A3
      real(DP) :: fact00,fact10,fact20,fact30,fact40
      real(DP) :: fact02,fact12,fact22,fact32,fact42
      real(DP) :: fact04,fact14,fact24,fact34,fact44
      real(DP) :: ptt,pdt,pdd

      kmin = ak(1)
      kmax = ak(nk) 

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

      pk_A11 = 0.0d0 
      pk_A12 = 0.0d0
      pk_A22 = 0.0d0
      pk_A23 = 0.0d0 
      pk_A33 = 0.0d0
      pk_tA11 = 0.0d0
      pk_tA12 = 0.0d0
      pk_tA22 = 0.0d0
      pk_tA23 = 0.0d0
      pk_tA33 = 0.0d0
      pk_aa11 = 0.0d0
      pk_aa12 = 0.0d0
      pk_aa22 = 0.0d0
      pk_aa23 = 0.0d0
      pk_aa33 = 0.0d0

      k = ak(ik)

      xmin = kmin / k
      xmax = kmax / k

!     ////// Gauss-Legendre integration //////  c

      if(k.lt.0.2) isub =200 
      if(k.ge.0.2) isub =0 
                
      call gaulegf(log(xmin),log(xmax),xx,ww,ixmax-isub)

      do ix=1, ixmax-isub
        xx(ix)= dexp(xx(ix))
        pk_B111 = pk_B111+ww(ix)*integ_fp_B(1,xx(ix),k,xmin,xmax)
        pk_B112 = pk_B112+ww(ix)*integ_fp_B(2,xx(ix),k,xmin,xmax)
        pk_B121 = pk_B121+ww(ix)*integ_fp_B(3,xx(ix),k,xmin,xmax)
        pk_B122 = pk_B122+ww(ix)*integ_fp_B(4,xx(ix),k,xmin,xmax)
        pk_B211 = pk_B211+ww(ix)*integ_fp_B(5,xx(ix),k,xmin,xmax)
        pk_B212 = pk_B212+ww(ix)*integ_fp_B(6,xx(ix),k,xmin,xmax)
        pk_B221 = pk_B221+ww(ix)*integ_fp_B(7,xx(ix),k,xmin,xmax)
        pk_B222 = pk_B222+ww(ix)*integ_fp_B(8,xx(ix),k,xmin,xmax)
        pk_B312 = pk_B312+ww(ix)*integ_fp_B(9,xx(ix),k,xmin,xmax)
        pk_B321 = pk_B321+ww(ix)*integ_fp_B(10,xx(ix),k,xmin,xmax)
        pk_B322 = pk_B322+ww(ix)*integ_fp_B(11,xx(ix),k,xmin,xmax)
        pk_B422 = pk_B422+ww(ix)*integ_fp_B(12,xx(ix),k,xmin,xmax)
        pk_A11 = pk_A11+ww(ix)*integ_fp_A(1,xx(ix),k,xmin,xmax)
        pk_A12 = pk_A12+ww(ix)*integ_fp_A(2,xx(ix),k,xmin,xmax)
        pk_A22 = pk_A22+ww(ix)*integ_fp_A(3,xx(ix),k,xmin,xmax)
        pk_A23 = pk_A23+ww(ix)*integ_fp_A(4,xx(ix),k,xmin,xmax)
        pk_A33 = pk_A33+ww(ix)*integ_fp_A(5,xx(ix),k,xmin,xmax)     
        pk_tA11 = pk_tA11+ww(ix)*integ_fp_A(6,xx(ix),k,xmin,xmax)
        pk_tA12 = pk_tA12+ww(ix)*integ_fp_A(7,xx(ix),k,xmin,xmax)
        pk_tA22 = pk_tA22+ww(ix)*integ_fp_A(8,xx(ix),k,xmin,xmax)
        pk_tA23 = pk_tA23+ww(ix)*integ_fp_A(9,xx(ix),k,xmin,xmax)
        pk_tA33 = pk_tA33+ww(ix)*integ_fp_A(10,xx(ix),k,xmin,xmax)
        pk_aa11 = pk_aa11+ww(ix)*integ_fp_A(11,xx(ix),k,xmin,xmax)
        pk_aa12 = pk_aa12+ww(ix)*integ_fp_A(12,xx(ix),k,xmin,xmax)
        pk_aa22 = pk_aa22+ww(ix)*integ_fp_A(13,xx(ix),k,xmin,xmax)
        pk_aa23 = pk_aa23+ww(ix)*integ_fp_A(14,xx(ix),k,xmin,xmax)
        pk_aa33 = pk_aa33+ww(ix)*integ_fp_A(15,xx(ix),k,xmin,xmax)
      enddo

      kfact = k**3 / (2.*pi)**2

      ! I am a little unsure about the factor of 2 here 

      pk_B111 = 2.d0 * pk_B111 * kfact
      pk_B112 = - 2.d0 * pk_B112 * kfact
      pk_B121 = - 2.d0 * pk_B121 * kfact
      pk_B122 = 2.d0 * pk_B122 * kfact
      pk_B211 = 2.d0 * pk_B211 * kfact
      pk_B212 = - 2.d0 * pk_B212 * kfact
      pk_B221 = - 2.d0 * pk_B221 * kfact
      pk_B222 = 2.d0 * pk_B222 * kfact
      pk_B312 = - 2.d0 * pk_B312 * kfact
      pk_B321 = - 2.d0 * pk_B321 * kfact
      pk_B322 = 2.d0 * pk_B322 * kfact
      pk_B422 = 2.d0 * pk_B422 * kfact
      pk_A11 = 2.d0 * pk_A11 * kfact
      pk_A12 = 2.d0 * pk_A12 * kfact
      pk_A22 = 2.d0 * pk_A22 * kfact
      pk_A23 = 2.d0 * pk_A23 * kfact
      pk_A33 = 2.d0 * pk_A33 * kfact
      pk_tA11 = 2.d0 * pk_tA11 * kfact
      pk_tA12 = 2.d0 * pk_tA12 * kfact
      pk_tA22 = 2.d0 * pk_tA22 * kfact
      pk_tA23 = 2.d0 * pk_tA23 * kfact
      pk_tA33 = 2.d0 * pk_tA33 * kfact
      pk_aa11 = 2.d0 * pk_aa11 * kfact
      pk_aa12 = 2.d0 * pk_aa12 * kfact
      pk_aa22 = 2.d0 * pk_aa22 * kfact
      pk_aa23 = 2.d0 * pk_aa23 * kfact
      pk_aa33 = 2.d0 * pk_aa33 * kfact

      !write(6,'(i4,1p4e18.10)') ik,k,pk(ik),pk_B111(ik),pk_B112(ik)
         
      pk_B1 = ff**2*pk_B111 + ff**3*pk_B112 + ff**3*pk_B121 + ff**4*pk_B122
      pk_B2 = ff**2*pk_B211 + ff**3*pk_B212 + ff**3*pk_B221 + ff**4*pk_B222
      pk_B3 =                 ff**3*pk_B312 + ff**3*pk_B321 + ff**4*pk_B322
      pk_B4 =                                                 ff**4*pk_B422

      pk_A1 = ff*(pk_A11+pk_tA11+pk_aa11) + ff**2*(pk_A12+pk_tA12+pk_aa12) 
      pk_A2 =                               ff**2*(pk_A22+pk_tA22+pk_aa22) + ff**3*(pk_A23+pk_tA23+pk_aa23) 
      pk_A3 =                                                                ff**3*(pk_A33+pk_tA33+pk_aa33) 

      alpha = (k*ff*sigmav)**2.0

      fact00 = fact(0,0,alpha) 
      fact10 = fact(1,0,alpha) 
      fact20 = fact(2,0,alpha)
      fact30 = fact(3,0,alpha)
      fact40 = fact(4,0,alpha)

      fact02 = fact(0,2,alpha) 
      fact12 = fact(1,2,alpha) 
      fact22 = fact(2,2,alpha)
      fact32 = fact(3,2,alpha)
      fact42 = fact(4,2,alpha)

      fact04 = fact(0,4,alpha) 
      fact14 = fact(1,4,alpha) 
      fact24 = fact(2,4,alpha)
      fact34 = fact(3,4,alpha)
      fact44 = fact(4,4,alpha)

      pk_val = find_pk(k)
      
      pdd = pk_val*b1**2
      pdt = pk_val*b1
      ptt = pk_val

      pk0dd(ik) = fact00*pdd
      pk2dd(ik) = fact02*pdd
      pk4dd(ik) = fact04*pdd

      pk0dt(ik) = 2*ff*fact10*pdt
      pk2dt(ik) = 2*ff*fact12*pdt
      pk4dt(ik) = 2*ff*fact14*pdt

      pk0tt(ik) = ff**2*fact20*ptt
      pk2tt(ik) = ff**2*fact22*ptt
      pk4tt(ik) = ff**2*fact24*ptt

      !write(*,*) k,pdd*(1.0+2.0/3.0*ff/b1+1.0/5.0*(ff/b1)**2),pk0dd(ik)+pk0dt(ik)+pk0tt(ik)
      !write(*,*) k,pdd*(4.0/3.0*ff/b1+4.0/7.0*(ff/b1)**2),pk2dd(ik)+pk2dt(ik)+pk2tt(ik)

      ! Note A(k,mu,f) ~ kmuf (eqn 19) and B(k,mu,f) ~ (kmuf)^2 (eqn 20 of Taruya 2010)
      ! Writing as A(k,mu,beta) so f->beta means B ~ b_1^2 and A ~ b_1

      pk0corr_B(ik) = (fact10 * pk_B1 + fact20 * pk_B2 + fact30 * pk_B3 + fact40 * pk_B4)*b1**2.0
      pk2corr_B(ik) = (fact12 * pk_B1 + fact22 * pk_B2 + fact32 * pk_B3 + fact42 * pk_B4)*b1**2.0
      pk4corr_B(ik) = (fact14 * pk_B1 + fact24 * pk_B2 + fact34 * pk_B3 + fact44 * pk_B4)*b1**2.0

      pk0corr_A(ik) = (fact10 * pk_A1 + fact20 * pk_A2 + fact30 * pk_A3)*b1**1.0
      pk2corr_A(ik) = (fact12 * pk_A1 + fact22 * pk_A2 + fact32 * pk_A3)*b1**1.0
      pk4corr_A(ik) = (fact14 * pk_A1 + fact24 * pk_A2 + fact34 * pk_A3)*b1**1.0
         
      !write(6,'(i4,1p100e11.3)') ik,k,pk0dd(ik),pk2dd(ik),pk4dd(ik),&
      ! pk0dt(ik),pk2dt(ik),pk4dt(ik),pk0tt(ik),pk2tt(ik),pk4tt(ik), &
      !  pk0corr_B(ik),pk2corr_B(ik),pk4corr_B(ik),pk0corr_A(ik),pk2corr_A(ik),pk4corr_A(ik)
         
    end subroutine calc_correction

! ******************************************************* 

    subroutine calc_pkred(filename)
  
!     Summing up all contributions to redshift P(k) in PT
!     and calculating monopole, quadrupole and hexadecapole
!     spectra

      character(len=200) filename
      integer ik
      
      !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDUlE(DYNAMIC) &
      !$OMP & PRIVATE(ik)
      do ik=1,nk
         call calc_correction(ik)
      end do
      !$OMP END PARAllEl DO
      
      open(unit=11,file=trim(filename)//'_l0.dat',status='unknown')
      do ik=1,nk
         write(11,'(1p100e11.3)') ak(ik),pk0dd(ik),pk0dt(ik),pk0tt(ik),pk0corr_A(ik),pk0corr_B(ik)
      end do
      close(11)

      open(unit=11,file=trim(filename)//'_l2.dat',status='unknown')
      do ik=1,nk
         write(11,'(1p100e11.3)') ak(ik),pk2dd(ik),pk2dt(ik),pk2tt(ik),pk2corr_A(ik),pk2corr_B(ik)
      end do
      close(11)

      open(unit=11,file=trim(filename)//'_l4.dat',status='unknown')
      do ik=1,nk
         write(11,'(1p100e11.3)') ak(ik),pk4dd(ik),pk4dt(ik),pk4tt(ik),pk4corr_A(ik),pk4corr_B(ik)
      end do
      close(11)


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


