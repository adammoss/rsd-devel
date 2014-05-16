c ********************************************************** c
c
      program read_stfile
c
c               Time-stamp: <2012-07-31 13:54:31 ataruya>    c
c ********************************************************** c
c
c     Reading the diagrams data created by RegPT code (default file name: 
c     st_pt.dat), this supplemental code provides RegPT, standard PT, and 
c     LRT results of matter power spectrum at any redshift (in flat wCDM 
c     model). 
c
c     How to setup and how to run:
c
c         Just compile this file with fortran compiler, ifort or gfortran:
c
c             ifort [or gfortran] read_stfile.f -o read_stfile.exe
c
c         To run the code, just type,  ./read_stfile.exe
c         The code will ask you to type the file name of diagram data, and 
c         you may just follow the instructions. The computed results are
c         saved to 'pk_total.dat' for power spectrum, and 'xi_total.dat' 
c         for correlation function (optional). 
c
c
c     Input files: 
c         PT diagram data (created by RegPT code with option -direct)
c
c     Output format for 'pk_total.dat':
c
c     k, P_no-wiggle(k), P_lin(k), P_RegPT1loop(k), P_RegPT2loop(k), 
c           P_SPT1loop(k), P_SPT2loop(k), P_LRT1loop(k), P_LRT2loop(k)
c     
c     Here, k is in units of [h/Mpc], and all of P(k) are in units of [(Mpc/h)^3].
c     
c     Output format for 'xi_total.dat' (optional):
c
c     r, Xi_lin(r), Xi_RegPT1loop(r), Xi_RegPT2loop(r), 
c           Xi_LRT1loop(r), Xi_LRT2loop(r)
c
c     where r is in units of [Mpc/h], and all of Xi(r) are dimensionless.
c
      implicit none
      integer ik, ikmax, ik_max
      integer itype, ir, ir_max, iparams
      parameter(ikmax=5000)
      real*8   ak(ikmax), pk_lin(ikmax), pk_EH(ikmax), dummy
      real*8   G1a_1loop(ikmax), G1a_2loop(ikmax)
      real*8   G1b_1loop(ikmax), G1b_2loop(ikmax)
      real*8   pkcorr_G2_tree_tree(ikmax), pkcorr_G2_tree_1loop(ikmax)
      real*8   pkcorr_G2_1loop_1loop(ikmax), pkcorr_G3_tree(ikmax)
      real*8   sigmav, exp_factor
      real*8   G1a_reg, G1b_reg, pkcorr_G1, pkcorr_G2, pkcorr_G3
      real*8   pkRegPT_1loop(ikmax), pkRegPT_2loop(ikmax)
      real*8   pk22, pk13, pk33, pk24, pk15
      real*8   pkSPT_1loop(ikmax), pkSPT_2loop(ikmax)
      real*8   A11, A22, A13, exponent_A, exponent_B
      real*8   pkLRT_1loop(ikmax), pkLRT_2loop(ikmax)
      real*8   ar(ikmax), xi_lin(ikmax)
      real*8   xiRegPT_1loop(ikmax), xiRegPT_2loop(ikmax)
      real*8   xiLRT_1loop(ikmax), xiLRT_2loop(ikmax)
      real*8   kmin, kmax, dd, pi
      common /pk_data/ ak, pk_lin, pkRegPT_1loop, pkRegPT_2loop, 
     &     pkLRT_1loop, pkLRT_2loop, ik_max
      common /growth_factor/ dd
      character istfile*60, infile*50, outfile*50
      pi = 4.d0 * datan(1.d0)
c     ---------------------------------------------------
c
c     ////// File name of output data for RegPTmain  //////
c
      write(6,*) ' Type file name of PT diagram data '
      read(5,*) istfile
c
c     ////// Read pkcorr_Gamma1.dat  //////
c
      open(9, file=istfile, status='unknown')
      do ik=1, ikmax
         read(9,*,END=5) ak(ik), pk_EH(ik), pk_lin(ik), 
     &        G1a_1loop(ik), G1a_2loop(ik), 
     &        pkcorr_G2_tree_tree(ik), pkcorr_G2_tree_1loop(ik), 
     &        pkcorr_G2_1loop_1loop(ik), pkcorr_G3_tree(ik)
         G1b_1loop(ik) = G1a_1loop(ik)
         G1b_2loop(ik) = G1a_2loop(ik)
      enddo
c
 5    ik_max = ik-1
      close(9)
      kmin = ak(1)
      kmax = ak(ik_max)
c
c     ///// linear P(k) of no-wiggles: Eisenstein & Hu (1998) ///// c
c
      call get_growth
c
      write(6,*)
      write(6,*) 'D(z)=',dd
      write(6,*)
c
c     ////// sigma_v  //////
c
ccc      call calc_sigmav(sigmav)
c
      call calc_A11_A22_A13(kmin, kmax, A11, A22, A13)
c
      A11 = A11 * dd**2
      A22 = A22 * dd**4
      A13 = A13 * dd**4
c
      do ik=1, ik_max
c
         call calc_running_sigmav(kmin, kmax, ak(ik), sigmav)
c
         exp_factor = 0.5d0 * (ak(ik)*sigmav*dd)**2 
c
c     ////// Power spectrum at two-loop order ////// 
c
         if(exp_factor.gt.100.d0) then
c
c     to avoid underflow
            pkRegPT_2loop(ik) = 0.d0
            pkRegPT_1loop(ik) = 0.d0
            pkLRT_2loop(ik) = 0.d0
            pkLRT_1loop(ik) = 0.d0
c
         else
c
c     ----- RegPT ----
            G1a_reg = 1.d0 + exp_factor + 0.5*exp_factor**2 + 
     &           dd**2*G1a_1loop(ik) * (1.d0 + exp_factor) + 
     &           dd**4*G1a_2loop(ik) 
            G1b_reg = 1.d0 + exp_factor + 0.5*exp_factor**2 + 
     &           dd**2*G1b_1loop(ik) * (1.d0 + exp_factor) + 
     &           dd**4*G1b_2loop(ik) 
            G1a_reg = G1a_reg * dd * dexp(-exp_factor)
            G1b_reg = G1b_reg * dd * dexp(-exp_factor)
c
            pkcorr_G1 = G1a_reg*G1b_reg * pk_lin(ik)
c
            pkcorr_G2 = pkcorr_G2_tree_tree(ik) * 
     &           (1.d0 + exp_factor)**2 + pkcorr_G2_tree_1loop(ik) * 
     &           dd**2*(1.d0 + exp_factor) + pkcorr_G2_1loop_1loop(ik) * 
     &           dd**4 
            pkcorr_G2 = pkcorr_G2 * dd**4 * dexp(-2.d0 * exp_factor)
c
            pkcorr_G3 = pkcorr_G3_tree(ik) * dd**6 
     &           * dexp(-2.d0 * exp_factor)
c
            pkRegPT_2loop(ik) = pkcorr_G1 + pkcorr_G2 + pkcorr_G3 ! RegPT
c
c     ----- Standard PT ----- 
            pk22 = pkcorr_G2_tree_tree(ik)
            pk13 = 2.d0*pk_lin(ik) * G1a_1loop(ik)
            pk33 = pk_lin(ik) * G1a_1loop(ik)**2 + pkcorr_G3_tree(ik)
            pk24 = pkcorr_G2_tree_1loop(ik)
            pk15 = 2.d0*pk_lin(ik) * G1a_2loop(ik)
c
            pkSPT_2loop(ik) = dd**2*pk_lin(ik)
     &           + dd**4*(pk22 + pk13) + dd**6*(pk33 + pk24 + pk15)
c
c     ----- Lagrangian Ressumation Theory (LRT) -----
            exponent_A = ak(ik)**2/(6.d0*pi**2) 
     &           * ( A11 + A22 + 2.d0*A13)
            exponent_B = ak(ik)**2/(6.d0*pi**2) * A11
c
            pkLRT_2loop(ik) = pkSPT_2loop(ik) + dd**2*pk_lin(ik) *
     &           (exponent_A + 0.5d0*exponent_B**2) + 
     &           exponent_B * dd**4*(pk22 + pk13) 
            pkLRT_2loop(ik) = pkLRT_2loop(ik) * dexp(-exponent_A)
c
c     ////// Power spectrum at one-loop order //////
c
c     ----- RegPT -----
            G1a_reg = dd * ( 1.d0 + exp_factor + dd**2*G1a_1loop(ik) )
     &           * dexp(-exp_factor)
            G1b_reg = dd * ( 1.d0 + exp_factor + dd**2*G1b_1loop(ik) )
     &           * dexp(-exp_factor)
c
            pkcorr_G1 = G1a_reg*G1b_reg * pk_lin(ik)
c
            pkcorr_G2 = pkcorr_G2_tree_tree(ik) * dd**4 
     &           * dexp(-2.d0 * exp_factor)
c
            pkRegPT_1loop(ik) = pkcorr_G1 + pkcorr_G2
c
c     ----- Standard PT ----- 
            pkSPT_1loop(ik) = dd**2*pk_lin(ik)
     &           + dd**4*(pk22 + pk13)
c
c     ----- LRT -----
            pkLRT_1loop(ik) = pkSPT_1loop(ik) +
     &           exponent_B * dd**2*pk_lin(ik)
            pkLRT_1loop(ik) = pkLRT_1loop(ik) * dexp(-exponent_B)
c
         endif
c
      enddo
c
      close(9)
c
c     //////// Save power spectrum results ////// c
c
      open(9, file='pk_total.dat', status='unknown')
c
      do ik=1, ik_max
         write(9,'(1p9e18.10)') ak(ik), 
     &        dd**2*pk_EH(ik), dd**2*pk_lin(ik), 
     &        pkRegPT_1loop(ik), pkRegPT_2loop(ik), 
     &        pkSPT_1loop(ik), pkSPT_2loop(ik), 
     &        pkLRT_1loop(ik), pkLRT_2loop(ik) 
      enddo
c
      close(9)
c
c
      write(6,*) 
      write(6,*) 'Output file for power spectrum --->  pk_total.dat '
      write(6,*) 
c
c     //////// Save correlation function results (optional) ////// c      
c
      write(6,*) 'Save Xi(r) results for RegPT and LRT ?  Yes[0], No[1]'
      read(5,*) itype
      if(itype.eq.1) stop
c
      if(kmax .lt. 5.d0) then
         write(6,*) 'Cutoff k_max is too small. Convergence of    '
         write(6,*) 'the Fourier integral might not be guaranteed.' 
      endif
c
      call calc_xi(ar, kmin, kmax, ik_max, 
     &     xi_lin, xiRegPT_1loop, xiRegPT_2loop, 
     &     xiLRT_1loop, xiLRT_2loop, ir_max)
c
      open(11, file='xi_total.dat', status='unknown')
c
      do ir = 1, ir_max
         write(11,'(1p6e18.10)') ar(ir), dd**2*xi_lin(ir), 
     &        xiRegPT_1loop(ir), xiRegPT_2loop(ir), 
     &        xiRegPT_1loop(ir), xiRegPT_2loop(ir) 
      enddo
c     
      close(11)
c     
      write(6,*) 
      write(6,*) 'Output file for correlation function'
      write(6,*) '                              --->  xi_total.dat '
      write(6,*) 
c
      end
c
c ******************************************************* c
c
      subroutine calc_running_sigmav(kmin, kmax, k, sigmav)
c
c ******************************************************* c
c
c     sigmav^2 = int_0^{k/2} dq P0(q)/(6*pi^2)
c
      implicit none
c
      integer ik, ik_max, ix, ixmax
      parameter(ixmax=500)
      real*8  qq(ixmax), ww(ixmax)
      real*8  k, q, kmin, kmax, k_cut, pklin, sigmav, sigma_v2, pi
      pi = 4.d0 * datan(1.d0)
c     --------------------------------------------------------
c
      sigma_v2 = 0.d0
      k_cut = k/2.d0
c
      if(k_cut.le.kmin) goto 1      ! out of k-range, just skip
      if(k_cut.ge.kmax) k_cut = kmax  
c
      call gauleg(log(kmin), log(k_cut), qq, ww, ixmax)
c
      do ix=1, ixmax
         q = dexp(qq(ix))
         call find_pk(0, q, pklin)
         sigma_v2 = sigma_v2 + ww(ix) * pklin * q 
      enddo
c
 1    sigma_v2 = sigma_v2 / (6.d0*pi**2)
      sigmav = dsqrt(sigma_v2)
c
      end
c
c ******************************************************* c
c
      subroutine calc_A11_A22_A13(kmin, kmax, A11, A22, A13)
c
c ******************************************************* c
c
      implicit none
      integer ip1, ip2, ixmax
      parameter(ixmax=600)
      real*8  pi
      real*8  A11, A22, A13
      real*8  kmin, kmax
      real*8  p1, p2, pklin1, pklin2
      real*8  ww(ixmax), kk(ixmax)
      real*8  int_ss, ss, kernel_LRT
      pi = 4.d0 * datan(1.d0)
c     ---------------------------------------------------
c
      call gauleg(log(kmin),log(kmax),kk,ww,ixmax)
c
      A11 = 0.d0
      ss = 0.d0
c
c     ////// Gauss-Legendre integration for A11, A22 //////  c
c
      do 10 ip1=1, ixmax
c
         p1 = dexp(kk(ip1))
         call find_pk(0, p1, pklin1)
         A11 = A11 + ww(ip1) * pklin1 * p1
c
         int_ss = 0.d0
c
         do ip2=1, ixmax
c
            p2 = dexp(kk(ip2))
            call find_pk(0, p2, pklin2) 
            int_ss = int_ss + ww(ip2)
     &           * kernel_LRT(p1,p2) * pklin2 * p2**2
         enddo
c
         ss = ss + int_ss * ww(ip1) * pklin1 * p1**2
c
 10   continue
c
      A22 = ss * 9.d0/ (392.d0*pi*pi)
      A13 = ss * 5.d0/ (84.d0*pi*pi)
c
      end
c
c ******************************************************* c
c
      function  kernel_LRT(p1, p2)
c
c ******************************************************* c
c
      implicit none
      real*8 kernel_LRT, y, s, p1, p2
c     -------------------------------------
c
      kernel_LRT = 0.0
c
      y = ( p1/p2 + p2/p1 ) / 2.d0
c
      if(y.lt.0.) then
         write(6,*) ' warning: kernel_LRT negative !! : y=',y
         kernel_LRT = 0.0
      elseif(y .lt. 5.d-2) then
         kernel_LRT = (8./3.) * (y - y**3 + y**5/5.)
      elseif(y .gt. 5.d2) then
         s = 1.d0 / y
         kernel_LRT = 8./315. *s*(21. + 3.*s*s + s**4)
      elseif(y.ge.0.99 .and. y.le.1.01) then
         if(abs(y-1.).le.1.d-3) then
            kernel_LRT = 2./3. - 4./3.*(y-1.) 
         else
            kernel_LRT = 2./3. - 4./3.*(y-1.) + 
     &           ( -3.+dlog(4.d0)-2.*dlog(dabs(y-1.)) )*(y-1.)**2
         endif
      else 
         kernel_LRT = (5./3.)*y - y**3 
     &        - 0.5*(y**2-1.)**2 * dlog(dabs((y-1.)/(y+1.)))
      endif
c
      end
c
c ************************************************************
c
      subroutine calc_xi(ar, kmin, kmax, ik_max, 
     &     xi_lin, xiRegPT_1loop, xiRegPT_2loop, 
     &     xiLRT_1loop, xiLRT_2loop, ir_max)
c
c ************************************************************
c
      implicit none
      integer ik, ir, ix, ikmax, ik_max, ixmax, irmax, ir_max
      parameter(ikmax=5000, ixmax=10000, irmax=501)
      real*8  kk(ixmax), wk(ixmax)
      real*8  ar(ikmax), xi_lin(ikmax)
      real*8  xiRegPT_1loop(ikmax), xiRegPT_2loop(ikmax)
      real*8  xiLRT_1loop(ikmax), xiLRT_2loop(ikmax)
      real*8  rmin, rmax, kmin, kmax, k1, k2
      real*8  pklin1,  pkRegPT_1loop1, pkRegPT_2loop1, sinc, pi
      real*8  pklin2,  pkRegPT_1loop2, pkRegPT_2loop2
      real*8  pkLRT_1loop1, pkLRT_2loop1
      real*8  pkLRT_1loop2, pkLRT_2loop2
      parameter(rmin=30.d0, rmax=200.d0)
      pi = 4.d0 * datan(1.d0)
c     ------------------------------------------------
c
      ir_max = irmax
c
      do ir = 1, irmax 
         ar(ir) = rmin * (rmax/rmin)**(dble(ir-1)/dble(irmax-1))
         xi_lin(ir) = 0.0d0
         xiRegPT_1loop(ir) = 0.0d0
         xiRegPT_2loop(ir) = 0.0d0
         xiLRT_1loop(ir) = 0.0d0
         xiLRT_2loop(ir) = 0.0d0
c
         k1 = kmin
         call find_pk(0, k1, pklin1)
         call find_pk(1, k1, pkRegPT_1loop1)
         call find_pk(2, k1, pkRegPT_2loop1)
         call find_pk(3, k1, pkLRT_1loop1)
         call find_pk(4, k1, pkLRT_2loop1)
c
         do ix = 1, ixmax-1
c
            k2 = kmin * (kmax/kmin) ** (dble(ix)/dble(ixmax-1))
c
            call find_pk(0, k2, pklin2)
            call find_pk(1, k2, pkRegPT_1loop2)
            call find_pk(2, k2, pkRegPT_2loop2)
            call find_pk(3, k2, pkLRT_1loop2)
            call find_pk(4, k2, pkLRT_2loop2)
c
            xi_lin(ir) = xi_lin(ir) + 
     &           ( k1*k1 * sinc(k1*ar(ir)) * pklin1 
     &           + k2*k2 * sinc(k2*ar(ir)) * pklin2 ) 
     &           * (k2 - k1)/2.d0
            xiRegPT_1loop(ir) = xiRegPT_1loop(ir) + 
     &           ( k1*k1 * sinc(k1*ar(ir)) * pkRegPT_1loop1 
     &           + k2*k2 * sinc(k2*ar(ir)) * pkRegPT_1loop2 ) 
     &           * (k2 - k1)/2.d0
            xiRegPT_2loop(ir) = xiRegPT_2loop(ir) + 
     &           ( k1*k1 * sinc(k1*ar(ir)) * pkRegPT_2loop1 
     &           + k2*k2 * sinc(k2*ar(ir)) * pkRegPT_2loop2 ) 
     &           * (k2 - k1)/2.d0
            xiLRT_1loop(ir) = xiLRT_1loop(ir) + 
     &           ( k1*k1 * sinc(k1*ar(ir)) * pkLRT_1loop1 
     &           + k2*k2 * sinc(k2*ar(ir)) * pkLRT_1loop2 ) 
     &           * (k2 - k1)/2.d0
            xiLRT_2loop(ir) = xiLRT_2loop(ir) + 
     &           ( k1*k1 * sinc(k1*ar(ir)) * pkLRT_2loop1 
     &           + k2*k2 * sinc(k2*ar(ir)) * pkLRT_2loop2 ) 
     &           * (k2 - k1)/2.d0
c
            pklin1 = pklin2
            pkRegPT_1loop1 = pkRegPT_1loop2
            pkRegPT_2loop1 = pkRegPT_2loop2
            pkLRT_1loop1 = pkLRT_1loop2
            pkLRT_2loop1 = pkLRT_2loop2
            k1 = k2
c
         enddo
c
         xi_lin(ir)   = xi_lin(ir)   / (2.d0*pi*pi)
         xiRegPT_1loop(ir) = xiRegPT_1loop(ir) / (2.d0*pi*pi)
         xiRegPT_2loop(ir) = xiRegPT_2loop(ir) / (2.d0*pi*pi)
         xiLRT_1loop(ir) = xiLRT_1loop(ir) / (2.d0*pi*pi)
         xiLRT_2loop(ir) = xiLRT_2loop(ir) / (2.d0*pi*pi)
c     
      enddo
c     
      end
c
c ******************************************************* c
c
      function sinc(x)
c
c ******************************************************* c
c
      real*8 x, sinc
c
      if(dabs(x).le.1.d-3) then 
         sinc = 1 - x**2/6.d0 + x**4/120.d0
      else
         sinc = sin(x)/x
      endif
c
      end
c
c ******************************************************* c
c
      subroutine find_pk(ipk, kk, pk)
c
c ******************************************************* c
c
      implicit none
      integer ik_max, ikmax, ipk
      integer j, jmin, jmax
      parameter(ikmax=5000)
      real*8 ak(ikmax), pk_lin(ikmax)
      real*8 pkRegPT_1loop(ikmax), pkRegPT_2loop(ikmax)
      real*8 pkLRT_1loop(ikmax), pkLRT_2loop(ikmax)
      real*8 kk, s, ds, pk
      common /pk_data/ ak, pk_lin, pkRegPT_1loop, pkRegPT_2loop,
     &     pkLRT_1loop, pkLRT_2loop, ik_max
c     -------------------------------------------
c
      call hunt(ak, ik_max, kk, j)
c
      jmin = j - 2
      jmax = j + 2
      if(jmin.lt.1) jmin = 1
      if(jmax.ge.ik_max) jmax = ik_max
c
      if(ipk.eq.0) 
     &     call polint(ak(jmin),pk_lin(jmin),jmax-jmin+1,kk,s,ds)
      if(ipk.eq.1) 
     &     call polint(ak(jmin),pkRegPT_1loop(jmin),jmax-jmin+1,kk,s,ds)
      if(ipk.eq.2) 
     &     call polint(ak(jmin),pkRegPT_2loop(jmin),jmax-jmin+1,kk,s,ds)
      if(ipk.eq.3) 
     &     call polint(ak(jmin),pkLRT_1loop(jmin),jmax-jmin+1,kk,s,ds)
      if(ipk.eq.4) 
     &     call polint(ak(jmin),pkLRT_2loop(jmin),jmax-jmin+1,kk,s,ds)
      pk = s
c      
      end
c
c ************************************************ c
c
      SUBROUTINE hunt(xx,n,x,jlo)
c
c ************************************************ c
c
      implicit none
      INTEGER jlo,n
      REAL*8 x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END
c
c ******************************************************* c
c
      function growth(zred, Omega_m, Omega_v, w_de)
c
c ******************************************************* c
c
c     Linear growth factor for flat cosmology
c
      implicit none
      real*8 growth, zred, Omega_m, Omega_v, a, b, c
      real*8 zred1, zi, zz, gz, g0
      real*8 w_de
c     --------------------------------------------
c
      a = -1.d0 / (3.d0 * w_de)
      b = (w_de - 1.d0)/ (2.d0 * w_de)
      c = 1.d0 - 5.d0 / (6.d0 * w_de)
c
      zred1 = 1.d0 + zred
      zi = - Omega_v / Omega_m
      zz = zi * zred1**(3.d0*w_de) 
c
      call HYGFX(a,b,c,zi, g0)  
      call HYGFX(a,b,c,zz, gz)  
c
      growth = (gz/g0) / zred1
c
      end
c
c ************************************************ c
c
      subroutine get_growth
c
c ************************************************ c
c
      implicit none
      integer iparams
      real*8  omega0, w_de, growth, dd, zred
      common /growth_factor/ dd
c     ----------------------------------
c
      write(6,*)
      write(6,*) 'Output redshift ?'
      read(5,*) zred
c
c ///// cosmological parameters (default, wmap5) ///// c
c
      omega0 = 0.279d0
      w_de = -1.d0
c
 5    write(6,*) 'Cosmological parameters used for growth factor'
      write(6,*)
      write(6,*) '[1] omega_m =',omega0
      write(6,*) '[2]    w_de =',w_de
c
      write(6,*)
      write(6,*) 'change cosmological parameter? [1-5] or n[0]'
      write(6,*)
      read(5,*)  iparams
      if(iparams.eq.0) goto 10
      if(iparams.eq.1) then
         write(6,*) 'type omega_m'
         read(5,*) omega0
      elseif(iparams.eq.2) then
         write(6,*) 'type w_de'
         read(5,*) w_de
      endif
      goto 5
c
 10   dd = growth(zred, omega0, 1.d0-omega0, w_de)
c
      end
c
c ************************************************ c
c
      function var(k) 
c
c ************************************************ c
c
      implicit none
      real*8  var, k, x, w_th, Pk_lin_EH
      real*8  pi, r_th
      common /R_tophat/ r_th
c     -----------------------------------
c
      pi = 4.d0 * datan(1.d0)
c
      x = k * r_th
      w_th = 3.d0 * (sin(x)-x*cos(x)) / x**3 
      var = k**2 * w_th**2 * Pk_lin_EH(k) / (2.d0*pi**2)
c
      end
c
c ************************************************ c
c
      function Pk_lin_EH(k)
c
c ************************************************ c
c
c     compute un-normalized linear P(k) 
c     based on eq.(29) of Eisenstein & Hu (1998)
c     (no-wiggle approximation)  
c
      implicit none
      real*8 k, ss, alpha_gam, theta_cmb
      real*8 gamma_eff, q, L0, C0, T_EH, Pk_lin_EH
      real*8 omegab, omega0, h, Tcmb, n_s
      common /cosm_param/ omegab, omega0, h, Tcmb, n_s
c     -----------------------------------------
c
c ///// fitting formula for no-wiggle P(k) (Eq.[29] of EH98)
c
      ss = 44.5 * h * dlog( 9.83 / (omega0*h*h) ) / 
     &     dsqrt( 1.d0 + 10.d0 * (omegab*h*h)**0.75 )
      alpha_gam = 1.d0 
     &     - 0.328 * dlog( 431. * omega0*h*h ) * omegab/omega0
     &     + 0.38 * dlog( 22.3 * omega0*h*h ) * (omegab/omega0)**2
      theta_cmb = Tcmb / 2.70 
      gamma_eff = omega0 * h * 
     &     ( alpha_gam + (1.d0 - alpha_gam) / (1.d0 + (0.43*k*ss)**4) )
c
      q = k * theta_cmb**2 / gamma_eff
      L0 = dlog( 2.d0 * dexp(1.d0) + 1.8 * q ) 
      C0 = 14.2 + 731.d0 / ( 1.d0 + 62.5 * q )
c
      T_EH = L0 / (L0 + C0*q*q )
c
      Pk_lin_EH = k ** n_s * T_EH**2 
c
      end
c
c ************************************************************
c
      SUBROUTINE qromo(func,a,b,ss,choose)
c
c ************************************************************
c
      implicit none
      INTEGER JMAX,JMAXP,K,KM
      REAL*8 a,b,func,ss,EPS
      EXTERNAL func,choose
      PARAMETER (EPS=1.d-6, JMAX=14, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        call choose(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=h(j)/9.
11    continue
      write(6,*) 'too many steps in qromo'
      stop
      END
c
c ************************************************************
c
      SUBROUTINE polint(xa,ya,n,x,y,dy)
c
c ************************************************************
c
      implicit none
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) then
             write(6,*) 'failure in polint'
             stop
          endif
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
c
c
C ************************************************************
C
      SUBROUTINE midpnt(func,a,b,s,n)
C
C ************************************************************
C
      implicit none
      INTEGER n
      REAL*8 a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL*8 ddel,del,sum,tnm,x
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END
C
C ************************************************************
C
      SUBROUTINE midinf(funk,aa,bb,s,n)
C
C ************************************************************
C
      implicit none
      INTEGER n
      REAL*8 aa,bb,s,funk
      EXTERNAL funk
      INTEGER it,j
      REAL*8 a,b,ddel,del,sum,tnm,func,x
      func(x)=funk(1./x)/x**2
      b=1./aa
      a=1./bb
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END
c
c
c ******************************************************* c
c
        SUBROUTINE HYGFX(A,B,C,X,HF)
c
c ******************************************************* c
C
C       ====================================================
C       Purpose: Compute hypergeometric function F(a,b,c,x)
C       Input :  a --- Parameter
C                b --- Parameter
C                c --- Parameter, c <> 0,-1,-2,...
C                x --- Argument   ( x < 1 )
C       Output:  HF --- F(a,b,c,x)
C       Routines called:
C            (1) GAMMAX for computing gamma function
C            (2) PSI for computing psi function
C       ====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        LOGICAL L0,L1,L2,L3,L4,L5
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        L0=C.EQ.INT(C).AND.C.LT.0.0
        L1=1.0D0-X.LT.1.0D-15.AND.C-A-B.LE.0.0
        L2=A.EQ.INT(A).AND.A.LT.0.0
        L3=B.EQ.INT(B).AND.B.LT.0.0
        L4=C-A.EQ.INT(C-A).AND.C-A.LE.0.0
        L5=C-B.EQ.INT(C-B).AND.C-B.LE.0.0
        IF (L0.OR.L1) THEN
           WRITE(*,*)'The hypergeometric series is divergent'
           RETURN
        ENDIF
        EPS=1.0D-15
        IF (X.GT.0.95) EPS=1.0D-8
        IF (X.EQ.0.0.OR.A.EQ.0.0.OR.B.EQ.0.0) THEN
           HF=1.0D0
           RETURN
        ELSE IF (1.0D0-X.EQ.EPS.AND.C-A-B.GT.0.0) THEN
           CALL GAMMAX(C,GC)
           CALL GAMMAX(C-A-B,GCAB)
           CALL GAMMAX(C-A,GCA)
           CALL GAMMAX(C-B,GCB)
           HF=GC*GCAB/(GCA*GCB)
           RETURN
        ELSE IF (1.0D0+X.LE.EPS.AND.DABS(C-A+B-1.0).LE.EPS) THEN
           G0=DSQRT(PI)*2.0D0**(-A)
           CALL GAMMAX(C,G1)
           CALL GAMMAX(1.0D0+A/2.0-B,G2)
           CALL GAMMAX(0.5D0+0.5*A,G3)
           HF=G0*G1/(G2*G3)
           RETURN
        ELSE IF (L2.OR.L3) THEN
           IF (L2) NM=INT(ABS(A))
           IF (L3) NM=INT(ABS(B))
           HF=1.0D0
           R=1.0D0
           DO 10 K=1,NM
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X
10            HF=HF+R
           RETURN
        ELSE IF (L4.OR.L5) THEN
           IF (L4) NM=INT(ABS(C-A))
           IF (L5) NM=INT(ABS(C-B))
           HF=1.0D0
           R=1.0D0
           DO 15 K=1,NM
              R=R*(C-A+K-1.0D0)*(C-B+K-1.0D0)/(K*(C+K-1.0D0))*X
15            HF=HF+R
           HF=(1.0D0-X)**(C-A-B)*HF
           RETURN
        ENDIF
        AA=A
        BB=B
        X1=X
        IF (X.LT.0.0D0) THEN
           X=X/(X-1.0D0)
           IF (C.GT.A.AND.B.LT.A.AND.B.GT.0.0) THEN
              A=BB
              B=AA
           ENDIF
           B=C-B
        ENDIF
        IF (X.GE.0.75D0) THEN
           GM=0.0D0
           IF (DABS(C-A-B-INT(C-A-B)).LT.1.0D-15) THEN
              M=INT(C-A-B)
              CALL GAMMAX(A,GA)
              CALL GAMMAX(B,GB)
              CALL GAMMAX(C,GC)
              CALL GAMMAX(A+M,GAM)
              CALL GAMMAX(B+M,GBM)
              CALL PSI(A,PA)
              CALL PSI(B,PB)
              IF (M.NE.0) GM=1.0D0
              DO 30 J=1,ABS(M)-1
30               GM=GM*J
              RM=1.0D0
              DO 35 J=1,ABS(M)
35               RM=RM*J
              F0=1.0D0
              R0=1.0D0
              R1=1.0D0
              SP0=0.D0
              SP=0.0D0
              IF (M.GE.0) THEN
                 C0=GM*GC/(GAM*GBM)
                 C1=-GC*(X-1.0D0)**M/(GA*GB*RM)
                 DO 40 K=1,M-1
                    R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(K-M))*(1.0-X)
40                  F0=F0+R0
                 DO 45 K=1,M
45                  SP0=SP0+1.0D0/(A+K-1.0)+1.0/(B+K-1.0)-1.0/K
                 F1=PA+PB+SP0+2.0D0*EL+DLOG(1.0D0-X)
                 DO 55 K=1,250
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0))
                    SM=0.0D0
                    DO 50 J=1,M
50                     SM=SM+(1.0D0-A)/((J+K)*(A+J+K-1.0))+1.0/
     &                    (B+J+K-1.0)
                    RP=PA+PB+2.0D0*EL+SP+SM+DLOG(1.0D0-X)
                    R1=R1*(A+M+K-1.0D0)*(B+M+K-1.0)/(K*(M+K))*(1.0-X)
                    F1=F1+R1*RP
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 60
55                  HW=F1
60               HF=F0*C0+F1*C1
              ELSE IF (M.LT.0) THEN
                 M=-M
                 C0=GM*GC/(GA*GB*(1.0D0-X)**M)
                 C1=-(-1)**M*GC/(GAM*GBM*RM)
                 DO 65 K=1,M-1
                    R0=R0*(A-M+K-1.0D0)*(B-M+K-1.0)/(K*(K-M))*(1.0-X)
65                  F0=F0+R0
                 DO 70 K=1,M
70                  SP0=SP0+1.0D0/K
                 F1=PA+PB-SP0+2.0D0*EL+DLOG(1.0D0-X)
                 DO 80 K=1,250
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0))
                    SM=0.0D0
                    DO 75 J=1,M
75                     SM=SM+1.0D0/(J+K)
                    RP=PA+PB+2.0D0*EL+SP-SM+DLOG(1.0D0-X)
                    R1=R1*(A+K-1.0D0)*(B+K-1.0)/(K*(M+K))*(1.0-X)
                    F1=F1+R1*RP
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 85
80                  HW=F1
85               HF=F0*C0+F1*C1
              ENDIF
           ELSE
              CALL GAMMAX(A,GA)
              CALL GAMMAX(B,GB)
              CALL GAMMAX(C,GC)
              CALL GAMMAX(C-A,GCA)
              CALL GAMMAX(C-B,GCB)
              CALL GAMMAX(C-A-B,GCAB)
              CALL GAMMAX(A+B-C,GABC)
              C0=GC*GCAB/(GCA*GCB)
              C1=GC*GABC/(GA*GB)*(1.0D0-X)**(C-A-B)
              HF=0.0D0
              R0=C0
              R1=C1
              DO 90 K=1,250
                 R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(A+B-C+K))*(1.0-X)
                 R1=R1*(C-A+K-1.0D0)*(C-B+K-1.0)/(K*(C-A-B+K))
     &              *(1.0-X)
                 HF=HF+R0+R1
                 IF (DABS(HF-HW).LT.DABS(HF)*EPS) GO TO 95
90               HW=HF
95            HF=HF+C0+C1
           ENDIF
        ELSE
           A0=1.0D0
           IF (C.GT.A.AND.C.LT.2.0D0*A.AND.
     &         C.GT.B.AND.C.LT.2.0D0*B) THEN
              A0=(1.0D0-X)**(C-A-B)
              A=C-A
              B=C-B
           ENDIF
           HF=1.0D0
           R=1.0D0
           DO 100 K=1,250
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X
              HF=HF+R
              IF (DABS(HF-HW).LE.DABS(HF)*EPS) GO TO 105
100           HW=HF
105        HF=A0*HF
        ENDIF
        IF (X1.LT.0.0D0) THEN
           X=X1
           C0=1.0D0/(1.0D0-X)**AA
           HF=C0*HF
        ENDIF
        A=AA
        B=BB
        IF (K.GT.120) WRITE(*,115)
115     FORMAT(1X,'Warning! You should check the accuracy')
        RETURN
        END
c
c ******************************************************* c
c
        SUBROUTINE GAMMAX(X,GA)
c
c ******************************************************* c
C
C       ==================================================
C       Purpose: Compute gamma function â(x)
C       Input :  x  --- Argument of â(x)
C                       ( x is not equal to 0,-1,-2,úúú)
C       Output:  GA --- â(x)
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,
     &          -0.6558780715202538D0, -0.420026350340952D-1,
     &          0.1665386113822915D0,-.421977345555443D-1,
     &          -.96219715278770D-2, .72189432466630D-2,
     &          -.11651675918591D-2, -.2152416741149D-3,
     &          .1280502823882D-3, -.201348547807D-4,
     &          -.12504934821D-5, .11330272320D-5,
     &          -.2056338417D-6, .61160950D-8,
     &          .50020075D-8, -.11812746D-8,
     &          .1043427D-9, .77823D-11,
     &          -.36968D-11, .51D-12,
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END
c
c ******************************************************* c
c
        SUBROUTINE PSI(X,PS)
c
c ******************************************************* c
C
C       ======================================
C       Purpose: Compute Psi function
C       Input :  x  --- Argument of psi(x)
C       Output:  PS --- psi(x)
C       ======================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XA=DABS(X)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        S=0.0D0
        IF (X.EQ.INT(X).AND.X.LE.0.0) THEN
           PS=1.0D+300
           RETURN
        ELSE IF (XA.EQ.INT(XA)) THEN
           N=XA
           DO 10 K=1 ,N-1
10            S=S+1.0D0/K
           PS=-EL+S
        ELSE IF (XA+.5.EQ.INT(XA+.5)) THEN
           N=XA-.5
           DO 20 K=1,N
20            S=S+1.0/(2.0D0*K-1.0D0)
           PS=-EL+2.0D0*S-1.386294361119891D0
        ELSE
           IF (XA.LT.10.0) THEN
              N=10-INT(XA)
              DO 30 K=0,N-1
30               S=S+1.0D0/(XA+K)
              XA=XA+N
           ENDIF
           X2=1.0D0/(XA*XA)
           A1=-.8333333333333D-01
           A2=.83333333333333333D-02
           A3=-.39682539682539683D-02
           A4=.41666666666666667D-02
           A5=-.75757575757575758D-02
           A6=.21092796092796093D-01
           A7=-.83333333333333333D-01
           A8=.4432598039215686D0
           PS=DLOG(XA)-.5D0/XA+X2*(((((((A8*X2+A7)*X2+
     &        A6)*X2+A5)*X2+A4)*X2+A3)*X2+A2)*X2+A1)
           PS=PS-S
        ENDIF
        IF (X.LT.0.0) PS=PS-PI*DCOS(PI*X)/DSIN(PI*X)-1.0D0/X
        RETURN
        END
c
c ************************************************************
c
      SUBROUTINE gauleg(x1,x2,x,w,n)
c
c ************************************************************
c
      INTEGER n
      REAL*8 x1,x2,x(n),w(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END

c
