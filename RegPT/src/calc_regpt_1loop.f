c ********************************************************** c
c
      subroutine calc_regpt_1loop(D_growth, frack)
c
c ********************************************************** c
c     This routine computes the aray pk_main at 1-loop order
c ********************************************************** c
c
      implicit none
c
      integer ik, iz, inum_z, inum_k,ikmax,ik_max,ik_max_fast
      integer frack, ff
      parameter(ikmax=2000)
      real*8  D_growth(100)
      real*8  G1a_1loop, G1b_1loop, pkcorr_G2_tree_tree
      real*8  sigmav2_running, pk_lin, exp_factor, G1a_reg, G1b_reg
      real*8  pkcorr_G1, pkcorr_G2
      real*8  k, ak(ikmax), pk(ikmax), dummy_a, dummy_b
      real*8  ak_fast(ikmax), pk_fast(ikmax,100), dpk_fast(ikmax,100)
c
      integer iverbose,ifast,icalc_xi,ifid
c
      common  /pk_data/ ak, pk, ik_max
      common  /dpk_fast/ ak_fast, pk_fast, dpk_fast, ik_max_fast
      common  /commandline2/ iverbose,ifast,inum_z,icalc_xi,ifid
c     -----------------------------------------------------
c
      if (iverbose.ge.1) 
     &     write(6,*) '... direct 1-loop computation running...'
      if (iverbose.ge.2) then
         write(6,*) ' > direct 1-loop results'
         write(6,*) ' > iz, k, pk_1loop'
      endif
c      
      ff=frack
      if (icalc_xi.eq.1) ff=1
c
      do ik=1, ik_max, ff
c
         k = ak(ik)
c
c     ///// RegPTmain calculations (at 1-loop order) ////
c
         call one_loop_Gamma1(1, k, G1a_1loop)
         G1b_1loop = G1a_1loop
c
         call calc_pkcorr_G2_tree_tree(1, 1, k, pkcorr_G2_tree_tree)
c
         ak_fast(ik) = k
c
         call calc_running_sigmav2(k, sigmav2_running)
c
         call find_pklin(2, k, pk_lin)
c
         do iz=1, inum_z
            exp_factor = 0.5d0 * (k*D_growth(iz))**2 * sigmav2_running
            if(exp_factor.ge.80.d0) goto 5
c
            G1a_reg = 1.d0 + exp_factor + D_growth(iz)**2*G1a_1loop 
            G1b_reg = 1.d0 + exp_factor + D_growth(iz)**2*G1b_1loop 
            G1a_reg = G1a_reg * D_growth(iz) * dexp(-exp_factor)
            G1b_reg = G1b_reg * D_growth(iz) * dexp(-exp_factor)
c     
            pkcorr_G1 = G1a_reg*G1b_reg * pk_lin
            pkcorr_G2 = pkcorr_G2_tree_tree * D_growth(iz)**4 * 
     &           dexp(-2.d0 * exp_factor)
c
            pk_fast(ik,iz) = pkcorr_G1
            dpk_fast(ik,iz) = pkcorr_G2
c
            if (iverbose.ge.2 .and. mod((ik-1)/ff,10).eq.0) then
               write(6,'(4x,I2,1p2e18.10)') 
     &              iz,k,pk_fast(ik, iz)+dpk_fast(ik, iz)
            endif
 5          continue
         enddo
c
      enddo
c
      ik_max_fast = ik_max
c
      end
c
c ********************************************************** c
c
      subroutine calc_pkcorr_G2_tree_tree(a, b, k, pkcorr_G2_tree_tree)
c
c ********************************************************** c
c
      implicit none
      integer a, b, ix, ixmax, imu, imu_max, ikmax, ik_max
      parameter(ikmax=2000)
      parameter(ixmax=200, imu_max=15)
      real*8  ak(ikmax), pk(ikmax)
      real*8  k, xmin, xmax, integ_pkcorr, pkcorr_G2_tree_tree
      real*8  x, xx(ixmax), wx(ixmax), pi
      real*8  mumin, mumax, mmu(imu_max), wmu(imu_max)
      real*8  q, kq, pk_q, pk_kq, Gamma2_stdPT
      common  /pk_data/ ak, pk, ik_max
      pi = 4.d0 * datan(1.d0)
c     -----------------------------------------------------
c
      xmin = ak(1) / k
      xmax = ak(ik_max) / k
c
      call gauleg(dlog(xmin),dlog(xmax),xx,wx,ixmax)
c
      pkcorr_G2_tree_tree = 0.d0
c
      do ix=1, ixmax
         x = dexp(xx(ix))
c
         mumin = max(-1.d0, (1.d0+x**2-xmax**2)/2.d0/x)
         mumax = min( 1.d0, (1.d0+x**2-xmin**2)/2.d0/x)
         if(x.ge.0.5d0) mumax= 0.5d0/x
c
         call gauleg(mumin, mumax, mmu, wmu, imu_max)
c
         integ_pkcorr = 0.d0
         do imu=1, imu_max
            q = k*x 
            kq = k*dsqrt(1.d0+x**2-2.d0*mmu(imu)*x)
c
            call find_pk(q, pk_q)
            call find_pk(kq, pk_kq)
c
            integ_pkcorr = integ_pkcorr + wmu(imu) * 2.d0 *pk_q*pk_kq
     &           * Gamma2_stdPT(a,0,kq,q,k)*Gamma2_stdPT(b,0,kq,q,k)
         enddo
c
         pkcorr_G2_tree_tree = pkcorr_G2_tree_tree + 
     &        wx(ix) * integ_pkcorr * x**3
c
      enddo
c
      pkcorr_G2_tree_tree = 2.d0 * pkcorr_G2_tree_tree * 
     *     k**3 / (2.d0*pi)**2
c
         end
