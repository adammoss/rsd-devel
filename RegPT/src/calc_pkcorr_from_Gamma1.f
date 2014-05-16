c ******************************************************* c
c
      subroutine calc_pkcorr_from_Gamma1(a, b, k, G1a_1loop, G1a_2loop,
     &     G1b_1loop, G1b_2loop)
c
c ******************************************************* c
c
c     Gamma_a^(1)(k) * Gamma_b^(1)(k) * P(k)
c
c     without exp. cutoff
c
c     labels a,  b = 1 : density
c                  = 2 : velocity divergence
c
      implicit none
c
      integer  a, b
      real*8   k, exp_factor
      real*8   G1a_1loop, G1a_2loop
      real*8   G1b_1loop, G1b_2loop
      real*8   G1a_oneloop, G1a_twoloop, G1b_oneloop, G1b_twoloop
      real*8   sigmav
      common /sigma_v/ sigmav 
c     --------------------------------------------------------
c
      exp_factor = 0.5d0 * (k*sigmav)**2
c
      if(k.lt.1.d-2 .or. exp_factor.gt.50.d0) then
         G1a_1loop = 0.d0
         G1a_2loop = 0.d0
         G1b_1loop = 0.d0
         G1b_2loop = 0.d0
      else
         if(a.eq.b) then
            call one_loop_Gamma1(a, k, G1a_1loop)
            call two_loop_Gamma1(a, k, G1a_2loop)
            G1a_2loop = G1a_2loop + G1a_1loop**2/2.d0 
c
            G1b_1loop = G1a_1loop
            G1b_2loop = G1a_2loop
         else
            call one_loop_Gamma1(a, k, G1a_1loop)
            call two_loop_Gamma1(a, k, G1a_2loop)
            call one_loop_Gamma1(b, k, G1b_1loop)
            call two_loop_Gamma1(b, k, G1b_2loop)
            G1a_2loop = G1a_2loop + G1a_1loop**2/2.d0 
            G1b_2loop = G1b_2loop + G1b_1loop**2/2.d0 
         endif
      endif
c
      end
c
