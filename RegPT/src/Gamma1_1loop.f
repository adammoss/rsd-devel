c ******************************************************* c
c
      subroutine one_loop_Gamma1(a, k, ss)
c
c ******************************************************* c
c
c     computing Gamma1 one-loop
c
c     label a = 1:   density
c           a = 2:   velociy divergence
c     
      implicit none
c
      integer  ik, ikmax, ik_max
      integer  ix, ixmax, a
      parameter(ikmax=2000)
      parameter(ixmax=600)
      real*8   ak(ikmax), pk(ikmax)
      real*8   ak2(ikmax), ff(ikmax), gg(ikmax)
      real*8   xx(ixmax), wx(ixmax)
      real*8   kernel_ff, kernel_gg
      real*8   pi, kmin, kmax, xmin, xmax
      real*8   pklin, x, k, ss
      character infile*50
      common /pk_data/ ak, pk, ik_max
      pi = 4.d0 * datan(1.d0)
c     -------------------------------------------------
c
      kmin = ak(1)
      kmax = ak(ik_max)
c
c
      xmin= kmin / k
      xmax= kmax / k
c
      ss = 0.0
c
      call gauleg(dlog(xmin), dlog(xmax), xx, wx, ixmax)
c
      do ix=1, ixmax
c
         x = dexp(xx(ix))
         call find_pk(k*x, pklin)

         if(a.eq.1) ss = ss + kernel_ff(x) *pklin * wx(ix)
         if(a.eq.2) ss = ss + kernel_gg(x) *pklin * wx(ix)
c     
      enddo
c
         ss =  ss *  k**3 / (2.d0*pi*pi) 
c
      end
c
c ******************************************************* c
c
      function kernel_ff(x)
c
c ******************************************************* c
c
      implicit none
      real*8 kernel_ff, x
c
      if((x.le.10.d0).and. abs(x-1.d0).ge.1.d-2) then
         kernel_ff = 6.d0/x/x - 79.d0 + 50.d0*x*x - 21.d0*x*x*x*x 
     &        + 0.75d0*(1.d0/x-x)**3 * (2.d0 + 7.d0*x*x) 
     &        * dlog(abs((1.d0-x)/(1.d0+x))**2)
         kernel_ff = kernel_ff / 504.d0
      elseif(abs(x-1.d0).lt.1.d-2) then 
         kernel_ff = - 11.d0/126.d0 + (x-1.d0)/126.d0 
     &        - 29.d0/252.d0 * (x-1.d0)**2
      elseif(x.gt.10.d0) then
         kernel_ff = - 61.d0/630.d0 + 2.d0/105.d0/x/x 
     &        - 10.d0/1323.d0 /x/x/x/x
      endif
c
      kernel_ff = kernel_ff * x
c
      end
c
c ******************************************************* c
c
      function kernel_gg(x)
c
c ******************************************************* c
c
      implicit none
      real*8 kernel_gg, x
c
      if((x.le.10.d0).and. abs(x-1.d0).ge.1.d-2) then
         kernel_gg = 6.d0/x/x - 41.d0 + 2.d0*x*x - 3.d0*x*x*x*x 
     &        + 0.75d0*(1.d0/x-x)**3 * (2.d0 + x*x) 
     &        * dlog(abs((1.d0-x)/(1.d0+x))**2)
         kernel_gg = kernel_gg / 168.d0
      elseif(abs(x-1.d0).lt.1.d-2) then 
         kernel_gg = - 3.d0/14.d0 - 5.d0/42.d0*(x-1.d0) 
     &        +  (x-1.d0)**2/84.d0
      elseif(x.gt.10.d0) then
         kernel_gg = - 3.d0/10.d0 + 26.d0/245.d0/x/x 
     &        - 38.d0/2205.d0 /x/x/x/x
      endif
c
      kernel_gg = kernel_gg * x
c
      end
