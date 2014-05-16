c ******************************************************* c
c
      function sigmaab(n, a, b)
c
c ******************************************************* c
c
      implicit none
      integer n, a, b
      real*8  sigmaab
c
      if(a.eq.1 .and. b.eq.1) then
         sigmaab = 2.d0*dble(n) + 1.d0
      elseif(a.eq.1 .and. b.eq.2) then
         sigmaab = 2.d0
      elseif(a.eq.2 .and. b.eq.1) then
         sigmaab = 3.d0
      elseif(a.eq.2 .and. b.eq.2) then
         sigmaab = 2.d0*dble(n)
      else
         sigmaab = 0.d0
      endif
c
      sigmaab = sigmaab / (2.d0*dble(n) + 3.d0) / (dble(n) - 1.d0)
c
      end
c
c ******************************************************* c
c
      function gam_matrix(a, b, c, p, q)
c
c ******************************************************* c
c
      implicit none
      integer a, b, c
      real*8  gam_matrix, p(3), q(3), pp, qq, pq
c
      pp =  p(1)**2 + p(2)**2 + p(3)**2 
      qq =  q(1)**2 + q(2)**2 + q(3)**2 
      pq =  p(1)*q(1) + p(2)*q(2) + p(3)*q(3)
c
      if(a.eq.1 .and. b.eq.1 .and. c.eq.2) then
         gam_matrix = 1.d0 + pq / qq 
      elseif(a.eq.1 .and. b.eq.2 .and. c.eq.1) then
         gam_matrix = 1.d0 + pq / pp 
      elseif(a.eq.2 .and. b.eq.2 .and. c.eq.2) then
         gam_matrix = pq * ( pp+qq+2.d0*pq ) /( pp*qq )
      else
         gam_matrix = 0.d0
      endif
c
      gam_matrix = gam_matrix / 2.d0
c
      end
c
c ******************************************************* c
c
      function F2_sym(a, p, q)
c
c ******************************************************* c
c
      implicit none
      integer a
      real*8 F2_sym, p(3), q(3), pp, qq, mu
      real*8 sigmaab, gam_matrix
c     -------------------------------------------
c
      pp = dsqrt( p(1)**2 + p(2)**2 + p(3)**2 )
      qq = dsqrt( q(1)**2 + q(2)**2 + q(3)**2 )
      mu = ( p(1)*q(1) + p(2)*q(2) + p(3)*q(3) )/( pp*qq )
c
      if(a.eq.1) F2_sym = 
     &     5.d0/7.d0 + mu/2.d0*(pp/qq + qq/pp) + 2.d0/7.d0*mu**2 
      if(a.eq.2) F2_sym = 
     &     3.d0/7.d0 + mu/2.d0*(pp/qq + qq/pp) + 4.d0/7.d0*mu**2 
c
cc      F2_sym = sigmaab(2, a, 1) * ( gam_matrix(1,1,2,p,q) + gam_matrix(1,2,1,p,q))
cc     &     + sigmaab(2, a, 2) * gam_matrix(2,2,2,p,q) 
c
      end
c
c ******************************************************* c
c
      function F3(a, p, q, r)
c
c ******************************************************* c
c
      implicit none
      integer a, i
      real*8 F3, F2_sym, p(3), q(3), r(3), qr(3), sigmaab, gam_matrix
c
      qr(1) = q(1) + r(1)
      qr(2) = q(2) + r(2)
      qr(3) = q(3) + r(3)
c
      if( (dabs(qr(1)).lt.1.d-30 .and. dabs(qr(2)).lt.1.d-30 .and. 
     &     dabs(qr(3)).lt.1.d-30)  .or. 
     &     (dabs(p(1)+qr(1)).lt.1.d-30 .and. dabs(p(2)+qr(2)).lt.1.d-30 
     &     .and. dabs(p(3)+qr(3)).lt.1.d-30) 
     &     ) then
         F3 = 0.d0
      else
         F3 = ( sigmaab(3,a,1) * gam_matrix(1,1,2,p,qr) 
     &        + sigmaab(3,a,2) * gam_matrix(2,2,2,p,qr) ) 
     &        * F2_sym(2,q,r) 
     &        + sigmaab(3,a,1) * gam_matrix(1,2,1,p,qr) 
     &        * F2_sym(1,q,r)
         F3 = 2.d0 * F3
      endif
c
ccc      write(6,*) 'F3=',F3, (p(i),i=1,3),(q(i),i=1,3),(r(i),i=1,3)
      end
c
c ******************************************************* c
c
      function F3_sym(a, p, q, r)
c
c ******************************************************* c
c
c     fully symmetric F3 :  p <--> q <--> r
c
      implicit none
      integer a
      real*8 F3_sym, F3, p(3), q(3), r(3)
c
      F3_sym = ( F3(a, p, q, r) + F3(a, r, p, q) + F3(a, q, r, p) )/3.d0
c
      end
c
c ******************************************************* c
c
      function F4a(a, p, q, r, s)
c
c ******************************************************* c
c
c     symmetric sub-part of F4:   q <--> r <--> s
c
      implicit none
      integer a, i
      real*8 F4a, F3_sym, p(3), q(3), r(3), s(3), qrs(3)
      real*8 sigmaab, gam_matrix
c
      qrs(1) = q(1) + r(1) + s(1)
      qrs(2) = q(2) + r(2) + s(2)
      qrs(3) = q(3) + r(3) + s(3)
c
      if( 
     &     (dabs(qrs(1)).lt.1.d-30 .and. dabs(qrs(2)).lt.1.d-30 .and. 
     &     dabs(qrs(3)).lt.1.d-30)  .or. 
     &     (dabs(p(1)+qrs(1)).lt.1.d-30 .and. dabs(p(2)+qrs(2)).lt.
     &     1.d-30 .and. dabs(p(3)+qrs(3)).lt.1.d-30) 
     &     ) then
         F4a = 0.d0
      else
         F4a = ( sigmaab(4,a,1) * gam_matrix(1,1,2,p,qrs) 
     &        + sigmaab(4,a,2) * gam_matrix(2,2,2,p,qrs) ) 
     &        * F3_sym(2,q,r,s) 
     &        + sigmaab(4,a,1) * gam_matrix(1,2,1,p,qrs) 
     &        * F3_sym(1,q,r,s)
         F4a = 2.d0 * F4a
      endif
c
      end
c
c ******************************************************* c
c
      function F4b(a, p, q, r, s)
c
c ******************************************************* c
c
c     symmetric sub-part of F4:   p <--> q,  r <--> s
c
      implicit none
      integer a, i
      real*8 F4b, F2_sym, p(3), q(3), r(3), s(3), pq(3), rs(3)
      real*8 sigmaab, gam_matrix
c
      pq(1) = p(1) + q(1)
      pq(2) = p(2) + q(2)
      pq(3) = p(3) + q(3)
      rs(1) = r(1) + s(1)
      rs(2) = r(2) + s(2)
      rs(3) = r(3) + s(3)
c
      if(
     &     (dabs(pq(1)).lt.1.d-30 .and. dabs(pq(2)).lt.1.d-30 .and. 
     &     dabs(pq(3)).lt.1.d-30) .or. 
     &     (dabs(rs(1)).lt.1.d-30 .and. dabs(rs(2)).lt.1.d-30 .and.
     &     dabs(rs(3)).lt.1.d-30) .or. 
     &     (dabs(pq(1)+rs(1)).lt.1.d-30 .and. dabs(pq(2)+rs(2)).lt.
     &     1.d-30 .and. dabs(pq(3)+rs(3)).lt.1.d-30) 
     &     ) then
         F4b = 0.d0
      else
         F4b =  sigmaab(4,a,1) * gam_matrix(1,1,2,pq,rs) * F2_sym(1,p,q) 
     &        * F2_sym(2,r,s)
     &        + sigmaab(4,a,1) * gam_matrix(1,2,1,pq,rs) * F2_sym(2,p,q) 
     &        * F2_sym(1,r,s)
     &        + sigmaab(4,a,2) * gam_matrix(2,2,2,pq,rs) * F2_sym(2,p,q) 
     &        * F2_sym(2,r,s)
      endif
c
      end
c
c ******************************************************* c
c
      function F4_sym(a, p, q, r, s)
c
c ******************************************************* c
c
c     fully symmetric F4:  p <--> q <--> r <--> s
c
      implicit none
      integer a, i
      real*8 F4a, F4b, F4_sym, p(3), q(3), r(3), s(3)
c
      F4_sym = ( F4a(a, p, q, r, s) + F4a(a, s, p, q, r) + 
     &     F4a(a, r, s, p, q) + F4a(a, q, r, s, p) )/4.d0 
c
      F4_sym = F4_sym + ( F4b(a, p, q, r, s) + F4b(a, p, r, q, s) 
     &     + F4b(a, p, s, q, r) + F4b(a, r, s, p, q) 
     &     + F4b(a, q, s, p, r) + F4b(a, q, r, p, s) )/6.d0 
c
      end
c
c ******************************************************* c
c
      function F5a(a, p, q, r, s, t)
c
c ******************************************************* c
c
c     symmetric sub-part of F5:   q <--> r <--> s <--> t
c
      implicit none
      integer a, i
      real*8 F5a, F4_sym, sigmaab, gam_matrix
      real*8 p(3), q(3), r(3), s(3), t(3), qrst(3) 
c
      qrst(1) = q(1) + r(1) + s(1) + t(1)
      qrst(2) = q(2) + r(2) + s(2) + t(2)
      qrst(3) = q(3) + r(3) + s(3) + t(3)
c
      if( 
     &     ( dabs(qrst(1)).lt.1.d-30 .and. dabs(qrst(2)).lt.1.d-30 
     &     .and. dabs(qrst(3)).lt.1.d-30)  .or. 
     &     ( dabs(p(1)+qrst(1)).lt.1.d-30 .and. dabs(p(2)+qrst(2)).lt.
     &     1.d-30 .and. dabs(p(3)+qrst(3)).lt.1.d-30)  
     &     ) then
         F5a = 0.d0
      else
         F5a = ( sigmaab(5,a,1) * gam_matrix(1,1,2,p,qrst) 
     &     + sigmaab(5,a,2) * gam_matrix(2,2,2,p,qrst) ) 
     &        * F4_sym(2,q,r,s,t) 
     &     + sigmaab(5,a,1) * gam_matrix(1,2,1,p,qrst) 
     &        * F4_sym(1,q,r,s,t)
         F5a = 2.d0 * F5a
      endif
c
      end
c
c ******************************************************* c
c
      function F5b(a, p, q, r, s, t)
c
c ******************************************************* c
c
c     symmetric sub-part of F5:   p <--> q,  r <--> s <--> t
c
      implicit none
      integer a, i
      real*8 F2_sym, F3_sym, sigmaab, gam_matrix, F5b
      real*8 p(3), q(3), r(3), s(3), t(3), pq(3), rst(3)
c
      pq(1) = p(1) + q(1)
      pq(2) = p(2) + q(2)
      pq(3) = p(3) + q(3)
      rst(1) = r(1) + s(1) + t(1)
      rst(2) = r(2) + s(2) + t(2)
      rst(3) = r(3) + s(3) + t(3)
c
      if( 
     &     ( dabs(pq(1)).lt.1.d-30 .and. dabs(pq(2)).lt.1.d-30 .and. 
     &     dabs(pq(3)).lt.1.d-30 )  .or. 
     &     ( dabs(rst(1)).lt.1.d-30 .and. dabs(rst(2)).lt.1.d-30 .and. 
     &     dabs(rst(3)).lt.1.d-30 )  .or. 
     &     ( dabs(pq(1)+rst(1)).lt.1.d-30 .and. dabs(pq(2)+rst(2)).lt.
     &     1.d-30 .and. dabs(pq(3)+rst(3)).lt.1.d-30 ) 
     &     ) then
         F5b = 0.d0
      else
         F5b =  sigmaab(5,a,1) * gam_matrix(1,1,2,pq,rst) * 
     &        F2_sym(1,p,q) * F3_sym(2,r,s,t)
     &        + sigmaab(5,a,1) * gam_matrix(1,2,1,pq,rst) * 
     &        F2_sym(2,p,q) * F3_sym(1,r,s,t)
     &        + sigmaab(5,a,2) * gam_matrix(2,2,2,pq,rst) * 
     &        F2_sym(2,p,q) * F3_sym(2,r,s,t)
         F5b = 2.d0 * F5b
      endif
c
      end
c
c ******************************************************* c
c
      function F5_sym(a, p, q, r, s, t)
c
c ******************************************************* c
c
c     fully symmetric F5:  p <--> q <--> r <--> s <--> t
c
      implicit none
      integer a, i
      real*8 F5a, F5b, F5_sym, p(3), q(3), r(3), s(3), t(3)
c
      F5_sym = ( F5a(a, p, q, r, s, t) + F5a(a, t, p, q, r, s) + 
     &     F5a(a, s, t, p, q, r) + F5a(a, r, s, t, p, q) + 
     &     F5a(a, q, r, s, t, p) ) / 5.d0
c
      F5_sym = F5_sym + ( F5b(a, p, q, r, s, t) + F5b(a, p, r, q, s, t) 
     &     + F5b(a, p, s, q, r, t) + F5b(a, p, t, q, r, s) 
     &     + F5b(a, q, r, p, s, t) + F5b(a, q, s, p, r, t) 
     &     + F5b(a, q, t, p, r, s) + F5b(a, r, s, p, q, t) 
     &     + F5b(a, r, t, p, q, s) + F5b(a, s, t, p, q, r) ) / 10.d0
c
ccc      write(6,*) 'in F5_sym', F5_sym
c
      end
c
