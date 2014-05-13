module am_routines
  implicit none

  INTERFACE bessi0
    MODULE PROCEDURE bessi0_s,bessi0_v
  END INTERFACE
  INTERFACE bessi1
    MODULE PROCEDURE bessi1_s,bessi1_v
  END INTERFACE
  INTERFACE bessk0
    MODULE PROCEDURE bessk0_s,bessk0_v
  END INTERFACE
  INTERFACE bessk1
    MODULE PROCEDURE bessk1_s,bessk1_v
  END INTERFACE
 INTERFACE bessk
    MODULE PROCEDURE bessk_s,bessk_v
  END INTERFACE

  contains

    FUNCTION select(k,arr)
      USE nrtype; USE nrutil, ONLY : assert,swap
      IMPLICIT NONE
      INTEGER(I4B), INTENT(IN) :: k
      REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
      REAL(SP) :: select
      INTEGER(I4B) :: i,r,j,l,n
      REAL(SP) :: a
      n=size(arr)
      call assert(k >= 1, k <= n, 'select args')
      l=1
      r=n
      do
         if (r-l <= 1) then
            if (r-l == 1) call swap(arr(l),arr(r),arr(l)>arr(r))
            select=arr(k)
            RETURN
         else
            i=(l+r)/2
            call swap(arr(i),arr(l+1))
            call swap(arr(l),arr(r),arr(l)>arr(r))
            call swap(arr(l+1),arr(r),arr(l+1)>arr(r))
            call swap(arr(l),arr(l+1),arr(l)>arr(l+1))
            i=l+1
            j=r
            a=arr(l+1)
            do
               do
                  i=i+1
                  if (arr(i) >= a) exit
               end do
               do
                  j=j-1
                  if (arr(j) <= a) exit
               end do
               if (j < i) exit
               call swap(arr(i),arr(j))
            end do
            arr(l+1)=arr(j)
            arr(j)=a
            if (j >= k) r=j-1
            if (j <= k) l=i
         end if
      end do
    END FUNCTION select

    SUBROUTINE medfit(x,y,a,b,abdev)
      USE nrtype; USE nrutil, ONLY : assert_eq
      IMPLICIT NONE
      REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
      REAL(SP), INTENT(OUT) :: a,b,abdev
      INTEGER(I4B) :: ndata
      REAL(SP) :: aa
      call medfit_private
    CONTAINS
      !BL
      SUBROUTINE medfit_private
	IMPLICIT NONE
 REAL(SP) :: b1,b2,bb,chisq,del,f,f1,f2,sigb,sx,sxx,sxy,sy
 REAL(SP), DIMENSION(size(x)) :: tmp
 ndata=assert_eq(size(x),size(y),'medfit')
 sx=sum(x)
 sy=sum(y)
 sxy=dot_product(x,y)
 sxx=dot_product(x,x)
 del=ndata*sxx-sx**2
 aa=(sxx*sy-sx*sxy)/del
 bb=(ndata*sxy-sx*sy)/del
 tmp(:)=y(:)-(aa+bb*x(:))
 chisq=dot_product(tmp,tmp)
 sigb=sqrt(chisq/del)
 b1=bb
 f1=rofunc(b1)
 b2=bb+sign(3.0_sp*sigb,f1)
 f2=rofunc(b2)
 if (b2 == b1) then
    a=aa
    b=bb
    RETURN
 endif
 do
    if (f1*f2 <= 0.0) exit
    bb=b2+1.6_sp*(b2-b1)
    b1=b2
    f1=f2
    b2=bb
    f2=rofunc(b2)
 end do
 sigb=0.01_sp*sigb
 do
    if (abs(b2-b1) <= sigb) exit
    bb=b1+0.5_sp*(b2-b1)
    if (bb == b1 .or. bb == b2) exit
    f=rofunc(bb)
    if (f*f1 >= 0.0) then
       f1=f
       b1=bb
    else
       f2=f
       b2=bb
    end if
 end do
 a=aa
 b=bb
 abdev=abdev/ndata
END SUBROUTINE medfit_private
!BL
FUNCTION rofunc(b)
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: b
  REAL(SP) :: rofunc
  REAL(SP), PARAMETER :: EPS=epsilon(b)
  INTEGER(I4B) :: j
  REAL(SP), DIMENSION(size(x)) :: arr,d
  arr(:)=y(:)-b*x(:)
  if (mod(ndata,2) == 0) then
     j=ndata/2
     aa=0.5_sp*(select(j,arr)+select(j+1,arr))
  else
     aa=select((ndata+1)/2,arr)
  end if
  d(:)=y(:)-(b*x(:)+aa)
  abdev=sum(abs(d))
  where (y(:) /= 0.0) d(:)=d(:)/abs(y(:))
  rofunc=sum(x(:)*sign(1.0_sp,d(:)), mask=(abs(d(:)) > EPS) )
END FUNCTION rofunc
END SUBROUTINE medfit
 
  FUNCTION bessi0_s(x)
	USE nrtype; USE nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: x
	REAL(DP) :: bessi0_s
	REAL(DP) :: ax
	REAL(DP), DIMENSION(7) :: p = (/1.0_dp,3.5156229_dp,&
		3.0899424_dp,1.2067492_dp,0.2659732_dp,0.360768e-1_dp,&
		0.45813e-2_dp/)
	REAL(DP), DIMENSION(9) :: q = (/0.39894228_dp,0.1328592e-1_dp,&
		0.225319e-2_dp,-0.157565e-2_dp,0.916281e-2_dp,&
		-0.2057706e-1_dp,0.2635537e-1_dp,-0.1647633e-1_dp,&
		0.392377e-2_dp/)
	ax=abs(x)
	if (ax < 3.75) then
		bessi0_s=poly(real((x/3.75_sp)**2,dp),p)
	else
		bessi0_s=(exp(ax)/sqrt(ax))*poly(real(3.75_sp/ax,dp),q)
	end if
	END FUNCTION bessi0_s


	FUNCTION bessi0_v(x)
	USE nrtype; USE nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: x
	REAL(DP), DIMENSION(size(x)) :: bessi0_v
	REAL(DP), DIMENSION(size(x)) :: ax
	REAL(DP), DIMENSION(size(x)) :: y
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	REAL(DP), DIMENSION(7) :: p = (/1.0_dp,3.5156229_dp,&
		3.0899424_dp,1.2067492_dp,0.2659732_dp,0.360768e-1_dp,&
		0.45813e-2_dp/)
	REAL(DP), DIMENSION(9) :: q = (/0.39894228_dp,0.1328592e-1_dp,&
		0.225319e-2_dp,-0.157565e-2_dp,0.916281e-2_dp,&
		-0.2057706e-1_dp,0.2635537e-1_dp,-0.1647633e-1_dp,&
		0.392377e-2_dp/)
	ax=abs(x)
	mask = (ax < 3.75)
	where (mask)
		bessi0_v=poly(real((x/3.75_sp)**2,dp),p,mask)
	elsewhere
		y=3.75_sp/ax
		bessi0_v=(exp(ax)/sqrt(ax))*poly(real(y,dp),q,.not. mask)
	end where
	END FUNCTION bessi0_v
  
     FUNCTION bessi1_s(x)
	USE nrtype; USE nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: x
	REAL(DP) :: bessi1_s
	REAL(DP) :: ax
	REAL(DP), DIMENSION(7) :: p = (/0.5_dp,0.87890594_dp,&
		0.51498869_dp,0.15084934_dp,0.2658733e-1_dp,&
		0.301532e-2_dp,0.32411e-3_dp/)
	REAL(DP), DIMENSION(9) :: q = (/0.39894228_dp,-0.3988024e-1_dp,&
		-0.362018e-2_dp,0.163801e-2_dp,-0.1031555e-1_dp,&
		0.2282967e-1_dp,-0.2895312e-1_dp,0.1787654e-1_dp,&
		-0.420059e-2_dp/)
	ax=abs(x)
	if (ax < 3.75) then
		bessi1_s=ax*poly(real((x/3.75_sp)**2,dp),p)
	else
		bessi1_s=(exp(ax)/sqrt(ax))*poly(real(3.75_sp/ax,dp),q)
	end if
	if (x < 0.0) bessi1_s=-bessi1_s
	END FUNCTION bessi1_s


	FUNCTION bessi1_v(x)
	USE nrtype; USE nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: x
	REAL(DP), DIMENSION(size(x)) :: bessi1_v
	REAL(DP), DIMENSION(size(x)) :: ax
	REAL(DP), DIMENSION(size(x)) :: y
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	REAL(DP), DIMENSION(7) :: p = (/0.5_dp,0.87890594_dp,&
		0.51498869_dp,0.15084934_dp,0.2658733e-1_dp,&
		0.301532e-2_dp,0.32411e-3_dp/)
	REAL(DP), DIMENSION(9) :: q = (/0.39894228_dp,-0.3988024e-1_dp,&
		-0.362018e-2_dp,0.163801e-2_dp,-0.1031555e-1_dp,&
		0.2282967e-1_dp,-0.2895312e-1_dp,0.1787654e-1_dp,&
		-0.420059e-2_dp/)
	ax=abs(x)
	mask = (ax < 3.75)
	where (mask)
		bessi1_v=ax*poly(real((x/3.75_sp)**2,dp),p,mask)
	elsewhere
		y=3.75_sp/ax
		bessi1_v=(exp(ax)/sqrt(ax))*poly(real(y,dp),q,.not. mask)
	end where
	where (x < 0.0) bessi1_v=-bessi1_v
	END FUNCTION bessi1_v  
  
    FUNCTION bessk0_s(x)
	USE nrtype; USE nrutil, ONLY : assert,poly
	!USE nr, ONLY : bessi0
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: x
	REAL(DP) :: bessk0_s
	REAL(DP) :: y
	REAL(DP), DIMENSION(7) :: p = (/-0.57721566_dp,0.42278420_dp,&
		0.23069756_dp,0.3488590e-1_dp,0.262698e-2_dp,0.10750e-3_dp,&
		0.74e-5_dp/)
	REAL(DP), DIMENSION(7) :: q = (/1.25331414_dp,-0.7832358e-1_dp,&
		0.2189568e-1_dp,-0.1062446e-1_dp,0.587872e-2_dp,&
		-0.251540e-2_dp,0.53208e-3_dp/)
	call assert(x > 0.0, 'bessk0_s arg')
	if (x <= 2.0) then
		y=x*x/4.0_sp
		bessk0_s=(-log(x/2.0_sp)*bessi0(x))+poly(y,p)
	else
		y=(2.0_sp/x)
		bessk0_s=(exp(-x)/sqrt(x))*poly(y,q)
	end if
	END FUNCTION bessk0_s


	FUNCTION bessk0_v(x)
	USE nrtype; USE nrutil, ONLY : assert,poly
	!USE nr, ONLY : bessi0
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: x
	REAL(DP), DIMENSION(size(x)) :: bessk0_v
	REAL(DP), DIMENSION(size(x)) :: y
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	REAL(DP), DIMENSION(7) :: p = (/-0.57721566_dp,0.42278420_dp,&
		0.23069756_dp,0.3488590e-1_dp,0.262698e-2_dp,0.10750e-3_dp,&
		0.74e-5_dp/)
	REAL(DP), DIMENSION(7) :: q = (/1.25331414_dp,-0.7832358e-1_dp,&
		0.2189568e-1_dp,-0.1062446e-1_dp,0.587872e-2_dp,&
		-0.251540e-2_dp,0.53208e-3_dp/)
	call assert(all(x > 0.0), 'bessk0_v arg')
	mask = (x <= 2.0)
	where (mask)
		y=x*x/4.0_sp
		bessk0_v=(-log(x/2.0_sp)*bessi0(x))+poly(y,p,mask)
	elsewhere
		y=(2.0_sp/x)
		bessk0_v=(exp(-x)/sqrt(x))*poly(y,q,.not. mask)
	end where
	END FUNCTION bessk0_v
  
   FUNCTION bessk1_s(x)
	USE nrtype; USE nrutil, ONLY : assert,poly
	!USE nr, ONLY : bessi1
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: x
	REAL(DP) :: bessk1_s
	REAL(DP) :: y
	REAL(DP), DIMENSION(7) :: p = (/1.0_dp,0.15443144_dp,&
		-0.67278579_dp,-0.18156897_dp,-0.1919402e-1_dp,&
		-0.110404e-2_dp,-0.4686e-4_dp/)
	REAL(DP), DIMENSION(7) :: q = (/1.25331414_dp,0.23498619_dp,&
		-0.3655620e-1_dp,0.1504268e-1_dp,-0.780353e-2_dp,&
		0.325614e-2_dp,-0.68245e-3_dp/)
	call assert(x > 0.0, 'bessk1_s arg')
	if (x <= 2.0) then
		y=x*x/4.0_sp
		bessk1_s=(log(x/2.0_sp)*bessi1(x))+(1.0_sp/x)*poly(y,p)
	else
		y=2.0_sp/x
		bessk1_s=(exp(-x)/sqrt(x))*poly(y,q)
	end if
	END FUNCTION bessk1_s


	FUNCTION bessk1_v(x)
	USE nrtype; USE nrutil, ONLY : assert,poly
	!USE nr, ONLY : bessi1
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: x
	REAL(DP), DIMENSION(size(x)) :: bessk1_v
	REAL(DP), DIMENSION(size(x)) :: y
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	REAL(DP), DIMENSION(7) :: p = (/1.0_dp,0.15443144_dp,&
		-0.67278579_dp,-0.18156897_dp,-0.1919402e-1_dp,&
		-0.110404e-2_dp,-0.4686e-4_dp/)
	REAL(DP), DIMENSION(7) :: q = (/1.25331414_dp,0.23498619_dp,&
		-0.3655620e-1_dp,0.1504268e-1_dp,-0.780353e-2_dp,&
		0.325614e-2_dp,-0.68245e-3_dp/)
	call assert(all(x > 0.0), 'bessk1_v arg')
	mask = (x <= 2.0)
	where (mask)
		y=x*x/4.0_sp
		bessk1_v=(log(x/2.0_sp)*bessi1(x))+(1.0_sp/x)*poly(y,p,mask)
	elsewhere
		y=2.0_sp/x
		bessk1_v=(exp(-x)/sqrt(x))*poly(y,q,.not. mask)
	end where
	END FUNCTION bessk1_v
	
	  FUNCTION bessk_s(n,x)
	USE nrtype; USE nrutil, ONLY : assert
	!USE nr, ONLY : bessk0,bessk1
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(DP), INTENT(IN) :: x
	REAL(DP) :: bessk_s
	INTEGER(I4B) :: j
	REAL(DP) :: bk,bkm,bkp,tox
	call assert(n >= 2, x > 0.0, 'bessk_s args')
	tox=2.0_sp/x
	bkm=bessk0(x)
	bk=bessk1(x)
	do j=1,n-1
		bkp=bkm+j*tox*bk
		bkm=bk
		bk=bkp
	end do
	bessk_s=bk
	END FUNCTION bessk_s


	FUNCTION bessk_v(n,x)
	USE nrtype; USE nrutil, ONLY : assert
	!USE nr, ONLY : bessk0,bessk1
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(DP), DIMENSION(:), INTENT(IN) :: x
	REAL(DP), DIMENSION(size(x)) :: bessk_v
	INTEGER(I4B) :: j
	REAL(DP), DIMENSION(size(x)) :: bk,bkm,bkp,tox
	call assert(n >= 2, all(x > 0.0), 'bessk_v args')
	tox=2.0_sp/x
	bkm=bessk0(x)
	bk=bessk1(x)
	do j=1,n-1
		bkp=bkm+j*tox*bk
		bkm=bk
		bk=bkp
	end do
	bessk_v=bk
	END FUNCTION bessk_v
	
        FUNCTION zbrent(func,x1,x2,tol,yval)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: x1,x2,tol
        REAL(DP), intent(in), optional :: yval ! Modified
	REAL(DP) :: zbrent
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: ITMAX=100
	REAL(DP), PARAMETER :: EPS=epsilon(x1)
	INTEGER(I4B) :: iter
	REAL(DP) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
	a=x1
	b=x2
        if(present(yval)) then
           fa=func(a)-yval
           fb=func(b)-yval
        else
           fa=func(a)
           fb=func(b)
        end if
	if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) &
		call nrerror('root must be bracketed for zbrent')
	c=b
	fc=fb
	do iter=1,ITMAX
		if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
			c=a
			fc=fa
			d=b-a
			e=d
		end if
		if (abs(fc) < abs(fb)) then
			a=b
			b=c
			c=a
			fa=fb
			fb=fc
			fc=fa
		end if
		tol1=2.0_sp*EPS*abs(b)+0.5_sp*tol
		xm=0.5_sp*(c-b)
		if (abs(xm) <= tol1 .or. fb == 0.0) then
			zbrent=b
			RETURN
		end if
		if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
			s=fb/fa
			if (a == c) then
				p=2.0_sp*xm*s
				q=1.0_sp-s
			else
				q=fa/fc
				r=fb/fc
				p=s*(2.0_sp*xm*q*(q-r)-(b-a)*(r-1.0_sp))
				q=(q-1.0_sp)*(r-1.0_sp)*(s-1.0_sp)
			end if
			if (p > 0.0) q=-q
			p=abs(p)
			if (2.0_sp*p  <  min(3.0_sp*xm*q-abs(tol1*q),abs(e*q))) then
				e=d
				d=p/q
			else
				d=xm
				e=d
			end if
		else
			d=xm
			e=d
		end if
		a=b
		fa=fb
		b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
                if(present(yval)) then
                   fb=func(b)-yval
                else
                   fb=func(b)
                end if
	end do
	call nrerror('zbrent: exceeded maximum iterations')
	zbrent=b
	END FUNCTION zbrent

	
	SUBROUTINE tridag(a,b,c,r,u)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: a,b,c,r
	REAL(DP), DIMENSION(:), INTENT(OUT) :: u
	REAL(DP), DIMENSION(size(b)) :: gam
	INTEGER(I4B) :: n,j
	REAL(DP) :: bet
	n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag')
	bet=b(1)
	if (bet == 0.0) call nrerror('tridag: Error at code stage 1')
	u(1)=r(1)/bet
	do j=2,n
		gam(j)=c(j-1)/bet
		bet=b(j)-a(j-1)*gam(j)
		if (bet == 0.0) &
			call nrerror('tridag: Error at code stage 2')
		u(j)=(r(j)-a(j-1)*u(j-1))/bet
	end do
	do j=n-1,1,-1
		u(j)=u(j)-gam(j+1)*u(j+1)
	end do
	END SUBROUTINE tridag

	RECURSIVE SUBROUTINE tridag_par(a,b,c,r,u)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	!USE nr, ONLY : tridag_ser
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: a,b,c,r
	REAL(DP), DIMENSION(:), INTENT(OUT) :: u
	INTEGER(I4B), PARAMETER :: NPAR_TRIDAG=4
	INTEGER(I4B) :: n,n2,nm,nx
	REAL(DP), DIMENSION(size(b)/2) :: y,q,piva
	REAL(DP), DIMENSION(size(b)/2-1) :: x,z
	REAL(DP), DIMENSION(size(a)/2) :: pivc
	n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_par')
	if (n < NPAR_TRIDAG) then
		call tridag(a,b,c,r,u)
	else
		if (maxval(abs(b(1:n))) == 0.0) &
			call nrerror('tridag_par: possible singular matrix')
		n2=size(y)
		nm=size(pivc)
		nx=size(x)
		piva = a(1:n-1:2)/b(1:n-1:2)
		pivc = c(2:n-1:2)/b(3:n:2)
		y(1:nm) = b(2:n-1:2)-piva(1:nm)*c(1:n-2:2)-pivc*a(2:n-1:2)
		q(1:nm) = r(2:n-1:2)-piva(1:nm)*r(1:n-2:2)-pivc*r(3:n:2)
		if (nm < n2) then
			y(n2) = b(n)-piva(n2)*c(n-1)
			q(n2) = r(n)-piva(n2)*r(n-1)
		end if
		x = -piva(2:n2)*a(2:n-2:2)
		z = -pivc(1:nx)*c(3:n-1:2)
		call tridag_par(x,y,z,q,u(2:n:2))
		u(1) = (r(1)-c(1)*u(2))/b(1)
		u(3:n-1:2) = (r(3:n-1:2)-a(2:n-2:2)*u(2:n-2:2) &
			-c(3:n-1:2)*u(4:n:2))/b(3:n-1:2)
		if (nm == n2) u(n)=(r(n)-a(n-1)*u(n-1))/b(n)
	end if
	END SUBROUTINE tridag_par
	
	
	SUBROUTINE spline_nr(x,y,yp1,ypn,y2)
	USE nrtype; USE nrutil, ONLY : assert_eq
	!USE nr, ONLY : tridag
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(DP), INTENT(IN) :: yp1,ypn
	REAL(DP), DIMENSION(:), INTENT(OUT) :: y2
	INTEGER(I4B) :: n
	REAL(DP), DIMENSION(size(x)) :: a,b,c,r
	n=assert_eq(size(x),size(y),size(y2),'spline')
	c(1:n-1)=x(2:n)-x(1:n-1)
	r(1:n-1)=6.0_sp*((y(2:n)-y(1:n-1))/c(1:n-1))
	r(2:n-1)=r(2:n-1)-r(1:n-2)
	a(2:n-1)=c(1:n-2)
	b(2:n-1)=2.0_sp*(c(2:n-1)+a(2:n-1))
	b(1)=1.0
	b(n)=1.0
	if (yp1 > 0.99e30_sp) then
		r(1)=0.0
		c(1)=0.0
	else
		r(1)=(3.0_sp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
		c(1)=0.5
	end if
	if (ypn > 0.99e30_sp) then
		r(n)=0.0
		a(n)=0.0
	else
		r(n)=(-3.0_sp/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
		a(n)=0.5
	end if
	call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
	END SUBROUTINE spline_nr	
  
   	FUNCTION splint_nr(xa,ya,y2a,x)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	!USE nr, ONLY: locate
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
	REAL(DP), INTENT(IN) :: x
	REAL(DP) :: splint_nr
	INTEGER(I4B) :: khi,klo,n
	REAL(DP) :: a,b,h
	n=assert_eq(size(xa),size(ya),size(y2a),'splint')
	klo=max(min(locate(xa,x),n-1),1)
	khi=klo+1
	h=xa(khi)-xa(klo)
	if (h == 0.0) call nrerror('bad xa input in splint')
	a=(xa(khi)-x)/h
	b=(x-xa(klo))/h
	splint_nr=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_sp
	END FUNCTION splint_nr

 	FUNCTION splint_deriv_nr(xa,ya,y2a,x)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	!USE nr, ONLY: locate
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
	REAL(DP), INTENT(IN) :: x
	REAL(DP) :: splint_deriv_nr
	INTEGER(I4B) :: khi,klo,n
	REAL(DP) :: a,b,h,dely
	n=assert_eq(size(xa),size(ya),size(y2a),'splint')
	klo=max(min(locate(xa,x),n-1),1)
	khi=klo+1
	h=xa(khi)-xa(klo)
	if (h == 0.0) call nrerror('bad xa input in splint')
	a=(xa(khi)-x)/h
	b=(x-xa(klo))/h
        dely=ya(khi)-ya(klo)
	splint_deriv_nr=dely/h-(3.d0*a**2-1.0d0)*h*y2a(klo)/6.0d0 + &
           (3.d0*b**2-1.0d0)*h*y2a(khi)/6.0d0 
	END FUNCTION splint_deriv_nr

	FUNCTION locate(xx,x)
	USE nrtype
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: xx
	REAL(DP), INTENT(IN) :: x
	INTEGER(I4B) :: locate
	INTEGER(I4B) :: n,jl,jm,ju
	LOGICAL :: ascnd
	n=size(xx)
	ascnd = (xx(n) >= xx(1))
	jl=0
	ju=n+1
	do
		if (ju-jl <= 1) exit
		jm=(ju+jl)/2
		if (ascnd .eqv. (x >= xx(jm))) then
			jl=jm
		else
			ju=jm
		end if
	end do
	if (x == xx(1)) then
		locate=1
	else if (x == xx(n)) then
		locate=n-1
	else
		locate=jl
	end if
	END FUNCTION locate

        SUBROUTINE trapzd(func,a,b,s,n)
          USE nrtype; USE nrutil, ONLY : arth
          IMPLICIT NONE
          REAL(DP), INTENT(IN) :: a,b
          REAL(DP), INTENT(INOUT) :: s
          INTEGER(I4B), INTENT(IN) :: n
          INTERFACE
             FUNCTION func(x)
               USE nrtype
               REAL(DP), DIMENSION(:), INTENT(IN) :: x
               REAL(DP), DIMENSION(size(x)) :: func
             END FUNCTION func
          END INTERFACE
          REAL(DP) :: del,fsum
          INTEGER(I4B) :: it
          if (n == 1) then
             s=0.5_sp*(b-a)*sum(func( (/ a,b /) ))
          else
             it=2**(n-2)
             del=(b-a)/it
             fsum=sum(func(arth(a+0.5_sp*del,del,it)))
             s=0.5_sp*(s+del*fsum)
          end if
	END SUBROUTINE trapzd
        
        FUNCTION qtrap(func,a,b)
          USE nrtype; USE nrutil, ONLY : nrerror
          IMPLICIT NONE
          REAL(DP), INTENT(IN) :: a,b
          REAL(DP) :: qtrap
          INTERFACE
             FUNCTION func(x)
               USE nrtype
               REAL(DP), DIMENSION(:), INTENT(IN) :: x
               REAL(DP), DIMENSION(size(x)) :: func
             END FUNCTION func
          END INTERFACE
          INTEGER(I4B), PARAMETER :: JMAX=20
          REAL(SP), PARAMETER :: EPS=1.0e-3_sp, UNLIKELY=-1.0e30_sp
          REAL(DP) :: olds
          INTEGER(I4B) :: j
          olds=UNLIKELY
          do j=1,JMAX
             call trapzd(func,a,b,qtrap,j)
             if (j > 5) then
                if (abs(qtrap-olds) < EPS*abs(olds) .or. &
                     (qtrap == 0.0 .and. olds == 0.0)) RETURN
             end if
             olds=qtrap
          end do
          call nrerror('qtrap: too many steps')
	END FUNCTION qtrap
        
        SUBROUTINE trapzd_f(func,a,b,s,n)
          INTEGER n
          double precision a,b,s,func
          EXTERNAL func
          INTEGER it,j
          double precision del,sum,tnm,x
          if (n.eq.1) then
             s=0.5*(b-a)*(func(a)+func(b))
          else
             it=2**(n-2)
             tnm=it
             del=(b-a)/tnm
             x=a+0.5*del
             sum=0.
             do j=1,it
                sum=sum+func(x)
                x=x+del
             end do
             s=0.5*(s+(b-a)*sum/tnm)
          endif
          return
        END SUBROUTINE trapzd_f

        SUBROUTINE QTRAP_F(FUNC,A,B,S)
          real, PARAMETER :: EPS=1.E-6
          integer, parameter :: JMAX=20
          integer j
          double precision olds,func,a,b,s
          external func
          OLDS=-1.E30
          DO J=1,JMAX
             CALL TRAPZD_F(FUNC,A,B,S,J)
             IF (ABS(S-OLDS).LT.EPS*ABS(OLDS)) RETURN
             OLDS=S
          end do
          stop 'Too many steps.'
        END SUBROUTINE QTRAP_F

       function ran1_f77(idum)
! RANDOM NUMBER GENERATOR FROM NUMERICAL RECIPES, 2ND EDITION 
      implicit double precision (a-h,o-z)     
      integer idum,ia,im,iq,ir,ntab,ndiv

      parameter (ia=16807,im=2147483647,am=1.d0/im,iq=127773,ir=2836)
      parameter (ntab=32,ndiv=1+(im-1)/ntab,eps=1.2d-7,rnmx=1.d0-eps)
      integer j,k,iv(ntab),iy
      save iv,iy
      data iv /ntab*0/, iy /0/

      if (idum.le.0.or.iy.eq.0) then
         idum=max(-idum,1)
         do j=ntab+8,1,-1
            k=idum/iq
            idum=ia*(idum-k*iq)-ir*k
            if (idum.lt.0) idum=idum+im
            if (j.le.ntab) iv(j)=idum
         enddo
         iy=iv(1)
      endif

      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum=idum+im
      j=1+iy/ndiv
      iy=iv(j)
      iv(j)=idum
      ran1_f77=min(am*iy,rnmx)
    end function ran1_f77

    SUBROUTINE polint(xa,ya,x,y,dy)
	USE nrtype; USE nrutil, ONLY : assert_eq,iminloc,nrerror
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya
	REAL(DP), INTENT(IN) :: x
	REAL(DP), INTENT(OUT) :: y,dy
	INTEGER(I4B) :: m,n,ns
	REAL(DP), DIMENSION(size(xa)) :: c,d,den,ho
	n=assert_eq(size(xa),size(ya),'polint')
	c=ya
	d=ya
	ho=xa-x
	ns=iminloc(abs(x-xa))
	y=ya(ns)
	ns=ns-1
	do m=1,n-1
		den(1:n-m)=ho(1:n-m)-ho(1+m:n)
		if (any(den(1:n-m) == 0.0)) &
			call nrerror('polint: calculation failure')
		den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
		d(1:n-m)=ho(1+m:n)*den(1:n-m)
		c(1:n-m)=ho(1:n-m)*den(1:n-m)
		if (2*ns < n-m) then
			dy=c(ns+1)
		else
			dy=d(ns)
			ns=ns-1
		end if
		y=y+dy
	end do
	END SUBROUTINE polint

 SUBROUTINE hunt(xx,x,jlo)
	USE nrtype
	IMPLICIT NONE
	INTEGER(I4B), INTENT(INOUT) :: jlo
	REAL(DP), INTENT(IN) :: x
	REAL(DP), DIMENSION(:), INTENT(IN) :: xx
	INTEGER(I4B) :: n,inc,jhi,jm
	LOGICAL :: ascnd
	n=size(xx)
	ascnd = (xx(n) >= xx(1))
	if (jlo <= 0 .or. jlo > n) then
		jlo=0
		jhi=n+1
	else
		inc=1
		if (x >= xx(jlo) .eqv. ascnd) then
			do
				jhi=jlo+inc
				if (jhi > n) then
					jhi=n+1
					exit
				else
					if (x < xx(jhi) .eqv. ascnd) exit
					jlo=jhi
					inc=inc+inc
				end if
			end do
		else
			jhi=jlo
			do
				jlo=jhi-inc
				if (jlo < 1) then
					jlo=0
					exit
				else
					if (x >= xx(jlo) .eqv. ascnd) exit
					jhi=jlo
					inc=inc+inc
				end if
			end do
		end if
	end if
	do
		if (jhi-jlo <= 1) then
			if (x == xx(n)) jlo=n-1
			if (x == xx(1)) jlo=1
			exit
		else
			jm=(jhi+jlo)/2
			if (x >= xx(jm) .eqv. ascnd) then
				jlo=jm
			else
				jhi=jm
			end if
		end if
	end do
	END SUBROUTINE hunt

SUBROUTINE gauleg(x1,x2,x,w)
	USE nrtype; USE nrutil, ONLY : arth,assert_eq,nrerror
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: x1,x2
	REAL(DP), DIMENSION(:), INTENT(OUT) :: x,w
	REAL(DP), PARAMETER :: EPS=3.0e-14_dp
	INTEGER(I4B) :: its,j,m,n
	INTEGER(I4B), PARAMETER :: MAXIT=10
	REAL(DP) :: xl,xm
	REAL(DP), DIMENSION((size(x)+1)/2) :: p1,p2,p3,pp,z,z1
	LOGICAL(LGT), DIMENSION((size(x)+1)/2) :: unfinished
	n=assert_eq(size(x),size(w),'gauleg')
	m=(n+1)/2
	xm=0.5_dp*(x2+x1)
	xl=0.5_dp*(x2-x1)
	z=cos(PI_D*(arth(1,1,m)-0.25_dp)/(n+0.5_dp))
	unfinished=.true.
	do its=1,MAXIT
		where (unfinished)
			p1=1.0
			p2=0.0
		end where
		do j=1,n
			where (unfinished)
				p3=p2
				p2=p1
				p1=((2.0_dp*j-1.0_dp)*z*p2-(j-1.0_dp)*p3)/j
			end where
		end do
		where (unfinished)
			pp=n*(z*p1-p2)/(z*z-1.0_dp)
			z1=z
			z=z1-p1/pp
			unfinished=(abs(z-z1) > EPS)
		end where
		if (.not. any(unfinished)) exit
	end do
	if (its == MAXIT+1) call nrerror('too many iterations in gauleg')
	x(1:m)=xm-xl*z
	x(n:n-m+1:-1)=xm+xl*z
	w(1:m)=2.0_dp*xl/((1.0_dp-z**2)*pp**2)
	w(n:n-m+1:-1)=w(1:m)
	END SUBROUTINE gauleg

subroutine gaulegf(x1, x2, x, w, n)
  implicit none
  integer, intent(in) :: n
  double precision, intent(in) :: x1, x2
  double precision, dimension(n), intent(out) :: x, w
  integer :: i, j, m
  double precision :: p1, p2, p3, pp, xl, xm, z, z1
  double precision, parameter :: eps=3.d-14
      
  m = (n+1)/2
  xm = 0.5d0*(x2+x1)
  xl = 0.5d0*(x2-x1)
  do i=1,m
    z = cos(3.141592654d0*(i-0.25d0)/(n+0.5d0))
    z1 = 0.0
    do while(abs(z-z1) .gt. eps)
      p1 = 1.0d0
      p2 = 0.0d0
      do j=1,n
        p3 = p2
        p2 = p1
        p1 = ((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j
      end do
      pp = n*(z*p1-p2)/(z*z-1.0d0)
      z1 = z
      z = z1 - p1/pp
    end do
    x(i) = xm - xl*z
    x(n+1-i) = xm + xl*z
    w(i) = (2.0d0*xl)/((1.0d0-z*z)*pp*pp)
    w(n+1-i) = w(i)
  end do
end subroutine gaulegf


! ************************************************ 
!
      function gamhalf(n)
!
! ************************************************ 
!
!     Gamma(n+1/2) up to n=6
!
      integer n
      real*8 gamhalf,pi
      pi = 4.d0 *datan(1.d0)

      if(n.eq.0) gamhalf = 1.d0
      if(n.eq.1) gamhalf = 0.5d0
      if(n.eq.2) gamhalf = 0.75d0
      if(n.eq.3) gamhalf = 1.875d0
      if(n.eq.4) gamhalf = 6.5625d0
      if(n.eq.5) gamhalf = 29.53125d0
      if(n.eq.6) gamhalf = 162.421875d0
      if(n.eq.7) gamhalf = 1055.7421875d0
      if(n.eq.8) gamhalf = 7918.06640625d0

      gamhalf = gamhalf * dsqrt(pi)

    end function gamhalf

! ************************************************************
!
      FUNCTION gammp(a,x)
!
! ************************************************************
!
      REAL*8 a,gammp,x
!CU    USES gcf,gser
      REAL*8 gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.)pause 'bad arguments in gammp'
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=1.-gammcf
      endif
      return
    END FUNCTION gammp

!
! ************************************************************
!
      SUBROUTINE gser(gamser,a,x,gln)
!
! ************************************************************
!
      INTEGER ITMAX
      REAL*8 a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.d-7)
!CU    USES gammln
      INTEGER n
      REAL*8 ap,del,sum
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)pause 'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
     end do
      pause 'a too large, ITMAX too small in gser'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
    END SUBROUTINE gser
!
! ************************************************************

      SUBROUTINE gcf(gammcf,a,x,gln)
!
! ************************************************************
!
      INTEGER ITMAX
      REAL*8 a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.d-7,FPMIN=1.d-30)
!CU    USES gammln
      INTEGER i
      REAL*8 an,b,c,d,del,h
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
     end do
      pause 'a too large, ITMAX too small in gcf'
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
    end SUBROUTINE gcf
!
! ************************************************************
!
      FUNCTION gammln(xx)
!
! ************************************************************
!
      REAL*8 gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0, &
     24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
     -.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
     end do
      gammln=tmp+log(stp*ser/x)
      return
    END FUNCTION gammln


end module am_routines





