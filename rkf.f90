	module types
	  integer,parameter::dp=selected_real_kind(15)
	end module types

	module vars
	use types
	implicit none
	Real(DP),parameter::pi=3.14159265359
	Real(DP),parameter::g=32.2, ht=200._dp, k=200._dp, tc=4._dp, fm=0.024, l=2000._dp, d=3.60901, db=1.5, qv0=10.027677, h0=168._dp
	Real(DP),parameter::u0=(4._dp*qv0)/(pi*(d)**2._dp)
	!Real(DP)::qv
	end module vars



	module check_realin
	  contains
	    subroutine realcheck(w, in_var)
	      use types
	      implicit none
	      character(*), intent(in)::w
	      real(dp), intent(out)::in_var
	      integer::ioerror
	      character(len=5):: fmt="(/,a)"
	    do
	      write(*,fmt,advance="no") w
	      read(*,*,iostat=ioerror) in_var
	      if(ioerror==0)exit
	      write(*,fmt) "INVALID INPUT!"
	    end do
	    end subroutine realcheck
	end module check_realin

	! module to change character to integer type
	module check_intgin
	  contains
	      subroutine intcheck(w, in_var)
	      use types
	      implicit none
	      character(*), intent(in)::w
	      integer, intent(out)::in_var
	      integer::ioerror
	      do
		  write(*,*) w
		  read(*,*,iostat=ioerror) in_var
		  if(ioerror==0)exit
		  write(*,*) "INVALID INPUT!"

	      end do
	    end subroutine intcheck
	end module check_intgin
!=============================================================
!=============================================================
	Program Surgetank

	use types
	use vars
	use check_realin
	use check_intgin
!	use openfile


	implicit none
	integer::ntimes, i, exitflag, neqn
	real(dp), allocatable, dimension(:)::y
	real(dp)::xstart, xend, h, hmin, hmax, epsi1, epsi2
	real(dp):: dx_report, x0, xfinish


    ! user input prompts. need error check modules
    character(40), parameter::w0="Number of equations?:",&
	    w1="Value of epsilon 1?:",&
	    w2="Value of epsilon 2?:",&
	    w3="Value for h minimum?:",&
	    w4="Value for h maximum?:",&
	    w5="Initial time value?:",&
	    w6="Final time value?:",&
	    w7="Number of intervals?:",&
	    w8="Value of h?:"
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! series of interfaces
    interface
      subroutine rkf(neqn, h, y, f, xstart, xend, hmin, hmax, epsi1, epsi2, exitflag)
	use types
	real(dp),dimension(:),intent(inout)::y
	real(dp),intent(in)::hmin, hmax, epsi1,epsi2
	real(dp),intent(inout)::h,xstart,xend
	integer,intent(out)::exitflag
	integer,intent(in)::neqn

	interface
	  function f(x,y,neqn)
	    use types
	    use vars
	    integer, intent(in)::neqn
	    real(dp),intent(in)::x
	    real(dp),dimension(neqn),intent(in)::y
    !be careful to dimension the function the same as the number of ODEs
	    real(dp),dimension(neqn)::f
	  end function f
	end interface
      end subroutine rkf

      function dy(x,y,neqn)
      use types
      use vars
	integer, intent(in)::neqn
	real(dp),intent(in)::x
	real(dp),dimension(neqn),intent(in)::y
	real(dp),dimension(neqn)::dy
      end function dy
    end interface
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !Call the user-entered values.
	call intcheck(w0,neqn)
	  allocate(y(neqn))
	call intcheck(w7,ntimes)
	call realcheck(w1,epsi1)
	call realcheck(w2,epsi2)
	call realcheck(w3,hmin)
	call realcheck(w4,hmax)
	call realcheck(w5,x0)
	call realcheck(w6,xfinish)

!     open(unit=19,file="rkfoutput.csv", status="replace")
    dx_report=(xfinish-x0)/ntimes


!call file_open(11,.false.,.false.)
    open(unit=19,file="rkfoutput.csv", status="replace")

	xstart=x0
	DO I=1,ntimes
	  xend=xstart+dx_report
	  ! call RKF subtroutine here, after writing the initial conditions
if(xstart==0)  then
y(1)=u0
y(2)=h0
end if
	call rkf(neqn,h,y,dy,xstart,xend,hmin,hmax,epsi1,epsi2,exitflag)


!write evaluated functions at end of reporting intervals. Won't report initial values
	    write(*,*)xend,y(2),y(1)
	    write(19,*)xend,y(2),y(1)

	    xstart=xend

	end do


    End Program Surgetank

!=================================================================================
!=================================================================================

    ! Subroutine for RKF method
    subroutine rkf(neqn,h,y,f,xstart,xend,hmin,hmax,epsi1,epsi2,exitflag)
      use types
      implicit none
      real(dp),dimension(:),intent(inout)::y
      real(dp),intent(in)::hmin,hmax,epsi1,epsi2,xstart
      real(dp),intent(inout)::h,xend
      integer,intent(out)::exitflag
      integer,intent(in)::neqn
      real(dp),dimension(neqn)::k1,k2,k3,k4,k5,k6,y4,ysave
      real(dp)::x,E,hsave
      logical::short 
      ! short is for when h is bigger than the reporting interval. 
      !allows a way to change the h size to fit in the interval

    ! These local vars have to go in the subroutine. Can't have them in a module
    real(dp),parameter::C2=1._dp/5._dp
    real(dp),parameter::C3=3._dp/10._dp, C3k1=3._dp/40._dp, C3k2=9._dp/40._dp
    real(dp),parameter::C4=3._dp/5._dp, C4k1=3._dp/10._dp, C4k2=-9._dp/10._dp, C4k3=6._dp/5._dp
    real(dp),parameter::C5k1=-11._dp/54._dp, C5k2=5._dp/2._dp, C5k3=-70._dp/27._dp, C5k4=35._dp/27._dp
    real(dp),parameter::C6=7._dp/8._dp, C6k1=1631._dp/55296._dp, C6k2=175._dp/512._dp, C6k3=575._dp/13824._dp
    real(dp),parameter::C6k4=44275._dp/110592._dp, C6k5=253._dp/4096._dp
    real(dp),parameter::Cy4k1=37._dp/378._dp, Cy4k3=250._dp/621._dp, Cy4k4=125._dp/594._dp, Cy4k6=512._dp/1771._dp
    real(dp),parameter::Cy5k1=2825._dp/27648._dp, Cy5k3=18575._dp/48384._dp, Cy5k4=13525._dp/55296._dp, Cy5k5=277._dp/14336._dp
    real(dp),parameter::Cy5k6=1._dp/4._dp

	!subroutine has interface to function to be evaluated from the RKF subroutine
	interface
	  function f(x,y,neqn)
	  use types
	  integer,intent(in)::neqn
	  real(dp),intent(in)::x
	  real(dp),dimension(:),intent(in)::y
	  real(dp),dimension(neqn)::f
	  end function f
	end interface
! Might want to error check here for input values or RKF intermediate values
      exitflag=0
      x=xstart
      if((x+h)>xend)h=xend-x ! adjusts initial h so it does not exceed the width of the interval
      short=.false.
      do
	      ysave=y(1:neqn)         ! store solutions in vector
	
	! compute k1,k2,k3,k4,k5,y,y4
	! six function references
	! These are the coefficients for the k values, k1 - k6. They INCLUDE THE SIGNS if negative
	! last two are for the 4th and 5th order estimates, respectively

	k1=f(x,y(1:neqn),neqn)
	k2=f(x+C2*h,y(1:neqn)+C2*k1*h,neqn)
	k3=f(x+C3*h,y(1:neqn)+C3k1*k1*h+C3k2*k2*h,neqn)
	k4=f(x+C4*h,y(1:neqn)+C4k1*k1*h+C4k2*k2*h+C4k3*k3*h,neqn)
	k5=f(x+h,y(1:neqn)+C5k1*k1*h+C5k2*k2*h+C5k3*k3*h+C5k4*k4*h,neqn)
	k6=f(x+C6*h,y(1:neqn)+C6k1*k1*h+C6k2*k2*h+C6k3*k3*h+C6k4*k4*h+C6k5*k5*h,neqn)
	y4=y(1:neqn)+h*(Cy4k1*k1+Cy4k3*k3+Cy4k4*k4+Cy4k6*k6)
	y(1:neqn)=y(1:neqn)+h*(Cy5k1*k1+Cy5k3*K3+Cy5k4*k4+Cy5k5*k5+Cy5k6*k6)

	! bookkeeping for h
	! compute E - estimate of error for each of the difeqns - this is the worst case truncation error estimate
	E=maxval(abs((y(1:neqn)-y4)/y(1:neqn)))
	! two parts - do it again or move forward
	if (E>epsi2) then
	!reduce step size and try again. or if we're at min step size, then quit
	  if(abs(h-hmin)<1E-6_dp) then 
	! or 1E-6_dp same as 1d-6
	    exitflag=2
	    Return
	  !!!check exitflag in the main program!!!!!
	  end if
	h=h/2._dp
	  if (h<hmin) h=hmin
	!gives it one more chance. don't need to close the if statement since it's a logical if.
	!now back up y 
	y(1:neqn)=ysave ! reset for repeat

	else
	  ! good answer. now we can update x, so it'll catch up to where y is
	  x=x+h  
	  ! do we want to lengthen the step size? based on truncation error criteria:
	    if(E<epsi1)then
	      h=h*2._dp
		if (h>hmax) h=hmax !don't let it get any bigger than the biggest step size that has been specified
	    end if
	    if(.not.short) hsave=h
	  if(abs(x-xend)<1E-6_dp) then
	    exitflag=1
	    h=hsave
	    return
	  end if
	! be sure the next step doesn't go off the end (past xend)
	  if(x+h>xend)then
	    short=.true.
	    h=xend-x
	  end if
	End if
	! maybe put write statements in here to check outputs, stepsize,etc
	end do
	return
	end subroutine RKF


      function dy(x,y,neqn)
	use types
	use vars
	implicit none

	integer,intent(in)::neqn
	real(dp),intent(in)::x
	real(dp),dimension(:),intent(in)::y
	real(dp),dimension(neqn)::dy
	integer::i


! functions go here in the following form, where dy(n) refers to the differential equation number, 
! which will reference another dependent variable that changes with time, also in the form y(n). 
! Be consistent when numeric assignments.
! Initial conditions go above, just before the call to the rkf subroutine, as logicals if statement.
! For this program, y(1)=u, y(2)=h, dy(1)=du/dt, dy(2)=dh/dt. 
! RKF will output solutions to dy(1) and dy(2) (integrated with respect to y and x(time) at each reporting time interval

    if (x>=0.and.tc>=x) then
    dy(2)=(y(1)*((d**2._dp)/(db**2._dp)))-((4._dp/(pi*db**2._dp))*(qv0*(1-(x/tc))**(y(2)/k)))
    dy(1)= ((1._dp/l)*(g*(ht-y(2))))-((1._dp/2._dp)*(fm/d)*y(1)*abs(y(1)))    !du/dt

    else
    dy(2)=(y(1)*d**2._dp)/(db**2._dp)
    dy(1)= ((1._dp/l)*(g*(ht-y(2))))-((1._dp/2._dp)*(fm/d)*y(1)*abs(y(1)))    !du/dt
    end if

    return
    end function

