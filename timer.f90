!***********************************************************************************************************************************
!
!                                                           T I M E R
!
!  Program:      TIMER
!
!  Programmer:   David G. Simpson
!                NASA Goddard Space Flight Center
!                Greenbelt, Maryland  20771
!
!  Date:         September 15, 2006
!
!  Language:     Fortran-90
!
!  Version:      1.00a
!
!  Description:  This module acts as a stopwatch.  It may be used to compute the total time required for a program to run.
!
!  Files:        Source file:
!
!                   timer.f90                   Timer module
!
!  Notes:        Three public subroutines are available:
!
!                   CALL TIMER_START                       starts timer
!                   CALL TIMER_STOP                        stops timer
!                   CALL TIMER_ELAPSED (timestr)           computes elapsed time, where timestr is a CHARACTER(LEN=8) string
!                                                             giving the elapsed time in the form hh:mm:ss.
!
!                The timer can measure elapsed times up to 100 hours.  After 100 hours, the clock will roll over to 0.
!
!  Example:
!                   CHARACTER(LEN=8) :: ELAPSED_TIME
!                   CALL TIMER_START
!                   ! program code goes here
!                   CALL TIMER_STOP
!                   CALL TIMER_ELAPSED(ELAPSED_TIME)
!                   PRINT *, ' TIME ELAPSED = '//ELAPSED_TIME
!
!***********************************************************************************************************************************

      MODULE TIMER

!
!     Global variable declarations.
!

      DOUBLE PRECISION, PRIVATE :: SYS_START_JD, SYS_STOP_JD

      CONTAINS

!-----------------------------------------------------------------------------------------------------------------------------------
!  TIMER_START
!
!  This subroutine starts the timer.
!-----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE TIMER_START

      IMPLICIT NONE

      CHARACTER(LEN=10) :: SYS_START_TIME, SYS_START_DATE, SYS_START_ZONE
      INTEGER, DIMENSION(8) :: SYS_START_DT


      CALL DATE_AND_TIME (SYS_START_DATE, SYS_START_TIME, SYS_START_ZONE, SYS_START_DT)
      SYS_START_JD = JULDAY(SYS_START_DT(2),SYS_START_DT(3),SYS_START_DT(1),SYS_START_DT(5),SYS_START_DT(6),SYS_START_DT(7))

      RETURN

      END SUBROUTINE TIMER_START





!-----------------------------------------------------------------------------------------------------------------------------------
!  TIMER_STOP
!
!  This subroutine stops the timer.
!-----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE TIMER_STOP

      IMPLICIT NONE

      CHARACTER(LEN=10) :: SYS_STOP_TIME, SYS_STOP_DATE, SYS_STOP_ZONE
      INTEGER, DIMENSION(8) :: SYS_STOP_DT

!
!     Stop timer and compute elapsed time.
!     Put this at the end of the program.
!

      CALL DATE_AND_TIME (SYS_STOP_DATE, SYS_STOP_TIME, SYS_STOP_ZONE, SYS_STOP_DT)
      SYS_STOP_JD = JULDAY(SYS_STOP_DT(2),SYS_STOP_DT(3),SYS_STOP_DT(1),SYS_STOP_DT(5),SYS_STOP_DT(6),SYS_STOP_DT(7))

      RETURN

      END SUBROUTINE TIMER_STOP





!-----------------------------------------------------------------------------------------------------------------------------------
!  TIMER_ELAPSED
!
!  This subroutine computes the elapsed time.  The result is returned as an 8-character string of the form hh:mm:ss.
!-----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE TIMER_ELAPSED (ELAPSED_STR)

      IMPLICIT NONE

      CHARACTER(LEN=8), INTENT(OUT) :: ELAPSED_STR

      INTEGER :: COMPUTE_HOURS, COMPUTE_MINUTES, COMPUTE_SECONDS
      DOUBLE PRECISION :: COMPUTE_TIME


!
!     Compute elapsed time and break down into hours, minutes, and seconds.
!

      COMPUTE_TIME = (SYS_STOP_JD - SYS_START_JD)*24.0D0                        ! hours
      COMPUTE_HOURS = INT(COMPUTE_TIME)                                         ! hours
      COMPUTE_TIME = COMPUTE_TIME - COMPUTE_HOURS                               ! hours
      COMPUTE_TIME = COMPUTE_TIME * 60.0D0                                      ! minutes
      COMPUTE_MINUTES = INT(COMPUTE_TIME)                                       ! minutes
      COMPUTE_TIME = COMPUTE_TIME - COMPUTE_MINUTES                             ! minutes
      COMPUTE_TIME = COMPUTE_TIME * 60.0D0                                      ! seconds
      COMPUTE_SECONDS = NINT(COMPUTE_TIME)                                      ! seconds

!
!     Handle rollover problems.
!

      IF (COMPUTE_SECONDS .EQ. 60) THEN
         COMPUTE_SECONDS = 0
         COMPUTE_MINUTES = COMPUTE_MINUTES + 1
      END IF

      IF (COMPUTE_MINUTES .EQ. 60) THEN
         COMPUTE_MINUTES = 0
         COMPUTE_HOURS = COMPUTE_HOURS + 1
      END IF

      IF (COMPUTE_HOURS .GE. 100) THEN
         COMPUTE_HOURS = MOD (COMPUTE_HOURS, 100)
      END IF

!
!     Print elapsed time to a string.
!

      WRITE (UNIT=ELAPSED_STR, FMT='(I2.2,2(1H:,I2.2))') COMPUTE_HOURS, COMPUTE_MINUTES, COMPUTE_SECONDS

      END SUBROUTINE TIMER_ELAPSED





!-----------------------------------------------------------------------------------------------------------------------------------
!  JULDAY - Compute Julian day from Gregorian date.
!-----------------------------------------------------------------------------------------------------------------------------------

      FUNCTION JULDAY (MONTH, DAY, YEAR, HOUR, MINUTE, SECOND) RESULT (JD)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: MONTH, DAY, YEAR, HOUR, MINUTE, SECOND
      DOUBLE PRECISION :: JD

      DOUBLE PRECISION :: DAYFRAC
      INTEGER :: Y, M, A, B


      DAYFRAC = DBLE(DAY) + HOUR/24.0D0 + MINUTE/1440.0D0 + SECOND/86400.0D0

      IF (MONTH .GT. 2) THEN
         Y = YEAR
         M = MONTH
      ELSE
         Y = YEAR - 1
         M = MONTH + 12
      END IF

      A = Y/100
      B = 2 - A + A/4

      JD = INT(365.25D0*(Y+4716)) + INT(30.6001D0*(M+1)) + DAYFRAC + B - 1524.5D0

      RETURN

      END FUNCTION JULDAY

      END MODULE TIMER

