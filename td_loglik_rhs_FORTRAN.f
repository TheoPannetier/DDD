
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! Example how to apply Fortran code with variable-length parameters
! ---------------------------------------------------------------------
! ------------------------------------------------------------------------

!==========================================================================
! Helper function: 
! fill vec with N elements from parms, starting at position ii
!==========================================================================

      SUBROUTINE fill1d (vec, DIMP, parms, II)
      IMPLICIT NONE
      INTEGER DIMP, II, I
      DOUBLE PRECISION vec(DIMP), parms(*)
      II = II
        DO I = 1, DIMP
          II = II + 1
          vec(I) = parms(II)
        ENDDO
        
      END SUBROUTINE fill1d

!==========================================================================
! module with declarations
!==========================================================================

      MODULE dimmod

      ! length of the vector -  decided in R-code
      INTEGER  :: N
      INTEGER  :: kk
      
      ! 1 parameter vectors with unknown length
      DOUBLE PRECISION, ALLOCATABLE  :: P(:)  
      
      ! Boolean: will become TRUE if the parameters have a value
      LOGICAL :: initialised = .FALSE.

      END MODULE dimmod

!==========================================================================
!==========================================================================
! Initialisation: name of this function as passed by "initfunc" argument
! Sets the fixed parameter vector, and allocates memory
!==========================================================================
!==========================================================================

      SUBROUTINE initmod (steadyparms)
      USE dimmod 

      IMPLICIT NONE
      EXTERNAL steadyparms

      INTEGER, PARAMETER :: nparsmall = 2  ! constant-length parameters
      
      DOUBLE PRECISION parms(nparsmall)
      COMMON /XCBPar/parms                 ! common block 

! Set the fixed parameters obtained from R
      CALL steadyparms(nparsmall, parms)

! first parameter has the length of the vector       
      N = INT(parms(1) + 1e-6)
      kk = INT(parms(2) + 1e-6)

! Allocate variable size arrays (state variables, derivatives and parameters)

      IF (ALLOCATED(P)) DEALLOCATE(P)  
      ALLOCATE(P(3 * (N + 2 + 2 * kk)))

      initialised = .FALSE.
       
      END SUBROUTINE initmod
      
!==========================================================================
!==========================================================================
! Dynamic routine: name of this function as passed by "func" argument
! variable parameter values are passed via yout
!==========================================================================
!==========================================================================
       
      SUBROUTINE runmod3 (neq, t, Conc, dConc, yout, ip)
      USE dimmod
      IMPLICIT NONE
! N is the length of probs and kk is the number of sigmas
!......................... declaration section.............................
      INTEGER           :: neq, ip(*), i, ii
      DOUBLE PRECISION  :: t, Conc(N - kk), dConc(N - kk), yout(*)
      DOUBLE PRECISION  :: V(N + 2)
      DOUBLE PRECISION  :: lavec((N - kk) + 2),muvec((N - kk) + 2)
      DOUBLE PRECISION  :: nn((N - kk) + 2)
      DOUBLE PRECISION  :: FF1, FF2, FF3

! parameters - named here
      DOUBLE PRECISION rn(2)
      COMMON /XCBPar/rn

! local variables
      CHARACTER(len=100) msg

!............................ statements ..................................

      IF (.NOT. Initialised) THEN
        ! check memory allocated to output variables
        IF (ip(1) < 1) CALL rexit("nout not large enough") 

        ! save parameter values in yout
        ii = ip(1)   ! Start of parameter values
        CALL fill1d(P, 3 * ((N - kk) + 2), yout, ii)   ! ii is updated in fill1d
        Initialised = .TRUE.          ! to prevent from initialising more than once
      ENDIF

! dynamics

 !  dx = lavec[(2:(lx+1))+kk-1] * nn[(2:(lx+1))+2*kk-1] * xx[(2:(lx+1))-1] + muvec[(2:(lx+1))+kk+1] 
    !  * nn[(2:(lx+1))+1] * xx[(2:(lx+1))+1] - (lavec[(2:(lx+1))+kk] + muvec[(2:(lx+1))+kk]) * nn[(2:(lx+1))+kk] * xx[2:(lx+1)]
 ! met kk = 0
 !  mutd = rep(mu,lrs)
 !  En = sum((0:(lx - 1)) * x[1:lx] )
 !  dsigdiv = mutd / En

      M = N - kk
      V(1) = 0
      DO I = 2, M + 1 
        V(I) = Conc(I - 1)
      ENDDO
      V(M + 2) = 0
      DO I = 1, M + 2 + 2 * 0
       lavec(I) = P(I)
       muvec(I) = P(I + M + 2)
       nn(I)    = P(I + 2 * (M + 2))
      ENDDO
      En = 0
      DO I = 1, M
       En       = En + (I - 1) * Conc[I]
      ENDDO

      DO I = 2, M + 1 
        FF1 = lavec(I + 0 - 1) * nn(I + 2 * 0 - 1) * V(I - 1)
        FF2 = muvec(I + 0 + 1) * nn(I + 1) * V(I + 1)
        FF3 = (lavec(I + 0) + muvec(I + 0)) * nn(I + 0) * V(I)
        dConc(I - 1) = FF1 + FF2 - FF3
      ENDDO
      
      DO I = 1,kk
        dConc(M + I) = muvec(1)/En
      ENDDO
  
      END SUBROUTINE runmod3
