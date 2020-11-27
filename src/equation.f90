!===============================================================================!
MODULE MOD_Equation
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE ExactFunction
  MODULE PROCEDURE ExactFunction
END INTERFACE

INTERFACE SourceTerms
  MODULE PROCEDURE SourceTerms
END INTERFACE

INTERFACE BoundaryConditions
  MODULE PROCEDURE BoundaryConditions
END INTERFACE

INTERFACE TimeStep
  MODULE PROCEDURE TimeStep
END INTERFACE

INTERFACE ConsToPrim
  MODULE PROCEDURE ConsToPrim
END INTERFACE

INTERFACE PrimToCons
  MODULE PROCEDURE PrimToCons
END INTERFACE

INTERFACE EvaluateFluxX
  MODULE PROCEDURE EvaluateFluxX
END INTERFACE

INTERFACE EvaluateFluxY
  MODULE PROCEDURE EvaluateFluxY
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: ExactFunction
PUBLIC :: SourceTerms
PUBLIC :: BoundaryConditions
PUBLIC :: TimeStep
PUBLIC :: ConsToPrim
PUBLIC :: PrimToCons
PUBLIC :: EvaluateFluxX
PUBLIC :: EvaluateFluxY
!-------------------------------------------------------------------------------!
!
!
!
!===============================================================================!
CONTAINS
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ExactFunction(WhichInitialCondition,t,x,Cons)
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: nDims
USE MOD_FiniteDifference2D_vars,ONLY: PI
USE MOD_FiniteDifference2D_vars,ONLY: MESH_X0
USE MOD_FiniteDifference2D_vars,ONLY: MESH_X1
USE MOD_FiniteDifference2D_vars,ONLY: MESH_SX
USE MOD_FiniteDifference2D_vars,ONLY: Gravity
USE MOD_FiniteDifference2D_vars,ONLY: PrimRefState1
USE MOD_FiniteDifference2D_vars,ONLY: PrimRefState2
USE MOD_FiniteDifference2D_vars,ONLY: PrimRefState3
USE MOD_FiniteDifference2D_vars,ONLY: PrimRefState4
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
INTEGER,INTENT(IN) :: WhichInitialCondition
REAL,INTENT(IN)    :: t
REAL,INTENT(IN)    :: x(1:nDims)
REAL,INTENT(OUT)   :: Cons(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: Prim(1:nVar)
REAL               :: xc(2), xm(2), r, r0, hl, hr
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

Cons = 0.0
Prim = 0.0
SELECT CASE(WhichInitialCondition)
  CASE(200)
    ! Constant State
    Prim(1) = 1.0
    Prim(2) = 0.0
    Prim(3) = 0.0
    CALL PrimToCons(Prim,Cons)
  
  CASE(211)
    ! Dam Break
    xm(1) = MESH_X0(1)+0.5*MESH_SX(1)
    xm(2) = MESH_X0(2)+0.5*MESH_SX(2)

    hl    = 2.5
    hr    = 0.5

    IF (x(1) .LE. xm(1)) THEN
      Prim(1) = hl
      Prim(2) = 0.0
      Prim(3) = 0.0
    ELSE
      Prim(1) = hr
      Prim(2) = 0.0
      Prim(3) = 0.0
    END IF

    CALL PrimToCons(Prim,Cons)
  CASE(212)
    ! Circular Dam Break
    xm(1) = MESH_X0(1)+0.5*MESH_SX(1)
    xm(2) = MESH_X0(2)+0.5*MESH_SX(2)
    xc(1) = x(1)-xm(1)
    xc(2) = x(2)-xm(2)
    r     = SQRT(xc(1)**2 + xc(2)**2)

    r0    = 2.5

    Prim(1) = 0.5
    Prim(2) = 0.0
    Prim(3) = 0.0

    IF (r .LE. r0) THEN
      Prim(1) = 2.5
      Prim(2) = 0.0
      Prim(3) = 0.0
    END IF

    CALL PrimToCons(Prim,Cons)
  CASE DEFAULT
    ErrorMessage = "Exact function not specified"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE ExactFunction
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE SourceTerms(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: S
USE MOD_FiniteDifference2D_vars,ONLY: MeshNodes
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER          :: ii, jj
!-------------------------------------------------------------------------------!

S = 0.0

!-------------------------------------------------------------------------------!
END SUBROUTINE SourceTerms
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE BoundaryConditions(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: nElemsX
USE MOD_FiniteDifference2D_vars,ONLY: nElemsY
USE MOD_FiniteDifference2D_vars,ONLY: nGhosts
USE MOD_FiniteDifference2D_vars,ONLY: MESH_X0
USE MOD_FiniteDifference2D_vars,ONLY: MESH_X1
USE MOD_FiniteDifference2D_vars,ONLY: MeshNodes
USE MOD_FiniteDifference2D_vars,ONLY: PrimRefState1
USE MOD_FiniteDifference2D_vars,ONLY: PrimRefState2
USE MOD_FiniteDifference2D_vars,ONLY: PrimRefState3
USE MOD_FiniteDifference2D_vars,ONLY: PrimRefState4
USE MOD_FiniteDifference2D_vars,ONLY: BoundaryConditionsType
USE MOD_FiniteDifference2D_vars,ONLY: U
USE MOD_FiniteDifference2D_vars,ONLY: V
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
INTEGER            :: idx_vx, idx_vy
REAL               :: x0, xc, xt
REAL               :: Prim_in(1:nVar), Prim_out(1:nVar)
REAL               :: Cons_in(1:nVar), Cons_out(1:nVar)
REAL               :: ConsRefState1(1:nVar), ConsRefState2(1:nVar)
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

idx_vx = 2
idx_vy = 3

!------------------------------!
! Left Boundary Conditions     !
!------------------------------!
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nElemsX-nGhosts+ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nGhosts-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState4(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = Cons_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState4(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = Cons_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nGhosts-ii+1,jj)
        U(idx_vx,-nGhosts+ii,jj) =-U(idx_vx,nGhosts-ii+1,jj)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Right Boundary Conditions    !
!------------------------------!
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY    
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,nElemsX-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState2(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = Cons_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState2(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = Cons_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,nElemsX-ii+1,jj)
        U(idx_vx,nElemsX+ii,jj) =-U(idx_vx,nElemsX-ii+1,jj)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Top Boundary Conditions      !
!------------------------------!
SELECT CASE(BoundaryConditionsType(3))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = U(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = U(1:nVar,ii,nElemsY-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState3(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = Cons_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState3(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = Cons_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = U(1:nVar,ii,nElemsY-jj+1)
        U(idx_vy,ii,nElemsY+jj) =-U(idx_vy,ii,nElemsY-jj+1)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Bottom Boundary Conditions   !
!------------------------------!
SELECT CASE(BoundaryConditionsType(1))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nElemsY-nGhosts+jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nGhosts-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState1(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = Cons_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState1(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = Cons_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nGhosts-jj+1)
        U(idx_vy,ii,-nGhosts+jj) =-U(idx_vy,ii,nGhosts-jj+1)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

DO jj=-nGhosts,nElemsY+nGhosts+1
  DO ii=-nGhosts,nElemsX+nGhosts+1
    CALL ConsToPrim(U(1:nVar,ii,jj),V(1:nVar,ii,jj))
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE BoundaryConditions
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION TimeStep()
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: U
USE MOD_FiniteDifference2D_vars,ONLY: CFL
USE MOD_FiniteDifference2D_vars,ONLY: MESH_DX
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: nElemsX
USE MOD_FiniteDifference2D_vars,ONLY: nElemsY
USE MOD_FiniteDifference2D_vars,ONLY: LambdaMaxX
USE MOD_FiniteDifference2D_vars,ONLY: LambdaMaxY
USE MOD_FiniteDifference2D_vars,ONLY: MIN_TIMESTEP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL    :: TimeStep
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL    :: FastestWaveX, FastestWaveY
REAL    :: Prim(1:nVar)
INTEGER :: ii, jj
!-------------------------------------------------------------------------------!

LambdaMaxX = 0.0
LambdaMaxY = 0.0
TimeStep = HUGE(1.0)

DO jj=1,nElemsY
  DO ii=1,nElemsX
    CALL ConsToPrim(U(1:nVar,ii,jj),Prim(1:nVar))
    CALL WaveSpeeds2D(Prim(1:nVar),FastestWaveX,FastestWaveY)
    LambdaMaxX = MAX(LambdaMaxX,ABS(FastestWaveX))
    LambdaMaxY = MAX(LambdaMaxY,ABS(FastestWaveY))
    TimeStep  = MIN(TimeStep,MESH_DX(1)/LambdaMaxX,MESH_DX(2)/LambdaMaxY)
  END DO
END DO

TimeStep = CFL*TimeStep

IF (TimeStep .LT. MIN_TIMESTEP) THEN
  TimeStep = MIN_TIMESTEP
END IF

!-------------------------------------------------------------------------------!
END FUNCTION TimeStep
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WaveSpeeds2D(Prim,fastestx,fastesty)
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: Gravity
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: fastestx
REAL,INTENT(OUT) :: fastesty
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: h, vx, vy
!-------------------------------------------------------------------------------!

h  = Prim(1)
vx = Prim(2)
vy = Prim(3)

fastestx = vx + SQRT(Gravity*h)
fastesty = vy + SQRT(Gravity*h)

!-------------------------------------------------------------------------------!
END SUBROUTINE WaveSpeeds2D
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ConsToPrim(Cons, Prim)
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: Gravity
USE MOD_FiniteDifference2D_vars,ONLY: MIN_DEPTH, MIN_SPEED
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Cons(1:nVar)
REAL,INTENT(OUT) :: Prim(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: h, hvx, hvy, ht
!-------------------------------------------------------------------------------!

h    = Cons(1)
hvx  = Cons(2)
hvy  = Cons(3)

IF (h .LT. MIN_DEPTH) THEN
  h = MIN_DEPTH
END IF
IF (ABS(hvx) .LT. MIN_SPEED) THEN
  hvx = 0.0
END IF
IF (ABS(hvy) .LT. MIN_SPEED) THEN
  hvy = 0.0
END IF

ht = h + MIN_DEPTH/h

Prim(1) = h
Prim(2) = hvx/ht
Prim(3) = hvy/ht

!-------------------------------------------------------------------------------!
END SUBROUTINE ConsToPrim
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE PrimToCons(Prim, Cons)
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: Gravity
USE MOD_FiniteDifference2D_vars,ONLY: MIN_DEPTH, MIN_SPEED
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: Cons(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: h, vx, vy
!-------------------------------------------------------------------------------!

h   = Prim(1)
vx  = Prim(2)
vy  = Prim(3)

IF (h .LT. MIN_DEPTH) THEN
  h = MIN_DEPTH
END IF
IF (ABS(vx) .LT. MIN_SPEED) THEN
  vx = 0.0
END IF
IF (ABS(vy) .LT. MIN_SPEED) THEN
  vy = 0.0
END IF

Cons(1) = h
Cons(2) = h*vx
Cons(3) = h*vy

!-------------------------------------------------------------------------------!
END SUBROUTINE PrimToCons
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE EvaluateFluxX(Prim,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: Gravity
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: h, vx, vy
!-------------------------------------------------------------------------------!

h  = Prim(1)
vx = Prim(2)
vy = Prim(3)

Flux(1) = h*vx
Flux(2) = h*vx**2 + 0.5*Gravity*h**2
Flux(3) = h*vx*vy

!-------------------------------------------------------------------------------!
END SUBROUTINE EvaluateFluxX
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE EvaluateFluxY(Prim,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: Gravity
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: h, vx, vy
!-------------------------------------------------------------------------------!

h  = Prim(1)
vx = Prim(2)
vy = Prim(3)

Flux(1) = h*vy
Flux(2) = h*vy*vx
Flux(3) = h*vy**2 + 0.5*Gravity*h**2

!-------------------------------------------------------------------------------!
END SUBROUTINE EvaluateFluxY
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_Equation
!-------------------------------------------------------------------------------!
