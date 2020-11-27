!===============================================================================!
MODULE MOD_Mesh
!===============================================================================!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE BuildMesh
  MODULE PROCEDURE BuildMesh
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: BuildMesh
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
SUBROUTINE BuildMesh()
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: MESH_SX
USE MOD_FiniteDifference2D_vars,ONLY: MESH_X0
USE MOD_FiniteDifference2D_vars,ONLY: MESH_X1
USE MOD_FiniteDifference2D_vars,ONLY: MESH_DX
USE MOD_FiniteDifference2D_vars,ONLY: nDims
USE MOD_FiniteDifference2D_vars,ONLY: nElemsX
USE MOD_FiniteDifference2D_vars,ONLY: nElemsY
USE MOD_FiniteDifference2D_vars,ONLY: MeshNodes
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj
!-------------------------------------------------------------------------------!

MeshNodes = 0.0

Mesh_SX    = ABS(Mesh_X1-Mesh_X0)
Mesh_DX(1) = ABS(Mesh_SX(1))/(REAL(nElemsX-1))
Mesh_DX(2) = ABS(Mesh_SX(2))/(REAL(nElemsY-1))

DO jj=1,nElemsY
  DO ii=1,nElemsX
    MeshNodes(1:nDims,ii,jj) = Mesh_X0(1:nDims) &
                             + (/REAL(ii-1),REAL(jj-1)/)*Mesh_DX(1:nDims)
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE BuildMesh
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_Mesh
!-------------------------------------------------------------------------------!
