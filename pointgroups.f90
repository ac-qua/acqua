MODULE pointgroups
 
 USE matmath
 USE prmat
 
 IMPLICIT NONE

 CONTAINS

  SUBROUTINE c2vz(NUMINAT,ELEM,POSMAT,NUMAT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calls: -matmultr                                                   !
  !        -matcomp                                                    !
  !                                                                    !
  ! variables: -NUMINAT (h)                                            !
  !            -ELEM(NUMINAT) (h)                                      !
  !            -POSMAT(3:NUMINAT) (h)                                  !
  !            -NUMAT: number of all atoms                             !
  !            -C2ZPOS: atoms by C2Z matrix                            !
  !                                                                    !
  ! transformation matrices:                                           !
  !  -C2Z: rotation about 180° around z axis                           !
  !              ( -1,  0,  0)                                         !
  !        C2Z = (  0, -1,  0)                                         !
  !              (  0,  0,  1)                                         !
  !  -SXZ: mirror plane along xz                                       !
  !              (  1,  0,  0)                                         !
  !        SXZ = (  0, -1,  0)                                         !
  !              (  0,  0,  1)                                         !
  !  -SYZ: mirror plane along yz                                       !
  !              ( -1,  0,  0)                                         !
  !        SYZ = (  0,  1,  0)                                         !
  !              (  0,  0,  1)                                         !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   INTEGER NUMINAT, ELEM(:), NUMAT,SMAT1,SMAT2,SMAT3
   REAL*8 POSMAT(:,:)
   REAL*8 C2Z(3,3), SXZ(3,3), SYZ(3,3)
   REAL*8, ALLOCATABLE :: C2ZPOS(:,:)
   REAL*8, ALLOCATABLE :: SXZPOS(:,:)
   REAL*8, ALLOCATABLE :: SYZPOS(:,:)
   REAL*8, ALLOCATABLE :: MAT1(:,:)
   REAL*8, ALLOCATABLE :: MAT2(:,:)
   REAL*8, ALLOCATABLE :: MAT3(:,:)

  !!! initialization !!!
  
   C2Z(1:3,1:3) = 0
   SXZ(1:3,1:3) = 0
   SYZ(1:3,1:3) = 0

  !!! construction of transformation matrices !!!
   
   C2Z(1,1) = -1
   C2Z(2,2) = -1
   C2Z(3,3) =  1
   SXZ(1,1) =  1
   SXZ(2,2) = -1
   SXZ(3,3) =  1
   SYZ(1,1) = -1
   SYZ(2,2) =  1
   SYZ(3,3) =  1

  !!! C2 rotation aroung z axis, mirror plane xz, mirror plane yz !!!
  
   CALL matmultr(C2Z,3,3,POSMAT,3,NUMINAT,C2ZPOS)
   CALL matmultr(SXZ,3,3,POSMAT,3,NUMINAT,SXZPOS)
   CALL matmultr(SYZ,3,3,POSMAT,3,NUMINAT,SYZPOS)
  !!! comparing pos matrices !!!

   CALL matcomp(POSMAT,3,NUMINAT,C2ZPOS,3,NUMINAT,MAT1,SMAT1)
   DEALLOCATE( C2ZPOS )
   CALL matcomp(MAT1,3,SMAT1,SXZPOS,3,NUMINAT,MAT2,SMAT2)
   DEALLOCATE( MAT1 )
   DEALLOCATE( SXZPOS )
   CALL matcomp(MAT2,3,SMAT2,SYZPOS,3,NUMINAT,MAT3,SMAT3)
   DEALLOCATE( MAT2 )
   DEALLOCATE( SYZPOS)
   CALL printmatr(MAT3,3,SMAT3)
  END SUBROUTINE c2vz

END MODULE
