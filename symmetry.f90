MODULE symmetry

 USE prmat
 USE pointgroups
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! contains following subroutines:                                      !
!  -pointgroup(): for SYSDIM=0; reads pointgroup and atom positions;   !
!  -readatomposca(NUMINAT,ELEM,POSMAT): reads number of inequivalent   ! 
!   atoms [NUMINAT], elements [ELEM()] and position of atoms           !
!   [POSMAT(3,NUMINAT)] as cartesian coordinates in Angstrom           !
!  -readposmatc(NUMINAT,ELEM,POSMAT): reads elements [ELEM(NUMINAT)]   !
!   and positions [POSMAT(3,NUMINAT)] of inequivalent atoms in         !
!   Cartesian inputs.                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 IMPLICIT NONE
 
 CONTAINS

  SUBROUTINE pointgroup()
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calls readatomposca()                                              !
  !       [pointgroup](NUMINAT,ELEM,POSMAT,NUMAT)                      !
  !                                                                    !
  ! variables: -PG: point group number                                 ! 
  !            -NUMINAT: number of inequivalent atoms                  !
  !            -NUMAT: number of all atoms                             !
  !            -CORDTYPE: information about type of inp. coord.        !
  !            -ELEM(:): vector containing the element numbers         !
  !            -POSMAT(3,:): matrix containing the atom positions      !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   INTEGER PG, NUMINAT, NUMAT
   CHARACTER*8 CORDTYPE
   INTEGER, ALLOCATABLE :: ELEM(:)
   REAL*8, ALLOCATABLE :: POSMAT(:,:)

   READ(*,*) PG
   READ(*,*) NUMINAT
   READ(*,*) CORDTYPE
   WRITE(*,*) "Point group:",PG
   WRITE(*,*) "Number of inequivalent atoms:",NUMINAT

   IF(CORDTYPE.EQ."CARTA") THEN
    WRITE(*,*) "Atom positions given in Cartesian coordinates [Ang]."
    CALL readatomposca(NUMINAT,ELEM,POSMAT)
   END IF

   IF(PG.EQ.15) THEN
    WRITE(*,*) "Point group C2v(z)."
    CALL c2vz(NUMINAT,ELEM,POSMAT,NUMAT)
   END IF
   
  END SUBROUTINE pointgroup

  SUBROUTINE readatomposca(NUMINAT,ELEM,POSMAT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calls readposmatc()                                                !
  !                                                                    !
  ! variables: -NUMINAT: number of inequivalent atoms                  !
  !            -ELEM(:): vector containing the element numbers         !
  !            -POSMAT(3,:): matrix containing the atom positions      !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   INTEGER NUMINAT
   INTEGER, ALLOCATABLE :: ELEM(:)
   REAL*8, ALLOCATABLE :: POSMAT(:,:)

   ALLOCATE( ELEM(NUMINAT) )
   ALLOCATE( POSMAT(3,NUMINAT) )
   
   POSMAT(1:3,1:NUMINAT) = 0
   
   CALL readposmatc(NUMINAT,ELEM,POSMAT)
   CALL printveci(ELEM,NUMINAT)
   CALL printmatr(POSMAT,3,NUMINAT)
   
  END SUBROUTINE readatomposca

  SUBROUTINE readposmatc(NUMINAT,ELEM,POSMAT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! reads cartesian geometry inputs                                    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   INTEGER NUMINAT, ELEM(:), I, E
   REAL*8 POSMAT(:,:), X, Y, Z

   DO I=1,NUMINAT
    READ(*,*) E, X, Y, Z
    ELEM(I) = E
    POSMAT(1,I) = X
    POSMAT(2,I) = Y
    POSMAT(3,I) = Z
   END DO
  
  END SUBROUTINE readposmatc
  
END MODULE
