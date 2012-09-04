MODULE prmat

 IMPLICIT NONE

 CONTAINS

  SUBROUTINE printmatr(A,Z,S)
  
   REAL*8 A(:,:)
   INTEGER Z,S,I

   DO I = 1,Z
    WRITE(*,*) A(I,1:S)
   END DO

  END SUBROUTINE printmatr

  SUBROUTINE printveci(V,D)

   INTEGER V(:)
   INTEGER D,I

   DO I = 1,D
    WRITE(*,*) V(I)
   END DO
  
  END SUBROUTINE printveci

END MODULE
   
