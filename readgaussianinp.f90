MODULE readgaussianinp

 IMPLICIT NONE

 CONTAINS
 
 SUBROUTINE gaussianinput()

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! variables: - A,i,j: variables being used for list-directed input    !
 !              PGCAR: Pointgroup (CHARACTER)  		               !
 !		LABEL: not important for calculations, just eye-candy  !
 !		CHRG: charge of the molecule			       !
 ! 		MULT: spin multiplicity				       !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Input Structure:						       !
 !    \\  link 0 section					       !
 !    \\  route section (specifies job type and model)		       !
 !    \\  <blank line>						       !
 !    \\  title section						       !
 !    \\  <blank line>						       !
 !    \\  molecule specification / input structure 		       !
 !    \\  <blank line>						       !
 !    \\  variables section (defines vars used in molecule section)    !
 !    \\  <blank line>						       !
 !    \\  <EOF>							       !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! useful help for list-directed input in fortran:
 ! http://objectmix.com/fortran/242381-reading-file-containing-blank-lines.html

  CHARACTER*8 A, PGCAR, LABEL
  INTEGER CHRG, MULT

 ! link 0 section incomplete
 
 10 READ(*,*) A
   IF( A .ne. "" ) THEN
     !get link 0 commands  
     GOTO 10
   END IF

   !!Moegliche Fehlerquelle: Wird sofort nach leerer Zeile fortgefahren?
   !!				und ist "" richtig?


  !! ROUTE SECTION FEHLT! WIE TRENN ICH DIE VON DER LINK0 SECTION AB?

 ! title section
  READ(*,*) PGCAR,LABEL  ! in der title section können angeblichauch mehrere lines sein


 ! molecule specification
 CHARACTER*8,ALLOCATABLE :: ATOM(:)
 INTEGER, ALLOCATABLE :: XYZ(:,:)
 INTEGER n
 
 READ(*,*) CHRG,MULT
 
 n=1
 11 READ(*,*) A
   IF( A .ne. "" ) THEN
     READ(A,*) ATOM(n),XYZ(n,1),XYZ(n,2),XYZ(n,3)     
     n=n+1   !am Ende dieser Schleife ist n= Anzahl Atome + 1
     GOTO 11
   END IF

 ALLOCATE( ATOM(n-1) )
 ALLOCATE( XYZ(n-1,3) )

  
 END SUBROUTINE gaussianinput

END MODULE readgaussianinp
