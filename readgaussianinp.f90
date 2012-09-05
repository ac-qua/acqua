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
  CHARACTER*16 A

  CHARACTER*8 X,METHOD,BASIS,TYPUS       !route

  CHARACTER*8 PGCAR, LABEL		!title
  INTEGER CHRG, MULT

  CHARACTER*8,ALLOCATABLE :: ATOM(:)	!molecule specs
  INTEGER, ALLOCATABLE :: XYZ(:,:)
  INTEGER n


 ! link 0 section incomplete

 10 READ(*,'(a)') A
    WRITE(*,*) "A:", A
   IF( A .ne. "" ) THEN
        IF ( SCAN(A,"%").eq.1 ) THEN
           WRITE(*,*) "get link 0 commands"
           !get link 0 commands  
        ELSE
           READ(A,*)X, METHOD,BASIS,TYPUS
              IF (X.ne."#") THEN 
                 WRITE(*,*) "ERROR: Something went terribly wrong, please check your input" 
              END IF
        END IF
     GOTO 10
   END IF


 WRITE(*,*) "X: ", X
 WRITE(*,*) "METHOD: ", METHOD
 WRITE(*,*) "BASIS: ", BASIS
 WRITE(*,*) "TYPUS: ", TYPUS

   !!Moegliche Fehlerquelle: Wird sofort nach leerer Zeile fortgefahren?
   !!				und ist "" richtig?


  !! ROUTE SECTION FEHLT! WIE TRENN ICH DIE VON DER LINK0 SECTION AB?



 ! title section
  READ(*,*) PGCAR,LABEL  ! in der title section können angeblichauch mehrere lines sein


 ! molecule specification
 
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
