subroutine scfalgorithm(number_of_atoms,noe,nao,slaterzeta,xyz,chrg,ncoeff,i,j,k)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! variables:	-noe, nao: number of electrons / atom orbitals  !
 !              -xyz, xyzauf: atom coordinates / aufpunkte      !
 !              -slaterzeta: zetas of AO's                      !
 !              -ncoeff: AO's of atoms                          !
 !              -limit: convergence limit                       !
 !              -z,c: primitive and contracted coefficients     !
 !              -kkenergy: energy of the nuclei (BO-Approx.)    !
 !              -edens: electrondensity              		!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 implicit none !real*8 (a-h,o-z)
 integer number_of_atoms, noe, nao  
 integer ncoeff,i,j,k,l,n,limit
 real*8 xyz,chrg,slaterzeta
 real*8 z(6,nao),c(6,nao),fockenergy,fockenergyalt,kkenergy,rAB,edens
 real*8, allocatable :: h00(:,:),s00(:,:),t00(:,:),v00(:,:),xyzauf(:,:),g(:,:,:,:),pac(:),eigenval(:),eigenvec(:,:),eigenvalmat(:,:),trans(:,:),tmp(:,:),fullG(:,:),densitymat(:,:),fock(:,:),fock2(:,:),cstrich(:,:),orbitalenergy(:),cohnestrich(:,:),abstandvec(:),&
& mulliken(:,:), mullikencharge(:),mulltmp(:), vecr(:)
 dimension slaterzeta(nao)
 dimension xyz(3,number_of_atoms) 
 dimension chrg(number_of_atoms) 
 dimension ncoeff(number_of_atoms) 

 allocate( s00(nao,nao) )
 allocate( t00(nao,nao) )
 allocate( v00(nao,nao) )
 allocate( h00(nao,nao) )
 allocate( xyzauf(3,nao) )
 allocate( abstandvec(3) )
 allocate( g(nao,nao,nao,nao) )
 allocate( pac(nao*(nao+1)/2) )
 allocate( eigenval(nao) )
 allocate( eigenvec(nao,nao) )
 allocate( eigenvalmat(nao,nao) )
 allocate( trans(nao,nao) )
 allocate( tmp(nao,nao) )
 allocate( fullG(nao,nao) )
 allocate( densitymat(nao,nao) )
 allocate( fock(nao,nao) )
 allocate( fock2(nao,nao) )
 allocate( cstrich(nao,nao) )
 allocate( orbitalenergy(nao) )
 allocate( cohnestrich(nao,nao) )
 allocate( mulliken(nao,nao) )
 allocate( mullikencharge(number_of_atoms) )
 allocate( mulltmp(number_of_atoms) )
 allocate( vecr(3) )
 
fockenergyalt=0
kkenergy=0

!write(*,*) 'numberofatoms:', number_of_atoms
!write(*,*) 'numberofelectrons:', noe
!write(*,*) 'numberofaos:', nao
!write(*,*) 'coordinates:', xyz
!write(*,*) 'charge:', chrg
!write(*,*) 'ncoeffs:', ncoeff
!write(*,*) 'zetas:', slaterzeta

do i=1,nao
    call slater(slaterzeta(i),z(1:6,i),c(1:6,i))
end do

!write(*,*) 'z:', z
!write(*,*) 'c:', c


!!!Aufpunkte produzieren
k=1
do i=1,number_of_atoms
  do j=1,ncoeff(i)
     xyzauf(1:3,k) = xyz(1:3,i)
     k=k+1
  end do
end do

!call prmat(6,xyzauf,3,nao,'Aufpunktmatrix:')

!!!calculate one-electron-integrals --> h00,s00

do i=1,nao
  do j=1,nao
    call stvint(6,6,number_of_atoms,xyz,chrg,xyzauf(1:3,i),xyzauf(1:3,j),z(1:6,i),z(1:6,j),c(1:6,i),c(1:6,j),s00(i,j),t00(i,j),v00(i,j))
    h00(i,j)=v00(i,j)+t00(i,j)
  end do
end do

!call prmat(6,v00,nao,nao,'v00:')
!call prmat(6,s00,nao,nao,'s00:')
!call prmat(6,h00,nao,nao,'h00:')
!call prmat(6,t00,nao,nao,'t00:')

!!!calculate two-electron-integrals --> g

do i=1,nao
  do j=1,nao
    do k=1,nao
      do l=1,nao
        call twoint(6,6,6,6,xyzauf(1:3,i),xyzauf(1:3,j),xyzauf(1:3,k),xyzauf(1:3,l),z(1:6,i),z(1:6,j),z(1:6,k),z(1:6,l),c(1:6,i),c(1:6,j),c(1:6,k),c(1:6,l),g(i,j,k,l))
      end do
    end do
  end do
end do

!! MAtrix packen
call matrixpac(s00,nao,pac)

!call rsp(pac,nao*(nao+1)/2,2,eigenval,eigenvec)
call rsp(pac,nao,2,eigenval,eigenvec)

!call prmat(6,eigenvec,nao,nao,'eigenvec:')
!write(*,*) 'Vektoren:', eigenvec
!write(*,*) 'Werte:', eigenval

!!Eigenvalmat erstellen; das ist hier schon hoch -0.5
do i=1,nao
  eigenvalmat(i,i)=(eigenval(i))**(-0.5)
end do

!call prmat(6,eigenvalmat,nao,nao,'Eigenvalmat:')

!! U * Eigenvalmat *  U

call matrixmulttrans(nao,tmp,eigenvalmat,eigenvec)
call matrixmult(nao,trans,eigenvec,tmp)

!call prmat(6,trans,nao,nao,'Transformationmatrix X:')

!! Xtransp * S * X = 1

!call matrixmulttrans(nao,tmp,
tmp = matmul(transpose(trans),matmul(s00,trans))

!call prmat(6,tmp,nao,nao,'XT*S*X=1?::')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! FOCK MATRIX BASTELN: !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
limit=50
densitymat=0
fockenergy=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Hier beginnt die Fock-IteratioN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do n=1,limit 

fullG=0

   do i=1,nao
     do j=1,nao
       do k=1,nao
         do l=1,nao
            fullG(i,j) = fullG(i,j) + densitymat(l,k)*(g(i,j,k,l)-0.5*g(i,l,k,j))
         end do
       end do
     end do
   end do

!call prmat(6,fullG,nao,nao,'FULLG:')

fock=h00+fullG

!call prmat(6,fock,nao,nao,'FOCK-MATRIX:')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Aus F F' machen und diagonalisieren um        !!! 
!!! C' und € zu erhalten +  Rücktrans von C' zu C !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call matrixmult(nao,tmp,fock,trans)
call matrixmulttrans2(nao,fock2,trans,tmp)

!call prmat(6,fock2,nao,nao,'FOCKSTRICH-MATRIX:')

call matrixpac(fock2,nao,pac)
call rsp(pac,nao,1,orbitalenergy,cstrich)

!call prmat(6,cstrich,nao,nao,'Cstrich:')
!write(*,*) 'Orbitalenergien:', orbitalenergy

call matrixmult(nao,cohnestrich,trans,cstrich)
!call prmat(6,cohnestrich,nao,nao,'C:')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Generiere neue P aus C                 !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   densitymat=0
   do k=1,(0.5*noe)
     do i=1,nao
       do j=1,nao
         densitymat(i,j)=densitymat(i,j)+2*cohnestrich(i,k)*cohnestrich(j,k)
       end do
     end do
   end do

   !call prmat(6,densitymat,nao,nao,'Densitymat:')

   !write(*,*) 'Orbitalenergien:', orbitalenergy


fockenergy=0
  do i=1,nao
    do j=1,nao
      fockenergy=fockenergy+0.5*densitymat(j,i)*(h00(i,j)+fock(i,j))
    end do
  end do


write(*,*) '----------------------------------------------------'
write(*,*) 'Iterationsschritt:', n
write(*,*) 'Fockenergie:', fockenergy
write(*,*) 'Änderung der Energie:', fockenergyalt-fockenergy

if(abs(fockenergyalt-fockenergy).gt.0.000001) then
fockenergyalt=fockenergy
else
exit
end if
end do

write(*,*) '----------------------------------------------------'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     ENDE FOCK-ITERATION      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Kern-Kern-Wecheselwirkung ausrechnen: Szabo (3.185)
kkenergy=0
if(number_of_atoms.gt.1) then

 do i=1,number_of_atoms
 do j=i+1,number_of_atoms
  if(i+1.le.number_of_atoms)  then
   abstandvec = xyz(1:3,i)-xyz(1:3,j)
!  write(*,*) 'Abstandsvektor:', abstandvec
  rAB = sqrt(abstandvec(1)**2+abstandvec(2)**2+abstandvec(3)**2)
!  write(*,*) 'rAB:', rAB
    kkenergy = kkenergy + chrg(i)*chrg(j) / rAB
  write(*,*) 'Kern-Kern-Energie:', kkenergy
  end if
 end do
 end do
fockenergy=fockenergy+kkenergy
end if

write(*,*) '----------------------------------------------------'
write(*,*) 'FOCK-ENERGY:', fockenergy
write(*,*) '----------------------------------------------------'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!          Mulliken-Populationsanalyse                         !!!!
!!!              siehe mullikenpop.f90                           !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


call mullikenpop(nao,mulliken,densitymat,s00,number_of_atoms,ncoeff,mulltmp,mullikencharge,chrg)
write(*,*) 'Mullikencharge:', mullikencharge
write(*,*) '----------------------------------------------------'


call edenscalc(nao,densitymat,z,xyzauf,edens)


deallocate( abstandvec )
deallocate( h00 )
deallocate( s00 )
deallocate( v00 )
deallocate( t00 )
deallocate( xyzauf )
deallocate( g )
deallocate( pac )
deallocate( eigenval )
deallocate( eigenvec )
deallocate( eigenvalmat )
deallocate( trans )
deallocate( tmp )
deallocate( fullG )
deallocate( densitymat )
deallocate( fock )
deallocate( fock2 )
deallocate( cstrich )
deallocate( orbitalenergy )
deallocate( cohnestrich )
deallocate( mulliken )
deallocate( mullikencharge )
deallocate( mulltmp )
deallocate( vecr )

end subroutine scfalgorithm

