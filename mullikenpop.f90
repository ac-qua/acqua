!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!          Mulliken-Populationsanalyse                         !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine mullikenpop(nao,mulliken,densitymat,s00,number_of_atoms,ncoeff,mulltmp,mullikencharge,chrg)
 implicit none

 integer nao,i,j,k,number_of_atoms,ncoeff
 real*8 mulliken,densitymat,s00,mulltmp,mullikencharge,chrg

dimension s00(nao,nao) 
dimension densitymat(nao,nao) 
dimension chrg(number_of_atoms)
dimension ncoeff(number_of_atoms)
dimension mulliken(nao,nao) 
dimension mullikencharge(number_of_atoms) 
dimension mulltmp(number_of_atoms) 


call matrixmult(nao,mulliken,densitymat,s00)
!call prmat(6,mulliken,nao,nao,'P*s00:')


k=1
mulltmp=0

do i=1,number_of_atoms
  do j=1,ncoeff(i)
     mulltmp(i) = mulltmp(i) + mulliken(k,k)
     k=k+1 
  end do
     mullikencharge(i)=chrg(i) - mulltmp(i)
end do

end subroutine mullikenpop
