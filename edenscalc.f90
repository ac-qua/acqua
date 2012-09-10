!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!          Elektronendichteberechnung                          !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine edenscalc(nao,densitymat,z,xyzauf,edens)
implicit none

integer i,j,k,nao
real*8 edens,densitymat,xyzauf,z
real*8, allocatable :: vecr(:)

dimension densitymat(nao,nao) 
dimension z(6,nao)
dimension xyzauf(3,nao) 

real*8 pi 


pi= 3.14159265358979d0

allocate( vecr(3) )

! Benutzereingabe vom Punkt r

!write(*,*) 'Geben sie den Punkt an, an dem sie die Elektronendichte wissen möchten:'
!do i=1,3
!  if(i.eq.1) then
!   write(*,*) 'X:'
!  else if(i.eq.2) then
!   write(*,*) 'Y:'
!  else if(i.eq.3) then
!   write(*,*) 'Z:'
!  read(*,*) vecr(i)
!  end if
!end do

vecr(1)=0
vecr(2)=0
vecr(3)=0



! Berechnung der Elektrondichte nach Gaußschen Produkttheorem; Szabo 3.207 & 3.203

edens=0


!Böser Fehler: Z ist eine Matrix
!Und Formel besser aufsplitten für übersichtlichkeit

do i=1,nao
  do j=1,nao      
    do k=1,3
   edens = edens + densitymat(i,j) * ((2*z(i)*z(j)/((z(i)+z(j))*pi))**(3/4) * exp( - z(i)*z(j) / (z(i) + z(j)) * abs(xyzauf(k,i) - xyzauf(k,j))**2)) * (2*(z(i)+z(j))/pi)**(3/4)*exp(-(z(i)+z(j))*abs(vecr(k)-(z(i)*xyzauf(k,i)+z(j)*xyzauf(k,j)/(z(i)+z(j))))**2) 
    end do
  end do           
end do



!do i=1,nao
!  do j=1,nao      
!   edens = edens + densitymat(i,j) * ((2*c(i)*c(j)/((c(i)+c(j))*pi))**(3/4) * exp( - c(i)*c(j) / (c(i) + c(j)) * abs(xyzauf(1:3,i) - xyzauf(1:3,j))**2)) * (2*(c(i)+c(j))/pi)**(3/4)*exp(-(c(i)+c(j))*abs(vecr(1:3)-(c(i)*xyzauf(1:3,i)+c(j)*xyzauf(1:3,j)/(c(i)+c(j))))**2) 
!    edens = edens + densitymat(i,j) * ((2*c(i)*c(j)/((c(i)+c(j))*pi))**(3/4) * exp( - c(i)*c(j) / (c(i) + c(j)) * abs(xyzauf(1:3,i) - xyzauf(1:3,j))**2)) * (2*(c(i)+c(j))/pi)**(3/4)*exp(-(c(i)+c(j))*abs(vecr(1:3)-(c(i)*xyzauf(1:3,i)+c(j)*xyzauf(1:3,j)/(c(i)+c(j))))**2) 
!    edens = edens + densitymat(i,j) * ((2*alpha*beta/((alpha+beta)*pi))**(3/4) * exp( - alpha*beta / (alpha + beta) * abs(R_A - R_B)**2)) * (2(alpha+beta)/pi)**(3/4) * exp(-(alpha+beta)*abs(r-R_A)**2)
!  end do           
!end do

deallocate( vecr )

end subroutine edenscalc

