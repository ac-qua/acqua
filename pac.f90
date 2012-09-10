subroutine matrixpac(matin,dim,vecout)
integer dim,i,j,k
real*8 :: matin,vecout

dimension matin(dim,dim)
dimension vecout(dim*(dim+1)/2)

k=1
do j=1,dim
do i=1,dim
  if (j.ge.i) then
    vecout(k)= matin(i,j)
    k=k+1
  end if
end do 
end do

end subroutine matrixpac
