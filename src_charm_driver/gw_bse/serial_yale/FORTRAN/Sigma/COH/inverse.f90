! epsilon matrix inverse
! Oct. 2014  minjung.kim@yale.edu

! Using LAPACK
! for more information, visit 
! http://www.netlib.org/lapack/explore-html/d0/db3/zgetri_8f.html


subroutine inverse( n, M )

   use constant
   implicit none
   
   external ZGESV
   external ZGETRF
   external ZGETRI
   
   integer, intent(in) :: n
   complex*16, intent(inout) :: M(n,n)
   
   integer :: INFO
   integer, allocatable, dimension(:) :: IPIV
   complex*16, allocatable, dimension(:,:) :: WORK
   
   integer :: i
   
   
   allocate( WORK(n,n), IPIV(n) )
   
   WORK = cmplx( 0.d0, 0.d0 )
   do i = 1, n
      WORK(i,i) = cmplx( 1.d0, 0.d0 )
   enddo
   
   call ZGESV(n, n, M, n, IPIV, WORK, n, INFO)
   
   M = WORK

   deallocate( WORK, IPIV )


end subroutine