! epsilon matrix inverse
! Oct. 2014  minjung.kim@yale.edu

! Using LAPACK


! n = number of G vectors in epsilon matrix
! Dimension of M : n x n
subroutine inverse( n, M )

   use constant
   implicit none
   
   external ZGESV
   
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





subroutine inverse_iterative( n, M )

   use constant
   implicit none
   
   integer, intent(in) :: n
   complex(dp), intent(inout) :: M(n,n)
   ! local variables
   integer :: i, j, iter, nmax
   real(dp), allocatable :: mv(:)
   real(dp) :: alpha, tmp
   complex(dp), allocatable :: X(:,:)
   integer, allocatable :: eye(:,:)
   
   allocate( mv(n), X(n,n), eye(n,n) )
   
   ! maximum number of iteration
   nmax = 30
   
   ! fill identity matrix
   eye = 0.d0
   do i = 1, n
      eye(i,i) = 1.d0
   enddo
   
   ! find alpha   
   X = M * transpose( conjg( M ) )
   
   alpha = 0.d0
   do i = 1, n
      do j= 1, n
         mv(j) = abs( X(i,j) )
      enddo
      tmp = sum(mv)
      if ( tmp .ge. alpha ) alpha = tmp
   enddo
   deallocate( mv )
   
   print*, 'alpha'
   alpha = 1/alpha
   X = alpha * transpose( conjg( M ) )
   
   do iter = 1, nmax
      X = MATMUL( X,  2*eye - MATMUL( M, X ) )
   enddo
   
   M(:,:) = X(:,:)
   
   deallocate(eye,X)
   
end subroutine