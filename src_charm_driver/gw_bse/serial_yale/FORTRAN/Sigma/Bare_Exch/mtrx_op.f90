! Matrix operation for 3x3 matrix
! Oct 2014. minjung.kim@yale.edu



! Inverse 3x3 matrix
! 3x3 matrix inversion

! double precision and real numbers

subroutine dmatinv( m , inv_m ) 

   use constant
   implicit none
   real(dp), intent(in) :: m(3,3)
   real(dp), intent(inout) :: inv_m(3,3)
   real(dp) :: det, tr ! determinant and trace
   
   ! work variables:
   real(dp) :: a(3,3), del, x
   integer :: i, j
   !---------------------------------------------------------------
   !
   ! compute matrix of cofactors
   !
   a(1,1) = m(2,2)*m(3,3) - m(2,3)*m(3,2)
   a(2,1) = -m(2,1)*m(3,3) + m(2,3)*m(3,1)
   a(3,1) = m(2,1)*m(3,2) - m(2,2)*m(3,1)
   a(1,2) = -m(1,2)*m(3,3) + m(1,3)*m(3,2)
   a(2,2) = m(1,1)*m(3,3) - m(1,3)*m(3,1)
   a(3,2) = -m(1,1)*m(3,2) + m(1,2)*m(3,1)
   a(1,3) = m(1,2)*m(2,3) - m(1,3)*m(2,2)
   a(2,3) = -m(1,1)*m(2,3) + m(1,3)*m(2,1)
   a(3,3) = m(1,1)*m(2,2) - m(1,2)*m(2,1)
   
   ! compute determinant
   det = m(1,1)*a(1,1) + m(1,2)*a(2,1) + m(1,3)*a(3,1)
   ! compute trace
   tr = m(1,1) + m(2,2) + m(3,3)
   del = 1.0d-5
   if (abs(det) < del ) then
      call print_error_exit( 'Determinant of A is zero. Program exits.' )
   endif
   
   
   do i = 1, 3
      do j = 1, 3
         x = a(i,j) / det
         inv_m(i,j) = x
      enddo
   enddo
   
end subroutine


subroutine ddet( m , det )
   use constant
   implicit none
   real(dp), intent(in) :: m(3,3)
   real(dp), intent(inout) :: det

   real(dp) :: a(3,3), del, tr

   a(1,1) = m(2,2)*m(3,3) - m(2,3)*m(3,2)
   a(2,1) = -m(2,1)*m(3,3) + m(2,3)*m(3,1)
   a(3,1) = m(2,1)*m(3,2) - m(2,2)*m(3,1)
   a(1,2) = -m(1,2)*m(3,3) + m(1,3)*m(3,2)
   a(2,2) = m(1,1)*m(3,3) - m(1,3)*m(3,1)
   a(3,2) = -m(1,1)*m(3,2) + m(1,2)*m(3,1)
   a(1,3) = m(1,2)*m(2,3) - m(1,3)*m(2,2)
   a(2,3) = -m(1,1)*m(2,3) + m(1,3)*m(2,1)
   a(3,3) = m(1,1)*m(2,2) - m(1,2)*m(2,1)
   
   ! compute determinant
   det = m(1,1)*a(1,1) + m(1,2)*a(2,1) + m(1,3)*a(3,1)
   ! compute trace
   tr = m(1,1) + m(2,2) + m(3,3)
   del = 1.0d-5
   if (abs(det) < del ) then
      call print_error_exit( 'Determinant of A is zero. Program exits.' )
   endif
end subroutine




! matrix vector multiplication: a =  M * b
! double precision and real numbers

subroutine dmatvec( a, b, M )

   use constant
   implicit none
   
   real(dp), intent(inout) :: a(3)
   real(dp), intent(in) :: M(3,3), b(3)
   integer :: i, j
   
   a = 0.d0
   do i = 1, 3
      do j = 1, 3
         a(i) = M(i,j)*b(j) + a(i)
      enddo
   enddo
   

end subroutine
