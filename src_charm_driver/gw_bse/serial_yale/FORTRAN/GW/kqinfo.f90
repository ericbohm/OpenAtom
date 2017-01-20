! Set q vectors
!
! Oct. 2014  minjung.kim@yale.edu
!

! This subroutine gives you an information of k points 
! and q vectors 


subroutine get_qvec( k, q )

   use constant
   use electronic_structure
   implicit none
   
   type(kptinfo), intent(inout) :: k, q
   
   ! local variables
   integer :: ik
   
   ! if k is not shifted, then {k} = {q}
   ! if k is shifted from origin, then we need to set {q} 
   
   allocate( q%vec(3,k%nk) )
   ! Let's set up q vectors
   do ik = 1, k%nk
      q%vec(:,ik) = k%vec(:,ik) - k%vec(:,1)
   enddo
   
   ! number of k points
   q%nk = k%nk
   
end subroutine






subroutine get_kq_index( k, q, ik, iq, ikq, Uvec )

   use constant
   use electronic_structure

   implicit none

   type(kptinfo),intent(in) :: k, q
   integer, intent(in) :: iq, ik
   integer, intent(inout) :: ikq
   
   ! work variables
   real(dp), dimension(3) :: v1, v2, kplusq
   integer :: i, j
   real(dp), dimension(:,:), allocatable :: kvtmp
   logical :: idxfound
   real(dp) :: test
   ! for umklapp process
   integer, dimension(3) :: Uvec


   ! assign kvec at ik
   v1 = k%vec(:,ik)
   ! assign qvec at iq
   v2 = q%vec(:,iq)

   ! calculate kvec(ik)+qvec(iq)
   v2 = v1 + v2

   ! save k+q
   kplusq = v2


   ! if v2 is outside of BZ, let's get it back to BZ
   do i =1, 3
      if ( v2(i) .ge. dble( one ) ) then
         v2(i) = v2(i) - int(v2(i))
      elseif ( v2(i) .lt. dble( zero )  ) then
         v2(i) = v2(i) + dble( one )
      endif
   enddo
   

   ! calculate umklapp vector

   do i = 1, 3
      Uvec(i) = int( kplusq(i) - v2(i) )
   enddo


   ! mapping to original k vectors
   idxfound = .false. 
   do i = 1, k%nk
      test = distsq( v2, k%vec(:,i) )
      if ( test .lt. 1.0d-10 ) then
         ikq = i
         idxfound = .true.
      endif
   enddo

   if ( idxfound .eqv. .false. ) then
      call print_error_exit( 'k+q index is not found.')
   endif


contains

   real(dp) function distsq( a, b )

      implicit none
      real(dp),dimension(3),intent(in) :: a, b
      integer :: i,j,k
      real(dp) :: tmp

      tmp = 0.d0
      do i = 1, 3
         tmp = ( a(i)-b(i) ) **2 + tmp
      enddo

      distsq = tmp

   end function

end subroutine






! Quantum Espresso's output k vectors are in cartesian basis, not the reciprocal lattice vector basis.
subroutine cartesian_to_crystal( sys, k, k_ )
   
   use constant
   use electronic_structure

   implicit none
   type(sysinfo), intent(in) :: sys
   type(kptinfo), intent(in) :: k
   type(kptinfo), intent(inout) :: k_
   real(dp) :: inv_b(3,3)
   
   integer :: ii, ik, nk
   
   call dmatinv( sys%bvec, inv_b  )

   nk = k%nk
   k_%nk = nk
   allocate( k_%vec(3,nk) ) 
   k_%vec = 0.d0
   
   do ik = 1, nk
      call dmatvec( k_%vec(:,ik), k%vec(:,ik), inv_b )  ! sol, inp, mat
   enddo


end subroutine
