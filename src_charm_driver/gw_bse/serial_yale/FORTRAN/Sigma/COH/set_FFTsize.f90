! This subroutine set the size of FFT

subroutine set_FFTsize( psi, FFTsize )

   use constant
   use electronic_structure

   implicit none
   type(wfstruc), intent(inout) :: psi
   integer, intent(inout) :: FFTsize(3)

   ! work variables
   integer :: ig, ik, j
   real(dp) :: cutoff, tmpboxsize(3)
   integer :: box(3)  ! FFT box size for {Gmax} (density energy cutoff)   
   integer, allocatable :: kbox(:,:)  ! FFT box size for each k vectors
   integer :: kboxmax(3) ! the largest box

   ! Find the number of maximum points in each axis
   do j = 1, 3
      box(j) = maxval( psi%gvec(j,:) ) - minval( psi%gvec(j,:) ) + 1
   enddo
   
   ! Same as above, but this one is related to k and wavefunction cutoff defined by DFT calculations
   allocate( kbox(3,psi%nkpt) )
   do ik = 1, psi%nkpt
      do j = 1, 3
         kbox(j,ik) = maxval( psi%wk(ik)%gvec(j,:) ) - minval( psi%wk(ik)%gvec(j,:) ) + 1
      enddo
   enddo
   
   do j = 1, 3
      kboxmax(j) = maxval( kbox(j,:) )
   enddo

   ! here is the multiplier that sets the final FFTsize
   ! later, this should be read from an input file or something like that.
   cutoff = 0.9d0
   tmpboxsize(:) = dble( kboxmax(:) ) * cutoff

   FFTsize = int( tmpboxsize ) + 1
  !FFTsize = 14
   do j = 1, 3
      ! make it an even number
      if ( MOD( FFTsize(j), 2 ) .ne. 0 ) FFTsize(j) = FFTsize(j) + 1
   enddo
      
   print*, 'Size of FFT box:', FFTsize

end subroutine
