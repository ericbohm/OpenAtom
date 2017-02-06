! This subroutine deals with the shifted k grid wavefunctions
! Relevant for q=0

! Nov. 2014 minjung.kim@yale.edu

subroutine shifted_k_wfn( FFTsize, inp, psi, psi_R, k, sys )

   use constant
   use electronic_structure
   use usrinput

   implicit none
   integer, intent(in) :: FFTsize(3)
   type( input ), intent(in) :: inp
   type( wfstruc ), intent(inout)  :: psi, psi_R
   type( kptinfo ), intent( inout ) :: k
   type( sysinfo ) :: sys 

   ! READ wavefunctions with shifted grid
   call read_wfn( inp%wfqname, sys, psi, k )

   ! 
   call FFT_all_bands( psi, psi_R, FFTsize, sys )


end subroutine
