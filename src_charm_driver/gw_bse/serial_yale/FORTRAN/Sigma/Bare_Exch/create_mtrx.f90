! Oct 2014 minjung.kim@yale.edu

subroutine create_Pmtrx_struc( pol, FFTsize )


   use constant
   use gw_structure

   implicit none
   type(polarizability), intent(inout) :: pol
   integer, intent(in) :: FFTsize(3)
   integer :: npt, nq, iq
      
   npt = FFTsize(1)*FFTsize(2)*FFTsize(3)
   nq = Pol%nq
   allocate( pol%mtrx(nq) )
   do iq = 1, nq
      allocate( pol%mtrx(iq)%C( npt, npt) )
      pol%mtrx(iq)%C(:,:) = 0.d0
   enddo
      
end subroutine

! for full frequency calculation
subroutine create_P_freq_struc( P_freq, ffreq, FFTsize )

   use constant
   use gw_structure
   type(full_polarizability), intent(inout) :: P_freq
   type(full_frequency), intent(in) :: ffreq
   integer, intent(in) :: FFTsize(3)
   integer :: npt, nq, iq
   
   npt = FFTsize(1)*FFTsize(2)*FFTsize(3)
   nq = P_freq%nq
   nE = ffreq%nEstep
   allocate( P_freq%mtrx(nE, nq ) )
   
   do iq = 1, nq
      do iE = 1, nE
         allocate( P_freq%mtrx(iE,iq)%C(npt,npt) )
         P_freq%mtrx(iE,iq)%C(:,:) = 0.d0
      enddo
   enddo
   
end subroutine
