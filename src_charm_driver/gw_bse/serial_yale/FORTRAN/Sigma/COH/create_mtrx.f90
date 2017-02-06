! Oct 2014 minjung.kim@yale.edu

subroutine create_Pmtrx_struc( pol, FFTsize )


   use constant
   use electronic_structure
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


