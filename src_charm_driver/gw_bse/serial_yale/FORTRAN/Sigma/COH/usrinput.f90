! User provided input values defined here

! Nov. 2014 minjung.kim@yale.edu

! input file name: in (needs to be changed)

module usrinput

   use constant
   implicit none

   type input

      ! wavefunction file names
      character(len=100) :: wfname
      character(len=100) :: wfqname
      
      ! Decide whether the system is crystal or not
      logical :: crystal

      ! fft cutoff 
      real(dp) :: fftcut 
      ! epsilon cutoff (Ry)
      real(dp) :: epscut

   end type


   contains
   
   subroutine read_input( inp )
      implicit none
      type(input), intent(inout) :: inp
      integer, parameter :: iunit = 30


      open(iunit,file='in',status='old',form='formatted')
      
      inp%wfname = 'wfn.dat'
      inp%wfqname = 'wfnq.dat'
      inp%crystal = .true.
      inp%epscut = 14.d0

      close(iunit)


   end subroutine

end module

