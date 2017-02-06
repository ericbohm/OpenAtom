! User provided input values defined here

! Nov. 2014 minjung.kim@yale.edu

! input file name: in (needs to be changed)

module usrinput

   use constant
   use gw_structure
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
      real(dp) :: Pcut
      real(dp) :: epscut
      
      ! GPP related variables
      logical :: gpp_flag_is_on
      character(len=100) :: rhoName
      
      ! Full-frequency calculations
      logical :: ffreq_flag_is_on

   end type


   contains
   
   subroutine read_input( inp )

      type(input), intent(inout) :: inp
      integer, parameter :: iunit = 30


      open(iunit,file='in',status='old',form='formatted')
      
      read(iunit,*) inp%wfname
      read(iunit,*) inp%wfqname
      read(iunit,*) inp%Pcut
      read(iunit,*) inp%epscut
      !read(iunit,*) inp%gpp_flag_is_on
      !read(iunit,*) inp%rhoName
      close(iunit)


      ! Let's initialize all of the input parameters here
      
      !inp%wfname = 'wfn.dat'
      !inp%wfqname = 'wfnq.dat'
      inp%crystal = .true.
      !np%epscut = 15.d0
      inp%gpp_flag_is_on = .true.
      inp%rhoName = 'rho.dat'
      inp%ffreq_flag_is_on = .false.

   end subroutine
   

   !*** WARNING: this is temporary subroutine !
   subroutine set_full_frequency ( ffreq )
   
      type(full_frequency), intent(inout) :: ffreq
      
      ffreq%is_on = .false.
      ! energy is in eV unit !
      ffreq%Emin = 0
      ffreq%Emax = 5
      ffreq%Estep = 1
      ffreq%Ebrdn = 0.2
      ffreq%nEstep = (ffreq%Emax - ffreq%Emin) / ffreq%Estep + 1
   
   end subroutine
   
   subroutine set_gpp_frequency( gpp )
      
      type(gpp_type), intent(inout) :: gpp
      
      gpp%Emin = 0
      gpp%Emax = 15
      gpp%Estep = 1
      gpp%Ebrdn = 0.2
      gpp%nEstep = (gpp%Emax-gpp%Emin) / gpp%Estep  + 1
      
   end subroutine
   
   ! complete this subroutine
   !subroutine initialize( inp )
   !end subroutine 


end module
