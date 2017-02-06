! STRUCTURE MODULE 
!                     
!    Sep. 2014   minjung.kim@yale.edu

!==========================================================
module electronic_structure

   use constant

   ! k points info
   type kptinfo

      integer :: nk ! number of k poitns
      real(dp), dimension(:,:), allocatable :: vec ! in reciprocal space
      real(dp), dimension(:), allocatable :: wt ! weight factor
      real(dp), dimension(3) :: shift   ! shift 

   end type

   ! simulation cell information 
   type, extends(kptinfo) :: sysinfo
   
      integer :: ng      ! number of total g vectors
      integer, dimension(:,:), allocatable :: gvec   ! g vectors
      
      real(dp) :: alat                   ! lattice constant
      real(dp), dimension(3,3) :: avec   ! lattice vectors
      real(dp), dimension(3,3) :: bvec   ! reciprocal lattice vectors
      real(dp) :: vol ! volume of the simulation cell
                  
   end type 

   ! g vectors 
   type gspace
      
      integer :: ng         ! total number of g 
      integer, dimension(:,:), allocatable :: gvec   ! g vectors
      integer, dimension(3) :: Nfft      ! # of FFT grid in each direction

   end type

   
   ! wavefunction data at k
   type, extends(gspace) :: wfn
      
      integer :: nb     ! number of bands
      integer :: nv     ! number of valence bands
      integer :: nc     ! number of conduction bands
      integer :: kidx   ! k index
 
      real(dp), dimension(:), allocatable :: eig  ! eigenvalue (energy)
      real(dp), dimension(:), allocatable :: occ  ! occupancy
      
      complex(dp), dimension(:,:), allocatable :: cg  ! Fourier coefficient, for all band index
      
   end type 


   ! total wave function structure (we need it because we save all data on memory) 
   type, extends(gspace) :: wfstruc

      integer :: nband  ! number of total bands
      integer :: nspin  ! number of spin
      integer :: nkpt   ! number of k vectors
      type(wfn), dimension(:), allocatable :: wk

   end type


end module electronic_structure
!==========================================================



!==========================================================
module gw_structure

   use constant
   use electronic_structure
   implicit none

   ! 1D array
   type onedim_array
   
      real(dp), dimension(:), allocatable :: R
      complex(dp), dimension(:), allocatable :: C

   end type
   
   ! rank 2 matrix
   type rank2_mtrx
   
      complex(dp),dimension(:,:),allocatable :: C
      real(dp),dimension(:,:),allocatable :: R
   
   end type
   
   
   type polarizability

      integer :: nq ! number of q vectors considered
      type(rank2_mtrx), dimension(:), allocatable :: mtrx
      
   end type 
   
   ! epsilon
   type epsilon_type
   
      integer :: nq   ! number of q vectors
      ! epsilon matrix and its inverse
      type(rank2_mtrx), dimension(:), allocatable :: mtrx
      type(gspace), dimension(:), allocatable :: gs
   
   end type 
   
   
   ! Generalized Plasmon-Pole model
   type gpp_type
   
      integer :: nq  ! number of q vectors
      integer, dimension(:), allocatable :: ng
      type(rank2_mtrx), dimension(:), allocatable :: eigvec
      type(onedim_array), dimension(:), allocatable :: eigval
      type(onedim_array), dimension(:), allocatable :: omsq
      type(rank2_mtrx),dimension(:,:), allocatable :: mtrx
      real(dp) :: Emin
      real(dp) :: Emax
      real(dp) :: Estep
      real(dp) :: Ebrdn
      integer  :: nEstep
      
   end type


   ! User provided value for the full frequency calculation
   type full_frequency
      logical  :: is_on  ! true if full frequency is on, default if false
      real(dp) :: Emin   ! The lowest evaluation energy, set 0 as default 
      real(dp) :: Emax   ! The highest evaluation energy
      real(dp) :: Estep  ! increment of the evaluation energy
      real(dp) :: Ebrdn  ! Lorentzian broadening, set 0.2 as default
      integer  :: nEstep !
   end type

   type full_polarizability
      integer :: nq
      type(rank2_mtrx), dimension(:,:), allocatable :: mtrx
   end type
   
   type epsilon_ffreq 
      integer :: nq
      type(rank2_mtrx), dimension(:,:), allocatable :: mtrx   
   end type



end module gw_structure
!==========================================================
