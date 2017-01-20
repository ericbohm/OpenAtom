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
   
   ! rank 2 matrix
   type rank2_mtrx
   
      complex(dp),dimension(:,:),allocatable :: C
      real(dp),dimension(:,:),allocatable :: R
      type(gspace), dimension(:), allocatable :: gs
   
   end type
   
   
   type polarizability

      integer :: nq ! number of q vectors considered
      type(rank2_mtrx), dimension(:), allocatable :: mtrx
      
   end type 
   
   ! epsilon
   type epsilon_type
   
      integer :: nq   ! q index

      type(rank2_mtrx), dimension(:), allocatable :: mtrx   ! epsilon matrix
      type(rank2_mtrx), dimension(:), allocatable :: invm   ! and its inverse
      type(gspace), dimension(:), allocatable :: gs
   
   end type 


end module gw_structure
!==========================================================
