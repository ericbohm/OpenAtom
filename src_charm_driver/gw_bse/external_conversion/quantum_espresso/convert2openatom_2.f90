!
! Sep. 2014   minjung.kim@yale.edu
!
! Last modified: Nov. 2015 by MK


Program converter

   use qexml_module

   implicit none

   character(len=5) :: codename = 'PW2OA'
   integer, parameter :: iunit = 10
   integer, parameter :: dp = kind(1.0d0)

   type cell_info
      real(dp) :: alat
      real(dp) :: a1(3), a2(3), a3(3)
      real(dp) :: b1(3), b2(3), b3(3)
      character(len=256) :: aunit, bunit
   end type
   
   type(cell_info) :: cell
  
   integer :: ngk_tot, nb
   integer, allocatable :: ipwk(:), igk_all(:)

   character(256) :: work_dir, prefix, dirname, filename, fsysname
   logical :: shift_flag ! whether wavefunctions are shifted or not
   integer :: nk ! it should be decided at first
   integer :: nspin

   ! total number of g vectors
   integer :: ngm
   ! master gvector list
   integer, allocatable :: master_gv(:,:)
   ! eigenvalues and its occupancies   
   real(dp), allocatable :: eig(:,:,:), occ(:,:,:)
   ! wavefunction
   complex(dp), allocatable :: wfn(:)
   ! k points and their weight
   real(dp), allocatable :: xk(:,:), wk(:)
   ! fftgrid for density (dense fft grid)
   integer :: fftsize(3)
   ! wavefunction cutoff
   real(dp) :: ecutwfc
   ! array to collect all g list in all k points
   integer, allocatable :: gkidx_allkpt(:)
   ! kmax_cp
   integer :: kmax_cp(3)
   ! number of g vectors
   integer :: ngdoublepack, ngkpt
   ! doublepack has only half sphere
   logical :: doublepack
   ! new g list that fits to OpenAtom
   integer, allocatable :: glist_doublepack(:,:), glist_kpt(:,:)
   ! number of planewaves at each k points (same numbers!)
   integer :: ncoeff 
   
   integer :: ib, ispin, ik, npwk
   integer, allocatable :: idxgk(:)

   ! GPP calculation related variables
   integer :: nr(3)
   real(dp), allocatable :: rho(:,:,:)
   logical :: gpp_flag
   
   ! Vxc variables
   logical :: Vxc_flag
   character(256) :: fname_vxc

   integer :: stdin = 5

   ! some working variables
   integer :: i, j
   integer :: ierr

   NAMELIST /INPUT/ prefix, work_dir, fsysname, doublepack, shift_flag, gpp_flag, Vxc_flag

!---------------------------------------------------------------------
! starts here
!---------------------------------------------------------------------

   ! default setting
   work_dir = './'
   doublepack = .false.
   shift_flag = .false.
   gpp_flag = .false.
   Vxc_flag = .false.
   fname_vxc = 'Vxcr.dat'
   
   read( stdin, INPUT, iostat=ierr )
   
   ! initialize QEXML library
   
   dirname = trim(work_dir) // '/' // trim(prefix) // '.save/'
   call qexml_init( iunit, Dir=dirname )

   filename = trim(dirname)//"data-file.xml"
   
   call qexml_openfile( filename, "read", IERR=ierr )


   ! read cell information
   call get_cell_info( cell )

   ! read number of k points = nk
   call qexml_read_bands_info( NUM_K_POINTS=nk, IERR=ierr)

   ! read BZ: k points and their weight factor
   allocate( xk(3,nk), wk(nk) )
   call qexml_read_bz( XK=xk, WK=wk, IERR=ierr )

   ! get master g vector
   call qexml_read_planewaves( NGM=ngm, IERR=ierr )
   allocate( master_gv(3,ngm) )
   call qexml_read_planewaves( IGV=master_gv, IERR=ierr )

   ! number of gvectors in each k
   allocate( ipwk( nk ) )
   call read_pw_index( nk, ipwk )

   ! Read number of bands and number of spins
   call qexml_read_bands_info( NBND=nb, NSPIN=nspin, IERR=ierr)

   ! read eigenvalues and occupancies
   allocate( eig(nb, nk, nspin), occ(nb, nk, nspin) )
   call read_eig_occ( nspin, nk, nb, eig, occ )

   ! read FFTsize (dense FFT grid for density)
   call qexml_read_planewaves( NR1=fftsize(1), NR2=fftsize(2), NR3=fftsize(3), IERR=ierr )

   ! read wavefunction cutoff
   call qexml_read_planewaves( ECUTWFC=ecutwfc, IERR=ierr )

   !---------------------------------------------------------------------------
   !   here it starts setting the same number of k points at each k points
   !---------------------------------------------------------------------------
   ! Adjust ecutwfc so that all g vectors at each k points are inside of the new ecutwfc
   ! of course, 0.5*|g+k|^2 < ecutwfc(original) but g itself can be outside of ecutwfc(original)
   call adjust_ecutwfc( cell, ecutwfc, nk, xk )

   ! set gkidx_allkpt and initialize
   allocate( gkidx_allkpt( ngm ) )
   gkidx_allkpt(:) = 0 ! 0 - do not include, 1 - include

   do ik = 1, nk
      ! read number of plane-waves at this k point
      call qexml_read_gk( ik, NPWK=npwk, IERR=ierr)
      allocate( idxgk( npwk ) )
      ! read plane-wave index to idxgk. idxgk connects to master_gv
      call qexml_read_gk( ik, index=idxgk, IERR=ierr )
      ! now, update gkidx_allkpt
      call set_gkidx_allkpt( npwk, idxgk, ngm, master_gv, gkidx_allkpt )

      deallocate( idxgk )
   enddo
   ! finally, we got the number of planewaves at each k point (each k point has the same number of planewaves)
   ncoeff = sum(gkidx_allkpt)
   ! print out
   print*, "(new) number of planewaves in state.out is:", ncoeff



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ THIS IS FOR KPT COMPARISON FOR WFNMINIMIZATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Here we set kmax_cp
   ! kmax_cp is the largest integer value in kx,ky,kz direction
   kmax_cp = 0 
   do ispin = 1, nspin
      do ik = 1, nk
         ! read number of plane-waves at this k point
         call qexml_read_gk( ik, NPWK=npwk, IERR=ierr )
         allocate( idxgk( npwk ) )
         ! read planewave index to idxgk. idxgk connects to master_gv
         call qexml_read_gk( ik, index=idxgk, IERR=ierr )
         ! now, we want to convert wavefunctions to openatom output format
         ! all k points should include the same number of g vectors
         call get_kmax( ngm, master_gv, npwk, idxgk, kmax_cp )
         deallocate( idxgk )
      enddo
   enddo
   ! to be safe, let's add 1 in each direction
   kmax_cp(:) = kmax_cp(:) + 1
   print*, 'kmax_cp is : ', kmax_cp(1), kmax_cp(2), kmax_cp(3)

   ! let's set the list of G vectors here
   ! warning: gamma point will contain only half sphere if doublepack flag is on
   if (doublepack) then
      call countkvec3d_sm( kmax_cp, ngdoublepack, ecutwfc, cell )
      allocate( glist_doublepack(3,ngdoublepack) )
      call setkvec3d_sm_simple( kmax_cp, ngdoublepack, ecutwfc, cell, glist_doublepack )
   endif
   
   call countkvec3d_sm_kpt( kmax_cp, ngkpt, ecutwfc, cell )
   allocate( glist_kpt(3,ngkpt) )
   call setkvec3d_sm_kpt( kmax_cp, ngkpt, ecutwfc, cell, glist_kpt)

   ! we need to setup gkidx_allkpt again for this particular test.
   ! it could be wrtten as a separate subroutine!!!
   ! initialized gkidx_allkpt AGAIN!!!
   gkidx_allkpt = 0
   do i = 1, ngm
      do j = 1, ngkpt ! loop over new g list set by setkvec3d_sm_kpt
         if( master_gv(1,i) == glist_kpt(1,j) .and. &
             master_gv(2,i) == glist_kpt(2,j) .and. &
             master_gv(3,i) == glist_kpt(3,j) ) then
            gkidx_allkpt(i) = 1
         endif
      enddo
   enddo
      
   print*, 'ngkpt:', ngkpt
   print*, 'sum gkidx_allkpt:', sum(gkidx_allkpt)
   print*, 'ecutwfc:', ecutwfc
   open(22,file='kptdata',status='replace',form='formatted')
   write(22,*) ngkpt
   do ik = 1, ngkpt
      write(22,*) glist_kpt(1:3,ik)
   enddo
   close(22)
   print*, 'done'


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END ROUTINES FOR KPT COMPARISON ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


   ! read and write wfc into state.out file
   do ispin = 1, nspin
      do ik = 1, nk
         ! get number of plane-waves 
         call qexml_read_gk( ik, NPWK=npwk, IERR=ierr)
         allocate( idxgk( npwk ) )
         ! read plane-wave index to idxgk. idxgk connects to master_gv
         call qexml_read_gk( ik, index=idxgk, IERR=ierr )
         do ib = 1, nb
            ! if k is at gamma point and doublepack flag is on
            ! it saves only half sphere
            if ( ik==1 .and. doublepack .eqv. .true. ) then
               call read_write_wfn( ib, ik, ispin, npwk, idxgk, master_gv, ngm,&
                    fftsize, shift_flag, gkidx_allkpt, doublepack)
            ! if not gamma point or not doublepack
            else
               call read_write_wfn( ib, ik, ispin, npwk, idxgk, master_gv, ngm,&
                    fftsize, shift_flag, gkidx_allkpt, doublepack)
            endif
         enddo! ib loop
         deallocate( idxgk )
      enddo! ik loop
   enddo! is loop
   ! wfn order: is-ik-ib-ig
   
   call qexml_closefile( "read", IERR=ierr )

   
   if ( shift_flag .eqv. .false.) then
      ! Write system information
      call write_system_info( cell, ncoeff, xk, wk, &
                       nspin, nk, nb, fftsize, fsysname )
   else
      continue
   endif
                          
   ! write eigenvalues and occupation numbers
   do ispin = 1, nspin
      do ik = 1, nk
         call write_eig_occ(ispin, ik, nb, eig(:,ik,ispin), occ(:,ik,ispin), shift_flag)
      enddo
   enddo   

   
! If RHO flag is on, create rho.dat file for GPP calculations
   if ( gpp_flag ) then
      call qexml_read_rho( NR1=nr(1), NR2=nr(2), NR3=nr(3), IERR=ierr )
      allocate( rho(nr(1),nr(2),nr(3)) )
      call write_rho( nr, rho )
   endif

! If Vxc_flag is on, create vxc.dat file for sigma calculations
!  if ( Vxc_flag ) then
!     call write_vxc( fname_vxc, codename  )
!  endif


   deallocate(xk, wk, eig, occ, master_gv )





!****************************************************************!
!
!                        SUBROUTINES
!
!****************************************************************!
contains


!************ subroutines ********************************************

!subroutine write_vxc( fname_vxc, codename )
!  ! adopted from pw2bgw.f90 
!  use environment,          only: environment_start, environment_end
!  use fft_base,             only: dfftp
!  use fft_interfaces,       only: fwfft
!  use ener,                 only: etxc, vtxc
!  use gvect,                only: ngm, ngm_g, ig_l2g, nl, mill
!  use lsda_mod,             only: nspin
!  use scf,                  only: rho, rho_core, rhog_core
!  use wavefunctions_module, only: psic

!  character(len=100), intent(in) :: fname_vxc
!  character(len=5), intent(in) :: codename
!  integer, parameter :: iu = 13

!  ! local variables
!  integer :: id, ig, is, ir
!  integer :: nr, ns, nd
!  integer, allocatable :: gvec(:,:)
!  
!  real(dp), allocatable :: vxcr(:,:)
!  complex(dp), allocatable :: vxcg(:,:)
!  
!  
!  ! assign variables
!  ns = nspin       ! number of spin
!  nr = dfftp%nnr   ! number of r grid
!  nd = 3           ! dimension


!  open(unit=iu,file=trim(fname_vxc),status='replace',form='formatted')

!  call environment_start( codename )
!  
!  call read_file()
!  
!  allocate( vxcr( nr, ns ) )
!  vxcr(:,:) = 0.d0
!  rho_core(:) = 0.d0
!  rhog_core(:) = 0.d0
!  call v_xc( rho, rho_core, rhog_core, etxc, vtxc, vxcr )
!  

!  
!  ! ngm is the number of g vectors in this process (with gamma tricks, only G>=0)
!  ! ngm_g is the total number of g vectors
!  ! if serial, ngm = ngm_g
!  allocate( gvec( 3, ngm ) )
!  gvec = 0.d0
!  do ig = 1, ngm
!     do id = 1, nd
!        ! ig_l2g converts the local G-vector index into the global index
!        ! mill is a miller index
!        gvec( id, ig_l2g(ig) ) = mill( id, ig )
!     enddo
!  enddo
!  
!  ! now we need to get vxcg
!  ! nl: fft index to G index
!  allocate( vxcg( ngm_g, ns ) )
!  do is = 1, ns
!     do ir = 1, nr
!        psic(ir) = complex( vxcr(ir,is), 0.d0 )
!     enddo
!     ! forward Fourier transform
!     call fwfft( 'Dense', psic, dfftp )
!     do ig = 1, ngm
!        vxcg( ig_l2g(ig), is ) = psic( nl(ig) )
!     enddo
!  enddo  
!  
!  do is = 1, ns
!     do ig = 1, ngm_g
!        write(iu,*) vxcg( ig, is ), ( gvec(i,ig), i=1,3 )
!     enddo
!  enddo
!  
!  deallocate( gvec, vxcr, vxcg )
!  
!  call environment_end( codename ) 
!  
!end subroutine


  
subroutine write_eig_occ( ispin, ik, nstate, eig, occ, shift_flag ) 

   integer, intent(in) :: ispin, ik, nstate
   real(dp), dimension(nstate), intent(in) :: eig, occ
   logical, intent(in) :: shift_flag
   character(len=100) :: fdir, fname
   
   integer :: iunit = 30
   
   ! directory name
   if (shift_flag .eqv. .false.) then
      if ( ik .le. 10 ) then
         write( fdir, '( "./STATES_IN/Spin.", I1, "_Kpt.", I1, "_Bead.0_Temper.0/" )' ) ispin-1, ik-1
      elseif ( ik .le. 100 ) then
         write( fdir, '( "./STATES_IN/Spin.", I1, "_Kpt.", I2, "_Bead.0_Temper.0/" )' ) ispin-1, ik-1   
      else
         print*, 'Cannot locate directory because you have more than 100 k points'
         stop
      endif
   elseif (shift_flag .eqv. .true. ) then
      if ( ik .le. 10 ) then
         write( fdir, '( "./STATES_IN/Spin.", I1, "_Kpt.0", I1, "_Bead.0_Temper.0/" )' ) ispin-1, ik-1
      elseif ( ik .le. 100 ) then
         write( fdir, '( "./STATES_IN/Spin.", I1, "_Kpt.0", I2, "_Bead.0_Temper.0/" )' ) ispin-1, ik-1   
      else
         print*, 'Cannot locate directory because you have more than 100 k points'
         stop
      endif
   endif   
   
   
   fname = 'eigenvalues.in'
   
   
   ! open file
   
   open(iunit, file=trim(fdir)//trim(fname), form='formatted', status='unknown')
   
   do i = 1, nstate
      write(iunit,*) eig(i)  !, occ(i) ! openatom won't read occupation 
   enddo
   close(iunit)

end subroutine





  
  


!------------------------------------------------------------------------------
! name of the output file: state\\ib\\.out.gz
subroutine read_write_wfn( ib, ik, ispin, npwk, idxgk, master_gv, ngm, &
                           fftsize, shift_flag, gkidx_allkpt, doublepack)

   integer, intent(in) :: ib, ik, ispin, npwk, fftsize(3), ngm
   integer, intent(in) :: master_gv(3,ngm)
   integer, intent(in) :: idxgk(npwk)
   logical, intent(in) :: shift_flag
   integer, intent(in) :: gkidx_allkpt(ngm)
   logical, intent(in) :: doublepack
   complex(dp), allocatable :: wfn(:,:)
   integer :: iunit = 30
   character (len=100) :: fplace, fname

   integer :: i, j
   ! total number of gvector to be written
   integer :: ngktot
   logical :: notfound
   real(dp) :: ZERO
   ZERO = dble(0)

   ! read wavefunction values
   allocate( wfn( npwk, 1 ) )
   if ( ispin .eq. 1 )  then
      call qexml_read_wfc ( IBNDS=ib, IBNDE=ib, IK=ik, WF=wfn, IERR=ierr)
   else
      call qexml_read_wfc ( IBNDS=ib, IBNDE=ib, IK=ik, ISPIN=ispin, WF=wfn, IERR=ierr)
   endif

   ! set up the directory where state files are stored
   if ( shift_flag .eqv. .false. ) then
      if ( ik .le. 10 ) then
         write( fplace, '( "./STATES_IN/Spin.", I1, "_Kpt.", I1, "_Bead.0_Temper.0/" )' ) ispin-1, ik-1
      elseif ( ik .le. 100 ) then
         write( fplace, '( "./STATES_IN/Spin.", I1, "_Kpt.", I2, "_Bead.0_Temper.0/" )' ) ispin-1, ik-1   
      else
         print*, 'Cannot locate directory because you have more than 100 k points'
         stop
      endif
   elseif( shift_flag .eqv. .true. ) then
      if ( ik .le. 10 ) then
         write( fplace, '( "./STATES_IN/Spin.", I1, "_Kpt.0", I1, "_Bead.0_Temper.0/" )' ) ispin-1, ik-1
      elseif ( ik .le. 100 ) then
         write( fplace, '( "./STATES_IN/Spin.", I1, "_Kpt.0", I2, "_Bead.0_Temper.0/" )' ) ispin-1, ik-1   
      else
         print*, 'Cannot locate directory because you have more than 100 k points'
         stop
      endif
   endif
   
   ! set file name state.out
   if ( ib .lt. 10) write( fname, '( "state", I1, ".out" )' )ib
   if ( ib .ge. 10 .and. ib .lt. 100) &
      write( fname, '( "state", I2, ".out" )' ) ib
   if ( ib .ge. 100 .and. ib .lt. 1000) &
      write( fname, '( "state", I3, ".out" )' ) ib
   if ( ib .ge. 1000 .and. ib .lt. 10000) &
      write( fname, '( "state", I4, ".out" )' ) ib
   if ( ib .ge. 10000 .and. ib .lt. 100000) &
      write( fname, '( "state", I2, ".out" )' ) ib
   
   !------------------------------------------
   ! writing starts here
   
   ! open state.out file
   open(iunit, file=trim(fplace)//trim(fname), form='formatted', status='replace')

   ! 1. write total number of g vectors and dense FFT grid size
   ngktot = sum(gkidx_allkpt)
   write(iunit,*) ngktot, fftsize(1:3)

   ! 2. write rest of them
   do i = 1, ngm
      notfound = .true.
      if ( gkidx_allkpt(i) == 0 ) then
         ! do nothing
      elseif ( gkidx_allkpt(i) == 1) then
         ! let's find if 
         do j = 1, npwk
            if( idxgk(j) == i ) then
               notfound = .false.
               write(iunit,*) realpart( wfn(j,1) ), imagpart( wfn(j,1) ),&
                    master_gv(1:3,idxgk(j))
               ! finish do j loop
               exit
            endif
         enddo
         if( notfound ) then
            write(iunit,*) ZERO, ZERO, master_gv(1:3,i)
         endif  
      else
         ! something wrong
         print*, '@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@'
         print*, 'gkdix_allkpt has illegal value for the point', i
         print*, '@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@'
         stop
      endif
   enddo

   close(iunit)
   deallocate( wfn )
   
end subroutine


!------------------------------------------------------------------------------
subroutine write_system_info( cell, ncoeff, xk, wk, &
                        nspin, nk, nb, fftsize,  outname )

   integer, parameter :: iunit = 20

   type(cell_info) :: cell
   integer, intent(in) :: ncoeff
   integer, intent(in) :: nspin, nk, nb
   integer, intent(in) :: fftsize(3)

   real(dp), intent(in) :: xk(3,nk), wk(nk)
   integer :: i, j, ik, ib, ispin
   character(len=256), intent(in) :: outname

   open(iunit,file=trim(outname),form='formatted', status='replace')
 
   ! cell lattice
   write(iunit,'(f12.6)') cell%alat
   
   ! lattice vectors
   write(iunit,'(3f12.6)') ( cell%a1(i), i=1,3 )
   write(iunit,'(3f12.6)') ( cell%a2(i), i=1,3 )
   write(iunit,'(3f12.6)') ( cell%a3(i), i=1,3 )


   ! reciprocal lattice vectors
   write(iunit,'(3f12.6)') ( cell%b1(i), i=1,3 )
   write(iunit,'(3f12.6)') ( cell%b2(i), i=1,3 )
   write(iunit,'(3f12.6)') ( cell%b3(i), i=1,3 )
   


   write(iunit,*) nspin, nk, nb
   
   ! number of planewaves at each k point -> will might need to be removed
   do ik = 1, nk
      write(iunit,*) ncoeff
   enddo

   ! k points and weights
   do ik = 1, nk
      write(iunit,'(4f12.6)') ( xk(i,ik), i=1,3 ), wk(ik)
   enddo
   
   ! dense fftsize
   write(iunit,*) fftsize(1), fftsize(2), fftsize(3)
   
   

   close(iunit)

end subroutine



!---------------------------------------------------------------------
! get information of simulation cell
subroutine get_cell_info( c )
   
   type(cell_info), intent(inout) :: c
   integer :: ierr
   
   call qexml_read_cell( ALAT=c%alat, &
                        A1=c%a1(:), A2=c%a2(:), A3=c%a3(:), &
                        B1=c%b1(:), B2=c%b2(:), B3=c%b3(:), &
                        A_UNITS=c%aunit, B_UNITS=c%bunit, IERR=ierr)
   
end subroutine




!-------------------------------------------------------------------
subroutine read_eig_occ( nspin, nk, nb, eigen, occp )

   integer, intent(in) :: nb, nk, nspin
   real(dp), intent(inout) :: eigen(nb,nk,nspin), occp(nb,nk,nspin)

   real(dp), allocatable :: eigtmp(:), occtmp(:)
   integer :: ik, ib, ispin

   ! work variables
   
   ! Comments:
   ! if the system does not have any spin index, including ISPIN in subroutine causes something nasty.
   ! have to check for spin polarized system if ISPIN works.
   ! same for qexml_read_wfc subroutine
 
   do ispin = 1, nspin
      do ik = 1, nk
         allocate( eigtmp(nb), occtmp(nb) )
         if ( nspin .eq. 1 ) then
            call qexml_read_bands( IK=ik, EIG=eigtmp, OCC=occtmp, IERR=ierr )
         else 
            call qexml_read_bands( IK=ik, ISPIN=ispin, EIG=eigtmp, OCC=occtmp, IERR=ierr)
         endif
         eigen(:,ik,ispin) = eigtmp(:)
         occp(:,ik,ispin) = occtmp(:)

         deallocate( eigtmp, occtmp )
      enddo
   enddo

end subroutine

!---------------------------------------------------------------------
! get number of plane waves at each k points
Subroutine read_pw_index( nk, ipwk )

   integer, intent(in) :: nk
   integer :: npwk, ik, ierr

   integer, intent(inout) :: ipwk(nk)

   do ik = 1, nk
      call qexml_read_gk( ik, NPWK=npwk, IERR=ierr )
      ipwk( ik ) = npwk
   enddo

end subroutine



!------------------------------------------------------------------------------
!  we need to save rho for GPP calculations
subroutine write_rho( nr, rho )
  
   integer, intent(in) :: nr(3)
   real(dp), intent(inout) :: rho(nr(1),nr(2),nr(3))
   integer,parameter :: iu = 31
   integer :: i,j,k
  
   call qexml_read_rho( RHO=rho, IERR=ierr )

   ! for FORTRAN CODE
  !open(iu,file='rho.dat',form='unformatted',status='unknown')

  !write(iu) nr(1:3)   
  !do k = 1, nr(3)
  !   write(iu) (rho(1:nr(1),j,k),j=1,nr(2))
  !enddo
  !
  !close(iu)


   ! for C++ CODE
   print*, "Rho is printed to a file. This rho is for C++ serial code"
   open(iu,file='rho.dat',form='formatted',status='unknown')
   write(iu,*) nr(1:3)
   do i = 1, nr(1)
      do j = 1, nr(2)
         write(iu,*) rho(i,j,1:nr(3))
      enddo
   enddo
   
 end subroutine write_rho





!-----------------------------------------------------------------------------
! find new ecutwfc. It will be slightly bigger than original ecutwfc
subroutine adjust_ecutwfc( cell, ecutwfc, nk, xk )

  ! cell information (i.e., lattice vectors and reciprocal lattice vectors)
  type(cell_info) :: cell
  ! wavefunction cutoff
  real(dp), intent(inout) :: ecutwfc
  ! number of k points
  integer, intent(in) :: nk
  ! k points in 2pi/alat inverse bohr unit
  real(dp), dimension(3,nk) ,intent(in) :: xk
  
  ! size of the k vector
  real(dp) :: trykmax, kmax
  ! define PI
  real(dp), parameter :: PI = 3.14159265359
  integer :: ik
  

  ! 1. find the largest k vector in inverse bohr unit
  ! initializie
  kmax = 0
  do ik = 1, nk
     trykmax = xk(1,ik)**2 + xk(2,ik)**2 + xk(3,ik)**2
     kmax = sqrt( kmax )
     if ( trykmax .ge. kmax ) then
        kmax = trykmax
     endif
  enddo

  print*, 'original ecutwfc: (Hartree)', ecutwfc
  ! new cutwfc
  ecutwfc = 0.5 * (sqrt(2*ecutwfc) + kmax*2*PI/cell%alat ) **2
  print*, 'new ecutwfc: (Hartree)', ecutwfc

end subroutine adjust_ecutwfc





!-----------------------------------------------------------------------------
! find new ecutwfc. It will be slightly bigger than original ecutwfc
  subroutine set_gkidx_allkpt( npwk, idxgk, ngm, master_gv, gkidx_allkpt )
     ! number of plane waves at this k point
     integer, intent(in) :: npwk
     ! connect to master_gv list
     integer, intent(in) :: idxgk(npwk)
     ! number of plane waves in the master g vector list
     integer, intent(in) :: ngm
     ! master g vector list
     integer, intent(in) :: master_gv(3,ngm)
     ! gkidx_allkpt. it keeps updating
     integer, intent(inout) :: gkidx_allkpt(ngm)

     integer :: i, this_gidx

     ! initially, gkidx_allkpt(:) is ZERO
     ! loop over npwk
     do i = 1, npwk
        this_gidx = idxgk(i)
        if( gkidx_allkpt( this_gidx ) .ne. 1 ) then
           gkidx_allkpt( this_gidx ) = 1
        endif
     enddo

     print*, 'Total g vector list up to this point:', sum( gkidx_allkpt ) 

  end subroutine set_gkidx_allkpt



!-----------------------------------------------------------------------------
! find universial g list
! kmax_cp is updated throughout spin/kpt loop
subroutine get_kmax( ngm, master_gv, npwk, idxgk, kmax_cp)
  ! number of g in the master gvec list
  integer, intent(in) :: ngm
  ! master gvector list
  integer, intent(in) :: master_gv(3,ngm)
  ! number of planewave at this k point
  integer, intent(in) :: npwk
  ! mapping to the master_gv
  integer, intent(in) :: idxgk(npwk)
  ! kmax_cp is used to count and set new g list
  integer, intent(inout) :: kmax_cp(3)

  integer :: ipw, i, tmpkmax(3)

  do ipw = 1, npwk
     do i = 1, 3
        tmpkmax(i) = master_gv(i,idxgk(ipw))
        ! if tmpkmax(i) is negative, make it positive
        if (tmpkmax(i) .lt. 0) then
           tmpkmax(i) = -1 * tmpkmax(i)
        endif
        ! compare if tmpkmax(i) is bigger than kmax_cp(i)
        if ( tmpkmax(i) .ge. kmax_cp(i) ) then
           ! if bigger, then we set new kmax_cp
           kmax_cp(i) = tmpkmax(i)
        endif
     enddo
  enddo
  
end subroutine get_kmax


 !-----------------------------------------------------------------------------
 ! Below subroutines come from OpenAtom program
 ! It is actually NOT called inside of the main program
 ! But I leave it here in case we may modify this converter later
 !-----------------------------------------------------------------------------
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BEGIN DO NOT USE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!------------------------------------------------------------------------------
!  count number of g vectors for doublepack (only half sphere is saved)
 subroutine countkvec3d_sm( kmax_cp, ng , ecut, cell )
   integer, intent(in) :: kmax_cp(3)
   integer, intent(inout) :: ng
   real(dp), intent(in) :: ecut
   type(cell_info) :: cell

   integer :: i1, i2, i3
   real(dp), parameter :: PI = 3.14159265359
   real(dp) :: xk, yk, zk, tryme
   real(dp) :: factor
   integer :: i, kamax, kbmax, kcmax, kbmin, kcmin, ka, kb, kc, icount
   real(dp) :: aka, akb, akc, g, gmin, gmax

   ! parameter set-up
   factor = 2*PI/cell%alat
   ! initialization
   gmin = 10000000
   gmax = 0
   icount = 0

! FIXME(?) Below doesn't work for non-cubic cell (like diamond structure), so I'll skip this for the moment
!  ! starts from b1 direction
!  i1 = kmax_cp(1)

!  ! expands to only b1 direction
!  do i = 1, i1
!     xk = i * cell%b1(1) * factor
!     yk = i * cell%b1(2) * factor
!     zk = i * cell%b1(3) * factor
!     tryme = (xk*xk + yk*yk + zk*zk)*0.5
!     print*, i, tryme
!     if (tryme .gt. ecut) exit
!  enddo

!  kamax = i-1  ! kamax =  maximum b1 
!  i1 = kamax

!  ! ka goes through 0 to the maximum b1  half sphere includes only ka >= 0
!  ! b1, b2, b3 are real reciprocal lattice vectors
!  do ka = 0, i1
!     aka = dble( ka )
!     kbmin = -kmax_cp(2)
!     if (ka==0) kbmin = 0
!     
!     ! 1. find the real minimum for b2 (kbmin) at this "ka(=aka)"
!     do i = kbmin, 0 ! kbmin = -kmax_cp[2]
!        akb = dble( i )
!        xk = ( aka * cell%b1(1) + akb * cell%b2(1) ) * factor
!        yk = ( aka * cell%b1(2) + akb * cell%b2(2) ) * factor
!        zk = ( aka * cell%b1(3) + akb * cell%b2(3) ) * factor
!        tryme = ( xk*xk + yk*yk + zk*zk )*0.5
!        if (tryme .le. ecut ) exit
!     enddo
!     kbmin = i ! this is the minimum for b2

!     ! 2. find the real maximum for b2 (kbmax) at this "ka(=aka)"
!     i2 = kmax_cp(2)
!     do i = 1, i2  ! 1 to kmax_cp[2]
!        akb = dble( i )
!        xk = ( aka * cell%b1(1) + akb * cell%b2(1) ) * factor
!        yk = ( aka * cell%b1(2) + akb * cell%b2(2) ) * factor
!        zk = ( aka * cell%b1(3) + akb * cell%b2(3) ) * factor
!        tryme = ( xk*xk + yk*yk + zk*zk ) * 0.5
!        if (tryme .gt. ecut ) exit
!     enddo
!     kbmax = i - 1

!     ! then loop over kbmin =< kb =< kbmax
!     i2 = kbmax
!     do kb = kbmin, i2
!        akb = dble( i )
!        kcmin = -kmax_cp(3)
!        if ( ka==0 .and. kb==0 ) kcmin = 1
!        ! 3. Find kcmin at this "ka(=aka)" and "kb(=akb)"
!        do i = kcmin, 0
!           akc = dble( i )
!           xk = ( aka * cell%b1(1) + akb * cell%b2(1) + akc * cell%b3(1) ) * factor
!           yk = ( aka * cell%b1(2) + akb * cell%b2(2) + akc * cell%b3(2) ) * factor
!           zk = ( aka * cell%b1(3) + akb * cell%b2(3) + akc * cell%b3(3) ) * factor
!           tryme = ( xk*xk + yk*yk + zk*zk ) * 0.5
!           if (tryme .le. ecut) exit
!        enddo
!        kcmin = i

!        ! 4. Fine kcmax at this "ka(=aka)" and "kb(=akb)"
!        i3 = kmax_cp(3)
!        do i = 1, i3
!           akc = dble( i )
!           xk = ( aka * cell%b1(1) + akb * cell%b2(1) + akc * cell%b3(1) ) * factor
!           yk = ( aka * cell%b1(2) + akb * cell%b2(2) + akc * cell%b3(2) ) * factor
!           zk = ( aka * cell%b1(3) + akb * cell%b2(3) + akc * cell%b3(3) ) * factor
!           tryme = ( xk*xk + yk*yk + zk*zk ) * 0.5
!           if (tryme .gt. ecut) exit
!        enddo
!        kcmax = i - 1

!        i3 = kcmax
!        ! then loop over kcmin =< kc =< kcmax
!        do kc = kcmin, i3
!           akc = dble( kc )
!           xk = ( aka * cell%b1(1) + akb * cell%b2(1) + akc * cell%b3(1) ) * factor
!           yk = ( aka * cell%b1(2) + akb * cell%b2(2) + akc * cell%b3(2) ) * factor
!           zk = ( aka * cell%b1(3) + akb * cell%b2(3) + akc * cell%b3(3) ) * factor
!           g = sqrt( xk*xk + yk*yk + zk*zk )
!           if( gmin .ge. g ) gmin = g
!           if( gmax .le. g ) gmax = g
!           icount = icount + 1
!        enddo! kc loop

!     enddo! kb loop
!  enddo! ka loop
!  
!  ! total number of g vectors except g = 0
!  ng = icount
!  ! add one to include g = 0
!  ng = ng + 1

   ! count again to check using kmax_cp
   icount = 0
   do ka = 0, kmax_cp(1)
      aka = dble( ka )
      kbmin = -kmax_cp(2)
      if( ka==0 ) kbmin = 0
      do kb = kbmin, kmax_cp(2)
         akb = dble( kb )
         kcmin = -kmax_cp(3)
         if( ka==0 .and. kb==0 ) kcmin = 1
         do kc = kcmin, kmax_cp(3)
            akc = dble( kc )
            xk = ( aka * cell%b1(1) + akb * cell%b2(1) + akc * cell%b3(1) ) * factor
            yk = ( aka * cell%b1(2) + akb * cell%b2(2) + akc * cell%b3(2) ) * factor
            zk = ( aka * cell%b1(3) + akb * cell%b2(3) + akc * cell%b3(3) ) * factor
            tryme = ( xk*xk + yk*yk + zk*zk ) * 0.5
            if( tryme .le. ecut ) icount = icount + 1
         enddo! kc
      enddo! kb
   enddo! ka

   ng = icount
   ng = ng + 1 ! to include g=0

   ! print results
   print*, 'gvector counts: half sphere without g=0 has', icount , 'vectors'
   print*, 'Ecutwfc:', ecut, 'kmax_cp:', kmax_cp(1), kmax_cp(2), kmax_cp(3)
   !print*, 'kamax, kbmax, kbmin, kcmax, kcmin:', kamax, kbmax, kbmin, kcmax, kcmin
   !print*, 'gmax:', gmax, 'gmin:', gmin
            
 end subroutine countkvec3d_sm

 subroutine setkvec3d_sm_simple( kmax_cp, ng, ecut, cell, glist )
   integer, intent(in) :: kmax_cp(3)
   integer, intent(in) :: ng
   real(dp), intent(in) :: ecut
   type(cell_info) :: cell
   integer, intent(inout) :: glist(3,ng)

   real(dp), parameter :: PI = 3.14159265359
   real(dp) :: xk, yk, zk, tryme
   real(dp) :: factor
   integer :: i, kamax, kbmax, kcmax, kbmin, kcmin, ka, kb, kc, icount
   real(dp) :: aka, akb, akc
   
   ! parameter set-up
   factor = 2*PI/cell%alat
   icount = 0
   do ka = 0, kmax_cp(1)
      aka = dble( ka )
      kbmin = -kmax_cp(2)
      if( ka==0 ) kbmin = 0
      do kb = kbmin, kmax_cp(2)
         akb = dble( kb )
         kcmin = -kmax_cp(3)
         if( ka==0 .and. kb==0 ) kcmin = 1
         do kc = kcmin, kmax_cp(3)
            akc = dble( kc )
            xk = ( aka * cell%b1(1) + akb * cell%b2(1) + akc * cell%b3(1) ) * factor
            yk = ( aka * cell%b1(2) + akb * cell%b2(2) + akc * cell%b3(2) ) * factor
            zk = ( aka * cell%b1(3) + akb * cell%b2(3) + akc * cell%b3(3) ) * factor
            tryme = ( xk*xk + yk*yk + zk*zk ) * 0.5
            if( tryme .le. ecut ) then
               icount = icount + 1
               glist(1,icount) = ka
               glist(2,icount) = kb
               glist(3,icount) = kc
            endif
            
         enddo! kc
      enddo! kb
   enddo! ka
   icount = icount + 1
   if ( ng .ne. icount ) then
      print*, '@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@'
      print*, ' number of g vectors do not match'
      print*, ' ng', ng , 'vs', icount
      print*, '@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@'
      stop
   endif
   ! g=0 is added
   glist(:,icount) = 0
      
 end subroutine setkvec3d_sm_simple
 



 subroutine countkvec3d_sm_kpt( kmax_cp, ng, ecut, cell )
   integer, intent(in) :: kmax_cp(3)
   integer, intent(inout) :: ng
   real(dp), intent(in) :: ecut
   type(cell_info) :: cell

   real(dp), parameter :: PI = 3.14159265359
   real(dp) :: factor
   integer :: icount, ka, kb, kc
   real(dp) :: xk, yk, zk, aka, akb, akc, g, gmin, gmax, tryme
   
   factor = 2 * PI / cell%alat
   ! initialization
   gmin = 10000000
   gmax = 0

   icount = 0
   do ka = -kmax_cp(1), kmax_cp(1)
      aka = dble( ka )
      do kb = -kmax_cp(2), kmax_cp(2)
         akb = dble( kb )
         do kc = -kmax_cp(3), kmax_cp(3)
            akc = dble( kc )
            xk = ( aka * cell%b1(1) + akb * cell%b2(1) + akc * cell%b3(1) ) * factor
            yk = ( aka * cell%b1(2) + akb * cell%b2(2) + akc * cell%b3(2) ) * factor
            zk = ( aka * cell%b1(3) + akb * cell%b2(3) + akc * cell%b3(3) ) * factor
            g = sqrt( xk*xk + yk*yk + zk*zk )
            if ( g .ge. gmax ) gmax = g
            if ( g .le. gmin ) gmin = g
            tryme = ( xk*xk + yk*yk + zk*zk ) * 0.5
            if( tryme .le. ecut ) icount = icount + 1
         enddo
      enddo
   enddo

   ! number of k points
   ng = icount
   
   ! printint results
   print*, 'gvector counts: full sphere has', icount ,'vectors'
   print*, 'ecutwfc:', ecut, 'kmax_cp:', kmax_cp(1), kmax_cp(2), kmax_cp(3)
   !print*, 'gmax:', gmax, 'gmin:', gmin

 end subroutine countkvec3d_sm_kpt

 subroutine setkvec3d_sm_kpt( kmax_cp, ng, ecut, cell, glist )
   integer, intent(in) :: kmax_cp(3)
   integer, intent(in) :: ng
   real(dp), intent(in) :: ecut
   type(cell_info) :: cell
   integer, intent(inout) :: glist(3,ng)

   real(dp), parameter :: PI = 3.14159265359
   real(dp) :: factor
   integer :: icount, ka, kb, kc
   real(dp) :: xk, yk, zk, aka, akb, akc, tryme
   
   factor = 2 * PI / cell%alat
   icount = 0

   do ka = -kmax_cp(1), kmax_cp(1)
      aka = dble( ka )
      do kb = -kmax_cp(2), kmax_cp(2)
         akb = dble( kb )
         do kc = -kmax_cp(3), kmax_cp(3)
            akc = dble( kc )
            xk = ( aka * cell%b1(1) + akb * cell%b2(1) + akc * cell%b3(1) ) * factor
            yk = ( aka * cell%b1(2) + akb * cell%b2(2) + akc * cell%b3(2) ) * factor
            zk = ( aka * cell%b1(3) + akb * cell%b2(3) + akc * cell%b3(3) ) * factor
            tryme = ( xk*xk + yk*yk + zk*zk ) * 0.5
            if( tryme .le. ecut ) then
               icount = icount + 1
               glist(1,icount) = ka
               glist(2,icount) = kb
               glist(3,icount) = kc
            endif
         enddo
      enddo
   enddo

   ! printint error
   if ( icount .ne. ng ) then
      print*, '@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@'
      print*, 'incorrect number of small kvectors in full sphere'
      print*, icount, 'vs', ng
      print*, '@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@'
   endif

 end subroutine setkvec3d_sm_kpt
 


end program










! we don't use this

!subroutine setkvec3d_sm( kmax_cp, ng , ecut, cell, glist )
!  integer, intent(in) :: kmax_cp(3)
!  integer, intent(in) :: ng
!  real(dp), intent(in) :: ecut
!  type(cell_info) :: cell
!  integer, intent(inout) :: glist(3,ng)

!  integer :: i1, i2, i3
!  real(dp), parameter :: PI = 3.14159265359
!  real(dp) :: xk, yk, zk, tryme
!  real(dp) :: factor
!  integer :: i, kamax, kbmax, kcmax, kbmin, kcmin, ka, kb, kc, icount
!  real(dp) :: aka, akb, akc, g, gmin, gmax

!  ! parameter set-up
!  factor = 2*PI/cell%alat
!  ! initialization
!  gmin = 10000000
!  gmax = 0
!  icount = 0

!  ! starts from b1 direction
!  i1 = kmax_cp(1)

!  ! expands to only b1 direction
!  do i = 1, i1
!     xk = i * cell%b1(1) * factor
!     yk = i * cell%b1(2) * factor
!     zk = i * cell%b1(3) * factor
!     tryme = (xk*xk + yk*yk + zk*zk)*0.5
!     if (tryme .gt. ecut) exit
!  enddo

!  kamax = i-1  ! kamax =  maximum b1 - 1
!  i1 = kamax

!  ! ka goes through 0 to the maximum b1  half sphere includes only ka >= 0
!  ! b1, b2, b3 are real reciprocal lattice vectors
!  do ka = 0, i1
!     aka = dble( ka )
!     kbmin = -kmax_cp(2)
!     if (ka==0) kbmin = 0
!     
!     ! 1. find the real minimum for b2 (kbmin) at this "ka(=aka)"
!     do i = kbmin, 0 ! kbmin = -kmax_cp[2]
!        akb = dble( i )
!        xk = ( aka * cell%b1(1) + akb * cell%b2(1) ) * factor
!        yk = ( aka * cell%b1(2) + akb * cell%b2(2) ) * factor
!        zk = ( aka * cell%b1(3) + akb * cell%b2(3) ) * factor
!        tryme = ( xk*xk + yk*yk + zk*zk )*0.5
!        if (tryme .le. ecut ) exit
!     enddo
!     kbmin = i ! this is the minimum for b2

!     ! 2. find the real maximum for b2 (kbmax) at this "ka(=aka)"
!     i2 = kmax_cp(2)
!     do i = 1, i2  ! 1 to kmax_cp[2]
!        akb = dble( i )
!        xk = ( aka * cell%b1(1) + akb * cell%b2(1) ) * factor
!        yk = ( aka * cell%b1(2) + akb * cell%b2(2) ) * factor
!        zk = ( aka * cell%b1(3) + akb * cell%b2(3) ) * factor
!        tryme = ( xk*xk + yk*yk + zk*zk ) * 0.5
!        if (tryme .gt. ecut ) exit
!     enddo
!     kbmax = i - 1

!     ! then loop over kbmin =< kb =< kbmax
!     i2 = kbmax
!     do kb = kbmin, i2
!        akb = dble( i )
!        kcmin = -kmax_cp(3)
!        if ( ka==0 .and. kb==0 ) kcmin = 1
!        ! 3. Find kcmin at this "ka(=aka)" and "kb(=akb)"
!        do i = kcmin, 0
!           akc = dble( i )
!           xk = ( aka * cell%b1(1) + akb * cell%b2(1) + akc * cell%b3(1) ) * factor
!           yk = ( aka * cell%b1(2) + akb * cell%b2(2) + akc * cell%b3(2) ) * factor
!           zk = ( aka * cell%b1(3) + akb * cell%b2(3) + akc * cell%b3(3) ) * factor
!           tryme = ( xk*xk + yk*yk + zk*zk ) * 0.5
!           if (tryme .le. ecut) exit
!        enddo
!        kcmin = i

!        ! 4. Fine kcmax at this "ka(=aka)" and "kb(=akb)"
!        i3 = kmax_cp(3)
!        do i = 1, i3
!           akc = dble( i )
!           xk = ( aka * cell%b1(1) + akb * cell%b2(1) + akc * cell%b3(1) ) * factor
!           yk = ( aka * cell%b1(2) + akb * cell%b2(2) + akc * cell%b3(2) ) * factor
!           zk = ( aka * cell%b1(3) + akb * cell%b2(3) + akc * cell%b3(3) ) * factor
!           tryme = ( xk*xk + yk*yk + zk*zk ) * 0.5
!           if (tryme .gt. ecut) exit
!        enddo
!        kcmax = i - 1

!        i3 = kcmax
!        ! then loop over kcmin =< kc =< kcmax
!        do kc = kcmin, i3
!           akc = dble( kc )
!           xk = ( aka * cell%b1(1) + akb * cell%b2(1) + akc * cell%b3(1) ) * factor
!           yk = ( aka * cell%b1(2) + akb * cell%b2(2) + akc * cell%b3(2) ) * factor
!           zk = ( aka * cell%b1(3) + akb * cell%b2(3) + akc * cell%b3(3) ) * factor
!           g = sqrt( xk*xk + yk*yk + zk*zk )
!           if( gmin .ge. g ) gmin = g
!           if( gmax .le. g ) gmax = g
!           icount = icount + 1
!           glist(1,icount) = ka
!           glist(2,icount) = kb
!           glist(3,icount) = kc
!        enddo! kc loop

!     enddo! kb loop
!  enddo! ka loop

!  if ( ng-1 .ne. icount ) then
!     print*, "@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@"
!     print*, "Mismatch number of small kvectors:"
!     print*, icount ,'vs', ng
!     print*, "@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@"
!  endif
!  
!  ! total number of g vectors
!  icount = icount + 1
!  glist(:,icount) = 0

!end subroutine setkvec3d_sm
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END DO NOT USE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
