!
! Sep. 2014   minjung.kim@yale.edu
!
! Last modified: Feb. 2015 by MK


Program save_wfn

   use qexml_module

   implicit none

   integer, parameter :: iunit = 10
   integer, parameter :: dp = kind(1.0d0)

   type cell_info
      real(dp) :: alat
      real(dp) :: avec(3,3)
      real(dp) :: bvec(3,3)
      character(len=256) :: aunit, bunit
   end type
   
   type(cell_info) :: cell
  
   integer :: ngk_tot, nb
   integer, allocatable :: ipwk(:), igk_all(:)

   integer :: i
   integer :: ierr
   
   character(256) :: work_dir, prefix, dirname, filename, wfn_name
   integer :: nk ! it should be decided at first
   integer :: nspin

   integer :: ngm
   integer, allocatable :: master_gv(:,:)
   real(dp), allocatable :: eig(:,:,:), occ(:,:,:)
   complex(dp), allocatable :: wfn(:)
   real(dp), allocatable :: xk(:,:), wk(:) 

   ! GPP calculation related variables
   integer :: nr(3)
   real(dp), allocatable :: rho(:,:,:)
   logical :: gpp_flag
   
   
   integer :: stdin = 5

   NAMELIST /INPUT/ prefix, work_dir, wfn_name, gpp_flag

!---------------------------------------------------------------------
! starts here
!---------------------------------------------------------------------

   ! default setting
   work_dir = './'
   gpp_flag = .false.
print*, 'read in file'   
   read( stdin, INPUT, iostat=ierr )
   
      
   ! initialize QEXML library
   
   dirname = trim(work_dir) // '/' // trim(prefix) // '.save/'
   call qexml_init( iunit, Dir=dirname )
print*, 'qexml_init'
   filename = trim(dirname)//"data-file.xml"
   
   call qexml_openfile( filename, "read", IERR=ierr )
print*, 'qexml_openfile'

   ! read cell information
   call get_cell_info( cell )
print*, 'get_cell_info'

   ! read number of k points = nk
   call qexml_read_bands_info( NUM_K_POINTS=nk, IERR=ierr)
print*, 'qexml_read_bands_info'
   ! read BZ: k points and their weight factor
   allocate( xk(3,nk), wk(nk) )
   call qexml_read_bz( XK=xk, WK=wk, IERR=ierr )
print*, 'qexml_read_bz'
   ! calculate total number of g vectors (from first to the last k points)
   ! ipwk: #pws in each k point
   allocate( ipwk( nk ) )
   
   
   
   ! initialization
   ipwk = 0
   call read_pw_index( nk, ipwk )






print*, 'read_pw_index'
   ! total number of g index (sum (npw_k) )
   ngk_tot = sum( ipwk(1:nk) )
print*, 'ipwk', ipwk
print*, 'nk, ngk_tot:', nk, ngk_tot
   allocate( igk_all( ngk_tot ) )
   ! get gk index at each k points (from 1st to nk k points)
   call get_igk_all( nk, ngk_tot, igk_all )
print*, 'get_igk_all'
   ! get master g vector
   call qexml_read_planewaves( NGM=ngm, IERR=ierr )
print*, 'qexml_read_planewaves'
   allocate( master_gv(3,ngm) )
   call qexml_read_planewaves( IGV=master_gv, IERR=ierr )
print*, 'qexml_read_planewaves'

   ! Read number of bands and number of spins
   call qexml_read_bands_info( NBND=nb, NSPIN=nspin, IERR=ierr)
print*, 'qexml_read_bands_info'

   ! read eigenvalues and occupancies
   allocate( eig(nb, nk, nspin), occ(nb, nk, nspin) )
   call read_eig_occ( nspin, nk, nb, eig, occ )
print*, 'read_eig_occ'

   ! Read wave function
   allocate( wfn( nspin*ngk_tot*nb ) )
   
   call read_wfn( nb, nk, nspin, ipwk, ngk_tot, wfn )
print*, 'read_wfn'

   ! wfn order: is-ik-ib-ig
   
   call qexml_closefile( "read", IERR=ierr )
print*, 'qexml_closefile'

! Write wfn data in one file

   call write_outfile( cell, ngm, master_gv, ipwk, ngk_tot, igk_all, xk, wk, &
                       nspin, nk, nb, eig, occ, wfn(:), wfn_name )
print*, 'write_outfile'


! If RHO flag is on, create rho.dat file for GPP calculations
   if ( gpp_flag ) then
      call qexml_read_rho( NR1=nr(1), NR2=nr(2), NR3=nr(3), IERR=ierr )
      allocate( rho(nr(1),nr(2),nr(3)) )
      call write_rho( nr, rho )
      deallocate(rho )
   endif
print*, 'finished'

deallocate(xk, wk, ipwk, igk_all, master_gv, eig, occ, wfn)

contains


!************ subroutines ********************************************


Subroutine write_outfile( cell, ngm, master_gv, npwk, sum_npwk, igk_all, xk, wk, &
                        nspin, nk, nb, eig, occ, wfn, outname)

   integer, parameter :: iunit = 20

   type(cell_info) :: cell
   integer, intent(in) :: ngm, nspin, nk, nb, sum_npwk
   integer, intent(in) :: master_gv(3,ngm)
   integer, intent(in) :: npwk(nk), igk_all(sum_npwk)

   real(dp), intent(in) :: eig(nb,nk,nspin), occ(nb,nk,nspin)
   complex(dp), intent(in) :: wfn( sum_npwk * nb * nspin )
   real(dp), intent(in) :: xk(3,nk), wk(nk)
   integer :: i, j, ik, ib, ispin
   character(len=256), intent(in) :: outname

   open(iunit, file=trim(outname), form='unformatted',status='unknown')
   
   write(iunit) cell%alat, (( cell%avec(i,j), i=1,3 ), j=1,3 )
   write(iunit) (( cell%bvec(i,j), i=1,3 ), j=1,3 )
   write(iunit) ngm, nspin, nk, nb, sum_npwk
   write(iunit) (( master_gv(i,j),i=1,3), j=1,ngm)
   write(iunit) npwk
   write(iunit) igk_all
   write(iunit) ((xk(i,ik),i=1,3), ik=1,nk)
   write(iunit) wk(1:nk)
   write(iunit) ((( eig(ib,ik,ispin),ib=1,nb), ik=1,nk), ispin=1,nspin)
   write(iunit) ((( occ(ib,ik,ispin),ib=1,nb), ik=1,nk), ispin=1,nspin)
   write(iunit)
   write(iunit) wfn

end subroutine



!---------------------------------------------------------------------
! get information of simulation cell
Subroutine get_cell_info( c )
   
   type(cell_info), intent(inout) :: c
   integer :: ierr
   
   call qexml_read_cell( ALAT=c%alat, &
                        A1=c%avec(:,1), A2=c%avec(:,2), A3=c%avec(:,3), &
                        B1=c%bvec(:,1), B2=c%bvec(:,2), B3=c%bvec(:,3), &
                        A_UNITS=c%aunit, B_UNITS=c%bunit, IERR=ierr)
   
end subroutine



!---------------------------------------------------------------------
! get number of plane waves at each k points
Subroutine read_pw_index( nk, ipwk )

   integer, intent(in) :: nk
   integer :: npwk, ik, ierr

   integer, intent(inout) :: ipwk(nk)

   do ik = 1, nk
      call qexml_read_gk( ik, NPWK=npwk, IERR=ierr )
      print*, 'npwk at :', ik, npwk
      ipwk( ik ) = npwk
   enddo

end subroutine



!---------------------------------------------------------------------
! get k specific gvector index and save it in 1-D array

Subroutine get_igk_all( nk, ngk_tot, igk_all )

   integer :: nk, ngk_tot
   integer, allocatable :: igk(:) ! g vector indexes at each k point
   integer, intent(inout) :: igk_all(ngk_tot) !  g indexes in all k points (dim=sum(ipwk))
   integer :: npwk, ierr
   
   ! work variables
   integer :: ik, counter, inpw
   
   ! store gkvec index in 1D array : igv_all
   counter = 0
   do ik = 1, nk
      ! call qexml
      call qexml_read_gk( ik, NPWK=npwk, IERR=ierr )
      allocate( igk( npwk ) )
      call qexml_read_gk( ik, index=igk, IERR=ierr )
      do inpw = 1, npwk
         counter = counter + 1
         igk_all( counter ) = igk( inpw )
      enddo
      deallocate( igk )
   enddo

   if ( counter .ne. ngk_tot ) then
     print*, 'Something is missing... Check get_igk_all subroutine!'
   endif 

end Subroutine 


!-------------------------------------------------------------------
Subroutine read_eig_occ( nspin, nk, nb, eigen, occp )

   integer, intent(in) :: nb, nk, nspin
   real(dp), intent(inout) :: eigen(nb,nk,nspin), occp(nb,nk,nspin)

   real(dp), allocatable :: eigtmp(:,:), occtmp(:,:)
   integer :: ik, ib, ispin

   ! work variables
   
   ! Comments:
   ! if the system does not have any spin index, including ISPIN in subroutine causes something nasty.
   ! have to check for spin polarized system if ISPIN works.
   ! same for qexml_read_wfc subroutine
 
 
   ! modification for qe-6.0 (Later than qe-5.3 is affected by this)
   logical :: lsda, lkpoint_dir
   character(len=256) :: filename
   integer :: nkstot, iks
   
   nkstot = nspin * nk
   if( nspin==1 ) then
      lsda = .false.
   elseif( nspin==2) then
      lsda = .true.
   endif
   
   lkpoint_dir = .true.
   filename = 'data-file.xml'
   allocate( eigtmp(nb, nkstot), occtmp(nb,nkstot) )
 
   call qexml_read_bands_pw( NUM_K_POINTS=nk, NBND=nb, NKSTOT=nkstot, LSDA=lsda,&
    LKPOINT_DIR=lkpoint_dir, FILENAME=filename, ET=eigtmp, WG=occtmp, IERR=ierr) 
   iks = 0
   do ispin = 1, nspin
      do ik = 1, nk
         iks = iks + 1
         eigen(:,ik,ispin) = eigtmp(:,iks)
         occp(:,ik,ispin) = occtmp(:,iks)
      enddo
   enddo
   
   deallocate (eigtmp, occtmp)

end subroutine


!------------------------------------------------------------------------------
Subroutine read_wfn( nb, nk, nspin, ipwk, ngk_tot, wfn )

   integer, intent(in) :: nk, nb, nspin, ipwk(nk)
   integer, intent(in) :: ngk_tot
   complex(dp), intent(inout) :: wfn(nspin*ngk_tot*nb)
   complex(dp), allocatable :: wftmp(:,:)
   integer :: ik, ib, npw, nn, ini, fin, ispin

   do ispin = 1, nspin   
      do ik = 1, nk
         npw = ipwk(ik)
         nn = sum( ipwk(1:ik) ) - npw
         do ib = 1, nb
            allocate( wftmp( npw, 1 ) )
            if ( ispin .eq. 1 )  then
               call qexml_read_wfc ( IBNDS=ib, IBNDE=ib, IK=ik, WF=wftmp, IERR=ierr)
            else
               call qexml_read_wfc ( IBNDS=ib, IBNDE=ib, IK=ik, ISPIN=ispin, WF=wftmp, IERR=ierr)
            endif
            ini = (ispin-1)*ngk_tot*nb + nn*nb + npw*(ib-1) + 1
            fin = (ispin-1)*ngk_tot*nb + nn*nb + npw*ib
            wfn( ini : fin ) = wftmp(1:npw,1)
            deallocate( wftmp )
         enddo
      enddo
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

   open(iu,file='rho.dat',form='unformatted',status='unknown')

   write(iu) nr(1:3)   
   do k = 1, nr(3)
      write(iu) (rho(1:nr(1),j,k),j=1,nr(2))
   enddo
   
   close(iu)
   
end subroutine

end program
