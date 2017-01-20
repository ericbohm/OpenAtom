! 
!   This subroutine reads wavefunction data  
!      and store to the type 'wfstruc'
!
!
!   Sep. 2014   minjung.kim@yale.edu


! WFN structure: (obtained from save_QEwfn.f90)
! alat, avec(3:3) (in atomic unit)
! bvec(3:3) ( in reciprocal lattie vectors )
! ngm, nspin, nk, nb, npwk_all
! gvec(3,ngm)
! npwk(nk)
! igk_all(sum(npwk))
! (empty line)
! eigenvalue(nb,nk,nspin)
! occupancy(nb,nk,nspin)
! empty line
! wfn( sum(npwk)*nb*nspin )   (order: is-ik-ib-ig)


subroutine read_wfn( fname, sys, psi, k )

   use constant
   use usrinput
   use electronic_structure
    
   implicit none
   character(len=100) :: fname
   type(sysinfo), intent(inout) :: sys
   type(wfstruc), intent(inout) :: psi
   type(kptinfo), intent(inout) :: k

   ! work variables
   integer, parameter :: iunit=20   ! input file unit
   integer :: ng, ns, nk, nb, ngk
   integer :: ii,jj,ib,ik,is,iks
   integer :: istart, iend, counter
   integer :: npwk_all    ! number 
   integer,allocatable :: npwk(:), igk_all(:)
   real(dp), allocatable :: eig(:,:,:), occ(:,:,:)
   complex(dp), allocatable :: wfntmp(:)
   
   integer :: ierr
   real(dp) :: avectmp(3,3), vol
   real(dp) :: scale


   
!------------------ READING BLOCK ---------------------!
   open(unit=iunit,file=trim(fname),form='unformatted',status='unknown',iostat=ierr)
   if ( ierr .gt. 0 ) then  
      call print_error_exit( 'Error opening wavefunction data file (wfn.dat)' )
   endif

   ! Start reading 

   ! cell information
   read(iunit) sys%alat, ((sys%avec(ii,jj),ii=1,3),jj=1,3) ! in cartesian unit
   read(iunit) ((sys%bvec(ii,jj),ii=1,3),jj=1,3) ! unit is 2pi/alat

   ! calculate volume of the system
   avectmp = transpose( sys%avec )
   call ddet( avectmp, sys%vol )
   
   print*, 'volume of the simulation cell:', sys%vol

   ! number of total gvec, # spin, # kpt, # bands 
   read(iunit) psi%ng, psi%nspin, psi%nkpt, psi%nband
   
   ng=psi%ng; ns=psi%nspin; nk=psi%nkpt; nb=psi%nband
   
   allocate( psi%gvec(3,ng) )
   allocate( npwk(nk) )
   
   ! read g vectors (ng depends on Rho energy cutoff)
   read(iunit) ( ( psi%gvec(ii,jj), ii=1,3 ), jj=1,ng )
   read(iunit) ( npwk(ii), ii=1,nk )

   npwk_all = sum( npwk )
   allocate( igk_all(npwk_all) ) 

   read(iunit) igk_all
   
   allocate( k%vec(3,nk), k%wt(nk) )
   read(iunit) ( ( k%vec(ii,ik), ii=1,3 ), ik=1,nk )
   read(iunit) k%wt
   k%nk = nk
      
   allocate( eig(nb,nk,ns), occ(nb,nk,ns) )
   
   read(iunit) (((eig(ib,ik,is),ib=1,nb),ik=1,nk),is=1,ns)
   read(iunit) (((occ(ib,ik,is),ib=1,nb),ik=1,nk),is=1,ns)
   read(iunit) ! empty line

   allocate( wfntmp(npwk_all*nb*ns) )
   read(iunit) wfntmp
   
   close(iunit)
!------------- READING BLOCK ENDS ---------------------!   

   ! assign values to psi data structure
   
   ! psi%wk(:) will save wavefunctions at each k points
   ! order: (1st spin) 1,2,3,....,nk, (2nd spin) 1,2,3,...,nk, ...etc

   allocate( psi%wk(nk*ns) )  

   do is = 1, ns
      do ik= 1, nk
         iks = nk*(is-1)+ik
         ! k index
         psi%wk(iks)%kidx = ik
         ! number of g vectors at ik
         psi%wk(iks)%ng = npwk(ik)
         ngk = npwk(ik)
    
         ! allocate arrays for eigenvalues, occupancies, and wfn data
         allocate( psi%wk(iks)%eig(nb) )
         allocate( psi%wk(iks)%occ(nb) )
         allocate( psi%wk(iks)%cg( psi%wk(iks)%ng , nb ) )
         allocate( psi%wk(iks)%gvec(3,psi%wk(iks)%ng) )
         
         psi%wk(iks)%eig(1:nb) = eig(1:nb,ik,is)
         psi%wk(iks)%occ(1:nb) = occ(1:nb,ik,is)
        
         ! let's assign g vectors
        
         istart = sum( npwk(1:ik) ) - npwk(ik) + 1
         iend = sum( npwk(1:ik) )

         counter = 0
         do ii = istart, iend
            counter = counter + 1
            psi%wk(iks)%gvec(1:3,counter) = psi%gvec( 1:3,igk_all(ii) )
         enddo
               
         ! wfn coefficient Band index!!! npwk_all*nb = total in one spin channel   npwk(ik)*nb
         do ib = 1, nb
            istart = ( (is-1)*npwk_all + sum( npwk(1:ik-1) ) )*nb + npwk(ik)*(ib-1) + 1
            iend = ( (is-1)*npwk_all + sum( npwk(1:ik-1) ) )*nb + npwk(ik)*ib
            psi%wk(iks)%cg( 1:npwk(ik), ib ) = wfntmp(istart:iend)
         enddo

      enddo
   enddo

   deallocate( npwk, igk_all, eig, occ, wfntmp )
   
end subroutine
