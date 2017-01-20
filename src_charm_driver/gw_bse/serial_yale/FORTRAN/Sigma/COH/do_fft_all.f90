!
! FFT for all bands
!

subroutine FFT_all_bands( psi_G , psi_R, fftsize, sys )

   use constant
   use electronic_structure
   implicit none
   
   type( wfstruc ), intent(inout) :: psi_G, psi_R
   !type( wfstruc ), allocatable :: psi_N
   integer, intent(in) :: fftsize(3)
   type( sysinfo ), intent(in) :: sys
   
   ! work variables
   integer :: nk ! number of k points
   integer :: nb ! number of bands
   integer :: ns ! number of spin
   integer :: ndata ! number of data (=size of the wave function (1-D) array)
   real(dp) :: scale
   
   integer :: ik

   nk = psi_G%nkpt
   nb = psi_G%nband
   ns = psi_G%nspin
   ndata = fftsize(1)*fftsize(2)*fftsize(3)
   
   scale = 1.d0 / sqrt( sys%vol )


   ! assign nk, nb, ns to psi_R structure
   psi_R%nkpt = nk
   psi_R%nband = nb
   psi_R%nspin = ns
   
   
   ! allocate memory to save wavefunction
   allocate( psi_R%wk(nk) )
   ! psi_R%wk(nk)= psi_N%wk(nk)
   do ik = 1, nk
     
      call allocate_psi_R( psi_G%wk(ik), psi_R%wk(ik), fftsize, nb )
      
      call psi_G_to_R( psi_G%wk(ik), psi_R%wk(ik), fftsize )

      ! divide by sqrt( cell_vol )
      psi_R%wk(ik)%cg = psi_R%wk(ik)%cg * scale
      
      deallocate( psi_G%wk(ik)%cg )
      deallocate( psi_G%wk(ik)%eig, psi_G%wk(ik)%occ )
    
   enddo
      
end subroutine


!-----------------------------------------------
! in this subroutine, the number of FFT performed = Nbnd.
! Input is the entire wave function ( 1st to the last band) for a specific k point
subroutine psi_G_to_R( wfn_G, wfn_R, FFTsize )
   
   use constant
   use electronic_structure
   implicit none
   
   type( wfn ), intent(inout) :: wfn_G, wfn_R
   integer, intent(in) :: FFTsize(3)
   
   ! work variables
   integer :: ib, nb ! band index and number of bands
   integer :: ng ! number of g vectors
   integer :: ig, gidx(3)
   integer, allocatable :: idx(:,:)
   complex(dp), allocatable :: fftbox(:,:,:)
   complex(dp), allocatable :: a_g(:)
   complex(dp), allocatable :: a_r(:)
   integer, parameter :: one = 1  ! FFT sign BACKWARD
   integer :: istart, iend
   integer :: ndata
   

! get this number from somewhere else... not here
   nb  = size(wfn_G%eig(:),1)

!-----------
   ndata = FFTsize(1)*FFTsize(2)*FFTsize(3)

   ng = wfn_G%ng

   allocate( idx(3,ng) )
   

   ! idx(3,ng) will be different at different k points. 
   ! idx(3,ng) is valid for all bands at the current k point.
   do ig = 1, ng
      
      gidx(1:3) = wfn_G%gvec(:,ig)
      
      call gidx_to_fft_index( gidx, FFTsize, idx(:,ig) ) ! this subroutine is correct
      
   enddo

   allocate ( a_g( ng ) )
   allocate( a_r( ndata ) )
   allocate( fftbox( 1:FFTsize(1), 1:FFTsize(2), 1:FFTsize(3) ) )


   do ib = 1, nb

      a_g = wfn_G%cg( 1:ng, ib )
      
      fftbox = 0.d0  
   
      call put_into_fftbox( ng, a_g, idx, FFTsize, fftbox ) 


      call do_fft( fftbox, FFTsize, one )

      a_r = 0
      call box_to_array( FFTsize, fftbox, a_r )
      
      wfn_R%cg( 1:ndata, ib ) = a_r( 1:ndata ) 
      
   enddo
   
   deallocate( idx, a_g, a_r )


end subroutine



!------------------------------------------------------
subroutine allocate_psi_R( wfn_G, wfn_R, fftsize, nband )

   use constant
   use electronic_structure
   implicit none
   
   type( wfn ), intent(inout) :: wfn_G, wfn_R
   integer, intent(in) :: fftsize(3)
   integer, intent(in) :: nband
   
   integer :: ndata
   
   ndata = fftsize(1)*fftsize(2)*fftsize(3)
   
   ! allocate memory for wavefunction data, eig, and occ
   allocate( wfn_R%cg(ndata,nband) )
   allocate( wfn_R%eig(nband) )
   allocate( wfn_R%occ(nband) )
   
   ! assign  values
   wfn_R%eig(1:nband) = wfn_G%eig(1:nband)
   wfn_R%occ(1:nband) = wfn_G%occ(1:nband)
   wfn_R%kidx = wfn_G%kidx

end subroutine
   
