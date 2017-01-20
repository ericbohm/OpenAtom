!  Oct. 2014    minjung.kim@yale.edu
!  This subroutine do FFTs for each bands and update Pmtrx
!
!  Spin is not yet implemented

subroutine calc_Prrp( psi_v, psi_c, Pmtrx, FFTsize, sys, Uvec )
   
   use constant
   use electronic_structure
   use gw_structure
   
   implicit none
   type(wfn), intent(in) :: psi_v, psi_c
   type(rank2_mtrx), intent(inout) :: Pmtrx
   integer, intent(in) :: FFTsize(3)
   type(sysinfo), intent(in) :: sys
   integer, intent(in) :: Uvec(3)

   ! work variables
   complex(dp), dimension(:), allocatable :: a_v, a_c
   integer :: ndata
   integer :: nv, nc, nb
   integer :: iv, ic
   integer :: utest
   
   ! Change it ------------------------------
   nv = int ( sum( psi_v%occ ) )
   !nb = int ( sizeof( psi_c%eig(:) ) )/ 8 
   nb = 52
   nc = nb - nv
   !-----------------------------------------

   ndata = fftsize(1)*fftsize(2)*fftsize(3)
   utest = abs( Uvec(1) ) + abs( Uvec(2) ) + abs( Uvec(3) )


   print*, 'number of bands, nv, nc:', nb, nv, nc
   allocate( a_v(ndata), a_c(ndata) )

   ! no spin

   ! loop for valence bands
   do iv = 1, nv

      a_v(1:ndata) = psi_v%cg(1:ndata,iv)
      
      ! change wavefunction for U-process
      if (utest .eq. 0 ) then 
         continue
      else
         call modify_Uproc_wfn ( ndata, a_v, Uvec, sys, FFTsize )
      endif

      ! loop for conduction bands
      do ic = 1, nc

         a_c(1:ndata) = psi_c%cg(1:ndata,nv+ic)
         
         call update_P( a_v, a_c, ndata, psi_v%eig(iv), psi_c%eig(nb-nc+ic), Pmtrx%C , sys%nk )    

      enddo

   enddo
   
   
   contains 
   subroutine update_P( a_v, a_c, ndata, Ev, Ec, P, nk )

      complex(dp), dimension(ndata), intent(inout) :: a_v, a_c 
      integer, intent(in) :: ndata
      real(dp), intent(in) :: Ev, Ec
      complex(dp), dimension(ndata,ndata), intent(inout) :: P
      integer, intent(in) :: nk
   
      ! work variables
      complex(dp), dimension(ndata) :: f
      integer :: i, j
      real(dp) :: fact


      f(1:ndata) = a_v(1:ndata) * conjg( a_c(1:ndata) )

      fact = 4.d0/(Ev-Ec)  !  for nspin = 1 
      ! If you want to compare it with BGW, factor should be 2/Nk for nspin=1 because 
      ! their energy unit isf in Rydberg

      do j = 1, ndata
         do i = 1, ndata

            p(i,j) = fact * f(i) * conjg( f(j) ) + p(i,j)
            
         enddo
      enddo
   end subroutine update_P
   
   
   subroutine modify_Uproc_wfn( ndata, a_c, Uvec, sys, FFTsize )

      integer, intent(in) :: ndata
      complex(dp), dimension(ndata), intent(inout) :: a_c
      integer, dimension(3), intent(in) :: Uvec
      type(sysinfo), intent(in) :: sys
      integer, dimension(3), intent(in) :: FFTsize
      
      ! work variables
      complex(dp), allocatable :: fftbox(:,:,:), fact(:,:,:)
      integer :: i, j, k, ii
      integer, allocatable :: idx(:,:)
      real(dp) :: a(3,3), b(3,3)
      real(dp) :: rijk(3), G0(3), phase
      
      allocate( idx(3,ndata) )
      allocate( fftbox(fftsize(1),fftsize(2),fftsize(3)) )
      allocate( fact(fftsize(1),fftsize(2),fftsize(3)) )
      
      ! set fftbox index
      call set_3Dbox_index( ndata, FFTsize, idx )
            
      call put_into_fftbox( ndata, a_c, idx, FFTsize, fftbox )
      
      
      a(1:3,1:3) = sys%avec(1:3,1:3)
      b(1:3,1:3) = sys%bvec(1:3,1:3)*2.d0*pi/sys%alat
   
      ! calculate factor to be multiplied to the fftbox
      do k = 1, FFTsize(3)
         do j = 1, FFTsize(2)
            do i = 1, FFTsize(1) 
               do ii = 1, 3
                  rijk(ii) = a(ii,1)*(i-1)/FFTsize(1) + a(ii,2)*(j-1)/FFTsize(2) + &
                             a(ii,3)*(k-1)/FFTsize(3)
                  G0(ii) = b(ii,1)*Uvec(1) + b(ii,2)*Uvec(2) + b(ii,3)*Uvec(3)
               enddo
               G0 = -1.d0*G0
               phase = dot_product( G0, rijk )
               fact(i,j,k) = cmplx( cos(phase), sin(phase) )
            enddo
         enddo
      enddo
      fftbox(1:fftsize(1),1:fftsize(2),1:fftsize(3)) =    &
            fact(1:fftsize(1),1:fftsize(2),1:fftsize(3))*    &
            fftbox(1:fftsize(1),1:fftsize(2),1:fftsize(3)) 
         
      
      call box_to_array( FFTsize, fftbox, a_c )
   
      deallocate( idx, fftbox, fact )
   end subroutine modify_Uproc_wfn
   
   
end subroutine! - End of Calc_Prrp subroutine



   
!###########################################################################
!
!                 FFT : P(r,r') -> P(G,G')
!                 #FFT = 2*ndata
!
subroutine P_r_to_g( iq, P, FFTsize )

   use constant
   use gw_structure
   implicit none
   integer, intent(in) :: iq
   integer, intent(in) :: FFTsize(3)
   type(rank2_mtrx), intent(inout) :: P
   
   ! work variables
   complex(dp), allocatable :: a_r(:), fftbox(:,:,:)
   integer :: icol, irow, colrow
   integer :: ndata
   integer, allocatable :: idx(:,:)
   
   ndata = FFTsize(1)*FFTsize(2)*FFTsize(3)
   allocate( idx(3, ndata) )
   allocate( a_r(ndata), fftbox(FFTsize(1),FFTsize(2),FFTsize(3)) )
   
   ! set index
   call set_3Dbox_index( ndata, FFTsize, idx )
   
   
   ! Time to FFT
   ! standard: -1, and +1
   do colrow = 1, 2
      
      do icol = 1, ndata
      
         a_r( 1:ndata ) = P%C( 1:ndata, icol )
         ! put 1D array to 3D box
         call put_into_fftbox( ndata, a_r, idx, FFTsize, fftbox ) 
         ! FFT WARNING: SIGN CHANGE
         if ( colrow .eq. 1 ) then
            if (iq .eq. 1) call do_fft( fftbox, FFTsize, one)
            if (iq .ne. 1) call do_fft( fftbox, FFTsize, minus_one)
         elseif ( colrow .eq. 2 ) then  
            if (iq .eq. 1) call do_fft( fftbox, FFTsize, minus_one)
            if (iq .ne. 1) call do_fft( fftbox, FFTsize, one)
         endif
         ! save FFTed values to a_r ( but its contents are g-space values )
         call box_to_array( FFTsize, fftbox, a_r )
          ! exchange column to FFTed values
         P%C( 1:ndata, icol ) = a_r( 1:ndata )
      enddo
      
      P%C(1:ndata,1:ndata) = Transpose( P%C(1:ndata,1:ndata) )
            
   enddo

end subroutine
