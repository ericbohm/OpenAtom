!  Oct. 2014    minjung.kim@yale.edu
!  This subroutine do FFTs for each bands and update Pmtrx
!
!  Spin is not yet implemented

! outer loop is in 'k'

subroutine calc_Prrp_ffreq( psi_v, psi_c, Pmtrx, FFTsize, sys, Uvec, ffreq, iE )
   
   use constant
   use electronic_structure
   use gw_structure
   
   implicit none
   type(wfn), intent(in) :: psi_v, psi_c
   type(rank2_mtrx), intent(inout) :: Pmtrx
   integer, intent(in) :: FFTsize(3)
   type(sysinfo), intent(in) :: sys
   integer, intent(in) :: Uvec(3)
   type(full_frequency), intent(inout) :: ffreq
   integer, intent(in) :: iE   ! used to calculate Eeval

   ! work variables
   complex(dp), dimension(:), allocatable :: a_v, a_c
   integer :: ndata
   integer :: nv, nc, nb
   integer :: iv, ic
   integer :: utest
   complex(dp) :: Eeval
   real(dp) :: Ereal
   real(dp),parameter :: Hartree = 27.211385
   
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

   ! Set energy value to evaluate P(q,w)
   Ereal = ffreq%Emin + dble( iE - 1 ) * ffreq%Estep
   Eeval = cmplx( Ereal/Hartree, ffreq%Ebrdn/Hartree )
print*, 'Evaluation energy:', Eeval

   ! no spin

   ! loop for valence bands
   do iv = 1, nv

      a_v(1:ndata) = psi_v%cg(1:ndata,iv)
      
      ! change wavefunction for U-process
      if (utest .eq. 0 ) then 
         continue
      else
         call modify_Uproc_wfn_ffreq ( ndata, a_v, Uvec, sys, FFTsize )
      endif

      ! loop for conduction bands
      do ic = 1, nc

         a_c(1:ndata) = psi_c%cg(1:ndata,nv+ic)

         call update_P_ffreq( a_v, a_c, ndata, psi_v%eig(iv), &
                              psi_c%eig(nb-nc+ic), Pmtrx%C , sys%nk, Eeval )    

      enddo

   enddo
   
   
   contains 
   ! ************************************
   subroutine update_P_ffreq( a_v, a_c, ndata, Ev, Ec, P, nk, Eeval )

      complex(dp), dimension(ndata), intent(inout) :: a_v, a_c 
      integer, intent(in) :: ndata
      real(dp), intent(in) :: Ev, Ec
      complex(dp), dimension(ndata,ndata), intent(inout) :: P
      integer, intent(in) :: nk
      complex(dp), intent(in) :: Eeval
   
      ! work variables
      complex(dp), dimension(ndata) :: f
      integer :: i, j
      complex(dp) :: fact


      f(1:ndata) = a_v(1:ndata) * conjg( a_c(1:ndata) )

      fact = 2.d0 * ( 1/( Ev-Ec-Eeval ) + 1/( Ev-Ec+Eeval ) )

      do j = 1, ndata
         do i = 1, ndata

            p(i,j) = fact * f(i) * conjg( f(j) ) + p(i,j)
            
         enddo
      enddo
   end subroutine update_P_ffreq
   
   
   subroutine modify_Uproc_wfn_ffreq( ndata, a_c, Uvec, sys, FFTsize )
      
      implicit none
      integer, intent(in) :: ndata
      complex(dp), dimension(ndata), intent(inout) :: a_c
      integer, dimension(3), intent(in) :: Uvec
      type(sysinfo), intent(in) :: sys
      integer, dimension(3), intent(in) :: FFTsize
      
      ! work variables
      complex(dp), allocatable :: fftbox(:,:,:), Ufac(:,:,:)
      integer :: i, j, k, ii
      integer, allocatable :: idx(:,:)
      real(dp) :: a(3,3), b(3,3)
      real(dp) :: rijk(3), G0(3), phase
      
      allocate( idx(3,ndata) )
      allocate( fftbox(fftsize(1),fftsize(2),fftsize(3)) )
      allocate( Ufac(fftsize(1),fftsize(2),fftsize(3)) )
      
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
               Ufac(i,j,k) = cmplx( cos(phase), sin(phase) )
            enddo
         enddo
      enddo
      fftbox(1:fftsize(1),1:fftsize(2),1:fftsize(3)) =    &
            Ufac(1:fftsize(1),1:fftsize(2),1:fftsize(3))*    &
            fftbox(1:fftsize(1),1:fftsize(2),1:fftsize(3)) 
         
      
      call box_to_array( FFTsize, fftbox, a_c )
   
      deallocate( idx, fftbox, Ufac )
      
   end subroutine modify_Uproc_wfn_ffreq 
   
end subroutine! - End of Calc_Prrp_ffreq subroutine



