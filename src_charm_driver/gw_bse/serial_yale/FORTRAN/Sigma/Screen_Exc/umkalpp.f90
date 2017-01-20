!  Oct. 2014    minjung.kim@yale.edu
!  This subroutine do FFTs for each bands and update Pmtrx
!
!  Spin is not yet implemented

!subroutine calc_v( psi_v, a_v, FFTsize, sys, Uvec, ndata, iv)
!   use constant
!   use electronic_structure
!   use gw_structure
   
!   implicit none
!   !type(wfn), intent(in) :: psi_v,  psi_N
!   type(wfn), intent(in) :: psi_v
!   complex(dp), intent(out) :: a_v(ndata)
!   !type(rank2_mtrx), intent(inout) :: Cmtrx
!   !type(rank2_mtrx)  :: C
!   integer, intent(in) :: FFTsize(3)
!   type(sysinfo), intent(in) :: sys
!   integer, intent(inout) :: Uvec(3), ndata
!
!   ! work variables
!   integer :: nv, nc, nb, nbmax, nbmin
!   integer :: iv, ic, ib
!   integer :: utest
!   
!   ! Change it ------------------------------
!   !nv = int ( sum( psi_v%occ ) )
!   nv = 4
!   !nb = int ( sizeof( psi_c%eig(:) ) )/ 8 
!   nb = 52
!   nc = nb - nv
!   nbmax = 1
!   nbmin = 1
!   !-----------------------------------------
!
!   ndata = fftsize(1)*fftsize(2)*fftsize(3)
!   utest = abs( Uvec(1) ) + abs( Uvec(2) ) + abs( Uvec(3) )
!
!   !print *,'ndata = ',ndata
!   !print*, 'number of bands, nv, nc, nbmax:', nb, nv, nc, nbmax
!     
!!   print *, 'psi_v', psi_v%cg(1,1)
!!   print *, 'psi_v', psi_v%cg(ndata,1)
!!   call flush(6)
!
!!   print *,'a_v(1)',a_v(1)
!!   do iv=1,ndata
!!      print *, 'av ',iv,' is',a_v(iv)
!!      call flush(6)
!!   end do
!   
!
!!     print*,'some norms'
!!     do iv=1,4
!!       a_v = 0.0d0
!!       a_v(1:ndata) = psi_v%cg(1:ndata,iv)
!!       print *,'iv = ',iv,sum(abs(a_v(:))**2)*270.10614592123204/ndata
!!     enddo
!
!     !Bmtrx%C = 0.0d0
!     do iv =  1,  nv
!     !do iv =  2,  2
!
!        !print *, 'doing iv=',iv
!        !print *, 'zeoring'
!        a_v = 0.0d0
!        print *,' considering iv=', iv
!        a_v(1:ndata) = psi_v%cg(1:ndata,iv)
!       ! print *, 'utest is ',utest
!        ! change wavefunction for U-process
!        if (utest .ne. 0 ) then 
!           call modify_Uproc_wfn ( ndata, a_v, Uvec, sys, FFTsize )
!        endif
!
!!        do ib =1, nbmax
!!         a_n(1:ndata) = psi_N%cg(1:ndata,ib)
!!         call update_f( a_v, a_n, ndata, f)
!
!!     enddo
!    enddo
!!   contains 
!   subroutine update_f( a_v, a_n, ndata, f )
!
!      complex(dp), dimension(ndata), intent(inout) :: a_v
!      complex(dp), dimension(ndata), intent(inout) :: a_n
!      integer, intent(in) :: ndata
!      !complex(dp), dimension(ndata,ndata), intent(inout) :: B
!      complex(dp), dimension(ndata), intent(inout) :: f 
!      ! work variables
!      integer :: i, j
!    !  real(dp) :: fact
!
!! for n=l 
!!for  n =B,
!!      f(1:ndata) = a_v(1:ndata) * conjg( a_n(1:ndata) )
!
!!      fact = 2.d0/dble(nk)/(Ev-Ec)
!!      diagonal element of the  X
!     
!    ! Bmtrx%C = 0.0d0
!      do i = 1, ndata
!
!         f(i) = a_v(i) * conjg( a_n(i) ) + f(i)
!            
!      enddo
!   end subroutine update_f
   



   subroutine modify_Uproc_wfn( ndata, a_c, Uvec, sys, FFTsize )

   use constant
   use electronic_structure
   use gw_structure
   
   implicit none
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
   
   
!end subroutine! - End of Calc_Crrp subroutine

