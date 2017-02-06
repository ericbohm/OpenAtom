! generalized plasmon-pole subroutines
! Feb.2015    minjung.kim@yale.edu

!**************************************************************
!                     construct S matrix

subroutine calc_Smtrx( ndata, mtrx, Vcoulb )

   use constant
   use gw_structure
   
   ! input variables
   integer ,intent(in) :: ndata 
   type(rank2_mtrx), intent(inout) :: mtrx
   real(dp),dimension(ndata), intent(inout) :: Vcoulb
   
   ! work variables
   integer :: i, j

   Vcoulb(:) = sqrt( Vcoulb(:) )
   
   do i = 1, ndata
      
      do j = 1, ndata

         if ( i .eq. j ) mtrx%C(i,j) = mtrx%C(i,j) - one
      
         mtrx%C(i,j) = Vcoulb(i) * mtrx%C(i,j) * Vcoulb(j)

      enddo
   
   enddo

end subroutine



!**************************************************************
!              eigen decomposition of S matrix
! LAPACK SUBROUTINE
! http://www.netlib.org/lapack/explore-html/df/db2/cheev_8f.html

subroutine eigen_decomposition( ndata, eigvec, eigval )

   use constant
   use gw_structure
   implicit none

   external ZHEEV
   
   ! input variables
   integer, intent(in) :: ndata
   type(rank2_mtrx), intent(inout) :: eigvec
   type(onedim_array), intent(inout) :: eigval
   ! LAPACK input variables
   character(len=1) :: JOBZ, UPLO
   complex(dp), dimension(:), allocatable :: WORK
   real(dp), dimension(:), allocatable :: RWORK
   integer :: INFO, LWORK, LWMAX
   
   JOBZ = 'V' ! find eigenvalues and eigenvectors
   UPLO = 'L' ! A contains lower triangle
   
   ! on exit
   ! eigval contains eigenvalues in ascending order 
   ! eigvec%C contains eigenvectors of S mtrx
   LWORK = -1; LWMAX = 250*ndata
   allocate( WORK(LWMAX), RWORK(3*ndata-2) )

   ! inquire optimal size of LWORK
   call ZHEEV(JOBZ, UPLO, ndata, eigvec%C, ndata, eigval%R, WORK, LWORK, RWORK, INFO )

   LWORK = MIN( LWMAX, INT( WORK(1) ) )
   call ZHEEV(JOBZ, UPLO, ndata, eigvec%C, ndata, eigval%R, WORK, LWORK, RWORK, INFO )

     
   if ( INFO .gt. 0 ) then
      print*, ''
      print*, '********** WARNING *************'
      print*, 'eigen decomposition algorithm failed to converge.'
      print*, INFO,'off-diagonal elements of an intermediate tridiagonal form did not converge to zero'
   elseif ( INFO .lt. 0 ) then
      print*, ''
      print*, '********** WARNING *************'
      print*, -INFO, 'th argument has illegal value'
   elseif ( INFO .eq. 0 ) then
      print*, 'eigen decomposition succeed!'
   endif
   
   deallocate( WORK, RWORK )
   
end subroutine








!**************************************************************
!              calculate Omega^2
! Let's calculate omegas
! first, you need density. This should be read from rho.dat file
 
subroutine calc_Omsq( inp, sys, ndata, eigvec, eigval, omsq, gvec, qvec )


   use constant
   use gw_structure
   use electronic_structure
   use usrinput
   
   implicit none
   
   type(input), intent(in) :: inp
   type(sysinfo), intent(in) :: sys
   integer, intent(in) :: ndata
   type(rank2_mtrx), intent(inout) :: eigvec
   type(onedim_array), intent(inout) :: eigval
   type(onedim_array), intent(inout) :: omsq
   integer, intent(in) :: gvec(3,ndata)
   real(dp), intent(in) :: qvec(3)
   
   ! work variables
   integer :: nr(3), nfft
   real(dp), allocatable :: rho_real(:,:,:)
   complex(dp), allocatable :: rho(:,:,:)
   integer, allocatable :: rhogidx(:,:)
   integer :: i, j, k, icol, irow
   complex(dp), allocatable :: tmpMat(:,:)
   real(dp) :: norm, qplusG(3), qplusGp(3)
   complex(dp) :: fact, rhoG_Gp, rhoZERO, Wpl2
   integer, parameter :: iu=101
   

   ! First, we need to find rho(G-G') 
   ! Read rho(R) from 'rho.dat' file, then do FFT
   
   if ( inp%gpp_flag_is_on ) then
      open(iu,file=trim(inp%rhoName),form='unformatted',status='unknown')
   else
      call print_error_exit('GPP flag is off. Check your input variables')
   endif
   
   ! read rho.dat file
   read(iu) nr(1:3)
   allocate( rho_real(nr(1),nr(2),nr(3)) )

   do k = 1, nr(3)
      read(iu) ( rho_real(1:nr(1),j,k), j=1,nr(2) )   
   enddo

   ! rho_real= |\psi|^2/unit_volume * nk    (unit_volume = vol / FFTsize)
   ! so we need to change the value here

   ! number of fft points
   nfft = nr(1)*nr(2)*nr(3)
   ! normalization (WARNING: rho_real comes from Quantum Espresso)
   rho_real = rho_real / nfft * sys%vol
   print*, 'rho sum', sum( rho_real)
   ! now, sum( rho_real ) = number of electrons

! ************** Maybe you need to change this ********************** ! 
   allocate( rho(nr(1),nr(2),nr(3)) )
   do k = 1, nr(3)
      do j = 1, nr(2) 
         do i = 1, nr(1)
            rho(i,j,k) = cmplx( rho_real(i,j,k), 0.d0 ) 
         enddo
      enddo
   enddo
   

   ! ********* FFTW **********
   ! rho(R) -> rho(G)

   call do_fft( rho, nr, minus_one ) 

   allocate( rhogidx(3,nfft) )

   ! get gidx
   call fftidx_to_gidx( nr, nfft, rhogidx )

   rho = rho / nfft  ! This scaling may not be true. Check.
   
   ! Omega^2 = U^-1 * MatGG' * U
   
   ! MatGG' needs q, G vector information
   ! q, G, density matrix
   allocate( tmpMat( ndata, ndata ) )
   
   tmpMat = 0.d0
   rhoZERO = rho(1,1,1)

   do icol = 1, ndata
      do irow = 1, ndata
      
         call cryst_to_cart( sys%bvec, sys%alat, qvec, gvec(:,irow), qplusG )
         call cryst_to_cart( sys%bvec, sys%alat, qvec, gvec(:,icol), qplusGp )
         call calc_rhoG_Gp( nfft, nr, rho, rhogidx, gvec(:,irow), gvec(:,icol),&
                            rhoG_Gp )
         call calc_mtxGGP( qplusG , qplusGp, rhoG_Gp, rhoZERO, &
                           tmpMat(irow,icol) )
         
      enddo
   enddo
   
   deallocate( rho, rhogidx )
   
   
   ! Do matrix multiplication \omega^2 = V^-1 * mtxGGp * V

   tmpMat = MATMUL( transpose( conjg(eigvec%C(:,:)) ), tmpMat(:,:) )
   tmpMat = MATMUL( tmpMat(:,:), eigvec%C(:,:) )
   
   ! \omega_{pl}^2 = 4\pi\rho(0)
   Wpl2 = 4.d0 * pi * sum(rho_real)/sys%vol

print*, "printing Omega_pl^2"
print*, Wpl2 

print*, "printing omsq value"
   do i =1, ndata

      omsq%R(i) = -Wpl2 * real( tmpMat(i,i) ) / eigval%R(i) ! -> this one seems to work
      !omsq%R(i) = (-4.d0 * pi * Wpl2) * real( tmpMat(i,i) ) &
      !           / ( eigval%R(i) * sys%vol * sys%nk )
print*, "omsq:", omsq%R(i)
   enddo

      
   deallocate( tmpMat )
   
   
   contains
   
   subroutine calc_rhoG_Gp( ndata, nr, rho, gidx, gv1, gv2, rhoG_Gp )
      
      integer, intent(in):: ndata
      integer, intent(in) :: nr(3)
      complex(dp), intent(in) :: rho(nr(1),nr(2),nr(3))
      integer, intent(in) :: gidx(3,ndata)
      integer, intent(in) :: gv1(3), gv2(3)
      complex(dp), intent(out) :: rhoG_Gp
      ! work variables
      integer :: i,j,k
      integer :: vecG_Gp(3), idxG_Gp, rem, ii, jj, kk
      
      
      vecG_Gp(1:3) = gv1(1:3) - gv2(1:3)



      do i = 1, ndata

         if (gidx(1,i) .eq. vecG_Gp(1) .and. &
             gidx(2,i) .eq. vecG_Gp(2) .and. &
             gidx(3,i) .eq. vecG_Gp(3) ) then
            
            
            ii = gidx(1,i) + 1
            jj = gidx(2,i) + 1
            kk = gidx(3,i) + 1
            
            if (ii .le. 0 ) ii = ii + nr(1)
            if (jj .le. 0 ) jj = jj + nr(2)
            if (kk .le. 0 ) kk = kk + nr(3)
                 
            rhoG_Gp = rho(ii,jj,kk)
            exit
         
         elseif( i .eq. ndata) then
            print*, " G-G' density element does not exist in density matrix"
            print*, "G  :", gv1
            print*, "G' :", gv2
            rhoG_Gp = 0.d0
         endif
      
      enddo

   end subroutine
   
   
   ! (q+G)(q+G')/(|q+G|^2|*|q+G'|^2) * \rho(G-G') / \rho(0)
   subroutine calc_mtxGGP( qG, qGp, rhoG_Gp, rhoZERO, mtxel )
      
      real(dp), intent(in) :: qG(3), qGp(3)
      complex(dp), intent(inout) :: mtxel, rhoG_Gp, rhoZERO

      mtxel = dot_product( qG, qGp )
      mtxel = mtxel / sqrt( dot_product(qG,qG) ) / sqrt( dot_product(qGp,qGp) )
      !mtxel = mtxel / ( dot_product(qG,qG) * dot_product(qGp,qGp) )
      mtxel = mtxel * rhoG_Gp / rhoZERO

   end subroutine

end subroutine





subroutine calc_Smtrx_w( ndata, eig, omsq, Eeval, V, Sw )

   use constant
   integer, intent(in) :: ndata
   real(dp), dimension(ndata), intent(in) :: eig, omsq
   complex(dp), intent(in) :: Eeval
   complex(dp), intent(in) :: V(ndata,ndata)
   complex(dp), intent(inout) :: Sw(ndata,ndata)
   ! local variables
   complex(dp), dimension(ndata) :: tmp

integer :: n(10), in

   
   do i = 1, ndata
      tmp(i) = eig(i) * omsq(i) / (omsq(i) - Eeval**2)
   enddo
 
 
print*, ''
print*, 'number of data:', ndata
print*, 'how many data you want to include (list 10 numbers:)'
!read(*,*) (n(i),i=1,10) 
n(1) = 1; n(2) = 5; n(3)=10; n(4)=15; n(5) = 20; n(6) = 25; n(7) = 30; n(8) = 35; 
n(9) = 40; n(10)=410
!do i = 1, 10
!n(i) = i*40
!enddo

do in  = 1, 10

   Sw = 0.d0
   do j = 1, ndata
      do i = 1, ndata
         !do k = 1, ndata
         do k = 1, n(in)
            Sw(i,j) = Sw(i,j) + V(i,k) * tmp(k) * conjg( V(j,k) )
         enddo
      enddo
   enddo

   write(200,*) n(in), REAL( Sw(1,1)+1 )

enddo



end subroutine



subroutine  get_epsinv_mtrx_from_gpp( nd, Vc, mtrx )

   use constant
   use gw_structure
   
   implicit none
   integer, intent(in) :: nd
   real(dp), intent(in) :: Vc(nd)
   complex(dp), intent(inout) :: mtrx(nd,nd)
   
   integer :: i, j

   do i = 1, nd
      do j = 1, nd
      
         !mtrx(i,j) = mtrx(i,j) / ( sqrt( Vc(i) ) * sqrt( Vc(j) ) ) 
         if ( i .eq. j ) then
            mtrx(i,j) = mtrx(i,j) + 1
         endif
      enddo
   enddo
   
   
   !call inverse( nd, mtrx )
end subroutine
