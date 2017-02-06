! 
! FFT subroutines  (uses FFTW2.x)
! Sep. 2014   
! Adapted from Sohrab's fft subroutines
!  modified by MK
!

subroutine do_fft( fftbox, Nfft, sign )

   use constant
   implicit none
   integer, intent(in) :: Nfft(3)
   integer, intent(in) :: sign
   double complex,intent(inout) :: fftbox( Nfft(1), Nfft(2), Nfft(3) )
   
   ! Get the FFTW constant definitions we need
   include 'fftw_f77.i'

   ! local vars
   integer, dimension(3), save :: Nfftold = 0
   integer*8, save:: plus_plan, minus_plan


   if (Nfftold(1)/=Nfft(1) .or. Nfftold(2)/=Nfft(2) .or. &
       Nfftold(3)/=Nfft(3) ) then
      call fftwnd_f77_create_plan(plus_plan,3,Nfft,FFTW_BACKWARD, &
          FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
      call fftwnd_f77_create_plan(minus_plan,3,Nfft,FFTW_FORWARD, &
          FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
      Nfftold(:) = Nfft(:)
   endif

   if (sign == 1) then
      call fftwnd_f77_one(plus_plan,fftbox,0)
   else if (sign == -1) then
      call fftwnd_f77_one(minus_plan,fftbox,0)
   else
      call print_error_exit('sign is not 1 or -1 in do_FFT')
   endif

end subroutine



! put a_g(:) into 3D box
subroutine put_into_fftbox( ndata, a_g, idx, FFTsize, fftbox )

   use constant
   implicit none
   
   integer, intent(in) :: ndata
   complex(dp), intent(in) :: a_g(ndata)
   integer, intent(in) :: idx(3,ndata)
   integer, intent(in) :: FFTsize(3)
   complex(dp), intent(inout) :: fftbox(FFTsize(1),FFTsize(2),FFTsize(3))

   integer :: i, test

   ! initialize fftbox
   fftbox = 0.d0
   
   do i = 1, ndata
      test = sum( idx )
      ! idx = -1 if gidx is not within FFT size
      ! so we reject a_g(i) if idx(:,i) is less than 0.
      if ( test .lt. 0 ) then
         continue  ! cycle may be used..?
      else
         fftbox(idx(1,i),idx(2,i),idx(3,i)) = a_g(i)
      endif
   enddo

end subroutine




! this is called in Psi_G_to_R subroutine !
! Note: In general, FFT box is bigger than max(gidx), but here,
! FFT box may be smaller than max(gidx). 
subroutine gidx_to_fft_index( gidx, FFTsize, idx )
   
   use constant
   implicit none
   
   integer, intent(in) :: gidx(3)
   integer, intent(in) :: FFTsize(3)
   integer, intent(out) :: idx(3)
   integer :: i 
   integer :: imax ! Reject if gvec is larger than this 
   integer :: include_data(3)   ! true if  |gidx| =< cutoff (related to FFT size)
   integer :: upper, lower
   
   
   ! Check if gidx(1:3) is in our cutoff (~ FFT size)
   include_data(1:3) = -1  ! initialization
   do i = 1, 3
      upper = FFTsize(i)/2
      lower = int(-1)*FFTsize(i)/2
      if (gidx(i) < upper .and. gidx(i) >= lower) then
         include_data( i ) = 0
      endif
   enddo
   
   ! Comments: Let's consider an example that FFTsize = 8
   ! Above block returns TRUE if gidx is within the range of [-4,-3,-2,-1,0,1,2,3]
   ! then, [-4,-3,-2,-1,0,1,2,3] -> [5,6,7,8,1,2,3,4] 
   if ( sum( include_data(1:3) ) .eq. 0 )  then  
      do i = 1, 3
         idx(i) = gidx(i) + 1
         if (idx(i) <= 0 ) then
            idx(i) = idx(i) + FFTsize(i)
         endif
      enddo
   else
      idx(:) = -1 ! We need to reject this one when put_into_fftbox is called. 
   endif
      
end subroutine
   

   
   

! fftbox -> a_r

subroutine box_to_array( FFTsize, fftbox, a_r )
   
   use constant
   implicit none
   
   integer, intent(in) :: FFTsize(3)
   complex(dp), intent(in) :: fftbox( FFTsize(1), FFTsize(2), FFTsize(3) )
   complex(dp), intent(inout) :: a_r( FFTsize(1)*FFTsize(2)*FFTsize(3) )

   ! work variables
   integer :: i, j, k
   integer :: ni, nj, nk
   integer :: istart, iend

   
   ni = FFTsize(1)
   nj = FFTsize(2)
   nk = FFTsize(3)

   do k = 1, nk
      do j = 1, nj
         istart = ni*nj*(k-1) + ni*(j-1) + 1
         iend = ni*nj*(k-1) + ni*j
         a_r( istart: iend ) = fftbox(:, j, k) 
      enddo
   enddo
   
end subroutine





subroutine fftidx_to_gidx( FFTsize, nfft, gidx )

   implicit none
   integer, intent(in) :: FFTsize(3)
   integer, intent(in) :: nfft
   integer, intent(inout) :: gidx(3,nfft)
   integer :: idx(3)
   integer :: i, j, k, ijk, test
   
   ! FFTidx - > gidx
   ! [1,2,3,4,5,6,.....N] -> [0,1,2,3,4,N/2-1,-N/2,-N/2+1, -N/2+2,...-1]
   
   ! Let's get 3d index from 1d index
   ijk = 0
   do k = 1, FFTsize(3)
      do j = 1, FFTsize(2)
         do i = 1, FFTsize(1) 
            ijk = ijk + 1
            test = i + (j-1)*FFTsize(1) + (k-1)*FFTsize(1)*FFTsize(2)
            if ( test .ne. ijk) then
               print*, 'Error in fftidx_to_gidx at :', i, j, k
            else
               gidx(1,ijk) = i
               gidx(2,ijk) = j
               gidx(3,ijk) = k
            endif
         enddo
      enddo
   enddo   
   
   do ijk = 1, nfft
      do i = 1, 3
         gidx(i,ijk) = gidx(i,ijk)-1
         if ( gidx(i,ijk) >= FFTsize(i)/2 ) then
            gidx(i,ijk) = gidx(i,ijk)- FFTsize(i)
         endif
      enddo    
   enddo
 
end subroutine




subroutine set_3Dbox_index( ndata, FFTsize, idx )
   implicit none
   integer, intent(in) :: ndata
   integer, intent(in) :: FFTsize(3)
   integer, intent(inout) :: idx(3, ndata)
   ! work variables
   integer :: i, j, k, ii

   ! set index
   ii = 0 
   do k = 1, FFTsize(3)
      do j = 1, FFTsize(2)
         do i = 1, FFTsize(1)
            ii = ii + 1
            idx(1,ii) = i
            idx(2,ii) = j
            idx(3,ii) = k
         enddo
      enddo
   enddo
end subroutine
