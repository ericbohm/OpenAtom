!  Sigma_x main program
!
!  Dec 2016  subhasish.mandal@yale.edu

program sigma

   use constant
   use usrinput
   use electronic_structure
   use gw_structure

   implicit none
   type(input) :: inp
   type(sysinfo) :: sys, sys_shft
   type(wfstruc) :: psi, psi_R, psi_shft, psi_R_shft
   type(wfn) :: psi_v,  psi_N 
   type(epsilon_type) :: eps, sig_mat
   ! k is in cartesian unit, k_ is in crystal coordinate
   type(kptinfo) :: k, q, k_, & 
                    k_shft
   real(dp) ::  Ha2ev
   real(dp) ::  Ry2ev
   type(gspace) :: g_pol     ! Polarizability g space
   integer :: FFTsize(3)
   integer :: nk, nq, iq, ik,  g
   integer :: icm, idm, iqg
   integer :: ikq ! kpt index for k+q vector
   integer :: nfft ! FFTsize(1)*FFTsize(2)*FFTsize(3)
   integer :: ndata ! number of g vector in epsilon matrix   
   integer :: ngdata ! number of g vector in epsilon matrix   
   integer, allocatable :: idx(:,:)
   integer, allocatable :: gidx(:) ! link G vectors between polarizability and epsilon
   integer :: Uvec(3) ! Umklapp vector
   complex(dp), allocatable :: fftbox(:,:,:)
   integer, allocatable :: gvec(:,:) ! g vectors
   integer:: i,j, idata, ii
   complex(dp), dimension(:), allocatable ::a_n
   complex(dp), dimension(:), allocatable ::a_v
   real(dp), dimension(:), allocatable ::coulb
   complex(dp), dimension(:), allocatable ::fr, fg, fgp, frp, f_r
   integer :: nv, nc, nb, nbmax, nbmin
   integer :: iv, ic, ib
   integer :: utest
   real(dp) ::  sig_x, sig_g, sig_g1, sig_g2
   real(dp) ::  sig_sum, sig_sum1, sig_sum2
!*****************************************
integer, allocatable :: gv(:,:), gpv(:,:)
integer :: ig, igp, test, kk
real(dp) :: g2, gp2
!*****************************************
!*******************
!******Band info*******
   nv = 4 ! This is for Silicon example
   nb = 52 ! This is for Silicon example
   nc = nb - nv
   nbmax = 1 ! At present not implemented 
   nbmin = 1 ! At present not implemented 
!********************
!********************
!Some Constants
   Ha2ev=27.211382543519
   Ry2ev=13.60569127
   utest = abs( Uvec(1) ) + abs( Uvec(2) ) + abs( Uvec(3) )
!********************
!********************
   ! Read input values
   call read_input( inp )
   ! Read wavefunction data and save into the memory
   call read_wfn( inp%wfname, sys, psi, k )

   ! check_dim is for testing.
   ! call check_dim( psi )

   ! This routine gives you the size of the FFT box
   call set_FFTsize( psi, FFTsize )

   nk = k%nk
   sys%nk = nk
   nq = nk  ! nq = nk
   q%nk = nq
   print*,'check Nq', nq, nk
   !stop

   ! retrieve k vectors in reciprocal lattice basis from cartesian basis
   ! (unit is 2*pi/a for both)
   ! relevant for QE output
   print *, 'Calling cartesian_to_crystal'
   call cartesian_to_crystal( sys, k, k_ )
   !call cartesian_to_crystal( sys, k_, k )
   ! Let's get k and q vector information
   print *, 'get_qvec'
   call get_qvec( k_, q )
   !call get_qvec( k, q )
   !print *, 
!   if (inp%crystal) then
!         q%vec(1:3,1) = k%vec(1:3,1) - k_shft%vec(1:3,1)
!   endif

   print *, 'norm of k=n=2 in g space is ',sum(abs(psi%wk(2)%cg(:,2))**2)


   !!print*, 'Calculate sigma at q:', q%vec(1:3,iq)
   ! Do FFT for all k points and bands
   ! psi (in gspace) will be removed from the memory
   !here psi==psi_G as input and psi_R as output
   print *, 'FFT_all_bands'
   call FFT_all_bands( psi, psi_R, FFTsize, sys )

   print *, 'Allocating a pile'

   ndata = FFTsize(1)*FFTsize(2)*FFTsize(3)
   allocate(a_n(ndata))
   print *, 'I am allocating a_v of size',ndata
   allocate(a_v(ndata))
   print *, 'Allocating a pile-----done-----'
   allocate(fr(ndata))
   allocate(frp(ndata))
   allocate(fg(ndata))
   allocate(fgp(ndata))
   allocate(coulb(ndata))
   allocate( gvec(3,ndata) )
   allocate( g_pol%gvec(3,ndata) , gidx(ndata) )
   allocate (sig_mat%mtrx(nq), sig_mat%gs(nq))
   g_pol%ng = ndata
!!
   print *, 'Allocating a pile-----done---2--'
!***************************************************************** 
!Zeroing bunch of variables
sig_sum=0.0
sig_g1=0.0
sig_g2=0.0
fr=0.0
frp=0.0
fg=0.0
fgp=0.0
! Sanity check 
   print *, 'norm of k=n=1 in R space is ',sum(abs(psi_R%wk(1)%cg(1:ndata,1))**2)*sys%vol / ndata

   !K-loop starts here
do ik = 1, nk

sig_g=0.0
sig_g1=0.0
sig_g2=0.0
sig_sum1=0.0
sig_sum2=0.0
sig_sum=0.0

fr=0.0
frp=0.0
fg=0.0
fgp=0.
!band(ib)-loop starts here
!
! These are the bands you want to compute for --this is user input but
!presently resides inside the code 
!
do  ib = 1, 2

sig_g=0.0
sig_g1=0.0
sig_g2=0.0
sig_sum1=0.0
sig_sum2=0.0
sig_sum=0.0

fr=0.0
frp=0.0
fg=0.0
fgp=0.0
!q-loop starts here
!At the end we need to sum up the sigma for all qs. 
      
do iq = 1, nq


fr=0.0
frp=0.0
fg=0.0
fgp=0.0

         call fftidx_to_gidx( FFTsize, ndata, gvec )
         print *, 'done fft_to_gidx'
         ! let's calculate which k vector is the same as k+q
         call get_kq_index( k_, q, ik, iq, ikq, Uvec )
         !print *, 'ik,iq,ikq,Uvec',ik,iq,ikq,Uvec
         ! psi(r,ib) - the state for which we are computing sigma
         a_n(1:ndata) = psi_R%wk(ik)%cg(1:ndata,ib)
         !print *, 'a_n(1) = ',a_n(1)*sys%vol/ndata
         !print *, 'a_n(1) = ',a_n(1)
         !Now calling valence bands i.ei psi(r, iv) 
         !print *, 'calling  a_v_r'
         !call calc_v(psi_R%wk(ik), psi_R%wk(ikq), f, FFTsize, sys, Uvec )
         !print *,'zeroing a_v'
         !call flush(6)
         !a_v = 0.0d0


         ! starts valence bands loop here
         do iv=1, nv
         !print *, 'a_v(1) = ', a_v(1), size(a_v)
         ! print *, 'a_v(1)*a_v(1) = ',a_v(1)*sys%vol/ndata*a_v(1)*sys%vol/ndata


          a_v = 0.0d0
          print *,' considering iv=', iv
          a_v(1:ndata) = psi_R%wk(ikq)%cg(1:ndata, iv)


            ! change wavefunction for U-process
         if (utest .ne. 0 ) then 
           call modify_Uproc_wfn ( ndata, a_v, Uvec, sys, FFTsize )

         endif
         ! Forming fr : fr= (psi_n)*(psi_l) at r
         fr(1:ndata) = a_v(1:ndata) * conjg( a_n(1:ndata) )*sys%vol/ndata
         ! Forming frp : fr= (psi_n)*(psi_l) at r'
         frp(1:ndata) = a_v(1:ndata) * conjg( a_n(1:ndata) )*sys%vol/ndata
         !print *, 'computed f_r', size(fr), ndata
         !f(1:ndata) = fr(1:ndata)
         !print *, 'computed f_r', fr(1) 
         !print *, 'computed f_rp', frp(1) 

         ! This routine FFTs fr to fg
         call f_r_to_g( -1, fr, fg, FFTsize, ndata )
         ! if you want to print some fg
!         print *,'some g,gvec,fg(g)'
!         do g=1,5
!           print *,g,gvec(:,g),fg(g)
!         enddo

         !frp(1:ndata) = a_v(1:ndata) * conjg( a_n(1:ndata) )*sys%vol/ndata
         call f_r_to_g( -1, frp, fgp, FFTsize, ndata )



! g-space reduction is not implemented here for Bare exachange 
!         g_pol%gvec(:,:) = gvec(:, :)
!         if (ik==1 .and. iv==1) then
!           call reduce_gspace( g_pol, sig_mat%gs(iq), inp, sys, q%vec(1:3,iq), gidx )
!         endif
!         ngdata = sig_mat%gs(iq)%ng
         ngdata = ndata
         print*, 'number of data in sigma calculations:', ngdata
         call calc_coulb( sys, ngdata, q%vec(:,iq), gvec, coulb, nq )

!print*, 'print out coulmb', ik
! ########PRINTING Results STARTS########
! !*****************************************
!print*, 'print out coulmb', ik
!!!
open(1,file='Sigma.out',form='formatted',status='unknown')
open(2,file='Sigma.log',form='formatted',status='unknown')
write(1,*) 'Which k-point you are on : ', ik, k_%vec(1:3, ik)
write(1,*) '******************************* '
write(1,*) '******************************* '
write(1,*)  'Which q-point you are on : ', q%vec(1:3, iq)
write(1,*) '******************************* '
write(1,*) '******************************* '
write(1,*) 'start for in =    ', iv
write(1,*) ' Vcoulb  ABS(D)  sum_sig_g   n1   ig   ga  gb  gc      '
! We always consider q==1 a special point. At present coulb(g) at q==1 is
! hardcoded to a particular value for Silicon
if (iq == 1) then
    print*, 'Which q-point you are on : ', q%vec(1:3, iq)
    print*, 'Which k-point you are on : ', k_%vec(1:3, ik)
    do g = 1, ngdata
    coulb(1)=0.38343474
    sig_g1=  sig_g1  + coulb(g)*real(abs(fg(g)*conjg(fgp(g)))) 
    write(1, '(3f15.8,610i8)') coulb(g),real(abs(fg(g)*conjg(fgp(g)))),sig_g1, iv, g, gvec(:,g)
    enddo !g ends
    write(1,*) 'accumulation of sigma at q : ', sig_g1
    write(1, '(1f15.8)') sig_g1+sig_g
elseif (iq .ne. 1 ) then 
print*, 'Which q-point you are on : ', q%vec(1:3, iq)
print*, 'Which k-point you are on : ', k_%vec(1:3, ik)


do g = 1, ngdata 

sig_g2=  sig_g2 + coulb(g)*real(abs(fg(g)*fgp(g))) 
sig_sum=sig_sum + sig_g
write(1, '(3f15.8,610i8)') coulb(g),real(abs(fg(g)*conjg(fgp(g)))), sig_g2, iv, g, gvec(:,g)
enddo
write(1,*) 'accumulation of sigma at q : '
write(1, '(1f15.8)') sig_g2
endif


enddo !end of iv

enddo !end of iq


do iq=1, nq
sig_sum= sig_g1 + sig_g2 + sig_g
enddo
print*, 'print sigma over all q at k', k_%vec(1:3, ik), '::', -sig_sum, -sig_sum*Ha2eV
write(1,*) 'band index ', ib, 'K-point ', k_%vec(1:3, ik)
write(2,*)  ' Dealing with K-point', k_%vec(1:3, ik)
write(2,*)  '****************'
write(2,*)  'Band index, Sigma'
write(1, '(2f15.8)' ) -sig_sum, -sig_sum*Ry2eV
write(2, '(i8,2f15.8)' ) ib,  -sig_sum*Ry2eV

enddo ! end of iband
!!
enddo
close(1)
close(2)



deallocate(a_n)

end program
