!  Sigma_x main program
!
! Last modified Oct, 2016
! subhasish.mandal@yale.edu

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
   type(rank2_mtrx) :: S
   type(rank2_mtrx) :: Eps_GW
   type(epsilon_type) :: eps, sig_mat
   type(kptinfo) :: k, q, k_, & 
                    k_shft
   real(dp) ::  Ha2ev, Ry2ev
   type(gspace) :: g_pol     ! Polarizability g space
   integer :: FFTsize(3)
   integer :: nk, nq, iq, ik,  g, gp
   integer :: icm, idm, iqg
   integer :: ikq ! kpt index for k+q vector
   integer :: nfft ! FFTsize(1)*FFTsize(2)*FFTsize(3)
   integer :: ndata ! number of g vector in epsilon matrix   
   integer :: ngdata ! number of g vector in epsilon matrix   
   integer, allocatable :: idx(:,:)
   !integer, allocatable :: ndata, nfft ! number of g vector in epsilon matrix   
   integer, allocatable :: gidx(:) ! link G vectors between polarizability and epsilon
   integer :: Uvec(3) ! Umklapp vector
   !integer :: gvec(3,:) ! g vectors
   complex(dp), allocatable :: fftbox(:,:,:)
   integer, allocatable :: gvec(:,:) ! g vectors
   integer, allocatable :: G_vec(:,:) ! g vectors
   integer, allocatable :: Gp_vec(:,:) ! g vectors
   integer:: i,j, idata, ii
   !complex(dp) :: Dgg
   complex(dp), dimension(:), allocatable ::a_n
   complex(dp), dimension(:), allocatable ::a_v
   real(dp), dimension(:), allocatable ::coulb
   complex(dp), dimension(:), allocatable ::fr, fg, fgp, frp, f_r
   complex(dp), dimension(:), allocatable ::Tg, Tg1
   integer :: nv, nc, nb, nbmax, nbmin
   integer :: iv, ic, ib, utest
   real(dp) ::  sig_x, sig_g, sig_g1, sig_g2
   real(dp) ::  sig_sum, sig_sum1, sig_sum2
!*****************************************
integer, allocatable :: gv(:,:), gpv(:,:)
integer :: ig, igp, test, kk
real(dp) :: g2, gp2
!*****************************************
!*******************
!******Band info*******
   nv = 4
   nb = 52
   nc = nb - nv
   nbmax = 1 ! Cross-diagonal terms not implemented yet
   nbmin = 1 ! Cross-diagonal terms not implemented yet
!********************
!********************
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

  ! CASE1: q = 0 
     iq = 1
     Uvec(1:3) = 0
     if (inp%crystal) then
      call shifted_k_wfn( FFTsize, inp, psi_shft, psi_R_shft, k_shft, sys_shft )


      else
      print*, 'This program does not support confined systems. Program exits.'

      stop

   endif
   !Allocating a pile
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
   allocate(S%C(ndata,ndata) )
   allocate(Eps_GW%C(ndata,ndata) )
   allocate(Tg(ndata) )
   allocate(Tg1(ndata) )
   allocate(coulb(ndata))
   allocate( gvec(3,ndata) )
   allocate( g_pol%gvec(3,ndata) , gidx(ndata) )
   allocate (sig_mat%mtrx(nq), sig_mat%gs(nq))
   g_pol%ng = ndata
   call fftidx_to_gidx( FFTsize, ndata, g_pol%gvec )
!!
!   print *, 'Allocating a pile-----done---2--'
!***************************************************************** 

sig_sum=0.0
sig_g1=0.0
sig_g2=0.0
!Brrp%C=0.0
fr=0.0
frp=0.0
fg=0.0
fgp=0.0
Tg=0.d0
Tg1=0.d0
!Sanity check 
   print *, 'norm of k=n=1 in R space is ',sum(abs(psi_R%wk(1)%cg(1:ndata,1))**2)*sys%vol / ndata

!do ik = 1, nk
do ik = 1, 2 ! Major loops for k-points

! Zeroing some variables
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
do  ib = 1, 2 ! Major loop for bands

sig_g=0.0
sig_g1=0.0
sig_g2=0.0
sig_sum1=0.0
sig_sum2=0.0
sig_sum=0.0


Tg=0.0
Tg1=0.0
fr=0.0
frp=0.0
fg=0.0
fgp=0.0

      
do iq = 1,  nq
!do iq = 2, 2
coulb = 0.0d0

fr=0.0
frp=0.0
fg=0.0
fgp=0.0
         call fftidx_to_gidx( FFTsize, ndata, gvec )
         print *, 'done fft_to_gidx'
         ! let's calculate which k vector is the same as k+q
         call get_kq_index( k_, q, ik, iq, ikq, Uvec )
         !print *, 'ik,iq,ikq,Uvec',ik,iq,ikq,Uvec
    ! To make fr, we first need psi(r,ib) - the state for which we are computing sigma
         a_n(1:ndata) = psi_R%wk(ik)%cg(1:ndata,ib)


deallocate( q%vec )
call get_qvec( k, q )

         g_pol%gvec(:,:) = gvec(:, :)
         ! Reducing g-space --depends on eps_cut
           call reduce_gspace( g_pol, sig_mat%gs(iq), inp, sys, q%vec(1:3,iq), gidx )
         ngdata = sig_mat%gs(iq)%ng
         !Calling Coulb
   call calc_coulb( sys, sig_mat%gs(iq)%ng, q%vec(:,iq), sig_mat%gs(iq)%gvec, coulb, nq )
   !Calling S-matrix
        call calc_S( S,Eps_GW, ngdata, nq, iq, q%vec(1:3,iq),gidx, sys )
     print*, "Checking the size of S and fg", size(fg), size(S%C(:, 1))
         print*, 'number of data in sigma calculations:', ngdata
deallocate( q%vec )

call get_qvec( k_, q ) 
! this was for reduce_gpsace which works only on k but
!not on k_

! Okay, Now starting all valence bands 
! Loop over it 
        do iv = 1, nv
          a_v = 0.0d0
          print *,' considering iv=', iv
         
          a_v(1:ndata) = psi_R%wk(ikq)%cg(1:ndata,iv)



         
          ! change wavefunction for

          ! U-process

          if (utest .ne. 0 ) then 

          call modify_Uproc_wfn ( ndata, a_v, Uvec,sys, FFTsize )

          endif
!Now contrcut fr 
         fr(1:ndata) = a_v(1:ndata) * conjg( a_n(1:ndata) )*sys%vol/ndata
!Now contrcut fr'
         frp(1:ndata) = a_v(1:ndata) * conjg( a_n(1:ndata) )*sys%vol/ndata
!FFT fr to fg
         call f_r_to_g( -1, fr, fg, FFTsize, ndata )
!FFT fr' to fg'

         call f_r_to_g( -1, frp, fgp, FFTsize, ndata )


! Print to DISK

open(10,file='Sigma_SX.out',form='formatted',status='unknown')
open(11,file='Sigma_SX_Details.out',form='formatted',status='unknown')
open(20,file='Sigma.log',form='formatted',status='unknown')
write(10,*) 'Which k-point you are on : ', ik, k_%vec(1:3, ik)
write(10,*) '******************************* '
write(10,*) '******************************* '
write(10,*) '******************************* '
write(10,*)  'Which q-point you are on : ', q%vec(1:3, iq)
write(10,*) '******************************* '
write(10,*) '******************************* '
write(10,*) 'start for in =    ', iv



write(11,*) 'Which k-point you are on : ', ik, k_%vec(1:3, ik)
write(11,*) '******************************* '
write(11,*) '******************************* '
write(11,*) '******************************* '
write(11,*)  'Which q-point you are on : ', q%vec(1:3, iq)
write(11,*) '******************************* '
write(11,*) '******************************* '
write(11,*) 'start for in =    ', iv

write(10,*)'  gp, sig_mat%gs(iq)%gvec(:,gp)  , Tg( gp ), real(Tg(gp))*coulb(gp), sig_g2'

 write(11,*) ' g,  gvec(:,g), gp, gvec(:,gp) , &
   coulb(g), coulb(gp),  Eps_GW%C( g, gp ) &
   , real(fg(gidx(g))*conjg(fg(gidx(gp)))), real(fg(gidx(g))*conjg(fg(gidx(gp))))*Eps_GW%C(g,gp), Tg(g)'
print*, 'Which q-point you are on : ', q%vec(1:3, iq)
print*, 'Which k-point you are on : ', k_%vec(1:3, ik)
print*, 'start for in =    ', iv


if (iq == 1) then
    print*, 'Which q-point you are on : ', q%vec(1:3, iq)
    print*, 'Which k-point you are on : ', k_%vec(1:3, ik)
    do gp=1, ngdata
    Tg=0 ! This is very important 
      do g=1, ngdata
    coulb(1)=0.38343474
  Tg = real(fg(gidx(g))*conjg(fg(gidx(gp))))*S%C(g,gp)+Tg
    write(11, '(8i3,12f15.8 )') g,  sig_mat%gs(iq)%gvec(:,g), gp,  sig_mat%gs(iq)%gvec(:,gp) , &
   coulb(g), coulb(gp),  S%C( g, gp ) &
   , real(fg(gidx(g))*conjg(fg(gidx(gp)))), real(fg(gidx(g))*conjg(fg(gidx(gp))))*S%C(g,gp), Tg(g)
    enddo
  sig_g1=real(Tg(gp)) + sig_g1
    write(10, '(4i3,12f15.8 )') gp, sig_mat%gs(iq)%gvec(:,gp)  , Tg( gp ), real(Tg(gp)), sig_g1
    enddo 

 print*,  'Total ScreenEx accumulated -', -sig_g1
 write(10,* ) 'Total EX-'
 write(10, '(1f15.8)')  -sig_g1*Ry2ev
 ! For iq ne 1 uncomment till endif 
    elseif ( iq .ne.1) then 


    print*, 'Which q-point you are on : ', q%vec(1:3, iq)
    print*, 'Which k-point you are on : ', k_%vec(1:3, ik)
    do gp=1, ngdata
    Tg=0 !This is important
      do g=1, ngdata
    Tg = real(fg(gidx(g))*conjg(fg(gidx(gp))))*S%C(g,gp)+Tg
    write(11, '(8i3,12f15.8 )') g,  sig_mat%gs(iq)%gvec(:,g), gp,  sig_mat%gs(iq)%gvec(:,gp) , &
   coulb(g), coulb(gp),  S%C( g, gp ) &
   , real(fg(gidx(g))*conjg(fg(gidx(gp)))), real(fg(gidx(g))*conjg(fg(gidx(gp))))*S%C(g,gp), Tg(g)
    enddo 

  sig_g2=real(Tg(gp)) + sig_g2
    write(10, '(4i3,12f15.8 )') gp, sig_mat%gs(iq)%gvec(:,gp)  , Tg( gp ), real(Tg(gp)), sig_g2
    enddo 

 write(10,* ) 'Total ScreenEx accumulated -'
 print*,  'Total EX-', -sig_g2
 write(10, '(1f15.8)')  -sig_g2*Ry2ev
endif

!
!
!
enddo !end of iv
enddo !end of iq

do iq=1, nq
sig_sum= sig_g1 + sig_g2 
enddo
print*, 'print sigma over all q at k', k_%vec(1:3, ik), ib, '::', -sig_sum, -sig_sum*Ry2eV
write(20,*) 'band index ', ib, 'K-point ', k_%vec(1:3, ik)
write(20,*)  ' Dealing with K-point', k_%vec(1:3, ik)
write(20,*)  '****************'
write(20,*)  'Band index, Sigma'
write(10, '(2f15.8)' ) -sig_sum, -sig_sum*Ry2eV
write(20, '(i8,2f15.8)' ) ib,  -sig_sum*Ry2eV
!
!


enddo ! end of iband
!!
enddo !end of ik
close(10)
close(20)



deallocate(a_n)
!130 continue
end program

