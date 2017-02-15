!  Sigma_x main program
!
!  Final Modification
!January, 2017
!subhasish.mandal@yale.edu

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
   type(gspace) :: g_coh     ! Polarizability g space
   integer :: FFTsize(3)
   integer :: nk, nq, iq, ik,  g, gp, gpp, grp
   integer :: icm, idm, iqg
   integer :: ikq ! kpt index for k+q vector
   integer :: nfft ! FFTsize(1)*FFTsize(2)*FFTsize(3)
   integer :: ndata ! number of g vector in epsilon matrix   
   integer :: ngdata ! number of g vector in epsilon matrix   
   integer :: ngppdata ! number of g vector to process COH matrix   
   integer, allocatable :: idx(:,:)
   integer, allocatable :: gidx(:) ! link G vectors between polarizability and epsilon
   integer, allocatable :: gpidx(:) ! link G vectors between polarizability and epsilon
   integer :: counter ! link G vectors between COH and epsilon
   integer :: Uvec(3) ! Umklapp vector
   complex(dp), allocatable :: fftbox(:,:,:)
   integer, allocatable :: gppvec(:,:) ! g vectors
   integer, allocatable :: gppvec_gp(:,:) ! gp vectors after matching g-g'
   integer, allocatable :: gppvec_g(:,:) ! g vectors after matching g-g'
   integer, allocatable :: grvec(:,:) !reduced g vectors for COH
   integer, allocatable :: gvec(:,:) ! g vectors
   !integer, allocatable :: G_vec(:,:) ! g vectors
   !integer, allocatable :: Gp_vec(:,:) ! g vectors
   integer:: i,j, idata, ii
   complex(dp), dimension(:), allocatable ::a_n
   complex(dp), dimension(:), allocatable ::a_v
   real(dp), dimension(:), allocatable ::coulb
   real(dp), dimension(:), allocatable ::coulb_gp
   real(dp), dimension(:), allocatable ::SR
   real(dp), dimension(:), allocatable ::E_spl
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
   nbmax = 1 !not implemented yet
   nbmin = 1 !not implemented yet
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
   print *, 'Allocating a pile'

   ndata = FFTsize(1)*FFTsize(2)*FFTsize(3)
   ngppdata = ndata*ndata
   allocate(a_n(ndata))
   print *, 'I am allocating a_v of size',ndata
   allocate(a_v(ndata))
   print *, 'Allocating a pile-----done-----'
   allocate(fr(ndata))
   allocate(frp(ndata))
   allocate(fg(ndata))
   allocate(fgp(ndata))
   allocate(S%C(ndata,ndata) )
   allocate(SR(ngppdata) )
   allocate(E_spl(ngppdata) )
 allocate(Tg(ndata) )
   allocate(Tg1(ndata) )
   allocate(coulb(ndata))
   allocate(coulb_gp(ngppdata))
   allocate( gvec(3,ndata) )
   allocate( gppvec(3,ngppdata) )
   allocate( gppvec_gp(3,ngppdata) )
   allocate( gppvec_g(3,ngppdata) )
   allocate( grvec(3,ngppdata) )
   allocate( g_pol%gvec(3,ndata) , gidx(ndata) )
   allocate (sig_mat%mtrx(nq), sig_mat%gs(nq))
   g_pol%ng = ndata
!   g_coh%ng = ngppdata
   call fftidx_to_gidx( FFTsize, ndata, g_pol%gvec )
   !call fftidx_to_gidx( FFTsize, ngppdata, g_coh%gppvec )
!!
   print *, 'Allocating a pile-----done---2--'
!***************************************************************** 
!Zeroing some  variable
sig_sum=0.0
sig_g1=0.0
sig_g2=0.0
fr=0.0
frp=0.0
fg=0.0
fgp=0.0
Tg=0.d0
Tg1=0.d0
   print *, 'norm of k=n=1 in R space is ',sum(abs(psi_R%wk(1)%cg(1:ndata,1))**2)*sys%vol / ndata

do ik = 1, nk ! Major loop for k

!Zeroing some  variable
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
do  ib = 1, 3 ! Major loop for bands (for user to change)

!Zeroing some  variable
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

      
do iq = 1,  nq ! Major Loop for q-points


!Zeroing some  variable

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
         ! Getting the psi(r,ib) - the state for which we are computing sigma
         a_n(1:ndata) = psi_R%wk(ik)%cg(1:ndata,ib)

! This deallocation is for reduce_gspace --starts here ---
deallocate( q%vec )
call get_qvec( k, q )

! Reducing g-space based on the cut-off 
           call reduce_gspace( g_pol, sig_mat%gs(iq), inp, sys, q%vec(1:3,iq), gidx )
          ngdata = sig_mat%gs(iq)%ng
g_pol%gvec(:,:) = gvec(:, :)
!Calling Coulmb  potential again in the main program 
   call calc_coulb( sys, sig_mat%gs(iq)%ng, q%vec(:,iq), sig_mat%gs(iq)%gvec, coulb, nq )
   !Calling  S-matrix 
        call calc_S( S, Eps_GW, ngdata, nq, iq, q%vec(1:3,iq),gidx, sys )
        !Sanity Check 
     print*, "Checking the size of S and fg", size(fg), size(S%C(:, 1))
         print*, 'number of data in sigma calculations:', ngdata

deallocate( q%vec )
call get_qvec( k_, q )

! This deallocation is for reduce_gspace --ends here ---
         
          ! change wavefunction for

          ! U-process

        !  if (utest .ne. 0 ) then 
!
     !     call modify_Uproc_wfn ( ndata, a_v, Uvec,sys, FFTsize )

        !  endif

       !Computing fr  ; Note there is no frp for this term
         fr(1:ndata) = a_n(1:ndata) * conjg( a_n(1:ndata) )*sys%vol/ndata
!FFT to g-space --from fr to fg
         call f_r_to_g( -1, fr, fg, FFTsize, ndata )



         !Writing to disk .....starts....
open(10,file='Sigma_COH.log_ecut14',form='formatted',status='unknown')
open(20,file='Sigma_COH_Details.out_ecut14',form='formatted',status='unknown')
write(20,*) 'Which k-point you are on : ', ik, k_%vec(1:3, ik)
write(20,*) '******************************* '
write(20,*) '******************************* '
write(20,*) '******************************* '
write(20,*)  'Which q-point you are on : ', q%vec(1:3, iq)
write(20,*) '******************************* '
write(20,*) '******************************* '

!write(20,*) ' old-g-index,  old-gvec ---, ----- & 
!new-g-index   ---new-grvec(:,grp)----,  --- &
!fg---,  --------  SR(grp)---- real(SR(grp)*fg(g)), sig_g1 '

! write(20, *)'   g,     gvec(:,gidx(g))     ,  gp,        gpvec(:,gidx(gp))  ,&
! Coulb(gp),      Eps_Inv,   BGW_eps,    fg_Real,   fg_Img,  schx=fg_Real*BGW_Eps, &
! accum(schx)   '

 write(20, *)'   g,     gvec(:,gidx(g))     ,  gp,        gpvec(:,gidx(gp))  ,Coulb(gp), &
   S_R ,   special_eps,    fg_Real,   fg_Img,  schx=fg_Real*special_eps, &
   accum(schx)'

! We always treat iq==0 as a special point since at q=0, coulb blows up 
if (iq == 1) then
    print*, 'Which q-point you are on : ', q%vec(1:3, iq)
    print*, 'Which k-point you are on : ', k_%vec(1:3, ik)
    counter=0
    igp=0
    do gp=1, ngdata
      do g=1, ngdata
      ! I am starting a counter to keep track of  g,g'; it makes a long list of
      ! gxg' elements
      counter=counter+1
      !For silicon we hardcode this coulb(1)
      coulb(1)=0.19171737
      !coulb(1)=0.38343474
      ! Building g''=g-g' 
      gppvec(:,counter)=gvec(:,gidx(g))-gvec(:,gidx(gp))
      ! Feeding SR vector with S-matrix element following counter
      SR(counter)=S%C(g,gp)
      ! Special S-matrix  treatment due to Physics-issue
      E_spl(counter)=S%C(g,gp)

      if ((g .eq. 1) .AND. (gp .ne. 1)  ) E_spl(counter)=0.000
      if ((gp .eq. 1) .AND. (g .ne. 1  )) E_spl(counter)=0.000

      ! Now  making  new coulb-list and g-vec  for list of length counter 
      coulb_gp(counter)=coulb(gp)
      gppvec_g(:,counter)=gvec(:,gidx(g))
      ! Now  making  new coulb-list and g'-vec  for list of length counter 
      gppvec_gp(:,counter)=gvec(:,gidx(gp))


      enddo ! g ends
      enddo !g' ends
! Now we will start computing Matrix elements 
igp=0
!grp=0
    do igp=1, counter
      do g=1, ngdata
      coulb(1)=0.38343474
      !Matching g-g' in the long list and picking up elements
      if( (gppvec(1,igp) .eq. gvec(1,gidx(g))) .AND.(gppvec(2,igp) .eq. gvec(2,gidx(g))) &
      .AND. (gppvec(3,igp) .eq. gvec(3,gidx(g)) )) then
! Summing over 
      sig_g1=  sig_g1 + real(E_spl(igp)*real(fg(gidx(g))))
      !Writting  Simga-matrix elemenet 
      write(20, '(8i8, 12f15.8 )') &  !g, gvec(:,gidx(g)), & ! igp,   gvec(:,igp), gppvec(:,igp), igp  ,&
      g, gppvec_g(:,igp), igp,  gppvec_gp(:,igp), coulb_gp(igp), SR(igp), &
      E_spl(igp), fg(gidx(g)), real(E_spl(igp)*conjg(fg(gidx(g)))), sig_g1
      endif

      enddo
    enddo
 print*, grp

 ! Writing to disk for the accumulated values 
       write(10,*) 'accumulation of sigma at q : ', iq, sig_g1, sig_g1*0.5, sig_g1*Ha2eV*0.5
       write(20,*) 'accumulation of sigma at q : ', iq,  sig_g1, sig_g1*0.5, sig_g1*Ha2eV*0.5

       print*, 'accumulation of sigma at q : ',iq,  sig_g1 !, sig_g2*0.5, sig_g2*Ry2eV*0.5
    
   
   
       !Here goes for rest of the q-vectors 
       !We do exactly same as above expect we compute Coulmb for g=g'=0
       elseif ( iq .ne.1) then 
!
    print*, 'Which q-point you are on : ', q%vec(1:3, iq)
    print*, 'Which k-point you are on : ', k_%vec(1:3, ik)
    counter=0
    do gp=1, ngdata
      do g=1, ngdata
      counter=counter+1
      gppvec(:,counter)=gvec(:,gidx(g))-gvec(:,gidx(gp))
      SR(counter)=S%C(g,gp)
      E_spl(counter)=S%C(g,gp)
      coulb_gp(counter)=coulb(gp)
      gppvec_g(:,counter)=gvec(:,gidx(g))
      gppvec_gp(:,counter)=gvec(:,gidx(gp))
     enddo
    enddo 

!
igp=0
grp=0
sig_g2=0.0
    do igp=1, counter
      do g=1, ngdata
      if( (gppvec(1,igp) .eq. gvec(1,gidx(g))) .AND.(gppvec(2,igp) .eq. gvec(2,gidx(g))) &
      .AND. (gppvec(3,igp) .eq. gvec(3,gidx(g)) )) then

      sig_g2=  sig_g2 + real(E_spl(igp)*real(fg(gidx(g))))
      write(20, '(8i8, 12f15.8 )') &  !g, gvec(:,gidx(g)), & ! igp,   gvec(:,igp), gppvec(:,igp), igp  ,&
      g, gppvec_g(:,igp), igp,  gppvec_gp(:,igp), coulb_gp(igp), SR(igp), &
      E_spl(igp), fg(gidx(g)), real(E_spl(igp)*conjg(fg(gidx(g)))), sig_g2


      endif
      enddo
    enddo
    sig_sum1=sig_sum1 + sig_g2

       write(10,*) 'accumulation of sigma at q : ', iq, sig_g2, sig_g2*0.5, sig_g2*Ry2eV*0.5
       write(20,*) 'accumulation of sigma at q : ', iq, sig_g2, sig_g2*0.5, sig_g2*Ry2eV*0.5
       print*, 'accumulation of sigma at q : ',iq,  sig_g2 !, sig_g2*0.5, sig_g2*Ry2eV*0.5


!
endif

enddo !end of iq
!Summing up contribution for iq=1 and iq>1
sig_sum= sig_g1 + sig_sum1
print*, 'print sigma over all q at k', k_%vec(1:3, ik), ib, '::', sig_sum, &
sig_sum*Ha2eV*0.5
write(20,*)  ' Dealing with K-point', k_%vec(1:3, ik)
write(20,*)  '****************'
write(20,*)  'Band index, Sigma'
write(20, '(i8,2f15.8)' ) ib,  -sig_sum*Ha2eV*0.5
write(10,*)  ' Dealing with K-point', k_%vec(1:3, ik)
write(10,*)  '****************'
write(10,*)  'Band index, Sigma'
write(10, '(i8, 2f15.8)' ) ib, sig_sum, sig_sum*Ha2eV*0.5

!
enddo ! end of iband
!!
enddo !end of ik
close(10)
close(20)



deallocate(a_n)

!130 continue
end program

