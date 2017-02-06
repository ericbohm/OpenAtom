! GW calculations - Generalized Plasmon Pole by diagonalization method

! June 2015  minjung.kim@yale.edu

! Full frequency calculation is also included!!!


program gw_gpp

   use constant
   use usrinput
   use electronic_structure
   use gw_structure

   implicit none
   type(input) :: inp
   type(sysinfo) :: sys, sys_shft
   type(wfstruc) :: psi, psi_R, psi_shft, psi_R_shft
   
   ! k is in cartesian unit, k_ is in crystal coordinate
   type(kptinfo) :: k, q, k_, & 
                    k_shft
   type(polarizability) :: Prrp
   type(epsilon_type) :: eps
   type(gspace) :: g_pol     ! Polarizability g space
   type(gpp_type) :: GPP
   type(full_frequency) :: ffreq
   type(full_polarizability) :: P_freq
   type(epsilon_ffreq) :: eps_freq
   
   real(dp), allocatable :: Vcoulb(:) ! coulomb potential

   ! determine FFT grid. This number is multiplied to the minimum cutoff radius
   integer :: FFTsize(3)
   integer :: nk, nq, iq, ik
   integer :: ikq ! kpt index for k+q vector
   integer :: nfft ! FFTsize(1)*FFTsize(2)*FFTsize(3)
   integer :: ndata ! number of g vector in epsilon matrix   
   integer, allocatable :: gidx(:) ! link G vectors between polarizability and epsilon
   integer :: Uvec(3) ! Umklapp vector
   
   integer:: i,j, idata, ii, iE

! **** Modification for GPP 
   integer :: irow
   complex(dp) :: Eeval
   real(dp), parameter :: Hartree = 27.211385


   ! Let's calculate only what I need. Set q point!
   complex(dp), allocatable :: eigval_gpp(:)
   integer, parameter :: qset = 1
   integer :: G1idx(3), G2idx(3), g1, g2
   character(len=100) :: fname


!*****************************************************************  
!    STEP 1:   Epsilon Inverse calculation
!*****************************************************************  


!*****************************************************************  
!                            READING
   ! Read input values
   call read_input( inp )
      
   
   ! Read wavefunction data and save into the memory
   call read_wfn( inp%wfname, sys, psi, k )



!*****************************************************************  
!                        INITIAL SETTING
   call set_FFTsize( psi, FFTsize )

   nk = k%nk
   sys%nk = nk
   nq = nk  ! nq = nk
   q%nk = nq

   ! retrieve k vectors in reciprocal lattice basis from cartesian basis
   ! (unit is 2*pi/a for both)
   ! relevant for QE output
   call cartesian_to_crystal( sys, k, k_ )

   ! Let's get k and q vector information
   call get_qvec( k_, q )
 
   ! Do FFT for all k points and bands
   ! psi (in gspace) will be removed from the memory
   call FFT_all_bands( psi, psi_R, FFTsize, sys )

   Prrp%nq = nq

   ! construct chi matrix
   call create_Pmtrx_struc( Prrp, FFTsize )


   ! If full-frequency flag is on:
   if (inp%ffreq_flag_is_on) then
      call set_full_frequency( ffreq )
      P_freq%nq = nq
      call create_P_freq_struc( P_freq, ffreq, FFTsize )
   endif
   
 
!*****************************************************************  
!                Polarizability matrix calculation


   ! CASE1: q = 0 
   Uvec(1:3) = 0
   if (inp%crystal) then
   
      call shifted_k_wfn( FFTsize, inp, psi_shft, psi_R_shft, k_shft, sys_shft )

      do ik = 1, nk

         call calc_Prrp( psi_R_shft%wk(ik), psi_R%wk(ik), Prrp%mtrx(iq),&
                          FFTsize, sys, Uvec )
         
         if ( ffreq%is_on ) then
            do iE = 1, ffreq%nEstep
               call calc_Prrp_ffreq( psi_R_shft%wk(ik), psi_R%wk(ik), P_freq%mtrx(iE,iq),&
                                    FFTsize, sys, Uvec , ffreq, iE )
            enddo
         endif
      enddo
   
   else
      print*, 'This program does not support confined systems. Program exits.'
      stop
   endif

   ! discard psi_R_shft 
   do ik = 1, nk 
      deallocate( psi_R_shft%wk(ik)%cg ) 
   enddo


   ! CASE2: q \= 0
   do iq = 2, nq
      do ik = 1, nk
         ! let's calculate which k vector is the same as k+q
         call get_kq_index( k_, q, ik, iq, ikq, Uvec )
         
         ! Prrp%mtrx keeps being updated here
         call calc_Prrp(psi_R%wk(ikq), psi_R%wk(ik), Prrp%mtrx(iq), FFTsize, sys, Uvec )

         if ( ffreq%is_on ) then
            do iE = 1, ffreq%nEstep
               call calc_Prrp_ffreq( psi_R%wk(ikq), psi_R%wk(ik), P_freq%mtrx(iE,iq),&
                                     FFTsize, sys, Uvec , ffreq, iE )
            enddo
         endif
      enddo !k-loop
   enddo

   ! remove psi_R 
   do ik = 1, nk 
      deallocate( psi_R%wk(ik)%cg ) 
   enddo
 
   print*, 'INVERT P MATRIX (Prrp -> Pggp)'

   ! dimension of the polarizability matrix: NFFT-by-NFFT
   nfft = FFTsize(1)*FFTsize(2)*FFTsize(3)

   do iq = 1, nq
      call P_r_to_g( iq, Prrp%mtrx(iq), FFTsize )
      Prrp%mtrx(iq)%C = Prrp%mtrx(iq)%C *sys%vol**2 / nfft**2
      
      if ( ffreq%is_on ) then
         do iE = 1, ffreq%nEstep
            call P_r_to_g( iq, P_freq%mtrx(iE,iq), FFTsize)
            P_freq%mtrx(iE,iq)%C = P_freq%mtrx(iE,iq)%C *sys%vol**2 / nfft**2
         enddo
      endif
     
   enddo


!*****************************************************************  
!           Epsilon and the inverse matrix calculations
      
   ! allocate for epsilon
   allocate( eps%mtrx(nq), eps%gs(nq) )
   
   if (ffreq%is_on) allocate( eps_freq%mtrx(ffreq%nEstep,nq) )
   
   ! find the dimension of epsilon matrix
   allocate( g_pol%gvec(3,nfft) , gidx(nfft) )

   g_pol%ng = nfft
   ! get polarization g index
   call fftidx_to_gidx( FFTsize, nfft, g_pol%gvec )



   ! Set q vectors in "Cartesian" coordinates

   deallocate( q%vec )
   ! k (from QE) is in cartesian coordinates
   call get_qvec( k, q )

   ! Set the q vector when q -> 0
   if (inp%crystal) then
      q%vec(1:3,1) = k%vec(1:3,1) - k_shft%vec(1:3,1)
   endif


   do iq = 1, nq
   
      call reduce_gspace( g_pol, eps%gs(iq), inp, sys, q%vec(1:3,iq), gidx )

      ndata = eps%gs(iq)%ng
      
      print*, ' '
      print*, 'number of data in epsilon calculations:', ndata

      allocate( eps%mtrx(iq)%C(ndata,ndata) )
      call calc_eps( Prrp%mtrx(iq), eps%mtrx(iq), ndata, nq, q%vec(:,iq), &
                     eps%gs(iq)%gvec, gidx, sys ) 
 
      deallocate( Prrp%mtrx(iq)%C ) 
      
      !call inverse( ndata, eps%mtrx(iq)%C )      
      
      call inverse_iterative( ndata, eps%mtrx(iq)%C )
 
      ! For the full-frequency calculations:
      if (ffreq%is_on) then
         do iE = 1, ffreq%nEstep
            allocate( eps_freq%mtrx(iE,iq)%C(ndata,ndata) )
            call calc_eps( P_freq%mtrx(iE,iq), eps_freq%mtrx(iE,iq), ndata, nq, &
                           q%vec(:,iq), eps%gs(iq)%gvec, gidx, sys )
            call inverse( ndata, eps_freq%mtrx(iE,iq)%C )               
         enddo
         deallocate( P_freq%mtrx )
      endif
      
   enddo




!*****************************************************************  
!   STEP 2: Generalized Plasmon-Pole calculations
!*****************************************************************  

!***************************************************************** 
!            Generalized Plasmon-Pole calculation 

   ! Epsilon inverse matrix : eps%mtrx(iq)%C

   call set_gpp_frequency( gpp ) 

   allocate( GPP%eigvec(nq), GPP%eigval(nq), GPP%ng(nq), GPP%omsq(nq) )
   allocate( GPP%mtrx(gpp%nEstep,nq) )

   
   do iq = 1, nq

      ndata = eps%gs(iq)%ng
      GPP%ng(iq) = ndata

      allocate( eigval_gpp( ndata ) )      
      allocate( Vcoulb( ndata ) )

      call calc_coulb( sys, ndata, q%vec(:,iq), eps%gs(iq)%gvec, Vcoulb, nq )

      ! eps%mtrx is going to be Smtrx
      call calc_Smtrx( ndata, eps%mtrx(iq), Vcoulb )  

      
      
      allocate( GPP%eigvec(iq)%C(ndata,ndata), GPP%eigval(iq)%R(ndata) )
      allocate( GPP%omsq(iq)%R(ndata) )
      
      
      GPP%eigvec(iq)%C(1:ndata,1:ndata) = eps%mtrx(iq)%C(1:ndata,1:ndata)

      
      ! remove eps matrix from the memory -> maybe you don't have to do this. Change your
      ! data structure. And keep it from epsilon calculations
      deallocate( eps%mtrx(iq)%C )

      call eigen_decomposition( ndata, GPP%eigvec(iq), GPP%eigval(iq) )
  
      call calc_Omsq( inp, sys, ndata, GPP%eigvec(iq), GPP%eigval(iq), GPP%omsq(iq), &
                       eps%gs(iq)%gvec, q%vec(1:3,iq) )


      
      ! calculate frequency dependent Smtrx(iq)
      do iE = 1, gpp%nEstep
         Eeval = CMPLX(( gpp%Emin + (iE-1)*gpp%Estep ) / Hartree, 0.2/Hartree )
      
         allocate( gpp%mtrx(iE,iq)%C(ndata,ndata) )

         call calc_Smtrx_w( ndata, GPP%eigval(iq)%R, GPP%omsq(iq)%R, Eeval, &
                            GPP%eigvec(iq)%C, GPP%mtrx(iE,iq)%C )

         ! calculate epsilon inverse
         call get_epsinv_mtrx_from_gpp( ndata, Vcoulb, gpp%mtrx(iE,iq)%C )
        
      enddo! END iE loop
      
      deallocate( Vcoulb, eigval_gpp )
      
      
   enddo  ! END iq loop
 
 
   



end program
