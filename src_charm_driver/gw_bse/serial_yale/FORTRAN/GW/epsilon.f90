!  Epsilon main program
!
!  Sep. 2014  minjung.kim@yale.edu
!

program epsilon

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

   ! determine FFT grid. This number is multiplied to the minimum cutoff radius
   integer :: FFTsize(3)
   integer :: nk, nq, iq, ik
   integer :: ikq ! kpt index for k+q vector
   integer :: nfft ! FFTsize(1)*FFTsize(2)*FFTsize(3)
   integer :: ndata ! number of g vector in epsilon matrix   
   integer, allocatable :: gidx(:) ! link G vectors between polarizability and epsilon
   integer :: Uvec(3) ! Umklapp vector
   
   integer:: i,j, idata, ii

!*****************************************
integer, allocatable :: gv(:,:), gpv(:,:)
integer :: ig, igp, test, kk
real(dp) :: g2, gp2
!*****************************************

   ! Read input values
   call read_input( inp )
      
   
   ! Read wavefunction data and save into the memory
   call read_wfn( inp%wfname, sys, psi, k )

  
   ! check_dim is for testing.
   ! call check_dim( psi )

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

   ! Polarizability matrix construction
   Prrp%nq = nq
   call create_Pmtrx_struc( Prrp, FFTsize )

 
 
!*****************************************************************  
!                  calculate Polarization matrix
   print*, 'Calculate Pmtrx'
   ! CASE1: q = 0 

   iq = 1
   Uvec(1:3) = 0
   if (inp%crystal) then
   
      call shifted_k_wfn( FFTsize, inp, psi_shft, psi_R_shft, k_shft, sys_shft )

      do ik = 1, nk

         call calc_Prrp( psi_R_shft%wk(ik), psi_R%wk(ik), Prrp%mtrx(iq), FFTsize , sys, Uvec )

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

      enddo
 
   enddo
   
   ! remove psi_R 
   do ik = 1, nk 
      deallocate( psi_R%wk(ik)%cg ) 
   enddo


   ! Converting Prrp to Pggp
   ! *** I am not going to make Pggp matrix 
   ! instead, Pggp is being stored in Prrp matrix

   print*, 'FFT (Prrp -> Pggp)'

   ! dimension of the polarizability matrix: NFFT-by-NFFT
   nfft = FFTsize(1)*FFTsize(2)*FFTsize(3)

   do iq = 1, nq
        
      call P_r_to_g( iq, Prrp%mtrx(iq), FFTsize )
      Prrp%mtrx(iq)%C = Prrp%mtrx(iq)%C * sys%vol**2 / nfft**2

   enddo



!*****************************************************************  
!            calculate Epsilon and inverse matrix 
      
   ! allocate for epsilon
   allocate( eps%mtrx(nq), eps%gs(nq) )
   
   ! find the dimension of epsilon matrix
   allocate( g_pol%gvec(3,nfft) , gidx(nfft) )

   g_pol%ng = nfft
   ! get polarization g index
   call fftidx_to_gidx( FFTsize, nfft, g_pol%gvec )


   ! Set q vectors in Cartesian coordinates

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
      print*, 'number of data in epsilon calculations:', ndata
      allocate( eps%mtrx(iq)%C(ndata,ndata) )

      print*, 'Calculate epsilon matrix at q:', q%vec(1:3,iq)
      call calc_eps( Prrp%mtrx(iq), eps%mtrx(iq), ndata, nq, q%vec(:,iq), &
                     eps%gs(iq)%gvec, gidx, sys )
                     
      ! remove Polarizability matrix from memory
      !deallocate( Prrp%mtrx(iq)%C )
      
      ! if you want to write epsilon matrix, this is the time to do that.
      ! eps%mtrx will be inverted in the next step
      print*, 'calculate epsilon inverse'
      
      
      !call inverse( ndata, eps%mtrx(iq)%C )
      
      call inverse_iterative( ndata, eps%mtrx(iq)%C )


!*****************************************
!print*, 'print out epsilon inverse matrix'
!if (iq.eq.1) open(190,file='EPS_INV_1.dat',form='formatted',status='unknown')
!if (iq.eq.2) open(190,file='EPS_INV_2.dat',form='formatted',status='unknown')


print*, "number of g points"
print*, eps%gs(1)%ng
print*, eps%gs(2)%ng

goto 130
write(190,*) '21609'
do j = 1,eps%gs(iq)%ng
   do i = 1, eps%gs(iq)%ng
      g2 = 0.d0
      gp2 = 0.d0
      do ii = 1, 3
         g2 = eps%gs(iq)%gvec(ii,i)**2 + g2
         gp2 = eps%gs(iq)%gvec(ii,j)**2 + gp2
      enddo
      if( g2 .le. 10.d0 .and. gp2 .le. 10.0d0 ) then 
         write(190,'(6i3,2f15.8)') eps%gs(iq)%gvec(:,i), eps%gs(iq)%gvec(:,j), eps%mtrx(iq)%C(i,j)
      endif 
   enddo
enddo
close(190)
!******************************************
130 continue

   enddo



end program
