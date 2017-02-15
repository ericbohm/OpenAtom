!
!  Feb 2015 subhasish.mandal@yale.edu

!subroutine calc_S( sggp, epsinv,  ndata, nq, qvec, gvec, gidx, sys )
!subroutine calc_S( S, ndata, nq, qvec, G_vec, gidx, sys )
!subroutine calc_S( S, ndata, nq, qvec, gvec, gidx, sys )

subroutine calc_S( S,Eps_BGW, ngdata, nq, iq, qvec, gidx, sys )
!subroutine calc_S(  ndata, nq, qvec, gvec, gidx, sys )

   use constant
   use usrinput
   use electronic_structure
   use gw_structure
   
   implicit none
   
   type(rank2_mtrx), intent(inout) ::S
   type(rank2_mtrx) ::  EpsInv
   !type(rank2_mtrx) ::  Eps_BGW
   type(rank2_mtrx), intent(out) ::  Eps_BGW
   !type(epsilon_type) :: eps_inv
   type(rank2_mtrx)  ::  E
   !type(rank2_mtrx), intent(inout) ::  sggp
   !type(rank2_mtrx) ::  sggp
   type(sysinfo), intent(in) :: sys
   !integer :: ngdata  ! number of fft points = total number of g points in P matrix
   integer :: ng  ! number of fft points = total number of g points in P matrix
   !integer, intent(in) :: ndata  ! number of fft points = total number of g points in P matrix
   integer, intent(in) :: ngdata  ! total number of q points in the system
   integer, intent(in) :: nq  ! total number of q points in the system
   integer, intent(in) :: iq  ! total number of q points in the system
   real(dp), intent(in) :: qvec(3) ! q vector
   integer, allocatable :: gvec(:,:) ! g vectors
   integer, allocatable  :: G_vec(:, :) ! g vectors
   !integer, allocatable :: Gp_vec(:, :) ! g vectors
   integer, intent(in) :: gidx(ngdata)
   integer :: i,j
   integer :: ierr

   
   ! work variables
   real(dp), dimension(ngdata) :: coulb
   real(dp) :: qplusg
   real(dp) :: vol
   complex(dp) :: Eggp
   
   integer :: g, gp

  allocate(gvec(3, ngdata))


   coulb=0.0d0

   !allocate (G_vec(3,ndata*ndata))

   allocate(EpsInv%C(ngdata,ngdata))
   allocate(Eps_BGW%C(ngdata,ngdata))
   !allocate(coulb(ndata))
   !allocate(eps_inv%mtrx(nq))
   !*** calculate coulomb potential at q+G
   ! call calc_coulb( sys, ndata, qvec, gvec(:,:), coulb, nq )
   
   ! to make symmetric epsilon matrix, 
   !coulb(:) = sqrt( coulb(:) ) 

   !do i =1, ndata
   !print*,"checking Coulomb", coulb(i), ndata, gvec(:,i)
   !enddo
!goto 120
 !  call read_eps(Gvec, Gpvec,eps_inv%mtrx)
   !print*, 'read g-vector from eps_inv matrix',  Gvec(:,:)

!ng=411
   !call calc_coulb( sys, ndata, qvec, gvec(:, g), coulb, nq )
   !call calc_coulb( sys, ndata, qvec, gvec(:, gp), coulb, nq )
S%C=0.d0
EPS_BGW%C=0.d0
!do iq =1, nq

  print*, ' reading EPS_INV for iq', iq 
print*, 'ndata for iq', ngdata


if (iq.eq.1) open(18,file='EPS_INV_1.dat',form='formatted',status='old',iostat=ierr)
if (iq.eq.2) open(18,file='EPS_INV_2.dat',form='formatted',status='old',iostat=ierr)
if (iq.eq.3) open(18,file='EPS_INV_3.dat',form='formatted',status='old',iostat=ierr)
if (iq.eq.4) open(18,file='EPS_INV_4.dat',form='formatted',status='old',iostat=ierr)
if (iq.eq.5) open(18,file='EPS_INV_5.dat',form='formatted',status='old',iostat=ierr)
if (iq.eq.6) open(18,file='EPS_INV_6.dat',form='formatted',status='old',iostat=ierr)
if (iq.eq.7) open(18,file='EPS_INV_7.dat',form='formatted',status='old',iostat=ierr)
if (iq.eq.8) open(18,file='EPS_INV_8.dat',form='formatted',status='old',iostat=ierr)
!open(unit=18,file='EPS_INV_allq.dat',form='formatted',status='old',iostat=ierr)
!open(1,file='Check-Coulmb-Check-EPS.out',form='formatted',status='unknown')
if (iq.eq.1) open(1,file='Check-Coulmb-Check-EPS_1.out',form='formatted',status='unknown')
if (iq.eq.2) open(1,file='Check-Coulmb-Check-EPS_2.out',form='formatted',status='unknown')
if (iq.eq.3) open(1,file='Check-Coulmb-Check-EPS_3.out',form='formatted',status='unknown')
if (iq.eq.4) open(1,file='Check-Coulmb-Check-EPS_4.out',form='formatted',status='unknown')
if (iq.eq.5) open(1,file='Check-Coulmb-Check-EPS_5.out',form='formatted',status='unknown')
if (iq.eq.6) open(1,file='Check-Coulmb-Check-EPS_6.out',form='formatted',status='unknown')
if (iq.eq.7) open(1,file='Check-Coulmb-Check-EPS_7.out',form='formatted',status='unknown')
if (iq.eq.8) open(1,file='Check-Coulmb-Check-EPS_8.out',form='formatted',status='unknown')
write(1,*) ' g       gp        coulb(g)      coulb(gp)     &
Re(EPS_inv) -1/0  Im(EPS_inv)-1/0    --Re(S)--  -Im(S)-- '
   do gp = 1, ngdata
   coulb=0.d0
      do g = 1, ngdata
   coulb=0.d0
!
  !read(18, '(6i3,2f15.8)' ) G_vec(:, g), G_vec(:,gp) , EpsInv%C(g, gp)
      read(18, '(6i3,2f15.8)' ) gvec(:, g), gvec(:,gp) , EpsInv%C(g, gp)
  !print*, ' reading EPS_INV for iq', iq 
  !check Coulomb',  gvec(:,g), gvec(:,gp) ,coulb(g) , coulb(gp)
   call calc_coulb( sys, ngdata, qvec, gvec(:,: ), coulb, nq )
   !call calc_coulb( sys, ndata, qvec, gvec(:, gp), coulb, nq )
   !coulb(:) = sqrt( coulb(:) ) 
    !coulb(1)=0.38343474
     if ( iq .eq. 1 ) coulb(1)=0.19171737
     !if ( iq .eq. 1 ) coulb(1)=0.38343474
    !coulb(1)=0.61922
       S%C( g, gp ) =sqrt( coulb(g)) * EpsInv%C(g, gp) * sqrt(coulb(gp))
       EPS_BGW%C( g, gp ) = EpsInv%C(g, gp)
       !if ( g .eq. gp ) S%C(g,gp) = (EpsInv%C(g,gp) - 1)*sqrt(coulb(g))*sqrt(coulb(gp))
       if ( g .eq. gp ) S%C(g,gp) = sqrt(coulb(g))*(EpsInv%C(g,gp) - 1)*sqrt(coulb(gp))
       !if ( g .eq. gp ) EPS_BGW%C(g,gp) = EpsInv%C(g,gp) - 1
!  print*, 'read g-vector from eps_inv matrix',  G_vec(:,g), G_vec(:,gp) ,coulb(g) , coulb(gp)
  write(1,'(6i3,6f15.8 )') gvec(:,g), gvec(:,gp),coulb(g) ,coulb(gp), &
   Eps_BGW%C(g, gp),  S%C( g, gp )
  !write(1,'(6i3,8f15.8 )') gvec(:,g), gvec(:,gp),coulb(g) ,coulb(gp),EpsInv%C(g, gp) &
 !, 1 - EpsInv%C(g, gp), real(S%C(g,gp))
!  print*, ' check Coulomb',  gvec(:,g), gvec(:,gp) , S%C(g, gp)
      enddo
   enddo




! elseif (iq .ne. 1 ) then
!
!   do gp = ndata*(iq-1), ndata*iq
!   coulb=0.d0
!      do g = ndata*(iq-1), ndata*iq
!   coulb=0.d0
!
  print*, ' reading EPS_INV for iq', iq 

!E%C=0.d0
!S%C=0.0d0
close(18)
close(1)
!enddo !iq ends
!120 continue
end subroutine






subroutine calc_coulb( sys, ngdata, qvec, gvec, coulb, nq )
   
   use constant
   use electronic_structure
   implicit none
   type( sysinfo ), intent(in) :: sys
   integer, intent(in) :: ngdata            ! number of g vectors = Nfft
   real(dp), intent(in) :: qvec(3)         ! q vectors
   integer, intent(in) :: gvec(3,ngdata)    ! g vectors
   real(dp), intent(inout) :: coulb(ngdata) ! coulomb operator
   integer, intent(in) :: nq               ! number of q points (=# k points)
   
   ! work variables
   real(dp) :: avec(3,3), bvec(3,3), vol, alat
   integer :: i, j
   real(dp) :: gqsq, avectmp(3,3)
   real(dp) :: gq(3)

   avec = sys%avec
   bvec = sys%bvec
   vol = sys%vol
   alat = sys%alat
   coulb = 0.0d0

   do i = 1, ngdata
      ! change g vector coordinates
      call cryst_to_cart( bvec, alat, qvec, gvec(:,i), gq )
      
  !    print*, "check g-vec", gvec(:,i)
      gqsq = dot_product( gq, gq )
      
      coulb(i) = ( 4.d0 * pi ) /( gqsq *vol*nq)
      !!coulb(i) = ( 8.d0 * pi ) / (vol * nq * gqsq)
      ! If you want to compare the results with BGW, 
      ! Vc = 8*PI / (G+q)^2 should be used
      
   enddo

end subroutine
  

subroutine cryst_to_cart( bvec, alat, qvec, gvec, gplusqcart )
   use constant
   real(dp), intent(in) :: bvec(3,3)
   real(dp), intent(in) :: alat
   real(dp), intent(in) :: qvec(3) ! q is in cartesian 
   integer, intent(in) :: gvec(3)
   real(dp), intent(inout) :: GplusQcart(3)

   real(dp) :: gcryst(3), gcart(3)

   gcryst(1:3) = dble( gvec(1:3) )
   
   ! transfer to cartesian coordinates
      
   gcart = MATMUL( bvec, gcryst )
      
   GplusQcart = qvec + gcart
   GplusQcart = 2.d0 * pi / alat * GplusQcart
end subroutine cryst_to_cart




