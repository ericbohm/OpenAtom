!
!  calculate sigma matrix
!  Oct 2014. minjung.kim@yale.edu
!  July 2015. subhasish.mandal@yale.edu

subroutine calc_sig( D, sig, ndata, qvec, gvec , sys)
!subroutine calc_sig( D, sig, ndata, sys, nq, coulb, g )

   use constant
   use electronic_structure
   use gw_structure
   implicit none
   type(rank2_mtrx), intent(inout) :: D
   real(dp), intent(inout) ::  sig
   integer, intent(in) :: ndata
   integer ::  nq
   !real(dp), intent(in) :: qvec(3) ! q vector
   real(dp) :: qvec(3) ! q vector
   !integer, intent(in) :: gvec(3,ndata) ! g vectors
   integer :: gvec(3,ndata) ! g vectors
   type(sysinfo), intent(in) :: sys

   ! work variables
   real(dp), allocatable :: coulb(:)
   real(dp) :: qplusg
   real(dp) :: vol
   integer :: g

   !*** calculate coulomb potential at q+G
   allocate(coulb(ndata))
   !allocate(g(ndata))
   call calc_coulb( sys, ndata, qvec, gvec, coulb, nq )
   sig=0.0d0
   do g = 1, ndata

   !      Dgg = D%C( gidx(g), gidx(g) )
         !sig%C( g, g ) = -coulb(g) * Dgg + sig%C( g, g )
         sig = -coulb(g) * D%C(g,g) + sig
!         sig = -coulb(g)

!         print*, 'check coulb:', g
   enddo
   deallocate(coulb)
   !print*, 'sigma summed over all g:',g ,ndata, sig
!   print*, 'sigma summed over all g:', coulb
end subroutine






subroutine calc_coulb( sys, ndata, qvec, gvec, coulb, nq )
   
   use constant
   use electronic_structure
   implicit none
   type( sysinfo ), intent(in) :: sys
   integer, intent(in) :: ndata, nq            ! number of g vectors = Nfft
   real(dp), intent(in) :: qvec(3)         ! q vectors
   integer, intent(in) :: gvec(3,ndata)    ! g vectors
   real(dp), intent(inout) :: coulb(ndata) ! coulomb operator
   
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
   do i = 1, ndata
      ! change g vector coordinates
      call cryst_to_cart( bvec, alat, qvec, gvec(:,i), gq )
      
      gqsq = dot_product( gq, gq )

      !coulb(i) = ( 4.d0 * pi )/( vol * nq * qplusg )
      
      ! Vc = 8*PI / (G+q)^2
      coulb(i) = ( 8.d0 * pi ) /( gqsq *vol*nq)
     !To compare with BGW
      !coulb(i) = ( 8.d0 * pi ) / gqsq
!      coulb(i) = sqrt( coulb(i) )

   enddo
   !print*, 'Coulb', coulb
   contains
   subroutine cryst_to_cart( bvec, alat, qvec, gvec, gplusqcart )
      real(dp), intent(in) :: bvec(3,3)
      real(dp), intent(in) :: alat
      real(dp), intent(in) :: qvec(3) ! q is in cartesian 
      integer, intent(in) :: gvec(3)
      real(dp), intent(inout) :: GplusQcart(3)

      real(dp) :: gcryst(3), gcart(3)

      gcryst(1:3) = dble( gvec(1:3) ) + dble( qvec(1:3) )
   
      ! transfer to cartesian coordinates
      
      gcart = MATMUL( bvec, gcryst )
      
      GplusQcart = gcart
      GplusQcart = 2.d0 * pi / alat * GplusQcart
   end subroutine cryst_to_cart

end subroutine


