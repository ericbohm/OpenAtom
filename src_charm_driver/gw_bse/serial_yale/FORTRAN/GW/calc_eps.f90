!
!  calculate epsilon matrix
!  Oct 2014. minjung.kim@yale.edu
!

subroutine calc_eps( P, eps, ndata, nq, qvec, gvec, gidx, sys )

   use constant
   use electronic_structure
   use gw_structure
   
   implicit none
   
   type(rank2_mtrx), intent(inout) :: P, eps
   type(sysinfo), intent(in) :: sys
   integer, intent(in) :: ndata  ! number of fft points = total number of g points in P matrix
   integer, intent(in) :: nq  ! total number of q points in the system
   real(dp), intent(in) :: qvec(3) ! q vector
   integer, intent(in) :: gvec(3,ndata) ! g vectors
   integer, intent(in) :: gidx(ndata)

   
   ! work variables
   real(dp), dimension(ndata) :: coulb
   real(dp) :: qplusg
   real(dp) :: vol
   complex(dp) :: Pggp
   
   integer :: g, gp

   !*** calculate coulomb potential at q+G
   call calc_coulb( sys, ndata, qvec, gvec, coulb, nq )
   
   ! to make symmetric epsilon matrix, 
   coulb(:) = sqrt( coulb(:) ) 

   do gp = 1, ndata
      do g = 1, ndata

         Pggp = P%C( gidx(g), gidx(gp) )

         eps%C( g, gp ) = -coulb(g) * Pggp * coulb(gp)

         if ( g .eq. gp ) eps%C(g,gp) = eps%C(g,gp) + 1

      enddo
   enddo


end subroutine






subroutine calc_coulb( sys, ndata, qvec, gvec, coulb, nq )
   
   use constant
   use electronic_structure
   implicit none
   type( sysinfo ), intent(in) :: sys
   integer, intent(in) :: ndata            ! number of g vectors = Nfft
   real(dp), intent(in) :: qvec(3)         ! q vectors
   integer, intent(in) :: gvec(3,ndata)    ! g vectors
   real(dp), intent(inout) :: coulb(ndata) ! coulomb operator
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

   do i = 1, ndata
      ! change g vector coordinates
      call cryst_to_cart( bvec, alat, qvec, gvec(:,i), gq )
      
      gqsq = dot_product( gq, gq )
      
      coulb(i) = ( 4.d0 * pi ) / (vol * nq * gqsq)
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




