! Creating G spaces

subroutine reduce_gspace( G_pol, G_eps, inp, sys, qvec, gidx )

   use constant
   use electronic_structure
   use usrinput
   implicit none
   
   type(gspace), intent(inout) :: G_pol, G_eps
   type(input), intent(in) :: inp
   type(sysinfo), intent(in) :: sys
   real(dp) :: qvec(3) ! this q vector is in 2pi/alat unit (cartesian)
   integer, intent(inout) :: gidx(G_pol%ng)
   
   !work variables
   real(dp) :: Ecut ! energy cutoff for epsilon matrix calculations
   integer :: ig, ng , i, j
   real(dp) :: ekin ! kinetic energy ( 1/2*(q+G)^2 )
   real(dp) :: tmpvec(3), Gcart(3)
   integer :: counter
   
   
   Ecut = inp%epscut ! in Ry unit
   
   ng = G_pol%ng ! number of G vectors in polarizability calculations

   gidx = 0; counter = 0
      
   do ig = 1, ng

      ekin = 0.d0
      tmpvec = 0.d0
      Gcart = 0.d0
      do i = 1, 3
         do j = 1, 3
            Gcart(i) = G_pol%gvec(j,ig) * sys%bvec(i,j) + Gcart(i)
         enddo
      enddo

! q in cartesian unit at this moment
      
      tmpvec(1:3) = 2*pi/sys%alat * (Gcart(1:3)+qvec(1:3))
      ekin = 0.5d0 * dot_product( tmpvec(1:3), tmpvec(1:3) )
      
      if ( ekin .le. Ecut ) then
         counter = counter + 1
         gidx( counter ) = ig
      endif
      
   enddo
   
   
   G_eps%ng = counter
   allocate( G_eps%gvec(3,counter) ) 
   
   do i = 1, counter
      G_eps%gvec(:,i) = G_pol%gvec(:,gidx(i) )
   enddo


   
end subroutine
