module constant

   implicit none
   
   ! double precision
   integer, parameter :: dp = kind(1.0d0)
   integer, parameter :: dpc = kind((1.0d0,1.0d0))
   
   ! mathematical constants
   real(dp), parameter :: pi = 3.1415926535897932d0
   real(dp), parameter :: ryd = 13.605826d0
   real(dp), parameter :: bohr = 0.529177d0

   ! integers
   integer, parameter :: zero = 0
   integer, parameter :: one = 1
   integer, parameter :: minus_one = -1
   integer, parameter :: two = 2
   integer, parameter :: minus_two = -2

end module
