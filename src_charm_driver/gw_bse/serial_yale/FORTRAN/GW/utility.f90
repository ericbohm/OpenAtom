! utilities

subroutine print_error_exit(str)

   character(len=*), intent(in) :: str
   
   write(*,'(a)') trim(str)
   stop
end subroutine
