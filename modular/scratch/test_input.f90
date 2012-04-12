program test_input

   use initialization

   implicit none

   call initialize_params()

   print*, 'top of domain is: ', pm
   print*, 'bottom of domain is: ', pb
   print*, 'theta range is: ', qm

end program test_input
