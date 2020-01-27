subroutine krome_initab()
   use krome_main
   use krome_user
   implicit none

   call krome_init()
   call krome_set_user_Av(2.9644d0)                                                    
   call krome_set_user_dgomega(5.d-1)
   call krome_set_user_zeta(1.d0)
   call krome_set_user_rad(1.d0)
   call krome_set_user_gRad(1.d-5)
   call krome_set_user_gArea(2.4d-22)
   call krome_set_user_fr(1.d0)
   call krome_set_user_desorb(1.d0)
   call krome_set_user_h2desorb(1.d0)
   call krome_set_user_crdesorb(1.d0)
   call krome_set_user_uvcr(1.d0)

end subroutine krome_initab
