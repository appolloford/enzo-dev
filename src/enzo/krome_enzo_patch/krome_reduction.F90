
!############### MODULE ##############
module krome_reduction
contains

  !**************************
  function fex_check(n,Tgas)
    use krome_commons
    use krome_tabs
    implicit none
    integer::i
    integer::r1,r2
    real*8::fex_check,n(nspec),k(nrea),rrmax,Tgas

    k(:) = coe_tab(n(:))
    rrmax = 0.d0
    n(idx_dummy) = 1.d0
    n(idx_g) = 1.d0
    n(idx_CR) = 1.d0
    do i=1,nrea
      r1 = arr_r1(i)
      r2 = arr_r2(i)
      arr_flux(i) = k(i)*n(r1)*n(r2)
      rrmax = max(rrmax, arr_flux(i))
    end do
    fex_check = rrmax

  end function fex_check

end module krome_reduction
