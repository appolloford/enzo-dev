
!############### MODULE ##############
module krome_tabs
contains

  !***********************+
  function coe_tab(n)
    !interface to tabs
    use krome_subs
    use krome_getphys
    use krome_phfuncs
    use krome_grfuncs
    use krome_constants
    use krome_commons
    use krome_user_commons
    implicit none
    integer::idx,j
    real*8::Tgas, coe_tab(nrea),n(nspec),small

    Tgas = max(n(idx_Tgas),phys_Tcmb)
    small = 0d0

    coe_tab(:) = coe(n(:))

    coe_tab(999)=knrco(n(:))

  end function coe_tab

end module krome_tabs
