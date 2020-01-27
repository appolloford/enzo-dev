  subroutine evaluate_tgas(d, e, ge, u, v, w, &
    krome_x,imethod,idual,idim,tgas,utem)
    use krome_user

    !     written by: KROME DEVELOPERS
    !     date: 2013
    !     
    !     PURPOSE:
    !     Evaluate temperature Tgas from the internal energy e
    !     or ge depending on the energy flag idual

    implicit none
    real*8::krome_x(krome_nmols),p2d
    real*8::d,ge,e,tgas,kmu,kgamma,krome_tiny
    real*8::u,v,w,temstart,utem,ntot
    integer::imethod,idim,idual

    krome_tiny = 1d-30
    kgamma = krome_get_gamma(krome_x(:), 3d2) !3d2 is a dummy
    kmu    = krome_get_mu(krome_x(:))         !evaluate mean molecular weight

    !Compute Pressure with various methods 
    if (imethod .eq. 2) then
       !Zeus - e() is really gas energy
       p2d = (kgamma - 1.d0)*e
    else
       if(idual.eq.1) then
          !PPM with dual energy -- use gas energy
          p2d = (kgamma - 1.d0)*ge
       else
          !PPM without dual energy -- use total energy
          p2d = e - 0.5d0*u**2
          if (idim.gt.1) p2d = p2d - 0.5d0 * v**2
          if (idim.gt.2) p2d = p2d - 0.5d0 * w**2
          p2d = max((kgamma - 1.d0) * p2d, krome_tiny)
       endif
    endif

    !compute temperature
    tgas = p2d * utem * kmu

  end subroutine evaluate_tgas
