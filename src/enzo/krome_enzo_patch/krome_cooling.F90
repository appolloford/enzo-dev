
!############### MODULE ##############
module KROME_cooling
  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2019-12-02 01:47:36
  !  Changeset xxxxxxx
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************
  integer,parameter::coolTab_n=int(1e2)
  integer,parameter::nZrate=0
  real*8::coolTab(nZrate,coolTab_n),coolTab_logTlow, coolTab_logTup
  real*8::coolTab_T(coolTab_n),inv_coolTab_T(coolTab_n-1),inv_coolTab_idx
contains

  !*******************
  function cooling(n,inTgas)
    use krome_commons
    implicit none
    real*8::n(:),inTgas,cooling,Tgas

    Tgas = inTgas
    cooling = sum(get_cooling_array(n(:),Tgas))

  end function cooling

  !*******************************
  function get_cooling_array(n, Tgas)
    use krome_commons
    implicit none
    real*8::n(:), Tgas
    real*8::get_cooling_array(ncools),cools(ncools)
    real*8::f1,f2,smooth

    f1 = 1d0
    f2 = 1d0

    !returns cooling in erg/cm3/s
    cools(:) = 0d0

    cools(idx_cool_custom) = cooling_custom(n(:),Tgas)

    get_cooling_array(:) = cools(:)

  end function get_cooling_array

  !*****************************
  function cooling_custom(n,Tgas)
    use krome_commons
    use krome_subs
    use krome_constants
    implicit none
    real*8::n(:),Tgas,cooling_custom

    cooling_custom = 0d0

  end function cooling_custom

  !**********************************
  function kpla(n,Tgas)
    !Planck opacity mean fit (Lenzuni+1996)
    !only denisity dependent (note that the
    ! fit provided by Lenzuni is wrong)
    ! valid for T<3e3 K
    !use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    real*8::kpla,rhogas,Tgas,n(:),y
    real*8::a0,a1,m(nspec)

    m(:) = get_mass()
    rhogas = sum(n(1:nmols)*m(1:nmols)) !g/cm3

    kpla = 0.d0
    !opacity is zero under 1e-12 g/cm3
    if(rhogas<1d-12) return

    !fit coefficients
    a0 = 1.000042d0
    a1 = 2.14989d0

    !log density cannot exceed 0.5 g/cm3
    y = log10(min(rhogas,0.5d0))

    kpla = 1d1**(a0*y + a1) !fit density only

  end function kpla

  !*****************************
  function coolingChem(n,Tgas)
    implicit none
    real*8::coolingChem,n(:),Tgas

    !note that this function is a dummy.
    ! For chemical cooling you should see
    ! heatingChem function in krome_heating.f90

    coolingChem = 0.d0

  end function coolingChem

  !***********************
  subroutine mylin2(a,b)
    !solve Ax=B analytically for a 2-levels system
    implicit none
    integer,parameter::n=2
    real*8::a(n,n),b(n),c(n),iab

    !uncomment this: safer but slower function
    !if(a(2,2)==a(2,1)) then
    !   print *,"ERROR: a22=a21 in mylin2"
    !   stop
    !end if
    iab = b(1)/(a(2,2)-a(2,1))
    c(1) = a(2,2) * iab
    c(2) = -a(2,1) * iab
    b(:) = c(:)

  end subroutine mylin2

  !************************
  subroutine mylin3(a,b)
    !solve Ax=B analytically for a 3-levels system
    implicit none
    integer,parameter::n=3
    real*8::iab,a(n,n),b(n),c(n)

    !uncomment this: safer but slower function
    !if(a(2,2)==a(2,3)) then
    !   print *,"ERROR: a22=a23 in mylin3"
    !   stop
    !end if

    !uncomment this: safer but slower
    !if(a(2,1)*a(3,2)+a(2,2)*a(3,3)+a(2,3)*a(3,1) == &
        !     a(2,1)*a(3,3)+a(2,2)*a(3,1)+a(2,3)*a(3,2)) then
    !   print *,"ERROR: division by zero in mylin3"
    !   stop
    !end if

    iab = b(1) / (a(2,1)*(a(3,3)-a(3,2)) + a(2,2)*(a(3,1)-a(3,3)) &
        + a(2,3)*(a(3,2)-a(3,1)))
    c(1) = (a(2,3)*a(3,2)-a(2,2)*a(3,3)) * iab
    c(2) = -(a(2,3)*a(3,1)-a(2,1)*a(3,3)) * iab
    c(3) = (a(3,1)*a(2,2)-a(2,1)*a(3,2)) * iab
    b(:) = c(:)

  end subroutine mylin3

  !************************************
  subroutine plot_cool(n)
    !routine to plot cooling at runtime
    real*8::n(:),Tgas,Tmin,Tmax
    real*8::cool_atomic,cool_H2,cool_HD,cool_tot, cool_totGP,cool_H2GP
    real*8::cool_dH,cool_Z
    integer::i,imax
    imax = 1000
    Tmin = log10(1d1)
    Tmax = log10(1d8)
    print *,"plotting cooling..."
    open(33,file="KROME_cooling_plot.dat",status="replace")
    do i=1,imax
      Tgas = 1d1**(i*(Tmax-Tmin)/imax+Tmin)
      cool_H2 = 0.d0
      cool_H2GP = 0.d0
      cool_HD = 0.d0
      cool_atomic = 0.d0
      cool_Z = 0.d0
      cool_dH = 0.d0
      cool_tot = cool_H2 + cool_atomic + cool_HD + cool_Z + cool_dH
      cool_totGP = cool_H2GP + cool_atomic + cool_HD + cool_Z + cool_dH
      write(33,'(99E12.3e3)') Tgas, cool_tot, cool_totGP, cool_H2, &
          cool_atomic, cool_HD, cool_H2GP, cool_Z, cool_dH
    end do
    close(33)
    print *,"done!"

  end subroutine plot_cool

  !***********************************
  !routine to dump cooling in unit nfile
  subroutine dump_cool(n,Tgas,nfile)
    use krome_commons
    implicit none
    real*8::Tgas,n(:),cools(ncools)
    integer::nfile

    cools(:) = get_cooling_array(n(:),Tgas)
    write(nfile,'(99E14.5e3)') Tgas, sum(cools), cools(:)

  end subroutine dump_cool

end module KROME_cooling
