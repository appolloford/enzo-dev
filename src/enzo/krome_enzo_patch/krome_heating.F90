
!############### MODULE ##############
module KROME_heating
contains

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

  !************************
  function heating(n,inTgas,k,nH2dust)
    implicit none
    real*8::n(:), Tgas, inTgas, k(:), nH2dust
    real*8::heating

    Tgas = inTgas
    heating = sum(get_heating_array(n(:),Tgas,k(:), nH2dust))

  end function heating

  !*******************************
  function get_heating_array(n, Tgas, k, nH2dust)
    use krome_commons
    implicit none
    real*8::n(:), Tgas, k(:), nH2dust
    real*8::get_heating_array(nheats),heats(nheats)
    real*8::smooth,f1,f2
    !returns heating in erg/cm3/s

    heats(:) = 0.d0

    f2 = 1.

    heats(idx_heat_custom) = heat_custom(n(:),Tgas)

    get_heating_array(:) = heats(:)

  end function get_heating_array

  !*************************
  function heat_custom(n,Tgas)
    use krome_commons
    use krome_subs
    use krome_constants
    implicit none
    real*8::n(:),Tgas,heat_custom

    heat_custom = 0d0

  end function heat_custom

end module KROME_heating
