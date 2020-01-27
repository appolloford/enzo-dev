
!############### MODULE ##############
module krome_subs
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
  !compute reaction rates cm^3(n-1)/s
  function coe(n)
    use krome_commons
    use krome_constants
    use krome_user_commons
    use krome_getphys
    use krome_grfuncs
    use krome_phfuncs
    use krome_fit
    implicit none
    real*8::coe(nrea),k(nrea),Tgas,n(nspec),kmax
    real*8::T32
    real*8::invT
    real*8::sqrTgas
    real*8::small,nmax
    integer::i
    real*8::mantle  !preproc from coevar
    real*8::Hnuclei  !preproc from coevar
    real*8::mantleabund  !preproc from coevar
    !Tgas is in K
    Tgas = max(n(idx_Tgas), phys_Tcmb)
    Tgas = min(Tgas,1d9)

    !maxn initialization can be removed and small can be
    ! replaced with a proper value according to the environment
    nmax = max(maxval(n(1:nmols)),1d0)
    small = 1d-40/(nmax*nmax*nmax*nmax)

    T32 = Tgas*0.0033333333333333335 !Tgas/(300 K) (#)
    invT = 1.d0/Tgas !inverse of T (1/K)
    sqrTgas = sqrt(Tgas) !Tgas rootsquare (K**0.5)

    Hnuclei = get_Hnuclei(n(:))
    mantle = get_mantle(n(:))
    mantleabund = mantle/Hnuclei

    k(:) = small !inizialize coefficients

    !CH + O -> HCO+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1) = small + (1.09e-11&
          *(T32)**(-2.19)*exp(-165.1*invT))
    end if

    !H+ + HNC -> HCN + H+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(2) = small + (1e-09)
    end if

    !H2 + CH -> C + H2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(3) = small + (6e-09&
          *exp(-40200.0*invT))
    end if

    !H2 + H2 -> H2 + H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(4) = small + (1e-08&
          *exp(-84100.0*invT))
    end if

    !H2 + H2O -> OH + H2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(5) = small + (5.8e-09&
          *exp(-52900.0*invT))
    end if

    !H2 + HOC+ -> HCO+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(6) = small + (3.8e-10)
    end if

    !H2 + O2 -> O + O + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(7) = small + (6e-09&
          *exp(-52300.0*invT))
    end if

    !H2 + OH -> O + H2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(8) = small + (6e-09&
          *exp(-50900.0*invT))
    end if

    !H2 + E -> H + H + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(9) = small + (3.22e-09&
          *(T32)**(0.35)*exp(-102000.0*invT))
    end if

    !H + CH -> C + H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(10) = small + (6e-09&
          *exp(-40200.0*invT))
    end if

    !H + H2 -> H + H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(11) = small + (4.67e-07&
          *(T32)**(-1.0)*exp(-55000.0*invT))
    end if

    !H + H2O -> OH + H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(12) = small + (5.8e-09&
          *exp(-52900.0*invT))
    end if

    !H + O2 -> O + O + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(13) = small + (6e-09&
          *exp(-52300.0*invT))
    end if

    !H + OH -> O + H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(14) = small + (6e-09&
          *exp(-50900.0*invT))
    end if

    !C+ + CH2 -> CH2+ + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(15) = small + (5.2e-10)
    end if

    !C+ + CH -> CH+ + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(16) = small + (3.8e-10&
          *(T32)**(-0.5))
    end if

    !C+ + H2CO -> H2CO+ + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(17) = small + (7.8e-10&
          *(T32)**(-0.5))
    end if

    !C+ + HCO -> HCO+ + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(18) = small + (4.8e-10&
          *(T32)**(-0.5))
    end if

    !C+ + MG -> MG+ + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(19) = small + (1.1e-09)
    end if

    !C+ + NH3 -> NH3+ + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(20) = small + (6.72e-10&
          *exp(+0.5*invT))
    end if

    !C+ + NO -> NO+ + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(21) = small + (7.05e-10&
          *(T32)**(-0.03)*exp(+16.7*invT))
    end if

    !C + CN+ -> CN + C+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(22) = small + (1.1e-10)
    end if

    !C + CO+ -> CO + C+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(23) = small + (1.1e-10)
    end if

    !C + N2+ -> N2 + C+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(24) = small + (1.1e-10)
    end if

    !C + O2+ -> O2 + C+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(25) = small + (5.2e-11)
    end if

    !CH+ + HCO -> HCO+ + CH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(26) = small + (4.6e-10&
          *(T32)**(-0.5))
    end if

    !CH+ + MG -> MG+ + CH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(27) = small + (3.6e-10)
    end if

    !CH+ + NH3 -> NH3+ + CH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(28) = small + (4.59e-10&
          *(T32)**(-0.5))
    end if

    !CH+ + NO -> NO+ + CH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(29) = small + (7.6e-10)
    end if

    !CH2+ + NO -> NO+ + CH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(30) = small + (4.2e-10)
    end if

    !CH2 + CN+ -> CN + CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(31) = small + (8.8e-10)
    end if

    !CH2 + CO+ -> CO + CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(32) = small + (4.3e-10)
    end if

    !CH2 + H2CO+ -> H2CO + CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(33) = small + (4.3e-10)
    end if

    !CH2 + H2O+ -> H2O + CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(34) = small + (4.7e-10)
    end if

    !CH2 + N2+ -> N2 + CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(35) = small + (8.7e-10)
    end if

    !CH2 + NH2+ -> NH2 + CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(36) = small + (4.9e-10)
    end if

    !CH2 + O+ -> O + CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(37) = small + (9.7e-10)
    end if

    !CH2 + O2+ -> O2 + CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(38) = small + (4.3e-10)
    end if

    !CH2 + OH+ -> OH + CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(39) = small + (4.8e-10)
    end if

    !CH3+ + HCO -> HCO+ + CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(40) = small + (4.4e-10&
          *(T32)**(-0.5))
    end if

    !CH3+ + MG -> MG+ + CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(41) = small + (3.5e-09)
    end if

    !CH3+ + NO -> NO+ + CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(42) = small + (1e-09)
    end if

    !CH4+ + H2CO -> H2CO+ + CH4
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(43) = small + (1.62e-09&
          *(T32)**(-0.5))
    end if

    !CH4+ + NH3 -> NH3+ + CH4
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(44) = small + (1.65e-09&
          *(T32)**(-0.5))
    end if

    !CH4+ + O2 -> O2+ + CH4
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(45) = small + (3.9e-10)
    end if

    !CH4 + CO+ -> CO + CH4+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(46) = small + (7.93e-10)
    end if

    !CH + CN+ -> CN + CH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(47) = small + (6.4e-10&
          *(T32)**(-0.5))
    end if

    !CH + CO+ -> CO + CH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(48) = small + (3.2e-10&
          *(T32)**(-0.5))
    end if

    !CH + H2CO+ -> H2CO + CH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(49) = small + (3.1e-10&
          *(T32)**(-0.5))
    end if

    !CH + H2O+ -> H2O + CH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(50) = small + (3.4e-10&
          *(T32)**(-0.5))
    end if

    !CH + N+ -> N + CH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(51) = small + (3.6e-10&
          *(T32)**(-0.5))
    end if

    !CH + N2+ -> N2 + CH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(52) = small + (6.3e-10&
          *(T32)**(-0.5))
    end if

    !CH + NH2+ -> NH2 + CH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(53) = small + (3.5e-10&
          *(T32)**(-0.5))
    end if

    !CH + O+ -> O + CH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(54) = small + (3.5e-10&
          *(T32)**(-0.5))
    end if

    !CH + O2+ -> O2 + CH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(55) = small + (3.1e-10&
          *(T32)**(-0.5))
    end if

    !CH + OH+ -> OH + CH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(56) = small + (3.5e-10&
          *(T32)**(-0.5))
    end if

    !CN+ + CO -> CO+ + CN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(57) = small + (6.3e-10)
    end if

    !CN+ + H2CO -> H2CO+ + CN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(58) = small + (5.2e-10&
          *(T32)**(-0.5))
    end if

    !CN+ + HCN -> HCN+ + CN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(59) = small + (1.79e-09&
          *(T32)**(-0.5))
    end if

    !CN+ + HCO -> HCO+ + CN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(60) = small + (3.7e-10&
          *(T32)**(-0.5))
    end if

    !CN+ + NO -> NO+ + CN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(61) = small + (5.7e-10)
    end if

    !CN+ + O2 -> O2+ + CN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(62) = small + (2.58e-10)
    end if

    !CN + N2+ -> N2 + CN+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(63) = small + (1e-10&
          *(T32)**(-0.5))
    end if

    !CO+ + H2CO -> H2CO+ + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(64) = small + (1.35e-09&
          *(T32)**(-0.5))
    end if

    !CO+ + HCO -> HCO+ + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(65) = small + (7.4e-10&
          *(T32)**(-0.5))
    end if

    !CO+ + NO -> NO+ + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(66) = small + (3.3e-10)
    end if

    !CO+ + O2 -> O2+ + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(67) = small + (1.2e-10)
    end if

    !CO + N2+ -> N2 + CO+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(68) = small + (7.4e-11)
    end if

    !H+ + CH2 -> CH2+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(69) = small + (1.4e-09)
    end if

    !H+ + CH3 -> CH3+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(70) = small + (3.4e-09)
    end if

    !H+ + CH4 -> CH4+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(71) = small + (1.5e-09)
    end if

    !H+ + CH -> CH+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(72) = small + (1.9e-09&
          *(T32)**(-0.5))
    end if

    !H+ + H2CO -> H2CO+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(73) = small + (2.96e-09&
          *(T32)**(-0.5))
    end if

    !H+ + H2O -> H2O+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(74) = small + (6.9e-09&
          *(T32)**(-0.5))
    end if

    !H+ + HCN -> HCN+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(75) = small + (1.05e-08&
          *(T32)**(-0.13))
    end if

    !H+ + HCO -> HCO+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(76) = small + (9.4e-10&
          *(T32)**(-0.5))
    end if

    !H+ + MG -> MG+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(77) = small + (1.1e-09)
    end if

    !H+ + NH2 -> NH2+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(78) = small + (2.9e-09&
          *(T32)**(-0.5))
    end if

    !H+ + NH3 -> NH3+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(79) = small + (3.7e-09&
          *(T32)**(-0.5))
    end if

    !H+ + NH -> NH+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(80) = small + (2.1e-09&
          *(T32)**(-0.5))
    end if

    !H+ + NO -> NO+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(81) = small + (2.9e-09)
    end if

    !H+ + O2 -> O2+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(82) = small + (2e-09)
    end if

    !H+ + O -> O+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(83) = small + (6.86e-10&
          *(T32)**(0.26)*exp(-224.3*invT))
    end if

    !H+ + OH -> OH+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(84) = small + (2.1e-09&
          *(T32)**(-0.5))
    end if

    !H2+ + CH2 -> CH2+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(85) = small + (1e-09)
    end if

    !H2+ + CH4 -> CH4+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(86) = small + (1.4e-09)
    end if

    !H2+ + CH -> CH+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(87) = small + (7.1e-10&
          *(T32)**(-0.5))
    end if

    !H2+ + CN -> CN+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(88) = small + (1.2e-09&
          *(T32)**(-0.5))
    end if

    !H2+ + CO -> CO+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(89) = small + (6.44e-10)
    end if

    !H2+ + H2CO -> H2CO+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(90) = small + (1.4e-09&
          *(T32)**(-0.5))
    end if

    !H2+ + H2O -> H2O+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(91) = small + (3.9e-09&
          *(T32)**(-0.5))
    end if

    !H2+ + HCN -> HCN+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(92) = small + (2.7e-09&
          *(T32)**(-0.5))
    end if

    !H2+ + HCO -> HCO+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(93) = small + (1e-09&
          *(T32)**(-0.5))
    end if

    !H2+ + NH2 -> NH2+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(94) = small + (2.1e-09&
          *(T32)**(-0.5))
    end if

    !H2+ + NH3 -> NH3+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(95) = small + (5.7e-09&
          *(T32)**(-0.5))
    end if

    !H2+ + NH -> NH+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(96) = small + (7.6e-10&
          *(T32)**(-0.5))
    end if

    !H2+ + NO -> NO+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(97) = small + (1.1e-09)
    end if

    !H2+ + O2 -> O2+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(98) = small + (8e-10)
    end if

    !H2+ + OH -> OH+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(99) = small + (7.6e-10&
          *(T32)**(-0.5))
    end if

    !H2 + HE+ -> HE + H2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(100) = small + (7.2e-15)
    end if

    !H2CO + O2+ -> O2 + H2CO+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(101) = small + (2.07e-09&
          *(T32)**(-0.5))
    end if

    !H2O+ + H2CO -> H2CO+ + H2O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(102) = small + (1.41e-09&
          *(T32)**(-0.5))
    end if

    !H2O+ + HCO -> HCO+ + H2O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(103) = small + (2.8e-10&
          *(T32)**(-0.5))
    end if

    !H2O+ + MG -> MG+ + H2O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(104) = small + (2.2e-09)
    end if

    !H2O+ + NO -> NO+ + H2O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(105) = small + (2.7e-10)
    end if

    !H2O+ + O2 -> O2+ + H2O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(106) = small + (4.6e-10)
    end if

    !H2O + CO+ -> CO + H2O+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(107) = small + (1.72e-09&
          *(T32)**(-0.5))
    end if

    !H2O + HCN+ -> HCN + H2O+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(108) = small + (1.8e-09&
          *(T32)**(-0.5))
    end if

    !H2O + N2+ -> N2 + H2O+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(109) = small + (2.3e-09&
          *(T32)**(-0.5))
    end if

    !H + CN+ -> CN + H+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(110) = small + (6.4e-10)
    end if

    !H + CO+ -> CO + H+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(111) = small + (7.5e-10)
    end if

    !H + H2+ -> H2 + H+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(112) = small + (6.4e-10)
    end if

    !H + HCN+ -> HCN + H+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(113) = small + (3.7e-11)
    end if

    !H + HE+ -> HE + H+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(114) = small + (1.2e-15&
          *(T32)**(0.25))
    end if

    !H + O+ -> O + H+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(115) = small + (5.66e-10&
          *(T32)**(0.36)*exp(+8.6*invT))
    end if

    !HCN+ + NO -> NO+ + HCN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(116) = small + (8.1e-10)
    end if

    !HCN+ + O2 -> O2+ + HCN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(117) = small + (3.2e-10)
    end if

    !HCN + CO+ -> CO + HCN+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(118) = small + (3.4e-09&
          *(T32)**(-0.5))
    end if

    !HCN + N2+ -> N2 + HCN+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(119) = small + (3.9e-10&
          *(T32)**(-0.5))
    end if

    !HCO + H2CO+ -> H2CO + HCO+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(120) = small + (3.6e-10&
          *(T32)**(-0.5))
    end if

    !HCO + O2+ -> O2 + HCO+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(121) = small + (3.6e-10&
          *(T32)**(-0.5))
    end if

    !HE+ + C -> C+ + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(122) = small + (6.3e-15&
          *(T32)**(0.75))
    end if

    !HE+ + CH4 -> CH4+ + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(123) = small + (5.1e-11)
    end if

    !HE+ + CH -> CH+ + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(124) = small + (5e-10&
          *(T32)**(-0.5))
    end if

    !HE+ + H2CO -> H2CO+ + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(125) = small + (9.69e-10&
          *(T32)**(-0.5))
    end if

    !HE+ + H2O -> H2O+ + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(126) = small + (6.05e-11&
          *(T32)**(-0.5))
    end if

    !HE+ + N2 -> N2+ + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(127) = small + (6.4e-10)
    end if

    !HE+ + NH3 -> NH3+ + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(128) = small + (2.64e-10&
          *(T32)**(-0.5))
    end if

    !HE+ + O2 -> O2+ + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(129) = small + (3.3e-11)
    end if

    !MG + H2CO+ -> H2CO + MG+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(130) = small + (2.9e-09)
    end if

    !MG + HCO+ -> HCO + MG+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(131) = small + (2.9e-09)
    end if

    !MG + N2+ -> N2 + MG+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(132) = small + (7e-10)
    end if

    !MG + NO+ -> NO + MG+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(133) = small + (8.1e-10)
    end if

    !MG + O2+ -> O2 + MG+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(134) = small + (1.2e-09)
    end if

    !N+ + CH2 -> CH2+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(135) = small + (1e-09)
    end if

    !N+ + CH4 -> CH4+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(136) = small + (2.8e-11)
    end if

    !N+ + CN -> CN+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(137) = small + (1.1e-09&
          *(T32)**(-0.5))
    end if

    !N+ + CO -> CO+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(138) = small + (8.25e-10)
    end if

    !N+ + H2CO -> H2CO+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(139) = small + (1.88e-09&
          *(T32)**(-0.5))
    end if

    !N+ + H2O -> H2O+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(140) = small + (2.8e-09&
          *(T32)**(-0.5))
    end if

    !N+ + HCN -> HCN+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(141) = small + (3.7e-09&
          *(T32)**(-0.5))
    end if

    !N+ + HCO -> HCO+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(142) = small + (4.5e-10&
          *(T32)**(-0.5))
    end if

    !N+ + MG -> MG+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(143) = small + (1.2e-09)
    end if

    !N+ + NH2 -> NH2+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(144) = small + (1e-09&
          *(T32)**(-0.5))
    end if

    !N+ + NH3 -> NH3+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(145) = small + (1.97e-09&
          *(T32)**(-0.5))
    end if

    !N+ + NH -> NH+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(146) = small + (3.7e-10&
          *(T32)**(-0.5))
    end if

    !N+ + NO -> NO+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(147) = small + (4.51e-10)
    end if

    !N+ + O2 -> O2+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(148) = small + (3.11e-10)
    end if

    !N+ + OH -> OH+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(149) = small + (3.7e-10&
          *(T32)**(-0.5))
    end if

    !N2+ + H2CO -> H2CO+ + N2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(150) = small + (3.77e-10&
          *(T32)**(-0.5))
    end if

    !N2+ + HCO -> HCO+ + N2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(151) = small + (3.7e-10&
          *(T32)**(-0.5))
    end if

    !N2+ + NO -> NO+ + N2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(152) = small + (4.4e-10)
    end if

    !N2+ + O2 -> O2+ + N2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(153) = small + (5e-11)
    end if

    !N + N2+ -> N2 + N+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(154) = small + (1e-11)
    end if

    !NH+ + H2CO -> H2CO+ + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(155) = small + (9.9e-10&
          *(T32)**(-0.5))
    end if

    !NH+ + H2O -> H2O+ + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(156) = small + (1.05e-09&
          *(T32)**(-0.5))
    end if

    !NH+ + NH3 -> NH3+ + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(157) = small + (1.8e-09&
          *(T32)**(-0.5))
    end if

    !NH+ + NO -> NO+ + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(158) = small + (7.12e-10)
    end if

    !NH+ + O2 -> O2+ + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(159) = small + (4.51e-10)
    end if

    !NH2+ + HCO -> HCO+ + NH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(160) = small + (4.3e-10&
          *(T32)**(-0.5))
    end if

    !NH2+ + NH3 -> NH3+ + NH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(161) = small + (6.9e-10&
          *(T32)**(-0.5))
    end if

    !NH2+ + NO -> NO+ + NH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(162) = small + (7e-10)
    end if

    !NH2 + CN+ -> CN + NH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(163) = small + (9.1e-10&
          *(T32)**(-0.5))
    end if

    !NH2 + CO+ -> CO + NH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(164) = small + (4.5e-10&
          *(T32)**(-0.5))
    end if

    !NH2 + H2O+ -> H2O + NH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(165) = small + (4.9e-10&
          *(T32)**(-0.5))
    end if

    !NH2 + N2+ -> N2 + NH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(166) = small + (8.9e-10&
          *(T32)**(-0.5))
    end if

    !NH2 + O2+ -> O2 + NH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(167) = small + (8.7e-10&
          *(T32)**(-0.5))
    end if

    !NH2 + OH+ -> OH + NH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(168) = small + (5e-10&
          *(T32)**(-0.5))
    end if

    !NH3+ + HCO -> HCO+ + NH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(169) = small + (4.2e-10&
          *(T32)**(-0.5))
    end if

    !NH3+ + MG -> MG+ + NH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(170) = small + (3.3e-09)
    end if

    !NH3+ + NO -> NO+ + NH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(171) = small + (7.2e-10)
    end if

    !NH3 + CO+ -> CO + NH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(172) = small + (2.02e-09&
          *(T32)**(-0.5))
    end if

    !NH3 + H2CO+ -> H2CO + NH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(173) = small + (4.25e-10&
          *(T32)**(-0.5))
    end if

    !NH3 + H2O+ -> H2O + NH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(174) = small + (2.21e-09&
          *(T32)**(-0.5))
    end if

    !NH3 + HCN+ -> HCN + NH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(175) = small + (1.68e-09&
          *(T32)**(-0.5))
    end if

    !NH3 + N2+ -> N2 + NH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(176) = small + (1.9e-09&
          *(T32)**(-0.5))
    end if

    !NH3 + O2+ -> O2 + NH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(177) = small + (2e-09&
          *(T32)**(-0.5))
    end if

    !NH + CN+ -> CN + NH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(178) = small + (6.5e-10&
          *(T32)**(-0.5))
    end if

    !NH + CO+ -> CO + NH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(179) = small + (3.2e-10&
          *(T32)**(-0.5))
    end if

    !NH + N2+ -> N2 + NH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(180) = small + (6.5e-10&
          *(T32)**(-0.5))
    end if

    !NH + O+ -> O + NH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(181) = small + (3.6e-10&
          *(T32)**(-0.5))
    end if

    !NO + H2CO+ -> H2CO + NO+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(182) = small + (7.8e-10)
    end if

    !NO + HNO+ -> HNO + NO+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(183) = small + (7e-10)
    end if

    !NO + O2+ -> O2 + NO+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(184) = small + (4.6e-10)
    end if

    !O+ + CH4 -> CH4+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(185) = small + (8.9e-10)
    end if

    !O+ + CO -> CO+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(186) = small + (4.9e-12&
          *(T32)**(0.5)*exp(-4580.0*invT))
    end if

    !O+ + H2CO -> H2CO+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(187) = small + (2.1e-09&
          *(T32)**(-0.5))
    end if

    !O+ + H2O -> H2O+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(188) = small + (3.2e-09&
          *(T32)**(-0.5))
    end if

    !O+ + HCO -> HCO+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(189) = small + (4.3e-10&
          *(T32)**(-0.5))
    end if

    !O+ + NH2 -> NH2+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(190) = small + (1e-09&
          *(T32)**(-0.5))
    end if

    !O+ + NH3 -> NH3+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(191) = small + (1.2e-09&
          *(T32)**(-0.5))
    end if

    !O+ + O2 -> O2+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(192) = small + (1.9e-11)
    end if

    !O+ + OH -> OH+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(193) = small + (3.6e-10&
          *(T32)**(-0.5))
    end if

    !O + CN+ -> CN + O+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(194) = small + (6.5e-11)
    end if

    !O + CO+ -> CO + O+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(195) = small + (1.4e-10)
    end if

    !O + N2+ -> N2 + O+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(196) = small + (1e-11)
    end if

    !OH+ + H2CO -> H2CO+ + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(197) = small + (7.44e-10&
          *(T32)**(-0.5))
    end if

    !OH+ + H2O -> H2O+ + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(198) = small + (1.59e-09&
          *(T32)**(-0.5))
    end if

    !OH+ + HCO -> HCO+ + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(199) = small + (2.8e-10&
          *(T32)**(-0.5))
    end if

    !OH+ + NH3 -> NH3+ + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(200) = small + (1.2e-09&
          *(T32)**(-0.5))
    end if

    !OH+ + NO -> NO+ + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(201) = small + (3.59e-10)
    end if

    !OH+ + O2 -> O2+ + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(202) = small + (5.9e-10)
    end if

    !OH + CN+ -> CN + OH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(203) = small + (6.4e-10&
          *(T32)**(-0.5))
    end if

    !OH + CO+ -> CO + OH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(204) = small + (3.1e-10&
          *(T32)**(-0.5))
    end if

    !OH + N2+ -> N2 + OH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(205) = small + (6.3e-10&
          *(T32)**(-0.5))
    end if

    !C -> C+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(206) = small + (2.3e-17&
          *user_zeta)
    end if

    !CO -> CO+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(207) = small + (3.9e-17&
          *user_zeta)
    end if

    !H2 -> H+ + H + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(208) = small + (2.86e-19&
          *user_zeta)
    end if

    !H2 -> H2+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(209) = small + (1.2e-17&
          *user_zeta)
    end if

    !H2 -> H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(210) = small + (1.3e-18&
          *user_zeta)
    end if

    !H -> H+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(211) = small + (5.98e-18&
          *user_zeta)
    end if

    !HE -> HE+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(212) = small + (6.5e-18&
          *user_zeta)
    end if

    !N -> N+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(213) = small + (2.7e-17&
          *user_zeta)
    end if

    !O -> O+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(214) = small + (3.4e-17&
          *user_zeta)
    end if

    !C -> C+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(215) = small + (1.3e-17&
          *255.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !CH+ -> C+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(216) = small + (1.3e-17&
          *88.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !CH2 -> CH2+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(217) = small + (1.3e-17&
          *250.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !CH2 -> CH + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(218) = small + (1.3e-17&
          *250.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !CH3 -> CH2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(219) = small + (1.3e-17&
          *250.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !CH3 -> CH3+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(220) = small + (1.3e-17&
          *250.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !CH3 -> CH + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(221) = small + (1.3e-17&
          *250.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !CH3OH -> H2CO + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(222) = small + (1.3e-17&
          *1584.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !CH3OH -> OH + CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(223) = small + (1.3e-17&
          *752.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !CH4 -> CH2 + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(224) = small + (1.3e-17&
          *1169.5*1.d0&
          /(1.d0-user_dgomega))
    end if

    !CH -> C + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(225) = small + (1.3e-17&
          *365.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !CN -> N + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(226) = small + (1.3e-17&
          *5290.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !CO2 -> CO + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(227) = small + (1.3e-17&
          *854.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !CO -> O + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(228) = small + (1.3e-17&
          *(T32)**(1.17)*105.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !H2CN -> HCN + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(229) = small + (1.3e-17&
          *750.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !H2CO -> CO + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(230) = small + (1.3e-17&
          *1329.5*1.d0&
          /(1.d0-user_dgomega))
    end if

    !H2O -> OH + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(231) = small + (1.3e-17&
          *485.5*1.d0&
          /(1.d0-user_dgomega))
    end if

    !H -> H+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(232) = small + (1.3e-17&
          *0.2*1.d0&
          /(1.d0-user_dgomega))
    end if

    !HCN -> CN + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(233) = small + (1.3e-17&
          *1557.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !HCO -> CO + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(234) = small + (1.3e-17&
          *210.5*1.d0&
          /(1.d0-user_dgomega))
    end if

    !HCO -> HCO+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(235) = small + (1.3e-17&
          *584.5*1.d0&
          /(1.d0-user_dgomega))
    end if

    !HNC -> CN + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(236) = small + (1.3e-17&
          *1500.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !HNCO -> NH + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(237) = small + (1.3e-17&
          *1500.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !HNO -> NO + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(238) = small + (1.3e-17&
          *500.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !HE -> HE+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(239) = small + (1.3e-17&
          *0.2*1.d0&
          /(1.d0-user_dgomega))
    end if

    !MG -> MG+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(240) = small + (1.3e-17&
          *66.5*1.d0&
          /(1.d0-user_dgomega))
    end if

    !N2 -> N + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(241) = small + (1.3e-17&
          *25.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !N -> N+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(242) = small + (1.3e-17&
          *1.1*1.d0&
          /(1.d0-user_dgomega))
    end if

    !NH2 -> NH2+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(243) = small + (1.3e-17&
          *324.5*1.d0&
          /(1.d0-user_dgomega))
    end if

    !NH2 -> NH + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(244) = small + (1.3e-17&
          *40.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !NH3 -> NH2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(245) = small + (1.3e-17&
          *657.5*1.d0&
          /(1.d0-user_dgomega))
    end if

    !NH3 -> NH3+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(246) = small + (1.3e-17&
          *288.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !NH3 -> NH + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(247) = small + (1.3e-17&
          *270.5*1.d0&
          /(1.d0-user_dgomega))
    end if

    !NH -> N + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(248) = small + (1.3e-17&
          *250.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !NH -> NH+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(249) = small + (1.3e-17&
          *250.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !NO2 -> NO + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(250) = small + (1.3e-17&
          *750.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !NO -> NO+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(251) = small + (1.3e-17&
          *247.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !NO -> O + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(252) = small + (1.3e-17&
          *231.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !O2 -> O2+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(253) = small + (1.3e-17&
          *58.5*1.d0&
          /(1.d0-user_dgomega))
    end if

    !O2 -> O + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(254) = small + (1.3e-17&
          *375.5*1.d0&
          /(1.d0-user_dgomega))
    end if

    !O2H -> O2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(255) = small + (1.3e-17&
          *375.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !O -> O+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(256) = small + (1.3e-17&
          *1.4*1.d0&
          /(1.d0-user_dgomega))
    end if

    !OCN -> CN + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(257) = small + (1.3e-17&
          *750.0*1.d0&
          /(1.d0-user_dgomega))
    end if

    !OH -> O + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(258) = small + (1.3e-17&
          *254.5*1.d0&
          /(1.d0-user_dgomega))
    end if

    !CH+ + E -> C + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(259) = small + (1.5e-07&
          *(T32)**(-0.42))
    end if

    !CH2+ + E -> C + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(260) = small + (7.68e-08&
          *(T32)**(-0.6))
    end if

    !CH2+ + E -> C + H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(261) = small + (4.03e-07&
          *(T32)**(-0.6))
    end if

    !CH2+ + E -> CH + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(262) = small + (1.6e-07&
          *(T32)**(-0.6))
    end if

    !CH3+ + E -> CH2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(263) = small + (7.75e-08&
          *(T32)**(-0.5))
    end if

    !CH3+ + E -> CH + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(264) = small + (1.95e-07&
          *(T32)**(-0.5))
    end if

    !CH3+ + E -> CH + H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(265) = small + (2e-07&
          *(T32)**(-0.4))
    end if

    !CH4+ + E -> CH2 + H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(266) = small + (1.75e-07&
          *(T32)**(-0.5))
    end if

    !CH4+ + E -> CH3 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(267) = small + (1.75e-07&
          *(T32)**(-0.5))
    end if

    !CN+ + E -> N + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(268) = small + (1.8e-07&
          *(T32)**(-0.5))
    end if

    !CO+ + E -> O + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(269) = small + (2e-07&
          *(T32)**(-0.48))
    end if

    !H2+ + E -> H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(270) = small + (1.6e-08&
          *(T32)**(-0.43))
    end if

    !H2CO+ + E -> CH2 + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(271) = small + (2.5e-08&
          *(T32)**(-0.7))
    end if

    !H2CO+ + E -> CO + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(272) = small + (7.5e-08&
          *(T32)**(-0.7))
    end if

    !H2CO+ + E -> CO + H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(273) = small + (2.5e-07&
          *(T32)**(-0.7))
    end if

    !H2CO+ + E -> HCO + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(274) = small + (1.6e-07&
          *(T32)**(-0.7))
    end if

    !H2NO+ + E -> HNO + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(275) = small + (1.5e-07&
          *(T32)**(-0.5))
    end if

    !H2NO+ + E -> NO + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(276) = small + (1.5e-07&
          *(T32)**(-0.5))
    end if

    !H2O+ + E -> O + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(277) = small + (3.9e-08&
          *(T32)**(-0.5))
    end if

    !H2O+ + E -> O + H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(278) = small + (3.05e-07&
          *(T32)**(-0.5))
    end if

    !H2O+ + E -> OH + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(279) = small + (8.6e-08&
          *(T32)**(-0.5))
    end if

    !H3+ + E -> H2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(280) = small + (2.34e-08&
          *(T32)**(-0.52))
    end if

    !H3+ + E -> H + H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(281) = small + (4.36e-08&
          *(T32)**(-0.52))
    end if

    !H3CO+ + E -> CH2 + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(282) = small + (4.2e-08&
          *(T32)**(-0.78))
    end if

    !H3CO+ + E -> CH + H2O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(283) = small + (1.4e-08&
          *(T32)**(-0.78))
    end if

    !H3CO+ + E -> CO + H2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(284) = small + (2.1e-07&
          *(T32)**(-0.78))
    end if

    !H3CO+ + E -> H2CO + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(285) = small + (2.17e-07&
          *(T32)**(-0.78))
    end if

    !H3CO+ + E -> HCO + H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(286) = small + (2.17e-07&
          *(T32)**(-0.78))
    end if

    !H3O+ + E -> H2O + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(287) = small + (7.09e-08&
          *(T32)**(-0.5))
    end if

    !H3O+ + E -> O + H2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(288) = small + (5.6e-09&
          *(T32)**(-0.5))
    end if

    !H3O+ + E -> OH + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(289) = small + (5.37e-08&
          *(T32)**(-0.5))
    end if

    !H3O+ + E -> OH + H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(290) = small + (3.05e-07&
          *(T32)**(-0.5))
    end if

    !HCN+ + E -> CN + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(291) = small + (2e-07&
          *(T32)**(-0.5))
    end if

    !HCNH+ + E -> CN + H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(292) = small + (9.3e-08&
          *(T32)**(-0.65))
    end if

    !HCNH+ + E -> HCN + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(293) = small + (9.5e-08&
          *(T32)**(-0.65))
    end if

    !HCNH+ + E -> HNC + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(294) = small + (9.5e-08&
          *(T32)**(-0.65))
    end if

    !HCO+ + E -> CO + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(295) = small + (2.4e-07&
          *(T32)**(-0.69))
    end if

    !HCO2+ + E -> CO2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(296) = small + (6e-08&
          *(T32)**(-0.64))
    end if

    !HCO2+ + E -> CO + O + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(297) = small + (8.1e-07&
          *(T32)**(-0.64))
    end if

    !HCO2+ + E -> CO + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(298) = small + (3.2e-07&
          *(T32)**(-0.64))
    end if

    !HNO+ + E -> NO + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(299) = small + (3e-07&
          *(T32)**(-0.5))
    end if

    !HOC+ + E -> CO + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(300) = small + (1.1e-07&
          *(T32)**(-1.0))
    end if

    !HEH+ + E -> HE + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(301) = small + (1e-08&
          *(T32)**(-0.6))
    end if

    !N2+ + E -> N + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(302) = small + (1.7e-07&
          *(T32)**(-0.3))
    end if

    !N2H+ + E -> N2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(303) = small + (2.77e-07&
          *(T32)**(-0.74))
    end if

    !N2H+ + E -> N + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(304) = small + (2.09e-08&
          *(T32)**(-0.74))
    end if

    !NH+ + E -> N + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(305) = small + (4.3e-08&
          *(T32)**(-0.5))
    end if

    !NH2+ + E -> N + H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(306) = small + (1.78e-07&
          *(T32)**(-0.8)*exp(-17.1*invT))
    end if

    !NH2+ + E -> NH + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(307) = small + (9.21e-08&
          *(T32)**(-0.79)*exp(-17.1*invT))
    end if

    !NH3+ + E -> NH2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(308) = small + (1.55e-07&
          *(T32)**(-0.5))
    end if

    !NH3+ + E -> NH + H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(309) = small + (1.55e-07&
          *(T32)**(-0.5))
    end if

    !NO+ + E -> O + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(310) = small + (4.3e-07&
          *(T32)**(-0.37))
    end if

    !O2+ + E -> O + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(311) = small + (1.95e-07&
          *(T32)**(-0.7))
    end if

    !O2H+ + E -> O2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(312) = small + (3e-07&
          *(T32)**(-0.5))
    end if

    !OH+ + E -> O + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(313) = small + (3.75e-08&
          *(T32)**(-0.5))
    end if

    !C+ + CH3OH -> H3CO+ + CH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(314) = small + (5.2e-10&
          *(T32)**(-0.5))
    end if

    !C+ + CH3OH -> HCO + CH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(315) = small + (2.08e-09&
          *(T32)**(-0.5))
    end if

    !C+ + CO2 -> CO+ + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(316) = small + (1.1e-09)
    end if

    !C+ + H2CO -> CO + CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(317) = small + (2.34e-09&
          *(T32)**(-0.5))
    end if

    !C+ + H2CO -> HCO+ + CH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(318) = small + (7.8e-10&
          *(T32)**(-0.5))
    end if

    !C+ + H2O -> HCO+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(319) = small + (9e-10&
          *(T32)**(-0.5))
    end if

    !C+ + H2O -> HOC+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(320) = small + (2.09e-09&
          *(T32)**(-0.5))
    end if

    !C+ + HCO -> CO + CH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(321) = small + (4.8e-10&
          *(T32)**(-0.5))
    end if

    !C+ + NH2 -> HCN+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(322) = small + (1.1e-09&
          *(T32)**(-0.5))
    end if

    !C+ + NH3 -> HCN+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(323) = small + (1.2e-10&
          *exp(+0.5*invT))
    end if

    !C+ + NH -> CN+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(324) = small + (7.8e-10&
          *(T32)**(-0.5))
    end if

    !C+ + O2 -> CO+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(325) = small + (3.42e-10)
    end if

    !C+ + O2 -> CO + O+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(326) = small + (4.54e-10)
    end if

    !C+ + OCN -> CO+ + CN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(327) = small + (3.8e-09)
    end if

    !C+ + OH -> CO+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(328) = small + (7.7e-10&
          *(T32)**(-0.5))
    end if

    !C + H2O+ -> OH + CH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(329) = small + (1.1e-09)
    end if

    !C + H3O+ -> HCO+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(330) = small + (1e-11)
    end if

    !C + HCN+ -> CN + CH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(331) = small + (1.1e-09)
    end if

    !C + HCO+ -> CO + CH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(332) = small + (1.1e-09)
    end if

    !C + HCO2+ -> CO2 + CH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(333) = small + (1e-09)
    end if

    !C + HNO+ -> NO + CH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(334) = small + (1e-09)
    end if

    !C + N2H+ -> N2 + CH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(335) = small + (1.1e-09)
    end if

    !C + NH+ -> N + CH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(336) = small + (1.6e-09)
    end if

    !C + O2+ -> CO+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(337) = small + (5.2e-11)
    end if

    !C + O2H+ -> O2 + CH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(338) = small + (1e-09)
    end if

    !C + OH+ -> O + CH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(339) = small + (1.2e-09)
    end if

    !CH+ + CH3OH -> H2CO + CH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(340) = small + (1.45e-09&
          *(T32)**(-0.5))
    end if

    !CH+ + CH3OH -> H3CO+ + CH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(341) = small + (2.9e-10&
          *(T32)**(-0.5))
    end if

    !CH+ + CO2 -> HCO+ + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(342) = small + (1.6e-09)
    end if

    !CH+ + H2CO -> CO + CH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(343) = small + (9.6e-10&
          *(T32)**(-0.5))
    end if

    !CH+ + H2CO -> H3CO+ + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(344) = small + (9.6e-10&
          *(T32)**(-0.5))
    end if

    !CH+ + H2CO -> HCO+ + CH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(345) = small + (9.6e-10&
          *(T32)**(-0.5))
    end if

    !CH+ + H2O -> H2CO+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(346) = small + (5.8e-10&
          *(T32)**(-0.5))
    end if

    !CH+ + H2O -> H3O+ + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(347) = small + (5.8e-10&
          *(T32)**(-0.5))
    end if

    !CH+ + H2O -> HCO+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(348) = small + (2.9e-09&
          *(T32)**(-0.5))
    end if

    !CH+ + HCN -> HCNH+ + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(349) = small + (1.8e-09&
          *(T32)**(-0.5))
    end if

    !CH+ + HCO -> CO + CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(350) = small + (4.6e-10&
          *(T32)**(-0.5))
    end if

    !CH+ + HNC -> HCNH+ + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(351) = small + (1.8e-09&
          *(T32)**(-0.5))
    end if

    !CH+ + N -> CN+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(352) = small + (1.9e-10)
    end if

    !CH+ + NH2 -> HCN+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(353) = small + (1.1e-09&
          *(T32)**(-0.5))
    end if

    !CH+ + NH -> CN+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(354) = small + (7.6e-10&
          *(T32)**(-0.5))
    end if

    !CH+ + O2 -> CO+ + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(355) = small + (1e-11)
    end if

    !CH+ + O2 -> HCO+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(356) = small + (9.7e-10)
    end if

    !CH+ + O2 -> HCO + O+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(357) = small + (1e-11)
    end if

    !CH+ + O -> CO+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(358) = small + (3.5e-10)
    end if

    !CH+ + OH -> CO+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(359) = small + (7.5e-10&
          *(T32)**(-0.5))
    end if

    !CH2+ + CO2 -> H2CO+ + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(360) = small + (1.6e-09)
    end if

    !CH2+ + H2CO -> HCO+ + CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(361) = small + (2.81e-09&
          *(T32)**(-0.5))
    end if

    !CH2+ + H2O -> H3CO+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(362) = small + (1.2e-09&
          *(T32)**(-0.5))
    end if

    !CH2+ + HCO -> CO + CH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(363) = small + (4.5e-10&
          *(T32)**(-0.5))
    end if

    !CH2+ + O2 -> HCO+ + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(364) = small + (9.1e-10)
    end if

    !CH2+ + O -> HCO+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(365) = small + (7.5e-10)
    end if

    !CH2 + CO+ -> HCO+ + CH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(366) = small + (4.3e-10)
    end if

    !CH2 + H2CO+ -> HCO + CH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(367) = small + (4.3e-10)
    end if

    !CH2 + H2O+ -> OH + CH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(368) = small + (4.7e-10)
    end if

    !CH2 + H3O+ -> H2O + CH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(369) = small + (9.4e-10)
    end if

    !CH2 + HCN+ -> CN + CH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(370) = small + (8.7e-10)
    end if

    !CH2 + HCNH+ -> HCN + CH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(371) = small + (4.35e-10)
    end if

    !CH2 + HCNH+ -> HNC + CH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(372) = small + (4.35e-10)
    end if

    !CH2 + HCO+ -> CO + CH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(373) = small + (8.6e-10)
    end if

    !CH2 + HNO+ -> NO + CH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(374) = small + (8.6e-10)
    end if

    !CH2 + N2H+ -> N2 + CH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(375) = small + (8.6e-10)
    end if

    !CH2 + NH+ -> CH3+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(376) = small + (1.4e-09)
    end if

    !CH2 + NH2+ -> CH3+ + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(377) = small + (4.9e-10)
    end if

    !CH2 + NH3+ -> NH2 + CH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(378) = small + (9.6e-10)
    end if

    !CH2 + O2+ -> H2CO+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(379) = small + (4.3e-10)
    end if

    !CH2 + O2H+ -> O2 + CH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(380) = small + (8.5e-10)
    end if

    !CH2 + OH+ -> O + CH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(381) = small + (4.8e-10)
    end if

    !CH3+ + CH3OH -> H3CO+ + CH4
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(382) = small + (2.3e-09&
          *(T32)**(-0.5))
    end if

    !CH3+ + H2CO -> HCO+ + CH4
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(383) = small + (1.6e-09&
          *(T32)**(-0.5))
    end if

    !CH3+ + HCO -> CO + CH4+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(384) = small + (4.4e-10&
          *(T32)**(-0.5))
    end if

    !CH3+ + O2 -> H3CO+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(385) = small + (5e-12)
    end if

    !CH3+ + O -> H2CO+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(386) = small + (4e-11)
    end if

    !CH3+ + O -> HCO+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(387) = small + (4e-10)
    end if

    !CH3+ + OH -> H2CO+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(388) = small + (7.2e-10&
          *(T32)**(-0.5))
    end if

    !CH4+ + CO2 -> HCO2+ + CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(389) = small + (1.2e-09)
    end if

    !CH4+ + CO -> HCO+ + CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(390) = small + (1.4e-09)
    end if

    !CH4+ + H2CO -> H3CO+ + CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(391) = small + (1.98e-09&
          *(T32)**(-0.5))
    end if

    !CH4+ + H2O -> H3O+ + CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(392) = small + (2.6e-09&
          *(T32)**(-0.5))
    end if

    !CH4 + CO+ -> HCO+ + CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(393) = small + (4.55e-10)
    end if

    !CH4 + H2CO+ -> H3CO+ + CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(394) = small + (9.35e-11)
    end if

    !CH4 + H2O+ -> H3O+ + CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(395) = small + (1.4e-09)
    end if

    !CH4 + HCN+ -> HCNH+ + CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(396) = small + (1.04e-09)
    end if

    !CH4 + N2+ -> N2 + CH2+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(397) = small + (7e-11)
    end if

    !CH4 + N2+ -> N2 + CH3+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(398) = small + (9.3e-10)
    end if

    !CH4 + OH+ -> H3O+ + CH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(399) = small + (1.31e-09)
    end if

    !CH + CO+ -> HCO+ + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(400) = small + (3.2e-10&
          *(T32)**(-0.5))
    end if

    !CH + H2CO+ -> HCO + CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(401) = small + (3.1e-10&
          *(T32)**(-0.5))
    end if

    !CH + H2O+ -> OH + CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(402) = small + (3.4e-10&
          *(T32)**(-0.5))
    end if

    !CH + H3CO+ -> H2CO + CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(403) = small + (6.2e-10&
          *(T32)**(-0.5))
    end if

    !CH + H3O+ -> H2O + CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(404) = small + (6.8e-10&
          *(T32)**(-0.5))
    end if

    !CH + HCN+ -> CN + CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(405) = small + (6.3e-10&
          *(T32)**(-0.5))
    end if

    !CH + HCNH+ -> HCN + CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(406) = small + (3.15e-10&
          *(T32)**(-0.5))
    end if

    !CH + HCNH+ -> HNC + CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(407) = small + (3.15e-10&
          *(T32)**(-0.5))
    end if

    !CH + HCO+ -> CO + CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(408) = small + (6.3e-10&
          *(T32)**(-0.5))
    end if

    !CH + HNO+ -> NO + CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(409) = small + (6.2e-10&
          *(T32)**(-0.5))
    end if

    !CH + N+ -> CN+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(410) = small + (3.6e-10&
          *(T32)**(-0.5))
    end if

    !CH + N2H+ -> N2 + CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(411) = small + (6.3e-10&
          *(T32)**(-0.5))
    end if

    !CH + NH+ -> CH2+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(412) = small + (9.9e-10&
          *(T32)**(-0.5))
    end if

    !CH + NH2+ -> NH + CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(413) = small + (3.5e-10&
          *(T32)**(-0.5))
    end if

    !CH + O+ -> CO+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(414) = small + (3.5e-10&
          *(T32)**(-0.5))
    end if

    !CH + O2+ -> HCO+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(415) = small + (3.1e-10&
          *(T32)**(-0.5))
    end if

    !CH + O2H+ -> O2 + CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(416) = small + (6.2e-10&
          *(T32)**(-0.5))
    end if

    !CH + OH+ -> O + CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(417) = small + (3.5e-10&
          *(T32)**(-0.5))
    end if

    !CN+ + H2CO -> HCO+ + HCN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(418) = small + (5.2e-10&
          *(T32)**(-0.5))
    end if

    !CN+ + HCO -> CO + HCN+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(419) = small + (3.7e-10&
          *(T32)**(-0.5))
    end if

    !CN+ + O2 -> NO+ + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(420) = small + (8.6e-11)
    end if

    !CN + HNO+ -> NO + HCN+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(421) = small + (8.7e-10&
          *(T32)**(-0.5))
    end if

    !CN + O2H+ -> O2 + HCN+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(422) = small + (8.6e-10&
          *(T32)**(-0.5))
    end if

    !CO+ + H2CO -> HCO+ + HCO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(423) = small + (1.65e-09&
          *(T32)**(-0.5))
    end if

    !CO + HCO2+ -> CO2 + HCO+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(424) = small + (7.8e-10)
    end if

    !CO + HNO+ -> NO + HCO+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(425) = small + (1e-10)
    end if

    !CO + N2H+ -> HCO+ + N2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(426) = small + (8.8e-10)
    end if

    !CO + O2H+ -> O2 + HCO+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(427) = small + (8.4e-10)
    end if

    !H+ + CH2 -> CH+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(428) = small + (1.4e-09)
    end if

    !H+ + CH3OH -> CH3+ + H2O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(429) = small + (5.9e-10&
          *(T32)**(-0.5))
    end if

    !H+ + CH3OH -> H3CO+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(430) = small + (3.84e-09&
          *(T32)**(-0.5))
    end if

    !H+ + CH3OH -> HCO+ + H2 + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(431) = small + (8.85e-10&
          *(T32)**(-0.5))
    end if

    !H+ + CH4 -> CH3+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(432) = small + (2.3e-09)
    end if

    !H+ + CO2 -> HCO+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(433) = small + (3.5e-09)
    end if

    !H+ + H2CO -> CO+ + H2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(434) = small + (1.06e-09&
          *(T32)**(-0.5))
    end if

    !H+ + H2CO -> HCO+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(435) = small + (3.57e-09&
          *(T32)**(-0.5))
    end if

    !H+ + HCO -> CO+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(436) = small + (9.4e-10&
          *(T32)**(-0.5))
    end if

    !H+ + HCO -> CO + H2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(437) = small + (9.4e-10&
          *(T32)**(-0.5))
    end if

    !H+ + HNCO -> NH2+ + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(438) = small + (7.94e-09&
          *(T32)**(-0.5))
    end if

    !H+ + HNO -> NO+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(439) = small + (4e-09&
          *(T32)**(-0.5))
    end if

    !H+ + NO2 -> NO+ + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(440) = small + (1.9e-09)
    end if

    !H2+ + C -> CH+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(441) = small + (2.4e-09)
    end if

    !H2+ + CH2 -> CH3+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(442) = small + (1e-09)
    end if

    !H2+ + CH4 -> CH3+ + H2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(443) = small + (2.3e-09)
    end if

    !H2+ + CH -> CH2+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(444) = small + (7.1e-10&
          *(T32)**(-0.5))
    end if

    !H2+ + CN -> HCN+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(445) = small + (1.2e-09&
          *(T32)**(-0.5))
    end if

    !H2+ + CO2 -> HCO2+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(446) = small + (2.35e-09)
    end if

    !H2+ + CO -> HCO+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(447) = small + (2.16e-09)
    end if

    !H2+ + H2 -> H3+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(448) = small + (2.08e-09)
    end if

    !H2+ + H2CO -> HCO+ + H2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(449) = small + (1.4e-09&
          *(T32)**(-0.5))
    end if

    !H2+ + H2O -> H3O+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(450) = small + (3.4e-09&
          *(T32)**(-0.5))
    end if

    !H2+ + HCO -> CO + H3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(451) = small + (1e-09&
          *(T32)**(-0.5))
    end if

    !H2+ + HE -> HEH+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(452) = small + (1.3e-10)
    end if

    !H2+ + N2 -> N2H+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(453) = small + (2e-09)
    end if

    !H2+ + N -> NH+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(454) = small + (1.9e-09)
    end if

    !H2+ + NH -> NH2+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(455) = small + (7.6e-10&
          *(T32)**(-0.5))
    end if

    !H2+ + NO -> HNO+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(456) = small + (1.1e-09)
    end if

    !H2+ + O2 -> O2H+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(457) = small + (1.9e-09)
    end if

    !H2+ + O -> OH+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(458) = small + (1.5e-09)
    end if

    !H2+ + OH -> H2O+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(459) = small + (7.6e-10&
          *(T32)**(-0.5))
    end if

    !H2 + C+ -> CH+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(460) = small + (1e-10&
          *exp(-4640.0*invT))
    end if

    !H2 + CH+ -> CH2+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(461) = small + (1.2e-09)
    end if

    !H2 + CH2+ -> CH3+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(462) = small + (1.6e-09)
    end if

    !H2 + CN+ -> HCN+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(463) = small + (1e-09)
    end if

    !H2 + CO+ -> HCO+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(464) = small + (7.5e-10)
    end if

    !H2 + CO+ -> HOC+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(465) = small + (7.5e-10)
    end if

    !H2 + H2O+ -> H3O+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(466) = small + (6.4e-10)
    end if

    !H2 + HCN+ -> HCNH+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(467) = small + (9e-10)
    end if

    !H2 + HE+ -> HE + H+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(468) = small + (3.7e-14&
          *exp(-35.0*invT))
    end if

    !H2 + HEH+ -> HE + H3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(469) = small + (1.5e-09)
    end if

    !H2 + N+ -> NH+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(470) = small + (1e-09&
          *exp(-85.0*invT))
    end if

    !H2 + N2+ -> N2H+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(471) = small + (2e-09)
    end if

    !H2 + NH+ -> N + H3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(472) = small + (2.25e-10)
    end if

    !H2 + NH+ -> NH2+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(473) = small + (1.28e-09)
    end if

    !H2 + NH2+ -> NH3+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(474) = small + (2.7e-10)
    end if

    !H2 + O+ -> OH+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(475) = small + (1.7e-09)
    end if

    !H2 + O2H+ -> O2 + H3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(476) = small + (6.4e-10)
    end if

    !H2 + OH+ -> H2O+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(477) = small + (1.01e-09)
    end if

    !H2CO+ + H2CO -> H3CO+ + HCO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(478) = small + (3.2e-09&
          *(T32)**(-0.5))
    end if

    !H2CO+ + O2 -> O2H + HCO+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(479) = small + (7.7e-11)
    end if

    !H2CO + HNO+ -> H3CO+ + NO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(480) = small + (1e-09&
          *(T32)**(-0.5))
    end if

    !H2CO + O2+ -> O2 + HCO+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(481) = small + (2.3e-10&
          *(T32)**(-0.5))
    end if

    !H2CO + O2H+ -> O2 + H3CO+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(482) = small + (9.8e-10&
          *(T32)**(-0.5))
    end if

    !H2O+ + CO -> HCO+ + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(483) = small + (5e-10)
    end if

    !H2O+ + H2CO -> H3CO+ + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(484) = small + (6.62e-10&
          *(T32)**(-0.5))
    end if

    !H2O+ + H2O -> H3O+ + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(485) = small + (2.1e-09&
          *(T32)**(-0.5))
    end if

    !H2O+ + HCN -> HCNH+ + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(486) = small + (2.1e-09&
          *(T32)**(-0.5))
    end if

    !H2O+ + HCO -> CO + H3O+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(487) = small + (2.8e-10&
          *(T32)**(-0.5))
    end if

    !H2O+ + HCO -> H2CO+ + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(488) = small + (2.8e-10&
          *(T32)**(-0.5))
    end if

    !H2O+ + HNC -> HCNH+ + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(489) = small + (1.1e-09&
          *(T32)**(-0.5))
    end if

    !H2O + CN+ -> HCN+ + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(490) = small + (1.6e-09&
          *(T32)**(-0.5))
    end if

    !H2O + CN+ -> HCO+ + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(491) = small + (1.6e-10&
          *(T32)**(-0.5))
    end if

    !H2O + CO+ -> HCO+ + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(492) = small + (8.84e-10&
          *(T32)**(-0.5))
    end if

    !H2O + H2CO+ -> HCO + H3O+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(493) = small + (2.6e-09&
          *(T32)**(-0.5))
    end if

    !H2O + H3CO+ -> H2CO + H3O+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(494) = small + (2.3e-10&
          *(T32)**(-0.5))
    end if

    !H2O + HCN+ -> CN + H3O+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(495) = small + (1.8e-09&
          *(T32)**(-0.5))
    end if

    !H2O + HCO+ -> CO + H3O+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(496) = small + (2.5e-09&
          *(T32)**(-0.5))
    end if

    !H2O + HCO2+ -> CO2 + H3O+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(497) = small + (2.3e-09&
          *(T32)**(-0.5))
    end if

    !H2O + HNO+ -> NO + H3O+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(498) = small + (2.3e-09&
          *(T32)**(-0.5))
    end if

    !H2O + N2+ -> N2H+ + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(499) = small + (5e-10&
          *(T32)**(-0.5))
    end if

    !H2O + N2H+ -> N2 + H3O+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(500) = small + (2.6e-09&
          *(T32)**(-0.5))
    end if

    !H2O + O2H+ -> O2 + H3O+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(501) = small + (8.2e-10&
          *(T32)**(-0.5))
    end if

    !H3+ + C -> CH+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(502) = small + (2e-09)
    end if

    !H3+ + CH2 -> CH3+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(503) = small + (1.7e-09)
    end if

    !H3+ + CH3 -> CH4+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(504) = small + (2.1e-09)
    end if

    !H3+ + CH3OH -> CH3+ + H2O + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(505) = small + (3.71e-09&
          *(T32)**(-0.5))
    end if

    !H3+ + CH -> CH2+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(506) = small + (1.2e-09&
          *(T32)**(-0.5))
    end if

    !H3+ + CN -> HCN+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(507) = small + (2e-09&
          *(T32)**(-0.5))
    end if

    !H3+ + CO2 -> HCO2+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(508) = small + (2e-09)
    end if

    !H3+ + CO -> HCO+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(509) = small + (1.36e-09&
          *(T32)**(-0.14)*exp(+3.4*invT))
    end if

    !H3+ + CO -> HOC+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(510) = small + (8.49e-10&
          *(T32)**(0.07)*exp(-5.2*invT))
    end if

    !H3+ + H2CO -> H3CO+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(511) = small + (6.3e-09&
          *(T32)**(-0.5))
    end if

    !H3+ + H2O -> H3O+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(512) = small + (5.9e-09&
          *(T32)**(-0.5))
    end if

    !H3+ + HCN -> HCNH+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(513) = small + (8.1e-09&
          *(T32)**(-0.5))
    end if

    !H3+ + HCO -> H2CO+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(514) = small + (1.7e-09&
          *(T32)**(-0.5))
    end if

    !H3+ + HNC -> HCNH+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(515) = small + (8.1e-09&
          *(T32)**(-0.5))
    end if

    !H3+ + HNO -> H2NO+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(516) = small + (2e-09&
          *(T32)**(-0.5))
    end if

    !H3+ + MG -> MG+ + H2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(517) = small + (1e-09)
    end if

    !H3+ + N2 -> N2H+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(518) = small + (1.8e-09)
    end if

    !H3+ + NH2 -> NH3+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(519) = small + (1.8e-09&
          *(T32)**(-0.5))
    end if

    !H3+ + NH -> NH2+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(520) = small + (1.3e-09&
          *(T32)**(-0.5))
    end if

    !H3+ + NO2 -> NO+ + OH + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(521) = small + (7e-10)
    end if

    !H3+ + NO -> HNO+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(522) = small + (1.1e-09)
    end if

    !H3+ + O2 -> O2H+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(523) = small + (9.3e-10&
          *exp(-100.0*invT))
    end if

    !H3+ + O -> H2O+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(524) = small + (3.42e-10&
          *(T32)**(-0.16)*exp(-1.4*invT))
    end if

    !H3+ + O -> OH+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(525) = small + (7.98e-10&
          *(T32)**(-0.16)*exp(-1.4*invT))
    end if

    !H3+ + OH -> H2O+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(526) = small + (1.3e-09&
          *(T32)**(-0.5))
    end if

    !H3O+ + H2CO -> H3CO+ + H2O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(527) = small + (3.4e-09&
          *(T32)**(-0.5))
    end if

    !H3O+ + HCN -> HCNH+ + H2O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(528) = small + (3.8e-09&
          *(T32)**(-0.5))
    end if

    !H3O+ + HNC -> HCNH+ + H2O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(529) = small + (4e-09&
          *(T32)**(-0.5))
    end if

    !H + CH+ -> C+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(530) = small + (9.06e-10&
          *(T32)**(-0.37)*exp(-29.1*invT))
    end if

    !H + CH2+ -> CH+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(531) = small + (1e-09&
          *exp(-7080.0*invT))
    end if

    !H + CH3+ -> CH2+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(532) = small + (7e-10&
          *exp(-10560.0*invT))
    end if

    !H + CH4+ -> CH3+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(533) = small + (1e-11)
    end if

    !H + HEH+ -> HE + H2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(534) = small + (9.1e-10)
    end if

    !HCN+ + CO2 -> HCO2+ + CN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(535) = small + (2.1e-10)
    end if

    !HCN+ + CO -> HCO+ + CN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(536) = small + (1.4e-10)
    end if

    !HCN+ + H2CO -> H3CO+ + CN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(537) = small + (1e-09&
          *(T32)**(-0.5))
    end if

    !HCN+ + HCN -> HCNH+ + CN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(538) = small + (1.6e-09&
          *(T32)**(-0.5))
    end if

    !HCN+ + HCO -> H2CO+ + CN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(539) = small + (3.7e-10&
          *(T32)**(-0.5))
    end if

    !HCN+ + HCO -> HCNH+ + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(540) = small + (3.7e-10&
          *(T32)**(-0.5))
    end if

    !HCN+ + HNC -> HCNH+ + CN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(541) = small + (1e-09&
          *(T32)**(-0.5))
    end if

    !HCN + H2CO+ -> HCO + HCNH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(542) = small + (1.4e-09&
          *(T32)**(-0.5))
    end if

    !HCN + H3CO+ -> H2CO + HCNH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(543) = small + (1.3e-09&
          *(T32)**(-0.5))
    end if

    !HCN + HCO+ -> HCNH+ + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(544) = small + (3.1e-09&
          *(T32)**(-0.5))
    end if

    !HCN + HNO+ -> NO + HCNH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(545) = small + (9.9e-10&
          *(T32)**(-0.5))
    end if

    !HCN + N2H+ -> HCNH+ + N2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(546) = small + (3.2e-09&
          *(T32)**(-0.5))
    end if

    !HCN + O2H+ -> O2 + HCNH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(547) = small + (9.7e-10&
          *(T32)**(-0.5))
    end if

    !HCNH+ + H2CO -> H3CO+ + HCN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(548) = small + (1.05e-09&
          *(T32)**(-0.5))
    end if

    !HCNH+ + H2CO -> H3CO+ + HNC
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(549) = small + (1.05e-09&
          *(T32)**(-0.5))
    end if

    !HCO+ + H2CO -> H3CO+ + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(550) = small + (3.3e-09&
          *(T32)**(-0.5))
    end if

    !HCO+ + HCO -> H2CO+ + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(551) = small + (7.3e-10&
          *(T32)**(-0.5))
    end if

    !HCO + H2CO+ -> H3CO+ + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(552) = small + (3.6e-10&
          *(T32)**(-0.5))
    end if

    !HCO + HNO+ -> H2CO+ + NO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(553) = small + (7.2e-10&
          *(T32)**(-0.5))
    end if

    !HCO + N2H+ -> H2CO+ + N2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(554) = small + (7.3e-10&
          *(T32)**(-0.5))
    end if

    !HCO + O2+ -> O2H+ + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(555) = small + (3.6e-10&
          *(T32)**(-0.5))
    end if

    !HCO + O2H+ -> O2 + H2CO+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(556) = small + (7.1e-10&
          *(T32)**(-0.5))
    end if

    !HNC + H2CO+ -> HCO + HCNH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(557) = small + (1.4e-09&
          *(T32)**(-0.5))
    end if

    !HNC + H3CO+ -> H2CO + HCNH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(558) = small + (1.3e-09&
          *(T32)**(-0.5))
    end if

    !HNC + HCO+ -> HCNH+ + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(559) = small + (3.1e-09&
          *(T32)**(-0.5))
    end if

    !HNC + HNO+ -> NO + HCNH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(560) = small + (9.9e-10&
          *(T32)**(-0.5))
    end if

    !HNC + N2H+ -> HCNH+ + N2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(561) = small + (3.2e-09&
          *(T32)**(-0.5))
    end if

    !HNC + O2H+ -> O2 + HCNH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(562) = small + (9.7e-10&
          *(T32)**(-0.5))
    end if

    !HNO+ + CO2 -> HCO2+ + NO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(563) = small + (1e-10)
    end if

    !HE+ + CH2 -> C+ + HE + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(564) = small + (7.5e-10)
    end if

    !HE+ + CH2 -> CH+ + HE + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(565) = small + (7.5e-10)
    end if

    !HE+ + CH3 -> CH+ + HE + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(566) = small + (1.8e-09)
    end if

    !HE+ + CH3OH -> OH+ + CH3 + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(567) = small + (1.1e-09&
          *(T32)**(-0.5))
    end if

    !HE+ + CH3OH -> OH + CH3+ + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(568) = small + (1.1e-09&
          *(T32)**(-0.5))
    end if

    !HE+ + CH4 -> CH+ + HE + H2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(569) = small + (2.4e-10)
    end if

    !HE+ + CH4 -> CH2+ + HE + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(570) = small + (9.5e-10)
    end if

    !HE+ + CH4 -> CH3+ + HE + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(571) = small + (8.5e-11)
    end if

    !HE+ + CH4 -> CH3 + HE + H+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(572) = small + (4.8e-10)
    end if

    !HE+ + CH -> C+ + HE + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(573) = small + (1.1e-09&
          *(T32)**(-0.5))
    end if

    !HE+ + CN -> N+ + C + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(574) = small + (8.8e-10&
          *(T32)**(-0.5))
    end if

    !HE+ + CN -> N + C+ + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(575) = small + (8.8e-10&
          *(T32)**(-0.5))
    end if

    !HE+ + CO2 -> CO+ + O + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(576) = small + (8.7e-10)
    end if

    !HE+ + CO2 -> CO + O+ + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(577) = small + (1e-10)
    end if

    !HE+ + CO2 -> O2+ + C + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(578) = small + (1.1e-11)
    end if

    !HE+ + CO2 -> O2 + C+ + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(579) = small + (4e-11)
    end if

    !HE+ + CO -> O + C+ + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(580) = small + (1.6e-09)
    end if

    !HE+ + H2CO -> CO+ + HE + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(581) = small + (1.88e-09&
          *(T32)**(-0.5))
    end if

    !HE+ + H2CO -> HCO+ + HE + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(582) = small + (1.14e-09&
          *(T32)**(-0.5))
    end if

    !HE+ + H2CO -> O + CH2+ + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(583) = small + (1.71e-09&
          *(T32)**(-0.5))
    end if

    !HE+ + H2O -> OH+ + HE + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(584) = small + (2.86e-10&
          *(T32)**(-0.5))
    end if

    !HE+ + H2O -> OH + HE + H+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(585) = small + (2.04e-10&
          *(T32)**(-0.5))
    end if

    !HE+ + HCN -> CN+ + HE + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(586) = small + (1.46e-09&
          *(T32)**(-0.5))
    end if

    !HE+ + HCN -> N+ + CH + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(587) = small + (2.17e-10&
          *(T32)**(-0.5))
    end if

    !HE+ + HCN -> N + C+ + HE + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(588) = small + (7.75e-10&
          *(T32)**(-0.5))
    end if

    !HE+ + HCN -> N + CH+ + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(589) = small + (6.51e-10&
          *(T32)**(-0.5))
    end if

    !HE+ + HCO -> CO+ + HE + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(590) = small + (4.9e-10&
          *(T32)**(-0.5))
    end if

    !HE+ + HCO -> CO + HEH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(591) = small + (3e-10&
          *(T32)**(-0.5))
    end if

    !HE+ + HCO -> O + CH+ + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(592) = small + (4.9e-10&
          *(T32)**(-0.5))
    end if

    !HE+ + HNC -> CN+ + HE + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(593) = small + (5e-10&
          *(T32)**(-0.5))
    end if

    !HE+ + HNC -> N + C+ + HE + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(594) = small + (5e-10&
          *(T32)**(-0.5))
    end if

    !HE+ + HNC -> NH+ + C + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(595) = small + (5e-10&
          *(T32)**(-0.5))
    end if

    !HE+ + HNO -> NO+ + HE + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(596) = small + (1e-09&
          *(T32)**(-0.5))
    end if

    !HE+ + HNO -> NO + HE + H+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(597) = small + (1e-09&
          *(T32)**(-0.5))
    end if

    !HE+ + N2 -> N+ + N + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(598) = small + (9.6e-10)
    end if

    !HE+ + NH2 -> N+ + HE + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(599) = small + (8e-10&
          *(T32)**(-0.5))
    end if

    !HE+ + NH2 -> NH+ + HE + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(600) = small + (8e-10&
          *(T32)**(-0.5))
    end if

    !HE+ + NH3 -> NH+ + HE + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(601) = small + (1.76e-10&
          *(T32)**(-0.5))
    end if

    !HE+ + NH3 -> NH2+ + HE + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(602) = small + (1.76e-09&
          *(T32)**(-0.5))
    end if

    !HE+ + NH -> N+ + HE + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(603) = small + (1.1e-09&
          *(T32)**(-0.5))
    end if

    !HE+ + NO -> O+ + N + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(604) = small + (2e-10)
    end if

    !HE+ + NO -> O + N+ + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(605) = small + (1.4e-09)
    end if

    !HE+ + O2 -> O+ + O + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(606) = small + (1.1e-09)
    end if

    !HE+ + OCN -> CN+ + O + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(607) = small + (3e-09)
    end if

    !HE+ + OCN -> CN + O+ + HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(608) = small + (3e-09)
    end if

    !HE+ + OH -> O+ + HE + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(609) = small + (1.1e-09&
          *(T32)**(-0.5))
    end if

    !N+ + CH3OH -> H2CO+ + NH + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(610) = small + (9.3e-10&
          *(T32)**(-0.5))
    end if

    !N+ + CH3OH -> H3CO+ + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(611) = small + (4.96e-10&
          *(T32)**(-0.5))
    end if

    !N+ + CH3OH -> NO+ + CH3 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(612) = small + (3.1e-10&
          *(T32)**(-0.5))
    end if

    !N+ + CH3OH -> NO + CH3+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(613) = small + (1.24e-10&
          *(T32)**(-0.5))
    end if

    !N+ + CH4 -> CH3+ + N + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(614) = small + (4.7e-10)
    end if

    !N+ + CH4 -> HCN+ + H2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(615) = small + (5.6e-11)
    end if

    !N+ + CH4 -> HCNH+ + H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(616) = small + (3.8e-10)
    end if

    !N+ + CO2 -> NO + CO+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(617) = small + (2.5e-10)
    end if

    !N+ + CO -> NO+ + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(618) = small + (1.45e-10)
    end if

    !N+ + H2CO -> HCO+ + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(619) = small + (7.25e-10&
          *(T32)**(-0.5))
    end if

    !N+ + H2CO -> NO+ + CH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(620) = small + (2.9e-10&
          *(T32)**(-0.5))
    end if

    !N+ + HCO -> CO + NH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(621) = small + (4.5e-10&
          *(T32)**(-0.5))
    end if

    !N+ + NH3 -> N2H+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(622) = small + (2.16e-10&
          *(T32)**(-0.5))
    end if

    !N+ + NH3 -> NH2+ + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(623) = small + (2.16e-10&
          *(T32)**(-0.5))
    end if

    !N+ + NH -> N2+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(624) = small + (3.7e-10&
          *(T32)**(-0.5))
    end if

    !N+ + NO -> N2+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(625) = small + (7.9e-11)
    end if

    !N+ + O2 -> NO+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(626) = small + (2.63e-10)
    end if

    !N+ + O2 -> NO + O+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(627) = small + (3.66e-11)
    end if

    !N2+ + H2CO -> HCO+ + N2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(628) = small + (2.52e-09&
          *(T32)**(-0.5))
    end if

    !N2+ + HCO -> N2H+ + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(629) = small + (3.7e-10&
          *(T32)**(-0.5))
    end if

    !N2 + HNO+ -> NO + N2H+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(630) = small + (1e-11)
    end if

    !N2 + O2H+ -> O2 + N2H+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(631) = small + (8e-10)
    end if

    !N2H+ + CO2 -> HCO2+ + N2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(632) = small + (9.8e-10)
    end if

    !N2H+ + H2CO -> H3CO+ + N2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(633) = small + (3.3e-09&
          *(T32)**(-0.5))
    end if

    !N + CH2+ -> HCN+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(634) = small + (2.2e-10)
    end if

    !N + CN+ -> N2+ + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(635) = small + (6.1e-10)
    end if

    !N + H2O+ -> HNO+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(636) = small + (1.12e-10)
    end if

    !N + H2O+ -> NO+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(637) = small + (2.8e-11)
    end if

    !N + NH+ -> N2+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(638) = small + (1.3e-09)
    end if

    !N + NH2+ -> N2H+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(639) = small + (9.1e-11)
    end if

    !N + O2+ -> NO+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(640) = small + (1.8e-10)
    end if

    !N + OH+ -> NO+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(641) = small + (8.9e-10)
    end if

    !NH+ + CN -> HCN+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(642) = small + (1.6e-09&
          *(T32)**(-0.5))
    end if

    !NH+ + CO2 -> HCO2+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(643) = small + (3.85e-10)
    end if

    !NH+ + CO2 -> HNO+ + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(644) = small + (3.85e-10)
    end if

    !NH+ + CO2 -> NO+ + HCO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(645) = small + (3.3e-10)
    end if

    !NH+ + CO -> HCO+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(646) = small + (4.41e-10)
    end if

    !NH+ + H2CO -> H3CO+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(647) = small + (4.95e-10&
          *(T32)**(-0.5))
    end if

    !NH+ + H2CO -> HCO+ + NH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(648) = small + (1.82e-09&
          *(T32)**(-0.5))
    end if

    !NH+ + H2O -> H3O+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(649) = small + (1.05e-09&
          *(T32)**(-0.5))
    end if

    !NH+ + H2O -> HNO+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(650) = small + (3.5e-10&
          *(T32)**(-0.5))
    end if

    !NH+ + H2O -> NH3+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(651) = small + (1.75e-10&
          *(T32)**(-0.5))
    end if

    !NH+ + H2O -> OH + NH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(652) = small + (8.75e-10&
          *(T32)**(-0.5))
    end if

    !NH+ + HCN -> HCNH+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(653) = small + (1.8e-09&
          *(T32)**(-0.5))
    end if

    !NH+ + HCO -> H2CO+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(654) = small + (1.3e-09&
          *(T32)**(-0.5))
    end if

    !NH+ + HNC -> HCNH+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(655) = small + (1.8e-09&
          *(T32)**(-0.5))
    end if

    !NH+ + N2 -> N2H+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(656) = small + (6.5e-10)
    end if

    !NH+ + NH2 -> NH3+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(657) = small + (1.5e-09&
          *(T32)**(-0.5))
    end if

    !NH+ + NH -> NH2+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(658) = small + (1e-09&
          *(T32)**(-0.5))
    end if

    !NH+ + NO -> N2H+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(659) = small + (1.78e-10)
    end if

    !NH+ + O2 -> NO+ + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(660) = small + (2.05e-10)
    end if

    !NH+ + O2 -> O2H+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(661) = small + (1.64e-10)
    end if

    !NH+ + O -> OH+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(662) = small + (1e-09)
    end if

    !NH+ + OH -> H2O+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(663) = small + (1e-09&
          *(T32)**(-0.5))
    end if

    !NH2+ + H2CO -> H3CO+ + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(664) = small + (2.24e-09&
          *(T32)**(-0.5))
    end if

    !NH2+ + H2CO -> HCO + NH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(665) = small + (5.6e-10&
          *(T32)**(-0.5))
    end if

    !NH2+ + H2O -> H3O+ + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(666) = small + (2.76e-09&
          *(T32)**(-0.5))
    end if

    !NH2+ + H2O -> NH3+ + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(667) = small + (1e-10&
          *(T32)**(-0.5))
    end if

    !NH2+ + HCN -> HCNH+ + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(668) = small + (1.2e-09&
          *(T32)**(-0.5))
    end if

    !NH2+ + HCO -> H2CO+ + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(669) = small + (4.3e-10&
          *(T32)**(-0.5))
    end if

    !NH2+ + HNC -> HCNH+ + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(670) = small + (1.2e-09&
          *(T32)**(-0.5))
    end if

    !NH2+ + NH2 -> NH3+ + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(671) = small + (1e-09&
          *(T32)**(-0.5))
    end if

    !NH2+ + O2 -> H2NO+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(672) = small + (1.19e-10)
    end if

    !NH2+ + O2 -> HNO+ + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(673) = small + (2.1e-11)
    end if

    !NH2 + CO+ -> HCO+ + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(674) = small + (4.5e-10&
          *(T32)**(-0.5))
    end if

    !NH2 + H2CO+ -> HCO + NH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(675) = small + (8.8e-10&
          *(T32)**(-0.5))
    end if

    !NH2 + H2O+ -> NH3+ + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(676) = small + (4.9e-10&
          *(T32)**(-0.5))
    end if

    !NH2 + H3CO+ -> H2CO + NH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(677) = small + (8.8e-10&
          *(T32)**(-0.5))
    end if

    !NH2 + H3O+ -> H2O + NH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(678) = small + (9.7e-10&
          *(T32)**(-0.5))
    end if

    !NH2 + HCN+ -> CN + NH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(679) = small + (9e-10&
          *(T32)**(-0.5))
    end if

    !NH2 + HCNH+ -> HCN + NH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(680) = small + (4.45e-10&
          *(T32)**(-0.5))
    end if

    !NH2 + HCNH+ -> HNC + NH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(681) = small + (4.45e-10&
          *(T32)**(-0.5))
    end if

    !NH2 + HCO+ -> CO + NH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(682) = small + (8.9e-10&
          *(T32)**(-0.5))
    end if

    !NH2 + HNO+ -> NO + NH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(683) = small + (8.8e-10&
          *(T32)**(-0.5))
    end if

    !NH2 + N2H+ -> N2 + NH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(684) = small + (8.9e-10&
          *(T32)**(-0.5))
    end if

    !NH2 + O2H+ -> O2 + NH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(685) = small + (8.7e-10&
          *(T32)**(-0.5))
    end if

    !NH2 + OH+ -> NH3+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(686) = small + (5e-10&
          *(T32)**(-0.5))
    end if

    !NH3 + CO+ -> HCO+ + NH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(687) = small + (4.08e-11&
          *(T32)**(-0.5))
    end if

    !NH3 + HCN+ -> HCNH+ + NH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(688) = small + (8.4e-10&
          *(T32)**(-0.5))
    end if

    !NH + CH3+ -> HCNH+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(689) = small + (7.4e-10&
          *(T32)**(-0.5))
    end if

    !NH + CO+ -> HCO+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(690) = small + (3.2e-10&
          *(T32)**(-0.5))
    end if

    !NH + H2CO+ -> H3CO+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(691) = small + (6.4e-10&
          *(T32)**(-0.5))
    end if

    !NH + H2O+ -> H3O+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(692) = small + (7.1e-10&
          *(T32)**(-0.5))
    end if

    !NH + HCN+ -> CN + NH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(693) = small + (6.5e-10&
          *(T32)**(-0.5))
    end if

    !NH + HCO+ -> CO + NH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(694) = small + (6.4e-10&
          *(T32)**(-0.5))
    end if

    !NH + HNO+ -> NO + NH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(695) = small + (6.3e-10&
          *(T32)**(-0.5))
    end if

    !NH + N2H+ -> N2 + NH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(696) = small + (6.4e-10&
          *(T32)**(-0.5))
    end if

    !NH + NH2+ -> NH3+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(697) = small + (7.3e-10&
          *(T32)**(-0.5))
    end if

    !NH + O+ -> NO+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(698) = small + (3.6e-10&
          *(T32)**(-0.5))
    end if

    !NH + O2+ -> HNO+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(699) = small + (3.2e-10&
          *(T32)**(-0.5))
    end if

    !NH + O2H+ -> O2 + NH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(700) = small + (6.3e-10&
          *(T32)**(-0.5))
    end if

    !NH + OH+ -> NH2+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(701) = small + (3.6e-10&
          *(T32)**(-0.5))
    end if

    !NO + O2H+ -> O2 + HNO+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(702) = small + (7.7e-10)
    end if

    !O+ + CH3OH -> H2CO+ + H2O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(703) = small + (9.5e-11&
          *(T32)**(-0.5))
    end if

    !O+ + CH3OH -> H3CO+ + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(704) = small + (1.33e-09&
          *(T32)**(-0.5))
    end if

    !O+ + CH4 -> OH + CH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(705) = small + (1.1e-10)
    end if

    !O+ + CN -> NO+ + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(706) = small + (1e-09&
          *(T32)**(-0.5))
    end if

    !O+ + CO2 -> O2+ + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(707) = small + (9.4e-10)
    end if

    !O+ + H2CO -> HCO+ + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(708) = small + (1.4e-09&
          *(T32)**(-0.5))
    end if

    !O+ + HCN -> HCO+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(709) = small + (1.2e-09&
          *(T32)**(-0.5))
    end if

    !O+ + HCN -> NO+ + CH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(710) = small + (1.2e-09&
          *(T32)**(-0.5))
    end if

    !O+ + HCO -> CO + OH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(711) = small + (4.3e-10&
          *(T32)**(-0.5))
    end if

    !O+ + N2 -> NO+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(712) = small + (2.42e-12&
          *(T32)**(-0.21)*exp(+44.0*invT))
    end if

    !O+ + NO2 -> O2 + NO+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(713) = small + (8.3e-10)
    end if

    !O+ + OH -> O2+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(714) = small + (3.6e-10&
          *(T32)**(-0.5))
    end if

    !O2+ + CH3OH -> O2 + H3CO+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(715) = small + (5e-10&
          *(T32)**(-0.5))
    end if

    !O2H+ + CO2 -> HCO2+ + O2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(716) = small + (1.1e-09)
    end if

    !O + CH4+ -> OH + CH3+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(717) = small + (1e-09)
    end if

    !O + H2O+ -> O2+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(718) = small + (4e-11)
    end if

    !O + HCO2+ -> O2 + HCO+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(719) = small + (1e-09)
    end if

    !O + N2+ -> NO+ + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(720) = small + (1.3e-10)
    end if

    !O + N2H+ -> N2 + OH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(721) = small + (1.4e-10)
    end if

    !O + NH2+ -> HNO+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(722) = small + (7.2e-11)
    end if

    !O + NH3+ -> HNO+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(723) = small + (1e-11)
    end if

    !O + O2H+ -> O2 + OH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(724) = small + (6.2e-10)
    end if

    !O + OH+ -> O2+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(725) = small + (7.1e-10)
    end if

    !OH+ + CN -> HCN+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(726) = small + (1e-09&
          *(T32)**(-0.5))
    end if

    !OH+ + CO2 -> HCO2+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(727) = small + (1.44e-09)
    end if

    !OH+ + CO -> HCO+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(728) = small + (1.05e-09)
    end if

    !OH+ + H2CO -> H3CO+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(729) = small + (1.12e-09&
          *(T32)**(-0.5))
    end if

    !OH+ + H2O -> H3O+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(730) = small + (1.3e-09&
          *(T32)**(-0.5))
    end if

    !OH+ + HCN -> HCNH+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(731) = small + (1.2e-09&
          *(T32)**(-0.5))
    end if

    !OH+ + HCO -> CO + H2O+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(732) = small + (2.8e-10&
          *(T32)**(-0.5))
    end if

    !OH+ + HCO -> H2CO+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(733) = small + (2.8e-10&
          *(T32)**(-0.5))
    end if

    !OH+ + HNC -> HCNH+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(734) = small + (1.2e-09&
          *(T32)**(-0.5))
    end if

    !OH+ + N2 -> N2H+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(735) = small + (3.6e-10)
    end if

    !OH+ + NO -> HNO+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(736) = small + (6.11e-10)
    end if

    !OH+ + OH -> H2O+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(737) = small + (7e-10&
          *(T32)**(-0.5))
    end if

    !OH + CO+ -> HCO+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(738) = small + (3.1e-10&
          *(T32)**(-0.5))
    end if

    !OH + H2O+ -> H3O+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(739) = small + (6.9e-10&
          *(T32)**(-0.5))
    end if

    !OH + HCN+ -> CN + H2O+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(740) = small + (6.3e-10&
          *(T32)**(-0.5))
    end if

    !OH + HCO+ -> CO + H2O+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(741) = small + (6.2e-10&
          *(T32)**(-0.5))
    end if

    !OH + HCO+ -> HCO2+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(742) = small + (1e-09&
          *(T32)**(-0.5))
    end if

    !OH + HNO+ -> NO + H2O+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(743) = small + (6.2e-10&
          *(T32)**(-0.5))
    end if

    !OH + N2H+ -> N2 + H2O+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(744) = small + (6.2e-10&
          *(T32)**(-0.5))
    end if

    !OH + O2H+ -> O2 + H2O+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(745) = small + (6.1e-10&
          *(T32)**(-0.5))
    end if

    !C + CH2 -> CH + CH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(746) = small + (2.69e-12&
          *exp(-23550.0*invT))
    end if

    !C + HCO -> CO + CH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(747) = small + (1e-10)
    end if

    !C + N2 -> CN + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(748) = small + (8.69e-11&
          *exp(-22600.0*invT))
    end if

    !C + NH2 -> HCN + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(749) = small + (3.26e-11&
          *(T32)**(-0.1)*exp(+9.0*invT))
    end if

    !C + NH2 -> HNC + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(750) = small + (3.26e-11&
          *(T32)**(-0.1)*exp(+9.0*invT))
    end if

    !C + NH2 -> NH + CH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(751) = small + (9.62e-13&
          *exp(-10517.0*invT))
    end if

    !C + NH -> CN + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(752) = small + (1.2e-10)
    end if

    !C + NH -> N + CH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(753) = small + (1.73e-11&
          *(T32)**(0.5)*exp(-4000.0*invT))
    end if

    !C + NO -> CN + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(754) = small + (6e-11&
          *(T32)**(-0.16))
    end if

    !C + NO -> CO + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(755) = small + (9e-11&
          *(T32)**(-0.16))
    end if

    !C + O2 -> CO + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(756) = small + (5.56e-11&
          *(T32)**(0.41)*exp(+26.9*invT))
    end if

    !C + OCN -> CO + CN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(757) = small + (1e-10)
    end if

    !C + OH -> CO + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(758) = small + (1e-10)
    end if

    !C + OH -> O + CH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(759) = small + (2.25e-11&
          *(T32)**(0.5)*exp(-14800.0*invT))
    end if

    !CH2 + CH2 -> CH3 + CH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(760) = small + (4e-10&
          *exp(-5000.0*invT))
    end if

    !CH2 + CH4 -> CH3 + CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(761) = small + (7.13e-12&
          *exp(-5050.0*invT))
    end if

    !CH2 + CN -> HCN + CH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(762) = small + (5.3e-12&
          *exp(-2500.0*invT))
    end if

    !CH2 + H2CO -> HCO + CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(763) = small + (3.3e-13&
          *exp(-3270.0*invT))
    end if

    !CH2 + HCO -> CO + CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(764) = small + (3e-11)
    end if

    !CH2 + HNO -> NO + CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(765) = small + (1.7e-11)
    end if

    !CH2 + N2 -> HCN + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(766) = small + (8e-12&
          *exp(-18000.0*invT))
    end if

    !CH2 + NO2 -> H2CO + NO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(767) = small + (6.91e-11)
    end if

    !CH2 + NO -> H2CO + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(768) = small + (2.7e-12&
          *exp(-3500.0*invT))
    end if

    !CH2 + NO -> HCN + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(769) = small + (3.65e-12)
    end if

    !CH2 + NO -> HNCO + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(770) = small + (3.65e-12)
    end if

    !CH2 + O2 -> CO2 + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(771) = small + (2.92e-11&
          *(T32)**(-3.3)*exp(-1443.0*invT))
    end if

    !CH2 + O2 -> CO2 + H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(772) = small + (3.65e-11&
          *(T32)**(-3.3)*exp(-1443.0*invT))
    end if

    !CH2 + O2 -> CO + H2O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(773) = small + (2.48e-10&
          *(T32)**(-3.3)*exp(-1443.0*invT))
    end if

    !CH2 + O2 -> H2CO + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(774) = small + (3.65e-11&
          *(T32)**(-3.3)*exp(-1443.0*invT))
    end if

    !CH2 + O2 -> HCO + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(775) = small + (4.1e-11&
          *exp(-750.0*invT))
    end if

    !CH2 + O -> CO + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(776) = small + (8e-11)
    end if

    !CH2 + O -> CO + H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(777) = small + (1.33e-10)
    end if

    !CH2 + O -> HCO + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(778) = small + (5.01e-11)
    end if

    !CH2 + O -> OH + CH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(779) = small + (4.98e-10&
          *exp(-6000.0*invT))
    end if

    !CH2 + OH -> H2CO + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(780) = small + (3e-11)
    end if

    !CH2 + OH -> H2O + CH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(781) = small + (1.44e-11&
          *(T32)**(0.5)*exp(-3000.0*invT))
    end if

    !CH2 + OH -> O + CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(782) = small + (1.44e-11&
          *(T32)**(0.5)*exp(-3000.0*invT))
    end if

    !CH3 + CH3 -> CH4 + CH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(783) = small + (7.13e-12&
          *exp(-5052.0*invT))
    end if

    !CH3 + CN -> HCN + CH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(784) = small + (9.21e-12&
          *(T32)**(0.7)*exp(-1500.0*invT))
    end if

    !CH3 + H2CO -> HCO + CH4
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(785) = small + (1.34e-15&
          *(T32)**(5.05)*exp(-1636.0*invT))
    end if

    !CH3 + H2O -> OH + CH4
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(786) = small + (2.3e-15&
          *(T32)**(3.47)*exp(-6681.0*invT))
    end if

    !CH3 + HCO -> CO + CH4
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(787) = small + (2e-10)
    end if

    !CH3 + HNO -> NO + CH4
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(788) = small + (3.32e-12)
    end if

    !CH3 + NH2 -> CH4 + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(789) = small + (4.76e-17&
          *(T32)**(5.77)*exp(+151.0*invT))
    end if

    !CH3 + NH3 -> CH4 + NH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(790) = small + (9.55e-14&
          *exp(-4890.0*invT))
    end if

    !CH3 + NO2 -> H2CO + HNO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(791) = small + (5.41e-12)
    end if

    !CH3 + NO -> HCN + H2O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(792) = small + (4e-12&
          *exp(-7900.0*invT))
    end if

    !CH3 + O2 -> H2CO + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(793) = small + (5.64e-13&
          *exp(-4500.0*invT))
    end if

    !CH3 + O2 -> HCO + H2O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(794) = small + (1.66e-12)
    end if

    !CH3 + O2 -> O2H + CH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(795) = small + (5.3e-12&
          *exp(-34975.0*invT))
    end if

    !CH3 + O2H -> O2 + CH4
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(796) = small + (6e-12)
    end if

    !CH3 + O -> CO + H2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(797) = small + (3.6e-11&
          *exp(-202.0*invT))
    end if

    !CH3 + O -> H2CO + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(798) = small + (1.3e-10)
    end if

    !CH3 + OH -> CH4 + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(799) = small + (3.27e-14&
          *(T32)**(2.2)*exp(-2240.0*invT))
    end if

    !CH3 + OH -> H2CO + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(800) = small + (1.7e-12)
    end if

    !CH3 + OH -> H2O + CH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(801) = small + (1.2e-10&
          *exp(-1400.0*invT))
    end if

    !CH4 + CN -> HCN + CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(802) = small + (3.14e-12&
          *(T32)**(1.53)*exp(-504.0*invT))
    end if

    !CH4 + O2 -> O2H + CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(803) = small + (6.7e-11&
          *exp(-28640.0*invT))
    end if

    !CH4 + OH -> H2O + CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(804) = small + (3.77e-13&
          *(T32)**(2.42)*exp(-1162.0*invT))
    end if

    !CH + CO2 -> HCO + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(805) = small + (2.94e-13&
          *(T32)**(0.5)*exp(-3000.0*invT))
    end if

    !CH + H2CO -> HCO + CH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(806) = small + (9.21e-12&
          *(T32)**(0.7)*exp(-2000.0*invT))
    end if

    !CH + HCO -> CO + CH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(807) = small + (2.87e-12&
          *(T32)**(0.7)*exp(-500.0*invT))
    end if

    !CH + HNO -> NO + CH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(808) = small + (1.73e-11)
    end if

    !CH + N2 -> HCN + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(809) = small + (5.6e-13&
          *(T32)**(0.88)*exp(-10128.0*invT))
    end if

    !CH + N -> CN + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(810) = small + (1.66e-10&
          *(T32)**(-0.09))
    end if

    !CH + N -> NH + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(811) = small + (3.03e-11&
          *(T32)**(0.65)*exp(-1207.0*invT))
    end if

    !CH + NO -> HCN + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(812) = small + (1.2e-10&
          *(T32)**(-0.13))
    end if

    !CH + NO -> HCO + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(813) = small + (1.16e-11&
          *(T32)**(-0.13))
    end if

    !CH + NO -> OCN + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(814) = small + (3.49e-11&
          *(T32)**(-0.13))
    end if

    !CH + O2 -> CO2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(815) = small + (1.14e-11&
          *(T32)**(-0.48))
    end if

    !CH + O2 -> CO + O + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(816) = small + (1.14e-11&
          *(T32)**(-0.48))
    end if

    !CH + O2 -> CO + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(817) = small + (7.6e-12&
          *(T32)**(-0.48))
    end if

    !CH + O2 -> HCO + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(818) = small + (7.6e-12&
          *(T32)**(-0.48))
    end if

    !CH + O2H -> HCO + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(819) = small + (1.44e-11&
          *(T32)**(0.5)*exp(-3000.0*invT))
    end if

    !CH + O2H -> O2 + CH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(820) = small + (2.94e-13&
          *(T32)**(0.5)*exp(-7550.0*invT))
    end if

    !CH + O -> CO + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(821) = small + (6.02e-11&
          *(T32)**(0.1)*exp(+4.5*invT))
    end if

    !CH + O -> OH + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(822) = small + (2.52e-11&
          *exp(-2381.0*invT))
    end if

    !CH + OH -> HCO + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(823) = small + (1.44e-11&
          *(T32)**(0.5)*exp(-5000.0*invT))
    end if

    !CN + H2CO -> HCO + HCN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(824) = small + (2.6e-10&
          *(T32)**(-0.47)*exp(-826.0*invT))
    end if

    !CN + HCO -> CO + HCN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(825) = small + (1e-10)
    end if

    !CN + HNO -> NO + HCN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(826) = small + (3e-11)
    end if

    !CN + NO2 -> NO + OCN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(827) = small + (7.02e-11&
          *(T32)**(-0.27)*exp(-8.3*invT))
    end if

    !CN + NO -> N2 + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(828) = small + (1.6e-13)
    end if

    !CN + NO -> OCN + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(829) = small + (1.62e-10&
          *exp(-21205.0*invT))
    end if

    !CN + O2 -> NO + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(830) = small + (5.12e-12&
          *(T32)**(-0.49)*exp(+5.2*invT))
    end if

    !CN + O2 -> OCN + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(831) = small + (2.02e-11&
          *(T32)**(-0.19)*exp(+31.9*invT))
    end if

    !CO + HNO -> CO2 + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(832) = small + (3.32e-12&
          *exp(-6170.0*invT))
    end if

    !CO + NO2 -> CO2 + NO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(833) = small + (1.48e-10&
          *exp(-17000.0*invT))
    end if

    !CO + O2 -> CO2 + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(834) = small + (5.99e-12&
          *exp(-24075.0*invT))
    end if

    !CO + O2H -> CO2 + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(835) = small + (5.6e-10&
          *exp(-12160.0*invT))
    end if

    !H2 + C -> CH + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(836) = small + (6.64e-10&
          *exp(-11700.0*invT))
    end if

    !H2 + CH2 -> CH3 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(837) = small + (5.18e-11&
          *(T32)**(0.17)*exp(-6400.0*invT))
    end if

    !H2 + CH3 -> CH4 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(838) = small + (6.86e-14&
          *(T32)**(2.74)*exp(-4740.0*invT))
    end if

    !H2 + CH -> CH2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(839) = small + (5.46e-10&
          *exp(-1943.0*invT))
    end if

    !H2 + CN -> HCN + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(840) = small + (4.04e-13&
          *(T32)**(2.87)*exp(-820.0*invT))
    end if

    !H2 + N -> NH + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(841) = small + (1.69e-09&
          *exp(-18095.0*invT))
    end if

    !H2 + NH2 -> NH3 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(842) = small + (2.05e-15&
          *(T32)**(3.89)*exp(-1400.0*invT))
    end if

    !H2 + NH -> NH2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(843) = small + (5.96e-11&
          *exp(-7782.0*invT))
    end if

    !H2 + O2 -> O2H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(844) = small + (2.4e-10&
          *exp(-28500.0*invT))
    end if

    !H2 + O2 -> OH + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(845) = small + (3.16e-10&
          *exp(-21890.0*invT))
    end if

    !H2 + O -> OH + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(846) = small + (3.14e-13&
          *(T32)**(2.7)*exp(-3150.0*invT))
    end if

    !H2 + OH -> H2O + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(847) = small + (2.05e-12&
          *(T32)**(1.52)*exp(-1736.0*invT))
    end if

    !H + CH2 -> CH + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(848) = small + (2.2e-10)
    end if

    !H + CH3 -> CH2 + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(849) = small + (1e-10&
          *exp(-7600.0*invT))
    end if

    !H + CH4 -> CH3 + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(850) = small + (5.94e-13&
          *(T32)**(3.0)*exp(-4045.0*invT))
    end if

    !H + CH -> C + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(851) = small + (1.31e-10&
          *exp(-80.0*invT))
    end if

    !H + CO2 -> CO + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(852) = small + (3.38e-10&
          *exp(-13163.0*invT))
    end if

    !H + CO -> OH + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(853) = small + (1.1e-10&
          *(T32)**(0.5)*exp(-77700.0*invT))
    end if

    !H + H2CN -> HCN + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(854) = small + (1e-10)
    end if

    !H + H2CO -> HCO + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(855) = small + (4.85e-12&
          *(T32)**(1.9)*exp(-1379.0*invT))
    end if

    !H + H2O -> OH + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(856) = small + (1.59e-11&
          *(T32)**(1.2)*exp(-9610.0*invT))
    end if

    !H + HCN -> CN + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(857) = small + (6.2e-10&
          *exp(-12500.0*invT))
    end if

    !H + HCO -> CO + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(858) = small + (1.5e-10)
    end if

    !H + HCO -> O + CH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(859) = small + (6.61e-11&
          *exp(-51598.0*invT))
    end if

    !H + HNC -> HCN + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(860) = small + (1.14e-13&
          *(T32)**(4.23)*exp(+114.6*invT))
    end if

    !H + HNO -> NH2 + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(861) = small + (1.05e-09&
          *(T32)**(-0.3)*exp(-14730.0*invT))
    end if

    !H + HNO -> NO + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(862) = small + (4.5e-11&
          *(T32)**(0.72)*exp(-329.0*invT))
    end if

    !H + HNO -> OH + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(863) = small + (2.4e-09&
          *(T32)**(-0.5)*exp(-9010.0*invT))
    end if

    !H + NH2 -> NH + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(864) = small + (4.56e-12&
          *(T32)**(1.02)*exp(-2161.0*invT))
    end if

    !H + NH3 -> NH2 + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(865) = small + (7.8e-13&
          *(T32)**(2.4)*exp(-4990.0*invT))
    end if

    !H + NH -> N + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(866) = small + (1.73e-11&
          *(T32)**(0.5)*exp(-2400.0*invT))
    end if

    !H + NO2 -> NO + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(867) = small + (1.4e-10&
          *exp(-740.0*invT))
    end if

    !H + NO -> O + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(868) = small + (9.29e-10&
          *(T32)**(-0.1)*exp(-35220.0*invT))
    end if

    !H + NO -> OH + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(869) = small + (3.6e-10&
          *exp(-24910.0*invT))
    end if

    !H + O2 -> OH + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(870) = small + (2.61e-10&
          *exp(-8156.0*invT))
    end if

    !H + O2H -> H2O + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(871) = small + (5e-11&
          *exp(-866.0*invT))
    end if

    !H + O2H -> O2 + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(872) = small + (2.06e-11&
          *(T32)**(0.84)*exp(-277.0*invT))
    end if

    !H + O2H -> OH + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(873) = small + (1.66e-10&
          *exp(-413.0*invT))
    end if

    !H + OCN -> HCN + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(874) = small + (1.87e-11&
          *(T32)**(0.9)*exp(-2924.0*invT))
    end if

    !H + OCN -> NH + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(875) = small + (1.26e-10&
          *exp(-515.0*invT))
    end if

    !H + OCN -> OH + CN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(876) = small + (1e-10)
    end if

    !H + OH -> O + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(877) = small + (6.99e-14&
          *(T32)**(2.8)*exp(-1950.0*invT))
    end if

    !HCO + HCO -> CO + CO + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(878) = small + (3.6e-11)
    end if

    !HCO + HCO -> H2CO + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(879) = small + (3e-11)
    end if

    !HCO + HNO -> H2CO + NO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(880) = small + (1e-12&
          *exp(-1000.0*invT))
    end if

    !HCO + NO -> HNO + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(881) = small + (1.2e-11)
    end if

    !HCO + O2 -> CO2 + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(882) = small + (7.6e-13)
    end if

    !HCO + O2 -> O2H + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(883) = small + (4.64e-12&
          *(T32)**(0.7)*exp(+25.6*invT))
    end if

    !HCO + O2H -> O2 + H2CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(884) = small + (5e-11)
    end if

    !HNCO + C -> CO + HNC
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(885) = small + (1e-12)
    end if

    !N + CH2 -> HCN + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(886) = small + (3.95e-11&
          *(T32)**(0.17))
    end if

    !N + CH2 -> HNC + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(887) = small + (3.95e-11&
          *(T32)**(0.17))
    end if

    !N + CH2 -> NH + CH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(888) = small + (9.96e-13&
          *exp(-20380.0*invT))
    end if

    !N + CH3 -> H2CN + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(889) = small + (7.4e-11&
          *(T32)**(0.26)*exp(-8.4*invT))
    end if

    !N + CH3 -> HCN + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(890) = small + (1.3e-11&
          *(T32)**(0.5))
    end if

    !N + CH3 -> HCN + H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(891) = small + (3.32e-13)
    end if

    !N + CN -> N2 + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(892) = small + (1e-10&
          *(T32)**(0.18))
    end if

    !N + CO2 -> NO + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(893) = small + (3.2e-13&
          *exp(-1710.0*invT))
    end if

    !N + H2CN -> HCN + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(894) = small + (1e-10&
          *exp(-200.0*invT))
    end if

    !N + HCO -> CO + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(895) = small + (5.71e-12&
          *(T32)**(0.5)*exp(-1000.0*invT))
    end if

    !N + HCO -> HCN + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(896) = small + (1.7e-10)
    end if

    !N + HCO -> OCN + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(897) = small + (1e-10)
    end if

    !N + HNO -> NO + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(898) = small + (2.94e-12&
          *(T32)**(0.5)*exp(-1000.0*invT))
    end if

    !N + NH -> N2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(899) = small + (4.98e-11)
    end if

    !N + NO2 -> N2 + O + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(900) = small + (2.41e-12)
    end if

    !N + NO2 -> NO + NO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(901) = small + (1e-12)
    end if

    !N + NO2 -> O2 + N2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(902) = small + (1e-12)
    end if

    !N + NO -> N2 + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(903) = small + (3.38e-11&
          *(T32)**(-0.17)*exp(+2.8*invT))
    end if

    !N + O2 -> NO + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(904) = small + (2.26e-12&
          *(T32)**(0.86)*exp(-3134.0*invT))
    end if

    !N + O2H -> O2 + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(905) = small + (1.7e-13)
    end if

    !N + OH -> NO + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(906) = small + (6.05e-11&
          *(T32)**(-0.23)*exp(-14.9*invT))
    end if

    !N + OH -> O + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(907) = small + (1.88e-11&
          *(T32)**(0.1)*exp(-10700.0*invT))
    end if

    !NH2 + CH4 -> CH3 + NH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(908) = small + (2.9e-13&
          *(T32)**(2.87)*exp(-5380.0*invT))
    end if

    !NH2 + NO -> N2 + H2O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(909) = small + (4.27e-11&
          *(T32)**(-2.5)*exp(-331.0*invT))
    end if

    !NH2 + NO -> N2 + OH + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(910) = small + (1.49e-12)
    end if

    !NH2 + OH -> H2O + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(911) = small + (1.35e-12&
          *(T32)**(1.25)*exp(+43.5*invT))
    end if

    !NH2 + OH -> NH3 + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(912) = small + (2.08e-13&
          *(T32)**(0.76)*exp(-262.0*invT))
    end if

    !NH3 + CN -> HCN + NH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(913) = small + (2.75e-11&
          *(T32)**(-1.14))
    end if

    !NH + CH4 -> CH3 + NH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(914) = small + (6.63e-16&
          *(T32)**(6.13)*exp(-5895.0*invT))
    end if

    !NH + CN -> HCN + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(915) = small + (2.94e-12&
          *(T32)**(0.5)*exp(-1000.0*invT))
    end if

    !NH + H2O -> OH + NH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(916) = small + (1.83e-12&
          *(T32)**(1.6)*exp(-14090.0*invT))
    end if

    !NH + NH3 -> NH2 + NH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(917) = small + (5.25e-10&
          *exp(-13470.0*invT))
    end if

    !NH + NH -> N2 + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(918) = small + (1.7e-11)
    end if

    !NH + NH -> N2 + H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(919) = small + (1.16e-09)
    end if

    !NH + NH -> NH2 + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(920) = small + (1.81e-13&
          *(T32)**(1.8)*exp(+70.0*invT))
    end if

    !NH + NO2 -> HNO + NO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(921) = small + (2.44e-11&
          *(T32)**(-1.94)*exp(-56.9*invT))
    end if

    !NH + NO -> N2 + O + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(922) = small + (7.4e-10&
          *exp(-10540.0*invT))
    end if

    !NH + NO -> N2 + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(923) = small + (1.33e-11&
          *(T32)**(-0.78)*exp(-40.0*invT))
    end if

    !NH + O2 -> HNO + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(924) = small + (6.88e-14&
          *(T32)**(2.07)*exp(-3281.0*invT))
    end if

    !NH + O2 -> NO + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(925) = small + (2.54e-14&
          *(T32)**(1.18)*exp(-312.0*invT))
    end if

    !NH + O -> NO + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(926) = small + (6.6e-11)
    end if

    !NH + O -> OH + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(927) = small + (1.16e-11)
    end if

    !NH + OH -> H2O + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(928) = small + (3.11e-12&
          *(T32)**(1.2))
    end if

    !NH + OH -> HNO + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(929) = small + (3.32e-11)
    end if

    !NH + OH -> NH2 + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(930) = small + (2.93e-12&
          *(T32)**(0.1)*exp(-5800.0*invT))
    end if

    !NO + NO -> O2 + N2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(931) = small + (2.51e-11&
          *exp(-30653.0*invT))
    end if

    !NO + O2 -> NO2 + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(932) = small + (2.8e-12&
          *exp(-23400.0*invT))
    end if

    !NO + OCN -> N2 + CO2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(933) = small + (4.55e-11&
          *(T32)**(-1.33)*exp(-242.0*invT))
    end if

    !O2 + OCN -> CO2 + NO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(934) = small + (1.32e-12)
    end if

    !O2 + OCN -> NO2 + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(935) = small + (8.1e-11&
          *exp(-773.0*invT))
    end if

    !O + CH4 -> OH + CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(936) = small + (2.29e-12&
          *(T32)**(2.2)*exp(-3820.0*invT))
    end if

    !O + CN -> CO + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(937) = small + (2.54e-11)
    end if

    !O + CN -> NO + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(938) = small + (5.37e-11&
          *exp(-13800.0*invT))
    end if

    !O + CO2 -> O2 + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(939) = small + (2.46e-11&
          *exp(-26567.0*invT))
    end if

    !O + H2CN -> OCN + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(940) = small + (1e-10)
    end if

    !O + H2CO -> HCO + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(941) = small + (1.07e-11&
          *(T32)**(1.17)*exp(-1242.0*invT))
    end if

    !O + H2O -> OH + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(942) = small + (1.85e-11&
          *(T32)**(0.95)*exp(-8571.0*invT))
    end if

    !O + HCN -> CN + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(943) = small + (6.21e-10&
          *exp(-12439.0*invT))
    end if

    !O + HCN -> CO + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(944) = small + (7.3e-13&
          *(T32)**(1.14)*exp(-3742.0*invT))
    end if

    !O + HCN -> OCN + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(945) = small + (1.36e-12&
          *(T32)**(1.38)*exp(-3693.0*invT))
    end if

    !O + HCO -> CO2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(946) = small + (5e-11)
    end if

    !O + HCO -> CO + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(947) = small + (5e-11)
    end if

    !O + HNO -> NO2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(948) = small + (1e-12)
    end if

    !O + HNO -> NO + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(949) = small + (3.8e-11&
          *(T32)**(-0.08))
    end if

    !O + HNO -> O2 + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(950) = small + (2.94e-12&
          *(T32)**(0.5)*exp(-3500.0*invT))
    end if

    !O + N2 -> NO + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(951) = small + (2.51e-10&
          *exp(-38602.0*invT))
    end if

    !O + NH2 -> HNO + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(952) = small + (6.3e-11&
          *(T32)**(-0.1))
    end if

    !O + NH2 -> OH + NH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(953) = small + (7e-12&
          *(T32)**(-0.1))
    end if

    !O + NH3 -> OH + NH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(954) = small + (1.89e-11&
          *exp(-4003.0*invT))
    end if

    !O + NO2 -> O2 + NO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(955) = small + (9.82e-12&
          *(T32)**(-0.21)*exp(-5.2*invT))
    end if

    !O + NO -> O2 + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(956) = small + (1.18e-11&
          *exp(-20413.0*invT))
    end if

    !O + O2H -> O2 + OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(957) = small + (5.76e-11&
          *(T32)**(-0.3)*exp(-7.5*invT))
    end if

    !O + OCN -> CO + NO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(958) = small + (1e-10)
    end if

    !O + OCN -> O2 + CN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(959) = small + (4.02e-10&
          *(T32)**(-1.43)*exp(-3501.0*invT))
    end if

    !O + OH -> O2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(960) = small + (3.69e-11&
          *(T32)**(-0.27)*exp(-12.9*invT))
    end if

    !OH + CN -> HCN + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(961) = small + (1e-11&
          *exp(-1000.0*invT))
    end if

    !OH + CN -> OCN + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(962) = small + (7e-11)
    end if

    !OH + CO -> CO2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(963) = small + (2.81e-13&
          *exp(-176.0*invT))
    end if

    !OH + H2CO -> HCO + H2O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(964) = small + (7.76e-12&
          *(T32)**(0.82)*exp(+30.6*invT))
    end if

    !OH + HCN -> CN + H2O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(965) = small + (1.87e-13&
          *(T32)**(1.5)*exp(-3887.0*invT))
    end if

    !OH + HCN -> CO + NH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(966) = small + (1.07e-13&
          *exp(-5892.0*invT))
    end if

    !OH + HCO -> CO + H2O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(967) = small + (1.7e-10)
    end if

    !OH + HNO -> NO + H2O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(968) = small + (6.17e-12&
          *(T32)**(1.23)*exp(+44.3*invT))
    end if

    !OH + NH3 -> H2O + NH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(969) = small + (1.47e-13&
          *(T32)**(2.05)*exp(-7.0*invT))
    end if

    !OH + NO -> NO2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(970) = small + (5.2e-12&
          *exp(-15100.0*invT))
    end if

    !OH + O2H -> O2 + H2O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(971) = small + (8.58e-11&
          *(T32)**(-0.56)*exp(-14.8*invT))
    end if

    !OH + OH -> H2O + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(972) = small + (1.65e-12&
          *(T32)**(1.14)*exp(-50.0*invT))
    end if

    !C -> C+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(973) = small + (3.1e-10&
          *exp(-3.3*user_Av)*user_rad&
          /1.7)
    end if

    !CH+ -> C + H+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(974) = small + (3.3e-10&
          *exp(-2.9*user_Av)*user_rad&
          /1.7)
    end if

    !CH2+ -> C+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(975) = small + (4.67e-11&
          *exp(-2.2*user_Av)*user_rad&
          /1.7)
    end if

    !CH2+ -> CH+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(976) = small + (4.67e-11&
          *exp(-2.2*user_Av)*user_rad&
          /1.7)
    end if

    !CH2+ -> CH + H+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(977) = small + (4.67e-11&
          *exp(-2.2*user_Av)*user_rad&
          /1.7)
    end if

    !CH2 -> CH2+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(978) = small + (1e-09&
          *exp(-2.3*user_Av)*user_rad&
          /1.7)
    end if

    !CH2 -> CH + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(979) = small + (5.8e-10&
          *exp(-2.0*user_Av)*user_rad&
          /1.7)
    end if

    !CH3+ -> CH+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(980) = small + (1e-09&
          *exp(-1.7*user_Av)*user_rad&
          /1.7)
    end if

    !CH3+ -> CH2+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(981) = small + (1e-09&
          *exp(-1.7*user_Av)*user_rad&
          /1.7)
    end if

    !CH3 -> CH2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(982) = small + (1.35e-10&
          *exp(-2.3*user_Av)*user_rad&
          /1.7)
    end if

    !CH3 -> CH3+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(983) = small + (1e-10&
          *exp(-2.1*user_Av)*user_rad&
          /1.7)
    end if

    !CH3 -> CH + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(984) = small + (1.35e-10&
          *exp(-2.3*user_Av)*user_rad&
          /1.7)
    end if

    !CH3OH -> H2CO + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(985) = small + (7e-10&
          *exp(-2.3*user_Av)*user_rad&
          /1.7)
    end if

    !CH3OH -> H3CO+ + H + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(986) = small + (1.3e-10&
          *exp(-2.6*user_Av)*user_rad&
          /1.7)
    end if

    !CH3OH -> OH + CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(987) = small + (7e-10&
          *exp(-2.3*user_Av)*user_rad&
          /1.7)
    end if

    !CH4+ -> CH2+ + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(988) = small + (2.27e-10&
          *exp(-2.7*user_Av)*user_rad&
          /1.7)
    end if

    !CH4+ -> CH3+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(989) = small + (5.33e-11&
          *exp(-2.7*user_Av)*user_rad&
          /1.7)
    end if

    !CH4 -> CH2 + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(990) = small + (9.8e-10&
          *exp(-2.6*user_Av)*user_rad&
          /1.7)
    end if

    !CH4 -> CH3 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(991) = small + (2.2e-10&
          *exp(-2.6*user_Av)*user_rad&
          /1.7)
    end if

    !CH4 -> CH4+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(992) = small + (6.8e-12&
          *exp(-3.9*user_Av)*user_rad&
          /1.7)
    end if

    !CH4 -> CH + H2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(993) = small + (2.2e-10&
          *exp(-2.6*user_Av)*user_rad&
          /1.7)
    end if

    !CH -> C + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(994) = small + (9.2e-10&
          *exp(-1.7*user_Av)*user_rad&
          /1.7)
    end if

    !CH -> CH+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(995) = small + (7.6e-10&
          *exp(-3.3*user_Av)*user_rad&
          /1.7)
    end if

    !CN -> N + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(996) = small + (2.9e-10&
          *exp(-3.5*user_Av)*user_rad&
          /1.7)
    end if

    !CO+ -> C+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(997) = small + (1e-10&
          *exp(-2.5*user_Av)*user_rad&
          /1.7)
    end if

    !CO2 -> CO + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(998) = small + (8.9e-10&
          *exp(-3.0*user_Av)*user_rad&
          /1.7)
    end if

    !CO -> O + C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(999) = small + (2e-10&
          *exp(-3.5*user_Av)*user_rad&
          /1.7)
    end if

    !H2+ -> H+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1000) = small + (5.7e-10&
          *exp(-2.4*user_Av)*user_rad&
          /1.7)
    end if

    !H2CN -> HCN + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1001) = small + (5.48e-10&
          *exp(-2.0*user_Av)*user_rad&
          /1.7)
    end if

    !H2CO -> CO + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1002) = small + (1e-09&
          *exp(-2.2*user_Av)*user_rad&
          /1.7)
    end if

    !H2CO -> CO + H + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1003) = small + (7e-10&
          *exp(-1.7*user_Av)*user_rad&
          /1.7)
    end if

    !H2CO -> H2CO+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1004) = small + (4.7e-10&
          *exp(-2.8*user_Av)*user_rad&
          /1.7)
    end if

    !H2CO -> HCO+ + H + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1005) = small + (1.4e-11&
          *exp(-3.1*user_Av)*user_rad&
          /1.7)
    end if

    !H2O+ -> OH+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1006) = small + (1e-12&
          *exp(-2.0*user_Av)*user_rad&
          /1.7)
    end if

    !H2O -> H2O+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1007) = small + (3.1e-11&
          *exp(-3.9*user_Av)*user_rad&
          /1.7)
    end if

    !H2O -> OH + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1008) = small + (8e-10&
          *exp(-2.2*user_Av)*user_rad&
          /1.7)
    end if

    !H3+ -> H2+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1009) = small + (5e-15&
          *exp(-2.3*user_Av)*user_rad&
          /1.7)
    end if

    !H3+ -> H2 + H+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1010) = small + (5e-15&
          *exp(-1.8*user_Av)*user_rad&
          /1.7)
    end if

    !HCN -> CN + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1011) = small + (1.6e-09&
          *exp(-2.7*user_Av)*user_rad&
          /1.7)
    end if

    !HCO+ -> CO+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1012) = small + (5.4e-12&
          *exp(-3.3*user_Av)*user_rad&
          /1.7)
    end if

    !HCO -> CO + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1013) = small + (1.1e-09&
          *exp(-1.1*user_Av)*user_rad&
          /1.7)
    end if

    !HCO -> HCO+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1014) = small + (5.6e-10&
          *exp(-2.1*user_Av)*user_rad&
          /1.7)
    end if

    !HNC -> CN + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1015) = small + (1.5e-09&
          *exp(-2.1*user_Av)*user_rad&
          /1.7)
    end if

    !HNCO -> NH + CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1016) = small + (1e-09&
          *exp(-1.7*user_Av)*user_rad&
          /1.7)
    end if

    !HNO -> NO + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1017) = small + (1.7e-10&
          *exp(-0.5*user_Av)*user_rad&
          /1.7)
    end if

    !MG -> MG+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1018) = small + (7.9e-11&
          *exp(-2.1*user_Av)*user_rad&
          /1.7)
    end if

    !N2 -> N + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1019) = small + (2.3e-10&
          *exp(-3.9*user_Av)*user_rad&
          /1.7)
    end if

    !NH+ -> N + H+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1020) = small + (5.4e-11&
          *exp(-1.6*user_Av)*user_rad&
          /1.7)
    end if

    !NH2 -> NH2+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1021) = small + (1.73e-10&
          *exp(-2.6*user_Av)*user_rad&
          /1.7)
    end if

    !NH2 -> NH + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1022) = small + (7.5e-10&
          *exp(-2.0*user_Av)*user_rad&
          /1.7)
    end if

    !NH3 -> NH2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1023) = small + (9.23e-10&
          *exp(-2.1*user_Av)*user_rad&
          /1.7)
    end if

    !NH3 -> NH3+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1024) = small + (2.8e-10&
          *exp(-3.1*user_Av)*user_rad&
          /1.7)
    end if

    !NH3 -> NH + H2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1025) = small + (2.76e-10&
          *exp(-2.1*user_Av)*user_rad&
          /1.7)
    end if

    !NH -> N + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1026) = small + (5e-10&
          *exp(-2.3*user_Av)*user_rad&
          /1.7)
    end if

    !NH -> NH+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1027) = small + (1e-11&
          *exp(-2.0*user_Av)*user_rad&
          /1.7)
    end if

    !NO2 -> NO + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1028) = small + (1.4e-09&
          *exp(-2.1*user_Av)*user_rad&
          /1.7)
    end if

    !NO -> NO+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1029) = small + (2.6e-10&
          *exp(-2.9*user_Av)*user_rad&
          /1.7)
    end if

    !NO -> O + N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1030) = small + (4.7e-10&
          *exp(-2.1*user_Av)*user_rad&
          /1.7)
    end if

    !O2+ -> O+ + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1031) = small + (3.5e-11&
          *exp(-2.0*user_Av)*user_rad&
          /1.7)
    end if

    !O2 -> O2+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1032) = small + (7.6e-11&
          *exp(-3.9*user_Av)*user_rad&
          /1.7)
    end if

    !O2 -> O + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1033) = small + (7.9e-10&
          *exp(-2.1*user_Av)*user_rad&
          /1.7)
    end if

    !O2H -> O2 + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1034) = small + (3.35e-10&
          *exp(-2.1*user_Av)*user_rad&
          /1.7)
    end if

    !O2H -> OH + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1035) = small + (3.35e-10&
          *exp(-2.1*user_Av)*user_rad&
          /1.7)
    end if

    !OCN -> CN + O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1036) = small + (1e-11&
          *exp(-2.0*user_Av)*user_rad&
          /1.7)
    end if

    !OH+ -> O+ + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1037) = small + (1.1e-11&
          *exp(-3.5*user_Av)*user_rad&
          /1.7)
    end if

    !OH -> O + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1038) = small + (3.9e-10&
          *exp(-2.2*user_Av)*user_rad&
          /1.7)
    end if

    !OH -> OH+ + E
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1039) = small + (1.6e-12&
          *exp(-3.1*user_Av)*user_rad&
          /1.7)
    end if

    !C+ + N -> CN+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1040) = small + (1.08e-18&
          *(T32)**(0.07)*exp(-57.5*invT))
    end if

    !C+ + O -> CO+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1041) = small + (3.14e-18&
          *(T32)**(-0.15)*exp(-68.0*invT))
    end if

    !C + N -> CN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1042) = small + (5.72e-19&
          *(T32)**(0.37)*exp(-51.0*invT))
    end if

    !C + O+ -> CO+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1043) = small + (5e-10&
          *(T32)**(-3.7)*exp(-800.0*invT))
    end if

    !C + O -> CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1044) = small + (4.69e-19&
          *(T32)**(1.52)*exp(+50.5*invT))
    end if

    !H+ + H -> H2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1045) = small + (1.15e-18&
          *(T32)**(1.49)*exp(-228.0*invT))
    end if

    !H+ + HE -> HEH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1046) = small + (5.26e-20&
          *(T32)**(-0.51))
    end if

    !H2 + C+ -> CH2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1047) = small + (2e-16&
          *(T32)**(-1.3)*exp(-23.0*invT))
    end if

    !H2 + C -> CH2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1048) = small + (1e-17)
    end if

    !H2 + CH -> CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1049) = small + (5.09e-18&
          *(T32)**(-0.71)*exp(-11.6*invT))
    end if

    !H + C+ -> CH+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1050) = small + (1.7e-17)
    end if

    !H + C -> CH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1051) = small + (1e-17)
    end if

    !H + O -> OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1052) = small + (9.9e-19&
          *(T32)**(-0.38))
    end if

    !H + OH -> H2O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1053) = small + (5.26e-18&
          *(T32)**(-5.22)*exp(-90.0*invT))
    end if

    !N+ + N -> N2+
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1054) = small + (3.71e-18&
          *(T32)**(0.24)*exp(-26.1*invT))
    end if

    !O + O -> O2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1055) = small + (4.9e-20&
          *(T32)**(1.58))
    end if

    !C+ + E -> C
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1056) = small + (2.36e-12&
          *(T32)**(-0.29)*exp(+17.6*invT))
    end if

    !CH3+ + E -> CH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1057) = small + (1.1e-10&
          *(T32)**(-0.5))
    end if

    !H+ + E -> H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1058) = small + (3.5e-12&
          *(T32)**(-0.75))
    end if

    !H2CO+ + E -> H2CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1059) = small + (1.1e-10&
          *(T32)**(-0.7))
    end if

    !HE+ + E -> HE
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1060) = small + (5.36e-12&
          *(T32)**(-0.5))
    end if

    !MG+ + E -> MG
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1061) = small + (2.78e-12&
          *(T32)**(-0.68))
    end if

    !N+ + E -> N
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1062) = small + (3.5e-12&
          *(T32)**(-0.53)*exp(+3.2*invT))
    end if

    !O+ + E -> O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1063) = small + (3.24e-12&
          *(T32)**(-0.66))
    end if

    !CH3OH -> CH3OH_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1064) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /32.d0)*Hnuclei)
    end if

    !HNCO -> HNCO_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1065) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /43.d0)*Hnuclei)
    end if

    !OCN -> HNCO_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1066) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /42.d0)*Hnuclei)
    end if

    !HOC+ -> H2CO_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1067) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /29.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !CO -> H2CO_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1068) = small + (4.57d4&
          *0.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /28.d0)*Hnuclei)
    end if

    !H3CO+ -> CH3OH_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1069) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /31.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !C -> CH4_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1070) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /12.d0)*Hnuclei)
    end if

    !CO -> CO_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1071) = small + (4.57d4&
          *0.9*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /28.d0)*Hnuclei)
    end if

    !H2CO -> H2CO_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1072) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /30.d0)*Hnuclei)
    end if

    !CH -> CH4_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1073) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /13.d0)*Hnuclei)
    end if

    !OH -> H2O_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1074) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /17.d0)*Hnuclei)
    end if

    !NO -> NO_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1075) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /30.d0)*Hnuclei)
    end if

    !CH2 -> CH4_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1076) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /14.d0)*Hnuclei)
    end if

    !H2O -> H2O_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1077) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /18.d0)*Hnuclei)
    end if

    !CO2 -> CO2_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1078) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /44.d0)*Hnuclei)
    end if

    !CH3 -> CH4_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1079) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /15.d0)*Hnuclei)
    end if

    !CH4 -> CH4_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1080) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /16.d0)*Hnuclei)
    end if

    !HCO -> H2CO_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1081) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /29.d0)*Hnuclei)
    end if

    !N2 -> N2_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1082) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /28.d0)*Hnuclei)
    end if

    !CN -> HCN_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1083) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /26.d0)*Hnuclei)
    end if

    !C+ -> CH4_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1084) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /12.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !NH -> NH3_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1085) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /15.d0)*Hnuclei)
    end if

    !HCN -> HCN_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1086) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /27.d0)*Hnuclei)
    end if

    !NH3 -> NH3_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1087) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /17.d0)*Hnuclei)
    end if

    !N+ -> NH3_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1088) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /14.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !O+ -> H2O_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1089) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /16.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !O2+ -> O2_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1090) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /32.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !N2+ -> N2_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1091) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /28.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !CH+ -> CH4_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1092) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /13.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !NH+ -> NH3_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1093) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /15.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !OH+ -> H2O_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1094) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /17.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !CO+ -> CO_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1095) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /28.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !CN+ -> HCN_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1096) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /26.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !NO+ -> NO_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1097) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /30.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !CH2+ -> CH4_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1098) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /14.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !NH2+ -> NH3_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1099) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /16.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !H2O+ -> H2O_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1100) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /18.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !HCO+ -> H2CO_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1101) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /29.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !HCN+ -> HCN_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1102) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /27.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !NH3+ -> NH3_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1103) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /17.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !H2CO+ -> H2CO_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1104) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /30.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !CH3+ -> CH4_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1105) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /15.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !H3O+ -> H2O_DUST + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1106) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /19.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !HCO2+ -> CO2_DUST + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1107) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /45.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !CH4+ -> CH4_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1108) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /16.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !NH2 -> NH3_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1109) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /16.d0)*Hnuclei)
    end if

    !N -> NH3_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1110) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /14.d0)*Hnuclei)
    end if

    !O -> H2O_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1111) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /16.d0)*Hnuclei)
    end if

    !O2 -> O2_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1112) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /32.d0)*Hnuclei)
    end if

    !NO2 -> NO2_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1113) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /46.d0)*Hnuclei)
    end if

    !HNO -> HNO_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1114) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /31.d0)*Hnuclei)
    end if

    !HNO+ -> HNO_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1115) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /31.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !H2NO+ -> HNO_DUST + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1116) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /32.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !O2H+ -> O2H_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1117) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /33.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !H2CN -> H2CN_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1118) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /28.d0)*Hnuclei)
    end if

    !MG+ -> MG_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1119) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /24.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !MG -> MG_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1120) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /24.d0)*Hnuclei)
    end if

    !HCNH+ -> HCN_DUST + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1121) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /28.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !N2H+ -> N2_DUST + H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1122) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /29.d0)&
          *(1.0+16.71d-4/(user_gRad*Tgas))*Hnuclei)
    end if

    !O2H -> O2H_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1123) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /33.d0)*Hnuclei)
    end if

    !CO -> CH3OH_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1124) = small + (4.57d4&
          *0.1*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /28.d0)*Hnuclei)
    end if

    !HNC -> HNC_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1125) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*sqrt(Tgas&
          /27.d0)*Hnuclei)
    end if

    !E -> E_DUST
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1126) = small + (4.57d4&
          *1.0*user_gArea*freeze(Tgas)*(1.0+16.71d-4&
          /(user_gRad*Tgas))&
          *Hnuclei)
    end if

    !CH4_DUST -> CH4
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1127) = small + (( desoh2(960.d0,mantleabund)&
          *(epsilon*1.d-17*sqrTgas&
          *n(idx_H)) + descr(960.d0,mantleabund)*(4*pi*user_zeta&
          *1.64d-4*user_gArea*phi) + deuvcr(960.d0,mantleabund)&
          *(user_gArea*uvy*4.875d3*exp(-1.8*user_Av)&
          *(1d0+user_rad&
          /uvcreff)) )  / mantleabund)
    end if

    !NH3_DUST -> NH3
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1128) = small + (( desoh2(5500.d0,mantleabund)&
          *(epsilon*1.d-17*sqrTgas&
          *n(idx_H)) + descr(5500.d0,mantleabund)*(4*pi*user_zeta&
          *1.64d-4*user_gArea*phi) + deuvcr(5500.d0,mantleabund)&
          *(user_gArea*uvy*4.875d3*exp(-1.8*user_Av)&
          *(1d0+user_rad&
          /uvcreff)) )  / mantleabund)
    end if

    !H2O_DUST -> H2O
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1129) = small + (( desoh2(5600.d0,mantleabund)&
          *(epsilon*1.d-17*sqrTgas&
          *n(idx_H)) + descr(5600.d0,mantleabund)*(4*pi*user_zeta&
          *1.64d-4*user_gArea*phi) + deuvcr(5600.d0,mantleabund)&
          *(user_gArea*uvy*4.875d3*exp(-1.8*user_Av)&
          *(1d0+user_rad&
          /uvcreff)) )  / mantleabund)
    end if

    !MG_DUST -> MG
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1130) = small + (( desoh2(5300.d0,mantleabund)&
          *(epsilon*1.d-17*sqrTgas&
          *n(idx_H)) + descr(5300.d0,mantleabund)*(4*pi*user_zeta&
          *1.64d-4*user_gArea*phi) + deuvcr(5300.d0,mantleabund)&
          *(user_gArea*uvy*4.875d3*exp(-1.8*user_Av)&
          *(1d0+user_rad&
          /uvcreff)) )  / mantleabund)
    end if

    !HCN_DUST -> HCN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1131) = small + (( desoh2(3700.d0,mantleabund)&
          *(epsilon*1.d-17*sqrTgas&
          *n(idx_H)) + descr(3700.d0,mantleabund)*(4*pi*user_zeta&
          *1.64d-4*user_gArea*phi) + deuvcr(3700.d0,mantleabund)&
          *(user_gArea*uvy*4.875d3*exp(-1.8*user_Av)&
          *(1d0+user_rad&
          /uvcreff)) )  / mantleabund)
    end if

    !HNC_DUST -> HNC
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1132) = small + (( desoh2(3800.d0,mantleabund)&
          *(epsilon*1.d-17*sqrTgas&
          *n(idx_H)) + descr(3800.d0,mantleabund)*(4*pi*user_zeta&
          *1.64d-4*user_gArea*phi) + deuvcr(3800.d0,mantleabund)&
          *(user_gArea*uvy*4.875d3*exp(-1.8*user_Av)&
          *(1d0+user_rad&
          /uvcreff)) )  / mantleabund)
    end if

    !CO_DUST -> CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1133) = small + (( desoh2(1300.d0,mantleabund)&
          *(epsilon*1.d-17*sqrTgas&
          *n(idx_H)) + descr(1300.d0,mantleabund)*(4*pi*user_zeta&
          *1.64d-4*user_gArea*phi) + deuvcr(1300.d0,mantleabund)&
          *(user_gArea*uvy*4.875d3*exp(-1.8*user_Av)&
          *(1d0+user_rad&
          /uvcreff)) )  / mantleabund)
    end if

    !N2_DUST -> N2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1134) = small + (( desoh2(1100.d0,mantleabund)&
          *(epsilon*1.d-17*sqrTgas&
          *n(idx_H)) + descr(1100.d0,mantleabund)*(4*pi*user_zeta&
          *1.64d-4*user_gArea*phi) + deuvcr(1100.d0,mantleabund)&
          *(user_gArea*uvy*4.875d3*exp(-1.8*user_Av)&
          *(1d0+user_rad&
          /uvcreff)) )  / mantleabund)
    end if

    !H2CN_DUST -> H2CN
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1135) = small + (( desoh2(2400.d0,mantleabund)&
          *(epsilon*1.d-17*sqrTgas&
          *n(idx_H)) + descr(2400.d0,mantleabund)*(4*pi*user_zeta&
          *1.64d-4*user_gArea*phi) + deuvcr(2400.d0,mantleabund)&
          *(user_gArea*uvy*4.875d3*exp(-1.8*user_Av)&
          *(1d0+user_rad&
          /uvcreff)) )  / mantleabund)
    end if

    !NO_DUST -> NO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1136) = small + (( desoh2(1600.d0,mantleabund)&
          *(epsilon*1.d-17*sqrTgas&
          *n(idx_H)) + descr(1600.d0,mantleabund)*(4*pi*user_zeta&
          *1.64d-4*user_gArea*phi) + deuvcr(1600.d0,mantleabund)&
          *(user_gArea*uvy*4.875d3*exp(-1.8*user_Av)&
          *(1d0+user_rad&
          /uvcreff)) )  / mantleabund)
    end if

    !H2CO_DUST -> H2CO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1137) = small + (( desoh2(4500.d0,mantleabund)&
          *(epsilon*1.d-17*sqrTgas&
          *n(idx_H)) + descr(4500.d0,mantleabund)*(4*pi*user_zeta&
          *1.64d-4*user_gArea*phi) + deuvcr(4500.d0,mantleabund)&
          *(user_gArea*uvy*4.875d3*exp(-1.8*user_Av)&
          *(1d0+user_rad&
          /uvcreff)) )  / mantleabund)
    end if

    !HNO_DUST -> HNO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1138) = small + (( desoh2(3000.d0,mantleabund)&
          *(epsilon*1.d-17*sqrTgas&
          *n(idx_H)) + descr(3000.d0,mantleabund)*(4*pi*user_zeta&
          *1.64d-4*user_gArea*phi) + deuvcr(3000.d0,mantleabund)&
          *(user_gArea*uvy*4.875d3*exp(-1.8*user_Av)&
          *(1d0+user_rad&
          /uvcreff)) )  / mantleabund)
    end if

    !O2_DUST -> O2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1139) = small + (( desoh2(1200.d0,mantleabund)&
          *(epsilon*1.d-17*sqrTgas&
          *n(idx_H)) + descr(1200.d0,mantleabund)*(4*pi*user_zeta&
          *1.64d-4*user_gArea*phi) + deuvcr(1200.d0,mantleabund)&
          *(user_gArea*uvy*4.875d3*exp(-1.8*user_Av)&
          *(1d0+user_rad&
          /uvcreff)) )  / mantleabund)
    end if

    !CH3OH_DUST -> CH3OH
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1140) = small + (( desoh2(5000.d0,mantleabund)&
          *(epsilon*1.d-17*sqrTgas&
          *n(idx_H)) + descr(5000.d0,mantleabund)*(4*pi*user_zeta&
          *1.64d-4*user_gArea*phi) + deuvcr(5000.d0,mantleabund)&
          *(user_gArea*uvy*4.875d3*exp(-1.8*user_Av)&
          *(1d0+user_rad&
          /uvcreff)) )  / mantleabund)
    end if

    !O2H_DUST -> O2H
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1141) = small + (( desoh2(5000.d0,mantleabund)&
          *(epsilon*1.d-17*sqrTgas&
          *n(idx_H)) + descr(5000.d0,mantleabund)*(4*pi*user_zeta&
          *1.64d-4*user_gArea*phi) + deuvcr(5000.d0,mantleabund)&
          *(user_gArea*uvy*4.875d3*exp(-1.8*user_Av)&
          *(1d0+user_rad&
          /uvcreff)) )  / mantleabund)
    end if

    !HNCO_DUST -> HNCO
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1142) = small + (( desoh2(4400.d0,mantleabund)&
          *(epsilon*1.d-17*sqrTgas&
          *n(idx_H)) + descr(4400.d0,mantleabund)*(4*pi*user_zeta&
          *1.64d-4*user_gArea*phi) + deuvcr(4400.d0,mantleabund)&
          *(user_gArea*uvy*4.875d3*exp(-1.8*user_Av)&
          *(1d0+user_rad&
          /uvcreff)) )  / mantleabund)
    end if

    !CO2_DUST -> CO2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1143) = small + (( desoh2(2600.d0,mantleabund)&
          *(epsilon*1.d-17*sqrTgas&
          *n(idx_H)) + descr(2600.d0,mantleabund)*(4*pi*user_zeta&
          *1.64d-4*user_gArea*phi) + deuvcr(2600.d0,mantleabund)&
          *(user_gArea*uvy*4.875d3*exp(-1.8*user_Av)&
          *(1d0+user_rad&
          /uvcreff)) )  / mantleabund)
    end if

    !NO2_DUST -> NO2
    if(Tgas.GE.0d0 .and. Tgas.LT.10000d0) then
      k(1144) = small + (( desoh2(2400.d0,mantleabund)&
          *(epsilon*1.d-17*sqrTgas&
          *n(idx_H)) + descr(2400.d0,mantleabund)*(4*pi*user_zeta&
          *1.64d-4*user_gArea*phi) + deuvcr(2400.d0,mantleabund)&
          *(user_gArea*uvy*4.875d3*exp(-1.8*user_Av)&
          *(1d0+user_rad&
          /uvcreff)) )  / mantleabund)
    end if

    coe(:) = k(:) !set coefficients to return variable

    !!uncomment below to check coefficient values
    !kmax = 1d0
    !if(maxval(k)>kmax.or.minval(k)<0d0) then
    !   print *,"***************"
    !   do i=1,size(k)
    !      if(k(i)<0d0.or.k(i)>kmax) print *,i,k(i)
    !   end do
    !end if
  end function coe

  !*************************
  subroutine loadReactionsVerbatim()
    use krome_commons
    implicit none
    character*50::fname,line
    integer::ios,i,nunit

    fname = "reactions_verbatim.dat"

    !verbatim reactions are loaded from file
    ! to increase compilation speed
    open(newunit=nunit,file=trim(fname),status="old",iostat=ios)
    if(ios/=0) then
      print *,"ERROR: "//trim(fname)//" file not present!"
      stop
    end if

    !load reactions from file
    do i=1,nrea
      read(nunit,'(a)',iostat=ios) line
      if(ios/=0) then
        print *,"ERROR: problem reading "//trim(fname)
        stop
      end if
      reactionNames(i) = trim(line)
    end do
    close(nunit)

  end subroutine loadReactionsVerbatim

  !*******************
  !The following functions compute the recombination rate
  ! on dust for H+, He+, C+, Si+, and O+. See Weingartner&Draine 2001
  ! dust2gas_ratio, D/D_sol, default is assumed equal to Z/Z_sol
  function H_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::H_recombination_on_dust

    H_recombination_on_dust = 0d0

    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    if(psi<=0) return

    H_recombination_on_dust =  1.225d-13*dust2gas_ratio &
        /(1.d0+8.074d-6*psi**(1.378)*(1.d0+5.087d2 &
        *Tgas**(0.01586)*psi**(-0.4723-1.102d-5*log(Tgas))))

  end function H_recombination_on_dust

  !******************
  function He_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::He_recombination_on_dust

    He_recombination_on_dust = 0d0
    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    if(psi<=0) return

    He_recombination_on_dust = 5.572d-14*dust2gas_ratio&
        /(1.d0+3.185d-7*psi**(1.512)*(1.d0+5.115d3&
        *Tgas**(3.903d-7)*psi**(-0.4956-5.494d-7*log(Tgas))))

  end function He_recombination_on_dust

  !*******************
  function C_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::C_recombination_on_dust

    C_recombination_on_dust = 0d0
    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    if(psi<=0) return

    C_recombination_on_dust = 4.558d-13*dust2gas_ratio&
        /(1.d0+6.089d-3*psi**(1.128)*(1.d0+4.331d2&
        *Tgas**(0.04845)*psi**(-0.8120-1.333d-4*log(Tgas))))

  end function C_recombination_on_dust

  !******************
  function Si_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::Si_recombination_on_dust

    Si_recombination_on_dust = 0d0
    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    if(psi<=0) return

    Si_recombination_on_dust = 2.166d-14*dust2gas_ratio&
        /(1.d0+5.678d-8*psi**(1.874)*(1.d0+4.375d4&
        *Tgas**(1.635d-6)*psi**(-0.8964-7.538d-5*log(Tgas))))

  end function Si_recombination_on_dust

  !********************
  function O_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,k_H
    real*8::O_recombination_on_dust

    k_H = H_recombination_on_dust(n(:),Tgas)
    O_recombination_on_dust = 0.25d0*k_H

  end function O_recombination_on_dust

  !*********************
  !This function returns the
  ! photorate of H2 occurring in the
  ! Lyman-Werner bands following the approximation
  ! provided by Glover&Jappsen 2007. Rate in 1/s.
  !Approximation valid at low-density, it assumes H2(nu = 0).
  !It also stores the rate as a common, needed for the photoheating
  function H2_solomonLW(myflux)
    use krome_commons
    use krome_constants
    implicit none
    real*8::H2_solomonLW,myflux

    !myflux is the radiation background at E = 12.87 eV
    !should be converted to erg
    H2_solomonLW = 1.38d9*myflux*eV_to_erg

  end function H2_solomonLW

  !****************************
  !tanh smoothing function that
  ! increses when xarg increases.
  ! xpos is the position of the transition point.
  ! slope is the steepness of the curve.
  function smooth_increase(xarg,xpos,slope)
    implicit none
    real*8::smooth_increase,xarg,xpos,slope

    smooth_increase = .5d0 * (tanh(slope * (xarg - xpos)) &
        + 1d0)

  end function smooth_increase

  !****************************
  !tanh smoothing function that
  ! decreses when xarg increases.
  ! xpos is the position of the transition point.
  ! slope is the steepness of the curve.
  function smooth_decrease(xarg,xpos,slope)
    implicit none
    real*8::smooth_decrease,xarg,xpos,slope

    smooth_decrease = .5d0 * (tanh(-slope * (xarg - xpos)) &
        + 1d0)

  end function smooth_decrease

  !*********************
  !sign: return 1d0 if x>=0d0,
  ! else return -1d0
  function get_sgn(x)
    implicit none
    real*8::x,get_sgn

    get_sgn = 1d0
    if(x==0d0) return
    get_sgn = x/abs(x)

  end function get_sgn

  !*********************
  function conserve(n,ni)
    use krome_commons
    implicit none
    real*8::conserve(nspec),n(nspec),ni(nspec),no(nspec)
    real*8::ntot,nitot,factor

    no(:) = n(:)

    conserve(:) = 0d0
    conserve(:) = no(:)

  end function conserve

  !*************************
  !this subroutine changes the x(:) mass fractions of the species
  ! to force conservation according to the reference ref(:)
  subroutine conserveLin_x(x,ref)
    use krome_commons
    use krome_getphys
    implicit none
    real*8::x(nmols),ref(natoms)
    real*8::A(natoms,natoms),B(natoms),m(nspec)

    m(:) = get_mass()
    A(:,:) = 0d0
    B(:) = ref(:)

    !charge conservation
    x(idx_E) = m(idx_E)*(+ 1d0*x(idx_HCOj) / m(idx_HCOj) &
        + 1d0*x(idx_Hj) / m(idx_Hj) &
        + 1d0*x(idx_HOCj) / m(idx_HOCj) &
        + 1d0*x(idx_Cj) / m(idx_Cj) &
        + 1d0*x(idx_CH2j) / m(idx_CH2j) &
        + 1d0*x(idx_CHj) / m(idx_CHj) &
        + 1d0*x(idx_H2COj) / m(idx_H2COj) &
        + 1d0*x(idx_MGj) / m(idx_MGj) &
        + 1d0*x(idx_NH3j) / m(idx_NH3j) &
        + 1d0*x(idx_NOj) / m(idx_NOj) &
        + 1d0*x(idx_CNj) / m(idx_CNj) &
        + 1d0*x(idx_COj) / m(idx_COj) &
        + 1d0*x(idx_N2j) / m(idx_N2j) &
        + 1d0*x(idx_O2j) / m(idx_O2j) &
        + 1d0*x(idx_H2Oj) / m(idx_H2Oj) &
        + 1d0*x(idx_NH2j) / m(idx_NH2j) &
        + 1d0*x(idx_Oj) / m(idx_Oj) &
        + 1d0*x(idx_OHj) / m(idx_OHj) &
        + 1d0*x(idx_CH3j) / m(idx_CH3j) &
        + 1d0*x(idx_CH4j) / m(idx_CH4j) &
        + 1d0*x(idx_Nj) / m(idx_Nj) &
        + 1d0*x(idx_HCNj) / m(idx_HCNj) &
        + 1d0*x(idx_NHj) / m(idx_NHj) &
        + 1d0*x(idx_H2j) / m(idx_H2j) &
        + 1d0*x(idx_HEj) / m(idx_HEj) &
        + 1d0*x(idx_HNOj) / m(idx_HNOj) &
        + 1d0*x(idx_H2NOj) / m(idx_H2NOj) &
        + 1d0*x(idx_H3j) / m(idx_H3j) &
        + 1d0*x(idx_H3COj) / m(idx_H3COj) &
        + 1d0*x(idx_H3Oj) / m(idx_H3Oj) &
        + 1d0*x(idx_HCNHj) / m(idx_HCNHj) &
        + 1d0*x(idx_HCO2j) / m(idx_HCO2j) &
        + 1d0*x(idx_HEHj) / m(idx_HEHj) &
        + 1d0*x(idx_N2Hj) / m(idx_N2Hj) &
        + 1d0*x(idx_O2Hj) / m(idx_O2Hj))
    !check if charge conservation goes wrong
    if(x(idx_E)<0d0) then
      print *,"ERROR in conserveLin, electrons < 0"
      stop
    end if

  end subroutine conserveLin_x

  !***************************
  !compute the total reference mass atom type by atom type
  function conserveLinGetRef_x(x)
    use krome_commons
    use krome_getphys
    implicit none
    real*8::conserveLinGetRef_x(natoms),x(nmols)
    real*8::m(nspec)

    m(:) = get_mass()
    conserveLinGetRef_x(:) = 0d0

  end function conserveLinGetRef_x

  !***************************
  !Ref: Sasaki & Takahara (1993)
  !This function evaluate the recombination rate
  ! for H+ + e --> H + gamma and the same
  ! for D+ + e --> D + gamma
  function elec_recomb_ST93(nabund,nelec,ntot,nucleiH,Trad)
    use krome_commons
    use krome_constants
    implicit none
    real*8::nabund,nelec,Trad
    real*8::nucleiH,elec_recomb_ST93
    real*8::al,ak,rc2,r2c
    real*8::a0,b0,c0,d0,e0
    real*8::a1,b1,c1,d1,e1,f1,g1,h1
    real*8::ntot,ratio

    al = 8.227d0
    ak = 22.06d0 / (hubble  *(1d0 + phys_zredshift) &
        * sqrt(1d0 + Omega0 * phys_zredshift))
    !Rc2 evaluation
    rc2 = 8.76d-11 * (1d0 + phys_zredshift)**(-0.58)
    !R2c evaluation
    r2c = (1.80d10 * Trad)**(1.5) &
        * exp(-3.9472d4 / Trad) * rc2

    !coefficients
    a0 = nucleiH * rc2
    b0 = ak * al * nucleiH
    c0 = ak * rc2 * nucleiH * nucleiH
    d0 = r2c * exp(-1.18416d5/Trad)
    e0 = ak * r2c * nucleiH

    !polynomial terms
    a1 = -d0 * (1d0 + b0)
    b1 = d0 * (1d0 + 2d0 * b0)
    c1 = a0 + b0 * (a0 - d0)
    d1 = -a0 * b0
    e1 = a0 * c0
    f1 = 1d0 + b0 + e0
    g1 = -(b0 + e0)
    h1 = c0

    ratio = nabund / ntot

    elec_recomb_ST93 = ntot*(a1 + b1*ratio + c1*ratio**2 + d1*ratio**3 &
        + e1*ratio**4) / (f1 + g1*ratio + h1*ratio**2)

    elec_recomb_ST93 = elec_recomb_ST93 / (nabund * nelec)

  end function elec_recomb_ST93

  !********************
  subroutine load_parts()
    use krome_commons
    implicit none

  end subroutine load_parts

  !*************************
  subroutine load_part(fname,array_part,min_part,dT_part)
    character(len=*)::fname
    integer::ios,icount,i,cv
    real*8,allocatable::array_part(:),emed(:)
    real*8::min_part,dT_part,Told,array_tmp(int(1e5)),rout(2)

    open(33,file=trim(fname),status="old",iostat=ios)
    if(ios.ne.0) then
      print *,"ERROR: partition function not found"
      print *," in file "//fname
      stop
    end if

    print *,"loading partition function from "//fname
    icount = 0
    min_part = 1d99
    Told = 0d0
    do
      read(33,*,iostat=ios) rout(:)
      if(ios<0) exit
      if(ios.ne.0) cycle
      icount = icount + 1
      min_part = min(min_part,rout(1))
      array_tmp(icount) = rout(2)
      dT_part = rout(1) - Told
      Told = rout(1)
    end do
    close(33)

    allocate(array_part(icount),emed(icount))
    array_part(:) = array_tmp(1:icount)

  end subroutine load_part

  !**********************
  function troe_falloff(k0,kinf,Fc,m)
    implicit none
    real*8::troe_falloff,k0,kinf,Fc,m,rm,xexp
    rm = k0*m/kinf
    xexp = 1d0/(1d0+log10(rm)**2)
    troe_falloff = k0*m/(1d0+rm)*Fc**xexp
  end function troe_falloff

  !*************************
  function k3body(k0,kinf,Fc,nM)
    implicit none
    real*8::k3body,k0,kinf,Fc,nM
    real*8::c,n,d,Pr,xexp,F

    c = -0.4d0-0.67d0*log10(Fc)
    n = 0.75d0-1.27d0*log10(Fc)
    d = 0.14d0
    Pr = k0*nM/kinf
    xexp = (log10(Pr)+c)/(n-d*(log10(Pr)+c))
    F = 1d1**(log10(Fc)/(1d0+xexp**2))
    k3body = kinf*(Pr/(1d0+Pr)) * F

  end function k3body

  !***********************
  !see http://kida.obs.u-bordeaux1.fr/help
  function KIDA3body(ka0,kb0,kc0,kaInf,kbInf,kcInf,kaFc,kbFc,&
        kcFc,kdFc,npart,Tgas,pmin,pmax)
    implicit none
    real*8::ka0,kb0,kc0,kaInf,kbInf,kcInf,kaFc,kbFc,kcFc,kdFc
    real*8::KIDA3body,kinf,p,f,npart,Tgas,fc,fexp,invT
    real*8::k0,cc,dd,nn,pmin,pmax

    KIDA3body = 0d0

    invT = 1d0/Tgas
    k0 = ka0*(Tgas/3d2)**kb0*exp(-kc0*invT)
    kinf = kainf*(Tgas/3d2)**kbinf*exp(-kcinf*invT)

    p = k0*npart/kinf
    if(p<pmin) return
    if(p>pmax) return

    fc = (1d0-kaFc)*exp(-Tgas/kbFc) + kaFc*exp(-Tgas/kbFc) &
        + exp(-kdFc*invT)

    cc = -0.4d0 - 0.67d0 *log10(fc)
    dd = 0.14d0
    nn = 0.75d0 - 1.27d0*log10(fc)
    fexp = 1d0 + ((log10(p)+cc)/(nn-dd*(log10(p)+cc)))**2

    f = fc**(1d0/fexp)

    KIDA3body = kinf*(p/(1d0+p))*f

  end function KIDA3body

  !******************************
  !collisional ionization rate from Verner+96
  ! unit: cm3/s
  function colion_v96(Tgas,dE,P,A,X,K)
    implicit none
    real*8::colion_v96,Tgas,dE,A,X,K,U,Te,P

    Te = Tgas * 8.621738d-5 !K to eV
    U = dE / Te
    colion_v96 = A * (1d0 + P*sqrt(U)) * U**K * exp(-U) / (X+U)

  end function colion_v96

  !****************************
  !radiative recombination rates from
  ! Verner routine, standard fit, cm3/s
  function recV96(Tgas,a,b)
    implicit none
    real*8::recV96,Tgas,a,b

    recV96 = a*(1d4/Tgas)**b

  end function recV96

  !****************************
  !radiative recombination rates from
  ! Verner routine, new fit, cm3/s
  function recNewV96(Tgas,r1,r2,r3,r4)
    implicit none
    real*8::recNewV96,Tgas,r1,r2,r3,r4,tt

    tt = sqrt(Tgas/r3)
    recNewV96 = r1/(tt*(tt + 1d0)**(1.-r2) &
        * (1d0 + sqrt(Tgas/r4))**(1.+r2))

  end function recNewV96

  !****************************
  !radiative recombination rates from
  ! Verner routine, iron only, cm3/s
  function recFeV96(Tgas,r1,r2,r3)
    implicit none
    real*8::recFeV96,Tgas,r1,r2,r3,tt

    tt = sqrt(Tgas*1d-4)
    recFeV96 = r1/tt**(r2 + r3 + log10(tt))

  end function recFeV96

  !******************************
  !radiative recombination rates from Verner+96
  ! unit: cm3/s
  function radrec_v96(Tgas,a,b,T0,T1)
    implicit none
    real*8::Tgas,a,b,T0,T1,radrec_v96,iT0

    iT0 = 1d0/T0
    radrec_v96 = a/(sqrt(Tgas*iT0) + (1d0*sqrt(Tgas*iT0))**(1.-b) &
        * (1d0+sqrt(Tgas/T1))**(1+b))

  end function radrec_v96

  !*******************************
  !radiative recombination rates low-temp fit, Verner+96
  ! unit: cm3/s
  function radrec_low_v96(Tgas,a,b,c,d,f)
    implicit none
    real*8::Tgas,a,b,c,d,f,radrec_low_v96,t,invt

    t = Tgas*1d-4
    invt = 1d0/t

    radrec_low_v96 = 1d-12 * (a*invt + b + c*t + d*t**2) &
        * t**(-1.5) * exp(-f*invt)

    radrec_low_v96 = max(0d0,radrec_low_v96)

  end function radrec_low_v96

  !***************************
  !Collisional dissociation rate (cm-3/s) by Martin et al. 1996
  ! H2+H->H+H+H
  !NOTE: the use of this rate is suggested
  ! for high-density regime and in the presence of UV backgrounds.
  ! if necessary it must be included in the reaction file as
  ! H2,H,,H,H,H,,NONE,NONE,dissH2_Martin96(n,Tgas)
  function dissH2_Martin96(n,Tgas)
    use krome_commons
    use krome_getphys
    integer::i
    real*8::n(nspec),Tgas,dissH2_Martin96
    real*8::CDrates,logTv(4),k_CIDm(21,2),k_CID,invT,logT,n_c1,n_c2,n_H
    real*8::logk_h1,logk_h2,logk_l1,logk_l2,logn_c1,logn_c2,p,logk_CID
    real*8::logT2,logT3

    !k_CID = collision-induced dissociation + dissociative tunneling

    !Collisional dissociation of H2
    k_CIDm(:,1) = (/-178.4239d0, -68.42243d0, 43.20243d0, -4.633167d0, &
        69.70086d0, 40870.38d0, -23705.70d0, 128.8953d0, -53.91334d0, &
        5.315517d0, -19.73427d0, 16780.95d0, -25786.11d0, 14.82123d0, &
        -4.890915d0, 0.4749030d0, -133.8283d0, -1.164408d0, 0.8227443d0,&
        0.5864073d0, -2.056313d0/)

    !Dissociative tunneling of H2
    k_CIDm(:,2) = (/-142.7664d0, 42.70741d0, -2.027365d0, -0.2582097d0, &
        21.36094d0, 27535.31d0, -21467.79d0, 60.34928d0, -27.43096d0, &
        2.676150d0, -11.28215d0, 14254.55d0, -23125.20d0, 9.305564d0, &
        -2.464009d0, 0.1985955d0, 743.0600d0, -1.174242d0, 0.7502286d0, &
        0.2358848d0, 2.937507d0/)

    n_H  = get_Hnuclei(n(:))
    logT = log10(Tgas)
    invT = 1.0d0/Tgas
    logT2 = logT*logT
    logT3 = logT2*logT
    logTv = (/1.d0, logT, logT2, logT3/)
    k_CID = 0.d0
    do i=1,2
      logk_h1 = k_CIDm(1,i)*logTv(1) + k_CIDm(2,i)*logTv(2) + &
          k_CIDm(3,i)*logTv(3) + k_CIDm(4,i)*logTv(4) + &
          k_CIDm(5,i)*log10(1.d0+k_CIDm(6,i)*invT)
      logk_h2 = k_CIDm(7,i)*invT
      logk_l1 = k_CIDm(8,i)*logTv(1) + k_CIDm(9,i)*logTv(2) + &
          k_CIDm(10,i)*logTv(3) + k_CIDm(11,i)*log10(1.d0+k_CIDm(12,i)*invT)
      logk_l2 = k_CIDm(13,i)*invT
      logn_c1 = k_CIDm(14,i)*logTv(1) + k_CIDm(15,i)*logTv(2) &
          + k_CIDm(16,i)*logTv(3) + k_CIDm(17,i)*invT
      logn_c2 = k_CIDm(18,i) + logn_c1
      p = k_CIDm(19,i) + k_CIDm(20,i)*exp(-Tgas/1.850d3) &
          + k_CIDm(21,i)*exp(-Tgas/4.40d2)
      n_c1 = 1d1**(logn_c1)
      n_c2 = 1d1**(logn_c2)
      logk_CID = logk_h1 - (logk_h1 - logk_l1) / (1.d0 + (n_H/n_c1)**p) &
          + logk_h2 - (logk_h2 - logk_l2) / (1.d0 + (n_H/n_c2)**p)
      k_CID = k_CID + 1.d1**logk_CID
    enddo

    dissH2_Martin96 = k_CID

  end function dissH2_Martin96

  !***********************************
  subroutine init_exp_table()
    use krome_commons
    implicit none
    integer::i
    real*8::a

    do i=1,exp_table_na
      a = (i-1)*(exp_table_aMax-exp_table_aMin)/(exp_table_na-1) + exp_table_aMin
      exp_table(i) = exp(-a)
    end do

  end subroutine init_exp_table

  !***************************
  !get the index of the specie name
  function get_index(name)
    use krome_commons
    use krome_getphys
    integer::get_index,i
    character*16::names(nspec)
    character*(*)::name
    names(:) = get_names()
    get_index = -1 !default index
    !loop on species to found the specie named name
    do i=1,nspec
      !when found store and break loop
      if(trim(names(i))== trim(name)) then
        get_index = i !store index
        exit
      end if
    end do

    !error if species not found
    if(get_index<0) then
      print *,"ERROR: can't find the index of ",name
      stop
    end if

  end function get_index

  !*****************************
  !computes revers kinetics from reaction and
  ! product indexes
  ! k_rev = k_for * revKc
  ! Note that reaction constant revKc is calculated with
  ! reactants and products from reverse reaction
  function revKc(Tgas,ridx,pidx)
    use krome_constants
    use krome_commons
    implicit none
    real*8::revKc,Tgas,dgibss,stoichiometricChange
    integer::ridx(:),pidx(:),i

    ! when considering forward reaction:
    ! Kc = (P°)**(p+p-r-r) * exp(-dGibss_forward°)
    ! where ° means at standard conditions of
    ! P° = 1 bar = (kb*T/1e6) dyn/cm^2 (cgs)
    ! when considering reverse:
    ! 1/Kc = revKc = (kb*T/1e6)**(p+p-r-r) * exp(-dGibss_reverse°)
    ! kb*T/1e6 is to go from 1 atm pressure to number density cm^-3
    ! When not at standard pressure this does not change:
    ! revKc = P**(p+p-r-r) *exp(-dGibss_reverse° - (p+p-r-r)*ln(P/P°))
    !       = (P°)**(p+p-r-r) * exp(-dGibss_reverse°)

    dgibss = 0.d0 ! Gibbs free energy/(R*T)
    stoichiometricChange = 0d0

    do i=1,size(pidx)
      dgibss = dgibss + revHS(Tgas,pidx(i))
      stoichiometricChange = stoichiometricChange + 1
    end do

    do i=1,size(ridx)
      dgibss = dgibss - revHS(Tgas,ridx(i))
      stoichiometricChange = stoichiometricChange - 1
    end do

    revKc = (boltzmann_erg * Tgas * 1e-6)**(-stoichiometricChange)&
        * exp(-dgibss)

  end function revKc

  !*****************************
  !compute H-S for species with index idx
  ! when temperature is Tgas
  function revHS(Tgas,idx)
    use krome_commons
    use krome_constants
    use krome_fit
    real*8::revHS,Tgas,Tgas2,Tgas3,Tgas4,invT,lnT,H,S
    real*8::Tnist,Tnist2,Tnist3,Tnist4,invTnist,invTnist2,lnTnist
    real*8::p1_nasa(92,7), p2_nasa(92,7), Tlim_nasa(92,3), p(7)
    real*8::p1_nist(92,7), p2_nist(92,7), Tlim_nist(92,3)
    integer::idx

    p(:) = 0.d0
    p1_nasa(:,:) = 0.d0
    p2_nasa(:,:) = 0.d0
    Tlim_nasa(:,:) = 0.d0
    p1_nist(:,:) = 0.d0
    p2_nist(:,:) = 0.d0
    Tlim_nist(:,:) = 0.d0
    Tgas2 = Tgas * Tgas
    Tgas3 = Tgas2 * Tgas
    Tgas4 = Tgas3 * Tgas
    invT = 1d0/Tgas
    lnT = log(Tgas)
    ! NIST polynomials are quite differernt
    ! it doesn't like easy stuff...
    Tnist = Tgas * 1.d-3
    Tnist2 = Tnist * Tnist
    Tnist3 = Tnist2 * Tnist
    Tnist4 = Tnist3 * Tnist2
    invTnist = 1d0/Tnist
    invTnist2 = invTnist * invTnist
    lnTnist = log(Tnist)

    p1_nasa(idx_CH,:)  = (/3.4897583d0,&
        0.0003243216d0,&
        -1.6899751d-06,&
        3.162842d-09,&
        -1.4061803d-12,&
        70660.755d0,&
        2.0842841d0/)
    p1_nasa(idx_O,:)  = (/3.1682671d0,&
        -0.00327931884d0,&
        6.64306396d-06,&
        -6.12806624d-09,&
        2.11265971d-12,&
        29122.2592d0,&
        2.05193346d0/)
    p1_nasa(idx_H2,:)  = (/2.34433112d0,&
        0.00798052075d0,&
        -1.9478151d-05,&
        2.01572094d-08,&
        -7.37611761d-12,&
        -917.935173d0,&
        0.683010238d0/)
    p1_nasa(idx_C,:)  = (/2.5542395d0,&
        -0.00032153772d0,&
        7.3379223d-07,&
        -7.3223487d-10,&
        2.6652144d-13,&
        85442.681d0,&
        4.5313085d0/)
    p1_nasa(idx_H,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        25473.66d0,&
        -0.44668285d0/)
    p1_nasa(idx_H2O,:)  = (/4.1986352d0,&
        -0.0020364017d0,&
        6.5203416d-06,&
        -5.4879269d-09,&
        1.771968d-12,&
        -30293.726d0,&
        -0.84900901d0/)
    p1_nasa(idx_OH,:)  = (/3.99198424d0,&
        -0.00240106655d0,&
        4.61664033d-06,&
        -3.87916306d-09,&
        1.36319502d-12,&
        3368.89836d0,&
        -0.103998477d0/)
    p1_nasa(idx_O2,:)  = (/3.78245636d0,&
        -0.00299673416d0,&
        9.84730201d-06,&
        -9.68129509d-09,&
        3.24372837d-12,&
        -1063.94356d0,&
        3.65767573d0/)
    p1_nasa(idx_CH2,:)  = (/3.84261832d0,&
        -7.36676871d-06,&
        6.16970693d-06,&
        -6.96689962d-09,&
        2.64620979d-12,&
        45863.1528d0,&
        1.2758447d0/)
    p1_nasa(idx_H2CO,:)  = (/4.65733258d0,&
        -0.00953742306d0,&
        4.04679152d-05,&
        -4.45317569d-08,&
        1.64761516d-11,&
        13861.5127d0,&
        1.97860732d0/)
    p1_nasa(idx_HCO,:)  = (/4.36380907d0,&
        -0.00535204137d0,&
        2.31954508d-05,&
        -2.6610904d-08,&
        1.02711962d-11,&
        25010.8717d0,&
        2.98106307d0/)
    p1_nasa(idx_NO,:)  = (/4.21859896d0,&
        -0.00463988124d0,&
        1.10443049d-05,&
        -9.34055507d-09,&
        2.80554874d-12,&
        9845.09964d0,&
        2.28061001d0/)
    p1_nasa(idx_CO,:)  = (/3.5795335d0,&
        -0.00061035369d0,&
        1.0168143d-06,&
        9.0700586d-10,&
        -9.0442449d-13,&
        -14344.086d0,&
        3.5084093d0/)
    p1_nasa(idx_N2,:)  = (/3.53100528d0,&
        -0.000123660988d0,&
        -5.02999433d-07,&
        2.43530612d-09,&
        -1.40881235d-12,&
        -1046.97628d0,&
        2.96747038d0/)
    p1_nasa(idx_CH3,:)  = (/3.6571797d0,&
        0.0021265979d0,&
        5.4583883d-06,&
        -6.6181003d-09,&
        2.4657074d-12,&
        16422.716d0,&
        1.6735354d0/)
    p1_nasa(idx_CH4,:)  = (/5.14825732d0,&
        -0.013700241d0,&
        4.93749414d-05,&
        -4.91952339d-08,&
        1.70097299d-11,&
        -10245.3222d0,&
        -4.63322726d0/)
    p1_nasa(idx_N,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        56104.638d0,&
        4.1939088d0/)
    p1_nasa(idx_HE,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        -745.375d0,&
        0.928723974d0/)
    p1_nasa(idx_CO2,:)  = (/2.356813d0,&
        0.0089841299d0,&
        -7.1220632d-06,&
        2.4573008d-09,&
        -1.4288548d-13,&
        -48371.971d0,&
        9.9009035d0/)
    p1_nasa(idx_NO2,:)  = (/3.9440312d0,&
        -0.001585429d0,&
        1.6657812d-05,&
        -2.0475426d-08,&
        7.8350564d-12,&
        2896.618d0,&
        6.3119919d0/)
    p1_nasa(idx_Hj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        184021.488d0,&
        -1.14064664d0/)
    p1_nasa(idx_Cj,:)  = (/2.61332254d0,&
        -0.000540148065d0,&
        1.03037233d-06,&
        -8.90092552d-10,&
        2.88500586d-13,&
        216862.274d0,&
        3.8345479d0/)
    p1_nasa(idx_COj,:)  = (/3.77061642d0,&
        -0.00201773246d0,&
        4.61081738d-06,&
        -2.99175463d-09,&
        6.06065045d-13,&
        149006.795d0,&
        3.38129783d0/)
    p1_nasa(idx_O2j,:)  = (/4.61017167d0,&
        -0.00635951952d0,&
        1.42425624d-05,&
        -1.20997923d-08,&
        3.70956878d-12,&
        139742.229d0,&
        -0.201326941d0/)
    p1_nasa(idx_H2Oj,:)  = (/4.02465912d0,&
        -0.00108851414d0,&
        5.13576558d-06,&
        -4.40027838d-09,&
        1.40726746d-12,&
        116895.616d0,&
        0.699968812d0/)
    p1_nasa(idx_Oj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        187935.284d0,&
        4.39337676d0/)
    p1_nasa(idx_OHj,:)  = (/3.50502572d0,&
        0.000241313747d0,&
        -1.42200948d-06,&
        2.64780232d-09,&
        -1.17038711d-12,&
        155210.676d0,&
        1.97907627d0/)
    p1_nasa(idx_H2j,:)  = (/3.77256072d0,&
        -0.0019574659d0,&
        4.54812047d-06,&
        -2.82152141d-09,&
        5.33969209d-13,&
        178694.654d0,&
        -3.96609192d0/)
    p1_nasa(idx_HEj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        285323.374d0,&
        1.62166556d0/)
    p1_nasa(idx_H3j,:)  = (/4.1795698d0,&
        -0.000868875627d0,&
        -1.09017371d-07,&
        4.13349766d-09,&
        -2.37877027d-12,&
        132635.537d0,&
        -5.838001d0/)
    p1_nasa(idx_H3Oj,:)  = (/3.79295251d0,&
        -0.000910852723d0,&
        1.16363521d-05,&
        -1.21364865d-08,&
        4.26159624d-12,&
        71402.7518d0,&
        1.47156927d0/)
    p2_nasa(idx_CH,:)  = (/2.5209369d0,&
        0.0017653639d0,&
        -4.614766d-07,&
        5.9289675d-11,&
        -3.3474501d-15,&
        70994.878d0,&
        7.4051829d0/)
    p2_nasa(idx_O,:)  = (/2.54363697d0,&
        -2.73162486d-05,&
        -4.1902952d-09,&
        4.95481845d-12,&
        -4.79553694d-16,&
        29226.012d0,&
        4.92229457d0/)
    p2_nasa(idx_H2,:)  = (/2.93286575d0,&
        0.000826608026d0,&
        -1.46402364d-07,&
        1.54100414d-11,&
        -6.888048d-16,&
        -813.065581d0,&
        -1.02432865d0/)
    p2_nasa(idx_C,:)  = (/2.605583d0,&
        -0.00019593434d0,&
        1.0673722d-07,&
        -1.642394d-11,&
        8.187058d-16,&
        85411.742d0,&
        4.1923868d0/)
    p2_nasa(idx_H,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        25473.66d0,&
        -0.44668285d0/)
    p2_nasa(idx_H2O,:)  = (/2.6770389d0,&
        0.0029731816d0,&
        -7.7376889d-07,&
        9.4433514d-11,&
        -4.2689991d-15,&
        -29885.894d0,&
        6.88255d0/)
    p2_nasa(idx_OH,:)  = (/2.83853033d0,&
        0.00110741289d0,&
        -2.94000209d-07,&
        4.20698729d-11,&
        -2.4228989d-15,&
        3697.80808d0,&
        5.84494652d0/)
    p2_nasa(idx_O2,:)  = (/3.66096065d0,&
        0.000656365811d0,&
        -1.41149627d-07,&
        2.05797935d-11,&
        -1.29913436d-15,&
        -1215.97718d0,&
        3.41536279d0/)
    p2_nasa(idx_CH2,:)  = (/3.11049513d0,&
        0.00373779517d0,&
        -1.37371977d-06,&
        2.23054839d-10,&
        -1.33567178d-14,&
        45971.5953d0,&
        4.62796405d0/)
    p2_nasa(idx_H2CO,:)  = (/3.65237205d0,&
        0.0055580706d0,&
        -1.97617181d-06,&
        3.16823378d-10,&
        -1.88747598d-14,&
        13553.6156d0,&
        4.2214084d0/)
    p2_nasa(idx_HCO,:)  = (/4.23892214d0,&
        0.0019657617d0,&
        -3.82075171d-07,&
        4.80137647d-11,&
        -3.11176347d-15,&
        24726.1645d0,&
        1.99698242d0/)
    p2_nasa(idx_NO,:)  = (/3.26071234d0,&
        0.00119101135d0,&
        -4.29122646d-07,&
        6.94481463d-11,&
        -4.03295681d-15,&
        9921.43132d0,&
        6.36900518d0/)
    p2_nasa(idx_CO,:)  = (/3.0484859d0,&
        0.0013517281d0,&
        -4.8579405d-07,&
        7.8853644d-11,&
        -4.6980746d-15,&
        -14266.117d0,&
        6.0170977d0/)
    p2_nasa(idx_N2,:)  = (/2.95257637d0,&
        0.0013969004d0,&
        -4.92631603d-07,&
        7.86010195d-11,&
        -4.60755204d-15,&
        -923.948688d0,&
        5.87188762d0/)
    p2_nasa(idx_CH3,:)  = (/2.9781206d0,&
        0.005797852d0,&
        -1.97558d-06,&
        3.072979d-10,&
        -1.7917416d-14,&
        16509.513d0,&
        4.7224799d0/)
    p2_nasa(idx_CH4,:)  = (/1.911786d0,&
        0.0096026796d0,&
        -3.38387841d-06,&
        5.3879724d-10,&
        -3.19306807d-14,&
        -10099.2136d0,&
        8.48241861d0/)
    p2_nasa(idx_N,:)  = (/2.4159429d0,&
        0.00017489065d0,&
        -1.1902369d-07,&
        3.0226244d-11,&
        -2.0360983d-15,&
        56133.775d0,&
        4.6496095d0/)
    p2_nasa(idx_HE,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        -745.375d0,&
        0.928723974d0/)
    p2_nasa(idx_CO2,:)  = (/4.6365111d0,&
        0.0027414569d0,&
        -9.9589759d-07,&
        1.6038666d-10,&
        -9.1619857d-15,&
        -49024.904d0,&
        -1.9348955d0/)
    p2_nasa(idx_NO2,:)  = (/4.884754d0,&
        0.0021723955d0,&
        -8.2806909d-07,&
        1.574751d-10,&
        -1.0510895d-14,&
        2316.4982d0,&
        -0.11741695d0/)
    p2_nasa(idx_Hj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        184021.488d0,&
        -1.14064664d0/)
    p2_nasa(idx_Cj,:)  = (/2.50827618d0,&
        -1.04354146d-05,&
        5.16160809d-09,&
        -1.14187475d-12,&
        9.43539946d-17,&
        216879.645d0,&
        4.3188599d0/)
    p2_nasa(idx_COj,:)  = (/2.93062935d0,&
        0.00156033262d0,&
        -6.16246355d-07,&
        1.09957336d-10,&
        -6.66119284d-15,&
        149147.222d0,&
        7.3384673d0/)
    p2_nasa(idx_O2j,:)  = (/3.31675922d0,&
        0.00111522244d0,&
        -3.83492556d-07,&
        5.72784687d-11,&
        -2.77648381d-15,&
        139876.823d0,&
        5.44726469d0/)
    p2_nasa(idx_H2Oj,:)  = (/3.31570445d0,&
        0.00210648746d0,&
        -3.76341515d-07,&
        3.47525972d-11,&
        -1.70335643d-15,&
        117017.475d0,&
        4.03220514d0/)
    p2_nasa(idx_Oj,:)  = (/2.48542028d0,&
        2.56978695d-05,&
        -1.28833378d-08,&
        1.65525487d-12,&
        1.09933344d-16,&
        187940.874d0,&
        4.47425446d0/)
    p2_nasa(idx_OHj,:)  = (/2.68358996d0,&
        0.00157006435d0,&
        -5.39972815d-07,&
        9.37643877d-11,&
        -5.70068067d-15,&
        155479.296d0,&
        6.44375894d0/)
    p2_nasa(idx_H2j,:)  = (/3.44204765d0,&
        0.000599083239d0,&
        6.69133685d-08,&
        -3.43574373d-11,&
        1.97626599d-15,&
        178650.236d0,&
        -2.79499055d0/)
    p2_nasa(idx_HEj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        285323.374d0,&
        1.62166556d0/)
    p2_nasa(idx_H3j,:)  = (/2.01435718d0,&
        0.00415925769d0,&
        -1.42664877d-06,&
        2.22372739d-10,&
        -1.29346518d-14,&
        133230.507d0,&
        5.46168967d0/)
    p2_nasa(idx_H3Oj,:)  = (/2.49647765d0,&
        0.0057284484d0,&
        -1.83953239d-06,&
        2.73577348d-10,&
        -1.54093917d-14,&
        71624.4227d0,&
        7.45850493d0/)
    Tlim_nasa(idx_CH,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_O,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H2,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_C,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H2O,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_OH,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_O2,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_CH2,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H2CO,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_HCO,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_NO,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_CO,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_N2,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_CH3,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_CH4,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_N,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_HE,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_CO2,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_NO2,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_Hj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_Cj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_COj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_O2j,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H2Oj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_Oj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_OHj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H2j,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_HEj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H3j,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H3Oj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)

    ! pick NASA data if present for species
    if (Tlim_nasa(idx,2) /= 0.d0) then
      !select set of NASA polynomials using temperature
      if(Tlim_nasa(idx,1).le.Tgas .and. Tgas.le.Tlim_nasa(idx,2)) then
        p(:) = p1_nasa(idx,:)

      else if(Tlim_nasa(idx,2)<Tgas .and. Tgas.le.Tlim_nasa(idx,3)) then
        p(:) = p2_nasa(idx,:)

        ! currently no option when Tgas not in Tlim range p(:) = 0
      end if

      !compute NASA polynomials for enthalpy and enthropy (unitless)
      H = p(1) + p(2)*0.5d0*Tgas + p(3)*Tgas2/3.d0 + p(4)*Tgas3*0.25d0 + &
          p(5)*Tgas4*0.2d0 + p(6)*invT
      S = p(1)*lnT + p(2)*Tgas + p(3)*Tgas2*0.5d0 + p(4)*Tgas3/3.d0 + &
          p(5)*Tgas4*0.25d0 + p(7)

      revHS = H - S

      ! else pick NIST data (if present)
    else if (Tlim_nist(idx,2) /= 0.d0) then
      if (Tlim_nist(idx,1) < Tgas .and. Tgas < Tlim_nist(idx,2)) then
        p(:) = p1_nist(idx,:)

      else if (Tlim_nist(idx,2) < Tgas .and. Tgas < Tlim_nist(idx,3)) then
        p(:) = p2_nist(idx,:)

        ! currently no option when Tgas not in Tlim range p(:) = 0
      end if

      !compute NIST polynomials for enthalpy and enthropy
      ! H in (kJ/mol)
      H = p(1)*Tnist + p(2)*0.5d0*Tnist2 + p(3)*Tnist3/3.d0 + p(4)*Tnist4*0.25d0&
          - p(5)*invTnist + p(6)
      !  Unitsless
      H = H / (Rgas_kJ * Tgas)

      ! S in (J/mol*K)
      S = p(1)*lnTnist + p(2)*Tnist + p(3)*Tnist2*0.5d0 + p(4)*Tnist3/3.d0&
          - p(5)*invTnist2*0.5d0 + p(7)
      !  Unitless. Note: do not use Tnist
      S = S / Rgas_J

      revHS = H - S

      ! return zero is no data exists
    else
      print *, "No thermochemical data of species index", idx
      revHS = 0.d0

    end if

  end function revHS

  !******************************
  subroutine print_best_flux(n,Tgas,nbestin)
    !print the first nbestin fluxes
    use krome_commons
    use krome_getphys
    implicit none
    real*8::n(nspec),Tgas,flux(nrea)
    integer::nbest,idx(nrea),i,nbestin
    character*50::name(nrea)

    nbest = min(nbestin,nrea) !cannot exceed the number of reactions

    flux(:) = get_flux(n(:),Tgas) !get fluxes
    name(:) = get_rnames() !get reaction names

    !call the sorting algorithm (bubblesort)
    idx(:) = idx_sort(flux(:))

    !print to screen
    print *,"***************"
    do i=1,nbest
      print '(I8,a1,a50,E17.8)',idx(i)," ",name(idx(i)),flux(idx(i))
    end do

  end subroutine print_best_flux

  !******************************
  subroutine print_best_flux_frac(n,Tgas,frac)
    !print the first nbestin fluxes
    use krome_commons
    use krome_getphys
    implicit none
    real*8::n(nspec),Tgas,flux(nrea),frac
    integer::idx(nrea),i
    character*50::name(nrea)

    if(frac>1d0) then
      print *,"ERROR: fraction in krome_print_best_flux should be <=1!"
      stop
    end if

    flux(:) = get_flux(n(:),Tgas) !get fluxes
    name(:) = get_rnames() !get reaction names

    !call the sorting algorithm (bubblesort)
    idx(:) = idx_sort(flux(:))

    !print to screen
    print *,"***************"
    do i=1,nrea
      if(flux(idx(i))<flux(idx(1))*frac) exit
      print '(I8,a1,a50,E17.8)',idx(i)," ",name(idx(i)),flux(idx(i))
    end do

  end subroutine print_best_flux_frac

  !******************************
  subroutine print_best_flux_spec(n,Tgas,nbestin,idx_found)
    !print the first nbestin fluxes for the reactions
    ! that contains the species with index idx_found
    use krome_commons
    use krome_getphys
    implicit none
    real*8::n(nspec),Tgas,flux(nrea),maxflux
    integer::nbest,idx(nrea),i,nbestin,idx_found
    character*50::name(nrea)
    logical::found

    nbest = min(nbestin,nrea) !cannot exceed the number of reactions
    maxflux = 0d0
    flux(:) = get_flux(n(:),Tgas) !get fluxes
    name(:) = get_rnames() !get reaction names
    do i=1,nrea
      found = .false.
      if(arr_r1(i) == idx_found) found = .true.
      if(arr_r2(i) == idx_found) found = .true.
      if(arr_p1(i) == idx_found) found = .true.
      if(arr_p2(i) == idx_found) found = .true.
      if(arr_p3(i) == idx_found) found = .true.
      if(arr_p4(i) == idx_found) found = .true.
      maxflux = max(maxflux,flux(i))
      if(.not.found) flux(i) = 0d0
    end do

    !call the sorting algorithm (bubblesort)
    idx(:) = idx_sort(flux(:))

    !print to screen
    print *,"***************"
    do i=1,nbest
      print '(I8,a1,a50,2E17.8)',idx(i)," ",name(idx(i)),flux(idx(i)),&
          flux(idx(i))/maxflux
    end do

  end subroutine print_best_flux_spec

  !*****************************
  function idx_sort(fin)
    !sorting algorithm: requires an array of real values fin
    ! and returns the sorted index list. descending.
    ! bubblesort: not very efficient, replace with what you prefer
    implicit none
    real*8::fin(:),f(size(fin)),ftmp
    integer::idx_sort(size(fin)),n,itmp,i
    logical::found

    f(:) = fin(:) !copy to local

    n = size(f)
    !init indexes
    do i=1,n
      idx_sort(i) = i
    end do

    !loop to sort
    do
      found = .false. !swapped something flag
      do i=2,n
        !> for descending, < for ascending
        if(f(i)>f(i-1)) then
          found = .true.
          !swap real value
          ftmp = f(i)
          f(i) = f(i-1)
          f(i-1) = ftmp
          !swap index
          itmp = idx_sort(i)
          idx_sort(i) = idx_sort(i-1)
          idx_sort(i-1) = itmp
        end if
      end do
      !if nothing swapped exit
      if(.not.found) exit
    end do

  end function idx_sort

  !******************************
  function get_flux(n,Tgas)
    !get the flux k*n*n*... of the rates
    use krome_commons
    implicit none
    integer::i
    integer::r1,r2
    real*8::get_flux(nrea),n(nspec),k(nrea),rrmax,Tgas

    k(:) = coe(n(:))
    rrmax = 0.d0
    n(idx_dummy) = 1.d0
    n(idx_g) = 1.d0
    n(idx_CR) = 1.d0
    do i=1,nrea
      r1 = arr_r1(i)
      r2 = arr_r2(i)
      arr_flux(i) = k(i)*n(r1)*n(r2)
    end do
    get_flux(:) = arr_flux(:)

  end function get_flux

  !*****************************
  subroutine load_arrays()
    !load the array containing reactants
    ! and product index
    use krome_commons

    arr_r1(1:1144) = (/2,55,6,6,6,6,6,6,6,8,8,8,8,8,57,57,57,57&
        ,57,57,57,7,7,7,7,59,59,59,59,58,12,12,12,12,12,12,12,12,12&
        ,72,72,72,73,73,73,23,2,2,2,2,2,2,2,2,2,2,64,64,64,64,64,64&
        ,18,65,65,65,65,19,55,55,55,55,55,55,55,55,55,55,55,55,55,55&
        ,55,55,77,77,77,77,77,77,77,77,77,77,77,77,77,77,77,6,13,68&
        ,68,68,68,68,9,9,9,8,8,8,8,8,8,75,75,5,5,14,14,78,78,78,78,78&
        ,78,78,78,15,15,15,15,15,74,74,74,74,74,74,74,74,74,74,74,74&
        ,74,74,74,66,66,66,66,24,76,76,76,76,76,69,69,69,21,21,21,21&
        ,21,21,62,62,62,16,16,16,16,16,16,25,25,25,25,17,17,17,70,70&
        ,70,70,70,70,70,70,70,3,3,3,71,71,71,71,71,71,10,10,10,7,19,6&
        ,6,6,8,26,24,3,7,59,12,12,22,22,22,28,28,23,2,18,29,19,30,13&
        ,9,8,5,14,14,4,31,27,26,15,20,24,21,21,16,16,16,25,25,32,17&
        ,17,11,11,33,3,34,10,59,58,58,58,72,72,72,73,73,64,65,77,60&
        ,60,60,60,80,80,68,68,68,81,81,82,82,82,82,82,83,83,83,83,75&
        ,84,84,84,54,85,85,85,79,56,86,66,87,87,76,69,69,62,62,63,67&
        ,88,71,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,7,7,7,7,7&
        ,7,7,7,7,7,7,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59&
        ,59,59,59,59,58,58,58,58,58,58,12,12,12,12,12,12,12,12,12,12&
        ,12,12,12,12,12,12,72,72,72,72,72,72,72,73,73,73,73,23,23,23&
        ,23,23,23,23,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,64,64,64,18&
        ,18,65,19,19,19,19,55,55,55,55,55,55,55,55,55,55,55,55,55,77&
        ,77,77,77,77,77,77,77,77,77,77,77,77,77,77,77,77,77,77,6,6,6&
        ,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,60,60,13,13,13,68,68,68,68,68&
        ,68,68,9,9,9,9,9,9,9,9,9,9,9,9,81,81,81,81,81,81,81,81,81,81&
        ,81,81,81,81,81,81,81,81,81,81,81,81,81,81,81,83,83,83,8,8,8&
        ,8,8,75,75,75,75,75,75,75,5,5,5,5,5,5,84,84,54,54,14,14,14,14&
        ,14,4,4,4,4,4,4,79,78,78,78,78,78,78,78,78,78,78,78,78,78,78&
        ,78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,78&
        ,78,78,78,78,78,78,78,78,78,78,78,78,74,74,74,74,74,74,74,74&
        ,74,74,74,74,74,74,74,74,74,74,66,66,20,20,87,87,24,24,24,24&
        ,24,24,24,24,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76&
        ,76,76,76,76,76,76,69,69,69,69,69,69,69,69,69,69,21,21,21,21&
        ,21,21,21,21,21,21,21,21,21,16,16,25,25,25,25,25,25,25,25,25&
        ,25,25,25,25,17,70,70,70,70,70,70,70,70,70,70,70,70,67,88,3,3&
        ,3,3,3,3,3,3,3,71,71,71,71,71,71,71,71,71,71,71,71,10,10,10&
        ,10,10,10,10,10,7,7,7,7,7,7,7,7,7,7,7,7,7,7,12,12,12,12,12,12&
        ,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,22,22,22&
        ,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,23,23,23,2,2&
        ,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,18,18,18,18,18,18,18,18,19&
        ,19,19,19,6,6,6,6,6,6,6,6,6,6,6,6,8,8,8,8,8,8,8,8,8,8,8,8,8,8&
        ,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,14,14,14,14,14,14,14,31,24&
        ,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24&
        ,24,21,21,21,21,21,16,25,25,25,25,25,25,25,25,25,25,25,25,25&
        ,25,25,25,25,17,17,17,11,11,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3&
        ,3,3,3,3,3,3,3,3,10,10,10,10,10,10,10,10,10,10,10,10,7,59,58&
        ,58,58,12,12,72,72,22,22,22,28,28,28,73,73,23,23,23,23,2,2,18&
        ,65,29,19,77,30,13,13,13,13,68,9,9,81,81,5,54,14,14,4,31,27&
        ,15,20,76,21,21,16,16,16,25,25,32,17,17,67,11,11,33,33,34,71&
        ,10,10,57,57,7,7,7,55,55,6,6,6,8,8,8,8,74,3,57,72,55,60,78,61&
        ,74,70,28,31,34,56,19,82,7,19,13,2,10,17,12,9,29,22,23,14,20&
        ,18,57,25,5,16,74,70,67,66,59,76,71,65,64,63,58,69,68,54,75&
        ,62,60,72,83,85,73,21,24,3,11,32,27,79,80,88,30,61,15,84,87&
        ,33,19,4,1,38,45,40,51,44,52,39,43,50,41,37,48,46,35,49,36,42&
        ,47/)
    arr_r2(1:1144) = (/3,4,2,6,9,56,11,10,1,2,6,9,11,10&
        ,12,2,13,14,15,16,17,64,65,66,67,14,15,16,17,17,64,65,60,68&
        ,66,69,70,67,71,14,15,17,13,16,11,65,64,65,60,68,74,66,69,70&
        ,67,71,19,13,5,14,17,11,66,13,14,17,11,66,12,22,23,2,13,9,5&
        ,14,15,21,16,25,17,11,3,10,12,23,2,18,19,13,9,5,14,21,16,25&
        ,17,11,10,78,67,13,14,15,17,11,65,75,66,64,65,77,75,78,70,17&
        ,11,65,66,60,67,7,23,2,13,9,20,16,11,60,54,66,63,67,12,23,18&
        ,19,13,9,5,14,15,21,16,25,17,11,10,13,14,17,11,66,13,9,16,17&
        ,11,14,16,17,64,65,68,66,67,71,14,15,17,65,60,68,75,66,67,64&
        ,65,66,70,60,79,67,23,19,13,9,14,21,16,11,10,64,65,66,13,9,14&
        ,16,17,11,64,65,66,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,1,1&
        ,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1&
        ,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,28,28,29,13,13&
        ,9,9,14,21,16,25,11,11,34,10,68,83,75,54,85,79,87,76,67,88,71&
        ,28,28,29,13,13,13,9,9,9,5,14,4,24,21,25,11,11,11,3,10,29,13&
        ,9,14,11,3,65,60,68,83,75,84,84,54,79,87,76,69,62,67,88,71,28&
        ,13,14,11,3,3,10,29,19,13,9,65,60,68,75,66,66,71,65,60,68,82&
        ,83,75,84,84,54,79,74,87,76,69,70,67,88,71,13,14,11,79,88,13&
        ,85,79,87,88,12,28,28,28,23,29,13,13,14,14,31,27,32,7,12,23,2&
        ,18,29,19,6,13,9,14,26,20,24,25,17,11,3,10,57,59,58,64,65,65&
        ,68,75,78,86,74,66,76,76,69,70,88,71,13,11,79,67,88,19,13,9,5&
        ,14,14,4,64,64,65,60,82,75,54,85,79,66,87,88,7,12,22,28,2,18&
        ,29,19,19,13,9,5,14,4,27,15,20,21,25,32,17,11,3,3,10,13,5,4&
        ,59,58,72,73,86,29,19,13,5,14,14,4,60,82,54,79,87,88,13,13,13&
        ,14,60,79,87,67,88,60,82,54,79,87,88,29,12,12,22,28,28,23,23&
        ,23,23,2,18,18,29,29,29,29,19,13,13,13,9,9,5,5,5,5,14,14,14,4&
        ,4,4,27,27,20,21,21,16,16,25,17,17,11,34,34,10,28,28,28,28,23&
        ,23,23,29,19,13,13,14,16,16,25,17,11,11,13,14,79,88,29,13,58&
        ,64,68,68,76,69,67,71,18,29,29,29,19,13,13,9,9,9,9,5,14,4,20&
        ,21,25,17,11,11,3,10,13,13,9,9,5,14,4,21,11,11,65,60,68,82,83&
        ,75,84,84,54,79,87,88,71,65,75,72,65,60,68,75,54,79,87,69,70&
        ,67,88,71,88,28,28,23,18,29,13,5,5,14,20,32,10,28,29,73,68,85&
        ,66,87,69,62,88,71,18,29,19,13,9,5,14,14,4,20,17,10,65,68,75&
        ,54,54,79,87,88,12,14,20,21,21,21,25,25,17,17,11,34,10,10,12&
        ,23,18,13,14,27,20,32,17,17,17,11,11,11,11,11,3,3,3,3,10,10&
        ,10,22,18,13,9,14,27,21,16,32,17,11,11,11,33,3,3,10,10,10,18&
        ,11,10,29,13,14,27,20,24,24,17,17,17,11,11,11,11,33,33,3,3,10&
        ,13,14,27,32,17,17,11,11,27,32,11,33,7,12,22,2,18,24,21,25,11&
        ,11,3,10,12,22,23,2,29,19,30,13,9,5,14,14,4,27,27,27,21,16,25&
        ,32,17,17,11,33,33,33,34,34,34,10,14,14,27,17,11,11,33,7,12&
        ,12,12,22,22,22,18,29,30,14,14,14,27,25,32,32,32,17,11,33,10&
        ,10,23,17,17,10,10,18,23,18,9,16,25,25,25,32,17,17,11,11,3,3&
        ,10,10,10,17,11,34,34,34,23,18,18,29,30,13,9,5,5,5,14,14,27&
        ,27,27,20,21,21,16,32,17,33,34,34,10,18,18,19,13,5,5,14,27,16&
        ,17,33,10,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,24,3,24,70,3,8,26,57,7,2,57,7&
        ,3,10,24,3,1,1,1,1,1,1,1,1,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92/)
    arr_p1(1:1144) = (/54,5,7&
        ,6,10,54,3,3,8,7,8,10,3,3,58,59,60,54,61,62,63,18,19,20,11,54&
        ,61,62,63,63,18,19,13,9,20,21,3,11,10,54,61,63,60,62,67,19,18&
        ,19,13,9,24,20,21,3,11,10,65,60,75,54,63,67,20,60,54,63,67,20&
        ,58,72,73,59,60,68,75,54,61,69,62,76,63,67,70,71,58,73,59,64&
        ,65,60,68,75,54,69,62,76,63,67,71,26,11,60,54,61,63,67,19,5&
        ,20,18,19,6,5,26,3,63,67,19,20,13,11,57,73,59,60,68,66,62,67&
        ,13,14,20,17,11,58,73,64,65,60,68,75,54,61,69,62,76,63,67,71&
        ,60,54,63,67,20,60,68,62,63,67,54,62,63,18,19,9,20,11,10,54&
        ,61,63,19,13,9,5,20,11,18,19,20,3,13,27,11,73,65,60,68,54,69&
        ,62,67,71,18,19,20,60,68,54,62,63,67,18,19,20,57,65,55,77,8&
        ,55,78,74,70,57,57,58,2,12,72,2,13,10,12,7,24,19,3,5,19,10,55&
        ,18,19,54,18,25,17,78,61,24,74,69,25,21,62,25,24,76,17,63,3&
        ,67,3,11,70,18,3,7,7,7,2,12,2,2,12,22,24,3,8,12,19,19,14,27&
        ,17,3,3,10,6,8,12,2,19,13,14,9,3,10,10,18,18,5,4,19,29,19,19&
        ,17,19,26,24,20,24,24,24,25,21,25,3,3,11,3,82,14,65,19,54,54&
        ,56,19,75,75,64,65,19,65,65,10,54,18,19,29,17,20,24,65,11,3&
        ,13,82,54,19,82,54,60,83,54,84,19,84,64,75,64,65,54,14,65,65&
        ,60,54,82,19,54,54,54,14,10,9,18,5,4,19,17,20,72,72,21,60,11&
        ,3,82,54,19,82,60,54,60,85,54,82,83,54,82,83,84,20,20,83,54&
        ,14,10,13,9,18,5,4,19,17,64,20,58,25,65,54,11,3,54,19,63,17&
        ,11,54,29,17,54,11,59,72,82,54,72,54,65,54,65,19,69,63,63,59&
        ,72,72,58,75,85,54,81,54,83,19,86,87,76,69,79,88,71,68,59,58&
        ,72,75,54,56,83,84,26,26,76,87,24,69,62,71,11,68,82,33,82,11&
        ,11,54,82,83,84,19,60,84,75,54,54,14,13,18,19,29,17,87,20,11&
        ,59,72,73,72,58,75,85,54,56,82,83,84,60,84,80,61,87,62,69,63&
        ,79,88,68,71,68,82,84,84,57,59,58,72,26,85,54,82,84,60,84,84&
        ,14,13,84,17,84,11,82,82,82,60,82,60,60,88,11,14,13,84,17,84&
        ,11,85,57,59,59,71,10,59,58,72,22,57,74,24,65,19,67,11,3,65&
        ,54,3,71,10,64,74,24,24,65,19,3,64,24,76,63,17,74,74,76,76,69&
        ,74,70,3,70,64,18,70,60,82,63,17,72,75,84,17,63,54,63,19,87&
        ,69,66,66,63,17,54,87,17,11,85,82,75,66,79,63,66,87,63,63,75&
        ,85,79,63,54,82,54,83,79,62,10,84,60,84,87,62,69,87,63,88,71&
        ,68,82,14,83,62,84,60,84,62,80,79,54,14,62,13,9,18,5,4,19,17&
        ,20,11,62,54,84,84,54,82,83,18,19,17,20,62,63,79,11,69,11,60&
        ,82,10,63,67,54,54,63,19,63,11,67,11,85,10,67,11,63,20,79,79&
        ,11,67,75,85,54,82,83,84,19,60,84,87,79,68,54,83,18,19,85,17&
        ,20,11,2,19,18,5,4,25,18,24,18,19,19,19,19,3,22,22,5,14,19,17&
        ,5,13,13,5,31,29,29,19,13,14,19,19,14,10,13,9,3,23,5,14,10,19&
        ,17,23,23,13,5,13,14,33,11,19,13,23,13,9,5,33,9,14,14,19,17,5&
        ,18,25,5,14,34,29,19,19,14,14,11,19,10,14,14,19,17,17,20,34&
        ,17,34,29,29,29,29,2,22,23,12,5,25,16,21,33,10,10,9,2,12,22,7&
        ,19,10,5,14,10,18,19,3,5,21,17,10,25,21,24,17,3,10,10,9,11,10&
        ,5,25,10,3,19,13,13,27,29,33,11,19,5,4,25,30,5,5,20,17,5,19,5&
        ,34,17,20,20,17,11,20,17,11,17,3,22,20,20,9,16,5,22,5,10,21&
        ,20,20,21,27,20,20,27,17,17,10,9,27,21,11,32,20,29,32,10,19&
        ,17,11,34,14,10,18,19,34,29,19,32,17,11,17,27,10,10,11,11,11&
        ,19,11,11,5,34,29,14,18,19,19,17,9,32,11,9,57,7,57,59,2,58,2&
        ,59,58,12,72,2,13,82,10,58,72,12,22,73,2,7,59,24,57,19,3,55,5&
        ,19,19,60,54,71,68,10,77,6,18,65,19,54,18,25,17,61,24,24,69&
        ,25,21,62,25,24,76,17,63,3,70,67,3,11,10,18,70,3,71,64,65,18&
        ,65,19,77,86,58,12,22,59,2,10,9,66,11,7,22,8,13,26,15,24,3,35&
        ,36,36,37,37,35,38,39,37,38,40,41,38,40,42,38,38,37,43,44,38&
        ,45,44,45,45,40,46,43,38,45,40,39,44,41,38,45,40,37,44,45,37&
        ,38,40,42,38,45,45,40,46,47,48,48,48,49,50,51,51,44,43,49,35&
        ,52,53,23,16,9,15,5,4,19,20,30,17,13,27,11,28,33,31,29&
        ,32/)
    arr_p2(1:1144) = (/1,55,6,8,6,6,3,6,8,8,8,8,3,8,7,7&
        ,7,7,7,7,7,57,57,57,57,2,2,2,2,12,58,58,58,58,58,58,58,58,58&
        ,22,22,22,23,23,23,73,59,59,59,59,59,59,59,59,59,59,18,18,18&
        ,18,18,18,64,19,19,19,19,65,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,6&
        ,6,6,6,6,6,6,6,6,6,6,6,6,6,6,77,60,9,9,9,9,9,68,68,68,55,55&
        ,55,55,55,55,5,5,75,75,54,54,26,26,26,26,26,26,26,26,61,61,61&
        ,61,61,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,20,20,20&
        ,20,74,25,25,25,25,25,21,21,21,69,69,69,69,69,69,16,16,16,62&
        ,62,62,62,62,62,76,76,76,76,63,63,63,3,3,3,3,3,3,3,3,3,70,70&
        ,70,10,10,10,10,10,10,71,71,71,1,1,8,1,8,1,1,1,1,1,8,1,8,8,1&
        ,6,6,22,6,8,7,3,7,8,6,8,1,8,8,1,8,19,8,1,1,24,1,1,8,8,1,6,8,1&
        ,3,1,24,1,3,8,1,3,8,8,6,8,8,8,6,8,8,8,7,7,8,3,6,8,8,8,6,6,8,8&
        ,8,8,10,9,6,8,8,8,6,6,8,8,8,8,8,8,8,3,10,8,8,8,24,8,25,8,8,8&
        ,8,8,24,3,8,8,2,72,19,58,2,8,8,59,8,6,8,3,70,18,8,59,6,59,59&
        ,59,59,59,59,3,59,59,72,12,19,72,7,12,8,7,6,7,58,7,8,6,6,10,3&
        ,70,8,6,19,22,8,72,10,8,2,72,72,72,72,72,72,72,72,72,24,25,72&
        ,3,72,72,23,23,73,3,8,6,6,22,22,22,22,22,22,22,22,58,72,12,7&
        ,58,58,58,58,58,58,58,58,58,8,58,24,58,8,3,58,58,5,75,19,75&
        ,75,14,54,54,20,54,6,9,6,6,6,3,6,6,6,77,19,6,10,8,8,6,8,8,8,8&
        ,8,6,8,81,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,55,81,8,8,81,8,8,8&
        ,81,8,14,54,17,54,82,10,10,10,10,83,10,10,10,25,10,83,83,83&
        ,83,83,83,10,83,83,6,6,6,9,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,10,6&
        ,6,8,6,6,9,9,9,6,6,6,6,77,18,18,18,18,18,19,18,84,84,19,84,20&
        ,84,5,4,19,19,19,17,20,19,60,84,84,19,84,20,84,17,26,26,26,22&
        ,72,26,26,26,26,26,7,57,3,70,7,57,57,26,26,58,26,26,26,2,57&
        ,59,26,86,59,26,57,7,26,26,24,26,26,26,26,26,24,74,3,3,70,26&
        ,25,25,22,72,24,6,8,65,7,25,12,76,6,25,8,3,3,70,20,19,87,87&
        ,20,20,8,7,8,6,8,8,3,8,24,24,19,14,24,24,21,24,6,3,69,24,24&
        ,24,24,24,24,3,10,24,24,24,25,62,25,10,25,25,25,25,3,10,25,62&
        ,10,62,62,62,62,62,62,62,62,62,3,21,21,6,24,24,24,69,69,69,69&
        ,24,8,3,69,3,79,9,10,72,7,19,10,24,2,71,24,63,8,82,11,72,6,54&
        ,24,71,8,6,71,8,3,3,3,3,3,3,68,3,3,3,3,3,3,3,68,68,8,68,68,68&
        ,2,2,24,8,8,2,8,2,3,24,3,18,8,2,2,22,2,22,22,22,25,17,24,10,8&
        ,6,8,9,3,10,6,8,8,2,8,2,22,12,12,23,23,23,23,25,21,27,9,10,9&
        ,12,23,6,8,3,6,12,22,22,22,19,12,12,12,24,8,7,3,24,8,8,3,10,3&
        ,10,12,8,7,8,5,5,5,34,19,24,19,3,25,17,3,10,8,8,8,8,8,8,8,8,8&
        ,10,8,8,6,6,6,6,10,7,6,6,6,6,6,12,8,3,6,25,6,6,6,10,25,24,3,3&
        ,6,10,3,19,18,6,19,19,17,19,10,19,13,4,8,8,2,8,6,8,7,19,25,25&
        ,3,8,25,8,3,17,20,3,3,25,8,25,16,9,10,25,3,21,21,24,21,21,6,8&
        ,24,17,3,10,3,10,8,24,24,8,3,20,3,29,17,19,22,24,7,19,6,10,10&
        ,10,25,8,8,10,8,10,25,24,8,25,21,17,24,10,17,18,8,3,8,8,9,9&
        ,21,9,9,21,8,9,3,1,55,6,8,55,1,8,6,8,8,1,6,6,8,22,6,8,6,8,1,6&
        ,8,1,7,3,3,7,8,8,6,8,1,8,8,1,8,8,55,8,8,8,1,8,19,8,1,24,55,1&
        ,8,8,1,6,8,1,3,1,24,3,1,3,8,3,3,8,8,1,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,8,8,92&
        ,92,92,92,92,92,92,92,8,92,92,92,92,8,8,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92/)
    arr_p3(1:1144) = (/92,92,8,8,8,92,6,8,1,8,8,8,8,8,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,1,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,8,92,92,92,8,8,92,92,92,92,92,92,8,92,92,92&
        ,92,8,92,92,8,92,92,8,92,8,92,8,92,8,92,8,92,92,92,92,8,92,92&
        ,92,92,92,92,92,92,8,92,92,8,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,6,8,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,6,92,92,8,92,92,92,92,92,92,92,92&
        ,8,92,92,92,92,92,8,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,8,92,92,92,92,92,92,92,92,92,92,92,92,8,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,6,92,92,92,92,92,92,92,92,92,92,92,8,92,92,92,6,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,6,8,6&
        ,26,26,6,6,8,55,8,26,26,26,26,26,26,26,6,8,26,8,55,8,26,26,26&
        ,8,92,26,8,26,26,8,55,26,6,8,6,8,8,26,26,26,26,26,8,8,92,8,8&
        ,8,8,8,92,92,92,92,92,92,92,92,92,92,92,8,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,8,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,8,92,92&
        ,92,92,8,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,8,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,8,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,6,92,92,92,92,92,92,92,92,92,92,92,92,8,92,92,92,92&
        ,92,92,92,92,3,92,92,92,92,92,92,92,92,92,8,92,92,92,92,92,92&
        ,92,92,8,92,92,8,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,1,92,92,92,92,92,92,8,92,92,92,92,92&
        ,92,92,92,92,8,92,1,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92/)
    arr_p4(1:1144) = (/92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,8,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,8&
        ,92,92,92,92,92,8,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92&
        ,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92/)

  end subroutine load_arrays

  ! ************************************
  ! solves linear least squares
  subroutine llsq(n, x, y, a, b)

    !****************************************************
    !
    !! LLSQ solves a linear least squares problem matching a line to data.
    !
    !  Discussion:
    !
    !    A formula for a line of the form Y = A * X + B is sought, which
    !    will minimize the root-mean-square error to N data points
    !    ( X(I), Y(I) );
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    07 March 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    In: N, the number of data values.
    !
    !    In: X(N), Y(N), the coordinates of the data points.
    !
    !    Out: A, B, the slope and Y-intercept of the
    !    least-squares approximant to the data.
    !
    implicit none
    integer,intent(in)::n
    real*8,intent(out)::a, b
    real*8,intent(in)::x(n), y(n)
    real*8::bot, top, xbar, ybar

    ! special case
    if(n == 1) then
      a = 0d0
      b = y(1)
      return
    end if

    ! average X and Y
    xbar = sum(x) / n
    ybar = sum(y) / n

    ! compute beta
    top = dot_product(x(:) - xbar, y(:) - ybar)
    bot = dot_product(x(:) - xbar, x(:) - xbar)

    ! if top is zero a is zero
    if(top==0d0) then
      a = 0d0
    else
      a = top / bot
    end if

    b = ybar - a * xbar

  end subroutine llsq

end module krome_subs
