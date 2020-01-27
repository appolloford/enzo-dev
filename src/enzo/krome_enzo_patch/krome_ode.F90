
!############### MODULE ##############
module krome_ode
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

  subroutine fex(neq,tt,nin,dn)
    use krome_commons
    use krome_constants
    use krome_subs
    use krome_cooling
    use krome_heating
    use krome_tabs
    use krome_photo
    use krome_gadiab
    use krome_getphys
    use krome_phfuncs
    use krome_fit
    use krome_user_commons
    implicit none
    integer::neq,idust
    real*8::tt,dn(neq),n(neq),k(nrea),krome_gamma
    real*8::gamma,Tgas,vgas,ntot,nH2dust,nd,nin(neq)
    real*8::rr
    integer::i,r1,r2,p1,p2,p3,p4

    n(:) = nin(:)

    nH2dust = 0.d0
    n(idx_CR) = 1.d0
    n(idx_g)  = 1.d0
    n(idx_dummy) = 1.d0

    dn(:) = 0.d0 !initialize differentials
    n(idx_Tgas) = max(n(idx_tgas),2.73d0)
    n(idx_Tgas) = min(n(idx_tgas),1d9)
    Tgas = n(idx_Tgas) !get temperature

    k(:) = coe_tab(n(:)) !compute coefficients

    n(idx_dummy) = 1.d0
    n(idx_g) = 1.d0
    n(idx_CR) = 1.d0
    do i=1,nrea
      r1 = arr_r1(i)
      r2 = arr_r2(i)
      p1 = arr_p1(i)
      p2 = arr_p2(i)
      p3 = arr_p3(i)
      p4 = arr_p4(i)
      rr = k(i)*n(r1)*n(r2)
      dn(r1) = dn(r1) - rr
      dn(r2) = dn(r2) - rr
      dn(p1) = dn(p1) + rr
      dn(p2) = dn(p2) + rr
      dn(p3) = dn(p3) + rr
      dn(p4) = dn(p4) + rr
    end do

    dn(idx_H) = dn(idx_H) - 2.0*( 1.d-17*sqrt(Tgas)*n(idx_H)*get_Hnuclei(n(:)) - h2d(n(:))*n(idx_H2) )
    dn(idx_H2) = dn(idx_H2) + 1.d-17*sqrt(Tgas)*n(idx_H)*get_Hnuclei(n(:)) - h2d(n(:))*n(idx_H2)

    last_coe(:) = k(:)

  end subroutine fex

  !***************************
  subroutine jes(neq, tt, n, j, ian, jan, pdj)
    use krome_commons
    use krome_subs
    use krome_tabs
    use krome_cooling
    use krome_heating
    use krome_constants
    use krome_gadiab
    use krome_getphys
    implicit none
    integer::neq, j, ian, jan, r1, r2, p1, p2, p3, i
    real*8::tt, n(neq), pdj(neq), dr1, dr2, kk,k(nrea),Tgas
    real*8::nn(neq),dn0,dn1,dnn,nH2dust,dn(neq),krome_gamma

    nH2dust = 0.d0
    Tgas = n(idx_Tgas)

    k(:) = last_coe(:) !get rate coefficients

    if(j==1) then
    elseif(j==1) then
      pdj(1) =  &
          -k(310)*n(idx_NOj)  &
          -k(1058)*n(idx_Hj)  &
          -k(1063)*n(idx_Oj)  &
          -k(291)*n(idx_HCNj)  &
          -k(299)*n(idx_HNOj)  &
          -k(286)*n(idx_H3COj)  &
          -k(292)*n(idx_HCNHj)  &
          -k(297)*n(idx_HCO2j)  &
          -k(270)*n(idx_H2j)  &
          -k(275)*n(idx_H2NOj)  &
          -k(308)*n(idx_NH3j)  &
          -k(263)*n(idx_CH3j)  &
          -k(1059)*n(idx_H2COj)  &
          -k(264)*n(idx_CH3j)  &
          -k(271)*n(idx_H2COj)  &
          -k(296)*n(idx_HCO2j)  &
          -k(301)*n(idx_HEHj)  &
          -k(302)*n(idx_N2j)  &
          -k(260)*n(idx_CH2j)  &
          -k(278)*n(idx_H2Oj)  &
          -k(1061)*n(idx_MGj)  &
          -k(273)*n(idx_H2COj)  &
          -k(284)*n(idx_H3COj)  &
          -k(287)*n(idx_H3Oj)  &
          -k(9)*n(idx_H2)  &
          -k(1062)*n(idx_Nj)  &
          -k(312)*n(idx_O2Hj)  &
          -k(289)*n(idx_H3Oj)  &
          -k(306)*n(idx_NH2j)  &
          -k(288)*n(idx_H3Oj)  &
          +k(9)*n(idx_H2)  &
          -k(1126)  &
          -k(309)*n(idx_NH3j)  &
          -k(262)*n(idx_CH2j)  &
          -k(277)*n(idx_H2Oj)  &
          -k(282)*n(idx_H3COj)  &
          -k(280)*n(idx_H3j)  &
          -k(313)*n(idx_OHj)  &
          -k(295)*n(idx_HCOj)  &
          -k(266)*n(idx_CH4j)  &
          -k(279)*n(idx_H2Oj)  &
          -k(285)*n(idx_H3COj)  &
          -k(1057)*n(idx_CH3j)  &
          -k(290)*n(idx_H3Oj)  &
          -k(305)*n(idx_NHj)  &
          -k(261)*n(idx_CH2j)  &
          -k(304)*n(idx_N2Hj)  &
          -k(300)*n(idx_HOCj)  &
          -k(293)*n(idx_HCNHj)  &
          -k(307)*n(idx_NH2j)  &
          -k(268)*n(idx_CNj)  &
          -k(1060)*n(idx_HEj)  &
          -k(276)*n(idx_H2NOj)  &
          -k(259)*n(idx_CHj)  &
          -k(283)*n(idx_H3COj)  &
          -k(267)*n(idx_CH4j)  &
          -k(269)*n(idx_COj)  &
          -k(274)*n(idx_H2COj)  &
          -k(294)*n(idx_HCNHj)  &
          -k(311)*n(idx_O2j)  &
          -k(1056)*n(idx_Cj)  &
          -k(272)*n(idx_H2COj)  &
          -k(281)*n(idx_H3j)  &
          -k(265)*n(idx_CH3j)  &
          -k(298)*n(idx_HCO2j)  &
          -k(303)*n(idx_N2Hj)
      pdj(2) =  &
          +k(262)*n(idx_CH2j)  &
          +k(265)*n(idx_CH3j)  &
          +k(264)*n(idx_CH3j)  &
          +k(283)*n(idx_H3COj)
      pdj(3) =  &
          +k(269)*n(idx_COj)  &
          +k(313)*n(idx_OHj)  &
          +k(297)*n(idx_HCO2j)  &
          +k(288)*n(idx_H3Oj)  &
          +k(1063)*n(idx_Oj)  &
          +k(271)*n(idx_H2COj)  &
          +k(277)*n(idx_H2Oj)  &
          +2.d0*k(311)*n(idx_O2j)  &
          +k(278)*n(idx_H2Oj)  &
          +k(310)*n(idx_NOj)
      pdj(4) =  &
          +k(294)*n(idx_HCNHj)
      pdj(5) =  &
          +k(293)*n(idx_HCNHj)
      pdj(6) =  &
          +k(276)*n(idx_H2NOj)  &
          +k(260)*n(idx_CH2j)  &
          +k(280)*n(idx_H3j)  &
          +k(288)*n(idx_H3Oj)  &
          +k(272)*n(idx_H2COj)  &
          -k(9)*n(idx_H2)  &
          +k(284)*n(idx_H3COj)  &
          +k(277)*n(idx_H2Oj)  &
          +k(289)*n(idx_H3Oj)  &
          +k(264)*n(idx_CH3j)
      pdj(7) =  &
          +k(260)*n(idx_CH2j)  &
          +k(269)*n(idx_COj)  &
          +k(1056)*n(idx_Cj)  &
          +k(261)*n(idx_CH2j)  &
          +k(268)*n(idx_CNj)  &
          +k(259)*n(idx_CHj)
      pdj(8) =  &
          +2.d0*k(265)*n(idx_CH3j)  &
          +k(1058)*n(idx_Hj)  &
          +2.d0*k(261)*n(idx_CH2j)  &
          +k(267)*n(idx_CH4j)  &
          +k(275)*n(idx_H2NOj)  &
          +k(279)*n(idx_H2Oj)  &
          +k(308)*n(idx_NH3j)  &
          +k(284)*n(idx_H3COj)  &
          +k(287)*n(idx_H3Oj)  &
          +k(300)*n(idx_HOCj)  &
          +2.d0*k(266)*n(idx_CH4j)  &
          +k(259)*n(idx_CHj)  &
          +k(312)*n(idx_O2Hj)  &
          +2.d0*k(273)*n(idx_H2COj)  &
          +2.d0*k(292)*n(idx_HCNHj)  &
          +k(299)*n(idx_HNOj)  &
          +2.d0*k(286)*n(idx_H3COj)  &
          +k(294)*n(idx_HCNHj)  &
          +3.d0*k(281)*n(idx_H3j)  &
          +k(285)*n(idx_H3COj)  &
          +2.d0*k(9)*n(idx_H2)  &
          +2.d0*k(306)*n(idx_NH2j)  &
          +k(301)*n(idx_HEHj)  &
          +k(263)*n(idx_CH3j)  &
          +k(297)*n(idx_HCO2j)  &
          +k(280)*n(idx_H3j)  &
          +2.d0*k(290)*n(idx_H3Oj)  &
          +k(288)*n(idx_H3Oj)  &
          +k(305)*n(idx_NHj)  &
          +k(274)*n(idx_H2COj)  &
          +k(296)*n(idx_HCO2j)  &
          +k(307)*n(idx_NH2j)  &
          +k(262)*n(idx_CH2j)  &
          +2.d0*k(309)*n(idx_NH3j)  &
          +k(303)*n(idx_N2Hj)  &
          +k(293)*n(idx_HCNHj)  &
          +k(313)*n(idx_OHj)  &
          +2.d0*k(270)*n(idx_H2j)  &
          +k(295)*n(idx_HCOj)  &
          +k(291)*n(idx_HCNj)  &
          +2.d0*k(278)*n(idx_H2Oj)
      pdj(9) =  &
          +k(287)*n(idx_H3Oj)  &
          +k(283)*n(idx_H3COj)
      pdj(10) =  &
          +k(279)*n(idx_H2Oj)  &
          +k(289)*n(idx_H3Oj)  &
          +k(290)*n(idx_H3Oj)  &
          +k(298)*n(idx_HCO2j)  &
          +k(282)*n(idx_H3COj)
      pdj(11) =  &
          +k(312)*n(idx_O2Hj)
      pdj(12) =  &
          +k(271)*n(idx_H2COj)  &
          +k(263)*n(idx_CH3j)  &
          +k(266)*n(idx_CH4j)  &
          +k(282)*n(idx_H3COj)
      pdj(13) =  &
          +k(1059)*n(idx_H2COj)  &
          +k(285)*n(idx_H3COj)
      pdj(14) =  &
          +k(286)*n(idx_H3COj)  &
          +k(274)*n(idx_H2COj)
      pdj(15) =  &
          +k(1061)*n(idx_MGj)
      pdj(17) =  &
          +k(276)*n(idx_H2NOj)  &
          +k(299)*n(idx_HNOj)
      pdj(18) =  &
          +k(291)*n(idx_HCNj)  &
          +k(292)*n(idx_HCNHj)
      pdj(19) =  &
          +k(298)*n(idx_HCO2j)  &
          +k(273)*n(idx_H2COj)  &
          +k(297)*n(idx_HCO2j)  &
          +k(272)*n(idx_H2COj)  &
          +k(295)*n(idx_HCOj)  &
          +k(284)*n(idx_H3COj)  &
          +k(300)*n(idx_HOCj)
      pdj(20) =  &
          +k(303)*n(idx_N2Hj)
      pdj(21) =  &
          +k(308)*n(idx_NH3j)
      pdj(22) =  &
          +k(1057)*n(idx_CH3j)  &
          +k(267)*n(idx_CH4j)
      pdj(24) =  &
          +k(1062)*n(idx_Nj)  &
          +k(305)*n(idx_NHj)  &
          +2.d0*k(302)*n(idx_N2j)  &
          +k(306)*n(idx_NH2j)  &
          +k(304)*n(idx_N2Hj)  &
          +k(268)*n(idx_CNj)  &
          +k(310)*n(idx_NOj)
      pdj(25) =  &
          +k(309)*n(idx_NH3j)  &
          +k(304)*n(idx_N2Hj)  &
          +k(307)*n(idx_NH2j)
      pdj(26) =  &
          +k(1060)*n(idx_HEj)  &
          +k(301)*n(idx_HEHj)
      pdj(27) =  &
          +k(275)*n(idx_H2NOj)
      pdj(29) =  &
          +k(296)*n(idx_HCO2j)
      pdj(53) =  &
          +k(1126)
      pdj(54) =  &
          -k(295)*n(idx_HCOj)
      pdj(55) =  &
          -k(1058)*n(idx_Hj)
      pdj(56) =  &
          -k(300)*n(idx_HOCj)
      pdj(57) =  &
          -k(1056)*n(idx_Cj)
      pdj(58) =  &
          -k(261)*n(idx_CH2j)  &
          -k(260)*n(idx_CH2j)  &
          -k(262)*n(idx_CH2j)
      pdj(59) =  &
          -k(259)*n(idx_CHj)
      pdj(60) =  &
          -k(273)*n(idx_H2COj)  &
          -k(272)*n(idx_H2COj)  &
          -k(274)*n(idx_H2COj)  &
          -k(1059)*n(idx_H2COj)  &
          -k(271)*n(idx_H2COj)
      pdj(61) =  &
          -k(1061)*n(idx_MGj)
      pdj(62) =  &
          -k(308)*n(idx_NH3j)  &
          -k(309)*n(idx_NH3j)
      pdj(63) =  &
          -k(310)*n(idx_NOj)
      pdj(64) =  &
          -k(268)*n(idx_CNj)
      pdj(65) =  &
          -k(269)*n(idx_COj)
      pdj(66) =  &
          -k(302)*n(idx_N2j)
      pdj(67) =  &
          -k(311)*n(idx_O2j)
      pdj(68) =  &
          -k(279)*n(idx_H2Oj)  &
          -k(277)*n(idx_H2Oj)  &
          -k(278)*n(idx_H2Oj)
      pdj(69) =  &
          -k(307)*n(idx_NH2j)  &
          -k(306)*n(idx_NH2j)
      pdj(70) =  &
          -k(1063)*n(idx_Oj)
      pdj(71) =  &
          -k(313)*n(idx_OHj)
      pdj(72) =  &
          -k(263)*n(idx_CH3j)  &
          -k(265)*n(idx_CH3j)  &
          -k(264)*n(idx_CH3j)  &
          -k(1057)*n(idx_CH3j)
      pdj(73) =  &
          -k(266)*n(idx_CH4j)  &
          -k(267)*n(idx_CH4j)
      pdj(74) =  &
          -k(1062)*n(idx_Nj)
      pdj(75) =  &
          -k(291)*n(idx_HCNj)
      pdj(76) =  &
          -k(305)*n(idx_NHj)
      pdj(77) =  &
          -k(270)*n(idx_H2j)
      pdj(78) =  &
          -k(1060)*n(idx_HEj)
      pdj(79) =  &
          -k(299)*n(idx_HNOj)
      pdj(80) =  &
          -k(276)*n(idx_H2NOj)  &
          -k(275)*n(idx_H2NOj)
      pdj(81) =  &
          -k(281)*n(idx_H3j)  &
          -k(280)*n(idx_H3j)
      pdj(82) =  &
          -k(286)*n(idx_H3COj)  &
          -k(283)*n(idx_H3COj)  &
          -k(282)*n(idx_H3COj)  &
          -k(285)*n(idx_H3COj)  &
          -k(284)*n(idx_H3COj)
      pdj(83) =  &
          -k(287)*n(idx_H3Oj)  &
          -k(288)*n(idx_H3Oj)  &
          -k(289)*n(idx_H3Oj)  &
          -k(290)*n(idx_H3Oj)
      pdj(84) =  &
          -k(292)*n(idx_HCNHj)  &
          -k(294)*n(idx_HCNHj)  &
          -k(293)*n(idx_HCNHj)
      pdj(85) =  &
          -k(298)*n(idx_HCO2j)  &
          -k(297)*n(idx_HCO2j)  &
          -k(296)*n(idx_HCO2j)
      pdj(86) =  &
          -k(301)*n(idx_HEHj)
      pdj(87) =  &
          -k(304)*n(idx_N2Hj)  &
          -k(303)*n(idx_N2Hj)
      pdj(88) =  &
          -k(312)*n(idx_O2Hj)
    elseif(j==2) then
      pdj(1) =  &
          +k(995)  &
          +k(1)*n(idx_O)
      pdj(2) =  &
          -k(410)*n(idx_Nj)  &
          -k(573)*n(idx_HEj)  &
          -k(809)*n(idx_N2)  &
          -k(406)*n(idx_HCNHj)  &
          -k(813)*n(idx_NO)  &
          -k(811)*n(idx_N)  &
          -k(414)*n(idx_Oj)  &
          -k(806)*n(idx_H2CO)  &
          -k(822)*n(idx_O)  &
          -k(53)*n(idx_NH2j)  &
          -k(823)*n(idx_OH)  &
          -k(52)*n(idx_N2j)  &
          -k(56)*n(idx_OHj)  &
          -k(407)*n(idx_HCNHj)  &
          -k(818)*n(idx_O2)  &
          -k(50)*n(idx_H2Oj)  &
          -k(400)*n(idx_COj)  &
          -k(839)*n(idx_H2)  &
          -k(416)*n(idx_O2Hj)  &
          -k(55)*n(idx_O2j)  &
          -k(413)*n(idx_NH2j)  &
          -k(51)*n(idx_Nj)  &
          -k(444)*n(idx_H2j)  &
          -k(805)*n(idx_CO2)  &
          -k(402)*n(idx_H2Oj)  &
          -k(409)*n(idx_HNOj)  &
          -k(821)*n(idx_O)  &
          -k(54)*n(idx_Oj)  &
          -k(995)  &
          -k(994)  &
          -k(1073)  &
          -k(812)*n(idx_NO)  &
          -k(851)*n(idx_H)  &
          -k(1049)*n(idx_H2)  &
          -k(807)*n(idx_HCO)  &
          -k(124)*n(idx_HEj)  &
          -k(808)*n(idx_HNO)  &
          -k(405)*n(idx_HCNj)  &
          -k(814)*n(idx_NO)  &
          -k(3)*n(idx_H2)  &
          -k(16)*n(idx_Cj)  &
          -k(47)*n(idx_CNj)  &
          -k(408)*n(idx_HCOj)  &
          -k(819)*n(idx_O2H)  &
          -k(506)*n(idx_H3j)  &
          -k(401)*n(idx_H2COj)  &
          -k(48)*n(idx_COj)  &
          -k(225)  &
          -k(49)*n(idx_H2COj)  &
          -k(815)*n(idx_O2)  &
          -k(810)*n(idx_N)  &
          -k(1)*n(idx_O)  &
          -k(820)*n(idx_O2H)  &
          -k(87)*n(idx_H2j)  &
          -k(417)*n(idx_OHj)  &
          -k(10)*n(idx_H)  &
          -k(415)*n(idx_O2j)  &
          -k(412)*n(idx_NHj)  &
          -k(403)*n(idx_H3COj)  &
          -k(404)*n(idx_H3Oj)  &
          -k(816)*n(idx_O2)  &
          -k(411)*n(idx_N2Hj)  &
          -k(817)*n(idx_O2)  &
          -k(72)*n(idx_Hj)
      pdj(3) =  &
          -k(1)*n(idx_O)  &
          +k(417)*n(idx_OHj)  &
          +k(54)*n(idx_Oj)  &
          -k(822)*n(idx_O)  &
          -k(821)*n(idx_O)  &
          +k(818)*n(idx_O2)  &
          +k(415)*n(idx_O2j)  &
          +k(812)*n(idx_NO)  &
          +k(816)*n(idx_O2)
      pdj(4) =  &
          +k(407)*n(idx_HCNHj)
      pdj(5) =  &
          +k(809)*n(idx_N2)  &
          +k(812)*n(idx_NO)  &
          +k(406)*n(idx_HCNHj)
      pdj(6) =  &
          +k(87)*n(idx_H2j)  &
          -k(3)*n(idx_H2)  &
          +k(3)*n(idx_H2)  &
          +k(506)*n(idx_H3j)  &
          +k(851)*n(idx_H)  &
          -k(839)*n(idx_H2)  &
          -k(1049)*n(idx_H2)
      pdj(7) =  &
          +k(10)*n(idx_H)  &
          +k(994)  &
          +k(400)*n(idx_COj)  &
          +k(16)*n(idx_Cj)  &
          +k(3)*n(idx_H2)  &
          +k(851)*n(idx_H)  &
          +k(225)  &
          +k(811)*n(idx_N)  &
          +k(822)*n(idx_O)
      pdj(8) =  &
          +k(573)*n(idx_HEj)  &
          +2.d0*k(10)*n(idx_H)  &
          +k(72)*n(idx_Hj)  &
          +k(821)*n(idx_O)  &
          +k(839)*n(idx_H2)  &
          +k(994)  &
          +k(815)*n(idx_O2)  &
          +k(444)*n(idx_H2j)  &
          +k(3)*n(idx_H2)  &
          -k(10)*n(idx_H)  &
          +k(816)*n(idx_O2)  &
          +k(410)*n(idx_Nj)  &
          +k(225)  &
          +k(814)*n(idx_NO)  &
          +k(823)*n(idx_OH)  &
          -k(851)*n(idx_H)  &
          +k(810)*n(idx_N)  &
          +k(414)*n(idx_Oj)
      pdj(9) =  &
          +k(50)*n(idx_H2Oj)  &
          +k(404)*n(idx_H3Oj)
      pdj(10) =  &
          +k(402)*n(idx_H2Oj)  &
          +k(819)*n(idx_O2H)  &
          -k(823)*n(idx_OH)  &
          +k(817)*n(idx_O2)  &
          +k(822)*n(idx_O)  &
          +k(56)*n(idx_OHj)
      pdj(11) =  &
          +k(416)*n(idx_O2Hj)  &
          +k(820)*n(idx_O2H)  &
          -k(815)*n(idx_O2)  &
          +k(55)*n(idx_O2j)  &
          -k(816)*n(idx_O2)  &
          -k(818)*n(idx_O2)  &
          -k(817)*n(idx_O2)
      pdj(12) =  &
          +k(808)*n(idx_HNO)  &
          +k(807)*n(idx_HCO)  &
          +k(820)*n(idx_O2H)  &
          +k(839)*n(idx_H2)  &
          +k(806)*n(idx_H2CO)
      pdj(13) =  &
          +k(403)*n(idx_H3COj)  &
          +k(49)*n(idx_H2COj)  &
          -k(806)*n(idx_H2CO)
      pdj(14) =  &
          -k(807)*n(idx_HCO)  &
          +k(823)*n(idx_OH)  &
          +k(819)*n(idx_O2H)  &
          +k(401)*n(idx_H2COj)  &
          +k(806)*n(idx_H2CO)  &
          +k(813)*n(idx_NO)  &
          +k(818)*n(idx_O2)  &
          +k(805)*n(idx_CO2)
      pdj(17) =  &
          -k(812)*n(idx_NO)  &
          +k(409)*n(idx_HNOj)  &
          -k(814)*n(idx_NO)  &
          -k(813)*n(idx_NO)  &
          +k(808)*n(idx_HNO)
      pdj(18) =  &
          +k(405)*n(idx_HCNj)  &
          +k(810)*n(idx_N)  &
          +k(47)*n(idx_CNj)
      pdj(19) =  &
          +k(821)*n(idx_O)  &
          +k(408)*n(idx_HCOj)  &
          +k(48)*n(idx_COj)  &
          +k(807)*n(idx_HCO)  &
          +k(816)*n(idx_O2)  &
          +k(817)*n(idx_O2)  &
          +k(805)*n(idx_CO2)
      pdj(20) =  &
          +k(411)*n(idx_N2Hj)  &
          +k(52)*n(idx_N2j)  &
          -k(809)*n(idx_N2)
      pdj(21) =  &
          +k(53)*n(idx_NH2j)
      pdj(22) =  &
          +k(1049)*n(idx_H2)
      pdj(24) =  &
          -k(811)*n(idx_N)  &
          -k(810)*n(idx_N)  &
          +k(51)*n(idx_Nj)  &
          +k(412)*n(idx_NHj)  &
          +k(813)*n(idx_NO)  &
          +k(809)*n(idx_N2)
      pdj(25) =  &
          +k(413)*n(idx_NH2j)  &
          +k(811)*n(idx_N)
      pdj(26) =  &
          +k(124)*n(idx_HEj)  &
          +k(573)*n(idx_HEj)
      pdj(27) =  &
          -k(808)*n(idx_HNO)
      pdj(29) =  &
          +k(815)*n(idx_O2)  &
          -k(805)*n(idx_CO2)
      pdj(33) =  &
          -k(820)*n(idx_O2H)  &
          -k(819)*n(idx_O2H)
      pdj(34) =  &
          +k(814)*n(idx_NO)
      pdj(38) =  &
          +k(1073)
      pdj(54) =  &
          +k(1)*n(idx_O)  &
          -k(408)*n(idx_HCOj)  &
          +k(415)*n(idx_O2j)  &
          +k(400)*n(idx_COj)
      pdj(55) =  &
          -k(72)*n(idx_Hj)
      pdj(57) =  &
          +k(573)*n(idx_HEj)  &
          -k(16)*n(idx_Cj)
      pdj(58) =  &
          +k(407)*n(idx_HCNHj)  &
          +k(409)*n(idx_HNOj)  &
          +k(406)*n(idx_HCNHj)  &
          +k(405)*n(idx_HCNj)  &
          +k(402)*n(idx_H2Oj)  &
          +k(417)*n(idx_OHj)  &
          +k(413)*n(idx_NH2j)  &
          +k(404)*n(idx_H3Oj)  &
          +k(401)*n(idx_H2COj)  &
          +k(506)*n(idx_H3j)  &
          +k(408)*n(idx_HCOj)  &
          +k(411)*n(idx_N2Hj)  &
          +k(416)*n(idx_O2Hj)  &
          +k(444)*n(idx_H2j)  &
          +k(412)*n(idx_NHj)  &
          +k(403)*n(idx_H3COj)
      pdj(59) =  &
          +k(50)*n(idx_H2Oj)  &
          +k(72)*n(idx_Hj)  &
          +k(55)*n(idx_O2j)  &
          +k(47)*n(idx_CNj)  &
          +k(995)  &
          +k(52)*n(idx_N2j)  &
          +k(51)*n(idx_Nj)  &
          +k(124)*n(idx_HEj)  &
          +k(53)*n(idx_NH2j)  &
          +k(48)*n(idx_COj)  &
          +k(87)*n(idx_H2j)  &
          +k(16)*n(idx_Cj)  &
          +k(49)*n(idx_H2COj)  &
          +k(54)*n(idx_Oj)  &
          +k(56)*n(idx_OHj)
      pdj(60) =  &
          -k(49)*n(idx_H2COj)  &
          -k(401)*n(idx_H2COj)
      pdj(64) =  &
          -k(47)*n(idx_CNj)  &
          +k(410)*n(idx_Nj)
      pdj(65) =  &
          +k(414)*n(idx_Oj)  &
          -k(400)*n(idx_COj)  &
          -k(48)*n(idx_COj)
      pdj(66) =  &
          -k(52)*n(idx_N2j)
      pdj(67) =  &
          -k(55)*n(idx_O2j)  &
          -k(415)*n(idx_O2j)
      pdj(68) =  &
          -k(50)*n(idx_H2Oj)  &
          -k(402)*n(idx_H2Oj)
      pdj(69) =  &
          -k(53)*n(idx_NH2j)  &
          -k(413)*n(idx_NH2j)
      pdj(70) =  &
          -k(54)*n(idx_Oj)  &
          -k(414)*n(idx_Oj)
      pdj(71) =  &
          -k(56)*n(idx_OHj)  &
          -k(417)*n(idx_OHj)
      pdj(74) =  &
          -k(410)*n(idx_Nj)  &
          -k(51)*n(idx_Nj)
      pdj(75) =  &
          -k(405)*n(idx_HCNj)
      pdj(76) =  &
          -k(412)*n(idx_NHj)
      pdj(77) =  &
          -k(444)*n(idx_H2j)  &
          -k(87)*n(idx_H2j)
      pdj(78) =  &
          -k(573)*n(idx_HEj)  &
          -k(124)*n(idx_HEj)
      pdj(79) =  &
          -k(409)*n(idx_HNOj)
      pdj(81) =  &
          -k(506)*n(idx_H3j)
      pdj(82) =  &
          -k(403)*n(idx_H3COj)
      pdj(83) =  &
          -k(404)*n(idx_H3Oj)
      pdj(84) =  &
          -k(406)*n(idx_HCNHj)  &
          -k(407)*n(idx_HCNHj)
      pdj(87) =  &
          -k(411)*n(idx_N2Hj)
      pdj(88) =  &
          -k(416)*n(idx_O2Hj)
    elseif(j==3) then
      pdj(1) =  &
          +k(214)  &
          +k(1)*n(idx_CH)  &
          +k(256)
      pdj(2) =  &
          -k(821)*n(idx_CH)  &
          -k(822)*n(idx_CH)  &
          +k(779)*n(idx_CH2)  &
          -k(1)*n(idx_CH)
      pdj(3) =  &
          -k(953)*n(idx_NH2)  &
          -k(957)*n(idx_O2H)  &
          -k(927)*n(idx_NH)  &
          -k(797)*n(idx_CH3)  &
          -k(942)*n(idx_H2O)  &
          -k(214)  &
          -k(725)*n(idx_OHj)  &
          -k(719)*n(idx_HCO2j)  &
          -k(940)*n(idx_H2CN)  &
          -k(722)*n(idx_NH2j)  &
          -k(777)*n(idx_CH2)  &
          -k(821)*n(idx_CH)  &
          -k(959)*n(idx_OCN)  &
          -k(946)*n(idx_HCO)  &
          -k(1052)*n(idx_H)  &
          -k(458)*n(idx_H2j)  &
          -k(718)*n(idx_H2Oj)  &
          -k(926)*n(idx_NH)  &
          -k(1111)  &
          -k(386)*n(idx_CH3j)  &
          -k(717)*n(idx_CH4j)  &
          -k(365)*n(idx_CH2j)  &
          -k(720)*n(idx_N2j)  &
          -k(936)*n(idx_CH4)  &
          -k(951)*n(idx_N2)  &
          -k(778)*n(idx_CH2)  &
          -k(952)*n(idx_NH2)  &
          -k(194)*n(idx_CNj)  &
          -k(960)*n(idx_OH)  &
          -k(1044)*n(idx_C)  &
          -k(938)*n(idx_CN)  &
          -k(941)*n(idx_H2CO)  &
          -k(776)*n(idx_CH2)  &
          -k(387)*n(idx_CH3j)  &
          -k(358)*n(idx_CHj)  &
          -k(958)*n(idx_OCN)  &
          -k(662)*n(idx_NHj)  &
          -k(83)*n(idx_Hj)  &
          -k(937)*n(idx_CN)  &
          -k(524)*n(idx_H3j)  &
          -k(949)*n(idx_HNO)  &
          -k(846)*n(idx_H2)  &
          -k(1041)*n(idx_Cj)  &
          -k(943)*n(idx_HCN)  &
          -k(947)*n(idx_HCO)  &
          -k(822)*n(idx_CH)  &
          -k(1)*n(idx_CH)  &
          -k(944)*n(idx_HCN)  &
          -k(195)*n(idx_COj)  &
          -k(955)*n(idx_NO2)  &
          -k(939)*n(idx_CO2)  &
          -k(721)*n(idx_N2Hj)  &
          -k(798)*n(idx_CH3)  &
          -k(196)*n(idx_N2j)  &
          -4.d0*k(1055)*n(idx_O)  &
          -k(779)*n(idx_CH2)  &
          -k(945)*n(idx_HCN)  &
          -k(948)*n(idx_HNO)  &
          -k(256)  &
          -k(950)*n(idx_HNO)  &
          -k(956)*n(idx_NO)  &
          -k(724)*n(idx_O2Hj)  &
          -k(954)*n(idx_NH3)  &
          -k(723)*n(idx_NH3j)  &
          -k(525)*n(idx_H3j)
      pdj(5) =  &
          -k(944)*n(idx_HCN)  &
          -k(945)*n(idx_HCN)  &
          -k(943)*n(idx_HCN)
      pdj(6) =  &
          +k(723)*n(idx_NH3j)  &
          +k(718)*n(idx_H2Oj)  &
          -k(846)*n(idx_H2)  &
          +k(387)*n(idx_CH3j)  &
          +k(940)*n(idx_H2CN)  &
          +k(525)*n(idx_H3j)  &
          +k(797)*n(idx_CH3)  &
          +k(776)*n(idx_CH2)
      pdj(7) =  &
          +k(822)*n(idx_CH)  &
          +k(938)*n(idx_CN)  &
          -k(1044)*n(idx_C)
      pdj(8) =  &
          +k(458)*n(idx_H2j)  &
          +2.d0*k(777)*n(idx_CH2)  &
          +k(948)*n(idx_HNO)  &
          +k(778)*n(idx_CH2)  &
          -k(1052)*n(idx_H)  &
          +k(365)*n(idx_CH2j)  &
          +k(846)*n(idx_H2)  &
          +k(725)*n(idx_OHj)  &
          +k(386)*n(idx_CH3j)  &
          +k(821)*n(idx_CH)  &
          +k(960)*n(idx_OH)  &
          +k(83)*n(idx_Hj)  &
          +k(798)*n(idx_CH3)  &
          +k(945)*n(idx_HCN)  &
          +k(524)*n(idx_H3j)  &
          +k(358)*n(idx_CHj)  &
          +k(926)*n(idx_NH)  &
          +k(946)*n(idx_HCO)  &
          +k(722)*n(idx_NH2j)  &
          +k(952)*n(idx_NH2)  &
          +k(797)*n(idx_CH3)
      pdj(9) =  &
          -k(942)*n(idx_H2O)
      pdj(10) =  &
          +k(822)*n(idx_CH)  &
          +k(941)*n(idx_H2CO)  &
          +k(954)*n(idx_NH3)  &
          +k(957)*n(idx_O2H)  &
          +k(846)*n(idx_H2)  &
          +2.d0*k(942)*n(idx_H2O)  &
          +k(953)*n(idx_NH2)  &
          +k(949)*n(idx_HNO)  &
          +k(936)*n(idx_CH4)  &
          +k(943)*n(idx_HCN)  &
          +k(947)*n(idx_HCO)  &
          +k(717)*n(idx_CH4j)  &
          +k(1052)*n(idx_H)  &
          +k(927)*n(idx_NH)  &
          +k(779)*n(idx_CH2)  &
          -k(960)*n(idx_OH)
      pdj(11) =  &
          +k(950)*n(idx_HNO)  &
          +k(956)*n(idx_NO)  &
          +k(959)*n(idx_OCN)  &
          +k(955)*n(idx_NO2)  &
          +k(939)*n(idx_CO2)  &
          +k(960)*n(idx_OH)  &
          +k(719)*n(idx_HCO2j)  &
          +2.d0*k(1055)*n(idx_O)  &
          +k(957)*n(idx_O2H)  &
          +k(724)*n(idx_O2Hj)
      pdj(12) =  &
          -k(779)*n(idx_CH2)  &
          -k(778)*n(idx_CH2)  &
          -k(777)*n(idx_CH2)  &
          -k(776)*n(idx_CH2)
      pdj(13) =  &
          -k(941)*n(idx_H2CO)  &
          +k(798)*n(idx_CH3)
      pdj(14) =  &
          -k(947)*n(idx_HCO)  &
          +k(941)*n(idx_H2CO)  &
          -k(946)*n(idx_HCO)  &
          +k(778)*n(idx_CH2)
      pdj(16) =  &
          -k(954)*n(idx_NH3)
      pdj(17) =  &
          +k(958)*n(idx_OCN)  &
          +k(955)*n(idx_NO2)  &
          +k(951)*n(idx_N2)  &
          +k(949)*n(idx_HNO)  &
          +k(938)*n(idx_CN)  &
          -k(956)*n(idx_NO)  &
          +k(926)*n(idx_NH)
      pdj(18) =  &
          -k(938)*n(idx_CN)  &
          +k(943)*n(idx_HCN)  &
          +k(194)*n(idx_CNj)  &
          +k(959)*n(idx_OCN)  &
          -k(937)*n(idx_CN)
      pdj(19) =  &
          +k(777)*n(idx_CH2)  &
          +k(944)*n(idx_HCN)  &
          +k(958)*n(idx_OCN)  &
          +k(939)*n(idx_CO2)  &
          +k(821)*n(idx_CH)  &
          +k(947)*n(idx_HCO)  &
          +k(797)*n(idx_CH3)  &
          +k(195)*n(idx_COj)  &
          +k(937)*n(idx_CN)  &
          +k(776)*n(idx_CH2)  &
          +k(1044)*n(idx_C)
      pdj(20) =  &
          -k(951)*n(idx_N2)  &
          +k(196)*n(idx_N2j)  &
          +k(721)*n(idx_N2Hj)
      pdj(21) =  &
          -k(953)*n(idx_NH2)  &
          +k(954)*n(idx_NH3)  &
          -k(952)*n(idx_NH2)
      pdj(22) =  &
          -k(798)*n(idx_CH3)  &
          +k(936)*n(idx_CH4)  &
          -k(797)*n(idx_CH3)
      pdj(23) =  &
          -k(936)*n(idx_CH4)
      pdj(24) =  &
          +k(956)*n(idx_NO)  &
          +k(720)*n(idx_N2j)  &
          +k(951)*n(idx_N2)  &
          +k(937)*n(idx_CN)  &
          +k(662)*n(idx_NHj)  &
          +k(927)*n(idx_NH)
      pdj(25) =  &
          -k(926)*n(idx_NH)  &
          +k(950)*n(idx_HNO)  &
          -k(927)*n(idx_NH)  &
          +k(944)*n(idx_HCN)  &
          +k(953)*n(idx_NH2)
      pdj(27) =  &
          -k(948)*n(idx_HNO)  &
          -k(950)*n(idx_HNO)  &
          +k(952)*n(idx_NH2)  &
          -k(949)*n(idx_HNO)
      pdj(29) =  &
          -k(939)*n(idx_CO2)  &
          +k(946)*n(idx_HCO)
      pdj(30) =  &
          -k(940)*n(idx_H2CN)
      pdj(32) =  &
          +k(948)*n(idx_HNO)  &
          -k(955)*n(idx_NO2)
      pdj(33) =  &
          -k(957)*n(idx_O2H)
      pdj(34) =  &
          -k(958)*n(idx_OCN)  &
          -k(959)*n(idx_OCN)  &
          +k(940)*n(idx_H2CN)  &
          +k(945)*n(idx_HCN)
      pdj(40) =  &
          +k(1111)
      pdj(54) =  &
          +k(365)*n(idx_CH2j)  &
          +k(387)*n(idx_CH3j)  &
          +k(1)*n(idx_CH)  &
          +k(719)*n(idx_HCO2j)
      pdj(55) =  &
          -k(83)*n(idx_Hj)
      pdj(57) =  &
          -k(1041)*n(idx_Cj)
      pdj(58) =  &
          -k(365)*n(idx_CH2j)
      pdj(59) =  &
          -k(358)*n(idx_CHj)
      pdj(60) =  &
          +k(386)*n(idx_CH3j)
      pdj(62) =  &
          -k(723)*n(idx_NH3j)
      pdj(63) =  &
          +k(720)*n(idx_N2j)
      pdj(64) =  &
          -k(194)*n(idx_CNj)
      pdj(65) =  &
          +k(358)*n(idx_CHj)  &
          +k(1041)*n(idx_Cj)  &
          -k(195)*n(idx_COj)
      pdj(66) =  &
          -k(196)*n(idx_N2j)  &
          -k(720)*n(idx_N2j)
      pdj(67) =  &
          +k(725)*n(idx_OHj)  &
          +k(718)*n(idx_H2Oj)
      pdj(68) =  &
          +k(524)*n(idx_H3j)  &
          -k(718)*n(idx_H2Oj)
      pdj(69) =  &
          -k(722)*n(idx_NH2j)
      pdj(70) =  &
          +k(256)  &
          +k(214)  &
          +k(194)*n(idx_CNj)  &
          +k(196)*n(idx_N2j)  &
          +k(83)*n(idx_Hj)  &
          +k(195)*n(idx_COj)
      pdj(71) =  &
          +k(458)*n(idx_H2j)  &
          +k(721)*n(idx_N2Hj)  &
          +k(525)*n(idx_H3j)  &
          +k(662)*n(idx_NHj)  &
          +k(724)*n(idx_O2Hj)  &
          -k(725)*n(idx_OHj)
      pdj(72) =  &
          -k(386)*n(idx_CH3j)  &
          -k(387)*n(idx_CH3j)  &
          +k(717)*n(idx_CH4j)
      pdj(73) =  &
          -k(717)*n(idx_CH4j)
      pdj(76) =  &
          -k(662)*n(idx_NHj)
      pdj(77) =  &
          -k(458)*n(idx_H2j)
      pdj(79) =  &
          +k(723)*n(idx_NH3j)  &
          +k(722)*n(idx_NH2j)
      pdj(81) =  &
          -k(524)*n(idx_H3j)  &
          -k(525)*n(idx_H3j)
      pdj(85) =  &
          -k(719)*n(idx_HCO2j)
      pdj(87) =  &
          -k(721)*n(idx_N2Hj)
      pdj(88) =  &
          -k(724)*n(idx_O2Hj)
    elseif(j==4) then
      pdj(3) =  &
          +k(734)*n(idx_OHj)
      pdj(4) =  &
          -k(529)*n(idx_H3Oj)  &
          -k(489)*n(idx_H2Oj)  &
          -k(655)*n(idx_NHj)  &
          -k(351)*n(idx_CHj)  &
          -k(595)*n(idx_HEj)  &
          -k(594)*n(idx_HEj)  &
          -k(734)*n(idx_OHj)  &
          -k(236)  &
          -k(562)*n(idx_O2Hj)  &
          -k(593)*n(idx_HEj)  &
          -k(561)*n(idx_N2Hj)  &
          -k(670)*n(idx_NH2j)  &
          -k(560)*n(idx_HNOj)  &
          -k(557)*n(idx_H2COj)  &
          -k(1015)  &
          -k(515)*n(idx_H3j)  &
          -k(2)*n(idx_Hj)  &
          -k(559)*n(idx_HCOj)  &
          -k(541)*n(idx_HCNj)  &
          -k(860)*n(idx_H)  &
          -k(1125)  &
          -k(558)*n(idx_H3COj)
      pdj(5) =  &
          +k(2)*n(idx_Hj)  &
          +k(860)*n(idx_H)
      pdj(6) =  &
          +k(515)*n(idx_H3j)
      pdj(7) =  &
          +k(351)*n(idx_CHj)  &
          +k(595)*n(idx_HEj)
      pdj(8) =  &
          +k(1015)  &
          +k(594)*n(idx_HEj)  &
          +k(593)*n(idx_HEj)  &
          -k(860)*n(idx_H)  &
          +k(236)  &
          +k(860)*n(idx_H)
      pdj(9) =  &
          +k(529)*n(idx_H3Oj)
      pdj(10) =  &
          +k(489)*n(idx_H2Oj)
      pdj(11) =  &
          +k(562)*n(idx_O2Hj)
      pdj(13) =  &
          +k(558)*n(idx_H3COj)
      pdj(14) =  &
          +k(557)*n(idx_H2COj)
      pdj(17) =  &
          +k(560)*n(idx_HNOj)
      pdj(18) =  &
          +k(541)*n(idx_HCNj)  &
          +k(1015)  &
          +k(236)
      pdj(19) =  &
          +k(559)*n(idx_HCOj)
      pdj(20) =  &
          +k(561)*n(idx_N2Hj)
      pdj(24) =  &
          +k(655)*n(idx_NHj)  &
          +k(594)*n(idx_HEj)
      pdj(25) =  &
          +k(670)*n(idx_NH2j)
      pdj(26) =  &
          +k(593)*n(idx_HEj)  &
          +k(595)*n(idx_HEj)  &
          +k(594)*n(idx_HEj)
      pdj(52) =  &
          +k(1125)
      pdj(54) =  &
          -k(559)*n(idx_HCOj)
      pdj(55) =  &
          +k(2)*n(idx_Hj)  &
          -k(2)*n(idx_Hj)
      pdj(57) =  &
          +k(594)*n(idx_HEj)
      pdj(59) =  &
          -k(351)*n(idx_CHj)
      pdj(60) =  &
          -k(557)*n(idx_H2COj)
      pdj(64) =  &
          +k(593)*n(idx_HEj)
      pdj(68) =  &
          -k(489)*n(idx_H2Oj)
      pdj(69) =  &
          -k(670)*n(idx_NH2j)
      pdj(71) =  &
          -k(734)*n(idx_OHj)
      pdj(75) =  &
          -k(541)*n(idx_HCNj)
      pdj(76) =  &
          -k(655)*n(idx_NHj)  &
          +k(595)*n(idx_HEj)
      pdj(78) =  &
          -k(593)*n(idx_HEj)  &
          -k(595)*n(idx_HEj)  &
          -k(594)*n(idx_HEj)
      pdj(79) =  &
          -k(560)*n(idx_HNOj)
      pdj(81) =  &
          -k(515)*n(idx_H3j)
      pdj(82) =  &
          -k(558)*n(idx_H3COj)
      pdj(83) =  &
          -k(529)*n(idx_H3Oj)
      pdj(84) =  &
          +k(541)*n(idx_HCNj)  &
          +k(655)*n(idx_NHj)  &
          +k(559)*n(idx_HCOj)  &
          +k(562)*n(idx_O2Hj)  &
          +k(561)*n(idx_N2Hj)  &
          +k(351)*n(idx_CHj)  &
          +k(734)*n(idx_OHj)  &
          +k(557)*n(idx_H2COj)  &
          +k(515)*n(idx_H3j)  &
          +k(489)*n(idx_H2Oj)  &
          +k(558)*n(idx_H3COj)  &
          +k(560)*n(idx_HNOj)  &
          +k(670)*n(idx_NH2j)  &
          +k(529)*n(idx_H3Oj)
      pdj(87) =  &
          -k(561)*n(idx_N2Hj)
      pdj(88) =  &
          -k(562)*n(idx_O2Hj)
    elseif(j==5) then
      pdj(2) =  &
          +k(710)*n(idx_Oj)  &
          +k(587)*n(idx_HEj)
      pdj(3) =  &
          -k(943)*n(idx_O)  &
          -k(945)*n(idx_O)  &
          -k(944)*n(idx_O)  &
          +k(731)*n(idx_OHj)
      pdj(5) =  &
          -k(943)*n(idx_O)  &
          -k(653)*n(idx_NHj)  &
          -k(709)*n(idx_Oj)  &
          -k(544)*n(idx_HCOj)  &
          -k(59)*n(idx_CNj)  &
          -k(118)*n(idx_COj)  &
          -k(944)*n(idx_O)  &
          -k(75)*n(idx_Hj)  &
          -k(1011)  &
          -k(349)*n(idx_CHj)  &
          -k(589)*n(idx_HEj)  &
          -k(141)*n(idx_Nj)  &
          -k(857)*n(idx_H)  &
          -k(92)*n(idx_H2j)  &
          -k(546)*n(idx_N2Hj)  &
          -k(545)*n(idx_HNOj)  &
          -k(710)*n(idx_Oj)  &
          -k(119)*n(idx_N2j)  &
          -k(588)*n(idx_HEj)  &
          -k(587)*n(idx_HEj)  &
          -k(513)*n(idx_H3j)  &
          -k(668)*n(idx_NH2j)  &
          -k(486)*n(idx_H2Oj)  &
          -k(945)*n(idx_O)  &
          -k(731)*n(idx_OHj)  &
          -k(1086)  &
          -k(965)*n(idx_OH)  &
          -k(547)*n(idx_O2Hj)  &
          -k(966)*n(idx_OH)  &
          -k(233)  &
          -k(542)*n(idx_H2COj)  &
          -k(543)*n(idx_H3COj)  &
          -k(586)*n(idx_HEj)  &
          -k(528)*n(idx_H3Oj)  &
          -k(538)*n(idx_HCNj)
      pdj(6) =  &
          +k(857)*n(idx_H)  &
          +k(92)*n(idx_H2j)  &
          +k(513)*n(idx_H3j)
      pdj(7) =  &
          +k(349)*n(idx_CHj)
      pdj(8) =  &
          +k(233)  &
          +k(945)*n(idx_O)  &
          -k(857)*n(idx_H)  &
          +k(1011)  &
          +k(588)*n(idx_HEj)  &
          +k(75)*n(idx_Hj)  &
          +k(586)*n(idx_HEj)
      pdj(9) =  &
          +k(965)*n(idx_OH)  &
          +k(528)*n(idx_H3Oj)
      pdj(10) =  &
          -k(965)*n(idx_OH)  &
          +k(943)*n(idx_O)  &
          +k(486)*n(idx_H2Oj)  &
          -k(966)*n(idx_OH)
      pdj(11) =  &
          +k(547)*n(idx_O2Hj)
      pdj(13) =  &
          +k(543)*n(idx_H3COj)
      pdj(14) =  &
          +k(542)*n(idx_H2COj)
      pdj(17) =  &
          +k(545)*n(idx_HNOj)
      pdj(18) =  &
          +k(233)  &
          +k(538)*n(idx_HCNj)  &
          +k(943)*n(idx_O)  &
          +k(1011)  &
          +k(965)*n(idx_OH)  &
          +k(59)*n(idx_CNj)  &
          +k(857)*n(idx_H)
      pdj(19) =  &
          +k(944)*n(idx_O)  &
          +k(118)*n(idx_COj)  &
          +k(544)*n(idx_HCOj)  &
          +k(966)*n(idx_OH)
      pdj(20) =  &
          +k(546)*n(idx_N2Hj)  &
          +k(119)*n(idx_N2j)
      pdj(21) =  &
          +k(966)*n(idx_OH)
      pdj(24) =  &
          +k(588)*n(idx_HEj)  &
          +k(589)*n(idx_HEj)  &
          +k(653)*n(idx_NHj)  &
          +k(709)*n(idx_Oj)  &
          +k(141)*n(idx_Nj)
      pdj(25) =  &
          +k(944)*n(idx_O)  &
          +k(668)*n(idx_NH2j)
      pdj(26) =  &
          +k(588)*n(idx_HEj)  &
          +k(589)*n(idx_HEj)  &
          +k(586)*n(idx_HEj)  &
          +k(587)*n(idx_HEj)
      pdj(34) =  &
          +k(945)*n(idx_O)
      pdj(44) =  &
          +k(1086)
      pdj(54) =  &
          -k(544)*n(idx_HCOj)  &
          +k(709)*n(idx_Oj)
      pdj(55) =  &
          -k(75)*n(idx_Hj)
      pdj(57) =  &
          +k(588)*n(idx_HEj)
      pdj(59) =  &
          -k(349)*n(idx_CHj)  &
          +k(589)*n(idx_HEj)
      pdj(60) =  &
          -k(542)*n(idx_H2COj)
      pdj(63) =  &
          +k(710)*n(idx_Oj)
      pdj(64) =  &
          -k(59)*n(idx_CNj)  &
          +k(586)*n(idx_HEj)
      pdj(65) =  &
          -k(118)*n(idx_COj)
      pdj(66) =  &
          -k(119)*n(idx_N2j)
      pdj(68) =  &
          -k(486)*n(idx_H2Oj)
      pdj(69) =  &
          -k(668)*n(idx_NH2j)
      pdj(70) =  &
          -k(710)*n(idx_Oj)  &
          -k(709)*n(idx_Oj)
      pdj(71) =  &
          -k(731)*n(idx_OHj)
      pdj(74) =  &
          -k(141)*n(idx_Nj)  &
          +k(587)*n(idx_HEj)
      pdj(75) =  &
          +k(92)*n(idx_H2j)  &
          +k(141)*n(idx_Nj)  &
          +k(118)*n(idx_COj)  &
          +k(59)*n(idx_CNj)  &
          +k(119)*n(idx_N2j)  &
          +k(75)*n(idx_Hj)  &
          -k(538)*n(idx_HCNj)
      pdj(76) =  &
          -k(653)*n(idx_NHj)
      pdj(77) =  &
          -k(92)*n(idx_H2j)
      pdj(78) =  &
          -k(586)*n(idx_HEj)  &
          -k(587)*n(idx_HEj)  &
          -k(588)*n(idx_HEj)  &
          -k(589)*n(idx_HEj)
      pdj(79) =  &
          -k(545)*n(idx_HNOj)
      pdj(81) =  &
          -k(513)*n(idx_H3j)
      pdj(82) =  &
          -k(543)*n(idx_H3COj)
      pdj(83) =  &
          -k(528)*n(idx_H3Oj)
      pdj(84) =  &
          +k(543)*n(idx_H3COj)  &
          +k(731)*n(idx_OHj)  &
          +k(513)*n(idx_H3j)  &
          +k(538)*n(idx_HCNj)  &
          +k(546)*n(idx_N2Hj)  &
          +k(542)*n(idx_H2COj)  &
          +k(349)*n(idx_CHj)  &
          +k(486)*n(idx_H2Oj)  &
          +k(545)*n(idx_HNOj)  &
          +k(653)*n(idx_NHj)  &
          +k(547)*n(idx_O2Hj)  &
          +k(528)*n(idx_H3Oj)  &
          +k(544)*n(idx_HCOj)  &
          +k(668)*n(idx_NH2j)
      pdj(87) =  &
          -k(546)*n(idx_N2Hj)
      pdj(88) =  &
          -k(547)*n(idx_O2Hj)
    elseif(j==6) then
      pdj(1) =  &
          +k(9)*n(idx_E)  &
          -k(9)*n(idx_E)  &
          +k(209)  &
          +k(208)
      pdj(2) =  &
          -k(3)*n(idx_CH)  &
          -k(839)*n(idx_CH)  &
          +k(836)*n(idx_C)  &
          -k(1049)*n(idx_CH)
      pdj(3) =  &
          +2.d0*k(7)*n(idx_O2)  &
          -k(846)*n(idx_O)  &
          +k(8)*n(idx_OH)
      pdj(5) =  &
          +k(840)*n(idx_CN)
      pdj(6) =  &
          -k(1049)*n(idx_CH)  &
          -k(9)*n(idx_E)  &
          -k(208)  &
          -k(8)*n(idx_OH)  &
          -k(476)*n(idx_O2Hj)  &
          -k(473)*n(idx_NHj)  &
          -k(469)*n(idx_HEHj)  &
          -k(463)*n(idx_CNj)  &
          +k(7)*n(idx_O2)  &
          -k(471)*n(idx_N2j)  &
          +k(5)*n(idx_H2O)  &
          -k(470)*n(idx_Nj)  &
          -k(475)*n(idx_Oj)  &
          -k(465)*n(idx_COj)  &
          -k(7)*n(idx_O2)  &
          -k(840)*n(idx_CN)  &
          -k(1047)*n(idx_Cj)  &
          -k(461)*n(idx_CHj)  &
          -k(1048)*n(idx_C)  &
          -k(847)*n(idx_OH)  &
          -k(846)*n(idx_O)  &
          -k(466)*n(idx_H2Oj)  &
          -k(210)  &
          -k(837)*n(idx_CH2)  &
          -k(836)*n(idx_C)  &
          -k(841)*n(idx_N)  &
          -k(209)  &
          -k(843)*n(idx_NH)  &
          -k(839)*n(idx_CH)  &
          -k(5)*n(idx_H2O)  &
          -k(838)*n(idx_CH3)  &
          -k(474)*n(idx_NH2j)  &
          +2.d0*k(4)*n(idx_H2)  &
          -k(460)*n(idx_Cj)  &
          +k(3)*n(idx_CH)  &
          -k(462)*n(idx_CH2j)  &
          -k(3)*n(idx_CH)  &
          -k(6)*n(idx_HOCj)  &
          -k(468)*n(idx_HEj)  &
          -4.d0*k(4)*n(idx_H2)  &
          -k(845)*n(idx_O2)  &
          -k(467)*n(idx_HCNj)  &
          -k(448)*n(idx_H2j)  &
          -k(477)*n(idx_OHj)  &
          -k(464)*n(idx_COj)  &
          -k(844)*n(idx_O2)  &
          -k(842)*n(idx_NH2)  &
          +k(8)*n(idx_OH)  &
          -k(100)*n(idx_HEj)  &
          -k(472)*n(idx_NHj)  &
          +k(6)*n(idx_HOCj)  &
          -k(11)*n(idx_H)
      pdj(7) =  &
          -k(836)*n(idx_C)  &
          -k(1048)*n(idx_C)  &
          +k(3)*n(idx_CH)
      pdj(8) =  &
          +k(474)*n(idx_NH2j)  &
          +k(5)*n(idx_H2O)  &
          +k(475)*n(idx_Oj)  &
          +4.d0*k(4)*n(idx_H2)  &
          +k(208)  &
          +k(839)*n(idx_CH)  &
          +k(838)*n(idx_CH3)  &
          +k(471)*n(idx_N2j)  &
          +k(477)*n(idx_OHj)  &
          +k(470)*n(idx_Nj)  &
          +k(842)*n(idx_NH2)  &
          +k(461)*n(idx_CHj)  &
          +k(465)*n(idx_COj)  &
          +k(462)*n(idx_CH2j)  &
          +k(837)*n(idx_CH2)  &
          +2.d0*k(9)*n(idx_E)  &
          +k(843)*n(idx_NH)  &
          -k(11)*n(idx_H)  &
          +k(467)*n(idx_HCNj)  &
          +k(3)*n(idx_CH)  &
          +k(460)*n(idx_Cj)  &
          +2.d0*k(210)  &
          +k(847)*n(idx_OH)  &
          +k(468)*n(idx_HEj)  &
          +k(473)*n(idx_NHj)  &
          +k(844)*n(idx_O2)  &
          +k(836)*n(idx_C)  &
          +k(841)*n(idx_N)  &
          +k(840)*n(idx_CN)  &
          +k(846)*n(idx_O)  &
          +3.d0*k(11)*n(idx_H)  &
          +k(448)*n(idx_H2j)  &
          +k(8)*n(idx_OH)  &
          +k(463)*n(idx_CNj)  &
          +k(464)*n(idx_COj)  &
          +k(466)*n(idx_H2Oj)
      pdj(9) =  &
          +k(847)*n(idx_OH)  &
          -k(5)*n(idx_H2O)
      pdj(10) =  &
          -k(8)*n(idx_OH)  &
          -k(847)*n(idx_OH)  &
          +2.d0*k(845)*n(idx_O2)  &
          +k(5)*n(idx_H2O)  &
          +k(846)*n(idx_O)
      pdj(11) =  &
          -k(844)*n(idx_O2)  &
          -k(7)*n(idx_O2)  &
          -k(845)*n(idx_O2)  &
          +k(476)*n(idx_O2Hj)
      pdj(12) =  &
          +k(839)*n(idx_CH)  &
          -k(837)*n(idx_CH2)  &
          +k(1048)*n(idx_C)
      pdj(16) =  &
          +k(842)*n(idx_NH2)
      pdj(18) =  &
          -k(840)*n(idx_CN)
      pdj(21) =  &
          -k(842)*n(idx_NH2)  &
          +k(843)*n(idx_NH)
      pdj(22) =  &
          +k(1049)*n(idx_CH)  &
          -k(838)*n(idx_CH3)  &
          +k(837)*n(idx_CH2)
      pdj(23) =  &
          +k(838)*n(idx_CH3)
      pdj(24) =  &
          +k(472)*n(idx_NHj)  &
          -k(841)*n(idx_N)
      pdj(25) =  &
          -k(843)*n(idx_NH)  &
          +k(841)*n(idx_N)
      pdj(26) =  &
          +k(469)*n(idx_HEHj)  &
          +k(468)*n(idx_HEj)  &
          +k(100)*n(idx_HEj)
      pdj(33) =  &
          +k(844)*n(idx_O2)
      pdj(54) =  &
          +k(464)*n(idx_COj)  &
          +k(6)*n(idx_HOCj)
      pdj(55) =  &
          +k(468)*n(idx_HEj)  &
          +k(208)
      pdj(56) =  &
          +k(465)*n(idx_COj)  &
          -k(6)*n(idx_HOCj)
      pdj(57) =  &
          -k(460)*n(idx_Cj)  &
          -k(1047)*n(idx_Cj)
      pdj(58) =  &
          -k(462)*n(idx_CH2j)  &
          +k(1047)*n(idx_Cj)  &
          +k(461)*n(idx_CHj)
      pdj(59) =  &
          -k(461)*n(idx_CHj)  &
          +k(460)*n(idx_Cj)
      pdj(62) =  &
          +k(474)*n(idx_NH2j)
      pdj(64) =  &
          -k(463)*n(idx_CNj)
      pdj(65) =  &
          -k(464)*n(idx_COj)  &
          -k(465)*n(idx_COj)
      pdj(66) =  &
          -k(471)*n(idx_N2j)
      pdj(68) =  &
          -k(466)*n(idx_H2Oj)  &
          +k(477)*n(idx_OHj)
      pdj(69) =  &
          +k(473)*n(idx_NHj)  &
          -k(474)*n(idx_NH2j)
      pdj(70) =  &
          -k(475)*n(idx_Oj)
      pdj(71) =  &
          -k(477)*n(idx_OHj)  &
          +k(475)*n(idx_Oj)
      pdj(72) =  &
          +k(462)*n(idx_CH2j)
      pdj(74) =  &
          -k(470)*n(idx_Nj)
      pdj(75) =  &
          -k(467)*n(idx_HCNj)  &
          +k(463)*n(idx_CNj)
      pdj(76) =  &
          -k(472)*n(idx_NHj)  &
          +k(470)*n(idx_Nj)  &
          -k(473)*n(idx_NHj)
      pdj(77) =  &
          -k(448)*n(idx_H2j)  &
          +k(100)*n(idx_HEj)  &
          +k(209)
      pdj(78) =  &
          -k(100)*n(idx_HEj)  &
          -k(468)*n(idx_HEj)
      pdj(81) =  &
          +k(469)*n(idx_HEHj)  &
          +k(476)*n(idx_O2Hj)  &
          +k(448)*n(idx_H2j)  &
          +k(472)*n(idx_NHj)
      pdj(83) =  &
          +k(466)*n(idx_H2Oj)
      pdj(84) =  &
          +k(467)*n(idx_HCNj)
      pdj(86) =  &
          -k(469)*n(idx_HEHj)
      pdj(87) =  &
          +k(471)*n(idx_N2j)
      pdj(88) =  &
          -k(476)*n(idx_O2Hj)
    elseif(j==7) then
      pdj(1) =  &
          +k(215)  &
          +k(206)  &
          +k(973)
      pdj(2) =  &
          +k(836)*n(idx_H2)  &
          +k(1051)*n(idx_H)  &
          +k(759)*n(idx_OH)  &
          +2.d0*k(746)*n(idx_CH2)  &
          +k(747)*n(idx_HCO)  &
          +k(751)*n(idx_NH2)  &
          +k(753)*n(idx_NH)
      pdj(3) =  &
          +k(756)*n(idx_O2)  &
          +k(754)*n(idx_NO)  &
          -k(1044)*n(idx_O)  &
          +k(337)*n(idx_O2j)  &
          +k(759)*n(idx_OH)  &
          +k(339)*n(idx_OHj)
      pdj(4) =  &
          +k(885)*n(idx_HNCO)  &
          +k(750)*n(idx_NH2)
      pdj(5) =  &
          +k(749)*n(idx_NH2)
      pdj(6) =  &
          +k(330)*n(idx_H3Oj)  &
          +k(502)*n(idx_H3j)  &
          -k(1048)*n(idx_H2)  &
          -k(836)*n(idx_H2)
      pdj(7) =  &
          -k(751)*n(idx_NH2)  &
          -k(338)*n(idx_O2Hj)  &
          -k(23)*n(idx_COj)  &
          -k(1043)*n(idx_Oj)  &
          -k(335)*n(idx_N2Hj)  &
          -k(750)*n(idx_NH2)  &
          -k(24)*n(idx_N2j)  &
          -k(757)*n(idx_OCN)  &
          -k(746)*n(idx_CH2)  &
          -k(1042)*n(idx_N)  &
          -k(755)*n(idx_NO)  &
          -k(206)  &
          -k(336)*n(idx_NHj)  &
          -k(973)  &
          -k(25)*n(idx_O2j)  &
          -k(748)*n(idx_N2)  &
          -k(333)*n(idx_HCO2j)  &
          -k(756)*n(idx_O2)  &
          -k(753)*n(idx_NH)  &
          -k(759)*n(idx_OH)  &
          -k(215)  &
          -k(332)*n(idx_HCOj)  &
          -k(122)*n(idx_HEj)  &
          -k(836)*n(idx_H2)  &
          -k(339)*n(idx_OHj)  &
          -k(758)*n(idx_OH)  &
          -k(331)*n(idx_HCNj)  &
          -k(1051)*n(idx_H)  &
          -k(502)*n(idx_H3j)  &
          -k(329)*n(idx_H2Oj)  &
          -k(441)*n(idx_H2j)  &
          -k(752)*n(idx_NH)  &
          -k(747)*n(idx_HCO)  &
          -k(337)*n(idx_O2j)  &
          -k(334)*n(idx_HNOj)  &
          -k(754)*n(idx_NO)  &
          -k(1048)*n(idx_H2)  &
          -k(1044)*n(idx_O)  &
          -k(22)*n(idx_CNj)  &
          -k(1070)  &
          -k(885)*n(idx_HNCO)  &
          -k(330)*n(idx_H3Oj)  &
          -k(749)*n(idx_NH2)
      pdj(8) =  &
          +k(749)*n(idx_NH2)  &
          +k(836)*n(idx_H2)  &
          -k(1051)*n(idx_H)  &
          +k(750)*n(idx_NH2)  &
          +k(441)*n(idx_H2j)  &
          +k(758)*n(idx_OH)  &
          +k(752)*n(idx_NH)
      pdj(10) =  &
          -k(758)*n(idx_OH)  &
          +k(329)*n(idx_H2Oj)  &
          -k(759)*n(idx_OH)
      pdj(11) =  &
          +k(338)*n(idx_O2Hj)  &
          -k(756)*n(idx_O2)  &
          +k(25)*n(idx_O2j)
      pdj(12) =  &
          -k(746)*n(idx_CH2)  &
          +k(1048)*n(idx_H2)
      pdj(14) =  &
          -k(747)*n(idx_HCO)
      pdj(17) =  &
          -k(754)*n(idx_NO)  &
          +k(334)*n(idx_HNOj)  &
          -k(755)*n(idx_NO)
      pdj(18) =  &
          +k(757)*n(idx_OCN)  &
          +k(22)*n(idx_CNj)  &
          +k(331)*n(idx_HCNj)  &
          +k(1042)*n(idx_N)  &
          +k(748)*n(idx_N2)  &
          +k(754)*n(idx_NO)  &
          +k(752)*n(idx_NH)
      pdj(19) =  &
          +k(757)*n(idx_OCN)  &
          +k(755)*n(idx_NO)  &
          +k(885)*n(idx_HNCO)  &
          +k(756)*n(idx_O2)  &
          +k(332)*n(idx_HCOj)  &
          +k(747)*n(idx_HCO)  &
          +k(23)*n(idx_COj)  &
          +k(1044)*n(idx_O)  &
          +k(758)*n(idx_OH)
      pdj(20) =  &
          -k(748)*n(idx_N2)  &
          +k(24)*n(idx_N2j)  &
          +k(335)*n(idx_N2Hj)
      pdj(21) =  &
          -k(751)*n(idx_NH2)  &
          -k(750)*n(idx_NH2)  &
          -k(749)*n(idx_NH2)
      pdj(24) =  &
          -k(1042)*n(idx_N)  &
          +k(753)*n(idx_NH)  &
          +k(748)*n(idx_N2)  &
          +k(336)*n(idx_NHj)  &
          +k(755)*n(idx_NO)
      pdj(25) =  &
          -k(753)*n(idx_NH)  &
          +k(751)*n(idx_NH2)  &
          -k(752)*n(idx_NH)
      pdj(26) =  &
          +k(122)*n(idx_HEj)
      pdj(29) =  &
          +k(333)*n(idx_HCO2j)
      pdj(31) =  &
          -k(885)*n(idx_HNCO)
      pdj(34) =  &
          -k(757)*n(idx_OCN)
      pdj(38) =  &
          +k(1070)
      pdj(54) =  &
          -k(332)*n(idx_HCOj)  &
          +k(330)*n(idx_H3Oj)
      pdj(57) =  &
          +k(122)*n(idx_HEj)  &
          +k(215)  &
          +k(206)  &
          +k(24)*n(idx_N2j)  &
          +k(973)  &
          +k(23)*n(idx_COj)  &
          +k(22)*n(idx_CNj)  &
          +k(25)*n(idx_O2j)
      pdj(59) =  &
          +k(331)*n(idx_HCNj)  &
          +k(332)*n(idx_HCOj)  &
          +k(335)*n(idx_N2Hj)  &
          +k(338)*n(idx_O2Hj)  &
          +k(333)*n(idx_HCO2j)  &
          +k(334)*n(idx_HNOj)  &
          +k(329)*n(idx_H2Oj)  &
          +k(339)*n(idx_OHj)  &
          +k(502)*n(idx_H3j)  &
          +k(336)*n(idx_NHj)  &
          +k(441)*n(idx_H2j)
      pdj(64) =  &
          -k(22)*n(idx_CNj)
      pdj(65) =  &
          +k(1043)*n(idx_Oj)  &
          +k(337)*n(idx_O2j)  &
          -k(23)*n(idx_COj)
      pdj(66) =  &
          -k(24)*n(idx_N2j)
      pdj(67) =  &
          -k(25)*n(idx_O2j)  &
          -k(337)*n(idx_O2j)
      pdj(68) =  &
          -k(329)*n(idx_H2Oj)
      pdj(70) =  &
          -k(1043)*n(idx_Oj)
      pdj(71) =  &
          -k(339)*n(idx_OHj)
      pdj(75) =  &
          -k(331)*n(idx_HCNj)
      pdj(76) =  &
          -k(336)*n(idx_NHj)
      pdj(77) =  &
          -k(441)*n(idx_H2j)
      pdj(78) =  &
          -k(122)*n(idx_HEj)
      pdj(79) =  &
          -k(334)*n(idx_HNOj)
      pdj(81) =  &
          -k(502)*n(idx_H3j)
      pdj(83) =  &
          -k(330)*n(idx_H3Oj)
      pdj(85) =  &
          -k(333)*n(idx_HCO2j)
      pdj(87) =  &
          -k(335)*n(idx_N2Hj)
      pdj(88) =  &
          -k(338)*n(idx_O2Hj)
    elseif(j==8) then
      pdj(1) =  &
          +k(211)  &
          +k(232)
      pdj(2) =  &
          +k(1051)*n(idx_C)  &
          -k(851)*n(idx_CH)  &
          +k(848)*n(idx_CH2)  &
          -k(10)*n(idx_CH)
      pdj(3) =  &
          -k(1052)*n(idx_O)  &
          +k(859)*n(idx_HCO)  &
          +k(874)*n(idx_OCN)  &
          +k(870)*n(idx_O2)  &
          +k(861)*n(idx_HNO)  &
          +k(14)*n(idx_OH)  &
          +k(877)*n(idx_OH)  &
          +2.d0*k(13)*n(idx_O2)  &
          +k(871)*n(idx_O2H)  &
          +k(868)*n(idx_NO)  &
          +k(115)*n(idx_Oj)
      pdj(4) =  &
          -k(860)*n(idx_HNC)
      pdj(5) =  &
          -k(857)*n(idx_HCN)  &
          +k(854)*n(idx_H2CN)  &
          +k(860)*n(idx_HNC)  &
          +k(874)*n(idx_OCN)  &
          +k(113)*n(idx_HCNj)
      pdj(6) =  &
          +k(858)*n(idx_HCO)  &
          +k(530)*n(idx_CHj)  &
          +k(855)*n(idx_H2CO)  &
          +k(857)*n(idx_HCN)  &
          +k(856)*n(idx_H2O)  &
          +k(865)*n(idx_NH3)  &
          +k(533)*n(idx_CH4j)  &
          -k(11)*n(idx_H2)  &
          +k(862)*n(idx_HNO)  &
          +k(872)*n(idx_O2H)  &
          +k(112)*n(idx_H2j)  &
          +k(850)*n(idx_CH4)  &
          +k(854)*n(idx_H2CN)  &
          +k(864)*n(idx_NH2)  &
          +k(532)*n(idx_CH3j)  &
          +k(877)*n(idx_OH)  &
          +k(531)*n(idx_CH2j)  &
          +k(849)*n(idx_CH3)  &
          +k(866)*n(idx_NH)  &
          +k(851)*n(idx_CH)  &
          +k(848)*n(idx_CH2)
      pdj(7) =  &
          +k(853)*n(idx_CO)  &
          +k(10)*n(idx_CH)  &
          +k(851)*n(idx_CH)  &
          -k(1051)*n(idx_C)
      pdj(8) =  &
          -k(849)*n(idx_CH3)  &
          -k(870)*n(idx_O2)  &
          -k(532)*n(idx_CH3j)  &
          -k(875)*n(idx_OCN)  &
          -k(866)*n(idx_NH)  &
          -k(534)*n(idx_HEHj)  &
          -k(848)*n(idx_CH2)  &
          +k(860)*n(idx_HNC)  &
          -k(12)*n(idx_H2O)  &
          -k(115)*n(idx_Oj)  &
          -k(865)*n(idx_NH3)  &
          -k(10)*n(idx_CH)  &
          -k(110)*n(idx_CNj)  &
          -k(860)*n(idx_HNC)  &
          -k(1050)*n(idx_Cj)  &
          +2.d0*k(10)*n(idx_CH)  &
          -k(531)*n(idx_CH2j)  &
          -k(868)*n(idx_NO)  &
          -k(855)*n(idx_H2CO)  &
          -k(862)*n(idx_HNO)  &
          -k(232)  &
          -k(876)*n(idx_OCN)  &
          -k(111)*n(idx_COj)  &
          +k(13)*n(idx_O2)  &
          -k(871)*n(idx_O2H)  &
          -k(863)*n(idx_HNO)  &
          -k(533)*n(idx_CH4j)  &
          -k(1053)*n(idx_OH)  &
          -k(856)*n(idx_H2O)  &
          -k(869)*n(idx_NO)  &
          -k(112)*n(idx_H2j)  &
          -k(13)*n(idx_O2)  &
          -k(852)*n(idx_CO2)  &
          -k(874)*n(idx_OCN)  &
          -k(858)*n(idx_HCO)  &
          -k(859)*n(idx_HCO)  &
          +3.d0*k(11)*n(idx_H2)  &
          -k(872)*n(idx_O2H)  &
          -k(877)*n(idx_OH)  &
          -k(854)*n(idx_H2CN)  &
          -k(1052)*n(idx_O)  &
          -k(114)*n(idx_HEj)  &
          -k(1045)*n(idx_Hj)  &
          -k(850)*n(idx_CH4)  &
          -k(14)*n(idx_OH)  &
          -k(11)*n(idx_H2)  &
          -k(857)*n(idx_HCN)  &
          -k(853)*n(idx_CO)  &
          -k(851)*n(idx_CH)  &
          -k(113)*n(idx_HCNj)  &
          -k(530)*n(idx_CHj)  &
          -k(861)*n(idx_HNO)  &
          -k(864)*n(idx_NH2)  &
          +2.d0*k(14)*n(idx_OH)  &
          -k(211)  &
          -k(873)*n(idx_O2H)  &
          +2.d0*k(12)*n(idx_H2O)  &
          -k(867)*n(idx_NO2)  &
          -k(1051)*n(idx_C)
      pdj(9) =  &
          -k(12)*n(idx_H2O)  &
          -k(856)*n(idx_H2O)  &
          +k(1053)*n(idx_OH)  &
          +k(871)*n(idx_O2H)
      pdj(10) =  &
          +k(856)*n(idx_H2O)  &
          +k(870)*n(idx_O2)  &
          -k(14)*n(idx_OH)  &
          +k(876)*n(idx_OCN)  &
          +k(1052)*n(idx_O)  &
          +k(852)*n(idx_CO2)  &
          -k(1053)*n(idx_OH)  &
          -k(877)*n(idx_OH)  &
          +k(867)*n(idx_NO2)  &
          +k(853)*n(idx_CO)  &
          +k(12)*n(idx_H2O)  &
          +k(863)*n(idx_HNO)  &
          +k(869)*n(idx_NO)  &
          +2.d0*k(873)*n(idx_O2H)
      pdj(11) =  &
          -k(13)*n(idx_O2)  &
          +k(872)*n(idx_O2H)  &
          -k(870)*n(idx_O2)
      pdj(12) =  &
          +k(859)*n(idx_HCO)  &
          +k(849)*n(idx_CH3)  &
          -k(848)*n(idx_CH2)
      pdj(13) =  &
          -k(855)*n(idx_H2CO)
      pdj(14) =  &
          -k(858)*n(idx_HCO)  &
          -k(859)*n(idx_HCO)  &
          +k(855)*n(idx_H2CO)
      pdj(16) =  &
          -k(865)*n(idx_NH3)
      pdj(17) =  &
          +k(867)*n(idx_NO2)  &
          +k(862)*n(idx_HNO)  &
          -k(868)*n(idx_NO)  &
          -k(869)*n(idx_NO)
      pdj(18) =  &
          +k(110)*n(idx_CNj)  &
          +k(857)*n(idx_HCN)  &
          +k(876)*n(idx_OCN)
      pdj(19) =  &
          +k(858)*n(idx_HCO)  &
          +k(875)*n(idx_OCN)  &
          +k(852)*n(idx_CO2)  &
          +k(111)*n(idx_COj)  &
          -k(853)*n(idx_CO)
      pdj(21) =  &
          +k(861)*n(idx_HNO)  &
          -k(864)*n(idx_NH2)  &
          +k(865)*n(idx_NH3)
      pdj(22) =  &
          +k(850)*n(idx_CH4)  &
          -k(849)*n(idx_CH3)
      pdj(23) =  &
          -k(850)*n(idx_CH4)
      pdj(24) =  &
          +k(866)*n(idx_NH)  &
          +k(869)*n(idx_NO)
      pdj(25) =  &
          +k(864)*n(idx_NH2)  &
          +k(868)*n(idx_NO)  &
          -k(866)*n(idx_NH)  &
          +k(863)*n(idx_HNO)  &
          +k(875)*n(idx_OCN)
      pdj(26) =  &
          +k(534)*n(idx_HEHj)  &
          +k(114)*n(idx_HEj)
      pdj(27) =  &
          -k(861)*n(idx_HNO)  &
          -k(862)*n(idx_HNO)  &
          -k(863)*n(idx_HNO)
      pdj(29) =  &
          -k(852)*n(idx_CO2)
      pdj(30) =  &
          -k(854)*n(idx_H2CN)
      pdj(32) =  &
          -k(867)*n(idx_NO2)
      pdj(33) =  &
          -k(873)*n(idx_O2H)  &
          -k(872)*n(idx_O2H)  &
          -k(871)*n(idx_O2H)
      pdj(34) =  &
          -k(876)*n(idx_OCN)  &
          -k(875)*n(idx_OCN)  &
          -k(874)*n(idx_OCN)
      pdj(55) =  &
          -k(1045)*n(idx_Hj)  &
          +k(114)*n(idx_HEj)  &
          +k(112)*n(idx_H2j)  &
          +k(113)*n(idx_HCNj)  &
          +k(110)*n(idx_CNj)  &
          +k(211)  &
          +k(232)  &
          +k(111)*n(idx_COj)  &
          +k(115)*n(idx_Oj)
      pdj(57) =  &
          +k(530)*n(idx_CHj)  &
          -k(1050)*n(idx_Cj)
      pdj(58) =  &
          +k(532)*n(idx_CH3j)  &
          -k(531)*n(idx_CH2j)
      pdj(59) =  &
          +k(531)*n(idx_CH2j)  &
          +k(1050)*n(idx_Cj)  &
          -k(530)*n(idx_CHj)
      pdj(64) =  &
          -k(110)*n(idx_CNj)
      pdj(65) =  &
          -k(111)*n(idx_COj)
      pdj(70) =  &
          -k(115)*n(idx_Oj)
      pdj(72) =  &
          +k(533)*n(idx_CH4j)  &
          -k(532)*n(idx_CH3j)
      pdj(73) =  &
          -k(533)*n(idx_CH4j)
      pdj(75) =  &
          -k(113)*n(idx_HCNj)
      pdj(77) =  &
          +k(534)*n(idx_HEHj)  &
          +k(1045)*n(idx_Hj)  &
          -k(112)*n(idx_H2j)
      pdj(78) =  &
          -k(114)*n(idx_HEj)
      pdj(86) =  &
          -k(534)*n(idx_HEHj)
    elseif(j==9) then
      pdj(1) =  &
          +k(1007)
      pdj(3) =  &
          +k(651)*n(idx_NHj)  &
          -k(942)*n(idx_O)  &
          +k(188)*n(idx_Oj)  &
          +k(730)*n(idx_OHj)
      pdj(5) =  &
          +k(108)*n(idx_HCNj)
      pdj(6) =  &
          +k(5)*n(idx_H2)  &
          -k(5)*n(idx_H2)  &
          +k(650)*n(idx_NHj)  &
          +k(856)*n(idx_H)  &
          +k(512)*n(idx_H3j)  &
          +k(348)*n(idx_CHj)  &
          +k(91)*n(idx_H2j)
      pdj(7) =  &
          +k(347)*n(idx_CHj)
      pdj(8) =  &
          +k(319)*n(idx_Cj)  &
          +k(231)  &
          +k(346)*n(idx_CHj)  &
          +k(320)*n(idx_Cj)  &
          +k(1008)  &
          +k(5)*n(idx_H2)  &
          +2.d0*k(12)*n(idx_H)  &
          +k(362)*n(idx_CH2j)  &
          +k(584)*n(idx_HEj)  &
          -k(12)*n(idx_H)  &
          +k(450)*n(idx_H2j)  &
          -k(856)*n(idx_H)  &
          +k(74)*n(idx_Hj)
      pdj(9) =  &
          -k(320)*n(idx_Cj)  &
          -k(198)*n(idx_OHj)  &
          -k(108)*n(idx_HCNj)  &
          -k(500)*n(idx_N2Hj)  &
          -k(1008)  &
          -k(650)*n(idx_NHj)  &
          -k(5)*n(idx_H2)  &
          -k(491)*n(idx_CNj)  &
          -k(392)*n(idx_CH4j)  &
          -k(512)*n(idx_H3j)  &
          -k(490)*n(idx_CNj)  &
          -k(652)*n(idx_NHj)  &
          -k(362)*n(idx_CH2j)  &
          -k(916)*n(idx_NH)  &
          -k(450)*n(idx_H2j)  &
          -k(651)*n(idx_NHj)  &
          -k(1077)  &
          -k(347)*n(idx_CHj)  &
          -k(107)*n(idx_COj)  &
          -k(126)*n(idx_HEj)  &
          -k(942)*n(idx_O)  &
          -k(156)*n(idx_NHj)  &
          -k(485)*n(idx_H2Oj)  &
          -k(493)*n(idx_H2COj)  &
          -k(492)*n(idx_COj)  &
          -k(188)*n(idx_Oj)  &
          -k(91)*n(idx_H2j)  &
          -k(585)*n(idx_HEj)  &
          -k(786)*n(idx_CH3)  &
          -k(856)*n(idx_H)  &
          -k(497)*n(idx_HCO2j)  &
          -k(730)*n(idx_OHj)  &
          -k(667)*n(idx_NH2j)  &
          -k(499)*n(idx_N2j)  &
          -k(12)*n(idx_H)  &
          -k(584)*n(idx_HEj)  &
          -k(140)*n(idx_Nj)  &
          -k(494)*n(idx_H3COj)  &
          -k(109)*n(idx_N2j)  &
          -k(348)*n(idx_CHj)  &
          -k(649)*n(idx_NHj)  &
          -k(1007)  &
          -k(74)*n(idx_Hj)  &
          -k(231)  &
          -k(666)*n(idx_NH2j)  &
          -k(501)*n(idx_O2Hj)  &
          -k(346)*n(idx_CHj)  &
          -k(495)*n(idx_HCNj)  &
          -k(319)*n(idx_Cj)  &
          -k(496)*n(idx_HCOj)  &
          -k(498)*n(idx_HNOj)
      pdj(10) =  &
          +k(198)*n(idx_OHj)  &
          +k(916)*n(idx_NH)  &
          +k(231)  &
          +k(1008)  &
          +2.d0*k(942)*n(idx_O)  &
          +k(5)*n(idx_H2)  &
          +k(12)*n(idx_H)  &
          +k(490)*n(idx_CNj)  &
          +k(485)*n(idx_H2Oj)  &
          +k(499)*n(idx_N2j)  &
          +k(786)*n(idx_CH3)  &
          +k(585)*n(idx_HEj)  &
          +k(667)*n(idx_NH2j)  &
          +k(492)*n(idx_COj)  &
          +k(652)*n(idx_NHj)  &
          +k(856)*n(idx_H)
      pdj(11) =  &
          +k(501)*n(idx_O2Hj)
      pdj(13) =  &
          +k(494)*n(idx_H3COj)
      pdj(14) =  &
          +k(493)*n(idx_H2COj)
      pdj(17) =  &
          +k(498)*n(idx_HNOj)
      pdj(18) =  &
          +k(495)*n(idx_HCNj)
      pdj(19) =  &
          +k(107)*n(idx_COj)  &
          +k(496)*n(idx_HCOj)
      pdj(20) =  &
          +k(109)*n(idx_N2j)  &
          +k(500)*n(idx_N2Hj)
      pdj(21) =  &
          +k(916)*n(idx_NH)
      pdj(22) =  &
          +k(392)*n(idx_CH4j)  &
          -k(786)*n(idx_CH3)
      pdj(23) =  &
          +k(786)*n(idx_CH3)
      pdj(24) =  &
          +k(140)*n(idx_Nj)  &
          +k(649)*n(idx_NHj)
      pdj(25) =  &
          +k(491)*n(idx_CNj)  &
          +k(156)*n(idx_NHj)  &
          -k(916)*n(idx_NH)  &
          +k(666)*n(idx_NH2j)
      pdj(26) =  &
          +k(126)*n(idx_HEj)  &
          +k(584)*n(idx_HEj)  &
          +k(585)*n(idx_HEj)
      pdj(29) =  &
          +k(497)*n(idx_HCO2j)
      pdj(40) =  &
          +k(1077)
      pdj(54) =  &
          +k(491)*n(idx_CNj)  &
          +k(348)*n(idx_CHj)  &
          +k(492)*n(idx_COj)  &
          -k(496)*n(idx_HCOj)  &
          +k(319)*n(idx_Cj)
      pdj(55) =  &
          -k(74)*n(idx_Hj)  &
          +k(585)*n(idx_HEj)
      pdj(56) =  &
          +k(320)*n(idx_Cj)
      pdj(57) =  &
          -k(320)*n(idx_Cj)  &
          -k(319)*n(idx_Cj)
      pdj(58) =  &
          -k(362)*n(idx_CH2j)
      pdj(59) =  &
          -k(348)*n(idx_CHj)  &
          -k(346)*n(idx_CHj)  &
          -k(347)*n(idx_CHj)
      pdj(60) =  &
          +k(346)*n(idx_CHj)  &
          -k(493)*n(idx_H2COj)
      pdj(62) =  &
          +k(651)*n(idx_NHj)  &
          +k(667)*n(idx_NH2j)
      pdj(64) =  &
          -k(490)*n(idx_CNj)  &
          -k(491)*n(idx_CNj)
      pdj(65) =  &
          -k(107)*n(idx_COj)  &
          -k(492)*n(idx_COj)
      pdj(66) =  &
          -k(499)*n(idx_N2j)  &
          -k(109)*n(idx_N2j)
      pdj(68) =  &
          +k(107)*n(idx_COj)  &
          +k(198)*n(idx_OHj)  &
          +k(108)*n(idx_HCNj)  &
          +k(156)*n(idx_NHj)  &
          +k(140)*n(idx_Nj)  &
          +k(188)*n(idx_Oj)  &
          +k(91)*n(idx_H2j)  &
          +k(109)*n(idx_N2j)  &
          +k(1007)  &
          +k(126)*n(idx_HEj)  &
          -k(485)*n(idx_H2Oj)  &
          +k(74)*n(idx_Hj)
      pdj(69) =  &
          -k(667)*n(idx_NH2j)  &
          +k(652)*n(idx_NHj)  &
          -k(666)*n(idx_NH2j)
      pdj(70) =  &
          -k(188)*n(idx_Oj)
      pdj(71) =  &
          -k(730)*n(idx_OHj)  &
          +k(584)*n(idx_HEj)  &
          -k(198)*n(idx_OHj)
      pdj(73) =  &
          -k(392)*n(idx_CH4j)
      pdj(74) =  &
          -k(140)*n(idx_Nj)
      pdj(75) =  &
          +k(490)*n(idx_CNj)  &
          -k(108)*n(idx_HCNj)  &
          -k(495)*n(idx_HCNj)
      pdj(76) =  &
          -k(156)*n(idx_NHj)  &
          -k(651)*n(idx_NHj)  &
          -k(649)*n(idx_NHj)  &
          -k(650)*n(idx_NHj)  &
          -k(652)*n(idx_NHj)
      pdj(77) =  &
          -k(91)*n(idx_H2j)  &
          -k(450)*n(idx_H2j)
      pdj(78) =  &
          -k(584)*n(idx_HEj)  &
          -k(585)*n(idx_HEj)  &
          -k(126)*n(idx_HEj)
      pdj(79) =  &
          -k(498)*n(idx_HNOj)  &
          +k(650)*n(idx_NHj)
      pdj(81) =  &
          -k(512)*n(idx_H3j)
      pdj(82) =  &
          -k(494)*n(idx_H3COj)  &
          +k(362)*n(idx_CH2j)
      pdj(83) =  &
          +k(500)*n(idx_N2Hj)  &
          +k(495)*n(idx_HCNj)  &
          +k(666)*n(idx_NH2j)  &
          +k(649)*n(idx_NHj)  &
          +k(485)*n(idx_H2Oj)  &
          +k(497)*n(idx_HCO2j)  &
          +k(347)*n(idx_CHj)  &
          +k(392)*n(idx_CH4j)  &
          +k(501)*n(idx_O2Hj)  &
          +k(512)*n(idx_H3j)  &
          +k(450)*n(idx_H2j)  &
          +k(493)*n(idx_H2COj)  &
          +k(498)*n(idx_HNOj)  &
          +k(494)*n(idx_H3COj)  &
          +k(496)*n(idx_HCOj)  &
          +k(730)*n(idx_OHj)
      pdj(85) =  &
          -k(497)*n(idx_HCO2j)
      pdj(87) =  &
          +k(499)*n(idx_N2j)  &
          -k(500)*n(idx_N2Hj)
      pdj(88) =  &
          -k(501)*n(idx_O2Hj)
    elseif(j==10) then
      pdj(1) =  &
          +k(1039)
      pdj(2) =  &
          +k(781)*n(idx_CH2)  &
          +k(759)*n(idx_C)  &
          -k(823)*n(idx_CH)
      pdj(3) =  &
          +k(961)*n(idx_CN)  &
          +k(739)*n(idx_H2Oj)  &
          +k(930)*n(idx_NH)  &
          +k(1038)  &
          -k(960)*n(idx_O)  &
          +k(799)*n(idx_CH3)  &
          +k(258)  &
          +k(782)*n(idx_CH2)  &
          +2.d0*k(972)*n(idx_OH)  &
          +k(738)*n(idx_COj)  &
          +k(912)*n(idx_NH2)  &
          +k(759)*n(idx_C)  &
          +k(193)*n(idx_Oj)  &
          +k(14)*n(idx_H)  &
          +k(907)*n(idx_N)  &
          +k(877)*n(idx_H)  &
          +k(8)*n(idx_H2)  &
          +k(737)*n(idx_OHj)
      pdj(5) =  &
          +k(961)*n(idx_CN)  &
          -k(966)*n(idx_HCN)  &
          -k(965)*n(idx_HCN)
      pdj(6) =  &
          +k(388)*n(idx_CH3j)  &
          -k(847)*n(idx_H2)  &
          +k(359)*n(idx_CHj)  &
          +k(526)*n(idx_H3j)  &
          -k(8)*n(idx_H2)  &
          +k(800)*n(idx_CH3)  &
          +k(99)*n(idx_H2j)  &
          +k(877)*n(idx_H)  &
          +k(8)*n(idx_H2)
      pdj(7) =  &
          -k(759)*n(idx_C)  &
          -k(758)*n(idx_C)
      pdj(8) =  &
          -k(877)*n(idx_H)  &
          +k(1038)  &
          +k(929)*n(idx_NH)  &
          -k(14)*n(idx_H)  &
          +k(328)*n(idx_Cj)  &
          +k(780)*n(idx_CH2)  &
          +k(258)  &
          +k(84)*n(idx_Hj)  &
          +k(459)*n(idx_H2j)  &
          +k(962)*n(idx_CN)  &
          +k(906)*n(idx_N)  &
          +2.d0*k(14)*n(idx_H)  &
          +k(758)*n(idx_C)  &
          +k(960)*n(idx_O)  &
          +k(847)*n(idx_H2)  &
          +k(8)*n(idx_H2)  &
          +k(714)*n(idx_Oj)  &
          -k(1053)*n(idx_H)  &
          +k(742)*n(idx_HCOj)  &
          +k(963)*n(idx_CO)  &
          +k(823)*n(idx_CH)  &
          +k(609)*n(idx_HEj)  &
          +k(970)*n(idx_NO)
      pdj(9) =  &
          +k(964)*n(idx_H2CO)  &
          +k(781)*n(idx_CH2)  &
          +k(965)*n(idx_HCN)  &
          +k(968)*n(idx_HNO)  &
          +k(969)*n(idx_NH3)  &
          +k(928)*n(idx_NH)  &
          +k(1053)*n(idx_H)  &
          +k(971)*n(idx_O2H)  &
          +2.d0*k(972)*n(idx_OH)  &
          +k(804)*n(idx_CH4)  &
          +k(847)*n(idx_H2)  &
          +k(801)*n(idx_CH3)  &
          +k(967)*n(idx_HCO)  &
          +k(911)*n(idx_NH2)
      pdj(10) =  &
          -k(877)*n(idx_H)  &
          -k(928)*n(idx_NH)  &
          -k(526)*n(idx_H3j)  &
          -k(969)*n(idx_NH3)  &
          -k(758)*n(idx_C)  &
          -k(971)*n(idx_O2H)  &
          -k(14)*n(idx_H)  &
          -k(911)*n(idx_NH2)  &
          -k(780)*n(idx_CH2)  &
          -k(963)*n(idx_CO)  &
          -k(804)*n(idx_CH4)  &
          -k(149)*n(idx_Nj)  &
          -k(960)*n(idx_O)  &
          -k(801)*n(idx_CH3)  &
          -k(799)*n(idx_CH3)  &
          -k(961)*n(idx_CN)  &
          -k(929)*n(idx_NH)  &
          -k(1038)  &
          -k(912)*n(idx_NH2)  &
          -k(99)*n(idx_H2j)  &
          -k(738)*n(idx_COj)  &
          -k(84)*n(idx_Hj)  &
          -k(388)*n(idx_CH3j)  &
          -k(741)*n(idx_HCOj)  &
          -k(823)*n(idx_CH)  &
          -k(359)*n(idx_CHj)  &
          -k(714)*n(idx_Oj)  &
          -k(742)*n(idx_HCOj)  &
          -k(205)*n(idx_N2j)  &
          -k(800)*n(idx_CH3)  &
          -k(737)*n(idx_OHj)  &
          -k(782)*n(idx_CH2)  &
          -k(930)*n(idx_NH)  &
          -4.d0*k(972)*n(idx_OH)  &
          -k(966)*n(idx_HCN)  &
          -k(459)*n(idx_H2j)  &
          -k(328)*n(idx_Cj)  &
          -k(743)*n(idx_HNOj)  &
          -k(8)*n(idx_H2)  &
          -k(258)  &
          -k(907)*n(idx_N)  &
          -k(193)*n(idx_Oj)  &
          -k(609)*n(idx_HEj)  &
          -k(745)*n(idx_O2Hj)  &
          -k(962)*n(idx_CN)  &
          -k(1053)*n(idx_H)  &
          -k(740)*n(idx_HCNj)  &
          -k(739)*n(idx_H2Oj)  &
          -k(964)*n(idx_H2CO)  &
          -k(663)*n(idx_NHj)  &
          -k(1039)  &
          -k(968)*n(idx_HNO)  &
          -k(759)*n(idx_C)  &
          -k(906)*n(idx_N)  &
          -k(847)*n(idx_H2)  &
          -k(967)*n(idx_HCO)  &
          -k(1074)  &
          -k(965)*n(idx_HCN)  &
          -k(203)*n(idx_CNj)  &
          -k(781)*n(idx_CH2)  &
          -k(204)*n(idx_COj)  &
          -k(744)*n(idx_N2Hj)  &
          -k(970)*n(idx_NO)
      pdj(11) =  &
          +k(971)*n(idx_O2H)  &
          +k(745)*n(idx_O2Hj)  &
          +k(960)*n(idx_O)
      pdj(12) =  &
          -k(780)*n(idx_CH2)  &
          -k(781)*n(idx_CH2)  &
          +k(801)*n(idx_CH3)  &
          -k(782)*n(idx_CH2)
      pdj(13) =  &
          +k(780)*n(idx_CH2)  &
          -k(964)*n(idx_H2CO)  &
          +k(800)*n(idx_CH3)
      pdj(14) =  &
          +k(964)*n(idx_H2CO)  &
          -k(967)*n(idx_HCO)  &
          +k(823)*n(idx_CH)
      pdj(16) =  &
          +k(912)*n(idx_NH2)  &
          -k(969)*n(idx_NH3)
      pdj(17) =  &
          +k(743)*n(idx_HNOj)  &
          +k(906)*n(idx_N)  &
          +k(968)*n(idx_HNO)  &
          -k(970)*n(idx_NO)
      pdj(18) =  &
          -k(961)*n(idx_CN)  &
          -k(962)*n(idx_CN)  &
          +k(740)*n(idx_HCNj)  &
          +k(203)*n(idx_CNj)  &
          +k(965)*n(idx_HCN)
      pdj(19) =  &
          -k(963)*n(idx_CO)  &
          +k(204)*n(idx_COj)  &
          +k(741)*n(idx_HCOj)  &
          +k(966)*n(idx_HCN)  &
          +k(758)*n(idx_C)  &
          +k(967)*n(idx_HCO)
      pdj(20) =  &
          +k(205)*n(idx_N2j)  &
          +k(744)*n(idx_N2Hj)
      pdj(21) =  &
          -k(912)*n(idx_NH2)  &
          +k(969)*n(idx_NH3)  &
          +k(930)*n(idx_NH)  &
          +k(966)*n(idx_HCN)  &
          -k(911)*n(idx_NH2)
      pdj(22) =  &
          -k(800)*n(idx_CH3)  &
          -k(799)*n(idx_CH3)  &
          +k(782)*n(idx_CH2)  &
          -k(801)*n(idx_CH3)  &
          +k(804)*n(idx_CH4)
      pdj(23) =  &
          +k(799)*n(idx_CH3)  &
          -k(804)*n(idx_CH4)
      pdj(24) =  &
          +k(928)*n(idx_NH)  &
          +k(663)*n(idx_NHj)  &
          -k(907)*n(idx_N)  &
          -k(906)*n(idx_N)  &
          +k(149)*n(idx_Nj)
      pdj(25) =  &
          -k(928)*n(idx_NH)  &
          -k(930)*n(idx_NH)  &
          +k(907)*n(idx_N)  &
          +k(911)*n(idx_NH2)  &
          -k(929)*n(idx_NH)
      pdj(26) =  &
          +k(609)*n(idx_HEj)
      pdj(27) =  &
          +k(929)*n(idx_NH)  &
          -k(968)*n(idx_HNO)
      pdj(29) =  &
          +k(963)*n(idx_CO)
      pdj(32) =  &
          +k(970)*n(idx_NO)
      pdj(33) =  &
          -k(971)*n(idx_O2H)
      pdj(34) =  &
          +k(962)*n(idx_CN)
      pdj(40) =  &
          +k(1074)
      pdj(54) =  &
          +k(738)*n(idx_COj)  &
          -k(741)*n(idx_HCOj)  &
          -k(742)*n(idx_HCOj)
      pdj(55) =  &
          -k(84)*n(idx_Hj)
      pdj(57) =  &
          -k(328)*n(idx_Cj)
      pdj(59) =  &
          -k(359)*n(idx_CHj)
      pdj(60) =  &
          +k(388)*n(idx_CH3j)
      pdj(64) =  &
          -k(203)*n(idx_CNj)
      pdj(65) =  &
          +k(359)*n(idx_CHj)  &
          +k(328)*n(idx_Cj)  &
          -k(204)*n(idx_COj)  &
          -k(738)*n(idx_COj)
      pdj(66) =  &
          -k(205)*n(idx_N2j)
      pdj(67) =  &
          +k(714)*n(idx_Oj)
      pdj(68) =  &
          +k(745)*n(idx_O2Hj)  &
          -k(739)*n(idx_H2Oj)  &
          +k(741)*n(idx_HCOj)  &
          +k(744)*n(idx_N2Hj)  &
          +k(459)*n(idx_H2j)  &
          +k(663)*n(idx_NHj)  &
          +k(526)*n(idx_H3j)  &
          +k(743)*n(idx_HNOj)  &
          +k(737)*n(idx_OHj)  &
          +k(740)*n(idx_HCNj)
      pdj(70) =  &
          -k(714)*n(idx_Oj)  &
          +k(609)*n(idx_HEj)  &
          -k(193)*n(idx_Oj)
      pdj(71) =  &
          +k(1039)  &
          +k(149)*n(idx_Nj)  &
          +k(204)*n(idx_COj)  &
          -k(737)*n(idx_OHj)  &
          +k(205)*n(idx_N2j)  &
          +k(84)*n(idx_Hj)  &
          +k(193)*n(idx_Oj)  &
          +k(99)*n(idx_H2j)  &
          +k(203)*n(idx_CNj)
      pdj(72) =  &
          -k(388)*n(idx_CH3j)
      pdj(74) =  &
          -k(149)*n(idx_Nj)
      pdj(75) =  &
          -k(740)*n(idx_HCNj)
      pdj(76) =  &
          -k(663)*n(idx_NHj)
      pdj(77) =  &
          -k(99)*n(idx_H2j)  &
          -k(459)*n(idx_H2j)
      pdj(78) =  &
          -k(609)*n(idx_HEj)
      pdj(79) =  &
          -k(743)*n(idx_HNOj)
      pdj(81) =  &
          -k(526)*n(idx_H3j)
      pdj(83) =  &
          +k(739)*n(idx_H2Oj)
      pdj(85) =  &
          +k(742)*n(idx_HCOj)
      pdj(87) =  &
          -k(744)*n(idx_N2Hj)
      pdj(88) =  &
          -k(745)*n(idx_O2Hj)
    elseif(j==11) then
      pdj(1) =  &
          +k(253)  &
          +k(1032)
      pdj(2) =  &
          -k(815)*n(idx_CH)  &
          -k(816)*n(idx_CH)  &
          -k(817)*n(idx_CH)  &
          -k(818)*n(idx_CH)
      pdj(3) =  &
          +2.d0*k(1033)  &
          +k(870)*n(idx_H)  &
          +k(626)*n(idx_Nj)  &
          +k(192)*n(idx_Oj)  &
          +k(756)*n(idx_C)  &
          +k(356)*n(idx_CHj)  &
          +k(924)*n(idx_NH)  &
          +k(834)*n(idx_CO)  &
          +k(774)*n(idx_CH2)  &
          +k(904)*n(idx_N)  &
          +k(672)*n(idx_NH2j)  &
          +k(818)*n(idx_CH)  &
          +k(831)*n(idx_CN)  &
          +k(325)*n(idx_Cj)  &
          +k(606)*n(idx_HEj)  &
          +k(385)*n(idx_CH3j)  &
          +2.d0*k(7)*n(idx_H2)  &
          +2.d0*k(13)*n(idx_H)  &
          +k(816)*n(idx_CH)  &
          +k(932)*n(idx_NO)  &
          +2.d0*k(254)
      pdj(5) =  &
          +k(117)*n(idx_HCNj)
      pdj(6) =  &
          -k(7)*n(idx_H2)  &
          -k(844)*n(idx_H2)  &
          -k(845)*n(idx_H2)  &
          +k(98)*n(idx_H2j)  &
          +k(771)*n(idx_CH2)  &
          +k(7)*n(idx_H2)  &
          +k(523)*n(idx_H3j)
      pdj(7) =  &
          -k(756)*n(idx_C)
      pdj(8) =  &
          +k(82)*n(idx_Hj)  &
          -k(870)*n(idx_H)  &
          +k(844)*n(idx_H2)  &
          +k(815)*n(idx_CH)  &
          +k(816)*n(idx_CH)  &
          +2.d0*k(772)*n(idx_CH2)  &
          -k(13)*n(idx_H)  &
          +k(457)*n(idx_H2j)  &
          +k(13)*n(idx_H)
      pdj(9) =  &
          +k(773)*n(idx_CH2)  &
          +k(794)*n(idx_CH3)  &
          +k(106)*n(idx_H2Oj)
      pdj(10) =  &
          +k(870)*n(idx_H)  &
          +k(775)*n(idx_CH2)  &
          +k(817)*n(idx_CH)  &
          +2.d0*k(845)*n(idx_H2)  &
          +k(660)*n(idx_NHj)  &
          +k(882)*n(idx_HCO)  &
          +k(925)*n(idx_NH)  &
          +k(364)*n(idx_CH2j)  &
          +k(673)*n(idx_NH2j)  &
          +k(355)*n(idx_CHj)  &
          +k(793)*n(idx_CH3)  &
          +k(202)*n(idx_OHj)
      pdj(11) =  &
          -k(7)*n(idx_H2)  &
          -k(816)*n(idx_CH)  &
          -k(795)*n(idx_CH3)  &
          -k(773)*n(idx_CH2)  &
          -k(818)*n(idx_CH)  &
          -k(771)*n(idx_CH2)  &
          -k(253)  &
          -k(883)*n(idx_HCO)  &
          -k(756)*n(idx_C)  &
          -k(98)*n(idx_H2j)  &
          -k(523)*n(idx_H3j)  &
          -k(803)*n(idx_CH4)  &
          -k(793)*n(idx_CH3)  &
          -k(831)*n(idx_CN)  &
          -k(606)*n(idx_HEj)  &
          -k(192)*n(idx_Oj)  &
          -k(932)*n(idx_NO)  &
          -k(129)*n(idx_HEj)  &
          -k(845)*n(idx_H2)  &
          -k(1112)  &
          -k(117)*n(idx_HCNj)  &
          -k(420)*n(idx_CNj)  &
          -k(45)*n(idx_CH4j)  &
          -k(934)*n(idx_OCN)  &
          -k(834)*n(idx_CO)  &
          -k(355)*n(idx_CHj)  &
          -k(357)*n(idx_CHj)  &
          -k(13)*n(idx_H)  &
          -k(67)*n(idx_COj)  &
          -k(1032)  &
          -k(62)*n(idx_CNj)  &
          -k(661)*n(idx_NHj)  &
          -k(844)*n(idx_H2)  &
          -k(148)*n(idx_Nj)  &
          -k(673)*n(idx_NH2j)  &
          -k(325)*n(idx_Cj)  &
          -k(925)*n(idx_NH)  &
          -k(385)*n(idx_CH3j)  &
          -k(870)*n(idx_H)  &
          -k(202)*n(idx_OHj)  &
          -k(660)*n(idx_NHj)  &
          -k(935)*n(idx_OCN)  &
          -k(364)*n(idx_CH2j)  &
          -k(775)*n(idx_CH2)  &
          -k(82)*n(idx_Hj)  &
          -k(772)*n(idx_CH2)  &
          -k(924)*n(idx_NH)  &
          -k(794)*n(idx_CH3)  &
          -k(153)*n(idx_N2j)  &
          -k(457)*n(idx_H2j)  &
          -k(904)*n(idx_N)  &
          -k(479)*n(idx_H2COj)  &
          -k(626)*n(idx_Nj)  &
          -k(356)*n(idx_CHj)  &
          -k(774)*n(idx_CH2)  &
          -k(882)*n(idx_HCO)  &
          -k(106)*n(idx_H2Oj)  &
          -k(672)*n(idx_NH2j)  &
          -k(817)*n(idx_CH)  &
          -k(326)*n(idx_Cj)  &
          -k(627)*n(idx_Nj)  &
          -k(815)*n(idx_CH)  &
          -k(1033)  &
          -k(159)*n(idx_NHj)  &
          -k(254)  &
          -k(830)*n(idx_CN)
      pdj(12) =  &
          -k(773)*n(idx_CH2)  &
          -k(771)*n(idx_CH2)  &
          -k(774)*n(idx_CH2)  &
          +k(795)*n(idx_CH3)  &
          -k(775)*n(idx_CH2)  &
          -k(772)*n(idx_CH2)
      pdj(13) =  &
          +k(793)*n(idx_CH3)  &
          +k(774)*n(idx_CH2)
      pdj(14) =  &
          +k(775)*n(idx_CH2)  &
          -k(883)*n(idx_HCO)  &
          -k(882)*n(idx_HCO)  &
          +k(818)*n(idx_CH)  &
          +k(357)*n(idx_CHj)  &
          +k(794)*n(idx_CH3)
      pdj(17) =  &
          +k(627)*n(idx_Nj)  &
          -k(932)*n(idx_NO)  &
          +k(925)*n(idx_NH)  &
          +k(830)*n(idx_CN)  &
          +k(904)*n(idx_N)  &
          +k(934)*n(idx_OCN)
      pdj(18) =  &
          +k(62)*n(idx_CNj)  &
          -k(831)*n(idx_CN)  &
          -k(830)*n(idx_CN)
      pdj(19) =  &
          +k(67)*n(idx_COj)  &
          +k(816)*n(idx_CH)  &
          +k(817)*n(idx_CH)  &
          +k(883)*n(idx_HCO)  &
          -k(834)*n(idx_CO)  &
          +k(420)*n(idx_CNj)  &
          +k(326)*n(idx_Cj)  &
          +k(935)*n(idx_OCN)  &
          +k(830)*n(idx_CN)  &
          +k(773)*n(idx_CH2)  &
          +k(756)*n(idx_C)
      pdj(20) =  &
          +k(153)*n(idx_N2j)
      pdj(22) =  &
          -k(794)*n(idx_CH3)  &
          +k(803)*n(idx_CH4)  &
          -k(793)*n(idx_CH3)  &
          -k(795)*n(idx_CH3)
      pdj(23) =  &
          +k(45)*n(idx_CH4j)  &
          -k(803)*n(idx_CH4)
      pdj(24) =  &
          -k(904)*n(idx_N)  &
          +k(661)*n(idx_NHj)  &
          +k(148)*n(idx_Nj)
      pdj(25) =  &
          +k(159)*n(idx_NHj)  &
          -k(924)*n(idx_NH)  &
          -k(925)*n(idx_NH)
      pdj(26) =  &
          +k(606)*n(idx_HEj)  &
          +k(129)*n(idx_HEj)
      pdj(27) =  &
          +k(924)*n(idx_NH)
      pdj(29) =  &
          +k(815)*n(idx_CH)  &
          +k(882)*n(idx_HCO)  &
          +k(772)*n(idx_CH2)  &
          +k(834)*n(idx_CO)  &
          +k(771)*n(idx_CH2)  &
          +k(934)*n(idx_OCN)
      pdj(32) =  &
          +k(932)*n(idx_NO)  &
          +k(935)*n(idx_OCN)
      pdj(33) =  &
          +k(803)*n(idx_CH4)  &
          +k(883)*n(idx_HCO)  &
          +k(844)*n(idx_H2)  &
          +k(479)*n(idx_H2COj)  &
          +k(795)*n(idx_CH3)
      pdj(34) =  &
          -k(934)*n(idx_OCN)  &
          -k(935)*n(idx_OCN)  &
          +k(831)*n(idx_CN)
      pdj(46) =  &
          +k(1112)
      pdj(54) =  &
          +k(364)*n(idx_CH2j)  &
          +k(356)*n(idx_CHj)  &
          +k(479)*n(idx_H2COj)
      pdj(55) =  &
          -k(82)*n(idx_Hj)
      pdj(57) =  &
          -k(325)*n(idx_Cj)  &
          -k(326)*n(idx_Cj)
      pdj(58) =  &
          -k(364)*n(idx_CH2j)
      pdj(59) =  &
          -k(355)*n(idx_CHj)  &
          -k(357)*n(idx_CHj)  &
          -k(356)*n(idx_CHj)
      pdj(60) =  &
          -k(479)*n(idx_H2COj)
      pdj(63) =  &
          +k(420)*n(idx_CNj)  &
          +k(626)*n(idx_Nj)  &
          +k(660)*n(idx_NHj)
      pdj(64) =  &
          -k(420)*n(idx_CNj)  &
          -k(62)*n(idx_CNj)
      pdj(65) =  &
          +k(355)*n(idx_CHj)  &
          +k(325)*n(idx_Cj)  &
          -k(67)*n(idx_COj)
      pdj(66) =  &
          -k(153)*n(idx_N2j)
      pdj(67) =  &
          +k(62)*n(idx_CNj)  &
          +k(67)*n(idx_COj)  &
          +k(82)*n(idx_Hj)  &
          +k(117)*n(idx_HCNj)  &
          +k(148)*n(idx_Nj)  &
          +k(192)*n(idx_Oj)  &
          +k(159)*n(idx_NHj)  &
          +k(153)*n(idx_N2j)  &
          +k(253)  &
          +k(1032)  &
          +k(98)*n(idx_H2j)  &
          +k(45)*n(idx_CH4j)  &
          +k(106)*n(idx_H2Oj)  &
          +k(202)*n(idx_OHj)  &
          +k(129)*n(idx_HEj)
      pdj(68) =  &
          -k(106)*n(idx_H2Oj)
      pdj(69) =  &
          -k(673)*n(idx_NH2j)  &
          -k(672)*n(idx_NH2j)
      pdj(70) =  &
          +k(627)*n(idx_Nj)  &
          +k(357)*n(idx_CHj)  &
          +k(326)*n(idx_Cj)  &
          -k(192)*n(idx_Oj)  &
          +k(606)*n(idx_HEj)
      pdj(71) =  &
          -k(202)*n(idx_OHj)
      pdj(72) =  &
          -k(385)*n(idx_CH3j)
      pdj(73) =  &
          -k(45)*n(idx_CH4j)
      pdj(74) =  &
          -k(627)*n(idx_Nj)  &
          -k(626)*n(idx_Nj)  &
          -k(148)*n(idx_Nj)
      pdj(75) =  &
          -k(117)*n(idx_HCNj)
      pdj(76) =  &
          -k(660)*n(idx_NHj)  &
          -k(159)*n(idx_NHj)  &
          -k(661)*n(idx_NHj)
      pdj(77) =  &
          -k(98)*n(idx_H2j)  &
          -k(457)*n(idx_H2j)
      pdj(78) =  &
          -k(129)*n(idx_HEj)  &
          -k(606)*n(idx_HEj)
      pdj(79) =  &
          +k(673)*n(idx_NH2j)
      pdj(80) =  &
          +k(672)*n(idx_NH2j)
      pdj(81) =  &
          -k(523)*n(idx_H3j)
      pdj(82) =  &
          +k(385)*n(idx_CH3j)
      pdj(88) =  &
          +k(661)*n(idx_NHj)  &
          +k(457)*n(idx_H2j)  &
          +k(523)*n(idx_H3j)
    elseif(j==12) then
      pdj(1) =  &
          +k(978)  &
          +k(217)
      pdj(2) =  &
          +k(888)*n(idx_N)  &
          +2.d0*k(746)*n(idx_C)  &
          +k(848)*n(idx_H)  &
          +k(366)*n(idx_COj)  &
          +k(781)*n(idx_OH)  &
          +k(979)  &
          +k(779)*n(idx_O)  &
          +k(762)*n(idx_CN)  &
          +2.d0*k(760)*n(idx_CH2)  &
          +k(218)
      pdj(3) =  &
          -k(776)*n(idx_O)  &
          +k(774)*n(idx_O2)  &
          +k(37)*n(idx_Oj)  &
          +k(381)*n(idx_OHj)  &
          +k(379)*n(idx_O2j)  &
          -k(777)*n(idx_O)  &
          -k(778)*n(idx_O)  &
          +k(782)*n(idx_OH)  &
          -k(779)*n(idx_O)
      pdj(4) =  &
          +k(372)*n(idx_HCNHj)  &
          +k(887)*n(idx_N)
      pdj(5) =  &
          +k(769)*n(idx_NO)  &
          +k(766)*n(idx_N2)  &
          +k(762)*n(idx_CN)  &
          +k(371)*n(idx_HCNHj)  &
          +k(886)*n(idx_N)
      pdj(6) =  &
          +k(85)*n(idx_H2j)  &
          +k(503)*n(idx_H3j)  &
          +k(848)*n(idx_H)  &
          +k(564)*n(idx_HEj)  &
          -k(837)*n(idx_H2)  &
          +k(428)*n(idx_Hj)  &
          +k(776)*n(idx_O)  &
          +k(771)*n(idx_O2)
      pdj(7) =  &
          +k(15)*n(idx_Cj)  &
          -k(746)*n(idx_C)
      pdj(8) =  &
          +k(770)*n(idx_NO)  &
          +k(565)*n(idx_HEj)  &
          +2.d0*k(777)*n(idx_O)  &
          +k(69)*n(idx_Hj)  &
          +k(837)*n(idx_H2)  &
          +k(886)*n(idx_N)  &
          +k(887)*n(idx_N)  &
          +2.d0*k(772)*n(idx_O2)  &
          +k(979)  &
          +k(778)*n(idx_O)  &
          +k(442)*n(idx_H2j)  &
          +k(780)*n(idx_OH)  &
          -k(848)*n(idx_H)  &
          +k(218)
      pdj(9) =  &
          +k(369)*n(idx_H3Oj)  &
          +k(34)*n(idx_H2Oj)  &
          +k(773)*n(idx_O2)  &
          +k(781)*n(idx_OH)
      pdj(10) =  &
          -k(782)*n(idx_OH)  &
          +k(39)*n(idx_OHj)  &
          +k(775)*n(idx_O2)  &
          -k(781)*n(idx_OH)  &
          +k(779)*n(idx_O)  &
          +k(368)*n(idx_H2Oj)  &
          +k(769)*n(idx_NO)  &
          -k(780)*n(idx_OH)
      pdj(11) =  &
          -k(775)*n(idx_O2)  &
          +k(380)*n(idx_O2Hj)  &
          -k(771)*n(idx_O2)  &
          -k(773)*n(idx_O2)  &
          +k(38)*n(idx_O2j)  &
          -k(774)*n(idx_O2)  &
          -k(772)*n(idx_O2)
      pdj(12) =  &
          -k(782)*n(idx_OH)  &
          -k(775)*n(idx_O2)  &
          -k(771)*n(idx_O2)  &
          -k(34)*n(idx_H2Oj)  &
          -k(366)*n(idx_COj)  &
          -k(375)*n(idx_N2Hj)  &
          -k(368)*n(idx_H2Oj)  &
          -k(32)*n(idx_COj)  &
          -k(565)*n(idx_HEj)  &
          -k(31)*n(idx_CNj)  &
          -k(764)*n(idx_HCO)  &
          -k(763)*n(idx_H2CO)  &
          -k(379)*n(idx_O2j)  &
          -k(373)*n(idx_HCOj)  &
          -k(746)*n(idx_C)  &
          -k(37)*n(idx_Oj)  &
          -k(35)*n(idx_N2j)  &
          -k(15)*n(idx_Cj)  &
          -k(780)*n(idx_OH)  &
          -k(376)*n(idx_NHj)  &
          -k(979)  &
          -k(377)*n(idx_NH2j)  &
          -k(371)*n(idx_HCNHj)  &
          -k(33)*n(idx_H2COj)  &
          -k(777)*n(idx_O)  &
          -k(762)*n(idx_CN)  &
          -k(773)*n(idx_O2)  &
          -k(369)*n(idx_H3Oj)  &
          -k(218)  &
          -k(778)*n(idx_O)  &
          -k(767)*n(idx_NO2)  &
          -k(769)*n(idx_NO)  &
          -k(38)*n(idx_O2j)  &
          -k(442)*n(idx_H2j)  &
          -k(887)*n(idx_N)  &
          -k(372)*n(idx_HCNHj)  &
          -k(69)*n(idx_Hj)  &
          -k(217)  &
          -k(978)  &
          -k(564)*n(idx_HEj)  &
          -k(370)*n(idx_HCNj)  &
          -k(135)*n(idx_Nj)  &
          -k(766)*n(idx_N2)  &
          -k(381)*n(idx_OHj)  &
          -4.d0*k(760)*n(idx_CH2)  &
          -k(39)*n(idx_OHj)  &
          -k(85)*n(idx_H2j)  &
          -k(848)*n(idx_H)  &
          -k(779)*n(idx_O)  &
          -k(367)*n(idx_H2COj)  &
          -k(428)*n(idx_Hj)  &
          -k(770)*n(idx_NO)  &
          -k(765)*n(idx_HNO)  &
          -k(1076)  &
          -k(374)*n(idx_HNOj)  &
          -k(503)*n(idx_H3j)  &
          -k(776)*n(idx_O)  &
          -k(774)*n(idx_O2)  &
          -k(781)*n(idx_OH)  &
          -k(36)*n(idx_NH2j)  &
          -k(378)*n(idx_NH3j)  &
          -k(380)*n(idx_O2Hj)  &
          -k(888)*n(idx_N)  &
          -k(886)*n(idx_N)  &
          -k(761)*n(idx_CH4)  &
          -k(772)*n(idx_O2)  &
          -k(768)*n(idx_NO)  &
          -k(837)*n(idx_H2)
      pdj(13) =  &
          -k(763)*n(idx_H2CO)  &
          +k(33)*n(idx_H2COj)  &
          +k(774)*n(idx_O2)  &
          +k(768)*n(idx_NO)  &
          +k(780)*n(idx_OH)  &
          +k(767)*n(idx_NO2)
      pdj(14) =  &
          +k(775)*n(idx_O2)  &
          +k(778)*n(idx_O)  &
          +k(763)*n(idx_H2CO)  &
          +k(367)*n(idx_H2COj)  &
          -k(764)*n(idx_HCO)
      pdj(17) =  &
          -k(770)*n(idx_NO)  &
          -k(768)*n(idx_NO)  &
          +k(765)*n(idx_HNO)  &
          +k(374)*n(idx_HNOj)  &
          -k(769)*n(idx_NO)  &
          +k(767)*n(idx_NO2)
      pdj(18) =  &
          -k(762)*n(idx_CN)  &
          +k(370)*n(idx_HCNj)  &
          +k(31)*n(idx_CNj)
      pdj(19) =  &
          +k(32)*n(idx_COj)  &
          +k(777)*n(idx_O)  &
          +k(764)*n(idx_HCO)  &
          +k(776)*n(idx_O)  &
          +k(373)*n(idx_HCOj)  &
          +k(773)*n(idx_O2)
      pdj(20) =  &
          -k(766)*n(idx_N2)  &
          +k(375)*n(idx_N2Hj)  &
          +k(35)*n(idx_N2j)
      pdj(21) =  &
          +k(378)*n(idx_NH3j)  &
          +k(36)*n(idx_NH2j)
      pdj(22) =  &
          +2.d0*k(761)*n(idx_CH4)  &
          +k(764)*n(idx_HCO)  &
          +k(837)*n(idx_H2)  &
          +k(765)*n(idx_HNO)  &
          +2.d0*k(760)*n(idx_CH2)  &
          +k(782)*n(idx_OH)  &
          +k(763)*n(idx_H2CO)
      pdj(23) =  &
          -k(761)*n(idx_CH4)
      pdj(24) =  &
          +k(376)*n(idx_NHj)  &
          +k(135)*n(idx_Nj)  &
          -k(886)*n(idx_N)  &
          -k(888)*n(idx_N)  &
          +k(768)*n(idx_NO)  &
          -k(887)*n(idx_N)
      pdj(25) =  &
          +k(377)*n(idx_NH2j)  &
          +k(766)*n(idx_N2)  &
          +k(888)*n(idx_N)
      pdj(26) =  &
          +k(564)*n(idx_HEj)  &
          +k(565)*n(idx_HEj)
      pdj(27) =  &
          -k(765)*n(idx_HNO)
      pdj(29) =  &
          +k(771)*n(idx_O2)  &
          +k(772)*n(idx_O2)
      pdj(31) =  &
          +k(770)*n(idx_NO)
      pdj(32) =  &
          -k(767)*n(idx_NO2)
      pdj(38) =  &
          +k(1076)
      pdj(54) =  &
          -k(373)*n(idx_HCOj)  &
          +k(366)*n(idx_COj)
      pdj(55) =  &
          -k(428)*n(idx_Hj)  &
          -k(69)*n(idx_Hj)
      pdj(57) =  &
          +k(564)*n(idx_HEj)  &
          -k(15)*n(idx_Cj)
      pdj(58) =  &
          +k(33)*n(idx_H2COj)  &
          +k(34)*n(idx_H2Oj)  &
          +k(39)*n(idx_OHj)  &
          +k(35)*n(idx_N2j)  &
          +k(135)*n(idx_Nj)  &
          +k(37)*n(idx_Oj)  &
          +k(38)*n(idx_O2j)  &
          +k(69)*n(idx_Hj)  &
          +k(978)  &
          +k(217)  &
          +k(15)*n(idx_Cj)  &
          +k(32)*n(idx_COj)  &
          +k(31)*n(idx_CNj)  &
          +k(36)*n(idx_NH2j)  &
          +k(85)*n(idx_H2j)
      pdj(59) =  &
          +k(565)*n(idx_HEj)  &
          +k(428)*n(idx_Hj)
      pdj(60) =  &
          -k(33)*n(idx_H2COj)  &
          +k(379)*n(idx_O2j)  &
          -k(367)*n(idx_H2COj)
      pdj(62) =  &
          -k(378)*n(idx_NH3j)
      pdj(64) =  &
          -k(31)*n(idx_CNj)
      pdj(65) =  &
          -k(32)*n(idx_COj)  &
          -k(366)*n(idx_COj)
      pdj(66) =  &
          -k(35)*n(idx_N2j)
      pdj(67) =  &
          -k(379)*n(idx_O2j)  &
          -k(38)*n(idx_O2j)
      pdj(68) =  &
          -k(34)*n(idx_H2Oj)  &
          -k(368)*n(idx_H2Oj)
      pdj(69) =  &
          -k(36)*n(idx_NH2j)  &
          -k(377)*n(idx_NH2j)
      pdj(70) =  &
          -k(37)*n(idx_Oj)
      pdj(71) =  &
          -k(381)*n(idx_OHj)  &
          -k(39)*n(idx_OHj)
      pdj(72) =  &
          +k(378)*n(idx_NH3j)  &
          +k(376)*n(idx_NHj)  &
          +k(380)*n(idx_O2Hj)  &
          +k(371)*n(idx_HCNHj)  &
          +k(503)*n(idx_H3j)  &
          +k(381)*n(idx_OHj)  &
          +k(369)*n(idx_H3Oj)  &
          +k(375)*n(idx_N2Hj)  &
          +k(372)*n(idx_HCNHj)  &
          +k(367)*n(idx_H2COj)  &
          +k(377)*n(idx_NH2j)  &
          +k(368)*n(idx_H2Oj)  &
          +k(373)*n(idx_HCOj)  &
          +k(374)*n(idx_HNOj)  &
          +k(370)*n(idx_HCNj)  &
          +k(442)*n(idx_H2j)
      pdj(74) =  &
          -k(135)*n(idx_Nj)
      pdj(75) =  &
          -k(370)*n(idx_HCNj)
      pdj(76) =  &
          -k(376)*n(idx_NHj)
      pdj(77) =  &
          -k(442)*n(idx_H2j)  &
          -k(85)*n(idx_H2j)
      pdj(78) =  &
          -k(564)*n(idx_HEj)  &
          -k(565)*n(idx_HEj)
      pdj(79) =  &
          -k(374)*n(idx_HNOj)
      pdj(81) =  &
          -k(503)*n(idx_H3j)
      pdj(83) =  &
          -k(369)*n(idx_H3Oj)
      pdj(84) =  &
          -k(372)*n(idx_HCNHj)  &
          -k(371)*n(idx_HCNHj)
      pdj(87) =  &
          -k(375)*n(idx_N2Hj)
      pdj(88) =  &
          -k(380)*n(idx_O2Hj)
    elseif(j==13) then
      pdj(1) =  &
          +k(1005)  &
          +k(1004)
      pdj(2) =  &
          -k(806)*n(idx_CH)  &
          +k(318)*n(idx_Cj)
      pdj(3) =  &
          +k(187)*n(idx_Oj)  &
          +k(729)*n(idx_OHj)  &
          +k(583)*n(idx_HEj)  &
          -k(941)*n(idx_O)
      pdj(4) =  &
          +k(549)*n(idx_HCNHj)
      pdj(5) =  &
          +k(418)*n(idx_CNj)  &
          +k(548)*n(idx_HCNHj)  &
          +k(824)*n(idx_CN)
      pdj(6) =  &
          +k(230)  &
          +k(449)*n(idx_H2j)  &
          +k(581)*n(idx_HEj)  &
          +k(434)*n(idx_Hj)  &
          +k(90)*n(idx_H2j)  &
          +k(855)*n(idx_H)  &
          +k(1002)  &
          +k(435)*n(idx_Hj)  &
          +k(511)*n(idx_H3j)
      pdj(7) =  &
          +k(344)*n(idx_CHj)  &
          +k(17)*n(idx_Cj)
      pdj(8) =  &
          +k(628)*n(idx_N2j)  &
          +k(481)*n(idx_O2j)  &
          +k(449)*n(idx_H2j)  &
          -k(855)*n(idx_H)  &
          +k(434)*n(idx_Hj)  &
          +k(582)*n(idx_HEj)  &
          +k(73)*n(idx_Hj)  &
          +2.d0*k(1003)  &
          +k(1005)
      pdj(9) =  &
          +k(964)*n(idx_OH)  &
          +k(527)*n(idx_H3Oj)  &
          +k(102)*n(idx_H2Oj)
      pdj(10) =  &
          +k(197)*n(idx_OHj)  &
          -k(964)*n(idx_OH)  &
          +k(941)*n(idx_O)  &
          +k(484)*n(idx_H2Oj)  &
          +k(708)*n(idx_Oj)
      pdj(11) =  &
          +k(481)*n(idx_O2j)  &
          +k(482)*n(idx_O2Hj)  &
          +k(101)*n(idx_O2j)
      pdj(12) =  &
          +k(806)*n(idx_CH)  &
          +k(620)*n(idx_Nj)  &
          -k(763)*n(idx_CH2)  &
          +k(345)*n(idx_CHj)
      pdj(13) =  &
          -k(17)*n(idx_Cj)  &
          -k(435)*n(idx_Hj)  &
          -k(125)*n(idx_HEj)  &
          -k(187)*n(idx_Oj)  &
          -k(824)*n(idx_CN)  &
          -k(785)*n(idx_CH3)  &
          -k(619)*n(idx_Nj)  &
          -k(511)*n(idx_H3j)  &
          -k(391)*n(idx_CH4j)  &
          -k(361)*n(idx_CH2j)  &
          -k(729)*n(idx_OHj)  &
          -k(434)*n(idx_Hj)  &
          -k(1002)  &
          -k(581)*n(idx_HEj)  &
          -k(423)*n(idx_COj)  &
          -k(550)*n(idx_HCOj)  &
          -k(583)*n(idx_HEj)  &
          -k(139)*n(idx_Nj)  &
          -k(763)*n(idx_CH2)  &
          -k(418)*n(idx_CNj)  &
          -k(230)  &
          -k(480)*n(idx_HNOj)  &
          -k(549)*n(idx_HCNHj)  &
          -k(647)*n(idx_NHj)  &
          -k(964)*n(idx_OH)  &
          -k(345)*n(idx_CHj)  &
          -k(43)*n(idx_CH4j)  &
          -k(478)*n(idx_H2COj)  &
          -k(482)*n(idx_O2Hj)  &
          -k(1004)  &
          -k(58)*n(idx_CNj)  &
          -k(648)*n(idx_NHj)  &
          -k(548)*n(idx_HCNHj)  &
          -k(343)*n(idx_CHj)  &
          -k(317)*n(idx_Cj)  &
          -k(537)*n(idx_HCNj)  &
          -k(582)*n(idx_HEj)  &
          -k(101)*n(idx_O2j)  &
          -k(90)*n(idx_H2j)  &
          -k(318)*n(idx_Cj)  &
          -k(484)*n(idx_H2Oj)  &
          -k(806)*n(idx_CH)  &
          -k(64)*n(idx_COj)  &
          -k(1003)  &
          -k(102)*n(idx_H2Oj)  &
          -k(708)*n(idx_Oj)  &
          -k(481)*n(idx_O2j)  &
          -k(1005)  &
          -k(383)*n(idx_CH3j)  &
          -k(665)*n(idx_NH2j)  &
          -k(941)*n(idx_O)  &
          -k(150)*n(idx_N2j)  &
          -k(664)*n(idx_NH2j)  &
          -k(633)*n(idx_N2Hj)  &
          -k(155)*n(idx_NHj)  &
          -k(620)*n(idx_Nj)  &
          -k(855)*n(idx_H)  &
          -k(344)*n(idx_CHj)  &
          -k(1072)  &
          -k(197)*n(idx_OHj)  &
          -k(527)*n(idx_H3Oj)  &
          -k(449)*n(idx_H2j)  &
          -k(73)*n(idx_Hj)  &
          -k(628)*n(idx_N2j)
      pdj(14) =  &
          +k(785)*n(idx_CH3)  &
          +k(941)*n(idx_O)  &
          +k(665)*n(idx_NH2j)  &
          +k(806)*n(idx_CH)  &
          +k(478)*n(idx_H2COj)  &
          +k(964)*n(idx_OH)  &
          +k(855)*n(idx_H)  &
          +k(824)*n(idx_CN)  &
          +k(423)*n(idx_COj)  &
          +k(763)*n(idx_CH2)
      pdj(17) =  &
          +k(480)*n(idx_HNOj)
      pdj(18) =  &
          +k(537)*n(idx_HCNj)  &
          -k(824)*n(idx_CN)  &
          +k(58)*n(idx_CNj)
      pdj(19) =  &
          +k(230)  &
          +k(1003)  &
          +k(343)*n(idx_CHj)  &
          +k(317)*n(idx_Cj)  &
          +k(1002)  &
          +k(550)*n(idx_HCOj)  &
          +k(64)*n(idx_COj)
      pdj(20) =  &
          +k(628)*n(idx_N2j)  &
          +k(150)*n(idx_N2j)  &
          +k(633)*n(idx_N2Hj)
      pdj(21) =  &
          +k(648)*n(idx_NHj)
      pdj(22) =  &
          -k(785)*n(idx_CH3)  &
          +k(391)*n(idx_CH4j)  &
          +k(361)*n(idx_CH2j)  &
          +k(763)*n(idx_CH2)
      pdj(23) =  &
          +k(43)*n(idx_CH4j)  &
          +k(785)*n(idx_CH3)  &
          +k(383)*n(idx_CH3j)
      pdj(24) =  &
          +k(647)*n(idx_NHj)  &
          +k(139)*n(idx_Nj)
      pdj(25) =  &
          +k(664)*n(idx_NH2j)  &
          +k(155)*n(idx_NHj)  &
          +k(619)*n(idx_Nj)
      pdj(26) =  &
          +k(581)*n(idx_HEj)  &
          +k(582)*n(idx_HEj)  &
          +k(125)*n(idx_HEj)  &
          +k(583)*n(idx_HEj)
      pdj(37) =  &
          +k(1072)
      pdj(54) =  &
          +k(628)*n(idx_N2j)  &
          +k(648)*n(idx_NHj)  &
          +k(481)*n(idx_O2j)  &
          +k(383)*n(idx_CH3j)  &
          +k(449)*n(idx_H2j)  &
          +k(582)*n(idx_HEj)  &
          +k(361)*n(idx_CH2j)  &
          +k(1005)  &
          +k(418)*n(idx_CNj)  &
          +k(619)*n(idx_Nj)  &
          +k(708)*n(idx_Oj)  &
          +k(435)*n(idx_Hj)  &
          +k(423)*n(idx_COj)  &
          -k(550)*n(idx_HCOj)  &
          +k(318)*n(idx_Cj)  &
          +k(345)*n(idx_CHj)
      pdj(55) =  &
          -k(434)*n(idx_Hj)  &
          -k(435)*n(idx_Hj)  &
          -k(73)*n(idx_Hj)
      pdj(57) =  &
          -k(17)*n(idx_Cj)  &
          -k(318)*n(idx_Cj)  &
          -k(317)*n(idx_Cj)
      pdj(58) =  &
          -k(361)*n(idx_CH2j)  &
          +k(583)*n(idx_HEj)  &
          +k(317)*n(idx_Cj)
      pdj(59) =  &
          -k(344)*n(idx_CHj)  &
          -k(345)*n(idx_CHj)  &
          -k(343)*n(idx_CHj)
      pdj(60) =  &
          +k(43)*n(idx_CH4j)  &
          +k(197)*n(idx_OHj)  &
          +k(125)*n(idx_HEj)  &
          +k(90)*n(idx_H2j)  &
          +k(58)*n(idx_CNj)  &
          +k(102)*n(idx_H2Oj)  &
          -k(478)*n(idx_H2COj)  &
          +k(73)*n(idx_Hj)  &
          +k(1004)  &
          +k(155)*n(idx_NHj)  &
          +k(101)*n(idx_O2j)  &
          +k(187)*n(idx_Oj)  &
          +k(17)*n(idx_Cj)  &
          +k(139)*n(idx_Nj)  &
          +k(150)*n(idx_N2j)  &
          +k(64)*n(idx_COj)
      pdj(62) =  &
          +k(665)*n(idx_NH2j)
      pdj(63) =  &
          +k(620)*n(idx_Nj)
      pdj(64) =  &
          -k(58)*n(idx_CNj)  &
          -k(418)*n(idx_CNj)
      pdj(65) =  &
          -k(64)*n(idx_COj)  &
          +k(434)*n(idx_Hj)  &
          +k(581)*n(idx_HEj)  &
          -k(423)*n(idx_COj)
      pdj(66) =  &
          -k(150)*n(idx_N2j)  &
          -k(628)*n(idx_N2j)
      pdj(67) =  &
          -k(481)*n(idx_O2j)  &
          -k(101)*n(idx_O2j)
      pdj(68) =  &
          -k(484)*n(idx_H2Oj)  &
          -k(102)*n(idx_H2Oj)
      pdj(69) =  &
          -k(665)*n(idx_NH2j)  &
          -k(664)*n(idx_NH2j)
      pdj(70) =  &
          -k(708)*n(idx_Oj)  &
          -k(187)*n(idx_Oj)
      pdj(71) =  &
          -k(729)*n(idx_OHj)  &
          -k(197)*n(idx_OHj)
      pdj(72) =  &
          -k(383)*n(idx_CH3j)  &
          +k(343)*n(idx_CHj)
      pdj(73) =  &
          -k(391)*n(idx_CH4j)  &
          -k(43)*n(idx_CH4j)
      pdj(74) =  &
          -k(139)*n(idx_Nj)  &
          -k(619)*n(idx_Nj)  &
          -k(620)*n(idx_Nj)
      pdj(75) =  &
          -k(537)*n(idx_HCNj)
      pdj(76) =  &
          -k(647)*n(idx_NHj)  &
          -k(648)*n(idx_NHj)  &
          -k(155)*n(idx_NHj)
      pdj(77) =  &
          -k(90)*n(idx_H2j)  &
          -k(449)*n(idx_H2j)
      pdj(78) =  &
          -k(583)*n(idx_HEj)  &
          -k(582)*n(idx_HEj)  &
          -k(125)*n(idx_HEj)  &
          -k(581)*n(idx_HEj)
      pdj(79) =  &
          -k(480)*n(idx_HNOj)
      pdj(81) =  &
          -k(511)*n(idx_H3j)
      pdj(82) =  &
          +k(537)*n(idx_HCNj)  &
          +k(391)*n(idx_CH4j)  &
          +k(633)*n(idx_N2Hj)  &
          +k(480)*n(idx_HNOj)  &
          +k(664)*n(idx_NH2j)  &
          +k(482)*n(idx_O2Hj)  &
          +k(484)*n(idx_H2Oj)  &
          +k(478)*n(idx_H2COj)  &
          +k(548)*n(idx_HCNHj)  &
          +k(647)*n(idx_NHj)  &
          +k(344)*n(idx_CHj)  &
          +k(729)*n(idx_OHj)  &
          +k(550)*n(idx_HCOj)  &
          +k(511)*n(idx_H3j)  &
          +k(527)*n(idx_H3Oj)  &
          +k(549)*n(idx_HCNHj)
      pdj(83) =  &
          -k(527)*n(idx_H3Oj)
      pdj(84) =  &
          -k(549)*n(idx_HCNHj)  &
          -k(548)*n(idx_HCNHj)
      pdj(87) =  &
          -k(633)*n(idx_N2Hj)
      pdj(88) =  &
          -k(482)*n(idx_O2Hj)
    elseif(j==14) then
      pdj(1) =  &
          +k(235)  &
          +k(1014)
      pdj(2) =  &
          +k(747)*n(idx_C)  &
          -k(807)*n(idx_CH)  &
          +k(26)*n(idx_CHj)
      pdj(3) =  &
          -k(947)*n(idx_O)  &
          +k(592)*n(idx_HEj)  &
          +k(859)*n(idx_H)  &
          +k(733)*n(idx_OHj)  &
          +k(189)*n(idx_Oj)  &
          -k(946)*n(idx_O)  &
          +k(896)*n(idx_N)
      pdj(5) =  &
          +k(896)*n(idx_N)  &
          +k(825)*n(idx_CN)
      pdj(6) =  &
          +k(514)*n(idx_H3j)  &
          +k(858)*n(idx_H)  &
          +k(436)*n(idx_Hj)  &
          +k(93)*n(idx_H2j)  &
          +2.d0*k(878)*n(idx_HCO)
      pdj(7) =  &
          -k(747)*n(idx_C)  &
          +k(18)*n(idx_Cj)
      pdj(8) =  &
          +k(76)*n(idx_Hj)  &
          +k(1013)  &
          -k(858)*n(idx_H)  &
          +k(897)*n(idx_N)  &
          -k(859)*n(idx_H)  &
          +k(234)  &
          +k(946)*n(idx_O)  &
          +k(590)*n(idx_HEj)
      pdj(9) =  &
          +k(967)*n(idx_OH)  &
          +k(103)*n(idx_H2Oj)
      pdj(10) =  &
          +k(199)*n(idx_OHj)  &
          +k(882)*n(idx_O2)  &
          +k(488)*n(idx_H2Oj)  &
          -k(967)*n(idx_OH)  &
          +k(947)*n(idx_O)
      pdj(11) =  &
          +k(121)*n(idx_O2j)  &
          -k(882)*n(idx_O2)  &
          +k(884)*n(idx_O2H)  &
          +k(556)*n(idx_O2Hj)  &
          -k(883)*n(idx_O2)
      pdj(12) =  &
          +k(807)*n(idx_CH)  &
          +k(859)*n(idx_H)  &
          -k(764)*n(idx_CH2)
      pdj(13) =  &
          +k(880)*n(idx_HNO)  &
          +2.d0*k(879)*n(idx_HCO)  &
          +k(120)*n(idx_H2COj)  &
          +k(884)*n(idx_O2H)
      pdj(14) =  &
          -k(189)*n(idx_Oj)  &
          -k(151)*n(idx_N2j)  &
          -k(621)*n(idx_Nj)  &
          -k(76)*n(idx_Hj)  &
          -k(40)*n(idx_CH3j)  &
          -k(947)*n(idx_O)  &
          -k(1081)  &
          -k(858)*n(idx_H)  &
          -k(160)*n(idx_NH2j)  &
          -k(807)*n(idx_CH)  &
          -k(18)*n(idx_Cj)  &
          -k(487)*n(idx_H2Oj)  &
          -k(591)*n(idx_HEj)  &
          -k(65)*n(idx_COj)  &
          -k(514)*n(idx_H3j)  &
          -k(946)*n(idx_O)  &
          -k(733)*n(idx_OHj)  &
          -k(1013)  &
          -k(554)*n(idx_N2Hj)  &
          -k(60)*n(idx_CNj)  &
          -k(93)*n(idx_H2j)  &
          -k(881)*n(idx_NO)  &
          -k(363)*n(idx_CH2j)  &
          -k(967)*n(idx_OH)  &
          -k(120)*n(idx_H2COj)  &
          -k(142)*n(idx_Nj)  &
          -k(437)*n(idx_Hj)  &
          -k(553)*n(idx_HNOj)  &
          -k(555)*n(idx_O2j)  &
          -4.d0*k(878)*n(idx_HCO)  &
          -k(880)*n(idx_HNO)  &
          -k(551)*n(idx_HCOj)  &
          -k(732)*n(idx_OHj)  &
          -k(825)*n(idx_CN)  &
          -k(436)*n(idx_Hj)  &
          -k(556)*n(idx_O2Hj)  &
          -k(884)*n(idx_O2H)  &
          -k(883)*n(idx_O2)  &
          -k(234)  &
          -k(321)*n(idx_Cj)  &
          -k(711)*n(idx_Oj)  &
          -k(654)*n(idx_NHj)  &
          -k(747)*n(idx_C)  &
          -k(199)*n(idx_OHj)  &
          -k(764)*n(idx_CH2)  &
          -k(896)*n(idx_N)  &
          -k(539)*n(idx_HCNj)  &
          -k(787)*n(idx_CH3)  &
          -k(590)*n(idx_HEj)  &
          -k(26)*n(idx_CHj)  &
          -k(103)*n(idx_H2Oj)  &
          -k(488)*n(idx_H2Oj)  &
          -k(882)*n(idx_O2)  &
          -k(895)*n(idx_N)  &
          -k(592)*n(idx_HEj)  &
          -k(552)*n(idx_H2COj)  &
          -k(669)*n(idx_NH2j)  &
          -k(540)*n(idx_HCNj)  &
          -k(629)*n(idx_N2j)  &
          -k(384)*n(idx_CH3j)  &
          -k(169)*n(idx_NH3j)  &
          -k(121)*n(idx_O2j)  &
          -k(859)*n(idx_H)  &
          -k(1014)  &
          -k(897)*n(idx_N)  &
          -k(451)*n(idx_H2j)  &
          -k(419)*n(idx_CNj)  &
          -k(235)  &
          -4.d0*k(879)*n(idx_HCO)  &
          -k(350)*n(idx_CHj)
      pdj(16) =  &
          +k(169)*n(idx_NH3j)
      pdj(17) =  &
          +k(880)*n(idx_HNO)  &
          +k(553)*n(idx_HNOj)  &
          -k(881)*n(idx_NO)
      pdj(18) =  &
          -k(825)*n(idx_CN)  &
          +k(60)*n(idx_CNj)  &
          +k(539)*n(idx_HCNj)
      pdj(19) =  &
          +k(967)*n(idx_OH)  &
          +k(732)*n(idx_OHj)  &
          +k(551)*n(idx_HCOj)  &
          +k(419)*n(idx_CNj)  &
          +k(825)*n(idx_CN)  &
          +k(711)*n(idx_Oj)  &
          +k(487)*n(idx_H2Oj)  &
          +k(363)*n(idx_CH2j)  &
          +4.d0*k(878)*n(idx_HCO)  &
          +k(787)*n(idx_CH3)  &
          +k(895)*n(idx_N)  &
          +k(591)*n(idx_HEj)  &
          +k(883)*n(idx_O2)  &
          +k(384)*n(idx_CH3j)  &
          +k(540)*n(idx_HCNj)  &
          +k(65)*n(idx_COj)  &
          +k(234)  &
          +k(552)*n(idx_H2COj)  &
          +k(555)*n(idx_O2j)  &
          +k(881)*n(idx_NO)  &
          +k(451)*n(idx_H2j)  &
          +k(747)*n(idx_C)  &
          +k(621)*n(idx_Nj)  &
          +k(807)*n(idx_CH)  &
          +k(321)*n(idx_Cj)  &
          +k(858)*n(idx_H)  &
          +k(764)*n(idx_CH2)  &
          +k(1013)  &
          +k(947)*n(idx_O)  &
          +2.d0*k(879)*n(idx_HCO)  &
          +k(350)*n(idx_CHj)  &
          +k(629)*n(idx_N2j)  &
          +k(437)*n(idx_Hj)
      pdj(20) =  &
          +k(151)*n(idx_N2j)  &
          +k(554)*n(idx_N2Hj)
      pdj(21) =  &
          +k(160)*n(idx_NH2j)
      pdj(22) =  &
          -k(787)*n(idx_CH3)  &
          +k(764)*n(idx_CH2)  &
          +k(40)*n(idx_CH3j)
      pdj(23) =  &
          +k(787)*n(idx_CH3)
      pdj(24) =  &
          -k(897)*n(idx_N)  &
          -k(896)*n(idx_N)  &
          +k(142)*n(idx_Nj)  &
          +k(654)*n(idx_NHj)  &
          -k(895)*n(idx_N)
      pdj(25) =  &
          +k(895)*n(idx_N)  &
          +k(669)*n(idx_NH2j)
      pdj(26) =  &
          +k(592)*n(idx_HEj)  &
          +k(590)*n(idx_HEj)
      pdj(27) =  &
          -k(880)*n(idx_HNO)  &
          +k(881)*n(idx_NO)
      pdj(29) =  &
          +k(946)*n(idx_O)  &
          +k(882)*n(idx_O2)
      pdj(33) =  &
          +k(883)*n(idx_O2)  &
          -k(884)*n(idx_O2H)
      pdj(34) =  &
          +k(897)*n(idx_N)
      pdj(37) =  &
          +k(1081)
      pdj(54) =  &
          +k(169)*n(idx_NH3j)  &
          +k(76)*n(idx_Hj)  &
          +k(103)*n(idx_H2Oj)  &
          +k(60)*n(idx_CNj)  &
          +k(40)*n(idx_CH3j)  &
          +k(121)*n(idx_O2j)  &
          +k(142)*n(idx_Nj)  &
          +k(151)*n(idx_N2j)  &
          -k(551)*n(idx_HCOj)  &
          +k(26)*n(idx_CHj)  &
          +k(160)*n(idx_NH2j)  &
          +k(235)  &
          +k(189)*n(idx_Oj)  &
          +k(18)*n(idx_Cj)  &
          +k(199)*n(idx_OHj)  &
          +k(1014)  &
          +k(120)*n(idx_H2COj)  &
          +k(65)*n(idx_COj)  &
          +k(93)*n(idx_H2j)
      pdj(55) =  &
          -k(437)*n(idx_Hj)  &
          -k(436)*n(idx_Hj)  &
          -k(76)*n(idx_Hj)
      pdj(57) =  &
          -k(321)*n(idx_Cj)  &
          -k(18)*n(idx_Cj)
      pdj(58) =  &
          +k(350)*n(idx_CHj)  &
          -k(363)*n(idx_CH2j)
      pdj(59) =  &
          -k(26)*n(idx_CHj)  &
          +k(592)*n(idx_HEj)  &
          +k(321)*n(idx_Cj)  &
          -k(350)*n(idx_CHj)
      pdj(60) =  &
          -k(552)*n(idx_H2COj)  &
          +k(556)*n(idx_O2Hj)  &
          +k(554)*n(idx_N2Hj)  &
          +k(551)*n(idx_HCOj)  &
          +k(514)*n(idx_H3j)  &
          +k(733)*n(idx_OHj)  &
          +k(488)*n(idx_H2Oj)  &
          +k(654)*n(idx_NHj)  &
          +k(553)*n(idx_HNOj)  &
          -k(120)*n(idx_H2COj)  &
          +k(669)*n(idx_NH2j)  &
          +k(539)*n(idx_HCNj)
      pdj(62) =  &
          -k(169)*n(idx_NH3j)
      pdj(64) =  &
          -k(60)*n(idx_CNj)  &
          -k(419)*n(idx_CNj)
      pdj(65) =  &
          -k(65)*n(idx_COj)  &
          +k(590)*n(idx_HEj)  &
          +k(436)*n(idx_Hj)
      pdj(66) =  &
          -k(151)*n(idx_N2j)  &
          -k(629)*n(idx_N2j)
      pdj(67) =  &
          -k(555)*n(idx_O2j)  &
          -k(121)*n(idx_O2j)
      pdj(68) =  &
          +k(732)*n(idx_OHj)  &
          -k(103)*n(idx_H2Oj)  &
          -k(488)*n(idx_H2Oj)  &
          -k(487)*n(idx_H2Oj)
      pdj(69) =  &
          -k(160)*n(idx_NH2j)  &
          -k(669)*n(idx_NH2j)
      pdj(70) =  &
          -k(189)*n(idx_Oj)  &
          -k(711)*n(idx_Oj)
      pdj(71) =  &
          -k(199)*n(idx_OHj)  &
          +k(711)*n(idx_Oj)  &
          -k(733)*n(idx_OHj)  &
          -k(732)*n(idx_OHj)
      pdj(72) =  &
          +k(363)*n(idx_CH2j)  &
          -k(40)*n(idx_CH3j)  &
          -k(384)*n(idx_CH3j)
      pdj(73) =  &
          +k(384)*n(idx_CH3j)
      pdj(74) =  &
          -k(142)*n(idx_Nj)  &
          -k(621)*n(idx_Nj)
      pdj(75) =  &
          +k(419)*n(idx_CNj)  &
          -k(539)*n(idx_HCNj)  &
          -k(540)*n(idx_HCNj)
      pdj(76) =  &
          +k(621)*n(idx_Nj)  &
          -k(654)*n(idx_NHj)
      pdj(77) =  &
          -k(451)*n(idx_H2j)  &
          -k(93)*n(idx_H2j)  &
          +k(437)*n(idx_Hj)
      pdj(78) =  &
          -k(591)*n(idx_HEj)  &
          -k(590)*n(idx_HEj)  &
          -k(592)*n(idx_HEj)
      pdj(79) =  &
          -k(553)*n(idx_HNOj)
      pdj(81) =  &
          +k(451)*n(idx_H2j)  &
          -k(514)*n(idx_H3j)
      pdj(82) =  &
          +k(552)*n(idx_H2COj)
      pdj(83) =  &
          +k(487)*n(idx_H2Oj)
      pdj(84) =  &
          +k(540)*n(idx_HCNj)
      pdj(86) =  &
          +k(591)*n(idx_HEj)
      pdj(87) =  &
          +k(629)*n(idx_N2j)  &
          -k(554)*n(idx_N2Hj)
      pdj(88) =  &
          -k(556)*n(idx_O2Hj)  &
          +k(555)*n(idx_O2j)
    elseif(j==15) then
      pdj(1) =  &
          +k(240)  &
          +k(1018)
      pdj(2) =  &
          +k(27)*n(idx_CHj)
      pdj(6) =  &
          +k(517)*n(idx_H3j)
      pdj(7) =  &
          +k(19)*n(idx_Cj)
      pdj(8) =  &
          +k(517)*n(idx_H3j)  &
          +k(77)*n(idx_Hj)
      pdj(9) =  &
          +k(104)*n(idx_H2Oj)
      pdj(11) =  &
          +k(134)*n(idx_O2j)
      pdj(13) =  &
          +k(130)*n(idx_H2COj)
      pdj(14) =  &
          +k(131)*n(idx_HCOj)
      pdj(15) =  &
          -k(130)*n(idx_H2COj)  &
          -k(27)*n(idx_CHj)  &
          -k(517)*n(idx_H3j)  &
          -k(41)*n(idx_CH3j)  &
          -k(131)*n(idx_HCOj)  &
          -k(1018)  &
          -k(19)*n(idx_Cj)  &
          -k(134)*n(idx_O2j)  &
          -k(104)*n(idx_H2Oj)  &
          -k(77)*n(idx_Hj)  &
          -k(240)  &
          -k(143)*n(idx_Nj)  &
          -k(132)*n(idx_N2j)  &
          -k(1120)  &
          -k(133)*n(idx_NOj)  &
          -k(170)*n(idx_NH3j)
      pdj(16) =  &
          +k(170)*n(idx_NH3j)
      pdj(17) =  &
          +k(133)*n(idx_NOj)
      pdj(20) =  &
          +k(132)*n(idx_N2j)
      pdj(22) =  &
          +k(41)*n(idx_CH3j)
      pdj(24) =  &
          +k(143)*n(idx_Nj)
      pdj(51) =  &
          +k(1120)
      pdj(54) =  &
          -k(131)*n(idx_HCOj)
      pdj(55) =  &
          -k(77)*n(idx_Hj)
      pdj(57) =  &
          -k(19)*n(idx_Cj)
      pdj(59) =  &
          -k(27)*n(idx_CHj)
      pdj(60) =  &
          -k(130)*n(idx_H2COj)
      pdj(61) =  &
          +k(170)*n(idx_NH3j)  &
          +k(131)*n(idx_HCOj)  &
          +k(517)*n(idx_H3j)  &
          +k(134)*n(idx_O2j)  &
          +k(104)*n(idx_H2Oj)  &
          +k(77)*n(idx_Hj)  &
          +k(143)*n(idx_Nj)  &
          +k(1018)  &
          +k(132)*n(idx_N2j)  &
          +k(27)*n(idx_CHj)  &
          +k(19)*n(idx_Cj)  &
          +k(133)*n(idx_NOj)  &
          +k(240)  &
          +k(41)*n(idx_CH3j)  &
          +k(130)*n(idx_H2COj)
      pdj(62) =  &
          -k(170)*n(idx_NH3j)
      pdj(63) =  &
          -k(133)*n(idx_NOj)
      pdj(66) =  &
          -k(132)*n(idx_N2j)
      pdj(67) =  &
          -k(134)*n(idx_O2j)
      pdj(68) =  &
          -k(104)*n(idx_H2Oj)
      pdj(72) =  &
          -k(41)*n(idx_CH3j)
      pdj(74) =  &
          -k(143)*n(idx_Nj)
      pdj(81) =  &
          -k(517)*n(idx_H3j)
    elseif(j==16) then
      pdj(1) =  &
          +k(246)  &
          +k(1024)
      pdj(2) =  &
          +k(28)*n(idx_CHj)
      pdj(3) =  &
          -k(954)*n(idx_O)  &
          +k(191)*n(idx_Oj)
      pdj(5) =  &
          +k(175)*n(idx_HCNj)  &
          +k(913)*n(idx_CN)
      pdj(6) =  &
          +k(323)*n(idx_Cj)  &
          +k(1025)  &
          +k(601)*n(idx_HEj)  &
          +k(865)*n(idx_H)  &
          +k(622)*n(idx_Nj)  &
          +k(95)*n(idx_H2j)  &
          +k(247)
      pdj(7) =  &
          +k(20)*n(idx_Cj)
      pdj(8) =  &
          +k(79)*n(idx_Hj)  &
          +k(602)*n(idx_HEj)  &
          +k(1023)  &
          +k(245)  &
          -k(865)*n(idx_H)
      pdj(9) =  &
          +k(174)*n(idx_H2Oj)  &
          +k(969)*n(idx_OH)
      pdj(10) =  &
          -k(969)*n(idx_OH)  &
          +k(200)*n(idx_OHj)  &
          +k(954)*n(idx_O)
      pdj(11) =  &
          +k(177)*n(idx_O2j)
      pdj(13) =  &
          +k(173)*n(idx_H2COj)
      pdj(16) =  &
          -k(623)*n(idx_Nj)  &
          -k(602)*n(idx_HEj)  &
          -k(954)*n(idx_O)  &
          -k(28)*n(idx_CHj)  &
          -k(245)  &
          -k(174)*n(idx_H2Oj)  &
          -k(79)*n(idx_Hj)  &
          -k(1087)  &
          -k(172)*n(idx_COj)  &
          -k(20)*n(idx_Cj)  &
          -k(917)*n(idx_NH)  &
          -k(95)*n(idx_H2j)  &
          -k(200)*n(idx_OHj)  &
          -k(173)*n(idx_H2COj)  &
          -k(247)  &
          -k(790)*n(idx_CH3)  &
          -k(323)*n(idx_Cj)  &
          -k(1024)  &
          -k(44)*n(idx_CH4j)  &
          -k(688)*n(idx_HCNj)  &
          -k(687)*n(idx_COj)  &
          -k(913)*n(idx_CN)  &
          -k(1023)  &
          -k(177)*n(idx_O2j)  &
          -k(1025)  &
          -k(157)*n(idx_NHj)  &
          -k(128)*n(idx_HEj)  &
          -k(191)*n(idx_Oj)  &
          -k(622)*n(idx_Nj)  &
          -k(969)*n(idx_OH)  &
          -k(865)*n(idx_H)  &
          -k(161)*n(idx_NH2j)  &
          -k(176)*n(idx_N2j)  &
          -k(145)*n(idx_Nj)  &
          -k(246)  &
          -k(601)*n(idx_HEj)  &
          -k(175)*n(idx_HCNj)
      pdj(18) =  &
          -k(913)*n(idx_CN)
      pdj(19) =  &
          +k(172)*n(idx_COj)
      pdj(20) =  &
          +k(176)*n(idx_N2j)
      pdj(21) =  &
          +2.d0*k(917)*n(idx_NH)  &
          +k(245)  &
          +k(969)*n(idx_OH)  &
          +k(954)*n(idx_O)  &
          +k(688)*n(idx_HCNj)  &
          +k(687)*n(idx_COj)  &
          +k(161)*n(idx_NH2j)  &
          +k(865)*n(idx_H)  &
          +k(1023)  &
          +k(790)*n(idx_CH3)  &
          +k(913)*n(idx_CN)
      pdj(22) =  &
          -k(790)*n(idx_CH3)
      pdj(23) =  &
          +k(44)*n(idx_CH4j)  &
          +k(790)*n(idx_CH3)
      pdj(24) =  &
          +k(145)*n(idx_Nj)
      pdj(25) =  &
          +k(247)  &
          -k(917)*n(idx_NH)  &
          +k(1025)  &
          +k(623)*n(idx_Nj)  &
          +k(157)*n(idx_NHj)
      pdj(26) =  &
          +k(602)*n(idx_HEj)  &
          +k(601)*n(idx_HEj)  &
          +k(128)*n(idx_HEj)
      pdj(45) =  &
          +k(1087)
      pdj(54) =  &
          +k(687)*n(idx_COj)
      pdj(55) =  &
          -k(79)*n(idx_Hj)
      pdj(57) =  &
          -k(20)*n(idx_Cj)  &
          -k(323)*n(idx_Cj)
      pdj(59) =  &
          -k(28)*n(idx_CHj)
      pdj(60) =  &
          -k(173)*n(idx_H2COj)
      pdj(62) =  &
          +k(246)  &
          +k(44)*n(idx_CH4j)  &
          +k(191)*n(idx_Oj)  &
          +k(28)*n(idx_CHj)  &
          +k(176)*n(idx_N2j)  &
          +k(157)*n(idx_NHj)  &
          +k(200)*n(idx_OHj)  &
          +k(1024)  &
          +k(175)*n(idx_HCNj)  &
          +k(172)*n(idx_COj)  &
          +k(145)*n(idx_Nj)  &
          +k(177)*n(idx_O2j)  &
          +k(161)*n(idx_NH2j)  &
          +k(128)*n(idx_HEj)  &
          +k(79)*n(idx_Hj)  &
          +k(20)*n(idx_Cj)  &
          +k(174)*n(idx_H2Oj)  &
          +k(95)*n(idx_H2j)  &
          +k(173)*n(idx_H2COj)
      pdj(65) =  &
          -k(687)*n(idx_COj)  &
          -k(172)*n(idx_COj)
      pdj(66) =  &
          -k(176)*n(idx_N2j)
      pdj(67) =  &
          -k(177)*n(idx_O2j)
      pdj(68) =  &
          -k(174)*n(idx_H2Oj)
      pdj(69) =  &
          +k(602)*n(idx_HEj)  &
          -k(161)*n(idx_NH2j)  &
          +k(623)*n(idx_Nj)
      pdj(70) =  &
          -k(191)*n(idx_Oj)
      pdj(71) =  &
          -k(200)*n(idx_OHj)
      pdj(73) =  &
          -k(44)*n(idx_CH4j)
      pdj(74) =  &
          -k(623)*n(idx_Nj)  &
          -k(145)*n(idx_Nj)  &
          -k(622)*n(idx_Nj)
      pdj(75) =  &
          +k(323)*n(idx_Cj)  &
          -k(175)*n(idx_HCNj)  &
          -k(688)*n(idx_HCNj)
      pdj(76) =  &
          +k(601)*n(idx_HEj)  &
          -k(157)*n(idx_NHj)
      pdj(77) =  &
          -k(95)*n(idx_H2j)
      pdj(78) =  &
          -k(128)*n(idx_HEj)  &
          -k(601)*n(idx_HEj)  &
          -k(602)*n(idx_HEj)
      pdj(84) =  &
          +k(688)*n(idx_HCNj)
      pdj(87) =  &
          +k(622)*n(idx_Nj)
    elseif(j==17) then
      pdj(1) =  &
          +k(251)  &
          +k(1029)
      pdj(2) =  &
          -k(812)*n(idx_CH)  &
          -k(813)*n(idx_CH)  &
          -k(814)*n(idx_CH)  &
          +k(29)*n(idx_CHj)
      pdj(3) =  &
          +k(903)*n(idx_N)  &
          +k(625)*n(idx_Nj)  &
          +k(754)*n(idx_C)  &
          +k(868)*n(idx_H)  &
          +k(605)*n(idx_HEj)  &
          -k(956)*n(idx_O)  &
          +k(932)*n(idx_O2)  &
          +k(1030)  &
          +k(252)  &
          +k(659)*n(idx_NHj)  &
          +k(736)*n(idx_OHj)  &
          +k(922)*n(idx_NH)  &
          +k(812)*n(idx_CH)
      pdj(5) =  &
          +k(769)*n(idx_CH2)  &
          +k(792)*n(idx_CH3)  &
          +k(116)*n(idx_HCNj)  &
          +k(812)*n(idx_CH)
      pdj(6) =  &
          +k(97)*n(idx_H2j)  &
          +k(522)*n(idx_H3j)
      pdj(7) =  &
          -k(754)*n(idx_C)  &
          -k(755)*n(idx_C)  &
          +k(21)*n(idx_Cj)
      pdj(8) =  &
          -k(868)*n(idx_H)  &
          +k(970)*n(idx_OH)  &
          +k(814)*n(idx_CH)  &
          +k(81)*n(idx_Hj)  &
          +k(910)*n(idx_NH2)  &
          +k(456)*n(idx_H2j)  &
          -k(869)*n(idx_H)  &
          +k(922)*n(idx_NH)  &
          +k(770)*n(idx_CH2)
      pdj(9) =  &
          +k(105)*n(idx_H2Oj)  &
          +k(792)*n(idx_CH3)  &
          +k(909)*n(idx_NH2)
      pdj(10) =  &
          +k(923)*n(idx_NH)  &
          +k(201)*n(idx_OHj)  &
          -k(970)*n(idx_OH)  &
          +k(869)*n(idx_H)  &
          +k(910)*n(idx_NH2)  &
          +k(769)*n(idx_CH2)
      pdj(11) =  &
          +k(702)*n(idx_O2Hj)  &
          +k(184)*n(idx_O2j)  &
          +2.d0*k(931)*n(idx_NO)  &
          -k(932)*n(idx_O2)  &
          +k(956)*n(idx_O)
      pdj(12) =  &
          +k(30)*n(idx_CH2j)  &
          -k(770)*n(idx_CH2)  &
          -k(768)*n(idx_CH2)  &
          -k(769)*n(idx_CH2)
      pdj(13) =  &
          +k(768)*n(idx_CH2)  &
          +k(182)*n(idx_H2COj)
      pdj(14) =  &
          +k(813)*n(idx_CH)  &
          -k(881)*n(idx_HCO)
      pdj(16) =  &
          +k(171)*n(idx_NH3j)
      pdj(17) =  &
          -k(868)*n(idx_H)  &
          -k(754)*n(idx_C)  &
          -k(755)*n(idx_C)  &
          -k(881)*n(idx_HCO)  &
          -4.d0*k(931)*n(idx_NO)  &
          -k(922)*n(idx_NH)  &
          -k(81)*n(idx_Hj)  &
          -k(158)*n(idx_NHj)  &
          -k(923)*n(idx_NH)  &
          -k(659)*n(idx_NHj)  &
          -k(933)*n(idx_OCN)  &
          -k(105)*n(idx_H2Oj)  &
          -k(522)*n(idx_H3j)  &
          -k(116)*n(idx_HCNj)  &
          -k(956)*n(idx_O)  &
          -k(768)*n(idx_CH2)  &
          -k(770)*n(idx_CH2)  &
          -k(970)*n(idx_OH)  &
          -k(814)*n(idx_CH)  &
          -k(909)*n(idx_NH2)  &
          -k(605)*n(idx_HEj)  &
          -k(812)*n(idx_CH)  &
          -k(792)*n(idx_CH3)  &
          -k(869)*n(idx_H)  &
          -k(1075)  &
          -k(251)  &
          -k(769)*n(idx_CH2)  &
          -k(42)*n(idx_CH3j)  &
          -k(66)*n(idx_COj)  &
          -k(184)*n(idx_O2j)  &
          -k(183)*n(idx_HNOj)  &
          -k(1030)  &
          -k(201)*n(idx_OHj)  &
          -k(29)*n(idx_CHj)  &
          -k(456)*n(idx_H2j)  &
          -k(21)*n(idx_Cj)  &
          -k(252)  &
          -k(828)*n(idx_CN)  &
          -k(97)*n(idx_H2j)  &
          -k(61)*n(idx_CNj)  &
          -k(932)*n(idx_O2)  &
          -k(152)*n(idx_N2j)  &
          -k(813)*n(idx_CH)  &
          -k(147)*n(idx_Nj)  &
          -k(910)*n(idx_NH2)  &
          -k(625)*n(idx_Nj)  &
          -k(604)*n(idx_HEj)  &
          -k(171)*n(idx_NH3j)  &
          -k(1029)  &
          -k(182)*n(idx_H2COj)  &
          -k(736)*n(idx_OHj)  &
          -k(162)*n(idx_NH2j)  &
          -k(30)*n(idx_CH2j)  &
          -k(702)*n(idx_O2Hj)  &
          -k(829)*n(idx_CN)  &
          -k(903)*n(idx_N)
      pdj(18) =  &
          -k(829)*n(idx_CN)  &
          +k(754)*n(idx_C)  &
          -k(828)*n(idx_CN)  &
          +k(61)*n(idx_CNj)
      pdj(19) =  &
          +k(755)*n(idx_C)  &
          +k(66)*n(idx_COj)  &
          +k(828)*n(idx_CN)  &
          +k(881)*n(idx_HCO)
      pdj(20) =  &
          +k(903)*n(idx_N)  &
          +k(923)*n(idx_NH)  &
          +k(933)*n(idx_OCN)  &
          +k(910)*n(idx_NH2)  &
          +2.d0*k(931)*n(idx_NO)  &
          +k(152)*n(idx_N2j)  &
          +k(828)*n(idx_CN)  &
          +k(909)*n(idx_NH2)  &
          +k(922)*n(idx_NH)
      pdj(21) =  &
          -k(909)*n(idx_NH2)  &
          -k(910)*n(idx_NH2)  &
          +k(162)*n(idx_NH2j)
      pdj(22) =  &
          +k(42)*n(idx_CH3j)  &
          -k(792)*n(idx_CH3)
      pdj(24) =  &
          +k(755)*n(idx_C)  &
          +k(147)*n(idx_Nj)  &
          +k(956)*n(idx_O)  &
          +k(813)*n(idx_CH)  &
          +k(869)*n(idx_H)  &
          +k(768)*n(idx_CH2)  &
          +k(604)*n(idx_HEj)  &
          +k(252)  &
          +k(1030)  &
          -k(903)*n(idx_N)  &
          +k(829)*n(idx_CN)
      pdj(25) =  &
          -k(922)*n(idx_NH)  &
          +k(158)*n(idx_NHj)  &
          +k(868)*n(idx_H)  &
          -k(923)*n(idx_NH)
      pdj(26) =  &
          +k(605)*n(idx_HEj)  &
          +k(604)*n(idx_HEj)
      pdj(27) =  &
          +k(183)*n(idx_HNOj)  &
          +k(881)*n(idx_HCO)
      pdj(29) =  &
          +k(933)*n(idx_OCN)
      pdj(31) =  &
          +k(770)*n(idx_CH2)
      pdj(32) =  &
          +k(932)*n(idx_O2)  &
          +k(970)*n(idx_OH)
      pdj(34) =  &
          -k(933)*n(idx_OCN)  &
          +k(814)*n(idx_CH)  &
          +k(829)*n(idx_CN)
      pdj(41) =  &
          +k(1075)
      pdj(55) =  &
          -k(81)*n(idx_Hj)
      pdj(57) =  &
          -k(21)*n(idx_Cj)
      pdj(58) =  &
          -k(30)*n(idx_CH2j)
      pdj(59) =  &
          -k(29)*n(idx_CHj)
      pdj(60) =  &
          -k(182)*n(idx_H2COj)
      pdj(62) =  &
          -k(171)*n(idx_NH3j)
      pdj(63) =  &
          +k(183)*n(idx_HNOj)  &
          +k(66)*n(idx_COj)  &
          +k(201)*n(idx_OHj)  &
          +k(171)*n(idx_NH3j)  &
          +k(21)*n(idx_Cj)  &
          +k(97)*n(idx_H2j)  &
          +k(30)*n(idx_CH2j)  &
          +k(251)  &
          +k(105)*n(idx_H2Oj)  &
          +k(1029)  &
          +k(158)*n(idx_NHj)  &
          +k(182)*n(idx_H2COj)  &
          +k(162)*n(idx_NH2j)  &
          +k(147)*n(idx_Nj)  &
          +k(81)*n(idx_Hj)  &
          +k(42)*n(idx_CH3j)  &
          +k(152)*n(idx_N2j)  &
          +k(116)*n(idx_HCNj)  &
          +k(184)*n(idx_O2j)  &
          +k(29)*n(idx_CHj)  &
          +k(61)*n(idx_CNj)
      pdj(64) =  &
          -k(61)*n(idx_CNj)
      pdj(65) =  &
          -k(66)*n(idx_COj)
      pdj(66) =  &
          +k(625)*n(idx_Nj)  &
          -k(152)*n(idx_N2j)
      pdj(67) =  &
          -k(184)*n(idx_O2j)
      pdj(68) =  &
          -k(105)*n(idx_H2Oj)
      pdj(69) =  &
          -k(162)*n(idx_NH2j)
      pdj(70) =  &
          +k(604)*n(idx_HEj)
      pdj(71) =  &
          -k(201)*n(idx_OHj)  &
          -k(736)*n(idx_OHj)
      pdj(72) =  &
          -k(42)*n(idx_CH3j)
      pdj(74) =  &
          -k(625)*n(idx_Nj)  &
          +k(605)*n(idx_HEj)  &
          -k(147)*n(idx_Nj)
      pdj(75) =  &
          -k(116)*n(idx_HCNj)
      pdj(76) =  &
          -k(659)*n(idx_NHj)  &
          -k(158)*n(idx_NHj)
      pdj(77) =  &
          -k(456)*n(idx_H2j)  &
          -k(97)*n(idx_H2j)
      pdj(78) =  &
          -k(605)*n(idx_HEj)  &
          -k(604)*n(idx_HEj)
      pdj(79) =  &
          +k(702)*n(idx_O2Hj)  &
          +k(736)*n(idx_OHj)  &
          +k(522)*n(idx_H3j)  &
          -k(183)*n(idx_HNOj)  &
          +k(456)*n(idx_H2j)
      pdj(81) =  &
          -k(522)*n(idx_H3j)
      pdj(87) =  &
          +k(659)*n(idx_NHj)
      pdj(88) =  &
          -k(702)*n(idx_O2Hj)
    elseif(j==18) then
      pdj(2) =  &
          +k(762)*n(idx_CH2)
      pdj(3) =  &
          +k(726)*n(idx_OHj)  &
          -k(938)*n(idx_O)  &
          -k(937)*n(idx_O)  &
          +k(831)*n(idx_O2)  &
          +k(961)*n(idx_OH)
      pdj(5) =  &
          +k(826)*n(idx_HNO)  &
          +k(824)*n(idx_H2CO)  &
          +k(802)*n(idx_CH4)  &
          +k(840)*n(idx_H2)  &
          +k(961)*n(idx_OH)  &
          +k(915)*n(idx_NH)  &
          +k(762)*n(idx_CH2)  &
          +k(825)*n(idx_HCO)  &
          +k(784)*n(idx_CH3)  &
          +k(913)*n(idx_NH3)
      pdj(6) =  &
          +k(88)*n(idx_H2j)  &
          -k(840)*n(idx_H2)  &
          +k(507)*n(idx_H3j)
      pdj(7) =  &
          +k(996)  &
          +k(938)*n(idx_O)  &
          +k(706)*n(idx_Oj)  &
          +k(892)*n(idx_N)  &
          +k(574)*n(idx_HEj)  &
          +k(226)
      pdj(8) =  &
          +k(840)*n(idx_H2)  &
          +k(445)*n(idx_H2j)  &
          +k(962)*n(idx_OH)
      pdj(10) =  &
          -k(961)*n(idx_OH)  &
          -k(962)*n(idx_OH)
      pdj(11) =  &
          -k(831)*n(idx_O2)  &
          +k(422)*n(idx_O2Hj)  &
          -k(830)*n(idx_O2)
      pdj(12) =  &
          +k(784)*n(idx_CH3)  &
          -k(762)*n(idx_CH2)
      pdj(13) =  &
          -k(824)*n(idx_H2CO)
      pdj(14) =  &
          +k(824)*n(idx_H2CO)  &
          -k(825)*n(idx_HCO)
      pdj(16) =  &
          -k(913)*n(idx_NH3)
      pdj(17) =  &
          +k(826)*n(idx_HNO)  &
          +k(421)*n(idx_HNOj)  &
          +k(938)*n(idx_O)  &
          -k(829)*n(idx_NO)  &
          +k(830)*n(idx_O2)  &
          -k(828)*n(idx_NO)  &
          +k(827)*n(idx_NO2)
      pdj(18) =  &
          -k(574)*n(idx_HEj)  &
          -k(913)*n(idx_NH3)  &
          -k(445)*n(idx_H2j)  &
          -k(726)*n(idx_OHj)  &
          -k(762)*n(idx_CH2)  &
          -k(962)*n(idx_OH)  &
          -k(802)*n(idx_CH4)  &
          -k(915)*n(idx_NH)  &
          -k(824)*n(idx_H2CO)  &
          -k(996)  &
          -k(226)  &
          -k(784)*n(idx_CH3)  &
          -k(706)*n(idx_Oj)  &
          -k(831)*n(idx_O2)  &
          -k(88)*n(idx_H2j)  &
          -k(507)*n(idx_H3j)  &
          -k(938)*n(idx_O)  &
          -k(828)*n(idx_NO)  &
          -k(892)*n(idx_N)  &
          -k(575)*n(idx_HEj)  &
          -k(421)*n(idx_HNOj)  &
          -k(830)*n(idx_O2)  &
          -k(137)*n(idx_Nj)  &
          -k(961)*n(idx_OH)  &
          -k(422)*n(idx_O2Hj)  &
          -k(1083)  &
          -k(829)*n(idx_NO)  &
          -k(937)*n(idx_O)  &
          -k(63)*n(idx_N2j)  &
          -k(642)*n(idx_NHj)  &
          -k(827)*n(idx_NO2)  &
          -k(826)*n(idx_HNO)  &
          -k(825)*n(idx_HCO)  &
          -k(840)*n(idx_H2)
      pdj(19) =  &
          +k(830)*n(idx_O2)  &
          +k(937)*n(idx_O)  &
          +k(828)*n(idx_NO)  &
          +k(825)*n(idx_HCO)
      pdj(20) =  &
          +k(63)*n(idx_N2j)  &
          +k(828)*n(idx_NO)  &
          +k(892)*n(idx_N)
      pdj(21) =  &
          +k(913)*n(idx_NH3)
      pdj(22) =  &
          +k(802)*n(idx_CH4)  &
          -k(784)*n(idx_CH3)
      pdj(23) =  &
          -k(802)*n(idx_CH4)
      pdj(24) =  &
          +k(937)*n(idx_O)  &
          +k(996)  &
          +k(915)*n(idx_NH)  &
          +k(642)*n(idx_NHj)  &
          +k(575)*n(idx_HEj)  &
          -k(892)*n(idx_N)  &
          +k(137)*n(idx_Nj)  &
          +k(829)*n(idx_NO)  &
          +k(226)
      pdj(25) =  &
          -k(915)*n(idx_NH)
      pdj(26) =  &
          +k(575)*n(idx_HEj)  &
          +k(574)*n(idx_HEj)
      pdj(27) =  &
          -k(826)*n(idx_HNO)
      pdj(32) =  &
          -k(827)*n(idx_NO2)
      pdj(34) =  &
          +k(829)*n(idx_NO)  &
          +k(831)*n(idx_O2)  &
          +k(962)*n(idx_OH)  &
          +k(827)*n(idx_NO2)
      pdj(44) =  &
          +k(1083)
      pdj(57) =  &
          +k(575)*n(idx_HEj)
      pdj(63) =  &
          +k(706)*n(idx_Oj)
      pdj(64) =  &
          +k(88)*n(idx_H2j)  &
          +k(63)*n(idx_N2j)  &
          +k(137)*n(idx_Nj)
      pdj(66) =  &
          -k(63)*n(idx_N2j)
      pdj(70) =  &
          -k(706)*n(idx_Oj)
      pdj(71) =  &
          -k(726)*n(idx_OHj)
      pdj(74) =  &
          -k(137)*n(idx_Nj)  &
          +k(574)*n(idx_HEj)
      pdj(75) =  &
          +k(726)*n(idx_OHj)  &
          +k(421)*n(idx_HNOj)  &
          +k(507)*n(idx_H3j)  &
          +k(422)*n(idx_O2Hj)  &
          +k(642)*n(idx_NHj)  &
          +k(445)*n(idx_H2j)
      pdj(76) =  &
          -k(642)*n(idx_NHj)
      pdj(77) =  &
          -k(88)*n(idx_H2j)  &
          -k(445)*n(idx_H2j)
      pdj(78) =  &
          -k(575)*n(idx_HEj)  &
          -k(574)*n(idx_HEj)
      pdj(79) =  &
          -k(421)*n(idx_HNOj)
      pdj(81) =  &
          -k(507)*n(idx_H3j)
      pdj(88) =  &
          -k(422)*n(idx_O2Hj)
    elseif(j==19) then
      pdj(1) =  &
          +k(207)
      pdj(3) =  &
          +k(580)*n(idx_HEj)  &
          +k(834)*n(idx_O2)  &
          +k(999)  &
          +k(728)*n(idx_OHj)  &
          +k(228)  &
          +k(186)*n(idx_Oj)
      pdj(6) =  &
          +k(89)*n(idx_H2j)  &
          +k(510)*n(idx_H3j)  &
          +k(509)*n(idx_H3j)
      pdj(7) =  &
          +k(853)*n(idx_H)  &
          +k(618)*n(idx_Nj)  &
          +k(228)  &
          +k(999)
      pdj(8) =  &
          +k(447)*n(idx_H2j)  &
          +k(963)*n(idx_OH)  &
          -k(853)*n(idx_H)
      pdj(10) =  &
          +k(835)*n(idx_O2H)  &
          -k(963)*n(idx_OH)  &
          +k(483)*n(idx_H2Oj)  &
          +k(853)*n(idx_H)
      pdj(11) =  &
          +k(427)*n(idx_O2Hj)  &
          -k(834)*n(idx_O2)
      pdj(17) =  &
          +k(833)*n(idx_NO2)  &
          +k(425)*n(idx_HNOj)
      pdj(18) =  &
          +k(536)*n(idx_HCNj)  &
          +k(57)*n(idx_CNj)
      pdj(19) =  &
          -k(832)*n(idx_HNO)  &
          -k(89)*n(idx_H2j)  &
          -k(833)*n(idx_NO2)  &
          -k(447)*n(idx_H2j)  &
          -k(68)*n(idx_N2j)  &
          -k(509)*n(idx_H3j)  &
          -k(424)*n(idx_HCO2j)  &
          -k(999)  &
          -k(728)*n(idx_OHj)  &
          -k(536)*n(idx_HCNj)  &
          -k(483)*n(idx_H2Oj)  &
          -k(427)*n(idx_O2Hj)  &
          -k(186)*n(idx_Oj)  &
          -k(580)*n(idx_HEj)  &
          -k(834)*n(idx_O2)  &
          -k(1124)  &
          -k(228)  &
          -k(390)*n(idx_CH4j)  &
          -k(1071)  &
          -k(510)*n(idx_H3j)  &
          -k(835)*n(idx_O2H)  &
          -k(426)*n(idx_N2Hj)  &
          -k(646)*n(idx_NHj)  &
          -k(1068)  &
          -k(138)*n(idx_Nj)  &
          -k(963)*n(idx_OH)  &
          -k(425)*n(idx_HNOj)  &
          -k(853)*n(idx_H)  &
          -k(207)  &
          -k(618)*n(idx_Nj)  &
          -k(57)*n(idx_CNj)
      pdj(20) =  &
          +k(68)*n(idx_N2j)  &
          +k(426)*n(idx_N2Hj)
      pdj(22) =  &
          +k(390)*n(idx_CH4j)
      pdj(24) =  &
          +k(646)*n(idx_NHj)  &
          +k(138)*n(idx_Nj)
      pdj(25) =  &
          +k(832)*n(idx_HNO)
      pdj(26) =  &
          +k(580)*n(idx_HEj)
      pdj(27) =  &
          -k(832)*n(idx_HNO)
      pdj(29) =  &
          +k(835)*n(idx_O2H)  &
          +k(834)*n(idx_O2)  &
          +k(832)*n(idx_HNO)  &
          +k(424)*n(idx_HCO2j)  &
          +k(833)*n(idx_NO2)  &
          +k(963)*n(idx_OH)
      pdj(32) =  &
          -k(833)*n(idx_NO2)
      pdj(33) =  &
          -k(835)*n(idx_O2H)
      pdj(35) =  &
          +k(1124)
      pdj(37) =  &
          +k(1068)
      pdj(39) =  &
          +k(1071)
      pdj(54) =  &
          +k(646)*n(idx_NHj)  &
          +k(425)*n(idx_HNOj)  &
          +k(390)*n(idx_CH4j)  &
          +k(483)*n(idx_H2Oj)  &
          +k(427)*n(idx_O2Hj)  &
          +k(424)*n(idx_HCO2j)  &
          +k(728)*n(idx_OHj)  &
          +k(536)*n(idx_HCNj)  &
          +k(426)*n(idx_N2Hj)  &
          +k(447)*n(idx_H2j)  &
          +k(509)*n(idx_H3j)
      pdj(56) =  &
          +k(510)*n(idx_H3j)
      pdj(57) =  &
          +k(580)*n(idx_HEj)
      pdj(63) =  &
          +k(618)*n(idx_Nj)
      pdj(64) =  &
          -k(57)*n(idx_CNj)
      pdj(65) =  &
          +k(207)  &
          +k(68)*n(idx_N2j)  &
          +k(57)*n(idx_CNj)  &
          +k(186)*n(idx_Oj)  &
          +k(89)*n(idx_H2j)  &
          +k(138)*n(idx_Nj)
      pdj(66) =  &
          -k(68)*n(idx_N2j)
      pdj(68) =  &
          -k(483)*n(idx_H2Oj)
      pdj(70) =  &
          -k(186)*n(idx_Oj)
      pdj(71) =  &
          -k(728)*n(idx_OHj)
      pdj(73) =  &
          -k(390)*n(idx_CH4j)
      pdj(74) =  &
          -k(618)*n(idx_Nj)  &
          -k(138)*n(idx_Nj)
      pdj(75) =  &
          -k(536)*n(idx_HCNj)
      pdj(76) =  &
          -k(646)*n(idx_NHj)
      pdj(77) =  &
          -k(89)*n(idx_H2j)  &
          -k(447)*n(idx_H2j)
      pdj(78) =  &
          -k(580)*n(idx_HEj)
      pdj(79) =  &
          -k(425)*n(idx_HNOj)
      pdj(81) =  &
          -k(509)*n(idx_H3j)  &
          -k(510)*n(idx_H3j)
      pdj(85) =  &
          -k(424)*n(idx_HCO2j)
      pdj(87) =  &
          -k(426)*n(idx_N2Hj)
      pdj(88) =  &
          -k(427)*n(idx_O2Hj)
    elseif(j==20) then
      pdj(2) =  &
          -k(809)*n(idx_CH)
      pdj(3) =  &
          -k(951)*n(idx_O)  &
          +k(735)*n(idx_OHj)
      pdj(5) =  &
          +k(766)*n(idx_CH2)  &
          +k(809)*n(idx_CH)
      pdj(6) =  &
          +k(518)*n(idx_H3j)
      pdj(7) =  &
          -k(748)*n(idx_C)
      pdj(8) =  &
          +k(453)*n(idx_H2j)
      pdj(11) =  &
          +k(631)*n(idx_O2Hj)
      pdj(12) =  &
          -k(766)*n(idx_CH2)
      pdj(17) =  &
          +k(951)*n(idx_O)  &
          +k(630)*n(idx_HNOj)
      pdj(18) =  &
          +k(748)*n(idx_C)
      pdj(20) =  &
          -k(951)*n(idx_O)  &
          -k(809)*n(idx_CH)  &
          -k(656)*n(idx_NHj)  &
          -k(735)*n(idx_OHj)  &
          -k(1082)  &
          -k(1019)  &
          -k(518)*n(idx_H3j)  &
          -k(766)*n(idx_CH2)  &
          -k(598)*n(idx_HEj)  &
          -k(127)*n(idx_HEj)  &
          -k(748)*n(idx_C)  &
          -k(453)*n(idx_H2j)  &
          -k(630)*n(idx_HNOj)  &
          -k(712)*n(idx_Oj)  &
          -k(241)  &
          -k(631)*n(idx_O2Hj)
      pdj(24) =  &
          +k(809)*n(idx_CH)  &
          +k(712)*n(idx_Oj)  &
          +k(951)*n(idx_O)  &
          +k(656)*n(idx_NHj)  &
          +k(598)*n(idx_HEj)  &
          +2.d0*k(1019)  &
          +2.d0*k(241)  &
          +k(748)*n(idx_C)
      pdj(25) =  &
          +k(766)*n(idx_CH2)
      pdj(26) =  &
          +k(127)*n(idx_HEj)  &
          +k(598)*n(idx_HEj)
      pdj(43) =  &
          +k(1082)
      pdj(63) =  &
          +k(712)*n(idx_Oj)
      pdj(66) =  &
          +k(127)*n(idx_HEj)
      pdj(70) =  &
          -k(712)*n(idx_Oj)
      pdj(71) =  &
          -k(735)*n(idx_OHj)
      pdj(74) =  &
          +k(598)*n(idx_HEj)
      pdj(76) =  &
          -k(656)*n(idx_NHj)
      pdj(77) =  &
          -k(453)*n(idx_H2j)
      pdj(78) =  &
          -k(598)*n(idx_HEj)  &
          -k(127)*n(idx_HEj)
      pdj(79) =  &
          -k(630)*n(idx_HNOj)
      pdj(81) =  &
          -k(518)*n(idx_H3j)
      pdj(87) =  &
          +k(735)*n(idx_OHj)  &
          +k(631)*n(idx_O2Hj)  &
          +k(630)*n(idx_HNOj)  &
          +k(656)*n(idx_NHj)  &
          +k(453)*n(idx_H2j)  &
          +k(518)*n(idx_H3j)
      pdj(88) =  &
          -k(631)*n(idx_O2Hj)
    elseif(j==21) then
      pdj(1) =  &
          +k(243)  &
          +k(1021)
      pdj(2) =  &
          +k(751)*n(idx_C)
      pdj(3) =  &
          +k(912)*n(idx_OH)  &
          -k(952)*n(idx_O)  &
          +k(190)*n(idx_Oj)  &
          +k(686)*n(idx_OHj)  &
          -k(953)*n(idx_O)
      pdj(4) =  &
          +k(750)*n(idx_C)  &
          +k(681)*n(idx_HCNHj)
      pdj(5) =  &
          +k(680)*n(idx_HCNHj)  &
          +k(749)*n(idx_C)
      pdj(6) =  &
          +k(94)*n(idx_H2j)  &
          -k(842)*n(idx_H2)  &
          +k(519)*n(idx_H3j)  &
          +k(353)*n(idx_CHj)  &
          +k(864)*n(idx_H)  &
          +k(599)*n(idx_HEj)
      pdj(7) =  &
          -k(749)*n(idx_C)  &
          -k(751)*n(idx_C)  &
          -k(750)*n(idx_C)
      pdj(8) =  &
          +k(750)*n(idx_C)  &
          +k(78)*n(idx_Hj)  &
          +k(749)*n(idx_C)  &
          +k(244)  &
          +k(600)*n(idx_HEj)  &
          +k(842)*n(idx_H2)  &
          +k(322)*n(idx_Cj)  &
          +k(910)*n(idx_NO)  &
          -k(864)*n(idx_H)  &
          +k(1022)  &
          +k(952)*n(idx_O)
      pdj(9) =  &
          +k(678)*n(idx_H3Oj)  &
          +k(909)*n(idx_NO)  &
          +k(165)*n(idx_H2Oj)  &
          +k(911)*n(idx_OH)
      pdj(10) =  &
          +k(168)*n(idx_OHj)  &
          -k(911)*n(idx_OH)  &
          -k(912)*n(idx_OH)  &
          +k(676)*n(idx_H2Oj)  &
          +k(953)*n(idx_O)  &
          +k(910)*n(idx_NO)
      pdj(11) =  &
          +k(685)*n(idx_O2Hj)  &
          +k(167)*n(idx_O2j)
      pdj(13) =  &
          +k(677)*n(idx_H3COj)
      pdj(14) =  &
          +k(675)*n(idx_H2COj)
      pdj(16) =  &
          +k(912)*n(idx_OH)  &
          +k(842)*n(idx_H2)  &
          +k(908)*n(idx_CH4)
      pdj(17) =  &
          +k(683)*n(idx_HNOj)  &
          -k(910)*n(idx_NO)  &
          -k(909)*n(idx_NO)
      pdj(18) =  &
          +k(163)*n(idx_CNj)  &
          +k(679)*n(idx_HCNj)
      pdj(19) =  &
          +k(682)*n(idx_HCOj)  &
          +k(164)*n(idx_COj)
      pdj(20) =  &
          +k(684)*n(idx_N2Hj)  &
          +k(910)*n(idx_NO)  &
          +k(909)*n(idx_NO)  &
          +k(166)*n(idx_N2j)
      pdj(21) =  &
          -k(789)*n(idx_CH3)  &
          -k(679)*n(idx_HCNj)  &
          -k(78)*n(idx_Hj)  &
          -k(144)*n(idx_Nj)  &
          -k(683)*n(idx_HNOj)  &
          -k(1022)  &
          -k(685)*n(idx_O2Hj)  &
          -k(842)*n(idx_H2)  &
          -k(519)*n(idx_H3j)  &
          -k(911)*n(idx_OH)  &
          -k(322)*n(idx_Cj)  &
          -k(1021)  &
          -k(681)*n(idx_HCNHj)  &
          -k(676)*n(idx_H2Oj)  &
          -k(909)*n(idx_NO)  &
          -k(353)*n(idx_CHj)  &
          -k(94)*n(idx_H2j)  &
          -k(680)*n(idx_HCNHj)  &
          -k(163)*n(idx_CNj)  &
          -k(671)*n(idx_NH2j)  &
          -k(600)*n(idx_HEj)  &
          -k(675)*n(idx_H2COj)  &
          -k(678)*n(idx_H3Oj)  &
          -k(243)  &
          -k(908)*n(idx_CH4)  &
          -k(953)*n(idx_O)  &
          -k(190)*n(idx_Oj)  &
          -k(164)*n(idx_COj)  &
          -k(165)*n(idx_H2Oj)  &
          -k(244)  &
          -k(751)*n(idx_C)  &
          -k(952)*n(idx_O)  &
          -k(674)*n(idx_COj)  &
          -k(682)*n(idx_HCOj)  &
          -k(166)*n(idx_N2j)  &
          -k(167)*n(idx_O2j)  &
          -k(684)*n(idx_N2Hj)  &
          -k(912)*n(idx_OH)  &
          -k(1109)  &
          -k(599)*n(idx_HEj)  &
          -k(910)*n(idx_NO)  &
          -k(677)*n(idx_H3COj)  &
          -k(749)*n(idx_C)  &
          -k(168)*n(idx_OHj)  &
          -k(657)*n(idx_NHj)  &
          -k(686)*n(idx_OHj)  &
          -k(864)*n(idx_H)  &
          -k(750)*n(idx_C)
      pdj(22) =  &
          +k(908)*n(idx_CH4)  &
          -k(789)*n(idx_CH3)
      pdj(23) =  &
          +k(789)*n(idx_CH3)  &
          -k(908)*n(idx_CH4)
      pdj(24) =  &
          +k(144)*n(idx_Nj)  &
          +k(657)*n(idx_NHj)
      pdj(25) =  &
          +k(671)*n(idx_NH2j)  &
          +k(244)  &
          +k(911)*n(idx_OH)  &
          +k(953)*n(idx_O)  &
          +k(789)*n(idx_CH3)  &
          +k(674)*n(idx_COj)  &
          +k(864)*n(idx_H)  &
          +k(751)*n(idx_C)  &
          +k(1022)
      pdj(26) =  &
          +k(599)*n(idx_HEj)  &
          +k(600)*n(idx_HEj)
      pdj(27) =  &
          +k(952)*n(idx_O)
      pdj(45) =  &
          +k(1109)
      pdj(54) =  &
          +k(674)*n(idx_COj)  &
          -k(682)*n(idx_HCOj)
      pdj(55) =  &
          -k(78)*n(idx_Hj)
      pdj(57) =  &
          -k(322)*n(idx_Cj)
      pdj(59) =  &
          -k(353)*n(idx_CHj)
      pdj(60) =  &
          -k(675)*n(idx_H2COj)
      pdj(62) =  &
          +k(675)*n(idx_H2COj)  &
          +k(671)*n(idx_NH2j)  &
          +k(677)*n(idx_H3COj)  &
          +k(681)*n(idx_HCNHj)  &
          +k(682)*n(idx_HCOj)  &
          +k(684)*n(idx_N2Hj)  &
          +k(678)*n(idx_H3Oj)  &
          +k(683)*n(idx_HNOj)  &
          +k(519)*n(idx_H3j)  &
          +k(676)*n(idx_H2Oj)  &
          +k(680)*n(idx_HCNHj)  &
          +k(685)*n(idx_O2Hj)  &
          +k(657)*n(idx_NHj)  &
          +k(679)*n(idx_HCNj)  &
          +k(686)*n(idx_OHj)
      pdj(64) =  &
          -k(163)*n(idx_CNj)
      pdj(65) =  &
          -k(674)*n(idx_COj)  &
          -k(164)*n(idx_COj)
      pdj(66) =  &
          -k(166)*n(idx_N2j)
      pdj(67) =  &
          -k(167)*n(idx_O2j)
      pdj(68) =  &
          -k(165)*n(idx_H2Oj)  &
          -k(676)*n(idx_H2Oj)
      pdj(69) =  &
          +k(144)*n(idx_Nj)  &
          +k(168)*n(idx_OHj)  &
          +k(167)*n(idx_O2j)  &
          +k(94)*n(idx_H2j)  &
          +k(78)*n(idx_Hj)  &
          +k(190)*n(idx_Oj)  &
          +k(164)*n(idx_COj)  &
          +k(163)*n(idx_CNj)  &
          +k(165)*n(idx_H2Oj)  &
          -k(671)*n(idx_NH2j)  &
          +k(1021)  &
          +k(166)*n(idx_N2j)  &
          +k(243)
      pdj(70) =  &
          -k(190)*n(idx_Oj)
      pdj(71) =  &
          -k(168)*n(idx_OHj)  &
          -k(686)*n(idx_OHj)
      pdj(74) =  &
          -k(144)*n(idx_Nj)  &
          +k(599)*n(idx_HEj)
      pdj(75) =  &
          -k(679)*n(idx_HCNj)  &
          +k(353)*n(idx_CHj)  &
          +k(322)*n(idx_Cj)
      pdj(76) =  &
          +k(600)*n(idx_HEj)  &
          -k(657)*n(idx_NHj)
      pdj(77) =  &
          -k(94)*n(idx_H2j)
      pdj(78) =  &
          -k(599)*n(idx_HEj)  &
          -k(600)*n(idx_HEj)
      pdj(79) =  &
          -k(683)*n(idx_HNOj)
      pdj(81) =  &
          -k(519)*n(idx_H3j)
      pdj(82) =  &
          -k(677)*n(idx_H3COj)
      pdj(83) =  &
          -k(678)*n(idx_H3Oj)
      pdj(84) =  &
          -k(681)*n(idx_HCNHj)  &
          -k(680)*n(idx_HCNHj)
      pdj(87) =  &
          -k(684)*n(idx_N2Hj)
      pdj(88) =  &
          -k(685)*n(idx_O2Hj)
    elseif(j==22) then
      pdj(1) =  &
          +k(220)  &
          +k(983)
      pdj(2) =  &
          +k(221)  &
          +k(984)
      pdj(3) =  &
          -k(797)*n(idx_O)  &
          +k(799)*n(idx_OH)  &
          -k(798)*n(idx_O)
      pdj(5) =  &
          +k(792)*n(idx_NO)  &
          +k(784)*n(idx_CN)  &
          +k(891)*n(idx_N)  &
          +k(890)*n(idx_N)
      pdj(6) =  &
          +k(984)  &
          +k(849)*n(idx_H)  &
          +k(797)*n(idx_O)  &
          +k(504)*n(idx_H3j)  &
          +k(221)  &
          +k(566)*n(idx_HEj)  &
          +k(890)*n(idx_N)  &
          +k(800)*n(idx_OH)  &
          -k(838)*n(idx_H2)
      pdj(8) =  &
          +k(219)  &
          +k(70)*n(idx_Hj)  &
          +k(797)*n(idx_O)  &
          +k(838)*n(idx_H2)  &
          -k(849)*n(idx_H)  &
          +k(798)*n(idx_O)  &
          +k(889)*n(idx_N)  &
          +2.d0*k(891)*n(idx_N)  &
          +k(982)
      pdj(9) =  &
          -k(786)*n(idx_H2O)  &
          +k(801)*n(idx_OH)  &
          +k(792)*n(idx_NO)  &
          +k(794)*n(idx_O2)
      pdj(10) =  &
          +k(786)*n(idx_H2O)  &
          -k(799)*n(idx_OH)  &
          +k(793)*n(idx_O2)  &
          -k(801)*n(idx_OH)  &
          -k(800)*n(idx_OH)
      pdj(11) =  &
          -k(794)*n(idx_O2)  &
          +k(796)*n(idx_O2H)  &
          -k(793)*n(idx_O2)  &
          -k(795)*n(idx_O2)
      pdj(12) =  &
          +k(219)  &
          +k(849)*n(idx_H)  &
          +k(801)*n(idx_OH)  &
          +k(795)*n(idx_O2)  &
          +k(784)*n(idx_CN)  &
          +k(982)  &
          +2.d0*k(783)*n(idx_CH3)
      pdj(13) =  &
          +k(800)*n(idx_OH)  &
          -k(785)*n(idx_H2CO)  &
          +k(793)*n(idx_O2)  &
          +k(791)*n(idx_NO2)  &
          +k(798)*n(idx_O)
      pdj(14) =  &
          +k(794)*n(idx_O2)  &
          +k(785)*n(idx_H2CO)  &
          -k(787)*n(idx_HCO)
      pdj(16) =  &
          -k(790)*n(idx_NH3)
      pdj(17) =  &
          -k(792)*n(idx_NO)  &
          +k(788)*n(idx_HNO)
      pdj(18) =  &
          -k(784)*n(idx_CN)
      pdj(19) =  &
          +k(787)*n(idx_HCO)  &
          +k(797)*n(idx_O)
      pdj(21) =  &
          -k(789)*n(idx_NH2)  &
          +k(790)*n(idx_NH3)
      pdj(22) =  &
          -k(70)*n(idx_Hj)  &
          -k(800)*n(idx_OH)  &
          -k(798)*n(idx_O)  &
          -k(785)*n(idx_H2CO)  &
          -k(787)*n(idx_HCO)  &
          -k(566)*n(idx_HEj)  &
          -k(789)*n(idx_NH2)  &
          -k(890)*n(idx_N)  &
          -k(788)*n(idx_HNO)  &
          -k(793)*n(idx_O2)  &
          -k(791)*n(idx_NO2)  &
          -k(801)*n(idx_OH)  &
          -k(784)*n(idx_CN)  &
          -k(1079)  &
          -k(984)  &
          -k(889)*n(idx_N)  &
          -k(220)  &
          -k(982)  &
          -k(796)*n(idx_O2H)  &
          -k(799)*n(idx_OH)  &
          -k(838)*n(idx_H2)  &
          -k(794)*n(idx_O2)  &
          -k(891)*n(idx_N)  &
          -k(790)*n(idx_NH3)  &
          -k(797)*n(idx_O)  &
          -k(786)*n(idx_H2O)  &
          -k(849)*n(idx_H)  &
          -k(504)*n(idx_H3j)  &
          -k(221)  &
          -k(219)  &
          -k(795)*n(idx_O2)  &
          -4.d0*k(783)*n(idx_CH3)  &
          -k(983)  &
          -k(792)*n(idx_NO)
      pdj(23) =  &
          +k(785)*n(idx_H2CO)  &
          +k(789)*n(idx_NH2)  &
          +k(788)*n(idx_HNO)  &
          +k(799)*n(idx_OH)  &
          +k(838)*n(idx_H2)  &
          +k(796)*n(idx_O2H)  &
          +k(790)*n(idx_NH3)  &
          +k(786)*n(idx_H2O)  &
          +2.d0*k(783)*n(idx_CH3)  &
          +k(787)*n(idx_HCO)
      pdj(24) =  &
          -k(889)*n(idx_N)  &
          -k(891)*n(idx_N)  &
          -k(890)*n(idx_N)
      pdj(25) =  &
          +k(789)*n(idx_NH2)
      pdj(26) =  &
          +k(566)*n(idx_HEj)
      pdj(27) =  &
          +k(791)*n(idx_NO2)  &
          -k(788)*n(idx_HNO)
      pdj(30) =  &
          +k(889)*n(idx_N)
      pdj(32) =  &
          -k(791)*n(idx_NO2)
      pdj(33) =  &
          +k(795)*n(idx_O2)  &
          -k(796)*n(idx_O2H)
      pdj(38) =  &
          +k(1079)
      pdj(55) =  &
          -k(70)*n(idx_Hj)
      pdj(59) =  &
          +k(566)*n(idx_HEj)
      pdj(72) =  &
          +k(70)*n(idx_Hj)  &
          +k(220)  &
          +k(983)
      pdj(73) =  &
          +k(504)*n(idx_H3j)
      pdj(78) =  &
          -k(566)*n(idx_HEj)
      pdj(81) =  &
          -k(504)*n(idx_H3j)
    elseif(j==23) then
      pdj(1) =  &
          +k(992)
      pdj(2) =  &
          +k(993)
      pdj(3) =  &
          -k(936)*n(idx_O)  &
          +k(185)*n(idx_Oj)
      pdj(5) =  &
          +k(802)*n(idx_CN)
      pdj(6) =  &
          +k(850)*n(idx_H)  &
          +k(569)*n(idx_HEj)  &
          +k(432)*n(idx_Hj)  &
          +k(570)*n(idx_HEj)  &
          +k(86)*n(idx_H2j)  &
          +k(993)  &
          +k(615)*n(idx_Nj)  &
          +k(990)  &
          +k(224)  &
          +k(443)*n(idx_H2j)  &
          +k(397)*n(idx_N2j)
      pdj(8) =  &
          +k(614)*n(idx_Nj)  &
          +k(398)*n(idx_N2j)  &
          +k(569)*n(idx_HEj)  &
          +k(571)*n(idx_HEj)  &
          -k(850)*n(idx_H)  &
          +2.d0*k(616)*n(idx_Nj)  &
          +k(993)  &
          +k(71)*n(idx_Hj)  &
          +k(615)*n(idx_Nj)  &
          +k(991)  &
          +k(443)*n(idx_H2j)
      pdj(9) =  &
          +k(804)*n(idx_OH)
      pdj(10) =  &
          -k(804)*n(idx_OH)  &
          +k(705)*n(idx_Oj)  &
          +k(936)*n(idx_O)
      pdj(11) =  &
          -k(803)*n(idx_O2)
      pdj(12) =  &
          +k(224)  &
          +k(399)*n(idx_OHj)  &
          -k(761)*n(idx_CH2)  &
          +k(990)
      pdj(16) =  &
          +k(908)*n(idx_NH2)
      pdj(18) =  &
          -k(802)*n(idx_CN)
      pdj(19) =  &
          +k(46)*n(idx_COj)
      pdj(20) =  &
          +k(398)*n(idx_N2j)  &
          +k(397)*n(idx_N2j)
      pdj(21) =  &
          -k(908)*n(idx_NH2)  &
          +k(914)*n(idx_NH)
      pdj(22) =  &
          +k(395)*n(idx_H2Oj)  &
          +k(803)*n(idx_O2)  &
          +2.d0*k(761)*n(idx_CH2)  &
          +k(908)*n(idx_NH2)  &
          +k(396)*n(idx_HCNj)  &
          +k(936)*n(idx_O)  &
          +k(393)*n(idx_COj)  &
          +k(394)*n(idx_H2COj)  &
          +k(850)*n(idx_H)  &
          +k(804)*n(idx_OH)  &
          +k(991)  &
          +k(914)*n(idx_NH)  &
          +k(572)*n(idx_HEj)  &
          +k(802)*n(idx_CN)
      pdj(23) =  &
          -k(803)*n(idx_O2)  &
          -k(394)*n(idx_H2COj)  &
          -k(936)*n(idx_O)  &
          -k(761)*n(idx_CH2)  &
          -k(46)*n(idx_COj)  &
          -k(914)*n(idx_NH)  &
          -k(908)*n(idx_NH2)  &
          -k(569)*n(idx_HEj)  &
          -k(571)*n(idx_HEj)  &
          -k(990)  &
          -k(224)  &
          -k(993)  &
          -k(992)  &
          -k(850)*n(idx_H)  &
          -k(802)*n(idx_CN)  &
          -k(123)*n(idx_HEj)  &
          -k(1080)  &
          -k(432)*n(idx_Hj)  &
          -k(396)*n(idx_HCNj)  &
          -k(615)*n(idx_Nj)  &
          -k(136)*n(idx_Nj)  &
          -k(804)*n(idx_OH)  &
          -k(393)*n(idx_COj)  &
          -k(572)*n(idx_HEj)  &
          -k(705)*n(idx_Oj)  &
          -k(71)*n(idx_Hj)  &
          -k(570)*n(idx_HEj)  &
          -k(399)*n(idx_OHj)  &
          -k(185)*n(idx_Oj)  &
          -k(614)*n(idx_Nj)  &
          -k(86)*n(idx_H2j)  &
          -k(443)*n(idx_H2j)  &
          -k(397)*n(idx_N2j)  &
          -k(616)*n(idx_Nj)  &
          -k(398)*n(idx_N2j)  &
          -k(395)*n(idx_H2Oj)  &
          -k(991)
      pdj(24) =  &
          +k(614)*n(idx_Nj)  &
          +k(136)*n(idx_Nj)
      pdj(25) =  &
          -k(914)*n(idx_NH)
      pdj(26) =  &
          +k(571)*n(idx_HEj)  &
          +k(570)*n(idx_HEj)  &
          +k(572)*n(idx_HEj)  &
          +k(123)*n(idx_HEj)  &
          +k(569)*n(idx_HEj)
      pdj(33) =  &
          +k(803)*n(idx_O2)
      pdj(38) =  &
          +k(1080)
      pdj(54) =  &
          +k(393)*n(idx_COj)
      pdj(55) =  &
          +k(572)*n(idx_HEj)  &
          -k(432)*n(idx_Hj)  &
          -k(71)*n(idx_Hj)
      pdj(58) =  &
          +k(570)*n(idx_HEj)  &
          +k(397)*n(idx_N2j)
      pdj(59) =  &
          +k(569)*n(idx_HEj)
      pdj(60) =  &
          -k(394)*n(idx_H2COj)
      pdj(65) =  &
          -k(46)*n(idx_COj)  &
          -k(393)*n(idx_COj)
      pdj(66) =  &
          -k(397)*n(idx_N2j)  &
          -k(398)*n(idx_N2j)
      pdj(68) =  &
          -k(395)*n(idx_H2Oj)
      pdj(70) =  &
          -k(705)*n(idx_Oj)  &
          -k(185)*n(idx_Oj)
      pdj(71) =  &
          -k(399)*n(idx_OHj)
      pdj(72) =  &
          +k(614)*n(idx_Nj)  &
          +k(398)*n(idx_N2j)  &
          +k(432)*n(idx_Hj)  &
          +k(571)*n(idx_HEj)  &
          +k(705)*n(idx_Oj)  &
          +k(443)*n(idx_H2j)
      pdj(73) =  &
          +k(992)  &
          +k(185)*n(idx_Oj)  &
          +k(123)*n(idx_HEj)  &
          +k(86)*n(idx_H2j)  &
          +k(136)*n(idx_Nj)  &
          +k(46)*n(idx_COj)  &
          +k(71)*n(idx_Hj)
      pdj(74) =  &
          -k(616)*n(idx_Nj)  &
          -k(615)*n(idx_Nj)  &
          -k(136)*n(idx_Nj)  &
          -k(614)*n(idx_Nj)
      pdj(75) =  &
          -k(396)*n(idx_HCNj)  &
          +k(615)*n(idx_Nj)
      pdj(77) =  &
          -k(443)*n(idx_H2j)  &
          -k(86)*n(idx_H2j)
      pdj(78) =  &
          -k(569)*n(idx_HEj)  &
          -k(571)*n(idx_HEj)  &
          -k(123)*n(idx_HEj)  &
          -k(570)*n(idx_HEj)  &
          -k(572)*n(idx_HEj)
      pdj(82) =  &
          +k(394)*n(idx_H2COj)
      pdj(83) =  &
          +k(399)*n(idx_OHj)  &
          +k(395)*n(idx_H2Oj)
      pdj(84) =  &
          +k(616)*n(idx_Nj)  &
          +k(396)*n(idx_HCNj)
    elseif(j==24) then
      pdj(1) =  &
          +k(242)  &
          +k(213)
      pdj(2) =  &
          -k(811)*n(idx_CH)  &
          +k(888)*n(idx_CH2)  &
          -k(810)*n(idx_CH)
      pdj(3) =  &
          +k(640)*n(idx_O2j)  &
          +2.d0*k(900)*n(idx_NO2)  &
          +k(896)*n(idx_HCO)  &
          +k(904)*n(idx_O2)  &
          +k(903)*n(idx_NO)  &
          +k(907)*n(idx_OH)
      pdj(4) =  &
          +k(887)*n(idx_CH2)
      pdj(5) =  &
          +k(891)*n(idx_CH3)  &
          +k(896)*n(idx_HCO)  &
          +k(886)*n(idx_CH2)  &
          +k(890)*n(idx_CH3)  &
          +k(894)*n(idx_H2CN)
      pdj(6) =  &
          -k(841)*n(idx_H2)  &
          +k(890)*n(idx_CH3)  &
          +k(637)*n(idx_H2Oj)
      pdj(7) =  &
          +k(635)*n(idx_CNj)  &
          +k(811)*n(idx_CH)  &
          -k(1042)*n(idx_C)  &
          +k(892)*n(idx_CN)
      pdj(8) =  &
          +k(906)*n(idx_OH)  &
          +k(638)*n(idx_NHj)  &
          +k(641)*n(idx_OHj)  &
          +k(889)*n(idx_CH3)  &
          +k(887)*n(idx_CH2)  &
          +k(639)*n(idx_NH2j)  &
          +k(899)*n(idx_NH)  &
          +2.d0*k(891)*n(idx_CH3)  &
          +k(636)*n(idx_H2Oj)  &
          +k(897)*n(idx_HCO)  &
          +k(886)*n(idx_CH2)  &
          +k(810)*n(idx_CH)  &
          +k(841)*n(idx_H2)  &
          +k(352)*n(idx_CHj)  &
          +k(454)*n(idx_H2j)  &
          +k(634)*n(idx_CH2j)
      pdj(10) =  &
          -k(907)*n(idx_OH)  &
          -k(906)*n(idx_OH)
      pdj(11) =  &
          -k(904)*n(idx_O2)  &
          +k(905)*n(idx_O2H)  &
          +k(902)*n(idx_NO2)
      pdj(12) =  &
          -k(887)*n(idx_CH2)  &
          -k(886)*n(idx_CH2)  &
          -k(888)*n(idx_CH2)
      pdj(14) =  &
          -k(897)*n(idx_HCO)  &
          -k(895)*n(idx_HCO)  &
          -k(896)*n(idx_HCO)
      pdj(17) =  &
          -k(903)*n(idx_NO)  &
          +k(893)*n(idx_CO2)  &
          +2.d0*k(901)*n(idx_NO2)  &
          +k(906)*n(idx_OH)  &
          +k(898)*n(idx_HNO)  &
          +k(904)*n(idx_O2)
      pdj(18) =  &
          +k(810)*n(idx_CH)  &
          +k(1042)*n(idx_C)  &
          -k(892)*n(idx_CN)
      pdj(19) =  &
          +k(895)*n(idx_HCO)  &
          +k(893)*n(idx_CO2)
      pdj(20) =  &
          +k(892)*n(idx_CN)  &
          +k(900)*n(idx_NO2)  &
          +k(899)*n(idx_NH)  &
          +k(154)*n(idx_N2j)  &
          +k(903)*n(idx_NO)  &
          +k(902)*n(idx_NO2)
      pdj(22) =  &
          -k(891)*n(idx_CH3)  &
          -k(889)*n(idx_CH3)  &
          -k(890)*n(idx_CH3)
      pdj(24) =  &
          -k(907)*n(idx_OH)  &
          -k(905)*n(idx_O2H)  &
          -k(891)*n(idx_CH3)  &
          -k(895)*n(idx_HCO)  &
          -k(841)*n(idx_H2)  &
          -k(897)*n(idx_HCO)  &
          -k(640)*n(idx_O2j)  &
          -k(899)*n(idx_NH)  &
          -k(810)*n(idx_CH)  &
          -k(906)*n(idx_OH)  &
          -k(902)*n(idx_NO2)  &
          -k(213)  &
          -k(892)*n(idx_CN)  &
          -k(888)*n(idx_CH2)  &
          -k(639)*n(idx_NH2j)  &
          -k(352)*n(idx_CHj)  &
          -k(637)*n(idx_H2Oj)  &
          -k(894)*n(idx_H2CN)  &
          -k(889)*n(idx_CH3)  &
          -k(890)*n(idx_CH3)  &
          -k(903)*n(idx_NO)  &
          -k(1054)*n(idx_Nj)  &
          -k(454)*n(idx_H2j)  &
          -k(898)*n(idx_HNO)  &
          -k(1110)  &
          -k(887)*n(idx_CH2)  &
          -k(634)*n(idx_CH2j)  &
          -k(893)*n(idx_CO2)  &
          -k(1040)*n(idx_Cj)  &
          -k(896)*n(idx_HCO)  &
          -k(886)*n(idx_CH2)  &
          -k(638)*n(idx_NHj)  &
          -k(154)*n(idx_N2j)  &
          -k(636)*n(idx_H2Oj)  &
          -k(641)*n(idx_OHj)  &
          -k(904)*n(idx_O2)  &
          -k(901)*n(idx_NO2)  &
          -k(811)*n(idx_CH)  &
          -k(900)*n(idx_NO2)  &
          -k(1042)*n(idx_C)  &
          -k(242)  &
          -k(635)*n(idx_CNj)
      pdj(25) =  &
          +k(888)*n(idx_CH2)  &
          +k(895)*n(idx_HCO)  &
          -k(899)*n(idx_NH)  &
          +k(905)*n(idx_O2H)  &
          +k(898)*n(idx_HNO)  &
          +k(894)*n(idx_H2CN)  &
          +k(841)*n(idx_H2)  &
          +k(811)*n(idx_CH)  &
          +k(907)*n(idx_OH)
      pdj(27) =  &
          -k(898)*n(idx_HNO)
      pdj(29) =  &
          -k(893)*n(idx_CO2)
      pdj(30) =  &
          -k(894)*n(idx_H2CN)  &
          +k(889)*n(idx_CH3)
      pdj(32) =  &
          -k(902)*n(idx_NO2)  &
          -k(900)*n(idx_NO2)  &
          -k(901)*n(idx_NO2)
      pdj(33) =  &
          -k(905)*n(idx_O2H)
      pdj(34) =  &
          +k(897)*n(idx_HCO)
      pdj(45) =  &
          +k(1110)
      pdj(57) =  &
          -k(1040)*n(idx_Cj)
      pdj(58) =  &
          -k(634)*n(idx_CH2j)
      pdj(59) =  &
          -k(352)*n(idx_CHj)
      pdj(63) =  &
          +k(640)*n(idx_O2j)  &
          +k(641)*n(idx_OHj)  &
          +k(637)*n(idx_H2Oj)
      pdj(64) =  &
          +k(1040)*n(idx_Cj)  &
          -k(635)*n(idx_CNj)  &
          +k(352)*n(idx_CHj)
      pdj(66) =  &
          -k(154)*n(idx_N2j)  &
          +k(635)*n(idx_CNj)  &
          +k(1054)*n(idx_Nj)  &
          +k(638)*n(idx_NHj)
      pdj(67) =  &
          -k(640)*n(idx_O2j)
      pdj(68) =  &
          -k(637)*n(idx_H2Oj)  &
          -k(636)*n(idx_H2Oj)
      pdj(69) =  &
          -k(639)*n(idx_NH2j)
      pdj(71) =  &
          -k(641)*n(idx_OHj)
      pdj(74) =  &
          +k(154)*n(idx_N2j)  &
          +k(242)  &
          -k(1054)*n(idx_Nj)  &
          +k(213)
      pdj(75) =  &
          +k(634)*n(idx_CH2j)
      pdj(76) =  &
          -k(638)*n(idx_NHj)  &
          +k(454)*n(idx_H2j)
      pdj(77) =  &
          -k(454)*n(idx_H2j)
      pdj(79) =  &
          +k(636)*n(idx_H2Oj)
      pdj(87) =  &
          +k(639)*n(idx_NH2j)
    elseif(j==25) then
      pdj(1) =  &
          +k(249)  &
          +k(1027)
      pdj(2) =  &
          +k(753)*n(idx_C)
      pdj(3) =  &
          -k(927)*n(idx_O)  &
          +k(701)*n(idx_OHj)  &
          +k(924)*n(idx_O2)  &
          +k(930)*n(idx_OH)  &
          +k(699)*n(idx_O2j)  &
          +k(181)*n(idx_Oj)  &
          +k(922)*n(idx_NO)  &
          -k(926)*n(idx_O)
      pdj(5) =  &
          +k(915)*n(idx_CN)
      pdj(6) =  &
          +k(866)*n(idx_H)  &
          +2.d0*k(918)*n(idx_NH)  &
          +k(520)*n(idx_H3j)  &
          +k(689)*n(idx_CH3j)  &
          +k(96)*n(idx_H2j)  &
          +k(354)*n(idx_CHj)  &
          -k(843)*n(idx_H2)
      pdj(7) =  &
          -k(752)*n(idx_C)  &
          -k(753)*n(idx_C)
      pdj(8) =  &
          +k(80)*n(idx_Hj)  &
          +k(698)*n(idx_Oj)  &
          +k(929)*n(idx_OH)  &
          +4.d0*k(919)*n(idx_NH)  &
          -k(866)*n(idx_H)  &
          +k(624)*n(idx_Nj)  &
          +k(899)*n(idx_N)  &
          +k(922)*n(idx_NO)  &
          +k(843)*n(idx_H2)  &
          +k(926)*n(idx_O)  &
          +k(248)  &
          +k(1026)  &
          +k(752)*n(idx_C)  &
          +k(603)*n(idx_HEj)  &
          +k(455)*n(idx_H2j)  &
          +k(324)*n(idx_Cj)
      pdj(9) =  &
          +k(928)*n(idx_OH)  &
          -k(916)*n(idx_H2O)
      pdj(10) =  &
          +k(923)*n(idx_NO)  &
          +k(927)*n(idx_O)  &
          +k(916)*n(idx_H2O)  &
          -k(928)*n(idx_OH)  &
          -k(929)*n(idx_OH)  &
          -k(930)*n(idx_OH)  &
          +k(925)*n(idx_O2)
      pdj(11) =  &
          -k(924)*n(idx_O2)  &
          -k(925)*n(idx_O2)  &
          +k(700)*n(idx_O2Hj)
      pdj(16) =  &
          -k(917)*n(idx_NH3)
      pdj(17) =  &
          -k(922)*n(idx_NO)  &
          +k(926)*n(idx_O)  &
          +k(695)*n(idx_HNOj)  &
          +k(921)*n(idx_NO2)  &
          -k(923)*n(idx_NO)  &
          +k(925)*n(idx_O2)
      pdj(18) =  &
          -k(915)*n(idx_CN)  &
          +k(693)*n(idx_HCNj)  &
          +k(752)*n(idx_C)  &
          +k(178)*n(idx_CNj)
      pdj(19) =  &
          +k(694)*n(idx_HCOj)  &
          +k(179)*n(idx_COj)
      pdj(20) =  &
          +2.d0*k(919)*n(idx_NH)  &
          +2.d0*k(918)*n(idx_NH)  &
          +k(923)*n(idx_NO)  &
          +k(696)*n(idx_N2Hj)  &
          +k(899)*n(idx_N)  &
          +k(180)*n(idx_N2j)  &
          +k(922)*n(idx_NO)
      pdj(21) =  &
          +2.d0*k(917)*n(idx_NH3)  &
          +k(914)*n(idx_CH4)  &
          +2.d0*k(920)*n(idx_NH)  &
          +k(930)*n(idx_OH)  &
          +k(843)*n(idx_H2)  &
          +k(916)*n(idx_H2O)
      pdj(22) =  &
          +k(914)*n(idx_CH4)
      pdj(23) =  &
          -k(914)*n(idx_CH4)
      pdj(24) =  &
          +k(866)*n(idx_H)  &
          +2.d0*k(920)*n(idx_NH)  &
          +k(927)*n(idx_O)  &
          -k(899)*n(idx_N)  &
          +k(915)*n(idx_CN)  &
          +k(697)*n(idx_NH2j)  &
          +k(928)*n(idx_OH)  &
          +k(692)*n(idx_H2Oj)  &
          +k(691)*n(idx_H2COj)  &
          +k(248)  &
          +k(1026)  &
          +k(146)*n(idx_Nj)  &
          +k(690)*n(idx_COj)  &
          +k(753)*n(idx_C)  &
          +k(658)*n(idx_NHj)
      pdj(25) =  &
          -k(1026)  &
          -k(690)*n(idx_COj)  &
          -k(930)*n(idx_OH)  &
          -k(752)*n(idx_C)  &
          -k(96)*n(idx_H2j)  &
          -k(520)*n(idx_H3j)  &
          -k(700)*n(idx_O2Hj)  &
          -k(916)*n(idx_H2O)  &
          -k(691)*n(idx_H2COj)  &
          -k(180)*n(idx_N2j)  &
          -k(692)*n(idx_H2Oj)  &
          -k(455)*n(idx_H2j)  &
          -k(249)  &
          -k(693)*n(idx_HCNj)  &
          -k(866)*n(idx_H)  &
          -k(927)*n(idx_O)  &
          -k(928)*n(idx_OH)  &
          -k(929)*n(idx_OH)  &
          -k(699)*n(idx_O2j)  &
          -k(1085)  &
          -k(843)*n(idx_H2)  &
          -k(922)*n(idx_NO)  &
          -k(80)*n(idx_Hj)  &
          -k(914)*n(idx_CH4)  &
          -k(248)  &
          -k(921)*n(idx_NO2)  &
          -k(701)*n(idx_OHj)  &
          -k(926)*n(idx_O)  &
          -k(181)*n(idx_Oj)  &
          -k(178)*n(idx_CNj)  &
          -4.d0*k(918)*n(idx_NH)  &
          -k(179)*n(idx_COj)  &
          -k(923)*n(idx_NO)  &
          -k(698)*n(idx_Oj)  &
          -k(695)*n(idx_HNOj)  &
          -k(753)*n(idx_C)  &
          -k(689)*n(idx_CH3j)  &
          -4.d0*k(920)*n(idx_NH)  &
          -k(603)*n(idx_HEj)  &
          -k(694)*n(idx_HCOj)  &
          -k(696)*n(idx_N2Hj)  &
          -k(925)*n(idx_O2)  &
          -k(899)*n(idx_N)  &
          -k(354)*n(idx_CHj)  &
          -k(146)*n(idx_Nj)  &
          -4.d0*k(919)*n(idx_NH)  &
          -k(915)*n(idx_CN)  &
          -k(658)*n(idx_NHj)  &
          -k(917)*n(idx_NH3)  &
          -k(324)*n(idx_Cj)  &
          -k(624)*n(idx_Nj)  &
          -k(924)*n(idx_O2)  &
          -k(697)*n(idx_NH2j)  &
          -k(1027)
      pdj(26) =  &
          +k(603)*n(idx_HEj)
      pdj(27) =  &
          +k(929)*n(idx_OH)  &
          +k(921)*n(idx_NO2)  &
          +k(924)*n(idx_O2)
      pdj(32) =  &
          -k(921)*n(idx_NO2)
      pdj(45) =  &
          +k(1085)
      pdj(54) =  &
          -k(694)*n(idx_HCOj)  &
          +k(690)*n(idx_COj)
      pdj(55) =  &
          -k(80)*n(idx_Hj)
      pdj(57) =  &
          -k(324)*n(idx_Cj)
      pdj(59) =  &
          -k(354)*n(idx_CHj)
      pdj(60) =  &
          -k(691)*n(idx_H2COj)
      pdj(62) =  &
          +k(697)*n(idx_NH2j)
      pdj(63) =  &
          +k(698)*n(idx_Oj)
      pdj(64) =  &
          -k(178)*n(idx_CNj)  &
          +k(324)*n(idx_Cj)  &
          +k(354)*n(idx_CHj)
      pdj(65) =  &
          -k(690)*n(idx_COj)  &
          -k(179)*n(idx_COj)
      pdj(66) =  &
          -k(180)*n(idx_N2j)  &
          +k(624)*n(idx_Nj)
      pdj(67) =  &
          -k(699)*n(idx_O2j)
      pdj(68) =  &
          -k(692)*n(idx_H2Oj)
      pdj(69) =  &
          +k(520)*n(idx_H3j)  &
          +k(700)*n(idx_O2Hj)  &
          +k(694)*n(idx_HCOj)  &
          +k(701)*n(idx_OHj)  &
          +k(695)*n(idx_HNOj)  &
          +k(658)*n(idx_NHj)  &
          +k(693)*n(idx_HCNj)  &
          -k(697)*n(idx_NH2j)  &
          +k(455)*n(idx_H2j)  &
          +k(696)*n(idx_N2Hj)
      pdj(70) =  &
          -k(181)*n(idx_Oj)  &
          -k(698)*n(idx_Oj)
      pdj(71) =  &
          -k(701)*n(idx_OHj)
      pdj(72) =  &
          -k(689)*n(idx_CH3j)
      pdj(74) =  &
          +k(603)*n(idx_HEj)  &
          -k(624)*n(idx_Nj)  &
          -k(146)*n(idx_Nj)
      pdj(75) =  &
          -k(693)*n(idx_HCNj)
      pdj(76) =  &
          +k(80)*n(idx_Hj)  &
          +k(1027)  &
          +k(96)*n(idx_H2j)  &
          +k(179)*n(idx_COj)  &
          +k(146)*n(idx_Nj)  &
          +k(249)  &
          +k(180)*n(idx_N2j)  &
          +k(181)*n(idx_Oj)  &
          +k(178)*n(idx_CNj)  &
          -k(658)*n(idx_NHj)
      pdj(77) =  &
          -k(455)*n(idx_H2j)  &
          -k(96)*n(idx_H2j)
      pdj(78) =  &
          -k(603)*n(idx_HEj)
      pdj(79) =  &
          -k(695)*n(idx_HNOj)  &
          +k(699)*n(idx_O2j)
      pdj(81) =  &
          -k(520)*n(idx_H3j)
      pdj(82) =  &
          +k(691)*n(idx_H2COj)
      pdj(83) =  &
          +k(692)*n(idx_H2Oj)
      pdj(84) =  &
          +k(689)*n(idx_CH3j)
      pdj(87) =  &
          -k(696)*n(idx_N2Hj)
      pdj(88) =  &
          -k(700)*n(idx_O2Hj)
    elseif(j==26) then
      pdj(1) =  &
          +k(212)  &
          +k(239)
      pdj(8) =  &
          +k(452)*n(idx_H2j)
      pdj(26) =  &
          -k(239)  &
          -k(452)*n(idx_H2j)  &
          -k(1046)*n(idx_Hj)  &
          -k(212)
      pdj(55) =  &
          -k(1046)*n(idx_Hj)
      pdj(77) =  &
          -k(452)*n(idx_H2j)
      pdj(78) =  &
          +k(212)  &
          +k(239)
      pdj(86) =  &
          +k(452)*n(idx_H2j)  &
          +k(1046)*n(idx_Hj)
    elseif(j==27) then
      pdj(2) =  &
          -k(808)*n(idx_CH)
      pdj(3) =  &
          -k(948)*n(idx_O)  &
          +k(861)*n(idx_H)  &
          -k(950)*n(idx_O)  &
          -k(949)*n(idx_O)
      pdj(5) =  &
          +k(826)*n(idx_CN)
      pdj(6) =  &
          +k(516)*n(idx_H3j)  &
          +k(862)*n(idx_H)  &
          +k(439)*n(idx_Hj)
      pdj(8) =  &
          +k(948)*n(idx_O)  &
          -k(862)*n(idx_H)  &
          -k(861)*n(idx_H)  &
          +k(596)*n(idx_HEj)  &
          +k(1017)  &
          +k(238)  &
          -k(863)*n(idx_H)
      pdj(9) =  &
          +k(968)*n(idx_OH)
      pdj(10) =  &
          +k(863)*n(idx_H)  &
          +k(949)*n(idx_O)  &
          -k(968)*n(idx_OH)
      pdj(11) =  &
          +k(950)*n(idx_O)
      pdj(12) =  &
          +k(808)*n(idx_CH)  &
          -k(765)*n(idx_CH2)
      pdj(13) =  &
          +k(880)*n(idx_HCO)
      pdj(14) =  &
          -k(880)*n(idx_HCO)
      pdj(17) =  &
          +k(949)*n(idx_O)  &
          +k(826)*n(idx_CN)  &
          +k(880)*n(idx_HCO)  &
          +k(788)*n(idx_CH3)  &
          +k(765)*n(idx_CH2)  &
          +k(597)*n(idx_HEj)  &
          +k(968)*n(idx_OH)  &
          +k(898)*n(idx_N)  &
          +k(862)*n(idx_H)  &
          +k(1017)  &
          +k(808)*n(idx_CH)  &
          +k(238)
      pdj(18) =  &
          -k(826)*n(idx_CN)
      pdj(19) =  &
          -k(832)*n(idx_CO)
      pdj(21) =  &
          +k(861)*n(idx_H)
      pdj(22) =  &
          -k(788)*n(idx_CH3)  &
          +k(765)*n(idx_CH2)
      pdj(23) =  &
          +k(788)*n(idx_CH3)
      pdj(24) =  &
          -k(898)*n(idx_N)
      pdj(25) =  &
          +k(863)*n(idx_H)  &
          +k(950)*n(idx_O)  &
          +k(832)*n(idx_CO)  &
          +k(898)*n(idx_N)
      pdj(26) =  &
          +k(597)*n(idx_HEj)  &
          +k(596)*n(idx_HEj)
      pdj(27) =  &
          -k(862)*n(idx_H)  &
          -k(238)  &
          -k(861)*n(idx_H)  &
          -k(597)*n(idx_HEj)  &
          -k(439)*n(idx_Hj)  &
          -k(1017)  &
          -k(788)*n(idx_CH3)  &
          -k(948)*n(idx_O)  &
          -k(880)*n(idx_HCO)  &
          -k(765)*n(idx_CH2)  &
          -k(826)*n(idx_CN)  &
          -k(516)*n(idx_H3j)  &
          -k(808)*n(idx_CH)  &
          -k(950)*n(idx_O)  &
          -k(1114)  &
          -k(898)*n(idx_N)  &
          -k(968)*n(idx_OH)  &
          -k(949)*n(idx_O)  &
          -k(596)*n(idx_HEj)  &
          -k(832)*n(idx_CO)  &
          -k(863)*n(idx_H)
      pdj(29) =  &
          +k(832)*n(idx_CO)
      pdj(32) =  &
          +k(948)*n(idx_O)
      pdj(48) =  &
          +k(1114)
      pdj(55) =  &
          +k(597)*n(idx_HEj)  &
          -k(439)*n(idx_Hj)
      pdj(63) =  &
          +k(439)*n(idx_Hj)  &
          +k(596)*n(idx_HEj)
      pdj(78) =  &
          -k(597)*n(idx_HEj)  &
          -k(596)*n(idx_HEj)
      pdj(80) =  &
          +k(516)*n(idx_H3j)
      pdj(81) =  &
          -k(516)*n(idx_H3j)
    elseif(j==28) then
      pdj(1) =  &
          +k(986)
      pdj(2) =  &
          +k(314)*n(idx_Cj)
      pdj(6) =  &
          +k(222)  &
          +k(505)*n(idx_H3j)  &
          +k(985)  &
          +2.d0*k(431)*n(idx_Hj)  &
          +k(430)*n(idx_Hj)
      pdj(8) =  &
          +k(715)*n(idx_O2j)  &
          +k(612)*n(idx_Nj)  &
          +k(610)*n(idx_Nj)  &
          +k(613)*n(idx_Nj)  &
          +k(986)
      pdj(9) =  &
          +k(429)*n(idx_Hj)  &
          +k(703)*n(idx_Oj)  &
          +k(505)*n(idx_H3j)
      pdj(10) =  &
          +k(223)  &
          +k(987)  &
          +k(704)*n(idx_Oj)  &
          +k(568)*n(idx_HEj)
      pdj(11) =  &
          +k(715)*n(idx_O2j)
      pdj(12) =  &
          +k(341)*n(idx_CHj)
      pdj(13) =  &
          +k(222)  &
          +k(985)  &
          +k(340)*n(idx_CHj)
      pdj(14) =  &
          +k(315)*n(idx_Cj)
      pdj(17) =  &
          +k(613)*n(idx_Nj)
      pdj(22) =  &
          +k(223)  &
          +k(612)*n(idx_Nj)  &
          +k(987)  &
          +k(567)*n(idx_HEj)
      pdj(23) =  &
          +k(382)*n(idx_CH3j)
      pdj(25) =  &
          +k(610)*n(idx_Nj)  &
          +k(611)*n(idx_Nj)
      pdj(26) =  &
          +k(567)*n(idx_HEj)  &
          +k(568)*n(idx_HEj)
      pdj(28) =  &
          -k(610)*n(idx_Nj)  &
          -k(715)*n(idx_O2j)  &
          -k(222)  &
          -k(1064)  &
          -k(315)*n(idx_Cj)  &
          -k(382)*n(idx_CH3j)  &
          -k(612)*n(idx_Nj)  &
          -k(986)  &
          -k(987)  &
          -k(505)*n(idx_H3j)  &
          -k(568)*n(idx_HEj)  &
          -k(613)*n(idx_Nj)  &
          -k(314)*n(idx_Cj)  &
          -k(985)  &
          -k(223)  &
          -k(340)*n(idx_CHj)  &
          -k(704)*n(idx_Oj)  &
          -k(430)*n(idx_Hj)  &
          -k(341)*n(idx_CHj)  &
          -k(429)*n(idx_Hj)  &
          -k(703)*n(idx_Oj)  &
          -k(431)*n(idx_Hj)  &
          -k(611)*n(idx_Nj)  &
          -k(567)*n(idx_HEj)
      pdj(35) =  &
          +k(1064)
      pdj(54) =  &
          +k(431)*n(idx_Hj)
      pdj(55) =  &
          -k(430)*n(idx_Hj)  &
          -k(429)*n(idx_Hj)  &
          -k(431)*n(idx_Hj)
      pdj(57) =  &
          -k(314)*n(idx_Cj)  &
          -k(315)*n(idx_Cj)
      pdj(59) =  &
          -k(340)*n(idx_CHj)  &
          -k(341)*n(idx_CHj)
      pdj(60) =  &
          +k(703)*n(idx_Oj)  &
          +k(610)*n(idx_Nj)
      pdj(63) =  &
          +k(612)*n(idx_Nj)
      pdj(67) =  &
          -k(715)*n(idx_O2j)
      pdj(70) =  &
          -k(704)*n(idx_Oj)  &
          -k(703)*n(idx_Oj)
      pdj(71) =  &
          +k(567)*n(idx_HEj)
      pdj(72) =  &
          +k(429)*n(idx_Hj)  &
          +k(315)*n(idx_Cj)  &
          +k(613)*n(idx_Nj)  &
          +k(340)*n(idx_CHj)  &
          -k(382)*n(idx_CH3j)  &
          +k(505)*n(idx_H3j)  &
          +k(568)*n(idx_HEj)
      pdj(74) =  &
          -k(610)*n(idx_Nj)  &
          -k(612)*n(idx_Nj)  &
          -k(611)*n(idx_Nj)  &
          -k(613)*n(idx_Nj)
      pdj(78) =  &
          -k(568)*n(idx_HEj)  &
          -k(567)*n(idx_HEj)
      pdj(81) =  &
          -k(505)*n(idx_H3j)
      pdj(82) =  &
          +k(314)*n(idx_Cj)  &
          +k(704)*n(idx_Oj)  &
          +k(430)*n(idx_Hj)  &
          +k(986)  &
          +k(341)*n(idx_CHj)  &
          +k(715)*n(idx_O2j)  &
          +k(611)*n(idx_Nj)  &
          +k(382)*n(idx_CH3j)
    elseif(j==29) then
      pdj(2) =  &
          -k(805)*n(idx_CH)
      pdj(3) =  &
          +k(576)*n(idx_HEj)  &
          +k(727)*n(idx_OHj)  &
          +k(998)  &
          -k(939)*n(idx_O)  &
          +k(433)*n(idx_Hj)  &
          +k(227)
      pdj(6) =  &
          +k(508)*n(idx_H3j)
      pdj(7) =  &
          +k(578)*n(idx_HEj)
      pdj(8) =  &
          +k(446)*n(idx_H2j)  &
          -k(852)*n(idx_H)
      pdj(10) =  &
          +k(852)*n(idx_H)
      pdj(11) =  &
          +k(579)*n(idx_HEj)  &
          +k(939)*n(idx_O)  &
          +k(716)*n(idx_O2Hj)
      pdj(14) =  &
          +k(645)*n(idx_NHj)  &
          +k(805)*n(idx_CH)
      pdj(17) =  &
          +k(617)*n(idx_Nj)  &
          +k(563)*n(idx_HNOj)  &
          +k(893)*n(idx_N)
      pdj(18) =  &
          +k(535)*n(idx_HCNj)
      pdj(19) =  &
          +k(998)  &
          +k(893)*n(idx_N)  &
          +k(316)*n(idx_Cj)  &
          +k(707)*n(idx_Oj)  &
          +k(342)*n(idx_CHj)  &
          +k(360)*n(idx_CH2j)  &
          +k(939)*n(idx_O)  &
          +k(852)*n(idx_H)  &
          +k(644)*n(idx_NHj)  &
          +k(227)  &
          +k(577)*n(idx_HEj)  &
          +k(805)*n(idx_CH)
      pdj(20) =  &
          +k(632)*n(idx_N2Hj)
      pdj(22) =  &
          +k(389)*n(idx_CH4j)
      pdj(24) =  &
          -k(893)*n(idx_N)  &
          +k(643)*n(idx_NHj)
      pdj(26) =  &
          +k(579)*n(idx_HEj)  &
          +k(578)*n(idx_HEj)  &
          +k(577)*n(idx_HEj)  &
          +k(576)*n(idx_HEj)
      pdj(29) =  &
          -k(707)*n(idx_Oj)  &
          -k(577)*n(idx_HEj)  &
          -k(716)*n(idx_O2Hj)  &
          -k(508)*n(idx_H3j)  &
          -k(535)*n(idx_HCNj)  &
          -k(643)*n(idx_NHj)  &
          -k(433)*n(idx_Hj)  &
          -k(805)*n(idx_CH)  &
          -k(939)*n(idx_O)  &
          -k(727)*n(idx_OHj)  &
          -k(617)*n(idx_Nj)  &
          -k(644)*n(idx_NHj)  &
          -k(578)*n(idx_HEj)  &
          -k(563)*n(idx_HNOj)  &
          -k(576)*n(idx_HEj)  &
          -k(342)*n(idx_CHj)  &
          -k(998)  &
          -k(446)*n(idx_H2j)  &
          -k(1078)  &
          -k(645)*n(idx_NHj)  &
          -k(389)*n(idx_CH4j)  &
          -k(360)*n(idx_CH2j)  &
          -k(579)*n(idx_HEj)  &
          -k(852)*n(idx_H)  &
          -k(893)*n(idx_N)  &
          -k(632)*n(idx_N2Hj)  &
          -k(316)*n(idx_Cj)  &
          -k(227)
      pdj(42) =  &
          +k(1078)
      pdj(54) =  &
          +k(342)*n(idx_CHj)  &
          +k(433)*n(idx_Hj)
      pdj(55) =  &
          -k(433)*n(idx_Hj)
      pdj(57) =  &
          +k(579)*n(idx_HEj)  &
          -k(316)*n(idx_Cj)
      pdj(58) =  &
          -k(360)*n(idx_CH2j)
      pdj(59) =  &
          -k(342)*n(idx_CHj)
      pdj(60) =  &
          +k(360)*n(idx_CH2j)
      pdj(63) =  &
          +k(645)*n(idx_NHj)
      pdj(65) =  &
          +k(316)*n(idx_Cj)  &
          +k(617)*n(idx_Nj)  &
          +k(576)*n(idx_HEj)
      pdj(67) =  &
          +k(578)*n(idx_HEj)  &
          +k(707)*n(idx_Oj)
      pdj(70) =  &
          -k(707)*n(idx_Oj)  &
          +k(577)*n(idx_HEj)
      pdj(71) =  &
          -k(727)*n(idx_OHj)
      pdj(73) =  &
          -k(389)*n(idx_CH4j)
      pdj(74) =  &
          -k(617)*n(idx_Nj)
      pdj(75) =  &
          -k(535)*n(idx_HCNj)
      pdj(76) =  &
          -k(643)*n(idx_NHj)  &
          -k(644)*n(idx_NHj)  &
          -k(645)*n(idx_NHj)
      pdj(77) =  &
          -k(446)*n(idx_H2j)
      pdj(78) =  &
          -k(579)*n(idx_HEj)  &
          -k(578)*n(idx_HEj)  &
          -k(577)*n(idx_HEj)  &
          -k(576)*n(idx_HEj)
      pdj(79) =  &
          +k(644)*n(idx_NHj)  &
          -k(563)*n(idx_HNOj)
      pdj(81) =  &
          -k(508)*n(idx_H3j)
      pdj(85) =  &
          +k(563)*n(idx_HNOj)  &
          +k(535)*n(idx_HCNj)  &
          +k(446)*n(idx_H2j)  &
          +k(727)*n(idx_OHj)  &
          +k(716)*n(idx_O2Hj)  &
          +k(643)*n(idx_NHj)  &
          +k(632)*n(idx_N2Hj)  &
          +k(389)*n(idx_CH4j)  &
          +k(508)*n(idx_H3j)
      pdj(87) =  &
          -k(632)*n(idx_N2Hj)
      pdj(88) =  &
          -k(716)*n(idx_O2Hj)
    elseif(j==30) then
      pdj(3) =  &
          -k(940)*n(idx_O)
      pdj(5) =  &
          +k(894)*n(idx_N)  &
          +k(229)  &
          +k(1001)  &
          +k(854)*n(idx_H)
      pdj(6) =  &
          +k(940)*n(idx_O)  &
          +k(854)*n(idx_H)
      pdj(8) =  &
          +k(1001)  &
          +k(229)  &
          -k(854)*n(idx_H)
      pdj(24) =  &
          -k(894)*n(idx_N)
      pdj(25) =  &
          +k(894)*n(idx_N)
      pdj(30) =  &
          -k(854)*n(idx_H)  &
          -k(1001)  &
          -k(229)  &
          -k(940)*n(idx_O)  &
          -k(894)*n(idx_N)  &
          -k(1118)
      pdj(34) =  &
          +k(940)*n(idx_O)
      pdj(50) =  &
          +k(1118)
    elseif(j==31) then
      pdj(4) =  &
          +k(885)*n(idx_C)
      pdj(7) =  &
          -k(885)*n(idx_C)
      pdj(19) =  &
          +k(885)*n(idx_C)  &
          +k(237)  &
          +k(438)*n(idx_Hj)  &
          +k(1016)
      pdj(25) =  &
          +k(1016)  &
          +k(237)
      pdj(31) =  &
          -k(438)*n(idx_Hj)  &
          -k(237)  &
          -k(1065)  &
          -k(885)*n(idx_C)  &
          -k(1016)
      pdj(36) =  &
          +k(1065)
      pdj(55) =  &
          -k(438)*n(idx_Hj)
      pdj(69) =  &
          +k(438)*n(idx_Hj)
    elseif(j==32) then
      pdj(3) =  &
          +k(1028)  &
          +2.d0*k(900)*n(idx_N)  &
          +k(250)  &
          -k(955)*n(idx_O)
      pdj(6) =  &
          +k(521)*n(idx_H3j)
      pdj(8) =  &
          -k(867)*n(idx_H)
      pdj(10) =  &
          +k(440)*n(idx_Hj)  &
          +k(867)*n(idx_H)  &
          +k(521)*n(idx_H3j)
      pdj(11) =  &
          +k(902)*n(idx_N)  &
          +k(955)*n(idx_O)  &
          +k(713)*n(idx_Oj)
      pdj(12) =  &
          -k(767)*n(idx_CH2)
      pdj(13) =  &
          +k(767)*n(idx_CH2)  &
          +k(791)*n(idx_CH3)
      pdj(17) =  &
          +k(921)*n(idx_NH)  &
          +k(1028)  &
          +2.d0*k(901)*n(idx_N)  &
          +k(250)  &
          +k(827)*n(idx_CN)  &
          +k(955)*n(idx_O)  &
          +k(767)*n(idx_CH2)  &
          +k(833)*n(idx_CO)  &
          +k(867)*n(idx_H)
      pdj(18) =  &
          -k(827)*n(idx_CN)
      pdj(19) =  &
          -k(833)*n(idx_CO)
      pdj(20) =  &
          +k(902)*n(idx_N)  &
          +k(900)*n(idx_N)
      pdj(22) =  &
          -k(791)*n(idx_CH3)
      pdj(24) =  &
          -k(901)*n(idx_N)  &
          -k(900)*n(idx_N)  &
          -k(902)*n(idx_N)
      pdj(25) =  &
          -k(921)*n(idx_NH)
      pdj(27) =  &
          +k(921)*n(idx_NH)  &
          +k(791)*n(idx_CH3)
      pdj(29) =  &
          +k(833)*n(idx_CO)
      pdj(32) =  &
          -k(767)*n(idx_CH2)  &
          -k(440)*n(idx_Hj)  &
          -k(901)*n(idx_N)  &
          -k(900)*n(idx_N)  &
          -k(1113)  &
          -k(791)*n(idx_CH3)  &
          -k(1028)  &
          -k(921)*n(idx_NH)  &
          -k(250)  &
          -k(521)*n(idx_H3j)  &
          -k(833)*n(idx_CO)  &
          -k(867)*n(idx_H)  &
          -k(827)*n(idx_CN)  &
          -k(955)*n(idx_O)  &
          -k(713)*n(idx_Oj)  &
          -k(902)*n(idx_N)
      pdj(34) =  &
          +k(827)*n(idx_CN)
      pdj(47) =  &
          +k(1113)
      pdj(55) =  &
          -k(440)*n(idx_Hj)
      pdj(63) =  &
          +k(713)*n(idx_Oj)  &
          +k(440)*n(idx_Hj)  &
          +k(521)*n(idx_H3j)
      pdj(70) =  &
          -k(713)*n(idx_Oj)
      pdj(81) =  &
          -k(521)*n(idx_H3j)
    elseif(j==33) then
      pdj(2) =  &
          -k(819)*n(idx_CH)  &
          -k(820)*n(idx_CH)
      pdj(3) =  &
          -k(957)*n(idx_O)  &
          +k(871)*n(idx_H)  &
          +k(1035)
      pdj(6) =  &
          +k(872)*n(idx_H)
      pdj(8) =  &
          +k(255)  &
          -k(872)*n(idx_H)  &
          -k(873)*n(idx_H)  &
          +k(1034)  &
          -k(871)*n(idx_H)
      pdj(9) =  &
          +k(871)*n(idx_H)  &
          +k(971)*n(idx_OH)
      pdj(10) =  &
          -k(971)*n(idx_OH)  &
          +k(819)*n(idx_CH)  &
          +k(1035)  &
          +k(957)*n(idx_O)  &
          +k(835)*n(idx_CO)  &
          +2.d0*k(873)*n(idx_H)
      pdj(11) =  &
          +k(255)  &
          +k(796)*n(idx_CH3)  &
          +k(820)*n(idx_CH)  &
          +k(1034)  &
          +k(872)*n(idx_H)  &
          +k(905)*n(idx_N)  &
          +k(971)*n(idx_OH)  &
          +k(884)*n(idx_HCO)  &
          +k(957)*n(idx_O)
      pdj(12) =  &
          +k(820)*n(idx_CH)
      pdj(13) =  &
          +k(884)*n(idx_HCO)
      pdj(14) =  &
          +k(819)*n(idx_CH)  &
          -k(884)*n(idx_HCO)
      pdj(19) =  &
          -k(835)*n(idx_CO)
      pdj(22) =  &
          -k(796)*n(idx_CH3)
      pdj(23) =  &
          +k(796)*n(idx_CH3)
      pdj(24) =  &
          -k(905)*n(idx_N)
      pdj(25) =  &
          +k(905)*n(idx_N)
      pdj(29) =  &
          +k(835)*n(idx_CO)
      pdj(33) =  &
          -k(872)*n(idx_H)  &
          -k(957)*n(idx_O)  &
          -k(905)*n(idx_N)  &
          -k(1123)  &
          -k(971)*n(idx_OH)  &
          -k(796)*n(idx_CH3)  &
          -k(1034)  &
          -k(1035)  &
          -k(819)*n(idx_CH)  &
          -k(873)*n(idx_H)  &
          -k(255)  &
          -k(835)*n(idx_CO)  &
          -k(820)*n(idx_CH)  &
          -k(884)*n(idx_HCO)  &
          -k(871)*n(idx_H)
      pdj(49) =  &
          +k(1123)
    elseif(j==34) then
      pdj(3) =  &
          +k(607)*n(idx_HEj)  &
          -k(958)*n(idx_O)  &
          +k(257)  &
          +k(1036)  &
          -k(959)*n(idx_O)  &
          +k(874)*n(idx_H)
      pdj(5) =  &
          +k(874)*n(idx_H)
      pdj(7) =  &
          -k(757)*n(idx_C)
      pdj(8) =  &
          -k(875)*n(idx_H)  &
          -k(876)*n(idx_H)  &
          -k(874)*n(idx_H)
      pdj(10) =  &
          +k(876)*n(idx_H)
      pdj(11) =  &
          -k(935)*n(idx_O2)  &
          -k(934)*n(idx_O2)  &
          +k(959)*n(idx_O)
      pdj(17) =  &
          +k(958)*n(idx_O)  &
          -k(933)*n(idx_NO)  &
          +k(934)*n(idx_O2)
      pdj(18) =  &
          +k(257)  &
          +k(1036)  &
          +k(757)*n(idx_C)  &
          +k(608)*n(idx_HEj)  &
          +k(959)*n(idx_O)  &
          +k(876)*n(idx_H)  &
          +k(327)*n(idx_Cj)
      pdj(19) =  &
          +k(875)*n(idx_H)  &
          +k(958)*n(idx_O)  &
          +k(935)*n(idx_O2)  &
          +k(757)*n(idx_C)
      pdj(20) =  &
          +k(933)*n(idx_NO)
      pdj(25) =  &
          +k(875)*n(idx_H)
      pdj(26) =  &
          +k(608)*n(idx_HEj)  &
          +k(607)*n(idx_HEj)
      pdj(29) =  &
          +k(933)*n(idx_NO)  &
          +k(934)*n(idx_O2)
      pdj(32) =  &
          +k(935)*n(idx_O2)
      pdj(34) =  &
          -k(327)*n(idx_Cj)  &
          -k(935)*n(idx_O2)  &
          -k(958)*n(idx_O)  &
          -k(757)*n(idx_C)  &
          -k(933)*n(idx_NO)  &
          -k(607)*n(idx_HEj)  &
          -k(876)*n(idx_H)  &
          -k(959)*n(idx_O)  &
          -k(1066)  &
          -k(874)*n(idx_H)  &
          -k(257)  &
          -k(875)*n(idx_H)  &
          -k(934)*n(idx_O2)  &
          -k(608)*n(idx_HEj)  &
          -k(1036)
      pdj(36) =  &
          +k(1066)
      pdj(57) =  &
          -k(327)*n(idx_Cj)
      pdj(64) =  &
          +k(607)*n(idx_HEj)
      pdj(65) =  &
          +k(327)*n(idx_Cj)
      pdj(70) =  &
          +k(608)*n(idx_HEj)
      pdj(78) =  &
          -k(607)*n(idx_HEj)  &
          -k(608)*n(idx_HEj)
    elseif(j==35) then
      pdj(28) =  &
          +k(1140)
      pdj(35) =  &
          -k(1140)
    elseif(j==36) then
      pdj(31) =  &
          +k(1142)
      pdj(36) =  &
          -k(1142)
    elseif(j==37) then
      pdj(13) =  &
          +k(1137)
      pdj(37) =  &
          -k(1137)
    elseif(j==38) then
      pdj(23) =  &
          +k(1127)
      pdj(38) =  &
          -k(1127)
    elseif(j==39) then
      pdj(19) =  &
          +k(1133)
      pdj(39) =  &
          -k(1133)
    elseif(j==40) then
      pdj(9) =  &
          +k(1129)
      pdj(40) =  &
          -k(1129)
    elseif(j==41) then
      pdj(17) =  &
          +k(1136)
      pdj(41) =  &
          -k(1136)
    elseif(j==42) then
      pdj(29) =  &
          +k(1143)
      pdj(42) =  &
          -k(1143)
    elseif(j==43) then
      pdj(20) =  &
          +k(1134)
      pdj(43) =  &
          -k(1134)
    elseif(j==44) then
      pdj(5) =  &
          +k(1131)
      pdj(44) =  &
          -k(1131)
    elseif(j==45) then
      pdj(16) =  &
          +k(1128)
      pdj(45) =  &
          -k(1128)
    elseif(j==46) then
      pdj(11) =  &
          +k(1139)
      pdj(46) =  &
          -k(1139)
    elseif(j==47) then
      pdj(32) =  &
          +k(1144)
      pdj(47) =  &
          -k(1144)
    elseif(j==48) then
      pdj(27) =  &
          +k(1138)
      pdj(48) =  &
          -k(1138)
    elseif(j==49) then
      pdj(33) =  &
          +k(1141)
      pdj(49) =  &
          -k(1141)
    elseif(j==50) then
      pdj(30) =  &
          +k(1135)
      pdj(50) =  &
          -k(1135)
    elseif(j==51) then
      pdj(15) =  &
          +k(1130)
      pdj(51) =  &
          -k(1130)
    elseif(j==52) then
      pdj(4) =  &
          +k(1132)
      pdj(52) =  &
          -k(1132)
    elseif(j==53) then
    elseif(j==54) then
      pdj(1) =  &
          -k(295)*n(idx_E)
      pdj(2) =  &
          -k(408)*n(idx_CH)
      pdj(4) =  &
          -k(559)*n(idx_HNC)
      pdj(5) =  &
          -k(544)*n(idx_HCN)
      pdj(7) =  &
          -k(332)*n(idx_C)
      pdj(8) =  &
          +k(295)*n(idx_E)  &
          +k(1012)  &
          +k(742)*n(idx_OH)
      pdj(9) =  &
          -k(496)*n(idx_H2O)
      pdj(10) =  &
          -k(741)*n(idx_OH)  &
          -k(742)*n(idx_OH)
      pdj(12) =  &
          -k(373)*n(idx_CH2)
      pdj(13) =  &
          -k(550)*n(idx_H2CO)
      pdj(14) =  &
          -k(551)*n(idx_HCO)  &
          +k(131)*n(idx_MG)
      pdj(15) =  &
          -k(131)*n(idx_MG)
      pdj(19) =  &
          +k(741)*n(idx_OH)  &
          +k(550)*n(idx_H2CO)  &
          +k(496)*n(idx_H2O)  &
          +k(408)*n(idx_CH)  &
          +k(332)*n(idx_C)  &
          +k(559)*n(idx_HNC)  &
          +k(544)*n(idx_HCN)  &
          +k(295)*n(idx_E)  &
          +k(373)*n(idx_CH2)  &
          +k(694)*n(idx_NH)  &
          +k(551)*n(idx_HCO)  &
          +k(682)*n(idx_NH2)
      pdj(21) =  &
          -k(682)*n(idx_NH2)
      pdj(25) =  &
          -k(694)*n(idx_NH)
      pdj(37) =  &
          +k(1101)
      pdj(54) =  &
          -k(682)*n(idx_NH2)  &
          -k(694)*n(idx_NH)  &
          -k(332)*n(idx_C)  &
          -k(295)*n(idx_E)  &
          -k(559)*n(idx_HNC)  &
          -k(551)*n(idx_HCO)  &
          -k(373)*n(idx_CH2)  &
          -k(741)*n(idx_OH)  &
          -k(742)*n(idx_OH)  &
          -k(131)*n(idx_MG)  &
          -k(1012)  &
          -k(544)*n(idx_HCN)  &
          -k(408)*n(idx_CH)  &
          -k(496)*n(idx_H2O)  &
          -k(1101)  &
          -k(550)*n(idx_H2CO)
      pdj(58) =  &
          +k(408)*n(idx_CH)
      pdj(59) =  &
          +k(332)*n(idx_C)
      pdj(60) =  &
          +k(551)*n(idx_HCO)
      pdj(61) =  &
          +k(131)*n(idx_MG)
      pdj(62) =  &
          +k(682)*n(idx_NH2)
      pdj(65) =  &
          +k(1012)
      pdj(68) =  &
          +k(741)*n(idx_OH)
      pdj(69) =  &
          +k(694)*n(idx_NH)
      pdj(72) =  &
          +k(373)*n(idx_CH2)
      pdj(82) =  &
          +k(550)*n(idx_H2CO)
      pdj(83) =  &
          +k(496)*n(idx_H2O)
      pdj(84) =  &
          +k(544)*n(idx_HCN)  &
          +k(559)*n(idx_HNC)
      pdj(85) =  &
          +k(742)*n(idx_OH)
    elseif(j==55) then
      pdj(1) =  &
          -k(1058)*n(idx_E)
      pdj(2) =  &
          -k(72)*n(idx_CH)
      pdj(3) =  &
          +k(433)*n(idx_CO2)  &
          -k(83)*n(idx_O)
      pdj(4) =  &
          -k(2)*n(idx_HNC)
      pdj(5) =  &
          +k(2)*n(idx_HNC)  &
          -k(75)*n(idx_HCN)
      pdj(6) =  &
          +k(432)*n(idx_CH4)  &
          +k(436)*n(idx_HCO)  &
          +k(434)*n(idx_H2CO)  &
          +k(439)*n(idx_HNO)  &
          +k(435)*n(idx_H2CO)  &
          +k(430)*n(idx_CH3OH)  &
          +k(428)*n(idx_CH2)  &
          +2.d0*k(431)*n(idx_CH3OH)
      pdj(8) =  &
          -k(1045)*n(idx_H)  &
          +k(79)*n(idx_NH3)  &
          +k(70)*n(idx_CH3)  &
          +k(434)*n(idx_H2CO)  &
          +k(80)*n(idx_NH)  &
          +k(69)*n(idx_CH2)  &
          +k(81)*n(idx_NO)  &
          +k(84)*n(idx_OH)  &
          +k(78)*n(idx_NH2)  &
          +k(83)*n(idx_O)  &
          +k(71)*n(idx_CH4)  &
          +k(73)*n(idx_H2CO)  &
          +k(82)*n(idx_O2)  &
          +k(72)*n(idx_CH)  &
          +k(75)*n(idx_HCN)  &
          +k(74)*n(idx_H2O)  &
          +k(76)*n(idx_HCO)  &
          +k(77)*n(idx_MG)  &
          +k(1058)*n(idx_E)
      pdj(9) =  &
          +k(429)*n(idx_CH3OH)  &
          -k(74)*n(idx_H2O)
      pdj(10) =  &
          +k(440)*n(idx_NO2)  &
          -k(84)*n(idx_OH)
      pdj(11) =  &
          -k(82)*n(idx_O2)
      pdj(12) =  &
          -k(428)*n(idx_CH2)  &
          -k(69)*n(idx_CH2)
      pdj(13) =  &
          -k(434)*n(idx_H2CO)  &
          -k(73)*n(idx_H2CO)  &
          -k(435)*n(idx_H2CO)
      pdj(14) =  &
          -k(437)*n(idx_HCO)  &
          -k(436)*n(idx_HCO)  &
          -k(76)*n(idx_HCO)
      pdj(15) =  &
          -k(77)*n(idx_MG)
      pdj(16) =  &
          -k(79)*n(idx_NH3)
      pdj(17) =  &
          -k(81)*n(idx_NO)
      pdj(19) =  &
          +k(438)*n(idx_HNCO)  &
          +k(437)*n(idx_HCO)
      pdj(21) =  &
          -k(78)*n(idx_NH2)
      pdj(22) =  &
          -k(70)*n(idx_CH3)
      pdj(23) =  &
          -k(432)*n(idx_CH4)  &
          -k(71)*n(idx_CH4)
      pdj(25) =  &
          -k(80)*n(idx_NH)
      pdj(26) =  &
          -k(1046)*n(idx_HE)
      pdj(27) =  &
          -k(439)*n(idx_HNO)
      pdj(28) =  &
          -k(430)*n(idx_CH3OH)  &
          -k(431)*n(idx_CH3OH)  &
          -k(429)*n(idx_CH3OH)
      pdj(29) =  &
          -k(433)*n(idx_CO2)
      pdj(31) =  &
          -k(438)*n(idx_HNCO)
      pdj(32) =  &
          -k(440)*n(idx_NO2)
      pdj(54) =  &
          +k(433)*n(idx_CO2)  &
          +k(76)*n(idx_HCO)  &
          +k(435)*n(idx_H2CO)  &
          +k(431)*n(idx_CH3OH)
      pdj(55) =  &
          -k(75)*n(idx_HCN)  &
          -k(434)*n(idx_H2CO)  &
          -k(71)*n(idx_CH4)  &
          -k(69)*n(idx_CH2)  &
          -k(433)*n(idx_CO2)  &
          -k(428)*n(idx_CH2)  &
          -k(440)*n(idx_NO2)  &
          -k(430)*n(idx_CH3OH)  &
          -k(74)*n(idx_H2O)  &
          -k(439)*n(idx_HNO)  &
          -k(2)*n(idx_HNC)  &
          -k(80)*n(idx_NH)  &
          -k(83)*n(idx_O)  &
          -k(1058)*n(idx_E)  &
          -k(70)*n(idx_CH3)  &
          -k(84)*n(idx_OH)  &
          -k(73)*n(idx_H2CO)  &
          -k(437)*n(idx_HCO)  &
          -k(431)*n(idx_CH3OH)  &
          -k(429)*n(idx_CH3OH)  &
          -k(81)*n(idx_NO)  &
          -k(435)*n(idx_H2CO)  &
          -k(72)*n(idx_CH)  &
          -k(78)*n(idx_NH2)  &
          -k(1046)*n(idx_HE)  &
          -k(432)*n(idx_CH4)  &
          -k(1045)*n(idx_H)  &
          -k(77)*n(idx_MG)  &
          -k(76)*n(idx_HCO)  &
          -k(79)*n(idx_NH3)  &
          +k(2)*n(idx_HNC)  &
          -k(436)*n(idx_HCO)  &
          -k(82)*n(idx_O2)  &
          -k(438)*n(idx_HNCO)
      pdj(58) =  &
          +k(69)*n(idx_CH2)
      pdj(59) =  &
          +k(72)*n(idx_CH)  &
          +k(428)*n(idx_CH2)
      pdj(60) =  &
          +k(73)*n(idx_H2CO)
      pdj(61) =  &
          +k(77)*n(idx_MG)
      pdj(62) =  &
          +k(79)*n(idx_NH3)
      pdj(63) =  &
          +k(439)*n(idx_HNO)  &
          +k(81)*n(idx_NO)  &
          +k(440)*n(idx_NO2)
      pdj(65) =  &
          +k(436)*n(idx_HCO)  &
          +k(434)*n(idx_H2CO)
      pdj(67) =  &
          +k(82)*n(idx_O2)
      pdj(68) =  &
          +k(74)*n(idx_H2O)
      pdj(69) =  &
          +k(78)*n(idx_NH2)  &
          +k(438)*n(idx_HNCO)
      pdj(70) =  &
          +k(83)*n(idx_O)
      pdj(71) =  &
          +k(84)*n(idx_OH)
      pdj(72) =  &
          +k(429)*n(idx_CH3OH)  &
          +k(432)*n(idx_CH4)  &
          +k(70)*n(idx_CH3)
      pdj(73) =  &
          +k(71)*n(idx_CH4)
      pdj(75) =  &
          +k(75)*n(idx_HCN)
      pdj(76) =  &
          +k(80)*n(idx_NH)
      pdj(77) =  &
          +k(437)*n(idx_HCO)  &
          +k(1045)*n(idx_H)
      pdj(82) =  &
          +k(430)*n(idx_CH3OH)
      pdj(86) =  &
          +k(1046)*n(idx_HE)
    elseif(j==56) then
      pdj(1) =  &
          -k(300)*n(idx_E)
      pdj(6) =  &
          -k(6)*n(idx_H2)  &
          +k(6)*n(idx_H2)
      pdj(8) =  &
          +k(300)*n(idx_E)
      pdj(19) =  &
          +k(300)*n(idx_E)
      pdj(37) =  &
          +k(1067)
      pdj(54) =  &
          +k(6)*n(idx_H2)
      pdj(56) =  &
          -k(300)*n(idx_E)  &
          -k(6)*n(idx_H2)  &
          -k(1067)
    elseif(j==57) then
      pdj(1) =  &
          -k(1056)*n(idx_E)
      pdj(2) =  &
          +k(318)*n(idx_H2CO)  &
          -k(16)*n(idx_CH)  &
          +k(314)*n(idx_CH3OH)
      pdj(3) =  &
          -k(1041)*n(idx_O)  &
          +k(325)*n(idx_O2)
      pdj(6) =  &
          -k(1047)*n(idx_H2)  &
          +k(323)*n(idx_NH3)  &
          -k(460)*n(idx_H2)
      pdj(7) =  &
          +k(20)*n(idx_NH3)  &
          +k(21)*n(idx_NO)  &
          +k(19)*n(idx_MG)  &
          +k(17)*n(idx_H2CO)  &
          +k(18)*n(idx_HCO)  &
          +k(15)*n(idx_CH2)  &
          +k(1056)*n(idx_E)  &
          +k(16)*n(idx_CH)
      pdj(8) =  &
          +k(320)*n(idx_H2O)  &
          +k(324)*n(idx_NH)  &
          +k(328)*n(idx_OH)  &
          +k(460)*n(idx_H2)  &
          +k(319)*n(idx_H2O)  &
          +k(322)*n(idx_NH2)  &
          -k(1050)*n(idx_H)
      pdj(9) =  &
          -k(319)*n(idx_H2O)  &
          -k(320)*n(idx_H2O)
      pdj(10) =  &
          -k(328)*n(idx_OH)
      pdj(11) =  &
          -k(325)*n(idx_O2)  &
          -k(326)*n(idx_O2)
      pdj(12) =  &
          -k(15)*n(idx_CH2)
      pdj(13) =  &
          -k(17)*n(idx_H2CO)  &
          -k(317)*n(idx_H2CO)  &
          -k(318)*n(idx_H2CO)
      pdj(14) =  &
          +k(315)*n(idx_CH3OH)  &
          -k(321)*n(idx_HCO)  &
          -k(18)*n(idx_HCO)
      pdj(15) =  &
          -k(19)*n(idx_MG)
      pdj(16) =  &
          -k(323)*n(idx_NH3)  &
          -k(20)*n(idx_NH3)
      pdj(17) =  &
          -k(21)*n(idx_NO)
      pdj(18) =  &
          +k(327)*n(idx_OCN)
      pdj(19) =  &
          +k(317)*n(idx_H2CO)  &
          +k(321)*n(idx_HCO)  &
          +k(316)*n(idx_CO2)  &
          +k(326)*n(idx_O2)
      pdj(21) =  &
          -k(322)*n(idx_NH2)
      pdj(24) =  &
          -k(1040)*n(idx_N)
      pdj(25) =  &
          -k(324)*n(idx_NH)
      pdj(28) =  &
          -k(315)*n(idx_CH3OH)  &
          -k(314)*n(idx_CH3OH)
      pdj(29) =  &
          -k(316)*n(idx_CO2)
      pdj(34) =  &
          -k(327)*n(idx_OCN)
      pdj(38) =  &
          +k(1084)
      pdj(54) =  &
          +k(318)*n(idx_H2CO)  &
          +k(319)*n(idx_H2O)  &
          +k(18)*n(idx_HCO)
      pdj(56) =  &
          +k(320)*n(idx_H2O)
      pdj(57) =  &
          -k(1040)*n(idx_N)  &
          -k(1056)*n(idx_E)  &
          -k(15)*n(idx_CH2)  &
          -k(17)*n(idx_H2CO)  &
          -k(21)*n(idx_NO)  &
          -k(1041)*n(idx_O)  &
          -k(317)*n(idx_H2CO)  &
          -k(16)*n(idx_CH)  &
          -k(314)*n(idx_CH3OH)  &
          -k(321)*n(idx_HCO)  &
          -k(19)*n(idx_MG)  &
          -k(316)*n(idx_CO2)  &
          -k(1047)*n(idx_H2)  &
          -k(18)*n(idx_HCO)  &
          -k(1050)*n(idx_H)  &
          -k(1084)  &
          -k(320)*n(idx_H2O)  &
          -k(328)*n(idx_OH)  &
          -k(315)*n(idx_CH3OH)  &
          -k(460)*n(idx_H2)  &
          -k(325)*n(idx_O2)  &
          -k(20)*n(idx_NH3)  &
          -k(322)*n(idx_NH2)  &
          -k(324)*n(idx_NH)  &
          -k(319)*n(idx_H2O)  &
          -k(323)*n(idx_NH3)  &
          -k(326)*n(idx_O2)  &
          -k(327)*n(idx_OCN)  &
          -k(318)*n(idx_H2CO)
      pdj(58) =  &
          +k(15)*n(idx_CH2)  &
          +k(317)*n(idx_H2CO)  &
          +k(1047)*n(idx_H2)
      pdj(59) =  &
          +k(460)*n(idx_H2)  &
          +k(321)*n(idx_HCO)  &
          +k(1050)*n(idx_H)  &
          +k(16)*n(idx_CH)
      pdj(60) =  &
          +k(17)*n(idx_H2CO)
      pdj(61) =  &
          +k(19)*n(idx_MG)
      pdj(62) =  &
          +k(20)*n(idx_NH3)
      pdj(63) =  &
          +k(21)*n(idx_NO)
      pdj(64) =  &
          +k(324)*n(idx_NH)  &
          +k(1040)*n(idx_N)
      pdj(65) =  &
          +k(327)*n(idx_OCN)  &
          +k(325)*n(idx_O2)  &
          +k(316)*n(idx_CO2)  &
          +k(1041)*n(idx_O)  &
          +k(328)*n(idx_OH)
      pdj(70) =  &
          +k(326)*n(idx_O2)
      pdj(72) =  &
          +k(315)*n(idx_CH3OH)
      pdj(75) =  &
          +k(323)*n(idx_NH3)  &
          +k(322)*n(idx_NH2)
      pdj(82) =  &
          +k(314)*n(idx_CH3OH)
    elseif(j==58) then
      pdj(1) =  &
          -k(262)*n(idx_E)  &
          -k(261)*n(idx_E)  &
          -k(260)*n(idx_E)
      pdj(2) =  &
          +k(262)*n(idx_E)  &
          +k(977)
      pdj(3) =  &
          -k(365)*n(idx_O)
      pdj(6) =  &
          +k(531)*n(idx_H)  &
          +k(975)  &
          -k(462)*n(idx_H2)  &
          +k(260)*n(idx_E)
      pdj(7) =  &
          +k(261)*n(idx_E)  &
          +k(260)*n(idx_E)
      pdj(8) =  &
          +k(262)*n(idx_E)  &
          +k(976)  &
          +k(365)*n(idx_O)  &
          +k(462)*n(idx_H2)  &
          +k(362)*n(idx_H2O)  &
          -k(531)*n(idx_H)  &
          +k(634)*n(idx_N)  &
          +2.d0*k(261)*n(idx_E)
      pdj(9) =  &
          -k(362)*n(idx_H2O)
      pdj(10) =  &
          +k(364)*n(idx_O2)
      pdj(11) =  &
          -k(364)*n(idx_O2)
      pdj(12) =  &
          +k(30)*n(idx_NO)
      pdj(13) =  &
          -k(361)*n(idx_H2CO)
      pdj(14) =  &
          -k(363)*n(idx_HCO)
      pdj(17) =  &
          -k(30)*n(idx_NO)
      pdj(19) =  &
          +k(363)*n(idx_HCO)  &
          +k(360)*n(idx_CO2)
      pdj(22) =  &
          +k(361)*n(idx_H2CO)
      pdj(24) =  &
          -k(634)*n(idx_N)
      pdj(29) =  &
          -k(360)*n(idx_CO2)
      pdj(38) =  &
          +k(1098)
      pdj(54) =  &
          +k(361)*n(idx_H2CO)  &
          +k(364)*n(idx_O2)  &
          +k(365)*n(idx_O)
      pdj(55) =  &
          +k(977)
      pdj(57) =  &
          +k(975)
      pdj(58) =  &
          -k(365)*n(idx_O)  &
          -k(363)*n(idx_HCO)  &
          -k(364)*n(idx_O2)  &
          -k(976)  &
          -k(634)*n(idx_N)  &
          -k(977)  &
          -k(261)*n(idx_E)  &
          -k(260)*n(idx_E)  &
          -k(462)*n(idx_H2)  &
          -k(30)*n(idx_NO)  &
          -k(975)  &
          -k(262)*n(idx_E)  &
          -k(531)*n(idx_H)  &
          -k(360)*n(idx_CO2)  &
          -k(1098)  &
          -k(362)*n(idx_H2O)  &
          -k(361)*n(idx_H2CO)
      pdj(59) =  &
          +k(531)*n(idx_H)  &
          +k(976)
      pdj(60) =  &
          +k(360)*n(idx_CO2)
      pdj(63) =  &
          +k(30)*n(idx_NO)
      pdj(72) =  &
          +k(462)*n(idx_H2)  &
          +k(363)*n(idx_HCO)
      pdj(75) =  &
          +k(634)*n(idx_N)
      pdj(82) =  &
          +k(362)*n(idx_H2O)
    elseif(j==59) then
      pdj(1) =  &
          -k(259)*n(idx_E)
      pdj(2) =  &
          +k(28)*n(idx_NH3)  &
          +k(26)*n(idx_HCO)  &
          +k(29)*n(idx_NO)  &
          +k(27)*n(idx_MG)
      pdj(3) =  &
          +k(356)*n(idx_O2)  &
          -k(358)*n(idx_O)
      pdj(4) =  &
          -k(351)*n(idx_HNC)
      pdj(5) =  &
          -k(349)*n(idx_HCN)
      pdj(6) =  &
          +k(359)*n(idx_OH)  &
          +k(353)*n(idx_NH2)  &
          -k(461)*n(idx_H2)  &
          +k(354)*n(idx_NH)  &
          +k(348)*n(idx_H2O)  &
          +k(530)*n(idx_H)
      pdj(7) =  &
          +k(344)*n(idx_H2CO)  &
          +k(349)*n(idx_HCN)  &
          +k(351)*n(idx_HNC)  &
          +k(347)*n(idx_H2O)  &
          +k(259)*n(idx_E)  &
          +k(974)
      pdj(8) =  &
          +k(352)*n(idx_N)  &
          -k(530)*n(idx_H)  &
          +k(358)*n(idx_O)  &
          +k(346)*n(idx_H2O)  &
          +k(216)  &
          +k(259)*n(idx_E)  &
          +k(461)*n(idx_H2)
      pdj(9) =  &
          -k(348)*n(idx_H2O)  &
          -k(346)*n(idx_H2O)  &
          -k(347)*n(idx_H2O)
      pdj(10) =  &
          +k(355)*n(idx_O2)  &
          -k(359)*n(idx_OH)
      pdj(11) =  &
          -k(357)*n(idx_O2)  &
          -k(356)*n(idx_O2)  &
          -k(355)*n(idx_O2)
      pdj(12) =  &
          +k(341)*n(idx_CH3OH)  &
          +k(345)*n(idx_H2CO)
      pdj(13) =  &
          -k(345)*n(idx_H2CO)  &
          -k(343)*n(idx_H2CO)  &
          +k(340)*n(idx_CH3OH)  &
          -k(344)*n(idx_H2CO)
      pdj(14) =  &
          -k(350)*n(idx_HCO)  &
          -k(26)*n(idx_HCO)  &
          +k(357)*n(idx_O2)
      pdj(15) =  &
          -k(27)*n(idx_MG)
      pdj(16) =  &
          -k(28)*n(idx_NH3)
      pdj(17) =  &
          -k(29)*n(idx_NO)
      pdj(19) =  &
          +k(343)*n(idx_H2CO)  &
          +k(350)*n(idx_HCO)  &
          +k(342)*n(idx_CO2)
      pdj(21) =  &
          -k(353)*n(idx_NH2)
      pdj(24) =  &
          -k(352)*n(idx_N)
      pdj(25) =  &
          -k(354)*n(idx_NH)
      pdj(28) =  &
          -k(340)*n(idx_CH3OH)  &
          -k(341)*n(idx_CH3OH)
      pdj(29) =  &
          -k(342)*n(idx_CO2)
      pdj(38) =  &
          +k(1092)
      pdj(54) =  &
          +k(356)*n(idx_O2)  &
          +k(348)*n(idx_H2O)  &
          +k(26)*n(idx_HCO)  &
          +k(345)*n(idx_H2CO)  &
          +k(342)*n(idx_CO2)
      pdj(55) =  &
          +k(974)
      pdj(57) =  &
          +k(216)  &
          +k(530)*n(idx_H)
      pdj(58) =  &
          +k(461)*n(idx_H2)  &
          +k(350)*n(idx_HCO)
      pdj(59) =  &
          -k(216)  &
          -k(530)*n(idx_H)  &
          -k(357)*n(idx_O2)  &
          -k(461)*n(idx_H2)  &
          -k(29)*n(idx_NO)  &
          -k(358)*n(idx_O)  &
          -k(345)*n(idx_H2CO)  &
          -k(350)*n(idx_HCO)  &
          -k(353)*n(idx_NH2)  &
          -k(341)*n(idx_CH3OH)  &
          -k(28)*n(idx_NH3)  &
          -k(359)*n(idx_OH)  &
          -k(1092)  &
          -k(27)*n(idx_MG)  &
          -k(340)*n(idx_CH3OH)  &
          -k(259)*n(idx_E)  &
          -k(355)*n(idx_O2)  &
          -k(351)*n(idx_HNC)  &
          -k(349)*n(idx_HCN)  &
          -k(26)*n(idx_HCO)  &
          -k(346)*n(idx_H2O)  &
          -k(344)*n(idx_H2CO)  &
          -k(348)*n(idx_H2O)  &
          -k(342)*n(idx_CO2)  &
          -k(343)*n(idx_H2CO)  &
          -k(354)*n(idx_NH)  &
          -k(356)*n(idx_O2)  &
          -k(974)  &
          -k(347)*n(idx_H2O)  &
          -k(352)*n(idx_N)
      pdj(60) =  &
          +k(346)*n(idx_H2O)
      pdj(61) =  &
          +k(27)*n(idx_MG)
      pdj(62) =  &
          +k(28)*n(idx_NH3)
      pdj(63) =  &
          +k(29)*n(idx_NO)
      pdj(64) =  &
          +k(352)*n(idx_N)  &
          +k(354)*n(idx_NH)
      pdj(65) =  &
          +k(355)*n(idx_O2)  &
          +k(359)*n(idx_OH)  &
          +k(358)*n(idx_O)
      pdj(70) =  &
          +k(357)*n(idx_O2)
      pdj(72) =  &
          +k(343)*n(idx_H2CO)  &
          +k(340)*n(idx_CH3OH)
      pdj(75) =  &
          +k(353)*n(idx_NH2)
      pdj(82) =  &
          +k(344)*n(idx_H2CO)  &
          +k(341)*n(idx_CH3OH)
      pdj(83) =  &
          +k(347)*n(idx_H2O)
      pdj(84) =  &
          +k(351)*n(idx_HNC)  &
          +k(349)*n(idx_HCN)
    elseif(j==60) then
      pdj(1) =  &
          -k(1059)*n(idx_E)  &
          -k(271)*n(idx_E)  &
          -k(272)*n(idx_E)  &
          -k(273)*n(idx_E)  &
          -k(274)*n(idx_E)
      pdj(2) =  &
          -k(49)*n(idx_CH)  &
          -k(401)*n(idx_CH)
      pdj(3) =  &
          +k(271)*n(idx_E)
      pdj(4) =  &
          -k(557)*n(idx_HNC)
      pdj(5) =  &
          -k(542)*n(idx_HCN)
      pdj(6) =  &
          +k(272)*n(idx_E)
      pdj(8) =  &
          +k(274)*n(idx_E)  &
          +2.d0*k(273)*n(idx_E)
      pdj(9) =  &
          -k(493)*n(idx_H2O)
      pdj(11) =  &
          -k(479)*n(idx_O2)
      pdj(12) =  &
          +k(271)*n(idx_E)  &
          -k(33)*n(idx_CH2)  &
          -k(367)*n(idx_CH2)
      pdj(13) =  &
          +k(33)*n(idx_CH2)  &
          +k(130)*n(idx_MG)  &
          +k(120)*n(idx_HCO)  &
          +k(182)*n(idx_NO)  &
          -k(478)*n(idx_H2CO)  &
          +k(49)*n(idx_CH)  &
          +k(1059)*n(idx_E)  &
          +k(173)*n(idx_NH3)
      pdj(14) =  &
          +k(274)*n(idx_E)  &
          +k(675)*n(idx_NH2)  &
          -k(120)*n(idx_HCO)  &
          +k(367)*n(idx_CH2)  &
          +k(542)*n(idx_HCN)  &
          +k(493)*n(idx_H2O)  &
          +k(478)*n(idx_H2CO)  &
          +k(401)*n(idx_CH)  &
          +k(557)*n(idx_HNC)  &
          -k(552)*n(idx_HCO)
      pdj(15) =  &
          -k(130)*n(idx_MG)
      pdj(16) =  &
          -k(173)*n(idx_NH3)
      pdj(17) =  &
          -k(182)*n(idx_NO)
      pdj(19) =  &
          +k(272)*n(idx_E)  &
          +k(273)*n(idx_E)  &
          +k(552)*n(idx_HCO)
      pdj(21) =  &
          -k(675)*n(idx_NH2)
      pdj(22) =  &
          +k(394)*n(idx_CH4)
      pdj(23) =  &
          -k(394)*n(idx_CH4)
      pdj(24) =  &
          +k(691)*n(idx_NH)
      pdj(25) =  &
          -k(691)*n(idx_NH)
      pdj(33) =  &
          +k(479)*n(idx_O2)
      pdj(37) =  &
          +k(1104)
      pdj(54) =  &
          +k(479)*n(idx_O2)  &
          +k(120)*n(idx_HCO)
      pdj(58) =  &
          +k(33)*n(idx_CH2)  &
          +k(401)*n(idx_CH)
      pdj(59) =  &
          +k(49)*n(idx_CH)
      pdj(60) =  &
          -k(273)*n(idx_E)  &
          -k(130)*n(idx_MG)  &
          -k(675)*n(idx_NH2)  &
          -k(401)*n(idx_CH)  &
          -k(182)*n(idx_NO)  &
          -k(120)*n(idx_HCO)  &
          -k(1104)  &
          -k(173)*n(idx_NH3)  &
          -k(493)*n(idx_H2O)  &
          -k(274)*n(idx_E)  &
          -k(272)*n(idx_E)  &
          -k(394)*n(idx_CH4)  &
          -k(691)*n(idx_NH)  &
          -k(542)*n(idx_HCN)  &
          -k(271)*n(idx_E)  &
          -k(478)*n(idx_H2CO)  &
          -k(33)*n(idx_CH2)  &
          -k(552)*n(idx_HCO)  &
          -k(49)*n(idx_CH)  &
          -k(479)*n(idx_O2)  &
          -k(1059)*n(idx_E)  &
          -k(557)*n(idx_HNC)  &
          -k(367)*n(idx_CH2)
      pdj(61) =  &
          +k(130)*n(idx_MG)
      pdj(62) =  &
          +k(173)*n(idx_NH3)  &
          +k(675)*n(idx_NH2)
      pdj(63) =  &
          +k(182)*n(idx_NO)
      pdj(72) =  &
          +k(367)*n(idx_CH2)
      pdj(82) =  &
          +k(394)*n(idx_CH4)  &
          +k(552)*n(idx_HCO)  &
          +k(691)*n(idx_NH)  &
          +k(478)*n(idx_H2CO)
      pdj(83) =  &
          +k(493)*n(idx_H2O)
      pdj(84) =  &
          +k(557)*n(idx_HNC)  &
          +k(542)*n(idx_HCN)
    elseif(j==61) then
      pdj(1) =  &
          -k(1061)*n(idx_E)
      pdj(15) =  &
          +k(1061)*n(idx_E)
      pdj(51) =  &
          +k(1119)
      pdj(61) =  &
          -k(1061)*n(idx_E)  &
          -k(1119)
    elseif(j==62) then
      pdj(1) =  &
          -k(308)*n(idx_E)  &
          -k(309)*n(idx_E)
      pdj(3) =  &
          -k(723)*n(idx_O)
      pdj(6) =  &
          +k(723)*n(idx_O)
      pdj(8) =  &
          +k(308)*n(idx_E)  &
          +2.d0*k(309)*n(idx_E)
      pdj(12) =  &
          -k(378)*n(idx_CH2)
      pdj(14) =  &
          -k(169)*n(idx_HCO)
      pdj(15) =  &
          -k(170)*n(idx_MG)
      pdj(16) =  &
          +k(170)*n(idx_MG)  &
          +k(169)*n(idx_HCO)  &
          +k(171)*n(idx_NO)
      pdj(17) =  &
          -k(171)*n(idx_NO)
      pdj(21) =  &
          +k(308)*n(idx_E)  &
          +k(378)*n(idx_CH2)
      pdj(25) =  &
          +k(309)*n(idx_E)
      pdj(45) =  &
          +k(1103)
      pdj(54) =  &
          +k(169)*n(idx_HCO)
      pdj(61) =  &
          +k(170)*n(idx_MG)
      pdj(62) =  &
          -k(308)*n(idx_E)  &
          -k(171)*n(idx_NO)  &
          -k(723)*n(idx_O)  &
          -k(309)*n(idx_E)  &
          -k(170)*n(idx_MG)  &
          -k(1103)  &
          -k(169)*n(idx_HCO)  &
          -k(378)*n(idx_CH2)
      pdj(63) =  &
          +k(171)*n(idx_NO)
      pdj(72) =  &
          +k(378)*n(idx_CH2)
      pdj(79) =  &
          +k(723)*n(idx_O)
    elseif(j==63) then
      pdj(1) =  &
          -k(310)*n(idx_E)
      pdj(3) =  &
          +k(310)*n(idx_E)
      pdj(15) =  &
          -k(133)*n(idx_MG)
      pdj(17) =  &
          +k(133)*n(idx_MG)
      pdj(24) =  &
          +k(310)*n(idx_E)
      pdj(41) =  &
          +k(1097)
      pdj(61) =  &
          +k(133)*n(idx_MG)
      pdj(63) =  &
          -k(1097)  &
          -k(133)*n(idx_MG)  &
          -k(310)*n(idx_E)
    elseif(j==64) then
      pdj(1) =  &
          -k(268)*n(idx_E)
      pdj(2) =  &
          -k(47)*n(idx_CH)
      pdj(3) =  &
          -k(194)*n(idx_O)
      pdj(5) =  &
          +k(418)*n(idx_H2CO)  &
          -k(59)*n(idx_HCN)
      pdj(6) =  &
          -k(463)*n(idx_H2)
      pdj(7) =  &
          +k(635)*n(idx_N)  &
          -k(22)*n(idx_C)  &
          +k(268)*n(idx_E)
      pdj(8) =  &
          +k(463)*n(idx_H2)  &
          -k(110)*n(idx_H)
      pdj(9) =  &
          -k(490)*n(idx_H2O)  &
          -k(491)*n(idx_H2O)
      pdj(10) =  &
          -k(203)*n(idx_OH)  &
          +k(490)*n(idx_H2O)
      pdj(11) =  &
          -k(62)*n(idx_O2)  &
          -k(420)*n(idx_O2)
      pdj(12) =  &
          -k(31)*n(idx_CH2)
      pdj(13) =  &
          -k(418)*n(idx_H2CO)  &
          -k(58)*n(idx_H2CO)
      pdj(14) =  &
          -k(60)*n(idx_HCO)  &
          -k(419)*n(idx_HCO)
      pdj(17) =  &
          -k(61)*n(idx_NO)
      pdj(18) =  &
          +k(47)*n(idx_CH)  &
          +k(22)*n(idx_C)  &
          +k(62)*n(idx_O2)  &
          +k(58)*n(idx_H2CO)  &
          +k(57)*n(idx_CO)  &
          +k(178)*n(idx_NH)  &
          +k(110)*n(idx_H)  &
          +k(59)*n(idx_HCN)  &
          +k(194)*n(idx_O)  &
          +k(61)*n(idx_NO)  &
          +k(31)*n(idx_CH2)  &
          +k(203)*n(idx_OH)  &
          +k(163)*n(idx_NH2)  &
          +k(60)*n(idx_HCO)
      pdj(19) =  &
          -k(57)*n(idx_CO)  &
          +k(420)*n(idx_O2)  &
          +k(419)*n(idx_HCO)
      pdj(21) =  &
          -k(163)*n(idx_NH2)
      pdj(24) =  &
          -k(635)*n(idx_N)  &
          +k(268)*n(idx_E)
      pdj(25) =  &
          -k(178)*n(idx_NH)  &
          +k(491)*n(idx_H2O)
      pdj(44) =  &
          +k(1096)
      pdj(54) =  &
          +k(418)*n(idx_H2CO)  &
          +k(491)*n(idx_H2O)  &
          +k(60)*n(idx_HCO)
      pdj(55) =  &
          +k(110)*n(idx_H)
      pdj(57) =  &
          +k(22)*n(idx_C)
      pdj(58) =  &
          +k(31)*n(idx_CH2)
      pdj(59) =  &
          +k(47)*n(idx_CH)
      pdj(60) =  &
          +k(58)*n(idx_H2CO)
      pdj(63) =  &
          +k(61)*n(idx_NO)  &
          +k(420)*n(idx_O2)
      pdj(64) =  &
          -k(62)*n(idx_O2)  &
          -k(31)*n(idx_CH2)  &
          -k(194)*n(idx_O)  &
          -k(635)*n(idx_N)  &
          -k(60)*n(idx_HCO)  &
          -k(59)*n(idx_HCN)  &
          -k(47)*n(idx_CH)  &
          -k(419)*n(idx_HCO)  &
          -k(22)*n(idx_C)  &
          -k(491)*n(idx_H2O)  &
          -k(463)*n(idx_H2)  &
          -k(61)*n(idx_NO)  &
          -k(58)*n(idx_H2CO)  &
          -k(418)*n(idx_H2CO)  &
          -k(490)*n(idx_H2O)  &
          -k(178)*n(idx_NH)  &
          -k(1096)  &
          -k(268)*n(idx_E)  &
          -k(420)*n(idx_O2)  &
          -k(57)*n(idx_CO)  &
          -k(163)*n(idx_NH2)  &
          -k(110)*n(idx_H)  &
          -k(203)*n(idx_OH)
      pdj(65) =  &
          +k(57)*n(idx_CO)
      pdj(66) =  &
          +k(635)*n(idx_N)
      pdj(67) =  &
          +k(62)*n(idx_O2)
      pdj(69) =  &
          +k(163)*n(idx_NH2)
      pdj(70) =  &
          +k(194)*n(idx_O)
      pdj(71) =  &
          +k(203)*n(idx_OH)
      pdj(75) =  &
          +k(59)*n(idx_HCN)  &
          +k(463)*n(idx_H2)  &
          +k(419)*n(idx_HCO)  &
          +k(490)*n(idx_H2O)
      pdj(76) =  &
          +k(178)*n(idx_NH)
    elseif(j==65) then
      pdj(1) =  &
          -k(269)*n(idx_E)
      pdj(2) =  &
          -k(400)*n(idx_CH)  &
          +k(366)*n(idx_CH2)  &
          -k(48)*n(idx_CH)
      pdj(3) =  &
          -k(195)*n(idx_O)  &
          +k(997)  &
          +k(738)*n(idx_OH)  &
          +k(269)*n(idx_E)
      pdj(5) =  &
          -k(118)*n(idx_HCN)
      pdj(6) =  &
          -k(464)*n(idx_H2)  &
          -k(465)*n(idx_H2)
      pdj(7) =  &
          +k(400)*n(idx_CH)  &
          -k(23)*n(idx_C)  &
          +k(269)*n(idx_E)
      pdj(8) =  &
          +k(465)*n(idx_H2)  &
          -k(111)*n(idx_H)  &
          +k(464)*n(idx_H2)
      pdj(9) =  &
          -k(107)*n(idx_H2O)  &
          -k(492)*n(idx_H2O)
      pdj(10) =  &
          -k(204)*n(idx_OH)  &
          +k(492)*n(idx_H2O)  &
          -k(738)*n(idx_OH)
      pdj(11) =  &
          -k(67)*n(idx_O2)
      pdj(12) =  &
          -k(366)*n(idx_CH2)  &
          -k(32)*n(idx_CH2)
      pdj(13) =  &
          -k(423)*n(idx_H2CO)  &
          -k(64)*n(idx_H2CO)
      pdj(14) =  &
          -k(65)*n(idx_HCO)  &
          +k(423)*n(idx_H2CO)
      pdj(16) =  &
          -k(172)*n(idx_NH3)  &
          -k(687)*n(idx_NH3)
      pdj(17) =  &
          -k(66)*n(idx_NO)
      pdj(19) =  &
          +k(46)*n(idx_CH4)  &
          +k(111)*n(idx_H)  &
          +k(118)*n(idx_HCN)  &
          +k(179)*n(idx_NH)  &
          +k(23)*n(idx_C)  &
          +k(164)*n(idx_NH2)  &
          +k(172)*n(idx_NH3)  &
          +k(67)*n(idx_O2)  &
          +k(48)*n(idx_CH)  &
          +k(195)*n(idx_O)  &
          +k(64)*n(idx_H2CO)  &
          +k(107)*n(idx_H2O)  &
          +k(65)*n(idx_HCO)  &
          +k(32)*n(idx_CH2)  &
          +k(204)*n(idx_OH)  &
          +k(66)*n(idx_NO)
      pdj(21) =  &
          -k(674)*n(idx_NH2)  &
          -k(164)*n(idx_NH2)  &
          +k(687)*n(idx_NH3)
      pdj(22) =  &
          +k(393)*n(idx_CH4)
      pdj(23) =  &
          -k(393)*n(idx_CH4)  &
          -k(46)*n(idx_CH4)
      pdj(24) =  &
          +k(690)*n(idx_NH)
      pdj(25) =  &
          -k(179)*n(idx_NH)  &
          +k(674)*n(idx_NH2)  &
          -k(690)*n(idx_NH)
      pdj(39) =  &
          +k(1095)
      pdj(54) =  &
          +k(400)*n(idx_CH)  &
          +k(674)*n(idx_NH2)  &
          +k(393)*n(idx_CH4)  &
          +k(366)*n(idx_CH2)  &
          +k(423)*n(idx_H2CO)  &
          +k(492)*n(idx_H2O)  &
          +k(738)*n(idx_OH)  &
          +k(464)*n(idx_H2)  &
          +k(687)*n(idx_NH3)  &
          +k(690)*n(idx_NH)  &
          +k(65)*n(idx_HCO)
      pdj(55) =  &
          +k(111)*n(idx_H)
      pdj(56) =  &
          +k(465)*n(idx_H2)
      pdj(57) =  &
          +k(997)  &
          +k(23)*n(idx_C)
      pdj(58) =  &
          +k(32)*n(idx_CH2)
      pdj(59) =  &
          +k(48)*n(idx_CH)
      pdj(60) =  &
          +k(64)*n(idx_H2CO)
      pdj(62) =  &
          +k(172)*n(idx_NH3)
      pdj(63) =  &
          +k(66)*n(idx_NO)
      pdj(65) =  &
          -k(465)*n(idx_H2)  &
          -k(269)*n(idx_E)  &
          -k(464)*n(idx_H2)  &
          -k(400)*n(idx_CH)  &
          -k(172)*n(idx_NH3)  &
          -k(1095)  &
          -k(64)*n(idx_H2CO)  &
          -k(204)*n(idx_OH)  &
          -k(366)*n(idx_CH2)  &
          -k(48)*n(idx_CH)  &
          -k(65)*n(idx_HCO)  &
          -k(423)*n(idx_H2CO)  &
          -k(107)*n(idx_H2O)  &
          -k(195)*n(idx_O)  &
          -k(179)*n(idx_NH)  &
          -k(111)*n(idx_H)  &
          -k(674)*n(idx_NH2)  &
          -k(492)*n(idx_H2O)  &
          -k(393)*n(idx_CH4)  &
          -k(690)*n(idx_NH)  &
          -k(46)*n(idx_CH4)  &
          -k(23)*n(idx_C)  &
          -k(66)*n(idx_NO)  &
          -k(687)*n(idx_NH3)  &
          -k(118)*n(idx_HCN)  &
          -k(67)*n(idx_O2)  &
          -k(997)  &
          -k(164)*n(idx_NH2)  &
          -k(738)*n(idx_OH)  &
          -k(32)*n(idx_CH2)
      pdj(67) =  &
          +k(67)*n(idx_O2)
      pdj(68) =  &
          +k(107)*n(idx_H2O)
      pdj(69) =  &
          +k(164)*n(idx_NH2)
      pdj(70) =  &
          +k(195)*n(idx_O)
      pdj(71) =  &
          +k(204)*n(idx_OH)
      pdj(73) =  &
          +k(46)*n(idx_CH4)
      pdj(75) =  &
          +k(118)*n(idx_HCN)
      pdj(76) =  &
          +k(179)*n(idx_NH)
    elseif(j==66) then
      pdj(1) =  &
          -k(302)*n(idx_E)
      pdj(2) =  &
          -k(52)*n(idx_CH)
      pdj(3) =  &
          -k(720)*n(idx_O)  &
          -k(196)*n(idx_O)
      pdj(5) =  &
          -k(119)*n(idx_HCN)
      pdj(6) =  &
          +k(397)*n(idx_CH4)  &
          -k(471)*n(idx_H2)
      pdj(7) =  &
          -k(24)*n(idx_C)
      pdj(8) =  &
          +k(471)*n(idx_H2)  &
          +k(628)*n(idx_H2CO)  &
          +k(398)*n(idx_CH4)
      pdj(9) =  &
          -k(499)*n(idx_H2O)  &
          -k(109)*n(idx_H2O)
      pdj(10) =  &
          +k(499)*n(idx_H2O)  &
          -k(205)*n(idx_OH)
      pdj(11) =  &
          -k(153)*n(idx_O2)
      pdj(12) =  &
          -k(35)*n(idx_CH2)
      pdj(13) =  &
          -k(150)*n(idx_H2CO)  &
          -k(628)*n(idx_H2CO)
      pdj(14) =  &
          -k(629)*n(idx_HCO)  &
          -k(151)*n(idx_HCO)
      pdj(15) =  &
          -k(132)*n(idx_MG)
      pdj(16) =  &
          -k(176)*n(idx_NH3)
      pdj(17) =  &
          -k(152)*n(idx_NO)
      pdj(18) =  &
          -k(63)*n(idx_CN)
      pdj(19) =  &
          +k(629)*n(idx_HCO)  &
          -k(68)*n(idx_CO)
      pdj(20) =  &
          +k(154)*n(idx_N)  &
          +k(152)*n(idx_NO)  &
          +k(35)*n(idx_CH2)  &
          +k(52)*n(idx_CH)  &
          +k(166)*n(idx_NH2)  &
          +k(132)*n(idx_MG)  &
          +k(628)*n(idx_H2CO)  &
          +k(196)*n(idx_O)  &
          +k(150)*n(idx_H2CO)  &
          +k(24)*n(idx_C)  &
          +k(397)*n(idx_CH4)  &
          +k(68)*n(idx_CO)  &
          +k(180)*n(idx_NH)  &
          +k(119)*n(idx_HCN)  &
          +k(109)*n(idx_H2O)  &
          +k(176)*n(idx_NH3)  &
          +k(63)*n(idx_CN)  &
          +k(153)*n(idx_O2)  &
          +k(151)*n(idx_HCO)  &
          +k(205)*n(idx_OH)  &
          +k(398)*n(idx_CH4)
      pdj(21) =  &
          -k(166)*n(idx_NH2)
      pdj(23) =  &
          -k(398)*n(idx_CH4)  &
          -k(397)*n(idx_CH4)
      pdj(24) =  &
          +k(720)*n(idx_O)  &
          -k(154)*n(idx_N)  &
          +2.d0*k(302)*n(idx_E)
      pdj(25) =  &
          -k(180)*n(idx_NH)
      pdj(43) =  &
          +k(1091)
      pdj(54) =  &
          +k(151)*n(idx_HCO)  &
          +k(628)*n(idx_H2CO)
      pdj(57) =  &
          +k(24)*n(idx_C)
      pdj(58) =  &
          +k(397)*n(idx_CH4)  &
          +k(35)*n(idx_CH2)
      pdj(59) =  &
          +k(52)*n(idx_CH)
      pdj(60) =  &
          +k(150)*n(idx_H2CO)
      pdj(61) =  &
          +k(132)*n(idx_MG)
      pdj(62) =  &
          +k(176)*n(idx_NH3)
      pdj(63) =  &
          +k(152)*n(idx_NO)  &
          +k(720)*n(idx_O)
      pdj(64) =  &
          +k(63)*n(idx_CN)
      pdj(65) =  &
          +k(68)*n(idx_CO)
      pdj(66) =  &
          -k(720)*n(idx_O)  &
          -k(151)*n(idx_HCO)  &
          -k(176)*n(idx_NH3)  &
          -k(302)*n(idx_E)  &
          -k(205)*n(idx_OH)  &
          -k(397)*n(idx_CH4)  &
          -k(629)*n(idx_HCO)  &
          -k(119)*n(idx_HCN)  &
          -k(109)*n(idx_H2O)  &
          -k(398)*n(idx_CH4)  &
          -k(132)*n(idx_MG)  &
          -k(150)*n(idx_H2CO)  &
          -k(154)*n(idx_N)  &
          -k(471)*n(idx_H2)  &
          -k(35)*n(idx_CH2)  &
          -k(24)*n(idx_C)  &
          -k(196)*n(idx_O)  &
          -k(180)*n(idx_NH)  &
          -k(63)*n(idx_CN)  &
          -k(499)*n(idx_H2O)  &
          -k(68)*n(idx_CO)  &
          -k(1091)  &
          -k(152)*n(idx_NO)  &
          -k(153)*n(idx_O2)  &
          -k(628)*n(idx_H2CO)  &
          -k(166)*n(idx_NH2)  &
          -k(52)*n(idx_CH)
      pdj(67) =  &
          +k(153)*n(idx_O2)
      pdj(68) =  &
          +k(109)*n(idx_H2O)
      pdj(69) =  &
          +k(166)*n(idx_NH2)
      pdj(70) =  &
          +k(196)*n(idx_O)
      pdj(71) =  &
          +k(205)*n(idx_OH)
      pdj(72) =  &
          +k(398)*n(idx_CH4)
      pdj(74) =  &
          +k(154)*n(idx_N)
      pdj(75) =  &
          +k(119)*n(idx_HCN)
      pdj(76) =  &
          +k(180)*n(idx_NH)
      pdj(87) =  &
          +k(471)*n(idx_H2)  &
          +k(499)*n(idx_H2O)  &
          +k(629)*n(idx_HCO)
    elseif(j==67) then
      pdj(1) =  &
          -k(311)*n(idx_E)
      pdj(2) =  &
          -k(415)*n(idx_CH)  &
          -k(55)*n(idx_CH)
      pdj(3) =  &
          +k(337)*n(idx_C)  &
          +k(699)*n(idx_NH)  &
          +k(640)*n(idx_N)  &
          +k(415)*n(idx_CH)  &
          +2.d0*k(311)*n(idx_E)  &
          +k(379)*n(idx_CH2)  &
          +k(1031)
      pdj(7) =  &
          -k(337)*n(idx_C)  &
          -k(25)*n(idx_C)
      pdj(8) =  &
          +k(715)*n(idx_CH3OH)  &
          +k(481)*n(idx_H2CO)
      pdj(11) =  &
          +k(715)*n(idx_CH3OH)  &
          +k(134)*n(idx_MG)  &
          +k(121)*n(idx_HCO)  &
          +k(25)*n(idx_C)  &
          +k(177)*n(idx_NH3)  &
          +k(55)*n(idx_CH)  &
          +k(184)*n(idx_NO)  &
          +k(101)*n(idx_H2CO)  &
          +k(167)*n(idx_NH2)  &
          +k(38)*n(idx_CH2)  &
          +k(481)*n(idx_H2CO)
      pdj(12) =  &
          -k(38)*n(idx_CH2)  &
          -k(379)*n(idx_CH2)
      pdj(13) =  &
          -k(481)*n(idx_H2CO)  &
          -k(101)*n(idx_H2CO)
      pdj(14) =  &
          -k(121)*n(idx_HCO)  &
          -k(555)*n(idx_HCO)
      pdj(15) =  &
          -k(134)*n(idx_MG)
      pdj(16) =  &
          -k(177)*n(idx_NH3)
      pdj(17) =  &
          -k(184)*n(idx_NO)
      pdj(19) =  &
          +k(555)*n(idx_HCO)
      pdj(21) =  &
          -k(167)*n(idx_NH2)
      pdj(24) =  &
          -k(640)*n(idx_N)
      pdj(25) =  &
          -k(699)*n(idx_NH)
      pdj(28) =  &
          -k(715)*n(idx_CH3OH)
      pdj(46) =  &
          +k(1090)
      pdj(54) =  &
          +k(121)*n(idx_HCO)  &
          +k(415)*n(idx_CH)  &
          +k(481)*n(idx_H2CO)
      pdj(57) =  &
          +k(25)*n(idx_C)
      pdj(58) =  &
          +k(38)*n(idx_CH2)
      pdj(59) =  &
          +k(55)*n(idx_CH)
      pdj(60) =  &
          +k(379)*n(idx_CH2)  &
          +k(101)*n(idx_H2CO)
      pdj(61) =  &
          +k(134)*n(idx_MG)
      pdj(62) =  &
          +k(177)*n(idx_NH3)
      pdj(63) =  &
          +k(640)*n(idx_N)  &
          +k(184)*n(idx_NO)
      pdj(65) =  &
          +k(337)*n(idx_C)
      pdj(67) =  &
          -k(1031)  &
          -k(55)*n(idx_CH)  &
          -k(337)*n(idx_C)  &
          -k(167)*n(idx_NH2)  &
          -k(1090)  &
          -k(101)*n(idx_H2CO)  &
          -k(699)*n(idx_NH)  &
          -k(311)*n(idx_E)  &
          -k(184)*n(idx_NO)  &
          -k(640)*n(idx_N)  &
          -k(555)*n(idx_HCO)  &
          -k(121)*n(idx_HCO)  &
          -k(38)*n(idx_CH2)  &
          -k(415)*n(idx_CH)  &
          -k(715)*n(idx_CH3OH)  &
          -k(481)*n(idx_H2CO)  &
          -k(177)*n(idx_NH3)  &
          -k(379)*n(idx_CH2)  &
          -k(25)*n(idx_C)  &
          -k(134)*n(idx_MG)
      pdj(69) =  &
          +k(167)*n(idx_NH2)
      pdj(70) =  &
          +k(1031)
      pdj(79) =  &
          +k(699)*n(idx_NH)
      pdj(82) =  &
          +k(715)*n(idx_CH3OH)
      pdj(88) =  &
          +k(555)*n(idx_HCO)
    elseif(j==68) then
      pdj(1) =  &
          -k(277)*n(idx_E)  &
          -k(278)*n(idx_E)  &
          -k(279)*n(idx_E)
      pdj(2) =  &
          -k(50)*n(idx_CH)  &
          -k(402)*n(idx_CH)
      pdj(3) =  &
          -k(718)*n(idx_O)  &
          +k(278)*n(idx_E)  &
          +k(739)*n(idx_OH)  &
          +k(277)*n(idx_E)
      pdj(4) =  &
          -k(489)*n(idx_HNC)
      pdj(5) =  &
          -k(486)*n(idx_HCN)
      pdj(6) =  &
          +k(718)*n(idx_O)  &
          +k(637)*n(idx_N)  &
          -k(466)*n(idx_H2)  &
          +k(277)*n(idx_E)
      pdj(7) =  &
          -k(329)*n(idx_C)
      pdj(8) =  &
          +k(1006)  &
          +k(636)*n(idx_N)  &
          +2.d0*k(278)*n(idx_E)  &
          +k(279)*n(idx_E)  &
          +k(466)*n(idx_H2)
      pdj(9) =  &
          +k(50)*n(idx_CH)  &
          +k(165)*n(idx_NH2)  &
          +k(103)*n(idx_HCO)  &
          +k(106)*n(idx_O2)  &
          +k(34)*n(idx_CH2)  &
          +k(105)*n(idx_NO)  &
          +k(102)*n(idx_H2CO)  &
          +k(174)*n(idx_NH3)  &
          -k(485)*n(idx_H2O)  &
          +k(104)*n(idx_MG)
      pdj(10) =  &
          +k(402)*n(idx_CH)  &
          +k(676)*n(idx_NH2)  &
          +k(489)*n(idx_HNC)  &
          +k(484)*n(idx_H2CO)  &
          +k(485)*n(idx_H2O)  &
          +k(483)*n(idx_CO)  &
          +k(279)*n(idx_E)  &
          +k(488)*n(idx_HCO)  &
          +k(329)*n(idx_C)  &
          -k(739)*n(idx_OH)  &
          +k(368)*n(idx_CH2)  &
          +k(486)*n(idx_HCN)
      pdj(11) =  &
          -k(106)*n(idx_O2)
      pdj(12) =  &
          -k(34)*n(idx_CH2)  &
          -k(368)*n(idx_CH2)
      pdj(13) =  &
          -k(102)*n(idx_H2CO)  &
          -k(484)*n(idx_H2CO)
      pdj(14) =  &
          -k(103)*n(idx_HCO)  &
          -k(488)*n(idx_HCO)  &
          -k(487)*n(idx_HCO)
      pdj(15) =  &
          -k(104)*n(idx_MG)
      pdj(16) =  &
          -k(174)*n(idx_NH3)
      pdj(17) =  &
          -k(105)*n(idx_NO)
      pdj(19) =  &
          -k(483)*n(idx_CO)  &
          +k(487)*n(idx_HCO)
      pdj(21) =  &
          -k(676)*n(idx_NH2)  &
          -k(165)*n(idx_NH2)
      pdj(22) =  &
          +k(395)*n(idx_CH4)
      pdj(23) =  &
          -k(395)*n(idx_CH4)
      pdj(24) =  &
          -k(637)*n(idx_N)  &
          -k(636)*n(idx_N)  &
          +k(692)*n(idx_NH)
      pdj(25) =  &
          -k(692)*n(idx_NH)
      pdj(40) =  &
          +k(1100)
      pdj(54) =  &
          +k(103)*n(idx_HCO)  &
          +k(483)*n(idx_CO)
      pdj(58) =  &
          +k(402)*n(idx_CH)  &
          +k(34)*n(idx_CH2)
      pdj(59) =  &
          +k(50)*n(idx_CH)  &
          +k(329)*n(idx_C)
      pdj(60) =  &
          +k(102)*n(idx_H2CO)  &
          +k(488)*n(idx_HCO)
      pdj(61) =  &
          +k(104)*n(idx_MG)
      pdj(62) =  &
          +k(174)*n(idx_NH3)  &
          +k(676)*n(idx_NH2)
      pdj(63) =  &
          +k(105)*n(idx_NO)  &
          +k(637)*n(idx_N)
      pdj(67) =  &
          +k(718)*n(idx_O)  &
          +k(106)*n(idx_O2)
      pdj(68) =  &
          -k(50)*n(idx_CH)  &
          -k(174)*n(idx_NH3)  &
          -k(486)*n(idx_HCN)  &
          -k(402)*n(idx_CH)  &
          -k(34)*n(idx_CH2)  &
          -k(165)*n(idx_NH2)  &
          -k(484)*n(idx_H2CO)  &
          -k(483)*n(idx_CO)  &
          -k(488)*n(idx_HCO)  &
          -k(102)*n(idx_H2CO)  &
          -k(487)*n(idx_HCO)  &
          -k(1006)  &
          -k(277)*n(idx_E)  &
          -k(329)*n(idx_C)  &
          -k(105)*n(idx_NO)  &
          -k(637)*n(idx_N)  &
          -k(485)*n(idx_H2O)  &
          -k(739)*n(idx_OH)  &
          -k(106)*n(idx_O2)  &
          -k(489)*n(idx_HNC)  &
          -k(104)*n(idx_MG)  &
          -k(103)*n(idx_HCO)  &
          -k(636)*n(idx_N)  &
          -k(718)*n(idx_O)  &
          -k(676)*n(idx_NH2)  &
          -k(278)*n(idx_E)  &
          -k(1100)  &
          -k(395)*n(idx_CH4)  &
          -k(279)*n(idx_E)  &
          -k(466)*n(idx_H2)  &
          -k(692)*n(idx_NH)  &
          -k(368)*n(idx_CH2)
      pdj(69) =  &
          +k(165)*n(idx_NH2)
      pdj(71) =  &
          +k(1006)
      pdj(72) =  &
          +k(368)*n(idx_CH2)
      pdj(79) =  &
          +k(636)*n(idx_N)
      pdj(82) =  &
          +k(484)*n(idx_H2CO)
      pdj(83) =  &
          +k(739)*n(idx_OH)  &
          +k(485)*n(idx_H2O)  &
          +k(692)*n(idx_NH)  &
          +k(395)*n(idx_CH4)  &
          +k(487)*n(idx_HCO)  &
          +k(466)*n(idx_H2)
      pdj(84) =  &
          +k(489)*n(idx_HNC)  &
          +k(486)*n(idx_HCN)
    elseif(j==69) then
      pdj(1) =  &
          -k(306)*n(idx_E)  &
          -k(307)*n(idx_E)
      pdj(2) =  &
          -k(413)*n(idx_CH)  &
          -k(53)*n(idx_CH)
      pdj(3) =  &
          -k(722)*n(idx_O)  &
          +k(672)*n(idx_O2)
      pdj(4) =  &
          -k(670)*n(idx_HNC)
      pdj(5) =  &
          -k(668)*n(idx_HCN)
      pdj(6) =  &
          -k(474)*n(idx_H2)
      pdj(8) =  &
          +2.d0*k(306)*n(idx_E)  &
          +k(307)*n(idx_E)  &
          +k(639)*n(idx_N)  &
          +k(722)*n(idx_O)  &
          +k(474)*n(idx_H2)
      pdj(9) =  &
          -k(666)*n(idx_H2O)  &
          -k(667)*n(idx_H2O)
      pdj(10) =  &
          +k(673)*n(idx_O2)  &
          +k(667)*n(idx_H2O)
      pdj(11) =  &
          -k(672)*n(idx_O2)  &
          -k(673)*n(idx_O2)
      pdj(12) =  &
          -k(36)*n(idx_CH2)  &
          -k(377)*n(idx_CH2)
      pdj(13) =  &
          -k(665)*n(idx_H2CO)  &
          -k(664)*n(idx_H2CO)
      pdj(14) =  &
          -k(160)*n(idx_HCO)  &
          +k(665)*n(idx_H2CO)  &
          -k(669)*n(idx_HCO)
      pdj(16) =  &
          -k(161)*n(idx_NH3)
      pdj(17) =  &
          -k(162)*n(idx_NO)
      pdj(21) =  &
          +k(160)*n(idx_HCO)  &
          +k(36)*n(idx_CH2)  &
          +k(162)*n(idx_NO)  &
          +k(53)*n(idx_CH)  &
          -k(671)*n(idx_NH2)  &
          +k(161)*n(idx_NH3)
      pdj(24) =  &
          +k(306)*n(idx_E)  &
          +k(697)*n(idx_NH)  &
          -k(639)*n(idx_N)
      pdj(25) =  &
          -k(697)*n(idx_NH)  &
          +k(668)*n(idx_HCN)  &
          +k(413)*n(idx_CH)  &
          +k(377)*n(idx_CH2)  &
          +k(307)*n(idx_E)  &
          +k(664)*n(idx_H2CO)  &
          +k(670)*n(idx_HNC)  &
          +k(671)*n(idx_NH2)  &
          +k(666)*n(idx_H2O)  &
          +k(669)*n(idx_HCO)
      pdj(45) =  &
          +k(1099)
      pdj(54) =  &
          +k(160)*n(idx_HCO)
      pdj(58) =  &
          +k(36)*n(idx_CH2)  &
          +k(413)*n(idx_CH)
      pdj(59) =  &
          +k(53)*n(idx_CH)
      pdj(60) =  &
          +k(669)*n(idx_HCO)
      pdj(62) =  &
          +k(665)*n(idx_H2CO)  &
          +k(697)*n(idx_NH)  &
          +k(667)*n(idx_H2O)  &
          +k(474)*n(idx_H2)  &
          +k(161)*n(idx_NH3)  &
          +k(671)*n(idx_NH2)
      pdj(63) =  &
          +k(162)*n(idx_NO)
      pdj(69) =  &
          -k(474)*n(idx_H2)  &
          -k(306)*n(idx_E)  &
          -k(1099)  &
          -k(670)*n(idx_HNC)  &
          -k(671)*n(idx_NH2)  &
          -k(36)*n(idx_CH2)  &
          -k(722)*n(idx_O)  &
          -k(53)*n(idx_CH)  &
          -k(697)*n(idx_NH)  &
          -k(666)*n(idx_H2O)  &
          -k(665)*n(idx_H2CO)  &
          -k(669)*n(idx_HCO)  &
          -k(160)*n(idx_HCO)  &
          -k(673)*n(idx_O2)  &
          -k(377)*n(idx_CH2)  &
          -k(672)*n(idx_O2)  &
          -k(639)*n(idx_N)  &
          -k(413)*n(idx_CH)  &
          -k(162)*n(idx_NO)  &
          -k(307)*n(idx_E)  &
          -k(664)*n(idx_H2CO)  &
          -k(667)*n(idx_H2O)  &
          -k(668)*n(idx_HCN)  &
          -k(161)*n(idx_NH3)
      pdj(72) =  &
          +k(377)*n(idx_CH2)
      pdj(79) =  &
          +k(722)*n(idx_O)  &
          +k(673)*n(idx_O2)
      pdj(80) =  &
          +k(672)*n(idx_O2)
      pdj(82) =  &
          +k(664)*n(idx_H2CO)
      pdj(83) =  &
          +k(666)*n(idx_H2O)
      pdj(84) =  &
          +k(668)*n(idx_HCN)  &
          +k(670)*n(idx_HNC)
      pdj(87) =  &
          +k(639)*n(idx_N)
    elseif(j==70) then
      pdj(1) =  &
          -k(1063)*n(idx_E)
      pdj(2) =  &
          +k(710)*n(idx_HCN)  &
          -k(54)*n(idx_CH)  &
          -k(414)*n(idx_CH)
      pdj(3) =  &
          +k(37)*n(idx_CH2)  &
          +k(186)*n(idx_CO)  &
          +k(115)*n(idx_H)  &
          +k(188)*n(idx_H2O)  &
          +k(181)*n(idx_NH)  &
          +k(54)*n(idx_CH)  &
          +k(185)*n(idx_CH4)  &
          +k(191)*n(idx_NH3)  &
          +k(189)*n(idx_HCO)  &
          +k(190)*n(idx_NH2)  &
          +k(1063)*n(idx_E)  &
          +k(187)*n(idx_H2CO)  &
          +k(192)*n(idx_O2)  &
          +k(193)*n(idx_OH)
      pdj(5) =  &
          -k(710)*n(idx_HCN)  &
          -k(709)*n(idx_HCN)
      pdj(6) =  &
          -k(475)*n(idx_H2)
      pdj(7) =  &
          -k(1043)*n(idx_C)  &
          +k(706)*n(idx_CN)
      pdj(8) =  &
          +k(414)*n(idx_CH)  &
          +k(698)*n(idx_NH)  &
          +k(714)*n(idx_OH)  &
          -k(115)*n(idx_H)  &
          +k(475)*n(idx_H2)
      pdj(9) =  &
          -k(188)*n(idx_H2O)  &
          +k(703)*n(idx_CH3OH)
      pdj(10) =  &
          +k(704)*n(idx_CH3OH)  &
          -k(193)*n(idx_OH)  &
          +k(708)*n(idx_H2CO)  &
          +k(705)*n(idx_CH4)  &
          -k(714)*n(idx_OH)
      pdj(11) =  &
          -k(192)*n(idx_O2)  &
          +k(713)*n(idx_NO2)
      pdj(12) =  &
          -k(37)*n(idx_CH2)
      pdj(13) =  &
          -k(708)*n(idx_H2CO)  &
          -k(187)*n(idx_H2CO)
      pdj(14) =  &
          -k(189)*n(idx_HCO)  &
          -k(711)*n(idx_HCO)
      pdj(16) =  &
          -k(191)*n(idx_NH3)
      pdj(18) =  &
          -k(706)*n(idx_CN)
      pdj(19) =  &
          +k(711)*n(idx_HCO)  &
          -k(186)*n(idx_CO)  &
          +k(707)*n(idx_CO2)
      pdj(20) =  &
          -k(712)*n(idx_N2)
      pdj(21) =  &
          -k(190)*n(idx_NH2)
      pdj(23) =  &
          -k(185)*n(idx_CH4)  &
          -k(705)*n(idx_CH4)
      pdj(24) =  &
          +k(709)*n(idx_HCN)  &
          +k(712)*n(idx_N2)
      pdj(25) =  &
          -k(698)*n(idx_NH)  &
          -k(181)*n(idx_NH)
      pdj(28) =  &
          -k(703)*n(idx_CH3OH)  &
          -k(704)*n(idx_CH3OH)
      pdj(29) =  &
          -k(707)*n(idx_CO2)
      pdj(32) =  &
          -k(713)*n(idx_NO2)
      pdj(40) =  &
          +k(1089)
      pdj(54) =  &
          +k(708)*n(idx_H2CO)  &
          +k(189)*n(idx_HCO)  &
          +k(709)*n(idx_HCN)
      pdj(55) =  &
          +k(115)*n(idx_H)
      pdj(58) =  &
          +k(37)*n(idx_CH2)
      pdj(59) =  &
          +k(54)*n(idx_CH)
      pdj(60) =  &
          +k(703)*n(idx_CH3OH)  &
          +k(187)*n(idx_H2CO)
      pdj(62) =  &
          +k(191)*n(idx_NH3)
      pdj(63) =  &
          +k(712)*n(idx_N2)  &
          +k(698)*n(idx_NH)  &
          +k(706)*n(idx_CN)  &
          +k(713)*n(idx_NO2)  &
          +k(710)*n(idx_HCN)
      pdj(65) =  &
          +k(1043)*n(idx_C)  &
          +k(414)*n(idx_CH)  &
          +k(186)*n(idx_CO)
      pdj(67) =  &
          +k(707)*n(idx_CO2)  &
          +k(714)*n(idx_OH)  &
          +k(192)*n(idx_O2)
      pdj(68) =  &
          +k(188)*n(idx_H2O)
      pdj(69) =  &
          +k(190)*n(idx_NH2)
      pdj(70) =  &
          -k(710)*n(idx_HCN)  &
          -k(185)*n(idx_CH4)  &
          -k(1089)  &
          -k(709)*n(idx_HCN)  &
          -k(189)*n(idx_HCO)  &
          -k(37)*n(idx_CH2)  &
          -k(191)*n(idx_NH3)  &
          -k(703)*n(idx_CH3OH)  &
          -k(704)*n(idx_CH3OH)  &
          -k(475)*n(idx_H2)  &
          -k(193)*n(idx_OH)  &
          -k(698)*n(idx_NH)  &
          -k(705)*n(idx_CH4)  &
          -k(188)*n(idx_H2O)  &
          -k(187)*n(idx_H2CO)  &
          -k(190)*n(idx_NH2)  &
          -k(192)*n(idx_O2)  &
          -k(708)*n(idx_H2CO)  &
          -k(186)*n(idx_CO)  &
          -k(712)*n(idx_N2)  &
          -k(713)*n(idx_NO2)  &
          -k(115)*n(idx_H)  &
          -k(1063)*n(idx_E)  &
          -k(714)*n(idx_OH)  &
          -k(181)*n(idx_NH)  &
          -k(707)*n(idx_CO2)  &
          -k(414)*n(idx_CH)  &
          -k(54)*n(idx_CH)  &
          -k(711)*n(idx_HCO)  &
          -k(1043)*n(idx_C)  &
          -k(706)*n(idx_CN)
      pdj(71) =  &
          +k(711)*n(idx_HCO)  &
          +k(475)*n(idx_H2)  &
          +k(193)*n(idx_OH)
      pdj(72) =  &
          +k(705)*n(idx_CH4)
      pdj(73) =  &
          +k(185)*n(idx_CH4)
      pdj(76) =  &
          +k(181)*n(idx_NH)
      pdj(82) =  &
          +k(704)*n(idx_CH3OH)
    elseif(j==71) then
      pdj(1) =  &
          -k(313)*n(idx_E)
      pdj(2) =  &
          -k(56)*n(idx_CH)  &
          -k(417)*n(idx_CH)
      pdj(3) =  &
          +k(729)*n(idx_H2CO)  &
          +k(735)*n(idx_N2)  &
          +k(730)*n(idx_H2O)  &
          +k(417)*n(idx_CH)  &
          +k(734)*n(idx_HNC)  &
          +k(736)*n(idx_NO)  &
          +k(733)*n(idx_HCO)  &
          +k(701)*n(idx_NH)  &
          -k(725)*n(idx_O)  &
          +k(313)*n(idx_E)  &
          +k(728)*n(idx_CO)  &
          +k(381)*n(idx_CH2)  &
          +k(731)*n(idx_HCN)  &
          +k(686)*n(idx_NH2)  &
          +k(726)*n(idx_CN)  &
          +k(727)*n(idx_CO2)  &
          +k(339)*n(idx_C)  &
          +k(737)*n(idx_OH)
      pdj(4) =  &
          -k(734)*n(idx_HNC)
      pdj(5) =  &
          -k(731)*n(idx_HCN)
      pdj(6) =  &
          -k(477)*n(idx_H2)
      pdj(7) =  &
          -k(339)*n(idx_C)
      pdj(8) =  &
          +k(1037)  &
          +k(641)*n(idx_N)  &
          +k(313)*n(idx_E)  &
          +k(725)*n(idx_O)  &
          +k(477)*n(idx_H2)
      pdj(9) =  &
          -k(730)*n(idx_H2O)  &
          -k(198)*n(idx_H2O)
      pdj(10) =  &
          +k(168)*n(idx_NH2)  &
          +k(197)*n(idx_H2CO)  &
          -k(737)*n(idx_OH)  &
          +k(198)*n(idx_H2O)  &
          +k(56)*n(idx_CH)  &
          +k(199)*n(idx_HCO)  &
          +k(39)*n(idx_CH2)  &
          +k(200)*n(idx_NH3)  &
          +k(201)*n(idx_NO)  &
          +k(202)*n(idx_O2)
      pdj(11) =  &
          -k(202)*n(idx_O2)
      pdj(12) =  &
          -k(39)*n(idx_CH2)  &
          +k(399)*n(idx_CH4)  &
          -k(381)*n(idx_CH2)
      pdj(13) =  &
          -k(729)*n(idx_H2CO)  &
          -k(197)*n(idx_H2CO)
      pdj(14) =  &
          -k(199)*n(idx_HCO)  &
          -k(733)*n(idx_HCO)  &
          -k(732)*n(idx_HCO)
      pdj(16) =  &
          -k(200)*n(idx_NH3)
      pdj(17) =  &
          -k(736)*n(idx_NO)  &
          -k(201)*n(idx_NO)
      pdj(18) =  &
          -k(726)*n(idx_CN)
      pdj(19) =  &
          -k(728)*n(idx_CO)  &
          +k(732)*n(idx_HCO)
      pdj(20) =  &
          -k(735)*n(idx_N2)
      pdj(21) =  &
          -k(686)*n(idx_NH2)  &
          -k(168)*n(idx_NH2)
      pdj(23) =  &
          -k(399)*n(idx_CH4)
      pdj(24) =  &
          -k(641)*n(idx_N)
      pdj(25) =  &
          -k(701)*n(idx_NH)
      pdj(29) =  &
          -k(727)*n(idx_CO2)
      pdj(40) =  &
          +k(1094)
      pdj(54) =  &
          +k(728)*n(idx_CO)  &
          +k(199)*n(idx_HCO)
      pdj(58) =  &
          +k(417)*n(idx_CH)  &
          +k(39)*n(idx_CH2)
      pdj(59) =  &
          +k(56)*n(idx_CH)  &
          +k(339)*n(idx_C)
      pdj(60) =  &
          +k(733)*n(idx_HCO)  &
          +k(197)*n(idx_H2CO)
      pdj(62) =  &
          +k(686)*n(idx_NH2)  &
          +k(200)*n(idx_NH3)
      pdj(63) =  &
          +k(641)*n(idx_N)  &
          +k(201)*n(idx_NO)
      pdj(67) =  &
          +k(725)*n(idx_O)  &
          +k(202)*n(idx_O2)
      pdj(68) =  &
          +k(732)*n(idx_HCO)  &
          +k(198)*n(idx_H2O)  &
          +k(737)*n(idx_OH)  &
          +k(477)*n(idx_H2)
      pdj(69) =  &
          +k(168)*n(idx_NH2)  &
          +k(701)*n(idx_NH)
      pdj(70) =  &
          +k(1037)
      pdj(71) =  &
          -k(728)*n(idx_CO)  &
          -k(726)*n(idx_CN)  &
          -k(477)*n(idx_H2)  &
          -k(725)*n(idx_O)  &
          -k(727)*n(idx_CO2)  &
          -k(736)*n(idx_NO)  &
          -k(381)*n(idx_CH2)  &
          -k(199)*n(idx_HCO)  &
          -k(313)*n(idx_E)  &
          -k(168)*n(idx_NH2)  &
          -k(731)*n(idx_HCN)  &
          -k(198)*n(idx_H2O)  &
          -k(735)*n(idx_N2)  &
          -k(200)*n(idx_NH3)  &
          -k(730)*n(idx_H2O)  &
          -k(399)*n(idx_CH4)  &
          -k(1037)  &
          -k(729)*n(idx_H2CO)  &
          -k(339)*n(idx_C)  &
          -k(734)*n(idx_HNC)  &
          -k(1094)  &
          -k(202)*n(idx_O2)  &
          -k(201)*n(idx_NO)  &
          -k(732)*n(idx_HCO)  &
          -k(39)*n(idx_CH2)  &
          -k(737)*n(idx_OH)  &
          -k(641)*n(idx_N)  &
          -k(56)*n(idx_CH)  &
          -k(701)*n(idx_NH)  &
          -k(417)*n(idx_CH)  &
          -k(733)*n(idx_HCO)  &
          -k(686)*n(idx_NH2)  &
          -k(197)*n(idx_H2CO)
      pdj(72) =  &
          +k(381)*n(idx_CH2)
      pdj(75) =  &
          +k(726)*n(idx_CN)
      pdj(79) =  &
          +k(736)*n(idx_NO)
      pdj(82) =  &
          +k(729)*n(idx_H2CO)
      pdj(83) =  &
          +k(399)*n(idx_CH4)  &
          +k(730)*n(idx_H2O)
      pdj(84) =  &
          +k(731)*n(idx_HCN)  &
          +k(734)*n(idx_HNC)
      pdj(85) =  &
          +k(727)*n(idx_CO2)
      pdj(87) =  &
          +k(735)*n(idx_N2)
    elseif(j==72) then
      pdj(1) =  &
          -k(265)*n(idx_E)  &
          -k(264)*n(idx_E)  &
          -k(263)*n(idx_E)  &
          -k(1057)*n(idx_E)
      pdj(2) =  &
          +k(265)*n(idx_E)  &
          +k(264)*n(idx_E)
      pdj(3) =  &
          +k(385)*n(idx_O2)  &
          -k(386)*n(idx_O)  &
          -k(387)*n(idx_O)
      pdj(6) =  &
          +k(532)*n(idx_H)  &
          +k(264)*n(idx_E)  &
          +k(387)*n(idx_O)  &
          +k(980)  &
          +k(689)*n(idx_NH)  &
          +k(388)*n(idx_OH)
      pdj(8) =  &
          +k(263)*n(idx_E)  &
          +k(386)*n(idx_O)  &
          +2.d0*k(265)*n(idx_E)  &
          +k(981)  &
          -k(532)*n(idx_H)
      pdj(10) =  &
          -k(388)*n(idx_OH)
      pdj(11) =  &
          -k(385)*n(idx_O2)
      pdj(12) =  &
          +k(263)*n(idx_E)
      pdj(13) =  &
          -k(383)*n(idx_H2CO)
      pdj(14) =  &
          -k(40)*n(idx_HCO)  &
          -k(384)*n(idx_HCO)
      pdj(15) =  &
          -k(41)*n(idx_MG)
      pdj(17) =  &
          -k(42)*n(idx_NO)
      pdj(19) =  &
          +k(384)*n(idx_HCO)
      pdj(22) =  &
          +k(41)*n(idx_MG)  &
          +k(42)*n(idx_NO)  &
          +k(1057)*n(idx_E)  &
          +k(40)*n(idx_HCO)
      pdj(23) =  &
          +k(382)*n(idx_CH3OH)  &
          +k(383)*n(idx_H2CO)
      pdj(25) =  &
          -k(689)*n(idx_NH)
      pdj(28) =  &
          -k(382)*n(idx_CH3OH)
      pdj(38) =  &
          +k(1105)
      pdj(54) =  &
          +k(387)*n(idx_O)  &
          +k(383)*n(idx_H2CO)  &
          +k(40)*n(idx_HCO)
      pdj(58) =  &
          +k(981)  &
          +k(532)*n(idx_H)
      pdj(59) =  &
          +k(980)
      pdj(60) =  &
          +k(386)*n(idx_O)  &
          +k(388)*n(idx_OH)
      pdj(61) =  &
          +k(41)*n(idx_MG)
      pdj(63) =  &
          +k(42)*n(idx_NO)
      pdj(72) =  &
          -k(1057)*n(idx_E)  &
          -k(689)*n(idx_NH)  &
          -k(265)*n(idx_E)  &
          -k(264)*n(idx_E)  &
          -k(263)*n(idx_E)  &
          -k(532)*n(idx_H)  &
          -k(387)*n(idx_O)  &
          -k(385)*n(idx_O2)  &
          -k(980)  &
          -k(1105)  &
          -k(386)*n(idx_O)  &
          -k(388)*n(idx_OH)  &
          -k(384)*n(idx_HCO)  &
          -k(40)*n(idx_HCO)  &
          -k(981)  &
          -k(41)*n(idx_MG)  &
          -k(382)*n(idx_CH3OH)  &
          -k(42)*n(idx_NO)  &
          -k(383)*n(idx_H2CO)
      pdj(73) =  &
          +k(384)*n(idx_HCO)
      pdj(82) =  &
          +k(385)*n(idx_O2)  &
          +k(382)*n(idx_CH3OH)
      pdj(84) =  &
          +k(689)*n(idx_NH)
    elseif(j==73) then
      pdj(1) =  &
          -k(267)*n(idx_E)  &
          -k(266)*n(idx_E)
      pdj(3) =  &
          -k(717)*n(idx_O)
      pdj(6) =  &
          +k(533)*n(idx_H)  &
          +k(988)
      pdj(8) =  &
          +k(989)  &
          -k(533)*n(idx_H)  &
          +k(267)*n(idx_E)  &
          +2.d0*k(266)*n(idx_E)
      pdj(9) =  &
          -k(392)*n(idx_H2O)
      pdj(10) =  &
          +k(717)*n(idx_O)
      pdj(11) =  &
          -k(45)*n(idx_O2)
      pdj(12) =  &
          +k(266)*n(idx_E)
      pdj(13) =  &
          -k(43)*n(idx_H2CO)  &
          -k(391)*n(idx_H2CO)
      pdj(16) =  &
          -k(44)*n(idx_NH3)
      pdj(19) =  &
          -k(390)*n(idx_CO)
      pdj(22) =  &
          +k(389)*n(idx_CO2)  &
          +k(267)*n(idx_E)  &
          +k(390)*n(idx_CO)  &
          +k(391)*n(idx_H2CO)  &
          +k(392)*n(idx_H2O)
      pdj(23) =  &
          +k(45)*n(idx_O2)  &
          +k(43)*n(idx_H2CO)  &
          +k(44)*n(idx_NH3)
      pdj(29) =  &
          -k(389)*n(idx_CO2)
      pdj(38) =  &
          +k(1108)
      pdj(54) =  &
          +k(390)*n(idx_CO)
      pdj(58) =  &
          +k(988)
      pdj(60) =  &
          +k(43)*n(idx_H2CO)
      pdj(62) =  &
          +k(44)*n(idx_NH3)
      pdj(67) =  &
          +k(45)*n(idx_O2)
      pdj(72) =  &
          +k(989)  &
          +k(533)*n(idx_H)  &
          +k(717)*n(idx_O)
      pdj(73) =  &
          -k(1108)  &
          -k(267)*n(idx_E)  &
          -k(717)*n(idx_O)  &
          -k(391)*n(idx_H2CO)  &
          -k(43)*n(idx_H2CO)  &
          -k(266)*n(idx_E)  &
          -k(988)  &
          -k(989)  &
          -k(390)*n(idx_CO)  &
          -k(45)*n(idx_O2)  &
          -k(533)*n(idx_H)  &
          -k(392)*n(idx_H2O)  &
          -k(44)*n(idx_NH3)  &
          -k(389)*n(idx_CO2)
      pdj(82) =  &
          +k(391)*n(idx_H2CO)
      pdj(83) =  &
          +k(392)*n(idx_H2O)
      pdj(85) =  &
          +k(389)*n(idx_CO2)
    elseif(j==74) then
      pdj(1) =  &
          -k(1062)*n(idx_E)
      pdj(2) =  &
          -k(410)*n(idx_CH)  &
          -k(51)*n(idx_CH)
      pdj(3) =  &
          +k(626)*n(idx_O2)  &
          +k(625)*n(idx_NO)
      pdj(5) =  &
          -k(141)*n(idx_HCN)
      pdj(6) =  &
          +k(622)*n(idx_NH3)  &
          +k(615)*n(idx_CH4)  &
          -k(470)*n(idx_H2)
      pdj(7) =  &
          +k(618)*n(idx_CO)
      pdj(8) =  &
          +k(410)*n(idx_CH)  &
          +k(610)*n(idx_CH3OH)  &
          +k(470)*n(idx_H2)  &
          +k(612)*n(idx_CH3OH)  &
          +k(613)*n(idx_CH3OH)  &
          +k(614)*n(idx_CH4)  &
          +2.d0*k(616)*n(idx_CH4)  &
          +k(615)*n(idx_CH4)  &
          +k(624)*n(idx_NH)
      pdj(9) =  &
          -k(140)*n(idx_H2O)
      pdj(10) =  &
          -k(149)*n(idx_OH)
      pdj(11) =  &
          -k(148)*n(idx_O2)  &
          -k(627)*n(idx_O2)  &
          -k(626)*n(idx_O2)
      pdj(12) =  &
          +k(620)*n(idx_H2CO)  &
          -k(135)*n(idx_CH2)
      pdj(13) =  &
          -k(620)*n(idx_H2CO)  &
          -k(139)*n(idx_H2CO)  &
          -k(619)*n(idx_H2CO)
      pdj(14) =  &
          -k(142)*n(idx_HCO)  &
          -k(621)*n(idx_HCO)
      pdj(15) =  &
          -k(143)*n(idx_MG)
      pdj(16) =  &
          -k(623)*n(idx_NH3)  &
          -k(622)*n(idx_NH3)  &
          -k(145)*n(idx_NH3)
      pdj(17) =  &
          -k(147)*n(idx_NO)  &
          +k(613)*n(idx_CH3OH)  &
          +k(617)*n(idx_CO2)  &
          -k(625)*n(idx_NO)  &
          +k(627)*n(idx_O2)
      pdj(18) =  &
          -k(137)*n(idx_CN)
      pdj(19) =  &
          +k(621)*n(idx_HCO)  &
          -k(138)*n(idx_CO)  &
          -k(618)*n(idx_CO)
      pdj(21) =  &
          -k(144)*n(idx_NH2)
      pdj(22) =  &
          +k(612)*n(idx_CH3OH)
      pdj(23) =  &
          -k(614)*n(idx_CH4)  &
          -k(136)*n(idx_CH4)  &
          -k(615)*n(idx_CH4)  &
          -k(616)*n(idx_CH4)
      pdj(24) =  &
          +k(148)*n(idx_O2)  &
          +k(141)*n(idx_HCN)  &
          +k(149)*n(idx_OH)  &
          +k(136)*n(idx_CH4)  &
          +k(614)*n(idx_CH4)  &
          +k(146)*n(idx_NH)  &
          -k(1054)*n(idx_N)  &
          +k(137)*n(idx_CN)  &
          +k(143)*n(idx_MG)  &
          +k(51)*n(idx_CH)  &
          +k(140)*n(idx_H2O)  &
          +k(138)*n(idx_CO)  &
          +k(139)*n(idx_H2CO)  &
          +k(135)*n(idx_CH2)  &
          +k(1062)*n(idx_E)  &
          +k(147)*n(idx_NO)  &
          +k(142)*n(idx_HCO)  &
          +k(145)*n(idx_NH3)  &
          +k(144)*n(idx_NH2)
      pdj(25) =  &
          +k(619)*n(idx_H2CO)  &
          -k(146)*n(idx_NH)  &
          -k(624)*n(idx_NH)  &
          +k(610)*n(idx_CH3OH)  &
          +k(611)*n(idx_CH3OH)  &
          +k(623)*n(idx_NH3)
      pdj(28) =  &
          -k(611)*n(idx_CH3OH)  &
          -k(613)*n(idx_CH3OH)  &
          -k(610)*n(idx_CH3OH)  &
          -k(612)*n(idx_CH3OH)
      pdj(29) =  &
          -k(617)*n(idx_CO2)
      pdj(45) =  &
          +k(1088)
      pdj(54) =  &
          +k(619)*n(idx_H2CO)  &
          +k(142)*n(idx_HCO)
      pdj(58) =  &
          +k(135)*n(idx_CH2)
      pdj(59) =  &
          +k(51)*n(idx_CH)
      pdj(60) =  &
          +k(610)*n(idx_CH3OH)  &
          +k(139)*n(idx_H2CO)
      pdj(61) =  &
          +k(143)*n(idx_MG)
      pdj(62) =  &
          +k(145)*n(idx_NH3)
      pdj(63) =  &
          +k(626)*n(idx_O2)  &
          +k(147)*n(idx_NO)  &
          +k(612)*n(idx_CH3OH)  &
          +k(620)*n(idx_H2CO)  &
          +k(618)*n(idx_CO)
      pdj(64) =  &
          +k(137)*n(idx_CN)  &
          +k(410)*n(idx_CH)
      pdj(65) =  &
          +k(617)*n(idx_CO2)  &
          +k(138)*n(idx_CO)
      pdj(66) =  &
          +k(625)*n(idx_NO)  &
          +k(624)*n(idx_NH)  &
          +k(1054)*n(idx_N)
      pdj(67) =  &
          +k(148)*n(idx_O2)
      pdj(68) =  &
          +k(140)*n(idx_H2O)
      pdj(69) =  &
          +k(623)*n(idx_NH3)  &
          +k(144)*n(idx_NH2)
      pdj(70) =  &
          +k(627)*n(idx_O2)
      pdj(71) =  &
          +k(149)*n(idx_OH)
      pdj(72) =  &
          +k(613)*n(idx_CH3OH)  &
          +k(614)*n(idx_CH4)
      pdj(73) =  &
          +k(136)*n(idx_CH4)
      pdj(74) =  &
          -k(146)*n(idx_NH)  &
          -k(148)*n(idx_O2)  &
          -k(611)*n(idx_CH3OH)  &
          -k(614)*n(idx_CH4)  &
          -k(1054)*n(idx_N)  &
          -k(470)*n(idx_H2)  &
          -k(622)*n(idx_NH3)  &
          -k(142)*n(idx_HCO)  &
          -k(136)*n(idx_CH4)  &
          -k(610)*n(idx_CH3OH)  &
          -k(617)*n(idx_CO2)  &
          -k(138)*n(idx_CO)  &
          -k(410)*n(idx_CH)  &
          -k(141)*n(idx_HCN)  &
          -k(612)*n(idx_CH3OH)  &
          -k(149)*n(idx_OH)  &
          -k(616)*n(idx_CH4)  &
          -k(625)*n(idx_NO)  &
          -k(144)*n(idx_NH2)  &
          -k(621)*n(idx_HCO)  &
          -k(137)*n(idx_CN)  &
          -k(626)*n(idx_O2)  &
          -k(623)*n(idx_NH3)  &
          -k(140)*n(idx_H2O)  &
          -k(1088)  &
          -k(619)*n(idx_H2CO)  &
          -k(613)*n(idx_CH3OH)  &
          -k(51)*n(idx_CH)  &
          -k(624)*n(idx_NH)  &
          -k(139)*n(idx_H2CO)  &
          -k(627)*n(idx_O2)  &
          -k(620)*n(idx_H2CO)  &
          -k(143)*n(idx_MG)  &
          -k(615)*n(idx_CH4)  &
          -k(135)*n(idx_CH2)  &
          -k(145)*n(idx_NH3)  &
          -k(147)*n(idx_NO)  &
          -k(1062)*n(idx_E)  &
          -k(618)*n(idx_CO)
      pdj(75) =  &
          +k(141)*n(idx_HCN)  &
          +k(615)*n(idx_CH4)
      pdj(76) =  &
          +k(146)*n(idx_NH)  &
          +k(621)*n(idx_HCO)  &
          +k(470)*n(idx_H2)
      pdj(82) =  &
          +k(611)*n(idx_CH3OH)
      pdj(84) =  &
          +k(616)*n(idx_CH4)
      pdj(87) =  &
          +k(622)*n(idx_NH3)
    elseif(j==75) then
      pdj(1) =  &
          -k(291)*n(idx_E)
      pdj(2) =  &
          -k(405)*n(idx_CH)
      pdj(4) =  &
          -k(541)*n(idx_HNC)
      pdj(5) =  &
          +k(117)*n(idx_O2)  &
          +k(175)*n(idx_NH3)  &
          +k(108)*n(idx_H2O)  &
          +k(113)*n(idx_H)  &
          +k(116)*n(idx_NO)  &
          -k(538)*n(idx_HCN)
      pdj(6) =  &
          -k(467)*n(idx_H2)
      pdj(7) =  &
          -k(331)*n(idx_C)
      pdj(8) =  &
          -k(113)*n(idx_H)  &
          +k(467)*n(idx_H2)  &
          +k(291)*n(idx_E)
      pdj(9) =  &
          -k(495)*n(idx_H2O)  &
          -k(108)*n(idx_H2O)
      pdj(10) =  &
          -k(740)*n(idx_OH)
      pdj(11) =  &
          -k(117)*n(idx_O2)
      pdj(12) =  &
          -k(370)*n(idx_CH2)
      pdj(13) =  &
          -k(537)*n(idx_H2CO)
      pdj(14) =  &
          -k(540)*n(idx_HCO)  &
          -k(539)*n(idx_HCO)
      pdj(16) =  &
          -k(688)*n(idx_NH3)  &
          -k(175)*n(idx_NH3)
      pdj(17) =  &
          -k(116)*n(idx_NO)
      pdj(18) =  &
          +k(538)*n(idx_HCN)  &
          +k(539)*n(idx_HCO)  &
          +k(370)*n(idx_CH2)  &
          +k(740)*n(idx_OH)  &
          +k(495)*n(idx_H2O)  &
          +k(537)*n(idx_H2CO)  &
          +k(541)*n(idx_HNC)  &
          +k(291)*n(idx_E)  &
          +k(536)*n(idx_CO)  &
          +k(405)*n(idx_CH)  &
          +k(331)*n(idx_C)  &
          +k(535)*n(idx_CO2)  &
          +k(693)*n(idx_NH)  &
          +k(679)*n(idx_NH2)
      pdj(19) =  &
          +k(540)*n(idx_HCO)  &
          -k(536)*n(idx_CO)
      pdj(21) =  &
          -k(679)*n(idx_NH2)  &
          +k(688)*n(idx_NH3)
      pdj(22) =  &
          +k(396)*n(idx_CH4)
      pdj(23) =  &
          -k(396)*n(idx_CH4)
      pdj(25) =  &
          -k(693)*n(idx_NH)
      pdj(29) =  &
          -k(535)*n(idx_CO2)
      pdj(44) =  &
          +k(1102)
      pdj(54) =  &
          +k(536)*n(idx_CO)
      pdj(55) =  &
          +k(113)*n(idx_H)
      pdj(58) =  &
          +k(405)*n(idx_CH)
      pdj(59) =  &
          +k(331)*n(idx_C)
      pdj(60) =  &
          +k(539)*n(idx_HCO)
      pdj(62) =  &
          +k(175)*n(idx_NH3)  &
          +k(679)*n(idx_NH2)
      pdj(63) =  &
          +k(116)*n(idx_NO)
      pdj(67) =  &
          +k(117)*n(idx_O2)
      pdj(68) =  &
          +k(108)*n(idx_H2O)  &
          +k(740)*n(idx_OH)
      pdj(69) =  &
          +k(693)*n(idx_NH)
      pdj(72) =  &
          +k(370)*n(idx_CH2)
      pdj(75) =  &
          -k(108)*n(idx_H2O)  &
          -k(740)*n(idx_OH)  &
          -k(331)*n(idx_C)  &
          -k(396)*n(idx_CH4)  &
          -k(116)*n(idx_NO)  &
          -k(693)*n(idx_NH)  &
          -k(405)*n(idx_CH)  &
          -k(113)*n(idx_H)  &
          -k(535)*n(idx_CO2)  &
          -k(1102)  &
          -k(117)*n(idx_O2)  &
          -k(688)*n(idx_NH3)  &
          -k(495)*n(idx_H2O)  &
          -k(370)*n(idx_CH2)  &
          -k(537)*n(idx_H2CO)  &
          -k(467)*n(idx_H2)  &
          -k(679)*n(idx_NH2)  &
          -k(291)*n(idx_E)  &
          -k(536)*n(idx_CO)  &
          -k(540)*n(idx_HCO)  &
          -k(539)*n(idx_HCO)  &
          -k(538)*n(idx_HCN)  &
          -k(175)*n(idx_NH3)  &
          -k(541)*n(idx_HNC)
      pdj(82) =  &
          +k(537)*n(idx_H2CO)
      pdj(83) =  &
          +k(495)*n(idx_H2O)
      pdj(84) =  &
          +k(396)*n(idx_CH4)  &
          +k(540)*n(idx_HCO)  &
          +k(541)*n(idx_HNC)  &
          +k(688)*n(idx_NH3)  &
          +k(467)*n(idx_H2)  &
          +k(538)*n(idx_HCN)
      pdj(85) =  &
          +k(535)*n(idx_CO2)
    elseif(j==76) then
      pdj(1) =  &
          -k(305)*n(idx_E)
      pdj(2) =  &
          -k(412)*n(idx_CH)
      pdj(3) =  &
          -k(662)*n(idx_O)  &
          +k(659)*n(idx_NO)  &
          +k(651)*n(idx_H2O)
      pdj(4) =  &
          -k(655)*n(idx_HNC)
      pdj(5) =  &
          -k(653)*n(idx_HCN)
      pdj(6) =  &
          +k(650)*n(idx_H2O)  &
          -k(473)*n(idx_H2)  &
          -k(472)*n(idx_H2)
      pdj(7) =  &
          -k(336)*n(idx_C)
      pdj(8) =  &
          +k(473)*n(idx_H2)  &
          +k(305)*n(idx_E)  &
          +k(638)*n(idx_N)
      pdj(9) =  &
          -k(652)*n(idx_H2O)  &
          -k(651)*n(idx_H2O)  &
          -k(156)*n(idx_H2O)  &
          -k(649)*n(idx_H2O)  &
          -k(650)*n(idx_H2O)
      pdj(10) =  &
          +k(660)*n(idx_O2)  &
          +k(652)*n(idx_H2O)  &
          -k(663)*n(idx_OH)
      pdj(11) =  &
          -k(159)*n(idx_O2)  &
          -k(661)*n(idx_O2)  &
          -k(660)*n(idx_O2)
      pdj(12) =  &
          -k(376)*n(idx_CH2)
      pdj(13) =  &
          -k(647)*n(idx_H2CO)  &
          -k(648)*n(idx_H2CO)  &
          -k(155)*n(idx_H2CO)
      pdj(14) =  &
          +k(645)*n(idx_CO2)  &
          -k(654)*n(idx_HCO)
      pdj(16) =  &
          -k(157)*n(idx_NH3)
      pdj(17) =  &
          -k(158)*n(idx_NO)  &
          -k(659)*n(idx_NO)
      pdj(18) =  &
          -k(642)*n(idx_CN)
      pdj(19) =  &
          +k(644)*n(idx_CO2)  &
          -k(646)*n(idx_CO)
      pdj(20) =  &
          -k(656)*n(idx_N2)
      pdj(21) =  &
          +k(648)*n(idx_H2CO)  &
          -k(657)*n(idx_NH2)
      pdj(24) =  &
          +k(1020)  &
          +k(376)*n(idx_CH2)  &
          +k(643)*n(idx_CO2)  &
          +k(657)*n(idx_NH2)  &
          +k(646)*n(idx_CO)  &
          +k(663)*n(idx_OH)  &
          +k(336)*n(idx_C)  &
          +k(412)*n(idx_CH)  &
          +k(653)*n(idx_HCN)  &
          +k(642)*n(idx_CN)  &
          -k(638)*n(idx_N)  &
          +k(661)*n(idx_O2)  &
          +k(647)*n(idx_H2CO)  &
          +k(658)*n(idx_NH)  &
          +k(655)*n(idx_HNC)  &
          +k(656)*n(idx_N2)  &
          +k(662)*n(idx_O)  &
          +k(472)*n(idx_H2)  &
          +k(305)*n(idx_E)  &
          +k(649)*n(idx_H2O)  &
          +k(654)*n(idx_HCO)
      pdj(25) =  &
          +k(156)*n(idx_H2O)  &
          -k(658)*n(idx_NH)  &
          +k(155)*n(idx_H2CO)  &
          +k(158)*n(idx_NO)  &
          +k(159)*n(idx_O2)  &
          +k(157)*n(idx_NH3)
      pdj(29) =  &
          -k(643)*n(idx_CO2)  &
          -k(644)*n(idx_CO2)  &
          -k(645)*n(idx_CO2)
      pdj(45) =  &
          +k(1093)
      pdj(54) =  &
          +k(648)*n(idx_H2CO)  &
          +k(646)*n(idx_CO)
      pdj(55) =  &
          +k(1020)
      pdj(58) =  &
          +k(412)*n(idx_CH)
      pdj(59) =  &
          +k(336)*n(idx_C)
      pdj(60) =  &
          +k(155)*n(idx_H2CO)  &
          +k(654)*n(idx_HCO)
      pdj(62) =  &
          +k(657)*n(idx_NH2)  &
          +k(157)*n(idx_NH3)  &
          +k(651)*n(idx_H2O)
      pdj(63) =  &
          +k(645)*n(idx_CO2)  &
          +k(660)*n(idx_O2)  &
          +k(158)*n(idx_NO)
      pdj(66) =  &
          +k(638)*n(idx_N)
      pdj(67) =  &
          +k(159)*n(idx_O2)
      pdj(68) =  &
          +k(156)*n(idx_H2O)  &
          +k(663)*n(idx_OH)
      pdj(69) =  &
          +k(473)*n(idx_H2)  &
          +k(658)*n(idx_NH)  &
          +k(652)*n(idx_H2O)
      pdj(71) =  &
          +k(662)*n(idx_O)
      pdj(72) =  &
          +k(376)*n(idx_CH2)
      pdj(75) =  &
          +k(642)*n(idx_CN)
      pdj(76) =  &
          -k(473)*n(idx_H2)  &
          -k(654)*n(idx_HCO)  &
          -k(1020)  &
          -k(157)*n(idx_NH3)  &
          -k(651)*n(idx_H2O)  &
          -k(642)*n(idx_CN)  &
          -k(305)*n(idx_E)  &
          -k(661)*n(idx_O2)  &
          -k(656)*n(idx_N2)  &
          -k(657)*n(idx_NH2)  &
          -k(663)*n(idx_OH)  &
          -k(376)*n(idx_CH2)  &
          -k(652)*n(idx_H2O)  &
          -k(638)*n(idx_N)  &
          -k(644)*n(idx_CO2)  &
          -k(660)*n(idx_O2)  &
          -k(653)*n(idx_HCN)  &
          -k(156)*n(idx_H2O)  &
          -k(649)*n(idx_H2O)  &
          -k(336)*n(idx_C)  &
          -k(158)*n(idx_NO)  &
          -k(472)*n(idx_H2)  &
          -k(645)*n(idx_CO2)  &
          -k(646)*n(idx_CO)  &
          -k(658)*n(idx_NH)  &
          -k(159)*n(idx_O2)  &
          -k(650)*n(idx_H2O)  &
          -k(655)*n(idx_HNC)  &
          -k(659)*n(idx_NO)  &
          -k(648)*n(idx_H2CO)  &
          -k(155)*n(idx_H2CO)  &
          -k(662)*n(idx_O)  &
          -k(1093)  &
          -k(412)*n(idx_CH)  &
          -k(647)*n(idx_H2CO)  &
          -k(643)*n(idx_CO2)
      pdj(79) =  &
          +k(644)*n(idx_CO2)  &
          +k(650)*n(idx_H2O)
      pdj(81) =  &
          +k(472)*n(idx_H2)
      pdj(82) =  &
          +k(647)*n(idx_H2CO)
      pdj(83) =  &
          +k(649)*n(idx_H2O)
      pdj(84) =  &
          +k(655)*n(idx_HNC)  &
          +k(653)*n(idx_HCN)
      pdj(85) =  &
          +k(643)*n(idx_CO2)
      pdj(87) =  &
          +k(656)*n(idx_N2)  &
          +k(659)*n(idx_NO)
      pdj(88) =  &
          +k(661)*n(idx_O2)
    elseif(j==77) then
      pdj(1) =  &
          -k(270)*n(idx_E)
      pdj(2) =  &
          -k(87)*n(idx_CH)  &
          -k(444)*n(idx_CH)
      pdj(3) =  &
          -k(458)*n(idx_O)
      pdj(5) =  &
          -k(92)*n(idx_HCN)
      pdj(6) =  &
          -k(448)*n(idx_H2)  &
          +k(90)*n(idx_H2CO)  &
          +k(88)*n(idx_CN)  &
          +k(112)*n(idx_H)  &
          +k(97)*n(idx_NO)  &
          +k(85)*n(idx_CH2)  &
          +k(92)*n(idx_HCN)  &
          +k(89)*n(idx_CO)  &
          +k(91)*n(idx_H2O)  &
          +k(443)*n(idx_CH4)  &
          +k(94)*n(idx_NH2)  &
          +k(96)*n(idx_NH)  &
          +k(87)*n(idx_CH)  &
          +k(86)*n(idx_CH4)  &
          +k(449)*n(idx_H2CO)  &
          +k(99)*n(idx_OH)  &
          +k(98)*n(idx_O2)  &
          +k(95)*n(idx_NH3)  &
          +k(93)*n(idx_HCO)
      pdj(7) =  &
          -k(441)*n(idx_C)
      pdj(8) =  &
          +k(446)*n(idx_CO2)  &
          +k(445)*n(idx_CN)  &
          +k(450)*n(idx_H2O)  &
          +k(455)*n(idx_NH)  &
          +2.d0*k(270)*n(idx_E)  &
          +k(448)*n(idx_H2)  &
          +k(454)*n(idx_N)  &
          +k(441)*n(idx_C)  &
          +k(453)*n(idx_N2)  &
          +k(1000)  &
          +k(458)*n(idx_O)  &
          +k(444)*n(idx_CH)  &
          +k(449)*n(idx_H2CO)  &
          +k(452)*n(idx_HE)  &
          +k(457)*n(idx_O2)  &
          +k(459)*n(idx_OH)  &
          +k(442)*n(idx_CH2)  &
          +k(443)*n(idx_CH4)  &
          +k(447)*n(idx_CO)  &
          -k(112)*n(idx_H)  &
          +k(456)*n(idx_NO)
      pdj(9) =  &
          -k(91)*n(idx_H2O)  &
          -k(450)*n(idx_H2O)
      pdj(10) =  &
          -k(459)*n(idx_OH)  &
          -k(99)*n(idx_OH)
      pdj(11) =  &
          -k(98)*n(idx_O2)  &
          -k(457)*n(idx_O2)
      pdj(12) =  &
          -k(85)*n(idx_CH2)  &
          -k(442)*n(idx_CH2)
      pdj(13) =  &
          -k(449)*n(idx_H2CO)  &
          -k(90)*n(idx_H2CO)
      pdj(14) =  &
          -k(451)*n(idx_HCO)  &
          -k(93)*n(idx_HCO)
      pdj(16) =  &
          -k(95)*n(idx_NH3)
      pdj(17) =  &
          -k(97)*n(idx_NO)  &
          -k(456)*n(idx_NO)
      pdj(18) =  &
          -k(88)*n(idx_CN)  &
          -k(445)*n(idx_CN)
      pdj(19) =  &
          -k(447)*n(idx_CO)  &
          -k(89)*n(idx_CO)  &
          +k(451)*n(idx_HCO)
      pdj(20) =  &
          -k(453)*n(idx_N2)
      pdj(21) =  &
          -k(94)*n(idx_NH2)
      pdj(23) =  &
          -k(86)*n(idx_CH4)  &
          -k(443)*n(idx_CH4)
      pdj(24) =  &
          -k(454)*n(idx_N)
      pdj(25) =  &
          -k(96)*n(idx_NH)  &
          -k(455)*n(idx_NH)
      pdj(26) =  &
          -k(452)*n(idx_HE)
      pdj(29) =  &
          -k(446)*n(idx_CO2)
      pdj(54) =  &
          +k(449)*n(idx_H2CO)  &
          +k(447)*n(idx_CO)  &
          +k(93)*n(idx_HCO)
      pdj(55) =  &
          +k(1000)  &
          +k(112)*n(idx_H)
      pdj(58) =  &
          +k(444)*n(idx_CH)  &
          +k(85)*n(idx_CH2)
      pdj(59) =  &
          +k(441)*n(idx_C)  &
          +k(87)*n(idx_CH)
      pdj(60) =  &
          +k(90)*n(idx_H2CO)
      pdj(62) =  &
          +k(95)*n(idx_NH3)
      pdj(63) =  &
          +k(97)*n(idx_NO)
      pdj(64) =  &
          +k(88)*n(idx_CN)
      pdj(65) =  &
          +k(89)*n(idx_CO)
      pdj(67) =  &
          +k(98)*n(idx_O2)
      pdj(68) =  &
          +k(91)*n(idx_H2O)  &
          +k(459)*n(idx_OH)
      pdj(69) =  &
          +k(94)*n(idx_NH2)  &
          +k(455)*n(idx_NH)
      pdj(71) =  &
          +k(99)*n(idx_OH)  &
          +k(458)*n(idx_O)
      pdj(72) =  &
          +k(443)*n(idx_CH4)  &
          +k(442)*n(idx_CH2)
      pdj(73) =  &
          +k(86)*n(idx_CH4)
      pdj(75) =  &
          +k(445)*n(idx_CN)  &
          +k(92)*n(idx_HCN)
      pdj(76) =  &
          +k(454)*n(idx_N)  &
          +k(96)*n(idx_NH)
      pdj(77) =  &
          -k(455)*n(idx_NH)  &
          -k(449)*n(idx_H2CO)  &
          -k(442)*n(idx_CH2)  &
          -k(444)*n(idx_CH)  &
          -k(90)*n(idx_H2CO)  &
          -k(452)*n(idx_HE)  &
          -k(91)*n(idx_H2O)  &
          -k(457)*n(idx_O2)  &
          -k(112)*n(idx_H)  &
          -k(445)*n(idx_CN)  &
          -k(86)*n(idx_CH4)  &
          -k(270)*n(idx_E)  &
          -k(92)*n(idx_HCN)  &
          -k(94)*n(idx_NH2)  &
          -k(447)*n(idx_CO)  &
          -k(1000)  &
          -k(99)*n(idx_OH)  &
          -k(451)*n(idx_HCO)  &
          -k(88)*n(idx_CN)  &
          -k(446)*n(idx_CO2)  &
          -k(441)*n(idx_C)  &
          -k(97)*n(idx_NO)  &
          -k(459)*n(idx_OH)  &
          -k(453)*n(idx_N2)  &
          -k(95)*n(idx_NH3)  &
          -k(454)*n(idx_N)  &
          -k(93)*n(idx_HCO)  &
          -k(458)*n(idx_O)  &
          -k(87)*n(idx_CH)  &
          -k(448)*n(idx_H2)  &
          -k(85)*n(idx_CH2)  &
          -k(98)*n(idx_O2)  &
          -k(443)*n(idx_CH4)  &
          -k(456)*n(idx_NO)  &
          -k(96)*n(idx_NH)  &
          -k(450)*n(idx_H2O)  &
          -k(89)*n(idx_CO)
      pdj(79) =  &
          +k(456)*n(idx_NO)
      pdj(81) =  &
          +k(448)*n(idx_H2)  &
          +k(451)*n(idx_HCO)
      pdj(83) =  &
          +k(450)*n(idx_H2O)
      pdj(85) =  &
          +k(446)*n(idx_CO2)
      pdj(86) =  &
          +k(452)*n(idx_HE)
      pdj(87) =  &
          +k(453)*n(idx_N2)
      pdj(88) =  &
          +k(457)*n(idx_O2)
    elseif(j==78) then
      pdj(1) =  &
          -k(1060)*n(idx_E)
      pdj(2) =  &
          -k(124)*n(idx_CH)  &
          +k(587)*n(idx_HCN)  &
          -k(573)*n(idx_CH)
      pdj(3) =  &
          +k(592)*n(idx_HCO)  &
          +k(576)*n(idx_CO2)  &
          +k(607)*n(idx_OCN)  &
          +k(580)*n(idx_CO)  &
          +k(605)*n(idx_NO)  &
          +k(606)*n(idx_O2)  &
          +k(583)*n(idx_H2CO)
      pdj(4) =  &
          -k(595)*n(idx_HNC)  &
          -k(593)*n(idx_HNC)  &
          -k(594)*n(idx_HNC)
      pdj(5) =  &
          -k(586)*n(idx_HCN)  &
          -k(587)*n(idx_HCN)  &
          -k(588)*n(idx_HCN)  &
          -k(589)*n(idx_HCN)
      pdj(6) =  &
          +k(564)*n(idx_CH2)  &
          +k(601)*n(idx_NH3)  &
          +k(566)*n(idx_CH3)  &
          -k(100)*n(idx_H2)  &
          +k(569)*n(idx_CH4)  &
          -k(468)*n(idx_H2)  &
          +k(599)*n(idx_NH2)  &
          +k(581)*n(idx_H2CO)  &
          +k(570)*n(idx_CH4)
      pdj(7) =  &
          -k(122)*n(idx_C)  &
          +k(595)*n(idx_HNC)  &
          +k(578)*n(idx_CO2)  &
          +k(574)*n(idx_CN)
      pdj(8) =  &
          -k(114)*n(idx_H)  &
          +k(596)*n(idx_HNO)  &
          +k(603)*n(idx_NH)  &
          +k(594)*n(idx_HNC)  &
          +k(609)*n(idx_OH)  &
          +k(602)*n(idx_NH3)  &
          +k(571)*n(idx_CH4)  &
          +k(600)*n(idx_NH2)  &
          +k(569)*n(idx_CH4)  &
          +k(584)*n(idx_H2O)  &
          +k(573)*n(idx_CH)  &
          +k(468)*n(idx_H2)  &
          +k(586)*n(idx_HCN)  &
          +k(593)*n(idx_HNC)  &
          +k(588)*n(idx_HCN)  &
          +k(565)*n(idx_CH2)  &
          +k(582)*n(idx_H2CO)  &
          +k(590)*n(idx_HCO)
      pdj(9) =  &
          -k(584)*n(idx_H2O)  &
          -k(126)*n(idx_H2O)  &
          -k(585)*n(idx_H2O)
      pdj(10) =  &
          +k(568)*n(idx_CH3OH)  &
          +k(585)*n(idx_H2O)  &
          -k(609)*n(idx_OH)
      pdj(11) =  &
          -k(606)*n(idx_O2)  &
          -k(129)*n(idx_O2)  &
          +k(579)*n(idx_CO2)
      pdj(12) =  &
          -k(564)*n(idx_CH2)  &
          -k(565)*n(idx_CH2)
      pdj(13) =  &
          -k(582)*n(idx_H2CO)  &
          -k(581)*n(idx_H2CO)  &
          -k(125)*n(idx_H2CO)  &
          -k(583)*n(idx_H2CO)
      pdj(14) =  &
          -k(590)*n(idx_HCO)  &
          -k(591)*n(idx_HCO)  &
          -k(592)*n(idx_HCO)
      pdj(16) =  &
          -k(128)*n(idx_NH3)  &
          -k(602)*n(idx_NH3)  &
          -k(601)*n(idx_NH3)
      pdj(17) =  &
          -k(604)*n(idx_NO)  &
          -k(605)*n(idx_NO)  &
          +k(597)*n(idx_HNO)
      pdj(18) =  &
          +k(608)*n(idx_OCN)  &
          -k(575)*n(idx_CN)  &
          -k(574)*n(idx_CN)
      pdj(19) =  &
          -k(580)*n(idx_CO)  &
          +k(591)*n(idx_HCO)  &
          +k(577)*n(idx_CO2)
      pdj(20) =  &
          -k(127)*n(idx_N2)  &
          -k(598)*n(idx_N2)
      pdj(21) =  &
          -k(600)*n(idx_NH2)  &
          -k(599)*n(idx_NH2)
      pdj(22) =  &
          -k(566)*n(idx_CH3)  &
          +k(567)*n(idx_CH3OH)  &
          +k(572)*n(idx_CH4)
      pdj(23) =  &
          -k(123)*n(idx_CH4)  &
          -k(571)*n(idx_CH4)  &
          -k(572)*n(idx_CH4)  &
          -k(570)*n(idx_CH4)  &
          -k(569)*n(idx_CH4)
      pdj(24) =  &
          +k(589)*n(idx_HCN)  &
          +k(575)*n(idx_CN)  &
          +k(588)*n(idx_HCN)  &
          +k(594)*n(idx_HNC)  &
          +k(598)*n(idx_N2)  &
          +k(604)*n(idx_NO)
      pdj(25) =  &
          -k(603)*n(idx_NH)
      pdj(26) =  &
          +k(596)*n(idx_HNO)  &
          +k(603)*n(idx_NH)  &
          +k(587)*n(idx_HCN)  &
          +k(574)*n(idx_CN)  &
          +k(122)*n(idx_C)  &
          +k(608)*n(idx_OCN)  &
          +k(576)*n(idx_CO2)  &
          +k(573)*n(idx_CH)  &
          +k(126)*n(idx_H2O)  &
          +k(575)*n(idx_CN)  &
          +k(599)*n(idx_NH2)  &
          +k(588)*n(idx_HCN)  &
          +k(577)*n(idx_CO2)  &
          +k(579)*n(idx_CO2)  &
          +k(592)*n(idx_HCO)  &
          +k(123)*n(idx_CH4)  &
          +k(601)*n(idx_NH3)  &
          +k(590)*n(idx_HCO)  &
          +k(602)*n(idx_NH3)  &
          +k(571)*n(idx_CH4)  &
          +k(589)*n(idx_HCN)  &
          +k(600)*n(idx_NH2)  &
          +k(568)*n(idx_CH3OH)  &
          +k(572)*n(idx_CH4)  &
          +k(468)*n(idx_H2)  &
          +k(604)*n(idx_NO)  &
          +k(583)*n(idx_H2CO)  &
          +k(129)*n(idx_O2)  &
          +k(580)*n(idx_CO)  &
          +k(570)*n(idx_CH4)  &
          +k(124)*n(idx_CH)  &
          +k(586)*n(idx_HCN)  &
          +k(114)*n(idx_H)  &
          +k(597)*n(idx_HNO)  &
          +k(125)*n(idx_H2CO)  &
          +k(609)*n(idx_OH)  &
          +k(607)*n(idx_OCN)  &
          +k(1060)*n(idx_E)  &
          +k(593)*n(idx_HNC)  &
          +k(100)*n(idx_H2)  &
          +k(127)*n(idx_N2)  &
          +k(565)*n(idx_CH2)  &
          +k(582)*n(idx_H2CO)  &
          +k(605)*n(idx_NO)  &
          +k(578)*n(idx_CO2)  &
          +k(564)*n(idx_CH2)  &
          +k(595)*n(idx_HNC)  &
          +k(566)*n(idx_CH3)  &
          +k(569)*n(idx_CH4)  &
          +k(584)*n(idx_H2O)  &
          +k(567)*n(idx_CH3OH)  &
          +k(585)*n(idx_H2O)  &
          +k(581)*n(idx_H2CO)  &
          +k(128)*n(idx_NH3)  &
          +k(594)*n(idx_HNC)  &
          +k(598)*n(idx_N2)  &
          +k(606)*n(idx_O2)
      pdj(27) =  &
          -k(597)*n(idx_HNO)  &
          -k(596)*n(idx_HNO)
      pdj(28) =  &
          -k(567)*n(idx_CH3OH)  &
          -k(568)*n(idx_CH3OH)
      pdj(29) =  &
          -k(578)*n(idx_CO2)  &
          -k(579)*n(idx_CO2)  &
          -k(576)*n(idx_CO2)  &
          -k(577)*n(idx_CO2)
      pdj(34) =  &
          -k(607)*n(idx_OCN)  &
          -k(608)*n(idx_OCN)
      pdj(54) =  &
          +k(582)*n(idx_H2CO)
      pdj(55) =  &
          +k(572)*n(idx_CH4)  &
          +k(114)*n(idx_H)  &
          +k(468)*n(idx_H2)  &
          +k(597)*n(idx_HNO)  &
          +k(585)*n(idx_H2O)
      pdj(57) =  &
          +k(564)*n(idx_CH2)  &
          +k(122)*n(idx_C)  &
          +k(573)*n(idx_CH)  &
          +k(575)*n(idx_CN)  &
          +k(588)*n(idx_HCN)  &
          +k(594)*n(idx_HNC)  &
          +k(580)*n(idx_CO)  &
          +k(579)*n(idx_CO2)
      pdj(58) =  &
          +k(583)*n(idx_H2CO)  &
          +k(570)*n(idx_CH4)
      pdj(59) =  &
          +k(592)*n(idx_HCO)  &
          +k(124)*n(idx_CH)  &
          +k(566)*n(idx_CH3)  &
          +k(589)*n(idx_HCN)  &
          +k(569)*n(idx_CH4)  &
          +k(565)*n(idx_CH2)
      pdj(60) =  &
          +k(125)*n(idx_H2CO)
      pdj(62) =  &
          +k(128)*n(idx_NH3)
      pdj(63) =  &
          +k(596)*n(idx_HNO)
      pdj(64) =  &
          +k(593)*n(idx_HNC)  &
          +k(586)*n(idx_HCN)  &
          +k(607)*n(idx_OCN)
      pdj(65) =  &
          +k(581)*n(idx_H2CO)  &
          +k(576)*n(idx_CO2)  &
          +k(590)*n(idx_HCO)
      pdj(66) =  &
          +k(127)*n(idx_N2)
      pdj(67) =  &
          +k(129)*n(idx_O2)  &
          +k(578)*n(idx_CO2)
      pdj(68) =  &
          +k(126)*n(idx_H2O)
      pdj(69) =  &
          +k(602)*n(idx_NH3)
      pdj(70) =  &
          +k(604)*n(idx_NO)  &
          +k(608)*n(idx_OCN)  &
          +k(577)*n(idx_CO2)  &
          +k(606)*n(idx_O2)  &
          +k(609)*n(idx_OH)
      pdj(71) =  &
          +k(584)*n(idx_H2O)  &
          +k(567)*n(idx_CH3OH)
      pdj(72) =  &
          +k(568)*n(idx_CH3OH)  &
          +k(571)*n(idx_CH4)
      pdj(73) =  &
          +k(123)*n(idx_CH4)
      pdj(74) =  &
          +k(603)*n(idx_NH)  &
          +k(587)*n(idx_HCN)  &
          +k(574)*n(idx_CN)  &
          +k(598)*n(idx_N2)  &
          +k(599)*n(idx_NH2)  &
          +k(605)*n(idx_NO)
      pdj(76) =  &
          +k(600)*n(idx_NH2)  &
          +k(601)*n(idx_NH3)  &
          +k(595)*n(idx_HNC)
      pdj(77) =  &
          +k(100)*n(idx_H2)
      pdj(78) =  &
          -k(114)*n(idx_H)  &
          -k(604)*n(idx_NO)  &
          -k(584)*n(idx_H2O)  &
          -k(591)*n(idx_HCO)  &
          -k(598)*n(idx_N2)  &
          -k(590)*n(idx_HCO)  &
          -k(576)*n(idx_CO2)  &
          -k(602)*n(idx_NH3)  &
          -k(583)*n(idx_H2CO)  &
          -k(125)*n(idx_H2CO)  &
          -k(129)*n(idx_O2)  &
          -k(124)*n(idx_CH)  &
          -k(593)*n(idx_HNC)  &
          -k(596)*n(idx_HNO)  &
          -k(580)*n(idx_CO)  &
          -k(570)*n(idx_CH4)  &
          -k(597)*n(idx_HNO)  &
          -k(577)*n(idx_CO2)  &
          -k(578)*n(idx_CO2)  &
          -k(123)*n(idx_CH4)  &
          -k(605)*n(idx_NO)  &
          -k(585)*n(idx_H2O)  &
          -k(564)*n(idx_CH2)  &
          -k(609)*n(idx_OH)  &
          -k(100)*n(idx_H2)  &
          -k(579)*n(idx_CO2)  &
          -k(581)*n(idx_H2CO)  &
          -k(1060)*n(idx_E)  &
          -k(468)*n(idx_H2)  &
          -k(608)*n(idx_OCN)  &
          -k(600)*n(idx_NH2)  &
          -k(565)*n(idx_CH2)  &
          -k(568)*n(idx_CH3OH)  &
          -k(572)*n(idx_CH4)  &
          -k(595)*n(idx_HNC)  &
          -k(588)*n(idx_HCN)  &
          -k(574)*n(idx_CN)  &
          -k(586)*n(idx_HCN)  &
          -k(126)*n(idx_H2O)  &
          -k(569)*n(idx_CH4)  &
          -k(587)*n(idx_HCN)  &
          -k(594)*n(idx_HNC)  &
          -k(582)*n(idx_H2CO)  &
          -k(573)*n(idx_CH)  &
          -k(589)*n(idx_HCN)  &
          -k(566)*n(idx_CH3)  &
          -k(122)*n(idx_C)  &
          -k(575)*n(idx_CN)  &
          -k(592)*n(idx_HCO)  &
          -k(128)*n(idx_NH3)  &
          -k(606)*n(idx_O2)  &
          -k(601)*n(idx_NH3)  &
          -k(599)*n(idx_NH2)  &
          -k(571)*n(idx_CH4)  &
          -k(607)*n(idx_OCN)  &
          -k(603)*n(idx_NH)  &
          -k(127)*n(idx_N2)  &
          -k(567)*n(idx_CH3OH)
      pdj(86) =  &
          +k(591)*n(idx_HCO)
    elseif(j==79) then
      pdj(1) =  &
          -k(299)*n(idx_E)
      pdj(2) =  &
          -k(409)*n(idx_CH)
      pdj(4) =  &
          -k(560)*n(idx_HNC)
      pdj(5) =  &
          -k(545)*n(idx_HCN)
      pdj(7) =  &
          -k(334)*n(idx_C)
      pdj(8) =  &
          +k(299)*n(idx_E)
      pdj(9) =  &
          -k(498)*n(idx_H2O)
      pdj(10) =  &
          -k(743)*n(idx_OH)
      pdj(12) =  &
          -k(374)*n(idx_CH2)
      pdj(13) =  &
          -k(480)*n(idx_H2CO)
      pdj(14) =  &
          -k(553)*n(idx_HCO)
      pdj(17) =  &
          -k(183)*n(idx_NO)  &
          +k(480)*n(idx_H2CO)  &
          +k(545)*n(idx_HCN)  &
          +k(498)*n(idx_H2O)  &
          +k(299)*n(idx_E)  &
          +k(553)*n(idx_HCO)  &
          +k(630)*n(idx_N2)  &
          +k(421)*n(idx_CN)  &
          +k(560)*n(idx_HNC)  &
          +k(374)*n(idx_CH2)  &
          +k(743)*n(idx_OH)  &
          +k(695)*n(idx_NH)  &
          +k(683)*n(idx_NH2)  &
          +k(425)*n(idx_CO)  &
          +k(563)*n(idx_CO2)  &
          +k(409)*n(idx_CH)  &
          +k(334)*n(idx_C)
      pdj(18) =  &
          -k(421)*n(idx_CN)
      pdj(19) =  &
          -k(425)*n(idx_CO)
      pdj(20) =  &
          -k(630)*n(idx_N2)
      pdj(21) =  &
          -k(683)*n(idx_NH2)
      pdj(25) =  &
          -k(695)*n(idx_NH)
      pdj(27) =  &
          +k(183)*n(idx_NO)
      pdj(29) =  &
          -k(563)*n(idx_CO2)
      pdj(48) =  &
          +k(1115)
      pdj(54) =  &
          +k(425)*n(idx_CO)
      pdj(58) =  &
          +k(409)*n(idx_CH)
      pdj(59) =  &
          +k(334)*n(idx_C)
      pdj(60) =  &
          +k(553)*n(idx_HCO)
      pdj(62) =  &
          +k(683)*n(idx_NH2)
      pdj(63) =  &
          +k(183)*n(idx_NO)
      pdj(68) =  &
          +k(743)*n(idx_OH)
      pdj(69) =  &
          +k(695)*n(idx_NH)
      pdj(72) =  &
          +k(374)*n(idx_CH2)
      pdj(75) =  &
          +k(421)*n(idx_CN)
      pdj(79) =  &
          -k(183)*n(idx_NO)  &
          -k(421)*n(idx_CN)  &
          -k(425)*n(idx_CO)  &
          -k(695)*n(idx_NH)  &
          -k(545)*n(idx_HCN)  &
          -k(299)*n(idx_E)  &
          -k(553)*n(idx_HCO)  &
          -k(560)*n(idx_HNC)  &
          -k(1115)  &
          -k(498)*n(idx_H2O)  &
          -k(409)*n(idx_CH)  &
          -k(334)*n(idx_C)  &
          -k(563)*n(idx_CO2)  &
          -k(374)*n(idx_CH2)  &
          -k(480)*n(idx_H2CO)  &
          -k(630)*n(idx_N2)  &
          -k(743)*n(idx_OH)  &
          -k(683)*n(idx_NH2)
      pdj(82) =  &
          +k(480)*n(idx_H2CO)
      pdj(83) =  &
          +k(498)*n(idx_H2O)
      pdj(84) =  &
          +k(545)*n(idx_HCN)  &
          +k(560)*n(idx_HNC)
      pdj(85) =  &
          +k(563)*n(idx_CO2)
      pdj(87) =  &
          +k(630)*n(idx_N2)
    elseif(j==80) then
      pdj(1) =  &
          -k(276)*n(idx_E)  &
          -k(275)*n(idx_E)
      pdj(6) =  &
          +k(276)*n(idx_E)
      pdj(8) =  &
          +k(275)*n(idx_E)  &
          +k(1116)
      pdj(17) =  &
          +k(276)*n(idx_E)
      pdj(27) =  &
          +k(275)*n(idx_E)
      pdj(48) =  &
          +k(1116)
      pdj(80) =  &
          -k(276)*n(idx_E)  &
          -k(1116)  &
          -k(275)*n(idx_E)
    elseif(j==81) then
      pdj(1) =  &
          -k(281)*n(idx_E)  &
          -k(280)*n(idx_E)
      pdj(2) =  &
          -k(506)*n(idx_CH)
      pdj(3) =  &
          -k(524)*n(idx_O)  &
          -k(525)*n(idx_O)
      pdj(4) =  &
          -k(515)*n(idx_HNC)
      pdj(5) =  &
          -k(513)*n(idx_HCN)
      pdj(6) =  &
          +k(512)*n(idx_H2O)  &
          +k(514)*n(idx_HCO)  &
          +k(506)*n(idx_CH)  &
          +k(508)*n(idx_CO2)  &
          +k(523)*n(idx_O2)  &
          +k(522)*n(idx_NO)  &
          +k(511)*n(idx_H2CO)  &
          +k(509)*n(idx_CO)  &
          +k(505)*n(idx_CH3OH)  &
          +k(516)*n(idx_HNO)  &
          +k(502)*n(idx_C)  &
          +k(503)*n(idx_CH2)  &
          +k(280)*n(idx_E)  &
          +k(521)*n(idx_NO2)  &
          +k(525)*n(idx_O)  &
          +k(520)*n(idx_NH)  &
          +k(510)*n(idx_CO)  &
          +k(507)*n(idx_CN)  &
          +k(517)*n(idx_MG)  &
          +k(519)*n(idx_NH2)  &
          +k(504)*n(idx_CH3)  &
          +k(1010)  &
          +k(515)*n(idx_HNC)  &
          +k(513)*n(idx_HCN)  &
          +k(526)*n(idx_OH)  &
          +k(518)*n(idx_N2)
      pdj(7) =  &
          -k(502)*n(idx_C)
      pdj(8) =  &
          +3.d0*k(281)*n(idx_E)  &
          +k(280)*n(idx_E)  &
          +k(524)*n(idx_O)  &
          +k(517)*n(idx_MG)  &
          +k(1009)
      pdj(9) =  &
          -k(512)*n(idx_H2O)  &
          +k(505)*n(idx_CH3OH)
      pdj(10) =  &
          -k(526)*n(idx_OH)  &
          +k(521)*n(idx_NO2)
      pdj(11) =  &
          -k(523)*n(idx_O2)
      pdj(12) =  &
          -k(503)*n(idx_CH2)
      pdj(13) =  &
          -k(511)*n(idx_H2CO)
      pdj(14) =  &
          -k(514)*n(idx_HCO)
      pdj(15) =  &
          -k(517)*n(idx_MG)
      pdj(17) =  &
          -k(522)*n(idx_NO)
      pdj(18) =  &
          -k(507)*n(idx_CN)
      pdj(19) =  &
          -k(510)*n(idx_CO)  &
          -k(509)*n(idx_CO)
      pdj(20) =  &
          -k(518)*n(idx_N2)
      pdj(21) =  &
          -k(519)*n(idx_NH2)
      pdj(22) =  &
          -k(504)*n(idx_CH3)
      pdj(25) =  &
          -k(520)*n(idx_NH)
      pdj(27) =  &
          -k(516)*n(idx_HNO)
      pdj(28) =  &
          -k(505)*n(idx_CH3OH)
      pdj(29) =  &
          -k(508)*n(idx_CO2)
      pdj(32) =  &
          -k(521)*n(idx_NO2)
      pdj(54) =  &
          +k(509)*n(idx_CO)
      pdj(55) =  &
          +k(1010)
      pdj(56) =  &
          +k(510)*n(idx_CO)
      pdj(58) =  &
          +k(506)*n(idx_CH)
      pdj(59) =  &
          +k(502)*n(idx_C)
      pdj(60) =  &
          +k(514)*n(idx_HCO)
      pdj(61) =  &
          +k(517)*n(idx_MG)
      pdj(62) =  &
          +k(519)*n(idx_NH2)
      pdj(63) =  &
          +k(521)*n(idx_NO2)
      pdj(68) =  &
          +k(524)*n(idx_O)  &
          +k(526)*n(idx_OH)
      pdj(69) =  &
          +k(520)*n(idx_NH)
      pdj(71) =  &
          +k(525)*n(idx_O)
      pdj(72) =  &
          +k(503)*n(idx_CH2)  &
          +k(505)*n(idx_CH3OH)
      pdj(73) =  &
          +k(504)*n(idx_CH3)
      pdj(75) =  &
          +k(507)*n(idx_CN)
      pdj(77) =  &
          +k(1009)
      pdj(79) =  &
          +k(522)*n(idx_NO)
      pdj(80) =  &
          +k(516)*n(idx_HNO)
      pdj(81) =  &
          -k(507)*n(idx_CN)  &
          -k(506)*n(idx_CH)  &
          -k(515)*n(idx_HNC)  &
          -k(513)*n(idx_HCN)  &
          -k(524)*n(idx_O)  &
          -k(280)*n(idx_E)  &
          -k(520)*n(idx_NH)  &
          -k(505)*n(idx_CH3OH)  &
          -k(512)*n(idx_H2O)  &
          -k(1010)  &
          -k(516)*n(idx_HNO)  &
          -k(523)*n(idx_O2)  &
          -k(281)*n(idx_E)  &
          -k(522)*n(idx_NO)  &
          -k(503)*n(idx_CH2)  &
          -k(1009)  &
          -k(511)*n(idx_H2CO)  &
          -k(509)*n(idx_CO)  &
          -k(521)*n(idx_NO2)  &
          -k(517)*n(idx_MG)  &
          -k(519)*n(idx_NH2)  &
          -k(514)*n(idx_HCO)  &
          -k(525)*n(idx_O)  &
          -k(526)*n(idx_OH)  &
          -k(510)*n(idx_CO)  &
          -k(508)*n(idx_CO2)  &
          -k(502)*n(idx_C)  &
          -k(504)*n(idx_CH3)  &
          -k(518)*n(idx_N2)
      pdj(82) =  &
          +k(511)*n(idx_H2CO)
      pdj(83) =  &
          +k(512)*n(idx_H2O)
      pdj(84) =  &
          +k(513)*n(idx_HCN)  &
          +k(515)*n(idx_HNC)
      pdj(85) =  &
          +k(508)*n(idx_CO2)
      pdj(87) =  &
          +k(518)*n(idx_N2)
      pdj(88) =  &
          +k(523)*n(idx_O2)
    elseif(j==82) then
      pdj(1) =  &
          -k(286)*n(idx_E)  &
          -k(283)*n(idx_E)  &
          -k(282)*n(idx_E)  &
          -k(285)*n(idx_E)  &
          -k(284)*n(idx_E)
      pdj(2) =  &
          -k(403)*n(idx_CH)  &
          +k(283)*n(idx_E)
      pdj(4) =  &
          -k(558)*n(idx_HNC)
      pdj(5) =  &
          -k(543)*n(idx_HCN)
      pdj(6) =  &
          +k(284)*n(idx_E)
      pdj(8) =  &
          +2.d0*k(286)*n(idx_E)  &
          +k(285)*n(idx_E)  &
          +k(284)*n(idx_E)
      pdj(9) =  &
          -k(494)*n(idx_H2O)  &
          +k(283)*n(idx_E)
      pdj(10) =  &
          +k(282)*n(idx_E)
      pdj(12) =  &
          +k(282)*n(idx_E)
      pdj(13) =  &
          +k(558)*n(idx_HNC)  &
          +k(543)*n(idx_HCN)  &
          +k(285)*n(idx_E)  &
          +k(403)*n(idx_CH)  &
          +k(494)*n(idx_H2O)  &
          +k(677)*n(idx_NH2)
      pdj(14) =  &
          +k(286)*n(idx_E)
      pdj(19) =  &
          +k(284)*n(idx_E)
      pdj(21) =  &
          -k(677)*n(idx_NH2)
      pdj(35) =  &
          +k(1069)
      pdj(58) =  &
          +k(403)*n(idx_CH)
      pdj(62) =  &
          +k(677)*n(idx_NH2)
      pdj(82) =  &
          -k(543)*n(idx_HCN)  &
          -k(677)*n(idx_NH2)  &
          -k(494)*n(idx_H2O)  &
          -k(403)*n(idx_CH)  &
          -k(558)*n(idx_HNC)  &
          -k(285)*n(idx_E)  &
          -k(1069)  &
          -k(286)*n(idx_E)  &
          -k(283)*n(idx_E)  &
          -k(282)*n(idx_E)  &
          -k(284)*n(idx_E)
      pdj(83) =  &
          +k(494)*n(idx_H2O)
      pdj(84) =  &
          +k(558)*n(idx_HNC)  &
          +k(543)*n(idx_HCN)
    elseif(j==83) then
      pdj(1) =  &
          -k(287)*n(idx_E)  &
          -k(290)*n(idx_E)  &
          -k(289)*n(idx_E)  &
          -k(288)*n(idx_E)
      pdj(2) =  &
          -k(404)*n(idx_CH)
      pdj(3) =  &
          +k(288)*n(idx_E)
      pdj(4) =  &
          -k(529)*n(idx_HNC)
      pdj(5) =  &
          -k(528)*n(idx_HCN)
      pdj(6) =  &
          +k(289)*n(idx_E)  &
          +k(288)*n(idx_E)  &
          +k(330)*n(idx_C)
      pdj(7) =  &
          -k(330)*n(idx_C)
      pdj(8) =  &
          +k(288)*n(idx_E)  &
          +k(287)*n(idx_E)  &
          +k(1106)  &
          +2.d0*k(290)*n(idx_E)
      pdj(9) =  &
          +k(528)*n(idx_HCN)  &
          +k(678)*n(idx_NH2)  &
          +k(287)*n(idx_E)  &
          +k(527)*n(idx_H2CO)  &
          +k(369)*n(idx_CH2)  &
          +k(529)*n(idx_HNC)  &
          +k(404)*n(idx_CH)
      pdj(10) =  &
          +k(289)*n(idx_E)  &
          +k(290)*n(idx_E)
      pdj(12) =  &
          -k(369)*n(idx_CH2)
      pdj(13) =  &
          -k(527)*n(idx_H2CO)
      pdj(21) =  &
          -k(678)*n(idx_NH2)
      pdj(40) =  &
          +k(1106)
      pdj(54) =  &
          +k(330)*n(idx_C)
      pdj(58) =  &
          +k(404)*n(idx_CH)
      pdj(62) =  &
          +k(678)*n(idx_NH2)
      pdj(72) =  &
          +k(369)*n(idx_CH2)
      pdj(82) =  &
          +k(527)*n(idx_H2CO)
      pdj(83) =  &
          -k(290)*n(idx_E)  &
          -k(529)*n(idx_HNC)  &
          -k(528)*n(idx_HCN)  &
          -k(369)*n(idx_CH2)  &
          -k(1106)  &
          -k(678)*n(idx_NH2)  &
          -k(330)*n(idx_C)  &
          -k(527)*n(idx_H2CO)  &
          -k(289)*n(idx_E)  &
          -k(287)*n(idx_E)  &
          -k(404)*n(idx_CH)  &
          -k(288)*n(idx_E)
      pdj(84) =  &
          +k(528)*n(idx_HCN)  &
          +k(529)*n(idx_HNC)
    elseif(j==84) then
      pdj(1) =  &
          -k(294)*n(idx_E)  &
          -k(292)*n(idx_E)  &
          -k(293)*n(idx_E)
      pdj(2) =  &
          -k(406)*n(idx_CH)  &
          -k(407)*n(idx_CH)
      pdj(4) =  &
          +k(549)*n(idx_H2CO)  &
          +k(681)*n(idx_NH2)  &
          +k(294)*n(idx_E)  &
          +k(372)*n(idx_CH2)  &
          +k(407)*n(idx_CH)
      pdj(5) =  &
          +k(293)*n(idx_E)  &
          +k(371)*n(idx_CH2)  &
          +k(680)*n(idx_NH2)  &
          +k(548)*n(idx_H2CO)  &
          +k(406)*n(idx_CH)
      pdj(8) =  &
          +2.d0*k(292)*n(idx_E)  &
          +k(293)*n(idx_E)  &
          +k(294)*n(idx_E)  &
          +k(1121)
      pdj(12) =  &
          -k(372)*n(idx_CH2)  &
          -k(371)*n(idx_CH2)
      pdj(13) =  &
          -k(548)*n(idx_H2CO)  &
          -k(549)*n(idx_H2CO)
      pdj(18) =  &
          +k(292)*n(idx_E)
      pdj(21) =  &
          -k(680)*n(idx_NH2)  &
          -k(681)*n(idx_NH2)
      pdj(44) =  &
          +k(1121)
      pdj(58) =  &
          +k(407)*n(idx_CH)  &
          +k(406)*n(idx_CH)
      pdj(62) =  &
          +k(681)*n(idx_NH2)  &
          +k(680)*n(idx_NH2)
      pdj(72) =  &
          +k(371)*n(idx_CH2)  &
          +k(372)*n(idx_CH2)
      pdj(82) =  &
          +k(548)*n(idx_H2CO)  &
          +k(549)*n(idx_H2CO)
      pdj(84) =  &
          -k(372)*n(idx_CH2)  &
          -k(294)*n(idx_E)  &
          -k(293)*n(idx_E)  &
          -k(1121)  &
          -k(407)*n(idx_CH)  &
          -k(371)*n(idx_CH2)  &
          -k(406)*n(idx_CH)  &
          -k(548)*n(idx_H2CO)  &
          -k(681)*n(idx_NH2)  &
          -k(549)*n(idx_H2CO)  &
          -k(680)*n(idx_NH2)  &
          -k(292)*n(idx_E)
    elseif(j==85) then
      pdj(1) =  &
          -k(298)*n(idx_E)  &
          -k(296)*n(idx_E)  &
          -k(297)*n(idx_E)
      pdj(3) =  &
          -k(719)*n(idx_O)  &
          +k(297)*n(idx_E)
      pdj(7) =  &
          -k(333)*n(idx_C)
      pdj(8) =  &
          +k(1107)  &
          +k(296)*n(idx_E)  &
          +k(297)*n(idx_E)
      pdj(9) =  &
          -k(497)*n(idx_H2O)
      pdj(10) =  &
          +k(298)*n(idx_E)
      pdj(11) =  &
          +k(719)*n(idx_O)
      pdj(19) =  &
          -k(424)*n(idx_CO)  &
          +k(297)*n(idx_E)  &
          +k(298)*n(idx_E)
      pdj(29) =  &
          +k(333)*n(idx_C)  &
          +k(497)*n(idx_H2O)  &
          +k(424)*n(idx_CO)  &
          +k(296)*n(idx_E)
      pdj(42) =  &
          +k(1107)
      pdj(54) =  &
          +k(424)*n(idx_CO)  &
          +k(719)*n(idx_O)
      pdj(59) =  &
          +k(333)*n(idx_C)
      pdj(83) =  &
          +k(497)*n(idx_H2O)
      pdj(85) =  &
          -k(296)*n(idx_E)  &
          -k(298)*n(idx_E)  &
          -k(424)*n(idx_CO)  &
          -k(1107)  &
          -k(297)*n(idx_E)  &
          -k(719)*n(idx_O)  &
          -k(333)*n(idx_C)  &
          -k(497)*n(idx_H2O)
    elseif(j==86) then
      pdj(1) =  &
          -k(301)*n(idx_E)
      pdj(6) =  &
          -k(469)*n(idx_H2)
      pdj(8) =  &
          +k(301)*n(idx_E)  &
          -k(534)*n(idx_H)
      pdj(26) =  &
          +k(534)*n(idx_H)  &
          +k(301)*n(idx_E)  &
          +k(469)*n(idx_H2)
      pdj(77) =  &
          +k(534)*n(idx_H)
      pdj(81) =  &
          +k(469)*n(idx_H2)
      pdj(86) =  &
          -k(534)*n(idx_H)  &
          -k(301)*n(idx_E)  &
          -k(469)*n(idx_H2)
    elseif(j==87) then
      pdj(1) =  &
          -k(303)*n(idx_E)  &
          -k(304)*n(idx_E)
      pdj(2) =  &
          -k(411)*n(idx_CH)
      pdj(3) =  &
          -k(721)*n(idx_O)
      pdj(4) =  &
          -k(561)*n(idx_HNC)
      pdj(5) =  &
          -k(546)*n(idx_HCN)
      pdj(7) =  &
          -k(335)*n(idx_C)
      pdj(8) =  &
          +k(1122)  &
          +k(303)*n(idx_E)
      pdj(9) =  &
          -k(500)*n(idx_H2O)
      pdj(10) =  &
          -k(744)*n(idx_OH)
      pdj(12) =  &
          -k(375)*n(idx_CH2)
      pdj(13) =  &
          -k(633)*n(idx_H2CO)
      pdj(14) =  &
          -k(554)*n(idx_HCO)
      pdj(19) =  &
          -k(426)*n(idx_CO)
      pdj(20) =  &
          +k(684)*n(idx_NH2)  &
          +k(633)*n(idx_H2CO)  &
          +k(335)*n(idx_C)  &
          +k(561)*n(idx_HNC)  &
          +k(375)*n(idx_CH2)  &
          +k(411)*n(idx_CH)  &
          +k(546)*n(idx_HCN)  &
          +k(500)*n(idx_H2O)  &
          +k(721)*n(idx_O)  &
          +k(744)*n(idx_OH)  &
          +k(632)*n(idx_CO2)  &
          +k(696)*n(idx_NH)  &
          +k(426)*n(idx_CO)  &
          +k(554)*n(idx_HCO)  &
          +k(303)*n(idx_E)
      pdj(21) =  &
          -k(684)*n(idx_NH2)
      pdj(24) =  &
          +k(304)*n(idx_E)
      pdj(25) =  &
          +k(304)*n(idx_E)  &
          -k(696)*n(idx_NH)
      pdj(29) =  &
          -k(632)*n(idx_CO2)
      pdj(43) =  &
          +k(1122)
      pdj(54) =  &
          +k(426)*n(idx_CO)
      pdj(58) =  &
          +k(411)*n(idx_CH)
      pdj(59) =  &
          +k(335)*n(idx_C)
      pdj(60) =  &
          +k(554)*n(idx_HCO)
      pdj(62) =  &
          +k(684)*n(idx_NH2)
      pdj(68) =  &
          +k(744)*n(idx_OH)
      pdj(69) =  &
          +k(696)*n(idx_NH)
      pdj(71) =  &
          +k(721)*n(idx_O)
      pdj(72) =  &
          +k(375)*n(idx_CH2)
      pdj(82) =  &
          +k(633)*n(idx_H2CO)
      pdj(83) =  &
          +k(500)*n(idx_H2O)
      pdj(84) =  &
          +k(546)*n(idx_HCN)  &
          +k(561)*n(idx_HNC)
      pdj(85) =  &
          +k(632)*n(idx_CO2)
      pdj(87) =  &
          -k(375)*n(idx_CH2)  &
          -k(554)*n(idx_HCO)  &
          -k(335)*n(idx_C)  &
          -k(304)*n(idx_E)  &
          -k(426)*n(idx_CO)  &
          -k(561)*n(idx_HNC)  &
          -k(633)*n(idx_H2CO)  &
          -k(411)*n(idx_CH)  &
          -k(303)*n(idx_E)  &
          -k(546)*n(idx_HCN)  &
          -k(721)*n(idx_O)  &
          -k(696)*n(idx_NH)  &
          -k(500)*n(idx_H2O)  &
          -k(632)*n(idx_CO2)  &
          -k(1122)  &
          -k(684)*n(idx_NH2)  &
          -k(744)*n(idx_OH)
    elseif(j==88) then
      pdj(1) =  &
          -k(312)*n(idx_E)
      pdj(2) =  &
          -k(416)*n(idx_CH)
      pdj(3) =  &
          -k(724)*n(idx_O)
      pdj(4) =  &
          -k(562)*n(idx_HNC)
      pdj(5) =  &
          -k(547)*n(idx_HCN)
      pdj(6) =  &
          -k(476)*n(idx_H2)
      pdj(7) =  &
          -k(338)*n(idx_C)
      pdj(8) =  &
          +k(312)*n(idx_E)
      pdj(9) =  &
          -k(501)*n(idx_H2O)
      pdj(10) =  &
          -k(745)*n(idx_OH)
      pdj(11) =  &
          +k(556)*n(idx_HCO)  &
          +k(501)*n(idx_H2O)  &
          +k(702)*n(idx_NO)  &
          +k(685)*n(idx_NH2)  &
          +k(422)*n(idx_CN)  &
          +k(482)*n(idx_H2CO)  &
          +k(476)*n(idx_H2)  &
          +k(338)*n(idx_C)  &
          +k(700)*n(idx_NH)  &
          +k(631)*n(idx_N2)  &
          +k(380)*n(idx_CH2)  &
          +k(562)*n(idx_HNC)  &
          +k(312)*n(idx_E)  &
          +k(745)*n(idx_OH)  &
          +k(716)*n(idx_CO2)  &
          +k(547)*n(idx_HCN)  &
          +k(724)*n(idx_O)  &
          +k(427)*n(idx_CO)  &
          +k(416)*n(idx_CH)
      pdj(12) =  &
          -k(380)*n(idx_CH2)
      pdj(13) =  &
          -k(482)*n(idx_H2CO)
      pdj(14) =  &
          -k(556)*n(idx_HCO)
      pdj(17) =  &
          -k(702)*n(idx_NO)
      pdj(18) =  &
          -k(422)*n(idx_CN)
      pdj(19) =  &
          -k(427)*n(idx_CO)
      pdj(20) =  &
          -k(631)*n(idx_N2)
      pdj(21) =  &
          -k(685)*n(idx_NH2)
      pdj(25) =  &
          -k(700)*n(idx_NH)
      pdj(29) =  &
          -k(716)*n(idx_CO2)
      pdj(49) =  &
          +k(1117)
      pdj(54) =  &
          +k(427)*n(idx_CO)
      pdj(58) =  &
          +k(416)*n(idx_CH)
      pdj(59) =  &
          +k(338)*n(idx_C)
      pdj(60) =  &
          +k(556)*n(idx_HCO)
      pdj(62) =  &
          +k(685)*n(idx_NH2)
      pdj(68) =  &
          +k(745)*n(idx_OH)
      pdj(69) =  &
          +k(700)*n(idx_NH)
      pdj(71) =  &
          +k(724)*n(idx_O)
      pdj(72) =  &
          +k(380)*n(idx_CH2)
      pdj(75) =  &
          +k(422)*n(idx_CN)
      pdj(79) =  &
          +k(702)*n(idx_NO)
      pdj(81) =  &
          +k(476)*n(idx_H2)
      pdj(82) =  &
          +k(482)*n(idx_H2CO)
      pdj(83) =  &
          +k(501)*n(idx_H2O)
      pdj(84) =  &
          +k(562)*n(idx_HNC)  &
          +k(547)*n(idx_HCN)
      pdj(85) =  &
          +k(716)*n(idx_CO2)
      pdj(87) =  &
          +k(631)*n(idx_N2)
      pdj(88) =  &
          -k(422)*n(idx_CN)  &
          -k(416)*n(idx_CH)  &
          -k(685)*n(idx_NH2)  &
          -k(702)*n(idx_NO)  &
          -k(700)*n(idx_NH)  &
          -k(501)*n(idx_H2O)  &
          -k(380)*n(idx_CH2)  &
          -k(476)*n(idx_H2)  &
          -k(745)*n(idx_OH)  &
          -k(562)*n(idx_HNC)  &
          -k(312)*n(idx_E)  &
          -k(338)*n(idx_C)  &
          -k(724)*n(idx_O)  &
          -k(427)*n(idx_CO)  &
          -k(631)*n(idx_N2)  &
          -k(1117)  &
          -k(547)*n(idx_HCN)  &
          -k(716)*n(idx_CO2)  &
          -k(556)*n(idx_HCO)  &
          -k(482)*n(idx_H2CO)
    elseif(j==89) then
    elseif(j==90) then
    elseif(j==91) then

    elseif(j==92) then
    end if

    return
  end subroutine jes

  !*************************
  subroutine jex(neq,t,n,ml,mu,pd,npd)
    use krome_commons
    use krome_tabs
    use krome_cooling
    use krome_heating
    use krome_constants
    use krome_subs
    use krome_gadiab
    implicit none
    real*8::n(neq),pd(neq,neq),t,k(nrea),dn0,dn1,dnn,Tgas
    real*8::krome_gamma,nn(neq),nH2dust
    integer::neq,ml,mu,npd

    Tgas = n(idx_Tgas)
    npd = neq
    k(:) = coe_tab(n(:))
    pd(:,:) = 0d0
    krome_gamma = gamma_index(n(:))

    !d[E_dot]/d[E]
    pd(1,1) =  &
        -k(310)*n(idx_NOj)  &
        -k(1058)*n(idx_Hj)  &
        -k(1063)*n(idx_Oj)  &
        -k(291)*n(idx_HCNj)  &
        -k(299)*n(idx_HNOj)  &
        -k(286)*n(idx_H3COj)  &
        -k(292)*n(idx_HCNHj)  &
        -k(297)*n(idx_HCO2j)  &
        -k(270)*n(idx_H2j)  &
        -k(275)*n(idx_H2NOj)  &
        -k(308)*n(idx_NH3j)  &
        -k(263)*n(idx_CH3j)  &
        -k(1059)*n(idx_H2COj)  &
        -k(264)*n(idx_CH3j)  &
        -k(271)*n(idx_H2COj)  &
        -k(296)*n(idx_HCO2j)  &
        -k(301)*n(idx_HEHj)  &
        -k(302)*n(idx_N2j)  &
        -k(260)*n(idx_CH2j)  &
        -k(278)*n(idx_H2Oj)  &
        -k(1061)*n(idx_MGj)  &
        -k(273)*n(idx_H2COj)  &
        -k(284)*n(idx_H3COj)  &
        -k(287)*n(idx_H3Oj)  &
        -k(9)*n(idx_H2)  &
        -k(1062)*n(idx_Nj)  &
        -k(312)*n(idx_O2Hj)  &
        -k(289)*n(idx_H3Oj)  &
        -k(306)*n(idx_NH2j)  &
        -k(288)*n(idx_H3Oj)  &
        +k(9)*n(idx_H2)  &
        -k(1126)  &
        -k(309)*n(idx_NH3j)  &
        -k(262)*n(idx_CH2j)  &
        -k(277)*n(idx_H2Oj)  &
        -k(282)*n(idx_H3COj)  &
        -k(280)*n(idx_H3j)  &
        -k(313)*n(idx_OHj)  &
        -k(295)*n(idx_HCOj)  &
        -k(266)*n(idx_CH4j)  &
        -k(279)*n(idx_H2Oj)  &
        -k(285)*n(idx_H3COj)  &
        -k(1057)*n(idx_CH3j)  &
        -k(290)*n(idx_H3Oj)  &
        -k(305)*n(idx_NHj)  &
        -k(261)*n(idx_CH2j)  &
        -k(304)*n(idx_N2Hj)  &
        -k(300)*n(idx_HOCj)  &
        -k(293)*n(idx_HCNHj)  &
        -k(307)*n(idx_NH2j)  &
        -k(268)*n(idx_CNj)  &
        -k(1060)*n(idx_HEj)  &
        -k(276)*n(idx_H2NOj)  &
        -k(259)*n(idx_CHj)  &
        -k(283)*n(idx_H3COj)  &
        -k(267)*n(idx_CH4j)  &
        -k(269)*n(idx_COj)  &
        -k(274)*n(idx_H2COj)  &
        -k(294)*n(idx_HCNHj)  &
        -k(311)*n(idx_O2j)  &
        -k(1056)*n(idx_Cj)  &
        -k(272)*n(idx_H2COj)  &
        -k(281)*n(idx_H3j)  &
        -k(265)*n(idx_CH3j)  &
        -k(298)*n(idx_HCO2j)  &
        -k(303)*n(idx_N2Hj)

    !d[CH_dot]/d[E]
    pd(2,1) =  &
        +k(262)*n(idx_CH2j)  &
        +k(265)*n(idx_CH3j)  &
        +k(264)*n(idx_CH3j)  &
        +k(283)*n(idx_H3COj)

    !d[O_dot]/d[E]
    pd(3,1) =  &
        +k(269)*n(idx_COj)  &
        +k(313)*n(idx_OHj)  &
        +k(297)*n(idx_HCO2j)  &
        +k(288)*n(idx_H3Oj)  &
        +k(1063)*n(idx_Oj)  &
        +k(271)*n(idx_H2COj)  &
        +k(277)*n(idx_H2Oj)  &
        +2.d0*k(311)*n(idx_O2j)  &
        +k(278)*n(idx_H2Oj)  &
        +k(310)*n(idx_NOj)

    !d[HNC_dot]/d[E]
    pd(4,1) =  &
        +k(294)*n(idx_HCNHj)

    !d[HCN_dot]/d[E]
    pd(5,1) =  &
        +k(293)*n(idx_HCNHj)

    !d[H2_dot]/d[E]
    pd(6,1) =  &
        +k(276)*n(idx_H2NOj)  &
        +k(260)*n(idx_CH2j)  &
        +k(280)*n(idx_H3j)  &
        +k(288)*n(idx_H3Oj)  &
        +k(272)*n(idx_H2COj)  &
        -k(9)*n(idx_H2)  &
        +k(284)*n(idx_H3COj)  &
        +k(277)*n(idx_H2Oj)  &
        +k(289)*n(idx_H3Oj)  &
        +k(264)*n(idx_CH3j)

    !d[C_dot]/d[E]
    pd(7,1) =  &
        +k(260)*n(idx_CH2j)  &
        +k(269)*n(idx_COj)  &
        +k(1056)*n(idx_Cj)  &
        +k(261)*n(idx_CH2j)  &
        +k(268)*n(idx_CNj)  &
        +k(259)*n(idx_CHj)

    !d[H_dot]/d[E]
    pd(8,1) =  &
        +2.d0*k(265)*n(idx_CH3j)  &
        +k(1058)*n(idx_Hj)  &
        +2.d0*k(261)*n(idx_CH2j)  &
        +k(267)*n(idx_CH4j)  &
        +k(275)*n(idx_H2NOj)  &
        +k(279)*n(idx_H2Oj)  &
        +k(308)*n(idx_NH3j)  &
        +k(284)*n(idx_H3COj)  &
        +k(287)*n(idx_H3Oj)  &
        +k(300)*n(idx_HOCj)  &
        +2.d0*k(266)*n(idx_CH4j)  &
        +k(259)*n(idx_CHj)  &
        +k(312)*n(idx_O2Hj)  &
        +2.d0*k(273)*n(idx_H2COj)  &
        +2.d0*k(292)*n(idx_HCNHj)  &
        +k(299)*n(idx_HNOj)  &
        +2.d0*k(286)*n(idx_H3COj)  &
        +k(294)*n(idx_HCNHj)  &
        +3.d0*k(281)*n(idx_H3j)  &
        +k(285)*n(idx_H3COj)  &
        +2.d0*k(9)*n(idx_H2)  &
        +2.d0*k(306)*n(idx_NH2j)  &
        +k(301)*n(idx_HEHj)  &
        +k(263)*n(idx_CH3j)  &
        +k(297)*n(idx_HCO2j)  &
        +k(280)*n(idx_H3j)  &
        +2.d0*k(290)*n(idx_H3Oj)  &
        +k(288)*n(idx_H3Oj)  &
        +k(305)*n(idx_NHj)  &
        +k(274)*n(idx_H2COj)  &
        +k(296)*n(idx_HCO2j)  &
        +k(307)*n(idx_NH2j)  &
        +k(262)*n(idx_CH2j)  &
        +2.d0*k(309)*n(idx_NH3j)  &
        +k(303)*n(idx_N2Hj)  &
        +k(293)*n(idx_HCNHj)  &
        +k(313)*n(idx_OHj)  &
        +2.d0*k(270)*n(idx_H2j)  &
        +k(295)*n(idx_HCOj)  &
        +k(291)*n(idx_HCNj)  &
        +2.d0*k(278)*n(idx_H2Oj)

    !d[H2O_dot]/d[E]
    pd(9,1) =  &
        +k(287)*n(idx_H3Oj)  &
        +k(283)*n(idx_H3COj)

    !d[OH_dot]/d[E]
    pd(10,1) =  &
        +k(279)*n(idx_H2Oj)  &
        +k(289)*n(idx_H3Oj)  &
        +k(290)*n(idx_H3Oj)  &
        +k(298)*n(idx_HCO2j)  &
        +k(282)*n(idx_H3COj)

    !d[O2_dot]/d[E]
    pd(11,1) =  &
        +k(312)*n(idx_O2Hj)

    !d[CH2_dot]/d[E]
    pd(12,1) =  &
        +k(271)*n(idx_H2COj)  &
        +k(263)*n(idx_CH3j)  &
        +k(266)*n(idx_CH4j)  &
        +k(282)*n(idx_H3COj)

    !d[H2CO_dot]/d[E]
    pd(13,1) =  &
        +k(1059)*n(idx_H2COj)  &
        +k(285)*n(idx_H3COj)

    !d[HCO_dot]/d[E]
    pd(14,1) =  &
        +k(286)*n(idx_H3COj)  &
        +k(274)*n(idx_H2COj)

    !d[MG_dot]/d[E]
    pd(15,1) =  &
        +k(1061)*n(idx_MGj)

    !d[NO_dot]/d[E]
    pd(17,1) =  &
        +k(276)*n(idx_H2NOj)  &
        +k(299)*n(idx_HNOj)

    !d[CN_dot]/d[E]
    pd(18,1) =  &
        +k(291)*n(idx_HCNj)  &
        +k(292)*n(idx_HCNHj)

    !d[CO_dot]/d[E]
    pd(19,1) =  &
        +k(298)*n(idx_HCO2j)  &
        +k(273)*n(idx_H2COj)  &
        +k(297)*n(idx_HCO2j)  &
        +k(272)*n(idx_H2COj)  &
        +k(295)*n(idx_HCOj)  &
        +k(284)*n(idx_H3COj)  &
        +k(300)*n(idx_HOCj)

    !d[N2_dot]/d[E]
    pd(20,1) =  &
        +k(303)*n(idx_N2Hj)

    !d[NH2_dot]/d[E]
    pd(21,1) =  &
        +k(308)*n(idx_NH3j)

    !d[CH3_dot]/d[E]
    pd(22,1) =  &
        +k(1057)*n(idx_CH3j)  &
        +k(267)*n(idx_CH4j)

    !d[N_dot]/d[E]
    pd(24,1) =  &
        +k(1062)*n(idx_Nj)  &
        +k(305)*n(idx_NHj)  &
        +2.d0*k(302)*n(idx_N2j)  &
        +k(306)*n(idx_NH2j)  &
        +k(304)*n(idx_N2Hj)  &
        +k(268)*n(idx_CNj)  &
        +k(310)*n(idx_NOj)

    !d[NH_dot]/d[E]
    pd(25,1) =  &
        +k(309)*n(idx_NH3j)  &
        +k(304)*n(idx_N2Hj)  &
        +k(307)*n(idx_NH2j)

    !d[HE_dot]/d[E]
    pd(26,1) =  &
        +k(1060)*n(idx_HEj)  &
        +k(301)*n(idx_HEHj)

    !d[HNO_dot]/d[E]
    pd(27,1) =  &
        +k(275)*n(idx_H2NOj)

    !d[CO2_dot]/d[E]
    pd(29,1) =  &
        +k(296)*n(idx_HCO2j)

    !d[E_DUST_dot]/d[E]
    pd(53,1) =  &
        +k(1126)

    !d[HCO+_dot]/d[E]
    pd(54,1) =  &
        -k(295)*n(idx_HCOj)

    !d[H+_dot]/d[E]
    pd(55,1) =  &
        -k(1058)*n(idx_Hj)

    !d[HOC+_dot]/d[E]
    pd(56,1) =  &
        -k(300)*n(idx_HOCj)

    !d[C+_dot]/d[E]
    pd(57,1) =  &
        -k(1056)*n(idx_Cj)

    !d[CH2+_dot]/d[E]
    pd(58,1) =  &
        -k(261)*n(idx_CH2j)  &
        -k(260)*n(idx_CH2j)  &
        -k(262)*n(idx_CH2j)

    !d[CH+_dot]/d[E]
    pd(59,1) =  &
        -k(259)*n(idx_CHj)

    !d[H2CO+_dot]/d[E]
    pd(60,1) =  &
        -k(273)*n(idx_H2COj)  &
        -k(272)*n(idx_H2COj)  &
        -k(274)*n(idx_H2COj)  &
        -k(1059)*n(idx_H2COj)  &
        -k(271)*n(idx_H2COj)

    !d[MG+_dot]/d[E]
    pd(61,1) =  &
        -k(1061)*n(idx_MGj)

    !d[NH3+_dot]/d[E]
    pd(62,1) =  &
        -k(308)*n(idx_NH3j)  &
        -k(309)*n(idx_NH3j)

    !d[NO+_dot]/d[E]
    pd(63,1) =  &
        -k(310)*n(idx_NOj)

    !d[CN+_dot]/d[E]
    pd(64,1) =  &
        -k(268)*n(idx_CNj)

    !d[CO+_dot]/d[E]
    pd(65,1) =  &
        -k(269)*n(idx_COj)

    !d[N2+_dot]/d[E]
    pd(66,1) =  &
        -k(302)*n(idx_N2j)

    !d[O2+_dot]/d[E]
    pd(67,1) =  &
        -k(311)*n(idx_O2j)

    !d[H2O+_dot]/d[E]
    pd(68,1) =  &
        -k(279)*n(idx_H2Oj)  &
        -k(277)*n(idx_H2Oj)  &
        -k(278)*n(idx_H2Oj)

    !d[NH2+_dot]/d[E]
    pd(69,1) =  &
        -k(307)*n(idx_NH2j)  &
        -k(306)*n(idx_NH2j)

    !d[O+_dot]/d[E]
    pd(70,1) =  &
        -k(1063)*n(idx_Oj)

    !d[OH+_dot]/d[E]
    pd(71,1) =  &
        -k(313)*n(idx_OHj)

    !d[CH3+_dot]/d[E]
    pd(72,1) =  &
        -k(263)*n(idx_CH3j)  &
        -k(265)*n(idx_CH3j)  &
        -k(264)*n(idx_CH3j)  &
        -k(1057)*n(idx_CH3j)

    !d[CH4+_dot]/d[E]
    pd(73,1) =  &
        -k(266)*n(idx_CH4j)  &
        -k(267)*n(idx_CH4j)

    !d[N+_dot]/d[E]
    pd(74,1) =  &
        -k(1062)*n(idx_Nj)

    !d[HCN+_dot]/d[E]
    pd(75,1) =  &
        -k(291)*n(idx_HCNj)

    !d[NH+_dot]/d[E]
    pd(76,1) =  &
        -k(305)*n(idx_NHj)

    !d[H2+_dot]/d[E]
    pd(77,1) =  &
        -k(270)*n(idx_H2j)

    !d[HE+_dot]/d[E]
    pd(78,1) =  &
        -k(1060)*n(idx_HEj)

    !d[HNO+_dot]/d[E]
    pd(79,1) =  &
        -k(299)*n(idx_HNOj)

    !d[H2NO+_dot]/d[E]
    pd(80,1) =  &
        -k(276)*n(idx_H2NOj)  &
        -k(275)*n(idx_H2NOj)

    !d[H3+_dot]/d[E]
    pd(81,1) =  &
        -k(281)*n(idx_H3j)  &
        -k(280)*n(idx_H3j)

    !d[H3CO+_dot]/d[E]
    pd(82,1) =  &
        -k(286)*n(idx_H3COj)  &
        -k(283)*n(idx_H3COj)  &
        -k(282)*n(idx_H3COj)  &
        -k(285)*n(idx_H3COj)  &
        -k(284)*n(idx_H3COj)

    !d[H3O+_dot]/d[E]
    pd(83,1) =  &
        -k(287)*n(idx_H3Oj)  &
        -k(288)*n(idx_H3Oj)  &
        -k(289)*n(idx_H3Oj)  &
        -k(290)*n(idx_H3Oj)

    !d[HCNH+_dot]/d[E]
    pd(84,1) =  &
        -k(292)*n(idx_HCNHj)  &
        -k(294)*n(idx_HCNHj)  &
        -k(293)*n(idx_HCNHj)

    !d[HCO2+_dot]/d[E]
    pd(85,1) =  &
        -k(298)*n(idx_HCO2j)  &
        -k(297)*n(idx_HCO2j)  &
        -k(296)*n(idx_HCO2j)

    !d[HEH+_dot]/d[E]
    pd(86,1) =  &
        -k(301)*n(idx_HEHj)

    !d[N2H+_dot]/d[E]
    pd(87,1) =  &
        -k(304)*n(idx_N2Hj)  &
        -k(303)*n(idx_N2Hj)

    !d[O2H+_dot]/d[E]
    pd(88,1) =  &
        -k(312)*n(idx_O2Hj)

    !d[E_dot]/d[CH]
    pd(1,2) =  &
        +k(995)  &
        +k(1)*n(idx_O)

    !d[CH_dot]/d[CH]
    pd(2,2) =  &
        -k(410)*n(idx_Nj)  &
        -k(573)*n(idx_HEj)  &
        -k(809)*n(idx_N2)  &
        -k(406)*n(idx_HCNHj)  &
        -k(813)*n(idx_NO)  &
        -k(811)*n(idx_N)  &
        -k(414)*n(idx_Oj)  &
        -k(806)*n(idx_H2CO)  &
        -k(822)*n(idx_O)  &
        -k(53)*n(idx_NH2j)  &
        -k(823)*n(idx_OH)  &
        -k(52)*n(idx_N2j)  &
        -k(56)*n(idx_OHj)  &
        -k(407)*n(idx_HCNHj)  &
        -k(818)*n(idx_O2)  &
        -k(50)*n(idx_H2Oj)  &
        -k(400)*n(idx_COj)  &
        -k(839)*n(idx_H2)  &
        -k(416)*n(idx_O2Hj)  &
        -k(55)*n(idx_O2j)  &
        -k(413)*n(idx_NH2j)  &
        -k(51)*n(idx_Nj)  &
        -k(444)*n(idx_H2j)  &
        -k(805)*n(idx_CO2)  &
        -k(402)*n(idx_H2Oj)  &
        -k(409)*n(idx_HNOj)  &
        -k(821)*n(idx_O)  &
        -k(54)*n(idx_Oj)  &
        -k(995)  &
        -k(994)  &
        -k(1073)  &
        -k(812)*n(idx_NO)  &
        -k(851)*n(idx_H)  &
        -k(1049)*n(idx_H2)  &
        -k(807)*n(idx_HCO)  &
        -k(124)*n(idx_HEj)  &
        -k(808)*n(idx_HNO)  &
        -k(405)*n(idx_HCNj)  &
        -k(814)*n(idx_NO)  &
        -k(3)*n(idx_H2)  &
        -k(16)*n(idx_Cj)  &
        -k(47)*n(idx_CNj)  &
        -k(408)*n(idx_HCOj)  &
        -k(819)*n(idx_O2H)  &
        -k(506)*n(idx_H3j)  &
        -k(401)*n(idx_H2COj)  &
        -k(48)*n(idx_COj)  &
        -k(225)  &
        -k(49)*n(idx_H2COj)  &
        -k(815)*n(idx_O2)  &
        -k(810)*n(idx_N)  &
        -k(1)*n(idx_O)  &
        -k(820)*n(idx_O2H)  &
        -k(87)*n(idx_H2j)  &
        -k(417)*n(idx_OHj)  &
        -k(10)*n(idx_H)  &
        -k(415)*n(idx_O2j)  &
        -k(412)*n(idx_NHj)  &
        -k(403)*n(idx_H3COj)  &
        -k(404)*n(idx_H3Oj)  &
        -k(816)*n(idx_O2)  &
        -k(411)*n(idx_N2Hj)  &
        -k(817)*n(idx_O2)  &
        -k(72)*n(idx_Hj)

    !d[O_dot]/d[CH]
    pd(3,2) =  &
        -k(1)*n(idx_O)  &
        +k(417)*n(idx_OHj)  &
        +k(54)*n(idx_Oj)  &
        -k(822)*n(idx_O)  &
        -k(821)*n(idx_O)  &
        +k(818)*n(idx_O2)  &
        +k(415)*n(idx_O2j)  &
        +k(812)*n(idx_NO)  &
        +k(816)*n(idx_O2)

    !d[HNC_dot]/d[CH]
    pd(4,2) =  &
        +k(407)*n(idx_HCNHj)

    !d[HCN_dot]/d[CH]
    pd(5,2) =  &
        +k(809)*n(idx_N2)  &
        +k(812)*n(idx_NO)  &
        +k(406)*n(idx_HCNHj)

    !d[H2_dot]/d[CH]
    pd(6,2) =  &
        +k(87)*n(idx_H2j)  &
        -k(3)*n(idx_H2)  &
        +k(3)*n(idx_H2)  &
        +k(506)*n(idx_H3j)  &
        +k(851)*n(idx_H)  &
        -k(839)*n(idx_H2)  &
        -k(1049)*n(idx_H2)

    !d[C_dot]/d[CH]
    pd(7,2) =  &
        +k(10)*n(idx_H)  &
        +k(994)  &
        +k(400)*n(idx_COj)  &
        +k(16)*n(idx_Cj)  &
        +k(3)*n(idx_H2)  &
        +k(851)*n(idx_H)  &
        +k(225)  &
        +k(811)*n(idx_N)  &
        +k(822)*n(idx_O)

    !d[H_dot]/d[CH]
    pd(8,2) =  &
        +k(573)*n(idx_HEj)  &
        +2.d0*k(10)*n(idx_H)  &
        +k(72)*n(idx_Hj)  &
        +k(821)*n(idx_O)  &
        +k(839)*n(idx_H2)  &
        +k(994)  &
        +k(815)*n(idx_O2)  &
        +k(444)*n(idx_H2j)  &
        +k(3)*n(idx_H2)  &
        -k(10)*n(idx_H)  &
        +k(816)*n(idx_O2)  &
        +k(410)*n(idx_Nj)  &
        +k(225)  &
        +k(814)*n(idx_NO)  &
        +k(823)*n(idx_OH)  &
        -k(851)*n(idx_H)  &
        +k(810)*n(idx_N)  &
        +k(414)*n(idx_Oj)

    !d[H2O_dot]/d[CH]
    pd(9,2) =  &
        +k(50)*n(idx_H2Oj)  &
        +k(404)*n(idx_H3Oj)

    !d[OH_dot]/d[CH]
    pd(10,2) =  &
        +k(402)*n(idx_H2Oj)  &
        +k(819)*n(idx_O2H)  &
        -k(823)*n(idx_OH)  &
        +k(817)*n(idx_O2)  &
        +k(822)*n(idx_O)  &
        +k(56)*n(idx_OHj)

    !d[O2_dot]/d[CH]
    pd(11,2) =  &
        +k(416)*n(idx_O2Hj)  &
        +k(820)*n(idx_O2H)  &
        -k(815)*n(idx_O2)  &
        +k(55)*n(idx_O2j)  &
        -k(816)*n(idx_O2)  &
        -k(818)*n(idx_O2)  &
        -k(817)*n(idx_O2)

    !d[CH2_dot]/d[CH]
    pd(12,2) =  &
        +k(808)*n(idx_HNO)  &
        +k(807)*n(idx_HCO)  &
        +k(820)*n(idx_O2H)  &
        +k(839)*n(idx_H2)  &
        +k(806)*n(idx_H2CO)

    !d[H2CO_dot]/d[CH]
    pd(13,2) =  &
        +k(403)*n(idx_H3COj)  &
        +k(49)*n(idx_H2COj)  &
        -k(806)*n(idx_H2CO)

    !d[HCO_dot]/d[CH]
    pd(14,2) =  &
        -k(807)*n(idx_HCO)  &
        +k(823)*n(idx_OH)  &
        +k(819)*n(idx_O2H)  &
        +k(401)*n(idx_H2COj)  &
        +k(806)*n(idx_H2CO)  &
        +k(813)*n(idx_NO)  &
        +k(818)*n(idx_O2)  &
        +k(805)*n(idx_CO2)

    !d[NO_dot]/d[CH]
    pd(17,2) =  &
        -k(812)*n(idx_NO)  &
        +k(409)*n(idx_HNOj)  &
        -k(814)*n(idx_NO)  &
        -k(813)*n(idx_NO)  &
        +k(808)*n(idx_HNO)

    !d[CN_dot]/d[CH]
    pd(18,2) =  &
        +k(405)*n(idx_HCNj)  &
        +k(810)*n(idx_N)  &
        +k(47)*n(idx_CNj)

    !d[CO_dot]/d[CH]
    pd(19,2) =  &
        +k(821)*n(idx_O)  &
        +k(408)*n(idx_HCOj)  &
        +k(48)*n(idx_COj)  &
        +k(807)*n(idx_HCO)  &
        +k(816)*n(idx_O2)  &
        +k(817)*n(idx_O2)  &
        +k(805)*n(idx_CO2)

    !d[N2_dot]/d[CH]
    pd(20,2) =  &
        +k(411)*n(idx_N2Hj)  &
        +k(52)*n(idx_N2j)  &
        -k(809)*n(idx_N2)

    !d[NH2_dot]/d[CH]
    pd(21,2) =  &
        +k(53)*n(idx_NH2j)

    !d[CH3_dot]/d[CH]
    pd(22,2) =  &
        +k(1049)*n(idx_H2)

    !d[N_dot]/d[CH]
    pd(24,2) =  &
        -k(811)*n(idx_N)  &
        -k(810)*n(idx_N)  &
        +k(51)*n(idx_Nj)  &
        +k(412)*n(idx_NHj)  &
        +k(813)*n(idx_NO)  &
        +k(809)*n(idx_N2)

    !d[NH_dot]/d[CH]
    pd(25,2) =  &
        +k(413)*n(idx_NH2j)  &
        +k(811)*n(idx_N)

    !d[HE_dot]/d[CH]
    pd(26,2) =  &
        +k(124)*n(idx_HEj)  &
        +k(573)*n(idx_HEj)

    !d[HNO_dot]/d[CH]
    pd(27,2) =  &
        -k(808)*n(idx_HNO)

    !d[CO2_dot]/d[CH]
    pd(29,2) =  &
        +k(815)*n(idx_O2)  &
        -k(805)*n(idx_CO2)

    !d[O2H_dot]/d[CH]
    pd(33,2) =  &
        -k(820)*n(idx_O2H)  &
        -k(819)*n(idx_O2H)

    !d[OCN_dot]/d[CH]
    pd(34,2) =  &
        +k(814)*n(idx_NO)

    !d[CH4_DUST_dot]/d[CH]
    pd(38,2) =  &
        +k(1073)

    !d[HCO+_dot]/d[CH]
    pd(54,2) =  &
        +k(1)*n(idx_O)  &
        -k(408)*n(idx_HCOj)  &
        +k(415)*n(idx_O2j)  &
        +k(400)*n(idx_COj)

    !d[H+_dot]/d[CH]
    pd(55,2) =  &
        -k(72)*n(idx_Hj)

    !d[C+_dot]/d[CH]
    pd(57,2) =  &
        +k(573)*n(idx_HEj)  &
        -k(16)*n(idx_Cj)

    !d[CH2+_dot]/d[CH]
    pd(58,2) =  &
        +k(407)*n(idx_HCNHj)  &
        +k(409)*n(idx_HNOj)  &
        +k(406)*n(idx_HCNHj)  &
        +k(405)*n(idx_HCNj)  &
        +k(402)*n(idx_H2Oj)  &
        +k(417)*n(idx_OHj)  &
        +k(413)*n(idx_NH2j)  &
        +k(404)*n(idx_H3Oj)  &
        +k(401)*n(idx_H2COj)  &
        +k(506)*n(idx_H3j)  &
        +k(408)*n(idx_HCOj)  &
        +k(411)*n(idx_N2Hj)  &
        +k(416)*n(idx_O2Hj)  &
        +k(444)*n(idx_H2j)  &
        +k(412)*n(idx_NHj)  &
        +k(403)*n(idx_H3COj)

    !d[CH+_dot]/d[CH]
    pd(59,2) =  &
        +k(50)*n(idx_H2Oj)  &
        +k(72)*n(idx_Hj)  &
        +k(55)*n(idx_O2j)  &
        +k(47)*n(idx_CNj)  &
        +k(995)  &
        +k(52)*n(idx_N2j)  &
        +k(51)*n(idx_Nj)  &
        +k(124)*n(idx_HEj)  &
        +k(53)*n(idx_NH2j)  &
        +k(48)*n(idx_COj)  &
        +k(87)*n(idx_H2j)  &
        +k(16)*n(idx_Cj)  &
        +k(49)*n(idx_H2COj)  &
        +k(54)*n(idx_Oj)  &
        +k(56)*n(idx_OHj)

    !d[H2CO+_dot]/d[CH]
    pd(60,2) =  &
        -k(49)*n(idx_H2COj)  &
        -k(401)*n(idx_H2COj)

    !d[CN+_dot]/d[CH]
    pd(64,2) =  &
        -k(47)*n(idx_CNj)  &
        +k(410)*n(idx_Nj)

    !d[CO+_dot]/d[CH]
    pd(65,2) =  &
        +k(414)*n(idx_Oj)  &
        -k(400)*n(idx_COj)  &
        -k(48)*n(idx_COj)

    !d[N2+_dot]/d[CH]
    pd(66,2) =  &
        -k(52)*n(idx_N2j)

    !d[O2+_dot]/d[CH]
    pd(67,2) =  &
        -k(55)*n(idx_O2j)  &
        -k(415)*n(idx_O2j)

    !d[H2O+_dot]/d[CH]
    pd(68,2) =  &
        -k(50)*n(idx_H2Oj)  &
        -k(402)*n(idx_H2Oj)

    !d[NH2+_dot]/d[CH]
    pd(69,2) =  &
        -k(53)*n(idx_NH2j)  &
        -k(413)*n(idx_NH2j)

    !d[O+_dot]/d[CH]
    pd(70,2) =  &
        -k(54)*n(idx_Oj)  &
        -k(414)*n(idx_Oj)

    !d[OH+_dot]/d[CH]
    pd(71,2) =  &
        -k(56)*n(idx_OHj)  &
        -k(417)*n(idx_OHj)

    !d[N+_dot]/d[CH]
    pd(74,2) =  &
        -k(410)*n(idx_Nj)  &
        -k(51)*n(idx_Nj)

    !d[HCN+_dot]/d[CH]
    pd(75,2) =  &
        -k(405)*n(idx_HCNj)

    !d[NH+_dot]/d[CH]
    pd(76,2) =  &
        -k(412)*n(idx_NHj)

    !d[H2+_dot]/d[CH]
    pd(77,2) =  &
        -k(444)*n(idx_H2j)  &
        -k(87)*n(idx_H2j)

    !d[HE+_dot]/d[CH]
    pd(78,2) =  &
        -k(573)*n(idx_HEj)  &
        -k(124)*n(idx_HEj)

    !d[HNO+_dot]/d[CH]
    pd(79,2) =  &
        -k(409)*n(idx_HNOj)

    !d[H3+_dot]/d[CH]
    pd(81,2) =  &
        -k(506)*n(idx_H3j)

    !d[H3CO+_dot]/d[CH]
    pd(82,2) =  &
        -k(403)*n(idx_H3COj)

    !d[H3O+_dot]/d[CH]
    pd(83,2) =  &
        -k(404)*n(idx_H3Oj)

    !d[HCNH+_dot]/d[CH]
    pd(84,2) =  &
        -k(406)*n(idx_HCNHj)  &
        -k(407)*n(idx_HCNHj)

    !d[N2H+_dot]/d[CH]
    pd(87,2) =  &
        -k(411)*n(idx_N2Hj)

    !d[O2H+_dot]/d[CH]
    pd(88,2) =  &
        -k(416)*n(idx_O2Hj)

    !d[E_dot]/d[O]
    pd(1,3) =  &
        +k(214)  &
        +k(1)*n(idx_CH)  &
        +k(256)

    !d[CH_dot]/d[O]
    pd(2,3) =  &
        -k(821)*n(idx_CH)  &
        -k(822)*n(idx_CH)  &
        +k(779)*n(idx_CH2)  &
        -k(1)*n(idx_CH)

    !d[O_dot]/d[O]
    pd(3,3) =  &
        -k(953)*n(idx_NH2)  &
        -k(957)*n(idx_O2H)  &
        -k(927)*n(idx_NH)  &
        -k(797)*n(idx_CH3)  &
        -k(942)*n(idx_H2O)  &
        -k(214)  &
        -k(725)*n(idx_OHj)  &
        -k(719)*n(idx_HCO2j)  &
        -k(940)*n(idx_H2CN)  &
        -k(722)*n(idx_NH2j)  &
        -k(777)*n(idx_CH2)  &
        -k(821)*n(idx_CH)  &
        -k(959)*n(idx_OCN)  &
        -k(946)*n(idx_HCO)  &
        -k(1052)*n(idx_H)  &
        -k(458)*n(idx_H2j)  &
        -k(718)*n(idx_H2Oj)  &
        -k(926)*n(idx_NH)  &
        -k(1111)  &
        -k(386)*n(idx_CH3j)  &
        -k(717)*n(idx_CH4j)  &
        -k(365)*n(idx_CH2j)  &
        -k(720)*n(idx_N2j)  &
        -k(936)*n(idx_CH4)  &
        -k(951)*n(idx_N2)  &
        -k(778)*n(idx_CH2)  &
        -k(952)*n(idx_NH2)  &
        -k(194)*n(idx_CNj)  &
        -k(960)*n(idx_OH)  &
        -k(1044)*n(idx_C)  &
        -k(938)*n(idx_CN)  &
        -k(941)*n(idx_H2CO)  &
        -k(776)*n(idx_CH2)  &
        -k(387)*n(idx_CH3j)  &
        -k(358)*n(idx_CHj)  &
        -k(958)*n(idx_OCN)  &
        -k(662)*n(idx_NHj)  &
        -k(83)*n(idx_Hj)  &
        -k(937)*n(idx_CN)  &
        -k(524)*n(idx_H3j)  &
        -k(949)*n(idx_HNO)  &
        -k(846)*n(idx_H2)  &
        -k(1041)*n(idx_Cj)  &
        -k(943)*n(idx_HCN)  &
        -k(947)*n(idx_HCO)  &
        -k(822)*n(idx_CH)  &
        -k(1)*n(idx_CH)  &
        -k(944)*n(idx_HCN)  &
        -k(195)*n(idx_COj)  &
        -k(955)*n(idx_NO2)  &
        -k(939)*n(idx_CO2)  &
        -k(721)*n(idx_N2Hj)  &
        -k(798)*n(idx_CH3)  &
        -k(196)*n(idx_N2j)  &
        -4.d0*k(1055)*n(idx_O)  &
        -k(779)*n(idx_CH2)  &
        -k(945)*n(idx_HCN)  &
        -k(948)*n(idx_HNO)  &
        -k(256)  &
        -k(950)*n(idx_HNO)  &
        -k(956)*n(idx_NO)  &
        -k(724)*n(idx_O2Hj)  &
        -k(954)*n(idx_NH3)  &
        -k(723)*n(idx_NH3j)  &
        -k(525)*n(idx_H3j)

    !d[HCN_dot]/d[O]
    pd(5,3) =  &
        -k(944)*n(idx_HCN)  &
        -k(945)*n(idx_HCN)  &
        -k(943)*n(idx_HCN)

    !d[H2_dot]/d[O]
    pd(6,3) =  &
        +k(723)*n(idx_NH3j)  &
        +k(718)*n(idx_H2Oj)  &
        -k(846)*n(idx_H2)  &
        +k(387)*n(idx_CH3j)  &
        +k(940)*n(idx_H2CN)  &
        +k(525)*n(idx_H3j)  &
        +k(797)*n(idx_CH3)  &
        +k(776)*n(idx_CH2)

    !d[C_dot]/d[O]
    pd(7,3) =  &
        +k(822)*n(idx_CH)  &
        +k(938)*n(idx_CN)  &
        -k(1044)*n(idx_C)

    !d[H_dot]/d[O]
    pd(8,3) =  &
        +k(458)*n(idx_H2j)  &
        +2.d0*k(777)*n(idx_CH2)  &
        +k(948)*n(idx_HNO)  &
        +k(778)*n(idx_CH2)  &
        -k(1052)*n(idx_H)  &
        +k(365)*n(idx_CH2j)  &
        +k(846)*n(idx_H2)  &
        +k(725)*n(idx_OHj)  &
        +k(386)*n(idx_CH3j)  &
        +k(821)*n(idx_CH)  &
        +k(960)*n(idx_OH)  &
        +k(83)*n(idx_Hj)  &
        +k(798)*n(idx_CH3)  &
        +k(945)*n(idx_HCN)  &
        +k(524)*n(idx_H3j)  &
        +k(358)*n(idx_CHj)  &
        +k(926)*n(idx_NH)  &
        +k(946)*n(idx_HCO)  &
        +k(722)*n(idx_NH2j)  &
        +k(952)*n(idx_NH2)  &
        +k(797)*n(idx_CH3)

    !d[H2O_dot]/d[O]
    pd(9,3) =  &
        -k(942)*n(idx_H2O)

    !d[OH_dot]/d[O]
    pd(10,3) =  &
        +k(822)*n(idx_CH)  &
        +k(941)*n(idx_H2CO)  &
        +k(954)*n(idx_NH3)  &
        +k(957)*n(idx_O2H)  &
        +k(846)*n(idx_H2)  &
        +2.d0*k(942)*n(idx_H2O)  &
        +k(953)*n(idx_NH2)  &
        +k(949)*n(idx_HNO)  &
        +k(936)*n(idx_CH4)  &
        +k(943)*n(idx_HCN)  &
        +k(947)*n(idx_HCO)  &
        +k(717)*n(idx_CH4j)  &
        +k(1052)*n(idx_H)  &
        +k(927)*n(idx_NH)  &
        +k(779)*n(idx_CH2)  &
        -k(960)*n(idx_OH)

    !d[O2_dot]/d[O]
    pd(11,3) =  &
        +k(950)*n(idx_HNO)  &
        +k(956)*n(idx_NO)  &
        +k(959)*n(idx_OCN)  &
        +k(955)*n(idx_NO2)  &
        +k(939)*n(idx_CO2)  &
        +k(960)*n(idx_OH)  &
        +k(719)*n(idx_HCO2j)  &
        +2.d0*k(1055)*n(idx_O)  &
        +k(957)*n(idx_O2H)  &
        +k(724)*n(idx_O2Hj)

    !d[CH2_dot]/d[O]
    pd(12,3) =  &
        -k(779)*n(idx_CH2)  &
        -k(778)*n(idx_CH2)  &
        -k(777)*n(idx_CH2)  &
        -k(776)*n(idx_CH2)

    !d[H2CO_dot]/d[O]
    pd(13,3) =  &
        -k(941)*n(idx_H2CO)  &
        +k(798)*n(idx_CH3)

    !d[HCO_dot]/d[O]
    pd(14,3) =  &
        -k(947)*n(idx_HCO)  &
        +k(941)*n(idx_H2CO)  &
        -k(946)*n(idx_HCO)  &
        +k(778)*n(idx_CH2)

    !d[NH3_dot]/d[O]
    pd(16,3) =  &
        -k(954)*n(idx_NH3)

    !d[NO_dot]/d[O]
    pd(17,3) =  &
        +k(958)*n(idx_OCN)  &
        +k(955)*n(idx_NO2)  &
        +k(951)*n(idx_N2)  &
        +k(949)*n(idx_HNO)  &
        +k(938)*n(idx_CN)  &
        -k(956)*n(idx_NO)  &
        +k(926)*n(idx_NH)

    !d[CN_dot]/d[O]
    pd(18,3) =  &
        -k(938)*n(idx_CN)  &
        +k(943)*n(idx_HCN)  &
        +k(194)*n(idx_CNj)  &
        +k(959)*n(idx_OCN)  &
        -k(937)*n(idx_CN)

    !d[CO_dot]/d[O]
    pd(19,3) =  &
        +k(777)*n(idx_CH2)  &
        +k(944)*n(idx_HCN)  &
        +k(958)*n(idx_OCN)  &
        +k(939)*n(idx_CO2)  &
        +k(821)*n(idx_CH)  &
        +k(947)*n(idx_HCO)  &
        +k(797)*n(idx_CH3)  &
        +k(195)*n(idx_COj)  &
        +k(937)*n(idx_CN)  &
        +k(776)*n(idx_CH2)  &
        +k(1044)*n(idx_C)

    !d[N2_dot]/d[O]
    pd(20,3) =  &
        -k(951)*n(idx_N2)  &
        +k(196)*n(idx_N2j)  &
        +k(721)*n(idx_N2Hj)

    !d[NH2_dot]/d[O]
    pd(21,3) =  &
        -k(953)*n(idx_NH2)  &
        +k(954)*n(idx_NH3)  &
        -k(952)*n(idx_NH2)

    !d[CH3_dot]/d[O]
    pd(22,3) =  &
        -k(798)*n(idx_CH3)  &
        +k(936)*n(idx_CH4)  &
        -k(797)*n(idx_CH3)

    !d[CH4_dot]/d[O]
    pd(23,3) =  &
        -k(936)*n(idx_CH4)

    !d[N_dot]/d[O]
    pd(24,3) =  &
        +k(956)*n(idx_NO)  &
        +k(720)*n(idx_N2j)  &
        +k(951)*n(idx_N2)  &
        +k(937)*n(idx_CN)  &
        +k(662)*n(idx_NHj)  &
        +k(927)*n(idx_NH)

    !d[NH_dot]/d[O]
    pd(25,3) =  &
        -k(926)*n(idx_NH)  &
        +k(950)*n(idx_HNO)  &
        -k(927)*n(idx_NH)  &
        +k(944)*n(idx_HCN)  &
        +k(953)*n(idx_NH2)

    !d[HNO_dot]/d[O]
    pd(27,3) =  &
        -k(948)*n(idx_HNO)  &
        -k(950)*n(idx_HNO)  &
        +k(952)*n(idx_NH2)  &
        -k(949)*n(idx_HNO)

    !d[CO2_dot]/d[O]
    pd(29,3) =  &
        -k(939)*n(idx_CO2)  &
        +k(946)*n(idx_HCO)

    !d[H2CN_dot]/d[O]
    pd(30,3) =  &
        -k(940)*n(idx_H2CN)

    !d[NO2_dot]/d[O]
    pd(32,3) =  &
        +k(948)*n(idx_HNO)  &
        -k(955)*n(idx_NO2)

    !d[O2H_dot]/d[O]
    pd(33,3) =  &
        -k(957)*n(idx_O2H)

    !d[OCN_dot]/d[O]
    pd(34,3) =  &
        -k(958)*n(idx_OCN)  &
        -k(959)*n(idx_OCN)  &
        +k(940)*n(idx_H2CN)  &
        +k(945)*n(idx_HCN)

    !d[H2O_DUST_dot]/d[O]
    pd(40,3) =  &
        +k(1111)

    !d[HCO+_dot]/d[O]
    pd(54,3) =  &
        +k(365)*n(idx_CH2j)  &
        +k(387)*n(idx_CH3j)  &
        +k(1)*n(idx_CH)  &
        +k(719)*n(idx_HCO2j)

    !d[H+_dot]/d[O]
    pd(55,3) =  &
        -k(83)*n(idx_Hj)

    !d[C+_dot]/d[O]
    pd(57,3) =  &
        -k(1041)*n(idx_Cj)

    !d[CH2+_dot]/d[O]
    pd(58,3) =  &
        -k(365)*n(idx_CH2j)

    !d[CH+_dot]/d[O]
    pd(59,3) =  &
        -k(358)*n(idx_CHj)

    !d[H2CO+_dot]/d[O]
    pd(60,3) =  &
        +k(386)*n(idx_CH3j)

    !d[NH3+_dot]/d[O]
    pd(62,3) =  &
        -k(723)*n(idx_NH3j)

    !d[NO+_dot]/d[O]
    pd(63,3) =  &
        +k(720)*n(idx_N2j)

    !d[CN+_dot]/d[O]
    pd(64,3) =  &
        -k(194)*n(idx_CNj)

    !d[CO+_dot]/d[O]
    pd(65,3) =  &
        +k(358)*n(idx_CHj)  &
        +k(1041)*n(idx_Cj)  &
        -k(195)*n(idx_COj)

    !d[N2+_dot]/d[O]
    pd(66,3) =  &
        -k(196)*n(idx_N2j)  &
        -k(720)*n(idx_N2j)

    !d[O2+_dot]/d[O]
    pd(67,3) =  &
        +k(725)*n(idx_OHj)  &
        +k(718)*n(idx_H2Oj)

    !d[H2O+_dot]/d[O]
    pd(68,3) =  &
        +k(524)*n(idx_H3j)  &
        -k(718)*n(idx_H2Oj)

    !d[NH2+_dot]/d[O]
    pd(69,3) =  &
        -k(722)*n(idx_NH2j)

    !d[O+_dot]/d[O]
    pd(70,3) =  &
        +k(256)  &
        +k(214)  &
        +k(194)*n(idx_CNj)  &
        +k(196)*n(idx_N2j)  &
        +k(83)*n(idx_Hj)  &
        +k(195)*n(idx_COj)

    !d[OH+_dot]/d[O]
    pd(71,3) =  &
        +k(458)*n(idx_H2j)  &
        +k(721)*n(idx_N2Hj)  &
        +k(525)*n(idx_H3j)  &
        +k(662)*n(idx_NHj)  &
        +k(724)*n(idx_O2Hj)  &
        -k(725)*n(idx_OHj)

    !d[CH3+_dot]/d[O]
    pd(72,3) =  &
        -k(386)*n(idx_CH3j)  &
        -k(387)*n(idx_CH3j)  &
        +k(717)*n(idx_CH4j)

    !d[CH4+_dot]/d[O]
    pd(73,3) =  &
        -k(717)*n(idx_CH4j)

    !d[NH+_dot]/d[O]
    pd(76,3) =  &
        -k(662)*n(idx_NHj)

    !d[H2+_dot]/d[O]
    pd(77,3) =  &
        -k(458)*n(idx_H2j)

    !d[HNO+_dot]/d[O]
    pd(79,3) =  &
        +k(723)*n(idx_NH3j)  &
        +k(722)*n(idx_NH2j)

    !d[H3+_dot]/d[O]
    pd(81,3) =  &
        -k(524)*n(idx_H3j)  &
        -k(525)*n(idx_H3j)

    !d[HCO2+_dot]/d[O]
    pd(85,3) =  &
        -k(719)*n(idx_HCO2j)

    !d[N2H+_dot]/d[O]
    pd(87,3) =  &
        -k(721)*n(idx_N2Hj)

    !d[O2H+_dot]/d[O]
    pd(88,3) =  &
        -k(724)*n(idx_O2Hj)

    !d[O_dot]/d[HNC]
    pd(3,4) =  &
        +k(734)*n(idx_OHj)

    !d[HNC_dot]/d[HNC]
    pd(4,4) =  &
        -k(529)*n(idx_H3Oj)  &
        -k(489)*n(idx_H2Oj)  &
        -k(655)*n(idx_NHj)  &
        -k(351)*n(idx_CHj)  &
        -k(595)*n(idx_HEj)  &
        -k(594)*n(idx_HEj)  &
        -k(734)*n(idx_OHj)  &
        -k(236)  &
        -k(562)*n(idx_O2Hj)  &
        -k(593)*n(idx_HEj)  &
        -k(561)*n(idx_N2Hj)  &
        -k(670)*n(idx_NH2j)  &
        -k(560)*n(idx_HNOj)  &
        -k(557)*n(idx_H2COj)  &
        -k(1015)  &
        -k(515)*n(idx_H3j)  &
        -k(2)*n(idx_Hj)  &
        -k(559)*n(idx_HCOj)  &
        -k(541)*n(idx_HCNj)  &
        -k(860)*n(idx_H)  &
        -k(1125)  &
        -k(558)*n(idx_H3COj)

    !d[HCN_dot]/d[HNC]
    pd(5,4) =  &
        +k(2)*n(idx_Hj)  &
        +k(860)*n(idx_H)

    !d[H2_dot]/d[HNC]
    pd(6,4) =  &
        +k(515)*n(idx_H3j)

    !d[C_dot]/d[HNC]
    pd(7,4) =  &
        +k(351)*n(idx_CHj)  &
        +k(595)*n(idx_HEj)

    !d[H_dot]/d[HNC]
    pd(8,4) =  &
        +k(1015)  &
        +k(594)*n(idx_HEj)  &
        +k(593)*n(idx_HEj)  &
        -k(860)*n(idx_H)  &
        +k(236)  &
        +k(860)*n(idx_H)

    !d[H2O_dot]/d[HNC]
    pd(9,4) =  &
        +k(529)*n(idx_H3Oj)

    !d[OH_dot]/d[HNC]
    pd(10,4) =  &
        +k(489)*n(idx_H2Oj)

    !d[O2_dot]/d[HNC]
    pd(11,4) =  &
        +k(562)*n(idx_O2Hj)

    !d[H2CO_dot]/d[HNC]
    pd(13,4) =  &
        +k(558)*n(idx_H3COj)

    !d[HCO_dot]/d[HNC]
    pd(14,4) =  &
        +k(557)*n(idx_H2COj)

    !d[NO_dot]/d[HNC]
    pd(17,4) =  &
        +k(560)*n(idx_HNOj)

    !d[CN_dot]/d[HNC]
    pd(18,4) =  &
        +k(541)*n(idx_HCNj)  &
        +k(1015)  &
        +k(236)

    !d[CO_dot]/d[HNC]
    pd(19,4) =  &
        +k(559)*n(idx_HCOj)

    !d[N2_dot]/d[HNC]
    pd(20,4) =  &
        +k(561)*n(idx_N2Hj)

    !d[N_dot]/d[HNC]
    pd(24,4) =  &
        +k(655)*n(idx_NHj)  &
        +k(594)*n(idx_HEj)

    !d[NH_dot]/d[HNC]
    pd(25,4) =  &
        +k(670)*n(idx_NH2j)

    !d[HE_dot]/d[HNC]
    pd(26,4) =  &
        +k(593)*n(idx_HEj)  &
        +k(595)*n(idx_HEj)  &
        +k(594)*n(idx_HEj)

    !d[HNC_DUST_dot]/d[HNC]
    pd(52,4) =  &
        +k(1125)

    !d[HCO+_dot]/d[HNC]
    pd(54,4) =  &
        -k(559)*n(idx_HCOj)

    !d[H+_dot]/d[HNC]
    pd(55,4) =  &
        +k(2)*n(idx_Hj)  &
        -k(2)*n(idx_Hj)

    !d[C+_dot]/d[HNC]
    pd(57,4) =  &
        +k(594)*n(idx_HEj)

    !d[CH+_dot]/d[HNC]
    pd(59,4) =  &
        -k(351)*n(idx_CHj)

    !d[H2CO+_dot]/d[HNC]
    pd(60,4) =  &
        -k(557)*n(idx_H2COj)

    !d[CN+_dot]/d[HNC]
    pd(64,4) =  &
        +k(593)*n(idx_HEj)

    !d[H2O+_dot]/d[HNC]
    pd(68,4) =  &
        -k(489)*n(idx_H2Oj)

    !d[NH2+_dot]/d[HNC]
    pd(69,4) =  &
        -k(670)*n(idx_NH2j)

    !d[OH+_dot]/d[HNC]
    pd(71,4) =  &
        -k(734)*n(idx_OHj)

    !d[HCN+_dot]/d[HNC]
    pd(75,4) =  &
        -k(541)*n(idx_HCNj)

    !d[NH+_dot]/d[HNC]
    pd(76,4) =  &
        -k(655)*n(idx_NHj)  &
        +k(595)*n(idx_HEj)

    !d[HE+_dot]/d[HNC]
    pd(78,4) =  &
        -k(593)*n(idx_HEj)  &
        -k(595)*n(idx_HEj)  &
        -k(594)*n(idx_HEj)

    !d[HNO+_dot]/d[HNC]
    pd(79,4) =  &
        -k(560)*n(idx_HNOj)

    !d[H3+_dot]/d[HNC]
    pd(81,4) =  &
        -k(515)*n(idx_H3j)

    !d[H3CO+_dot]/d[HNC]
    pd(82,4) =  &
        -k(558)*n(idx_H3COj)

    !d[H3O+_dot]/d[HNC]
    pd(83,4) =  &
        -k(529)*n(idx_H3Oj)

    !d[HCNH+_dot]/d[HNC]
    pd(84,4) =  &
        +k(541)*n(idx_HCNj)  &
        +k(655)*n(idx_NHj)  &
        +k(559)*n(idx_HCOj)  &
        +k(562)*n(idx_O2Hj)  &
        +k(561)*n(idx_N2Hj)  &
        +k(351)*n(idx_CHj)  &
        +k(734)*n(idx_OHj)  &
        +k(557)*n(idx_H2COj)  &
        +k(515)*n(idx_H3j)  &
        +k(489)*n(idx_H2Oj)  &
        +k(558)*n(idx_H3COj)  &
        +k(560)*n(idx_HNOj)  &
        +k(670)*n(idx_NH2j)  &
        +k(529)*n(idx_H3Oj)

    !d[N2H+_dot]/d[HNC]
    pd(87,4) =  &
        -k(561)*n(idx_N2Hj)

    !d[O2H+_dot]/d[HNC]
    pd(88,4) =  &
        -k(562)*n(idx_O2Hj)

    !d[CH_dot]/d[HCN]
    pd(2,5) =  &
        +k(710)*n(idx_Oj)  &
        +k(587)*n(idx_HEj)

    !d[O_dot]/d[HCN]
    pd(3,5) =  &
        -k(943)*n(idx_O)  &
        -k(945)*n(idx_O)  &
        -k(944)*n(idx_O)  &
        +k(731)*n(idx_OHj)

    !d[HCN_dot]/d[HCN]
    pd(5,5) =  &
        -k(943)*n(idx_O)  &
        -k(653)*n(idx_NHj)  &
        -k(709)*n(idx_Oj)  &
        -k(544)*n(idx_HCOj)  &
        -k(59)*n(idx_CNj)  &
        -k(118)*n(idx_COj)  &
        -k(944)*n(idx_O)  &
        -k(75)*n(idx_Hj)  &
        -k(1011)  &
        -k(349)*n(idx_CHj)  &
        -k(589)*n(idx_HEj)  &
        -k(141)*n(idx_Nj)  &
        -k(857)*n(idx_H)  &
        -k(92)*n(idx_H2j)  &
        -k(546)*n(idx_N2Hj)  &
        -k(545)*n(idx_HNOj)  &
        -k(710)*n(idx_Oj)  &
        -k(119)*n(idx_N2j)  &
        -k(588)*n(idx_HEj)  &
        -k(587)*n(idx_HEj)  &
        -k(513)*n(idx_H3j)  &
        -k(668)*n(idx_NH2j)  &
        -k(486)*n(idx_H2Oj)  &
        -k(945)*n(idx_O)  &
        -k(731)*n(idx_OHj)  &
        -k(1086)  &
        -k(965)*n(idx_OH)  &
        -k(547)*n(idx_O2Hj)  &
        -k(966)*n(idx_OH)  &
        -k(233)  &
        -k(542)*n(idx_H2COj)  &
        -k(543)*n(idx_H3COj)  &
        -k(586)*n(idx_HEj)  &
        -k(528)*n(idx_H3Oj)  &
        -k(538)*n(idx_HCNj)

    !d[H2_dot]/d[HCN]
    pd(6,5) =  &
        +k(857)*n(idx_H)  &
        +k(92)*n(idx_H2j)  &
        +k(513)*n(idx_H3j)

    !d[C_dot]/d[HCN]
    pd(7,5) =  &
        +k(349)*n(idx_CHj)

    !d[H_dot]/d[HCN]
    pd(8,5) =  &
        +k(233)  &
        +k(945)*n(idx_O)  &
        -k(857)*n(idx_H)  &
        +k(1011)  &
        +k(588)*n(idx_HEj)  &
        +k(75)*n(idx_Hj)  &
        +k(586)*n(idx_HEj)

    !d[H2O_dot]/d[HCN]
    pd(9,5) =  &
        +k(965)*n(idx_OH)  &
        +k(528)*n(idx_H3Oj)

    !d[OH_dot]/d[HCN]
    pd(10,5) =  &
        -k(965)*n(idx_OH)  &
        +k(943)*n(idx_O)  &
        +k(486)*n(idx_H2Oj)  &
        -k(966)*n(idx_OH)

    !d[O2_dot]/d[HCN]
    pd(11,5) =  &
        +k(547)*n(idx_O2Hj)

    !d[H2CO_dot]/d[HCN]
    pd(13,5) =  &
        +k(543)*n(idx_H3COj)

    !d[HCO_dot]/d[HCN]
    pd(14,5) =  &
        +k(542)*n(idx_H2COj)

    !d[NO_dot]/d[HCN]
    pd(17,5) =  &
        +k(545)*n(idx_HNOj)

    !d[CN_dot]/d[HCN]
    pd(18,5) =  &
        +k(233)  &
        +k(538)*n(idx_HCNj)  &
        +k(943)*n(idx_O)  &
        +k(1011)  &
        +k(965)*n(idx_OH)  &
        +k(59)*n(idx_CNj)  &
        +k(857)*n(idx_H)

    !d[CO_dot]/d[HCN]
    pd(19,5) =  &
        +k(944)*n(idx_O)  &
        +k(118)*n(idx_COj)  &
        +k(544)*n(idx_HCOj)  &
        +k(966)*n(idx_OH)

    !d[N2_dot]/d[HCN]
    pd(20,5) =  &
        +k(546)*n(idx_N2Hj)  &
        +k(119)*n(idx_N2j)

    !d[NH2_dot]/d[HCN]
    pd(21,5) =  &
        +k(966)*n(idx_OH)

    !d[N_dot]/d[HCN]
    pd(24,5) =  &
        +k(588)*n(idx_HEj)  &
        +k(589)*n(idx_HEj)  &
        +k(653)*n(idx_NHj)  &
        +k(709)*n(idx_Oj)  &
        +k(141)*n(idx_Nj)

    !d[NH_dot]/d[HCN]
    pd(25,5) =  &
        +k(944)*n(idx_O)  &
        +k(668)*n(idx_NH2j)

    !d[HE_dot]/d[HCN]
    pd(26,5) =  &
        +k(588)*n(idx_HEj)  &
        +k(589)*n(idx_HEj)  &
        +k(586)*n(idx_HEj)  &
        +k(587)*n(idx_HEj)

    !d[OCN_dot]/d[HCN]
    pd(34,5) =  &
        +k(945)*n(idx_O)

    !d[HCN_DUST_dot]/d[HCN]
    pd(44,5) =  &
        +k(1086)

    !d[HCO+_dot]/d[HCN]
    pd(54,5) =  &
        -k(544)*n(idx_HCOj)  &
        +k(709)*n(idx_Oj)

    !d[H+_dot]/d[HCN]
    pd(55,5) =  &
        -k(75)*n(idx_Hj)

    !d[C+_dot]/d[HCN]
    pd(57,5) =  &
        +k(588)*n(idx_HEj)

    !d[CH+_dot]/d[HCN]
    pd(59,5) =  &
        -k(349)*n(idx_CHj)  &
        +k(589)*n(idx_HEj)

    !d[H2CO+_dot]/d[HCN]
    pd(60,5) =  &
        -k(542)*n(idx_H2COj)

    !d[NO+_dot]/d[HCN]
    pd(63,5) =  &
        +k(710)*n(idx_Oj)

    !d[CN+_dot]/d[HCN]
    pd(64,5) =  &
        -k(59)*n(idx_CNj)  &
        +k(586)*n(idx_HEj)

    !d[CO+_dot]/d[HCN]
    pd(65,5) =  &
        -k(118)*n(idx_COj)

    !d[N2+_dot]/d[HCN]
    pd(66,5) =  &
        -k(119)*n(idx_N2j)

    !d[H2O+_dot]/d[HCN]
    pd(68,5) =  &
        -k(486)*n(idx_H2Oj)

    !d[NH2+_dot]/d[HCN]
    pd(69,5) =  &
        -k(668)*n(idx_NH2j)

    !d[O+_dot]/d[HCN]
    pd(70,5) =  &
        -k(710)*n(idx_Oj)  &
        -k(709)*n(idx_Oj)

    !d[OH+_dot]/d[HCN]
    pd(71,5) =  &
        -k(731)*n(idx_OHj)

    !d[N+_dot]/d[HCN]
    pd(74,5) =  &
        -k(141)*n(idx_Nj)  &
        +k(587)*n(idx_HEj)

    !d[HCN+_dot]/d[HCN]
    pd(75,5) =  &
        +k(92)*n(idx_H2j)  &
        +k(141)*n(idx_Nj)  &
        +k(118)*n(idx_COj)  &
        +k(59)*n(idx_CNj)  &
        +k(119)*n(idx_N2j)  &
        +k(75)*n(idx_Hj)  &
        -k(538)*n(idx_HCNj)

    !d[NH+_dot]/d[HCN]
    pd(76,5) =  &
        -k(653)*n(idx_NHj)

    !d[H2+_dot]/d[HCN]
    pd(77,5) =  &
        -k(92)*n(idx_H2j)

    !d[HE+_dot]/d[HCN]
    pd(78,5) =  &
        -k(586)*n(idx_HEj)  &
        -k(587)*n(idx_HEj)  &
        -k(588)*n(idx_HEj)  &
        -k(589)*n(idx_HEj)

    !d[HNO+_dot]/d[HCN]
    pd(79,5) =  &
        -k(545)*n(idx_HNOj)

    !d[H3+_dot]/d[HCN]
    pd(81,5) =  &
        -k(513)*n(idx_H3j)

    !d[H3CO+_dot]/d[HCN]
    pd(82,5) =  &
        -k(543)*n(idx_H3COj)

    !d[H3O+_dot]/d[HCN]
    pd(83,5) =  &
        -k(528)*n(idx_H3Oj)

    !d[HCNH+_dot]/d[HCN]
    pd(84,5) =  &
        +k(543)*n(idx_H3COj)  &
        +k(731)*n(idx_OHj)  &
        +k(513)*n(idx_H3j)  &
        +k(538)*n(idx_HCNj)  &
        +k(546)*n(idx_N2Hj)  &
        +k(542)*n(idx_H2COj)  &
        +k(349)*n(idx_CHj)  &
        +k(486)*n(idx_H2Oj)  &
        +k(545)*n(idx_HNOj)  &
        +k(653)*n(idx_NHj)  &
        +k(547)*n(idx_O2Hj)  &
        +k(528)*n(idx_H3Oj)  &
        +k(544)*n(idx_HCOj)  &
        +k(668)*n(idx_NH2j)

    !d[N2H+_dot]/d[HCN]
    pd(87,5) =  &
        -k(546)*n(idx_N2Hj)

    !d[O2H+_dot]/d[HCN]
    pd(88,5) =  &
        -k(547)*n(idx_O2Hj)

    !d[E_dot]/d[H2]
    pd(1,6) =  &
        +k(9)*n(idx_E)  &
        -k(9)*n(idx_E)  &
        +k(209)  &
        +k(208)

    !d[CH_dot]/d[H2]
    pd(2,6) =  &
        -k(3)*n(idx_CH)  &
        -k(839)*n(idx_CH)  &
        +k(836)*n(idx_C)  &
        -k(1049)*n(idx_CH)

    !d[O_dot]/d[H2]
    pd(3,6) =  &
        +2.d0*k(7)*n(idx_O2)  &
        -k(846)*n(idx_O)  &
        +k(8)*n(idx_OH)

    !d[HCN_dot]/d[H2]
    pd(5,6) =  &
        +k(840)*n(idx_CN)

    !d[H2_dot]/d[H2]
    pd(6,6) =  &
        -k(1049)*n(idx_CH)  &
        -k(9)*n(idx_E)  &
        -k(208)  &
        -k(8)*n(idx_OH)  &
        -k(476)*n(idx_O2Hj)  &
        -k(473)*n(idx_NHj)  &
        -k(469)*n(idx_HEHj)  &
        -k(463)*n(idx_CNj)  &
        +k(7)*n(idx_O2)  &
        -k(471)*n(idx_N2j)  &
        +k(5)*n(idx_H2O)  &
        -k(470)*n(idx_Nj)  &
        -k(475)*n(idx_Oj)  &
        -k(465)*n(idx_COj)  &
        -k(7)*n(idx_O2)  &
        -k(840)*n(idx_CN)  &
        -k(1047)*n(idx_Cj)  &
        -k(461)*n(idx_CHj)  &
        -k(1048)*n(idx_C)  &
        -k(847)*n(idx_OH)  &
        -k(846)*n(idx_O)  &
        -k(466)*n(idx_H2Oj)  &
        -k(210)  &
        -k(837)*n(idx_CH2)  &
        -k(836)*n(idx_C)  &
        -k(841)*n(idx_N)  &
        -k(209)  &
        -k(843)*n(idx_NH)  &
        -k(839)*n(idx_CH)  &
        -k(5)*n(idx_H2O)  &
        -k(838)*n(idx_CH3)  &
        -k(474)*n(idx_NH2j)  &
        +2.d0*k(4)*n(idx_H2)  &
        -k(460)*n(idx_Cj)  &
        +k(3)*n(idx_CH)  &
        -k(462)*n(idx_CH2j)  &
        -k(3)*n(idx_CH)  &
        -k(6)*n(idx_HOCj)  &
        -k(468)*n(idx_HEj)  &
        -4.d0*k(4)*n(idx_H2)  &
        -k(845)*n(idx_O2)  &
        -k(467)*n(idx_HCNj)  &
        -k(448)*n(idx_H2j)  &
        -k(477)*n(idx_OHj)  &
        -k(464)*n(idx_COj)  &
        -k(844)*n(idx_O2)  &
        -k(842)*n(idx_NH2)  &
        +k(8)*n(idx_OH)  &
        -k(100)*n(idx_HEj)  &
        -k(472)*n(idx_NHj)  &
        +k(6)*n(idx_HOCj)  &
        -k(11)*n(idx_H)

    !d[C_dot]/d[H2]
    pd(7,6) =  &
        -k(836)*n(idx_C)  &
        -k(1048)*n(idx_C)  &
        +k(3)*n(idx_CH)

    !d[H_dot]/d[H2]
    pd(8,6) =  &
        +k(474)*n(idx_NH2j)  &
        +k(5)*n(idx_H2O)  &
        +k(475)*n(idx_Oj)  &
        +4.d0*k(4)*n(idx_H2)  &
        +k(208)  &
        +k(839)*n(idx_CH)  &
        +k(838)*n(idx_CH3)  &
        +k(471)*n(idx_N2j)  &
        +k(477)*n(idx_OHj)  &
        +k(470)*n(idx_Nj)  &
        +k(842)*n(idx_NH2)  &
        +k(461)*n(idx_CHj)  &
        +k(465)*n(idx_COj)  &
        +k(462)*n(idx_CH2j)  &
        +k(837)*n(idx_CH2)  &
        +2.d0*k(9)*n(idx_E)  &
        +k(843)*n(idx_NH)  &
        -k(11)*n(idx_H)  &
        +k(467)*n(idx_HCNj)  &
        +k(3)*n(idx_CH)  &
        +k(460)*n(idx_Cj)  &
        +2.d0*k(210)  &
        +k(847)*n(idx_OH)  &
        +k(468)*n(idx_HEj)  &
        +k(473)*n(idx_NHj)  &
        +k(844)*n(idx_O2)  &
        +k(836)*n(idx_C)  &
        +k(841)*n(idx_N)  &
        +k(840)*n(idx_CN)  &
        +k(846)*n(idx_O)  &
        +3.d0*k(11)*n(idx_H)  &
        +k(448)*n(idx_H2j)  &
        +k(8)*n(idx_OH)  &
        +k(463)*n(idx_CNj)  &
        +k(464)*n(idx_COj)  &
        +k(466)*n(idx_H2Oj)

    !d[H2O_dot]/d[H2]
    pd(9,6) =  &
        +k(847)*n(idx_OH)  &
        -k(5)*n(idx_H2O)

    !d[OH_dot]/d[H2]
    pd(10,6) =  &
        -k(8)*n(idx_OH)  &
        -k(847)*n(idx_OH)  &
        +2.d0*k(845)*n(idx_O2)  &
        +k(5)*n(idx_H2O)  &
        +k(846)*n(idx_O)

    !d[O2_dot]/d[H2]
    pd(11,6) =  &
        -k(844)*n(idx_O2)  &
        -k(7)*n(idx_O2)  &
        -k(845)*n(idx_O2)  &
        +k(476)*n(idx_O2Hj)

    !d[CH2_dot]/d[H2]
    pd(12,6) =  &
        +k(839)*n(idx_CH)  &
        -k(837)*n(idx_CH2)  &
        +k(1048)*n(idx_C)

    !d[NH3_dot]/d[H2]
    pd(16,6) =  &
        +k(842)*n(idx_NH2)

    !d[CN_dot]/d[H2]
    pd(18,6) =  &
        -k(840)*n(idx_CN)

    !d[NH2_dot]/d[H2]
    pd(21,6) =  &
        -k(842)*n(idx_NH2)  &
        +k(843)*n(idx_NH)

    !d[CH3_dot]/d[H2]
    pd(22,6) =  &
        +k(1049)*n(idx_CH)  &
        -k(838)*n(idx_CH3)  &
        +k(837)*n(idx_CH2)

    !d[CH4_dot]/d[H2]
    pd(23,6) =  &
        +k(838)*n(idx_CH3)

    !d[N_dot]/d[H2]
    pd(24,6) =  &
        +k(472)*n(idx_NHj)  &
        -k(841)*n(idx_N)

    !d[NH_dot]/d[H2]
    pd(25,6) =  &
        -k(843)*n(idx_NH)  &
        +k(841)*n(idx_N)

    !d[HE_dot]/d[H2]
    pd(26,6) =  &
        +k(469)*n(idx_HEHj)  &
        +k(468)*n(idx_HEj)  &
        +k(100)*n(idx_HEj)

    !d[O2H_dot]/d[H2]
    pd(33,6) =  &
        +k(844)*n(idx_O2)

    !d[HCO+_dot]/d[H2]
    pd(54,6) =  &
        +k(464)*n(idx_COj)  &
        +k(6)*n(idx_HOCj)

    !d[H+_dot]/d[H2]
    pd(55,6) =  &
        +k(468)*n(idx_HEj)  &
        +k(208)

    !d[HOC+_dot]/d[H2]
    pd(56,6) =  &
        +k(465)*n(idx_COj)  &
        -k(6)*n(idx_HOCj)

    !d[C+_dot]/d[H2]
    pd(57,6) =  &
        -k(460)*n(idx_Cj)  &
        -k(1047)*n(idx_Cj)

    !d[CH2+_dot]/d[H2]
    pd(58,6) =  &
        -k(462)*n(idx_CH2j)  &
        +k(1047)*n(idx_Cj)  &
        +k(461)*n(idx_CHj)

    !d[CH+_dot]/d[H2]
    pd(59,6) =  &
        -k(461)*n(idx_CHj)  &
        +k(460)*n(idx_Cj)

    !d[NH3+_dot]/d[H2]
    pd(62,6) =  &
        +k(474)*n(idx_NH2j)

    !d[CN+_dot]/d[H2]
    pd(64,6) =  &
        -k(463)*n(idx_CNj)

    !d[CO+_dot]/d[H2]
    pd(65,6) =  &
        -k(464)*n(idx_COj)  &
        -k(465)*n(idx_COj)

    !d[N2+_dot]/d[H2]
    pd(66,6) =  &
        -k(471)*n(idx_N2j)

    !d[H2O+_dot]/d[H2]
    pd(68,6) =  &
        -k(466)*n(idx_H2Oj)  &
        +k(477)*n(idx_OHj)

    !d[NH2+_dot]/d[H2]
    pd(69,6) =  &
        +k(473)*n(idx_NHj)  &
        -k(474)*n(idx_NH2j)

    !d[O+_dot]/d[H2]
    pd(70,6) =  &
        -k(475)*n(idx_Oj)

    !d[OH+_dot]/d[H2]
    pd(71,6) =  &
        -k(477)*n(idx_OHj)  &
        +k(475)*n(idx_Oj)

    !d[CH3+_dot]/d[H2]
    pd(72,6) =  &
        +k(462)*n(idx_CH2j)

    !d[N+_dot]/d[H2]
    pd(74,6) =  &
        -k(470)*n(idx_Nj)

    !d[HCN+_dot]/d[H2]
    pd(75,6) =  &
        -k(467)*n(idx_HCNj)  &
        +k(463)*n(idx_CNj)

    !d[NH+_dot]/d[H2]
    pd(76,6) =  &
        -k(472)*n(idx_NHj)  &
        +k(470)*n(idx_Nj)  &
        -k(473)*n(idx_NHj)

    !d[H2+_dot]/d[H2]
    pd(77,6) =  &
        -k(448)*n(idx_H2j)  &
        +k(100)*n(idx_HEj)  &
        +k(209)

    !d[HE+_dot]/d[H2]
    pd(78,6) =  &
        -k(100)*n(idx_HEj)  &
        -k(468)*n(idx_HEj)

    !d[H3+_dot]/d[H2]
    pd(81,6) =  &
        +k(469)*n(idx_HEHj)  &
        +k(476)*n(idx_O2Hj)  &
        +k(448)*n(idx_H2j)  &
        +k(472)*n(idx_NHj)

    !d[H3O+_dot]/d[H2]
    pd(83,6) =  &
        +k(466)*n(idx_H2Oj)

    !d[HCNH+_dot]/d[H2]
    pd(84,6) =  &
        +k(467)*n(idx_HCNj)

    !d[HEH+_dot]/d[H2]
    pd(86,6) =  &
        -k(469)*n(idx_HEHj)

    !d[N2H+_dot]/d[H2]
    pd(87,6) =  &
        +k(471)*n(idx_N2j)

    !d[O2H+_dot]/d[H2]
    pd(88,6) =  &
        -k(476)*n(idx_O2Hj)

    !d[E_dot]/d[C]
    pd(1,7) =  &
        +k(215)  &
        +k(206)  &
        +k(973)

    !d[CH_dot]/d[C]
    pd(2,7) =  &
        +k(836)*n(idx_H2)  &
        +k(1051)*n(idx_H)  &
        +k(759)*n(idx_OH)  &
        +2.d0*k(746)*n(idx_CH2)  &
        +k(747)*n(idx_HCO)  &
        +k(751)*n(idx_NH2)  &
        +k(753)*n(idx_NH)

    !d[O_dot]/d[C]
    pd(3,7) =  &
        +k(756)*n(idx_O2)  &
        +k(754)*n(idx_NO)  &
        -k(1044)*n(idx_O)  &
        +k(337)*n(idx_O2j)  &
        +k(759)*n(idx_OH)  &
        +k(339)*n(idx_OHj)

    !d[HNC_dot]/d[C]
    pd(4,7) =  &
        +k(885)*n(idx_HNCO)  &
        +k(750)*n(idx_NH2)

    !d[HCN_dot]/d[C]
    pd(5,7) =  &
        +k(749)*n(idx_NH2)

    !d[H2_dot]/d[C]
    pd(6,7) =  &
        +k(330)*n(idx_H3Oj)  &
        +k(502)*n(idx_H3j)  &
        -k(1048)*n(idx_H2)  &
        -k(836)*n(idx_H2)

    !d[C_dot]/d[C]
    pd(7,7) =  &
        -k(751)*n(idx_NH2)  &
        -k(338)*n(idx_O2Hj)  &
        -k(23)*n(idx_COj)  &
        -k(1043)*n(idx_Oj)  &
        -k(335)*n(idx_N2Hj)  &
        -k(750)*n(idx_NH2)  &
        -k(24)*n(idx_N2j)  &
        -k(757)*n(idx_OCN)  &
        -k(746)*n(idx_CH2)  &
        -k(1042)*n(idx_N)  &
        -k(755)*n(idx_NO)  &
        -k(206)  &
        -k(336)*n(idx_NHj)  &
        -k(973)  &
        -k(25)*n(idx_O2j)  &
        -k(748)*n(idx_N2)  &
        -k(333)*n(idx_HCO2j)  &
        -k(756)*n(idx_O2)  &
        -k(753)*n(idx_NH)  &
        -k(759)*n(idx_OH)  &
        -k(215)  &
        -k(332)*n(idx_HCOj)  &
        -k(122)*n(idx_HEj)  &
        -k(836)*n(idx_H2)  &
        -k(339)*n(idx_OHj)  &
        -k(758)*n(idx_OH)  &
        -k(331)*n(idx_HCNj)  &
        -k(1051)*n(idx_H)  &
        -k(502)*n(idx_H3j)  &
        -k(329)*n(idx_H2Oj)  &
        -k(441)*n(idx_H2j)  &
        -k(752)*n(idx_NH)  &
        -k(747)*n(idx_HCO)  &
        -k(337)*n(idx_O2j)  &
        -k(334)*n(idx_HNOj)  &
        -k(754)*n(idx_NO)  &
        -k(1048)*n(idx_H2)  &
        -k(1044)*n(idx_O)  &
        -k(22)*n(idx_CNj)  &
        -k(1070)  &
        -k(885)*n(idx_HNCO)  &
        -k(330)*n(idx_H3Oj)  &
        -k(749)*n(idx_NH2)

    !d[H_dot]/d[C]
    pd(8,7) =  &
        +k(749)*n(idx_NH2)  &
        +k(836)*n(idx_H2)  &
        -k(1051)*n(idx_H)  &
        +k(750)*n(idx_NH2)  &
        +k(441)*n(idx_H2j)  &
        +k(758)*n(idx_OH)  &
        +k(752)*n(idx_NH)

    !d[OH_dot]/d[C]
    pd(10,7) =  &
        -k(758)*n(idx_OH)  &
        +k(329)*n(idx_H2Oj)  &
        -k(759)*n(idx_OH)

    !d[O2_dot]/d[C]
    pd(11,7) =  &
        +k(338)*n(idx_O2Hj)  &
        -k(756)*n(idx_O2)  &
        +k(25)*n(idx_O2j)

    !d[CH2_dot]/d[C]
    pd(12,7) =  &
        -k(746)*n(idx_CH2)  &
        +k(1048)*n(idx_H2)

    !d[HCO_dot]/d[C]
    pd(14,7) =  &
        -k(747)*n(idx_HCO)

    !d[NO_dot]/d[C]
    pd(17,7) =  &
        -k(754)*n(idx_NO)  &
        +k(334)*n(idx_HNOj)  &
        -k(755)*n(idx_NO)

    !d[CN_dot]/d[C]
    pd(18,7) =  &
        +k(757)*n(idx_OCN)  &
        +k(22)*n(idx_CNj)  &
        +k(331)*n(idx_HCNj)  &
        +k(1042)*n(idx_N)  &
        +k(748)*n(idx_N2)  &
        +k(754)*n(idx_NO)  &
        +k(752)*n(idx_NH)

    !d[CO_dot]/d[C]
    pd(19,7) =  &
        +k(757)*n(idx_OCN)  &
        +k(755)*n(idx_NO)  &
        +k(885)*n(idx_HNCO)  &
        +k(756)*n(idx_O2)  &
        +k(332)*n(idx_HCOj)  &
        +k(747)*n(idx_HCO)  &
        +k(23)*n(idx_COj)  &
        +k(1044)*n(idx_O)  &
        +k(758)*n(idx_OH)

    !d[N2_dot]/d[C]
    pd(20,7) =  &
        -k(748)*n(idx_N2)  &
        +k(24)*n(idx_N2j)  &
        +k(335)*n(idx_N2Hj)

    !d[NH2_dot]/d[C]
    pd(21,7) =  &
        -k(751)*n(idx_NH2)  &
        -k(750)*n(idx_NH2)  &
        -k(749)*n(idx_NH2)

    !d[N_dot]/d[C]
    pd(24,7) =  &
        -k(1042)*n(idx_N)  &
        +k(753)*n(idx_NH)  &
        +k(748)*n(idx_N2)  &
        +k(336)*n(idx_NHj)  &
        +k(755)*n(idx_NO)

    !d[NH_dot]/d[C]
    pd(25,7) =  &
        -k(753)*n(idx_NH)  &
        +k(751)*n(idx_NH2)  &
        -k(752)*n(idx_NH)

    !d[HE_dot]/d[C]
    pd(26,7) =  &
        +k(122)*n(idx_HEj)

    !d[CO2_dot]/d[C]
    pd(29,7) =  &
        +k(333)*n(idx_HCO2j)

    !d[HNCO_dot]/d[C]
    pd(31,7) =  &
        -k(885)*n(idx_HNCO)

    !d[OCN_dot]/d[C]
    pd(34,7) =  &
        -k(757)*n(idx_OCN)

    !d[CH4_DUST_dot]/d[C]
    pd(38,7) =  &
        +k(1070)

    !d[HCO+_dot]/d[C]
    pd(54,7) =  &
        -k(332)*n(idx_HCOj)  &
        +k(330)*n(idx_H3Oj)

    !d[C+_dot]/d[C]
    pd(57,7) =  &
        +k(122)*n(idx_HEj)  &
        +k(215)  &
        +k(206)  &
        +k(24)*n(idx_N2j)  &
        +k(973)  &
        +k(23)*n(idx_COj)  &
        +k(22)*n(idx_CNj)  &
        +k(25)*n(idx_O2j)

    !d[CH+_dot]/d[C]
    pd(59,7) =  &
        +k(331)*n(idx_HCNj)  &
        +k(332)*n(idx_HCOj)  &
        +k(335)*n(idx_N2Hj)  &
        +k(338)*n(idx_O2Hj)  &
        +k(333)*n(idx_HCO2j)  &
        +k(334)*n(idx_HNOj)  &
        +k(329)*n(idx_H2Oj)  &
        +k(339)*n(idx_OHj)  &
        +k(502)*n(idx_H3j)  &
        +k(336)*n(idx_NHj)  &
        +k(441)*n(idx_H2j)

    !d[CN+_dot]/d[C]
    pd(64,7) =  &
        -k(22)*n(idx_CNj)

    !d[CO+_dot]/d[C]
    pd(65,7) =  &
        +k(1043)*n(idx_Oj)  &
        +k(337)*n(idx_O2j)  &
        -k(23)*n(idx_COj)

    !d[N2+_dot]/d[C]
    pd(66,7) =  &
        -k(24)*n(idx_N2j)

    !d[O2+_dot]/d[C]
    pd(67,7) =  &
        -k(25)*n(idx_O2j)  &
        -k(337)*n(idx_O2j)

    !d[H2O+_dot]/d[C]
    pd(68,7) =  &
        -k(329)*n(idx_H2Oj)

    !d[O+_dot]/d[C]
    pd(70,7) =  &
        -k(1043)*n(idx_Oj)

    !d[OH+_dot]/d[C]
    pd(71,7) =  &
        -k(339)*n(idx_OHj)

    !d[HCN+_dot]/d[C]
    pd(75,7) =  &
        -k(331)*n(idx_HCNj)

    !d[NH+_dot]/d[C]
    pd(76,7) =  &
        -k(336)*n(idx_NHj)

    !d[H2+_dot]/d[C]
    pd(77,7) =  &
        -k(441)*n(idx_H2j)

    !d[HE+_dot]/d[C]
    pd(78,7) =  &
        -k(122)*n(idx_HEj)

    !d[HNO+_dot]/d[C]
    pd(79,7) =  &
        -k(334)*n(idx_HNOj)

    !d[H3+_dot]/d[C]
    pd(81,7) =  &
        -k(502)*n(idx_H3j)

    !d[H3O+_dot]/d[C]
    pd(83,7) =  &
        -k(330)*n(idx_H3Oj)

    !d[HCO2+_dot]/d[C]
    pd(85,7) =  &
        -k(333)*n(idx_HCO2j)

    !d[N2H+_dot]/d[C]
    pd(87,7) =  &
        -k(335)*n(idx_N2Hj)

    !d[O2H+_dot]/d[C]
    pd(88,7) =  &
        -k(338)*n(idx_O2Hj)

    !d[E_dot]/d[H]
    pd(1,8) =  &
        +k(211)  &
        +k(232)

    !d[CH_dot]/d[H]
    pd(2,8) =  &
        +k(1051)*n(idx_C)  &
        -k(851)*n(idx_CH)  &
        +k(848)*n(idx_CH2)  &
        -k(10)*n(idx_CH)

    !d[O_dot]/d[H]
    pd(3,8) =  &
        -k(1052)*n(idx_O)  &
        +k(859)*n(idx_HCO)  &
        +k(874)*n(idx_OCN)  &
        +k(870)*n(idx_O2)  &
        +k(861)*n(idx_HNO)  &
        +k(14)*n(idx_OH)  &
        +k(877)*n(idx_OH)  &
        +2.d0*k(13)*n(idx_O2)  &
        +k(871)*n(idx_O2H)  &
        +k(868)*n(idx_NO)  &
        +k(115)*n(idx_Oj)

    !d[HNC_dot]/d[H]
    pd(4,8) =  &
        -k(860)*n(idx_HNC)

    !d[HCN_dot]/d[H]
    pd(5,8) =  &
        -k(857)*n(idx_HCN)  &
        +k(854)*n(idx_H2CN)  &
        +k(860)*n(idx_HNC)  &
        +k(874)*n(idx_OCN)  &
        +k(113)*n(idx_HCNj)

    !d[H2_dot]/d[H]
    pd(6,8) =  &
        +k(858)*n(idx_HCO)  &
        +k(530)*n(idx_CHj)  &
        +k(855)*n(idx_H2CO)  &
        +k(857)*n(idx_HCN)  &
        +k(856)*n(idx_H2O)  &
        +k(865)*n(idx_NH3)  &
        +k(533)*n(idx_CH4j)  &
        -k(11)*n(idx_H2)  &
        +k(862)*n(idx_HNO)  &
        +k(872)*n(idx_O2H)  &
        +k(112)*n(idx_H2j)  &
        +k(850)*n(idx_CH4)  &
        +k(854)*n(idx_H2CN)  &
        +k(864)*n(idx_NH2)  &
        +k(532)*n(idx_CH3j)  &
        +k(877)*n(idx_OH)  &
        +k(531)*n(idx_CH2j)  &
        +k(849)*n(idx_CH3)  &
        +k(866)*n(idx_NH)  &
        +k(851)*n(idx_CH)  &
        +k(848)*n(idx_CH2)

    !d[C_dot]/d[H]
    pd(7,8) =  &
        +k(853)*n(idx_CO)  &
        +k(10)*n(idx_CH)  &
        +k(851)*n(idx_CH)  &
        -k(1051)*n(idx_C)

    !d[H_dot]/d[H]
    pd(8,8) =  &
        -k(849)*n(idx_CH3)  &
        -k(870)*n(idx_O2)  &
        -k(532)*n(idx_CH3j)  &
        -k(875)*n(idx_OCN)  &
        -k(866)*n(idx_NH)  &
        -k(534)*n(idx_HEHj)  &
        -k(848)*n(idx_CH2)  &
        +k(860)*n(idx_HNC)  &
        -k(12)*n(idx_H2O)  &
        -k(115)*n(idx_Oj)  &
        -k(865)*n(idx_NH3)  &
        -k(10)*n(idx_CH)  &
        -k(110)*n(idx_CNj)  &
        -k(860)*n(idx_HNC)  &
        -k(1050)*n(idx_Cj)  &
        +2.d0*k(10)*n(idx_CH)  &
        -k(531)*n(idx_CH2j)  &
        -k(868)*n(idx_NO)  &
        -k(855)*n(idx_H2CO)  &
        -k(862)*n(idx_HNO)  &
        -k(232)  &
        -k(876)*n(idx_OCN)  &
        -k(111)*n(idx_COj)  &
        +k(13)*n(idx_O2)  &
        -k(871)*n(idx_O2H)  &
        -k(863)*n(idx_HNO)  &
        -k(533)*n(idx_CH4j)  &
        -k(1053)*n(idx_OH)  &
        -k(856)*n(idx_H2O)  &
        -k(869)*n(idx_NO)  &
        -k(112)*n(idx_H2j)  &
        -k(13)*n(idx_O2)  &
        -k(852)*n(idx_CO2)  &
        -k(874)*n(idx_OCN)  &
        -k(858)*n(idx_HCO)  &
        -k(859)*n(idx_HCO)  &
        +3.d0*k(11)*n(idx_H2)  &
        -k(872)*n(idx_O2H)  &
        -k(877)*n(idx_OH)  &
        -k(854)*n(idx_H2CN)  &
        -k(1052)*n(idx_O)  &
        -k(114)*n(idx_HEj)  &
        -k(1045)*n(idx_Hj)  &
        -k(850)*n(idx_CH4)  &
        -k(14)*n(idx_OH)  &
        -k(11)*n(idx_H2)  &
        -k(857)*n(idx_HCN)  &
        -k(853)*n(idx_CO)  &
        -k(851)*n(idx_CH)  &
        -k(113)*n(idx_HCNj)  &
        -k(530)*n(idx_CHj)  &
        -k(861)*n(idx_HNO)  &
        -k(864)*n(idx_NH2)  &
        +2.d0*k(14)*n(idx_OH)  &
        -k(211)  &
        -k(873)*n(idx_O2H)  &
        +2.d0*k(12)*n(idx_H2O)  &
        -k(867)*n(idx_NO2)  &
        -k(1051)*n(idx_C)

    !d[H2O_dot]/d[H]
    pd(9,8) =  &
        -k(12)*n(idx_H2O)  &
        -k(856)*n(idx_H2O)  &
        +k(1053)*n(idx_OH)  &
        +k(871)*n(idx_O2H)

    !d[OH_dot]/d[H]
    pd(10,8) =  &
        +k(856)*n(idx_H2O)  &
        +k(870)*n(idx_O2)  &
        -k(14)*n(idx_OH)  &
        +k(876)*n(idx_OCN)  &
        +k(1052)*n(idx_O)  &
        +k(852)*n(idx_CO2)  &
        -k(1053)*n(idx_OH)  &
        -k(877)*n(idx_OH)  &
        +k(867)*n(idx_NO2)  &
        +k(853)*n(idx_CO)  &
        +k(12)*n(idx_H2O)  &
        +k(863)*n(idx_HNO)  &
        +k(869)*n(idx_NO)  &
        +2.d0*k(873)*n(idx_O2H)

    !d[O2_dot]/d[H]
    pd(11,8) =  &
        -k(13)*n(idx_O2)  &
        +k(872)*n(idx_O2H)  &
        -k(870)*n(idx_O2)

    !d[CH2_dot]/d[H]
    pd(12,8) =  &
        +k(859)*n(idx_HCO)  &
        +k(849)*n(idx_CH3)  &
        -k(848)*n(idx_CH2)

    !d[H2CO_dot]/d[H]
    pd(13,8) =  &
        -k(855)*n(idx_H2CO)

    !d[HCO_dot]/d[H]
    pd(14,8) =  &
        -k(858)*n(idx_HCO)  &
        -k(859)*n(idx_HCO)  &
        +k(855)*n(idx_H2CO)

    !d[NH3_dot]/d[H]
    pd(16,8) =  &
        -k(865)*n(idx_NH3)

    !d[NO_dot]/d[H]
    pd(17,8) =  &
        +k(867)*n(idx_NO2)  &
        +k(862)*n(idx_HNO)  &
        -k(868)*n(idx_NO)  &
        -k(869)*n(idx_NO)

    !d[CN_dot]/d[H]
    pd(18,8) =  &
        +k(110)*n(idx_CNj)  &
        +k(857)*n(idx_HCN)  &
        +k(876)*n(idx_OCN)

    !d[CO_dot]/d[H]
    pd(19,8) =  &
        +k(858)*n(idx_HCO)  &
        +k(875)*n(idx_OCN)  &
        +k(852)*n(idx_CO2)  &
        +k(111)*n(idx_COj)  &
        -k(853)*n(idx_CO)

    !d[NH2_dot]/d[H]
    pd(21,8) =  &
        +k(861)*n(idx_HNO)  &
        -k(864)*n(idx_NH2)  &
        +k(865)*n(idx_NH3)

    !d[CH3_dot]/d[H]
    pd(22,8) =  &
        +k(850)*n(idx_CH4)  &
        -k(849)*n(idx_CH3)

    !d[CH4_dot]/d[H]
    pd(23,8) =  &
        -k(850)*n(idx_CH4)

    !d[N_dot]/d[H]
    pd(24,8) =  &
        +k(866)*n(idx_NH)  &
        +k(869)*n(idx_NO)

    !d[NH_dot]/d[H]
    pd(25,8) =  &
        +k(864)*n(idx_NH2)  &
        +k(868)*n(idx_NO)  &
        -k(866)*n(idx_NH)  &
        +k(863)*n(idx_HNO)  &
        +k(875)*n(idx_OCN)

    !d[HE_dot]/d[H]
    pd(26,8) =  &
        +k(534)*n(idx_HEHj)  &
        +k(114)*n(idx_HEj)

    !d[HNO_dot]/d[H]
    pd(27,8) =  &
        -k(861)*n(idx_HNO)  &
        -k(862)*n(idx_HNO)  &
        -k(863)*n(idx_HNO)

    !d[CO2_dot]/d[H]
    pd(29,8) =  &
        -k(852)*n(idx_CO2)

    !d[H2CN_dot]/d[H]
    pd(30,8) =  &
        -k(854)*n(idx_H2CN)

    !d[NO2_dot]/d[H]
    pd(32,8) =  &
        -k(867)*n(idx_NO2)

    !d[O2H_dot]/d[H]
    pd(33,8) =  &
        -k(873)*n(idx_O2H)  &
        -k(872)*n(idx_O2H)  &
        -k(871)*n(idx_O2H)

    !d[OCN_dot]/d[H]
    pd(34,8) =  &
        -k(876)*n(idx_OCN)  &
        -k(875)*n(idx_OCN)  &
        -k(874)*n(idx_OCN)

    !d[H+_dot]/d[H]
    pd(55,8) =  &
        -k(1045)*n(idx_Hj)  &
        +k(114)*n(idx_HEj)  &
        +k(112)*n(idx_H2j)  &
        +k(113)*n(idx_HCNj)  &
        +k(110)*n(idx_CNj)  &
        +k(211)  &
        +k(232)  &
        +k(111)*n(idx_COj)  &
        +k(115)*n(idx_Oj)

    !d[C+_dot]/d[H]
    pd(57,8) =  &
        +k(530)*n(idx_CHj)  &
        -k(1050)*n(idx_Cj)

    !d[CH2+_dot]/d[H]
    pd(58,8) =  &
        +k(532)*n(idx_CH3j)  &
        -k(531)*n(idx_CH2j)

    !d[CH+_dot]/d[H]
    pd(59,8) =  &
        +k(531)*n(idx_CH2j)  &
        +k(1050)*n(idx_Cj)  &
        -k(530)*n(idx_CHj)

    !d[CN+_dot]/d[H]
    pd(64,8) =  &
        -k(110)*n(idx_CNj)

    !d[CO+_dot]/d[H]
    pd(65,8) =  &
        -k(111)*n(idx_COj)

    !d[O+_dot]/d[H]
    pd(70,8) =  &
        -k(115)*n(idx_Oj)

    !d[CH3+_dot]/d[H]
    pd(72,8) =  &
        +k(533)*n(idx_CH4j)  &
        -k(532)*n(idx_CH3j)

    !d[CH4+_dot]/d[H]
    pd(73,8) =  &
        -k(533)*n(idx_CH4j)

    !d[HCN+_dot]/d[H]
    pd(75,8) =  &
        -k(113)*n(idx_HCNj)

    !d[H2+_dot]/d[H]
    pd(77,8) =  &
        +k(534)*n(idx_HEHj)  &
        +k(1045)*n(idx_Hj)  &
        -k(112)*n(idx_H2j)

    !d[HE+_dot]/d[H]
    pd(78,8) =  &
        -k(114)*n(idx_HEj)

    !d[HEH+_dot]/d[H]
    pd(86,8) =  &
        -k(534)*n(idx_HEHj)

    !d[E_dot]/d[H2O]
    pd(1,9) =  &
        +k(1007)

    !d[O_dot]/d[H2O]
    pd(3,9) =  &
        +k(651)*n(idx_NHj)  &
        -k(942)*n(idx_O)  &
        +k(188)*n(idx_Oj)  &
        +k(730)*n(idx_OHj)

    !d[HCN_dot]/d[H2O]
    pd(5,9) =  &
        +k(108)*n(idx_HCNj)

    !d[H2_dot]/d[H2O]
    pd(6,9) =  &
        +k(5)*n(idx_H2)  &
        -k(5)*n(idx_H2)  &
        +k(650)*n(idx_NHj)  &
        +k(856)*n(idx_H)  &
        +k(512)*n(idx_H3j)  &
        +k(348)*n(idx_CHj)  &
        +k(91)*n(idx_H2j)

    !d[C_dot]/d[H2O]
    pd(7,9) =  &
        +k(347)*n(idx_CHj)

    !d[H_dot]/d[H2O]
    pd(8,9) =  &
        +k(319)*n(idx_Cj)  &
        +k(231)  &
        +k(346)*n(idx_CHj)  &
        +k(320)*n(idx_Cj)  &
        +k(1008)  &
        +k(5)*n(idx_H2)  &
        +2.d0*k(12)*n(idx_H)  &
        +k(362)*n(idx_CH2j)  &
        +k(584)*n(idx_HEj)  &
        -k(12)*n(idx_H)  &
        +k(450)*n(idx_H2j)  &
        -k(856)*n(idx_H)  &
        +k(74)*n(idx_Hj)

    !d[H2O_dot]/d[H2O]
    pd(9,9) =  &
        -k(320)*n(idx_Cj)  &
        -k(198)*n(idx_OHj)  &
        -k(108)*n(idx_HCNj)  &
        -k(500)*n(idx_N2Hj)  &
        -k(1008)  &
        -k(650)*n(idx_NHj)  &
        -k(5)*n(idx_H2)  &
        -k(491)*n(idx_CNj)  &
        -k(392)*n(idx_CH4j)  &
        -k(512)*n(idx_H3j)  &
        -k(490)*n(idx_CNj)  &
        -k(652)*n(idx_NHj)  &
        -k(362)*n(idx_CH2j)  &
        -k(916)*n(idx_NH)  &
        -k(450)*n(idx_H2j)  &
        -k(651)*n(idx_NHj)  &
        -k(1077)  &
        -k(347)*n(idx_CHj)  &
        -k(107)*n(idx_COj)  &
        -k(126)*n(idx_HEj)  &
        -k(942)*n(idx_O)  &
        -k(156)*n(idx_NHj)  &
        -k(485)*n(idx_H2Oj)  &
        -k(493)*n(idx_H2COj)  &
        -k(492)*n(idx_COj)  &
        -k(188)*n(idx_Oj)  &
        -k(91)*n(idx_H2j)  &
        -k(585)*n(idx_HEj)  &
        -k(786)*n(idx_CH3)  &
        -k(856)*n(idx_H)  &
        -k(497)*n(idx_HCO2j)  &
        -k(730)*n(idx_OHj)  &
        -k(667)*n(idx_NH2j)  &
        -k(499)*n(idx_N2j)  &
        -k(12)*n(idx_H)  &
        -k(584)*n(idx_HEj)  &
        -k(140)*n(idx_Nj)  &
        -k(494)*n(idx_H3COj)  &
        -k(109)*n(idx_N2j)  &
        -k(348)*n(idx_CHj)  &
        -k(649)*n(idx_NHj)  &
        -k(1007)  &
        -k(74)*n(idx_Hj)  &
        -k(231)  &
        -k(666)*n(idx_NH2j)  &
        -k(501)*n(idx_O2Hj)  &
        -k(346)*n(idx_CHj)  &
        -k(495)*n(idx_HCNj)  &
        -k(319)*n(idx_Cj)  &
        -k(496)*n(idx_HCOj)  &
        -k(498)*n(idx_HNOj)

    !d[OH_dot]/d[H2O]
    pd(10,9) =  &
        +k(198)*n(idx_OHj)  &
        +k(916)*n(idx_NH)  &
        +k(231)  &
        +k(1008)  &
        +2.d0*k(942)*n(idx_O)  &
        +k(5)*n(idx_H2)  &
        +k(12)*n(idx_H)  &
        +k(490)*n(idx_CNj)  &
        +k(485)*n(idx_H2Oj)  &
        +k(499)*n(idx_N2j)  &
        +k(786)*n(idx_CH3)  &
        +k(585)*n(idx_HEj)  &
        +k(667)*n(idx_NH2j)  &
        +k(492)*n(idx_COj)  &
        +k(652)*n(idx_NHj)  &
        +k(856)*n(idx_H)

    !d[O2_dot]/d[H2O]
    pd(11,9) =  &
        +k(501)*n(idx_O2Hj)

    !d[H2CO_dot]/d[H2O]
    pd(13,9) =  &
        +k(494)*n(idx_H3COj)

    !d[HCO_dot]/d[H2O]
    pd(14,9) =  &
        +k(493)*n(idx_H2COj)

    !d[NO_dot]/d[H2O]
    pd(17,9) =  &
        +k(498)*n(idx_HNOj)

    !d[CN_dot]/d[H2O]
    pd(18,9) =  &
        +k(495)*n(idx_HCNj)

    !d[CO_dot]/d[H2O]
    pd(19,9) =  &
        +k(107)*n(idx_COj)  &
        +k(496)*n(idx_HCOj)

    !d[N2_dot]/d[H2O]
    pd(20,9) =  &
        +k(109)*n(idx_N2j)  &
        +k(500)*n(idx_N2Hj)

    !d[NH2_dot]/d[H2O]
    pd(21,9) =  &
        +k(916)*n(idx_NH)

    !d[CH3_dot]/d[H2O]
    pd(22,9) =  &
        +k(392)*n(idx_CH4j)  &
        -k(786)*n(idx_CH3)

    !d[CH4_dot]/d[H2O]
    pd(23,9) =  &
        +k(786)*n(idx_CH3)

    !d[N_dot]/d[H2O]
    pd(24,9) =  &
        +k(140)*n(idx_Nj)  &
        +k(649)*n(idx_NHj)

    !d[NH_dot]/d[H2O]
    pd(25,9) =  &
        +k(491)*n(idx_CNj)  &
        +k(156)*n(idx_NHj)  &
        -k(916)*n(idx_NH)  &
        +k(666)*n(idx_NH2j)

    !d[HE_dot]/d[H2O]
    pd(26,9) =  &
        +k(126)*n(idx_HEj)  &
        +k(584)*n(idx_HEj)  &
        +k(585)*n(idx_HEj)

    !d[CO2_dot]/d[H2O]
    pd(29,9) =  &
        +k(497)*n(idx_HCO2j)

    !d[H2O_DUST_dot]/d[H2O]
    pd(40,9) =  &
        +k(1077)

    !d[HCO+_dot]/d[H2O]
    pd(54,9) =  &
        +k(491)*n(idx_CNj)  &
        +k(348)*n(idx_CHj)  &
        +k(492)*n(idx_COj)  &
        -k(496)*n(idx_HCOj)  &
        +k(319)*n(idx_Cj)

    !d[H+_dot]/d[H2O]
    pd(55,9) =  &
        -k(74)*n(idx_Hj)  &
        +k(585)*n(idx_HEj)

    !d[HOC+_dot]/d[H2O]
    pd(56,9) =  &
        +k(320)*n(idx_Cj)

    !d[C+_dot]/d[H2O]
    pd(57,9) =  &
        -k(320)*n(idx_Cj)  &
        -k(319)*n(idx_Cj)

    !d[CH2+_dot]/d[H2O]
    pd(58,9) =  &
        -k(362)*n(idx_CH2j)

    !d[CH+_dot]/d[H2O]
    pd(59,9) =  &
        -k(348)*n(idx_CHj)  &
        -k(346)*n(idx_CHj)  &
        -k(347)*n(idx_CHj)

    !d[H2CO+_dot]/d[H2O]
    pd(60,9) =  &
        +k(346)*n(idx_CHj)  &
        -k(493)*n(idx_H2COj)

    !d[NH3+_dot]/d[H2O]
    pd(62,9) =  &
        +k(651)*n(idx_NHj)  &
        +k(667)*n(idx_NH2j)

    !d[CN+_dot]/d[H2O]
    pd(64,9) =  &
        -k(490)*n(idx_CNj)  &
        -k(491)*n(idx_CNj)

    !d[CO+_dot]/d[H2O]
    pd(65,9) =  &
        -k(107)*n(idx_COj)  &
        -k(492)*n(idx_COj)

    !d[N2+_dot]/d[H2O]
    pd(66,9) =  &
        -k(499)*n(idx_N2j)  &
        -k(109)*n(idx_N2j)

    !d[H2O+_dot]/d[H2O]
    pd(68,9) =  &
        +k(107)*n(idx_COj)  &
        +k(198)*n(idx_OHj)  &
        +k(108)*n(idx_HCNj)  &
        +k(156)*n(idx_NHj)  &
        +k(140)*n(idx_Nj)  &
        +k(188)*n(idx_Oj)  &
        +k(91)*n(idx_H2j)  &
        +k(109)*n(idx_N2j)  &
        +k(1007)  &
        +k(126)*n(idx_HEj)  &
        -k(485)*n(idx_H2Oj)  &
        +k(74)*n(idx_Hj)

    !d[NH2+_dot]/d[H2O]
    pd(69,9) =  &
        -k(667)*n(idx_NH2j)  &
        +k(652)*n(idx_NHj)  &
        -k(666)*n(idx_NH2j)

    !d[O+_dot]/d[H2O]
    pd(70,9) =  &
        -k(188)*n(idx_Oj)

    !d[OH+_dot]/d[H2O]
    pd(71,9) =  &
        -k(730)*n(idx_OHj)  &
        +k(584)*n(idx_HEj)  &
        -k(198)*n(idx_OHj)

    !d[CH4+_dot]/d[H2O]
    pd(73,9) =  &
        -k(392)*n(idx_CH4j)

    !d[N+_dot]/d[H2O]
    pd(74,9) =  &
        -k(140)*n(idx_Nj)

    !d[HCN+_dot]/d[H2O]
    pd(75,9) =  &
        +k(490)*n(idx_CNj)  &
        -k(108)*n(idx_HCNj)  &
        -k(495)*n(idx_HCNj)

    !d[NH+_dot]/d[H2O]
    pd(76,9) =  &
        -k(156)*n(idx_NHj)  &
        -k(651)*n(idx_NHj)  &
        -k(649)*n(idx_NHj)  &
        -k(650)*n(idx_NHj)  &
        -k(652)*n(idx_NHj)

    !d[H2+_dot]/d[H2O]
    pd(77,9) =  &
        -k(91)*n(idx_H2j)  &
        -k(450)*n(idx_H2j)

    !d[HE+_dot]/d[H2O]
    pd(78,9) =  &
        -k(584)*n(idx_HEj)  &
        -k(585)*n(idx_HEj)  &
        -k(126)*n(idx_HEj)

    !d[HNO+_dot]/d[H2O]
    pd(79,9) =  &
        -k(498)*n(idx_HNOj)  &
        +k(650)*n(idx_NHj)

    !d[H3+_dot]/d[H2O]
    pd(81,9) =  &
        -k(512)*n(idx_H3j)

    !d[H3CO+_dot]/d[H2O]
    pd(82,9) =  &
        -k(494)*n(idx_H3COj)  &
        +k(362)*n(idx_CH2j)

    !d[H3O+_dot]/d[H2O]
    pd(83,9) =  &
        +k(500)*n(idx_N2Hj)  &
        +k(495)*n(idx_HCNj)  &
        +k(666)*n(idx_NH2j)  &
        +k(649)*n(idx_NHj)  &
        +k(485)*n(idx_H2Oj)  &
        +k(497)*n(idx_HCO2j)  &
        +k(347)*n(idx_CHj)  &
        +k(392)*n(idx_CH4j)  &
        +k(501)*n(idx_O2Hj)  &
        +k(512)*n(idx_H3j)  &
        +k(450)*n(idx_H2j)  &
        +k(493)*n(idx_H2COj)  &
        +k(498)*n(idx_HNOj)  &
        +k(494)*n(idx_H3COj)  &
        +k(496)*n(idx_HCOj)  &
        +k(730)*n(idx_OHj)

    !d[HCO2+_dot]/d[H2O]
    pd(85,9) =  &
        -k(497)*n(idx_HCO2j)

    !d[N2H+_dot]/d[H2O]
    pd(87,9) =  &
        +k(499)*n(idx_N2j)  &
        -k(500)*n(idx_N2Hj)

    !d[O2H+_dot]/d[H2O]
    pd(88,9) =  &
        -k(501)*n(idx_O2Hj)

    !d[E_dot]/d[OH]
    pd(1,10) =  &
        +k(1039)

    !d[CH_dot]/d[OH]
    pd(2,10) =  &
        +k(781)*n(idx_CH2)  &
        +k(759)*n(idx_C)  &
        -k(823)*n(idx_CH)

    !d[O_dot]/d[OH]
    pd(3,10) =  &
        +k(961)*n(idx_CN)  &
        +k(739)*n(idx_H2Oj)  &
        +k(930)*n(idx_NH)  &
        +k(1038)  &
        -k(960)*n(idx_O)  &
        +k(799)*n(idx_CH3)  &
        +k(258)  &
        +k(782)*n(idx_CH2)  &
        +2.d0*k(972)*n(idx_OH)  &
        +k(738)*n(idx_COj)  &
        +k(912)*n(idx_NH2)  &
        +k(759)*n(idx_C)  &
        +k(193)*n(idx_Oj)  &
        +k(14)*n(idx_H)  &
        +k(907)*n(idx_N)  &
        +k(877)*n(idx_H)  &
        +k(8)*n(idx_H2)  &
        +k(737)*n(idx_OHj)

    !d[HCN_dot]/d[OH]
    pd(5,10) =  &
        +k(961)*n(idx_CN)  &
        -k(966)*n(idx_HCN)  &
        -k(965)*n(idx_HCN)

    !d[H2_dot]/d[OH]
    pd(6,10) =  &
        +k(388)*n(idx_CH3j)  &
        -k(847)*n(idx_H2)  &
        +k(359)*n(idx_CHj)  &
        +k(526)*n(idx_H3j)  &
        -k(8)*n(idx_H2)  &
        +k(800)*n(idx_CH3)  &
        +k(99)*n(idx_H2j)  &
        +k(877)*n(idx_H)  &
        +k(8)*n(idx_H2)

    !d[C_dot]/d[OH]
    pd(7,10) =  &
        -k(759)*n(idx_C)  &
        -k(758)*n(idx_C)

    !d[H_dot]/d[OH]
    pd(8,10) =  &
        -k(877)*n(idx_H)  &
        +k(1038)  &
        +k(929)*n(idx_NH)  &
        -k(14)*n(idx_H)  &
        +k(328)*n(idx_Cj)  &
        +k(780)*n(idx_CH2)  &
        +k(258)  &
        +k(84)*n(idx_Hj)  &
        +k(459)*n(idx_H2j)  &
        +k(962)*n(idx_CN)  &
        +k(906)*n(idx_N)  &
        +2.d0*k(14)*n(idx_H)  &
        +k(758)*n(idx_C)  &
        +k(960)*n(idx_O)  &
        +k(847)*n(idx_H2)  &
        +k(8)*n(idx_H2)  &
        +k(714)*n(idx_Oj)  &
        -k(1053)*n(idx_H)  &
        +k(742)*n(idx_HCOj)  &
        +k(963)*n(idx_CO)  &
        +k(823)*n(idx_CH)  &
        +k(609)*n(idx_HEj)  &
        +k(970)*n(idx_NO)

    !d[H2O_dot]/d[OH]
    pd(9,10) =  &
        +k(964)*n(idx_H2CO)  &
        +k(781)*n(idx_CH2)  &
        +k(965)*n(idx_HCN)  &
        +k(968)*n(idx_HNO)  &
        +k(969)*n(idx_NH3)  &
        +k(928)*n(idx_NH)  &
        +k(1053)*n(idx_H)  &
        +k(971)*n(idx_O2H)  &
        +2.d0*k(972)*n(idx_OH)  &
        +k(804)*n(idx_CH4)  &
        +k(847)*n(idx_H2)  &
        +k(801)*n(idx_CH3)  &
        +k(967)*n(idx_HCO)  &
        +k(911)*n(idx_NH2)

    !d[OH_dot]/d[OH]
    pd(10,10) =  &
        -k(877)*n(idx_H)  &
        -k(928)*n(idx_NH)  &
        -k(526)*n(idx_H3j)  &
        -k(969)*n(idx_NH3)  &
        -k(758)*n(idx_C)  &
        -k(971)*n(idx_O2H)  &
        -k(14)*n(idx_H)  &
        -k(911)*n(idx_NH2)  &
        -k(780)*n(idx_CH2)  &
        -k(963)*n(idx_CO)  &
        -k(804)*n(idx_CH4)  &
        -k(149)*n(idx_Nj)  &
        -k(960)*n(idx_O)  &
        -k(801)*n(idx_CH3)  &
        -k(799)*n(idx_CH3)  &
        -k(961)*n(idx_CN)  &
        -k(929)*n(idx_NH)  &
        -k(1038)  &
        -k(912)*n(idx_NH2)  &
        -k(99)*n(idx_H2j)  &
        -k(738)*n(idx_COj)  &
        -k(84)*n(idx_Hj)  &
        -k(388)*n(idx_CH3j)  &
        -k(741)*n(idx_HCOj)  &
        -k(823)*n(idx_CH)  &
        -k(359)*n(idx_CHj)  &
        -k(714)*n(idx_Oj)  &
        -k(742)*n(idx_HCOj)  &
        -k(205)*n(idx_N2j)  &
        -k(800)*n(idx_CH3)  &
        -k(737)*n(idx_OHj)  &
        -k(782)*n(idx_CH2)  &
        -k(930)*n(idx_NH)  &
        -4.d0*k(972)*n(idx_OH)  &
        -k(966)*n(idx_HCN)  &
        -k(459)*n(idx_H2j)  &
        -k(328)*n(idx_Cj)  &
        -k(743)*n(idx_HNOj)  &
        -k(8)*n(idx_H2)  &
        -k(258)  &
        -k(907)*n(idx_N)  &
        -k(193)*n(idx_Oj)  &
        -k(609)*n(idx_HEj)  &
        -k(745)*n(idx_O2Hj)  &
        -k(962)*n(idx_CN)  &
        -k(1053)*n(idx_H)  &
        -k(740)*n(idx_HCNj)  &
        -k(739)*n(idx_H2Oj)  &
        -k(964)*n(idx_H2CO)  &
        -k(663)*n(idx_NHj)  &
        -k(1039)  &
        -k(968)*n(idx_HNO)  &
        -k(759)*n(idx_C)  &
        -k(906)*n(idx_N)  &
        -k(847)*n(idx_H2)  &
        -k(967)*n(idx_HCO)  &
        -k(1074)  &
        -k(965)*n(idx_HCN)  &
        -k(203)*n(idx_CNj)  &
        -k(781)*n(idx_CH2)  &
        -k(204)*n(idx_COj)  &
        -k(744)*n(idx_N2Hj)  &
        -k(970)*n(idx_NO)

    !d[O2_dot]/d[OH]
    pd(11,10) =  &
        +k(971)*n(idx_O2H)  &
        +k(745)*n(idx_O2Hj)  &
        +k(960)*n(idx_O)

    !d[CH2_dot]/d[OH]
    pd(12,10) =  &
        -k(780)*n(idx_CH2)  &
        -k(781)*n(idx_CH2)  &
        +k(801)*n(idx_CH3)  &
        -k(782)*n(idx_CH2)

    !d[H2CO_dot]/d[OH]
    pd(13,10) =  &
        +k(780)*n(idx_CH2)  &
        -k(964)*n(idx_H2CO)  &
        +k(800)*n(idx_CH3)

    !d[HCO_dot]/d[OH]
    pd(14,10) =  &
        +k(964)*n(idx_H2CO)  &
        -k(967)*n(idx_HCO)  &
        +k(823)*n(idx_CH)

    !d[NH3_dot]/d[OH]
    pd(16,10) =  &
        +k(912)*n(idx_NH2)  &
        -k(969)*n(idx_NH3)

    !d[NO_dot]/d[OH]
    pd(17,10) =  &
        +k(743)*n(idx_HNOj)  &
        +k(906)*n(idx_N)  &
        +k(968)*n(idx_HNO)  &
        -k(970)*n(idx_NO)

    !d[CN_dot]/d[OH]
    pd(18,10) =  &
        -k(961)*n(idx_CN)  &
        -k(962)*n(idx_CN)  &
        +k(740)*n(idx_HCNj)  &
        +k(203)*n(idx_CNj)  &
        +k(965)*n(idx_HCN)

    !d[CO_dot]/d[OH]
    pd(19,10) =  &
        -k(963)*n(idx_CO)  &
        +k(204)*n(idx_COj)  &
        +k(741)*n(idx_HCOj)  &
        +k(966)*n(idx_HCN)  &
        +k(758)*n(idx_C)  &
        +k(967)*n(idx_HCO)

    !d[N2_dot]/d[OH]
    pd(20,10) =  &
        +k(205)*n(idx_N2j)  &
        +k(744)*n(idx_N2Hj)

    !d[NH2_dot]/d[OH]
    pd(21,10) =  &
        -k(912)*n(idx_NH2)  &
        +k(969)*n(idx_NH3)  &
        +k(930)*n(idx_NH)  &
        +k(966)*n(idx_HCN)  &
        -k(911)*n(idx_NH2)

    !d[CH3_dot]/d[OH]
    pd(22,10) =  &
        -k(800)*n(idx_CH3)  &
        -k(799)*n(idx_CH3)  &
        +k(782)*n(idx_CH2)  &
        -k(801)*n(idx_CH3)  &
        +k(804)*n(idx_CH4)

    !d[CH4_dot]/d[OH]
    pd(23,10) =  &
        +k(799)*n(idx_CH3)  &
        -k(804)*n(idx_CH4)

    !d[N_dot]/d[OH]
    pd(24,10) =  &
        +k(928)*n(idx_NH)  &
        +k(663)*n(idx_NHj)  &
        -k(907)*n(idx_N)  &
        -k(906)*n(idx_N)  &
        +k(149)*n(idx_Nj)

    !d[NH_dot]/d[OH]
    pd(25,10) =  &
        -k(928)*n(idx_NH)  &
        -k(930)*n(idx_NH)  &
        +k(907)*n(idx_N)  &
        +k(911)*n(idx_NH2)  &
        -k(929)*n(idx_NH)

    !d[HE_dot]/d[OH]
    pd(26,10) =  &
        +k(609)*n(idx_HEj)

    !d[HNO_dot]/d[OH]
    pd(27,10) =  &
        +k(929)*n(idx_NH)  &
        -k(968)*n(idx_HNO)

    !d[CO2_dot]/d[OH]
    pd(29,10) =  &
        +k(963)*n(idx_CO)

    !d[NO2_dot]/d[OH]
    pd(32,10) =  &
        +k(970)*n(idx_NO)

    !d[O2H_dot]/d[OH]
    pd(33,10) =  &
        -k(971)*n(idx_O2H)

    !d[OCN_dot]/d[OH]
    pd(34,10) =  &
        +k(962)*n(idx_CN)

    !d[H2O_DUST_dot]/d[OH]
    pd(40,10) =  &
        +k(1074)

    !d[HCO+_dot]/d[OH]
    pd(54,10) =  &
        +k(738)*n(idx_COj)  &
        -k(741)*n(idx_HCOj)  &
        -k(742)*n(idx_HCOj)

    !d[H+_dot]/d[OH]
    pd(55,10) =  &
        -k(84)*n(idx_Hj)

    !d[C+_dot]/d[OH]
    pd(57,10) =  &
        -k(328)*n(idx_Cj)

    !d[CH+_dot]/d[OH]
    pd(59,10) =  &
        -k(359)*n(idx_CHj)

    !d[H2CO+_dot]/d[OH]
    pd(60,10) =  &
        +k(388)*n(idx_CH3j)

    !d[CN+_dot]/d[OH]
    pd(64,10) =  &
        -k(203)*n(idx_CNj)

    !d[CO+_dot]/d[OH]
    pd(65,10) =  &
        +k(359)*n(idx_CHj)  &
        +k(328)*n(idx_Cj)  &
        -k(204)*n(idx_COj)  &
        -k(738)*n(idx_COj)

    !d[N2+_dot]/d[OH]
    pd(66,10) =  &
        -k(205)*n(idx_N2j)

    !d[O2+_dot]/d[OH]
    pd(67,10) =  &
        +k(714)*n(idx_Oj)

    !d[H2O+_dot]/d[OH]
    pd(68,10) =  &
        +k(745)*n(idx_O2Hj)  &
        -k(739)*n(idx_H2Oj)  &
        +k(741)*n(idx_HCOj)  &
        +k(744)*n(idx_N2Hj)  &
        +k(459)*n(idx_H2j)  &
        +k(663)*n(idx_NHj)  &
        +k(526)*n(idx_H3j)  &
        +k(743)*n(idx_HNOj)  &
        +k(737)*n(idx_OHj)  &
        +k(740)*n(idx_HCNj)

    !d[O+_dot]/d[OH]
    pd(70,10) =  &
        -k(714)*n(idx_Oj)  &
        +k(609)*n(idx_HEj)  &
        -k(193)*n(idx_Oj)

    !d[OH+_dot]/d[OH]
    pd(71,10) =  &
        +k(1039)  &
        +k(149)*n(idx_Nj)  &
        +k(204)*n(idx_COj)  &
        -k(737)*n(idx_OHj)  &
        +k(205)*n(idx_N2j)  &
        +k(84)*n(idx_Hj)  &
        +k(193)*n(idx_Oj)  &
        +k(99)*n(idx_H2j)  &
        +k(203)*n(idx_CNj)

    !d[CH3+_dot]/d[OH]
    pd(72,10) =  &
        -k(388)*n(idx_CH3j)

    !d[N+_dot]/d[OH]
    pd(74,10) =  &
        -k(149)*n(idx_Nj)

    !d[HCN+_dot]/d[OH]
    pd(75,10) =  &
        -k(740)*n(idx_HCNj)

    !d[NH+_dot]/d[OH]
    pd(76,10) =  &
        -k(663)*n(idx_NHj)

    !d[H2+_dot]/d[OH]
    pd(77,10) =  &
        -k(99)*n(idx_H2j)  &
        -k(459)*n(idx_H2j)

    !d[HE+_dot]/d[OH]
    pd(78,10) =  &
        -k(609)*n(idx_HEj)

    !d[HNO+_dot]/d[OH]
    pd(79,10) =  &
        -k(743)*n(idx_HNOj)

    !d[H3+_dot]/d[OH]
    pd(81,10) =  &
        -k(526)*n(idx_H3j)

    !d[H3O+_dot]/d[OH]
    pd(83,10) =  &
        +k(739)*n(idx_H2Oj)

    !d[HCO2+_dot]/d[OH]
    pd(85,10) =  &
        +k(742)*n(idx_HCOj)

    !d[N2H+_dot]/d[OH]
    pd(87,10) =  &
        -k(744)*n(idx_N2Hj)

    !d[O2H+_dot]/d[OH]
    pd(88,10) =  &
        -k(745)*n(idx_O2Hj)

    !d[E_dot]/d[O2]
    pd(1,11) =  &
        +k(253)  &
        +k(1032)

    !d[CH_dot]/d[O2]
    pd(2,11) =  &
        -k(815)*n(idx_CH)  &
        -k(816)*n(idx_CH)  &
        -k(817)*n(idx_CH)  &
        -k(818)*n(idx_CH)

    !d[O_dot]/d[O2]
    pd(3,11) =  &
        +2.d0*k(1033)  &
        +k(870)*n(idx_H)  &
        +k(626)*n(idx_Nj)  &
        +k(192)*n(idx_Oj)  &
        +k(756)*n(idx_C)  &
        +k(356)*n(idx_CHj)  &
        +k(924)*n(idx_NH)  &
        +k(834)*n(idx_CO)  &
        +k(774)*n(idx_CH2)  &
        +k(904)*n(idx_N)  &
        +k(672)*n(idx_NH2j)  &
        +k(818)*n(idx_CH)  &
        +k(831)*n(idx_CN)  &
        +k(325)*n(idx_Cj)  &
        +k(606)*n(idx_HEj)  &
        +k(385)*n(idx_CH3j)  &
        +2.d0*k(7)*n(idx_H2)  &
        +2.d0*k(13)*n(idx_H)  &
        +k(816)*n(idx_CH)  &
        +k(932)*n(idx_NO)  &
        +2.d0*k(254)

    !d[HCN_dot]/d[O2]
    pd(5,11) =  &
        +k(117)*n(idx_HCNj)

    !d[H2_dot]/d[O2]
    pd(6,11) =  &
        -k(7)*n(idx_H2)  &
        -k(844)*n(idx_H2)  &
        -k(845)*n(idx_H2)  &
        +k(98)*n(idx_H2j)  &
        +k(771)*n(idx_CH2)  &
        +k(7)*n(idx_H2)  &
        +k(523)*n(idx_H3j)

    !d[C_dot]/d[O2]
    pd(7,11) =  &
        -k(756)*n(idx_C)

    !d[H_dot]/d[O2]
    pd(8,11) =  &
        +k(82)*n(idx_Hj)  &
        -k(870)*n(idx_H)  &
        +k(844)*n(idx_H2)  &
        +k(815)*n(idx_CH)  &
        +k(816)*n(idx_CH)  &
        +2.d0*k(772)*n(idx_CH2)  &
        -k(13)*n(idx_H)  &
        +k(457)*n(idx_H2j)  &
        +k(13)*n(idx_H)

    !d[H2O_dot]/d[O2]
    pd(9,11) =  &
        +k(773)*n(idx_CH2)  &
        +k(794)*n(idx_CH3)  &
        +k(106)*n(idx_H2Oj)

    !d[OH_dot]/d[O2]
    pd(10,11) =  &
        +k(870)*n(idx_H)  &
        +k(775)*n(idx_CH2)  &
        +k(817)*n(idx_CH)  &
        +2.d0*k(845)*n(idx_H2)  &
        +k(660)*n(idx_NHj)  &
        +k(882)*n(idx_HCO)  &
        +k(925)*n(idx_NH)  &
        +k(364)*n(idx_CH2j)  &
        +k(673)*n(idx_NH2j)  &
        +k(355)*n(idx_CHj)  &
        +k(793)*n(idx_CH3)  &
        +k(202)*n(idx_OHj)

    !d[O2_dot]/d[O2]
    pd(11,11) =  &
        -k(7)*n(idx_H2)  &
        -k(816)*n(idx_CH)  &
        -k(795)*n(idx_CH3)  &
        -k(773)*n(idx_CH2)  &
        -k(818)*n(idx_CH)  &
        -k(771)*n(idx_CH2)  &
        -k(253)  &
        -k(883)*n(idx_HCO)  &
        -k(756)*n(idx_C)  &
        -k(98)*n(idx_H2j)  &
        -k(523)*n(idx_H3j)  &
        -k(803)*n(idx_CH4)  &
        -k(793)*n(idx_CH3)  &
        -k(831)*n(idx_CN)  &
        -k(606)*n(idx_HEj)  &
        -k(192)*n(idx_Oj)  &
        -k(932)*n(idx_NO)  &
        -k(129)*n(idx_HEj)  &
        -k(845)*n(idx_H2)  &
        -k(1112)  &
        -k(117)*n(idx_HCNj)  &
        -k(420)*n(idx_CNj)  &
        -k(45)*n(idx_CH4j)  &
        -k(934)*n(idx_OCN)  &
        -k(834)*n(idx_CO)  &
        -k(355)*n(idx_CHj)  &
        -k(357)*n(idx_CHj)  &
        -k(13)*n(idx_H)  &
        -k(67)*n(idx_COj)  &
        -k(1032)  &
        -k(62)*n(idx_CNj)  &
        -k(661)*n(idx_NHj)  &
        -k(844)*n(idx_H2)  &
        -k(148)*n(idx_Nj)  &
        -k(673)*n(idx_NH2j)  &
        -k(325)*n(idx_Cj)  &
        -k(925)*n(idx_NH)  &
        -k(385)*n(idx_CH3j)  &
        -k(870)*n(idx_H)  &
        -k(202)*n(idx_OHj)  &
        -k(660)*n(idx_NHj)  &
        -k(935)*n(idx_OCN)  &
        -k(364)*n(idx_CH2j)  &
        -k(775)*n(idx_CH2)  &
        -k(82)*n(idx_Hj)  &
        -k(772)*n(idx_CH2)  &
        -k(924)*n(idx_NH)  &
        -k(794)*n(idx_CH3)  &
        -k(153)*n(idx_N2j)  &
        -k(457)*n(idx_H2j)  &
        -k(904)*n(idx_N)  &
        -k(479)*n(idx_H2COj)  &
        -k(626)*n(idx_Nj)  &
        -k(356)*n(idx_CHj)  &
        -k(774)*n(idx_CH2)  &
        -k(882)*n(idx_HCO)  &
        -k(106)*n(idx_H2Oj)  &
        -k(672)*n(idx_NH2j)  &
        -k(817)*n(idx_CH)  &
        -k(326)*n(idx_Cj)  &
        -k(627)*n(idx_Nj)  &
        -k(815)*n(idx_CH)  &
        -k(1033)  &
        -k(159)*n(idx_NHj)  &
        -k(254)  &
        -k(830)*n(idx_CN)

    !d[CH2_dot]/d[O2]
    pd(12,11) =  &
        -k(773)*n(idx_CH2)  &
        -k(771)*n(idx_CH2)  &
        -k(774)*n(idx_CH2)  &
        +k(795)*n(idx_CH3)  &
        -k(775)*n(idx_CH2)  &
        -k(772)*n(idx_CH2)

    !d[H2CO_dot]/d[O2]
    pd(13,11) =  &
        +k(793)*n(idx_CH3)  &
        +k(774)*n(idx_CH2)

    !d[HCO_dot]/d[O2]
    pd(14,11) =  &
        +k(775)*n(idx_CH2)  &
        -k(883)*n(idx_HCO)  &
        -k(882)*n(idx_HCO)  &
        +k(818)*n(idx_CH)  &
        +k(357)*n(idx_CHj)  &
        +k(794)*n(idx_CH3)

    !d[NO_dot]/d[O2]
    pd(17,11) =  &
        +k(627)*n(idx_Nj)  &
        -k(932)*n(idx_NO)  &
        +k(925)*n(idx_NH)  &
        +k(830)*n(idx_CN)  &
        +k(904)*n(idx_N)  &
        +k(934)*n(idx_OCN)

    !d[CN_dot]/d[O2]
    pd(18,11) =  &
        +k(62)*n(idx_CNj)  &
        -k(831)*n(idx_CN)  &
        -k(830)*n(idx_CN)

    !d[CO_dot]/d[O2]
    pd(19,11) =  &
        +k(67)*n(idx_COj)  &
        +k(816)*n(idx_CH)  &
        +k(817)*n(idx_CH)  &
        +k(883)*n(idx_HCO)  &
        -k(834)*n(idx_CO)  &
        +k(420)*n(idx_CNj)  &
        +k(326)*n(idx_Cj)  &
        +k(935)*n(idx_OCN)  &
        +k(830)*n(idx_CN)  &
        +k(773)*n(idx_CH2)  &
        +k(756)*n(idx_C)

    !d[N2_dot]/d[O2]
    pd(20,11) =  &
        +k(153)*n(idx_N2j)

    !d[CH3_dot]/d[O2]
    pd(22,11) =  &
        -k(794)*n(idx_CH3)  &
        +k(803)*n(idx_CH4)  &
        -k(793)*n(idx_CH3)  &
        -k(795)*n(idx_CH3)

    !d[CH4_dot]/d[O2]
    pd(23,11) =  &
        +k(45)*n(idx_CH4j)  &
        -k(803)*n(idx_CH4)

    !d[N_dot]/d[O2]
    pd(24,11) =  &
        -k(904)*n(idx_N)  &
        +k(661)*n(idx_NHj)  &
        +k(148)*n(idx_Nj)

    !d[NH_dot]/d[O2]
    pd(25,11) =  &
        +k(159)*n(idx_NHj)  &
        -k(924)*n(idx_NH)  &
        -k(925)*n(idx_NH)

    !d[HE_dot]/d[O2]
    pd(26,11) =  &
        +k(606)*n(idx_HEj)  &
        +k(129)*n(idx_HEj)

    !d[HNO_dot]/d[O2]
    pd(27,11) =  &
        +k(924)*n(idx_NH)

    !d[CO2_dot]/d[O2]
    pd(29,11) =  &
        +k(815)*n(idx_CH)  &
        +k(882)*n(idx_HCO)  &
        +k(772)*n(idx_CH2)  &
        +k(834)*n(idx_CO)  &
        +k(771)*n(idx_CH2)  &
        +k(934)*n(idx_OCN)

    !d[NO2_dot]/d[O2]
    pd(32,11) =  &
        +k(932)*n(idx_NO)  &
        +k(935)*n(idx_OCN)

    !d[O2H_dot]/d[O2]
    pd(33,11) =  &
        +k(803)*n(idx_CH4)  &
        +k(883)*n(idx_HCO)  &
        +k(844)*n(idx_H2)  &
        +k(479)*n(idx_H2COj)  &
        +k(795)*n(idx_CH3)

    !d[OCN_dot]/d[O2]
    pd(34,11) =  &
        -k(934)*n(idx_OCN)  &
        -k(935)*n(idx_OCN)  &
        +k(831)*n(idx_CN)

    !d[O2_DUST_dot]/d[O2]
    pd(46,11) =  &
        +k(1112)

    !d[HCO+_dot]/d[O2]
    pd(54,11) =  &
        +k(364)*n(idx_CH2j)  &
        +k(356)*n(idx_CHj)  &
        +k(479)*n(idx_H2COj)

    !d[H+_dot]/d[O2]
    pd(55,11) =  &
        -k(82)*n(idx_Hj)

    !d[C+_dot]/d[O2]
    pd(57,11) =  &
        -k(325)*n(idx_Cj)  &
        -k(326)*n(idx_Cj)

    !d[CH2+_dot]/d[O2]
    pd(58,11) =  &
        -k(364)*n(idx_CH2j)

    !d[CH+_dot]/d[O2]
    pd(59,11) =  &
        -k(355)*n(idx_CHj)  &
        -k(357)*n(idx_CHj)  &
        -k(356)*n(idx_CHj)

    !d[H2CO+_dot]/d[O2]
    pd(60,11) =  &
        -k(479)*n(idx_H2COj)

    !d[NO+_dot]/d[O2]
    pd(63,11) =  &
        +k(420)*n(idx_CNj)  &
        +k(626)*n(idx_Nj)  &
        +k(660)*n(idx_NHj)

    !d[CN+_dot]/d[O2]
    pd(64,11) =  &
        -k(420)*n(idx_CNj)  &
        -k(62)*n(idx_CNj)

    !d[CO+_dot]/d[O2]
    pd(65,11) =  &
        +k(355)*n(idx_CHj)  &
        +k(325)*n(idx_Cj)  &
        -k(67)*n(idx_COj)

    !d[N2+_dot]/d[O2]
    pd(66,11) =  &
        -k(153)*n(idx_N2j)

    !d[O2+_dot]/d[O2]
    pd(67,11) =  &
        +k(62)*n(idx_CNj)  &
        +k(67)*n(idx_COj)  &
        +k(82)*n(idx_Hj)  &
        +k(117)*n(idx_HCNj)  &
        +k(148)*n(idx_Nj)  &
        +k(192)*n(idx_Oj)  &
        +k(159)*n(idx_NHj)  &
        +k(153)*n(idx_N2j)  &
        +k(253)  &
        +k(1032)  &
        +k(98)*n(idx_H2j)  &
        +k(45)*n(idx_CH4j)  &
        +k(106)*n(idx_H2Oj)  &
        +k(202)*n(idx_OHj)  &
        +k(129)*n(idx_HEj)

    !d[H2O+_dot]/d[O2]
    pd(68,11) =  &
        -k(106)*n(idx_H2Oj)

    !d[NH2+_dot]/d[O2]
    pd(69,11) =  &
        -k(673)*n(idx_NH2j)  &
        -k(672)*n(idx_NH2j)

    !d[O+_dot]/d[O2]
    pd(70,11) =  &
        +k(627)*n(idx_Nj)  &
        +k(357)*n(idx_CHj)  &
        +k(326)*n(idx_Cj)  &
        -k(192)*n(idx_Oj)  &
        +k(606)*n(idx_HEj)

    !d[OH+_dot]/d[O2]
    pd(71,11) =  &
        -k(202)*n(idx_OHj)

    !d[CH3+_dot]/d[O2]
    pd(72,11) =  &
        -k(385)*n(idx_CH3j)

    !d[CH4+_dot]/d[O2]
    pd(73,11) =  &
        -k(45)*n(idx_CH4j)

    !d[N+_dot]/d[O2]
    pd(74,11) =  &
        -k(627)*n(idx_Nj)  &
        -k(626)*n(idx_Nj)  &
        -k(148)*n(idx_Nj)

    !d[HCN+_dot]/d[O2]
    pd(75,11) =  &
        -k(117)*n(idx_HCNj)

    !d[NH+_dot]/d[O2]
    pd(76,11) =  &
        -k(660)*n(idx_NHj)  &
        -k(159)*n(idx_NHj)  &
        -k(661)*n(idx_NHj)

    !d[H2+_dot]/d[O2]
    pd(77,11) =  &
        -k(98)*n(idx_H2j)  &
        -k(457)*n(idx_H2j)

    !d[HE+_dot]/d[O2]
    pd(78,11) =  &
        -k(129)*n(idx_HEj)  &
        -k(606)*n(idx_HEj)

    !d[HNO+_dot]/d[O2]
    pd(79,11) =  &
        +k(673)*n(idx_NH2j)

    !d[H2NO+_dot]/d[O2]
    pd(80,11) =  &
        +k(672)*n(idx_NH2j)

    !d[H3+_dot]/d[O2]
    pd(81,11) =  &
        -k(523)*n(idx_H3j)

    !d[H3CO+_dot]/d[O2]
    pd(82,11) =  &
        +k(385)*n(idx_CH3j)

    !d[O2H+_dot]/d[O2]
    pd(88,11) =  &
        +k(661)*n(idx_NHj)  &
        +k(457)*n(idx_H2j)  &
        +k(523)*n(idx_H3j)

    !d[E_dot]/d[CH2]
    pd(1,12) =  &
        +k(978)  &
        +k(217)

    !d[CH_dot]/d[CH2]
    pd(2,12) =  &
        +k(888)*n(idx_N)  &
        +2.d0*k(746)*n(idx_C)  &
        +k(848)*n(idx_H)  &
        +k(366)*n(idx_COj)  &
        +k(781)*n(idx_OH)  &
        +k(979)  &
        +k(779)*n(idx_O)  &
        +k(762)*n(idx_CN)  &
        +2.d0*k(760)*n(idx_CH2)  &
        +k(218)

    !d[O_dot]/d[CH2]
    pd(3,12) =  &
        -k(776)*n(idx_O)  &
        +k(774)*n(idx_O2)  &
        +k(37)*n(idx_Oj)  &
        +k(381)*n(idx_OHj)  &
        +k(379)*n(idx_O2j)  &
        -k(777)*n(idx_O)  &
        -k(778)*n(idx_O)  &
        +k(782)*n(idx_OH)  &
        -k(779)*n(idx_O)

    !d[HNC_dot]/d[CH2]
    pd(4,12) =  &
        +k(372)*n(idx_HCNHj)  &
        +k(887)*n(idx_N)

    !d[HCN_dot]/d[CH2]
    pd(5,12) =  &
        +k(769)*n(idx_NO)  &
        +k(766)*n(idx_N2)  &
        +k(762)*n(idx_CN)  &
        +k(371)*n(idx_HCNHj)  &
        +k(886)*n(idx_N)

    !d[H2_dot]/d[CH2]
    pd(6,12) =  &
        +k(85)*n(idx_H2j)  &
        +k(503)*n(idx_H3j)  &
        +k(848)*n(idx_H)  &
        +k(564)*n(idx_HEj)  &
        -k(837)*n(idx_H2)  &
        +k(428)*n(idx_Hj)  &
        +k(776)*n(idx_O)  &
        +k(771)*n(idx_O2)

    !d[C_dot]/d[CH2]
    pd(7,12) =  &
        +k(15)*n(idx_Cj)  &
        -k(746)*n(idx_C)

    !d[H_dot]/d[CH2]
    pd(8,12) =  &
        +k(770)*n(idx_NO)  &
        +k(565)*n(idx_HEj)  &
        +2.d0*k(777)*n(idx_O)  &
        +k(69)*n(idx_Hj)  &
        +k(837)*n(idx_H2)  &
        +k(886)*n(idx_N)  &
        +k(887)*n(idx_N)  &
        +2.d0*k(772)*n(idx_O2)  &
        +k(979)  &
        +k(778)*n(idx_O)  &
        +k(442)*n(idx_H2j)  &
        +k(780)*n(idx_OH)  &
        -k(848)*n(idx_H)  &
        +k(218)

    !d[H2O_dot]/d[CH2]
    pd(9,12) =  &
        +k(369)*n(idx_H3Oj)  &
        +k(34)*n(idx_H2Oj)  &
        +k(773)*n(idx_O2)  &
        +k(781)*n(idx_OH)

    !d[OH_dot]/d[CH2]
    pd(10,12) =  &
        -k(782)*n(idx_OH)  &
        +k(39)*n(idx_OHj)  &
        +k(775)*n(idx_O2)  &
        -k(781)*n(idx_OH)  &
        +k(779)*n(idx_O)  &
        +k(368)*n(idx_H2Oj)  &
        +k(769)*n(idx_NO)  &
        -k(780)*n(idx_OH)

    !d[O2_dot]/d[CH2]
    pd(11,12) =  &
        -k(775)*n(idx_O2)  &
        +k(380)*n(idx_O2Hj)  &
        -k(771)*n(idx_O2)  &
        -k(773)*n(idx_O2)  &
        +k(38)*n(idx_O2j)  &
        -k(774)*n(idx_O2)  &
        -k(772)*n(idx_O2)

    !d[CH2_dot]/d[CH2]
    pd(12,12) =  &
        -k(782)*n(idx_OH)  &
        -k(775)*n(idx_O2)  &
        -k(771)*n(idx_O2)  &
        -k(34)*n(idx_H2Oj)  &
        -k(366)*n(idx_COj)  &
        -k(375)*n(idx_N2Hj)  &
        -k(368)*n(idx_H2Oj)  &
        -k(32)*n(idx_COj)  &
        -k(565)*n(idx_HEj)  &
        -k(31)*n(idx_CNj)  &
        -k(764)*n(idx_HCO)  &
        -k(763)*n(idx_H2CO)  &
        -k(379)*n(idx_O2j)  &
        -k(373)*n(idx_HCOj)  &
        -k(746)*n(idx_C)  &
        -k(37)*n(idx_Oj)  &
        -k(35)*n(idx_N2j)  &
        -k(15)*n(idx_Cj)  &
        -k(780)*n(idx_OH)  &
        -k(376)*n(idx_NHj)  &
        -k(979)  &
        -k(377)*n(idx_NH2j)  &
        -k(371)*n(idx_HCNHj)  &
        -k(33)*n(idx_H2COj)  &
        -k(777)*n(idx_O)  &
        -k(762)*n(idx_CN)  &
        -k(773)*n(idx_O2)  &
        -k(369)*n(idx_H3Oj)  &
        -k(218)  &
        -k(778)*n(idx_O)  &
        -k(767)*n(idx_NO2)  &
        -k(769)*n(idx_NO)  &
        -k(38)*n(idx_O2j)  &
        -k(442)*n(idx_H2j)  &
        -k(887)*n(idx_N)  &
        -k(372)*n(idx_HCNHj)  &
        -k(69)*n(idx_Hj)  &
        -k(217)  &
        -k(978)  &
        -k(564)*n(idx_HEj)  &
        -k(370)*n(idx_HCNj)  &
        -k(135)*n(idx_Nj)  &
        -k(766)*n(idx_N2)  &
        -k(381)*n(idx_OHj)  &
        -4.d0*k(760)*n(idx_CH2)  &
        -k(39)*n(idx_OHj)  &
        -k(85)*n(idx_H2j)  &
        -k(848)*n(idx_H)  &
        -k(779)*n(idx_O)  &
        -k(367)*n(idx_H2COj)  &
        -k(428)*n(idx_Hj)  &
        -k(770)*n(idx_NO)  &
        -k(765)*n(idx_HNO)  &
        -k(1076)  &
        -k(374)*n(idx_HNOj)  &
        -k(503)*n(idx_H3j)  &
        -k(776)*n(idx_O)  &
        -k(774)*n(idx_O2)  &
        -k(781)*n(idx_OH)  &
        -k(36)*n(idx_NH2j)  &
        -k(378)*n(idx_NH3j)  &
        -k(380)*n(idx_O2Hj)  &
        -k(888)*n(idx_N)  &
        -k(886)*n(idx_N)  &
        -k(761)*n(idx_CH4)  &
        -k(772)*n(idx_O2)  &
        -k(768)*n(idx_NO)  &
        -k(837)*n(idx_H2)

    !d[H2CO_dot]/d[CH2]
    pd(13,12) =  &
        -k(763)*n(idx_H2CO)  &
        +k(33)*n(idx_H2COj)  &
        +k(774)*n(idx_O2)  &
        +k(768)*n(idx_NO)  &
        +k(780)*n(idx_OH)  &
        +k(767)*n(idx_NO2)

    !d[HCO_dot]/d[CH2]
    pd(14,12) =  &
        +k(775)*n(idx_O2)  &
        +k(778)*n(idx_O)  &
        +k(763)*n(idx_H2CO)  &
        +k(367)*n(idx_H2COj)  &
        -k(764)*n(idx_HCO)

    !d[NO_dot]/d[CH2]
    pd(17,12) =  &
        -k(770)*n(idx_NO)  &
        -k(768)*n(idx_NO)  &
        +k(765)*n(idx_HNO)  &
        +k(374)*n(idx_HNOj)  &
        -k(769)*n(idx_NO)  &
        +k(767)*n(idx_NO2)

    !d[CN_dot]/d[CH2]
    pd(18,12) =  &
        -k(762)*n(idx_CN)  &
        +k(370)*n(idx_HCNj)  &
        +k(31)*n(idx_CNj)

    !d[CO_dot]/d[CH2]
    pd(19,12) =  &
        +k(32)*n(idx_COj)  &
        +k(777)*n(idx_O)  &
        +k(764)*n(idx_HCO)  &
        +k(776)*n(idx_O)  &
        +k(373)*n(idx_HCOj)  &
        +k(773)*n(idx_O2)

    !d[N2_dot]/d[CH2]
    pd(20,12) =  &
        -k(766)*n(idx_N2)  &
        +k(375)*n(idx_N2Hj)  &
        +k(35)*n(idx_N2j)

    !d[NH2_dot]/d[CH2]
    pd(21,12) =  &
        +k(378)*n(idx_NH3j)  &
        +k(36)*n(idx_NH2j)

    !d[CH3_dot]/d[CH2]
    pd(22,12) =  &
        +2.d0*k(761)*n(idx_CH4)  &
        +k(764)*n(idx_HCO)  &
        +k(837)*n(idx_H2)  &
        +k(765)*n(idx_HNO)  &
        +2.d0*k(760)*n(idx_CH2)  &
        +k(782)*n(idx_OH)  &
        +k(763)*n(idx_H2CO)

    !d[CH4_dot]/d[CH2]
    pd(23,12) =  &
        -k(761)*n(idx_CH4)

    !d[N_dot]/d[CH2]
    pd(24,12) =  &
        +k(376)*n(idx_NHj)  &
        +k(135)*n(idx_Nj)  &
        -k(886)*n(idx_N)  &
        -k(888)*n(idx_N)  &
        +k(768)*n(idx_NO)  &
        -k(887)*n(idx_N)

    !d[NH_dot]/d[CH2]
    pd(25,12) =  &
        +k(377)*n(idx_NH2j)  &
        +k(766)*n(idx_N2)  &
        +k(888)*n(idx_N)

    !d[HE_dot]/d[CH2]
    pd(26,12) =  &
        +k(564)*n(idx_HEj)  &
        +k(565)*n(idx_HEj)

    !d[HNO_dot]/d[CH2]
    pd(27,12) =  &
        -k(765)*n(idx_HNO)

    !d[CO2_dot]/d[CH2]
    pd(29,12) =  &
        +k(771)*n(idx_O2)  &
        +k(772)*n(idx_O2)

    !d[HNCO_dot]/d[CH2]
    pd(31,12) =  &
        +k(770)*n(idx_NO)

    !d[NO2_dot]/d[CH2]
    pd(32,12) =  &
        -k(767)*n(idx_NO2)

    !d[CH4_DUST_dot]/d[CH2]
    pd(38,12) =  &
        +k(1076)

    !d[HCO+_dot]/d[CH2]
    pd(54,12) =  &
        -k(373)*n(idx_HCOj)  &
        +k(366)*n(idx_COj)

    !d[H+_dot]/d[CH2]
    pd(55,12) =  &
        -k(428)*n(idx_Hj)  &
        -k(69)*n(idx_Hj)

    !d[C+_dot]/d[CH2]
    pd(57,12) =  &
        +k(564)*n(idx_HEj)  &
        -k(15)*n(idx_Cj)

    !d[CH2+_dot]/d[CH2]
    pd(58,12) =  &
        +k(33)*n(idx_H2COj)  &
        +k(34)*n(idx_H2Oj)  &
        +k(39)*n(idx_OHj)  &
        +k(35)*n(idx_N2j)  &
        +k(135)*n(idx_Nj)  &
        +k(37)*n(idx_Oj)  &
        +k(38)*n(idx_O2j)  &
        +k(69)*n(idx_Hj)  &
        +k(978)  &
        +k(217)  &
        +k(15)*n(idx_Cj)  &
        +k(32)*n(idx_COj)  &
        +k(31)*n(idx_CNj)  &
        +k(36)*n(idx_NH2j)  &
        +k(85)*n(idx_H2j)

    !d[CH+_dot]/d[CH2]
    pd(59,12) =  &
        +k(565)*n(idx_HEj)  &
        +k(428)*n(idx_Hj)

    !d[H2CO+_dot]/d[CH2]
    pd(60,12) =  &
        -k(33)*n(idx_H2COj)  &
        +k(379)*n(idx_O2j)  &
        -k(367)*n(idx_H2COj)

    !d[NH3+_dot]/d[CH2]
    pd(62,12) =  &
        -k(378)*n(idx_NH3j)

    !d[CN+_dot]/d[CH2]
    pd(64,12) =  &
        -k(31)*n(idx_CNj)

    !d[CO+_dot]/d[CH2]
    pd(65,12) =  &
        -k(32)*n(idx_COj)  &
        -k(366)*n(idx_COj)

    !d[N2+_dot]/d[CH2]
    pd(66,12) =  &
        -k(35)*n(idx_N2j)

    !d[O2+_dot]/d[CH2]
    pd(67,12) =  &
        -k(379)*n(idx_O2j)  &
        -k(38)*n(idx_O2j)

    !d[H2O+_dot]/d[CH2]
    pd(68,12) =  &
        -k(34)*n(idx_H2Oj)  &
        -k(368)*n(idx_H2Oj)

    !d[NH2+_dot]/d[CH2]
    pd(69,12) =  &
        -k(36)*n(idx_NH2j)  &
        -k(377)*n(idx_NH2j)

    !d[O+_dot]/d[CH2]
    pd(70,12) =  &
        -k(37)*n(idx_Oj)

    !d[OH+_dot]/d[CH2]
    pd(71,12) =  &
        -k(381)*n(idx_OHj)  &
        -k(39)*n(idx_OHj)

    !d[CH3+_dot]/d[CH2]
    pd(72,12) =  &
        +k(378)*n(idx_NH3j)  &
        +k(376)*n(idx_NHj)  &
        +k(380)*n(idx_O2Hj)  &
        +k(371)*n(idx_HCNHj)  &
        +k(503)*n(idx_H3j)  &
        +k(381)*n(idx_OHj)  &
        +k(369)*n(idx_H3Oj)  &
        +k(375)*n(idx_N2Hj)  &
        +k(372)*n(idx_HCNHj)  &
        +k(367)*n(idx_H2COj)  &
        +k(377)*n(idx_NH2j)  &
        +k(368)*n(idx_H2Oj)  &
        +k(373)*n(idx_HCOj)  &
        +k(374)*n(idx_HNOj)  &
        +k(370)*n(idx_HCNj)  &
        +k(442)*n(idx_H2j)

    !d[N+_dot]/d[CH2]
    pd(74,12) =  &
        -k(135)*n(idx_Nj)

    !d[HCN+_dot]/d[CH2]
    pd(75,12) =  &
        -k(370)*n(idx_HCNj)

    !d[NH+_dot]/d[CH2]
    pd(76,12) =  &
        -k(376)*n(idx_NHj)

    !d[H2+_dot]/d[CH2]
    pd(77,12) =  &
        -k(442)*n(idx_H2j)  &
        -k(85)*n(idx_H2j)

    !d[HE+_dot]/d[CH2]
    pd(78,12) =  &
        -k(564)*n(idx_HEj)  &
        -k(565)*n(idx_HEj)

    !d[HNO+_dot]/d[CH2]
    pd(79,12) =  &
        -k(374)*n(idx_HNOj)

    !d[H3+_dot]/d[CH2]
    pd(81,12) =  &
        -k(503)*n(idx_H3j)

    !d[H3O+_dot]/d[CH2]
    pd(83,12) =  &
        -k(369)*n(idx_H3Oj)

    !d[HCNH+_dot]/d[CH2]
    pd(84,12) =  &
        -k(372)*n(idx_HCNHj)  &
        -k(371)*n(idx_HCNHj)

    !d[N2H+_dot]/d[CH2]
    pd(87,12) =  &
        -k(375)*n(idx_N2Hj)

    !d[O2H+_dot]/d[CH2]
    pd(88,12) =  &
        -k(380)*n(idx_O2Hj)

    !d[E_dot]/d[H2CO]
    pd(1,13) =  &
        +k(1005)  &
        +k(1004)

    !d[CH_dot]/d[H2CO]
    pd(2,13) =  &
        -k(806)*n(idx_CH)  &
        +k(318)*n(idx_Cj)

    !d[O_dot]/d[H2CO]
    pd(3,13) =  &
        +k(187)*n(idx_Oj)  &
        +k(729)*n(idx_OHj)  &
        +k(583)*n(idx_HEj)  &
        -k(941)*n(idx_O)

    !d[HNC_dot]/d[H2CO]
    pd(4,13) =  &
        +k(549)*n(idx_HCNHj)

    !d[HCN_dot]/d[H2CO]
    pd(5,13) =  &
        +k(418)*n(idx_CNj)  &
        +k(548)*n(idx_HCNHj)  &
        +k(824)*n(idx_CN)

    !d[H2_dot]/d[H2CO]
    pd(6,13) =  &
        +k(230)  &
        +k(449)*n(idx_H2j)  &
        +k(581)*n(idx_HEj)  &
        +k(434)*n(idx_Hj)  &
        +k(90)*n(idx_H2j)  &
        +k(855)*n(idx_H)  &
        +k(1002)  &
        +k(435)*n(idx_Hj)  &
        +k(511)*n(idx_H3j)

    !d[C_dot]/d[H2CO]
    pd(7,13) =  &
        +k(344)*n(idx_CHj)  &
        +k(17)*n(idx_Cj)

    !d[H_dot]/d[H2CO]
    pd(8,13) =  &
        +k(628)*n(idx_N2j)  &
        +k(481)*n(idx_O2j)  &
        +k(449)*n(idx_H2j)  &
        -k(855)*n(idx_H)  &
        +k(434)*n(idx_Hj)  &
        +k(582)*n(idx_HEj)  &
        +k(73)*n(idx_Hj)  &
        +2.d0*k(1003)  &
        +k(1005)

    !d[H2O_dot]/d[H2CO]
    pd(9,13) =  &
        +k(964)*n(idx_OH)  &
        +k(527)*n(idx_H3Oj)  &
        +k(102)*n(idx_H2Oj)

    !d[OH_dot]/d[H2CO]
    pd(10,13) =  &
        +k(197)*n(idx_OHj)  &
        -k(964)*n(idx_OH)  &
        +k(941)*n(idx_O)  &
        +k(484)*n(idx_H2Oj)  &
        +k(708)*n(idx_Oj)

    !d[O2_dot]/d[H2CO]
    pd(11,13) =  &
        +k(481)*n(idx_O2j)  &
        +k(482)*n(idx_O2Hj)  &
        +k(101)*n(idx_O2j)

    !d[CH2_dot]/d[H2CO]
    pd(12,13) =  &
        +k(806)*n(idx_CH)  &
        +k(620)*n(idx_Nj)  &
        -k(763)*n(idx_CH2)  &
        +k(345)*n(idx_CHj)

    !d[H2CO_dot]/d[H2CO]
    pd(13,13) =  &
        -k(17)*n(idx_Cj)  &
        -k(435)*n(idx_Hj)  &
        -k(125)*n(idx_HEj)  &
        -k(187)*n(idx_Oj)  &
        -k(824)*n(idx_CN)  &
        -k(785)*n(idx_CH3)  &
        -k(619)*n(idx_Nj)  &
        -k(511)*n(idx_H3j)  &
        -k(391)*n(idx_CH4j)  &
        -k(361)*n(idx_CH2j)  &
        -k(729)*n(idx_OHj)  &
        -k(434)*n(idx_Hj)  &
        -k(1002)  &
        -k(581)*n(idx_HEj)  &
        -k(423)*n(idx_COj)  &
        -k(550)*n(idx_HCOj)  &
        -k(583)*n(idx_HEj)  &
        -k(139)*n(idx_Nj)  &
        -k(763)*n(idx_CH2)  &
        -k(418)*n(idx_CNj)  &
        -k(230)  &
        -k(480)*n(idx_HNOj)  &
        -k(549)*n(idx_HCNHj)  &
        -k(647)*n(idx_NHj)  &
        -k(964)*n(idx_OH)  &
        -k(345)*n(idx_CHj)  &
        -k(43)*n(idx_CH4j)  &
        -k(478)*n(idx_H2COj)  &
        -k(482)*n(idx_O2Hj)  &
        -k(1004)  &
        -k(58)*n(idx_CNj)  &
        -k(648)*n(idx_NHj)  &
        -k(548)*n(idx_HCNHj)  &
        -k(343)*n(idx_CHj)  &
        -k(317)*n(idx_Cj)  &
        -k(537)*n(idx_HCNj)  &
        -k(582)*n(idx_HEj)  &
        -k(101)*n(idx_O2j)  &
        -k(90)*n(idx_H2j)  &
        -k(318)*n(idx_Cj)  &
        -k(484)*n(idx_H2Oj)  &
        -k(806)*n(idx_CH)  &
        -k(64)*n(idx_COj)  &
        -k(1003)  &
        -k(102)*n(idx_H2Oj)  &
        -k(708)*n(idx_Oj)  &
        -k(481)*n(idx_O2j)  &
        -k(1005)  &
        -k(383)*n(idx_CH3j)  &
        -k(665)*n(idx_NH2j)  &
        -k(941)*n(idx_O)  &
        -k(150)*n(idx_N2j)  &
        -k(664)*n(idx_NH2j)  &
        -k(633)*n(idx_N2Hj)  &
        -k(155)*n(idx_NHj)  &
        -k(620)*n(idx_Nj)  &
        -k(855)*n(idx_H)  &
        -k(344)*n(idx_CHj)  &
        -k(1072)  &
        -k(197)*n(idx_OHj)  &
        -k(527)*n(idx_H3Oj)  &
        -k(449)*n(idx_H2j)  &
        -k(73)*n(idx_Hj)  &
        -k(628)*n(idx_N2j)

    !d[HCO_dot]/d[H2CO]
    pd(14,13) =  &
        +k(785)*n(idx_CH3)  &
        +k(941)*n(idx_O)  &
        +k(665)*n(idx_NH2j)  &
        +k(806)*n(idx_CH)  &
        +k(478)*n(idx_H2COj)  &
        +k(964)*n(idx_OH)  &
        +k(855)*n(idx_H)  &
        +k(824)*n(idx_CN)  &
        +k(423)*n(idx_COj)  &
        +k(763)*n(idx_CH2)

    !d[NO_dot]/d[H2CO]
    pd(17,13) =  &
        +k(480)*n(idx_HNOj)

    !d[CN_dot]/d[H2CO]
    pd(18,13) =  &
        +k(537)*n(idx_HCNj)  &
        -k(824)*n(idx_CN)  &
        +k(58)*n(idx_CNj)

    !d[CO_dot]/d[H2CO]
    pd(19,13) =  &
        +k(230)  &
        +k(1003)  &
        +k(343)*n(idx_CHj)  &
        +k(317)*n(idx_Cj)  &
        +k(1002)  &
        +k(550)*n(idx_HCOj)  &
        +k(64)*n(idx_COj)

    !d[N2_dot]/d[H2CO]
    pd(20,13) =  &
        +k(628)*n(idx_N2j)  &
        +k(150)*n(idx_N2j)  &
        +k(633)*n(idx_N2Hj)

    !d[NH2_dot]/d[H2CO]
    pd(21,13) =  &
        +k(648)*n(idx_NHj)

    !d[CH3_dot]/d[H2CO]
    pd(22,13) =  &
        -k(785)*n(idx_CH3)  &
        +k(391)*n(idx_CH4j)  &
        +k(361)*n(idx_CH2j)  &
        +k(763)*n(idx_CH2)

    !d[CH4_dot]/d[H2CO]
    pd(23,13) =  &
        +k(43)*n(idx_CH4j)  &
        +k(785)*n(idx_CH3)  &
        +k(383)*n(idx_CH3j)

    !d[N_dot]/d[H2CO]
    pd(24,13) =  &
        +k(647)*n(idx_NHj)  &
        +k(139)*n(idx_Nj)

    !d[NH_dot]/d[H2CO]
    pd(25,13) =  &
        +k(664)*n(idx_NH2j)  &
        +k(155)*n(idx_NHj)  &
        +k(619)*n(idx_Nj)

    !d[HE_dot]/d[H2CO]
    pd(26,13) =  &
        +k(581)*n(idx_HEj)  &
        +k(582)*n(idx_HEj)  &
        +k(125)*n(idx_HEj)  &
        +k(583)*n(idx_HEj)

    !d[H2CO_DUST_dot]/d[H2CO]
    pd(37,13) =  &
        +k(1072)

    !d[HCO+_dot]/d[H2CO]
    pd(54,13) =  &
        +k(628)*n(idx_N2j)  &
        +k(648)*n(idx_NHj)  &
        +k(481)*n(idx_O2j)  &
        +k(383)*n(idx_CH3j)  &
        +k(449)*n(idx_H2j)  &
        +k(582)*n(idx_HEj)  &
        +k(361)*n(idx_CH2j)  &
        +k(1005)  &
        +k(418)*n(idx_CNj)  &
        +k(619)*n(idx_Nj)  &
        +k(708)*n(idx_Oj)  &
        +k(435)*n(idx_Hj)  &
        +k(423)*n(idx_COj)  &
        -k(550)*n(idx_HCOj)  &
        +k(318)*n(idx_Cj)  &
        +k(345)*n(idx_CHj)

    !d[H+_dot]/d[H2CO]
    pd(55,13) =  &
        -k(434)*n(idx_Hj)  &
        -k(435)*n(idx_Hj)  &
        -k(73)*n(idx_Hj)

    !d[C+_dot]/d[H2CO]
    pd(57,13) =  &
        -k(17)*n(idx_Cj)  &
        -k(318)*n(idx_Cj)  &
        -k(317)*n(idx_Cj)

    !d[CH2+_dot]/d[H2CO]
    pd(58,13) =  &
        -k(361)*n(idx_CH2j)  &
        +k(583)*n(idx_HEj)  &
        +k(317)*n(idx_Cj)

    !d[CH+_dot]/d[H2CO]
    pd(59,13) =  &
        -k(344)*n(idx_CHj)  &
        -k(345)*n(idx_CHj)  &
        -k(343)*n(idx_CHj)

    !d[H2CO+_dot]/d[H2CO]
    pd(60,13) =  &
        +k(43)*n(idx_CH4j)  &
        +k(197)*n(idx_OHj)  &
        +k(125)*n(idx_HEj)  &
        +k(90)*n(idx_H2j)  &
        +k(58)*n(idx_CNj)  &
        +k(102)*n(idx_H2Oj)  &
        -k(478)*n(idx_H2COj)  &
        +k(73)*n(idx_Hj)  &
        +k(1004)  &
        +k(155)*n(idx_NHj)  &
        +k(101)*n(idx_O2j)  &
        +k(187)*n(idx_Oj)  &
        +k(17)*n(idx_Cj)  &
        +k(139)*n(idx_Nj)  &
        +k(150)*n(idx_N2j)  &
        +k(64)*n(idx_COj)

    !d[NH3+_dot]/d[H2CO]
    pd(62,13) =  &
        +k(665)*n(idx_NH2j)

    !d[NO+_dot]/d[H2CO]
    pd(63,13) =  &
        +k(620)*n(idx_Nj)

    !d[CN+_dot]/d[H2CO]
    pd(64,13) =  &
        -k(58)*n(idx_CNj)  &
        -k(418)*n(idx_CNj)

    !d[CO+_dot]/d[H2CO]
    pd(65,13) =  &
        -k(64)*n(idx_COj)  &
        +k(434)*n(idx_Hj)  &
        +k(581)*n(idx_HEj)  &
        -k(423)*n(idx_COj)

    !d[N2+_dot]/d[H2CO]
    pd(66,13) =  &
        -k(150)*n(idx_N2j)  &
        -k(628)*n(idx_N2j)

    !d[O2+_dot]/d[H2CO]
    pd(67,13) =  &
        -k(481)*n(idx_O2j)  &
        -k(101)*n(idx_O2j)

    !d[H2O+_dot]/d[H2CO]
    pd(68,13) =  &
        -k(484)*n(idx_H2Oj)  &
        -k(102)*n(idx_H2Oj)

    !d[NH2+_dot]/d[H2CO]
    pd(69,13) =  &
        -k(665)*n(idx_NH2j)  &
        -k(664)*n(idx_NH2j)

    !d[O+_dot]/d[H2CO]
    pd(70,13) =  &
        -k(708)*n(idx_Oj)  &
        -k(187)*n(idx_Oj)

    !d[OH+_dot]/d[H2CO]
    pd(71,13) =  &
        -k(729)*n(idx_OHj)  &
        -k(197)*n(idx_OHj)

    !d[CH3+_dot]/d[H2CO]
    pd(72,13) =  &
        -k(383)*n(idx_CH3j)  &
        +k(343)*n(idx_CHj)

    !d[CH4+_dot]/d[H2CO]
    pd(73,13) =  &
        -k(391)*n(idx_CH4j)  &
        -k(43)*n(idx_CH4j)

    !d[N+_dot]/d[H2CO]
    pd(74,13) =  &
        -k(139)*n(idx_Nj)  &
        -k(619)*n(idx_Nj)  &
        -k(620)*n(idx_Nj)

    !d[HCN+_dot]/d[H2CO]
    pd(75,13) =  &
        -k(537)*n(idx_HCNj)

    !d[NH+_dot]/d[H2CO]
    pd(76,13) =  &
        -k(647)*n(idx_NHj)  &
        -k(648)*n(idx_NHj)  &
        -k(155)*n(idx_NHj)

    !d[H2+_dot]/d[H2CO]
    pd(77,13) =  &
        -k(90)*n(idx_H2j)  &
        -k(449)*n(idx_H2j)

    !d[HE+_dot]/d[H2CO]
    pd(78,13) =  &
        -k(583)*n(idx_HEj)  &
        -k(582)*n(idx_HEj)  &
        -k(125)*n(idx_HEj)  &
        -k(581)*n(idx_HEj)

    !d[HNO+_dot]/d[H2CO]
    pd(79,13) =  &
        -k(480)*n(idx_HNOj)

    !d[H3+_dot]/d[H2CO]
    pd(81,13) =  &
        -k(511)*n(idx_H3j)

    !d[H3CO+_dot]/d[H2CO]
    pd(82,13) =  &
        +k(537)*n(idx_HCNj)  &
        +k(391)*n(idx_CH4j)  &
        +k(633)*n(idx_N2Hj)  &
        +k(480)*n(idx_HNOj)  &
        +k(664)*n(idx_NH2j)  &
        +k(482)*n(idx_O2Hj)  &
        +k(484)*n(idx_H2Oj)  &
        +k(478)*n(idx_H2COj)  &
        +k(548)*n(idx_HCNHj)  &
        +k(647)*n(idx_NHj)  &
        +k(344)*n(idx_CHj)  &
        +k(729)*n(idx_OHj)  &
        +k(550)*n(idx_HCOj)  &
        +k(511)*n(idx_H3j)  &
        +k(527)*n(idx_H3Oj)  &
        +k(549)*n(idx_HCNHj)

    !d[H3O+_dot]/d[H2CO]
    pd(83,13) =  &
        -k(527)*n(idx_H3Oj)

    !d[HCNH+_dot]/d[H2CO]
    pd(84,13) =  &
        -k(549)*n(idx_HCNHj)  &
        -k(548)*n(idx_HCNHj)

    !d[N2H+_dot]/d[H2CO]
    pd(87,13) =  &
        -k(633)*n(idx_N2Hj)

    !d[O2H+_dot]/d[H2CO]
    pd(88,13) =  &
        -k(482)*n(idx_O2Hj)

    !d[E_dot]/d[HCO]
    pd(1,14) =  &
        +k(235)  &
        +k(1014)

    !d[CH_dot]/d[HCO]
    pd(2,14) =  &
        +k(747)*n(idx_C)  &
        -k(807)*n(idx_CH)  &
        +k(26)*n(idx_CHj)

    !d[O_dot]/d[HCO]
    pd(3,14) =  &
        -k(947)*n(idx_O)  &
        +k(592)*n(idx_HEj)  &
        +k(859)*n(idx_H)  &
        +k(733)*n(idx_OHj)  &
        +k(189)*n(idx_Oj)  &
        -k(946)*n(idx_O)  &
        +k(896)*n(idx_N)

    !d[HCN_dot]/d[HCO]
    pd(5,14) =  &
        +k(896)*n(idx_N)  &
        +k(825)*n(idx_CN)

    !d[H2_dot]/d[HCO]
    pd(6,14) =  &
        +k(514)*n(idx_H3j)  &
        +k(858)*n(idx_H)  &
        +k(436)*n(idx_Hj)  &
        +k(93)*n(idx_H2j)  &
        +2.d0*k(878)*n(idx_HCO)

    !d[C_dot]/d[HCO]
    pd(7,14) =  &
        -k(747)*n(idx_C)  &
        +k(18)*n(idx_Cj)

    !d[H_dot]/d[HCO]
    pd(8,14) =  &
        +k(76)*n(idx_Hj)  &
        +k(1013)  &
        -k(858)*n(idx_H)  &
        +k(897)*n(idx_N)  &
        -k(859)*n(idx_H)  &
        +k(234)  &
        +k(946)*n(idx_O)  &
        +k(590)*n(idx_HEj)

    !d[H2O_dot]/d[HCO]
    pd(9,14) =  &
        +k(967)*n(idx_OH)  &
        +k(103)*n(idx_H2Oj)

    !d[OH_dot]/d[HCO]
    pd(10,14) =  &
        +k(199)*n(idx_OHj)  &
        +k(882)*n(idx_O2)  &
        +k(488)*n(idx_H2Oj)  &
        -k(967)*n(idx_OH)  &
        +k(947)*n(idx_O)

    !d[O2_dot]/d[HCO]
    pd(11,14) =  &
        +k(121)*n(idx_O2j)  &
        -k(882)*n(idx_O2)  &
        +k(884)*n(idx_O2H)  &
        +k(556)*n(idx_O2Hj)  &
        -k(883)*n(idx_O2)

    !d[CH2_dot]/d[HCO]
    pd(12,14) =  &
        +k(807)*n(idx_CH)  &
        +k(859)*n(idx_H)  &
        -k(764)*n(idx_CH2)

    !d[H2CO_dot]/d[HCO]
    pd(13,14) =  &
        +k(880)*n(idx_HNO)  &
        +2.d0*k(879)*n(idx_HCO)  &
        +k(120)*n(idx_H2COj)  &
        +k(884)*n(idx_O2H)

    !d[HCO_dot]/d[HCO]
    pd(14,14) =  &
        -k(189)*n(idx_Oj)  &
        -k(151)*n(idx_N2j)  &
        -k(621)*n(idx_Nj)  &
        -k(76)*n(idx_Hj)  &
        -k(40)*n(idx_CH3j)  &
        -k(947)*n(idx_O)  &
        -k(1081)  &
        -k(858)*n(idx_H)  &
        -k(160)*n(idx_NH2j)  &
        -k(807)*n(idx_CH)  &
        -k(18)*n(idx_Cj)  &
        -k(487)*n(idx_H2Oj)  &
        -k(591)*n(idx_HEj)  &
        -k(65)*n(idx_COj)  &
        -k(514)*n(idx_H3j)  &
        -k(946)*n(idx_O)  &
        -k(733)*n(idx_OHj)  &
        -k(1013)  &
        -k(554)*n(idx_N2Hj)  &
        -k(60)*n(idx_CNj)  &
        -k(93)*n(idx_H2j)  &
        -k(881)*n(idx_NO)  &
        -k(363)*n(idx_CH2j)  &
        -k(967)*n(idx_OH)  &
        -k(120)*n(idx_H2COj)  &
        -k(142)*n(idx_Nj)  &
        -k(437)*n(idx_Hj)  &
        -k(553)*n(idx_HNOj)  &
        -k(555)*n(idx_O2j)  &
        -4.d0*k(878)*n(idx_HCO)  &
        -k(880)*n(idx_HNO)  &
        -k(551)*n(idx_HCOj)  &
        -k(732)*n(idx_OHj)  &
        -k(825)*n(idx_CN)  &
        -k(436)*n(idx_Hj)  &
        -k(556)*n(idx_O2Hj)  &
        -k(884)*n(idx_O2H)  &
        -k(883)*n(idx_O2)  &
        -k(234)  &
        -k(321)*n(idx_Cj)  &
        -k(711)*n(idx_Oj)  &
        -k(654)*n(idx_NHj)  &
        -k(747)*n(idx_C)  &
        -k(199)*n(idx_OHj)  &
        -k(764)*n(idx_CH2)  &
        -k(896)*n(idx_N)  &
        -k(539)*n(idx_HCNj)  &
        -k(787)*n(idx_CH3)  &
        -k(590)*n(idx_HEj)  &
        -k(26)*n(idx_CHj)  &
        -k(103)*n(idx_H2Oj)  &
        -k(488)*n(idx_H2Oj)  &
        -k(882)*n(idx_O2)  &
        -k(895)*n(idx_N)  &
        -k(592)*n(idx_HEj)  &
        -k(552)*n(idx_H2COj)  &
        -k(669)*n(idx_NH2j)  &
        -k(540)*n(idx_HCNj)  &
        -k(629)*n(idx_N2j)  &
        -k(384)*n(idx_CH3j)  &
        -k(169)*n(idx_NH3j)  &
        -k(121)*n(idx_O2j)  &
        -k(859)*n(idx_H)  &
        -k(1014)  &
        -k(897)*n(idx_N)  &
        -k(451)*n(idx_H2j)  &
        -k(419)*n(idx_CNj)  &
        -k(235)  &
        -4.d0*k(879)*n(idx_HCO)  &
        -k(350)*n(idx_CHj)

    !d[NH3_dot]/d[HCO]
    pd(16,14) =  &
        +k(169)*n(idx_NH3j)

    !d[NO_dot]/d[HCO]
    pd(17,14) =  &
        +k(880)*n(idx_HNO)  &
        +k(553)*n(idx_HNOj)  &
        -k(881)*n(idx_NO)

    !d[CN_dot]/d[HCO]
    pd(18,14) =  &
        -k(825)*n(idx_CN)  &
        +k(60)*n(idx_CNj)  &
        +k(539)*n(idx_HCNj)

    !d[CO_dot]/d[HCO]
    pd(19,14) =  &
        +k(967)*n(idx_OH)  &
        +k(732)*n(idx_OHj)  &
        +k(551)*n(idx_HCOj)  &
        +k(419)*n(idx_CNj)  &
        +k(825)*n(idx_CN)  &
        +k(711)*n(idx_Oj)  &
        +k(487)*n(idx_H2Oj)  &
        +k(363)*n(idx_CH2j)  &
        +4.d0*k(878)*n(idx_HCO)  &
        +k(787)*n(idx_CH3)  &
        +k(895)*n(idx_N)  &
        +k(591)*n(idx_HEj)  &
        +k(883)*n(idx_O2)  &
        +k(384)*n(idx_CH3j)  &
        +k(540)*n(idx_HCNj)  &
        +k(65)*n(idx_COj)  &
        +k(234)  &
        +k(552)*n(idx_H2COj)  &
        +k(555)*n(idx_O2j)  &
        +k(881)*n(idx_NO)  &
        +k(451)*n(idx_H2j)  &
        +k(747)*n(idx_C)  &
        +k(621)*n(idx_Nj)  &
        +k(807)*n(idx_CH)  &
        +k(321)*n(idx_Cj)  &
        +k(858)*n(idx_H)  &
        +k(764)*n(idx_CH2)  &
        +k(1013)  &
        +k(947)*n(idx_O)  &
        +2.d0*k(879)*n(idx_HCO)  &
        +k(350)*n(idx_CHj)  &
        +k(629)*n(idx_N2j)  &
        +k(437)*n(idx_Hj)

    !d[N2_dot]/d[HCO]
    pd(20,14) =  &
        +k(151)*n(idx_N2j)  &
        +k(554)*n(idx_N2Hj)

    !d[NH2_dot]/d[HCO]
    pd(21,14) =  &
        +k(160)*n(idx_NH2j)

    !d[CH3_dot]/d[HCO]
    pd(22,14) =  &
        -k(787)*n(idx_CH3)  &
        +k(764)*n(idx_CH2)  &
        +k(40)*n(idx_CH3j)

    !d[CH4_dot]/d[HCO]
    pd(23,14) =  &
        +k(787)*n(idx_CH3)

    !d[N_dot]/d[HCO]
    pd(24,14) =  &
        -k(897)*n(idx_N)  &
        -k(896)*n(idx_N)  &
        +k(142)*n(idx_Nj)  &
        +k(654)*n(idx_NHj)  &
        -k(895)*n(idx_N)

    !d[NH_dot]/d[HCO]
    pd(25,14) =  &
        +k(895)*n(idx_N)  &
        +k(669)*n(idx_NH2j)

    !d[HE_dot]/d[HCO]
    pd(26,14) =  &
        +k(592)*n(idx_HEj)  &
        +k(590)*n(idx_HEj)

    !d[HNO_dot]/d[HCO]
    pd(27,14) =  &
        -k(880)*n(idx_HNO)  &
        +k(881)*n(idx_NO)

    !d[CO2_dot]/d[HCO]
    pd(29,14) =  &
        +k(946)*n(idx_O)  &
        +k(882)*n(idx_O2)

    !d[O2H_dot]/d[HCO]
    pd(33,14) =  &
        +k(883)*n(idx_O2)  &
        -k(884)*n(idx_O2H)

    !d[OCN_dot]/d[HCO]
    pd(34,14) =  &
        +k(897)*n(idx_N)

    !d[H2CO_DUST_dot]/d[HCO]
    pd(37,14) =  &
        +k(1081)

    !d[HCO+_dot]/d[HCO]
    pd(54,14) =  &
        +k(169)*n(idx_NH3j)  &
        +k(76)*n(idx_Hj)  &
        +k(103)*n(idx_H2Oj)  &
        +k(60)*n(idx_CNj)  &
        +k(40)*n(idx_CH3j)  &
        +k(121)*n(idx_O2j)  &
        +k(142)*n(idx_Nj)  &
        +k(151)*n(idx_N2j)  &
        -k(551)*n(idx_HCOj)  &
        +k(26)*n(idx_CHj)  &
        +k(160)*n(idx_NH2j)  &
        +k(235)  &
        +k(189)*n(idx_Oj)  &
        +k(18)*n(idx_Cj)  &
        +k(199)*n(idx_OHj)  &
        +k(1014)  &
        +k(120)*n(idx_H2COj)  &
        +k(65)*n(idx_COj)  &
        +k(93)*n(idx_H2j)

    !d[H+_dot]/d[HCO]
    pd(55,14) =  &
        -k(437)*n(idx_Hj)  &
        -k(436)*n(idx_Hj)  &
        -k(76)*n(idx_Hj)

    !d[C+_dot]/d[HCO]
    pd(57,14) =  &
        -k(321)*n(idx_Cj)  &
        -k(18)*n(idx_Cj)

    !d[CH2+_dot]/d[HCO]
    pd(58,14) =  &
        +k(350)*n(idx_CHj)  &
        -k(363)*n(idx_CH2j)

    !d[CH+_dot]/d[HCO]
    pd(59,14) =  &
        -k(26)*n(idx_CHj)  &
        +k(592)*n(idx_HEj)  &
        +k(321)*n(idx_Cj)  &
        -k(350)*n(idx_CHj)

    !d[H2CO+_dot]/d[HCO]
    pd(60,14) =  &
        -k(552)*n(idx_H2COj)  &
        +k(556)*n(idx_O2Hj)  &
        +k(554)*n(idx_N2Hj)  &
        +k(551)*n(idx_HCOj)  &
        +k(514)*n(idx_H3j)  &
        +k(733)*n(idx_OHj)  &
        +k(488)*n(idx_H2Oj)  &
        +k(654)*n(idx_NHj)  &
        +k(553)*n(idx_HNOj)  &
        -k(120)*n(idx_H2COj)  &
        +k(669)*n(idx_NH2j)  &
        +k(539)*n(idx_HCNj)

    !d[NH3+_dot]/d[HCO]
    pd(62,14) =  &
        -k(169)*n(idx_NH3j)

    !d[CN+_dot]/d[HCO]
    pd(64,14) =  &
        -k(60)*n(idx_CNj)  &
        -k(419)*n(idx_CNj)

    !d[CO+_dot]/d[HCO]
    pd(65,14) =  &
        -k(65)*n(idx_COj)  &
        +k(590)*n(idx_HEj)  &
        +k(436)*n(idx_Hj)

    !d[N2+_dot]/d[HCO]
    pd(66,14) =  &
        -k(151)*n(idx_N2j)  &
        -k(629)*n(idx_N2j)

    !d[O2+_dot]/d[HCO]
    pd(67,14) =  &
        -k(555)*n(idx_O2j)  &
        -k(121)*n(idx_O2j)

    !d[H2O+_dot]/d[HCO]
    pd(68,14) =  &
        +k(732)*n(idx_OHj)  &
        -k(103)*n(idx_H2Oj)  &
        -k(488)*n(idx_H2Oj)  &
        -k(487)*n(idx_H2Oj)

    !d[NH2+_dot]/d[HCO]
    pd(69,14) =  &
        -k(160)*n(idx_NH2j)  &
        -k(669)*n(idx_NH2j)

    !d[O+_dot]/d[HCO]
    pd(70,14) =  &
        -k(189)*n(idx_Oj)  &
        -k(711)*n(idx_Oj)

    !d[OH+_dot]/d[HCO]
    pd(71,14) =  &
        -k(199)*n(idx_OHj)  &
        +k(711)*n(idx_Oj)  &
        -k(733)*n(idx_OHj)  &
        -k(732)*n(idx_OHj)

    !d[CH3+_dot]/d[HCO]
    pd(72,14) =  &
        +k(363)*n(idx_CH2j)  &
        -k(40)*n(idx_CH3j)  &
        -k(384)*n(idx_CH3j)

    !d[CH4+_dot]/d[HCO]
    pd(73,14) =  &
        +k(384)*n(idx_CH3j)

    !d[N+_dot]/d[HCO]
    pd(74,14) =  &
        -k(142)*n(idx_Nj)  &
        -k(621)*n(idx_Nj)

    !d[HCN+_dot]/d[HCO]
    pd(75,14) =  &
        +k(419)*n(idx_CNj)  &
        -k(539)*n(idx_HCNj)  &
        -k(540)*n(idx_HCNj)

    !d[NH+_dot]/d[HCO]
    pd(76,14) =  &
        +k(621)*n(idx_Nj)  &
        -k(654)*n(idx_NHj)

    !d[H2+_dot]/d[HCO]
    pd(77,14) =  &
        -k(451)*n(idx_H2j)  &
        -k(93)*n(idx_H2j)  &
        +k(437)*n(idx_Hj)

    !d[HE+_dot]/d[HCO]
    pd(78,14) =  &
        -k(591)*n(idx_HEj)  &
        -k(590)*n(idx_HEj)  &
        -k(592)*n(idx_HEj)

    !d[HNO+_dot]/d[HCO]
    pd(79,14) =  &
        -k(553)*n(idx_HNOj)

    !d[H3+_dot]/d[HCO]
    pd(81,14) =  &
        +k(451)*n(idx_H2j)  &
        -k(514)*n(idx_H3j)

    !d[H3CO+_dot]/d[HCO]
    pd(82,14) =  &
        +k(552)*n(idx_H2COj)

    !d[H3O+_dot]/d[HCO]
    pd(83,14) =  &
        +k(487)*n(idx_H2Oj)

    !d[HCNH+_dot]/d[HCO]
    pd(84,14) =  &
        +k(540)*n(idx_HCNj)

    !d[HEH+_dot]/d[HCO]
    pd(86,14) =  &
        +k(591)*n(idx_HEj)

    !d[N2H+_dot]/d[HCO]
    pd(87,14) =  &
        +k(629)*n(idx_N2j)  &
        -k(554)*n(idx_N2Hj)

    !d[O2H+_dot]/d[HCO]
    pd(88,14) =  &
        -k(556)*n(idx_O2Hj)  &
        +k(555)*n(idx_O2j)

    !d[E_dot]/d[MG]
    pd(1,15) =  &
        +k(240)  &
        +k(1018)

    !d[CH_dot]/d[MG]
    pd(2,15) =  &
        +k(27)*n(idx_CHj)

    !d[H2_dot]/d[MG]
    pd(6,15) =  &
        +k(517)*n(idx_H3j)

    !d[C_dot]/d[MG]
    pd(7,15) =  &
        +k(19)*n(idx_Cj)

    !d[H_dot]/d[MG]
    pd(8,15) =  &
        +k(517)*n(idx_H3j)  &
        +k(77)*n(idx_Hj)

    !d[H2O_dot]/d[MG]
    pd(9,15) =  &
        +k(104)*n(idx_H2Oj)

    !d[O2_dot]/d[MG]
    pd(11,15) =  &
        +k(134)*n(idx_O2j)

    !d[H2CO_dot]/d[MG]
    pd(13,15) =  &
        +k(130)*n(idx_H2COj)

    !d[HCO_dot]/d[MG]
    pd(14,15) =  &
        +k(131)*n(idx_HCOj)

    !d[MG_dot]/d[MG]
    pd(15,15) =  &
        -k(130)*n(idx_H2COj)  &
        -k(27)*n(idx_CHj)  &
        -k(517)*n(idx_H3j)  &
        -k(41)*n(idx_CH3j)  &
        -k(131)*n(idx_HCOj)  &
        -k(1018)  &
        -k(19)*n(idx_Cj)  &
        -k(134)*n(idx_O2j)  &
        -k(104)*n(idx_H2Oj)  &
        -k(77)*n(idx_Hj)  &
        -k(240)  &
        -k(143)*n(idx_Nj)  &
        -k(132)*n(idx_N2j)  &
        -k(1120)  &
        -k(133)*n(idx_NOj)  &
        -k(170)*n(idx_NH3j)

    !d[NH3_dot]/d[MG]
    pd(16,15) =  &
        +k(170)*n(idx_NH3j)

    !d[NO_dot]/d[MG]
    pd(17,15) =  &
        +k(133)*n(idx_NOj)

    !d[N2_dot]/d[MG]
    pd(20,15) =  &
        +k(132)*n(idx_N2j)

    !d[CH3_dot]/d[MG]
    pd(22,15) =  &
        +k(41)*n(idx_CH3j)

    !d[N_dot]/d[MG]
    pd(24,15) =  &
        +k(143)*n(idx_Nj)

    !d[MG_DUST_dot]/d[MG]
    pd(51,15) =  &
        +k(1120)

    !d[HCO+_dot]/d[MG]
    pd(54,15) =  &
        -k(131)*n(idx_HCOj)

    !d[H+_dot]/d[MG]
    pd(55,15) =  &
        -k(77)*n(idx_Hj)

    !d[C+_dot]/d[MG]
    pd(57,15) =  &
        -k(19)*n(idx_Cj)

    !d[CH+_dot]/d[MG]
    pd(59,15) =  &
        -k(27)*n(idx_CHj)

    !d[H2CO+_dot]/d[MG]
    pd(60,15) =  &
        -k(130)*n(idx_H2COj)

    !d[MG+_dot]/d[MG]
    pd(61,15) =  &
        +k(170)*n(idx_NH3j)  &
        +k(131)*n(idx_HCOj)  &
        +k(517)*n(idx_H3j)  &
        +k(134)*n(idx_O2j)  &
        +k(104)*n(idx_H2Oj)  &
        +k(77)*n(idx_Hj)  &
        +k(143)*n(idx_Nj)  &
        +k(1018)  &
        +k(132)*n(idx_N2j)  &
        +k(27)*n(idx_CHj)  &
        +k(19)*n(idx_Cj)  &
        +k(133)*n(idx_NOj)  &
        +k(240)  &
        +k(41)*n(idx_CH3j)  &
        +k(130)*n(idx_H2COj)

    !d[NH3+_dot]/d[MG]
    pd(62,15) =  &
        -k(170)*n(idx_NH3j)

    !d[NO+_dot]/d[MG]
    pd(63,15) =  &
        -k(133)*n(idx_NOj)

    !d[N2+_dot]/d[MG]
    pd(66,15) =  &
        -k(132)*n(idx_N2j)

    !d[O2+_dot]/d[MG]
    pd(67,15) =  &
        -k(134)*n(idx_O2j)

    !d[H2O+_dot]/d[MG]
    pd(68,15) =  &
        -k(104)*n(idx_H2Oj)

    !d[CH3+_dot]/d[MG]
    pd(72,15) =  &
        -k(41)*n(idx_CH3j)

    !d[N+_dot]/d[MG]
    pd(74,15) =  &
        -k(143)*n(idx_Nj)

    !d[H3+_dot]/d[MG]
    pd(81,15) =  &
        -k(517)*n(idx_H3j)

    !d[E_dot]/d[NH3]
    pd(1,16) =  &
        +k(246)  &
        +k(1024)

    !d[CH_dot]/d[NH3]
    pd(2,16) =  &
        +k(28)*n(idx_CHj)

    !d[O_dot]/d[NH3]
    pd(3,16) =  &
        -k(954)*n(idx_O)  &
        +k(191)*n(idx_Oj)

    !d[HCN_dot]/d[NH3]
    pd(5,16) =  &
        +k(175)*n(idx_HCNj)  &
        +k(913)*n(idx_CN)

    !d[H2_dot]/d[NH3]
    pd(6,16) =  &
        +k(323)*n(idx_Cj)  &
        +k(1025)  &
        +k(601)*n(idx_HEj)  &
        +k(865)*n(idx_H)  &
        +k(622)*n(idx_Nj)  &
        +k(95)*n(idx_H2j)  &
        +k(247)

    !d[C_dot]/d[NH3]
    pd(7,16) =  &
        +k(20)*n(idx_Cj)

    !d[H_dot]/d[NH3]
    pd(8,16) =  &
        +k(79)*n(idx_Hj)  &
        +k(602)*n(idx_HEj)  &
        +k(1023)  &
        +k(245)  &
        -k(865)*n(idx_H)

    !d[H2O_dot]/d[NH3]
    pd(9,16) =  &
        +k(174)*n(idx_H2Oj)  &
        +k(969)*n(idx_OH)

    !d[OH_dot]/d[NH3]
    pd(10,16) =  &
        -k(969)*n(idx_OH)  &
        +k(200)*n(idx_OHj)  &
        +k(954)*n(idx_O)

    !d[O2_dot]/d[NH3]
    pd(11,16) =  &
        +k(177)*n(idx_O2j)

    !d[H2CO_dot]/d[NH3]
    pd(13,16) =  &
        +k(173)*n(idx_H2COj)

    !d[NH3_dot]/d[NH3]
    pd(16,16) =  &
        -k(623)*n(idx_Nj)  &
        -k(602)*n(idx_HEj)  &
        -k(954)*n(idx_O)  &
        -k(28)*n(idx_CHj)  &
        -k(245)  &
        -k(174)*n(idx_H2Oj)  &
        -k(79)*n(idx_Hj)  &
        -k(1087)  &
        -k(172)*n(idx_COj)  &
        -k(20)*n(idx_Cj)  &
        -k(917)*n(idx_NH)  &
        -k(95)*n(idx_H2j)  &
        -k(200)*n(idx_OHj)  &
        -k(173)*n(idx_H2COj)  &
        -k(247)  &
        -k(790)*n(idx_CH3)  &
        -k(323)*n(idx_Cj)  &
        -k(1024)  &
        -k(44)*n(idx_CH4j)  &
        -k(688)*n(idx_HCNj)  &
        -k(687)*n(idx_COj)  &
        -k(913)*n(idx_CN)  &
        -k(1023)  &
        -k(177)*n(idx_O2j)  &
        -k(1025)  &
        -k(157)*n(idx_NHj)  &
        -k(128)*n(idx_HEj)  &
        -k(191)*n(idx_Oj)  &
        -k(622)*n(idx_Nj)  &
        -k(969)*n(idx_OH)  &
        -k(865)*n(idx_H)  &
        -k(161)*n(idx_NH2j)  &
        -k(176)*n(idx_N2j)  &
        -k(145)*n(idx_Nj)  &
        -k(246)  &
        -k(601)*n(idx_HEj)  &
        -k(175)*n(idx_HCNj)

    !d[CN_dot]/d[NH3]
    pd(18,16) =  &
        -k(913)*n(idx_CN)

    !d[CO_dot]/d[NH3]
    pd(19,16) =  &
        +k(172)*n(idx_COj)

    !d[N2_dot]/d[NH3]
    pd(20,16) =  &
        +k(176)*n(idx_N2j)

    !d[NH2_dot]/d[NH3]
    pd(21,16) =  &
        +2.d0*k(917)*n(idx_NH)  &
        +k(245)  &
        +k(969)*n(idx_OH)  &
        +k(954)*n(idx_O)  &
        +k(688)*n(idx_HCNj)  &
        +k(687)*n(idx_COj)  &
        +k(161)*n(idx_NH2j)  &
        +k(865)*n(idx_H)  &
        +k(1023)  &
        +k(790)*n(idx_CH3)  &
        +k(913)*n(idx_CN)

    !d[CH3_dot]/d[NH3]
    pd(22,16) =  &
        -k(790)*n(idx_CH3)

    !d[CH4_dot]/d[NH3]
    pd(23,16) =  &
        +k(44)*n(idx_CH4j)  &
        +k(790)*n(idx_CH3)

    !d[N_dot]/d[NH3]
    pd(24,16) =  &
        +k(145)*n(idx_Nj)

    !d[NH_dot]/d[NH3]
    pd(25,16) =  &
        +k(247)  &
        -k(917)*n(idx_NH)  &
        +k(1025)  &
        +k(623)*n(idx_Nj)  &
        +k(157)*n(idx_NHj)

    !d[HE_dot]/d[NH3]
    pd(26,16) =  &
        +k(602)*n(idx_HEj)  &
        +k(601)*n(idx_HEj)  &
        +k(128)*n(idx_HEj)

    !d[NH3_DUST_dot]/d[NH3]
    pd(45,16) =  &
        +k(1087)

    !d[HCO+_dot]/d[NH3]
    pd(54,16) =  &
        +k(687)*n(idx_COj)

    !d[H+_dot]/d[NH3]
    pd(55,16) =  &
        -k(79)*n(idx_Hj)

    !d[C+_dot]/d[NH3]
    pd(57,16) =  &
        -k(20)*n(idx_Cj)  &
        -k(323)*n(idx_Cj)

    !d[CH+_dot]/d[NH3]
    pd(59,16) =  &
        -k(28)*n(idx_CHj)

    !d[H2CO+_dot]/d[NH3]
    pd(60,16) =  &
        -k(173)*n(idx_H2COj)

    !d[NH3+_dot]/d[NH3]
    pd(62,16) =  &
        +k(246)  &
        +k(44)*n(idx_CH4j)  &
        +k(191)*n(idx_Oj)  &
        +k(28)*n(idx_CHj)  &
        +k(176)*n(idx_N2j)  &
        +k(157)*n(idx_NHj)  &
        +k(200)*n(idx_OHj)  &
        +k(1024)  &
        +k(175)*n(idx_HCNj)  &
        +k(172)*n(idx_COj)  &
        +k(145)*n(idx_Nj)  &
        +k(177)*n(idx_O2j)  &
        +k(161)*n(idx_NH2j)  &
        +k(128)*n(idx_HEj)  &
        +k(79)*n(idx_Hj)  &
        +k(20)*n(idx_Cj)  &
        +k(174)*n(idx_H2Oj)  &
        +k(95)*n(idx_H2j)  &
        +k(173)*n(idx_H2COj)

    !d[CO+_dot]/d[NH3]
    pd(65,16) =  &
        -k(687)*n(idx_COj)  &
        -k(172)*n(idx_COj)

    !d[N2+_dot]/d[NH3]
    pd(66,16) =  &
        -k(176)*n(idx_N2j)

    !d[O2+_dot]/d[NH3]
    pd(67,16) =  &
        -k(177)*n(idx_O2j)

    !d[H2O+_dot]/d[NH3]
    pd(68,16) =  &
        -k(174)*n(idx_H2Oj)

    !d[NH2+_dot]/d[NH3]
    pd(69,16) =  &
        +k(602)*n(idx_HEj)  &
        -k(161)*n(idx_NH2j)  &
        +k(623)*n(idx_Nj)

    !d[O+_dot]/d[NH3]
    pd(70,16) =  &
        -k(191)*n(idx_Oj)

    !d[OH+_dot]/d[NH3]
    pd(71,16) =  &
        -k(200)*n(idx_OHj)

    !d[CH4+_dot]/d[NH3]
    pd(73,16) =  &
        -k(44)*n(idx_CH4j)

    !d[N+_dot]/d[NH3]
    pd(74,16) =  &
        -k(623)*n(idx_Nj)  &
        -k(145)*n(idx_Nj)  &
        -k(622)*n(idx_Nj)

    !d[HCN+_dot]/d[NH3]
    pd(75,16) =  &
        +k(323)*n(idx_Cj)  &
        -k(175)*n(idx_HCNj)  &
        -k(688)*n(idx_HCNj)

    !d[NH+_dot]/d[NH3]
    pd(76,16) =  &
        +k(601)*n(idx_HEj)  &
        -k(157)*n(idx_NHj)

    !d[H2+_dot]/d[NH3]
    pd(77,16) =  &
        -k(95)*n(idx_H2j)

    !d[HE+_dot]/d[NH3]
    pd(78,16) =  &
        -k(128)*n(idx_HEj)  &
        -k(601)*n(idx_HEj)  &
        -k(602)*n(idx_HEj)

    !d[HCNH+_dot]/d[NH3]
    pd(84,16) =  &
        +k(688)*n(idx_HCNj)

    !d[N2H+_dot]/d[NH3]
    pd(87,16) =  &
        +k(622)*n(idx_Nj)

    !d[E_dot]/d[NO]
    pd(1,17) =  &
        +k(251)  &
        +k(1029)

    !d[CH_dot]/d[NO]
    pd(2,17) =  &
        -k(812)*n(idx_CH)  &
        -k(813)*n(idx_CH)  &
        -k(814)*n(idx_CH)  &
        +k(29)*n(idx_CHj)

    !d[O_dot]/d[NO]
    pd(3,17) =  &
        +k(903)*n(idx_N)  &
        +k(625)*n(idx_Nj)  &
        +k(754)*n(idx_C)  &
        +k(868)*n(idx_H)  &
        +k(605)*n(idx_HEj)  &
        -k(956)*n(idx_O)  &
        +k(932)*n(idx_O2)  &
        +k(1030)  &
        +k(252)  &
        +k(659)*n(idx_NHj)  &
        +k(736)*n(idx_OHj)  &
        +k(922)*n(idx_NH)  &
        +k(812)*n(idx_CH)

    !d[HCN_dot]/d[NO]
    pd(5,17) =  &
        +k(769)*n(idx_CH2)  &
        +k(792)*n(idx_CH3)  &
        +k(116)*n(idx_HCNj)  &
        +k(812)*n(idx_CH)

    !d[H2_dot]/d[NO]
    pd(6,17) =  &
        +k(97)*n(idx_H2j)  &
        +k(522)*n(idx_H3j)

    !d[C_dot]/d[NO]
    pd(7,17) =  &
        -k(754)*n(idx_C)  &
        -k(755)*n(idx_C)  &
        +k(21)*n(idx_Cj)

    !d[H_dot]/d[NO]
    pd(8,17) =  &
        -k(868)*n(idx_H)  &
        +k(970)*n(idx_OH)  &
        +k(814)*n(idx_CH)  &
        +k(81)*n(idx_Hj)  &
        +k(910)*n(idx_NH2)  &
        +k(456)*n(idx_H2j)  &
        -k(869)*n(idx_H)  &
        +k(922)*n(idx_NH)  &
        +k(770)*n(idx_CH2)

    !d[H2O_dot]/d[NO]
    pd(9,17) =  &
        +k(105)*n(idx_H2Oj)  &
        +k(792)*n(idx_CH3)  &
        +k(909)*n(idx_NH2)

    !d[OH_dot]/d[NO]
    pd(10,17) =  &
        +k(923)*n(idx_NH)  &
        +k(201)*n(idx_OHj)  &
        -k(970)*n(idx_OH)  &
        +k(869)*n(idx_H)  &
        +k(910)*n(idx_NH2)  &
        +k(769)*n(idx_CH2)

    !d[O2_dot]/d[NO]
    pd(11,17) =  &
        +k(702)*n(idx_O2Hj)  &
        +k(184)*n(idx_O2j)  &
        +2.d0*k(931)*n(idx_NO)  &
        -k(932)*n(idx_O2)  &
        +k(956)*n(idx_O)

    !d[CH2_dot]/d[NO]
    pd(12,17) =  &
        +k(30)*n(idx_CH2j)  &
        -k(770)*n(idx_CH2)  &
        -k(768)*n(idx_CH2)  &
        -k(769)*n(idx_CH2)

    !d[H2CO_dot]/d[NO]
    pd(13,17) =  &
        +k(768)*n(idx_CH2)  &
        +k(182)*n(idx_H2COj)

    !d[HCO_dot]/d[NO]
    pd(14,17) =  &
        +k(813)*n(idx_CH)  &
        -k(881)*n(idx_HCO)

    !d[NH3_dot]/d[NO]
    pd(16,17) =  &
        +k(171)*n(idx_NH3j)

    !d[NO_dot]/d[NO]
    pd(17,17) =  &
        -k(868)*n(idx_H)  &
        -k(754)*n(idx_C)  &
        -k(755)*n(idx_C)  &
        -k(881)*n(idx_HCO)  &
        -4.d0*k(931)*n(idx_NO)  &
        -k(922)*n(idx_NH)  &
        -k(81)*n(idx_Hj)  &
        -k(158)*n(idx_NHj)  &
        -k(923)*n(idx_NH)  &
        -k(659)*n(idx_NHj)  &
        -k(933)*n(idx_OCN)  &
        -k(105)*n(idx_H2Oj)  &
        -k(522)*n(idx_H3j)  &
        -k(116)*n(idx_HCNj)  &
        -k(956)*n(idx_O)  &
        -k(768)*n(idx_CH2)  &
        -k(770)*n(idx_CH2)  &
        -k(970)*n(idx_OH)  &
        -k(814)*n(idx_CH)  &
        -k(909)*n(idx_NH2)  &
        -k(605)*n(idx_HEj)  &
        -k(812)*n(idx_CH)  &
        -k(792)*n(idx_CH3)  &
        -k(869)*n(idx_H)  &
        -k(1075)  &
        -k(251)  &
        -k(769)*n(idx_CH2)  &
        -k(42)*n(idx_CH3j)  &
        -k(66)*n(idx_COj)  &
        -k(184)*n(idx_O2j)  &
        -k(183)*n(idx_HNOj)  &
        -k(1030)  &
        -k(201)*n(idx_OHj)  &
        -k(29)*n(idx_CHj)  &
        -k(456)*n(idx_H2j)  &
        -k(21)*n(idx_Cj)  &
        -k(252)  &
        -k(828)*n(idx_CN)  &
        -k(97)*n(idx_H2j)  &
        -k(61)*n(idx_CNj)  &
        -k(932)*n(idx_O2)  &
        -k(152)*n(idx_N2j)  &
        -k(813)*n(idx_CH)  &
        -k(147)*n(idx_Nj)  &
        -k(910)*n(idx_NH2)  &
        -k(625)*n(idx_Nj)  &
        -k(604)*n(idx_HEj)  &
        -k(171)*n(idx_NH3j)  &
        -k(1029)  &
        -k(182)*n(idx_H2COj)  &
        -k(736)*n(idx_OHj)  &
        -k(162)*n(idx_NH2j)  &
        -k(30)*n(idx_CH2j)  &
        -k(702)*n(idx_O2Hj)  &
        -k(829)*n(idx_CN)  &
        -k(903)*n(idx_N)

    !d[CN_dot]/d[NO]
    pd(18,17) =  &
        -k(829)*n(idx_CN)  &
        +k(754)*n(idx_C)  &
        -k(828)*n(idx_CN)  &
        +k(61)*n(idx_CNj)

    !d[CO_dot]/d[NO]
    pd(19,17) =  &
        +k(755)*n(idx_C)  &
        +k(66)*n(idx_COj)  &
        +k(828)*n(idx_CN)  &
        +k(881)*n(idx_HCO)

    !d[N2_dot]/d[NO]
    pd(20,17) =  &
        +k(903)*n(idx_N)  &
        +k(923)*n(idx_NH)  &
        +k(933)*n(idx_OCN)  &
        +k(910)*n(idx_NH2)  &
        +2.d0*k(931)*n(idx_NO)  &
        +k(152)*n(idx_N2j)  &
        +k(828)*n(idx_CN)  &
        +k(909)*n(idx_NH2)  &
        +k(922)*n(idx_NH)

    !d[NH2_dot]/d[NO]
    pd(21,17) =  &
        -k(909)*n(idx_NH2)  &
        -k(910)*n(idx_NH2)  &
        +k(162)*n(idx_NH2j)

    !d[CH3_dot]/d[NO]
    pd(22,17) =  &
        +k(42)*n(idx_CH3j)  &
        -k(792)*n(idx_CH3)

    !d[N_dot]/d[NO]
    pd(24,17) =  &
        +k(755)*n(idx_C)  &
        +k(147)*n(idx_Nj)  &
        +k(956)*n(idx_O)  &
        +k(813)*n(idx_CH)  &
        +k(869)*n(idx_H)  &
        +k(768)*n(idx_CH2)  &
        +k(604)*n(idx_HEj)  &
        +k(252)  &
        +k(1030)  &
        -k(903)*n(idx_N)  &
        +k(829)*n(idx_CN)

    !d[NH_dot]/d[NO]
    pd(25,17) =  &
        -k(922)*n(idx_NH)  &
        +k(158)*n(idx_NHj)  &
        +k(868)*n(idx_H)  &
        -k(923)*n(idx_NH)

    !d[HE_dot]/d[NO]
    pd(26,17) =  &
        +k(605)*n(idx_HEj)  &
        +k(604)*n(idx_HEj)

    !d[HNO_dot]/d[NO]
    pd(27,17) =  &
        +k(183)*n(idx_HNOj)  &
        +k(881)*n(idx_HCO)

    !d[CO2_dot]/d[NO]
    pd(29,17) =  &
        +k(933)*n(idx_OCN)

    !d[HNCO_dot]/d[NO]
    pd(31,17) =  &
        +k(770)*n(idx_CH2)

    !d[NO2_dot]/d[NO]
    pd(32,17) =  &
        +k(932)*n(idx_O2)  &
        +k(970)*n(idx_OH)

    !d[OCN_dot]/d[NO]
    pd(34,17) =  &
        -k(933)*n(idx_OCN)  &
        +k(814)*n(idx_CH)  &
        +k(829)*n(idx_CN)

    !d[NO_DUST_dot]/d[NO]
    pd(41,17) =  &
        +k(1075)

    !d[H+_dot]/d[NO]
    pd(55,17) =  &
        -k(81)*n(idx_Hj)

    !d[C+_dot]/d[NO]
    pd(57,17) =  &
        -k(21)*n(idx_Cj)

    !d[CH2+_dot]/d[NO]
    pd(58,17) =  &
        -k(30)*n(idx_CH2j)

    !d[CH+_dot]/d[NO]
    pd(59,17) =  &
        -k(29)*n(idx_CHj)

    !d[H2CO+_dot]/d[NO]
    pd(60,17) =  &
        -k(182)*n(idx_H2COj)

    !d[NH3+_dot]/d[NO]
    pd(62,17) =  &
        -k(171)*n(idx_NH3j)

    !d[NO+_dot]/d[NO]
    pd(63,17) =  &
        +k(183)*n(idx_HNOj)  &
        +k(66)*n(idx_COj)  &
        +k(201)*n(idx_OHj)  &
        +k(171)*n(idx_NH3j)  &
        +k(21)*n(idx_Cj)  &
        +k(97)*n(idx_H2j)  &
        +k(30)*n(idx_CH2j)  &
        +k(251)  &
        +k(105)*n(idx_H2Oj)  &
        +k(1029)  &
        +k(158)*n(idx_NHj)  &
        +k(182)*n(idx_H2COj)  &
        +k(162)*n(idx_NH2j)  &
        +k(147)*n(idx_Nj)  &
        +k(81)*n(idx_Hj)  &
        +k(42)*n(idx_CH3j)  &
        +k(152)*n(idx_N2j)  &
        +k(116)*n(idx_HCNj)  &
        +k(184)*n(idx_O2j)  &
        +k(29)*n(idx_CHj)  &
        +k(61)*n(idx_CNj)

    !d[CN+_dot]/d[NO]
    pd(64,17) =  &
        -k(61)*n(idx_CNj)

    !d[CO+_dot]/d[NO]
    pd(65,17) =  &
        -k(66)*n(idx_COj)

    !d[N2+_dot]/d[NO]
    pd(66,17) =  &
        +k(625)*n(idx_Nj)  &
        -k(152)*n(idx_N2j)

    !d[O2+_dot]/d[NO]
    pd(67,17) =  &
        -k(184)*n(idx_O2j)

    !d[H2O+_dot]/d[NO]
    pd(68,17) =  &
        -k(105)*n(idx_H2Oj)

    !d[NH2+_dot]/d[NO]
    pd(69,17) =  &
        -k(162)*n(idx_NH2j)

    !d[O+_dot]/d[NO]
    pd(70,17) =  &
        +k(604)*n(idx_HEj)

    !d[OH+_dot]/d[NO]
    pd(71,17) =  &
        -k(201)*n(idx_OHj)  &
        -k(736)*n(idx_OHj)

    !d[CH3+_dot]/d[NO]
    pd(72,17) =  &
        -k(42)*n(idx_CH3j)

    !d[N+_dot]/d[NO]
    pd(74,17) =  &
        -k(625)*n(idx_Nj)  &
        +k(605)*n(idx_HEj)  &
        -k(147)*n(idx_Nj)

    !d[HCN+_dot]/d[NO]
    pd(75,17) =  &
        -k(116)*n(idx_HCNj)

    !d[NH+_dot]/d[NO]
    pd(76,17) =  &
        -k(659)*n(idx_NHj)  &
        -k(158)*n(idx_NHj)

    !d[H2+_dot]/d[NO]
    pd(77,17) =  &
        -k(456)*n(idx_H2j)  &
        -k(97)*n(idx_H2j)

    !d[HE+_dot]/d[NO]
    pd(78,17) =  &
        -k(605)*n(idx_HEj)  &
        -k(604)*n(idx_HEj)

    !d[HNO+_dot]/d[NO]
    pd(79,17) =  &
        +k(702)*n(idx_O2Hj)  &
        +k(736)*n(idx_OHj)  &
        +k(522)*n(idx_H3j)  &
        -k(183)*n(idx_HNOj)  &
        +k(456)*n(idx_H2j)

    !d[H3+_dot]/d[NO]
    pd(81,17) =  &
        -k(522)*n(idx_H3j)

    !d[N2H+_dot]/d[NO]
    pd(87,17) =  &
        +k(659)*n(idx_NHj)

    !d[O2H+_dot]/d[NO]
    pd(88,17) =  &
        -k(702)*n(idx_O2Hj)

    !d[CH_dot]/d[CN]
    pd(2,18) =  &
        +k(762)*n(idx_CH2)

    !d[O_dot]/d[CN]
    pd(3,18) =  &
        +k(726)*n(idx_OHj)  &
        -k(938)*n(idx_O)  &
        -k(937)*n(idx_O)  &
        +k(831)*n(idx_O2)  &
        +k(961)*n(idx_OH)

    !d[HCN_dot]/d[CN]
    pd(5,18) =  &
        +k(826)*n(idx_HNO)  &
        +k(824)*n(idx_H2CO)  &
        +k(802)*n(idx_CH4)  &
        +k(840)*n(idx_H2)  &
        +k(961)*n(idx_OH)  &
        +k(915)*n(idx_NH)  &
        +k(762)*n(idx_CH2)  &
        +k(825)*n(idx_HCO)  &
        +k(784)*n(idx_CH3)  &
        +k(913)*n(idx_NH3)

    !d[H2_dot]/d[CN]
    pd(6,18) =  &
        +k(88)*n(idx_H2j)  &
        -k(840)*n(idx_H2)  &
        +k(507)*n(idx_H3j)

    !d[C_dot]/d[CN]
    pd(7,18) =  &
        +k(996)  &
        +k(938)*n(idx_O)  &
        +k(706)*n(idx_Oj)  &
        +k(892)*n(idx_N)  &
        +k(574)*n(idx_HEj)  &
        +k(226)

    !d[H_dot]/d[CN]
    pd(8,18) =  &
        +k(840)*n(idx_H2)  &
        +k(445)*n(idx_H2j)  &
        +k(962)*n(idx_OH)

    !d[OH_dot]/d[CN]
    pd(10,18) =  &
        -k(961)*n(idx_OH)  &
        -k(962)*n(idx_OH)

    !d[O2_dot]/d[CN]
    pd(11,18) =  &
        -k(831)*n(idx_O2)  &
        +k(422)*n(idx_O2Hj)  &
        -k(830)*n(idx_O2)

    !d[CH2_dot]/d[CN]
    pd(12,18) =  &
        +k(784)*n(idx_CH3)  &
        -k(762)*n(idx_CH2)

    !d[H2CO_dot]/d[CN]
    pd(13,18) =  &
        -k(824)*n(idx_H2CO)

    !d[HCO_dot]/d[CN]
    pd(14,18) =  &
        +k(824)*n(idx_H2CO)  &
        -k(825)*n(idx_HCO)

    !d[NH3_dot]/d[CN]
    pd(16,18) =  &
        -k(913)*n(idx_NH3)

    !d[NO_dot]/d[CN]
    pd(17,18) =  &
        +k(826)*n(idx_HNO)  &
        +k(421)*n(idx_HNOj)  &
        +k(938)*n(idx_O)  &
        -k(829)*n(idx_NO)  &
        +k(830)*n(idx_O2)  &
        -k(828)*n(idx_NO)  &
        +k(827)*n(idx_NO2)

    !d[CN_dot]/d[CN]
    pd(18,18) =  &
        -k(574)*n(idx_HEj)  &
        -k(913)*n(idx_NH3)  &
        -k(445)*n(idx_H2j)  &
        -k(726)*n(idx_OHj)  &
        -k(762)*n(idx_CH2)  &
        -k(962)*n(idx_OH)  &
        -k(802)*n(idx_CH4)  &
        -k(915)*n(idx_NH)  &
        -k(824)*n(idx_H2CO)  &
        -k(996)  &
        -k(226)  &
        -k(784)*n(idx_CH3)  &
        -k(706)*n(idx_Oj)  &
        -k(831)*n(idx_O2)  &
        -k(88)*n(idx_H2j)  &
        -k(507)*n(idx_H3j)  &
        -k(938)*n(idx_O)  &
        -k(828)*n(idx_NO)  &
        -k(892)*n(idx_N)  &
        -k(575)*n(idx_HEj)  &
        -k(421)*n(idx_HNOj)  &
        -k(830)*n(idx_O2)  &
        -k(137)*n(idx_Nj)  &
        -k(961)*n(idx_OH)  &
        -k(422)*n(idx_O2Hj)  &
        -k(1083)  &
        -k(829)*n(idx_NO)  &
        -k(937)*n(idx_O)  &
        -k(63)*n(idx_N2j)  &
        -k(642)*n(idx_NHj)  &
        -k(827)*n(idx_NO2)  &
        -k(826)*n(idx_HNO)  &
        -k(825)*n(idx_HCO)  &
        -k(840)*n(idx_H2)

    !d[CO_dot]/d[CN]
    pd(19,18) =  &
        +k(830)*n(idx_O2)  &
        +k(937)*n(idx_O)  &
        +k(828)*n(idx_NO)  &
        +k(825)*n(idx_HCO)

    !d[N2_dot]/d[CN]
    pd(20,18) =  &
        +k(63)*n(idx_N2j)  &
        +k(828)*n(idx_NO)  &
        +k(892)*n(idx_N)

    !d[NH2_dot]/d[CN]
    pd(21,18) =  &
        +k(913)*n(idx_NH3)

    !d[CH3_dot]/d[CN]
    pd(22,18) =  &
        +k(802)*n(idx_CH4)  &
        -k(784)*n(idx_CH3)

    !d[CH4_dot]/d[CN]
    pd(23,18) =  &
        -k(802)*n(idx_CH4)

    !d[N_dot]/d[CN]
    pd(24,18) =  &
        +k(937)*n(idx_O)  &
        +k(996)  &
        +k(915)*n(idx_NH)  &
        +k(642)*n(idx_NHj)  &
        +k(575)*n(idx_HEj)  &
        -k(892)*n(idx_N)  &
        +k(137)*n(idx_Nj)  &
        +k(829)*n(idx_NO)  &
        +k(226)

    !d[NH_dot]/d[CN]
    pd(25,18) =  &
        -k(915)*n(idx_NH)

    !d[HE_dot]/d[CN]
    pd(26,18) =  &
        +k(575)*n(idx_HEj)  &
        +k(574)*n(idx_HEj)

    !d[HNO_dot]/d[CN]
    pd(27,18) =  &
        -k(826)*n(idx_HNO)

    !d[NO2_dot]/d[CN]
    pd(32,18) =  &
        -k(827)*n(idx_NO2)

    !d[OCN_dot]/d[CN]
    pd(34,18) =  &
        +k(829)*n(idx_NO)  &
        +k(831)*n(idx_O2)  &
        +k(962)*n(idx_OH)  &
        +k(827)*n(idx_NO2)

    !d[HCN_DUST_dot]/d[CN]
    pd(44,18) =  &
        +k(1083)

    !d[C+_dot]/d[CN]
    pd(57,18) =  &
        +k(575)*n(idx_HEj)

    !d[NO+_dot]/d[CN]
    pd(63,18) =  &
        +k(706)*n(idx_Oj)

    !d[CN+_dot]/d[CN]
    pd(64,18) =  &
        +k(88)*n(idx_H2j)  &
        +k(63)*n(idx_N2j)  &
        +k(137)*n(idx_Nj)

    !d[N2+_dot]/d[CN]
    pd(66,18) =  &
        -k(63)*n(idx_N2j)

    !d[O+_dot]/d[CN]
    pd(70,18) =  &
        -k(706)*n(idx_Oj)

    !d[OH+_dot]/d[CN]
    pd(71,18) =  &
        -k(726)*n(idx_OHj)

    !d[N+_dot]/d[CN]
    pd(74,18) =  &
        -k(137)*n(idx_Nj)  &
        +k(574)*n(idx_HEj)

    !d[HCN+_dot]/d[CN]
    pd(75,18) =  &
        +k(726)*n(idx_OHj)  &
        +k(421)*n(idx_HNOj)  &
        +k(507)*n(idx_H3j)  &
        +k(422)*n(idx_O2Hj)  &
        +k(642)*n(idx_NHj)  &
        +k(445)*n(idx_H2j)

    !d[NH+_dot]/d[CN]
    pd(76,18) =  &
        -k(642)*n(idx_NHj)

    !d[H2+_dot]/d[CN]
    pd(77,18) =  &
        -k(88)*n(idx_H2j)  &
        -k(445)*n(idx_H2j)

    !d[HE+_dot]/d[CN]
    pd(78,18) =  &
        -k(575)*n(idx_HEj)  &
        -k(574)*n(idx_HEj)

    !d[HNO+_dot]/d[CN]
    pd(79,18) =  &
        -k(421)*n(idx_HNOj)

    !d[H3+_dot]/d[CN]
    pd(81,18) =  &
        -k(507)*n(idx_H3j)

    !d[O2H+_dot]/d[CN]
    pd(88,18) =  &
        -k(422)*n(idx_O2Hj)

    !d[E_dot]/d[CO]
    pd(1,19) =  &
        +k(207)

    !d[O_dot]/d[CO]
    pd(3,19) =  &
        +k(580)*n(idx_HEj)  &
        +k(834)*n(idx_O2)  &
        +k(999)  &
        +k(728)*n(idx_OHj)  &
        +k(228)  &
        +k(186)*n(idx_Oj)

    !d[H2_dot]/d[CO]
    pd(6,19) =  &
        +k(89)*n(idx_H2j)  &
        +k(510)*n(idx_H3j)  &
        +k(509)*n(idx_H3j)

    !d[C_dot]/d[CO]
    pd(7,19) =  &
        +k(853)*n(idx_H)  &
        +k(618)*n(idx_Nj)  &
        +k(228)  &
        +k(999)

    !d[H_dot]/d[CO]
    pd(8,19) =  &
        +k(447)*n(idx_H2j)  &
        +k(963)*n(idx_OH)  &
        -k(853)*n(idx_H)

    !d[OH_dot]/d[CO]
    pd(10,19) =  &
        +k(835)*n(idx_O2H)  &
        -k(963)*n(idx_OH)  &
        +k(483)*n(idx_H2Oj)  &
        +k(853)*n(idx_H)

    !d[O2_dot]/d[CO]
    pd(11,19) =  &
        +k(427)*n(idx_O2Hj)  &
        -k(834)*n(idx_O2)

    !d[NO_dot]/d[CO]
    pd(17,19) =  &
        +k(833)*n(idx_NO2)  &
        +k(425)*n(idx_HNOj)

    !d[CN_dot]/d[CO]
    pd(18,19) =  &
        +k(536)*n(idx_HCNj)  &
        +k(57)*n(idx_CNj)

    !d[CO_dot]/d[CO]
    pd(19,19) =  &
        -k(832)*n(idx_HNO)  &
        -k(89)*n(idx_H2j)  &
        -k(833)*n(idx_NO2)  &
        -k(447)*n(idx_H2j)  &
        -k(68)*n(idx_N2j)  &
        -k(509)*n(idx_H3j)  &
        -k(424)*n(idx_HCO2j)  &
        -k(999)  &
        -k(728)*n(idx_OHj)  &
        -k(536)*n(idx_HCNj)  &
        -k(483)*n(idx_H2Oj)  &
        -k(427)*n(idx_O2Hj)  &
        -k(186)*n(idx_Oj)  &
        -k(580)*n(idx_HEj)  &
        -k(834)*n(idx_O2)  &
        -k(1124)  &
        -k(228)  &
        -k(390)*n(idx_CH4j)  &
        -k(1071)  &
        -k(510)*n(idx_H3j)  &
        -k(835)*n(idx_O2H)  &
        -k(426)*n(idx_N2Hj)  &
        -k(646)*n(idx_NHj)  &
        -k(1068)  &
        -k(138)*n(idx_Nj)  &
        -k(963)*n(idx_OH)  &
        -k(425)*n(idx_HNOj)  &
        -k(853)*n(idx_H)  &
        -k(207)  &
        -k(618)*n(idx_Nj)  &
        -k(57)*n(idx_CNj)

    !d[N2_dot]/d[CO]
    pd(20,19) =  &
        +k(68)*n(idx_N2j)  &
        +k(426)*n(idx_N2Hj)

    !d[CH3_dot]/d[CO]
    pd(22,19) =  &
        +k(390)*n(idx_CH4j)

    !d[N_dot]/d[CO]
    pd(24,19) =  &
        +k(646)*n(idx_NHj)  &
        +k(138)*n(idx_Nj)

    !d[NH_dot]/d[CO]
    pd(25,19) =  &
        +k(832)*n(idx_HNO)

    !d[HE_dot]/d[CO]
    pd(26,19) =  &
        +k(580)*n(idx_HEj)

    !d[HNO_dot]/d[CO]
    pd(27,19) =  &
        -k(832)*n(idx_HNO)

    !d[CO2_dot]/d[CO]
    pd(29,19) =  &
        +k(835)*n(idx_O2H)  &
        +k(834)*n(idx_O2)  &
        +k(832)*n(idx_HNO)  &
        +k(424)*n(idx_HCO2j)  &
        +k(833)*n(idx_NO2)  &
        +k(963)*n(idx_OH)

    !d[NO2_dot]/d[CO]
    pd(32,19) =  &
        -k(833)*n(idx_NO2)

    !d[O2H_dot]/d[CO]
    pd(33,19) =  &
        -k(835)*n(idx_O2H)

    !d[CH3OH_DUST_dot]/d[CO]
    pd(35,19) =  &
        +k(1124)

    !d[H2CO_DUST_dot]/d[CO]
    pd(37,19) =  &
        +k(1068)

    !d[CO_DUST_dot]/d[CO]
    pd(39,19) =  &
        +k(1071)

    !d[HCO+_dot]/d[CO]
    pd(54,19) =  &
        +k(646)*n(idx_NHj)  &
        +k(425)*n(idx_HNOj)  &
        +k(390)*n(idx_CH4j)  &
        +k(483)*n(idx_H2Oj)  &
        +k(427)*n(idx_O2Hj)  &
        +k(424)*n(idx_HCO2j)  &
        +k(728)*n(idx_OHj)  &
        +k(536)*n(idx_HCNj)  &
        +k(426)*n(idx_N2Hj)  &
        +k(447)*n(idx_H2j)  &
        +k(509)*n(idx_H3j)

    !d[HOC+_dot]/d[CO]
    pd(56,19) =  &
        +k(510)*n(idx_H3j)

    !d[C+_dot]/d[CO]
    pd(57,19) =  &
        +k(580)*n(idx_HEj)

    !d[NO+_dot]/d[CO]
    pd(63,19) =  &
        +k(618)*n(idx_Nj)

    !d[CN+_dot]/d[CO]
    pd(64,19) =  &
        -k(57)*n(idx_CNj)

    !d[CO+_dot]/d[CO]
    pd(65,19) =  &
        +k(207)  &
        +k(68)*n(idx_N2j)  &
        +k(57)*n(idx_CNj)  &
        +k(186)*n(idx_Oj)  &
        +k(89)*n(idx_H2j)  &
        +k(138)*n(idx_Nj)

    !d[N2+_dot]/d[CO]
    pd(66,19) =  &
        -k(68)*n(idx_N2j)

    !d[H2O+_dot]/d[CO]
    pd(68,19) =  &
        -k(483)*n(idx_H2Oj)

    !d[O+_dot]/d[CO]
    pd(70,19) =  &
        -k(186)*n(idx_Oj)

    !d[OH+_dot]/d[CO]
    pd(71,19) =  &
        -k(728)*n(idx_OHj)

    !d[CH4+_dot]/d[CO]
    pd(73,19) =  &
        -k(390)*n(idx_CH4j)

    !d[N+_dot]/d[CO]
    pd(74,19) =  &
        -k(618)*n(idx_Nj)  &
        -k(138)*n(idx_Nj)

    !d[HCN+_dot]/d[CO]
    pd(75,19) =  &
        -k(536)*n(idx_HCNj)

    !d[NH+_dot]/d[CO]
    pd(76,19) =  &
        -k(646)*n(idx_NHj)

    !d[H2+_dot]/d[CO]
    pd(77,19) =  &
        -k(89)*n(idx_H2j)  &
        -k(447)*n(idx_H2j)

    !d[HE+_dot]/d[CO]
    pd(78,19) =  &
        -k(580)*n(idx_HEj)

    !d[HNO+_dot]/d[CO]
    pd(79,19) =  &
        -k(425)*n(idx_HNOj)

    !d[H3+_dot]/d[CO]
    pd(81,19) =  &
        -k(509)*n(idx_H3j)  &
        -k(510)*n(idx_H3j)

    !d[HCO2+_dot]/d[CO]
    pd(85,19) =  &
        -k(424)*n(idx_HCO2j)

    !d[N2H+_dot]/d[CO]
    pd(87,19) =  &
        -k(426)*n(idx_N2Hj)

    !d[O2H+_dot]/d[CO]
    pd(88,19) =  &
        -k(427)*n(idx_O2Hj)

    !d[CH_dot]/d[N2]
    pd(2,20) =  &
        -k(809)*n(idx_CH)

    !d[O_dot]/d[N2]
    pd(3,20) =  &
        -k(951)*n(idx_O)  &
        +k(735)*n(idx_OHj)

    !d[HCN_dot]/d[N2]
    pd(5,20) =  &
        +k(766)*n(idx_CH2)  &
        +k(809)*n(idx_CH)

    !d[H2_dot]/d[N2]
    pd(6,20) =  &
        +k(518)*n(idx_H3j)

    !d[C_dot]/d[N2]
    pd(7,20) =  &
        -k(748)*n(idx_C)

    !d[H_dot]/d[N2]
    pd(8,20) =  &
        +k(453)*n(idx_H2j)

    !d[O2_dot]/d[N2]
    pd(11,20) =  &
        +k(631)*n(idx_O2Hj)

    !d[CH2_dot]/d[N2]
    pd(12,20) =  &
        -k(766)*n(idx_CH2)

    !d[NO_dot]/d[N2]
    pd(17,20) =  &
        +k(951)*n(idx_O)  &
        +k(630)*n(idx_HNOj)

    !d[CN_dot]/d[N2]
    pd(18,20) =  &
        +k(748)*n(idx_C)

    !d[N2_dot]/d[N2]
    pd(20,20) =  &
        -k(951)*n(idx_O)  &
        -k(809)*n(idx_CH)  &
        -k(656)*n(idx_NHj)  &
        -k(735)*n(idx_OHj)  &
        -k(1082)  &
        -k(1019)  &
        -k(518)*n(idx_H3j)  &
        -k(766)*n(idx_CH2)  &
        -k(598)*n(idx_HEj)  &
        -k(127)*n(idx_HEj)  &
        -k(748)*n(idx_C)  &
        -k(453)*n(idx_H2j)  &
        -k(630)*n(idx_HNOj)  &
        -k(712)*n(idx_Oj)  &
        -k(241)  &
        -k(631)*n(idx_O2Hj)

    !d[N_dot]/d[N2]
    pd(24,20) =  &
        +k(809)*n(idx_CH)  &
        +k(712)*n(idx_Oj)  &
        +k(951)*n(idx_O)  &
        +k(656)*n(idx_NHj)  &
        +k(598)*n(idx_HEj)  &
        +2.d0*k(1019)  &
        +2.d0*k(241)  &
        +k(748)*n(idx_C)

    !d[NH_dot]/d[N2]
    pd(25,20) =  &
        +k(766)*n(idx_CH2)

    !d[HE_dot]/d[N2]
    pd(26,20) =  &
        +k(127)*n(idx_HEj)  &
        +k(598)*n(idx_HEj)

    !d[N2_DUST_dot]/d[N2]
    pd(43,20) =  &
        +k(1082)

    !d[NO+_dot]/d[N2]
    pd(63,20) =  &
        +k(712)*n(idx_Oj)

    !d[N2+_dot]/d[N2]
    pd(66,20) =  &
        +k(127)*n(idx_HEj)

    !d[O+_dot]/d[N2]
    pd(70,20) =  &
        -k(712)*n(idx_Oj)

    !d[OH+_dot]/d[N2]
    pd(71,20) =  &
        -k(735)*n(idx_OHj)

    !d[N+_dot]/d[N2]
    pd(74,20) =  &
        +k(598)*n(idx_HEj)

    !d[NH+_dot]/d[N2]
    pd(76,20) =  &
        -k(656)*n(idx_NHj)

    !d[H2+_dot]/d[N2]
    pd(77,20) =  &
        -k(453)*n(idx_H2j)

    !d[HE+_dot]/d[N2]
    pd(78,20) =  &
        -k(598)*n(idx_HEj)  &
        -k(127)*n(idx_HEj)

    !d[HNO+_dot]/d[N2]
    pd(79,20) =  &
        -k(630)*n(idx_HNOj)

    !d[H3+_dot]/d[N2]
    pd(81,20) =  &
        -k(518)*n(idx_H3j)

    !d[N2H+_dot]/d[N2]
    pd(87,20) =  &
        +k(735)*n(idx_OHj)  &
        +k(631)*n(idx_O2Hj)  &
        +k(630)*n(idx_HNOj)  &
        +k(656)*n(idx_NHj)  &
        +k(453)*n(idx_H2j)  &
        +k(518)*n(idx_H3j)

    !d[O2H+_dot]/d[N2]
    pd(88,20) =  &
        -k(631)*n(idx_O2Hj)

    !d[E_dot]/d[NH2]
    pd(1,21) =  &
        +k(243)  &
        +k(1021)

    !d[CH_dot]/d[NH2]
    pd(2,21) =  &
        +k(751)*n(idx_C)

    !d[O_dot]/d[NH2]
    pd(3,21) =  &
        +k(912)*n(idx_OH)  &
        -k(952)*n(idx_O)  &
        +k(190)*n(idx_Oj)  &
        +k(686)*n(idx_OHj)  &
        -k(953)*n(idx_O)

    !d[HNC_dot]/d[NH2]
    pd(4,21) =  &
        +k(750)*n(idx_C)  &
        +k(681)*n(idx_HCNHj)

    !d[HCN_dot]/d[NH2]
    pd(5,21) =  &
        +k(680)*n(idx_HCNHj)  &
        +k(749)*n(idx_C)

    !d[H2_dot]/d[NH2]
    pd(6,21) =  &
        +k(94)*n(idx_H2j)  &
        -k(842)*n(idx_H2)  &
        +k(519)*n(idx_H3j)  &
        +k(353)*n(idx_CHj)  &
        +k(864)*n(idx_H)  &
        +k(599)*n(idx_HEj)

    !d[C_dot]/d[NH2]
    pd(7,21) =  &
        -k(749)*n(idx_C)  &
        -k(751)*n(idx_C)  &
        -k(750)*n(idx_C)

    !d[H_dot]/d[NH2]
    pd(8,21) =  &
        +k(750)*n(idx_C)  &
        +k(78)*n(idx_Hj)  &
        +k(749)*n(idx_C)  &
        +k(244)  &
        +k(600)*n(idx_HEj)  &
        +k(842)*n(idx_H2)  &
        +k(322)*n(idx_Cj)  &
        +k(910)*n(idx_NO)  &
        -k(864)*n(idx_H)  &
        +k(1022)  &
        +k(952)*n(idx_O)

    !d[H2O_dot]/d[NH2]
    pd(9,21) =  &
        +k(678)*n(idx_H3Oj)  &
        +k(909)*n(idx_NO)  &
        +k(165)*n(idx_H2Oj)  &
        +k(911)*n(idx_OH)

    !d[OH_dot]/d[NH2]
    pd(10,21) =  &
        +k(168)*n(idx_OHj)  &
        -k(911)*n(idx_OH)  &
        -k(912)*n(idx_OH)  &
        +k(676)*n(idx_H2Oj)  &
        +k(953)*n(idx_O)  &
        +k(910)*n(idx_NO)

    !d[O2_dot]/d[NH2]
    pd(11,21) =  &
        +k(685)*n(idx_O2Hj)  &
        +k(167)*n(idx_O2j)

    !d[H2CO_dot]/d[NH2]
    pd(13,21) =  &
        +k(677)*n(idx_H3COj)

    !d[HCO_dot]/d[NH2]
    pd(14,21) =  &
        +k(675)*n(idx_H2COj)

    !d[NH3_dot]/d[NH2]
    pd(16,21) =  &
        +k(912)*n(idx_OH)  &
        +k(842)*n(idx_H2)  &
        +k(908)*n(idx_CH4)

    !d[NO_dot]/d[NH2]
    pd(17,21) =  &
        +k(683)*n(idx_HNOj)  &
        -k(910)*n(idx_NO)  &
        -k(909)*n(idx_NO)

    !d[CN_dot]/d[NH2]
    pd(18,21) =  &
        +k(163)*n(idx_CNj)  &
        +k(679)*n(idx_HCNj)

    !d[CO_dot]/d[NH2]
    pd(19,21) =  &
        +k(682)*n(idx_HCOj)  &
        +k(164)*n(idx_COj)

    !d[N2_dot]/d[NH2]
    pd(20,21) =  &
        +k(684)*n(idx_N2Hj)  &
        +k(910)*n(idx_NO)  &
        +k(909)*n(idx_NO)  &
        +k(166)*n(idx_N2j)

    !d[NH2_dot]/d[NH2]
    pd(21,21) =  &
        -k(789)*n(idx_CH3)  &
        -k(679)*n(idx_HCNj)  &
        -k(78)*n(idx_Hj)  &
        -k(144)*n(idx_Nj)  &
        -k(683)*n(idx_HNOj)  &
        -k(1022)  &
        -k(685)*n(idx_O2Hj)  &
        -k(842)*n(idx_H2)  &
        -k(519)*n(idx_H3j)  &
        -k(911)*n(idx_OH)  &
        -k(322)*n(idx_Cj)  &
        -k(1021)  &
        -k(681)*n(idx_HCNHj)  &
        -k(676)*n(idx_H2Oj)  &
        -k(909)*n(idx_NO)  &
        -k(353)*n(idx_CHj)  &
        -k(94)*n(idx_H2j)  &
        -k(680)*n(idx_HCNHj)  &
        -k(163)*n(idx_CNj)  &
        -k(671)*n(idx_NH2j)  &
        -k(600)*n(idx_HEj)  &
        -k(675)*n(idx_H2COj)  &
        -k(678)*n(idx_H3Oj)  &
        -k(243)  &
        -k(908)*n(idx_CH4)  &
        -k(953)*n(idx_O)  &
        -k(190)*n(idx_Oj)  &
        -k(164)*n(idx_COj)  &
        -k(165)*n(idx_H2Oj)  &
        -k(244)  &
        -k(751)*n(idx_C)  &
        -k(952)*n(idx_O)  &
        -k(674)*n(idx_COj)  &
        -k(682)*n(idx_HCOj)  &
        -k(166)*n(idx_N2j)  &
        -k(167)*n(idx_O2j)  &
        -k(684)*n(idx_N2Hj)  &
        -k(912)*n(idx_OH)  &
        -k(1109)  &
        -k(599)*n(idx_HEj)  &
        -k(910)*n(idx_NO)  &
        -k(677)*n(idx_H3COj)  &
        -k(749)*n(idx_C)  &
        -k(168)*n(idx_OHj)  &
        -k(657)*n(idx_NHj)  &
        -k(686)*n(idx_OHj)  &
        -k(864)*n(idx_H)  &
        -k(750)*n(idx_C)

    !d[CH3_dot]/d[NH2]
    pd(22,21) =  &
        +k(908)*n(idx_CH4)  &
        -k(789)*n(idx_CH3)

    !d[CH4_dot]/d[NH2]
    pd(23,21) =  &
        +k(789)*n(idx_CH3)  &
        -k(908)*n(idx_CH4)

    !d[N_dot]/d[NH2]
    pd(24,21) =  &
        +k(144)*n(idx_Nj)  &
        +k(657)*n(idx_NHj)

    !d[NH_dot]/d[NH2]
    pd(25,21) =  &
        +k(671)*n(idx_NH2j)  &
        +k(244)  &
        +k(911)*n(idx_OH)  &
        +k(953)*n(idx_O)  &
        +k(789)*n(idx_CH3)  &
        +k(674)*n(idx_COj)  &
        +k(864)*n(idx_H)  &
        +k(751)*n(idx_C)  &
        +k(1022)

    !d[HE_dot]/d[NH2]
    pd(26,21) =  &
        +k(599)*n(idx_HEj)  &
        +k(600)*n(idx_HEj)

    !d[HNO_dot]/d[NH2]
    pd(27,21) =  &
        +k(952)*n(idx_O)

    !d[NH3_DUST_dot]/d[NH2]
    pd(45,21) =  &
        +k(1109)

    !d[HCO+_dot]/d[NH2]
    pd(54,21) =  &
        +k(674)*n(idx_COj)  &
        -k(682)*n(idx_HCOj)

    !d[H+_dot]/d[NH2]
    pd(55,21) =  &
        -k(78)*n(idx_Hj)

    !d[C+_dot]/d[NH2]
    pd(57,21) =  &
        -k(322)*n(idx_Cj)

    !d[CH+_dot]/d[NH2]
    pd(59,21) =  &
        -k(353)*n(idx_CHj)

    !d[H2CO+_dot]/d[NH2]
    pd(60,21) =  &
        -k(675)*n(idx_H2COj)

    !d[NH3+_dot]/d[NH2]
    pd(62,21) =  &
        +k(675)*n(idx_H2COj)  &
        +k(671)*n(idx_NH2j)  &
        +k(677)*n(idx_H3COj)  &
        +k(681)*n(idx_HCNHj)  &
        +k(682)*n(idx_HCOj)  &
        +k(684)*n(idx_N2Hj)  &
        +k(678)*n(idx_H3Oj)  &
        +k(683)*n(idx_HNOj)  &
        +k(519)*n(idx_H3j)  &
        +k(676)*n(idx_H2Oj)  &
        +k(680)*n(idx_HCNHj)  &
        +k(685)*n(idx_O2Hj)  &
        +k(657)*n(idx_NHj)  &
        +k(679)*n(idx_HCNj)  &
        +k(686)*n(idx_OHj)

    !d[CN+_dot]/d[NH2]
    pd(64,21) =  &
        -k(163)*n(idx_CNj)

    !d[CO+_dot]/d[NH2]
    pd(65,21) =  &
        -k(674)*n(idx_COj)  &
        -k(164)*n(idx_COj)

    !d[N2+_dot]/d[NH2]
    pd(66,21) =  &
        -k(166)*n(idx_N2j)

    !d[O2+_dot]/d[NH2]
    pd(67,21) =  &
        -k(167)*n(idx_O2j)

    !d[H2O+_dot]/d[NH2]
    pd(68,21) =  &
        -k(165)*n(idx_H2Oj)  &
        -k(676)*n(idx_H2Oj)

    !d[NH2+_dot]/d[NH2]
    pd(69,21) =  &
        +k(144)*n(idx_Nj)  &
        +k(168)*n(idx_OHj)  &
        +k(167)*n(idx_O2j)  &
        +k(94)*n(idx_H2j)  &
        +k(78)*n(idx_Hj)  &
        +k(190)*n(idx_Oj)  &
        +k(164)*n(idx_COj)  &
        +k(163)*n(idx_CNj)  &
        +k(165)*n(idx_H2Oj)  &
        -k(671)*n(idx_NH2j)  &
        +k(1021)  &
        +k(166)*n(idx_N2j)  &
        +k(243)

    !d[O+_dot]/d[NH2]
    pd(70,21) =  &
        -k(190)*n(idx_Oj)

    !d[OH+_dot]/d[NH2]
    pd(71,21) =  &
        -k(168)*n(idx_OHj)  &
        -k(686)*n(idx_OHj)

    !d[N+_dot]/d[NH2]
    pd(74,21) =  &
        -k(144)*n(idx_Nj)  &
        +k(599)*n(idx_HEj)

    !d[HCN+_dot]/d[NH2]
    pd(75,21) =  &
        -k(679)*n(idx_HCNj)  &
        +k(353)*n(idx_CHj)  &
        +k(322)*n(idx_Cj)

    !d[NH+_dot]/d[NH2]
    pd(76,21) =  &
        +k(600)*n(idx_HEj)  &
        -k(657)*n(idx_NHj)

    !d[H2+_dot]/d[NH2]
    pd(77,21) =  &
        -k(94)*n(idx_H2j)

    !d[HE+_dot]/d[NH2]
    pd(78,21) =  &
        -k(599)*n(idx_HEj)  &
        -k(600)*n(idx_HEj)

    !d[HNO+_dot]/d[NH2]
    pd(79,21) =  &
        -k(683)*n(idx_HNOj)

    !d[H3+_dot]/d[NH2]
    pd(81,21) =  &
        -k(519)*n(idx_H3j)

    !d[H3CO+_dot]/d[NH2]
    pd(82,21) =  &
        -k(677)*n(idx_H3COj)

    !d[H3O+_dot]/d[NH2]
    pd(83,21) =  &
        -k(678)*n(idx_H3Oj)

    !d[HCNH+_dot]/d[NH2]
    pd(84,21) =  &
        -k(681)*n(idx_HCNHj)  &
        -k(680)*n(idx_HCNHj)

    !d[N2H+_dot]/d[NH2]
    pd(87,21) =  &
        -k(684)*n(idx_N2Hj)

    !d[O2H+_dot]/d[NH2]
    pd(88,21) =  &
        -k(685)*n(idx_O2Hj)

    !d[E_dot]/d[CH3]
    pd(1,22) =  &
        +k(220)  &
        +k(983)

    !d[CH_dot]/d[CH3]
    pd(2,22) =  &
        +k(221)  &
        +k(984)

    !d[O_dot]/d[CH3]
    pd(3,22) =  &
        -k(797)*n(idx_O)  &
        +k(799)*n(idx_OH)  &
        -k(798)*n(idx_O)

    !d[HCN_dot]/d[CH3]
    pd(5,22) =  &
        +k(792)*n(idx_NO)  &
        +k(784)*n(idx_CN)  &
        +k(891)*n(idx_N)  &
        +k(890)*n(idx_N)

    !d[H2_dot]/d[CH3]
    pd(6,22) =  &
        +k(984)  &
        +k(849)*n(idx_H)  &
        +k(797)*n(idx_O)  &
        +k(504)*n(idx_H3j)  &
        +k(221)  &
        +k(566)*n(idx_HEj)  &
        +k(890)*n(idx_N)  &
        +k(800)*n(idx_OH)  &
        -k(838)*n(idx_H2)

    !d[H_dot]/d[CH3]
    pd(8,22) =  &
        +k(219)  &
        +k(70)*n(idx_Hj)  &
        +k(797)*n(idx_O)  &
        +k(838)*n(idx_H2)  &
        -k(849)*n(idx_H)  &
        +k(798)*n(idx_O)  &
        +k(889)*n(idx_N)  &
        +2.d0*k(891)*n(idx_N)  &
        +k(982)

    !d[H2O_dot]/d[CH3]
    pd(9,22) =  &
        -k(786)*n(idx_H2O)  &
        +k(801)*n(idx_OH)  &
        +k(792)*n(idx_NO)  &
        +k(794)*n(idx_O2)

    !d[OH_dot]/d[CH3]
    pd(10,22) =  &
        +k(786)*n(idx_H2O)  &
        -k(799)*n(idx_OH)  &
        +k(793)*n(idx_O2)  &
        -k(801)*n(idx_OH)  &
        -k(800)*n(idx_OH)

    !d[O2_dot]/d[CH3]
    pd(11,22) =  &
        -k(794)*n(idx_O2)  &
        +k(796)*n(idx_O2H)  &
        -k(793)*n(idx_O2)  &
        -k(795)*n(idx_O2)

    !d[CH2_dot]/d[CH3]
    pd(12,22) =  &
        +k(219)  &
        +k(849)*n(idx_H)  &
        +k(801)*n(idx_OH)  &
        +k(795)*n(idx_O2)  &
        +k(784)*n(idx_CN)  &
        +k(982)  &
        +2.d0*k(783)*n(idx_CH3)

    !d[H2CO_dot]/d[CH3]
    pd(13,22) =  &
        +k(800)*n(idx_OH)  &
        -k(785)*n(idx_H2CO)  &
        +k(793)*n(idx_O2)  &
        +k(791)*n(idx_NO2)  &
        +k(798)*n(idx_O)

    !d[HCO_dot]/d[CH3]
    pd(14,22) =  &
        +k(794)*n(idx_O2)  &
        +k(785)*n(idx_H2CO)  &
        -k(787)*n(idx_HCO)

    !d[NH3_dot]/d[CH3]
    pd(16,22) =  &
        -k(790)*n(idx_NH3)

    !d[NO_dot]/d[CH3]
    pd(17,22) =  &
        -k(792)*n(idx_NO)  &
        +k(788)*n(idx_HNO)

    !d[CN_dot]/d[CH3]
    pd(18,22) =  &
        -k(784)*n(idx_CN)

    !d[CO_dot]/d[CH3]
    pd(19,22) =  &
        +k(787)*n(idx_HCO)  &
        +k(797)*n(idx_O)

    !d[NH2_dot]/d[CH3]
    pd(21,22) =  &
        -k(789)*n(idx_NH2)  &
        +k(790)*n(idx_NH3)

    !d[CH3_dot]/d[CH3]
    pd(22,22) =  &
        -k(70)*n(idx_Hj)  &
        -k(800)*n(idx_OH)  &
        -k(798)*n(idx_O)  &
        -k(785)*n(idx_H2CO)  &
        -k(787)*n(idx_HCO)  &
        -k(566)*n(idx_HEj)  &
        -k(789)*n(idx_NH2)  &
        -k(890)*n(idx_N)  &
        -k(788)*n(idx_HNO)  &
        -k(793)*n(idx_O2)  &
        -k(791)*n(idx_NO2)  &
        -k(801)*n(idx_OH)  &
        -k(784)*n(idx_CN)  &
        -k(1079)  &
        -k(984)  &
        -k(889)*n(idx_N)  &
        -k(220)  &
        -k(982)  &
        -k(796)*n(idx_O2H)  &
        -k(799)*n(idx_OH)  &
        -k(838)*n(idx_H2)  &
        -k(794)*n(idx_O2)  &
        -k(891)*n(idx_N)  &
        -k(790)*n(idx_NH3)  &
        -k(797)*n(idx_O)  &
        -k(786)*n(idx_H2O)  &
        -k(849)*n(idx_H)  &
        -k(504)*n(idx_H3j)  &
        -k(221)  &
        -k(219)  &
        -k(795)*n(idx_O2)  &
        -4.d0*k(783)*n(idx_CH3)  &
        -k(983)  &
        -k(792)*n(idx_NO)

    !d[CH4_dot]/d[CH3]
    pd(23,22) =  &
        +k(785)*n(idx_H2CO)  &
        +k(789)*n(idx_NH2)  &
        +k(788)*n(idx_HNO)  &
        +k(799)*n(idx_OH)  &
        +k(838)*n(idx_H2)  &
        +k(796)*n(idx_O2H)  &
        +k(790)*n(idx_NH3)  &
        +k(786)*n(idx_H2O)  &
        +2.d0*k(783)*n(idx_CH3)  &
        +k(787)*n(idx_HCO)

    !d[N_dot]/d[CH3]
    pd(24,22) =  &
        -k(889)*n(idx_N)  &
        -k(891)*n(idx_N)  &
        -k(890)*n(idx_N)

    !d[NH_dot]/d[CH3]
    pd(25,22) =  &
        +k(789)*n(idx_NH2)

    !d[HE_dot]/d[CH3]
    pd(26,22) =  &
        +k(566)*n(idx_HEj)

    !d[HNO_dot]/d[CH3]
    pd(27,22) =  &
        +k(791)*n(idx_NO2)  &
        -k(788)*n(idx_HNO)

    !d[H2CN_dot]/d[CH3]
    pd(30,22) =  &
        +k(889)*n(idx_N)

    !d[NO2_dot]/d[CH3]
    pd(32,22) =  &
        -k(791)*n(idx_NO2)

    !d[O2H_dot]/d[CH3]
    pd(33,22) =  &
        +k(795)*n(idx_O2)  &
        -k(796)*n(idx_O2H)

    !d[CH4_DUST_dot]/d[CH3]
    pd(38,22) =  &
        +k(1079)

    !d[H+_dot]/d[CH3]
    pd(55,22) =  &
        -k(70)*n(idx_Hj)

    !d[CH+_dot]/d[CH3]
    pd(59,22) =  &
        +k(566)*n(idx_HEj)

    !d[CH3+_dot]/d[CH3]
    pd(72,22) =  &
        +k(70)*n(idx_Hj)  &
        +k(220)  &
        +k(983)

    !d[CH4+_dot]/d[CH3]
    pd(73,22) =  &
        +k(504)*n(idx_H3j)

    !d[HE+_dot]/d[CH3]
    pd(78,22) =  &
        -k(566)*n(idx_HEj)

    !d[H3+_dot]/d[CH3]
    pd(81,22) =  &
        -k(504)*n(idx_H3j)

    !d[E_dot]/d[CH4]
    pd(1,23) =  &
        +k(992)

    !d[CH_dot]/d[CH4]
    pd(2,23) =  &
        +k(993)

    !d[O_dot]/d[CH4]
    pd(3,23) =  &
        -k(936)*n(idx_O)  &
        +k(185)*n(idx_Oj)

    !d[HCN_dot]/d[CH4]
    pd(5,23) =  &
        +k(802)*n(idx_CN)

    !d[H2_dot]/d[CH4]
    pd(6,23) =  &
        +k(850)*n(idx_H)  &
        +k(569)*n(idx_HEj)  &
        +k(432)*n(idx_Hj)  &
        +k(570)*n(idx_HEj)  &
        +k(86)*n(idx_H2j)  &
        +k(993)  &
        +k(615)*n(idx_Nj)  &
        +k(990)  &
        +k(224)  &
        +k(443)*n(idx_H2j)  &
        +k(397)*n(idx_N2j)

    !d[H_dot]/d[CH4]
    pd(8,23) =  &
        +k(614)*n(idx_Nj)  &
        +k(398)*n(idx_N2j)  &
        +k(569)*n(idx_HEj)  &
        +k(571)*n(idx_HEj)  &
        -k(850)*n(idx_H)  &
        +2.d0*k(616)*n(idx_Nj)  &
        +k(993)  &
        +k(71)*n(idx_Hj)  &
        +k(615)*n(idx_Nj)  &
        +k(991)  &
        +k(443)*n(idx_H2j)

    !d[H2O_dot]/d[CH4]
    pd(9,23) =  &
        +k(804)*n(idx_OH)

    !d[OH_dot]/d[CH4]
    pd(10,23) =  &
        -k(804)*n(idx_OH)  &
        +k(705)*n(idx_Oj)  &
        +k(936)*n(idx_O)

    !d[O2_dot]/d[CH4]
    pd(11,23) =  &
        -k(803)*n(idx_O2)

    !d[CH2_dot]/d[CH4]
    pd(12,23) =  &
        +k(224)  &
        +k(399)*n(idx_OHj)  &
        -k(761)*n(idx_CH2)  &
        +k(990)

    !d[NH3_dot]/d[CH4]
    pd(16,23) =  &
        +k(908)*n(idx_NH2)

    !d[CN_dot]/d[CH4]
    pd(18,23) =  &
        -k(802)*n(idx_CN)

    !d[CO_dot]/d[CH4]
    pd(19,23) =  &
        +k(46)*n(idx_COj)

    !d[N2_dot]/d[CH4]
    pd(20,23) =  &
        +k(398)*n(idx_N2j)  &
        +k(397)*n(idx_N2j)

    !d[NH2_dot]/d[CH4]
    pd(21,23) =  &
        -k(908)*n(idx_NH2)  &
        +k(914)*n(idx_NH)

    !d[CH3_dot]/d[CH4]
    pd(22,23) =  &
        +k(395)*n(idx_H2Oj)  &
        +k(803)*n(idx_O2)  &
        +2.d0*k(761)*n(idx_CH2)  &
        +k(908)*n(idx_NH2)  &
        +k(396)*n(idx_HCNj)  &
        +k(936)*n(idx_O)  &
        +k(393)*n(idx_COj)  &
        +k(394)*n(idx_H2COj)  &
        +k(850)*n(idx_H)  &
        +k(804)*n(idx_OH)  &
        +k(991)  &
        +k(914)*n(idx_NH)  &
        +k(572)*n(idx_HEj)  &
        +k(802)*n(idx_CN)

    !d[CH4_dot]/d[CH4]
    pd(23,23) =  &
        -k(803)*n(idx_O2)  &
        -k(394)*n(idx_H2COj)  &
        -k(936)*n(idx_O)  &
        -k(761)*n(idx_CH2)  &
        -k(46)*n(idx_COj)  &
        -k(914)*n(idx_NH)  &
        -k(908)*n(idx_NH2)  &
        -k(569)*n(idx_HEj)  &
        -k(571)*n(idx_HEj)  &
        -k(990)  &
        -k(224)  &
        -k(993)  &
        -k(992)  &
        -k(850)*n(idx_H)  &
        -k(802)*n(idx_CN)  &
        -k(123)*n(idx_HEj)  &
        -k(1080)  &
        -k(432)*n(idx_Hj)  &
        -k(396)*n(idx_HCNj)  &
        -k(615)*n(idx_Nj)  &
        -k(136)*n(idx_Nj)  &
        -k(804)*n(idx_OH)  &
        -k(393)*n(idx_COj)  &
        -k(572)*n(idx_HEj)  &
        -k(705)*n(idx_Oj)  &
        -k(71)*n(idx_Hj)  &
        -k(570)*n(idx_HEj)  &
        -k(399)*n(idx_OHj)  &
        -k(185)*n(idx_Oj)  &
        -k(614)*n(idx_Nj)  &
        -k(86)*n(idx_H2j)  &
        -k(443)*n(idx_H2j)  &
        -k(397)*n(idx_N2j)  &
        -k(616)*n(idx_Nj)  &
        -k(398)*n(idx_N2j)  &
        -k(395)*n(idx_H2Oj)  &
        -k(991)

    !d[N_dot]/d[CH4]
    pd(24,23) =  &
        +k(614)*n(idx_Nj)  &
        +k(136)*n(idx_Nj)

    !d[NH_dot]/d[CH4]
    pd(25,23) =  &
        -k(914)*n(idx_NH)

    !d[HE_dot]/d[CH4]
    pd(26,23) =  &
        +k(571)*n(idx_HEj)  &
        +k(570)*n(idx_HEj)  &
        +k(572)*n(idx_HEj)  &
        +k(123)*n(idx_HEj)  &
        +k(569)*n(idx_HEj)

    !d[O2H_dot]/d[CH4]
    pd(33,23) =  &
        +k(803)*n(idx_O2)

    !d[CH4_DUST_dot]/d[CH4]
    pd(38,23) =  &
        +k(1080)

    !d[HCO+_dot]/d[CH4]
    pd(54,23) =  &
        +k(393)*n(idx_COj)

    !d[H+_dot]/d[CH4]
    pd(55,23) =  &
        +k(572)*n(idx_HEj)  &
        -k(432)*n(idx_Hj)  &
        -k(71)*n(idx_Hj)

    !d[CH2+_dot]/d[CH4]
    pd(58,23) =  &
        +k(570)*n(idx_HEj)  &
        +k(397)*n(idx_N2j)

    !d[CH+_dot]/d[CH4]
    pd(59,23) =  &
        +k(569)*n(idx_HEj)

    !d[H2CO+_dot]/d[CH4]
    pd(60,23) =  &
        -k(394)*n(idx_H2COj)

    !d[CO+_dot]/d[CH4]
    pd(65,23) =  &
        -k(46)*n(idx_COj)  &
        -k(393)*n(idx_COj)

    !d[N2+_dot]/d[CH4]
    pd(66,23) =  &
        -k(397)*n(idx_N2j)  &
        -k(398)*n(idx_N2j)

    !d[H2O+_dot]/d[CH4]
    pd(68,23) =  &
        -k(395)*n(idx_H2Oj)

    !d[O+_dot]/d[CH4]
    pd(70,23) =  &
        -k(705)*n(idx_Oj)  &
        -k(185)*n(idx_Oj)

    !d[OH+_dot]/d[CH4]
    pd(71,23) =  &
        -k(399)*n(idx_OHj)

    !d[CH3+_dot]/d[CH4]
    pd(72,23) =  &
        +k(614)*n(idx_Nj)  &
        +k(398)*n(idx_N2j)  &
        +k(432)*n(idx_Hj)  &
        +k(571)*n(idx_HEj)  &
        +k(705)*n(idx_Oj)  &
        +k(443)*n(idx_H2j)

    !d[CH4+_dot]/d[CH4]
    pd(73,23) =  &
        +k(992)  &
        +k(185)*n(idx_Oj)  &
        +k(123)*n(idx_HEj)  &
        +k(86)*n(idx_H2j)  &
        +k(136)*n(idx_Nj)  &
        +k(46)*n(idx_COj)  &
        +k(71)*n(idx_Hj)

    !d[N+_dot]/d[CH4]
    pd(74,23) =  &
        -k(616)*n(idx_Nj)  &
        -k(615)*n(idx_Nj)  &
        -k(136)*n(idx_Nj)  &
        -k(614)*n(idx_Nj)

    !d[HCN+_dot]/d[CH4]
    pd(75,23) =  &
        -k(396)*n(idx_HCNj)  &
        +k(615)*n(idx_Nj)

    !d[H2+_dot]/d[CH4]
    pd(77,23) =  &
        -k(443)*n(idx_H2j)  &
        -k(86)*n(idx_H2j)

    !d[HE+_dot]/d[CH4]
    pd(78,23) =  &
        -k(569)*n(idx_HEj)  &
        -k(571)*n(idx_HEj)  &
        -k(123)*n(idx_HEj)  &
        -k(570)*n(idx_HEj)  &
        -k(572)*n(idx_HEj)

    !d[H3CO+_dot]/d[CH4]
    pd(82,23) =  &
        +k(394)*n(idx_H2COj)

    !d[H3O+_dot]/d[CH4]
    pd(83,23) =  &
        +k(399)*n(idx_OHj)  &
        +k(395)*n(idx_H2Oj)

    !d[HCNH+_dot]/d[CH4]
    pd(84,23) =  &
        +k(616)*n(idx_Nj)  &
        +k(396)*n(idx_HCNj)

    !d[E_dot]/d[N]
    pd(1,24) =  &
        +k(242)  &
        +k(213)

    !d[CH_dot]/d[N]
    pd(2,24) =  &
        -k(811)*n(idx_CH)  &
        +k(888)*n(idx_CH2)  &
        -k(810)*n(idx_CH)

    !d[O_dot]/d[N]
    pd(3,24) =  &
        +k(640)*n(idx_O2j)  &
        +2.d0*k(900)*n(idx_NO2)  &
        +k(896)*n(idx_HCO)  &
        +k(904)*n(idx_O2)  &
        +k(903)*n(idx_NO)  &
        +k(907)*n(idx_OH)

    !d[HNC_dot]/d[N]
    pd(4,24) =  &
        +k(887)*n(idx_CH2)

    !d[HCN_dot]/d[N]
    pd(5,24) =  &
        +k(891)*n(idx_CH3)  &
        +k(896)*n(idx_HCO)  &
        +k(886)*n(idx_CH2)  &
        +k(890)*n(idx_CH3)  &
        +k(894)*n(idx_H2CN)

    !d[H2_dot]/d[N]
    pd(6,24) =  &
        -k(841)*n(idx_H2)  &
        +k(890)*n(idx_CH3)  &
        +k(637)*n(idx_H2Oj)

    !d[C_dot]/d[N]
    pd(7,24) =  &
        +k(635)*n(idx_CNj)  &
        +k(811)*n(idx_CH)  &
        -k(1042)*n(idx_C)  &
        +k(892)*n(idx_CN)

    !d[H_dot]/d[N]
    pd(8,24) =  &
        +k(906)*n(idx_OH)  &
        +k(638)*n(idx_NHj)  &
        +k(641)*n(idx_OHj)  &
        +k(889)*n(idx_CH3)  &
        +k(887)*n(idx_CH2)  &
        +k(639)*n(idx_NH2j)  &
        +k(899)*n(idx_NH)  &
        +2.d0*k(891)*n(idx_CH3)  &
        +k(636)*n(idx_H2Oj)  &
        +k(897)*n(idx_HCO)  &
        +k(886)*n(idx_CH2)  &
        +k(810)*n(idx_CH)  &
        +k(841)*n(idx_H2)  &
        +k(352)*n(idx_CHj)  &
        +k(454)*n(idx_H2j)  &
        +k(634)*n(idx_CH2j)

    !d[OH_dot]/d[N]
    pd(10,24) =  &
        -k(907)*n(idx_OH)  &
        -k(906)*n(idx_OH)

    !d[O2_dot]/d[N]
    pd(11,24) =  &
        -k(904)*n(idx_O2)  &
        +k(905)*n(idx_O2H)  &
        +k(902)*n(idx_NO2)

    !d[CH2_dot]/d[N]
    pd(12,24) =  &
        -k(887)*n(idx_CH2)  &
        -k(886)*n(idx_CH2)  &
        -k(888)*n(idx_CH2)

    !d[HCO_dot]/d[N]
    pd(14,24) =  &
        -k(897)*n(idx_HCO)  &
        -k(895)*n(idx_HCO)  &
        -k(896)*n(idx_HCO)

    !d[NO_dot]/d[N]
    pd(17,24) =  &
        -k(903)*n(idx_NO)  &
        +k(893)*n(idx_CO2)  &
        +2.d0*k(901)*n(idx_NO2)  &
        +k(906)*n(idx_OH)  &
        +k(898)*n(idx_HNO)  &
        +k(904)*n(idx_O2)

    !d[CN_dot]/d[N]
    pd(18,24) =  &
        +k(810)*n(idx_CH)  &
        +k(1042)*n(idx_C)  &
        -k(892)*n(idx_CN)

    !d[CO_dot]/d[N]
    pd(19,24) =  &
        +k(895)*n(idx_HCO)  &
        +k(893)*n(idx_CO2)

    !d[N2_dot]/d[N]
    pd(20,24) =  &
        +k(892)*n(idx_CN)  &
        +k(900)*n(idx_NO2)  &
        +k(899)*n(idx_NH)  &
        +k(154)*n(idx_N2j)  &
        +k(903)*n(idx_NO)  &
        +k(902)*n(idx_NO2)

    !d[CH3_dot]/d[N]
    pd(22,24) =  &
        -k(891)*n(idx_CH3)  &
        -k(889)*n(idx_CH3)  &
        -k(890)*n(idx_CH3)

    !d[N_dot]/d[N]
    pd(24,24) =  &
        -k(907)*n(idx_OH)  &
        -k(905)*n(idx_O2H)  &
        -k(891)*n(idx_CH3)  &
        -k(895)*n(idx_HCO)  &
        -k(841)*n(idx_H2)  &
        -k(897)*n(idx_HCO)  &
        -k(640)*n(idx_O2j)  &
        -k(899)*n(idx_NH)  &
        -k(810)*n(idx_CH)  &
        -k(906)*n(idx_OH)  &
        -k(902)*n(idx_NO2)  &
        -k(213)  &
        -k(892)*n(idx_CN)  &
        -k(888)*n(idx_CH2)  &
        -k(639)*n(idx_NH2j)  &
        -k(352)*n(idx_CHj)  &
        -k(637)*n(idx_H2Oj)  &
        -k(894)*n(idx_H2CN)  &
        -k(889)*n(idx_CH3)  &
        -k(890)*n(idx_CH3)  &
        -k(903)*n(idx_NO)  &
        -k(1054)*n(idx_Nj)  &
        -k(454)*n(idx_H2j)  &
        -k(898)*n(idx_HNO)  &
        -k(1110)  &
        -k(887)*n(idx_CH2)  &
        -k(634)*n(idx_CH2j)  &
        -k(893)*n(idx_CO2)  &
        -k(1040)*n(idx_Cj)  &
        -k(896)*n(idx_HCO)  &
        -k(886)*n(idx_CH2)  &
        -k(638)*n(idx_NHj)  &
        -k(154)*n(idx_N2j)  &
        -k(636)*n(idx_H2Oj)  &
        -k(641)*n(idx_OHj)  &
        -k(904)*n(idx_O2)  &
        -k(901)*n(idx_NO2)  &
        -k(811)*n(idx_CH)  &
        -k(900)*n(idx_NO2)  &
        -k(1042)*n(idx_C)  &
        -k(242)  &
        -k(635)*n(idx_CNj)

    !d[NH_dot]/d[N]
    pd(25,24) =  &
        +k(888)*n(idx_CH2)  &
        +k(895)*n(idx_HCO)  &
        -k(899)*n(idx_NH)  &
        +k(905)*n(idx_O2H)  &
        +k(898)*n(idx_HNO)  &
        +k(894)*n(idx_H2CN)  &
        +k(841)*n(idx_H2)  &
        +k(811)*n(idx_CH)  &
        +k(907)*n(idx_OH)

    !d[HNO_dot]/d[N]
    pd(27,24) =  &
        -k(898)*n(idx_HNO)

    !d[CO2_dot]/d[N]
    pd(29,24) =  &
        -k(893)*n(idx_CO2)

    !d[H2CN_dot]/d[N]
    pd(30,24) =  &
        -k(894)*n(idx_H2CN)  &
        +k(889)*n(idx_CH3)

    !d[NO2_dot]/d[N]
    pd(32,24) =  &
        -k(902)*n(idx_NO2)  &
        -k(900)*n(idx_NO2)  &
        -k(901)*n(idx_NO2)

    !d[O2H_dot]/d[N]
    pd(33,24) =  &
        -k(905)*n(idx_O2H)

    !d[OCN_dot]/d[N]
    pd(34,24) =  &
        +k(897)*n(idx_HCO)

    !d[NH3_DUST_dot]/d[N]
    pd(45,24) =  &
        +k(1110)

    !d[C+_dot]/d[N]
    pd(57,24) =  &
        -k(1040)*n(idx_Cj)

    !d[CH2+_dot]/d[N]
    pd(58,24) =  &
        -k(634)*n(idx_CH2j)

    !d[CH+_dot]/d[N]
    pd(59,24) =  &
        -k(352)*n(idx_CHj)

    !d[NO+_dot]/d[N]
    pd(63,24) =  &
        +k(640)*n(idx_O2j)  &
        +k(641)*n(idx_OHj)  &
        +k(637)*n(idx_H2Oj)

    !d[CN+_dot]/d[N]
    pd(64,24) =  &
        +k(1040)*n(idx_Cj)  &
        -k(635)*n(idx_CNj)  &
        +k(352)*n(idx_CHj)

    !d[N2+_dot]/d[N]
    pd(66,24) =  &
        -k(154)*n(idx_N2j)  &
        +k(635)*n(idx_CNj)  &
        +k(1054)*n(idx_Nj)  &
        +k(638)*n(idx_NHj)

    !d[O2+_dot]/d[N]
    pd(67,24) =  &
        -k(640)*n(idx_O2j)

    !d[H2O+_dot]/d[N]
    pd(68,24) =  &
        -k(637)*n(idx_H2Oj)  &
        -k(636)*n(idx_H2Oj)

    !d[NH2+_dot]/d[N]
    pd(69,24) =  &
        -k(639)*n(idx_NH2j)

    !d[OH+_dot]/d[N]
    pd(71,24) =  &
        -k(641)*n(idx_OHj)

    !d[N+_dot]/d[N]
    pd(74,24) =  &
        +k(154)*n(idx_N2j)  &
        +k(242)  &
        -k(1054)*n(idx_Nj)  &
        +k(213)

    !d[HCN+_dot]/d[N]
    pd(75,24) =  &
        +k(634)*n(idx_CH2j)

    !d[NH+_dot]/d[N]
    pd(76,24) =  &
        -k(638)*n(idx_NHj)  &
        +k(454)*n(idx_H2j)

    !d[H2+_dot]/d[N]
    pd(77,24) =  &
        -k(454)*n(idx_H2j)

    !d[HNO+_dot]/d[N]
    pd(79,24) =  &
        +k(636)*n(idx_H2Oj)

    !d[N2H+_dot]/d[N]
    pd(87,24) =  &
        +k(639)*n(idx_NH2j)

    !d[E_dot]/d[NH]
    pd(1,25) =  &
        +k(249)  &
        +k(1027)

    !d[CH_dot]/d[NH]
    pd(2,25) =  &
        +k(753)*n(idx_C)

    !d[O_dot]/d[NH]
    pd(3,25) =  &
        -k(927)*n(idx_O)  &
        +k(701)*n(idx_OHj)  &
        +k(924)*n(idx_O2)  &
        +k(930)*n(idx_OH)  &
        +k(699)*n(idx_O2j)  &
        +k(181)*n(idx_Oj)  &
        +k(922)*n(idx_NO)  &
        -k(926)*n(idx_O)

    !d[HCN_dot]/d[NH]
    pd(5,25) =  &
        +k(915)*n(idx_CN)

    !d[H2_dot]/d[NH]
    pd(6,25) =  &
        +k(866)*n(idx_H)  &
        +2.d0*k(918)*n(idx_NH)  &
        +k(520)*n(idx_H3j)  &
        +k(689)*n(idx_CH3j)  &
        +k(96)*n(idx_H2j)  &
        +k(354)*n(idx_CHj)  &
        -k(843)*n(idx_H2)

    !d[C_dot]/d[NH]
    pd(7,25) =  &
        -k(752)*n(idx_C)  &
        -k(753)*n(idx_C)

    !d[H_dot]/d[NH]
    pd(8,25) =  &
        +k(80)*n(idx_Hj)  &
        +k(698)*n(idx_Oj)  &
        +k(929)*n(idx_OH)  &
        +4.d0*k(919)*n(idx_NH)  &
        -k(866)*n(idx_H)  &
        +k(624)*n(idx_Nj)  &
        +k(899)*n(idx_N)  &
        +k(922)*n(idx_NO)  &
        +k(843)*n(idx_H2)  &
        +k(926)*n(idx_O)  &
        +k(248)  &
        +k(1026)  &
        +k(752)*n(idx_C)  &
        +k(603)*n(idx_HEj)  &
        +k(455)*n(idx_H2j)  &
        +k(324)*n(idx_Cj)

    !d[H2O_dot]/d[NH]
    pd(9,25) =  &
        +k(928)*n(idx_OH)  &
        -k(916)*n(idx_H2O)

    !d[OH_dot]/d[NH]
    pd(10,25) =  &
        +k(923)*n(idx_NO)  &
        +k(927)*n(idx_O)  &
        +k(916)*n(idx_H2O)  &
        -k(928)*n(idx_OH)  &
        -k(929)*n(idx_OH)  &
        -k(930)*n(idx_OH)  &
        +k(925)*n(idx_O2)

    !d[O2_dot]/d[NH]
    pd(11,25) =  &
        -k(924)*n(idx_O2)  &
        -k(925)*n(idx_O2)  &
        +k(700)*n(idx_O2Hj)

    !d[NH3_dot]/d[NH]
    pd(16,25) =  &
        -k(917)*n(idx_NH3)

    !d[NO_dot]/d[NH]
    pd(17,25) =  &
        -k(922)*n(idx_NO)  &
        +k(926)*n(idx_O)  &
        +k(695)*n(idx_HNOj)  &
        +k(921)*n(idx_NO2)  &
        -k(923)*n(idx_NO)  &
        +k(925)*n(idx_O2)

    !d[CN_dot]/d[NH]
    pd(18,25) =  &
        -k(915)*n(idx_CN)  &
        +k(693)*n(idx_HCNj)  &
        +k(752)*n(idx_C)  &
        +k(178)*n(idx_CNj)

    !d[CO_dot]/d[NH]
    pd(19,25) =  &
        +k(694)*n(idx_HCOj)  &
        +k(179)*n(idx_COj)

    !d[N2_dot]/d[NH]
    pd(20,25) =  &
        +2.d0*k(919)*n(idx_NH)  &
        +2.d0*k(918)*n(idx_NH)  &
        +k(923)*n(idx_NO)  &
        +k(696)*n(idx_N2Hj)  &
        +k(899)*n(idx_N)  &
        +k(180)*n(idx_N2j)  &
        +k(922)*n(idx_NO)

    !d[NH2_dot]/d[NH]
    pd(21,25) =  &
        +2.d0*k(917)*n(idx_NH3)  &
        +k(914)*n(idx_CH4)  &
        +2.d0*k(920)*n(idx_NH)  &
        +k(930)*n(idx_OH)  &
        +k(843)*n(idx_H2)  &
        +k(916)*n(idx_H2O)

    !d[CH3_dot]/d[NH]
    pd(22,25) =  &
        +k(914)*n(idx_CH4)

    !d[CH4_dot]/d[NH]
    pd(23,25) =  &
        -k(914)*n(idx_CH4)

    !d[N_dot]/d[NH]
    pd(24,25) =  &
        +k(866)*n(idx_H)  &
        +2.d0*k(920)*n(idx_NH)  &
        +k(927)*n(idx_O)  &
        -k(899)*n(idx_N)  &
        +k(915)*n(idx_CN)  &
        +k(697)*n(idx_NH2j)  &
        +k(928)*n(idx_OH)  &
        +k(692)*n(idx_H2Oj)  &
        +k(691)*n(idx_H2COj)  &
        +k(248)  &
        +k(1026)  &
        +k(146)*n(idx_Nj)  &
        +k(690)*n(idx_COj)  &
        +k(753)*n(idx_C)  &
        +k(658)*n(idx_NHj)

    !d[NH_dot]/d[NH]
    pd(25,25) =  &
        -k(1026)  &
        -k(690)*n(idx_COj)  &
        -k(930)*n(idx_OH)  &
        -k(752)*n(idx_C)  &
        -k(96)*n(idx_H2j)  &
        -k(520)*n(idx_H3j)  &
        -k(700)*n(idx_O2Hj)  &
        -k(916)*n(idx_H2O)  &
        -k(691)*n(idx_H2COj)  &
        -k(180)*n(idx_N2j)  &
        -k(692)*n(idx_H2Oj)  &
        -k(455)*n(idx_H2j)  &
        -k(249)  &
        -k(693)*n(idx_HCNj)  &
        -k(866)*n(idx_H)  &
        -k(927)*n(idx_O)  &
        -k(928)*n(idx_OH)  &
        -k(929)*n(idx_OH)  &
        -k(699)*n(idx_O2j)  &
        -k(1085)  &
        -k(843)*n(idx_H2)  &
        -k(922)*n(idx_NO)  &
        -k(80)*n(idx_Hj)  &
        -k(914)*n(idx_CH4)  &
        -k(248)  &
        -k(921)*n(idx_NO2)  &
        -k(701)*n(idx_OHj)  &
        -k(926)*n(idx_O)  &
        -k(181)*n(idx_Oj)  &
        -k(178)*n(idx_CNj)  &
        -4.d0*k(918)*n(idx_NH)  &
        -k(179)*n(idx_COj)  &
        -k(923)*n(idx_NO)  &
        -k(698)*n(idx_Oj)  &
        -k(695)*n(idx_HNOj)  &
        -k(753)*n(idx_C)  &
        -k(689)*n(idx_CH3j)  &
        -4.d0*k(920)*n(idx_NH)  &
        -k(603)*n(idx_HEj)  &
        -k(694)*n(idx_HCOj)  &
        -k(696)*n(idx_N2Hj)  &
        -k(925)*n(idx_O2)  &
        -k(899)*n(idx_N)  &
        -k(354)*n(idx_CHj)  &
        -k(146)*n(idx_Nj)  &
        -4.d0*k(919)*n(idx_NH)  &
        -k(915)*n(idx_CN)  &
        -k(658)*n(idx_NHj)  &
        -k(917)*n(idx_NH3)  &
        -k(324)*n(idx_Cj)  &
        -k(624)*n(idx_Nj)  &
        -k(924)*n(idx_O2)  &
        -k(697)*n(idx_NH2j)  &
        -k(1027)

    !d[HE_dot]/d[NH]
    pd(26,25) =  &
        +k(603)*n(idx_HEj)

    !d[HNO_dot]/d[NH]
    pd(27,25) =  &
        +k(929)*n(idx_OH)  &
        +k(921)*n(idx_NO2)  &
        +k(924)*n(idx_O2)

    !d[NO2_dot]/d[NH]
    pd(32,25) =  &
        -k(921)*n(idx_NO2)

    !d[NH3_DUST_dot]/d[NH]
    pd(45,25) =  &
        +k(1085)

    !d[HCO+_dot]/d[NH]
    pd(54,25) =  &
        -k(694)*n(idx_HCOj)  &
        +k(690)*n(idx_COj)

    !d[H+_dot]/d[NH]
    pd(55,25) =  &
        -k(80)*n(idx_Hj)

    !d[C+_dot]/d[NH]
    pd(57,25) =  &
        -k(324)*n(idx_Cj)

    !d[CH+_dot]/d[NH]
    pd(59,25) =  &
        -k(354)*n(idx_CHj)

    !d[H2CO+_dot]/d[NH]
    pd(60,25) =  &
        -k(691)*n(idx_H2COj)

    !d[NH3+_dot]/d[NH]
    pd(62,25) =  &
        +k(697)*n(idx_NH2j)

    !d[NO+_dot]/d[NH]
    pd(63,25) =  &
        +k(698)*n(idx_Oj)

    !d[CN+_dot]/d[NH]
    pd(64,25) =  &
        -k(178)*n(idx_CNj)  &
        +k(324)*n(idx_Cj)  &
        +k(354)*n(idx_CHj)

    !d[CO+_dot]/d[NH]
    pd(65,25) =  &
        -k(690)*n(idx_COj)  &
        -k(179)*n(idx_COj)

    !d[N2+_dot]/d[NH]
    pd(66,25) =  &
        -k(180)*n(idx_N2j)  &
        +k(624)*n(idx_Nj)

    !d[O2+_dot]/d[NH]
    pd(67,25) =  &
        -k(699)*n(idx_O2j)

    !d[H2O+_dot]/d[NH]
    pd(68,25) =  &
        -k(692)*n(idx_H2Oj)

    !d[NH2+_dot]/d[NH]
    pd(69,25) =  &
        +k(520)*n(idx_H3j)  &
        +k(700)*n(idx_O2Hj)  &
        +k(694)*n(idx_HCOj)  &
        +k(701)*n(idx_OHj)  &
        +k(695)*n(idx_HNOj)  &
        +k(658)*n(idx_NHj)  &
        +k(693)*n(idx_HCNj)  &
        -k(697)*n(idx_NH2j)  &
        +k(455)*n(idx_H2j)  &
        +k(696)*n(idx_N2Hj)

    !d[O+_dot]/d[NH]
    pd(70,25) =  &
        -k(181)*n(idx_Oj)  &
        -k(698)*n(idx_Oj)

    !d[OH+_dot]/d[NH]
    pd(71,25) =  &
        -k(701)*n(idx_OHj)

    !d[CH3+_dot]/d[NH]
    pd(72,25) =  &
        -k(689)*n(idx_CH3j)

    !d[N+_dot]/d[NH]
    pd(74,25) =  &
        +k(603)*n(idx_HEj)  &
        -k(624)*n(idx_Nj)  &
        -k(146)*n(idx_Nj)

    !d[HCN+_dot]/d[NH]
    pd(75,25) =  &
        -k(693)*n(idx_HCNj)

    !d[NH+_dot]/d[NH]
    pd(76,25) =  &
        +k(80)*n(idx_Hj)  &
        +k(1027)  &
        +k(96)*n(idx_H2j)  &
        +k(179)*n(idx_COj)  &
        +k(146)*n(idx_Nj)  &
        +k(249)  &
        +k(180)*n(idx_N2j)  &
        +k(181)*n(idx_Oj)  &
        +k(178)*n(idx_CNj)  &
        -k(658)*n(idx_NHj)

    !d[H2+_dot]/d[NH]
    pd(77,25) =  &
        -k(455)*n(idx_H2j)  &
        -k(96)*n(idx_H2j)

    !d[HE+_dot]/d[NH]
    pd(78,25) =  &
        -k(603)*n(idx_HEj)

    !d[HNO+_dot]/d[NH]
    pd(79,25) =  &
        -k(695)*n(idx_HNOj)  &
        +k(699)*n(idx_O2j)

    !d[H3+_dot]/d[NH]
    pd(81,25) =  &
        -k(520)*n(idx_H3j)

    !d[H3CO+_dot]/d[NH]
    pd(82,25) =  &
        +k(691)*n(idx_H2COj)

    !d[H3O+_dot]/d[NH]
    pd(83,25) =  &
        +k(692)*n(idx_H2Oj)

    !d[HCNH+_dot]/d[NH]
    pd(84,25) =  &
        +k(689)*n(idx_CH3j)

    !d[N2H+_dot]/d[NH]
    pd(87,25) =  &
        -k(696)*n(idx_N2Hj)

    !d[O2H+_dot]/d[NH]
    pd(88,25) =  &
        -k(700)*n(idx_O2Hj)

    !d[E_dot]/d[HE]
    pd(1,26) =  &
        +k(212)  &
        +k(239)

    !d[H_dot]/d[HE]
    pd(8,26) =  &
        +k(452)*n(idx_H2j)

    !d[HE_dot]/d[HE]
    pd(26,26) =  &
        -k(239)  &
        -k(452)*n(idx_H2j)  &
        -k(1046)*n(idx_Hj)  &
        -k(212)

    !d[H+_dot]/d[HE]
    pd(55,26) =  &
        -k(1046)*n(idx_Hj)

    !d[H2+_dot]/d[HE]
    pd(77,26) =  &
        -k(452)*n(idx_H2j)

    !d[HE+_dot]/d[HE]
    pd(78,26) =  &
        +k(212)  &
        +k(239)

    !d[HEH+_dot]/d[HE]
    pd(86,26) =  &
        +k(452)*n(idx_H2j)  &
        +k(1046)*n(idx_Hj)

    !d[CH_dot]/d[HNO]
    pd(2,27) =  &
        -k(808)*n(idx_CH)

    !d[O_dot]/d[HNO]
    pd(3,27) =  &
        -k(948)*n(idx_O)  &
        +k(861)*n(idx_H)  &
        -k(950)*n(idx_O)  &
        -k(949)*n(idx_O)

    !d[HCN_dot]/d[HNO]
    pd(5,27) =  &
        +k(826)*n(idx_CN)

    !d[H2_dot]/d[HNO]
    pd(6,27) =  &
        +k(516)*n(idx_H3j)  &
        +k(862)*n(idx_H)  &
        +k(439)*n(idx_Hj)

    !d[H_dot]/d[HNO]
    pd(8,27) =  &
        +k(948)*n(idx_O)  &
        -k(862)*n(idx_H)  &
        -k(861)*n(idx_H)  &
        +k(596)*n(idx_HEj)  &
        +k(1017)  &
        +k(238)  &
        -k(863)*n(idx_H)

    !d[H2O_dot]/d[HNO]
    pd(9,27) =  &
        +k(968)*n(idx_OH)

    !d[OH_dot]/d[HNO]
    pd(10,27) =  &
        +k(863)*n(idx_H)  &
        +k(949)*n(idx_O)  &
        -k(968)*n(idx_OH)

    !d[O2_dot]/d[HNO]
    pd(11,27) =  &
        +k(950)*n(idx_O)

    !d[CH2_dot]/d[HNO]
    pd(12,27) =  &
        +k(808)*n(idx_CH)  &
        -k(765)*n(idx_CH2)

    !d[H2CO_dot]/d[HNO]
    pd(13,27) =  &
        +k(880)*n(idx_HCO)

    !d[HCO_dot]/d[HNO]
    pd(14,27) =  &
        -k(880)*n(idx_HCO)

    !d[NO_dot]/d[HNO]
    pd(17,27) =  &
        +k(949)*n(idx_O)  &
        +k(826)*n(idx_CN)  &
        +k(880)*n(idx_HCO)  &
        +k(788)*n(idx_CH3)  &
        +k(765)*n(idx_CH2)  &
        +k(597)*n(idx_HEj)  &
        +k(968)*n(idx_OH)  &
        +k(898)*n(idx_N)  &
        +k(862)*n(idx_H)  &
        +k(1017)  &
        +k(808)*n(idx_CH)  &
        +k(238)

    !d[CN_dot]/d[HNO]
    pd(18,27) =  &
        -k(826)*n(idx_CN)

    !d[CO_dot]/d[HNO]
    pd(19,27) =  &
        -k(832)*n(idx_CO)

    !d[NH2_dot]/d[HNO]
    pd(21,27) =  &
        +k(861)*n(idx_H)

    !d[CH3_dot]/d[HNO]
    pd(22,27) =  &
        -k(788)*n(idx_CH3)  &
        +k(765)*n(idx_CH2)

    !d[CH4_dot]/d[HNO]
    pd(23,27) =  &
        +k(788)*n(idx_CH3)

    !d[N_dot]/d[HNO]
    pd(24,27) =  &
        -k(898)*n(idx_N)

    !d[NH_dot]/d[HNO]
    pd(25,27) =  &
        +k(863)*n(idx_H)  &
        +k(950)*n(idx_O)  &
        +k(832)*n(idx_CO)  &
        +k(898)*n(idx_N)

    !d[HE_dot]/d[HNO]
    pd(26,27) =  &
        +k(597)*n(idx_HEj)  &
        +k(596)*n(idx_HEj)

    !d[HNO_dot]/d[HNO]
    pd(27,27) =  &
        -k(862)*n(idx_H)  &
        -k(238)  &
        -k(861)*n(idx_H)  &
        -k(597)*n(idx_HEj)  &
        -k(439)*n(idx_Hj)  &
        -k(1017)  &
        -k(788)*n(idx_CH3)  &
        -k(948)*n(idx_O)  &
        -k(880)*n(idx_HCO)  &
        -k(765)*n(idx_CH2)  &
        -k(826)*n(idx_CN)  &
        -k(516)*n(idx_H3j)  &
        -k(808)*n(idx_CH)  &
        -k(950)*n(idx_O)  &
        -k(1114)  &
        -k(898)*n(idx_N)  &
        -k(968)*n(idx_OH)  &
        -k(949)*n(idx_O)  &
        -k(596)*n(idx_HEj)  &
        -k(832)*n(idx_CO)  &
        -k(863)*n(idx_H)

    !d[CO2_dot]/d[HNO]
    pd(29,27) =  &
        +k(832)*n(idx_CO)

    !d[NO2_dot]/d[HNO]
    pd(32,27) =  &
        +k(948)*n(idx_O)

    !d[HNO_DUST_dot]/d[HNO]
    pd(48,27) =  &
        +k(1114)

    !d[H+_dot]/d[HNO]
    pd(55,27) =  &
        +k(597)*n(idx_HEj)  &
        -k(439)*n(idx_Hj)

    !d[NO+_dot]/d[HNO]
    pd(63,27) =  &
        +k(439)*n(idx_Hj)  &
        +k(596)*n(idx_HEj)

    !d[HE+_dot]/d[HNO]
    pd(78,27) =  &
        -k(597)*n(idx_HEj)  &
        -k(596)*n(idx_HEj)

    !d[H2NO+_dot]/d[HNO]
    pd(80,27) =  &
        +k(516)*n(idx_H3j)

    !d[H3+_dot]/d[HNO]
    pd(81,27) =  &
        -k(516)*n(idx_H3j)

    !d[E_dot]/d[CH3OH]
    pd(1,28) =  &
        +k(986)

    !d[CH_dot]/d[CH3OH]
    pd(2,28) =  &
        +k(314)*n(idx_Cj)

    !d[H2_dot]/d[CH3OH]
    pd(6,28) =  &
        +k(222)  &
        +k(505)*n(idx_H3j)  &
        +k(985)  &
        +2.d0*k(431)*n(idx_Hj)  &
        +k(430)*n(idx_Hj)

    !d[H_dot]/d[CH3OH]
    pd(8,28) =  &
        +k(715)*n(idx_O2j)  &
        +k(612)*n(idx_Nj)  &
        +k(610)*n(idx_Nj)  &
        +k(613)*n(idx_Nj)  &
        +k(986)

    !d[H2O_dot]/d[CH3OH]
    pd(9,28) =  &
        +k(429)*n(idx_Hj)  &
        +k(703)*n(idx_Oj)  &
        +k(505)*n(idx_H3j)

    !d[OH_dot]/d[CH3OH]
    pd(10,28) =  &
        +k(223)  &
        +k(987)  &
        +k(704)*n(idx_Oj)  &
        +k(568)*n(idx_HEj)

    !d[O2_dot]/d[CH3OH]
    pd(11,28) =  &
        +k(715)*n(idx_O2j)

    !d[CH2_dot]/d[CH3OH]
    pd(12,28) =  &
        +k(341)*n(idx_CHj)

    !d[H2CO_dot]/d[CH3OH]
    pd(13,28) =  &
        +k(222)  &
        +k(985)  &
        +k(340)*n(idx_CHj)

    !d[HCO_dot]/d[CH3OH]
    pd(14,28) =  &
        +k(315)*n(idx_Cj)

    !d[NO_dot]/d[CH3OH]
    pd(17,28) =  &
        +k(613)*n(idx_Nj)

    !d[CH3_dot]/d[CH3OH]
    pd(22,28) =  &
        +k(223)  &
        +k(612)*n(idx_Nj)  &
        +k(987)  &
        +k(567)*n(idx_HEj)

    !d[CH4_dot]/d[CH3OH]
    pd(23,28) =  &
        +k(382)*n(idx_CH3j)

    !d[NH_dot]/d[CH3OH]
    pd(25,28) =  &
        +k(610)*n(idx_Nj)  &
        +k(611)*n(idx_Nj)

    !d[HE_dot]/d[CH3OH]
    pd(26,28) =  &
        +k(567)*n(idx_HEj)  &
        +k(568)*n(idx_HEj)

    !d[CH3OH_dot]/d[CH3OH]
    pd(28,28) =  &
        -k(610)*n(idx_Nj)  &
        -k(715)*n(idx_O2j)  &
        -k(222)  &
        -k(1064)  &
        -k(315)*n(idx_Cj)  &
        -k(382)*n(idx_CH3j)  &
        -k(612)*n(idx_Nj)  &
        -k(986)  &
        -k(987)  &
        -k(505)*n(idx_H3j)  &
        -k(568)*n(idx_HEj)  &
        -k(613)*n(idx_Nj)  &
        -k(314)*n(idx_Cj)  &
        -k(985)  &
        -k(223)  &
        -k(340)*n(idx_CHj)  &
        -k(704)*n(idx_Oj)  &
        -k(430)*n(idx_Hj)  &
        -k(341)*n(idx_CHj)  &
        -k(429)*n(idx_Hj)  &
        -k(703)*n(idx_Oj)  &
        -k(431)*n(idx_Hj)  &
        -k(611)*n(idx_Nj)  &
        -k(567)*n(idx_HEj)

    !d[CH3OH_DUST_dot]/d[CH3OH]
    pd(35,28) =  &
        +k(1064)

    !d[HCO+_dot]/d[CH3OH]
    pd(54,28) =  &
        +k(431)*n(idx_Hj)

    !d[H+_dot]/d[CH3OH]
    pd(55,28) =  &
        -k(430)*n(idx_Hj)  &
        -k(429)*n(idx_Hj)  &
        -k(431)*n(idx_Hj)

    !d[C+_dot]/d[CH3OH]
    pd(57,28) =  &
        -k(314)*n(idx_Cj)  &
        -k(315)*n(idx_Cj)

    !d[CH+_dot]/d[CH3OH]
    pd(59,28) =  &
        -k(340)*n(idx_CHj)  &
        -k(341)*n(idx_CHj)

    !d[H2CO+_dot]/d[CH3OH]
    pd(60,28) =  &
        +k(703)*n(idx_Oj)  &
        +k(610)*n(idx_Nj)

    !d[NO+_dot]/d[CH3OH]
    pd(63,28) =  &
        +k(612)*n(idx_Nj)

    !d[O2+_dot]/d[CH3OH]
    pd(67,28) =  &
        -k(715)*n(idx_O2j)

    !d[O+_dot]/d[CH3OH]
    pd(70,28) =  &
        -k(704)*n(idx_Oj)  &
        -k(703)*n(idx_Oj)

    !d[OH+_dot]/d[CH3OH]
    pd(71,28) =  &
        +k(567)*n(idx_HEj)

    !d[CH3+_dot]/d[CH3OH]
    pd(72,28) =  &
        +k(429)*n(idx_Hj)  &
        +k(315)*n(idx_Cj)  &
        +k(613)*n(idx_Nj)  &
        +k(340)*n(idx_CHj)  &
        -k(382)*n(idx_CH3j)  &
        +k(505)*n(idx_H3j)  &
        +k(568)*n(idx_HEj)

    !d[N+_dot]/d[CH3OH]
    pd(74,28) =  &
        -k(610)*n(idx_Nj)  &
        -k(612)*n(idx_Nj)  &
        -k(611)*n(idx_Nj)  &
        -k(613)*n(idx_Nj)

    !d[HE+_dot]/d[CH3OH]
    pd(78,28) =  &
        -k(568)*n(idx_HEj)  &
        -k(567)*n(idx_HEj)

    !d[H3+_dot]/d[CH3OH]
    pd(81,28) =  &
        -k(505)*n(idx_H3j)

    !d[H3CO+_dot]/d[CH3OH]
    pd(82,28) =  &
        +k(314)*n(idx_Cj)  &
        +k(704)*n(idx_Oj)  &
        +k(430)*n(idx_Hj)  &
        +k(986)  &
        +k(341)*n(idx_CHj)  &
        +k(715)*n(idx_O2j)  &
        +k(611)*n(idx_Nj)  &
        +k(382)*n(idx_CH3j)

    !d[CH_dot]/d[CO2]
    pd(2,29) =  &
        -k(805)*n(idx_CH)

    !d[O_dot]/d[CO2]
    pd(3,29) =  &
        +k(576)*n(idx_HEj)  &
        +k(727)*n(idx_OHj)  &
        +k(998)  &
        -k(939)*n(idx_O)  &
        +k(433)*n(idx_Hj)  &
        +k(227)

    !d[H2_dot]/d[CO2]
    pd(6,29) =  &
        +k(508)*n(idx_H3j)

    !d[C_dot]/d[CO2]
    pd(7,29) =  &
        +k(578)*n(idx_HEj)

    !d[H_dot]/d[CO2]
    pd(8,29) =  &
        +k(446)*n(idx_H2j)  &
        -k(852)*n(idx_H)

    !d[OH_dot]/d[CO2]
    pd(10,29) =  &
        +k(852)*n(idx_H)

    !d[O2_dot]/d[CO2]
    pd(11,29) =  &
        +k(579)*n(idx_HEj)  &
        +k(939)*n(idx_O)  &
        +k(716)*n(idx_O2Hj)

    !d[HCO_dot]/d[CO2]
    pd(14,29) =  &
        +k(645)*n(idx_NHj)  &
        +k(805)*n(idx_CH)

    !d[NO_dot]/d[CO2]
    pd(17,29) =  &
        +k(617)*n(idx_Nj)  &
        +k(563)*n(idx_HNOj)  &
        +k(893)*n(idx_N)

    !d[CN_dot]/d[CO2]
    pd(18,29) =  &
        +k(535)*n(idx_HCNj)

    !d[CO_dot]/d[CO2]
    pd(19,29) =  &
        +k(998)  &
        +k(893)*n(idx_N)  &
        +k(316)*n(idx_Cj)  &
        +k(707)*n(idx_Oj)  &
        +k(342)*n(idx_CHj)  &
        +k(360)*n(idx_CH2j)  &
        +k(939)*n(idx_O)  &
        +k(852)*n(idx_H)  &
        +k(644)*n(idx_NHj)  &
        +k(227)  &
        +k(577)*n(idx_HEj)  &
        +k(805)*n(idx_CH)

    !d[N2_dot]/d[CO2]
    pd(20,29) =  &
        +k(632)*n(idx_N2Hj)

    !d[CH3_dot]/d[CO2]
    pd(22,29) =  &
        +k(389)*n(idx_CH4j)

    !d[N_dot]/d[CO2]
    pd(24,29) =  &
        -k(893)*n(idx_N)  &
        +k(643)*n(idx_NHj)

    !d[HE_dot]/d[CO2]
    pd(26,29) =  &
        +k(579)*n(idx_HEj)  &
        +k(578)*n(idx_HEj)  &
        +k(577)*n(idx_HEj)  &
        +k(576)*n(idx_HEj)

    !d[CO2_dot]/d[CO2]
    pd(29,29) =  &
        -k(707)*n(idx_Oj)  &
        -k(577)*n(idx_HEj)  &
        -k(716)*n(idx_O2Hj)  &
        -k(508)*n(idx_H3j)  &
        -k(535)*n(idx_HCNj)  &
        -k(643)*n(idx_NHj)  &
        -k(433)*n(idx_Hj)  &
        -k(805)*n(idx_CH)  &
        -k(939)*n(idx_O)  &
        -k(727)*n(idx_OHj)  &
        -k(617)*n(idx_Nj)  &
        -k(644)*n(idx_NHj)  &
        -k(578)*n(idx_HEj)  &
        -k(563)*n(idx_HNOj)  &
        -k(576)*n(idx_HEj)  &
        -k(342)*n(idx_CHj)  &
        -k(998)  &
        -k(446)*n(idx_H2j)  &
        -k(1078)  &
        -k(645)*n(idx_NHj)  &
        -k(389)*n(idx_CH4j)  &
        -k(360)*n(idx_CH2j)  &
        -k(579)*n(idx_HEj)  &
        -k(852)*n(idx_H)  &
        -k(893)*n(idx_N)  &
        -k(632)*n(idx_N2Hj)  &
        -k(316)*n(idx_Cj)  &
        -k(227)

    !d[CO2_DUST_dot]/d[CO2]
    pd(42,29) =  &
        +k(1078)

    !d[HCO+_dot]/d[CO2]
    pd(54,29) =  &
        +k(342)*n(idx_CHj)  &
        +k(433)*n(idx_Hj)

    !d[H+_dot]/d[CO2]
    pd(55,29) =  &
        -k(433)*n(idx_Hj)

    !d[C+_dot]/d[CO2]
    pd(57,29) =  &
        +k(579)*n(idx_HEj)  &
        -k(316)*n(idx_Cj)

    !d[CH2+_dot]/d[CO2]
    pd(58,29) =  &
        -k(360)*n(idx_CH2j)

    !d[CH+_dot]/d[CO2]
    pd(59,29) =  &
        -k(342)*n(idx_CHj)

    !d[H2CO+_dot]/d[CO2]
    pd(60,29) =  &
        +k(360)*n(idx_CH2j)

    !d[NO+_dot]/d[CO2]
    pd(63,29) =  &
        +k(645)*n(idx_NHj)

    !d[CO+_dot]/d[CO2]
    pd(65,29) =  &
        +k(316)*n(idx_Cj)  &
        +k(617)*n(idx_Nj)  &
        +k(576)*n(idx_HEj)

    !d[O2+_dot]/d[CO2]
    pd(67,29) =  &
        +k(578)*n(idx_HEj)  &
        +k(707)*n(idx_Oj)

    !d[O+_dot]/d[CO2]
    pd(70,29) =  &
        -k(707)*n(idx_Oj)  &
        +k(577)*n(idx_HEj)

    !d[OH+_dot]/d[CO2]
    pd(71,29) =  &
        -k(727)*n(idx_OHj)

    !d[CH4+_dot]/d[CO2]
    pd(73,29) =  &
        -k(389)*n(idx_CH4j)

    !d[N+_dot]/d[CO2]
    pd(74,29) =  &
        -k(617)*n(idx_Nj)

    !d[HCN+_dot]/d[CO2]
    pd(75,29) =  &
        -k(535)*n(idx_HCNj)

    !d[NH+_dot]/d[CO2]
    pd(76,29) =  &
        -k(643)*n(idx_NHj)  &
        -k(644)*n(idx_NHj)  &
        -k(645)*n(idx_NHj)

    !d[H2+_dot]/d[CO2]
    pd(77,29) =  &
        -k(446)*n(idx_H2j)

    !d[HE+_dot]/d[CO2]
    pd(78,29) =  &
        -k(579)*n(idx_HEj)  &
        -k(578)*n(idx_HEj)  &
        -k(577)*n(idx_HEj)  &
        -k(576)*n(idx_HEj)

    !d[HNO+_dot]/d[CO2]
    pd(79,29) =  &
        +k(644)*n(idx_NHj)  &
        -k(563)*n(idx_HNOj)

    !d[H3+_dot]/d[CO2]
    pd(81,29) =  &
        -k(508)*n(idx_H3j)

    !d[HCO2+_dot]/d[CO2]
    pd(85,29) =  &
        +k(563)*n(idx_HNOj)  &
        +k(535)*n(idx_HCNj)  &
        +k(446)*n(idx_H2j)  &
        +k(727)*n(idx_OHj)  &
        +k(716)*n(idx_O2Hj)  &
        +k(643)*n(idx_NHj)  &
        +k(632)*n(idx_N2Hj)  &
        +k(389)*n(idx_CH4j)  &
        +k(508)*n(idx_H3j)

    !d[N2H+_dot]/d[CO2]
    pd(87,29) =  &
        -k(632)*n(idx_N2Hj)

    !d[O2H+_dot]/d[CO2]
    pd(88,29) =  &
        -k(716)*n(idx_O2Hj)

    !d[O_dot]/d[H2CN]
    pd(3,30) =  &
        -k(940)*n(idx_O)

    !d[HCN_dot]/d[H2CN]
    pd(5,30) =  &
        +k(894)*n(idx_N)  &
        +k(229)  &
        +k(1001)  &
        +k(854)*n(idx_H)

    !d[H2_dot]/d[H2CN]
    pd(6,30) =  &
        +k(940)*n(idx_O)  &
        +k(854)*n(idx_H)

    !d[H_dot]/d[H2CN]
    pd(8,30) =  &
        +k(1001)  &
        +k(229)  &
        -k(854)*n(idx_H)

    !d[N_dot]/d[H2CN]
    pd(24,30) =  &
        -k(894)*n(idx_N)

    !d[NH_dot]/d[H2CN]
    pd(25,30) =  &
        +k(894)*n(idx_N)

    !d[H2CN_dot]/d[H2CN]
    pd(30,30) =  &
        -k(854)*n(idx_H)  &
        -k(1001)  &
        -k(229)  &
        -k(940)*n(idx_O)  &
        -k(894)*n(idx_N)  &
        -k(1118)

    !d[OCN_dot]/d[H2CN]
    pd(34,30) =  &
        +k(940)*n(idx_O)

    !d[H2CN_DUST_dot]/d[H2CN]
    pd(50,30) =  &
        +k(1118)

    !d[HNC_dot]/d[HNCO]
    pd(4,31) =  &
        +k(885)*n(idx_C)

    !d[C_dot]/d[HNCO]
    pd(7,31) =  &
        -k(885)*n(idx_C)

    !d[CO_dot]/d[HNCO]
    pd(19,31) =  &
        +k(885)*n(idx_C)  &
        +k(237)  &
        +k(438)*n(idx_Hj)  &
        +k(1016)

    !d[NH_dot]/d[HNCO]
    pd(25,31) =  &
        +k(1016)  &
        +k(237)

    !d[HNCO_dot]/d[HNCO]
    pd(31,31) =  &
        -k(438)*n(idx_Hj)  &
        -k(237)  &
        -k(1065)  &
        -k(885)*n(idx_C)  &
        -k(1016)

    !d[HNCO_DUST_dot]/d[HNCO]
    pd(36,31) =  &
        +k(1065)

    !d[H+_dot]/d[HNCO]
    pd(55,31) =  &
        -k(438)*n(idx_Hj)

    !d[NH2+_dot]/d[HNCO]
    pd(69,31) =  &
        +k(438)*n(idx_Hj)

    !d[O_dot]/d[NO2]
    pd(3,32) =  &
        +k(1028)  &
        +2.d0*k(900)*n(idx_N)  &
        +k(250)  &
        -k(955)*n(idx_O)

    !d[H2_dot]/d[NO2]
    pd(6,32) =  &
        +k(521)*n(idx_H3j)

    !d[H_dot]/d[NO2]
    pd(8,32) =  &
        -k(867)*n(idx_H)

    !d[OH_dot]/d[NO2]
    pd(10,32) =  &
        +k(440)*n(idx_Hj)  &
        +k(867)*n(idx_H)  &
        +k(521)*n(idx_H3j)

    !d[O2_dot]/d[NO2]
    pd(11,32) =  &
        +k(902)*n(idx_N)  &
        +k(955)*n(idx_O)  &
        +k(713)*n(idx_Oj)

    !d[CH2_dot]/d[NO2]
    pd(12,32) =  &
        -k(767)*n(idx_CH2)

    !d[H2CO_dot]/d[NO2]
    pd(13,32) =  &
        +k(767)*n(idx_CH2)  &
        +k(791)*n(idx_CH3)

    !d[NO_dot]/d[NO2]
    pd(17,32) =  &
        +k(921)*n(idx_NH)  &
        +k(1028)  &
        +2.d0*k(901)*n(idx_N)  &
        +k(250)  &
        +k(827)*n(idx_CN)  &
        +k(955)*n(idx_O)  &
        +k(767)*n(idx_CH2)  &
        +k(833)*n(idx_CO)  &
        +k(867)*n(idx_H)

    !d[CN_dot]/d[NO2]
    pd(18,32) =  &
        -k(827)*n(idx_CN)

    !d[CO_dot]/d[NO2]
    pd(19,32) =  &
        -k(833)*n(idx_CO)

    !d[N2_dot]/d[NO2]
    pd(20,32) =  &
        +k(902)*n(idx_N)  &
        +k(900)*n(idx_N)

    !d[CH3_dot]/d[NO2]
    pd(22,32) =  &
        -k(791)*n(idx_CH3)

    !d[N_dot]/d[NO2]
    pd(24,32) =  &
        -k(901)*n(idx_N)  &
        -k(900)*n(idx_N)  &
        -k(902)*n(idx_N)

    !d[NH_dot]/d[NO2]
    pd(25,32) =  &
        -k(921)*n(idx_NH)

    !d[HNO_dot]/d[NO2]
    pd(27,32) =  &
        +k(921)*n(idx_NH)  &
        +k(791)*n(idx_CH3)

    !d[CO2_dot]/d[NO2]
    pd(29,32) =  &
        +k(833)*n(idx_CO)

    !d[NO2_dot]/d[NO2]
    pd(32,32) =  &
        -k(767)*n(idx_CH2)  &
        -k(440)*n(idx_Hj)  &
        -k(901)*n(idx_N)  &
        -k(900)*n(idx_N)  &
        -k(1113)  &
        -k(791)*n(idx_CH3)  &
        -k(1028)  &
        -k(921)*n(idx_NH)  &
        -k(250)  &
        -k(521)*n(idx_H3j)  &
        -k(833)*n(idx_CO)  &
        -k(867)*n(idx_H)  &
        -k(827)*n(idx_CN)  &
        -k(955)*n(idx_O)  &
        -k(713)*n(idx_Oj)  &
        -k(902)*n(idx_N)

    !d[OCN_dot]/d[NO2]
    pd(34,32) =  &
        +k(827)*n(idx_CN)

    !d[NO2_DUST_dot]/d[NO2]
    pd(47,32) =  &
        +k(1113)

    !d[H+_dot]/d[NO2]
    pd(55,32) =  &
        -k(440)*n(idx_Hj)

    !d[NO+_dot]/d[NO2]
    pd(63,32) =  &
        +k(713)*n(idx_Oj)  &
        +k(440)*n(idx_Hj)  &
        +k(521)*n(idx_H3j)

    !d[O+_dot]/d[NO2]
    pd(70,32) =  &
        -k(713)*n(idx_Oj)

    !d[H3+_dot]/d[NO2]
    pd(81,32) =  &
        -k(521)*n(idx_H3j)

    !d[CH_dot]/d[O2H]
    pd(2,33) =  &
        -k(819)*n(idx_CH)  &
        -k(820)*n(idx_CH)

    !d[O_dot]/d[O2H]
    pd(3,33) =  &
        -k(957)*n(idx_O)  &
        +k(871)*n(idx_H)  &
        +k(1035)

    !d[H2_dot]/d[O2H]
    pd(6,33) =  &
        +k(872)*n(idx_H)

    !d[H_dot]/d[O2H]
    pd(8,33) =  &
        +k(255)  &
        -k(872)*n(idx_H)  &
        -k(873)*n(idx_H)  &
        +k(1034)  &
        -k(871)*n(idx_H)

    !d[H2O_dot]/d[O2H]
    pd(9,33) =  &
        +k(871)*n(idx_H)  &
        +k(971)*n(idx_OH)

    !d[OH_dot]/d[O2H]
    pd(10,33) =  &
        -k(971)*n(idx_OH)  &
        +k(819)*n(idx_CH)  &
        +k(1035)  &
        +k(957)*n(idx_O)  &
        +k(835)*n(idx_CO)  &
        +2.d0*k(873)*n(idx_H)

    !d[O2_dot]/d[O2H]
    pd(11,33) =  &
        +k(255)  &
        +k(796)*n(idx_CH3)  &
        +k(820)*n(idx_CH)  &
        +k(1034)  &
        +k(872)*n(idx_H)  &
        +k(905)*n(idx_N)  &
        +k(971)*n(idx_OH)  &
        +k(884)*n(idx_HCO)  &
        +k(957)*n(idx_O)

    !d[CH2_dot]/d[O2H]
    pd(12,33) =  &
        +k(820)*n(idx_CH)

    !d[H2CO_dot]/d[O2H]
    pd(13,33) =  &
        +k(884)*n(idx_HCO)

    !d[HCO_dot]/d[O2H]
    pd(14,33) =  &
        +k(819)*n(idx_CH)  &
        -k(884)*n(idx_HCO)

    !d[CO_dot]/d[O2H]
    pd(19,33) =  &
        -k(835)*n(idx_CO)

    !d[CH3_dot]/d[O2H]
    pd(22,33) =  &
        -k(796)*n(idx_CH3)

    !d[CH4_dot]/d[O2H]
    pd(23,33) =  &
        +k(796)*n(idx_CH3)

    !d[N_dot]/d[O2H]
    pd(24,33) =  &
        -k(905)*n(idx_N)

    !d[NH_dot]/d[O2H]
    pd(25,33) =  &
        +k(905)*n(idx_N)

    !d[CO2_dot]/d[O2H]
    pd(29,33) =  &
        +k(835)*n(idx_CO)

    !d[O2H_dot]/d[O2H]
    pd(33,33) =  &
        -k(872)*n(idx_H)  &
        -k(957)*n(idx_O)  &
        -k(905)*n(idx_N)  &
        -k(1123)  &
        -k(971)*n(idx_OH)  &
        -k(796)*n(idx_CH3)  &
        -k(1034)  &
        -k(1035)  &
        -k(819)*n(idx_CH)  &
        -k(873)*n(idx_H)  &
        -k(255)  &
        -k(835)*n(idx_CO)  &
        -k(820)*n(idx_CH)  &
        -k(884)*n(idx_HCO)  &
        -k(871)*n(idx_H)

    !d[O2H_DUST_dot]/d[O2H]
    pd(49,33) =  &
        +k(1123)

    !d[O_dot]/d[OCN]
    pd(3,34) =  &
        +k(607)*n(idx_HEj)  &
        -k(958)*n(idx_O)  &
        +k(257)  &
        +k(1036)  &
        -k(959)*n(idx_O)  &
        +k(874)*n(idx_H)

    !d[HCN_dot]/d[OCN]
    pd(5,34) =  &
        +k(874)*n(idx_H)

    !d[C_dot]/d[OCN]
    pd(7,34) =  &
        -k(757)*n(idx_C)

    !d[H_dot]/d[OCN]
    pd(8,34) =  &
        -k(875)*n(idx_H)  &
        -k(876)*n(idx_H)  &
        -k(874)*n(idx_H)

    !d[OH_dot]/d[OCN]
    pd(10,34) =  &
        +k(876)*n(idx_H)

    !d[O2_dot]/d[OCN]
    pd(11,34) =  &
        -k(935)*n(idx_O2)  &
        -k(934)*n(idx_O2)  &
        +k(959)*n(idx_O)

    !d[NO_dot]/d[OCN]
    pd(17,34) =  &
        +k(958)*n(idx_O)  &
        -k(933)*n(idx_NO)  &
        +k(934)*n(idx_O2)

    !d[CN_dot]/d[OCN]
    pd(18,34) =  &
        +k(257)  &
        +k(1036)  &
        +k(757)*n(idx_C)  &
        +k(608)*n(idx_HEj)  &
        +k(959)*n(idx_O)  &
        +k(876)*n(idx_H)  &
        +k(327)*n(idx_Cj)

    !d[CO_dot]/d[OCN]
    pd(19,34) =  &
        +k(875)*n(idx_H)  &
        +k(958)*n(idx_O)  &
        +k(935)*n(idx_O2)  &
        +k(757)*n(idx_C)

    !d[N2_dot]/d[OCN]
    pd(20,34) =  &
        +k(933)*n(idx_NO)

    !d[NH_dot]/d[OCN]
    pd(25,34) =  &
        +k(875)*n(idx_H)

    !d[HE_dot]/d[OCN]
    pd(26,34) =  &
        +k(608)*n(idx_HEj)  &
        +k(607)*n(idx_HEj)

    !d[CO2_dot]/d[OCN]
    pd(29,34) =  &
        +k(933)*n(idx_NO)  &
        +k(934)*n(idx_O2)

    !d[NO2_dot]/d[OCN]
    pd(32,34) =  &
        +k(935)*n(idx_O2)

    !d[OCN_dot]/d[OCN]
    pd(34,34) =  &
        -k(327)*n(idx_Cj)  &
        -k(935)*n(idx_O2)  &
        -k(958)*n(idx_O)  &
        -k(757)*n(idx_C)  &
        -k(933)*n(idx_NO)  &
        -k(607)*n(idx_HEj)  &
        -k(876)*n(idx_H)  &
        -k(959)*n(idx_O)  &
        -k(1066)  &
        -k(874)*n(idx_H)  &
        -k(257)  &
        -k(875)*n(idx_H)  &
        -k(934)*n(idx_O2)  &
        -k(608)*n(idx_HEj)  &
        -k(1036)

    !d[HNCO_DUST_dot]/d[OCN]
    pd(36,34) =  &
        +k(1066)

    !d[C+_dot]/d[OCN]
    pd(57,34) =  &
        -k(327)*n(idx_Cj)

    !d[CN+_dot]/d[OCN]
    pd(64,34) =  &
        +k(607)*n(idx_HEj)

    !d[CO+_dot]/d[OCN]
    pd(65,34) =  &
        +k(327)*n(idx_Cj)

    !d[O+_dot]/d[OCN]
    pd(70,34) =  &
        +k(608)*n(idx_HEj)

    !d[HE+_dot]/d[OCN]
    pd(78,34) =  &
        -k(607)*n(idx_HEj)  &
        -k(608)*n(idx_HEj)

    !d[CH3OH_dot]/d[CH3OH_DUST]
    pd(28,35) =  &
        +k(1140)

    !d[CH3OH_DUST_dot]/d[CH3OH_DUST]
    pd(35,35) =  &
        -k(1140)

    !d[HNCO_dot]/d[HNCO_DUST]
    pd(31,36) =  &
        +k(1142)

    !d[HNCO_DUST_dot]/d[HNCO_DUST]
    pd(36,36) =  &
        -k(1142)

    !d[H2CO_dot]/d[H2CO_DUST]
    pd(13,37) =  &
        +k(1137)

    !d[H2CO_DUST_dot]/d[H2CO_DUST]
    pd(37,37) =  &
        -k(1137)

    !d[CH4_dot]/d[CH4_DUST]
    pd(23,38) =  &
        +k(1127)

    !d[CH4_DUST_dot]/d[CH4_DUST]
    pd(38,38) =  &
        -k(1127)

    !d[CO_dot]/d[CO_DUST]
    pd(19,39) =  &
        +k(1133)

    !d[CO_DUST_dot]/d[CO_DUST]
    pd(39,39) =  &
        -k(1133)

    !d[H2O_dot]/d[H2O_DUST]
    pd(9,40) =  &
        +k(1129)

    !d[H2O_DUST_dot]/d[H2O_DUST]
    pd(40,40) =  &
        -k(1129)

    !d[NO_dot]/d[NO_DUST]
    pd(17,41) =  &
        +k(1136)

    !d[NO_DUST_dot]/d[NO_DUST]
    pd(41,41) =  &
        -k(1136)

    !d[CO2_dot]/d[CO2_DUST]
    pd(29,42) =  &
        +k(1143)

    !d[CO2_DUST_dot]/d[CO2_DUST]
    pd(42,42) =  &
        -k(1143)

    !d[N2_dot]/d[N2_DUST]
    pd(20,43) =  &
        +k(1134)

    !d[N2_DUST_dot]/d[N2_DUST]
    pd(43,43) =  &
        -k(1134)

    !d[HCN_dot]/d[HCN_DUST]
    pd(5,44) =  &
        +k(1131)

    !d[HCN_DUST_dot]/d[HCN_DUST]
    pd(44,44) =  &
        -k(1131)

    !d[NH3_dot]/d[NH3_DUST]
    pd(16,45) =  &
        +k(1128)

    !d[NH3_DUST_dot]/d[NH3_DUST]
    pd(45,45) =  &
        -k(1128)

    !d[O2_dot]/d[O2_DUST]
    pd(11,46) =  &
        +k(1139)

    !d[O2_DUST_dot]/d[O2_DUST]
    pd(46,46) =  &
        -k(1139)

    !d[NO2_dot]/d[NO2_DUST]
    pd(32,47) =  &
        +k(1144)

    !d[NO2_DUST_dot]/d[NO2_DUST]
    pd(47,47) =  &
        -k(1144)

    !d[HNO_dot]/d[HNO_DUST]
    pd(27,48) =  &
        +k(1138)

    !d[HNO_DUST_dot]/d[HNO_DUST]
    pd(48,48) =  &
        -k(1138)

    !d[O2H_dot]/d[O2H_DUST]
    pd(33,49) =  &
        +k(1141)

    !d[O2H_DUST_dot]/d[O2H_DUST]
    pd(49,49) =  &
        -k(1141)

    !d[H2CN_dot]/d[H2CN_DUST]
    pd(30,50) =  &
        +k(1135)

    !d[H2CN_DUST_dot]/d[H2CN_DUST]
    pd(50,50) =  &
        -k(1135)

    !d[MG_dot]/d[MG_DUST]
    pd(15,51) =  &
        +k(1130)

    !d[MG_DUST_dot]/d[MG_DUST]
    pd(51,51) =  &
        -k(1130)

    !d[HNC_dot]/d[HNC_DUST]
    pd(4,52) =  &
        +k(1132)

    !d[HNC_DUST_dot]/d[HNC_DUST]
    pd(52,52) =  &
        -k(1132)

    !d[E_dot]/d[HCO+]
    pd(1,54) =  &
        -k(295)*n(idx_E)

    !d[CH_dot]/d[HCO+]
    pd(2,54) =  &
        -k(408)*n(idx_CH)

    !d[HNC_dot]/d[HCO+]
    pd(4,54) =  &
        -k(559)*n(idx_HNC)

    !d[HCN_dot]/d[HCO+]
    pd(5,54) =  &
        -k(544)*n(idx_HCN)

    !d[C_dot]/d[HCO+]
    pd(7,54) =  &
        -k(332)*n(idx_C)

    !d[H_dot]/d[HCO+]
    pd(8,54) =  &
        +k(295)*n(idx_E)  &
        +k(1012)  &
        +k(742)*n(idx_OH)

    !d[H2O_dot]/d[HCO+]
    pd(9,54) =  &
        -k(496)*n(idx_H2O)

    !d[OH_dot]/d[HCO+]
    pd(10,54) =  &
        -k(741)*n(idx_OH)  &
        -k(742)*n(idx_OH)

    !d[CH2_dot]/d[HCO+]
    pd(12,54) =  &
        -k(373)*n(idx_CH2)

    !d[H2CO_dot]/d[HCO+]
    pd(13,54) =  &
        -k(550)*n(idx_H2CO)

    !d[HCO_dot]/d[HCO+]
    pd(14,54) =  &
        -k(551)*n(idx_HCO)  &
        +k(131)*n(idx_MG)

    !d[MG_dot]/d[HCO+]
    pd(15,54) =  &
        -k(131)*n(idx_MG)

    !d[CO_dot]/d[HCO+]
    pd(19,54) =  &
        +k(741)*n(idx_OH)  &
        +k(550)*n(idx_H2CO)  &
        +k(496)*n(idx_H2O)  &
        +k(408)*n(idx_CH)  &
        +k(332)*n(idx_C)  &
        +k(559)*n(idx_HNC)  &
        +k(544)*n(idx_HCN)  &
        +k(295)*n(idx_E)  &
        +k(373)*n(idx_CH2)  &
        +k(694)*n(idx_NH)  &
        +k(551)*n(idx_HCO)  &
        +k(682)*n(idx_NH2)

    !d[NH2_dot]/d[HCO+]
    pd(21,54) =  &
        -k(682)*n(idx_NH2)

    !d[NH_dot]/d[HCO+]
    pd(25,54) =  &
        -k(694)*n(idx_NH)

    !d[H2CO_DUST_dot]/d[HCO+]
    pd(37,54) =  &
        +k(1101)

    !d[HCO+_dot]/d[HCO+]
    pd(54,54) =  &
        -k(682)*n(idx_NH2)  &
        -k(694)*n(idx_NH)  &
        -k(332)*n(idx_C)  &
        -k(295)*n(idx_E)  &
        -k(559)*n(idx_HNC)  &
        -k(551)*n(idx_HCO)  &
        -k(373)*n(idx_CH2)  &
        -k(741)*n(idx_OH)  &
        -k(742)*n(idx_OH)  &
        -k(131)*n(idx_MG)  &
        -k(1012)  &
        -k(544)*n(idx_HCN)  &
        -k(408)*n(idx_CH)  &
        -k(496)*n(idx_H2O)  &
        -k(1101)  &
        -k(550)*n(idx_H2CO)

    !d[CH2+_dot]/d[HCO+]
    pd(58,54) =  &
        +k(408)*n(idx_CH)

    !d[CH+_dot]/d[HCO+]
    pd(59,54) =  &
        +k(332)*n(idx_C)

    !d[H2CO+_dot]/d[HCO+]
    pd(60,54) =  &
        +k(551)*n(idx_HCO)

    !d[MG+_dot]/d[HCO+]
    pd(61,54) =  &
        +k(131)*n(idx_MG)

    !d[NH3+_dot]/d[HCO+]
    pd(62,54) =  &
        +k(682)*n(idx_NH2)

    !d[CO+_dot]/d[HCO+]
    pd(65,54) =  &
        +k(1012)

    !d[H2O+_dot]/d[HCO+]
    pd(68,54) =  &
        +k(741)*n(idx_OH)

    !d[NH2+_dot]/d[HCO+]
    pd(69,54) =  &
        +k(694)*n(idx_NH)

    !d[CH3+_dot]/d[HCO+]
    pd(72,54) =  &
        +k(373)*n(idx_CH2)

    !d[H3CO+_dot]/d[HCO+]
    pd(82,54) =  &
        +k(550)*n(idx_H2CO)

    !d[H3O+_dot]/d[HCO+]
    pd(83,54) =  &
        +k(496)*n(idx_H2O)

    !d[HCNH+_dot]/d[HCO+]
    pd(84,54) =  &
        +k(544)*n(idx_HCN)  &
        +k(559)*n(idx_HNC)

    !d[HCO2+_dot]/d[HCO+]
    pd(85,54) =  &
        +k(742)*n(idx_OH)

    !d[E_dot]/d[H+]
    pd(1,55) =  &
        -k(1058)*n(idx_E)

    !d[CH_dot]/d[H+]
    pd(2,55) =  &
        -k(72)*n(idx_CH)

    !d[O_dot]/d[H+]
    pd(3,55) =  &
        +k(433)*n(idx_CO2)  &
        -k(83)*n(idx_O)

    !d[HNC_dot]/d[H+]
    pd(4,55) =  &
        -k(2)*n(idx_HNC)

    !d[HCN_dot]/d[H+]
    pd(5,55) =  &
        +k(2)*n(idx_HNC)  &
        -k(75)*n(idx_HCN)

    !d[H2_dot]/d[H+]
    pd(6,55) =  &
        +k(432)*n(idx_CH4)  &
        +k(436)*n(idx_HCO)  &
        +k(434)*n(idx_H2CO)  &
        +k(439)*n(idx_HNO)  &
        +k(435)*n(idx_H2CO)  &
        +k(430)*n(idx_CH3OH)  &
        +k(428)*n(idx_CH2)  &
        +2.d0*k(431)*n(idx_CH3OH)

    !d[H_dot]/d[H+]
    pd(8,55) =  &
        -k(1045)*n(idx_H)  &
        +k(79)*n(idx_NH3)  &
        +k(70)*n(idx_CH3)  &
        +k(434)*n(idx_H2CO)  &
        +k(80)*n(idx_NH)  &
        +k(69)*n(idx_CH2)  &
        +k(81)*n(idx_NO)  &
        +k(84)*n(idx_OH)  &
        +k(78)*n(idx_NH2)  &
        +k(83)*n(idx_O)  &
        +k(71)*n(idx_CH4)  &
        +k(73)*n(idx_H2CO)  &
        +k(82)*n(idx_O2)  &
        +k(72)*n(idx_CH)  &
        +k(75)*n(idx_HCN)  &
        +k(74)*n(idx_H2O)  &
        +k(76)*n(idx_HCO)  &
        +k(77)*n(idx_MG)  &
        +k(1058)*n(idx_E)

    !d[H2O_dot]/d[H+]
    pd(9,55) =  &
        +k(429)*n(idx_CH3OH)  &
        -k(74)*n(idx_H2O)

    !d[OH_dot]/d[H+]
    pd(10,55) =  &
        +k(440)*n(idx_NO2)  &
        -k(84)*n(idx_OH)

    !d[O2_dot]/d[H+]
    pd(11,55) =  &
        -k(82)*n(idx_O2)

    !d[CH2_dot]/d[H+]
    pd(12,55) =  &
        -k(428)*n(idx_CH2)  &
        -k(69)*n(idx_CH2)

    !d[H2CO_dot]/d[H+]
    pd(13,55) =  &
        -k(434)*n(idx_H2CO)  &
        -k(73)*n(idx_H2CO)  &
        -k(435)*n(idx_H2CO)

    !d[HCO_dot]/d[H+]
    pd(14,55) =  &
        -k(437)*n(idx_HCO)  &
        -k(436)*n(idx_HCO)  &
        -k(76)*n(idx_HCO)

    !d[MG_dot]/d[H+]
    pd(15,55) =  &
        -k(77)*n(idx_MG)

    !d[NH3_dot]/d[H+]
    pd(16,55) =  &
        -k(79)*n(idx_NH3)

    !d[NO_dot]/d[H+]
    pd(17,55) =  &
        -k(81)*n(idx_NO)

    !d[CO_dot]/d[H+]
    pd(19,55) =  &
        +k(438)*n(idx_HNCO)  &
        +k(437)*n(idx_HCO)

    !d[NH2_dot]/d[H+]
    pd(21,55) =  &
        -k(78)*n(idx_NH2)

    !d[CH3_dot]/d[H+]
    pd(22,55) =  &
        -k(70)*n(idx_CH3)

    !d[CH4_dot]/d[H+]
    pd(23,55) =  &
        -k(432)*n(idx_CH4)  &
        -k(71)*n(idx_CH4)

    !d[NH_dot]/d[H+]
    pd(25,55) =  &
        -k(80)*n(idx_NH)

    !d[HE_dot]/d[H+]
    pd(26,55) =  &
        -k(1046)*n(idx_HE)

    !d[HNO_dot]/d[H+]
    pd(27,55) =  &
        -k(439)*n(idx_HNO)

    !d[CH3OH_dot]/d[H+]
    pd(28,55) =  &
        -k(430)*n(idx_CH3OH)  &
        -k(431)*n(idx_CH3OH)  &
        -k(429)*n(idx_CH3OH)

    !d[CO2_dot]/d[H+]
    pd(29,55) =  &
        -k(433)*n(idx_CO2)

    !d[HNCO_dot]/d[H+]
    pd(31,55) =  &
        -k(438)*n(idx_HNCO)

    !d[NO2_dot]/d[H+]
    pd(32,55) =  &
        -k(440)*n(idx_NO2)

    !d[HCO+_dot]/d[H+]
    pd(54,55) =  &
        +k(433)*n(idx_CO2)  &
        +k(76)*n(idx_HCO)  &
        +k(435)*n(idx_H2CO)  &
        +k(431)*n(idx_CH3OH)

    !d[H+_dot]/d[H+]
    pd(55,55) =  &
        -k(75)*n(idx_HCN)  &
        -k(434)*n(idx_H2CO)  &
        -k(71)*n(idx_CH4)  &
        -k(69)*n(idx_CH2)  &
        -k(433)*n(idx_CO2)  &
        -k(428)*n(idx_CH2)  &
        -k(440)*n(idx_NO2)  &
        -k(430)*n(idx_CH3OH)  &
        -k(74)*n(idx_H2O)  &
        -k(439)*n(idx_HNO)  &
        -k(2)*n(idx_HNC)  &
        -k(80)*n(idx_NH)  &
        -k(83)*n(idx_O)  &
        -k(1058)*n(idx_E)  &
        -k(70)*n(idx_CH3)  &
        -k(84)*n(idx_OH)  &
        -k(73)*n(idx_H2CO)  &
        -k(437)*n(idx_HCO)  &
        -k(431)*n(idx_CH3OH)  &
        -k(429)*n(idx_CH3OH)  &
        -k(81)*n(idx_NO)  &
        -k(435)*n(idx_H2CO)  &
        -k(72)*n(idx_CH)  &
        -k(78)*n(idx_NH2)  &
        -k(1046)*n(idx_HE)  &
        -k(432)*n(idx_CH4)  &
        -k(1045)*n(idx_H)  &
        -k(77)*n(idx_MG)  &
        -k(76)*n(idx_HCO)  &
        -k(79)*n(idx_NH3)  &
        +k(2)*n(idx_HNC)  &
        -k(436)*n(idx_HCO)  &
        -k(82)*n(idx_O2)  &
        -k(438)*n(idx_HNCO)

    !d[CH2+_dot]/d[H+]
    pd(58,55) =  &
        +k(69)*n(idx_CH2)

    !d[CH+_dot]/d[H+]
    pd(59,55) =  &
        +k(72)*n(idx_CH)  &
        +k(428)*n(idx_CH2)

    !d[H2CO+_dot]/d[H+]
    pd(60,55) =  &
        +k(73)*n(idx_H2CO)

    !d[MG+_dot]/d[H+]
    pd(61,55) =  &
        +k(77)*n(idx_MG)

    !d[NH3+_dot]/d[H+]
    pd(62,55) =  &
        +k(79)*n(idx_NH3)

    !d[NO+_dot]/d[H+]
    pd(63,55) =  &
        +k(439)*n(idx_HNO)  &
        +k(81)*n(idx_NO)  &
        +k(440)*n(idx_NO2)

    !d[CO+_dot]/d[H+]
    pd(65,55) =  &
        +k(436)*n(idx_HCO)  &
        +k(434)*n(idx_H2CO)

    !d[O2+_dot]/d[H+]
    pd(67,55) =  &
        +k(82)*n(idx_O2)

    !d[H2O+_dot]/d[H+]
    pd(68,55) =  &
        +k(74)*n(idx_H2O)

    !d[NH2+_dot]/d[H+]
    pd(69,55) =  &
        +k(78)*n(idx_NH2)  &
        +k(438)*n(idx_HNCO)

    !d[O+_dot]/d[H+]
    pd(70,55) =  &
        +k(83)*n(idx_O)

    !d[OH+_dot]/d[H+]
    pd(71,55) =  &
        +k(84)*n(idx_OH)

    !d[CH3+_dot]/d[H+]
    pd(72,55) =  &
        +k(429)*n(idx_CH3OH)  &
        +k(432)*n(idx_CH4)  &
        +k(70)*n(idx_CH3)

    !d[CH4+_dot]/d[H+]
    pd(73,55) =  &
        +k(71)*n(idx_CH4)

    !d[HCN+_dot]/d[H+]
    pd(75,55) =  &
        +k(75)*n(idx_HCN)

    !d[NH+_dot]/d[H+]
    pd(76,55) =  &
        +k(80)*n(idx_NH)

    !d[H2+_dot]/d[H+]
    pd(77,55) =  &
        +k(437)*n(idx_HCO)  &
        +k(1045)*n(idx_H)

    !d[H3CO+_dot]/d[H+]
    pd(82,55) =  &
        +k(430)*n(idx_CH3OH)

    !d[HEH+_dot]/d[H+]
    pd(86,55) =  &
        +k(1046)*n(idx_HE)

    !d[E_dot]/d[HOC+]
    pd(1,56) =  &
        -k(300)*n(idx_E)

    !d[H2_dot]/d[HOC+]
    pd(6,56) =  &
        -k(6)*n(idx_H2)  &
        +k(6)*n(idx_H2)

    !d[H_dot]/d[HOC+]
    pd(8,56) =  &
        +k(300)*n(idx_E)

    !d[CO_dot]/d[HOC+]
    pd(19,56) =  &
        +k(300)*n(idx_E)

    !d[H2CO_DUST_dot]/d[HOC+]
    pd(37,56) =  &
        +k(1067)

    !d[HCO+_dot]/d[HOC+]
    pd(54,56) =  &
        +k(6)*n(idx_H2)

    !d[HOC+_dot]/d[HOC+]
    pd(56,56) =  &
        -k(300)*n(idx_E)  &
        -k(6)*n(idx_H2)  &
        -k(1067)

    !d[E_dot]/d[C+]
    pd(1,57) =  &
        -k(1056)*n(idx_E)

    !d[CH_dot]/d[C+]
    pd(2,57) =  &
        +k(318)*n(idx_H2CO)  &
        -k(16)*n(idx_CH)  &
        +k(314)*n(idx_CH3OH)

    !d[O_dot]/d[C+]
    pd(3,57) =  &
        -k(1041)*n(idx_O)  &
        +k(325)*n(idx_O2)

    !d[H2_dot]/d[C+]
    pd(6,57) =  &
        -k(1047)*n(idx_H2)  &
        +k(323)*n(idx_NH3)  &
        -k(460)*n(idx_H2)

    !d[C_dot]/d[C+]
    pd(7,57) =  &
        +k(20)*n(idx_NH3)  &
        +k(21)*n(idx_NO)  &
        +k(19)*n(idx_MG)  &
        +k(17)*n(idx_H2CO)  &
        +k(18)*n(idx_HCO)  &
        +k(15)*n(idx_CH2)  &
        +k(1056)*n(idx_E)  &
        +k(16)*n(idx_CH)

    !d[H_dot]/d[C+]
    pd(8,57) =  &
        +k(320)*n(idx_H2O)  &
        +k(324)*n(idx_NH)  &
        +k(328)*n(idx_OH)  &
        +k(460)*n(idx_H2)  &
        +k(319)*n(idx_H2O)  &
        +k(322)*n(idx_NH2)  &
        -k(1050)*n(idx_H)

    !d[H2O_dot]/d[C+]
    pd(9,57) =  &
        -k(319)*n(idx_H2O)  &
        -k(320)*n(idx_H2O)

    !d[OH_dot]/d[C+]
    pd(10,57) =  &
        -k(328)*n(idx_OH)

    !d[O2_dot]/d[C+]
    pd(11,57) =  &
        -k(325)*n(idx_O2)  &
        -k(326)*n(idx_O2)

    !d[CH2_dot]/d[C+]
    pd(12,57) =  &
        -k(15)*n(idx_CH2)

    !d[H2CO_dot]/d[C+]
    pd(13,57) =  &
        -k(17)*n(idx_H2CO)  &
        -k(317)*n(idx_H2CO)  &
        -k(318)*n(idx_H2CO)

    !d[HCO_dot]/d[C+]
    pd(14,57) =  &
        +k(315)*n(idx_CH3OH)  &
        -k(321)*n(idx_HCO)  &
        -k(18)*n(idx_HCO)

    !d[MG_dot]/d[C+]
    pd(15,57) =  &
        -k(19)*n(idx_MG)

    !d[NH3_dot]/d[C+]
    pd(16,57) =  &
        -k(323)*n(idx_NH3)  &
        -k(20)*n(idx_NH3)

    !d[NO_dot]/d[C+]
    pd(17,57) =  &
        -k(21)*n(idx_NO)

    !d[CN_dot]/d[C+]
    pd(18,57) =  &
        +k(327)*n(idx_OCN)

    !d[CO_dot]/d[C+]
    pd(19,57) =  &
        +k(317)*n(idx_H2CO)  &
        +k(321)*n(idx_HCO)  &
        +k(316)*n(idx_CO2)  &
        +k(326)*n(idx_O2)

    !d[NH2_dot]/d[C+]
    pd(21,57) =  &
        -k(322)*n(idx_NH2)

    !d[N_dot]/d[C+]
    pd(24,57) =  &
        -k(1040)*n(idx_N)

    !d[NH_dot]/d[C+]
    pd(25,57) =  &
        -k(324)*n(idx_NH)

    !d[CH3OH_dot]/d[C+]
    pd(28,57) =  &
        -k(315)*n(idx_CH3OH)  &
        -k(314)*n(idx_CH3OH)

    !d[CO2_dot]/d[C+]
    pd(29,57) =  &
        -k(316)*n(idx_CO2)

    !d[OCN_dot]/d[C+]
    pd(34,57) =  &
        -k(327)*n(idx_OCN)

    !d[CH4_DUST_dot]/d[C+]
    pd(38,57) =  &
        +k(1084)

    !d[HCO+_dot]/d[C+]
    pd(54,57) =  &
        +k(318)*n(idx_H2CO)  &
        +k(319)*n(idx_H2O)  &
        +k(18)*n(idx_HCO)

    !d[HOC+_dot]/d[C+]
    pd(56,57) =  &
        +k(320)*n(idx_H2O)

    !d[C+_dot]/d[C+]
    pd(57,57) =  &
        -k(1040)*n(idx_N)  &
        -k(1056)*n(idx_E)  &
        -k(15)*n(idx_CH2)  &
        -k(17)*n(idx_H2CO)  &
        -k(21)*n(idx_NO)  &
        -k(1041)*n(idx_O)  &
        -k(317)*n(idx_H2CO)  &
        -k(16)*n(idx_CH)  &
        -k(314)*n(idx_CH3OH)  &
        -k(321)*n(idx_HCO)  &
        -k(19)*n(idx_MG)  &
        -k(316)*n(idx_CO2)  &
        -k(1047)*n(idx_H2)  &
        -k(18)*n(idx_HCO)  &
        -k(1050)*n(idx_H)  &
        -k(1084)  &
        -k(320)*n(idx_H2O)  &
        -k(328)*n(idx_OH)  &
        -k(315)*n(idx_CH3OH)  &
        -k(460)*n(idx_H2)  &
        -k(325)*n(idx_O2)  &
        -k(20)*n(idx_NH3)  &
        -k(322)*n(idx_NH2)  &
        -k(324)*n(idx_NH)  &
        -k(319)*n(idx_H2O)  &
        -k(323)*n(idx_NH3)  &
        -k(326)*n(idx_O2)  &
        -k(327)*n(idx_OCN)  &
        -k(318)*n(idx_H2CO)

    !d[CH2+_dot]/d[C+]
    pd(58,57) =  &
        +k(15)*n(idx_CH2)  &
        +k(317)*n(idx_H2CO)  &
        +k(1047)*n(idx_H2)

    !d[CH+_dot]/d[C+]
    pd(59,57) =  &
        +k(460)*n(idx_H2)  &
        +k(321)*n(idx_HCO)  &
        +k(1050)*n(idx_H)  &
        +k(16)*n(idx_CH)

    !d[H2CO+_dot]/d[C+]
    pd(60,57) =  &
        +k(17)*n(idx_H2CO)

    !d[MG+_dot]/d[C+]
    pd(61,57) =  &
        +k(19)*n(idx_MG)

    !d[NH3+_dot]/d[C+]
    pd(62,57) =  &
        +k(20)*n(idx_NH3)

    !d[NO+_dot]/d[C+]
    pd(63,57) =  &
        +k(21)*n(idx_NO)

    !d[CN+_dot]/d[C+]
    pd(64,57) =  &
        +k(324)*n(idx_NH)  &
        +k(1040)*n(idx_N)

    !d[CO+_dot]/d[C+]
    pd(65,57) =  &
        +k(327)*n(idx_OCN)  &
        +k(325)*n(idx_O2)  &
        +k(316)*n(idx_CO2)  &
        +k(1041)*n(idx_O)  &
        +k(328)*n(idx_OH)

    !d[O+_dot]/d[C+]
    pd(70,57) =  &
        +k(326)*n(idx_O2)

    !d[CH3+_dot]/d[C+]
    pd(72,57) =  &
        +k(315)*n(idx_CH3OH)

    !d[HCN+_dot]/d[C+]
    pd(75,57) =  &
        +k(323)*n(idx_NH3)  &
        +k(322)*n(idx_NH2)

    !d[H3CO+_dot]/d[C+]
    pd(82,57) =  &
        +k(314)*n(idx_CH3OH)

    !d[E_dot]/d[CH2+]
    pd(1,58) =  &
        -k(262)*n(idx_E)  &
        -k(261)*n(idx_E)  &
        -k(260)*n(idx_E)

    !d[CH_dot]/d[CH2+]
    pd(2,58) =  &
        +k(262)*n(idx_E)  &
        +k(977)

    !d[O_dot]/d[CH2+]
    pd(3,58) =  &
        -k(365)*n(idx_O)

    !d[H2_dot]/d[CH2+]
    pd(6,58) =  &
        +k(531)*n(idx_H)  &
        +k(975)  &
        -k(462)*n(idx_H2)  &
        +k(260)*n(idx_E)

    !d[C_dot]/d[CH2+]
    pd(7,58) =  &
        +k(261)*n(idx_E)  &
        +k(260)*n(idx_E)

    !d[H_dot]/d[CH2+]
    pd(8,58) =  &
        +k(262)*n(idx_E)  &
        +k(976)  &
        +k(365)*n(idx_O)  &
        +k(462)*n(idx_H2)  &
        +k(362)*n(idx_H2O)  &
        -k(531)*n(idx_H)  &
        +k(634)*n(idx_N)  &
        +2.d0*k(261)*n(idx_E)

    !d[H2O_dot]/d[CH2+]
    pd(9,58) =  &
        -k(362)*n(idx_H2O)

    !d[OH_dot]/d[CH2+]
    pd(10,58) =  &
        +k(364)*n(idx_O2)

    !d[O2_dot]/d[CH2+]
    pd(11,58) =  &
        -k(364)*n(idx_O2)

    !d[CH2_dot]/d[CH2+]
    pd(12,58) =  &
        +k(30)*n(idx_NO)

    !d[H2CO_dot]/d[CH2+]
    pd(13,58) =  &
        -k(361)*n(idx_H2CO)

    !d[HCO_dot]/d[CH2+]
    pd(14,58) =  &
        -k(363)*n(idx_HCO)

    !d[NO_dot]/d[CH2+]
    pd(17,58) =  &
        -k(30)*n(idx_NO)

    !d[CO_dot]/d[CH2+]
    pd(19,58) =  &
        +k(363)*n(idx_HCO)  &
        +k(360)*n(idx_CO2)

    !d[CH3_dot]/d[CH2+]
    pd(22,58) =  &
        +k(361)*n(idx_H2CO)

    !d[N_dot]/d[CH2+]
    pd(24,58) =  &
        -k(634)*n(idx_N)

    !d[CO2_dot]/d[CH2+]
    pd(29,58) =  &
        -k(360)*n(idx_CO2)

    !d[CH4_DUST_dot]/d[CH2+]
    pd(38,58) =  &
        +k(1098)

    !d[HCO+_dot]/d[CH2+]
    pd(54,58) =  &
        +k(361)*n(idx_H2CO)  &
        +k(364)*n(idx_O2)  &
        +k(365)*n(idx_O)

    !d[H+_dot]/d[CH2+]
    pd(55,58) =  &
        +k(977)

    !d[C+_dot]/d[CH2+]
    pd(57,58) =  &
        +k(975)

    !d[CH2+_dot]/d[CH2+]
    pd(58,58) =  &
        -k(365)*n(idx_O)  &
        -k(363)*n(idx_HCO)  &
        -k(364)*n(idx_O2)  &
        -k(976)  &
        -k(634)*n(idx_N)  &
        -k(977)  &
        -k(261)*n(idx_E)  &
        -k(260)*n(idx_E)  &
        -k(462)*n(idx_H2)  &
        -k(30)*n(idx_NO)  &
        -k(975)  &
        -k(262)*n(idx_E)  &
        -k(531)*n(idx_H)  &
        -k(360)*n(idx_CO2)  &
        -k(1098)  &
        -k(362)*n(idx_H2O)  &
        -k(361)*n(idx_H2CO)

    !d[CH+_dot]/d[CH2+]
    pd(59,58) =  &
        +k(531)*n(idx_H)  &
        +k(976)

    !d[H2CO+_dot]/d[CH2+]
    pd(60,58) =  &
        +k(360)*n(idx_CO2)

    !d[NO+_dot]/d[CH2+]
    pd(63,58) =  &
        +k(30)*n(idx_NO)

    !d[CH3+_dot]/d[CH2+]
    pd(72,58) =  &
        +k(462)*n(idx_H2)  &
        +k(363)*n(idx_HCO)

    !d[HCN+_dot]/d[CH2+]
    pd(75,58) =  &
        +k(634)*n(idx_N)

    !d[H3CO+_dot]/d[CH2+]
    pd(82,58) =  &
        +k(362)*n(idx_H2O)

    !d[E_dot]/d[CH+]
    pd(1,59) =  &
        -k(259)*n(idx_E)

    !d[CH_dot]/d[CH+]
    pd(2,59) =  &
        +k(28)*n(idx_NH3)  &
        +k(26)*n(idx_HCO)  &
        +k(29)*n(idx_NO)  &
        +k(27)*n(idx_MG)

    !d[O_dot]/d[CH+]
    pd(3,59) =  &
        +k(356)*n(idx_O2)  &
        -k(358)*n(idx_O)

    !d[HNC_dot]/d[CH+]
    pd(4,59) =  &
        -k(351)*n(idx_HNC)

    !d[HCN_dot]/d[CH+]
    pd(5,59) =  &
        -k(349)*n(idx_HCN)

    !d[H2_dot]/d[CH+]
    pd(6,59) =  &
        +k(359)*n(idx_OH)  &
        +k(353)*n(idx_NH2)  &
        -k(461)*n(idx_H2)  &
        +k(354)*n(idx_NH)  &
        +k(348)*n(idx_H2O)  &
        +k(530)*n(idx_H)

    !d[C_dot]/d[CH+]
    pd(7,59) =  &
        +k(344)*n(idx_H2CO)  &
        +k(349)*n(idx_HCN)  &
        +k(351)*n(idx_HNC)  &
        +k(347)*n(idx_H2O)  &
        +k(259)*n(idx_E)  &
        +k(974)

    !d[H_dot]/d[CH+]
    pd(8,59) =  &
        +k(352)*n(idx_N)  &
        -k(530)*n(idx_H)  &
        +k(358)*n(idx_O)  &
        +k(346)*n(idx_H2O)  &
        +k(216)  &
        +k(259)*n(idx_E)  &
        +k(461)*n(idx_H2)

    !d[H2O_dot]/d[CH+]
    pd(9,59) =  &
        -k(348)*n(idx_H2O)  &
        -k(346)*n(idx_H2O)  &
        -k(347)*n(idx_H2O)

    !d[OH_dot]/d[CH+]
    pd(10,59) =  &
        +k(355)*n(idx_O2)  &
        -k(359)*n(idx_OH)

    !d[O2_dot]/d[CH+]
    pd(11,59) =  &
        -k(357)*n(idx_O2)  &
        -k(356)*n(idx_O2)  &
        -k(355)*n(idx_O2)

    !d[CH2_dot]/d[CH+]
    pd(12,59) =  &
        +k(341)*n(idx_CH3OH)  &
        +k(345)*n(idx_H2CO)

    !d[H2CO_dot]/d[CH+]
    pd(13,59) =  &
        -k(345)*n(idx_H2CO)  &
        -k(343)*n(idx_H2CO)  &
        +k(340)*n(idx_CH3OH)  &
        -k(344)*n(idx_H2CO)

    !d[HCO_dot]/d[CH+]
    pd(14,59) =  &
        -k(350)*n(idx_HCO)  &
        -k(26)*n(idx_HCO)  &
        +k(357)*n(idx_O2)

    !d[MG_dot]/d[CH+]
    pd(15,59) =  &
        -k(27)*n(idx_MG)

    !d[NH3_dot]/d[CH+]
    pd(16,59) =  &
        -k(28)*n(idx_NH3)

    !d[NO_dot]/d[CH+]
    pd(17,59) =  &
        -k(29)*n(idx_NO)

    !d[CO_dot]/d[CH+]
    pd(19,59) =  &
        +k(343)*n(idx_H2CO)  &
        +k(350)*n(idx_HCO)  &
        +k(342)*n(idx_CO2)

    !d[NH2_dot]/d[CH+]
    pd(21,59) =  &
        -k(353)*n(idx_NH2)

    !d[N_dot]/d[CH+]
    pd(24,59) =  &
        -k(352)*n(idx_N)

    !d[NH_dot]/d[CH+]
    pd(25,59) =  &
        -k(354)*n(idx_NH)

    !d[CH3OH_dot]/d[CH+]
    pd(28,59) =  &
        -k(340)*n(idx_CH3OH)  &
        -k(341)*n(idx_CH3OH)

    !d[CO2_dot]/d[CH+]
    pd(29,59) =  &
        -k(342)*n(idx_CO2)

    !d[CH4_DUST_dot]/d[CH+]
    pd(38,59) =  &
        +k(1092)

    !d[HCO+_dot]/d[CH+]
    pd(54,59) =  &
        +k(356)*n(idx_O2)  &
        +k(348)*n(idx_H2O)  &
        +k(26)*n(idx_HCO)  &
        +k(345)*n(idx_H2CO)  &
        +k(342)*n(idx_CO2)

    !d[H+_dot]/d[CH+]
    pd(55,59) =  &
        +k(974)

    !d[C+_dot]/d[CH+]
    pd(57,59) =  &
        +k(216)  &
        +k(530)*n(idx_H)

    !d[CH2+_dot]/d[CH+]
    pd(58,59) =  &
        +k(461)*n(idx_H2)  &
        +k(350)*n(idx_HCO)

    !d[CH+_dot]/d[CH+]
    pd(59,59) =  &
        -k(216)  &
        -k(530)*n(idx_H)  &
        -k(357)*n(idx_O2)  &
        -k(461)*n(idx_H2)  &
        -k(29)*n(idx_NO)  &
        -k(358)*n(idx_O)  &
        -k(345)*n(idx_H2CO)  &
        -k(350)*n(idx_HCO)  &
        -k(353)*n(idx_NH2)  &
        -k(341)*n(idx_CH3OH)  &
        -k(28)*n(idx_NH3)  &
        -k(359)*n(idx_OH)  &
        -k(1092)  &
        -k(27)*n(idx_MG)  &
        -k(340)*n(idx_CH3OH)  &
        -k(259)*n(idx_E)  &
        -k(355)*n(idx_O2)  &
        -k(351)*n(idx_HNC)  &
        -k(349)*n(idx_HCN)  &
        -k(26)*n(idx_HCO)  &
        -k(346)*n(idx_H2O)  &
        -k(344)*n(idx_H2CO)  &
        -k(348)*n(idx_H2O)  &
        -k(342)*n(idx_CO2)  &
        -k(343)*n(idx_H2CO)  &
        -k(354)*n(idx_NH)  &
        -k(356)*n(idx_O2)  &
        -k(974)  &
        -k(347)*n(idx_H2O)  &
        -k(352)*n(idx_N)

    !d[H2CO+_dot]/d[CH+]
    pd(60,59) =  &
        +k(346)*n(idx_H2O)

    !d[MG+_dot]/d[CH+]
    pd(61,59) =  &
        +k(27)*n(idx_MG)

    !d[NH3+_dot]/d[CH+]
    pd(62,59) =  &
        +k(28)*n(idx_NH3)

    !d[NO+_dot]/d[CH+]
    pd(63,59) =  &
        +k(29)*n(idx_NO)

    !d[CN+_dot]/d[CH+]
    pd(64,59) =  &
        +k(352)*n(idx_N)  &
        +k(354)*n(idx_NH)

    !d[CO+_dot]/d[CH+]
    pd(65,59) =  &
        +k(355)*n(idx_O2)  &
        +k(359)*n(idx_OH)  &
        +k(358)*n(idx_O)

    !d[O+_dot]/d[CH+]
    pd(70,59) =  &
        +k(357)*n(idx_O2)

    !d[CH3+_dot]/d[CH+]
    pd(72,59) =  &
        +k(343)*n(idx_H2CO)  &
        +k(340)*n(idx_CH3OH)

    !d[HCN+_dot]/d[CH+]
    pd(75,59) =  &
        +k(353)*n(idx_NH2)

    !d[H3CO+_dot]/d[CH+]
    pd(82,59) =  &
        +k(344)*n(idx_H2CO)  &
        +k(341)*n(idx_CH3OH)

    !d[H3O+_dot]/d[CH+]
    pd(83,59) =  &
        +k(347)*n(idx_H2O)

    !d[HCNH+_dot]/d[CH+]
    pd(84,59) =  &
        +k(351)*n(idx_HNC)  &
        +k(349)*n(idx_HCN)

    !d[E_dot]/d[H2CO+]
    pd(1,60) =  &
        -k(1059)*n(idx_E)  &
        -k(271)*n(idx_E)  &
        -k(272)*n(idx_E)  &
        -k(273)*n(idx_E)  &
        -k(274)*n(idx_E)

    !d[CH_dot]/d[H2CO+]
    pd(2,60) =  &
        -k(49)*n(idx_CH)  &
        -k(401)*n(idx_CH)

    !d[O_dot]/d[H2CO+]
    pd(3,60) =  &
        +k(271)*n(idx_E)

    !d[HNC_dot]/d[H2CO+]
    pd(4,60) =  &
        -k(557)*n(idx_HNC)

    !d[HCN_dot]/d[H2CO+]
    pd(5,60) =  &
        -k(542)*n(idx_HCN)

    !d[H2_dot]/d[H2CO+]
    pd(6,60) =  &
        +k(272)*n(idx_E)

    !d[H_dot]/d[H2CO+]
    pd(8,60) =  &
        +k(274)*n(idx_E)  &
        +2.d0*k(273)*n(idx_E)

    !d[H2O_dot]/d[H2CO+]
    pd(9,60) =  &
        -k(493)*n(idx_H2O)

    !d[O2_dot]/d[H2CO+]
    pd(11,60) =  &
        -k(479)*n(idx_O2)

    !d[CH2_dot]/d[H2CO+]
    pd(12,60) =  &
        +k(271)*n(idx_E)  &
        -k(33)*n(idx_CH2)  &
        -k(367)*n(idx_CH2)

    !d[H2CO_dot]/d[H2CO+]
    pd(13,60) =  &
        +k(33)*n(idx_CH2)  &
        +k(130)*n(idx_MG)  &
        +k(120)*n(idx_HCO)  &
        +k(182)*n(idx_NO)  &
        -k(478)*n(idx_H2CO)  &
        +k(49)*n(idx_CH)  &
        +k(1059)*n(idx_E)  &
        +k(173)*n(idx_NH3)

    !d[HCO_dot]/d[H2CO+]
    pd(14,60) =  &
        +k(274)*n(idx_E)  &
        +k(675)*n(idx_NH2)  &
        -k(120)*n(idx_HCO)  &
        +k(367)*n(idx_CH2)  &
        +k(542)*n(idx_HCN)  &
        +k(493)*n(idx_H2O)  &
        +k(478)*n(idx_H2CO)  &
        +k(401)*n(idx_CH)  &
        +k(557)*n(idx_HNC)  &
        -k(552)*n(idx_HCO)

    !d[MG_dot]/d[H2CO+]
    pd(15,60) =  &
        -k(130)*n(idx_MG)

    !d[NH3_dot]/d[H2CO+]
    pd(16,60) =  &
        -k(173)*n(idx_NH3)

    !d[NO_dot]/d[H2CO+]
    pd(17,60) =  &
        -k(182)*n(idx_NO)

    !d[CO_dot]/d[H2CO+]
    pd(19,60) =  &
        +k(272)*n(idx_E)  &
        +k(273)*n(idx_E)  &
        +k(552)*n(idx_HCO)

    !d[NH2_dot]/d[H2CO+]
    pd(21,60) =  &
        -k(675)*n(idx_NH2)

    !d[CH3_dot]/d[H2CO+]
    pd(22,60) =  &
        +k(394)*n(idx_CH4)

    !d[CH4_dot]/d[H2CO+]
    pd(23,60) =  &
        -k(394)*n(idx_CH4)

    !d[N_dot]/d[H2CO+]
    pd(24,60) =  &
        +k(691)*n(idx_NH)

    !d[NH_dot]/d[H2CO+]
    pd(25,60) =  &
        -k(691)*n(idx_NH)

    !d[O2H_dot]/d[H2CO+]
    pd(33,60) =  &
        +k(479)*n(idx_O2)

    !d[H2CO_DUST_dot]/d[H2CO+]
    pd(37,60) =  &
        +k(1104)

    !d[HCO+_dot]/d[H2CO+]
    pd(54,60) =  &
        +k(479)*n(idx_O2)  &
        +k(120)*n(idx_HCO)

    !d[CH2+_dot]/d[H2CO+]
    pd(58,60) =  &
        +k(33)*n(idx_CH2)  &
        +k(401)*n(idx_CH)

    !d[CH+_dot]/d[H2CO+]
    pd(59,60) =  &
        +k(49)*n(idx_CH)

    !d[H2CO+_dot]/d[H2CO+]
    pd(60,60) =  &
        -k(273)*n(idx_E)  &
        -k(130)*n(idx_MG)  &
        -k(675)*n(idx_NH2)  &
        -k(401)*n(idx_CH)  &
        -k(182)*n(idx_NO)  &
        -k(120)*n(idx_HCO)  &
        -k(1104)  &
        -k(173)*n(idx_NH3)  &
        -k(493)*n(idx_H2O)  &
        -k(274)*n(idx_E)  &
        -k(272)*n(idx_E)  &
        -k(394)*n(idx_CH4)  &
        -k(691)*n(idx_NH)  &
        -k(542)*n(idx_HCN)  &
        -k(271)*n(idx_E)  &
        -k(478)*n(idx_H2CO)  &
        -k(33)*n(idx_CH2)  &
        -k(552)*n(idx_HCO)  &
        -k(49)*n(idx_CH)  &
        -k(479)*n(idx_O2)  &
        -k(1059)*n(idx_E)  &
        -k(557)*n(idx_HNC)  &
        -k(367)*n(idx_CH2)

    !d[MG+_dot]/d[H2CO+]
    pd(61,60) =  &
        +k(130)*n(idx_MG)

    !d[NH3+_dot]/d[H2CO+]
    pd(62,60) =  &
        +k(173)*n(idx_NH3)  &
        +k(675)*n(idx_NH2)

    !d[NO+_dot]/d[H2CO+]
    pd(63,60) =  &
        +k(182)*n(idx_NO)

    !d[CH3+_dot]/d[H2CO+]
    pd(72,60) =  &
        +k(367)*n(idx_CH2)

    !d[H3CO+_dot]/d[H2CO+]
    pd(82,60) =  &
        +k(394)*n(idx_CH4)  &
        +k(552)*n(idx_HCO)  &
        +k(691)*n(idx_NH)  &
        +k(478)*n(idx_H2CO)

    !d[H3O+_dot]/d[H2CO+]
    pd(83,60) =  &
        +k(493)*n(idx_H2O)

    !d[HCNH+_dot]/d[H2CO+]
    pd(84,60) =  &
        +k(557)*n(idx_HNC)  &
        +k(542)*n(idx_HCN)

    !d[E_dot]/d[MG+]
    pd(1,61) =  &
        -k(1061)*n(idx_E)

    !d[MG_dot]/d[MG+]
    pd(15,61) =  &
        +k(1061)*n(idx_E)

    !d[MG_DUST_dot]/d[MG+]
    pd(51,61) =  &
        +k(1119)

    !d[MG+_dot]/d[MG+]
    pd(61,61) =  &
        -k(1061)*n(idx_E)  &
        -k(1119)

    !d[E_dot]/d[NH3+]
    pd(1,62) =  &
        -k(308)*n(idx_E)  &
        -k(309)*n(idx_E)

    !d[O_dot]/d[NH3+]
    pd(3,62) =  &
        -k(723)*n(idx_O)

    !d[H2_dot]/d[NH3+]
    pd(6,62) =  &
        +k(723)*n(idx_O)

    !d[H_dot]/d[NH3+]
    pd(8,62) =  &
        +k(308)*n(idx_E)  &
        +2.d0*k(309)*n(idx_E)

    !d[CH2_dot]/d[NH3+]
    pd(12,62) =  &
        -k(378)*n(idx_CH2)

    !d[HCO_dot]/d[NH3+]
    pd(14,62) =  &
        -k(169)*n(idx_HCO)

    !d[MG_dot]/d[NH3+]
    pd(15,62) =  &
        -k(170)*n(idx_MG)

    !d[NH3_dot]/d[NH3+]
    pd(16,62) =  &
        +k(170)*n(idx_MG)  &
        +k(169)*n(idx_HCO)  &
        +k(171)*n(idx_NO)

    !d[NO_dot]/d[NH3+]
    pd(17,62) =  &
        -k(171)*n(idx_NO)

    !d[NH2_dot]/d[NH3+]
    pd(21,62) =  &
        +k(308)*n(idx_E)  &
        +k(378)*n(idx_CH2)

    !d[NH_dot]/d[NH3+]
    pd(25,62) =  &
        +k(309)*n(idx_E)

    !d[NH3_DUST_dot]/d[NH3+]
    pd(45,62) =  &
        +k(1103)

    !d[HCO+_dot]/d[NH3+]
    pd(54,62) =  &
        +k(169)*n(idx_HCO)

    !d[MG+_dot]/d[NH3+]
    pd(61,62) =  &
        +k(170)*n(idx_MG)

    !d[NH3+_dot]/d[NH3+]
    pd(62,62) =  &
        -k(308)*n(idx_E)  &
        -k(171)*n(idx_NO)  &
        -k(723)*n(idx_O)  &
        -k(309)*n(idx_E)  &
        -k(170)*n(idx_MG)  &
        -k(1103)  &
        -k(169)*n(idx_HCO)  &
        -k(378)*n(idx_CH2)

    !d[NO+_dot]/d[NH3+]
    pd(63,62) =  &
        +k(171)*n(idx_NO)

    !d[CH3+_dot]/d[NH3+]
    pd(72,62) =  &
        +k(378)*n(idx_CH2)

    !d[HNO+_dot]/d[NH3+]
    pd(79,62) =  &
        +k(723)*n(idx_O)

    !d[E_dot]/d[NO+]
    pd(1,63) =  &
        -k(310)*n(idx_E)

    !d[O_dot]/d[NO+]
    pd(3,63) =  &
        +k(310)*n(idx_E)

    !d[MG_dot]/d[NO+]
    pd(15,63) =  &
        -k(133)*n(idx_MG)

    !d[NO_dot]/d[NO+]
    pd(17,63) =  &
        +k(133)*n(idx_MG)

    !d[N_dot]/d[NO+]
    pd(24,63) =  &
        +k(310)*n(idx_E)

    !d[NO_DUST_dot]/d[NO+]
    pd(41,63) =  &
        +k(1097)

    !d[MG+_dot]/d[NO+]
    pd(61,63) =  &
        +k(133)*n(idx_MG)

    !d[NO+_dot]/d[NO+]
    pd(63,63) =  &
        -k(1097)  &
        -k(133)*n(idx_MG)  &
        -k(310)*n(idx_E)

    !d[E_dot]/d[CN+]
    pd(1,64) =  &
        -k(268)*n(idx_E)

    !d[CH_dot]/d[CN+]
    pd(2,64) =  &
        -k(47)*n(idx_CH)

    !d[O_dot]/d[CN+]
    pd(3,64) =  &
        -k(194)*n(idx_O)

    !d[HCN_dot]/d[CN+]
    pd(5,64) =  &
        +k(418)*n(idx_H2CO)  &
        -k(59)*n(idx_HCN)

    !d[H2_dot]/d[CN+]
    pd(6,64) =  &
        -k(463)*n(idx_H2)

    !d[C_dot]/d[CN+]
    pd(7,64) =  &
        +k(635)*n(idx_N)  &
        -k(22)*n(idx_C)  &
        +k(268)*n(idx_E)

    !d[H_dot]/d[CN+]
    pd(8,64) =  &
        +k(463)*n(idx_H2)  &
        -k(110)*n(idx_H)

    !d[H2O_dot]/d[CN+]
    pd(9,64) =  &
        -k(490)*n(idx_H2O)  &
        -k(491)*n(idx_H2O)

    !d[OH_dot]/d[CN+]
    pd(10,64) =  &
        -k(203)*n(idx_OH)  &
        +k(490)*n(idx_H2O)

    !d[O2_dot]/d[CN+]
    pd(11,64) =  &
        -k(62)*n(idx_O2)  &
        -k(420)*n(idx_O2)

    !d[CH2_dot]/d[CN+]
    pd(12,64) =  &
        -k(31)*n(idx_CH2)

    !d[H2CO_dot]/d[CN+]
    pd(13,64) =  &
        -k(418)*n(idx_H2CO)  &
        -k(58)*n(idx_H2CO)

    !d[HCO_dot]/d[CN+]
    pd(14,64) =  &
        -k(60)*n(idx_HCO)  &
        -k(419)*n(idx_HCO)

    !d[NO_dot]/d[CN+]
    pd(17,64) =  &
        -k(61)*n(idx_NO)

    !d[CN_dot]/d[CN+]
    pd(18,64) =  &
        +k(47)*n(idx_CH)  &
        +k(22)*n(idx_C)  &
        +k(62)*n(idx_O2)  &
        +k(58)*n(idx_H2CO)  &
        +k(57)*n(idx_CO)  &
        +k(178)*n(idx_NH)  &
        +k(110)*n(idx_H)  &
        +k(59)*n(idx_HCN)  &
        +k(194)*n(idx_O)  &
        +k(61)*n(idx_NO)  &
        +k(31)*n(idx_CH2)  &
        +k(203)*n(idx_OH)  &
        +k(163)*n(idx_NH2)  &
        +k(60)*n(idx_HCO)

    !d[CO_dot]/d[CN+]
    pd(19,64) =  &
        -k(57)*n(idx_CO)  &
        +k(420)*n(idx_O2)  &
        +k(419)*n(idx_HCO)

    !d[NH2_dot]/d[CN+]
    pd(21,64) =  &
        -k(163)*n(idx_NH2)

    !d[N_dot]/d[CN+]
    pd(24,64) =  &
        -k(635)*n(idx_N)  &
        +k(268)*n(idx_E)

    !d[NH_dot]/d[CN+]
    pd(25,64) =  &
        -k(178)*n(idx_NH)  &
        +k(491)*n(idx_H2O)

    !d[HCN_DUST_dot]/d[CN+]
    pd(44,64) =  &
        +k(1096)

    !d[HCO+_dot]/d[CN+]
    pd(54,64) =  &
        +k(418)*n(idx_H2CO)  &
        +k(491)*n(idx_H2O)  &
        +k(60)*n(idx_HCO)

    !d[H+_dot]/d[CN+]
    pd(55,64) =  &
        +k(110)*n(idx_H)

    !d[C+_dot]/d[CN+]
    pd(57,64) =  &
        +k(22)*n(idx_C)

    !d[CH2+_dot]/d[CN+]
    pd(58,64) =  &
        +k(31)*n(idx_CH2)

    !d[CH+_dot]/d[CN+]
    pd(59,64) =  &
        +k(47)*n(idx_CH)

    !d[H2CO+_dot]/d[CN+]
    pd(60,64) =  &
        +k(58)*n(idx_H2CO)

    !d[NO+_dot]/d[CN+]
    pd(63,64) =  &
        +k(61)*n(idx_NO)  &
        +k(420)*n(idx_O2)

    !d[CN+_dot]/d[CN+]
    pd(64,64) =  &
        -k(62)*n(idx_O2)  &
        -k(31)*n(idx_CH2)  &
        -k(194)*n(idx_O)  &
        -k(635)*n(idx_N)  &
        -k(60)*n(idx_HCO)  &
        -k(59)*n(idx_HCN)  &
        -k(47)*n(idx_CH)  &
        -k(419)*n(idx_HCO)  &
        -k(22)*n(idx_C)  &
        -k(491)*n(idx_H2O)  &
        -k(463)*n(idx_H2)  &
        -k(61)*n(idx_NO)  &
        -k(58)*n(idx_H2CO)  &
        -k(418)*n(idx_H2CO)  &
        -k(490)*n(idx_H2O)  &
        -k(178)*n(idx_NH)  &
        -k(1096)  &
        -k(268)*n(idx_E)  &
        -k(420)*n(idx_O2)  &
        -k(57)*n(idx_CO)  &
        -k(163)*n(idx_NH2)  &
        -k(110)*n(idx_H)  &
        -k(203)*n(idx_OH)

    !d[CO+_dot]/d[CN+]
    pd(65,64) =  &
        +k(57)*n(idx_CO)

    !d[N2+_dot]/d[CN+]
    pd(66,64) =  &
        +k(635)*n(idx_N)

    !d[O2+_dot]/d[CN+]
    pd(67,64) =  &
        +k(62)*n(idx_O2)

    !d[NH2+_dot]/d[CN+]
    pd(69,64) =  &
        +k(163)*n(idx_NH2)

    !d[O+_dot]/d[CN+]
    pd(70,64) =  &
        +k(194)*n(idx_O)

    !d[OH+_dot]/d[CN+]
    pd(71,64) =  &
        +k(203)*n(idx_OH)

    !d[HCN+_dot]/d[CN+]
    pd(75,64) =  &
        +k(59)*n(idx_HCN)  &
        +k(463)*n(idx_H2)  &
        +k(419)*n(idx_HCO)  &
        +k(490)*n(idx_H2O)

    !d[NH+_dot]/d[CN+]
    pd(76,64) =  &
        +k(178)*n(idx_NH)

    !d[E_dot]/d[CO+]
    pd(1,65) =  &
        -k(269)*n(idx_E)

    !d[CH_dot]/d[CO+]
    pd(2,65) =  &
        -k(400)*n(idx_CH)  &
        +k(366)*n(idx_CH2)  &
        -k(48)*n(idx_CH)

    !d[O_dot]/d[CO+]
    pd(3,65) =  &
        -k(195)*n(idx_O)  &
        +k(997)  &
        +k(738)*n(idx_OH)  &
        +k(269)*n(idx_E)

    !d[HCN_dot]/d[CO+]
    pd(5,65) =  &
        -k(118)*n(idx_HCN)

    !d[H2_dot]/d[CO+]
    pd(6,65) =  &
        -k(464)*n(idx_H2)  &
        -k(465)*n(idx_H2)

    !d[C_dot]/d[CO+]
    pd(7,65) =  &
        +k(400)*n(idx_CH)  &
        -k(23)*n(idx_C)  &
        +k(269)*n(idx_E)

    !d[H_dot]/d[CO+]
    pd(8,65) =  &
        +k(465)*n(idx_H2)  &
        -k(111)*n(idx_H)  &
        +k(464)*n(idx_H2)

    !d[H2O_dot]/d[CO+]
    pd(9,65) =  &
        -k(107)*n(idx_H2O)  &
        -k(492)*n(idx_H2O)

    !d[OH_dot]/d[CO+]
    pd(10,65) =  &
        -k(204)*n(idx_OH)  &
        +k(492)*n(idx_H2O)  &
        -k(738)*n(idx_OH)

    !d[O2_dot]/d[CO+]
    pd(11,65) =  &
        -k(67)*n(idx_O2)

    !d[CH2_dot]/d[CO+]
    pd(12,65) =  &
        -k(366)*n(idx_CH2)  &
        -k(32)*n(idx_CH2)

    !d[H2CO_dot]/d[CO+]
    pd(13,65) =  &
        -k(423)*n(idx_H2CO)  &
        -k(64)*n(idx_H2CO)

    !d[HCO_dot]/d[CO+]
    pd(14,65) =  &
        -k(65)*n(idx_HCO)  &
        +k(423)*n(idx_H2CO)

    !d[NH3_dot]/d[CO+]
    pd(16,65) =  &
        -k(172)*n(idx_NH3)  &
        -k(687)*n(idx_NH3)

    !d[NO_dot]/d[CO+]
    pd(17,65) =  &
        -k(66)*n(idx_NO)

    !d[CO_dot]/d[CO+]
    pd(19,65) =  &
        +k(46)*n(idx_CH4)  &
        +k(111)*n(idx_H)  &
        +k(118)*n(idx_HCN)  &
        +k(179)*n(idx_NH)  &
        +k(23)*n(idx_C)  &
        +k(164)*n(idx_NH2)  &
        +k(172)*n(idx_NH3)  &
        +k(67)*n(idx_O2)  &
        +k(48)*n(idx_CH)  &
        +k(195)*n(idx_O)  &
        +k(64)*n(idx_H2CO)  &
        +k(107)*n(idx_H2O)  &
        +k(65)*n(idx_HCO)  &
        +k(32)*n(idx_CH2)  &
        +k(204)*n(idx_OH)  &
        +k(66)*n(idx_NO)

    !d[NH2_dot]/d[CO+]
    pd(21,65) =  &
        -k(674)*n(idx_NH2)  &
        -k(164)*n(idx_NH2)  &
        +k(687)*n(idx_NH3)

    !d[CH3_dot]/d[CO+]
    pd(22,65) =  &
        +k(393)*n(idx_CH4)

    !d[CH4_dot]/d[CO+]
    pd(23,65) =  &
        -k(393)*n(idx_CH4)  &
        -k(46)*n(idx_CH4)

    !d[N_dot]/d[CO+]
    pd(24,65) =  &
        +k(690)*n(idx_NH)

    !d[NH_dot]/d[CO+]
    pd(25,65) =  &
        -k(179)*n(idx_NH)  &
        +k(674)*n(idx_NH2)  &
        -k(690)*n(idx_NH)

    !d[CO_DUST_dot]/d[CO+]
    pd(39,65) =  &
        +k(1095)

    !d[HCO+_dot]/d[CO+]
    pd(54,65) =  &
        +k(400)*n(idx_CH)  &
        +k(674)*n(idx_NH2)  &
        +k(393)*n(idx_CH4)  &
        +k(366)*n(idx_CH2)  &
        +k(423)*n(idx_H2CO)  &
        +k(492)*n(idx_H2O)  &
        +k(738)*n(idx_OH)  &
        +k(464)*n(idx_H2)  &
        +k(687)*n(idx_NH3)  &
        +k(690)*n(idx_NH)  &
        +k(65)*n(idx_HCO)

    !d[H+_dot]/d[CO+]
    pd(55,65) =  &
        +k(111)*n(idx_H)

    !d[HOC+_dot]/d[CO+]
    pd(56,65) =  &
        +k(465)*n(idx_H2)

    !d[C+_dot]/d[CO+]
    pd(57,65) =  &
        +k(997)  &
        +k(23)*n(idx_C)

    !d[CH2+_dot]/d[CO+]
    pd(58,65) =  &
        +k(32)*n(idx_CH2)

    !d[CH+_dot]/d[CO+]
    pd(59,65) =  &
        +k(48)*n(idx_CH)

    !d[H2CO+_dot]/d[CO+]
    pd(60,65) =  &
        +k(64)*n(idx_H2CO)

    !d[NH3+_dot]/d[CO+]
    pd(62,65) =  &
        +k(172)*n(idx_NH3)

    !d[NO+_dot]/d[CO+]
    pd(63,65) =  &
        +k(66)*n(idx_NO)

    !d[CO+_dot]/d[CO+]
    pd(65,65) =  &
        -k(465)*n(idx_H2)  &
        -k(269)*n(idx_E)  &
        -k(464)*n(idx_H2)  &
        -k(400)*n(idx_CH)  &
        -k(172)*n(idx_NH3)  &
        -k(1095)  &
        -k(64)*n(idx_H2CO)  &
        -k(204)*n(idx_OH)  &
        -k(366)*n(idx_CH2)  &
        -k(48)*n(idx_CH)  &
        -k(65)*n(idx_HCO)  &
        -k(423)*n(idx_H2CO)  &
        -k(107)*n(idx_H2O)  &
        -k(195)*n(idx_O)  &
        -k(179)*n(idx_NH)  &
        -k(111)*n(idx_H)  &
        -k(674)*n(idx_NH2)  &
        -k(492)*n(idx_H2O)  &
        -k(393)*n(idx_CH4)  &
        -k(690)*n(idx_NH)  &
        -k(46)*n(idx_CH4)  &
        -k(23)*n(idx_C)  &
        -k(66)*n(idx_NO)  &
        -k(687)*n(idx_NH3)  &
        -k(118)*n(idx_HCN)  &
        -k(67)*n(idx_O2)  &
        -k(997)  &
        -k(164)*n(idx_NH2)  &
        -k(738)*n(idx_OH)  &
        -k(32)*n(idx_CH2)

    !d[O2+_dot]/d[CO+]
    pd(67,65) =  &
        +k(67)*n(idx_O2)

    !d[H2O+_dot]/d[CO+]
    pd(68,65) =  &
        +k(107)*n(idx_H2O)

    !d[NH2+_dot]/d[CO+]
    pd(69,65) =  &
        +k(164)*n(idx_NH2)

    !d[O+_dot]/d[CO+]
    pd(70,65) =  &
        +k(195)*n(idx_O)

    !d[OH+_dot]/d[CO+]
    pd(71,65) =  &
        +k(204)*n(idx_OH)

    !d[CH4+_dot]/d[CO+]
    pd(73,65) =  &
        +k(46)*n(idx_CH4)

    !d[HCN+_dot]/d[CO+]
    pd(75,65) =  &
        +k(118)*n(idx_HCN)

    !d[NH+_dot]/d[CO+]
    pd(76,65) =  &
        +k(179)*n(idx_NH)

    !d[E_dot]/d[N2+]
    pd(1,66) =  &
        -k(302)*n(idx_E)

    !d[CH_dot]/d[N2+]
    pd(2,66) =  &
        -k(52)*n(idx_CH)

    !d[O_dot]/d[N2+]
    pd(3,66) =  &
        -k(720)*n(idx_O)  &
        -k(196)*n(idx_O)

    !d[HCN_dot]/d[N2+]
    pd(5,66) =  &
        -k(119)*n(idx_HCN)

    !d[H2_dot]/d[N2+]
    pd(6,66) =  &
        +k(397)*n(idx_CH4)  &
        -k(471)*n(idx_H2)

    !d[C_dot]/d[N2+]
    pd(7,66) =  &
        -k(24)*n(idx_C)

    !d[H_dot]/d[N2+]
    pd(8,66) =  &
        +k(471)*n(idx_H2)  &
        +k(628)*n(idx_H2CO)  &
        +k(398)*n(idx_CH4)

    !d[H2O_dot]/d[N2+]
    pd(9,66) =  &
        -k(499)*n(idx_H2O)  &
        -k(109)*n(idx_H2O)

    !d[OH_dot]/d[N2+]
    pd(10,66) =  &
        +k(499)*n(idx_H2O)  &
        -k(205)*n(idx_OH)

    !d[O2_dot]/d[N2+]
    pd(11,66) =  &
        -k(153)*n(idx_O2)

    !d[CH2_dot]/d[N2+]
    pd(12,66) =  &
        -k(35)*n(idx_CH2)

    !d[H2CO_dot]/d[N2+]
    pd(13,66) =  &
        -k(150)*n(idx_H2CO)  &
        -k(628)*n(idx_H2CO)

    !d[HCO_dot]/d[N2+]
    pd(14,66) =  &
        -k(629)*n(idx_HCO)  &
        -k(151)*n(idx_HCO)

    !d[MG_dot]/d[N2+]
    pd(15,66) =  &
        -k(132)*n(idx_MG)

    !d[NH3_dot]/d[N2+]
    pd(16,66) =  &
        -k(176)*n(idx_NH3)

    !d[NO_dot]/d[N2+]
    pd(17,66) =  &
        -k(152)*n(idx_NO)

    !d[CN_dot]/d[N2+]
    pd(18,66) =  &
        -k(63)*n(idx_CN)

    !d[CO_dot]/d[N2+]
    pd(19,66) =  &
        +k(629)*n(idx_HCO)  &
        -k(68)*n(idx_CO)

    !d[N2_dot]/d[N2+]
    pd(20,66) =  &
        +k(154)*n(idx_N)  &
        +k(152)*n(idx_NO)  &
        +k(35)*n(idx_CH2)  &
        +k(52)*n(idx_CH)  &
        +k(166)*n(idx_NH2)  &
        +k(132)*n(idx_MG)  &
        +k(628)*n(idx_H2CO)  &
        +k(196)*n(idx_O)  &
        +k(150)*n(idx_H2CO)  &
        +k(24)*n(idx_C)  &
        +k(397)*n(idx_CH4)  &
        +k(68)*n(idx_CO)  &
        +k(180)*n(idx_NH)  &
        +k(119)*n(idx_HCN)  &
        +k(109)*n(idx_H2O)  &
        +k(176)*n(idx_NH3)  &
        +k(63)*n(idx_CN)  &
        +k(153)*n(idx_O2)  &
        +k(151)*n(idx_HCO)  &
        +k(205)*n(idx_OH)  &
        +k(398)*n(idx_CH4)

    !d[NH2_dot]/d[N2+]
    pd(21,66) =  &
        -k(166)*n(idx_NH2)

    !d[CH4_dot]/d[N2+]
    pd(23,66) =  &
        -k(398)*n(idx_CH4)  &
        -k(397)*n(idx_CH4)

    !d[N_dot]/d[N2+]
    pd(24,66) =  &
        +k(720)*n(idx_O)  &
        -k(154)*n(idx_N)  &
        +2.d0*k(302)*n(idx_E)

    !d[NH_dot]/d[N2+]
    pd(25,66) =  &
        -k(180)*n(idx_NH)

    !d[N2_DUST_dot]/d[N2+]
    pd(43,66) =  &
        +k(1091)

    !d[HCO+_dot]/d[N2+]
    pd(54,66) =  &
        +k(151)*n(idx_HCO)  &
        +k(628)*n(idx_H2CO)

    !d[C+_dot]/d[N2+]
    pd(57,66) =  &
        +k(24)*n(idx_C)

    !d[CH2+_dot]/d[N2+]
    pd(58,66) =  &
        +k(397)*n(idx_CH4)  &
        +k(35)*n(idx_CH2)

    !d[CH+_dot]/d[N2+]
    pd(59,66) =  &
        +k(52)*n(idx_CH)

    !d[H2CO+_dot]/d[N2+]
    pd(60,66) =  &
        +k(150)*n(idx_H2CO)

    !d[MG+_dot]/d[N2+]
    pd(61,66) =  &
        +k(132)*n(idx_MG)

    !d[NH3+_dot]/d[N2+]
    pd(62,66) =  &
        +k(176)*n(idx_NH3)

    !d[NO+_dot]/d[N2+]
    pd(63,66) =  &
        +k(152)*n(idx_NO)  &
        +k(720)*n(idx_O)

    !d[CN+_dot]/d[N2+]
    pd(64,66) =  &
        +k(63)*n(idx_CN)

    !d[CO+_dot]/d[N2+]
    pd(65,66) =  &
        +k(68)*n(idx_CO)

    !d[N2+_dot]/d[N2+]
    pd(66,66) =  &
        -k(720)*n(idx_O)  &
        -k(151)*n(idx_HCO)  &
        -k(176)*n(idx_NH3)  &
        -k(302)*n(idx_E)  &
        -k(205)*n(idx_OH)  &
        -k(397)*n(idx_CH4)  &
        -k(629)*n(idx_HCO)  &
        -k(119)*n(idx_HCN)  &
        -k(109)*n(idx_H2O)  &
        -k(398)*n(idx_CH4)  &
        -k(132)*n(idx_MG)  &
        -k(150)*n(idx_H2CO)  &
        -k(154)*n(idx_N)  &
        -k(471)*n(idx_H2)  &
        -k(35)*n(idx_CH2)  &
        -k(24)*n(idx_C)  &
        -k(196)*n(idx_O)  &
        -k(180)*n(idx_NH)  &
        -k(63)*n(idx_CN)  &
        -k(499)*n(idx_H2O)  &
        -k(68)*n(idx_CO)  &
        -k(1091)  &
        -k(152)*n(idx_NO)  &
        -k(153)*n(idx_O2)  &
        -k(628)*n(idx_H2CO)  &
        -k(166)*n(idx_NH2)  &
        -k(52)*n(idx_CH)

    !d[O2+_dot]/d[N2+]
    pd(67,66) =  &
        +k(153)*n(idx_O2)

    !d[H2O+_dot]/d[N2+]
    pd(68,66) =  &
        +k(109)*n(idx_H2O)

    !d[NH2+_dot]/d[N2+]
    pd(69,66) =  &
        +k(166)*n(idx_NH2)

    !d[O+_dot]/d[N2+]
    pd(70,66) =  &
        +k(196)*n(idx_O)

    !d[OH+_dot]/d[N2+]
    pd(71,66) =  &
        +k(205)*n(idx_OH)

    !d[CH3+_dot]/d[N2+]
    pd(72,66) =  &
        +k(398)*n(idx_CH4)

    !d[N+_dot]/d[N2+]
    pd(74,66) =  &
        +k(154)*n(idx_N)

    !d[HCN+_dot]/d[N2+]
    pd(75,66) =  &
        +k(119)*n(idx_HCN)

    !d[NH+_dot]/d[N2+]
    pd(76,66) =  &
        +k(180)*n(idx_NH)

    !d[N2H+_dot]/d[N2+]
    pd(87,66) =  &
        +k(471)*n(idx_H2)  &
        +k(499)*n(idx_H2O)  &
        +k(629)*n(idx_HCO)

    !d[E_dot]/d[O2+]
    pd(1,67) =  &
        -k(311)*n(idx_E)

    !d[CH_dot]/d[O2+]
    pd(2,67) =  &
        -k(415)*n(idx_CH)  &
        -k(55)*n(idx_CH)

    !d[O_dot]/d[O2+]
    pd(3,67) =  &
        +k(337)*n(idx_C)  &
        +k(699)*n(idx_NH)  &
        +k(640)*n(idx_N)  &
        +k(415)*n(idx_CH)  &
        +2.d0*k(311)*n(idx_E)  &
        +k(379)*n(idx_CH2)  &
        +k(1031)

    !d[C_dot]/d[O2+]
    pd(7,67) =  &
        -k(337)*n(idx_C)  &
        -k(25)*n(idx_C)

    !d[H_dot]/d[O2+]
    pd(8,67) =  &
        +k(715)*n(idx_CH3OH)  &
        +k(481)*n(idx_H2CO)

    !d[O2_dot]/d[O2+]
    pd(11,67) =  &
        +k(715)*n(idx_CH3OH)  &
        +k(134)*n(idx_MG)  &
        +k(121)*n(idx_HCO)  &
        +k(25)*n(idx_C)  &
        +k(177)*n(idx_NH3)  &
        +k(55)*n(idx_CH)  &
        +k(184)*n(idx_NO)  &
        +k(101)*n(idx_H2CO)  &
        +k(167)*n(idx_NH2)  &
        +k(38)*n(idx_CH2)  &
        +k(481)*n(idx_H2CO)

    !d[CH2_dot]/d[O2+]
    pd(12,67) =  &
        -k(38)*n(idx_CH2)  &
        -k(379)*n(idx_CH2)

    !d[H2CO_dot]/d[O2+]
    pd(13,67) =  &
        -k(481)*n(idx_H2CO)  &
        -k(101)*n(idx_H2CO)

    !d[HCO_dot]/d[O2+]
    pd(14,67) =  &
        -k(121)*n(idx_HCO)  &
        -k(555)*n(idx_HCO)

    !d[MG_dot]/d[O2+]
    pd(15,67) =  &
        -k(134)*n(idx_MG)

    !d[NH3_dot]/d[O2+]
    pd(16,67) =  &
        -k(177)*n(idx_NH3)

    !d[NO_dot]/d[O2+]
    pd(17,67) =  &
        -k(184)*n(idx_NO)

    !d[CO_dot]/d[O2+]
    pd(19,67) =  &
        +k(555)*n(idx_HCO)

    !d[NH2_dot]/d[O2+]
    pd(21,67) =  &
        -k(167)*n(idx_NH2)

    !d[N_dot]/d[O2+]
    pd(24,67) =  &
        -k(640)*n(idx_N)

    !d[NH_dot]/d[O2+]
    pd(25,67) =  &
        -k(699)*n(idx_NH)

    !d[CH3OH_dot]/d[O2+]
    pd(28,67) =  &
        -k(715)*n(idx_CH3OH)

    !d[O2_DUST_dot]/d[O2+]
    pd(46,67) =  &
        +k(1090)

    !d[HCO+_dot]/d[O2+]
    pd(54,67) =  &
        +k(121)*n(idx_HCO)  &
        +k(415)*n(idx_CH)  &
        +k(481)*n(idx_H2CO)

    !d[C+_dot]/d[O2+]
    pd(57,67) =  &
        +k(25)*n(idx_C)

    !d[CH2+_dot]/d[O2+]
    pd(58,67) =  &
        +k(38)*n(idx_CH2)

    !d[CH+_dot]/d[O2+]
    pd(59,67) =  &
        +k(55)*n(idx_CH)

    !d[H2CO+_dot]/d[O2+]
    pd(60,67) =  &
        +k(379)*n(idx_CH2)  &
        +k(101)*n(idx_H2CO)

    !d[MG+_dot]/d[O2+]
    pd(61,67) =  &
        +k(134)*n(idx_MG)

    !d[NH3+_dot]/d[O2+]
    pd(62,67) =  &
        +k(177)*n(idx_NH3)

    !d[NO+_dot]/d[O2+]
    pd(63,67) =  &
        +k(640)*n(idx_N)  &
        +k(184)*n(idx_NO)

    !d[CO+_dot]/d[O2+]
    pd(65,67) =  &
        +k(337)*n(idx_C)

    !d[O2+_dot]/d[O2+]
    pd(67,67) =  &
        -k(1031)  &
        -k(55)*n(idx_CH)  &
        -k(337)*n(idx_C)  &
        -k(167)*n(idx_NH2)  &
        -k(1090)  &
        -k(101)*n(idx_H2CO)  &
        -k(699)*n(idx_NH)  &
        -k(311)*n(idx_E)  &
        -k(184)*n(idx_NO)  &
        -k(640)*n(idx_N)  &
        -k(555)*n(idx_HCO)  &
        -k(121)*n(idx_HCO)  &
        -k(38)*n(idx_CH2)  &
        -k(415)*n(idx_CH)  &
        -k(715)*n(idx_CH3OH)  &
        -k(481)*n(idx_H2CO)  &
        -k(177)*n(idx_NH3)  &
        -k(379)*n(idx_CH2)  &
        -k(25)*n(idx_C)  &
        -k(134)*n(idx_MG)

    !d[NH2+_dot]/d[O2+]
    pd(69,67) =  &
        +k(167)*n(idx_NH2)

    !d[O+_dot]/d[O2+]
    pd(70,67) =  &
        +k(1031)

    !d[HNO+_dot]/d[O2+]
    pd(79,67) =  &
        +k(699)*n(idx_NH)

    !d[H3CO+_dot]/d[O2+]
    pd(82,67) =  &
        +k(715)*n(idx_CH3OH)

    !d[O2H+_dot]/d[O2+]
    pd(88,67) =  &
        +k(555)*n(idx_HCO)

    !d[E_dot]/d[H2O+]
    pd(1,68) =  &
        -k(277)*n(idx_E)  &
        -k(278)*n(idx_E)  &
        -k(279)*n(idx_E)

    !d[CH_dot]/d[H2O+]
    pd(2,68) =  &
        -k(50)*n(idx_CH)  &
        -k(402)*n(idx_CH)

    !d[O_dot]/d[H2O+]
    pd(3,68) =  &
        -k(718)*n(idx_O)  &
        +k(278)*n(idx_E)  &
        +k(739)*n(idx_OH)  &
        +k(277)*n(idx_E)

    !d[HNC_dot]/d[H2O+]
    pd(4,68) =  &
        -k(489)*n(idx_HNC)

    !d[HCN_dot]/d[H2O+]
    pd(5,68) =  &
        -k(486)*n(idx_HCN)

    !d[H2_dot]/d[H2O+]
    pd(6,68) =  &
        +k(718)*n(idx_O)  &
        +k(637)*n(idx_N)  &
        -k(466)*n(idx_H2)  &
        +k(277)*n(idx_E)

    !d[C_dot]/d[H2O+]
    pd(7,68) =  &
        -k(329)*n(idx_C)

    !d[H_dot]/d[H2O+]
    pd(8,68) =  &
        +k(1006)  &
        +k(636)*n(idx_N)  &
        +2.d0*k(278)*n(idx_E)  &
        +k(279)*n(idx_E)  &
        +k(466)*n(idx_H2)

    !d[H2O_dot]/d[H2O+]
    pd(9,68) =  &
        +k(50)*n(idx_CH)  &
        +k(165)*n(idx_NH2)  &
        +k(103)*n(idx_HCO)  &
        +k(106)*n(idx_O2)  &
        +k(34)*n(idx_CH2)  &
        +k(105)*n(idx_NO)  &
        +k(102)*n(idx_H2CO)  &
        +k(174)*n(idx_NH3)  &
        -k(485)*n(idx_H2O)  &
        +k(104)*n(idx_MG)

    !d[OH_dot]/d[H2O+]
    pd(10,68) =  &
        +k(402)*n(idx_CH)  &
        +k(676)*n(idx_NH2)  &
        +k(489)*n(idx_HNC)  &
        +k(484)*n(idx_H2CO)  &
        +k(485)*n(idx_H2O)  &
        +k(483)*n(idx_CO)  &
        +k(279)*n(idx_E)  &
        +k(488)*n(idx_HCO)  &
        +k(329)*n(idx_C)  &
        -k(739)*n(idx_OH)  &
        +k(368)*n(idx_CH2)  &
        +k(486)*n(idx_HCN)

    !d[O2_dot]/d[H2O+]
    pd(11,68) =  &
        -k(106)*n(idx_O2)

    !d[CH2_dot]/d[H2O+]
    pd(12,68) =  &
        -k(34)*n(idx_CH2)  &
        -k(368)*n(idx_CH2)

    !d[H2CO_dot]/d[H2O+]
    pd(13,68) =  &
        -k(102)*n(idx_H2CO)  &
        -k(484)*n(idx_H2CO)

    !d[HCO_dot]/d[H2O+]
    pd(14,68) =  &
        -k(103)*n(idx_HCO)  &
        -k(488)*n(idx_HCO)  &
        -k(487)*n(idx_HCO)

    !d[MG_dot]/d[H2O+]
    pd(15,68) =  &
        -k(104)*n(idx_MG)

    !d[NH3_dot]/d[H2O+]
    pd(16,68) =  &
        -k(174)*n(idx_NH3)

    !d[NO_dot]/d[H2O+]
    pd(17,68) =  &
        -k(105)*n(idx_NO)

    !d[CO_dot]/d[H2O+]
    pd(19,68) =  &
        -k(483)*n(idx_CO)  &
        +k(487)*n(idx_HCO)

    !d[NH2_dot]/d[H2O+]
    pd(21,68) =  &
        -k(676)*n(idx_NH2)  &
        -k(165)*n(idx_NH2)

    !d[CH3_dot]/d[H2O+]
    pd(22,68) =  &
        +k(395)*n(idx_CH4)

    !d[CH4_dot]/d[H2O+]
    pd(23,68) =  &
        -k(395)*n(idx_CH4)

    !d[N_dot]/d[H2O+]
    pd(24,68) =  &
        -k(637)*n(idx_N)  &
        -k(636)*n(idx_N)  &
        +k(692)*n(idx_NH)

    !d[NH_dot]/d[H2O+]
    pd(25,68) =  &
        -k(692)*n(idx_NH)

    !d[H2O_DUST_dot]/d[H2O+]
    pd(40,68) =  &
        +k(1100)

    !d[HCO+_dot]/d[H2O+]
    pd(54,68) =  &
        +k(103)*n(idx_HCO)  &
        +k(483)*n(idx_CO)

    !d[CH2+_dot]/d[H2O+]
    pd(58,68) =  &
        +k(402)*n(idx_CH)  &
        +k(34)*n(idx_CH2)

    !d[CH+_dot]/d[H2O+]
    pd(59,68) =  &
        +k(50)*n(idx_CH)  &
        +k(329)*n(idx_C)

    !d[H2CO+_dot]/d[H2O+]
    pd(60,68) =  &
        +k(102)*n(idx_H2CO)  &
        +k(488)*n(idx_HCO)

    !d[MG+_dot]/d[H2O+]
    pd(61,68) =  &
        +k(104)*n(idx_MG)

    !d[NH3+_dot]/d[H2O+]
    pd(62,68) =  &
        +k(174)*n(idx_NH3)  &
        +k(676)*n(idx_NH2)

    !d[NO+_dot]/d[H2O+]
    pd(63,68) =  &
        +k(105)*n(idx_NO)  &
        +k(637)*n(idx_N)

    !d[O2+_dot]/d[H2O+]
    pd(67,68) =  &
        +k(718)*n(idx_O)  &
        +k(106)*n(idx_O2)

    !d[H2O+_dot]/d[H2O+]
    pd(68,68) =  &
        -k(50)*n(idx_CH)  &
        -k(174)*n(idx_NH3)  &
        -k(486)*n(idx_HCN)  &
        -k(402)*n(idx_CH)  &
        -k(34)*n(idx_CH2)  &
        -k(165)*n(idx_NH2)  &
        -k(484)*n(idx_H2CO)  &
        -k(483)*n(idx_CO)  &
        -k(488)*n(idx_HCO)  &
        -k(102)*n(idx_H2CO)  &
        -k(487)*n(idx_HCO)  &
        -k(1006)  &
        -k(277)*n(idx_E)  &
        -k(329)*n(idx_C)  &
        -k(105)*n(idx_NO)  &
        -k(637)*n(idx_N)  &
        -k(485)*n(idx_H2O)  &
        -k(739)*n(idx_OH)  &
        -k(106)*n(idx_O2)  &
        -k(489)*n(idx_HNC)  &
        -k(104)*n(idx_MG)  &
        -k(103)*n(idx_HCO)  &
        -k(636)*n(idx_N)  &
        -k(718)*n(idx_O)  &
        -k(676)*n(idx_NH2)  &
        -k(278)*n(idx_E)  &
        -k(1100)  &
        -k(395)*n(idx_CH4)  &
        -k(279)*n(idx_E)  &
        -k(466)*n(idx_H2)  &
        -k(692)*n(idx_NH)  &
        -k(368)*n(idx_CH2)

    !d[NH2+_dot]/d[H2O+]
    pd(69,68) =  &
        +k(165)*n(idx_NH2)

    !d[OH+_dot]/d[H2O+]
    pd(71,68) =  &
        +k(1006)

    !d[CH3+_dot]/d[H2O+]
    pd(72,68) =  &
        +k(368)*n(idx_CH2)

    !d[HNO+_dot]/d[H2O+]
    pd(79,68) =  &
        +k(636)*n(idx_N)

    !d[H3CO+_dot]/d[H2O+]
    pd(82,68) =  &
        +k(484)*n(idx_H2CO)

    !d[H3O+_dot]/d[H2O+]
    pd(83,68) =  &
        +k(739)*n(idx_OH)  &
        +k(485)*n(idx_H2O)  &
        +k(692)*n(idx_NH)  &
        +k(395)*n(idx_CH4)  &
        +k(487)*n(idx_HCO)  &
        +k(466)*n(idx_H2)

    !d[HCNH+_dot]/d[H2O+]
    pd(84,68) =  &
        +k(489)*n(idx_HNC)  &
        +k(486)*n(idx_HCN)

    !d[E_dot]/d[NH2+]
    pd(1,69) =  &
        -k(306)*n(idx_E)  &
        -k(307)*n(idx_E)

    !d[CH_dot]/d[NH2+]
    pd(2,69) =  &
        -k(413)*n(idx_CH)  &
        -k(53)*n(idx_CH)

    !d[O_dot]/d[NH2+]
    pd(3,69) =  &
        -k(722)*n(idx_O)  &
        +k(672)*n(idx_O2)

    !d[HNC_dot]/d[NH2+]
    pd(4,69) =  &
        -k(670)*n(idx_HNC)

    !d[HCN_dot]/d[NH2+]
    pd(5,69) =  &
        -k(668)*n(idx_HCN)

    !d[H2_dot]/d[NH2+]
    pd(6,69) =  &
        -k(474)*n(idx_H2)

    !d[H_dot]/d[NH2+]
    pd(8,69) =  &
        +2.d0*k(306)*n(idx_E)  &
        +k(307)*n(idx_E)  &
        +k(639)*n(idx_N)  &
        +k(722)*n(idx_O)  &
        +k(474)*n(idx_H2)

    !d[H2O_dot]/d[NH2+]
    pd(9,69) =  &
        -k(666)*n(idx_H2O)  &
        -k(667)*n(idx_H2O)

    !d[OH_dot]/d[NH2+]
    pd(10,69) =  &
        +k(673)*n(idx_O2)  &
        +k(667)*n(idx_H2O)

    !d[O2_dot]/d[NH2+]
    pd(11,69) =  &
        -k(672)*n(idx_O2)  &
        -k(673)*n(idx_O2)

    !d[CH2_dot]/d[NH2+]
    pd(12,69) =  &
        -k(36)*n(idx_CH2)  &
        -k(377)*n(idx_CH2)

    !d[H2CO_dot]/d[NH2+]
    pd(13,69) =  &
        -k(665)*n(idx_H2CO)  &
        -k(664)*n(idx_H2CO)

    !d[HCO_dot]/d[NH2+]
    pd(14,69) =  &
        -k(160)*n(idx_HCO)  &
        +k(665)*n(idx_H2CO)  &
        -k(669)*n(idx_HCO)

    !d[NH3_dot]/d[NH2+]
    pd(16,69) =  &
        -k(161)*n(idx_NH3)

    !d[NO_dot]/d[NH2+]
    pd(17,69) =  &
        -k(162)*n(idx_NO)

    !d[NH2_dot]/d[NH2+]
    pd(21,69) =  &
        +k(160)*n(idx_HCO)  &
        +k(36)*n(idx_CH2)  &
        +k(162)*n(idx_NO)  &
        +k(53)*n(idx_CH)  &
        -k(671)*n(idx_NH2)  &
        +k(161)*n(idx_NH3)

    !d[N_dot]/d[NH2+]
    pd(24,69) =  &
        +k(306)*n(idx_E)  &
        +k(697)*n(idx_NH)  &
        -k(639)*n(idx_N)

    !d[NH_dot]/d[NH2+]
    pd(25,69) =  &
        -k(697)*n(idx_NH)  &
        +k(668)*n(idx_HCN)  &
        +k(413)*n(idx_CH)  &
        +k(377)*n(idx_CH2)  &
        +k(307)*n(idx_E)  &
        +k(664)*n(idx_H2CO)  &
        +k(670)*n(idx_HNC)  &
        +k(671)*n(idx_NH2)  &
        +k(666)*n(idx_H2O)  &
        +k(669)*n(idx_HCO)

    !d[NH3_DUST_dot]/d[NH2+]
    pd(45,69) =  &
        +k(1099)

    !d[HCO+_dot]/d[NH2+]
    pd(54,69) =  &
        +k(160)*n(idx_HCO)

    !d[CH2+_dot]/d[NH2+]
    pd(58,69) =  &
        +k(36)*n(idx_CH2)  &
        +k(413)*n(idx_CH)

    !d[CH+_dot]/d[NH2+]
    pd(59,69) =  &
        +k(53)*n(idx_CH)

    !d[H2CO+_dot]/d[NH2+]
    pd(60,69) =  &
        +k(669)*n(idx_HCO)

    !d[NH3+_dot]/d[NH2+]
    pd(62,69) =  &
        +k(665)*n(idx_H2CO)  &
        +k(697)*n(idx_NH)  &
        +k(667)*n(idx_H2O)  &
        +k(474)*n(idx_H2)  &
        +k(161)*n(idx_NH3)  &
        +k(671)*n(idx_NH2)

    !d[NO+_dot]/d[NH2+]
    pd(63,69) =  &
        +k(162)*n(idx_NO)

    !d[NH2+_dot]/d[NH2+]
    pd(69,69) =  &
        -k(474)*n(idx_H2)  &
        -k(306)*n(idx_E)  &
        -k(1099)  &
        -k(670)*n(idx_HNC)  &
        -k(671)*n(idx_NH2)  &
        -k(36)*n(idx_CH2)  &
        -k(722)*n(idx_O)  &
        -k(53)*n(idx_CH)  &
        -k(697)*n(idx_NH)  &
        -k(666)*n(idx_H2O)  &
        -k(665)*n(idx_H2CO)  &
        -k(669)*n(idx_HCO)  &
        -k(160)*n(idx_HCO)  &
        -k(673)*n(idx_O2)  &
        -k(377)*n(idx_CH2)  &
        -k(672)*n(idx_O2)  &
        -k(639)*n(idx_N)  &
        -k(413)*n(idx_CH)  &
        -k(162)*n(idx_NO)  &
        -k(307)*n(idx_E)  &
        -k(664)*n(idx_H2CO)  &
        -k(667)*n(idx_H2O)  &
        -k(668)*n(idx_HCN)  &
        -k(161)*n(idx_NH3)

    !d[CH3+_dot]/d[NH2+]
    pd(72,69) =  &
        +k(377)*n(idx_CH2)

    !d[HNO+_dot]/d[NH2+]
    pd(79,69) =  &
        +k(722)*n(idx_O)  &
        +k(673)*n(idx_O2)

    !d[H2NO+_dot]/d[NH2+]
    pd(80,69) =  &
        +k(672)*n(idx_O2)

    !d[H3CO+_dot]/d[NH2+]
    pd(82,69) =  &
        +k(664)*n(idx_H2CO)

    !d[H3O+_dot]/d[NH2+]
    pd(83,69) =  &
        +k(666)*n(idx_H2O)

    !d[HCNH+_dot]/d[NH2+]
    pd(84,69) =  &
        +k(668)*n(idx_HCN)  &
        +k(670)*n(idx_HNC)

    !d[N2H+_dot]/d[NH2+]
    pd(87,69) =  &
        +k(639)*n(idx_N)

    !d[E_dot]/d[O+]
    pd(1,70) =  &
        -k(1063)*n(idx_E)

    !d[CH_dot]/d[O+]
    pd(2,70) =  &
        +k(710)*n(idx_HCN)  &
        -k(54)*n(idx_CH)  &
        -k(414)*n(idx_CH)

    !d[O_dot]/d[O+]
    pd(3,70) =  &
        +k(37)*n(idx_CH2)  &
        +k(186)*n(idx_CO)  &
        +k(115)*n(idx_H)  &
        +k(188)*n(idx_H2O)  &
        +k(181)*n(idx_NH)  &
        +k(54)*n(idx_CH)  &
        +k(185)*n(idx_CH4)  &
        +k(191)*n(idx_NH3)  &
        +k(189)*n(idx_HCO)  &
        +k(190)*n(idx_NH2)  &
        +k(1063)*n(idx_E)  &
        +k(187)*n(idx_H2CO)  &
        +k(192)*n(idx_O2)  &
        +k(193)*n(idx_OH)

    !d[HCN_dot]/d[O+]
    pd(5,70) =  &
        -k(710)*n(idx_HCN)  &
        -k(709)*n(idx_HCN)

    !d[H2_dot]/d[O+]
    pd(6,70) =  &
        -k(475)*n(idx_H2)

    !d[C_dot]/d[O+]
    pd(7,70) =  &
        -k(1043)*n(idx_C)  &
        +k(706)*n(idx_CN)

    !d[H_dot]/d[O+]
    pd(8,70) =  &
        +k(414)*n(idx_CH)  &
        +k(698)*n(idx_NH)  &
        +k(714)*n(idx_OH)  &
        -k(115)*n(idx_H)  &
        +k(475)*n(idx_H2)

    !d[H2O_dot]/d[O+]
    pd(9,70) =  &
        -k(188)*n(idx_H2O)  &
        +k(703)*n(idx_CH3OH)

    !d[OH_dot]/d[O+]
    pd(10,70) =  &
        +k(704)*n(idx_CH3OH)  &
        -k(193)*n(idx_OH)  &
        +k(708)*n(idx_H2CO)  &
        +k(705)*n(idx_CH4)  &
        -k(714)*n(idx_OH)

    !d[O2_dot]/d[O+]
    pd(11,70) =  &
        -k(192)*n(idx_O2)  &
        +k(713)*n(idx_NO2)

    !d[CH2_dot]/d[O+]
    pd(12,70) =  &
        -k(37)*n(idx_CH2)

    !d[H2CO_dot]/d[O+]
    pd(13,70) =  &
        -k(708)*n(idx_H2CO)  &
        -k(187)*n(idx_H2CO)

    !d[HCO_dot]/d[O+]
    pd(14,70) =  &
        -k(189)*n(idx_HCO)  &
        -k(711)*n(idx_HCO)

    !d[NH3_dot]/d[O+]
    pd(16,70) =  &
        -k(191)*n(idx_NH3)

    !d[CN_dot]/d[O+]
    pd(18,70) =  &
        -k(706)*n(idx_CN)

    !d[CO_dot]/d[O+]
    pd(19,70) =  &
        +k(711)*n(idx_HCO)  &
        -k(186)*n(idx_CO)  &
        +k(707)*n(idx_CO2)

    !d[N2_dot]/d[O+]
    pd(20,70) =  &
        -k(712)*n(idx_N2)

    !d[NH2_dot]/d[O+]
    pd(21,70) =  &
        -k(190)*n(idx_NH2)

    !d[CH4_dot]/d[O+]
    pd(23,70) =  &
        -k(185)*n(idx_CH4)  &
        -k(705)*n(idx_CH4)

    !d[N_dot]/d[O+]
    pd(24,70) =  &
        +k(709)*n(idx_HCN)  &
        +k(712)*n(idx_N2)

    !d[NH_dot]/d[O+]
    pd(25,70) =  &
        -k(698)*n(idx_NH)  &
        -k(181)*n(idx_NH)

    !d[CH3OH_dot]/d[O+]
    pd(28,70) =  &
        -k(703)*n(idx_CH3OH)  &
        -k(704)*n(idx_CH3OH)

    !d[CO2_dot]/d[O+]
    pd(29,70) =  &
        -k(707)*n(idx_CO2)

    !d[NO2_dot]/d[O+]
    pd(32,70) =  &
        -k(713)*n(idx_NO2)

    !d[H2O_DUST_dot]/d[O+]
    pd(40,70) =  &
        +k(1089)

    !d[HCO+_dot]/d[O+]
    pd(54,70) =  &
        +k(708)*n(idx_H2CO)  &
        +k(189)*n(idx_HCO)  &
        +k(709)*n(idx_HCN)

    !d[H+_dot]/d[O+]
    pd(55,70) =  &
        +k(115)*n(idx_H)

    !d[CH2+_dot]/d[O+]
    pd(58,70) =  &
        +k(37)*n(idx_CH2)

    !d[CH+_dot]/d[O+]
    pd(59,70) =  &
        +k(54)*n(idx_CH)

    !d[H2CO+_dot]/d[O+]
    pd(60,70) =  &
        +k(703)*n(idx_CH3OH)  &
        +k(187)*n(idx_H2CO)

    !d[NH3+_dot]/d[O+]
    pd(62,70) =  &
        +k(191)*n(idx_NH3)

    !d[NO+_dot]/d[O+]
    pd(63,70) =  &
        +k(712)*n(idx_N2)  &
        +k(698)*n(idx_NH)  &
        +k(706)*n(idx_CN)  &
        +k(713)*n(idx_NO2)  &
        +k(710)*n(idx_HCN)

    !d[CO+_dot]/d[O+]
    pd(65,70) =  &
        +k(1043)*n(idx_C)  &
        +k(414)*n(idx_CH)  &
        +k(186)*n(idx_CO)

    !d[O2+_dot]/d[O+]
    pd(67,70) =  &
        +k(707)*n(idx_CO2)  &
        +k(714)*n(idx_OH)  &
        +k(192)*n(idx_O2)

    !d[H2O+_dot]/d[O+]
    pd(68,70) =  &
        +k(188)*n(idx_H2O)

    !d[NH2+_dot]/d[O+]
    pd(69,70) =  &
        +k(190)*n(idx_NH2)

    !d[O+_dot]/d[O+]
    pd(70,70) =  &
        -k(710)*n(idx_HCN)  &
        -k(185)*n(idx_CH4)  &
        -k(1089)  &
        -k(709)*n(idx_HCN)  &
        -k(189)*n(idx_HCO)  &
        -k(37)*n(idx_CH2)  &
        -k(191)*n(idx_NH3)  &
        -k(703)*n(idx_CH3OH)  &
        -k(704)*n(idx_CH3OH)  &
        -k(475)*n(idx_H2)  &
        -k(193)*n(idx_OH)  &
        -k(698)*n(idx_NH)  &
        -k(705)*n(idx_CH4)  &
        -k(188)*n(idx_H2O)  &
        -k(187)*n(idx_H2CO)  &
        -k(190)*n(idx_NH2)  &
        -k(192)*n(idx_O2)  &
        -k(708)*n(idx_H2CO)  &
        -k(186)*n(idx_CO)  &
        -k(712)*n(idx_N2)  &
        -k(713)*n(idx_NO2)  &
        -k(115)*n(idx_H)  &
        -k(1063)*n(idx_E)  &
        -k(714)*n(idx_OH)  &
        -k(181)*n(idx_NH)  &
        -k(707)*n(idx_CO2)  &
        -k(414)*n(idx_CH)  &
        -k(54)*n(idx_CH)  &
        -k(711)*n(idx_HCO)  &
        -k(1043)*n(idx_C)  &
        -k(706)*n(idx_CN)

    !d[OH+_dot]/d[O+]
    pd(71,70) =  &
        +k(711)*n(idx_HCO)  &
        +k(475)*n(idx_H2)  &
        +k(193)*n(idx_OH)

    !d[CH3+_dot]/d[O+]
    pd(72,70) =  &
        +k(705)*n(idx_CH4)

    !d[CH4+_dot]/d[O+]
    pd(73,70) =  &
        +k(185)*n(idx_CH4)

    !d[NH+_dot]/d[O+]
    pd(76,70) =  &
        +k(181)*n(idx_NH)

    !d[H3CO+_dot]/d[O+]
    pd(82,70) =  &
        +k(704)*n(idx_CH3OH)

    !d[E_dot]/d[OH+]
    pd(1,71) =  &
        -k(313)*n(idx_E)

    !d[CH_dot]/d[OH+]
    pd(2,71) =  &
        -k(56)*n(idx_CH)  &
        -k(417)*n(idx_CH)

    !d[O_dot]/d[OH+]
    pd(3,71) =  &
        +k(729)*n(idx_H2CO)  &
        +k(735)*n(idx_N2)  &
        +k(730)*n(idx_H2O)  &
        +k(417)*n(idx_CH)  &
        +k(734)*n(idx_HNC)  &
        +k(736)*n(idx_NO)  &
        +k(733)*n(idx_HCO)  &
        +k(701)*n(idx_NH)  &
        -k(725)*n(idx_O)  &
        +k(313)*n(idx_E)  &
        +k(728)*n(idx_CO)  &
        +k(381)*n(idx_CH2)  &
        +k(731)*n(idx_HCN)  &
        +k(686)*n(idx_NH2)  &
        +k(726)*n(idx_CN)  &
        +k(727)*n(idx_CO2)  &
        +k(339)*n(idx_C)  &
        +k(737)*n(idx_OH)

    !d[HNC_dot]/d[OH+]
    pd(4,71) =  &
        -k(734)*n(idx_HNC)

    !d[HCN_dot]/d[OH+]
    pd(5,71) =  &
        -k(731)*n(idx_HCN)

    !d[H2_dot]/d[OH+]
    pd(6,71) =  &
        -k(477)*n(idx_H2)

    !d[C_dot]/d[OH+]
    pd(7,71) =  &
        -k(339)*n(idx_C)

    !d[H_dot]/d[OH+]
    pd(8,71) =  &
        +k(1037)  &
        +k(641)*n(idx_N)  &
        +k(313)*n(idx_E)  &
        +k(725)*n(idx_O)  &
        +k(477)*n(idx_H2)

    !d[H2O_dot]/d[OH+]
    pd(9,71) =  &
        -k(730)*n(idx_H2O)  &
        -k(198)*n(idx_H2O)

    !d[OH_dot]/d[OH+]
    pd(10,71) =  &
        +k(168)*n(idx_NH2)  &
        +k(197)*n(idx_H2CO)  &
        -k(737)*n(idx_OH)  &
        +k(198)*n(idx_H2O)  &
        +k(56)*n(idx_CH)  &
        +k(199)*n(idx_HCO)  &
        +k(39)*n(idx_CH2)  &
        +k(200)*n(idx_NH3)  &
        +k(201)*n(idx_NO)  &
        +k(202)*n(idx_O2)

    !d[O2_dot]/d[OH+]
    pd(11,71) =  &
        -k(202)*n(idx_O2)

    !d[CH2_dot]/d[OH+]
    pd(12,71) =  &
        -k(39)*n(idx_CH2)  &
        +k(399)*n(idx_CH4)  &
        -k(381)*n(idx_CH2)

    !d[H2CO_dot]/d[OH+]
    pd(13,71) =  &
        -k(729)*n(idx_H2CO)  &
        -k(197)*n(idx_H2CO)

    !d[HCO_dot]/d[OH+]
    pd(14,71) =  &
        -k(199)*n(idx_HCO)  &
        -k(733)*n(idx_HCO)  &
        -k(732)*n(idx_HCO)

    !d[NH3_dot]/d[OH+]
    pd(16,71) =  &
        -k(200)*n(idx_NH3)

    !d[NO_dot]/d[OH+]
    pd(17,71) =  &
        -k(736)*n(idx_NO)  &
        -k(201)*n(idx_NO)

    !d[CN_dot]/d[OH+]
    pd(18,71) =  &
        -k(726)*n(idx_CN)

    !d[CO_dot]/d[OH+]
    pd(19,71) =  &
        -k(728)*n(idx_CO)  &
        +k(732)*n(idx_HCO)

    !d[N2_dot]/d[OH+]
    pd(20,71) =  &
        -k(735)*n(idx_N2)

    !d[NH2_dot]/d[OH+]
    pd(21,71) =  &
        -k(686)*n(idx_NH2)  &
        -k(168)*n(idx_NH2)

    !d[CH4_dot]/d[OH+]
    pd(23,71) =  &
        -k(399)*n(idx_CH4)

    !d[N_dot]/d[OH+]
    pd(24,71) =  &
        -k(641)*n(idx_N)

    !d[NH_dot]/d[OH+]
    pd(25,71) =  &
        -k(701)*n(idx_NH)

    !d[CO2_dot]/d[OH+]
    pd(29,71) =  &
        -k(727)*n(idx_CO2)

    !d[H2O_DUST_dot]/d[OH+]
    pd(40,71) =  &
        +k(1094)

    !d[HCO+_dot]/d[OH+]
    pd(54,71) =  &
        +k(728)*n(idx_CO)  &
        +k(199)*n(idx_HCO)

    !d[CH2+_dot]/d[OH+]
    pd(58,71) =  &
        +k(417)*n(idx_CH)  &
        +k(39)*n(idx_CH2)

    !d[CH+_dot]/d[OH+]
    pd(59,71) =  &
        +k(56)*n(idx_CH)  &
        +k(339)*n(idx_C)

    !d[H2CO+_dot]/d[OH+]
    pd(60,71) =  &
        +k(733)*n(idx_HCO)  &
        +k(197)*n(idx_H2CO)

    !d[NH3+_dot]/d[OH+]
    pd(62,71) =  &
        +k(686)*n(idx_NH2)  &
        +k(200)*n(idx_NH3)

    !d[NO+_dot]/d[OH+]
    pd(63,71) =  &
        +k(641)*n(idx_N)  &
        +k(201)*n(idx_NO)

    !d[O2+_dot]/d[OH+]
    pd(67,71) =  &
        +k(725)*n(idx_O)  &
        +k(202)*n(idx_O2)

    !d[H2O+_dot]/d[OH+]
    pd(68,71) =  &
        +k(732)*n(idx_HCO)  &
        +k(198)*n(idx_H2O)  &
        +k(737)*n(idx_OH)  &
        +k(477)*n(idx_H2)

    !d[NH2+_dot]/d[OH+]
    pd(69,71) =  &
        +k(168)*n(idx_NH2)  &
        +k(701)*n(idx_NH)

    !d[O+_dot]/d[OH+]
    pd(70,71) =  &
        +k(1037)

    !d[OH+_dot]/d[OH+]
    pd(71,71) =  &
        -k(728)*n(idx_CO)  &
        -k(726)*n(idx_CN)  &
        -k(477)*n(idx_H2)  &
        -k(725)*n(idx_O)  &
        -k(727)*n(idx_CO2)  &
        -k(736)*n(idx_NO)  &
        -k(381)*n(idx_CH2)  &
        -k(199)*n(idx_HCO)  &
        -k(313)*n(idx_E)  &
        -k(168)*n(idx_NH2)  &
        -k(731)*n(idx_HCN)  &
        -k(198)*n(idx_H2O)  &
        -k(735)*n(idx_N2)  &
        -k(200)*n(idx_NH3)  &
        -k(730)*n(idx_H2O)  &
        -k(399)*n(idx_CH4)  &
        -k(1037)  &
        -k(729)*n(idx_H2CO)  &
        -k(339)*n(idx_C)  &
        -k(734)*n(idx_HNC)  &
        -k(1094)  &
        -k(202)*n(idx_O2)  &
        -k(201)*n(idx_NO)  &
        -k(732)*n(idx_HCO)  &
        -k(39)*n(idx_CH2)  &
        -k(737)*n(idx_OH)  &
        -k(641)*n(idx_N)  &
        -k(56)*n(idx_CH)  &
        -k(701)*n(idx_NH)  &
        -k(417)*n(idx_CH)  &
        -k(733)*n(idx_HCO)  &
        -k(686)*n(idx_NH2)  &
        -k(197)*n(idx_H2CO)

    !d[CH3+_dot]/d[OH+]
    pd(72,71) =  &
        +k(381)*n(idx_CH2)

    !d[HCN+_dot]/d[OH+]
    pd(75,71) =  &
        +k(726)*n(idx_CN)

    !d[HNO+_dot]/d[OH+]
    pd(79,71) =  &
        +k(736)*n(idx_NO)

    !d[H3CO+_dot]/d[OH+]
    pd(82,71) =  &
        +k(729)*n(idx_H2CO)

    !d[H3O+_dot]/d[OH+]
    pd(83,71) =  &
        +k(399)*n(idx_CH4)  &
        +k(730)*n(idx_H2O)

    !d[HCNH+_dot]/d[OH+]
    pd(84,71) =  &
        +k(731)*n(idx_HCN)  &
        +k(734)*n(idx_HNC)

    !d[HCO2+_dot]/d[OH+]
    pd(85,71) =  &
        +k(727)*n(idx_CO2)

    !d[N2H+_dot]/d[OH+]
    pd(87,71) =  &
        +k(735)*n(idx_N2)

    !d[E_dot]/d[CH3+]
    pd(1,72) =  &
        -k(265)*n(idx_E)  &
        -k(264)*n(idx_E)  &
        -k(263)*n(idx_E)  &
        -k(1057)*n(idx_E)

    !d[CH_dot]/d[CH3+]
    pd(2,72) =  &
        +k(265)*n(idx_E)  &
        +k(264)*n(idx_E)

    !d[O_dot]/d[CH3+]
    pd(3,72) =  &
        +k(385)*n(idx_O2)  &
        -k(386)*n(idx_O)  &
        -k(387)*n(idx_O)

    !d[H2_dot]/d[CH3+]
    pd(6,72) =  &
        +k(532)*n(idx_H)  &
        +k(264)*n(idx_E)  &
        +k(387)*n(idx_O)  &
        +k(980)  &
        +k(689)*n(idx_NH)  &
        +k(388)*n(idx_OH)

    !d[H_dot]/d[CH3+]
    pd(8,72) =  &
        +k(263)*n(idx_E)  &
        +k(386)*n(idx_O)  &
        +2.d0*k(265)*n(idx_E)  &
        +k(981)  &
        -k(532)*n(idx_H)

    !d[OH_dot]/d[CH3+]
    pd(10,72) =  &
        -k(388)*n(idx_OH)

    !d[O2_dot]/d[CH3+]
    pd(11,72) =  &
        -k(385)*n(idx_O2)

    !d[CH2_dot]/d[CH3+]
    pd(12,72) =  &
        +k(263)*n(idx_E)

    !d[H2CO_dot]/d[CH3+]
    pd(13,72) =  &
        -k(383)*n(idx_H2CO)

    !d[HCO_dot]/d[CH3+]
    pd(14,72) =  &
        -k(40)*n(idx_HCO)  &
        -k(384)*n(idx_HCO)

    !d[MG_dot]/d[CH3+]
    pd(15,72) =  &
        -k(41)*n(idx_MG)

    !d[NO_dot]/d[CH3+]
    pd(17,72) =  &
        -k(42)*n(idx_NO)

    !d[CO_dot]/d[CH3+]
    pd(19,72) =  &
        +k(384)*n(idx_HCO)

    !d[CH3_dot]/d[CH3+]
    pd(22,72) =  &
        +k(41)*n(idx_MG)  &
        +k(42)*n(idx_NO)  &
        +k(1057)*n(idx_E)  &
        +k(40)*n(idx_HCO)

    !d[CH4_dot]/d[CH3+]
    pd(23,72) =  &
        +k(382)*n(idx_CH3OH)  &
        +k(383)*n(idx_H2CO)

    !d[NH_dot]/d[CH3+]
    pd(25,72) =  &
        -k(689)*n(idx_NH)

    !d[CH3OH_dot]/d[CH3+]
    pd(28,72) =  &
        -k(382)*n(idx_CH3OH)

    !d[CH4_DUST_dot]/d[CH3+]
    pd(38,72) =  &
        +k(1105)

    !d[HCO+_dot]/d[CH3+]
    pd(54,72) =  &
        +k(387)*n(idx_O)  &
        +k(383)*n(idx_H2CO)  &
        +k(40)*n(idx_HCO)

    !d[CH2+_dot]/d[CH3+]
    pd(58,72) =  &
        +k(981)  &
        +k(532)*n(idx_H)

    !d[CH+_dot]/d[CH3+]
    pd(59,72) =  &
        +k(980)

    !d[H2CO+_dot]/d[CH3+]
    pd(60,72) =  &
        +k(386)*n(idx_O)  &
        +k(388)*n(idx_OH)

    !d[MG+_dot]/d[CH3+]
    pd(61,72) =  &
        +k(41)*n(idx_MG)

    !d[NO+_dot]/d[CH3+]
    pd(63,72) =  &
        +k(42)*n(idx_NO)

    !d[CH3+_dot]/d[CH3+]
    pd(72,72) =  &
        -k(1057)*n(idx_E)  &
        -k(689)*n(idx_NH)  &
        -k(265)*n(idx_E)  &
        -k(264)*n(idx_E)  &
        -k(263)*n(idx_E)  &
        -k(532)*n(idx_H)  &
        -k(387)*n(idx_O)  &
        -k(385)*n(idx_O2)  &
        -k(980)  &
        -k(1105)  &
        -k(386)*n(idx_O)  &
        -k(388)*n(idx_OH)  &
        -k(384)*n(idx_HCO)  &
        -k(40)*n(idx_HCO)  &
        -k(981)  &
        -k(41)*n(idx_MG)  &
        -k(382)*n(idx_CH3OH)  &
        -k(42)*n(idx_NO)  &
        -k(383)*n(idx_H2CO)

    !d[CH4+_dot]/d[CH3+]
    pd(73,72) =  &
        +k(384)*n(idx_HCO)

    !d[H3CO+_dot]/d[CH3+]
    pd(82,72) =  &
        +k(385)*n(idx_O2)  &
        +k(382)*n(idx_CH3OH)

    !d[HCNH+_dot]/d[CH3+]
    pd(84,72) =  &
        +k(689)*n(idx_NH)

    !d[E_dot]/d[CH4+]
    pd(1,73) =  &
        -k(267)*n(idx_E)  &
        -k(266)*n(idx_E)

    !d[O_dot]/d[CH4+]
    pd(3,73) =  &
        -k(717)*n(idx_O)

    !d[H2_dot]/d[CH4+]
    pd(6,73) =  &
        +k(533)*n(idx_H)  &
        +k(988)

    !d[H_dot]/d[CH4+]
    pd(8,73) =  &
        +k(989)  &
        -k(533)*n(idx_H)  &
        +k(267)*n(idx_E)  &
        +2.d0*k(266)*n(idx_E)

    !d[H2O_dot]/d[CH4+]
    pd(9,73) =  &
        -k(392)*n(idx_H2O)

    !d[OH_dot]/d[CH4+]
    pd(10,73) =  &
        +k(717)*n(idx_O)

    !d[O2_dot]/d[CH4+]
    pd(11,73) =  &
        -k(45)*n(idx_O2)

    !d[CH2_dot]/d[CH4+]
    pd(12,73) =  &
        +k(266)*n(idx_E)

    !d[H2CO_dot]/d[CH4+]
    pd(13,73) =  &
        -k(43)*n(idx_H2CO)  &
        -k(391)*n(idx_H2CO)

    !d[NH3_dot]/d[CH4+]
    pd(16,73) =  &
        -k(44)*n(idx_NH3)

    !d[CO_dot]/d[CH4+]
    pd(19,73) =  &
        -k(390)*n(idx_CO)

    !d[CH3_dot]/d[CH4+]
    pd(22,73) =  &
        +k(389)*n(idx_CO2)  &
        +k(267)*n(idx_E)  &
        +k(390)*n(idx_CO)  &
        +k(391)*n(idx_H2CO)  &
        +k(392)*n(idx_H2O)

    !d[CH4_dot]/d[CH4+]
    pd(23,73) =  &
        +k(45)*n(idx_O2)  &
        +k(43)*n(idx_H2CO)  &
        +k(44)*n(idx_NH3)

    !d[CO2_dot]/d[CH4+]
    pd(29,73) =  &
        -k(389)*n(idx_CO2)

    !d[CH4_DUST_dot]/d[CH4+]
    pd(38,73) =  &
        +k(1108)

    !d[HCO+_dot]/d[CH4+]
    pd(54,73) =  &
        +k(390)*n(idx_CO)

    !d[CH2+_dot]/d[CH4+]
    pd(58,73) =  &
        +k(988)

    !d[H2CO+_dot]/d[CH4+]
    pd(60,73) =  &
        +k(43)*n(idx_H2CO)

    !d[NH3+_dot]/d[CH4+]
    pd(62,73) =  &
        +k(44)*n(idx_NH3)

    !d[O2+_dot]/d[CH4+]
    pd(67,73) =  &
        +k(45)*n(idx_O2)

    !d[CH3+_dot]/d[CH4+]
    pd(72,73) =  &
        +k(989)  &
        +k(533)*n(idx_H)  &
        +k(717)*n(idx_O)

    !d[CH4+_dot]/d[CH4+]
    pd(73,73) =  &
        -k(1108)  &
        -k(267)*n(idx_E)  &
        -k(717)*n(idx_O)  &
        -k(391)*n(idx_H2CO)  &
        -k(43)*n(idx_H2CO)  &
        -k(266)*n(idx_E)  &
        -k(988)  &
        -k(989)  &
        -k(390)*n(idx_CO)  &
        -k(45)*n(idx_O2)  &
        -k(533)*n(idx_H)  &
        -k(392)*n(idx_H2O)  &
        -k(44)*n(idx_NH3)  &
        -k(389)*n(idx_CO2)

    !d[H3CO+_dot]/d[CH4+]
    pd(82,73) =  &
        +k(391)*n(idx_H2CO)

    !d[H3O+_dot]/d[CH4+]
    pd(83,73) =  &
        +k(392)*n(idx_H2O)

    !d[HCO2+_dot]/d[CH4+]
    pd(85,73) =  &
        +k(389)*n(idx_CO2)

    !d[E_dot]/d[N+]
    pd(1,74) =  &
        -k(1062)*n(idx_E)

    !d[CH_dot]/d[N+]
    pd(2,74) =  &
        -k(410)*n(idx_CH)  &
        -k(51)*n(idx_CH)

    !d[O_dot]/d[N+]
    pd(3,74) =  &
        +k(626)*n(idx_O2)  &
        +k(625)*n(idx_NO)

    !d[HCN_dot]/d[N+]
    pd(5,74) =  &
        -k(141)*n(idx_HCN)

    !d[H2_dot]/d[N+]
    pd(6,74) =  &
        +k(622)*n(idx_NH3)  &
        +k(615)*n(idx_CH4)  &
        -k(470)*n(idx_H2)

    !d[C_dot]/d[N+]
    pd(7,74) =  &
        +k(618)*n(idx_CO)

    !d[H_dot]/d[N+]
    pd(8,74) =  &
        +k(410)*n(idx_CH)  &
        +k(610)*n(idx_CH3OH)  &
        +k(470)*n(idx_H2)  &
        +k(612)*n(idx_CH3OH)  &
        +k(613)*n(idx_CH3OH)  &
        +k(614)*n(idx_CH4)  &
        +2.d0*k(616)*n(idx_CH4)  &
        +k(615)*n(idx_CH4)  &
        +k(624)*n(idx_NH)

    !d[H2O_dot]/d[N+]
    pd(9,74) =  &
        -k(140)*n(idx_H2O)

    !d[OH_dot]/d[N+]
    pd(10,74) =  &
        -k(149)*n(idx_OH)

    !d[O2_dot]/d[N+]
    pd(11,74) =  &
        -k(148)*n(idx_O2)  &
        -k(627)*n(idx_O2)  &
        -k(626)*n(idx_O2)

    !d[CH2_dot]/d[N+]
    pd(12,74) =  &
        +k(620)*n(idx_H2CO)  &
        -k(135)*n(idx_CH2)

    !d[H2CO_dot]/d[N+]
    pd(13,74) =  &
        -k(620)*n(idx_H2CO)  &
        -k(139)*n(idx_H2CO)  &
        -k(619)*n(idx_H2CO)

    !d[HCO_dot]/d[N+]
    pd(14,74) =  &
        -k(142)*n(idx_HCO)  &
        -k(621)*n(idx_HCO)

    !d[MG_dot]/d[N+]
    pd(15,74) =  &
        -k(143)*n(idx_MG)

    !d[NH3_dot]/d[N+]
    pd(16,74) =  &
        -k(623)*n(idx_NH3)  &
        -k(622)*n(idx_NH3)  &
        -k(145)*n(idx_NH3)

    !d[NO_dot]/d[N+]
    pd(17,74) =  &
        -k(147)*n(idx_NO)  &
        +k(613)*n(idx_CH3OH)  &
        +k(617)*n(idx_CO2)  &
        -k(625)*n(idx_NO)  &
        +k(627)*n(idx_O2)

    !d[CN_dot]/d[N+]
    pd(18,74) =  &
        -k(137)*n(idx_CN)

    !d[CO_dot]/d[N+]
    pd(19,74) =  &
        +k(621)*n(idx_HCO)  &
        -k(138)*n(idx_CO)  &
        -k(618)*n(idx_CO)

    !d[NH2_dot]/d[N+]
    pd(21,74) =  &
        -k(144)*n(idx_NH2)

    !d[CH3_dot]/d[N+]
    pd(22,74) =  &
        +k(612)*n(idx_CH3OH)

    !d[CH4_dot]/d[N+]
    pd(23,74) =  &
        -k(614)*n(idx_CH4)  &
        -k(136)*n(idx_CH4)  &
        -k(615)*n(idx_CH4)  &
        -k(616)*n(idx_CH4)

    !d[N_dot]/d[N+]
    pd(24,74) =  &
        +k(148)*n(idx_O2)  &
        +k(141)*n(idx_HCN)  &
        +k(149)*n(idx_OH)  &
        +k(136)*n(idx_CH4)  &
        +k(614)*n(idx_CH4)  &
        +k(146)*n(idx_NH)  &
        -k(1054)*n(idx_N)  &
        +k(137)*n(idx_CN)  &
        +k(143)*n(idx_MG)  &
        +k(51)*n(idx_CH)  &
        +k(140)*n(idx_H2O)  &
        +k(138)*n(idx_CO)  &
        +k(139)*n(idx_H2CO)  &
        +k(135)*n(idx_CH2)  &
        +k(1062)*n(idx_E)  &
        +k(147)*n(idx_NO)  &
        +k(142)*n(idx_HCO)  &
        +k(145)*n(idx_NH3)  &
        +k(144)*n(idx_NH2)

    !d[NH_dot]/d[N+]
    pd(25,74) =  &
        +k(619)*n(idx_H2CO)  &
        -k(146)*n(idx_NH)  &
        -k(624)*n(idx_NH)  &
        +k(610)*n(idx_CH3OH)  &
        +k(611)*n(idx_CH3OH)  &
        +k(623)*n(idx_NH3)

    !d[CH3OH_dot]/d[N+]
    pd(28,74) =  &
        -k(611)*n(idx_CH3OH)  &
        -k(613)*n(idx_CH3OH)  &
        -k(610)*n(idx_CH3OH)  &
        -k(612)*n(idx_CH3OH)

    !d[CO2_dot]/d[N+]
    pd(29,74) =  &
        -k(617)*n(idx_CO2)

    !d[NH3_DUST_dot]/d[N+]
    pd(45,74) =  &
        +k(1088)

    !d[HCO+_dot]/d[N+]
    pd(54,74) =  &
        +k(619)*n(idx_H2CO)  &
        +k(142)*n(idx_HCO)

    !d[CH2+_dot]/d[N+]
    pd(58,74) =  &
        +k(135)*n(idx_CH2)

    !d[CH+_dot]/d[N+]
    pd(59,74) =  &
        +k(51)*n(idx_CH)

    !d[H2CO+_dot]/d[N+]
    pd(60,74) =  &
        +k(610)*n(idx_CH3OH)  &
        +k(139)*n(idx_H2CO)

    !d[MG+_dot]/d[N+]
    pd(61,74) =  &
        +k(143)*n(idx_MG)

    !d[NH3+_dot]/d[N+]
    pd(62,74) =  &
        +k(145)*n(idx_NH3)

    !d[NO+_dot]/d[N+]
    pd(63,74) =  &
        +k(626)*n(idx_O2)  &
        +k(147)*n(idx_NO)  &
        +k(612)*n(idx_CH3OH)  &
        +k(620)*n(idx_H2CO)  &
        +k(618)*n(idx_CO)

    !d[CN+_dot]/d[N+]
    pd(64,74) =  &
        +k(137)*n(idx_CN)  &
        +k(410)*n(idx_CH)

    !d[CO+_dot]/d[N+]
    pd(65,74) =  &
        +k(617)*n(idx_CO2)  &
        +k(138)*n(idx_CO)

    !d[N2+_dot]/d[N+]
    pd(66,74) =  &
        +k(625)*n(idx_NO)  &
        +k(624)*n(idx_NH)  &
        +k(1054)*n(idx_N)

    !d[O2+_dot]/d[N+]
    pd(67,74) =  &
        +k(148)*n(idx_O2)

    !d[H2O+_dot]/d[N+]
    pd(68,74) =  &
        +k(140)*n(idx_H2O)

    !d[NH2+_dot]/d[N+]
    pd(69,74) =  &
        +k(623)*n(idx_NH3)  &
        +k(144)*n(idx_NH2)

    !d[O+_dot]/d[N+]
    pd(70,74) =  &
        +k(627)*n(idx_O2)

    !d[OH+_dot]/d[N+]
    pd(71,74) =  &
        +k(149)*n(idx_OH)

    !d[CH3+_dot]/d[N+]
    pd(72,74) =  &
        +k(613)*n(idx_CH3OH)  &
        +k(614)*n(idx_CH4)

    !d[CH4+_dot]/d[N+]
    pd(73,74) =  &
        +k(136)*n(idx_CH4)

    !d[N+_dot]/d[N+]
    pd(74,74) =  &
        -k(146)*n(idx_NH)  &
        -k(148)*n(idx_O2)  &
        -k(611)*n(idx_CH3OH)  &
        -k(614)*n(idx_CH4)  &
        -k(1054)*n(idx_N)  &
        -k(470)*n(idx_H2)  &
        -k(622)*n(idx_NH3)  &
        -k(142)*n(idx_HCO)  &
        -k(136)*n(idx_CH4)  &
        -k(610)*n(idx_CH3OH)  &
        -k(617)*n(idx_CO2)  &
        -k(138)*n(idx_CO)  &
        -k(410)*n(idx_CH)  &
        -k(141)*n(idx_HCN)  &
        -k(612)*n(idx_CH3OH)  &
        -k(149)*n(idx_OH)  &
        -k(616)*n(idx_CH4)  &
        -k(625)*n(idx_NO)  &
        -k(144)*n(idx_NH2)  &
        -k(621)*n(idx_HCO)  &
        -k(137)*n(idx_CN)  &
        -k(626)*n(idx_O2)  &
        -k(623)*n(idx_NH3)  &
        -k(140)*n(idx_H2O)  &
        -k(1088)  &
        -k(619)*n(idx_H2CO)  &
        -k(613)*n(idx_CH3OH)  &
        -k(51)*n(idx_CH)  &
        -k(624)*n(idx_NH)  &
        -k(139)*n(idx_H2CO)  &
        -k(627)*n(idx_O2)  &
        -k(620)*n(idx_H2CO)  &
        -k(143)*n(idx_MG)  &
        -k(615)*n(idx_CH4)  &
        -k(135)*n(idx_CH2)  &
        -k(145)*n(idx_NH3)  &
        -k(147)*n(idx_NO)  &
        -k(1062)*n(idx_E)  &
        -k(618)*n(idx_CO)

    !d[HCN+_dot]/d[N+]
    pd(75,74) =  &
        +k(141)*n(idx_HCN)  &
        +k(615)*n(idx_CH4)

    !d[NH+_dot]/d[N+]
    pd(76,74) =  &
        +k(146)*n(idx_NH)  &
        +k(621)*n(idx_HCO)  &
        +k(470)*n(idx_H2)

    !d[H3CO+_dot]/d[N+]
    pd(82,74) =  &
        +k(611)*n(idx_CH3OH)

    !d[HCNH+_dot]/d[N+]
    pd(84,74) =  &
        +k(616)*n(idx_CH4)

    !d[N2H+_dot]/d[N+]
    pd(87,74) =  &
        +k(622)*n(idx_NH3)

    !d[E_dot]/d[HCN+]
    pd(1,75) =  &
        -k(291)*n(idx_E)

    !d[CH_dot]/d[HCN+]
    pd(2,75) =  &
        -k(405)*n(idx_CH)

    !d[HNC_dot]/d[HCN+]
    pd(4,75) =  &
        -k(541)*n(idx_HNC)

    !d[HCN_dot]/d[HCN+]
    pd(5,75) =  &
        +k(117)*n(idx_O2)  &
        +k(175)*n(idx_NH3)  &
        +k(108)*n(idx_H2O)  &
        +k(113)*n(idx_H)  &
        +k(116)*n(idx_NO)  &
        -k(538)*n(idx_HCN)

    !d[H2_dot]/d[HCN+]
    pd(6,75) =  &
        -k(467)*n(idx_H2)

    !d[C_dot]/d[HCN+]
    pd(7,75) =  &
        -k(331)*n(idx_C)

    !d[H_dot]/d[HCN+]
    pd(8,75) =  &
        -k(113)*n(idx_H)  &
        +k(467)*n(idx_H2)  &
        +k(291)*n(idx_E)

    !d[H2O_dot]/d[HCN+]
    pd(9,75) =  &
        -k(495)*n(idx_H2O)  &
        -k(108)*n(idx_H2O)

    !d[OH_dot]/d[HCN+]
    pd(10,75) =  &
        -k(740)*n(idx_OH)

    !d[O2_dot]/d[HCN+]
    pd(11,75) =  &
        -k(117)*n(idx_O2)

    !d[CH2_dot]/d[HCN+]
    pd(12,75) =  &
        -k(370)*n(idx_CH2)

    !d[H2CO_dot]/d[HCN+]
    pd(13,75) =  &
        -k(537)*n(idx_H2CO)

    !d[HCO_dot]/d[HCN+]
    pd(14,75) =  &
        -k(540)*n(idx_HCO)  &
        -k(539)*n(idx_HCO)

    !d[NH3_dot]/d[HCN+]
    pd(16,75) =  &
        -k(688)*n(idx_NH3)  &
        -k(175)*n(idx_NH3)

    !d[NO_dot]/d[HCN+]
    pd(17,75) =  &
        -k(116)*n(idx_NO)

    !d[CN_dot]/d[HCN+]
    pd(18,75) =  &
        +k(538)*n(idx_HCN)  &
        +k(539)*n(idx_HCO)  &
        +k(370)*n(idx_CH2)  &
        +k(740)*n(idx_OH)  &
        +k(495)*n(idx_H2O)  &
        +k(537)*n(idx_H2CO)  &
        +k(541)*n(idx_HNC)  &
        +k(291)*n(idx_E)  &
        +k(536)*n(idx_CO)  &
        +k(405)*n(idx_CH)  &
        +k(331)*n(idx_C)  &
        +k(535)*n(idx_CO2)  &
        +k(693)*n(idx_NH)  &
        +k(679)*n(idx_NH2)

    !d[CO_dot]/d[HCN+]
    pd(19,75) =  &
        +k(540)*n(idx_HCO)  &
        -k(536)*n(idx_CO)

    !d[NH2_dot]/d[HCN+]
    pd(21,75) =  &
        -k(679)*n(idx_NH2)  &
        +k(688)*n(idx_NH3)

    !d[CH3_dot]/d[HCN+]
    pd(22,75) =  &
        +k(396)*n(idx_CH4)

    !d[CH4_dot]/d[HCN+]
    pd(23,75) =  &
        -k(396)*n(idx_CH4)

    !d[NH_dot]/d[HCN+]
    pd(25,75) =  &
        -k(693)*n(idx_NH)

    !d[CO2_dot]/d[HCN+]
    pd(29,75) =  &
        -k(535)*n(idx_CO2)

    !d[HCN_DUST_dot]/d[HCN+]
    pd(44,75) =  &
        +k(1102)

    !d[HCO+_dot]/d[HCN+]
    pd(54,75) =  &
        +k(536)*n(idx_CO)

    !d[H+_dot]/d[HCN+]
    pd(55,75) =  &
        +k(113)*n(idx_H)

    !d[CH2+_dot]/d[HCN+]
    pd(58,75) =  &
        +k(405)*n(idx_CH)

    !d[CH+_dot]/d[HCN+]
    pd(59,75) =  &
        +k(331)*n(idx_C)

    !d[H2CO+_dot]/d[HCN+]
    pd(60,75) =  &
        +k(539)*n(idx_HCO)

    !d[NH3+_dot]/d[HCN+]
    pd(62,75) =  &
        +k(175)*n(idx_NH3)  &
        +k(679)*n(idx_NH2)

    !d[NO+_dot]/d[HCN+]
    pd(63,75) =  &
        +k(116)*n(idx_NO)

    !d[O2+_dot]/d[HCN+]
    pd(67,75) =  &
        +k(117)*n(idx_O2)

    !d[H2O+_dot]/d[HCN+]
    pd(68,75) =  &
        +k(108)*n(idx_H2O)  &
        +k(740)*n(idx_OH)

    !d[NH2+_dot]/d[HCN+]
    pd(69,75) =  &
        +k(693)*n(idx_NH)

    !d[CH3+_dot]/d[HCN+]
    pd(72,75) =  &
        +k(370)*n(idx_CH2)

    !d[HCN+_dot]/d[HCN+]
    pd(75,75) =  &
        -k(108)*n(idx_H2O)  &
        -k(740)*n(idx_OH)  &
        -k(331)*n(idx_C)  &
        -k(396)*n(idx_CH4)  &
        -k(116)*n(idx_NO)  &
        -k(693)*n(idx_NH)  &
        -k(405)*n(idx_CH)  &
        -k(113)*n(idx_H)  &
        -k(535)*n(idx_CO2)  &
        -k(1102)  &
        -k(117)*n(idx_O2)  &
        -k(688)*n(idx_NH3)  &
        -k(495)*n(idx_H2O)  &
        -k(370)*n(idx_CH2)  &
        -k(537)*n(idx_H2CO)  &
        -k(467)*n(idx_H2)  &
        -k(679)*n(idx_NH2)  &
        -k(291)*n(idx_E)  &
        -k(536)*n(idx_CO)  &
        -k(540)*n(idx_HCO)  &
        -k(539)*n(idx_HCO)  &
        -k(538)*n(idx_HCN)  &
        -k(175)*n(idx_NH3)  &
        -k(541)*n(idx_HNC)

    !d[H3CO+_dot]/d[HCN+]
    pd(82,75) =  &
        +k(537)*n(idx_H2CO)

    !d[H3O+_dot]/d[HCN+]
    pd(83,75) =  &
        +k(495)*n(idx_H2O)

    !d[HCNH+_dot]/d[HCN+]
    pd(84,75) =  &
        +k(396)*n(idx_CH4)  &
        +k(540)*n(idx_HCO)  &
        +k(541)*n(idx_HNC)  &
        +k(688)*n(idx_NH3)  &
        +k(467)*n(idx_H2)  &
        +k(538)*n(idx_HCN)

    !d[HCO2+_dot]/d[HCN+]
    pd(85,75) =  &
        +k(535)*n(idx_CO2)

    !d[E_dot]/d[NH+]
    pd(1,76) =  &
        -k(305)*n(idx_E)

    !d[CH_dot]/d[NH+]
    pd(2,76) =  &
        -k(412)*n(idx_CH)

    !d[O_dot]/d[NH+]
    pd(3,76) =  &
        -k(662)*n(idx_O)  &
        +k(659)*n(idx_NO)  &
        +k(651)*n(idx_H2O)

    !d[HNC_dot]/d[NH+]
    pd(4,76) =  &
        -k(655)*n(idx_HNC)

    !d[HCN_dot]/d[NH+]
    pd(5,76) =  &
        -k(653)*n(idx_HCN)

    !d[H2_dot]/d[NH+]
    pd(6,76) =  &
        +k(650)*n(idx_H2O)  &
        -k(473)*n(idx_H2)  &
        -k(472)*n(idx_H2)

    !d[C_dot]/d[NH+]
    pd(7,76) =  &
        -k(336)*n(idx_C)

    !d[H_dot]/d[NH+]
    pd(8,76) =  &
        +k(473)*n(idx_H2)  &
        +k(305)*n(idx_E)  &
        +k(638)*n(idx_N)

    !d[H2O_dot]/d[NH+]
    pd(9,76) =  &
        -k(652)*n(idx_H2O)  &
        -k(651)*n(idx_H2O)  &
        -k(156)*n(idx_H2O)  &
        -k(649)*n(idx_H2O)  &
        -k(650)*n(idx_H2O)

    !d[OH_dot]/d[NH+]
    pd(10,76) =  &
        +k(660)*n(idx_O2)  &
        +k(652)*n(idx_H2O)  &
        -k(663)*n(idx_OH)

    !d[O2_dot]/d[NH+]
    pd(11,76) =  &
        -k(159)*n(idx_O2)  &
        -k(661)*n(idx_O2)  &
        -k(660)*n(idx_O2)

    !d[CH2_dot]/d[NH+]
    pd(12,76) =  &
        -k(376)*n(idx_CH2)

    !d[H2CO_dot]/d[NH+]
    pd(13,76) =  &
        -k(647)*n(idx_H2CO)  &
        -k(648)*n(idx_H2CO)  &
        -k(155)*n(idx_H2CO)

    !d[HCO_dot]/d[NH+]
    pd(14,76) =  &
        +k(645)*n(idx_CO2)  &
        -k(654)*n(idx_HCO)

    !d[NH3_dot]/d[NH+]
    pd(16,76) =  &
        -k(157)*n(idx_NH3)

    !d[NO_dot]/d[NH+]
    pd(17,76) =  &
        -k(158)*n(idx_NO)  &
        -k(659)*n(idx_NO)

    !d[CN_dot]/d[NH+]
    pd(18,76) =  &
        -k(642)*n(idx_CN)

    !d[CO_dot]/d[NH+]
    pd(19,76) =  &
        +k(644)*n(idx_CO2)  &
        -k(646)*n(idx_CO)

    !d[N2_dot]/d[NH+]
    pd(20,76) =  &
        -k(656)*n(idx_N2)

    !d[NH2_dot]/d[NH+]
    pd(21,76) =  &
        +k(648)*n(idx_H2CO)  &
        -k(657)*n(idx_NH2)

    !d[N_dot]/d[NH+]
    pd(24,76) =  &
        +k(1020)  &
        +k(376)*n(idx_CH2)  &
        +k(643)*n(idx_CO2)  &
        +k(657)*n(idx_NH2)  &
        +k(646)*n(idx_CO)  &
        +k(663)*n(idx_OH)  &
        +k(336)*n(idx_C)  &
        +k(412)*n(idx_CH)  &
        +k(653)*n(idx_HCN)  &
        +k(642)*n(idx_CN)  &
        -k(638)*n(idx_N)  &
        +k(661)*n(idx_O2)  &
        +k(647)*n(idx_H2CO)  &
        +k(658)*n(idx_NH)  &
        +k(655)*n(idx_HNC)  &
        +k(656)*n(idx_N2)  &
        +k(662)*n(idx_O)  &
        +k(472)*n(idx_H2)  &
        +k(305)*n(idx_E)  &
        +k(649)*n(idx_H2O)  &
        +k(654)*n(idx_HCO)

    !d[NH_dot]/d[NH+]
    pd(25,76) =  &
        +k(156)*n(idx_H2O)  &
        -k(658)*n(idx_NH)  &
        +k(155)*n(idx_H2CO)  &
        +k(158)*n(idx_NO)  &
        +k(159)*n(idx_O2)  &
        +k(157)*n(idx_NH3)

    !d[CO2_dot]/d[NH+]
    pd(29,76) =  &
        -k(643)*n(idx_CO2)  &
        -k(644)*n(idx_CO2)  &
        -k(645)*n(idx_CO2)

    !d[NH3_DUST_dot]/d[NH+]
    pd(45,76) =  &
        +k(1093)

    !d[HCO+_dot]/d[NH+]
    pd(54,76) =  &
        +k(648)*n(idx_H2CO)  &
        +k(646)*n(idx_CO)

    !d[H+_dot]/d[NH+]
    pd(55,76) =  &
        +k(1020)

    !d[CH2+_dot]/d[NH+]
    pd(58,76) =  &
        +k(412)*n(idx_CH)

    !d[CH+_dot]/d[NH+]
    pd(59,76) =  &
        +k(336)*n(idx_C)

    !d[H2CO+_dot]/d[NH+]
    pd(60,76) =  &
        +k(155)*n(idx_H2CO)  &
        +k(654)*n(idx_HCO)

    !d[NH3+_dot]/d[NH+]
    pd(62,76) =  &
        +k(657)*n(idx_NH2)  &
        +k(157)*n(idx_NH3)  &
        +k(651)*n(idx_H2O)

    !d[NO+_dot]/d[NH+]
    pd(63,76) =  &
        +k(645)*n(idx_CO2)  &
        +k(660)*n(idx_O2)  &
        +k(158)*n(idx_NO)

    !d[N2+_dot]/d[NH+]
    pd(66,76) =  &
        +k(638)*n(idx_N)

    !d[O2+_dot]/d[NH+]
    pd(67,76) =  &
        +k(159)*n(idx_O2)

    !d[H2O+_dot]/d[NH+]
    pd(68,76) =  &
        +k(156)*n(idx_H2O)  &
        +k(663)*n(idx_OH)

    !d[NH2+_dot]/d[NH+]
    pd(69,76) =  &
        +k(473)*n(idx_H2)  &
        +k(658)*n(idx_NH)  &
        +k(652)*n(idx_H2O)

    !d[OH+_dot]/d[NH+]
    pd(71,76) =  &
        +k(662)*n(idx_O)

    !d[CH3+_dot]/d[NH+]
    pd(72,76) =  &
        +k(376)*n(idx_CH2)

    !d[HCN+_dot]/d[NH+]
    pd(75,76) =  &
        +k(642)*n(idx_CN)

    !d[NH+_dot]/d[NH+]
    pd(76,76) =  &
        -k(473)*n(idx_H2)  &
        -k(654)*n(idx_HCO)  &
        -k(1020)  &
        -k(157)*n(idx_NH3)  &
        -k(651)*n(idx_H2O)  &
        -k(642)*n(idx_CN)  &
        -k(305)*n(idx_E)  &
        -k(661)*n(idx_O2)  &
        -k(656)*n(idx_N2)  &
        -k(657)*n(idx_NH2)  &
        -k(663)*n(idx_OH)  &
        -k(376)*n(idx_CH2)  &
        -k(652)*n(idx_H2O)  &
        -k(638)*n(idx_N)  &
        -k(644)*n(idx_CO2)  &
        -k(660)*n(idx_O2)  &
        -k(653)*n(idx_HCN)  &
        -k(156)*n(idx_H2O)  &
        -k(649)*n(idx_H2O)  &
        -k(336)*n(idx_C)  &
        -k(158)*n(idx_NO)  &
        -k(472)*n(idx_H2)  &
        -k(645)*n(idx_CO2)  &
        -k(646)*n(idx_CO)  &
        -k(658)*n(idx_NH)  &
        -k(159)*n(idx_O2)  &
        -k(650)*n(idx_H2O)  &
        -k(655)*n(idx_HNC)  &
        -k(659)*n(idx_NO)  &
        -k(648)*n(idx_H2CO)  &
        -k(155)*n(idx_H2CO)  &
        -k(662)*n(idx_O)  &
        -k(1093)  &
        -k(412)*n(idx_CH)  &
        -k(647)*n(idx_H2CO)  &
        -k(643)*n(idx_CO2)

    !d[HNO+_dot]/d[NH+]
    pd(79,76) =  &
        +k(644)*n(idx_CO2)  &
        +k(650)*n(idx_H2O)

    !d[H3+_dot]/d[NH+]
    pd(81,76) =  &
        +k(472)*n(idx_H2)

    !d[H3CO+_dot]/d[NH+]
    pd(82,76) =  &
        +k(647)*n(idx_H2CO)

    !d[H3O+_dot]/d[NH+]
    pd(83,76) =  &
        +k(649)*n(idx_H2O)

    !d[HCNH+_dot]/d[NH+]
    pd(84,76) =  &
        +k(655)*n(idx_HNC)  &
        +k(653)*n(idx_HCN)

    !d[HCO2+_dot]/d[NH+]
    pd(85,76) =  &
        +k(643)*n(idx_CO2)

    !d[N2H+_dot]/d[NH+]
    pd(87,76) =  &
        +k(656)*n(idx_N2)  &
        +k(659)*n(idx_NO)

    !d[O2H+_dot]/d[NH+]
    pd(88,76) =  &
        +k(661)*n(idx_O2)

    !d[E_dot]/d[H2+]
    pd(1,77) =  &
        -k(270)*n(idx_E)

    !d[CH_dot]/d[H2+]
    pd(2,77) =  &
        -k(87)*n(idx_CH)  &
        -k(444)*n(idx_CH)

    !d[O_dot]/d[H2+]
    pd(3,77) =  &
        -k(458)*n(idx_O)

    !d[HCN_dot]/d[H2+]
    pd(5,77) =  &
        -k(92)*n(idx_HCN)

    !d[H2_dot]/d[H2+]
    pd(6,77) =  &
        -k(448)*n(idx_H2)  &
        +k(90)*n(idx_H2CO)  &
        +k(88)*n(idx_CN)  &
        +k(112)*n(idx_H)  &
        +k(97)*n(idx_NO)  &
        +k(85)*n(idx_CH2)  &
        +k(92)*n(idx_HCN)  &
        +k(89)*n(idx_CO)  &
        +k(91)*n(idx_H2O)  &
        +k(443)*n(idx_CH4)  &
        +k(94)*n(idx_NH2)  &
        +k(96)*n(idx_NH)  &
        +k(87)*n(idx_CH)  &
        +k(86)*n(idx_CH4)  &
        +k(449)*n(idx_H2CO)  &
        +k(99)*n(idx_OH)  &
        +k(98)*n(idx_O2)  &
        +k(95)*n(idx_NH3)  &
        +k(93)*n(idx_HCO)

    !d[C_dot]/d[H2+]
    pd(7,77) =  &
        -k(441)*n(idx_C)

    !d[H_dot]/d[H2+]
    pd(8,77) =  &
        +k(446)*n(idx_CO2)  &
        +k(445)*n(idx_CN)  &
        +k(450)*n(idx_H2O)  &
        +k(455)*n(idx_NH)  &
        +2.d0*k(270)*n(idx_E)  &
        +k(448)*n(idx_H2)  &
        +k(454)*n(idx_N)  &
        +k(441)*n(idx_C)  &
        +k(453)*n(idx_N2)  &
        +k(1000)  &
        +k(458)*n(idx_O)  &
        +k(444)*n(idx_CH)  &
        +k(449)*n(idx_H2CO)  &
        +k(452)*n(idx_HE)  &
        +k(457)*n(idx_O2)  &
        +k(459)*n(idx_OH)  &
        +k(442)*n(idx_CH2)  &
        +k(443)*n(idx_CH4)  &
        +k(447)*n(idx_CO)  &
        -k(112)*n(idx_H)  &
        +k(456)*n(idx_NO)

    !d[H2O_dot]/d[H2+]
    pd(9,77) =  &
        -k(91)*n(idx_H2O)  &
        -k(450)*n(idx_H2O)

    !d[OH_dot]/d[H2+]
    pd(10,77) =  &
        -k(459)*n(idx_OH)  &
        -k(99)*n(idx_OH)

    !d[O2_dot]/d[H2+]
    pd(11,77) =  &
        -k(98)*n(idx_O2)  &
        -k(457)*n(idx_O2)

    !d[CH2_dot]/d[H2+]
    pd(12,77) =  &
        -k(85)*n(idx_CH2)  &
        -k(442)*n(idx_CH2)

    !d[H2CO_dot]/d[H2+]
    pd(13,77) =  &
        -k(449)*n(idx_H2CO)  &
        -k(90)*n(idx_H2CO)

    !d[HCO_dot]/d[H2+]
    pd(14,77) =  &
        -k(451)*n(idx_HCO)  &
        -k(93)*n(idx_HCO)

    !d[NH3_dot]/d[H2+]
    pd(16,77) =  &
        -k(95)*n(idx_NH3)

    !d[NO_dot]/d[H2+]
    pd(17,77) =  &
        -k(97)*n(idx_NO)  &
        -k(456)*n(idx_NO)

    !d[CN_dot]/d[H2+]
    pd(18,77) =  &
        -k(88)*n(idx_CN)  &
        -k(445)*n(idx_CN)

    !d[CO_dot]/d[H2+]
    pd(19,77) =  &
        -k(447)*n(idx_CO)  &
        -k(89)*n(idx_CO)  &
        +k(451)*n(idx_HCO)

    !d[N2_dot]/d[H2+]
    pd(20,77) =  &
        -k(453)*n(idx_N2)

    !d[NH2_dot]/d[H2+]
    pd(21,77) =  &
        -k(94)*n(idx_NH2)

    !d[CH4_dot]/d[H2+]
    pd(23,77) =  &
        -k(86)*n(idx_CH4)  &
        -k(443)*n(idx_CH4)

    !d[N_dot]/d[H2+]
    pd(24,77) =  &
        -k(454)*n(idx_N)

    !d[NH_dot]/d[H2+]
    pd(25,77) =  &
        -k(96)*n(idx_NH)  &
        -k(455)*n(idx_NH)

    !d[HE_dot]/d[H2+]
    pd(26,77) =  &
        -k(452)*n(idx_HE)

    !d[CO2_dot]/d[H2+]
    pd(29,77) =  &
        -k(446)*n(idx_CO2)

    !d[HCO+_dot]/d[H2+]
    pd(54,77) =  &
        +k(449)*n(idx_H2CO)  &
        +k(447)*n(idx_CO)  &
        +k(93)*n(idx_HCO)

    !d[H+_dot]/d[H2+]
    pd(55,77) =  &
        +k(1000)  &
        +k(112)*n(idx_H)

    !d[CH2+_dot]/d[H2+]
    pd(58,77) =  &
        +k(444)*n(idx_CH)  &
        +k(85)*n(idx_CH2)

    !d[CH+_dot]/d[H2+]
    pd(59,77) =  &
        +k(441)*n(idx_C)  &
        +k(87)*n(idx_CH)

    !d[H2CO+_dot]/d[H2+]
    pd(60,77) =  &
        +k(90)*n(idx_H2CO)

    !d[NH3+_dot]/d[H2+]
    pd(62,77) =  &
        +k(95)*n(idx_NH3)

    !d[NO+_dot]/d[H2+]
    pd(63,77) =  &
        +k(97)*n(idx_NO)

    !d[CN+_dot]/d[H2+]
    pd(64,77) =  &
        +k(88)*n(idx_CN)

    !d[CO+_dot]/d[H2+]
    pd(65,77) =  &
        +k(89)*n(idx_CO)

    !d[O2+_dot]/d[H2+]
    pd(67,77) =  &
        +k(98)*n(idx_O2)

    !d[H2O+_dot]/d[H2+]
    pd(68,77) =  &
        +k(91)*n(idx_H2O)  &
        +k(459)*n(idx_OH)

    !d[NH2+_dot]/d[H2+]
    pd(69,77) =  &
        +k(94)*n(idx_NH2)  &
        +k(455)*n(idx_NH)

    !d[OH+_dot]/d[H2+]
    pd(71,77) =  &
        +k(99)*n(idx_OH)  &
        +k(458)*n(idx_O)

    !d[CH3+_dot]/d[H2+]
    pd(72,77) =  &
        +k(443)*n(idx_CH4)  &
        +k(442)*n(idx_CH2)

    !d[CH4+_dot]/d[H2+]
    pd(73,77) =  &
        +k(86)*n(idx_CH4)

    !d[HCN+_dot]/d[H2+]
    pd(75,77) =  &
        +k(445)*n(idx_CN)  &
        +k(92)*n(idx_HCN)

    !d[NH+_dot]/d[H2+]
    pd(76,77) =  &
        +k(454)*n(idx_N)  &
        +k(96)*n(idx_NH)

    !d[H2+_dot]/d[H2+]
    pd(77,77) =  &
        -k(455)*n(idx_NH)  &
        -k(449)*n(idx_H2CO)  &
        -k(442)*n(idx_CH2)  &
        -k(444)*n(idx_CH)  &
        -k(90)*n(idx_H2CO)  &
        -k(452)*n(idx_HE)  &
        -k(91)*n(idx_H2O)  &
        -k(457)*n(idx_O2)  &
        -k(112)*n(idx_H)  &
        -k(445)*n(idx_CN)  &
        -k(86)*n(idx_CH4)  &
        -k(270)*n(idx_E)  &
        -k(92)*n(idx_HCN)  &
        -k(94)*n(idx_NH2)  &
        -k(447)*n(idx_CO)  &
        -k(1000)  &
        -k(99)*n(idx_OH)  &
        -k(451)*n(idx_HCO)  &
        -k(88)*n(idx_CN)  &
        -k(446)*n(idx_CO2)  &
        -k(441)*n(idx_C)  &
        -k(97)*n(idx_NO)  &
        -k(459)*n(idx_OH)  &
        -k(453)*n(idx_N2)  &
        -k(95)*n(idx_NH3)  &
        -k(454)*n(idx_N)  &
        -k(93)*n(idx_HCO)  &
        -k(458)*n(idx_O)  &
        -k(87)*n(idx_CH)  &
        -k(448)*n(idx_H2)  &
        -k(85)*n(idx_CH2)  &
        -k(98)*n(idx_O2)  &
        -k(443)*n(idx_CH4)  &
        -k(456)*n(idx_NO)  &
        -k(96)*n(idx_NH)  &
        -k(450)*n(idx_H2O)  &
        -k(89)*n(idx_CO)

    !d[HNO+_dot]/d[H2+]
    pd(79,77) =  &
        +k(456)*n(idx_NO)

    !d[H3+_dot]/d[H2+]
    pd(81,77) =  &
        +k(448)*n(idx_H2)  &
        +k(451)*n(idx_HCO)

    !d[H3O+_dot]/d[H2+]
    pd(83,77) =  &
        +k(450)*n(idx_H2O)

    !d[HCO2+_dot]/d[H2+]
    pd(85,77) =  &
        +k(446)*n(idx_CO2)

    !d[HEH+_dot]/d[H2+]
    pd(86,77) =  &
        +k(452)*n(idx_HE)

    !d[N2H+_dot]/d[H2+]
    pd(87,77) =  &
        +k(453)*n(idx_N2)

    !d[O2H+_dot]/d[H2+]
    pd(88,77) =  &
        +k(457)*n(idx_O2)

    !d[E_dot]/d[HE+]
    pd(1,78) =  &
        -k(1060)*n(idx_E)

    !d[CH_dot]/d[HE+]
    pd(2,78) =  &
        -k(124)*n(idx_CH)  &
        +k(587)*n(idx_HCN)  &
        -k(573)*n(idx_CH)

    !d[O_dot]/d[HE+]
    pd(3,78) =  &
        +k(592)*n(idx_HCO)  &
        +k(576)*n(idx_CO2)  &
        +k(607)*n(idx_OCN)  &
        +k(580)*n(idx_CO)  &
        +k(605)*n(idx_NO)  &
        +k(606)*n(idx_O2)  &
        +k(583)*n(idx_H2CO)

    !d[HNC_dot]/d[HE+]
    pd(4,78) =  &
        -k(595)*n(idx_HNC)  &
        -k(593)*n(idx_HNC)  &
        -k(594)*n(idx_HNC)

    !d[HCN_dot]/d[HE+]
    pd(5,78) =  &
        -k(586)*n(idx_HCN)  &
        -k(587)*n(idx_HCN)  &
        -k(588)*n(idx_HCN)  &
        -k(589)*n(idx_HCN)

    !d[H2_dot]/d[HE+]
    pd(6,78) =  &
        +k(564)*n(idx_CH2)  &
        +k(601)*n(idx_NH3)  &
        +k(566)*n(idx_CH3)  &
        -k(100)*n(idx_H2)  &
        +k(569)*n(idx_CH4)  &
        -k(468)*n(idx_H2)  &
        +k(599)*n(idx_NH2)  &
        +k(581)*n(idx_H2CO)  &
        +k(570)*n(idx_CH4)

    !d[C_dot]/d[HE+]
    pd(7,78) =  &
        -k(122)*n(idx_C)  &
        +k(595)*n(idx_HNC)  &
        +k(578)*n(idx_CO2)  &
        +k(574)*n(idx_CN)

    !d[H_dot]/d[HE+]
    pd(8,78) =  &
        -k(114)*n(idx_H)  &
        +k(596)*n(idx_HNO)  &
        +k(603)*n(idx_NH)  &
        +k(594)*n(idx_HNC)  &
        +k(609)*n(idx_OH)  &
        +k(602)*n(idx_NH3)  &
        +k(571)*n(idx_CH4)  &
        +k(600)*n(idx_NH2)  &
        +k(569)*n(idx_CH4)  &
        +k(584)*n(idx_H2O)  &
        +k(573)*n(idx_CH)  &
        +k(468)*n(idx_H2)  &
        +k(586)*n(idx_HCN)  &
        +k(593)*n(idx_HNC)  &
        +k(588)*n(idx_HCN)  &
        +k(565)*n(idx_CH2)  &
        +k(582)*n(idx_H2CO)  &
        +k(590)*n(idx_HCO)

    !d[H2O_dot]/d[HE+]
    pd(9,78) =  &
        -k(584)*n(idx_H2O)  &
        -k(126)*n(idx_H2O)  &
        -k(585)*n(idx_H2O)

    !d[OH_dot]/d[HE+]
    pd(10,78) =  &
        +k(568)*n(idx_CH3OH)  &
        +k(585)*n(idx_H2O)  &
        -k(609)*n(idx_OH)

    !d[O2_dot]/d[HE+]
    pd(11,78) =  &
        -k(606)*n(idx_O2)  &
        -k(129)*n(idx_O2)  &
        +k(579)*n(idx_CO2)

    !d[CH2_dot]/d[HE+]
    pd(12,78) =  &
        -k(564)*n(idx_CH2)  &
        -k(565)*n(idx_CH2)

    !d[H2CO_dot]/d[HE+]
    pd(13,78) =  &
        -k(582)*n(idx_H2CO)  &
        -k(581)*n(idx_H2CO)  &
        -k(125)*n(idx_H2CO)  &
        -k(583)*n(idx_H2CO)

    !d[HCO_dot]/d[HE+]
    pd(14,78) =  &
        -k(590)*n(idx_HCO)  &
        -k(591)*n(idx_HCO)  &
        -k(592)*n(idx_HCO)

    !d[NH3_dot]/d[HE+]
    pd(16,78) =  &
        -k(128)*n(idx_NH3)  &
        -k(602)*n(idx_NH3)  &
        -k(601)*n(idx_NH3)

    !d[NO_dot]/d[HE+]
    pd(17,78) =  &
        -k(604)*n(idx_NO)  &
        -k(605)*n(idx_NO)  &
        +k(597)*n(idx_HNO)

    !d[CN_dot]/d[HE+]
    pd(18,78) =  &
        +k(608)*n(idx_OCN)  &
        -k(575)*n(idx_CN)  &
        -k(574)*n(idx_CN)

    !d[CO_dot]/d[HE+]
    pd(19,78) =  &
        -k(580)*n(idx_CO)  &
        +k(591)*n(idx_HCO)  &
        +k(577)*n(idx_CO2)

    !d[N2_dot]/d[HE+]
    pd(20,78) =  &
        -k(127)*n(idx_N2)  &
        -k(598)*n(idx_N2)

    !d[NH2_dot]/d[HE+]
    pd(21,78) =  &
        -k(600)*n(idx_NH2)  &
        -k(599)*n(idx_NH2)

    !d[CH3_dot]/d[HE+]
    pd(22,78) =  &
        -k(566)*n(idx_CH3)  &
        +k(567)*n(idx_CH3OH)  &
        +k(572)*n(idx_CH4)

    !d[CH4_dot]/d[HE+]
    pd(23,78) =  &
        -k(123)*n(idx_CH4)  &
        -k(571)*n(idx_CH4)  &
        -k(572)*n(idx_CH4)  &
        -k(570)*n(idx_CH4)  &
        -k(569)*n(idx_CH4)

    !d[N_dot]/d[HE+]
    pd(24,78) =  &
        +k(589)*n(idx_HCN)  &
        +k(575)*n(idx_CN)  &
        +k(588)*n(idx_HCN)  &
        +k(594)*n(idx_HNC)  &
        +k(598)*n(idx_N2)  &
        +k(604)*n(idx_NO)

    !d[NH_dot]/d[HE+]
    pd(25,78) =  &
        -k(603)*n(idx_NH)

    !d[HE_dot]/d[HE+]
    pd(26,78) =  &
        +k(596)*n(idx_HNO)  &
        +k(603)*n(idx_NH)  &
        +k(587)*n(idx_HCN)  &
        +k(574)*n(idx_CN)  &
        +k(122)*n(idx_C)  &
        +k(608)*n(idx_OCN)  &
        +k(576)*n(idx_CO2)  &
        +k(573)*n(idx_CH)  &
        +k(126)*n(idx_H2O)  &
        +k(575)*n(idx_CN)  &
        +k(599)*n(idx_NH2)  &
        +k(588)*n(idx_HCN)  &
        +k(577)*n(idx_CO2)  &
        +k(579)*n(idx_CO2)  &
        +k(592)*n(idx_HCO)  &
        +k(123)*n(idx_CH4)  &
        +k(601)*n(idx_NH3)  &
        +k(590)*n(idx_HCO)  &
        +k(602)*n(idx_NH3)  &
        +k(571)*n(idx_CH4)  &
        +k(589)*n(idx_HCN)  &
        +k(600)*n(idx_NH2)  &
        +k(568)*n(idx_CH3OH)  &
        +k(572)*n(idx_CH4)  &
        +k(468)*n(idx_H2)  &
        +k(604)*n(idx_NO)  &
        +k(583)*n(idx_H2CO)  &
        +k(129)*n(idx_O2)  &
        +k(580)*n(idx_CO)  &
        +k(570)*n(idx_CH4)  &
        +k(124)*n(idx_CH)  &
        +k(586)*n(idx_HCN)  &
        +k(114)*n(idx_H)  &
        +k(597)*n(idx_HNO)  &
        +k(125)*n(idx_H2CO)  &
        +k(609)*n(idx_OH)  &
        +k(607)*n(idx_OCN)  &
        +k(1060)*n(idx_E)  &
        +k(593)*n(idx_HNC)  &
        +k(100)*n(idx_H2)  &
        +k(127)*n(idx_N2)  &
        +k(565)*n(idx_CH2)  &
        +k(582)*n(idx_H2CO)  &
        +k(605)*n(idx_NO)  &
        +k(578)*n(idx_CO2)  &
        +k(564)*n(idx_CH2)  &
        +k(595)*n(idx_HNC)  &
        +k(566)*n(idx_CH3)  &
        +k(569)*n(idx_CH4)  &
        +k(584)*n(idx_H2O)  &
        +k(567)*n(idx_CH3OH)  &
        +k(585)*n(idx_H2O)  &
        +k(581)*n(idx_H2CO)  &
        +k(128)*n(idx_NH3)  &
        +k(594)*n(idx_HNC)  &
        +k(598)*n(idx_N2)  &
        +k(606)*n(idx_O2)

    !d[HNO_dot]/d[HE+]
    pd(27,78) =  &
        -k(597)*n(idx_HNO)  &
        -k(596)*n(idx_HNO)

    !d[CH3OH_dot]/d[HE+]
    pd(28,78) =  &
        -k(567)*n(idx_CH3OH)  &
        -k(568)*n(idx_CH3OH)

    !d[CO2_dot]/d[HE+]
    pd(29,78) =  &
        -k(578)*n(idx_CO2)  &
        -k(579)*n(idx_CO2)  &
        -k(576)*n(idx_CO2)  &
        -k(577)*n(idx_CO2)

    !d[OCN_dot]/d[HE+]
    pd(34,78) =  &
        -k(607)*n(idx_OCN)  &
        -k(608)*n(idx_OCN)

    !d[HCO+_dot]/d[HE+]
    pd(54,78) =  &
        +k(582)*n(idx_H2CO)

    !d[H+_dot]/d[HE+]
    pd(55,78) =  &
        +k(572)*n(idx_CH4)  &
        +k(114)*n(idx_H)  &
        +k(468)*n(idx_H2)  &
        +k(597)*n(idx_HNO)  &
        +k(585)*n(idx_H2O)

    !d[C+_dot]/d[HE+]
    pd(57,78) =  &
        +k(564)*n(idx_CH2)  &
        +k(122)*n(idx_C)  &
        +k(573)*n(idx_CH)  &
        +k(575)*n(idx_CN)  &
        +k(588)*n(idx_HCN)  &
        +k(594)*n(idx_HNC)  &
        +k(580)*n(idx_CO)  &
        +k(579)*n(idx_CO2)

    !d[CH2+_dot]/d[HE+]
    pd(58,78) =  &
        +k(583)*n(idx_H2CO)  &
        +k(570)*n(idx_CH4)

    !d[CH+_dot]/d[HE+]
    pd(59,78) =  &
        +k(592)*n(idx_HCO)  &
        +k(124)*n(idx_CH)  &
        +k(566)*n(idx_CH3)  &
        +k(589)*n(idx_HCN)  &
        +k(569)*n(idx_CH4)  &
        +k(565)*n(idx_CH2)

    !d[H2CO+_dot]/d[HE+]
    pd(60,78) =  &
        +k(125)*n(idx_H2CO)

    !d[NH3+_dot]/d[HE+]
    pd(62,78) =  &
        +k(128)*n(idx_NH3)

    !d[NO+_dot]/d[HE+]
    pd(63,78) =  &
        +k(596)*n(idx_HNO)

    !d[CN+_dot]/d[HE+]
    pd(64,78) =  &
        +k(593)*n(idx_HNC)  &
        +k(586)*n(idx_HCN)  &
        +k(607)*n(idx_OCN)

    !d[CO+_dot]/d[HE+]
    pd(65,78) =  &
        +k(581)*n(idx_H2CO)  &
        +k(576)*n(idx_CO2)  &
        +k(590)*n(idx_HCO)

    !d[N2+_dot]/d[HE+]
    pd(66,78) =  &
        +k(127)*n(idx_N2)

    !d[O2+_dot]/d[HE+]
    pd(67,78) =  &
        +k(129)*n(idx_O2)  &
        +k(578)*n(idx_CO2)

    !d[H2O+_dot]/d[HE+]
    pd(68,78) =  &
        +k(126)*n(idx_H2O)

    !d[NH2+_dot]/d[HE+]
    pd(69,78) =  &
        +k(602)*n(idx_NH3)

    !d[O+_dot]/d[HE+]
    pd(70,78) =  &
        +k(604)*n(idx_NO)  &
        +k(608)*n(idx_OCN)  &
        +k(577)*n(idx_CO2)  &
        +k(606)*n(idx_O2)  &
        +k(609)*n(idx_OH)

    !d[OH+_dot]/d[HE+]
    pd(71,78) =  &
        +k(584)*n(idx_H2O)  &
        +k(567)*n(idx_CH3OH)

    !d[CH3+_dot]/d[HE+]
    pd(72,78) =  &
        +k(568)*n(idx_CH3OH)  &
        +k(571)*n(idx_CH4)

    !d[CH4+_dot]/d[HE+]
    pd(73,78) =  &
        +k(123)*n(idx_CH4)

    !d[N+_dot]/d[HE+]
    pd(74,78) =  &
        +k(603)*n(idx_NH)  &
        +k(587)*n(idx_HCN)  &
        +k(574)*n(idx_CN)  &
        +k(598)*n(idx_N2)  &
        +k(599)*n(idx_NH2)  &
        +k(605)*n(idx_NO)

    !d[NH+_dot]/d[HE+]
    pd(76,78) =  &
        +k(600)*n(idx_NH2)  &
        +k(601)*n(idx_NH3)  &
        +k(595)*n(idx_HNC)

    !d[H2+_dot]/d[HE+]
    pd(77,78) =  &
        +k(100)*n(idx_H2)

    !d[HE+_dot]/d[HE+]
    pd(78,78) =  &
        -k(114)*n(idx_H)  &
        -k(604)*n(idx_NO)  &
        -k(584)*n(idx_H2O)  &
        -k(591)*n(idx_HCO)  &
        -k(598)*n(idx_N2)  &
        -k(590)*n(idx_HCO)  &
        -k(576)*n(idx_CO2)  &
        -k(602)*n(idx_NH3)  &
        -k(583)*n(idx_H2CO)  &
        -k(125)*n(idx_H2CO)  &
        -k(129)*n(idx_O2)  &
        -k(124)*n(idx_CH)  &
        -k(593)*n(idx_HNC)  &
        -k(596)*n(idx_HNO)  &
        -k(580)*n(idx_CO)  &
        -k(570)*n(idx_CH4)  &
        -k(597)*n(idx_HNO)  &
        -k(577)*n(idx_CO2)  &
        -k(578)*n(idx_CO2)  &
        -k(123)*n(idx_CH4)  &
        -k(605)*n(idx_NO)  &
        -k(585)*n(idx_H2O)  &
        -k(564)*n(idx_CH2)  &
        -k(609)*n(idx_OH)  &
        -k(100)*n(idx_H2)  &
        -k(579)*n(idx_CO2)  &
        -k(581)*n(idx_H2CO)  &
        -k(1060)*n(idx_E)  &
        -k(468)*n(idx_H2)  &
        -k(608)*n(idx_OCN)  &
        -k(600)*n(idx_NH2)  &
        -k(565)*n(idx_CH2)  &
        -k(568)*n(idx_CH3OH)  &
        -k(572)*n(idx_CH4)  &
        -k(595)*n(idx_HNC)  &
        -k(588)*n(idx_HCN)  &
        -k(574)*n(idx_CN)  &
        -k(586)*n(idx_HCN)  &
        -k(126)*n(idx_H2O)  &
        -k(569)*n(idx_CH4)  &
        -k(587)*n(idx_HCN)  &
        -k(594)*n(idx_HNC)  &
        -k(582)*n(idx_H2CO)  &
        -k(573)*n(idx_CH)  &
        -k(589)*n(idx_HCN)  &
        -k(566)*n(idx_CH3)  &
        -k(122)*n(idx_C)  &
        -k(575)*n(idx_CN)  &
        -k(592)*n(idx_HCO)  &
        -k(128)*n(idx_NH3)  &
        -k(606)*n(idx_O2)  &
        -k(601)*n(idx_NH3)  &
        -k(599)*n(idx_NH2)  &
        -k(571)*n(idx_CH4)  &
        -k(607)*n(idx_OCN)  &
        -k(603)*n(idx_NH)  &
        -k(127)*n(idx_N2)  &
        -k(567)*n(idx_CH3OH)

    !d[HEH+_dot]/d[HE+]
    pd(86,78) =  &
        +k(591)*n(idx_HCO)

    !d[E_dot]/d[HNO+]
    pd(1,79) =  &
        -k(299)*n(idx_E)

    !d[CH_dot]/d[HNO+]
    pd(2,79) =  &
        -k(409)*n(idx_CH)

    !d[HNC_dot]/d[HNO+]
    pd(4,79) =  &
        -k(560)*n(idx_HNC)

    !d[HCN_dot]/d[HNO+]
    pd(5,79) =  &
        -k(545)*n(idx_HCN)

    !d[C_dot]/d[HNO+]
    pd(7,79) =  &
        -k(334)*n(idx_C)

    !d[H_dot]/d[HNO+]
    pd(8,79) =  &
        +k(299)*n(idx_E)

    !d[H2O_dot]/d[HNO+]
    pd(9,79) =  &
        -k(498)*n(idx_H2O)

    !d[OH_dot]/d[HNO+]
    pd(10,79) =  &
        -k(743)*n(idx_OH)

    !d[CH2_dot]/d[HNO+]
    pd(12,79) =  &
        -k(374)*n(idx_CH2)

    !d[H2CO_dot]/d[HNO+]
    pd(13,79) =  &
        -k(480)*n(idx_H2CO)

    !d[HCO_dot]/d[HNO+]
    pd(14,79) =  &
        -k(553)*n(idx_HCO)

    !d[NO_dot]/d[HNO+]
    pd(17,79) =  &
        -k(183)*n(idx_NO)  &
        +k(480)*n(idx_H2CO)  &
        +k(545)*n(idx_HCN)  &
        +k(498)*n(idx_H2O)  &
        +k(299)*n(idx_E)  &
        +k(553)*n(idx_HCO)  &
        +k(630)*n(idx_N2)  &
        +k(421)*n(idx_CN)  &
        +k(560)*n(idx_HNC)  &
        +k(374)*n(idx_CH2)  &
        +k(743)*n(idx_OH)  &
        +k(695)*n(idx_NH)  &
        +k(683)*n(idx_NH2)  &
        +k(425)*n(idx_CO)  &
        +k(563)*n(idx_CO2)  &
        +k(409)*n(idx_CH)  &
        +k(334)*n(idx_C)

    !d[CN_dot]/d[HNO+]
    pd(18,79) =  &
        -k(421)*n(idx_CN)

    !d[CO_dot]/d[HNO+]
    pd(19,79) =  &
        -k(425)*n(idx_CO)

    !d[N2_dot]/d[HNO+]
    pd(20,79) =  &
        -k(630)*n(idx_N2)

    !d[NH2_dot]/d[HNO+]
    pd(21,79) =  &
        -k(683)*n(idx_NH2)

    !d[NH_dot]/d[HNO+]
    pd(25,79) =  &
        -k(695)*n(idx_NH)

    !d[HNO_dot]/d[HNO+]
    pd(27,79) =  &
        +k(183)*n(idx_NO)

    !d[CO2_dot]/d[HNO+]
    pd(29,79) =  &
        -k(563)*n(idx_CO2)

    !d[HNO_DUST_dot]/d[HNO+]
    pd(48,79) =  &
        +k(1115)

    !d[HCO+_dot]/d[HNO+]
    pd(54,79) =  &
        +k(425)*n(idx_CO)

    !d[CH2+_dot]/d[HNO+]
    pd(58,79) =  &
        +k(409)*n(idx_CH)

    !d[CH+_dot]/d[HNO+]
    pd(59,79) =  &
        +k(334)*n(idx_C)

    !d[H2CO+_dot]/d[HNO+]
    pd(60,79) =  &
        +k(553)*n(idx_HCO)

    !d[NH3+_dot]/d[HNO+]
    pd(62,79) =  &
        +k(683)*n(idx_NH2)

    !d[NO+_dot]/d[HNO+]
    pd(63,79) =  &
        +k(183)*n(idx_NO)

    !d[H2O+_dot]/d[HNO+]
    pd(68,79) =  &
        +k(743)*n(idx_OH)

    !d[NH2+_dot]/d[HNO+]
    pd(69,79) =  &
        +k(695)*n(idx_NH)

    !d[CH3+_dot]/d[HNO+]
    pd(72,79) =  &
        +k(374)*n(idx_CH2)

    !d[HCN+_dot]/d[HNO+]
    pd(75,79) =  &
        +k(421)*n(idx_CN)

    !d[HNO+_dot]/d[HNO+]
    pd(79,79) =  &
        -k(183)*n(idx_NO)  &
        -k(421)*n(idx_CN)  &
        -k(425)*n(idx_CO)  &
        -k(695)*n(idx_NH)  &
        -k(545)*n(idx_HCN)  &
        -k(299)*n(idx_E)  &
        -k(553)*n(idx_HCO)  &
        -k(560)*n(idx_HNC)  &
        -k(1115)  &
        -k(498)*n(idx_H2O)  &
        -k(409)*n(idx_CH)  &
        -k(334)*n(idx_C)  &
        -k(563)*n(idx_CO2)  &
        -k(374)*n(idx_CH2)  &
        -k(480)*n(idx_H2CO)  &
        -k(630)*n(idx_N2)  &
        -k(743)*n(idx_OH)  &
        -k(683)*n(idx_NH2)

    !d[H3CO+_dot]/d[HNO+]
    pd(82,79) =  &
        +k(480)*n(idx_H2CO)

    !d[H3O+_dot]/d[HNO+]
    pd(83,79) =  &
        +k(498)*n(idx_H2O)

    !d[HCNH+_dot]/d[HNO+]
    pd(84,79) =  &
        +k(545)*n(idx_HCN)  &
        +k(560)*n(idx_HNC)

    !d[HCO2+_dot]/d[HNO+]
    pd(85,79) =  &
        +k(563)*n(idx_CO2)

    !d[N2H+_dot]/d[HNO+]
    pd(87,79) =  &
        +k(630)*n(idx_N2)

    !d[E_dot]/d[H2NO+]
    pd(1,80) =  &
        -k(276)*n(idx_E)  &
        -k(275)*n(idx_E)

    !d[H2_dot]/d[H2NO+]
    pd(6,80) =  &
        +k(276)*n(idx_E)

    !d[H_dot]/d[H2NO+]
    pd(8,80) =  &
        +k(275)*n(idx_E)  &
        +k(1116)

    !d[NO_dot]/d[H2NO+]
    pd(17,80) =  &
        +k(276)*n(idx_E)

    !d[HNO_dot]/d[H2NO+]
    pd(27,80) =  &
        +k(275)*n(idx_E)

    !d[HNO_DUST_dot]/d[H2NO+]
    pd(48,80) =  &
        +k(1116)

    !d[H2NO+_dot]/d[H2NO+]
    pd(80,80) =  &
        -k(276)*n(idx_E)  &
        -k(1116)  &
        -k(275)*n(idx_E)

    !d[E_dot]/d[H3+]
    pd(1,81) =  &
        -k(281)*n(idx_E)  &
        -k(280)*n(idx_E)

    !d[CH_dot]/d[H3+]
    pd(2,81) =  &
        -k(506)*n(idx_CH)

    !d[O_dot]/d[H3+]
    pd(3,81) =  &
        -k(524)*n(idx_O)  &
        -k(525)*n(idx_O)

    !d[HNC_dot]/d[H3+]
    pd(4,81) =  &
        -k(515)*n(idx_HNC)

    !d[HCN_dot]/d[H3+]
    pd(5,81) =  &
        -k(513)*n(idx_HCN)

    !d[H2_dot]/d[H3+]
    pd(6,81) =  &
        +k(512)*n(idx_H2O)  &
        +k(514)*n(idx_HCO)  &
        +k(506)*n(idx_CH)  &
        +k(508)*n(idx_CO2)  &
        +k(523)*n(idx_O2)  &
        +k(522)*n(idx_NO)  &
        +k(511)*n(idx_H2CO)  &
        +k(509)*n(idx_CO)  &
        +k(505)*n(idx_CH3OH)  &
        +k(516)*n(idx_HNO)  &
        +k(502)*n(idx_C)  &
        +k(503)*n(idx_CH2)  &
        +k(280)*n(idx_E)  &
        +k(521)*n(idx_NO2)  &
        +k(525)*n(idx_O)  &
        +k(520)*n(idx_NH)  &
        +k(510)*n(idx_CO)  &
        +k(507)*n(idx_CN)  &
        +k(517)*n(idx_MG)  &
        +k(519)*n(idx_NH2)  &
        +k(504)*n(idx_CH3)  &
        +k(1010)  &
        +k(515)*n(idx_HNC)  &
        +k(513)*n(idx_HCN)  &
        +k(526)*n(idx_OH)  &
        +k(518)*n(idx_N2)

    !d[C_dot]/d[H3+]
    pd(7,81) =  &
        -k(502)*n(idx_C)

    !d[H_dot]/d[H3+]
    pd(8,81) =  &
        +3.d0*k(281)*n(idx_E)  &
        +k(280)*n(idx_E)  &
        +k(524)*n(idx_O)  &
        +k(517)*n(idx_MG)  &
        +k(1009)

    !d[H2O_dot]/d[H3+]
    pd(9,81) =  &
        -k(512)*n(idx_H2O)  &
        +k(505)*n(idx_CH3OH)

    !d[OH_dot]/d[H3+]
    pd(10,81) =  &
        -k(526)*n(idx_OH)  &
        +k(521)*n(idx_NO2)

    !d[O2_dot]/d[H3+]
    pd(11,81) =  &
        -k(523)*n(idx_O2)

    !d[CH2_dot]/d[H3+]
    pd(12,81) =  &
        -k(503)*n(idx_CH2)

    !d[H2CO_dot]/d[H3+]
    pd(13,81) =  &
        -k(511)*n(idx_H2CO)

    !d[HCO_dot]/d[H3+]
    pd(14,81) =  &
        -k(514)*n(idx_HCO)

    !d[MG_dot]/d[H3+]
    pd(15,81) =  &
        -k(517)*n(idx_MG)

    !d[NO_dot]/d[H3+]
    pd(17,81) =  &
        -k(522)*n(idx_NO)

    !d[CN_dot]/d[H3+]
    pd(18,81) =  &
        -k(507)*n(idx_CN)

    !d[CO_dot]/d[H3+]
    pd(19,81) =  &
        -k(510)*n(idx_CO)  &
        -k(509)*n(idx_CO)

    !d[N2_dot]/d[H3+]
    pd(20,81) =  &
        -k(518)*n(idx_N2)

    !d[NH2_dot]/d[H3+]
    pd(21,81) =  &
        -k(519)*n(idx_NH2)

    !d[CH3_dot]/d[H3+]
    pd(22,81) =  &
        -k(504)*n(idx_CH3)

    !d[NH_dot]/d[H3+]
    pd(25,81) =  &
        -k(520)*n(idx_NH)

    !d[HNO_dot]/d[H3+]
    pd(27,81) =  &
        -k(516)*n(idx_HNO)

    !d[CH3OH_dot]/d[H3+]
    pd(28,81) =  &
        -k(505)*n(idx_CH3OH)

    !d[CO2_dot]/d[H3+]
    pd(29,81) =  &
        -k(508)*n(idx_CO2)

    !d[NO2_dot]/d[H3+]
    pd(32,81) =  &
        -k(521)*n(idx_NO2)

    !d[HCO+_dot]/d[H3+]
    pd(54,81) =  &
        +k(509)*n(idx_CO)

    !d[H+_dot]/d[H3+]
    pd(55,81) =  &
        +k(1010)

    !d[HOC+_dot]/d[H3+]
    pd(56,81) =  &
        +k(510)*n(idx_CO)

    !d[CH2+_dot]/d[H3+]
    pd(58,81) =  &
        +k(506)*n(idx_CH)

    !d[CH+_dot]/d[H3+]
    pd(59,81) =  &
        +k(502)*n(idx_C)

    !d[H2CO+_dot]/d[H3+]
    pd(60,81) =  &
        +k(514)*n(idx_HCO)

    !d[MG+_dot]/d[H3+]
    pd(61,81) =  &
        +k(517)*n(idx_MG)

    !d[NH3+_dot]/d[H3+]
    pd(62,81) =  &
        +k(519)*n(idx_NH2)

    !d[NO+_dot]/d[H3+]
    pd(63,81) =  &
        +k(521)*n(idx_NO2)

    !d[H2O+_dot]/d[H3+]
    pd(68,81) =  &
        +k(524)*n(idx_O)  &
        +k(526)*n(idx_OH)

    !d[NH2+_dot]/d[H3+]
    pd(69,81) =  &
        +k(520)*n(idx_NH)

    !d[OH+_dot]/d[H3+]
    pd(71,81) =  &
        +k(525)*n(idx_O)

    !d[CH3+_dot]/d[H3+]
    pd(72,81) =  &
        +k(503)*n(idx_CH2)  &
        +k(505)*n(idx_CH3OH)

    !d[CH4+_dot]/d[H3+]
    pd(73,81) =  &
        +k(504)*n(idx_CH3)

    !d[HCN+_dot]/d[H3+]
    pd(75,81) =  &
        +k(507)*n(idx_CN)

    !d[H2+_dot]/d[H3+]
    pd(77,81) =  &
        +k(1009)

    !d[HNO+_dot]/d[H3+]
    pd(79,81) =  &
        +k(522)*n(idx_NO)

    !d[H2NO+_dot]/d[H3+]
    pd(80,81) =  &
        +k(516)*n(idx_HNO)

    !d[H3+_dot]/d[H3+]
    pd(81,81) =  &
        -k(507)*n(idx_CN)  &
        -k(506)*n(idx_CH)  &
        -k(515)*n(idx_HNC)  &
        -k(513)*n(idx_HCN)  &
        -k(524)*n(idx_O)  &
        -k(280)*n(idx_E)  &
        -k(520)*n(idx_NH)  &
        -k(505)*n(idx_CH3OH)  &
        -k(512)*n(idx_H2O)  &
        -k(1010)  &
        -k(516)*n(idx_HNO)  &
        -k(523)*n(idx_O2)  &
        -k(281)*n(idx_E)  &
        -k(522)*n(idx_NO)  &
        -k(503)*n(idx_CH2)  &
        -k(1009)  &
        -k(511)*n(idx_H2CO)  &
        -k(509)*n(idx_CO)  &
        -k(521)*n(idx_NO2)  &
        -k(517)*n(idx_MG)  &
        -k(519)*n(idx_NH2)  &
        -k(514)*n(idx_HCO)  &
        -k(525)*n(idx_O)  &
        -k(526)*n(idx_OH)  &
        -k(510)*n(idx_CO)  &
        -k(508)*n(idx_CO2)  &
        -k(502)*n(idx_C)  &
        -k(504)*n(idx_CH3)  &
        -k(518)*n(idx_N2)

    !d[H3CO+_dot]/d[H3+]
    pd(82,81) =  &
        +k(511)*n(idx_H2CO)

    !d[H3O+_dot]/d[H3+]
    pd(83,81) =  &
        +k(512)*n(idx_H2O)

    !d[HCNH+_dot]/d[H3+]
    pd(84,81) =  &
        +k(513)*n(idx_HCN)  &
        +k(515)*n(idx_HNC)

    !d[HCO2+_dot]/d[H3+]
    pd(85,81) =  &
        +k(508)*n(idx_CO2)

    !d[N2H+_dot]/d[H3+]
    pd(87,81) =  &
        +k(518)*n(idx_N2)

    !d[O2H+_dot]/d[H3+]
    pd(88,81) =  &
        +k(523)*n(idx_O2)

    !d[E_dot]/d[H3CO+]
    pd(1,82) =  &
        -k(286)*n(idx_E)  &
        -k(283)*n(idx_E)  &
        -k(282)*n(idx_E)  &
        -k(285)*n(idx_E)  &
        -k(284)*n(idx_E)

    !d[CH_dot]/d[H3CO+]
    pd(2,82) =  &
        -k(403)*n(idx_CH)  &
        +k(283)*n(idx_E)

    !d[HNC_dot]/d[H3CO+]
    pd(4,82) =  &
        -k(558)*n(idx_HNC)

    !d[HCN_dot]/d[H3CO+]
    pd(5,82) =  &
        -k(543)*n(idx_HCN)

    !d[H2_dot]/d[H3CO+]
    pd(6,82) =  &
        +k(284)*n(idx_E)

    !d[H_dot]/d[H3CO+]
    pd(8,82) =  &
        +2.d0*k(286)*n(idx_E)  &
        +k(285)*n(idx_E)  &
        +k(284)*n(idx_E)

    !d[H2O_dot]/d[H3CO+]
    pd(9,82) =  &
        -k(494)*n(idx_H2O)  &
        +k(283)*n(idx_E)

    !d[OH_dot]/d[H3CO+]
    pd(10,82) =  &
        +k(282)*n(idx_E)

    !d[CH2_dot]/d[H3CO+]
    pd(12,82) =  &
        +k(282)*n(idx_E)

    !d[H2CO_dot]/d[H3CO+]
    pd(13,82) =  &
        +k(558)*n(idx_HNC)  &
        +k(543)*n(idx_HCN)  &
        +k(285)*n(idx_E)  &
        +k(403)*n(idx_CH)  &
        +k(494)*n(idx_H2O)  &
        +k(677)*n(idx_NH2)

    !d[HCO_dot]/d[H3CO+]
    pd(14,82) =  &
        +k(286)*n(idx_E)

    !d[CO_dot]/d[H3CO+]
    pd(19,82) =  &
        +k(284)*n(idx_E)

    !d[NH2_dot]/d[H3CO+]
    pd(21,82) =  &
        -k(677)*n(idx_NH2)

    !d[CH3OH_DUST_dot]/d[H3CO+]
    pd(35,82) =  &
        +k(1069)

    !d[CH2+_dot]/d[H3CO+]
    pd(58,82) =  &
        +k(403)*n(idx_CH)

    !d[NH3+_dot]/d[H3CO+]
    pd(62,82) =  &
        +k(677)*n(idx_NH2)

    !d[H3CO+_dot]/d[H3CO+]
    pd(82,82) =  &
        -k(543)*n(idx_HCN)  &
        -k(677)*n(idx_NH2)  &
        -k(494)*n(idx_H2O)  &
        -k(403)*n(idx_CH)  &
        -k(558)*n(idx_HNC)  &
        -k(285)*n(idx_E)  &
        -k(1069)  &
        -k(286)*n(idx_E)  &
        -k(283)*n(idx_E)  &
        -k(282)*n(idx_E)  &
        -k(284)*n(idx_E)

    !d[H3O+_dot]/d[H3CO+]
    pd(83,82) =  &
        +k(494)*n(idx_H2O)

    !d[HCNH+_dot]/d[H3CO+]
    pd(84,82) =  &
        +k(558)*n(idx_HNC)  &
        +k(543)*n(idx_HCN)

    !d[E_dot]/d[H3O+]
    pd(1,83) =  &
        -k(287)*n(idx_E)  &
        -k(290)*n(idx_E)  &
        -k(289)*n(idx_E)  &
        -k(288)*n(idx_E)

    !d[CH_dot]/d[H3O+]
    pd(2,83) =  &
        -k(404)*n(idx_CH)

    !d[O_dot]/d[H3O+]
    pd(3,83) =  &
        +k(288)*n(idx_E)

    !d[HNC_dot]/d[H3O+]
    pd(4,83) =  &
        -k(529)*n(idx_HNC)

    !d[HCN_dot]/d[H3O+]
    pd(5,83) =  &
        -k(528)*n(idx_HCN)

    !d[H2_dot]/d[H3O+]
    pd(6,83) =  &
        +k(289)*n(idx_E)  &
        +k(288)*n(idx_E)  &
        +k(330)*n(idx_C)

    !d[C_dot]/d[H3O+]
    pd(7,83) =  &
        -k(330)*n(idx_C)

    !d[H_dot]/d[H3O+]
    pd(8,83) =  &
        +k(288)*n(idx_E)  &
        +k(287)*n(idx_E)  &
        +k(1106)  &
        +2.d0*k(290)*n(idx_E)

    !d[H2O_dot]/d[H3O+]
    pd(9,83) =  &
        +k(528)*n(idx_HCN)  &
        +k(678)*n(idx_NH2)  &
        +k(287)*n(idx_E)  &
        +k(527)*n(idx_H2CO)  &
        +k(369)*n(idx_CH2)  &
        +k(529)*n(idx_HNC)  &
        +k(404)*n(idx_CH)

    !d[OH_dot]/d[H3O+]
    pd(10,83) =  &
        +k(289)*n(idx_E)  &
        +k(290)*n(idx_E)

    !d[CH2_dot]/d[H3O+]
    pd(12,83) =  &
        -k(369)*n(idx_CH2)

    !d[H2CO_dot]/d[H3O+]
    pd(13,83) =  &
        -k(527)*n(idx_H2CO)

    !d[NH2_dot]/d[H3O+]
    pd(21,83) =  &
        -k(678)*n(idx_NH2)

    !d[H2O_DUST_dot]/d[H3O+]
    pd(40,83) =  &
        +k(1106)

    !d[HCO+_dot]/d[H3O+]
    pd(54,83) =  &
        +k(330)*n(idx_C)

    !d[CH2+_dot]/d[H3O+]
    pd(58,83) =  &
        +k(404)*n(idx_CH)

    !d[NH3+_dot]/d[H3O+]
    pd(62,83) =  &
        +k(678)*n(idx_NH2)

    !d[CH3+_dot]/d[H3O+]
    pd(72,83) =  &
        +k(369)*n(idx_CH2)

    !d[H3CO+_dot]/d[H3O+]
    pd(82,83) =  &
        +k(527)*n(idx_H2CO)

    !d[H3O+_dot]/d[H3O+]
    pd(83,83) =  &
        -k(290)*n(idx_E)  &
        -k(529)*n(idx_HNC)  &
        -k(528)*n(idx_HCN)  &
        -k(369)*n(idx_CH2)  &
        -k(1106)  &
        -k(678)*n(idx_NH2)  &
        -k(330)*n(idx_C)  &
        -k(527)*n(idx_H2CO)  &
        -k(289)*n(idx_E)  &
        -k(287)*n(idx_E)  &
        -k(404)*n(idx_CH)  &
        -k(288)*n(idx_E)

    !d[HCNH+_dot]/d[H3O+]
    pd(84,83) =  &
        +k(528)*n(idx_HCN)  &
        +k(529)*n(idx_HNC)

    !d[E_dot]/d[HCNH+]
    pd(1,84) =  &
        -k(294)*n(idx_E)  &
        -k(292)*n(idx_E)  &
        -k(293)*n(idx_E)

    !d[CH_dot]/d[HCNH+]
    pd(2,84) =  &
        -k(406)*n(idx_CH)  &
        -k(407)*n(idx_CH)

    !d[HNC_dot]/d[HCNH+]
    pd(4,84) =  &
        +k(549)*n(idx_H2CO)  &
        +k(681)*n(idx_NH2)  &
        +k(294)*n(idx_E)  &
        +k(372)*n(idx_CH2)  &
        +k(407)*n(idx_CH)

    !d[HCN_dot]/d[HCNH+]
    pd(5,84) =  &
        +k(293)*n(idx_E)  &
        +k(371)*n(idx_CH2)  &
        +k(680)*n(idx_NH2)  &
        +k(548)*n(idx_H2CO)  &
        +k(406)*n(idx_CH)

    !d[H_dot]/d[HCNH+]
    pd(8,84) =  &
        +2.d0*k(292)*n(idx_E)  &
        +k(293)*n(idx_E)  &
        +k(294)*n(idx_E)  &
        +k(1121)

    !d[CH2_dot]/d[HCNH+]
    pd(12,84) =  &
        -k(372)*n(idx_CH2)  &
        -k(371)*n(idx_CH2)

    !d[H2CO_dot]/d[HCNH+]
    pd(13,84) =  &
        -k(548)*n(idx_H2CO)  &
        -k(549)*n(idx_H2CO)

    !d[CN_dot]/d[HCNH+]
    pd(18,84) =  &
        +k(292)*n(idx_E)

    !d[NH2_dot]/d[HCNH+]
    pd(21,84) =  &
        -k(680)*n(idx_NH2)  &
        -k(681)*n(idx_NH2)

    !d[HCN_DUST_dot]/d[HCNH+]
    pd(44,84) =  &
        +k(1121)

    !d[CH2+_dot]/d[HCNH+]
    pd(58,84) =  &
        +k(407)*n(idx_CH)  &
        +k(406)*n(idx_CH)

    !d[NH3+_dot]/d[HCNH+]
    pd(62,84) =  &
        +k(681)*n(idx_NH2)  &
        +k(680)*n(idx_NH2)

    !d[CH3+_dot]/d[HCNH+]
    pd(72,84) =  &
        +k(371)*n(idx_CH2)  &
        +k(372)*n(idx_CH2)

    !d[H3CO+_dot]/d[HCNH+]
    pd(82,84) =  &
        +k(548)*n(idx_H2CO)  &
        +k(549)*n(idx_H2CO)

    !d[HCNH+_dot]/d[HCNH+]
    pd(84,84) =  &
        -k(372)*n(idx_CH2)  &
        -k(294)*n(idx_E)  &
        -k(293)*n(idx_E)  &
        -k(1121)  &
        -k(407)*n(idx_CH)  &
        -k(371)*n(idx_CH2)  &
        -k(406)*n(idx_CH)  &
        -k(548)*n(idx_H2CO)  &
        -k(681)*n(idx_NH2)  &
        -k(549)*n(idx_H2CO)  &
        -k(680)*n(idx_NH2)  &
        -k(292)*n(idx_E)

    !d[E_dot]/d[HCO2+]
    pd(1,85) =  &
        -k(298)*n(idx_E)  &
        -k(296)*n(idx_E)  &
        -k(297)*n(idx_E)

    !d[O_dot]/d[HCO2+]
    pd(3,85) =  &
        -k(719)*n(idx_O)  &
        +k(297)*n(idx_E)

    !d[C_dot]/d[HCO2+]
    pd(7,85) =  &
        -k(333)*n(idx_C)

    !d[H_dot]/d[HCO2+]
    pd(8,85) =  &
        +k(1107)  &
        +k(296)*n(idx_E)  &
        +k(297)*n(idx_E)

    !d[H2O_dot]/d[HCO2+]
    pd(9,85) =  &
        -k(497)*n(idx_H2O)

    !d[OH_dot]/d[HCO2+]
    pd(10,85) =  &
        +k(298)*n(idx_E)

    !d[O2_dot]/d[HCO2+]
    pd(11,85) =  &
        +k(719)*n(idx_O)

    !d[CO_dot]/d[HCO2+]
    pd(19,85) =  &
        -k(424)*n(idx_CO)  &
        +k(297)*n(idx_E)  &
        +k(298)*n(idx_E)

    !d[CO2_dot]/d[HCO2+]
    pd(29,85) =  &
        +k(333)*n(idx_C)  &
        +k(497)*n(idx_H2O)  &
        +k(424)*n(idx_CO)  &
        +k(296)*n(idx_E)

    !d[CO2_DUST_dot]/d[HCO2+]
    pd(42,85) =  &
        +k(1107)

    !d[HCO+_dot]/d[HCO2+]
    pd(54,85) =  &
        +k(424)*n(idx_CO)  &
        +k(719)*n(idx_O)

    !d[CH+_dot]/d[HCO2+]
    pd(59,85) =  &
        +k(333)*n(idx_C)

    !d[H3O+_dot]/d[HCO2+]
    pd(83,85) =  &
        +k(497)*n(idx_H2O)

    !d[HCO2+_dot]/d[HCO2+]
    pd(85,85) =  &
        -k(296)*n(idx_E)  &
        -k(298)*n(idx_E)  &
        -k(424)*n(idx_CO)  &
        -k(1107)  &
        -k(297)*n(idx_E)  &
        -k(719)*n(idx_O)  &
        -k(333)*n(idx_C)  &
        -k(497)*n(idx_H2O)

    !d[E_dot]/d[HEH+]
    pd(1,86) =  &
        -k(301)*n(idx_E)

    !d[H2_dot]/d[HEH+]
    pd(6,86) =  &
        -k(469)*n(idx_H2)

    !d[H_dot]/d[HEH+]
    pd(8,86) =  &
        +k(301)*n(idx_E)  &
        -k(534)*n(idx_H)

    !d[HE_dot]/d[HEH+]
    pd(26,86) =  &
        +k(534)*n(idx_H)  &
        +k(301)*n(idx_E)  &
        +k(469)*n(idx_H2)

    !d[H2+_dot]/d[HEH+]
    pd(77,86) =  &
        +k(534)*n(idx_H)

    !d[H3+_dot]/d[HEH+]
    pd(81,86) =  &
        +k(469)*n(idx_H2)

    !d[HEH+_dot]/d[HEH+]
    pd(86,86) =  &
        -k(534)*n(idx_H)  &
        -k(301)*n(idx_E)  &
        -k(469)*n(idx_H2)

    !d[E_dot]/d[N2H+]
    pd(1,87) =  &
        -k(303)*n(idx_E)  &
        -k(304)*n(idx_E)

    !d[CH_dot]/d[N2H+]
    pd(2,87) =  &
        -k(411)*n(idx_CH)

    !d[O_dot]/d[N2H+]
    pd(3,87) =  &
        -k(721)*n(idx_O)

    !d[HNC_dot]/d[N2H+]
    pd(4,87) =  &
        -k(561)*n(idx_HNC)

    !d[HCN_dot]/d[N2H+]
    pd(5,87) =  &
        -k(546)*n(idx_HCN)

    !d[C_dot]/d[N2H+]
    pd(7,87) =  &
        -k(335)*n(idx_C)

    !d[H_dot]/d[N2H+]
    pd(8,87) =  &
        +k(1122)  &
        +k(303)*n(idx_E)

    !d[H2O_dot]/d[N2H+]
    pd(9,87) =  &
        -k(500)*n(idx_H2O)

    !d[OH_dot]/d[N2H+]
    pd(10,87) =  &
        -k(744)*n(idx_OH)

    !d[CH2_dot]/d[N2H+]
    pd(12,87) =  &
        -k(375)*n(idx_CH2)

    !d[H2CO_dot]/d[N2H+]
    pd(13,87) =  &
        -k(633)*n(idx_H2CO)

    !d[HCO_dot]/d[N2H+]
    pd(14,87) =  &
        -k(554)*n(idx_HCO)

    !d[CO_dot]/d[N2H+]
    pd(19,87) =  &
        -k(426)*n(idx_CO)

    !d[N2_dot]/d[N2H+]
    pd(20,87) =  &
        +k(684)*n(idx_NH2)  &
        +k(633)*n(idx_H2CO)  &
        +k(335)*n(idx_C)  &
        +k(561)*n(idx_HNC)  &
        +k(375)*n(idx_CH2)  &
        +k(411)*n(idx_CH)  &
        +k(546)*n(idx_HCN)  &
        +k(500)*n(idx_H2O)  &
        +k(721)*n(idx_O)  &
        +k(744)*n(idx_OH)  &
        +k(632)*n(idx_CO2)  &
        +k(696)*n(idx_NH)  &
        +k(426)*n(idx_CO)  &
        +k(554)*n(idx_HCO)  &
        +k(303)*n(idx_E)

    !d[NH2_dot]/d[N2H+]
    pd(21,87) =  &
        -k(684)*n(idx_NH2)

    !d[N_dot]/d[N2H+]
    pd(24,87) =  &
        +k(304)*n(idx_E)

    !d[NH_dot]/d[N2H+]
    pd(25,87) =  &
        +k(304)*n(idx_E)  &
        -k(696)*n(idx_NH)

    !d[CO2_dot]/d[N2H+]
    pd(29,87) =  &
        -k(632)*n(idx_CO2)

    !d[N2_DUST_dot]/d[N2H+]
    pd(43,87) =  &
        +k(1122)

    !d[HCO+_dot]/d[N2H+]
    pd(54,87) =  &
        +k(426)*n(idx_CO)

    !d[CH2+_dot]/d[N2H+]
    pd(58,87) =  &
        +k(411)*n(idx_CH)

    !d[CH+_dot]/d[N2H+]
    pd(59,87) =  &
        +k(335)*n(idx_C)

    !d[H2CO+_dot]/d[N2H+]
    pd(60,87) =  &
        +k(554)*n(idx_HCO)

    !d[NH3+_dot]/d[N2H+]
    pd(62,87) =  &
        +k(684)*n(idx_NH2)

    !d[H2O+_dot]/d[N2H+]
    pd(68,87) =  &
        +k(744)*n(idx_OH)

    !d[NH2+_dot]/d[N2H+]
    pd(69,87) =  &
        +k(696)*n(idx_NH)

    !d[OH+_dot]/d[N2H+]
    pd(71,87) =  &
        +k(721)*n(idx_O)

    !d[CH3+_dot]/d[N2H+]
    pd(72,87) =  &
        +k(375)*n(idx_CH2)

    !d[H3CO+_dot]/d[N2H+]
    pd(82,87) =  &
        +k(633)*n(idx_H2CO)

    !d[H3O+_dot]/d[N2H+]
    pd(83,87) =  &
        +k(500)*n(idx_H2O)

    !d[HCNH+_dot]/d[N2H+]
    pd(84,87) =  &
        +k(546)*n(idx_HCN)  &
        +k(561)*n(idx_HNC)

    !d[HCO2+_dot]/d[N2H+]
    pd(85,87) =  &
        +k(632)*n(idx_CO2)

    !d[N2H+_dot]/d[N2H+]
    pd(87,87) =  &
        -k(375)*n(idx_CH2)  &
        -k(554)*n(idx_HCO)  &
        -k(335)*n(idx_C)  &
        -k(304)*n(idx_E)  &
        -k(426)*n(idx_CO)  &
        -k(561)*n(idx_HNC)  &
        -k(633)*n(idx_H2CO)  &
        -k(411)*n(idx_CH)  &
        -k(303)*n(idx_E)  &
        -k(546)*n(idx_HCN)  &
        -k(721)*n(idx_O)  &
        -k(696)*n(idx_NH)  &
        -k(500)*n(idx_H2O)  &
        -k(632)*n(idx_CO2)  &
        -k(1122)  &
        -k(684)*n(idx_NH2)  &
        -k(744)*n(idx_OH)

    !d[E_dot]/d[O2H+]
    pd(1,88) =  &
        -k(312)*n(idx_E)

    !d[CH_dot]/d[O2H+]
    pd(2,88) =  &
        -k(416)*n(idx_CH)

    !d[O_dot]/d[O2H+]
    pd(3,88) =  &
        -k(724)*n(idx_O)

    !d[HNC_dot]/d[O2H+]
    pd(4,88) =  &
        -k(562)*n(idx_HNC)

    !d[HCN_dot]/d[O2H+]
    pd(5,88) =  &
        -k(547)*n(idx_HCN)

    !d[H2_dot]/d[O2H+]
    pd(6,88) =  &
        -k(476)*n(idx_H2)

    !d[C_dot]/d[O2H+]
    pd(7,88) =  &
        -k(338)*n(idx_C)

    !d[H_dot]/d[O2H+]
    pd(8,88) =  &
        +k(312)*n(idx_E)

    !d[H2O_dot]/d[O2H+]
    pd(9,88) =  &
        -k(501)*n(idx_H2O)

    !d[OH_dot]/d[O2H+]
    pd(10,88) =  &
        -k(745)*n(idx_OH)

    !d[O2_dot]/d[O2H+]
    pd(11,88) =  &
        +k(556)*n(idx_HCO)  &
        +k(501)*n(idx_H2O)  &
        +k(702)*n(idx_NO)  &
        +k(685)*n(idx_NH2)  &
        +k(422)*n(idx_CN)  &
        +k(482)*n(idx_H2CO)  &
        +k(476)*n(idx_H2)  &
        +k(338)*n(idx_C)  &
        +k(700)*n(idx_NH)  &
        +k(631)*n(idx_N2)  &
        +k(380)*n(idx_CH2)  &
        +k(562)*n(idx_HNC)  &
        +k(312)*n(idx_E)  &
        +k(745)*n(idx_OH)  &
        +k(716)*n(idx_CO2)  &
        +k(547)*n(idx_HCN)  &
        +k(724)*n(idx_O)  &
        +k(427)*n(idx_CO)  &
        +k(416)*n(idx_CH)

    !d[CH2_dot]/d[O2H+]
    pd(12,88) =  &
        -k(380)*n(idx_CH2)

    !d[H2CO_dot]/d[O2H+]
    pd(13,88) =  &
        -k(482)*n(idx_H2CO)

    !d[HCO_dot]/d[O2H+]
    pd(14,88) =  &
        -k(556)*n(idx_HCO)

    !d[NO_dot]/d[O2H+]
    pd(17,88) =  &
        -k(702)*n(idx_NO)

    !d[CN_dot]/d[O2H+]
    pd(18,88) =  &
        -k(422)*n(idx_CN)

    !d[CO_dot]/d[O2H+]
    pd(19,88) =  &
        -k(427)*n(idx_CO)

    !d[N2_dot]/d[O2H+]
    pd(20,88) =  &
        -k(631)*n(idx_N2)

    !d[NH2_dot]/d[O2H+]
    pd(21,88) =  &
        -k(685)*n(idx_NH2)

    !d[NH_dot]/d[O2H+]
    pd(25,88) =  &
        -k(700)*n(idx_NH)

    !d[CO2_dot]/d[O2H+]
    pd(29,88) =  &
        -k(716)*n(idx_CO2)

    !d[O2H_DUST_dot]/d[O2H+]
    pd(49,88) =  &
        +k(1117)

    !d[HCO+_dot]/d[O2H+]
    pd(54,88) =  &
        +k(427)*n(idx_CO)

    !d[CH2+_dot]/d[O2H+]
    pd(58,88) =  &
        +k(416)*n(idx_CH)

    !d[CH+_dot]/d[O2H+]
    pd(59,88) =  &
        +k(338)*n(idx_C)

    !d[H2CO+_dot]/d[O2H+]
    pd(60,88) =  &
        +k(556)*n(idx_HCO)

    !d[NH3+_dot]/d[O2H+]
    pd(62,88) =  &
        +k(685)*n(idx_NH2)

    !d[H2O+_dot]/d[O2H+]
    pd(68,88) =  &
        +k(745)*n(idx_OH)

    !d[NH2+_dot]/d[O2H+]
    pd(69,88) =  &
        +k(700)*n(idx_NH)

    !d[OH+_dot]/d[O2H+]
    pd(71,88) =  &
        +k(724)*n(idx_O)

    !d[CH3+_dot]/d[O2H+]
    pd(72,88) =  &
        +k(380)*n(idx_CH2)

    !d[HCN+_dot]/d[O2H+]
    pd(75,88) =  &
        +k(422)*n(idx_CN)

    !d[HNO+_dot]/d[O2H+]
    pd(79,88) =  &
        +k(702)*n(idx_NO)

    !d[H3+_dot]/d[O2H+]
    pd(81,88) =  &
        +k(476)*n(idx_H2)

    !d[H3CO+_dot]/d[O2H+]
    pd(82,88) =  &
        +k(482)*n(idx_H2CO)

    !d[H3O+_dot]/d[O2H+]
    pd(83,88) =  &
        +k(501)*n(idx_H2O)

    !d[HCNH+_dot]/d[O2H+]
    pd(84,88) =  &
        +k(562)*n(idx_HNC)  &
        +k(547)*n(idx_HCN)

    !d[HCO2+_dot]/d[O2H+]
    pd(85,88) =  &
        +k(716)*n(idx_CO2)

    !d[N2H+_dot]/d[O2H+]
    pd(87,88) =  &
        +k(631)*n(idx_N2)

    !d[O2H+_dot]/d[O2H+]
    pd(88,88) =  &
        -k(422)*n(idx_CN)  &
        -k(416)*n(idx_CH)  &
        -k(685)*n(idx_NH2)  &
        -k(702)*n(idx_NO)  &
        -k(700)*n(idx_NH)  &
        -k(501)*n(idx_H2O)  &
        -k(380)*n(idx_CH2)  &
        -k(476)*n(idx_H2)  &
        -k(745)*n(idx_OH)  &
        -k(562)*n(idx_HNC)  &
        -k(312)*n(idx_E)  &
        -k(338)*n(idx_C)  &
        -k(724)*n(idx_O)  &
        -k(427)*n(idx_CO)  &
        -k(631)*n(idx_N2)  &
        -k(1117)  &
        -k(547)*n(idx_HCN)  &
        -k(716)*n(idx_CO2)  &
        -k(556)*n(idx_HCO)  &
        -k(482)*n(idx_H2CO)

  end subroutine jex

end module krome_ode
