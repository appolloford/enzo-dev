!This module contains useful routines to get physical
! quantities, like mean molecular weight, mass density,
! mass, jeans length, etc. etc.

!############### MODULE ##############
module krome_getphys
contains

  !*****************************
  !get the mean molecular weight
  function get_mu(n)
    use krome_commons
    use krome_constants
    implicit none
    real*8::n(:),get_mu,m(nspec)
    m(:) = get_mass()

    !ip_mass is 1/proton_mass_in_g
    get_mu = max(sum(n(1:nmols)*m(1:nmols)),1d-40) &
        / max(sum(n(1:nmols)),1d-40) * ip_mass

  end function get_mu

  !***************************
  !get mean molecular weight
  function get_mu_rho(n,rhogas)
    use krome_commons
    use krome_constants
    implicit none
    real*8::get_mu_rho,rhogas,n(:)

    !ip_mass is 1/proton_mass_in_g
    get_mu_rho = rhogas / max(sum(n(1:nmols)),1d-40) * ip_mass

  end function get_mu_rho

  !************************
  !get species masses (g)
  function get_mass()
    use krome_commons
    implicit none
    real*8::get_mass(nspec)

    get_mass(1) = 9.10938188d-28	!E
    get_mass(2) = 2.17497276273d-23	!CH
    get_mass(3) = 2.67682601455d-23	!O
    get_mass(4) = 4.51719552546d-23	!HNC
    get_mass(5) = 4.51719552546d-23	!HCN
    get_mass(6) = 3.34706503638d-24	!H2
    get_mass(7) = 2.00761951091d-23	!C
    get_mass(8) = 1.67353251819d-24	!H
    get_mass(9) = 3.01153251819d-23	!H2O
    get_mass(10) = 2.84417926637d-23	!OH
    get_mass(11) = 5.3536520291d-23	!O2
    get_mass(12) = 2.34232601455d-23	!CH2
    get_mass(13) = 5.0191520291d-23	!H2CO
    get_mass(14) = 4.85179877728d-23	!HCO
    get_mass(15) = 4.01523902183d-23	!MG
    get_mass(16) = 2.84428251819d-23	!NH3
    get_mass(17) = 5.01904877728d-23	!NO
    get_mass(18) = 4.34984227364d-23	!CN
    get_mass(19) = 4.68444552546d-23	!CO
    get_mass(20) = 4.68444552546d-23	!N2
    get_mass(21) = 2.67692926637d-23	!NH2
    get_mass(22) = 2.50967926637d-23	!CH3
    get_mass(23) = 2.67703251819d-23	!CH4
    get_mass(24) = 2.34222276273d-23	!N
    get_mass(25) = 2.50957601455d-23	!NH
    get_mass(26) = 6.69206503638d-24	!HE
    get_mass(27) = 5.1864020291d-23	!HNO
    get_mass(28) = 5.35385853274d-23	!CH3OH
    get_mass(29) = 7.36127154001d-23	!CO2
    get_mass(30) = 4.68454877728d-23	!H2CN
    get_mass(31) = 7.19402154001d-23	!HNCO
    get_mass(32) = 7.69587479183d-23	!NO2
    get_mass(33) = 5.52100528092d-23	!O2H
    get_mass(34) = 7.02666828819d-23	!OCN
    get_mass(35) = 5.35385853274d-23	!CH3OH_DUST
    get_mass(36) = 7.19402154001d-23	!HNCO_DUST
    get_mass(37) = 5.0191520291d-23	!H2CO_DUST
    get_mass(38) = 2.67703251819d-23	!CH4_DUST
    get_mass(39) = 4.68444552546d-23	!CO_DUST
    get_mass(40) = 3.01153251819d-23	!H2O_DUST
    get_mass(41) = 5.01904877728d-23	!NO_DUST
    get_mass(42) = 7.36127154001d-23	!CO2_DUST
    get_mass(43) = 4.68444552546d-23	!N2_DUST
    get_mass(44) = 4.51719552546d-23	!HCN_DUST
    get_mass(45) = 2.84428251819d-23	!NH3_DUST
    get_mass(46) = 5.3536520291d-23	!O2_DUST
    get_mass(47) = 7.69587479183d-23	!NO2_DUST
    get_mass(48) = 5.1864020291d-23	!HNO_DUST
    get_mass(49) = 5.52100528092d-23	!O2H_DUST
    get_mass(50) = 4.68454877728d-23	!H2CN_DUST
    get_mass(51) = 4.01523902183d-23	!MG_DUST
    get_mass(52) = 4.51719552546d-23	!HNC_DUST
    get_mass(53) = 9.10938188d-28	!E_DUST
    get_mass(54) = 4.85170768346d-23	!HCO+
    get_mass(55) = 1.67262158d-24	!H+
    get_mass(56) = 4.85170768346d-23	!HOC+
    get_mass(57) = 2.00752841709d-23	!C+
    get_mass(58) = 2.34223492073d-23	!CH2+
    get_mass(59) = 2.17488166891d-23	!CH+
    get_mass(60) = 5.01906093528d-23	!H2CO+
    get_mass(61) = 4.01514792801d-23	!MG+
    get_mass(62) = 2.84419142437d-23	!NH3+
    get_mass(63) = 5.01895768346d-23	!NO+
    get_mass(64) = 4.34975117983d-23	!CN+
    get_mass(65) = 4.68435443164d-23	!CO+
    get_mass(66) = 4.68435443164d-23	!N2+
    get_mass(67) = 5.35356093528d-23	!O2+
    get_mass(68) = 3.01144142437d-23	!H2O+
    get_mass(69) = 2.67683817255d-23	!NH2+
    get_mass(70) = 2.67673492073d-23	!O+
    get_mass(71) = 2.84408817255d-23	!OH+
    get_mass(72) = 2.50958817255d-23	!CH3+
    get_mass(73) = 2.67694142437d-23	!CH4+
    get_mass(74) = 2.34213166891d-23	!N+
    get_mass(75) = 4.51710443164d-23	!HCN+
    get_mass(76) = 2.50948492073d-23	!NH+
    get_mass(77) = 3.34615409819d-24	!H2+
    get_mass(78) = 6.69115409819d-24	!HE+
    get_mass(79) = 5.18631093528d-23	!HNO+
    get_mass(80) = 5.3536641871d-23	!H2NO+
    get_mass(81) = 5.01968661638d-24	!H3+
    get_mass(82) = 5.1864141871d-23	!H3CO+
    get_mass(83) = 3.17879467619d-23	!H3O+
    get_mass(84) = 4.68445768346d-23	!HCNH+
    get_mass(85) = 7.52853369801d-23	!HCO2+
    get_mass(86) = 8.36468661638d-24	!HEH+
    get_mass(87) = 4.85170768346d-23	!N2H+
    get_mass(88) = 5.5209141871d-23	!O2H+
    get_mass(89) = 0.d0	!CR
    get_mass(90) = 0.d0	!g
    get_mass(91) = 0.d0	!Tgas
    get_mass(92) = 0.d0	!dummy

  end function get_mass

  !************************
  !get sqrt of the inverse of the masses (1/sqrt(g))
  function get_imass_sqrt()
    use krome_commons
    implicit none
    real*8::get_imass_sqrt(nspec)

    get_imass_sqrt(1) = 3.31326021505d+13	!E
    get_imass_sqrt(2) = 2.14423849574d+11	!CH
    get_imass_sqrt(3) = 1.93281339991d+11	!O
    get_imass_sqrt(4) = 1.48787194664d+11	!HNC
    get_imass_sqrt(5) = 1.48787194664d+11	!HCN
    get_imass_sqrt(6) = 5.46597856701d+11	!H2
    get_imass_sqrt(7) = 2.23182067346d+11	!C
    get_imass_sqrt(8) = 7.73006102111d+11	!H
    get_imass_sqrt(9) = 1.82224271009d+11	!H2O
    get_imass_sqrt(10) = 1.87508740611d+11	!OH
    get_imass_sqrt(11) = 1.36670546184d+11	!O2
    get_imass_sqrt(12) = 2.06621889668d+11	!CH2
    get_imass_sqrt(13) = 1.4115128127d+11	!H2CO
    get_imass_sqrt(14) = 1.43565011358d+11	!HCO
    get_imass_sqrt(15) = 1.57813553259d+11	!MG
    get_imass_sqrt(16) = 1.87505337153d+11	!NH3
    get_imass_sqrt(17) = 1.41152733144d+11	!NO
    get_imass_sqrt(18) = 1.51622357573d+11	!CN
    get_imass_sqrt(19) = 1.46106959624d+11	!CO
    get_imass_sqrt(20) = 1.46106959624d+11	!N2
    get_imass_sqrt(21) = 1.93277612428d+11	!NH2
    get_imass_sqrt(22) = 1.99613949989d+11	!CH3
    get_imass_sqrt(23) = 1.93273885081d+11	!CH4
    get_imass_sqrt(24) = 2.06626443857d+11	!N
    get_imass_sqrt(25) = 1.99618056318d+11	!NH
    get_imass_sqrt(26) = 3.86562679981d+11	!HE
    get_imass_sqrt(27) = 1.38856722679d+11	!HNO
    get_imass_sqrt(28) = 1.36667910399d+11	!CH3OH
    get_imass_sqrt(29) = 1.16553033405d+11	!CO2
    get_imass_sqrt(30) = 1.46105349448d+11	!H2CN
    get_imass_sqrt(31) = 1.17900089041d+11	!HNCO
    get_imass_sqrt(32) = 1.13991115426d+11	!NO2
    get_imass_sqrt(33) = 1.34583221186d+11	!O2H
    get_imass_sqrt(34) = 1.19295832983d+11	!OCN
    get_imass_sqrt(35) = 1.36667910399d+11	!CH3OH_DUST
    get_imass_sqrt(36) = 1.17900089041d+11	!HNCO_DUST
    get_imass_sqrt(37) = 1.4115128127d+11	!H2CO_DUST
    get_imass_sqrt(38) = 1.93273885081d+11	!CH4_DUST
    get_imass_sqrt(39) = 1.46106959624d+11	!CO_DUST
    get_imass_sqrt(40) = 1.82224271009d+11	!H2O_DUST
    get_imass_sqrt(41) = 1.41152733144d+11	!NO_DUST
    get_imass_sqrt(42) = 1.16553033405d+11	!CO2_DUST
    get_imass_sqrt(43) = 1.46106959624d+11	!N2_DUST
    get_imass_sqrt(44) = 1.48787194664d+11	!HCN_DUST
    get_imass_sqrt(45) = 1.87505337153d+11	!NH3_DUST
    get_imass_sqrt(46) = 1.36670546184d+11	!O2_DUST
    get_imass_sqrt(47) = 1.13991115426d+11	!NO2_DUST
    get_imass_sqrt(48) = 1.38856722679d+11	!HNO_DUST
    get_imass_sqrt(49) = 1.34583221186d+11	!O2H_DUST
    get_imass_sqrt(50) = 1.46105349448d+11	!H2CN_DUST
    get_imass_sqrt(51) = 1.57813553259d+11	!MG_DUST
    get_imass_sqrt(52) = 1.48787194664d+11	!HNC_DUST
    get_imass_sqrt(53) = 3.31326021505d+13	!E_DUST
    get_imass_sqrt(54) = 1.43566359113d+11	!HCO+
    get_imass_sqrt(55) = 7.732165696d+11	!H+
    get_imass_sqrt(56) = 1.43566359113d+11	!HOC+
    get_imass_sqrt(57) = 2.23187130855d+11	!C+
    get_imass_sqrt(58) = 2.06625907582d+11	!CH2+
    get_imass_sqrt(59) = 2.14428340044d+11	!CH+
    get_imass_sqrt(60) = 1.41152562182d+11	!H2CO+
    get_imass_sqrt(61) = 1.5781534345d+11	!MG+
    get_imass_sqrt(62) = 1.87508339841d+11	!NH3+
    get_imass_sqrt(63) = 1.41154014095d+11	!NO+
    get_imass_sqrt(64) = 1.51623945226d+11	!CN+
    get_imass_sqrt(65) = 1.46108380244d+11	!CO+
    get_imass_sqrt(66) = 1.46108380244d+11	!N2+
    get_imass_sqrt(67) = 1.36671708942d+11	!O2+
    get_imass_sqrt(68) = 1.82227027061d+11	!H2O+
    get_imass_sqrt(69) = 1.93280901055d+11	!NH2+
    get_imass_sqrt(70) = 1.93284628808d+11	!O+
    get_imass_sqrt(71) = 1.87511743463d+11	!OH+
    get_imass_sqrt(72) = 1.99617572781d+11	!CH3+
    get_imass_sqrt(73) = 1.93277173518d+11	!CH4+
    get_imass_sqrt(74) = 2.06630462037d+11	!N+
    get_imass_sqrt(75) = 1.48788694909d+11	!HCN+
    get_imass_sqrt(76) = 1.99621679333d+11	!NH+
    get_imass_sqrt(77) = 5.46672253003d+11	!H2+
    get_imass_sqrt(78) = 3.86588992536d+11	!HE+
    get_imass_sqrt(79) = 1.38857942133d+11	!HNO+
    get_imass_sqrt(80) = 1.36670390997d+11	!H2NO+
    get_imass_sqrt(81) = 4.463357746d+11	!H3+
    get_imass_sqrt(82) = 1.38856559925d+11	!H3CO+
    get_imass_sqrt(83) = 1.77365342346d+11	!H3O+
    get_imass_sqrt(84) = 1.46106770022d+11	!HCNH+
    get_imass_sqrt(85) = 1.15251026098d+11	!HCO2+
    get_imass_sqrt(86) = 3.45760328884d+11	!HEH+
    get_imass_sqrt(87) = 1.43566359113d+11	!N2H+
    get_imass_sqrt(88) = 1.34584331478d+11	!O2H+
    get_imass_sqrt(89) = 0.d0	!CR
    get_imass_sqrt(90) = 0.d0	!g
    get_imass_sqrt(91) = 0.d0	!Tgas
    get_imass_sqrt(92) = 0.d0	!dummy

  end function get_imass_sqrt

  !************************
  !get inverse of the species masses (1/g)
  function get_imass()
    use krome_commons
    implicit none
    real*8::get_imass(nspec)

    get_imass(1) = 1.09776932527d+27	!E
    get_imass(2) = 4.59775872662d+22	!CH
    get_imass(3) = 3.73576763885d+22	!O
    get_imass(4) = 2.21376292959d+22	!HNC
    get_imass(5) = 2.21376292959d+22	!HCN
    get_imass(6) = 2.9876921695d+23	!H2
    get_imass(7) = 4.98102351847d+22	!C
    get_imass(8) = 5.97538433901d+23	!H
    get_imass(9) = 3.32056849448d+22	!H2O
    get_imass(10) = 3.51595278056d+22	!OH
    get_imass(11) = 1.86788381943d+22	!O2
    get_imass(12) = 4.26926052901d+22	!CH2
    get_imass(13) = 1.99236842041d+22	!H2CO
    get_imass(14) = 2.06109124864d+22	!HCO
    get_imass(15) = 2.49051175924d+22	!MG
    get_imass(16) = 3.51582514608d+22	!NH3
    get_imass(17) = 1.99240940739d+22	!NO
    get_imass(18) = 2.2989339316d+22	!CN
    get_imass(19) = 2.13472436506d+22	!CO
    get_imass(20) = 2.13472436506d+22	!N2
    get_imass(21) = 3.73562354659d+22	!NH2
    get_imass(22) = 3.984572903d+22	!CH3
    get_imass(23) = 3.73547946544d+22	!CH4
    get_imass(24) = 4.26944873012d+22	!N
    get_imass(25) = 3.98473684081d+22	!NH
    get_imass(26) = 1.49430705554d+23	!HE
    get_imass(27) = 1.92811894332d+22	!HNO
    get_imass(28) = 1.86781177329d+22	!CH3OH
    get_imass(29) = 1.35846095958d+22	!CO2
    get_imass(30) = 2.13467731375d+22	!H2CN
    get_imass(31) = 1.39004309959d+22	!HNCO
    get_imass(32) = 1.2993974396d+22	!NO2
    get_imass(33) = 1.81126434248d+22	!O2H
    get_imass(34) = 1.42314957671d+22	!OCN
    get_imass(35) = 1.86781177329d+22	!CH3OH_DUST
    get_imass(36) = 1.39004309959d+22	!HNCO_DUST
    get_imass(37) = 1.99236842041d+22	!H2CO_DUST
    get_imass(38) = 3.73547946544d+22	!CH4_DUST
    get_imass(39) = 2.13472436506d+22	!CO_DUST
    get_imass(40) = 3.32056849448d+22	!H2O_DUST
    get_imass(41) = 1.99240940739d+22	!NO_DUST
    get_imass(42) = 1.35846095958d+22	!CO2_DUST
    get_imass(43) = 2.13472436506d+22	!N2_DUST
    get_imass(44) = 2.21376292959d+22	!HCN_DUST
    get_imass(45) = 3.51582514608d+22	!NH3_DUST
    get_imass(46) = 1.86788381943d+22	!O2_DUST
    get_imass(47) = 1.2993974396d+22	!NO2_DUST
    get_imass(48) = 1.92811894332d+22	!HNO_DUST
    get_imass(49) = 1.81126434248d+22	!O2H_DUST
    get_imass(50) = 2.13467731375d+22	!H2CN_DUST
    get_imass(51) = 2.49051175924d+22	!MG_DUST
    get_imass(52) = 2.21376292959d+22	!HNC_DUST
    get_imass(53) = 1.09776932527d+27	!E_DUST
    get_imass(54) = 2.0611299469d+22	!HCO+
    get_imass(55) = 5.97863863505d+23	!H+
    get_imass(56) = 2.0611299469d+22	!HOC+
    get_imass(57) = 4.98124953791d+22	!C+
    get_imass(58) = 4.2694265684d+22	!CH2+
    get_imass(59) = 4.59795130141d+22	!CH+
    get_imass(60) = 1.99240458105d+22	!H2CO+
    get_imass(61) = 2.49056826281d+22	!MG+
    get_imass(62) = 3.515937751d+22	!NH3+
    get_imass(63) = 1.99244556952d+22	!NO+
    get_imass(64) = 2.29898207658d+22	!CN+
    get_imass(65) = 2.13476587776d+22	!CO+
    get_imass(66) = 2.13476587776d+22	!N2+
    get_imass(67) = 1.86791560251d+22	!O2+
    get_imass(68) = 3.32066893916d+22	!H2O+
    get_imass(69) = 3.73575067128d+22	!NH2+
    get_imass(70) = 3.73589477335d+22	!O+
    get_imass(71) = 3.51606539365d+22	!OH+
    get_imass(72) = 3.98471753628d+22	!CH3+
    get_imass(73) = 3.73560658032d+22	!CH4+
    get_imass(74) = 4.26961478414d+22	!N+
    get_imass(75) = 2.21380757326d+22	!HCN+
    get_imass(76) = 3.98488148599d+22	!NH+
    get_imass(77) = 2.98850552203d+23	!H2+
    get_imass(78) = 1.4945104915d+23	!HE+
    get_imass(79) = 1.92815280934d+22	!HNO+
    get_imass(80) = 1.86787957752d+22	!H2NO+
    get_imass(81) = 1.99215623688d+23	!H3+
    get_imass(82) = 1.92811442342d+22	!H3CO+
    get_imass(83) = 3.14584646656d+22	!H3O+
    get_imass(84) = 2.13471882461d+22	!HCNH+
    get_imass(85) = 1.32827990165d+22	!HCO2+
    get_imass(86) = 1.1955020503d+23	!HEH+
    get_imass(87) = 2.0611299469d+22	!N2H+
    get_imass(88) = 1.81129422793d+22	!O2H+
    get_imass(89) = 0.d0	!CR
    get_imass(90) = 0.d0	!g
    get_imass(91) = 0.d0	!Tgas
    get_imass(92) = 0.d0	!dummy

  end function get_imass

  !************************
  !species binding energies (surface=BARE), K
  function get_EbindBare()
    use krome_commons
    implicit none
    real*8::get_EbindBare(nspec)

    get_EbindBare(:) = 1d99

    get_EbindBare(idx_O) = 1700.0d0
    get_EbindBare(idx_H2) = 300.0d0
    get_EbindBare(idx_H) = 500.0d0
    get_EbindBare(idx_H2O) = 4800.0d0
    get_EbindBare(idx_OH) = 1360.0d0
    get_EbindBare(idx_O2) = 1250.0d0
    get_EbindBare(idx_H2CO) = 1100.0d0
    get_EbindBare(idx_HCO) = 1100.0d0
    get_EbindBare(idx_CO) = 1100.0d0
    get_EbindBare(idx_CH3OH) = 1100.0d0
    get_EbindBare(idx_CO2) = 2300.0d0

  end function get_EbindBare

  !************************
  !species binding energies (surface=ICE), K
  function get_EbindIce()
    use krome_commons
    implicit none
    real*8::get_EbindIce(nspec)

    get_EbindIce(:) = 1d99

    get_EbindIce(idx_O) = 1700.0d0
    get_EbindIce(idx_H2) = 300.0d0
    get_EbindIce(idx_H) = 650.0d0
    get_EbindIce(idx_H2O) = 4800.0d0
    get_EbindIce(idx_OH) = 3500.0d0
    get_EbindIce(idx_O2) = 900.0d0
    get_EbindIce(idx_H2CO) = 3100.0d0
    get_EbindIce(idx_HCO) = 3100.0d0
    get_EbindIce(idx_CO) = 1300.0d0
    get_EbindIce(idx_CH3OH) = 3100.0d0
    get_EbindIce(idx_CO2) = 2300.0d0

  end function get_EbindIce

  !************************
  function get_kevap70()
    use krome_commons
    implicit none
    real*8::get_kevap70(nspec)

    get_kevap70(idx_E) = 0d0
    get_kevap70(idx_CH) = 0d0
    get_kevap70(idx_O) = 28.3692788833
    get_kevap70(idx_HNC) = 0d0
    get_kevap70(idx_HCN) = 0d0
    get_kevap70(idx_H2) = 13763786733.1
    get_kevap70(idx_C) = 0d0
    get_kevap70(idx_H) = 790490323.12
    get_kevap70(idx_H2O) = 1.65884938156e-18
    get_kevap70(idx_OH) = 3649.88043081
    get_kevap70(idx_O2) = 17568.7715065
    get_kevap70(idx_CH2) = 0d0
    get_kevap70(idx_H2CO) = 149751.929641
    get_kevap70(idx_HCO) = 149751.929641
    get_kevap70(idx_MG) = 0d0
    get_kevap70(idx_NH3) = 0d0
    get_kevap70(idx_NO) = 0d0
    get_kevap70(idx_CN) = 0d0
    get_kevap70(idx_CO) = 149751.929641
    get_kevap70(idx_N2) = 0d0
    get_kevap70(idx_NH2) = 0d0
    get_kevap70(idx_CH3) = 0d0
    get_kevap70(idx_CH4) = 0d0
    get_kevap70(idx_N) = 0d0
    get_kevap70(idx_NH) = 0d0
    get_kevap70(idx_HE) = 0d0
    get_kevap70(idx_HNO) = 0d0
    get_kevap70(idx_CH3OH) = 149751.929641
    get_kevap70(idx_CO2) = 0.00537432797219
    get_kevap70(idx_H2CN) = 0d0
    get_kevap70(idx_HNCO) = 0d0
    get_kevap70(idx_NO2) = 0d0
    get_kevap70(idx_O2H) = 0d0
    get_kevap70(idx_OCN) = 0d0
    get_kevap70(idx_CH3OH_DUST) = 0d0
    get_kevap70(idx_HNCO_DUST) = 0d0
    get_kevap70(idx_H2CO_DUST) = 0d0
    get_kevap70(idx_CH4_DUST) = 0d0
    get_kevap70(idx_CO_DUST) = 0d0
    get_kevap70(idx_H2O_DUST) = 0d0
    get_kevap70(idx_NO_DUST) = 0d0
    get_kevap70(idx_CO2_DUST) = 0d0
    get_kevap70(idx_N2_DUST) = 0d0
    get_kevap70(idx_HCN_DUST) = 0d0
    get_kevap70(idx_NH3_DUST) = 0d0
    get_kevap70(idx_O2_DUST) = 0d0
    get_kevap70(idx_NO2_DUST) = 0d0
    get_kevap70(idx_HNO_DUST) = 0d0
    get_kevap70(idx_O2H_DUST) = 0d0
    get_kevap70(idx_H2CN_DUST) = 0d0
    get_kevap70(idx_MG_DUST) = 0d0
    get_kevap70(idx_HNC_DUST) = 0d0
    get_kevap70(idx_E_DUST) = 0d0
    get_kevap70(idx_HCOj) = 0d0
    get_kevap70(idx_Hj) = 0d0
    get_kevap70(idx_HOCj) = 0d0
    get_kevap70(idx_Cj) = 0d0
    get_kevap70(idx_CH2j) = 0d0
    get_kevap70(idx_CHj) = 0d0
    get_kevap70(idx_H2COj) = 0d0
    get_kevap70(idx_MGj) = 0d0
    get_kevap70(idx_NH3j) = 0d0
    get_kevap70(idx_NOj) = 0d0
    get_kevap70(idx_CNj) = 0d0
    get_kevap70(idx_COj) = 0d0
    get_kevap70(idx_N2j) = 0d0
    get_kevap70(idx_O2j) = 0d0
    get_kevap70(idx_H2Oj) = 0d0
    get_kevap70(idx_NH2j) = 0d0
    get_kevap70(idx_Oj) = 0d0
    get_kevap70(idx_OHj) = 0d0
    get_kevap70(idx_CH3j) = 0d0
    get_kevap70(idx_CH4j) = 0d0
    get_kevap70(idx_Nj) = 0d0
    get_kevap70(idx_HCNj) = 0d0
    get_kevap70(idx_NHj) = 0d0
    get_kevap70(idx_H2j) = 0d0
    get_kevap70(idx_HEj) = 0d0
    get_kevap70(idx_HNOj) = 0d0
    get_kevap70(idx_H2NOj) = 0d0
    get_kevap70(idx_H3j) = 0d0
    get_kevap70(idx_H3COj) = 0d0
    get_kevap70(idx_H3Oj) = 0d0
    get_kevap70(idx_HCNHj) = 0d0
    get_kevap70(idx_HCO2j) = 0d0
    get_kevap70(idx_HEHj) = 0d0
    get_kevap70(idx_N2Hj) = 0d0
    get_kevap70(idx_O2Hj) = 0d0
    get_kevap70(idx_CR) = 0d0
    get_kevap70(idx_g) = 0d0
    get_kevap70(idx_Tgas) = 0d0
    get_kevap70(idx_dummy) = 0d0

  end function get_kevap70

  !************************
  !get verbatim reaction names
  function get_rnames()
    use krome_commons
    implicit none
    character*50::get_rnames(nrea)

    !reaction names are loaded from file
    get_rnames(:) = reactionNames(:)

  end function get_rnames

  !************************
  !get species names
  function get_names()
    use krome_commons
    implicit none
    character*16::get_names(nspec)

    get_names(1) = "E"
    get_names(2) = "CH"
    get_names(3) = "O"
    get_names(4) = "HNC"
    get_names(5) = "HCN"
    get_names(6) = "H2"
    get_names(7) = "C"
    get_names(8) = "H"
    get_names(9) = "H2O"
    get_names(10) = "OH"
    get_names(11) = "O2"
    get_names(12) = "CH2"
    get_names(13) = "H2CO"
    get_names(14) = "HCO"
    get_names(15) = "MG"
    get_names(16) = "NH3"
    get_names(17) = "NO"
    get_names(18) = "CN"
    get_names(19) = "CO"
    get_names(20) = "N2"
    get_names(21) = "NH2"
    get_names(22) = "CH3"
    get_names(23) = "CH4"
    get_names(24) = "N"
    get_names(25) = "NH"
    get_names(26) = "HE"
    get_names(27) = "HNO"
    get_names(28) = "CH3OH"
    get_names(29) = "CO2"
    get_names(30) = "H2CN"
    get_names(31) = "HNCO"
    get_names(32) = "NO2"
    get_names(33) = "O2H"
    get_names(34) = "OCN"
    get_names(35) = "CH3OH_DUST"
    get_names(36) = "HNCO_DUST"
    get_names(37) = "H2CO_DUST"
    get_names(38) = "CH4_DUST"
    get_names(39) = "CO_DUST"
    get_names(40) = "H2O_DUST"
    get_names(41) = "NO_DUST"
    get_names(42) = "CO2_DUST"
    get_names(43) = "N2_DUST"
    get_names(44) = "HCN_DUST"
    get_names(45) = "NH3_DUST"
    get_names(46) = "O2_DUST"
    get_names(47) = "NO2_DUST"
    get_names(48) = "HNO_DUST"
    get_names(49) = "O2H_DUST"
    get_names(50) = "H2CN_DUST"
    get_names(51) = "MG_DUST"
    get_names(52) = "HNC_DUST"
    get_names(53) = "E_DUST"
    get_names(54) = "HCO+"
    get_names(55) = "H+"
    get_names(56) = "HOC+"
    get_names(57) = "C+"
    get_names(58) = "CH2+"
    get_names(59) = "CH+"
    get_names(60) = "H2CO+"
    get_names(61) = "MG+"
    get_names(62) = "NH3+"
    get_names(63) = "NO+"
    get_names(64) = "CN+"
    get_names(65) = "CO+"
    get_names(66) = "N2+"
    get_names(67) = "O2+"
    get_names(68) = "H2O+"
    get_names(69) = "NH2+"
    get_names(70) = "O+"
    get_names(71) = "OH+"
    get_names(72) = "CH3+"
    get_names(73) = "CH4+"
    get_names(74) = "N+"
    get_names(75) = "HCN+"
    get_names(76) = "NH+"
    get_names(77) = "H2+"
    get_names(78) = "HE+"
    get_names(79) = "HNO+"
    get_names(80) = "H2NO+"
    get_names(81) = "H3+"
    get_names(82) = "H3CO+"
    get_names(83) = "H3O+"
    get_names(84) = "HCNH+"
    get_names(85) = "HCO2+"
    get_names(86) = "HEH+"
    get_names(87) = "N2H+"
    get_names(88) = "O2H+"
    get_names(89) = "CR"
    get_names(90) = "g"
    get_names(91) = "Tgas"
    get_names(92) = "dummy"

  end function get_names

  !************************
  !get cooling names list (empty element if cooling not present)
  function get_cooling_names()
    use krome_commons
    implicit none
    character*16::get_cooling_names(ncools)

    get_cooling_names(:) = ""

    get_cooling_names(idx_cool_h2) = "H2"
    get_cooling_names(idx_cool_h2gp) = "H2GP"
    get_cooling_names(idx_cool_atomic) = "ATOMIC"
    get_cooling_names(idx_cool_cen) = "CEN"
    get_cooling_names(idx_cool_hd) = "HD"
    get_cooling_names(idx_cool_metal) = "METAL"
    get_cooling_names(idx_cool_z) = "Z"
    get_cooling_names(idx_cool_dh) = "DH"
    get_cooling_names(idx_cool_enthalpic) = "ENTHALPIC"
    get_cooling_names(idx_cool_dust) = "DUST"
    get_cooling_names(idx_cool_compton) = "COMPTON"
    get_cooling_names(idx_cool_cie) = "CIE"
    get_cooling_names(idx_cool_cont) = "CONT"
    get_cooling_names(idx_cool_continuum) = "CONTINUUM"
    get_cooling_names(idx_cool_expansion) = "EXPANSION"
    get_cooling_names(idx_cool_exp) = "EXP"
    get_cooling_names(idx_cool_ff) = "FF"
    get_cooling_names(idx_cool_bss) = "BSS"
    get_cooling_names(idx_cool_custom) = "CUSTOM"
    get_cooling_names(idx_cool_co) = "CO"
    get_cooling_names(idx_cool_zcie) = "ZCIE"
    get_cooling_names(idx_cool_zcienouv) = "ZCIENOUV"
    get_cooling_names(idx_cool_zextend) = "ZEXTEND"
    get_cooling_names(idx_cool_gh) = "GH"

  end function get_cooling_names

  !************************
  !get heating names list (empty element if heating not present)
  function get_heating_names()
    use krome_commons
    implicit none
    character*16::get_heating_names(nheats)

    get_heating_names(:) = ""

    get_heating_names(idx_heat_chem) = "CHEM"
    get_heating_names(idx_heat_compress) = "COMPRESS"
    get_heating_names(idx_heat_compr) = "COMPR"
    get_heating_names(idx_heat_photo) = "PHOTO"
    get_heating_names(idx_heat_dh) = "DH"
    get_heating_names(idx_heat_enthalpic) = "ENTHALPIC"
    get_heating_names(idx_heat_av) = "AV"
    get_heating_names(idx_heat_photoav) = "PHOTOAV"
    get_heating_names(idx_heat_cr) = "CR"
    get_heating_names(idx_heat_dust) = "DUST"
    get_heating_names(idx_heat_xray) = "XRAY"
    get_heating_names(idx_heat_viscous) = "VISCOUS"
    get_heating_names(idx_heat_visc) = "VISC"
    get_heating_names(idx_heat_custom) = "CUSTOM"
    get_heating_names(idx_heat_zcie) = "ZCIE"

  end function get_heating_names

  !******************************
  !get the total number of H nuclei
  function get_Hnuclei(n)
    use krome_commons
    real*8::n(:),get_Hnuclei,nH

    nH = n(idx_CH) + &
        n(idx_HNC) + &
        n(idx_HCN) + &
        n(idx_H2)*2d0 + &
        n(idx_H) + &
        n(idx_H2O)*2d0 + &
        n(idx_OH) + &
        n(idx_CH2)*2d0 + &
        n(idx_H2CO)*2d0 + &
        n(idx_HCO) + &
        n(idx_NH3)*3d0 + &
        n(idx_NH2)*2d0 + &
        n(idx_CH3)*3d0 + &
        n(idx_CH4)*4d0 + &
        n(idx_NH) + &
        n(idx_HNO) + &
        n(idx_CH3OH)*4d0 + &
        n(idx_H2CN)*2d0 + &
        n(idx_HNCO) + &
        n(idx_O2H) + &
        n(idx_CH3OH_DUST)*4d0 + &
        n(idx_HNCO_DUST) + &
        n(idx_H2CO_DUST)*2d0 + &
        n(idx_CH4_DUST)*4d0 + &
        n(idx_H2O_DUST)*2d0 + &
        n(idx_HCN_DUST) + &
        n(idx_NH3_DUST)*3d0 + &
        n(idx_HNO_DUST) + &
        n(idx_O2H_DUST) + &
        n(idx_H2CN_DUST)*2d0 + &
        n(idx_HNC_DUST) + &
        n(idx_HCOj) + &
        n(idx_Hj) + &
        n(idx_HOCj) + &
        n(idx_CH2j)*2d0 + &
        n(idx_CHj) + &
        n(idx_H2COj)*2d0 + &
        n(idx_NH3j)*3d0 + &
        n(idx_H2Oj)*2d0 + &
        n(idx_NH2j)*2d0 + &
        n(idx_OHj) + &
        n(idx_CH3j)*3d0 + &
        n(idx_CH4j)*4d0 + &
        n(idx_HCNj) + &
        n(idx_NHj) + &
        n(idx_H2j)*2d0 + &
        n(idx_HNOj) + &
        n(idx_H2NOj)*2d0 + &
        n(idx_H3j)*3d0 + &
        n(idx_H3COj)*3d0 + &
        n(idx_H3Oj)*3d0 + &
        n(idx_HCNHj)*2d0 + &
        n(idx_HCO2j) + &
        n(idx_HEHj) + &
        n(idx_N2Hj) + &
        n(idx_O2Hj)
    get_Hnuclei = nH

  end function get_Hnuclei

  !******************************
  !get the total number of mantles
  function get_mantle(n)
    use krome_commons
    real*8::n(:),get_mantle,mantle

    mantle = n(idx_CH3OH_DUST) + &
        n(idx_HNCO_DUST) + &
        n(idx_H2CO_DUST) + &
        n(idx_CH4_DUST) + &
        n(idx_CO_DUST) + &
        n(idx_H2O_DUST) + &
        n(idx_NO_DUST) + &
        n(idx_CO2_DUST) + &
        n(idx_N2_DUST) + &
        n(idx_HCN_DUST) + &
        n(idx_NH3_DUST) + &
        n(idx_O2_DUST) + &
        n(idx_NO2_DUST) + &
        n(idx_HNO_DUST) + &
        n(idx_O2H_DUST) + &
        n(idx_H2CN_DUST) + &
        n(idx_MG_DUST) + &
        n(idx_HNC_DUST) + &
        n(idx_E_DUST)
    get_mantle = mantle

  end function get_mantle

  !***************************
  function get_zatoms()
    use krome_commons
    implicit none
    integer::get_zatoms(nspec)

    get_zatoms(1) = 0	!E
    get_zatoms(2) = 7	!CH
    get_zatoms(3) = 8	!O
    get_zatoms(4) = 14	!HNC
    get_zatoms(5) = 14	!HCN
    get_zatoms(6) = 2	!H2
    get_zatoms(7) = 6	!C
    get_zatoms(8) = 1	!H
    get_zatoms(9) = 10	!H2O
    get_zatoms(10) = 9	!OH
    get_zatoms(11) = 16	!O2
    get_zatoms(12) = 8	!CH2
    get_zatoms(13) = 16	!H2CO
    get_zatoms(14) = 15	!HCO
    get_zatoms(15) = 12	!MG
    get_zatoms(16) = 10	!NH3
    get_zatoms(17) = 15	!NO
    get_zatoms(18) = 13	!CN
    get_zatoms(19) = 14	!CO
    get_zatoms(20) = 14	!N2
    get_zatoms(21) = 9	!NH2
    get_zatoms(22) = 9	!CH3
    get_zatoms(23) = 10	!CH4
    get_zatoms(24) = 7	!N
    get_zatoms(25) = 8	!NH
    get_zatoms(26) = 2	!HE
    get_zatoms(27) = 16	!HNO
    get_zatoms(28) = 18	!CH3OH
    get_zatoms(29) = 22	!CO2
    get_zatoms(30) = 15	!H2CN
    get_zatoms(31) = 22	!HNCO
    get_zatoms(32) = 23	!NO2
    get_zatoms(33) = 17	!O2H
    get_zatoms(34) = 21	!OCN
    get_zatoms(35) = 18	!CH3OH_DUST
    get_zatoms(36) = 22	!HNCO_DUST
    get_zatoms(37) = 16	!H2CO_DUST
    get_zatoms(38) = 10	!CH4_DUST
    get_zatoms(39) = 14	!CO_DUST
    get_zatoms(40) = 10	!H2O_DUST
    get_zatoms(41) = 15	!NO_DUST
    get_zatoms(42) = 22	!CO2_DUST
    get_zatoms(43) = 14	!N2_DUST
    get_zatoms(44) = 14	!HCN_DUST
    get_zatoms(45) = 10	!NH3_DUST
    get_zatoms(46) = 16	!O2_DUST
    get_zatoms(47) = 23	!NO2_DUST
    get_zatoms(48) = 16	!HNO_DUST
    get_zatoms(49) = 17	!O2H_DUST
    get_zatoms(50) = 15	!H2CN_DUST
    get_zatoms(51) = 12	!MG_DUST
    get_zatoms(52) = 14	!HNC_DUST
    get_zatoms(53) = 0	!E_DUST
    get_zatoms(54) = 15	!HCO+
    get_zatoms(55) = 1	!H+
    get_zatoms(56) = 15	!HOC+
    get_zatoms(57) = 6	!C+
    get_zatoms(58) = 8	!CH2+
    get_zatoms(59) = 7	!CH+
    get_zatoms(60) = 16	!H2CO+
    get_zatoms(61) = 12	!MG+
    get_zatoms(62) = 10	!NH3+
    get_zatoms(63) = 15	!NO+
    get_zatoms(64) = 13	!CN+
    get_zatoms(65) = 14	!CO+
    get_zatoms(66) = 14	!N2+
    get_zatoms(67) = 16	!O2+
    get_zatoms(68) = 10	!H2O+
    get_zatoms(69) = 9	!NH2+
    get_zatoms(70) = 8	!O+
    get_zatoms(71) = 9	!OH+
    get_zatoms(72) = 9	!CH3+
    get_zatoms(73) = 10	!CH4+
    get_zatoms(74) = 7	!N+
    get_zatoms(75) = 14	!HCN+
    get_zatoms(76) = 8	!NH+
    get_zatoms(77) = 2	!H2+
    get_zatoms(78) = 2	!HE+
    get_zatoms(79) = 16	!HNO+
    get_zatoms(80) = 17	!H2NO+
    get_zatoms(81) = 3	!H3+
    get_zatoms(82) = 17	!H3CO+
    get_zatoms(83) = 11	!H3O+
    get_zatoms(84) = 15	!HCNH+
    get_zatoms(85) = 23	!HCO2+
    get_zatoms(86) = 3	!HEH+
    get_zatoms(87) = 15	!N2H+
    get_zatoms(88) = 17	!O2H+
    get_zatoms(89) = 0	!CR
    get_zatoms(90) = 0	!g
    get_zatoms(91) = 0	!Tgas
    get_zatoms(92) = 0	!dummy

  end function get_zatoms

  !******************************
  function get_qeff()
    use krome_commons
    implicit none
    real*8::get_qeff(nrea)

    get_qeff(:) = 0e0

  end function get_qeff

  !**************************
  function get_free_fall_time(n)
    use krome_constants
    use krome_commons
    implicit none
    real*8::n(:),m(nspec)
    real*8::rhogas,get_free_fall_time

    m(:) = get_mass()
    rhogas = sum(n(1:nmols)*m(1:nmols))
    get_free_fall_time = sqrt(3d0*pi/32d0/gravity/rhogas)

  end function get_free_fall_time

  !**************************
  function get_free_fall_time_rho(rhogas)
    use krome_constants
    implicit none
    real*8::rhogas,get_free_fall_time_rho

    get_free_fall_time_rho = sqrt(3d0*pi/32d0/gravity/rhogas)

  end function get_free_fall_time_rho

  !********************************
  function get_jeans_length(n,Tgas)
    !get jeans length in cm
    use krome_constants
    use krome_commons
    implicit none
    real*8::n(:),Tgas,mu,rhogas
    real*8::m(nspec),get_jeans_length
    m(:) = get_mass()
    rhogas = max(sum(n(1:nmols)*m(1:nmols)),1d-40)
    mu = get_mu_rho(n(:),rhogas)
    get_jeans_length = sqrt(pi*boltzmann_erg*Tgas/rhogas&
        /p_mass/gravity/mu)

  end function get_jeans_length

  !********************************
  function get_jeans_length_rho(n,Tgas,rhogas)
    !get jeans length in cm
    use krome_constants
    use krome_commons
    implicit none
    real*8::n(:),Tgas,mu,rhogas
    real*8::get_jeans_length_rho

    mu = get_mu_rho(n(:),rhogas)
    get_jeans_length_rho = sqrt(pi*boltzmann_erg*Tgas/rhogas&
        /p_mass/gravity/mu)

  end function get_jeans_length_rho

  !***************************
  !number density to column density conversion
  function num2col(ncalc,n)
    use krome_commons
    implicit none
    real*8::num2col,ncalc,n(:),Tgas
    Tgas = max(n(idx_Tgas),phys_Tcmb)

    num2col = 1.87d21*(max(ncalc,1d-40)*1d-3)**(2./3.)

  end function num2col

  !***********************
  !column density to number density conversion
  function col2num(ncalc,n)
    use krome_commons
    implicit none
    real*8::col2num,ncalc,n(:),Tgas
    Tgas = max(n(idx_Tgas),phys_Tcmb)

    col2num = 1d3 * (max(ncalc,1d-40)/1.87d21)**1.5

  end function col2num

  !************************
  !get electrons by balancing charges
  function get_electrons(n)
    use krome_commons
    implicit none
    real*8::get_electrons,n(nspec)

    get_electrons =  + n(idx_HCOj) &
        + n(idx_Hj) &
        + n(idx_HOCj) &
        + n(idx_Cj) &
        + n(idx_CH2j) &
        + n(idx_CHj) &
        + n(idx_H2COj) &
        + n(idx_MGj) &
        + n(idx_NH3j) &
        + n(idx_NOj) &
        + n(idx_CNj) &
        + n(idx_COj) &
        + n(idx_N2j) &
        + n(idx_O2j) &
        + n(idx_H2Oj) &
        + n(idx_NH2j) &
        + n(idx_Oj) &
        + n(idx_OHj) &
        + n(idx_CH3j) &
        + n(idx_CH4j) &
        + n(idx_Nj) &
        + n(idx_HCNj) &
        + n(idx_NHj) &
        + n(idx_H2j) &
        + n(idx_HEj) &
        + n(idx_HNOj) &
        + n(idx_H2NOj) &
        + n(idx_H3j) &
        + n(idx_H3COj) &
        + n(idx_H3Oj) &
        + n(idx_HCNHj) &
        + n(idx_HCO2j) &
        + n(idx_HEHj) &
        + n(idx_N2Hj) &
        + n(idx_O2Hj)
    get_electrons = max(get_electrons,0d0)

  end function get_electrons

  !************************
  !get species charges
  function get_charges()
    use krome_commons
    implicit none
    integer::get_charges(nspec)

    get_charges(1) = -1.d0 	!E
    get_charges(2) = 0.d0 	!CH
    get_charges(3) = 0.d0 	!O
    get_charges(4) = 0.d0 	!HNC
    get_charges(5) = 0.d0 	!HCN
    get_charges(6) = 0.d0 	!H2
    get_charges(7) = 0.d0 	!C
    get_charges(8) = 0.d0 	!H
    get_charges(9) = 0.d0 	!H2O
    get_charges(10) = 0.d0 	!OH
    get_charges(11) = 0.d0 	!O2
    get_charges(12) = 0.d0 	!CH2
    get_charges(13) = 0.d0 	!H2CO
    get_charges(14) = 0.d0 	!HCO
    get_charges(15) = 0.d0 	!MG
    get_charges(16) = 0.d0 	!NH3
    get_charges(17) = 0.d0 	!NO
    get_charges(18) = 0.d0 	!CN
    get_charges(19) = 0.d0 	!CO
    get_charges(20) = 0.d0 	!N2
    get_charges(21) = 0.d0 	!NH2
    get_charges(22) = 0.d0 	!CH3
    get_charges(23) = 0.d0 	!CH4
    get_charges(24) = 0.d0 	!N
    get_charges(25) = 0.d0 	!NH
    get_charges(26) = 0.d0 	!HE
    get_charges(27) = 0.d0 	!HNO
    get_charges(28) = 0.d0 	!CH3OH
    get_charges(29) = 0.d0 	!CO2
    get_charges(30) = 0.d0 	!H2CN
    get_charges(31) = 0.d0 	!HNCO
    get_charges(32) = 0.d0 	!NO2
    get_charges(33) = 0.d0 	!O2H
    get_charges(34) = 0.d0 	!OCN
    get_charges(35) = 0.d0 	!CH3OH_DUST
    get_charges(36) = 0.d0 	!HNCO_DUST
    get_charges(37) = 0.d0 	!H2CO_DUST
    get_charges(38) = 0.d0 	!CH4_DUST
    get_charges(39) = 0.d0 	!CO_DUST
    get_charges(40) = 0.d0 	!H2O_DUST
    get_charges(41) = 0.d0 	!NO_DUST
    get_charges(42) = 0.d0 	!CO2_DUST
    get_charges(43) = 0.d0 	!N2_DUST
    get_charges(44) = 0.d0 	!HCN_DUST
    get_charges(45) = 0.d0 	!NH3_DUST
    get_charges(46) = 0.d0 	!O2_DUST
    get_charges(47) = 0.d0 	!NO2_DUST
    get_charges(48) = 0.d0 	!HNO_DUST
    get_charges(49) = 0.d0 	!O2H_DUST
    get_charges(50) = 0.d0 	!H2CN_DUST
    get_charges(51) = 0.d0 	!MG_DUST
    get_charges(52) = 0.d0 	!HNC_DUST
    get_charges(53) = 0.d0 	!E_DUST
    get_charges(54) = 1.d0 	!HCO+
    get_charges(55) = 1.d0 	!H+
    get_charges(56) = 1.d0 	!HOC+
    get_charges(57) = 1.d0 	!C+
    get_charges(58) = 1.d0 	!CH2+
    get_charges(59) = 1.d0 	!CH+
    get_charges(60) = 1.d0 	!H2CO+
    get_charges(61) = 1.d0 	!MG+
    get_charges(62) = 1.d0 	!NH3+
    get_charges(63) = 1.d0 	!NO+
    get_charges(64) = 1.d0 	!CN+
    get_charges(65) = 1.d0 	!CO+
    get_charges(66) = 1.d0 	!N2+
    get_charges(67) = 1.d0 	!O2+
    get_charges(68) = 1.d0 	!H2O+
    get_charges(69) = 1.d0 	!NH2+
    get_charges(70) = 1.d0 	!O+
    get_charges(71) = 1.d0 	!OH+
    get_charges(72) = 1.d0 	!CH3+
    get_charges(73) = 1.d0 	!CH4+
    get_charges(74) = 1.d0 	!N+
    get_charges(75) = 1.d0 	!HCN+
    get_charges(76) = 1.d0 	!NH+
    get_charges(77) = 1.d0 	!H2+
    get_charges(78) = 1.d0 	!HE+
    get_charges(79) = 1.d0 	!HNO+
    get_charges(80) = 1.d0 	!H2NO+
    get_charges(81) = 1.d0 	!H3+
    get_charges(82) = 1.d0 	!H3CO+
    get_charges(83) = 1.d0 	!H3O+
    get_charges(84) = 1.d0 	!HCNH+
    get_charges(85) = 1.d0 	!HCO2+
    get_charges(86) = 1.d0 	!HEH+
    get_charges(87) = 1.d0 	!N2H+
    get_charges(88) = 1.d0 	!O2H+
    get_charges(89) = 0.d0 	!CR
    get_charges(90) = 0.d0 	!g
    get_charges(91) = 0.d0 	!Tgas
    get_charges(92) = 0.d0 	!dummy

  end function get_charges

  !*****************************
  ! get metallicity using C as reference
  function get_metallicityC(n)
    use krome_commons
    implicit none
    real*8::n(:),get_metallicityC,zC,nH

    nH = get_Hnuclei(n(:))

    zC = n(idx_CH) &
        + n(idx_HNC) &
        + n(idx_HCN) &
        + n(idx_C) &
        + n(idx_CH2) &
        + n(idx_H2CO) &
        + n(idx_HCO) &
        + n(idx_CN) &
        + n(idx_CO) &
        + n(idx_CH3) &
        + n(idx_CH4) &
        + n(idx_CH3OH) &
        + n(idx_CO2) &
        + n(idx_H2CN) &
        + n(idx_HNCO) &
        + n(idx_OCN) &
        + n(idx_CH3OH_DUST) &
        + n(idx_HNCO_DUST) &
        + n(idx_H2CO_DUST) &
        + n(idx_CH4_DUST) &
        + n(idx_CO_DUST) &
        + n(idx_CO2_DUST) &
        + n(idx_HCN_DUST) &
        + n(idx_H2CN_DUST) &
        + n(idx_HNC_DUST) &
        + n(idx_HCOj) &
        + n(idx_HOCj) &
        + n(idx_Cj) &
        + n(idx_CH2j) &
        + n(idx_CHj) &
        + n(idx_H2COj) &
        + n(idx_CNj) &
        + n(idx_COj) &
        + n(idx_CH3j) &
        + n(idx_CH4j) &
        + n(idx_HCNj) &
        + n(idx_H3COj) &
        + n(idx_HCNHj) &
        + n(idx_HCO2j)

    zC = max(zC, 0d0)

    get_metallicityC = log10(zC/nH+1d-40) - (-3.57)

    phys_metallicity = get_metallicityC

  end function get_metallicityC

  !*****************************
  ! get metallicity using Mg as reference
  function get_metallicityMg(n)
    use krome_commons
    implicit none
    real*8::n(:),get_metallicityMg,zMg,nH

    nH = get_Hnuclei(n(:))

    zMg = n(idx_MG) &
        + n(idx_MG_DUST) &
        + n(idx_MGj)

    zMg = max(zMg, 0d0)

    get_metallicityMg = log10(zMg/nH+1d-40) - (-4.4)

    phys_metallicity = get_metallicityMg

  end function get_metallicityMg

  !*****************************
  ! get metallicity using O as reference
  function get_metallicityO(n)
    use krome_commons
    implicit none
    real*8::n(:),get_metallicityO,zO,nH

    nH = get_Hnuclei(n(:))

    zO = n(idx_O) &
        + n(idx_H2O) &
        + n(idx_OH) &
        + 2d0*n(idx_O2) &
        + n(idx_H2CO) &
        + n(idx_HCO) &
        + n(idx_NO) &
        + n(idx_CO) &
        + n(idx_HNO) &
        + n(idx_CH3OH) &
        + 2d0*n(idx_CO2) &
        + n(idx_HNCO) &
        + 2d0*n(idx_NO2) &
        + 2d0*n(idx_O2H) &
        + n(idx_OCN) &
        + n(idx_CH3OH_DUST) &
        + n(idx_HNCO_DUST) &
        + n(idx_H2CO_DUST) &
        + n(idx_CO_DUST) &
        + n(idx_H2O_DUST) &
        + n(idx_NO_DUST) &
        + 2d0*n(idx_CO2_DUST) &
        + 2d0*n(idx_O2_DUST) &
        + 2d0*n(idx_NO2_DUST) &
        + n(idx_HNO_DUST) &
        + 2d0*n(idx_O2H_DUST) &
        + n(idx_HCOj) &
        + n(idx_HOCj) &
        + n(idx_H2COj) &
        + n(idx_NOj) &
        + n(idx_COj) &
        + 2d0*n(idx_O2j) &
        + n(idx_H2Oj) &
        + n(idx_Oj) &
        + n(idx_OHj) &
        + n(idx_HNOj) &
        + n(idx_H2NOj) &
        + n(idx_H3COj) &
        + n(idx_H3Oj) &
        + 2d0*n(idx_HCO2j) &
        + 2d0*n(idx_O2Hj)

    zO = max(zO, 0d0)

    get_metallicityO = log10(zO/nH+1d-40) - (-3.31)

    phys_metallicity = get_metallicityO

  end function get_metallicityO

  !*****************************
  ! get metallicity using N as reference
  function get_metallicityN(n)
    use krome_commons
    implicit none
    real*8::n(:),get_metallicityN,zN,nH

    nH = get_Hnuclei(n(:))

    zN = n(idx_HNC) &
        + n(idx_HCN) &
        + n(idx_NH3) &
        + n(idx_NO) &
        + n(idx_CN) &
        + 2d0*n(idx_N2) &
        + n(idx_NH2) &
        + n(idx_N) &
        + n(idx_NH) &
        + n(idx_HNO) &
        + n(idx_H2CN) &
        + n(idx_HNCO) &
        + n(idx_NO2) &
        + n(idx_OCN) &
        + n(idx_HNCO_DUST) &
        + n(idx_NO_DUST) &
        + 2d0*n(idx_N2_DUST) &
        + n(idx_HCN_DUST) &
        + n(idx_NH3_DUST) &
        + n(idx_NO2_DUST) &
        + n(idx_HNO_DUST) &
        + n(idx_H2CN_DUST) &
        + n(idx_HNC_DUST) &
        + n(idx_NH3j) &
        + n(idx_NOj) &
        + n(idx_CNj) &
        + 2d0*n(idx_N2j) &
        + n(idx_NH2j) &
        + n(idx_Nj) &
        + n(idx_HCNj) &
        + n(idx_NHj) &
        + n(idx_HNOj) &
        + n(idx_H2NOj) &
        + n(idx_HCNHj) &
        + 2d0*n(idx_N2Hj)

    zN = max(zN, 0d0)

    get_metallicityN = log10(zN/nH+1d-40) - (-4.17)

    phys_metallicity = get_metallicityN

  end function get_metallicityN

end module krome_getphys
