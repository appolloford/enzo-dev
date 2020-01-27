!KROME_DRIVER

subroutine krome_driver(d, e, ge, u, v, w, &
      De, CHI, OI, HNCI, &
      HCNI, H2I, CI, HI, &
      H2OI, OHI, O2I, CH2I, &
      H2COI, HCOI, MGI, NH3I, &
      NOI, CNI, COI, N2I, &
      NH2I, CH3I, CH4I, NI, &
      NHI, HeI, HNOI, CH3OHI, &
      CO2I, H2CNI, HNCOI, NO2I, &
      O2HI, OCNI, CH3OH_DUSTI, HNCO_DUSTI, &
      H2CO_DUSTI, CH4_DUSTI, CO_DUSTI, H2O_DUSTI, &
      NO_DUSTI, CO2_DUSTI, N2_DUSTI, HCN_DUSTI, &
      NH3_DUSTI, O2_DUSTI, NO2_DUSTI, HNO_DUSTI, &
      O2H_DUSTI, H2CN_DUSTI, MG_DUSTI, HNC_DUSTI, &
      E_DUSTI, HCOII, HII, HOCII, &
      CII, CH2II, CHII, H2COII, &
      MGII, NH3II, NOII, CNII, &
      COII, N2II, O2II, H2OII, &
      NH2II, OII, OHII, CH3II, &
      CH4II, NII, HCNII, NHII, &
      H2II, HeII, HNOII, H2NOII, &
      H3II, H3COII, H3OII, HCNHII, &
      HCO2II, HeHII, N2HII, O2HII, &
      in, jn, kn, imethod, &
      idual, idim, &
      is, js, ks, ie, je, ke, &
      dt, aye, &
      utem, uxyz, uaye, urho, utim, &
      gamma, fh, dtoh)

  !     SOLVE MULTI-SPECIES RATE EQUATIONS AND RADIATIVE COOLING
  !
  !     2014, KROME DEVELOPERS to interface the package with ENZO
  !
  !     PURPOSE:
  !     Solve the multi-species rate and cool equations via KROME.
  !
  !     INPUTS:
  !     in,jn,kn - dimensions of 3D fields
  !
  !     d        - total density field
  !
  !     is,ie    - start and end indices of active region (zero based)
  !     idual    - dual energy formalism flag (0 = off, 1 = on)
  !     idim     - dimensionality (rank) of problem
  !     imethod  - Hydro method (0 = PPMDE, 2 = ZEUS-type)
  !
  !     fh       - Hydrogen mass fraction (typically 0.76)
  !     dtoh     - Deuterium to H mass ratio
  !     dt       - timestep to integrate over
  !     aye      - expansion factor (in code units)
  !
  !     utim     - time units (i.e. code units to CGS conversion factor)
  !     uaye     - expansion factor conversion factor (uaye = 1/(1+zinit))
  !     urho     - density units
  !     uxyz     - length units
  !     utem     - temperature(-like) units
  !
  !     OUTPUTS:
  !     update chemical abundances densities (HI, HII, etc) and energy
  !
  !     PARAMETERS:
  !     mh      - H mass in cgs units
  !
  !-----------------------------------------------------------------------
  !     USE KROME
  use krome_main
  use krome_user
  use krome_constants

  implicit none

  real*8,parameter::mh=p_mass !mass_h
  real*8::dom,factor,tgas,tgasold,krome_x(krome_nmols)
  real*8::dt,dt_hydro,idom,edot,krome_tiny
  real*8::d(in,jn,kn),e(in,jn,kn),ge(in,jn,kn)
  real*8::u(in,jn,kn),v(in,jn,kn),w(in,jn,kn)
  real*8::aye,utem,uxyz,uaye,urho,utim,gamma,fh,dtoh
  integer::in,jn,kn,imethod,idual,is,js,ks,ie,je,ke,idim
  integer::i,j,k

  real*8::De(in,jn,kn)
  real*8::CHI(in,jn,kn)
  real*8::OI(in,jn,kn)
  real*8::HNCI(in,jn,kn)
  real*8::HCNI(in,jn,kn)
  real*8::H2I(in,jn,kn)
  real*8::CI(in,jn,kn)
  real*8::HI(in,jn,kn)
  real*8::H2OI(in,jn,kn)
  real*8::OHI(in,jn,kn)
  real*8::O2I(in,jn,kn)
  real*8::CH2I(in,jn,kn)
  real*8::H2COI(in,jn,kn)
  real*8::HCOI(in,jn,kn)
  real*8::MGI(in,jn,kn)
  real*8::NH3I(in,jn,kn)
  real*8::NOI(in,jn,kn)
  real*8::CNI(in,jn,kn)
  real*8::COI(in,jn,kn)
  real*8::N2I(in,jn,kn)
  real*8::NH2I(in,jn,kn)
  real*8::CH3I(in,jn,kn)
  real*8::CH4I(in,jn,kn)
  real*8::NI(in,jn,kn)
  real*8::NHI(in,jn,kn)
  real*8::HeI(in,jn,kn)
  real*8::HNOI(in,jn,kn)
  real*8::CH3OHI(in,jn,kn)
  real*8::CO2I(in,jn,kn)
  real*8::H2CNI(in,jn,kn)
  real*8::HNCOI(in,jn,kn)
  real*8::NO2I(in,jn,kn)
  real*8::O2HI(in,jn,kn)
  real*8::OCNI(in,jn,kn)
  real*8::CH3OH_DUSTI(in,jn,kn)
  real*8::HNCO_DUSTI(in,jn,kn)
  real*8::H2CO_DUSTI(in,jn,kn)
  real*8::CH4_DUSTI(in,jn,kn)
  real*8::CO_DUSTI(in,jn,kn)
  real*8::H2O_DUSTI(in,jn,kn)
  real*8::NO_DUSTI(in,jn,kn)
  real*8::CO2_DUSTI(in,jn,kn)
  real*8::N2_DUSTI(in,jn,kn)
  real*8::HCN_DUSTI(in,jn,kn)
  real*8::NH3_DUSTI(in,jn,kn)
  real*8::O2_DUSTI(in,jn,kn)
  real*8::NO2_DUSTI(in,jn,kn)
  real*8::HNO_DUSTI(in,jn,kn)
  real*8::O2H_DUSTI(in,jn,kn)
  real*8::H2CN_DUSTI(in,jn,kn)
  real*8::MG_DUSTI(in,jn,kn)
  real*8::HNC_DUSTI(in,jn,kn)
  real*8::E_DUSTI(in,jn,kn)
  real*8::HCOII(in,jn,kn)
  real*8::HII(in,jn,kn)
  real*8::HOCII(in,jn,kn)
  real*8::CII(in,jn,kn)
  real*8::CH2II(in,jn,kn)
  real*8::CHII(in,jn,kn)
  real*8::H2COII(in,jn,kn)
  real*8::MGII(in,jn,kn)
  real*8::NH3II(in,jn,kn)
  real*8::NOII(in,jn,kn)
  real*8::CNII(in,jn,kn)
  real*8::COII(in,jn,kn)
  real*8::N2II(in,jn,kn)
  real*8::O2II(in,jn,kn)
  real*8::H2OII(in,jn,kn)
  real*8::NH2II(in,jn,kn)
  real*8::OII(in,jn,kn)
  real*8::OHII(in,jn,kn)
  real*8::CH3II(in,jn,kn)
  real*8::CH4II(in,jn,kn)
  real*8::NII(in,jn,kn)
  real*8::HCNII(in,jn,kn)
  real*8::NHII(in,jn,kn)
  real*8::H2II(in,jn,kn)
  real*8::HeII(in,jn,kn)
  real*8::HNOII(in,jn,kn)
  real*8::H2NOII(in,jn,kn)
  real*8::H3II(in,jn,kn)
  real*8::H3COII(in,jn,kn)
  real*8::H3OII(in,jn,kn)
  real*8::HCNHII(in,jn,kn)
  real*8::HCO2II(in,jn,kn)
  real*8::HeHII(in,jn,kn)
  real*8::N2HII(in,jn,kn)
  real*8::O2HII(in,jn,kn)

  !******************************

  !set units
  dom = urho*(aye**3)/mh

  !scaling factor for comoving->proper
  factor = aye**(-3)

  !check minimal value and comoving->proper
  krome_tiny = 1d-30
  do k = ks+1, ke+1
    do j = js+1, je+1
      do i = is+1, ie+1
        d(i,j,k) = d(i,j,k) * factor
        !scale comoving->proper
        De(i,j,k) = De(i,j,k) * factor
        CHI(i,j,k) = CHI(i,j,k) * factor
        OI(i,j,k) = OI(i,j,k) * factor
        HNCI(i,j,k) = HNCI(i,j,k) * factor
        HCNI(i,j,k) = HCNI(i,j,k) * factor
        H2I(i,j,k) = H2I(i,j,k) * factor
        CI(i,j,k) = CI(i,j,k) * factor
        HI(i,j,k) = HI(i,j,k) * factor
        H2OI(i,j,k) = H2OI(i,j,k) * factor
        OHI(i,j,k) = OHI(i,j,k) * factor
        O2I(i,j,k) = O2I(i,j,k) * factor
        CH2I(i,j,k) = CH2I(i,j,k) * factor
        H2COI(i,j,k) = H2COI(i,j,k) * factor
        HCOI(i,j,k) = HCOI(i,j,k) * factor
        MGI(i,j,k) = MGI(i,j,k) * factor
        NH3I(i,j,k) = NH3I(i,j,k) * factor
        NOI(i,j,k) = NOI(i,j,k) * factor
        CNI(i,j,k) = CNI(i,j,k) * factor
        COI(i,j,k) = COI(i,j,k) * factor
        N2I(i,j,k) = N2I(i,j,k) * factor
        NH2I(i,j,k) = NH2I(i,j,k) * factor
        CH3I(i,j,k) = CH3I(i,j,k) * factor
        CH4I(i,j,k) = CH4I(i,j,k) * factor
        NI(i,j,k) = NI(i,j,k) * factor
        NHI(i,j,k) = NHI(i,j,k) * factor
        HeI(i,j,k) = HeI(i,j,k) * factor
        HNOI(i,j,k) = HNOI(i,j,k) * factor
        CH3OHI(i,j,k) = CH3OHI(i,j,k) * factor
        CO2I(i,j,k) = CO2I(i,j,k) * factor
        H2CNI(i,j,k) = H2CNI(i,j,k) * factor
        HNCOI(i,j,k) = HNCOI(i,j,k) * factor
        NO2I(i,j,k) = NO2I(i,j,k) * factor
        O2HI(i,j,k) = O2HI(i,j,k) * factor
        OCNI(i,j,k) = OCNI(i,j,k) * factor
        CH3OH_DUSTI(i,j,k) = CH3OH_DUSTI(i,j,k) * factor
        HNCO_DUSTI(i,j,k) = HNCO_DUSTI(i,j,k) * factor
        H2CO_DUSTI(i,j,k) = H2CO_DUSTI(i,j,k) * factor
        CH4_DUSTI(i,j,k) = CH4_DUSTI(i,j,k) * factor
        CO_DUSTI(i,j,k) = CO_DUSTI(i,j,k) * factor
        H2O_DUSTI(i,j,k) = H2O_DUSTI(i,j,k) * factor
        NO_DUSTI(i,j,k) = NO_DUSTI(i,j,k) * factor
        CO2_DUSTI(i,j,k) = CO2_DUSTI(i,j,k) * factor
        N2_DUSTI(i,j,k) = N2_DUSTI(i,j,k) * factor
        HCN_DUSTI(i,j,k) = HCN_DUSTI(i,j,k) * factor
        NH3_DUSTI(i,j,k) = NH3_DUSTI(i,j,k) * factor
        O2_DUSTI(i,j,k) = O2_DUSTI(i,j,k) * factor
        NO2_DUSTI(i,j,k) = NO2_DUSTI(i,j,k) * factor
        HNO_DUSTI(i,j,k) = HNO_DUSTI(i,j,k) * factor
        O2H_DUSTI(i,j,k) = O2H_DUSTI(i,j,k) * factor
        H2CN_DUSTI(i,j,k) = H2CN_DUSTI(i,j,k) * factor
        MG_DUSTI(i,j,k) = MG_DUSTI(i,j,k) * factor
        HNC_DUSTI(i,j,k) = HNC_DUSTI(i,j,k) * factor
        E_DUSTI(i,j,k) = E_DUSTI(i,j,k) * factor
        HCOII(i,j,k) = HCOII(i,j,k) * factor
        HII(i,j,k) = HII(i,j,k) * factor
        HOCII(i,j,k) = HOCII(i,j,k) * factor
        CII(i,j,k) = CII(i,j,k) * factor
        CH2II(i,j,k) = CH2II(i,j,k) * factor
        CHII(i,j,k) = CHII(i,j,k) * factor
        H2COII(i,j,k) = H2COII(i,j,k) * factor
        MGII(i,j,k) = MGII(i,j,k) * factor
        NH3II(i,j,k) = NH3II(i,j,k) * factor
        NOII(i,j,k) = NOII(i,j,k) * factor
        CNII(i,j,k) = CNII(i,j,k) * factor
        COII(i,j,k) = COII(i,j,k) * factor
        N2II(i,j,k) = N2II(i,j,k) * factor
        O2II(i,j,k) = O2II(i,j,k) * factor
        H2OII(i,j,k) = H2OII(i,j,k) * factor
        NH2II(i,j,k) = NH2II(i,j,k) * factor
        OII(i,j,k) = OII(i,j,k) * factor
        OHII(i,j,k) = OHII(i,j,k) * factor
        CH3II(i,j,k) = CH3II(i,j,k) * factor
        CH4II(i,j,k) = CH4II(i,j,k) * factor
        NII(i,j,k) = NII(i,j,k) * factor
        HCNII(i,j,k) = HCNII(i,j,k) * factor
        NHII(i,j,k) = NHII(i,j,k) * factor
        H2II(i,j,k) = H2II(i,j,k) * factor
        HeII(i,j,k) = HeII(i,j,k) * factor
        HNOII(i,j,k) = HNOII(i,j,k) * factor
        H2NOII(i,j,k) = H2NOII(i,j,k) * factor
        H3II(i,j,k) = H3II(i,j,k) * factor
        H3COII(i,j,k) = H3COII(i,j,k) * factor
        H3OII(i,j,k) = H3OII(i,j,k) * factor
        HCNHII(i,j,k) = HCNHII(i,j,k) * factor
        HCO2II(i,j,k) = HCO2II(i,j,k) * factor
        HeHII(i,j,k) = HeHII(i,j,k) * factor
        N2HII(i,j,k) = N2HII(i,j,k) * factor
        O2HII(i,j,k) = O2HII(i,j,k) * factor

        !mimimal value check
        De(i,j,k) = max(De(i,j,k), krome_tiny)
        CHI(i,j,k) = max(CHI(i,j,k), krome_tiny)
        OI(i,j,k) = max(OI(i,j,k), krome_tiny)
        HNCI(i,j,k) = max(HNCI(i,j,k), krome_tiny)
        HCNI(i,j,k) = max(HCNI(i,j,k), krome_tiny)
        H2I(i,j,k) = max(H2I(i,j,k), krome_tiny)
        CI(i,j,k) = max(CI(i,j,k), krome_tiny)
        HI(i,j,k) = max(HI(i,j,k), krome_tiny)
        H2OI(i,j,k) = max(H2OI(i,j,k), krome_tiny)
        OHI(i,j,k) = max(OHI(i,j,k), krome_tiny)
        O2I(i,j,k) = max(O2I(i,j,k), krome_tiny)
        CH2I(i,j,k) = max(CH2I(i,j,k), krome_tiny)
        H2COI(i,j,k) = max(H2COI(i,j,k), krome_tiny)
        HCOI(i,j,k) = max(HCOI(i,j,k), krome_tiny)
        MGI(i,j,k) = max(MGI(i,j,k), krome_tiny)
        NH3I(i,j,k) = max(NH3I(i,j,k), krome_tiny)
        NOI(i,j,k) = max(NOI(i,j,k), krome_tiny)
        CNI(i,j,k) = max(CNI(i,j,k), krome_tiny)
        COI(i,j,k) = max(COI(i,j,k), krome_tiny)
        N2I(i,j,k) = max(N2I(i,j,k), krome_tiny)
        NH2I(i,j,k) = max(NH2I(i,j,k), krome_tiny)
        CH3I(i,j,k) = max(CH3I(i,j,k), krome_tiny)
        CH4I(i,j,k) = max(CH4I(i,j,k), krome_tiny)
        NI(i,j,k) = max(NI(i,j,k), krome_tiny)
        NHI(i,j,k) = max(NHI(i,j,k), krome_tiny)
        HeI(i,j,k) = max(HeI(i,j,k), krome_tiny)
        HNOI(i,j,k) = max(HNOI(i,j,k), krome_tiny)
        CH3OHI(i,j,k) = max(CH3OHI(i,j,k), krome_tiny)
        CO2I(i,j,k) = max(CO2I(i,j,k), krome_tiny)
        H2CNI(i,j,k) = max(H2CNI(i,j,k), krome_tiny)
        HNCOI(i,j,k) = max(HNCOI(i,j,k), krome_tiny)
        NO2I(i,j,k) = max(NO2I(i,j,k), krome_tiny)
        O2HI(i,j,k) = max(O2HI(i,j,k), krome_tiny)
        OCNI(i,j,k) = max(OCNI(i,j,k), krome_tiny)
        CH3OH_DUSTI(i,j,k) = max(CH3OH_DUSTI(i,j,k), krome_tiny)
        HNCO_DUSTI(i,j,k) = max(HNCO_DUSTI(i,j,k), krome_tiny)
        H2CO_DUSTI(i,j,k) = max(H2CO_DUSTI(i,j,k), krome_tiny)
        CH4_DUSTI(i,j,k) = max(CH4_DUSTI(i,j,k), krome_tiny)
        CO_DUSTI(i,j,k) = max(CO_DUSTI(i,j,k), krome_tiny)
        H2O_DUSTI(i,j,k) = max(H2O_DUSTI(i,j,k), krome_tiny)
        NO_DUSTI(i,j,k) = max(NO_DUSTI(i,j,k), krome_tiny)
        CO2_DUSTI(i,j,k) = max(CO2_DUSTI(i,j,k), krome_tiny)
        N2_DUSTI(i,j,k) = max(N2_DUSTI(i,j,k), krome_tiny)
        HCN_DUSTI(i,j,k) = max(HCN_DUSTI(i,j,k), krome_tiny)
        NH3_DUSTI(i,j,k) = max(NH3_DUSTI(i,j,k), krome_tiny)
        O2_DUSTI(i,j,k) = max(O2_DUSTI(i,j,k), krome_tiny)
        NO2_DUSTI(i,j,k) = max(NO2_DUSTI(i,j,k), krome_tiny)
        HNO_DUSTI(i,j,k) = max(HNO_DUSTI(i,j,k), krome_tiny)
        O2H_DUSTI(i,j,k) = max(O2H_DUSTI(i,j,k), krome_tiny)
        H2CN_DUSTI(i,j,k) = max(H2CN_DUSTI(i,j,k), krome_tiny)
        MG_DUSTI(i,j,k) = max(MG_DUSTI(i,j,k), krome_tiny)
        HNC_DUSTI(i,j,k) = max(HNC_DUSTI(i,j,k), krome_tiny)
        E_DUSTI(i,j,k) = max(E_DUSTI(i,j,k), krome_tiny)
        HCOII(i,j,k) = max(HCOII(i,j,k), krome_tiny)
        HII(i,j,k) = max(HII(i,j,k), krome_tiny)
        HOCII(i,j,k) = max(HOCII(i,j,k), krome_tiny)
        CII(i,j,k) = max(CII(i,j,k), krome_tiny)
        CH2II(i,j,k) = max(CH2II(i,j,k), krome_tiny)
        CHII(i,j,k) = max(CHII(i,j,k), krome_tiny)
        H2COII(i,j,k) = max(H2COII(i,j,k), krome_tiny)
        MGII(i,j,k) = max(MGII(i,j,k), krome_tiny)
        NH3II(i,j,k) = max(NH3II(i,j,k), krome_tiny)
        NOII(i,j,k) = max(NOII(i,j,k), krome_tiny)
        CNII(i,j,k) = max(CNII(i,j,k), krome_tiny)
        COII(i,j,k) = max(COII(i,j,k), krome_tiny)
        N2II(i,j,k) = max(N2II(i,j,k), krome_tiny)
        O2II(i,j,k) = max(O2II(i,j,k), krome_tiny)
        H2OII(i,j,k) = max(H2OII(i,j,k), krome_tiny)
        NH2II(i,j,k) = max(NH2II(i,j,k), krome_tiny)
        OII(i,j,k) = max(OII(i,j,k), krome_tiny)
        OHII(i,j,k) = max(OHII(i,j,k), krome_tiny)
        CH3II(i,j,k) = max(CH3II(i,j,k), krome_tiny)
        CH4II(i,j,k) = max(CH4II(i,j,k), krome_tiny)
        NII(i,j,k) = max(NII(i,j,k), krome_tiny)
        HCNII(i,j,k) = max(HCNII(i,j,k), krome_tiny)
        NHII(i,j,k) = max(NHII(i,j,k), krome_tiny)
        H2II(i,j,k) = max(H2II(i,j,k), krome_tiny)
        HeII(i,j,k) = max(HeII(i,j,k), krome_tiny)
        HNOII(i,j,k) = max(HNOII(i,j,k), krome_tiny)
        H2NOII(i,j,k) = max(H2NOII(i,j,k), krome_tiny)
        H3II(i,j,k) = max(H3II(i,j,k), krome_tiny)
        H3COII(i,j,k) = max(H3COII(i,j,k), krome_tiny)
        H3OII(i,j,k) = max(H3OII(i,j,k), krome_tiny)
        HCNHII(i,j,k) = max(HCNHII(i,j,k), krome_tiny)
        HCO2II(i,j,k) = max(HCO2II(i,j,k), krome_tiny)
        HeHII(i,j,k) = max(HeHII(i,j,k), krome_tiny)
        N2HII(i,j,k) = max(N2HII(i,j,k), krome_tiny)
        O2HII(i,j,k) = max(O2HII(i,j,k), krome_tiny)

      end do
    end do
  end do

  !loop over zones
  do k = ks+1, ke+1
    do j = js+1, je+1
      do i = is+1, ie+1

        !rhogas = #KROME_sum to be removed

        !convert to number densities
        krome_x(krome_idx_E) = De(i,j,k) * dom
        krome_x(krome_idx_CH) = CHI(i,j,k) * dom * 0.0769230769231d0
        krome_x(krome_idx_O) = OI(i,j,k) * dom * 0.0625d0
        krome_x(krome_idx_HNC) = HNCI(i,j,k) * dom * 0.037037037037d0
        krome_x(krome_idx_HCN) = HCNI(i,j,k) * dom * 0.037037037037d0
        krome_x(krome_idx_H2) = H2I(i,j,k) * dom * 0.5d0
        krome_x(krome_idx_C) = CI(i,j,k) * dom * 0.0833333333333d0
        krome_x(krome_idx_H) = HI(i,j,k) * dom
        krome_x(krome_idx_H2O) = H2OI(i,j,k) * dom * 0.0555555555556d0
        krome_x(krome_idx_OH) = OHI(i,j,k) * dom * 0.0588235294118d0
        krome_x(krome_idx_O2) = O2I(i,j,k) * dom * 0.03125d0
        krome_x(krome_idx_CH2) = CH2I(i,j,k) * dom * 0.0714285714286d0
        krome_x(krome_idx_H2CO) = H2COI(i,j,k) * dom * 0.0333333333333d0
        krome_x(krome_idx_HCO) = HCOI(i,j,k) * dom * 0.0344827586207d0
        krome_x(krome_idx_MG) = MGI(i,j,k) * dom * 0.0416666666667d0
        krome_x(krome_idx_NH3) = NH3I(i,j,k) * dom * 0.0588235294118d0
        krome_x(krome_idx_NO) = NOI(i,j,k) * dom * 0.0333333333333d0
        krome_x(krome_idx_CN) = CNI(i,j,k) * dom * 0.0384615384615d0
        krome_x(krome_idx_CO) = COI(i,j,k) * dom * 0.0357142857143d0
        krome_x(krome_idx_N2) = N2I(i,j,k) * dom * 0.0357142857143d0
        krome_x(krome_idx_NH2) = NH2I(i,j,k) * dom * 0.0625d0
        krome_x(krome_idx_CH3) = CH3I(i,j,k) * dom * 0.0666666666667d0
        krome_x(krome_idx_CH4) = CH4I(i,j,k) * dom * 0.0625d0
        krome_x(krome_idx_N) = NI(i,j,k) * dom * 0.0714285714286d0
        krome_x(krome_idx_NH) = NHI(i,j,k) * dom * 0.0666666666667d0
        krome_x(krome_idx_HE) = HeI(i,j,k) * dom * 0.25d0
        krome_x(krome_idx_HNO) = HNOI(i,j,k) * dom * 0.0322580645161d0
        krome_x(krome_idx_CH3OH) = CH3OHI(i,j,k) * dom * 0.03125d0
        krome_x(krome_idx_CO2) = CO2I(i,j,k) * dom * 0.0227272727273d0
        krome_x(krome_idx_H2CN) = H2CNI(i,j,k) * dom * 0.0357142857143d0
        krome_x(krome_idx_HNCO) = HNCOI(i,j,k) * dom * 0.0232558139535d0
        krome_x(krome_idx_NO2) = NO2I(i,j,k) * dom * 0.0217391304348d0
        krome_x(krome_idx_O2H) = O2HI(i,j,k) * dom * 0.030303030303d0
        krome_x(krome_idx_OCN) = OCNI(i,j,k) * dom * 0.0238095238095d0
        krome_x(krome_idx_CH3OH_DUST) = CH3OH_DUSTI(i,j,k) * dom * 0.03125d0
        krome_x(krome_idx_HNCO_DUST) = HNCO_DUSTI(i,j,k) * dom * 0.0232558139535d0
        krome_x(krome_idx_H2CO_DUST) = H2CO_DUSTI(i,j,k) * dom * 0.0333333333333d0
        krome_x(krome_idx_CH4_DUST) = CH4_DUSTI(i,j,k) * dom * 0.0625d0
        krome_x(krome_idx_CO_DUST) = CO_DUSTI(i,j,k) * dom * 0.0357142857143d0
        krome_x(krome_idx_H2O_DUST) = H2O_DUSTI(i,j,k) * dom * 0.0555555555556d0
        krome_x(krome_idx_NO_DUST) = NO_DUSTI(i,j,k) * dom * 0.0333333333333d0
        krome_x(krome_idx_CO2_DUST) = CO2_DUSTI(i,j,k) * dom * 0.0227272727273d0
        krome_x(krome_idx_N2_DUST) = N2_DUSTI(i,j,k) * dom * 0.0357142857143d0
        krome_x(krome_idx_HCN_DUST) = HCN_DUSTI(i,j,k) * dom * 0.037037037037d0
        krome_x(krome_idx_NH3_DUST) = NH3_DUSTI(i,j,k) * dom * 0.0588235294118d0
        krome_x(krome_idx_O2_DUST) = O2_DUSTI(i,j,k) * dom * 0.03125d0
        krome_x(krome_idx_NO2_DUST) = NO2_DUSTI(i,j,k) * dom * 0.0217391304348d0
        krome_x(krome_idx_HNO_DUST) = HNO_DUSTI(i,j,k) * dom * 0.0322580645161d0
        krome_x(krome_idx_O2H_DUST) = O2H_DUSTI(i,j,k) * dom * 0.030303030303d0
        krome_x(krome_idx_H2CN_DUST) = H2CN_DUSTI(i,j,k) * dom * 0.0357142857143d0
        krome_x(krome_idx_MG_DUST) = MG_DUSTI(i,j,k) * dom * 0.0416666666667d0
        krome_x(krome_idx_HNC_DUST) = HNC_DUSTI(i,j,k) * dom * 0.037037037037d0
        krome_x(krome_idx_E_DUST) = E_DUSTI(i,j,k) * dom
        krome_x(krome_idx_HCOj) = HCOII(i,j,k) * dom * 0.0344827586207d0
        krome_x(krome_idx_Hj) = HII(i,j,k) * dom
        krome_x(krome_idx_HOCj) = HOCII(i,j,k) * dom * 0.0344827586207d0
        krome_x(krome_idx_Cj) = CII(i,j,k) * dom * 0.0833333333333d0
        krome_x(krome_idx_CH2j) = CH2II(i,j,k) * dom * 0.0714285714286d0
        krome_x(krome_idx_CHj) = CHII(i,j,k) * dom * 0.0769230769231d0
        krome_x(krome_idx_H2COj) = H2COII(i,j,k) * dom * 0.0333333333333d0
        krome_x(krome_idx_MGj) = MGII(i,j,k) * dom * 0.0416666666667d0
        krome_x(krome_idx_NH3j) = NH3II(i,j,k) * dom * 0.0588235294118d0
        krome_x(krome_idx_NOj) = NOII(i,j,k) * dom * 0.0333333333333d0
        krome_x(krome_idx_CNj) = CNII(i,j,k) * dom * 0.0384615384615d0
        krome_x(krome_idx_COj) = COII(i,j,k) * dom * 0.0357142857143d0
        krome_x(krome_idx_N2j) = N2II(i,j,k) * dom * 0.0357142857143d0
        krome_x(krome_idx_O2j) = O2II(i,j,k) * dom * 0.03125d0
        krome_x(krome_idx_H2Oj) = H2OII(i,j,k) * dom * 0.0555555555556d0
        krome_x(krome_idx_NH2j) = NH2II(i,j,k) * dom * 0.0625d0
        krome_x(krome_idx_Oj) = OII(i,j,k) * dom * 0.0625d0
        krome_x(krome_idx_OHj) = OHII(i,j,k) * dom * 0.0588235294118d0
        krome_x(krome_idx_CH3j) = CH3II(i,j,k) * dom * 0.0666666666667d0
        krome_x(krome_idx_CH4j) = CH4II(i,j,k) * dom * 0.0625d0
        krome_x(krome_idx_Nj) = NII(i,j,k) * dom * 0.0714285714286d0
        krome_x(krome_idx_HCNj) = HCNII(i,j,k) * dom * 0.037037037037d0
        krome_x(krome_idx_NHj) = NHII(i,j,k) * dom * 0.0666666666667d0
        krome_x(krome_idx_H2j) = H2II(i,j,k) * dom * 0.5d0
        krome_x(krome_idx_HEj) = HeII(i,j,k) * dom * 0.25d0
        krome_x(krome_idx_HNOj) = HNOII(i,j,k) * dom * 0.0322580645161d0
        krome_x(krome_idx_H2NOj) = H2NOII(i,j,k) * dom * 0.03125d0
        krome_x(krome_idx_H3j) = H3II(i,j,k) * dom * 0.333333333333d0
        krome_x(krome_idx_H3COj) = H3COII(i,j,k) * dom * 0.0322580645161d0
        krome_x(krome_idx_H3Oj) = H3OII(i,j,k) * dom * 0.0526315789474d0
        krome_x(krome_idx_HCNHj) = HCNHII(i,j,k) * dom * 0.0357142857143d0
        krome_x(krome_idx_HCO2j) = HCO2II(i,j,k) * dom * 0.0222222222222d0
        krome_x(krome_idx_HEHj) = HeHII(i,j,k) * dom * 0.2d0
        krome_x(krome_idx_N2Hj) = N2HII(i,j,k) * dom * 0.0344827586207d0
        krome_x(krome_idx_O2Hj) = O2HII(i,j,k) * dom * 0.030303030303d0

        call evaluate_tgas(d(i,j,k), e(i,j,k), ge(i,j,k),&
            u(i,j,k), v(i,j,k), w(i,j,k),&
            krome_x(:),imethod,idual,idim,tgas,&
            utem)

        !store old tgas
        tgasold = tgas

        dt_hydro = utim*dt !dt*time_conversion

        !call KROME solver
        call krome(krome_x(:),tgas,dt_hydro)

        idom = 1.d0/dom
        !convert back to code units
        De(i,j,k) = krome_x(krome_idx_E) * idom
        CHI(i,j,k) = krome_x(krome_idx_CH) * idom * 13d0
        OI(i,j,k) = krome_x(krome_idx_O) * idom * 16d0
        HNCI(i,j,k) = krome_x(krome_idx_HNC) * idom * 27d0
        HCNI(i,j,k) = krome_x(krome_idx_HCN) * idom * 27d0
        H2I(i,j,k) = krome_x(krome_idx_H2) * idom * 2d0
        CI(i,j,k) = krome_x(krome_idx_C) * idom * 12d0
        HI(i,j,k) = krome_x(krome_idx_H) * idom
        H2OI(i,j,k) = krome_x(krome_idx_H2O) * idom * 18d0
        OHI(i,j,k) = krome_x(krome_idx_OH) * idom * 17d0
        O2I(i,j,k) = krome_x(krome_idx_O2) * idom * 32d0
        CH2I(i,j,k) = krome_x(krome_idx_CH2) * idom * 14d0
        H2COI(i,j,k) = krome_x(krome_idx_H2CO) * idom * 30d0
        HCOI(i,j,k) = krome_x(krome_idx_HCO) * idom * 29d0
        MGI(i,j,k) = krome_x(krome_idx_MG) * idom * 24d0
        NH3I(i,j,k) = krome_x(krome_idx_NH3) * idom * 17d0
        NOI(i,j,k) = krome_x(krome_idx_NO) * idom * 30d0
        CNI(i,j,k) = krome_x(krome_idx_CN) * idom * 26d0
        COI(i,j,k) = krome_x(krome_idx_CO) * idom * 28d0
        N2I(i,j,k) = krome_x(krome_idx_N2) * idom * 28d0
        NH2I(i,j,k) = krome_x(krome_idx_NH2) * idom * 16d0
        CH3I(i,j,k) = krome_x(krome_idx_CH3) * idom * 15d0
        CH4I(i,j,k) = krome_x(krome_idx_CH4) * idom * 16d0
        NI(i,j,k) = krome_x(krome_idx_N) * idom * 14d0
        NHI(i,j,k) = krome_x(krome_idx_NH) * idom * 15d0
        HeI(i,j,k) = krome_x(krome_idx_HE) * idom * 4d0
        HNOI(i,j,k) = krome_x(krome_idx_HNO) * idom * 31d0
        CH3OHI(i,j,k) = krome_x(krome_idx_CH3OH) * idom * 32d0
        CO2I(i,j,k) = krome_x(krome_idx_CO2) * idom * 44d0
        H2CNI(i,j,k) = krome_x(krome_idx_H2CN) * idom * 28d0
        HNCOI(i,j,k) = krome_x(krome_idx_HNCO) * idom * 43d0
        NO2I(i,j,k) = krome_x(krome_idx_NO2) * idom * 46d0
        O2HI(i,j,k) = krome_x(krome_idx_O2H) * idom * 33d0
        OCNI(i,j,k) = krome_x(krome_idx_OCN) * idom * 42d0
        CH3OH_DUSTI(i,j,k) = krome_x(krome_idx_CH3OH_DUST) * idom * 32d0
        HNCO_DUSTI(i,j,k) = krome_x(krome_idx_HNCO_DUST) * idom * 43d0
        H2CO_DUSTI(i,j,k) = krome_x(krome_idx_H2CO_DUST) * idom * 30d0
        CH4_DUSTI(i,j,k) = krome_x(krome_idx_CH4_DUST) * idom * 16d0
        CO_DUSTI(i,j,k) = krome_x(krome_idx_CO_DUST) * idom * 28d0
        H2O_DUSTI(i,j,k) = krome_x(krome_idx_H2O_DUST) * idom * 18d0
        NO_DUSTI(i,j,k) = krome_x(krome_idx_NO_DUST) * idom * 30d0
        CO2_DUSTI(i,j,k) = krome_x(krome_idx_CO2_DUST) * idom * 44d0
        N2_DUSTI(i,j,k) = krome_x(krome_idx_N2_DUST) * idom * 28d0
        HCN_DUSTI(i,j,k) = krome_x(krome_idx_HCN_DUST) * idom * 27d0
        NH3_DUSTI(i,j,k) = krome_x(krome_idx_NH3_DUST) * idom * 17d0
        O2_DUSTI(i,j,k) = krome_x(krome_idx_O2_DUST) * idom * 32d0
        NO2_DUSTI(i,j,k) = krome_x(krome_idx_NO2_DUST) * idom * 46d0
        HNO_DUSTI(i,j,k) = krome_x(krome_idx_HNO_DUST) * idom * 31d0
        O2H_DUSTI(i,j,k) = krome_x(krome_idx_O2H_DUST) * idom * 33d0
        H2CN_DUSTI(i,j,k) = krome_x(krome_idx_H2CN_DUST) * idom * 28d0
        MG_DUSTI(i,j,k) = krome_x(krome_idx_MG_DUST) * idom * 24d0
        HNC_DUSTI(i,j,k) = krome_x(krome_idx_HNC_DUST) * idom * 27d0
        E_DUSTI(i,j,k) = krome_x(krome_idx_E_DUST) * idom
        HCOII(i,j,k) = krome_x(krome_idx_HCOj) * idom * 29d0
        HII(i,j,k) = krome_x(krome_idx_Hj) * idom
        HOCII(i,j,k) = krome_x(krome_idx_HOCj) * idom * 29d0
        CII(i,j,k) = krome_x(krome_idx_Cj) * idom * 12d0
        CH2II(i,j,k) = krome_x(krome_idx_CH2j) * idom * 14d0
        CHII(i,j,k) = krome_x(krome_idx_CHj) * idom * 13d0
        H2COII(i,j,k) = krome_x(krome_idx_H2COj) * idom * 30d0
        MGII(i,j,k) = krome_x(krome_idx_MGj) * idom * 24d0
        NH3II(i,j,k) = krome_x(krome_idx_NH3j) * idom * 17d0
        NOII(i,j,k) = krome_x(krome_idx_NOj) * idom * 30d0
        CNII(i,j,k) = krome_x(krome_idx_CNj) * idom * 26d0
        COII(i,j,k) = krome_x(krome_idx_COj) * idom * 28d0
        N2II(i,j,k) = krome_x(krome_idx_N2j) * idom * 28d0
        O2II(i,j,k) = krome_x(krome_idx_O2j) * idom * 32d0
        H2OII(i,j,k) = krome_x(krome_idx_H2Oj) * idom * 18d0
        NH2II(i,j,k) = krome_x(krome_idx_NH2j) * idom * 16d0
        OII(i,j,k) = krome_x(krome_idx_Oj) * idom * 16d0
        OHII(i,j,k) = krome_x(krome_idx_OHj) * idom * 17d0
        CH3II(i,j,k) = krome_x(krome_idx_CH3j) * idom * 15d0
        CH4II(i,j,k) = krome_x(krome_idx_CH4j) * idom * 16d0
        NII(i,j,k) = krome_x(krome_idx_Nj) * idom * 14d0
        HCNII(i,j,k) = krome_x(krome_idx_HCNj) * idom * 27d0
        NHII(i,j,k) = krome_x(krome_idx_NHj) * idom * 15d0
        H2II(i,j,k) = krome_x(krome_idx_H2j) * idom * 2d0
        HeII(i,j,k) = krome_x(krome_idx_HEj) * idom * 4d0
        HNOII(i,j,k) = krome_x(krome_idx_HNOj) * idom * 31d0
        H2NOII(i,j,k) = krome_x(krome_idx_H2NOj) * idom * 32d0
        H3II(i,j,k) = krome_x(krome_idx_H3j) * idom * 3d0
        H3COII(i,j,k) = krome_x(krome_idx_H3COj) * idom * 31d0
        H3OII(i,j,k) = krome_x(krome_idx_H3Oj) * idom * 19d0
        HCNHII(i,j,k) = krome_x(krome_idx_HCNHj) * idom * 28d0
        HCO2II(i,j,k) = krome_x(krome_idx_HCO2j) * idom * 45d0
        HeHII(i,j,k) = krome_x(krome_idx_HEHj) * idom * 5d0
        N2HII(i,j,k) = krome_x(krome_idx_N2Hj) * idom * 29d0
        O2HII(i,j,k) = krome_x(krome_idx_O2Hj) * idom * 33d0

        !evaluate energy from temperature difference
        edot = (tgas - tgasold) * d(i,j,k) &
            / ((gamma - 1.d0) * utem * dt)

        !update internal energy
        e(i,j,k)  = e(i,j,k) + edot / d(i,j,k) * dt
        !when using dual
        if (idual .eq. 1) ge(i,j,k) = ge(i,j,k)+ edot / d(i,j,k) * dt
      end do
    end do
  end do

  !scale comoving<-proper
  factor = aye**3
  do k = ks+1, ke+1
    do j = js+1, je+1
      do i = is+1, ie+1
        d(i,j,k) = d(i,j,k) * factor
        !scale comoving->proper
        De(i,j,k) = De(i,j,k) * factor
        CHI(i,j,k) = CHI(i,j,k) * factor
        OI(i,j,k) = OI(i,j,k) * factor
        HNCI(i,j,k) = HNCI(i,j,k) * factor
        HCNI(i,j,k) = HCNI(i,j,k) * factor
        H2I(i,j,k) = H2I(i,j,k) * factor
        CI(i,j,k) = CI(i,j,k) * factor
        HI(i,j,k) = HI(i,j,k) * factor
        H2OI(i,j,k) = H2OI(i,j,k) * factor
        OHI(i,j,k) = OHI(i,j,k) * factor
        O2I(i,j,k) = O2I(i,j,k) * factor
        CH2I(i,j,k) = CH2I(i,j,k) * factor
        H2COI(i,j,k) = H2COI(i,j,k) * factor
        HCOI(i,j,k) = HCOI(i,j,k) * factor
        MGI(i,j,k) = MGI(i,j,k) * factor
        NH3I(i,j,k) = NH3I(i,j,k) * factor
        NOI(i,j,k) = NOI(i,j,k) * factor
        CNI(i,j,k) = CNI(i,j,k) * factor
        COI(i,j,k) = COI(i,j,k) * factor
        N2I(i,j,k) = N2I(i,j,k) * factor
        NH2I(i,j,k) = NH2I(i,j,k) * factor
        CH3I(i,j,k) = CH3I(i,j,k) * factor
        CH4I(i,j,k) = CH4I(i,j,k) * factor
        NI(i,j,k) = NI(i,j,k) * factor
        NHI(i,j,k) = NHI(i,j,k) * factor
        HeI(i,j,k) = HeI(i,j,k) * factor
        HNOI(i,j,k) = HNOI(i,j,k) * factor
        CH3OHI(i,j,k) = CH3OHI(i,j,k) * factor
        CO2I(i,j,k) = CO2I(i,j,k) * factor
        H2CNI(i,j,k) = H2CNI(i,j,k) * factor
        HNCOI(i,j,k) = HNCOI(i,j,k) * factor
        NO2I(i,j,k) = NO2I(i,j,k) * factor
        O2HI(i,j,k) = O2HI(i,j,k) * factor
        OCNI(i,j,k) = OCNI(i,j,k) * factor
        CH3OH_DUSTI(i,j,k) = CH3OH_DUSTI(i,j,k) * factor
        HNCO_DUSTI(i,j,k) = HNCO_DUSTI(i,j,k) * factor
        H2CO_DUSTI(i,j,k) = H2CO_DUSTI(i,j,k) * factor
        CH4_DUSTI(i,j,k) = CH4_DUSTI(i,j,k) * factor
        CO_DUSTI(i,j,k) = CO_DUSTI(i,j,k) * factor
        H2O_DUSTI(i,j,k) = H2O_DUSTI(i,j,k) * factor
        NO_DUSTI(i,j,k) = NO_DUSTI(i,j,k) * factor
        CO2_DUSTI(i,j,k) = CO2_DUSTI(i,j,k) * factor
        N2_DUSTI(i,j,k) = N2_DUSTI(i,j,k) * factor
        HCN_DUSTI(i,j,k) = HCN_DUSTI(i,j,k) * factor
        NH3_DUSTI(i,j,k) = NH3_DUSTI(i,j,k) * factor
        O2_DUSTI(i,j,k) = O2_DUSTI(i,j,k) * factor
        NO2_DUSTI(i,j,k) = NO2_DUSTI(i,j,k) * factor
        HNO_DUSTI(i,j,k) = HNO_DUSTI(i,j,k) * factor
        O2H_DUSTI(i,j,k) = O2H_DUSTI(i,j,k) * factor
        H2CN_DUSTI(i,j,k) = H2CN_DUSTI(i,j,k) * factor
        MG_DUSTI(i,j,k) = MG_DUSTI(i,j,k) * factor
        HNC_DUSTI(i,j,k) = HNC_DUSTI(i,j,k) * factor
        E_DUSTI(i,j,k) = E_DUSTI(i,j,k) * factor
        HCOII(i,j,k) = HCOII(i,j,k) * factor
        HII(i,j,k) = HII(i,j,k) * factor
        HOCII(i,j,k) = HOCII(i,j,k) * factor
        CII(i,j,k) = CII(i,j,k) * factor
        CH2II(i,j,k) = CH2II(i,j,k) * factor
        CHII(i,j,k) = CHII(i,j,k) * factor
        H2COII(i,j,k) = H2COII(i,j,k) * factor
        MGII(i,j,k) = MGII(i,j,k) * factor
        NH3II(i,j,k) = NH3II(i,j,k) * factor
        NOII(i,j,k) = NOII(i,j,k) * factor
        CNII(i,j,k) = CNII(i,j,k) * factor
        COII(i,j,k) = COII(i,j,k) * factor
        N2II(i,j,k) = N2II(i,j,k) * factor
        O2II(i,j,k) = O2II(i,j,k) * factor
        H2OII(i,j,k) = H2OII(i,j,k) * factor
        NH2II(i,j,k) = NH2II(i,j,k) * factor
        OII(i,j,k) = OII(i,j,k) * factor
        OHII(i,j,k) = OHII(i,j,k) * factor
        CH3II(i,j,k) = CH3II(i,j,k) * factor
        CH4II(i,j,k) = CH4II(i,j,k) * factor
        NII(i,j,k) = NII(i,j,k) * factor
        HCNII(i,j,k) = HCNII(i,j,k) * factor
        NHII(i,j,k) = NHII(i,j,k) * factor
        H2II(i,j,k) = H2II(i,j,k) * factor
        HeII(i,j,k) = HeII(i,j,k) * factor
        HNOII(i,j,k) = HNOII(i,j,k) * factor
        H2NOII(i,j,k) = H2NOII(i,j,k) * factor
        H3II(i,j,k) = H3II(i,j,k) * factor
        H3COII(i,j,k) = H3COII(i,j,k) * factor
        H3OII(i,j,k) = H3OII(i,j,k) * factor
        HCNHII(i,j,k) = HCNHII(i,j,k) * factor
        HCO2II(i,j,k) = HCO2II(i,j,k) * factor
        HeHII(i,j,k) = HeHII(i,j,k) * factor
        N2HII(i,j,k) = N2HII(i,j,k) * factor
        O2HII(i,j,k) = O2HII(i,j,k) * factor

        !mimimal value check
        De(i,j,k) = max(De(i,j,k), krome_tiny)
        CHI(i,j,k) = max(CHI(i,j,k), krome_tiny)
        OI(i,j,k) = max(OI(i,j,k), krome_tiny)
        HNCI(i,j,k) = max(HNCI(i,j,k), krome_tiny)
        HCNI(i,j,k) = max(HCNI(i,j,k), krome_tiny)
        H2I(i,j,k) = max(H2I(i,j,k), krome_tiny)
        CI(i,j,k) = max(CI(i,j,k), krome_tiny)
        HI(i,j,k) = max(HI(i,j,k), krome_tiny)
        H2OI(i,j,k) = max(H2OI(i,j,k), krome_tiny)
        OHI(i,j,k) = max(OHI(i,j,k), krome_tiny)
        O2I(i,j,k) = max(O2I(i,j,k), krome_tiny)
        CH2I(i,j,k) = max(CH2I(i,j,k), krome_tiny)
        H2COI(i,j,k) = max(H2COI(i,j,k), krome_tiny)
        HCOI(i,j,k) = max(HCOI(i,j,k), krome_tiny)
        MGI(i,j,k) = max(MGI(i,j,k), krome_tiny)
        NH3I(i,j,k) = max(NH3I(i,j,k), krome_tiny)
        NOI(i,j,k) = max(NOI(i,j,k), krome_tiny)
        CNI(i,j,k) = max(CNI(i,j,k), krome_tiny)
        COI(i,j,k) = max(COI(i,j,k), krome_tiny)
        N2I(i,j,k) = max(N2I(i,j,k), krome_tiny)
        NH2I(i,j,k) = max(NH2I(i,j,k), krome_tiny)
        CH3I(i,j,k) = max(CH3I(i,j,k), krome_tiny)
        CH4I(i,j,k) = max(CH4I(i,j,k), krome_tiny)
        NI(i,j,k) = max(NI(i,j,k), krome_tiny)
        NHI(i,j,k) = max(NHI(i,j,k), krome_tiny)
        HeI(i,j,k) = max(HeI(i,j,k), krome_tiny)
        HNOI(i,j,k) = max(HNOI(i,j,k), krome_tiny)
        CH3OHI(i,j,k) = max(CH3OHI(i,j,k), krome_tiny)
        CO2I(i,j,k) = max(CO2I(i,j,k), krome_tiny)
        H2CNI(i,j,k) = max(H2CNI(i,j,k), krome_tiny)
        HNCOI(i,j,k) = max(HNCOI(i,j,k), krome_tiny)
        NO2I(i,j,k) = max(NO2I(i,j,k), krome_tiny)
        O2HI(i,j,k) = max(O2HI(i,j,k), krome_tiny)
        OCNI(i,j,k) = max(OCNI(i,j,k), krome_tiny)
        CH3OH_DUSTI(i,j,k) = max(CH3OH_DUSTI(i,j,k), krome_tiny)
        HNCO_DUSTI(i,j,k) = max(HNCO_DUSTI(i,j,k), krome_tiny)
        H2CO_DUSTI(i,j,k) = max(H2CO_DUSTI(i,j,k), krome_tiny)
        CH4_DUSTI(i,j,k) = max(CH4_DUSTI(i,j,k), krome_tiny)
        CO_DUSTI(i,j,k) = max(CO_DUSTI(i,j,k), krome_tiny)
        H2O_DUSTI(i,j,k) = max(H2O_DUSTI(i,j,k), krome_tiny)
        NO_DUSTI(i,j,k) = max(NO_DUSTI(i,j,k), krome_tiny)
        CO2_DUSTI(i,j,k) = max(CO2_DUSTI(i,j,k), krome_tiny)
        N2_DUSTI(i,j,k) = max(N2_DUSTI(i,j,k), krome_tiny)
        HCN_DUSTI(i,j,k) = max(HCN_DUSTI(i,j,k), krome_tiny)
        NH3_DUSTI(i,j,k) = max(NH3_DUSTI(i,j,k), krome_tiny)
        O2_DUSTI(i,j,k) = max(O2_DUSTI(i,j,k), krome_tiny)
        NO2_DUSTI(i,j,k) = max(NO2_DUSTI(i,j,k), krome_tiny)
        HNO_DUSTI(i,j,k) = max(HNO_DUSTI(i,j,k), krome_tiny)
        O2H_DUSTI(i,j,k) = max(O2H_DUSTI(i,j,k), krome_tiny)
        H2CN_DUSTI(i,j,k) = max(H2CN_DUSTI(i,j,k), krome_tiny)
        MG_DUSTI(i,j,k) = max(MG_DUSTI(i,j,k), krome_tiny)
        HNC_DUSTI(i,j,k) = max(HNC_DUSTI(i,j,k), krome_tiny)
        E_DUSTI(i,j,k) = max(E_DUSTI(i,j,k), krome_tiny)
        HCOII(i,j,k) = max(HCOII(i,j,k), krome_tiny)
        HII(i,j,k) = max(HII(i,j,k), krome_tiny)
        HOCII(i,j,k) = max(HOCII(i,j,k), krome_tiny)
        CII(i,j,k) = max(CII(i,j,k), krome_tiny)
        CH2II(i,j,k) = max(CH2II(i,j,k), krome_tiny)
        CHII(i,j,k) = max(CHII(i,j,k), krome_tiny)
        H2COII(i,j,k) = max(H2COII(i,j,k), krome_tiny)
        MGII(i,j,k) = max(MGII(i,j,k), krome_tiny)
        NH3II(i,j,k) = max(NH3II(i,j,k), krome_tiny)
        NOII(i,j,k) = max(NOII(i,j,k), krome_tiny)
        CNII(i,j,k) = max(CNII(i,j,k), krome_tiny)
        COII(i,j,k) = max(COII(i,j,k), krome_tiny)
        N2II(i,j,k) = max(N2II(i,j,k), krome_tiny)
        O2II(i,j,k) = max(O2II(i,j,k), krome_tiny)
        H2OII(i,j,k) = max(H2OII(i,j,k), krome_tiny)
        NH2II(i,j,k) = max(NH2II(i,j,k), krome_tiny)
        OII(i,j,k) = max(OII(i,j,k), krome_tiny)
        OHII(i,j,k) = max(OHII(i,j,k), krome_tiny)
        CH3II(i,j,k) = max(CH3II(i,j,k), krome_tiny)
        CH4II(i,j,k) = max(CH4II(i,j,k), krome_tiny)
        NII(i,j,k) = max(NII(i,j,k), krome_tiny)
        HCNII(i,j,k) = max(HCNII(i,j,k), krome_tiny)
        NHII(i,j,k) = max(NHII(i,j,k), krome_tiny)
        H2II(i,j,k) = max(H2II(i,j,k), krome_tiny)
        HeII(i,j,k) = max(HeII(i,j,k), krome_tiny)
        HNOII(i,j,k) = max(HNOII(i,j,k), krome_tiny)
        H2NOII(i,j,k) = max(H2NOII(i,j,k), krome_tiny)
        H3II(i,j,k) = max(H3II(i,j,k), krome_tiny)
        H3COII(i,j,k) = max(H3COII(i,j,k), krome_tiny)
        H3OII(i,j,k) = max(H3OII(i,j,k), krome_tiny)
        HCNHII(i,j,k) = max(HCNHII(i,j,k), krome_tiny)
        HCO2II(i,j,k) = max(HCO2II(i,j,k), krome_tiny)
        HeHII(i,j,k) = max(HeHII(i,j,k), krome_tiny)
        N2HII(i,j,k) = max(N2HII(i,j,k), krome_tiny)
        O2HII(i,j,k) = max(O2HII(i,j,k), krome_tiny)

      end do
    end do
  end do

end subroutine krome_driver
