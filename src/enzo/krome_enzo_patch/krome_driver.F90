!KROME_DRIVER

subroutine krome_driver(d, e, ge, u, v, w, &
      De, CHI, OI, HNCI, &
      HCNI, H2I, CI, HI, &
      H2OI, OHI, O2I, CH2I, &
      H2COI, HCOI, MGI, NH3I, &
      NOI, SII, SIC2I, SIC3I, &
      SICI, SIH2I, SIH3I, CNI, &
      COI, N2I, NH2I, CH3I, &
      CH4I, NI, NHI, SIH4I, &
      SIHI, SIOI, HeI, HNOI, &
      CH3OHI, CO2I, H2CNI, H2SIOI, &
      HNCOI, NO2I, O2HI, OCNI, &
      CH3OH_DUSTI, HNCO_DUSTI, H2CO_DUSTI, SIH4_DUSTI, &
      H2SIO_DUSTI, SIC_DUSTI, SIC2_DUSTI, SIC3_DUSTI, &
      CH4_DUSTI, CO_DUSTI, H2O_DUSTI, NO_DUSTI, &
      CO2_DUSTI, N2_DUSTI, HCN_DUSTI, NH3_DUSTI, &
      O2_DUSTI, NO2_DUSTI, HNO_DUSTI, O2H_DUSTI, &
      H2CN_DUSTI, MG_DUSTI, HNC_DUSTI, E_DUSTI, &
      SIO_DUSTI, HCOII, HII, HOCII, &
      CII, CH2II, CHII, H2COII, &
      MGII, NH3II, NOII, SIII, &
      SIC2II, SIC3II, SICII, SIH2II, &
      SIH3II, CNII, COII, N2II, &
      O2II, H2OII, NH2II, OII, &
      OHII, CH3II, CH4II, NII, &
      HCNII, NHII, SIH4II, SIHII, &
      SIOII, H2II, HeII, HNOII, &
      H2NOII, H3II, H3COII, H3OII, &
      HCNHII, HCO2II, HeHII, N2HII, &
      O2HII, SIH5II, SIOHII, &
      in, jn, kn, imethod, &
      idual, idim, &
      is, js, ks, ie, je, ke, &
      dt, aye, &
      utem, uxyz, uaye, urho, utim, &
      gamma, fh, dtoh, cellsize, opt_spu)

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
  use krome_user_commons
  use krome_constants

  implicit none

  integer :: projectiles(6)=(/krome_idx_H2, krome_idx_HE, krome_idx_C, krome_idx_O, krome_idx_SI, krome_idx_CO/)
  ! real*8, parameter :: grainrad = 1d-5

  real*8,parameter::mh=p_mass !mass_h
  real*8::dom,factor,tgas,tgasold,krome_x(krome_nmols)
  real*8::dt,dt_hydro,idom,edot,krome_tiny
  real*8::d(in,jn,kn),e(in,jn,kn),ge(in,jn,kn)
  real*8::u(in,jn,kn),v(in,jn,kn),w(in,jn,kn)
  real*8::numdens(in,jn,kn)
  real*8::ul,ur,vl,vr,wl,wr
  real*8::ulyields,uryields,vlyields,vryields,wlyields,wryields
  real*8::ulflux,urflux,vlflux,vrflux,wlflux,wrflux
  real*8::sputtering
  real*8::aye,utem,uxyz,uaye,urho,utim,gamma,fh,dtoh
  real*8::cellsize,vel,uvel
  integer::in,jn,kn,imethod,idual,is,js,ks,ie,je,ke,idim
  integer::i,j,k
  integer::projidx, specidx
  integer::opt_spu

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
  real*8::SII(in,jn,kn)
  real*8::SIC2I(in,jn,kn)
  real*8::SIC3I(in,jn,kn)
  real*8::SICI(in,jn,kn)
  real*8::SIH2I(in,jn,kn)
  real*8::SIH3I(in,jn,kn)
  real*8::CNI(in,jn,kn)
  real*8::COI(in,jn,kn)
  real*8::N2I(in,jn,kn)
  real*8::NH2I(in,jn,kn)
  real*8::CH3I(in,jn,kn)
  real*8::CH4I(in,jn,kn)
  real*8::NI(in,jn,kn)
  real*8::NHI(in,jn,kn)
  real*8::SIH4I(in,jn,kn)
  real*8::SIHI(in,jn,kn)
  real*8::SIOI(in,jn,kn)
  real*8::HeI(in,jn,kn)
  real*8::HNOI(in,jn,kn)
  real*8::CH3OHI(in,jn,kn)
  real*8::CO2I(in,jn,kn)
  real*8::H2CNI(in,jn,kn)
  real*8::H2SIOI(in,jn,kn)
  real*8::HNCOI(in,jn,kn)
  real*8::NO2I(in,jn,kn)
  real*8::O2HI(in,jn,kn)
  real*8::OCNI(in,jn,kn)
  real*8::CH3OH_DUSTI(in,jn,kn)
  real*8::HNCO_DUSTI(in,jn,kn)
  real*8::H2CO_DUSTI(in,jn,kn)
  real*8::SIH4_DUSTI(in,jn,kn)
  real*8::H2SIO_DUSTI(in,jn,kn)
  real*8::SIC_DUSTI(in,jn,kn)
  real*8::SIC2_DUSTI(in,jn,kn)
  real*8::SIC3_DUSTI(in,jn,kn)
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
  real*8::SIO_DUSTI(in,jn,kn)
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
  real*8::SIII(in,jn,kn)
  real*8::SIC2II(in,jn,kn)
  real*8::SIC3II(in,jn,kn)
  real*8::SICII(in,jn,kn)
  real*8::SIH2II(in,jn,kn)
  real*8::SIH3II(in,jn,kn)
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
  real*8::SIH4II(in,jn,kn)
  real*8::SIHII(in,jn,kn)
  real*8::SIOII(in,jn,kn)
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
  real*8::SIH5II(in,jn,kn)
  real*8::SIOHII(in,jn,kn)

  !******************************

  !set units
  dom = urho*(aye**3)/mh
  uvel = uxyz/utim

  !scaling factor for comoving->proper
  factor = aye**(-3)

  !check minimal value and comoving->proper
  d = d * factor
  De = max(De * factor, krome_tiny)
  CHI = max(CHI * factor, krome_tiny)
  OI = max(OI * factor, krome_tiny)
  HNCI = max(HNCI * factor, krome_tiny)
  HCNI = max(HCNI * factor, krome_tiny)
  H2I = max(H2I * factor, krome_tiny)
  CI = max(CI * factor, krome_tiny)
  HI = max(HI * factor, krome_tiny)
  H2OI = max(H2OI * factor, krome_tiny)
  OHI = max(OHI * factor, krome_tiny)
  O2I = max(O2I * factor, krome_tiny)
  CH2I = max(CH2I * factor, krome_tiny)
  H2COI = max(H2COI * factor, krome_tiny)
  HCOI = max(HCOI * factor, krome_tiny)
  MGI = max(MGI * factor, krome_tiny)
  NH3I = max(NH3I * factor, krome_tiny)
  NOI = max(NOI * factor, krome_tiny)
  SII = max(SII * factor, krome_tiny)
  SIC2I = max(SIC2I * factor, krome_tiny)
  SIC3I = max(SIC3I * factor, krome_tiny)
  SICI = max(SICI * factor, krome_tiny)
  SIH2I = max(SIH2I * factor, krome_tiny)
  SIH3I = max(SIH3I * factor, krome_tiny)
  CNI = max(CNI * factor, krome_tiny)
  COI = max(COI * factor, krome_tiny)
  N2I = max(N2I * factor, krome_tiny)
  NH2I = max(NH2I * factor, krome_tiny)
  CH3I = max(CH3I * factor, krome_tiny)
  CH4I = max(CH4I * factor, krome_tiny)
  NI = max(NI * factor, krome_tiny)
  NHI = max(NHI * factor, krome_tiny)
  SIH4I = max(SIH4I * factor, krome_tiny)
  SIHI = max(SIHI * factor, krome_tiny)
  SIOI = max(SIOI * factor, krome_tiny)
  HeI = max(HeI * factor, krome_tiny)
  HNOI = max(HNOI * factor, krome_tiny)
  CH3OHI = max(CH3OHI * factor, krome_tiny)
  CO2I = max(CO2I * factor, krome_tiny)
  H2CNI = max(H2CNI * factor, krome_tiny)
  H2SIOI = max(H2SIOI * factor, krome_tiny)
  HNCOI = max(HNCOI * factor, krome_tiny)
  NO2I = max(NO2I * factor, krome_tiny)
  O2HI = max(O2HI * factor, krome_tiny)
  OCNI = max(OCNI * factor, krome_tiny)
  CH3OH_DUSTI = max(CH3OH_DUSTI * factor, krome_tiny)
  HNCO_DUSTI = max(HNCO_DUSTI * factor, krome_tiny)
  H2CO_DUSTI = max(H2CO_DUSTI * factor, krome_tiny)
  SIH4_DUSTI = max(SIH4_DUSTI * factor, krome_tiny)
  H2SIO_DUSTI = max(H2SIO_DUSTI * factor, krome_tiny)
  SIC_DUSTI = max(SIC_DUSTI * factor, krome_tiny)
  SIC2_DUSTI = max(SIC2_DUSTI * factor, krome_tiny)
  SIC3_DUSTI = max(SIC3_DUSTI * factor, krome_tiny)
  CH4_DUSTI = max(CH4_DUSTI * factor, krome_tiny)
  CO_DUSTI = max(CO_DUSTI * factor, krome_tiny)
  H2O_DUSTI = max(H2O_DUSTI * factor, krome_tiny)
  NO_DUSTI = max(NO_DUSTI * factor, krome_tiny)
  CO2_DUSTI = max(CO2_DUSTI * factor, krome_tiny)
  N2_DUSTI = max(N2_DUSTI * factor, krome_tiny)
  HCN_DUSTI = max(HCN_DUSTI * factor, krome_tiny)
  NH3_DUSTI = max(NH3_DUSTI * factor, krome_tiny)
  O2_DUSTI = max(O2_DUSTI * factor, krome_tiny)
  NO2_DUSTI = max(NO2_DUSTI * factor, krome_tiny)
  HNO_DUSTI = max(HNO_DUSTI * factor, krome_tiny)
  O2H_DUSTI = max(O2H_DUSTI * factor, krome_tiny)
  H2CN_DUSTI = max(H2CN_DUSTI * factor, krome_tiny)
  MG_DUSTI = max(MG_DUSTI * factor, krome_tiny)
  HNC_DUSTI = max(HNC_DUSTI * factor, krome_tiny)
  E_DUSTI = max(E_DUSTI * factor, krome_tiny)
  SIO_DUSTI = max(SIO_DUSTI * factor, krome_tiny)
  HCOII = max(HCOII * factor, krome_tiny)
  HII = max(HII * factor, krome_tiny)
  HOCII = max(HOCII * factor, krome_tiny)
  CII = max(CII * factor, krome_tiny)
  CH2II = max(CH2II * factor, krome_tiny)
  CHII = max(CHII * factor, krome_tiny)
  H2COII = max(H2COII * factor, krome_tiny)
  MGII = max(MGII * factor, krome_tiny)
  NH3II = max(NH3II * factor, krome_tiny)
  NOII = max(NOII * factor, krome_tiny)
  SIII = max(SIII * factor, krome_tiny)
  SIC2II = max(SIC2II * factor, krome_tiny)
  SIC3II = max(SIC3II * factor, krome_tiny)
  SICII = max(SICII * factor, krome_tiny)
  SIH2II = max(SIH2II * factor, krome_tiny)
  SIH3II = max(SIH3II * factor, krome_tiny)
  CNII = max(CNII * factor, krome_tiny)
  COII = max(COII * factor, krome_tiny)
  N2II = max(N2II * factor, krome_tiny)
  O2II = max(O2II * factor, krome_tiny)
  H2OII = max(H2OII * factor, krome_tiny)
  NH2II = max(NH2II * factor, krome_tiny)
  OII = max(OII * factor, krome_tiny)
  OHII = max(OHII * factor, krome_tiny)
  CH3II = max(CH3II * factor, krome_tiny)
  CH4II = max(CH4II * factor, krome_tiny)
  NII = max(NII * factor, krome_tiny)
  HCNII = max(HCNII * factor, krome_tiny)
  NHII = max(NHII * factor, krome_tiny)
  SIH4II = max(SIH4II * factor, krome_tiny)
  SIHII = max(SIHII * factor, krome_tiny)
  SIOII = max(SIOII * factor, krome_tiny)
  H2II = max(H2II * factor, krome_tiny)
  HeII = max(HeII * factor, krome_tiny)
  HNOII = max(HNOII * factor, krome_tiny)
  H2NOII = max(H2NOII * factor, krome_tiny)
  H3II = max(H3II * factor, krome_tiny)
  H3COII = max(H3COII * factor, krome_tiny)
  H3OII = max(H3OII * factor, krome_tiny)
  HCNHII = max(HCNHII * factor, krome_tiny)
  HCO2II = max(HCO2II * factor, krome_tiny)
  HeHII = max(HeHII * factor, krome_tiny)
  N2HII = max(N2HII * factor, krome_tiny)
  O2HII = max(O2HII * factor, krome_tiny)
  SIH5II = max(SIH5II * factor, krome_tiny)
  SIOHII = max(SIOHII * factor, krome_tiny)

        !convert to number densities
  De = De * dom
  CHI = CHI * dom * 0.0769230769231d0
  OI = OI * dom * 0.0625d0
  HNCI = HNCI * dom * 0.037037037037d0
  HCNI = HCNI * dom * 0.037037037037d0
  H2I = H2I * dom * 0.5d0
  CI = CI * dom * 0.0833333333333d0
  HI = HI * dom
  H2OI = H2OI * dom * 0.0555555555556d0
  OHI = OHI * dom * 0.0588235294118d0
  O2I = O2I * dom * 0.03125d0
  CH2I = CH2I * dom * 0.0714285714286d0
  H2COI = H2COI * dom * 0.0333333333333d0
  HCOI = HCOI * dom * 0.0344827586207d0
  MGI = MGI * dom * 0.0416666666667d0
  NH3I = NH3I * dom * 0.0588235294118d0
  NOI = NOI * dom * 0.0333333333333d0
  SII = SII * dom * 0.0357142857143d0
  SIC2I = SIC2I * dom * 0.0192307692308d0
  SIC3I = SIC3I * dom * 0.015625d0
  SICI = SICI * dom * 0.025d0
  SIH2I = SIH2I * dom * 0.0333333333333d0
  SIH3I = SIH3I * dom * 0.0322580645161d0
  CNI = CNI * dom * 0.0384615384615d0
  COI = COI * dom * 0.0357142857143d0
  N2I = N2I * dom * 0.0357142857143d0
  NH2I = NH2I * dom * 0.0625d0
  CH3I = CH3I * dom * 0.0666666666667d0
  CH4I = CH4I * dom * 0.0625d0
  NI = NI * dom * 0.0714285714286d0
  NHI = NHI * dom * 0.0666666666667d0
  SIH4I = SIH4I * dom * 0.03125d0
  SIHI = SIHI * dom * 0.0344827586207d0
  SIOI = SIOI * dom * 0.0227272727273d0
  HeI = HeI * dom * 0.25d0
  HNOI = HNOI * dom * 0.0322580645161d0
  CH3OHI = CH3OHI * dom * 0.03125d0
  CO2I = CO2I * dom * 0.0227272727273d0
  H2CNI = H2CNI * dom * 0.0357142857143d0
  H2SIOI = H2SIOI * dom * 0.0217391304348d0
  HNCOI = HNCOI * dom * 0.0232558139535d0
  NO2I = NO2I * dom * 0.0217391304348d0
  O2HI = O2HI * dom * 0.030303030303d0
  OCNI = OCNI * dom * 0.0238095238095d0
  CH3OH_DUSTI = CH3OH_DUSTI * dom * 0.03125d0
  HNCO_DUSTI = HNCO_DUSTI * dom * 0.0232558139535d0
  H2CO_DUSTI = H2CO_DUSTI * dom * 0.0333333333333d0
  SIH4_DUSTI = SIH4_DUSTI * dom * 0.03125d0
  H2SIO_DUSTI = H2SIO_DUSTI * dom * 0.0217391304348d0
  SIC_DUSTI = SIC_DUSTI * dom * 0.025d0
  SIC2_DUSTI = SIC2_DUSTI * dom * 0.0192307692308d0
  SIC3_DUSTI = SIC3_DUSTI * dom * 0.015625d0
  CH4_DUSTI = CH4_DUSTI * dom * 0.0625d0
  CO_DUSTI = CO_DUSTI * dom * 0.0357142857143d0
  H2O_DUSTI = H2O_DUSTI * dom * 0.0555555555556d0
  NO_DUSTI = NO_DUSTI * dom * 0.0333333333333d0
  CO2_DUSTI = CO2_DUSTI * dom * 0.0227272727273d0
  N2_DUSTI = N2_DUSTI * dom * 0.0357142857143d0
  HCN_DUSTI = HCN_DUSTI * dom * 0.037037037037d0
  NH3_DUSTI = NH3_DUSTI * dom * 0.0588235294118d0
  O2_DUSTI = O2_DUSTI * dom * 0.03125d0
  NO2_DUSTI = NO2_DUSTI * dom * 0.0217391304348d0
  HNO_DUSTI = HNO_DUSTI * dom * 0.0322580645161d0
  O2H_DUSTI = O2H_DUSTI * dom * 0.030303030303d0
  H2CN_DUSTI = H2CN_DUSTI * dom * 0.0357142857143d0
  MG_DUSTI = MG_DUSTI * dom * 0.0416666666667d0
  HNC_DUSTI = HNC_DUSTI * dom * 0.037037037037d0
  E_DUSTI = E_DUSTI * dom
  SIO_DUSTI = SIO_DUSTI * dom * 0.0227272727273d0
  HCOII = HCOII * dom * 0.0344827586207d0
  HII = HII * dom
  HOCII = HOCII * dom * 0.0344827586207d0
  CII = CII * dom * 0.0833333333333d0
  CH2II = CH2II * dom * 0.0714285714286d0
  CHII = CHII * dom * 0.0769230769231d0
  H2COII = H2COII * dom * 0.0333333333333d0
  MGII = MGII * dom * 0.0416666666667d0
  NH3II = NH3II * dom * 0.0588235294118d0
  NOII = NOII * dom * 0.0333333333333d0
  SIII = SIII * dom * 0.0357142857143d0
  SIC2II = SIC2II * dom * 0.0192307692308d0
  SIC3II = SIC3II * dom * 0.015625d0
  SICII = SICII * dom * 0.025d0
  SIH2II = SIH2II * dom * 0.0333333333333d0
  SIH3II = SIH3II * dom * 0.0322580645161d0
  CNII = CNII * dom * 0.0384615384615d0
  COII = COII * dom * 0.0357142857143d0
  N2II = N2II * dom * 0.0357142857143d0
  O2II = O2II * dom * 0.03125d0
  H2OII = H2OII * dom * 0.0555555555556d0
  NH2II = NH2II * dom * 0.0625d0
  OII = OII * dom * 0.0625d0
  OHII = OHII * dom * 0.0588235294118d0
  CH3II = CH3II * dom * 0.0666666666667d0
  CH4II = CH4II * dom * 0.0625d0
  NII = NII * dom * 0.0714285714286d0
  HCNII = HCNII * dom * 0.037037037037d0
  NHII = NHII * dom * 0.0666666666667d0
  SIH4II = SIH4II * dom * 0.03125d0
  SIHII = SIHII * dom * 0.0344827586207d0
  SIOII = SIOII * dom * 0.0227272727273d0
  H2II = H2II * dom * 0.5d0
  HeII = HeII * dom * 0.25d0
  HNOII = HNOII * dom * 0.0322580645161d0
  H2NOII = H2NOII * dom * 0.03125d0
  H3II = H3II * dom * 0.333333333333d0
  H3COII = H3COII * dom * 0.0322580645161d0
  H3OII = H3OII * dom * 0.0526315789474d0
  HCNHII = HCNHII * dom * 0.0357142857143d0
  HCO2II = HCO2II * dom * 0.0222222222222d0
  HeHII = HeHII * dom * 0.2d0
  N2HII = N2HII * dom * 0.0344827586207d0
  O2HII = O2HII * dom * 0.030303030303d0
  SIH5II = SIH5II * dom * 0.030303030303d0
  SIOHII = SIOHII * dom * 0.0222222222222d0

  ! uflux = 0.0
  ! uflux(1:in-1, :, :) = max(d(1:in-1, :, :) * u(1:in-1, :, :) - d(2:in, :, :) * u(2:in, :, :), 0.0)
  ! uflux(2:in, :, :) -= min(d(2:in, :, :) * u(2:in, :, :) - d(1:in-1, :, :) * u(1:in-1, :, :), 0.0)
  
  ! vflux = 0.0
  ! vflux(1:in-1, :, :) = max(d(1:in-1, :, :) * v(1:in-1, :, :) - d(2:in, :, :) * v(2:in, :, :), 0.0)
  ! vflux(2:in, :, :) -= min(d(2:in, :, :) * v(2:in, :, :) - d(1:in-1, :, :) * v(1:in-1, :, :), 0.0)

  ! wflux = 0.0
  ! wflux(1:in-1, :, :) = max(d(1:in-1, :, :) * w(1:in-1, :, :) - d(2:in, :, :) * w(2:in, :, :), 0.0)
  ! wflux(2:in, :, :) -= min(d(2:in, :, :) * w(2:in, :, :) - d(1:in-1, :, :) * w(1:in-1, :, :), 0.0)

  ! totalinflux = uflux + vflux + wflux

  !loop over zones
  do k = ks+1, ke+1
    do j = js+1, je+1
      do i = is+1, ie+1

        !rhogas = #KROME_sum to be removed

        !convert to number densities
        krome_x(krome_idx_E) = De(i,j,k) 
        krome_x(krome_idx_CH) = CHI(i,j,k) 
        krome_x(krome_idx_O) = OI(i,j,k) 
        krome_x(krome_idx_HNC) = HNCI(i,j,k) 
        krome_x(krome_idx_HCN) = HCNI(i,j,k) 
        krome_x(krome_idx_H2) = H2I(i,j,k) 
        krome_x(krome_idx_C) = CI(i,j,k) 
        krome_x(krome_idx_H) = HI(i,j,k) 
        krome_x(krome_idx_H2O) = H2OI(i,j,k) 
        krome_x(krome_idx_OH) = OHI(i,j,k) 
        krome_x(krome_idx_O2) = O2I(i,j,k) 
        krome_x(krome_idx_CH2) = CH2I(i,j,k) 
        krome_x(krome_idx_H2CO) = H2COI(i,j,k) 
        krome_x(krome_idx_HCO) = HCOI(i,j,k) 
        krome_x(krome_idx_MG) = MGI(i,j,k) 
        krome_x(krome_idx_NH3) = NH3I(i,j,k) 
        krome_x(krome_idx_NO) = NOI(i,j,k) 
        krome_x(krome_idx_SI) = SII(i,j,k) 
        krome_x(krome_idx_SIC2) = SIC2I(i,j,k) 
        krome_x(krome_idx_SIC3) = SIC3I(i,j,k) 
        krome_x(krome_idx_SIC) = SICI(i,j,k) 
        krome_x(krome_idx_SIH2) = SIH2I(i,j,k) 
        krome_x(krome_idx_SIH3) = SIH3I(i,j,k) 
        krome_x(krome_idx_CN) = CNI(i,j,k) 
        krome_x(krome_idx_CO) = COI(i,j,k) 
        krome_x(krome_idx_N2) = N2I(i,j,k) 
        krome_x(krome_idx_NH2) = NH2I(i,j,k) 
        krome_x(krome_idx_CH3) = CH3I(i,j,k) 
        krome_x(krome_idx_CH4) = CH4I(i,j,k) 
        krome_x(krome_idx_N) = NI(i,j,k) 
        krome_x(krome_idx_NH) = NHI(i,j,k) 
        krome_x(krome_idx_SIH4) = SIH4I(i,j,k) 
        krome_x(krome_idx_SIH) = SIHI(i,j,k) 
        krome_x(krome_idx_SIO) = SIOI(i,j,k) 
        krome_x(krome_idx_HE) = HeI(i,j,k) 
        krome_x(krome_idx_HNO) = HNOI(i,j,k) 
        krome_x(krome_idx_CH3OH) = CH3OHI(i,j,k) 
        krome_x(krome_idx_CO2) = CO2I(i,j,k) 
        krome_x(krome_idx_H2CN) = H2CNI(i,j,k) 
        krome_x(krome_idx_H2SIO) = H2SIOI(i,j,k) 
        krome_x(krome_idx_HNCO) = HNCOI(i,j,k) 
        krome_x(krome_idx_NO2) = NO2I(i,j,k) 
        krome_x(krome_idx_O2H) = O2HI(i,j,k) 
        krome_x(krome_idx_OCN) = OCNI(i,j,k) 
        krome_x(krome_idx_CH3OH_DUST) = CH3OH_DUSTI(i,j,k) 
        krome_x(krome_idx_HNCO_DUST) = HNCO_DUSTI(i,j,k) 
        krome_x(krome_idx_H2CO_DUST) = H2CO_DUSTI(i,j,k) 
        krome_x(krome_idx_SIH4_DUST) = SIH4_DUSTI(i,j,k) 
        krome_x(krome_idx_H2SIO_DUST) = H2SIO_DUSTI(i,j,k) 
        krome_x(krome_idx_SIC_DUST) = SIC_DUSTI(i,j,k) 
        krome_x(krome_idx_SIC2_DUST) = SIC2_DUSTI(i,j,k) 
        krome_x(krome_idx_SIC3_DUST) = SIC3_DUSTI(i,j,k) 
        krome_x(krome_idx_CH4_DUST) = CH4_DUSTI(i,j,k) 
        krome_x(krome_idx_CO_DUST) = CO_DUSTI(i,j,k) 
        krome_x(krome_idx_H2O_DUST) = H2O_DUSTI(i,j,k) 
        krome_x(krome_idx_NO_DUST) = NO_DUSTI(i,j,k) 
        krome_x(krome_idx_CO2_DUST) = CO2_DUSTI(i,j,k) 
        krome_x(krome_idx_N2_DUST) = N2_DUSTI(i,j,k) 
        krome_x(krome_idx_HCN_DUST) = HCN_DUSTI(i,j,k) 
        krome_x(krome_idx_NH3_DUST) = NH3_DUSTI(i,j,k) 
        krome_x(krome_idx_O2_DUST) = O2_DUSTI(i,j,k) 
        krome_x(krome_idx_NO2_DUST) = NO2_DUSTI(i,j,k) 
        krome_x(krome_idx_HNO_DUST) = HNO_DUSTI(i,j,k) 
        krome_x(krome_idx_O2H_DUST) = O2H_DUSTI(i,j,k) 
        krome_x(krome_idx_H2CN_DUST) = H2CN_DUSTI(i,j,k) 
        krome_x(krome_idx_MG_DUST) = MG_DUSTI(i,j,k) 
        krome_x(krome_idx_HNC_DUST) = HNC_DUSTI(i,j,k) 
        krome_x(krome_idx_E_DUST) = E_DUSTI(i,j,k) 
        krome_x(krome_idx_SIO_DUST) = SIO_DUSTI(i,j,k) 
        krome_x(krome_idx_HCOj) = HCOII(i,j,k) 
        krome_x(krome_idx_Hj) = HII(i,j,k) 
        krome_x(krome_idx_HOCj) = HOCII(i,j,k) 
        krome_x(krome_idx_Cj) = CII(i,j,k) 
        krome_x(krome_idx_CH2j) = CH2II(i,j,k) 
        krome_x(krome_idx_CHj) = CHII(i,j,k) 
        krome_x(krome_idx_H2COj) = H2COII(i,j,k) 
        krome_x(krome_idx_MGj) = MGII(i,j,k) 
        krome_x(krome_idx_NH3j) = NH3II(i,j,k) 
        krome_x(krome_idx_NOj) = NOII(i,j,k) 
        krome_x(krome_idx_SIj) = SIII(i,j,k) 
        krome_x(krome_idx_SIC2j) = SIC2II(i,j,k) 
        krome_x(krome_idx_SIC3j) = SIC3II(i,j,k) 
        krome_x(krome_idx_SICj) = SICII(i,j,k) 
        krome_x(krome_idx_SIH2j) = SIH2II(i,j,k) 
        krome_x(krome_idx_SIH3j) = SIH3II(i,j,k) 
        krome_x(krome_idx_CNj) = CNII(i,j,k) 
        krome_x(krome_idx_COj) = COII(i,j,k) 
        krome_x(krome_idx_N2j) = N2II(i,j,k) 
        krome_x(krome_idx_O2j) = O2II(i,j,k) 
        krome_x(krome_idx_H2Oj) = H2OII(i,j,k) 
        krome_x(krome_idx_NH2j) = NH2II(i,j,k) 
        krome_x(krome_idx_Oj) = OII(i,j,k) 
        krome_x(krome_idx_OHj) = OHII(i,j,k) 
        krome_x(krome_idx_CH3j) = CH3II(i,j,k) 
        krome_x(krome_idx_CH4j) = CH4II(i,j,k) 
        krome_x(krome_idx_Nj) = NII(i,j,k) 
        krome_x(krome_idx_HCNj) = HCNII(i,j,k) 
        krome_x(krome_idx_NHj) = NHII(i,j,k) 
        krome_x(krome_idx_SIH4j) = SIH4II(i,j,k) 
        krome_x(krome_idx_SIHj) = SIHII(i,j,k) 
        krome_x(krome_idx_SIOj) = SIOII(i,j,k) 
        krome_x(krome_idx_H2j) = H2II(i,j,k) 
        krome_x(krome_idx_HEj) = HeII(i,j,k) 
        krome_x(krome_idx_HNOj) = HNOII(i,j,k) 
        krome_x(krome_idx_H2NOj) = H2NOII(i,j,k) 
        krome_x(krome_idx_H3j) = H3II(i,j,k) 
        krome_x(krome_idx_H3COj) = H3COII(i,j,k) 
        krome_x(krome_idx_H3Oj) = H3OII(i,j,k) 
        krome_x(krome_idx_HCNHj) = HCNHII(i,j,k) 
        krome_x(krome_idx_HCO2j) = HCO2II(i,j,k) 
        krome_x(krome_idx_HEHj) = HeHII(i,j,k) 
        krome_x(krome_idx_N2Hj) = N2HII(i,j,k) 
        krome_x(krome_idx_O2Hj) = O2HII(i,j,k) 
        krome_x(krome_idx_SIH5j) = SIH5II(i,j,k) 
        krome_x(krome_idx_SIOHj) = SIOHII(i,j,k) 

        call evaluate_tgas(d(i,j,k), e(i,j,k), ge(i,j,k),&
            u(i,j,k), v(i,j,k), w(i,j,k),&
            krome_x(:),imethod,idual,idim,tgas,&
            utem)

        !store old tgas
        tgasold = tgas

        sputtering = 0.0
        ! >>>>>>>>>>>>>> calculate the sputtering rate <<<<<<<<<<<<<<<
        if (opt_spu > 0) then
          ul = min(u(i, j, k) - u(i-1, j, k), 0.d0)
          ur = max(u(i, j, k) - u(i+1, j, k), 0.d0)
          vl = min(v(i, j, k) - v(i, j-1, k), 0.d0)
          vr = max(v(i, j, k) - v(i, j+1, k), 0.d0)
          wl = min(w(i, j, k) - w(i, j, k-1), 0.d0)
          wr = max(w(i, j, k) - w(i, j, k+1), 0.d0)

          do projidx = 1, 6
            specidx = projectiles(projidx)

            if (projidx == 1) then
              numdens = H2I
            else if (projidx == 2) then
              numdens = HEI
            else if (projidx == 3) then
              numdens = CI
            else if (projidx == 4) then
              numdens = OI
            else if (projidx == 5) then
              numdens = SiI
            else if (projidx == 6) then
              numdens = COI
            end if

            ulyields = iceYield(specidx, ul)
            uryields = iceYield(specidx, ur)
            vlyields = iceYield(specidx, vl)
            vryields = iceYield(specidx, vr)
            wlyields = iceYield(specidx, wl)
            wryields = iceYield(specidx, wr)
            ulflux = min(numdens(i, j, k) * u(i, j, k) - numdens(i-1, j, k) * u(i-1, j, k), 0.d0)
            urflux = max(numdens(i, j, k) * u(i, j, k) - numdens(i+1, j, k) * u(i+1, j, k), 0.d0)
            vlflux = min(numdens(i, j, k) * v(i, j, k) - numdens(i, j-1, k) * v(i, j-1, k), 0.d0)
            vrflux = max(numdens(i, j, k) * v(i, j, k) - numdens(i, j+1, k) * v(i, j+1, k), 0.d0)
            wlflux = min(numdens(i, j, k) * w(i, j, k) - numdens(i, j, k-1) * w(i, j, k-1), 0.d0)
            wrflux = max(numdens(i, j, k) * w(i, j, k) - numdens(i, j, k+1) * w(i, j, k+1), 0.d0)

            sputtering = sputtering &
                    & + abs(ulflux) * ulyields + abs(urflux) * uryields &
                    & + abs(vlflux) * vlyields + abs(vrflux) * vryields &
                    & + abs(wlflux) * wlyields + abs(wrflux) * wryields 
          
          end do
        end if
        ! >>>>>>>>>>>>>> end of sputtering rate calculation <<<<<<<<<<<<<<<

        dt_hydro = utim*dt !dt*time_conversion

        !call KROME solver
        call krome(krome_x(:),tgas,dt_hydro,cellsize*uxyz,sputtering)

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
        SII(i,j,k) = krome_x(krome_idx_SI) * idom * 28d0
        SIC2I(i,j,k) = krome_x(krome_idx_SIC2) * idom * 52d0
        SIC3I(i,j,k) = krome_x(krome_idx_SIC3) * idom * 64d0
        SICI(i,j,k) = krome_x(krome_idx_SIC) * idom * 40d0
        SIH2I(i,j,k) = krome_x(krome_idx_SIH2) * idom * 30d0
        SIH3I(i,j,k) = krome_x(krome_idx_SIH3) * idom * 31d0
        CNI(i,j,k) = krome_x(krome_idx_CN) * idom * 26d0
        COI(i,j,k) = krome_x(krome_idx_CO) * idom * 28d0
        N2I(i,j,k) = krome_x(krome_idx_N2) * idom * 28d0
        NH2I(i,j,k) = krome_x(krome_idx_NH2) * idom * 16d0
        CH3I(i,j,k) = krome_x(krome_idx_CH3) * idom * 15d0
        CH4I(i,j,k) = krome_x(krome_idx_CH4) * idom * 16d0
        NI(i,j,k) = krome_x(krome_idx_N) * idom * 14d0
        NHI(i,j,k) = krome_x(krome_idx_NH) * idom * 15d0
        SIH4I(i,j,k) = krome_x(krome_idx_SIH4) * idom * 32d0
        SIHI(i,j,k) = krome_x(krome_idx_SIH) * idom * 29d0
        SIOI(i,j,k) = krome_x(krome_idx_SIO) * idom * 44d0
        HeI(i,j,k) = krome_x(krome_idx_HE) * idom * 4d0
        HNOI(i,j,k) = krome_x(krome_idx_HNO) * idom * 31d0
        CH3OHI(i,j,k) = krome_x(krome_idx_CH3OH) * idom * 32d0
        CO2I(i,j,k) = krome_x(krome_idx_CO2) * idom * 44d0
        H2CNI(i,j,k) = krome_x(krome_idx_H2CN) * idom * 28d0
        H2SIOI(i,j,k) = krome_x(krome_idx_H2SIO) * idom * 46d0
        HNCOI(i,j,k) = krome_x(krome_idx_HNCO) * idom * 43d0
        NO2I(i,j,k) = krome_x(krome_idx_NO2) * idom * 46d0
        O2HI(i,j,k) = krome_x(krome_idx_O2H) * idom * 33d0
        OCNI(i,j,k) = krome_x(krome_idx_OCN) * idom * 42d0
        CH3OH_DUSTI(i,j,k) = krome_x(krome_idx_CH3OH_DUST) * idom * 32d0
        HNCO_DUSTI(i,j,k) = krome_x(krome_idx_HNCO_DUST) * idom * 43d0
        H2CO_DUSTI(i,j,k) = krome_x(krome_idx_H2CO_DUST) * idom * 30d0
        SIH4_DUSTI(i,j,k) = krome_x(krome_idx_SIH4_DUST) * idom * 32d0
        H2SIO_DUSTI(i,j,k) = krome_x(krome_idx_H2SIO_DUST) * idom * 46d0
        SIC_DUSTI(i,j,k) = krome_x(krome_idx_SIC_DUST) * idom * 40d0
        SIC2_DUSTI(i,j,k) = krome_x(krome_idx_SIC2_DUST) * idom * 52d0
        SIC3_DUSTI(i,j,k) = krome_x(krome_idx_SIC3_DUST) * idom * 64d0
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
        SIO_DUSTI(i,j,k) = krome_x(krome_idx_SIO_DUST) * idom * 44d0
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
        SIII(i,j,k) = krome_x(krome_idx_SIj) * idom * 28d0
        SIC2II(i,j,k) = krome_x(krome_idx_SIC2j) * idom * 52d0
        SIC3II(i,j,k) = krome_x(krome_idx_SIC3j) * idom * 64d0
        SICII(i,j,k) = krome_x(krome_idx_SICj) * idom * 40d0
        SIH2II(i,j,k) = krome_x(krome_idx_SIH2j) * idom * 30d0
        SIH3II(i,j,k) = krome_x(krome_idx_SIH3j) * idom * 31d0
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
        SIH4II(i,j,k) = krome_x(krome_idx_SIH4j) * idom * 32d0
        SIHII(i,j,k) = krome_x(krome_idx_SIHj) * idom * 29d0
        SIOII(i,j,k) = krome_x(krome_idx_SIOj) * idom * 44d0
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
        SIH5II(i,j,k) = krome_x(krome_idx_SIH5j) * idom * 33d0
        SIOHII(i,j,k) = krome_x(krome_idx_SIOHj) * idom * 45d0

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
        SII(i,j,k) = SII(i,j,k) * factor
        SIC2I(i,j,k) = SIC2I(i,j,k) * factor
        SIC3I(i,j,k) = SIC3I(i,j,k) * factor
        SICI(i,j,k) = SICI(i,j,k) * factor
        SIH2I(i,j,k) = SIH2I(i,j,k) * factor
        SIH3I(i,j,k) = SIH3I(i,j,k) * factor
        CNI(i,j,k) = CNI(i,j,k) * factor
        COI(i,j,k) = COI(i,j,k) * factor
        N2I(i,j,k) = N2I(i,j,k) * factor
        NH2I(i,j,k) = NH2I(i,j,k) * factor
        CH3I(i,j,k) = CH3I(i,j,k) * factor
        CH4I(i,j,k) = CH4I(i,j,k) * factor
        NI(i,j,k) = NI(i,j,k) * factor
        NHI(i,j,k) = NHI(i,j,k) * factor
        SIH4I(i,j,k) = SIH4I(i,j,k) * factor
        SIHI(i,j,k) = SIHI(i,j,k) * factor
        SIOI(i,j,k) = SIOI(i,j,k) * factor
        HeI(i,j,k) = HeI(i,j,k) * factor
        HNOI(i,j,k) = HNOI(i,j,k) * factor
        CH3OHI(i,j,k) = CH3OHI(i,j,k) * factor
        CO2I(i,j,k) = CO2I(i,j,k) * factor
        H2CNI(i,j,k) = H2CNI(i,j,k) * factor
        H2SIOI(i,j,k) = H2SIOI(i,j,k) * factor
        HNCOI(i,j,k) = HNCOI(i,j,k) * factor
        NO2I(i,j,k) = NO2I(i,j,k) * factor
        O2HI(i,j,k) = O2HI(i,j,k) * factor
        OCNI(i,j,k) = OCNI(i,j,k) * factor
        CH3OH_DUSTI(i,j,k) = CH3OH_DUSTI(i,j,k) * factor
        HNCO_DUSTI(i,j,k) = HNCO_DUSTI(i,j,k) * factor
        H2CO_DUSTI(i,j,k) = H2CO_DUSTI(i,j,k) * factor
        SIH4_DUSTI(i,j,k) = SIH4_DUSTI(i,j,k) * factor
        H2SIO_DUSTI(i,j,k) = H2SIO_DUSTI(i,j,k) * factor
        SIC_DUSTI(i,j,k) = SIC_DUSTI(i,j,k) * factor
        SIC2_DUSTI(i,j,k) = SIC2_DUSTI(i,j,k) * factor
        SIC3_DUSTI(i,j,k) = SIC3_DUSTI(i,j,k) * factor
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
        SIO_DUSTI(i,j,k) = SIO_DUSTI(i,j,k) * factor
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
        SIII(i,j,k) = SIII(i,j,k) * factor
        SIC2II(i,j,k) = SIC2II(i,j,k) * factor
        SIC3II(i,j,k) = SIC3II(i,j,k) * factor
        SICII(i,j,k) = SICII(i,j,k) * factor
        SIH2II(i,j,k) = SIH2II(i,j,k) * factor
        SIH3II(i,j,k) = SIH3II(i,j,k) * factor
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
        SIH4II(i,j,k) = SIH4II(i,j,k) * factor
        SIHII(i,j,k) = SIHII(i,j,k) * factor
        SIOII(i,j,k) = SIOII(i,j,k) * factor
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
        SIH5II(i,j,k) = SIH5II(i,j,k) * factor
        SIOHII(i,j,k) = SIOHII(i,j,k) * factor

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
        SII(i,j,k) = max(SII(i,j,k), krome_tiny)
        SIC2I(i,j,k) = max(SIC2I(i,j,k), krome_tiny)
        SIC3I(i,j,k) = max(SIC3I(i,j,k), krome_tiny)
        SICI(i,j,k) = max(SICI(i,j,k), krome_tiny)
        SIH2I(i,j,k) = max(SIH2I(i,j,k), krome_tiny)
        SIH3I(i,j,k) = max(SIH3I(i,j,k), krome_tiny)
        CNI(i,j,k) = max(CNI(i,j,k), krome_tiny)
        COI(i,j,k) = max(COI(i,j,k), krome_tiny)
        N2I(i,j,k) = max(N2I(i,j,k), krome_tiny)
        NH2I(i,j,k) = max(NH2I(i,j,k), krome_tiny)
        CH3I(i,j,k) = max(CH3I(i,j,k), krome_tiny)
        CH4I(i,j,k) = max(CH4I(i,j,k), krome_tiny)
        NI(i,j,k) = max(NI(i,j,k), krome_tiny)
        NHI(i,j,k) = max(NHI(i,j,k), krome_tiny)
        SIH4I(i,j,k) = max(SIH4I(i,j,k), krome_tiny)
        SIHI(i,j,k) = max(SIHI(i,j,k), krome_tiny)
        SIOI(i,j,k) = max(SIOI(i,j,k), krome_tiny)
        HeI(i,j,k) = max(HeI(i,j,k), krome_tiny)
        HNOI(i,j,k) = max(HNOI(i,j,k), krome_tiny)
        CH3OHI(i,j,k) = max(CH3OHI(i,j,k), krome_tiny)
        CO2I(i,j,k) = max(CO2I(i,j,k), krome_tiny)
        H2CNI(i,j,k) = max(H2CNI(i,j,k), krome_tiny)
        H2SIOI(i,j,k) = max(H2SIOI(i,j,k), krome_tiny)
        HNCOI(i,j,k) = max(HNCOI(i,j,k), krome_tiny)
        NO2I(i,j,k) = max(NO2I(i,j,k), krome_tiny)
        O2HI(i,j,k) = max(O2HI(i,j,k), krome_tiny)
        OCNI(i,j,k) = max(OCNI(i,j,k), krome_tiny)
        CH3OH_DUSTI(i,j,k) = max(CH3OH_DUSTI(i,j,k), krome_tiny)
        HNCO_DUSTI(i,j,k) = max(HNCO_DUSTI(i,j,k), krome_tiny)
        H2CO_DUSTI(i,j,k) = max(H2CO_DUSTI(i,j,k), krome_tiny)
        SIH4_DUSTI(i,j,k) = max(SIH4_DUSTI(i,j,k), krome_tiny)
        H2SIO_DUSTI(i,j,k) = max(H2SIO_DUSTI(i,j,k), krome_tiny)
        SIC_DUSTI(i,j,k) = max(SIC_DUSTI(i,j,k), krome_tiny)
        SIC2_DUSTI(i,j,k) = max(SIC2_DUSTI(i,j,k), krome_tiny)
        SIC3_DUSTI(i,j,k) = max(SIC3_DUSTI(i,j,k), krome_tiny)
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
        SIO_DUSTI(i,j,k) = max(SIO_DUSTI(i,j,k), krome_tiny)
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
        SIII(i,j,k) = max(SIII(i,j,k), krome_tiny)
        SIC2II(i,j,k) = max(SIC2II(i,j,k), krome_tiny)
        SIC3II(i,j,k) = max(SIC3II(i,j,k), krome_tiny)
        SICII(i,j,k) = max(SICII(i,j,k), krome_tiny)
        SIH2II(i,j,k) = max(SIH2II(i,j,k), krome_tiny)
        SIH3II(i,j,k) = max(SIH3II(i,j,k), krome_tiny)
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
        SIH4II(i,j,k) = max(SIH4II(i,j,k), krome_tiny)
        SIHII(i,j,k) = max(SIHII(i,j,k), krome_tiny)
        SIOII(i,j,k) = max(SIOII(i,j,k), krome_tiny)
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
        SIH5II(i,j,k) = max(SIH5II(i,j,k), krome_tiny)
        SIOHII(i,j,k) = max(SIOHII(i,j,k), krome_tiny)

      end do
    end do
  end do

end subroutine krome_driver
