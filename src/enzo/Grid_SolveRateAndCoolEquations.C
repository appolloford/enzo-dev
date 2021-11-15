/***********************************************************************
/
/  GRID CLASS (SOLVE THE COOLING/HEATING AND RATE EQUATIONS)
/
/  written by: Greg Bryan
/  date:       October, 1996
/  modified1:  July, 2005 to solve cool and rate equations simultaneously
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"
#include "Gadget.h"

/* This parameter controls whether the cooling function recomputes
   the metal cooling rates.  It is reset by RadiationFieldUpdate. */

extern int RadiationFieldRecomputeMetalRates;

/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);
int RadiationFieldCalculateRates(FLOAT Time);
int FindField(int field, int farray[], int numfields);
double ReturnWallTime();
extern "C" void FORTRAN_NAME(solve_rate_cool)(
    float *d, float *e, float *ge, float *u, float *v, float *w, float *de,
    float *HI, float *HII, float *HeI, float *HeII, float *HeIII,
    int *in, int *jn, int *kn, int *nratec, int *iexpand,
    hydro_method *imethod,
    int *idual, int *ispecies, int *imetal, int *imcool, int *idust, int *idim,
    int *is, int *js, int *ks, int *ie, int *je, int *ke, int *ih2co,
    int *ipiht, int *igammah,
    FLOAT *dx, float *dt, float *aye, float *redshift, float *temstart, float *temend,
    float *utem, float *uxyz, float *uaye, float *urho, float *utim,
    float *eta1, float *eta2, float *gamma, float *fh, float *dtoh,
    float *z_solar,
    float *k1a, float *k2a, float *k3a, float *k4a, float *k5a,
    float *k6a, float *k7a, float *k8a, float *k9a, float *k10a,
    float *k11a, float *k12a, float *k13a, float *k13dda, float *k14a,
    float *k15a,
    float *k16a, float *k17a, float *k18a, float *k19a, float *k22a,
    float *k24, float *k25, float *k26, float *k27, float *k28, float *k29,
    float *k30, float *k31,
    float *k50a, float *k51a, float *k52a, float *k53a, float *k54a,
    float *k55a, float *k56a,
    int *ndratec, float *dtemstart, float *dtemend, float *h2dusta,
    float *ncrna, float *ncrd1a, float *ncrd2a,
    float *ceHIa, float *ceHeIa, float *ceHeIIa, float *ciHIa,
    float *ciHeIa,
    float *ciHeISa, float *ciHeIIa, float *reHIIa, float *reHeII1a,
    float *reHeII2a, float *reHeIIIa, float *brema, float *compa, float *gammaha,
    float *comp_xraya, float *comp_temp,
    float *piHI, float *piHeI, float *piHeII,
    float *HM, float *H2I, float *H2II, float *DI, float *DII, float *HDI,
    float *metal,
    float *hyd01ka, float *h2k01a, float *vibha, float *rotha, float *rotla,
    float *gpldl, float *gphdl, float *HDltea, float *HDlowa,
    float *gaHIa, float *gaH2a, float *gaHea, float *gaHpa, float *gaela,
    float *gasgra, float *metala, int *n_xe, float *xe_start, float *xe_end,
    float *inutot, int *iradtype, int *nfreq, int *imetalregen,
    int *iradshield, float *avgsighp, float *avgsighep, float *avgsighe2p,
    int *iradtrans, int *iradcoupled, int *iradstep, int *ierr,
    int *irt_honly,
    float *kphHI, float *kphHeI, float *kphHeII,
    float *kdissH2I, float *photogamma,
    int *ih2optical, int *iciecool, int *ithreebody, float *ciecoa,
    int *icmbTfloor, int *iClHeat,
    float *clEleFra, int *clGridRank, int *clGridDim,
    float *clPar1, float *clPar2, float *clPar3, float *clPar4, float *clPar5,
    int *clDataSize, float *clCooling, float *clHeating);

extern "C" void FORTRAN_NAME(krome_driver)(
    float *d, float *e, float *ge, float *u, float *v, float *w,
    float *De, float *CHI, float *OI, float *HNCI,
    float *HCNI, float *H2I, float *CI, float *HI,
    float *H2OI, float *OHI, float *O2I, float *CH2I,
    float *H2COI, float *HCOI, float *MGI, float *NH3I,
    float *NOI, float *SII, float *SIC2I, float *SIC3I,
    float *SICI, float *SIH2I, float *SIH3I,
    float *CNI, float *COI, float *N2I, float *NH2I,
    float *CH3I, float *CH4I, float *NI, float *NHI,
    float *SIH4I, float *SIHI, float *SIOI, float *HeI,
    float *HNOI, float *CH3OHI, float *CO2I,
    float *H2CNI, float *H2SIOI, float *HNCOI,
    float *NO2I, float *O2HI, float *OCNI, float *CH3OH_DUSTI,
    float *HNCO_DUSTI, float *H2CO_DUSTI, float *SIH4_DUSTI,
    float *H2SIO_DUSTI, float *SIC_DUSTI, float *SIC2_DUSTI,
    float *SIC3_DUSTI, float *CH4_DUSTI, float *CO_DUSTI,
    float *H2O_DUSTI, float *NO_DUSTI, float *CO2_DUSTI,
    float *N2_DUSTI, float *HCN_DUSTI, float *NH3_DUSTI,
    float *O2_DUSTI, float *NO2_DUSTI, float *HNO_DUSTI,
    float *O2H_DUSTI, float *H2CN_DUSTI, float *MG_DUSTI,
    float *HNC_DUSTI, float *E_DUSTI, float *SIO_DUSTI,
    float *HCOII, float *HII, float *HOCII, float *CII,
    float *CH2II, float *CHII, float *H2COII,
    float *MGII, float *NH3II, float *NOII, float *SIII,
    float *SIC2II, float *SIC3II, float *SICII,
    float *SIH2II, float *SIH3II, float *CNII,
    float *COII, float *N2II, float *O2II, float *H2OII,
    float *NH2II, float *OII, float *OHII, float *CH3II,
    float *CH4II, float *NII, float *HCNII, float *NHII,
    float *SIH4II, float *SIHII, float *SIOII,
    float *H2II, float *HeII, float *HNOII, float *H2NOII,
    float *H3II, float *H3COII, float *H3OII,
    float *HCNHII, float *HCO2II, float *HeHII,
    float *N2HII, float *O2HII, float *SIH5II,
    float *SIOHII, int *in, int *jn, int *kn,
    hydro_method *imethod,
    int *idual, int *idim,
    int *is, int *js, int *ks, int *ie, int *je, int *ke,
    float *dt, float *aye,
    float *utem, float *uxyz, float *uaye, float *urho, float *utim,
    float *gamma, float *fh, float *dtoh, float *gridsize, int *opt_spu);

int grid::SolveRateAndCoolEquations(int RTCoupledSolverIntermediateStep)
{
  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

    /* Return if this doesn't concern us. */
#ifdef USE_KROME
  if (!use_krome)
    return SUCCESS;

  // Nothing to do in the first step if use_kromestep == 3
  if (use_krome && use_kromestep == 3 && TopGridCycle == 0)
    return SUCCESS;

  if (MultiSpecies != KROMESPECIES)
  {
    DebugCheck("Warning: species don't match with krome patches. Skip solving chemistry.\n");
    return SUCCESS;
  }
  if (debug && MyProcessorNumber == ROOT_PROCESSOR)
  {
    if (use_kromestep == 1 && Time > KromeTime + KromeDt)
      printf("SolveRateAndCool, KromeTime:%lf, KromeDt: %lf\n", KromeTime, KromeDt);
    else if (use_kromestep == 2)
      printf("SolveRateAndCool: use kromecycle\n");
    // else if (use_kromestep == 3 && TopGridCycle == KromeCycle)
    //   printf("SolveRateAndCool: use kromestep == 3, TopGridCycle:%d, KromeCycle:%d, KromeCycleSkip:%d\n",
    //           TopGridCycle, KromeCycle, KromeCycleSkip);
  }
#else
  if (!(MultiSpecies && RadiativeCooling))
    return SUCCESS;
#endif

  this->DebugCheck("SolveRadiativeCooling");

  /* Declarations */

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  int CHINum, OINum, HNCINum, HCNINum, CINum, H2OINum,
      OHINum, O2INum, CH2INum, H2COINum, HCOINum,
      MGINum, NH3INum, NOINum, SIINum, SIC2INum,
      SIC3INum, SICINum, SIH2INum, SIH3INum, CNINum,
      COINum, N2INum, NH2INum, CH3INum, CH4INum,
      NINum, NHINum, SIH4INum, SIHINum, SIOINum,
      HNOINum, CH3OHINum, CO2INum, H2CNINum, H2SIOINum,
      HNCOINum, NO2INum, O2HINum, OCNINum, CH3OH_DUSTINum,
      HNCO_DUSTINum, H2CO_DUSTINum, SIH4_DUSTINum,
      H2SIO_DUSTINum, SIC_DUSTINum, SIC2_DUSTINum,
      SIC3_DUSTINum, CH4_DUSTINum, CO_DUSTINum,
      H2O_DUSTINum, NO_DUSTINum, CO2_DUSTINum,
      N2_DUSTINum, HCN_DUSTINum, NH3_DUSTINum,
      O2_DUSTINum, NO2_DUSTINum, HNO_DUSTINum,
      O2H_DUSTINum, H2CN_DUSTINum, MG_DUSTINum,
      HNC_DUSTINum, E_DUSTINum, SIO_DUSTINum, HCOIINum,
      HOCIINum, CIINum, CH2IINum, CHIINum, H2COIINum,
      MGIINum, NH3IINum, NOIINum, SIIINum, SIC2IINum,
      SIC3IINum, SICIINum, SIH2IINum, SIH3IINum,
      CNIINum, COIINum, N2IINum, O2IINum, H2OIINum,
      NH2IINum, OIINum, OHIINum, CH3IINum, CH4IINum,
      NIINum, HCNIINum, NHIINum, SIH4IINum, SIHIINum,
      SIOIINum, HNOIINum, H2NOIINum, H3IINum, H3COIINum,
      H3OIINum, HCNHIINum, HCO2IINum, HeHIINum,
      N2HIINum, O2HIINum, SIH5IINum, SIOHIINum;

  FLOAT a = 1.0, dadt;

  /* Find fields: density, total energy, velocity1-3. */

  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL)
  {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  /* Find Multi-species fields. New routine from KROME */

  if (MultiSpecies)
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                              HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL)
    {
      ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
    }

  if (MultiSpecies)
    if (IdentifySpeciesFieldsKrome(
            DeNum, CHINum, OINum, HNCINum, HCNINum, H2INum,
            CINum, HINum, H2OINum, OHINum, O2INum, CH2INum,
            H2COINum, HCOINum, MGINum, NH3INum, NOINum,
            SIINum, SIC2INum, SIC3INum, SICINum, SIH2INum,
            SIH3INum, CNINum, COINum, N2INum, NH2INum,
            CH3INum, CH4INum, NINum, NHINum, SIH4INum,
            SIHINum, SIOINum, HeINum, HNOINum, CH3OHINum,
            CO2INum, H2CNINum, H2SIOINum, HNCOINum, NO2INum,
            O2HINum, OCNINum, CH3OH_DUSTINum, HNCO_DUSTINum,
            H2CO_DUSTINum, SIH4_DUSTINum, H2SIO_DUSTINum,
            SIC_DUSTINum, SIC2_DUSTINum, SIC3_DUSTINum,
            CH4_DUSTINum, CO_DUSTINum, H2O_DUSTINum,
            NO_DUSTINum, CO2_DUSTINum, N2_DUSTINum, HCN_DUSTINum,
            NH3_DUSTINum, O2_DUSTINum, NO2_DUSTINum,
            HNO_DUSTINum, O2H_DUSTINum, H2CN_DUSTINum,
            MG_DUSTINum, HNC_DUSTINum, E_DUSTINum, SIO_DUSTINum,
            HCOIINum, HIINum, HOCIINum, CIINum, CH2IINum,
            CHIINum, H2COIINum, MGIINum, NH3IINum, NOIINum,
            SIIINum, SIC2IINum, SIC3IINum, SICIINum,
            SIH2IINum, SIH3IINum, CNIINum, COIINum, N2IINum,
            O2IINum, H2OIINum, NH2IINum, OIINum, OHIINum,
            CH3IINum, CH4IINum, NIINum, HCNIINum, NHIINum,
            SIH4IINum, SIHIINum, SIOIINum, H2IINum, HeIINum,
            HNOIINum, H2NOIINum, H3IINum, H3COIINum,
            H3OIINum, HCNHIINum, HCO2IINum, HeHIINum,
            N2HIINum, O2HIINum, SIH5IINum, SIOHIINum) == FAIL)
    {
      ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
    }

  /* Find photo-ionization fields */

  int kphHINum, kphHeINum, kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum;
  int gammaNum;
  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum,
                                  kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum);

  /* Compute size of the current grid. */

  int i, dim, size = 1;
  for (dim = 0; dim < GridRank; dim++)
  {
    size *= GridDimension[dim];
  }

  /* Get easy to handle pointers for each variable. */

  float *density = BaryonField[DensNum];
  float *totalenergy = BaryonField[TENum];
  float *gasenergy = BaryonField[GENum];
  float *velocity1 = BaryonField[Vel1Num];
  float *velocity2 = BaryonField[Vel2Num];
  float *velocity3 = BaryonField[Vel3Num];

  /* Compute total gas energy if using MHD */
  if (UseMHD)
  {
    totalenergy = new float[size];
    float B2;
    for (int n = 0; n < size; n++)
    {
      B2 = pow(BaryonField[B1Num][n], 2) + pow(BaryonField[B2Num][n], 2) + pow(BaryonField[B3Num][n], 2);
      totalenergy[n] = BaryonField[TENum][n] - 0.5 * B2 / BaryonField[DensNum][n];
    }
  }
  else
  {
    totalenergy = BaryonField[TENum];
  }

  /* If using cosmology, compute the expansion factor and get units. */

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
        VelocityUnits = 1, TimeUnits = 1, aUnits = 1;

  if (ComovingCoordinates)
  {

    if (CosmologyComputeExpansionFactor(Time + 0.5 * dtFixed, &a, &dadt) == FAIL)
    {
      ENZO_FAIL("Error in CosmologyComputeExpansionFactors.");
    }

    aUnits = 1.0 / (1.0 + InitialRedshift);
  }

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, Time) == FAIL)
  {
    ENZO_FAIL("Error in GetUnits.");
  }

  float afloat = float(a);

  /* Metal cooling codes. */

  int MetalNum = 0, SNColourNum = 0;
  int MetalFieldPresent = FALSE;

  // First see if there's a metal field (so we can conserve species in
  // the solver)
  MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields);
  SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields);
  MetalFieldPresent = (MetalNum != -1 || SNColourNum != -1);

  // Double check if there's a metal field when we have metal cooling
  if (MetalCooling && MetalFieldPresent == FALSE)
  {
    if (debug)
      fprintf(stderr, "Warning: No metal field found.  Turning OFF MetalCooling.\n");
    MetalCooling = FALSE;
    MetalNum = 0;
  }

  /* If both metal fields (Pop I/II and III) exist, create a field
     that contains their sum */

  float *MetalPointer;
  float *TotalMetals = NULL;

  if (MetalNum != -1 && SNColourNum != -1)
  {
    TotalMetals = new float[size];
    for (i = 0; i < size; i++)
      TotalMetals[i] = BaryonField[MetalNum][i] + BaryonField[SNColourNum][i];
    MetalPointer = TotalMetals;
  } // ENDIF both metal types
  else
  {
    if (MetalNum != -1)
      MetalPointer = BaryonField[MetalNum];
    else if (SNColourNum != -1)
      MetalPointer = BaryonField[SNColourNum];
  } // ENDELSE both metal types

  /* Calculate the rates due to the radiation field. */

  if (!GadgetEquilibriumCooling)
  {
    if (RadiationFieldCalculateRates(Time + 0.5 * dtFixed) == FAIL)
    {
      ENZO_FAIL("Error in RadiationFieldCalculateRates.");
    }
  }

  /* Set up information for rates which depend on the radiation field. 
     Precompute factors for self shielding (this is the cross section * dx). */

  float HIShieldFactor = RadiationData.HIAveragePhotoHeatingCrossSection *
                         double(LengthUnits) * CellWidth[0][0];
  float HeIShieldFactor = RadiationData.HeIAveragePhotoHeatingCrossSection *
                          double(LengthUnits) * CellWidth[0][0];
  float HeIIShieldFactor = RadiationData.HeIIAveragePhotoHeatingCrossSection *
                           double(LengthUnits) * CellWidth[0][0];

  float dtCool = dtFixed;
#ifdef TRANSFER
  dtCool = (RTCoupledSolverIntermediateStep == TRUE) ? dtPhoton : dtFixed;
#endif

  /* Call the fortran routine to solve cooling equations. */

#ifdef USE_KROME

  if (KromeDt > 0.0)
    dtCool = KromeDt;

  if (use_kromestep == 3)
  {
    dtCool = Time - KromeTime;
    if (debug && MyProcessorNumber == ROOT_PROCESSOR)
    {
      printf("SolveRateAndCool, KromeDt = %13.7e\n", dtCool * TimeUnits);
    }
  }

  if (use_kromestep == 0 || use_kromestep == 2 ||
      (use_kromestep == 3 && TopGridCycle == KromeCycle) ||
      (use_kromestep == 1 && Time > KromeTime + KromeDt))
  {

    FORTRAN_NAME(krome_driver)
    (
        density, totalenergy, gasenergy, velocity1, velocity2, velocity3,
        BaryonField[DeNum], BaryonField[CHINum], BaryonField[OINum],
        BaryonField[HNCINum], BaryonField[HCNINum],
        BaryonField[H2INum], BaryonField[CINum],
        BaryonField[HINum], BaryonField[H2OINum],
        BaryonField[OHINum], BaryonField[O2INum],
        BaryonField[CH2INum], BaryonField[H2COINum],
        BaryonField[HCOINum], BaryonField[MGINum],
        BaryonField[NH3INum], BaryonField[NOINum],
        BaryonField[SIINum], BaryonField[SIC2INum],
        BaryonField[SIC3INum], BaryonField[SICINum],
        BaryonField[SIH2INum], BaryonField[SIH3INum],
        BaryonField[CNINum], BaryonField[COINum],
        BaryonField[N2INum], BaryonField[NH2INum],
        BaryonField[CH3INum], BaryonField[CH4INum],
        BaryonField[NINum], BaryonField[NHINum],
        BaryonField[SIH4INum], BaryonField[SIHINum],
        BaryonField[SIOINum], BaryonField[HeINum],
        BaryonField[HNOINum], BaryonField[CH3OHINum],
        BaryonField[CO2INum], BaryonField[H2CNINum],
        BaryonField[H2SIOINum], BaryonField[HNCOINum],
        BaryonField[NO2INum], BaryonField[O2HINum],
        BaryonField[OCNINum], BaryonField[CH3OH_DUSTINum],
        BaryonField[HNCO_DUSTINum], BaryonField[H2CO_DUSTINum],
        BaryonField[SIH4_DUSTINum], BaryonField[H2SIO_DUSTINum],
        BaryonField[SIC_DUSTINum], BaryonField[SIC2_DUSTINum],
        BaryonField[SIC3_DUSTINum], BaryonField[CH4_DUSTINum],
        BaryonField[CO_DUSTINum], BaryonField[H2O_DUSTINum],
        BaryonField[NO_DUSTINum], BaryonField[CO2_DUSTINum],
        BaryonField[N2_DUSTINum], BaryonField[HCN_DUSTINum],
        BaryonField[NH3_DUSTINum], BaryonField[O2_DUSTINum],
        BaryonField[NO2_DUSTINum], BaryonField[HNO_DUSTINum],
        BaryonField[O2H_DUSTINum], BaryonField[H2CN_DUSTINum],
        BaryonField[MG_DUSTINum], BaryonField[HNC_DUSTINum],
        BaryonField[E_DUSTINum], BaryonField[SIO_DUSTINum],
        BaryonField[HCOIINum], BaryonField[HIINum],
        BaryonField[HOCIINum], BaryonField[CIINum],
        BaryonField[CH2IINum], BaryonField[CHIINum],
        BaryonField[H2COIINum], BaryonField[MGIINum],
        BaryonField[NH3IINum], BaryonField[NOIINum],
        BaryonField[SIIINum], BaryonField[SIC2IINum],
        BaryonField[SIC3IINum], BaryonField[SICIINum],
        BaryonField[SIH2IINum], BaryonField[SIH3IINum],
        BaryonField[CNIINum], BaryonField[COIINum],
        BaryonField[N2IINum], BaryonField[O2IINum],
        BaryonField[H2OIINum], BaryonField[NH2IINum],
        BaryonField[OIINum], BaryonField[OHIINum],
        BaryonField[CH3IINum], BaryonField[CH4IINum],
        BaryonField[NIINum], BaryonField[HCNIINum],
        BaryonField[NHIINum], BaryonField[SIH4IINum],
        BaryonField[SIHIINum], BaryonField[SIOIINum],
        BaryonField[H2IINum], BaryonField[HeIINum],
        BaryonField[HNOIINum], BaryonField[H2NOIINum],
        BaryonField[H3IINum], BaryonField[H3COIINum],
        BaryonField[H3OIINum], BaryonField[HCNHIINum],
        BaryonField[HCO2IINum], BaryonField[HeHIINum],
        BaryonField[N2HIINum], BaryonField[O2HIINum],
        BaryonField[SIH5IINum], BaryonField[SIOHIINum],

        GridDimension, GridDimension + 1, GridDimension + 2,
        &HydroMethod,
        &DualEnergyFormalism,
        &GridRank, GridStartIndex, GridStartIndex + 1, GridStartIndex + 2,
        GridEndIndex, GridEndIndex + 1, GridEndIndex + 2,
        &dtCool, &afloat,
        &TemperatureUnits, &LengthUnits, &aUnits, &DensityUnits, &TimeUnits,
        &Gamma,
        &CoolData.HydrogenFractionByMass, &CoolData.DeuteriumToHydrogenRatio, 
        CellWidth[0], &use_sputtering);
  }
#endif // USE_KROME

  /*
  int ierr = 0;
  int addRT = (RadiativeTransfer) || (RadiativeTransferFLD);
  int RTcoupled = RadiativeTransferCoupledRateSolver;
  if ((RadiativeTransferFLD) && (RadiativeTransfer==0))
    RTcoupled = 0;    // disable if using FLD and not ray-tracing

  FORTRAN_NAME(solve_rate_cool)(
    density, totalenergy, gasenergy, velocity1, velocity2, velocity3,
    BaryonField[DeNum], BaryonField[HINum], BaryonField[HIINum], 
    BaryonField[HeINum], BaryonField[HeIINum], BaryonField[HeIIINum], 
    GridDimension, GridDimension+1, GridDimension+2,
    &CoolData.NumberOfTemperatureBins, &ComovingCoordinates, &HydroMethod, 
    &DualEnergyFormalism, &MultiSpecies, &MetalFieldPresent, &MetalCooling, 
    &H2FormationOnDust, 
    &GridRank, GridStartIndex, GridStartIndex+1, GridStartIndex+2, 
    GridEndIndex, GridEndIndex+1, GridEndIndex+2,
    &CoolData.ih2co, &CoolData.ipiht, &PhotoelectricHeating,
    CellWidth[0], &dtCool, &afloat, &RadiationFieldRedshift, 
    &CoolData.TemperatureStart, &CoolData.TemperatureEnd,
    &TemperatureUnits, &LengthUnits, &aUnits, &DensityUnits, &TimeUnits,
    &DualEnergyFormalismEta1, &DualEnergyFormalismEta2, &Gamma,
    &CoolData.HydrogenFractionByMass, &CoolData.DeuteriumToHydrogenRatio,
    &CoolData.SolarMetalFractionByMass,
    RateData.k1, RateData.k2, RateData.k3, RateData.k4, RateData.k5, 
    RateData.k6, RateData.k7, RateData.k8, RateData.k9, RateData.k10,
    RateData.k11, RateData.k12, RateData.k13, RateData.k13dd, RateData.k14, 
    RateData.k15, RateData.k16,
    RateData.k17, RateData.k18, RateData.k19, RateData.k22,
    &RateData.k24, &RateData.k25, &RateData.k26, &RateData.k27,
    &RateData.k28, &RateData.k29, &RateData.k30, &RateData.k31,
    RateData.k50, RateData.k51, RateData.k52, RateData.k53,
    RateData.k54, RateData.k55, RateData.k56,
    &RateData.NumberOfDustTemperatureBins, &RateData.DustTemperatureStart, 
    &RateData.DustTemperatureEnd, RateData.h2dust, 
    RateData.n_cr_n, RateData.n_cr_d1, RateData.n_cr_d2,
    CoolData.ceHI, CoolData.ceHeI, CoolData.ceHeII, CoolData.ciHI,
    CoolData.ciHeI, 
    CoolData.ciHeIS, CoolData.ciHeII, CoolData.reHII, CoolData.reHeII1, 
    CoolData.reHeII2, CoolData.reHeIII, CoolData.brem, &CoolData.comp, &CoolData.gammah,
    &CoolData.comp_xray, &CoolData.temp_xray,
    &CoolData.piHI, &CoolData.piHeI, &CoolData.piHeII,
    BaryonField[HMNum], BaryonField[H2INum], BaryonField[H2IINum],
    BaryonField[DINum], BaryonField[DIINum], BaryonField[HDINum],
    MetalPointer,
    CoolData.hyd01k, CoolData.h2k01, CoolData.vibh, CoolData.roth,CoolData.rotl,
    CoolData.GP99LowDensityLimit, CoolData.GP99HighDensityLimit, 
    CoolData.HDlte, CoolData.HDlow,
    CoolData.GAHI, CoolData.GAH2, CoolData.GAHe, CoolData.GAHp,
    CoolData.GAel, CoolData.gas_grain, 
    CoolData.metals, &CoolData.NumberOfElectronFracBins, 
    &CoolData.ElectronFracStart, &CoolData.ElectronFracEnd,
    RadiationData.Spectrum[0], &RadiationFieldType, 
    &RadiationData.NumberOfFrequencyBins, 
    &RadiationFieldRecomputeMetalRates,
    &RadiationData.RadiationShield, &HIShieldFactor, &HeIShieldFactor, &HeIIShieldFactor,
    &addRT, &RTcoupled,
    &RTCoupledSolverIntermediateStep, &ierr,
    &RadiativeTransferHydrogenOnly,
    BaryonField[kphHINum], BaryonField[kphHeINum], BaryonField[kphHeIINum], 
    BaryonField[kdissH2INum], BaryonField[gammaNum],
    &H2OpticalDepthApproximation, &CIECooling, &ThreeBodyRate, CoolData.cieco,
    &CloudyCoolingData.CMBTemperatureFloor,
    &CloudyCoolingData.IncludeCloudyHeating,
    &CloudyCoolingData.CloudyElectronFractionFactor,
    &CloudyCoolingData.CloudyCoolingGridRank,
    CloudyCoolingData.CloudyCoolingGridDimension,
    CloudyCoolingData.CloudyCoolingGridParameters[0],
    CloudyCoolingData.CloudyCoolingGridParameters[1],
    CloudyCoolingData.CloudyCoolingGridParameters[2],
    CloudyCoolingData.CloudyCoolingGridParameters[3],
    CloudyCoolingData.CloudyCoolingGridParameters[4],
    &CloudyCoolingData.CloudyDataSize,
    CloudyCoolingData.CloudyCooling, CloudyCoolingData.CloudyHeating);

  if (ierr) {
      fprintf(stdout, "GridLeftEdge = %"FSYM" %"FSYM" %"FSYM"\n",
	      GridLeftEdge[0], GridLeftEdge[1], GridLeftEdge[2]);
      fprintf(stdout, "GridRightEdge = %"FSYM" %"FSYM" %"FSYM"\n",
	      GridRightEdge[0], GridRightEdge[1], GridRightEdge[2]);
      fprintf(stdout, "GridDimension = %"ISYM" %"ISYM" %"ISYM"\n",
	      GridDimension[0], GridDimension[1], GridDimension[2]);
      ENZO_FAIL("Error in FORTRAN rate/cool solver!\n");
  }
*/

  if (UseMHD)
  {
    float B2, v2;
    for (int n = 0; n < size; n++)
    {
      B2 = pow(BaryonField[B1Num][n], 2) + pow(BaryonField[B2Num][n], 2) + pow(BaryonField[B3Num][n], 2);

      /* Always trust gas energy in cooling routine */
      if (DualEnergyFormalism)
      {

        v2 = pow(BaryonField[Vel1Num][n], 2) +
             pow(BaryonField[Vel2Num][n], 2) + pow(BaryonField[Vel3Num][n], 2);
        BaryonField[TENum][n] = gasenergy[n] + 0.5 * v2 + 0.5 * B2 / BaryonField[DensNum][n];
      }
      else
      {
        BaryonField[TENum][n] = totalenergy[n] + 0.5 * B2 / BaryonField[DensNum][n];
      }
    }

    delete totalenergy;
  }

  delete[] TotalMetals;

  return SUCCESS;
}
