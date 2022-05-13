/***********************************************************************
/
/  GRID CLASS (WRAP THE NAUNET CHEMISTRY SOLVER)
/
/  written by: Chia-Jung Hsu
/  date:       2021
/  modified1:
/
/  PURPOSE: Solve chemistry with naunet.
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
// clang-format off
#include "preincludes.h"
#ifdef USE_NAUNET
#include "naunet_enzo.h"
#endif
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "phys_constants.h"

// function prototypes

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);

double AvTableX[91] = {
  -3.0, -2.9, -2.8, -2.7, -2.6, -2.5, -2.4,
  -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7,
  -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0,
  -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3,
  -0.2, -0.1,  0.0,  0.1,  0.2,  0.3,  0.4,
  0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.1,
  1.2,  1.3,  1.4,  1.5,  1.6,  1.7,  1.8,
  1.9,  2.0,  2.1,  2.2,  2.3,  2.4,  2.5,
  2.6,  2.7,  2.8,  2.9,  3.0,  3.1,  3.2,
  3.3,  3.4,  3.5,  3.6,  3.7,  3.8,  3.9,
  4.0,  4.1,  4.2,  4.3,  4.4,  4.5,  4.6,
  4.7,  4.8,  4.9,  5.0,  5.1,  5.2,  5.3,
  5.4,  5.5,  5.6,  5.7,  5.8,  5.9,  6.0
};

double AvTableY[91] = {
  3.85e-03, 3.90e-03, 3.97e-03, 4.06e-03, 4.18e-03, 4.32e-03, 4.48e-03, 4.69e-03,
  4.96e-03, 5.29e-03, 5.71e-03, 6.22e-03, 6.85e-03, 7.65e-03, 8.65e-03, 9.91e-03,
  1.13e-02, 1.30e-02, 1.51e-02, 1.78e-02, 2.12e-02, 2.33e-02, 2.60e-02, 2.93e-02,
  3.36e-02, 3.89e-02, 4.41e-02, 5.07e-02, 5.90e-02, 6.95e-02, 8.27e-02, 9.47e-02,
  1.10e-01, 1.29e-01, 1.53e-01, 1.83e-01, 2.02e-01, 2.25e-01, 2.55e-01, 2.92e-01,
  3.39e-01, 3.76e-01, 4.22e-01, 4.81e-01, 5.55e-01, 6.48e-01, 7.06e-01, 7.78e-01,
  8.68e-01, 9.83e-01, 1.13e+00, 1.21e+00, 1.31e+00, 1.44e+00, 1.60e+00, 1.80e+00,
  1.93e+00, 2.08e+00, 2.28e+00, 2.54e+00, 2.85e+00, 3.04e+00, 3.28e+00, 3.58e+00,
  3.95e+00, 4.42e+00, 4.70e+00, 5.05e+00, 5.48e+00, 6.03e+00, 6.72e+00, 7.12e+00,
  7.61e+00, 8.24e+00, 9.03e+00, 1.00e+01, 1.06e+01, 1.13e+01, 1.22e+01, 1.33e+01,
  1.47e+01, 1.54e+01, 1.64e+01, 1.76e+01, 1.92e+01, 2.11e+01, 2.22e+01, 2.35e+01,
  2.52e+01, 2.73e+01, 3.00e+01
};

double GetAv(double nh) {

  double lognh = log10(nh);

  if (lognh >= 6.0) {
    return 30.0;
  }
  else if (lognh < -3.0) {
    return 3.85e-03;
  }

  for (int i=0; i<91; i++) {
    if (lognh > AvTableX[i] && lognh <= AvTableX[i+1]) {
      return (AvTableY[i+1] * (lognh - AvTableX[i]) 
              + AvTableY[i] * (AvTableX[i+1] - lognh)) 
              / (AvTableX[i+1] - AvTableX[i]);
    }
  }

}

// double GetHNuclei(double *y) {
//     return 4.0e+00*y[IDX_GCH3OHI] + 4.0e+00*y[IDX_GCH4I] + 2.0e+00*y[IDX_GH2CNI] +
//         2.0e+00*y[IDX_GH2COI] + 2.0e+00*y[IDX_GH2OI] + 2.0e+00*y[IDX_GH2SiOI] +
//         1.0e+00*y[IDX_GHCNI] + 1.0e+00*y[IDX_GHNCI] + 1.0e+00*y[IDX_GHNCOI] +
//         1.0e+00*y[IDX_GHNOI] + 3.0e+00*y[IDX_GNH3I] + 1.0e+00*y[IDX_GO2HI] +
//         4.0e+00*y[IDX_GSiH4I] + 1.0e+00*y[IDX_CHI] + 1.0e+00*y[IDX_CHII] +
//         2.0e+00*y[IDX_CH2I] + 2.0e+00*y[IDX_CH2II] + 3.0e+00*y[IDX_CH3I] +
//         3.0e+00*y[IDX_CH3II] + 4.0e+00*y[IDX_CH3OHI] + 4.0e+00*y[IDX_CH4I] +
//         4.0e+00*y[IDX_CH4II] + 1.0e+00*y[IDX_HI] + 1.0e+00*y[IDX_HII] +
//         2.0e+00*y[IDX_H2I] + 2.0e+00*y[IDX_H2II] + 2.0e+00*y[IDX_H2CNI] +
//         2.0e+00*y[IDX_H2COI] + 2.0e+00*y[IDX_H2COII] + 2.0e+00*y[IDX_H2NOII] +
//         2.0e+00*y[IDX_H2OI] + 2.0e+00*y[IDX_H2OII] + 2.0e+00*y[IDX_H2SiOI] +
//         3.0e+00*y[IDX_H3II] + 3.0e+00*y[IDX_H3COII] + 3.0e+00*y[IDX_H3OII] +
//         1.0e+00*y[IDX_HCNI] + 1.0e+00*y[IDX_HCNII] + 2.0e+00*y[IDX_HCNHII] +
//         1.0e+00*y[IDX_HCOI] + 1.0e+00*y[IDX_HCOII] + 1.0e+00*y[IDX_HCO2II] +
//         1.0e+00*y[IDX_HeHII] + 1.0e+00*y[IDX_HNCI] + 1.0e+00*y[IDX_HNCOI] +
//         1.0e+00*y[IDX_HNOI] + 1.0e+00*y[IDX_HNOII] + 1.0e+00*y[IDX_HOCII] +
//         1.0e+00*y[IDX_N2HII] + 1.0e+00*y[IDX_NHI] + 1.0e+00*y[IDX_NHII] +
//         2.0e+00*y[IDX_NH2I] + 2.0e+00*y[IDX_NH2II] + 3.0e+00*y[IDX_NH3I] +
//         3.0e+00*y[IDX_NH3II] + 1.0e+00*y[IDX_O2HI] + 1.0e+00*y[IDX_O2HII] +
//         1.0e+00*y[IDX_OHI] + 1.0e+00*y[IDX_OHII] + 1.0e+00*y[IDX_SiHI] +
//         1.0e+00*y[IDX_SiHII] + 2.0e+00*y[IDX_SiH2I] + 2.0e+00*y[IDX_SiH2II] +
//         3.0e+00*y[IDX_SiH3I] + 3.0e+00*y[IDX_SiH3II] + 4.0e+00*y[IDX_SiH4I] +
//         4.0e+00*y[IDX_SiH4II] + 5.0e+00*y[IDX_SiH5II] + 1.0e+00*y[IDX_SiOHII];
// }

int grid::NaunetWrapper()
{

#ifdef USE_NAUNET

  if (use_naunet == FALSE)
    return SUCCESS;

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (use_naunetstep == 1 && (TopGridCycle != NaunetCycle || TopGridCycle == 0))
    return SUCCESS;

  if (MultiSpecies != NAUNET_SPECIES) {
    printf("NaunetWrapper Warning: MultiSpecies = %d isn't valid for naunet. \
            Skip solving chemistry.\n", MultiSpecies);
    return SUCCESS;
  }

  if (debug && MyProcessorNumber == ROOT_PROCESSOR) {
    if (use_naunetstep == 1 && TopGridCycle == NaunetCycle) {
      printf("NaunetWrapper: GridLevel=%d, TopGridCycle=%d, NaunetCycle=%d, NaunetCycleSkip:%d\n",
              GridLevel, TopGridCycle, NaunetCycle, NaunetCycleSkip);
    }
  }

  LCAPERF_START("grid_NaunetWrapper");

  int GCH3OHINum, GCH4INum, GCOINum, GCO2INum, GH2CNINum, GH2COINum, GH2OINum,
      GH2SiOINum, GHCNINum, GHNCINum, GHNCOINum, GHNOINum, GMgINum, GN2INum,
      GNH3INum, GNOINum, GNO2INum, GO2INum, GO2HINum, GSiCINum, GSiC2INum,
      GSiC3INum, GSiH4INum, GSiOINum, CINum, CIINum, CHINum, CHIINum, CH2INum,
      CH2IINum, CH3INum, CH3IINum, CH3OHINum, CH4INum, CH4IINum, CNINum,
      CNIINum, COINum, COIINum, CO2INum, DeNum, HINum, HIINum, H2INum, H2IINum,
      H2CNINum, H2COINum, H2COIINum, H2NOIINum, H2OINum, H2OIINum, H2SiOINum,
      H3IINum, H3COIINum, H3OIINum, HCNINum, HCNIINum, HCNHIINum, HCOINum,
      HCOIINum, HCO2IINum, HeINum, HeIINum, HeHIINum, HNCINum, HNCOINum,
      HNOINum, HNOIINum, HOCIINum, MgINum, MgIINum, NINum, NIINum, N2INum,
      N2IINum, N2HIINum, NHINum, NHIINum, NH2INum, NH2IINum, NH3INum, NH3IINum,
      NOINum, NOIINum, NO2INum, OINum, OIINum, O2INum, O2IINum, O2HINum,
      O2HIINum, OCNINum, OHINum, OHIINum, SiINum, SiIINum, SiCINum, SiCIINum,
      SiC2INum, SiC2IINum, SiC3INum, SiC3IINum, SiHINum, SiHIINum, SiH2INum,
      SiH2IINum, SiH3INum, SiH3IINum, SiH4INum, SiH4IINum, SiH5IINum, SiOINum,
      SiOIINum, SiOHIINum;

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;

  double dt_chem = dtFixed;

  if (use_naunetstep == 1) dt_chem = Time - NaunetTime;

  if (MyProcessorNumber == ROOT_PROCESSOR) printf("NaunetWrapper: dt_chem = %13.7e\n", dt_chem);
  
  // Compute the size of the fields.
 
  int i, j, k, igrid;
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  int activesize = 1;
  for (int dim = 0; dim < GridRank; dim++)
    activesize *= (GridDimension[dim] - 2*NumberOfGhostZones);

  Eint32 *g_grid_dimension, *g_grid_start, *g_grid_end;
  g_grid_dimension = new Eint32[GridRank];
  g_grid_start = new Eint32[GridRank];
  g_grid_end = new Eint32[GridRank];
  for (i = 0; i < GridRank; i++) {
    g_grid_dimension[i] = (Eint32) GridDimension[i];
    g_grid_start[i] = (Eint32) GridStartIndex[i];
    g_grid_end[i] = (Eint32) GridEndIndex[i];
  }
 
  // Find fields: density, total energy, velocity1-3.
 
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  // Find Multi-species fields.

  GCH3OHINum = GCH4INum = GCOINum = GCO2INum = GH2CNINum = GH2COINum =
    GH2OINum = GH2SiOINum = GHCNINum = GHNCINum = GHNCOINum = GHNOINum = GMgINum
    = GN2INum = GNH3INum = GNOINum = GNO2INum = GO2INum = GO2HINum = GSiCINum =
    GSiC2INum = GSiC3INum = GSiH4INum = GSiOINum = CINum = CIINum = CHINum =
    CHIINum = CH2INum = CH2IINum = CH3INum = CH3IINum = CH3OHINum = CH4INum =
    CH4IINum = CNINum = CNIINum = COINum = COIINum = CO2INum = DeNum = HINum =
    HIINum = H2INum = H2IINum = H2CNINum = H2COINum = H2COIINum = H2NOIINum =
    H2OINum = H2OIINum = H2SiOINum = H3IINum = H3COIINum = H3OIINum = HCNINum =
    HCNIINum = HCNHIINum = HCOINum = HCOIINum = HCO2IINum = HeINum = HeIINum =
    HeHIINum = HNCINum = HNCOINum = HNOINum = HNOIINum = HOCIINum = MgINum =
    MgIINum = NINum = NIINum = N2INum = N2IINum = N2HIINum = NHINum = NHIINum =
    NH2INum = NH2IINum = NH3INum = NH3IINum = NOINum = NOIINum = NO2INum = OINum
    = OIINum = O2INum = O2IINum = O2HINum = O2HIINum = OCNINum = OHINum =
    OHIINum = SiINum = SiIINum = SiCINum = SiCIINum = SiC2INum = SiC2IINum =
    SiC3INum = SiC3IINum = SiHINum = SiHIINum = SiH2INum = SiH2IINum = SiH3INum
    = SiH3IINum = SiH4INum = SiH4IINum = SiH5IINum = SiOINum = SiOIINum =
    SiOHIINum = 0;
 
  if (MultiSpecies == NAUNET_SPECIES)
    if (IdentifyNaunetSpeciesFields(GCH3OHINum, GCH4INum, GCOINum, GCO2INum,
                                    GH2CNINum, GH2COINum, GH2OINum, GH2SiOINum,
                                    GHCNINum, GHNCINum, GHNCOINum, GHNOINum,
                                    GMgINum, GN2INum, GNH3INum, GNOINum,
                                    GNO2INum, GO2INum, GO2HINum, GSiCINum,
                                    GSiC2INum, GSiC3INum, GSiH4INum, GSiOINum,
                                    CINum, CIINum, CHINum, CHIINum, CH2INum,
                                    CH2IINum, CH3INum, CH3IINum, CH3OHINum,
                                    CH4INum, CH4IINum, CNINum, CNIINum, COINum,
                                    COIINum, CO2INum, DeNum, HINum, HIINum,
                                    H2INum, H2IINum, H2CNINum, H2COINum,
                                    H2COIINum, H2NOIINum, H2OINum, H2OIINum,
                                    H2SiOINum, H3IINum, H3COIINum, H3OIINum,
                                    HCNINum, HCNIINum, HCNHIINum, HCOINum,
                                    HCOIINum, HCO2IINum, HeINum, HeIINum,
                                    HeHIINum, HNCINum, HNCOINum, HNOINum,
                                    HNOIINum, HOCIINum, MgINum, MgIINum, NINum,
                                    NIINum, N2INum, N2IINum, N2HIINum, NHINum,
                                    NHIINum, NH2INum, NH2IINum, NH3INum,
                                    NH3IINum, NOINum, NOIINum, NO2INum, OINum,
                                    OIINum, O2INum, O2IINum, O2HINum, O2HIINum,
                                    OCNINum, OHINum, OHIINum, SiINum, SiIINum,
                                    SiCINum, SiCIINum, SiC2INum, SiC2IINum,
                                    SiC3INum, SiC3IINum, SiHINum, SiHIINum,
                                    SiH2INum, SiH2IINum, SiH3INum, SiH3IINum,
                                    SiH4INum, SiH4IINum, SiH5IINum, SiOINum,
                                    SiOIINum, SiOHIINum) == FAIL) {
      ENZO_FAIL("Error in grid->IdentifyNaunetSpeciesFields.\n");
    }
 
  // Get easy to handle pointers for each variable.
 
  float *density     = BaryonField[DensNum];
  float *totalenergy = BaryonField[TENum];
  float *gasenergy   = BaryonField[GENum];
  float *velocity1   = BaryonField[Vel1Num];
  float *velocity2   = BaryonField[Vel2Num];
  float *velocity3   = BaryonField[Vel3Num];

  // Compute the cooling time.

  FLOAT a = 1.0, dadt;
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
        VelocityUnits = 1, TimeUnits = 1, aUnits = 1;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
           &TimeUnits, &VelocityUnits, Time);

  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(Time+0.5*dt_chem, &a, &dadt);
    aUnits = 1.0/(1.0 + InitialRedshift);
  } 
  else if (RadiationFieldRedshift > -1){
    a        = 1.0 / (1.0 + RadiationFieldRedshift);
    aUnits   = 1.0;
  }
  float afloat = float(a);

  float *temperature = new float[size]; 
  if (this->ComputeTemperatureField(temperature) == FAIL){
    ENZO_FAIL("Error in grid->ComputeTemperatureField.");
  }

  float *sputtering = new float[size];
  if (opt_sputtering) {
    if (this->ComputeSputteringRate(sputtering) == FAIL){
      ENZO_FAIL("Error in grid->ComputeSputteringRate.");
    }
  }
  else {
    for (i = 0; i < size; i++) {
      sputtering[i] = 0.0;
    }
  }

  float NumberDensityUnits = DensityUnits / mh;

  float atol = 1e-30, smallx = 1e-40;

  Naunet naunet;
  naunet.Init(1, atol, 1e-5, 1000);

  // TODO: comoving, heating/cooling
  
  // Set your parameters here
  NaunetData data;

  realtype y[NAUNET_NEQUATIONS], y_init[NAUNET_NEQUATIONS];

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      igrid = (k * GridDimension[1] + j) * GridDimension[0] + GridStartIndex[0];
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, igrid++) {

        double nH       = BaryonField[iden][igrid] * DensityUnits / (1.4 * mh);
        // printf("nH %13.7e, temperature: %13.7e, time: %13.7e\n", data.nH, temperature[igrid], dt_chem * TimeUnits / 86400.0 / 365.0);
        // data.Tgas     = 15.0;
        data.nH       = nH;
        // data.Tgas     = min(temperature[igrid], 300.0);
        data.Tgas     = temperature[igrid];
        data.zeta     = 1.0e-16;
        data.Av       = GetAv(nH);
        data.omega    = 0.5;
        data.G0       = 4.0;
        data.rG       = 1.0e-5;
        data.gdens    = 7.6394373e-13 * data.nH;
        data.sites    = 1.5e15;
        data.fr       = opt_freeze;
        data.opt_thd  = opt_thd;
        data.opt_crd  = opt_crd;
        data.opt_uvd  = opt_uvd;
        data.opt_h2d  = opt_h2d;
        data.eb_crd   = 2.00e3;
        data.eb_h2d   = 1.21e3;
        data.eb_uvd   = 1.00e4;
        data.crdeseff = 1.0e8;
        data.h2deseff = 1.0e-2;
        data.uvcreff  = 1.0e-3;
        data.ksp      = sputtering[igrid];
    
        y[IDX_GCH3OHI] = max(BaryonField[GCH3OHINum][igrid], smallx) * NumberDensityUnits / 32.0;
        y[IDX_GCH4I] = max(BaryonField[GCH4INum][igrid], smallx) * NumberDensityUnits / 16.0;
        y[IDX_GCOI] = max(BaryonField[GCOINum][igrid], smallx) * NumberDensityUnits / 28.0;
        y[IDX_GCO2I] = max(BaryonField[GCO2INum][igrid], smallx) * NumberDensityUnits / 44.0;
        y[IDX_GH2CNI] = max(BaryonField[GH2CNINum][igrid], smallx) * NumberDensityUnits / 28.0;
        y[IDX_GH2COI] = max(BaryonField[GH2COINum][igrid], smallx) * NumberDensityUnits / 30.0;
        y[IDX_GH2OI] = max(BaryonField[GH2OINum][igrid], smallx) * NumberDensityUnits / 18.0;
        y[IDX_GH2SiOI] = max(BaryonField[GH2SiOINum][igrid], smallx) * NumberDensityUnits / 46.0;
        y[IDX_GHCNI] = max(BaryonField[GHCNINum][igrid], smallx) * NumberDensityUnits / 27.0;
        y[IDX_GHNCI] = max(BaryonField[GHNCINum][igrid], smallx) * NumberDensityUnits / 27.0;
        y[IDX_GHNCOI] = max(BaryonField[GHNCOINum][igrid], smallx) * NumberDensityUnits / 43.0;
        y[IDX_GHNOI] = max(BaryonField[GHNOINum][igrid], smallx) * NumberDensityUnits / 31.0;
        y[IDX_GMgI] = max(BaryonField[GMgINum][igrid], smallx) * NumberDensityUnits / 24.0;
        y[IDX_GN2I] = max(BaryonField[GN2INum][igrid], smallx) * NumberDensityUnits / 28.0;
        y[IDX_GNH3I] = max(BaryonField[GNH3INum][igrid], smallx) * NumberDensityUnits / 17.0;
        y[IDX_GNOI] = max(BaryonField[GNOINum][igrid], smallx) * NumberDensityUnits / 30.0;
        y[IDX_GNO2I] = max(BaryonField[GNO2INum][igrid], smallx) * NumberDensityUnits / 46.0;
        y[IDX_GO2I] = max(BaryonField[GO2INum][igrid], smallx) * NumberDensityUnits / 32.0;
        y[IDX_GO2HI] = max(BaryonField[GO2HINum][igrid], smallx) * NumberDensityUnits / 33.0;
        y[IDX_GSiCI] = max(BaryonField[GSiCINum][igrid], smallx) * NumberDensityUnits / 40.0;
        y[IDX_GSiC2I] = max(BaryonField[GSiC2INum][igrid], smallx) * NumberDensityUnits / 52.0;
        y[IDX_GSiC3I] = max(BaryonField[GSiC3INum][igrid], smallx) * NumberDensityUnits / 64.0;
        y[IDX_GSiH4I] = max(BaryonField[GSiH4INum][igrid], smallx) * NumberDensityUnits / 32.0;
        y[IDX_GSiOI] = max(BaryonField[GSiOINum][igrid], smallx) * NumberDensityUnits / 44.0;
        y[IDX_CI] = max(BaryonField[CINum][igrid], smallx) * NumberDensityUnits / 12.0;
        y[IDX_CII] = max(BaryonField[CIINum][igrid], smallx) * NumberDensityUnits / 12.0;
        y[IDX_CHI] = max(BaryonField[CHINum][igrid], smallx) * NumberDensityUnits / 13.0;
        y[IDX_CHII] = max(BaryonField[CHIINum][igrid], smallx) * NumberDensityUnits / 13.0;
        y[IDX_CH2I] = max(BaryonField[CH2INum][igrid], smallx) * NumberDensityUnits / 14.0;
        y[IDX_CH2II] = max(BaryonField[CH2IINum][igrid], smallx) * NumberDensityUnits / 14.0;
        y[IDX_CH3I] = max(BaryonField[CH3INum][igrid], smallx) * NumberDensityUnits / 15.0;
        y[IDX_CH3II] = max(BaryonField[CH3IINum][igrid], smallx) * NumberDensityUnits / 15.0;
        y[IDX_CH3OHI] = max(BaryonField[CH3OHINum][igrid], smallx) * NumberDensityUnits / 32.0;
        y[IDX_CH4I] = max(BaryonField[CH4INum][igrid], smallx) * NumberDensityUnits / 16.0;
        y[IDX_CH4II] = max(BaryonField[CH4IINum][igrid], smallx) * NumberDensityUnits / 16.0;
        y[IDX_CNI] = max(BaryonField[CNINum][igrid], smallx) * NumberDensityUnits / 26.0;
        y[IDX_CNII] = max(BaryonField[CNIINum][igrid], smallx) * NumberDensityUnits / 26.0;
        y[IDX_COI] = max(BaryonField[COINum][igrid], smallx) * NumberDensityUnits / 28.0;
        y[IDX_COII] = max(BaryonField[COIINum][igrid], smallx) * NumberDensityUnits / 28.0;
        y[IDX_CO2I] = max(BaryonField[CO2INum][igrid], smallx) * NumberDensityUnits / 44.0;
        y[IDX_EM] = max(BaryonField[DeNum][igrid], smallx) * NumberDensityUnits / 1.0;
        y[IDX_HI] = max(BaryonField[HINum][igrid], smallx) * NumberDensityUnits / 1.0;
        y[IDX_HII] = max(BaryonField[HIINum][igrid], smallx) * NumberDensityUnits / 1.0;
        y[IDX_H2I] = max(BaryonField[H2INum][igrid], smallx) * NumberDensityUnits / 2.0;
        y[IDX_H2II] = max(BaryonField[H2IINum][igrid], smallx) * NumberDensityUnits / 2.0;
        y[IDX_H2CNI] = max(BaryonField[H2CNINum][igrid], smallx) * NumberDensityUnits / 28.0;
        y[IDX_H2COI] = max(BaryonField[H2COINum][igrid], smallx) * NumberDensityUnits / 30.0;
        y[IDX_H2COII] = max(BaryonField[H2COIINum][igrid], smallx) * NumberDensityUnits / 30.0;
        y[IDX_H2NOII] = max(BaryonField[H2NOIINum][igrid], smallx) * NumberDensityUnits / 32.0;
        y[IDX_H2OI] = max(BaryonField[H2OINum][igrid], smallx) * NumberDensityUnits / 18.0;
        y[IDX_H2OII] = max(BaryonField[H2OIINum][igrid], smallx) * NumberDensityUnits / 18.0;
        y[IDX_H2SiOI] = max(BaryonField[H2SiOINum][igrid], smallx) * NumberDensityUnits / 46.0;
        y[IDX_H3II] = max(BaryonField[H3IINum][igrid], smallx) * NumberDensityUnits / 3.0;
        y[IDX_H3COII] = max(BaryonField[H3COIINum][igrid], smallx) * NumberDensityUnits / 31.0;
        y[IDX_H3OII] = max(BaryonField[H3OIINum][igrid], smallx) * NumberDensityUnits / 19.0;
        y[IDX_HCNI] = max(BaryonField[HCNINum][igrid], smallx) * NumberDensityUnits / 27.0;
        y[IDX_HCNII] = max(BaryonField[HCNIINum][igrid], smallx) * NumberDensityUnits / 27.0;
        y[IDX_HCNHII] = max(BaryonField[HCNHIINum][igrid], smallx) * NumberDensityUnits / 28.0;
        y[IDX_HCOI] = max(BaryonField[HCOINum][igrid], smallx) * NumberDensityUnits / 29.0;
        y[IDX_HCOII] = max(BaryonField[HCOIINum][igrid], smallx) * NumberDensityUnits / 29.0;
        y[IDX_HCO2II] = max(BaryonField[HCO2IINum][igrid], smallx) * NumberDensityUnits / 45.0;
        y[IDX_HeI] = max(BaryonField[HeINum][igrid], smallx) * NumberDensityUnits / 4.0;
        y[IDX_HeII] = max(BaryonField[HeIINum][igrid], smallx) * NumberDensityUnits / 4.0;
        y[IDX_HeHII] = max(BaryonField[HeHIINum][igrid], smallx) * NumberDensityUnits / 5.0;
        y[IDX_HNCI] = max(BaryonField[HNCINum][igrid], smallx) * NumberDensityUnits / 27.0;
        y[IDX_HNCOI] = max(BaryonField[HNCOINum][igrid], smallx) * NumberDensityUnits / 43.0;
        y[IDX_HNOI] = max(BaryonField[HNOINum][igrid], smallx) * NumberDensityUnits / 31.0;
        y[IDX_HNOII] = max(BaryonField[HNOIINum][igrid], smallx) * NumberDensityUnits / 31.0;
        y[IDX_HOCII] = max(BaryonField[HOCIINum][igrid], smallx) * NumberDensityUnits / 29.0;
        y[IDX_MgI] = max(BaryonField[MgINum][igrid], smallx) * NumberDensityUnits / 24.0;
        y[IDX_MgII] = max(BaryonField[MgIINum][igrid], smallx) * NumberDensityUnits / 24.0;
        y[IDX_NI] = max(BaryonField[NINum][igrid], smallx) * NumberDensityUnits / 14.0;
        y[IDX_NII] = max(BaryonField[NIINum][igrid], smallx) * NumberDensityUnits / 14.0;
        y[IDX_N2I] = max(BaryonField[N2INum][igrid], smallx) * NumberDensityUnits / 28.0;
        y[IDX_N2II] = max(BaryonField[N2IINum][igrid], smallx) * NumberDensityUnits / 28.0;
        y[IDX_N2HII] = max(BaryonField[N2HIINum][igrid], smallx) * NumberDensityUnits / 29.0;
        y[IDX_NHI] = max(BaryonField[NHINum][igrid], smallx) * NumberDensityUnits / 15.0;
        y[IDX_NHII] = max(BaryonField[NHIINum][igrid], smallx) * NumberDensityUnits / 15.0;
        y[IDX_NH2I] = max(BaryonField[NH2INum][igrid], smallx) * NumberDensityUnits / 16.0;
        y[IDX_NH2II] = max(BaryonField[NH2IINum][igrid], smallx) * NumberDensityUnits / 16.0;
        y[IDX_NH3I] = max(BaryonField[NH3INum][igrid], smallx) * NumberDensityUnits / 17.0;
        y[IDX_NH3II] = max(BaryonField[NH3IINum][igrid], smallx) * NumberDensityUnits / 17.0;
        y[IDX_NOI] = max(BaryonField[NOINum][igrid], smallx) * NumberDensityUnits / 30.0;
        y[IDX_NOII] = max(BaryonField[NOIINum][igrid], smallx) * NumberDensityUnits / 30.0;
        y[IDX_NO2I] = max(BaryonField[NO2INum][igrid], smallx) * NumberDensityUnits / 46.0;
        y[IDX_OI] = max(BaryonField[OINum][igrid], smallx) * NumberDensityUnits / 16.0;
        y[IDX_OII] = max(BaryonField[OIINum][igrid], smallx) * NumberDensityUnits / 16.0;
        y[IDX_O2I] = max(BaryonField[O2INum][igrid], smallx) * NumberDensityUnits / 32.0;
        y[IDX_O2II] = max(BaryonField[O2IINum][igrid], smallx) * NumberDensityUnits / 32.0;
        y[IDX_O2HI] = max(BaryonField[O2HINum][igrid], smallx) * NumberDensityUnits / 33.0;
        y[IDX_O2HII] = max(BaryonField[O2HIINum][igrid], smallx) * NumberDensityUnits / 33.0;
        y[IDX_OCNI] = max(BaryonField[OCNINum][igrid], smallx) * NumberDensityUnits / 42.0;
        y[IDX_OHI] = max(BaryonField[OHINum][igrid], smallx) * NumberDensityUnits / 17.0;
        y[IDX_OHII] = max(BaryonField[OHIINum][igrid], smallx) * NumberDensityUnits / 17.0;
        y[IDX_SiI] = max(BaryonField[SiINum][igrid], smallx) * NumberDensityUnits / 28.0;
        y[IDX_SiII] = max(BaryonField[SiIINum][igrid], smallx) * NumberDensityUnits / 28.0;
        y[IDX_SiCI] = max(BaryonField[SiCINum][igrid], smallx) * NumberDensityUnits / 40.0;
        y[IDX_SiCII] = max(BaryonField[SiCIINum][igrid], smallx) * NumberDensityUnits / 40.0;
        y[IDX_SiC2I] = max(BaryonField[SiC2INum][igrid], smallx) * NumberDensityUnits / 52.0;
        y[IDX_SiC2II] = max(BaryonField[SiC2IINum][igrid], smallx) * NumberDensityUnits / 52.0;
        y[IDX_SiC3I] = max(BaryonField[SiC3INum][igrid], smallx) * NumberDensityUnits / 64.0;
        y[IDX_SiC3II] = max(BaryonField[SiC3IINum][igrid], smallx) * NumberDensityUnits / 64.0;
        y[IDX_SiHI] = max(BaryonField[SiHINum][igrid], smallx) * NumberDensityUnits / 29.0;
        y[IDX_SiHII] = max(BaryonField[SiHIINum][igrid], smallx) * NumberDensityUnits / 29.0;
        y[IDX_SiH2I] = max(BaryonField[SiH2INum][igrid], smallx) * NumberDensityUnits / 30.0;
        y[IDX_SiH2II] = max(BaryonField[SiH2IINum][igrid], smallx) * NumberDensityUnits / 30.0;
        y[IDX_SiH3I] = max(BaryonField[SiH3INum][igrid], smallx) * NumberDensityUnits / 31.0;
        y[IDX_SiH3II] = max(BaryonField[SiH3IINum][igrid], smallx) * NumberDensityUnits / 31.0;
        y[IDX_SiH4I] = max(BaryonField[SiH4INum][igrid], smallx) * NumberDensityUnits / 32.0;
        y[IDX_SiH4II] = max(BaryonField[SiH4IINum][igrid], smallx) * NumberDensityUnits / 32.0;
        y[IDX_SiH5II] = max(BaryonField[SiH5IINum][igrid], smallx) * NumberDensityUnits / 33.0;
        y[IDX_SiOI] = max(BaryonField[SiOINum][igrid], smallx) * NumberDensityUnits / 44.0;
        y[IDX_SiOII] = max(BaryonField[SiOIINum][igrid], smallx) * NumberDensityUnits / 44.0;
        y[IDX_SiOHII] = max(BaryonField[SiOHIINum][igrid], smallx) * NumberDensityUnits / 45.0;
    
        for (int idx = IDX_GCH3OHI; idx <= IDX_SiOHII; idx++) {
          y_init[idx] = y[idx];
        }
    
        int flag = naunet.Solve(y, dt_chem * TimeUnits, &data);
    
        if (flag == NAUNET_FAIL) {

          // Try varying atol to fix the problem
          for (int ntrial = 1; ntrial < 6; ntrial ++) {
            naunet.Reset(1, pow(10.0, ntrial) * atol, 1e-5, 1000);

            for (int idx = IDX_GCH3OHI; idx <= IDX_SiOHII; idx++) {
              y[idx] = y_init[idx];
            }

            flag = naunet.Solve(y, dt_chem * TimeUnits, &data);

            if (flag == NAUNET_SUCCESS) {
              naunet.Reset(1, atol, 1e-5, 1000);
              break;
            }
          }

          // Fail to fix up
          if (flag == NAUNET_FAIL) {
            naunet.Finalize();

            ENZO_FAIL("Naunet failed in NaunetWrapper.C!");
          }
        }

        BaryonField[GCH3OHINum][igrid] = max(y[IDX_GCH3OHI] * 32.0 / NumberDensityUnits, smallx);
        BaryonField[GCH4INum][igrid] = max(y[IDX_GCH4I] * 16.0 / NumberDensityUnits, smallx);
        BaryonField[GCOINum][igrid] = max(y[IDX_GCOI] * 28.0 / NumberDensityUnits, smallx);
        BaryonField[GCO2INum][igrid] = max(y[IDX_GCO2I] * 44.0 / NumberDensityUnits, smallx);
        BaryonField[GH2CNINum][igrid] = max(y[IDX_GH2CNI] * 28.0 / NumberDensityUnits, smallx);
        BaryonField[GH2COINum][igrid] = max(y[IDX_GH2COI] * 30.0 / NumberDensityUnits, smallx);
        BaryonField[GH2OINum][igrid] = max(y[IDX_GH2OI] * 18.0 / NumberDensityUnits, smallx);
        BaryonField[GH2SiOINum][igrid] = max(y[IDX_GH2SiOI] * 46.0 / NumberDensityUnits, smallx);
        BaryonField[GHCNINum][igrid] = max(y[IDX_GHCNI] * 27.0 / NumberDensityUnits, smallx);
        BaryonField[GHNCINum][igrid] = max(y[IDX_GHNCI] * 27.0 / NumberDensityUnits, smallx);
        BaryonField[GHNCOINum][igrid] = max(y[IDX_GHNCOI] * 43.0 / NumberDensityUnits, smallx);
        BaryonField[GHNOINum][igrid] = max(y[IDX_GHNOI] * 31.0 / NumberDensityUnits, smallx);
        BaryonField[GMgINum][igrid] = max(y[IDX_GMgI] * 24.0 / NumberDensityUnits, smallx);
        BaryonField[GN2INum][igrid] = max(y[IDX_GN2I] * 28.0 / NumberDensityUnits, smallx);
        BaryonField[GNH3INum][igrid] = max(y[IDX_GNH3I] * 17.0 / NumberDensityUnits, smallx);
        BaryonField[GNOINum][igrid] = max(y[IDX_GNOI] * 30.0 / NumberDensityUnits, smallx);
        BaryonField[GNO2INum][igrid] = max(y[IDX_GNO2I] * 46.0 / NumberDensityUnits, smallx);
        BaryonField[GO2INum][igrid] = max(y[IDX_GO2I] * 32.0 / NumberDensityUnits, smallx);
        BaryonField[GO2HINum][igrid] = max(y[IDX_GO2HI] * 33.0 / NumberDensityUnits, smallx);
        BaryonField[GSiCINum][igrid] = max(y[IDX_GSiCI] * 40.0 / NumberDensityUnits, smallx);
        BaryonField[GSiC2INum][igrid] = max(y[IDX_GSiC2I] * 52.0 / NumberDensityUnits, smallx);
        BaryonField[GSiC3INum][igrid] = max(y[IDX_GSiC3I] * 64.0 / NumberDensityUnits, smallx);
        BaryonField[GSiH4INum][igrid] = max(y[IDX_GSiH4I] * 32.0 / NumberDensityUnits, smallx);
        BaryonField[GSiOINum][igrid] = max(y[IDX_GSiOI] * 44.0 / NumberDensityUnits, smallx);
        BaryonField[CINum][igrid] = max(y[IDX_CI] * 12.0 / NumberDensityUnits, smallx);
        BaryonField[CIINum][igrid] = max(y[IDX_CII] * 12.0 / NumberDensityUnits, smallx);
        BaryonField[CHINum][igrid] = max(y[IDX_CHI] * 13.0 / NumberDensityUnits, smallx);
        BaryonField[CHIINum][igrid] = max(y[IDX_CHII] * 13.0 / NumberDensityUnits, smallx);
        BaryonField[CH2INum][igrid] = max(y[IDX_CH2I] * 14.0 / NumberDensityUnits, smallx);
        BaryonField[CH2IINum][igrid] = max(y[IDX_CH2II] * 14.0 / NumberDensityUnits, smallx);
        BaryonField[CH3INum][igrid] = max(y[IDX_CH3I] * 15.0 / NumberDensityUnits, smallx);
        BaryonField[CH3IINum][igrid] = max(y[IDX_CH3II] * 15.0 / NumberDensityUnits, smallx);
        BaryonField[CH3OHINum][igrid] = max(y[IDX_CH3OHI] * 32.0 / NumberDensityUnits, smallx);
        BaryonField[CH4INum][igrid] = max(y[IDX_CH4I] * 16.0 / NumberDensityUnits, smallx);
        BaryonField[CH4IINum][igrid] = max(y[IDX_CH4II] * 16.0 / NumberDensityUnits, smallx);
        BaryonField[CNINum][igrid] = max(y[IDX_CNI] * 26.0 / NumberDensityUnits, smallx);
        BaryonField[CNIINum][igrid] = max(y[IDX_CNII] * 26.0 / NumberDensityUnits, smallx);
        BaryonField[COINum][igrid] = max(y[IDX_COI] * 28.0 / NumberDensityUnits, smallx);
        BaryonField[COIINum][igrid] = max(y[IDX_COII] * 28.0 / NumberDensityUnits, smallx);
        BaryonField[CO2INum][igrid] = max(y[IDX_CO2I] * 44.0 / NumberDensityUnits, smallx);
        BaryonField[DeNum][igrid] = max(y[IDX_EM] * 1.0 / NumberDensityUnits, smallx);
        BaryonField[HINum][igrid] = max(y[IDX_HI] * 1.0 / NumberDensityUnits, smallx);
        BaryonField[HIINum][igrid] = max(y[IDX_HII] * 1.0 / NumberDensityUnits, smallx);
        BaryonField[H2INum][igrid] = max(y[IDX_H2I] * 2.0 / NumberDensityUnits, smallx);
        BaryonField[H2IINum][igrid] = max(y[IDX_H2II] * 2.0 / NumberDensityUnits, smallx);
        BaryonField[H2CNINum][igrid] = max(y[IDX_H2CNI] * 28.0 / NumberDensityUnits, smallx);
        BaryonField[H2COINum][igrid] = max(y[IDX_H2COI] * 30.0 / NumberDensityUnits, smallx);
        BaryonField[H2COIINum][igrid] = max(y[IDX_H2COII] * 30.0 / NumberDensityUnits, smallx);
        BaryonField[H2NOIINum][igrid] = max(y[IDX_H2NOII] * 32.0 / NumberDensityUnits, smallx);
        BaryonField[H2OINum][igrid] = max(y[IDX_H2OI] * 18.0 / NumberDensityUnits, smallx);
        BaryonField[H2OIINum][igrid] = max(y[IDX_H2OII] * 18.0 / NumberDensityUnits, smallx);
        BaryonField[H2SiOINum][igrid] = max(y[IDX_H2SiOI] * 46.0 / NumberDensityUnits, smallx);
        BaryonField[H3IINum][igrid] = max(y[IDX_H3II] * 3.0 / NumberDensityUnits, smallx);
        BaryonField[H3COIINum][igrid] = max(y[IDX_H3COII] * 31.0 / NumberDensityUnits, smallx);
        BaryonField[H3OIINum][igrid] = max(y[IDX_H3OII] * 19.0 / NumberDensityUnits, smallx);
        BaryonField[HCNINum][igrid] = max(y[IDX_HCNI] * 27.0 / NumberDensityUnits, smallx);
        BaryonField[HCNIINum][igrid] = max(y[IDX_HCNII] * 27.0 / NumberDensityUnits, smallx);
        BaryonField[HCNHIINum][igrid] = max(y[IDX_HCNHII] * 28.0 / NumberDensityUnits, smallx);
        BaryonField[HCOINum][igrid] = max(y[IDX_HCOI] * 29.0 / NumberDensityUnits, smallx);
        BaryonField[HCOIINum][igrid] = max(y[IDX_HCOII] * 29.0 / NumberDensityUnits, smallx);
        BaryonField[HCO2IINum][igrid] = max(y[IDX_HCO2II] * 45.0 / NumberDensityUnits, smallx);
        BaryonField[HeINum][igrid] = max(y[IDX_HeI] * 4.0 / NumberDensityUnits, smallx);
        BaryonField[HeIINum][igrid] = max(y[IDX_HeII] * 4.0 / NumberDensityUnits, smallx);
        BaryonField[HeHIINum][igrid] = max(y[IDX_HeHII] * 5.0 / NumberDensityUnits, smallx);
        BaryonField[HNCINum][igrid] = max(y[IDX_HNCI] * 27.0 / NumberDensityUnits, smallx);
        BaryonField[HNCOINum][igrid] = max(y[IDX_HNCOI] * 43.0 / NumberDensityUnits, smallx);
        BaryonField[HNOINum][igrid] = max(y[IDX_HNOI] * 31.0 / NumberDensityUnits, smallx);
        BaryonField[HNOIINum][igrid] = max(y[IDX_HNOII] * 31.0 / NumberDensityUnits, smallx);
        BaryonField[HOCIINum][igrid] = max(y[IDX_HOCII] * 29.0 / NumberDensityUnits, smallx);
        BaryonField[MgINum][igrid] = max(y[IDX_MgI] * 24.0 / NumberDensityUnits, smallx);
        BaryonField[MgIINum][igrid] = max(y[IDX_MgII] * 24.0 / NumberDensityUnits, smallx);
        BaryonField[NINum][igrid] = max(y[IDX_NI] * 14.0 / NumberDensityUnits, smallx);
        BaryonField[NIINum][igrid] = max(y[IDX_NII] * 14.0 / NumberDensityUnits, smallx);
        BaryonField[N2INum][igrid] = max(y[IDX_N2I] * 28.0 / NumberDensityUnits, smallx);
        BaryonField[N2IINum][igrid] = max(y[IDX_N2II] * 28.0 / NumberDensityUnits, smallx);
        BaryonField[N2HIINum][igrid] = max(y[IDX_N2HII] * 29.0 / NumberDensityUnits, smallx);
        BaryonField[NHINum][igrid] = max(y[IDX_NHI] * 15.0 / NumberDensityUnits, smallx);
        BaryonField[NHIINum][igrid] = max(y[IDX_NHII] * 15.0 / NumberDensityUnits, smallx);
        BaryonField[NH2INum][igrid] = max(y[IDX_NH2I] * 16.0 / NumberDensityUnits, smallx);
        BaryonField[NH2IINum][igrid] = max(y[IDX_NH2II] * 16.0 / NumberDensityUnits, smallx);
        BaryonField[NH3INum][igrid] = max(y[IDX_NH3I] * 17.0 / NumberDensityUnits, smallx);
        BaryonField[NH3IINum][igrid] = max(y[IDX_NH3II] * 17.0 / NumberDensityUnits, smallx);
        BaryonField[NOINum][igrid] = max(y[IDX_NOI] * 30.0 / NumberDensityUnits, smallx);
        BaryonField[NOIINum][igrid] = max(y[IDX_NOII] * 30.0 / NumberDensityUnits, smallx);
        BaryonField[NO2INum][igrid] = max(y[IDX_NO2I] * 46.0 / NumberDensityUnits, smallx);
        BaryonField[OINum][igrid] = max(y[IDX_OI] * 16.0 / NumberDensityUnits, smallx);
        BaryonField[OIINum][igrid] = max(y[IDX_OII] * 16.0 / NumberDensityUnits, smallx);
        BaryonField[O2INum][igrid] = max(y[IDX_O2I] * 32.0 / NumberDensityUnits, smallx);
        BaryonField[O2IINum][igrid] = max(y[IDX_O2II] * 32.0 / NumberDensityUnits, smallx);
        BaryonField[O2HINum][igrid] = max(y[IDX_O2HI] * 33.0 / NumberDensityUnits, smallx);
        BaryonField[O2HIINum][igrid] = max(y[IDX_O2HII] * 33.0 / NumberDensityUnits, smallx);
        BaryonField[OCNINum][igrid] = max(y[IDX_OCNI] * 42.0 / NumberDensityUnits, smallx);
        BaryonField[OHINum][igrid] = max(y[IDX_OHI] * 17.0 / NumberDensityUnits, smallx);
        BaryonField[OHIINum][igrid] = max(y[IDX_OHII] * 17.0 / NumberDensityUnits, smallx);
        BaryonField[SiINum][igrid] = max(y[IDX_SiI] * 28.0 / NumberDensityUnits, smallx);
        BaryonField[SiIINum][igrid] = max(y[IDX_SiII] * 28.0 / NumberDensityUnits, smallx);
        BaryonField[SiCINum][igrid] = max(y[IDX_SiCI] * 40.0 / NumberDensityUnits, smallx);
        BaryonField[SiCIINum][igrid] = max(y[IDX_SiCII] * 40.0 / NumberDensityUnits, smallx);
        BaryonField[SiC2INum][igrid] = max(y[IDX_SiC2I] * 52.0 / NumberDensityUnits, smallx);
        BaryonField[SiC2IINum][igrid] = max(y[IDX_SiC2II] * 52.0 / NumberDensityUnits, smallx);
        BaryonField[SiC3INum][igrid] = max(y[IDX_SiC3I] * 64.0 / NumberDensityUnits, smallx);
        BaryonField[SiC3IINum][igrid] = max(y[IDX_SiC3II] * 64.0 / NumberDensityUnits, smallx);
        BaryonField[SiHINum][igrid] = max(y[IDX_SiHI] * 29.0 / NumberDensityUnits, smallx);
        BaryonField[SiHIINum][igrid] = max(y[IDX_SiHII] * 29.0 / NumberDensityUnits, smallx);
        BaryonField[SiH2INum][igrid] = max(y[IDX_SiH2I] * 30.0 / NumberDensityUnits, smallx);
        BaryonField[SiH2IINum][igrid] = max(y[IDX_SiH2II] * 30.0 / NumberDensityUnits, smallx);
        BaryonField[SiH3INum][igrid] = max(y[IDX_SiH3I] * 31.0 / NumberDensityUnits, smallx);
        BaryonField[SiH3IINum][igrid] = max(y[IDX_SiH3II] * 31.0 / NumberDensityUnits, smallx);
        BaryonField[SiH4INum][igrid] = max(y[IDX_SiH4I] * 32.0 / NumberDensityUnits, smallx);
        BaryonField[SiH4IINum][igrid] = max(y[IDX_SiH4II] * 32.0 / NumberDensityUnits, smallx);
        BaryonField[SiH5IINum][igrid] = max(y[IDX_SiH5II] * 33.0 / NumberDensityUnits, smallx);
        BaryonField[SiOINum][igrid] = max(y[IDX_SiOI] * 44.0 / NumberDensityUnits, smallx);
        BaryonField[SiOIINum][igrid] = max(y[IDX_SiOII] * 44.0 / NumberDensityUnits, smallx);
        BaryonField[SiOHIINum][igrid] = max(y[IDX_SiOHII] * 45.0 / NumberDensityUnits, smallx);
        
      }
    }
  }
  
  
  naunet.Finalize();

  delete [] temperature;
  delete [] sputtering;

  delete [] g_grid_dimension;
  delete [] g_grid_start;
  delete [] g_grid_end;

  LCAPERF_STOP("grid_NaunetWrapper");
#endif

  return SUCCESS;
}
