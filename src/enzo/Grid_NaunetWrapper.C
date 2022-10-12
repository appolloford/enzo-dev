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

double GetAv(double nh) {

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

  // Eint32 *g_grid_dimension, *g_grid_start, *g_grid_end;
  // g_grid_dimension = new Eint32[GridRank];
  // g_grid_start = new Eint32[GridRank];
  // g_grid_end = new Eint32[GridRank];
  // for (i = 0; i < GridRank; i++) {
  //   g_grid_dimension[i] = (Eint32) GridDimension[i];
  //   g_grid_start[i] = (Eint32) GridStartIndex[i];
  //   g_grid_end[i] = (Eint32) GridEndIndex[i];
  // }
 
  // Find fields: density, total energy, velocity1-3.
 
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  // Find Multi-species fields.



  int specnum[NSPECIES] = {0};
  if (MultiSpecies == NAUNET_SPECIES) {
    IdentifyNaunetSpeciesFields(specnum);
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

  realtype y[NEQUATIONS] = {0.0};
  realtype y_init[NEQUATIONS];

  y[IDX_H2I] = 0.5;
  y[IDX_HI]  = 5.0e-5;
  y[IDX_HeI] = 9.75e-2;
  y[IDX_NI]  = 7.5e-5;
  y[IDX_OI]  = 1.8e-4;
  y[IDX_COI] = 1.4e-4;
  y[IDX_MgI] = 7.0e-9;
  y[IDX_SiI] = 8.0e-9;

  naunet.SetReferenceAbund(y, 1);

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
        data.zeta     = var_zeta;
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

        for (int sidx = 0; sidx < NSPECIES; sidx++) {
          int snum = specnum[sidx];
          y[sidx] = BaryonField[snum][igrid] * NumberDensityUnits / A_Table[sidx];
          y_init[sidx] = y[sidx];
        }

    
        int flag = naunet.Solve(y, dt_chem * TimeUnits, &data);
    
        if (flag == NAUNET_FAIL) {

          // Try varying atol to fix the problem
          for (int ntrial = 1; ntrial < 6; ntrial ++) {
            naunet.Reset(1, pow(10.0, ntrial) * atol, 1e-5, 1000);

            for (int idx = IDX_GH2CNI; idx <= IDX_HI; idx++) {
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

        // Renormalize to ensure element abundance conserved
        // for (int sidx = 0; sidx < NSPECIES; sidx++) {
        //   y[sidx] = max(y[sidx], 1.e-30 * density[igrid] * NumberDensityUnits / A_Table[sidx]);
        // }

        // if (naunet.Renorm(y) == NAUNET_FAIL) {
        //   ENZO_FAIL("Naunet renorm failed in Grid_NaunetWrapper!");
        // }

        for (int sidx = 0; sidx < NSPECIES; sidx++) {
          int snum = specnum[sidx];
          BaryonField[snum][igrid] = y[sidx] * A_Table[sidx] / NumberDensityUnits;
        }
        
      }
    }
  }
  
  
  naunet.Finalize();

  delete [] temperature;
  delete [] sputtering;

  // delete [] g_grid_dimension;
  // delete [] g_grid_start;
  // delete [] g_grid_end;

  LCAPERF_STOP("grid_NaunetWrapper");
#endif

  return SUCCESS;
}
