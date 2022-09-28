/***********************************************************************
/
/  GRID CLASS (WRAP THE NAUNET SOLVER FOR POSTPROCESSING)
/
/  written by: Chia-Jung Hsu
/  date:       2021
/  modified1:
/
/  PURPOSE: Renorm species with naunet.
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
// clang-format off
#include "preincludes.h"
#ifdef USE_NAUNET
#include "naunet_enzo.h"
#include "naunet_physics.h"
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

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);
double GetAv(double nh);

int grid::NaunetPostProcess(FLOAT TargetTime)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

#ifdef USE_NAUNET
  
  // Compute the size of the fields.
 
  int i, j, k, igrid;
  int size = 1;
  int activesize = 1;
  for (int dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
    activesize *= (GridDimension[dim] - 2*NumberOfGhostZones);
  }

  // Find fields: density, total energy, velocity1-3.

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
 
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  // Get easy to handle pointers for each variable.

  float *density     = BaryonField[DensNum];
  // float *totalenergy = BaryonField[TENum];
  // float *gasenergy   = BaryonField[GENum];
  // float *velocity1   = BaryonField[Vel1Num];
  // float *velocity2   = BaryonField[Vel2Num];
  // float *velocity3   = BaryonField[Vel3Num];

  // Compute the cooling time.

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
        VelocityUnits = 1, TimeUnits = 1, aUnits = 1;

  float NumberDensityUnits = DensityUnits / mh;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
           &TimeUnits, &VelocityUnits, Time);

  float *temperature = new float[size]; 
  if (this->ComputeTemperatureField(temperature) == FAIL){
    ENZO_FAIL("Error in grid->ComputeTemperatureField.");
  }

  int specnum[NSPECIES] = {0};
  if (MultiSpecies == NAUNET_SPECIES) {
    IdentifyNaunetSpeciesFields(specnum);
  }

  float atol = 1e-30, smallx = 1e-40;
  float dt_target = TargetTime * TimeUnits;
  float logtend = log10(dt_target);
  float logtfirst = logtend - 4.0;

  float *tout = new float[41];
  for (int step = 0; step < 41; step ++) {
    // tout[step] = pow(10, log10(dt_target) * (step+1) / 40.0);
    tout[step] = pow(10.0, logtfirst + step * 0.1);
  }

  Naunet naunet;
  NaunetData data;

  naunet.Init(1, atol, 1e-5, 1000);

  realtype y[NAUNET_NEQUATIONS], y_init[NAUNET_NEQUATIONS];

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      igrid = (k * GridDimension[1] + j) * GridDimension[0] + GridStartIndex[0];
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, igrid++) {

        double nH     = BaryonField[iden][igrid] * DensityUnits / (1.4 * mh);
        data.nH       = nH;
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

        // for (int sidx = 0; sidx < NSPECIES; sidx++) {
        //   int snum = specnum[sidx];
        //   y[sidx] = max(BaryonField[snum][igrid], smallx) * NumberDensityUnits / A_Table[sidx];
        //   y_init[sidx] = y[sidx];
        // }

        // for (int sidx = 0; sidx < NSPECIES; sidx ++) {
        //   y[sidx] = smallx * NumberDensityUnits / A_Table[sidx];
        // }

        // y[IDX_H2I] = density[igrid] * NumberDensityUnits * 1.00e+0 / 1.4;
        // y[IDX_HI]  = density[igrid] * NumberDensityUnits * 5.00e-5 / 1.4;
        // y[IDX_HeI] = density[igrid] * NumberDensityUnits * 3.90e-1 / 1.4;
        // y[IDX_NI ] = density[igrid] * NumberDensityUnits * 7.50e-5 * 14.0 / 1.4;
        // y[IDX_OI]  = density[igrid] * NumberDensityUnits * 1.80e-4 * 16.0 / 1.4;
        // y[IDX_COI] = density[igrid] * NumberDensityUnits * 1.40e-4 * 28.0 / 1.4;
        // y[IDX_MgI] = density[igrid] * NumberDensityUnits * 7.00e-9 * 24.0 / 1.4;
        // y[IDX_SiI] = density[igrid] * NumberDensityUnits * 8.00e-9 * 28.0 / 1.4;

        for (int sidx = 0; sidx < NSPECIES; sidx ++) {
          y[sidx] = 1e-30;
        }

        y[IDX_H2I] = nH * 1.00e+0;
        y[IDX_HI]  = nH * 5.00e-5;
        y[IDX_HeI] = nH * 9.75e-2;
        y[IDX_NI ] = nH * 7.50e-5;
        y[IDX_OI]  = nH * 1.80e-4;
        y[IDX_COI] = nH * 1.40e-4;
        y[IDX_MgI] = nH * 7.00e-9;
        y[IDX_SiI] = nH * 8.00e-9;

        for (int sidx = 0; sidx < NSPECIES; sidx++) {
          y_init[sidx] = y[sidx];
        }

        int flag = NAUNET_SUCCESS;
        float tcur = 0.0;

        for (int ntrial = 0; ntrial < 6; ntrial ++) {
          if (ntrial) {
            naunet.Reset(1, pow(10.0, ntrial) * atol, 1e-5, 1000);
            flag = NAUNET_SUCCESS;
            tcur = 0.0;
            for (int sidx = 0; sidx < NSPECIES; sidx++) {
              y[sidx] = y_init[sidx];
            }
          }

          for (int step = 0; step < 41; step ++) {
            float dt = tout[step] - tcur;
            tcur += dt;
            if (naunet.Solve(y, dt, &data) == NAUNET_FAIL) {
              flag = NAUNET_FAIL;
              break;
            }
          }

          if (flag == NAUNET_SUCCESS) {
            if (ntrial) {
              naunet.Reset(1, atol, 1e-5, 1000);
            }
            break;
          }
        }

        // Fail in all trials
        if (flag == NAUNET_FAIL) {
          naunet.Finalize();
          ENZO_FAIL("Naunet failed in NaunetPostProcess.C!");
        }

        // for (int step = 0; step < 41; step ++) {

        //   float dt = tout[step] - tcur;
        //   tcur += dt;

        //   int flag = naunet.Solve(y, dt, &data);

        //   if (flag == NAUNET_FAIL) {

        //     // Try varying atol to fix the problem
        //     for (int ntrial = 1; ntrial < 6; ntrial ++) {
        //       naunet.Reset(1, pow(10.0, ntrial) * atol, 1e-5, 1000);

        //       for (int idx = IDX_GCH3OHI; idx <= IDX_SiOHII; idx++) {
        //         y[idx] = y_init[idx];
        //       }

        //       flag = naunet.Solve(y, dt, &data);

        //       if (flag == NAUNET_SUCCESS) {
        //         naunet.Reset(1, atol, 1e-5, 1000);
        //         break;
        //       }
        //     }

        //     // Fail to fix up
        //     if (flag == NAUNET_FAIL) {
        //       naunet.Finalize();

        //       ENZO_FAIL("Naunet failed in NaunetPostProcess.C!");
        //     }
        //   }
        // }

        for (int sidx = 0; sidx < NSPECIES; sidx++) {
          int snum = specnum[sidx];
          BaryonField[snum][igrid] = max(y[sidx] * A_Table[sidx] / NumberDensityUnits, smallx);
        }
      }
    }
  }

  naunet.Finalize();

  delete [] temperature;
  delete [] tout;

#endif

  return SUCCESS;
}
