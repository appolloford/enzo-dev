/***********************************************************************
/
/  GRID CLASS (WRAP THE NAUNET RENORMALIZAZTION )
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

int grid::NaunetRenormWrapper()
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

#ifdef USE_NAUNET

  if (!use_naunetrenorm) return SUCCESS;

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
        VelocityUnits = 1, TimeUnits = 1, aUnits = 1;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
          &TimeUnits, &VelocityUnits, Time);

  float NumberDensityUnits = DensityUnits / mh;

  int specnum[NSPECIES] = {0};
  if (MultiSpecies == NAUNET_SPECIES) {
    IdentifyNaunetSpeciesFields(specnum);
  }

  // Initial abundances
  float yab[NAUNET_NEQUATIONS] = {0.0};
  yab[IDX_H2I] = 0.5;
  yab[IDX_HI]  = 5.0e-5;
  yab[IDX_HeI] = 9.75e-2;
  yab[IDX_NI]  = 7.5e-5;
  yab[IDX_OI]  = 1.8e-4;
  yab[IDX_COI] = 1.4e-4;
  yab[IDX_MgI] = 7.0e-9;
  yab[IDX_SiI] = 8.0e-9;

  Naunet naunet;
  naunet.SetReferenceAbund(yab, 1);

  // for (int igrid = 0; igrid < size; igrid++) {
  for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      int igrid = (k * GridDimension[1] + j) * GridDimension[0] + GridStartIndex[0];
      for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, igrid++) {


        // printf("%d\n", igrid);
        for (int sidx = 0; sidx < NSPECIES; sidx++) {
          int snum = specnum[sidx];
          yab[sidx] = BaryonField[snum][igrid] * NumberDensityUnits / A_Table[sidx];
        }

        int flag = naunet.Renorm(yab);

        // Fail to Renormalize
        if (flag == NAUNET_FAIL) {
          ENZO_FAIL("Naunet renorm failed!");
        }

        for (int sidx = 0; sidx < NSPECIES; sidx++) {
          int snum = specnum[sidx];
          BaryonField[snum][igrid] = yab[sidx] * A_Table[sidx] / NumberDensityUnits;
        }
      }
    }
  }

  LCAPERF_STOP("grid_NaunetRenormWrapper");
#endif

  return SUCCESS;
}
