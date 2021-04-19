/***********************************************************************
/
/  GRID CLASS (SET THE BACKGROUND BULK VELOCITY FIELDS)
/            
/
/  written by: Chia-Jung Hsu
/  date:       March, 2021
/  modified1:
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

int GetUnits(float *DensityUnits, float *LengthUnits, float *TemperatureUnits, 
             float *TimeUnits, float *VelocityUnits, FLOAT Time);

int grid::SetBulkVelocities(float RelativeVelocity)
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
      return SUCCESS;

  if (RelativeVelocity == 0.) 
    return SUCCESS; // wouldn't call this a normalization

  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  fprintf(stderr, "Grid_NormalizeVelocities: \n");
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "GPRFN: Error in IdentifyPhysicalQuantities.\n");
    ENZO_FAIL("");
  }

  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1.0, TimeUnits = 1.0,
    VelocityUnits = 1.0;
  if (UsePhysicalUnit)
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, Time);

  int i, j, k, index, dim;
  float x, y, z;
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {

        x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
        y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
        z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

        index = i + j*GridDimension[0] + k*GridDimension[0]*GridDimension[1];
        BaryonField[TENum][index] -= 0.5*(pow(BaryonField[Vel1Num][index], 2) +
                                          pow(BaryonField[Vel2Num][index], 2) +
                                          pow(BaryonField[Vel3Num][index], 2) );
        /* add bulk velocity */
        if (x < 0.5) {
          BaryonField[Vel1Num][index] += 0.5 * RelativeVelocity;
        } else {
          BaryonField[Vel1Num][index] -= 0.5 * RelativeVelocity;
        }

        BaryonField[TENum][index] += 0.5*(pow(BaryonField[Vel1Num][index], 2) +
                                          pow(BaryonField[Vel2Num][index], 2) +
                                          pow(BaryonField[Vel3Num][index], 2) );

      }

  return SUCCESS;
}

