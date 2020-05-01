/***********************************************************************
/
/  GRID CLASS (PREPARE NORMALIZATION FOR VELOCITY FIELDS)
/            
/
/  written by: Tom Abel
/  date:       September, 2009
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
 
int grid::NormalizeVelocities(Eflt fac)
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
      return SUCCESS;

  if (fac == 0.) 
    return SUCCESS; // wouldn't call this a normalization

  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  fprintf(stderr, "Grid_NormalizeVelocities: \n");
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "GPRFN: Error in IdentifyPhysicalQuantities.\n");
    ENZO_FAIL("");
  }

  int i, j, k, index, dim;
  float x, y, z;//
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {

        x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i]; //
        y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j]; //
        z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k]; //

        index = i + j*GridDimension[0] + k*GridDimension[0]*GridDimension[1];
        BaryonField[TENum][index] -= 0.5*(pow(BaryonField[Vel1Num][index], 2) +
                                          pow(BaryonField[Vel2Num][index], 2) +
                                          pow(BaryonField[Vel3Num][index], 2) );
        for (dim = 0; dim < GridRank; dim++) {
          int vel = Vel1Num + dim;
          BaryonField[vel][index] *= fac;
        }
        /* add bulk velocity */
        if (x < 0.5) {
          //BaryonField[Vel1Num][index] += 5.0; /* 20 km/s */ //to changevel
          BaryonField[Vel1Num][index] += 2.5; /* 10 km/s velocity: 1=2km/s. (prev:vel=2.0) */
          //BaryonField[Vel1Num][index] += 1.25; /* 5 km/s */
          //BaryonField[Vel1Num][index] += 0.0; 
        } else {
          //BaryonField[Vel1Num][index] += -5.0;
          BaryonField[Vel1Num][index] += -2.5;
          //BaryonField[Vel1Num][index] += -1.25;
          //BaryonField[Vel1Num][index] += 0.0;
        }

        BaryonField[TENum][index] += 0.5*(pow(BaryonField[Vel1Num][index], 2) +
                                          pow(BaryonField[Vel2Num][index], 2) +
                                          pow(BaryonField[Vel3Num][index], 2) );

      }

  return SUCCESS;
}

