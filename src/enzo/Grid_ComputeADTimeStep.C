/***********************************************************************
/
/  GRID CLASS (Compute thermal conduction time-scale)
/
/  written by:  Duncan Christie (based on Grid_ComputeConductionTimeStep.C)
/  date:        June 2016
/
/  PURPOSE:  Calculates the shortest time scale for ambipolar diffusion on a
/  grid
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <math.h> 
#include <stdio.h>
#include <stdlib.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "phys_constants.h"
#include "CosmologyParameters.h"


// Function prototypes
int GetUnits (float *DensityUnits, float *LengthUnits,
	      float *TemperatureUnits, float *TimeUnits,
	      float *VelocityUnits, double *MassUnits, FLOAT Time);
float GetResistivity(float rho, float B, float T, float Time);

// Member functions
int grid::ComputeADTimeStep (float &dt) {

  dt = huge_number;

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
  if (NumberOfBaryonFields == 0)
    return SUCCESS;
  this->DebugCheck("ComputeADTimeStep");

  // Some locals
  float rho,Bx,By,Bz,B2,eta,Bmag;
  int i,j,k;

  if (!UseAmbipolarDiffusion) {
    return SUCCESS;  
  }

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num,
    B1Num, B2Num, B3Num, PhiNum;
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num,
				   TENum, B1Num, B2Num, B3Num, PhiNum);

  /* Compute the field size. */
  int size = 1;
  int dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];


  float *Temp = new float[size];
  float T;
  /* The ambipolar diffusion module requires the temperature*/
  if (this->ComputeTemperatureField(Temp) == FAIL) {
    ENZO_FAIL("ComputeADTimeStep(): Error when calling ComputeTemperatureField()\n");
  }

  /*Grid sizes*/
  FLOAT dxinv = 1.0 / CellWidth[0][0];
  FLOAT dyinv = (GridRank > 1) ? 1.0 / CellWidth[1][0] : 0.0;
  FLOAT dzinv = (GridRank > 2) ? 1.0 / CellWidth[2][0] : 0.0;

  dxinv = max(dxinv,max(dyinv,dzinv));

  float dtinv_ad = 0.;
  int n = 0;
  for (k = 0; k < GridDimension[2]; k++) {
    for (j = 0; j < GridDimension[1]; j++) {
      for (i = 0; i < GridDimension[0]; i++, n++) {
	rho = BaryonField[DensNum][n];
	Bx  = BaryonField[B1Num][n];
	By  = BaryonField[B2Num][n];
	Bz  = BaryonField[B3Num][n];

	B2 = Bx*Bx + By*By + Bz*Bz;

	Bmag = sqrt(B2);
	T = Temp[n];

	eta = GetResistivity(rho,Bmag,T,Time);
	dtinv_ad = max(dtinv_ad,10.e0*eta*dxinv*dxinv);

      }
    }
  }

  delete [] Temp;

  if (dtinv_ad != 0) 
    dt = 1.0/dtinv_ad;


  if (dt == 0.) 
    fprintf(stderr,"AD Error: %f\n",dtinv_ad);

  return SUCCESS;
 
}
