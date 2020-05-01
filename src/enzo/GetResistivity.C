/***********************************************************************
/
/  GetResistivity()
/
/  written by: Duncan Christie
/  date:       April, 2016
/  PURPOSE: Returns the ambipolar diffusion resistivity coefficient.
/
/  NOTE:  The value returned includes all the factors of clight and 4*Pi
/  absorbed into the coefficient.  Some papers (see any papers from
/  the Mouschovias group, for example) keep those constants out of the
/  resistivity.  This makes more sense aesthetically, but it is preferential
/  to not have all those factors of clight floating around the code.
/ 
/            
/  INPUT:
/    rho - density
/    T - temperature
/    B - magnetic field
/    Time - time
/  RETURNS:
/    
/    The resistivity (in code units)
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
#include "phys_constants.h"

#define SMALL_LOG_VALUE -99.0

/**************************** Functions Prototypes ******************************/

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float Time);

float ResistivityTableLookup(float D,float T, float B);

float GetResistivity(float rho, float B, float T,float Time)
{
  float eta_perp = 0.;
  float rhoi = 0.;
  float D = 0.;
  float alpha;
  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1,
    TimeUnits = 1.0, VelocityUnits = 1.0;


  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
           &TimeUnits, &VelocityUnits, Time);
  float MagneticUnits = sqrt(DensityUnits*4.0*M_PI)*VelocityUnits;

  float clight_code = clight/VelocityUnits;
  

  rho *= DensityUnits;
  B *= MagneticUnits;

  //T *= TemperatureUnits;

  if (debug) {
    if (T != T) {
      fprintf(stderr, "GetResistivity(): T is NaN!\n");
      fprintf(stderr, "GetResistivity(): D = %e, T = %e, B= %e\n",D,T,B);
    }
  }

  switch (ADResistivityType) {
  case 1:
    return ADResistivityScale*TimeUnits/(LengthUnits*LengthUnits);
  case 2:
    eta_perp = 0.;
    eta_perp = ADResistivityScale*B*B/(4.*pi*pow(rho,1.5));
    eta_perp *= TimeUnits/(LengthUnits*LengthUnits);
    return eta_perp;
  case 3:
    eta_perp = clight_code*clight_code*ResistivityTableLookup(rho,T,B)/(4.*pi);
    eta_perp /= TimeUnits;
    return eta_perp;
  case 4:
    float ni,nn;
    nn = rho/(2.33*1.6733e-24);
    ni = 1e-3*pow(nn/1e5,0.5) + 1e-3*pow(nn/1e3,-2);
    rhoi = 29.0*1.6733e-24*ni;
    alpha = 3.7e13;
    eta_perp = B*B/(4.0*pi*alpha*rhoi*rho);
    eta_perp *= TimeUnits/(LengthUnits*LengthUnits);
    return eta_perp;
  default:
    return 0.;
  }

}

