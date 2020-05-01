/***********************************************************************
/
/  ResistivityTableLookup()
/
/  written by: Duncan Christie
/  date:       April, 2016
/  PURPOSE:  Looks up the resistivity from the table read in
/            
/  INPUT:
/    D - density
/    T - temperature
/    B - magnetic field
/  RETURNS:
/    
/    SUCCESS or FAIL
/
************************************************************************/

#include <math.h>
#include "hdf5.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#define SMALL_LOG_VALUE -99.0

/**************************** Functions Prototypes ******************************/

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float Time);

// Initialize Cloudy Cooling
float ResistivityTableLookup(float D,float T, float B)
{
  float Time;
  int p,q;
  int Didx,Bidx,Tidx;
  int index;

  float slope;
  float value[2],value2[2],value3;
  float eta_perp;

  float V_D, V_B,V_T;

  if (ADUseLogInterp) {
    V_D = log10(D);
    V_B = log10(B);
    V_T = log10(T);

    Tidx = max(0,min(int((V_T-ResistivityData.TempParam[0])/ResistivityData.dlogT),ResistivityData.nT-2));
    Bidx = max(0,min(int((V_B-ResistivityData.MagParam[0])/ResistivityData.dlogB),ResistivityData.nB-2));
    Didx = max(0,min(int((V_D-ResistivityData.DensParam[0])/ResistivityData.dlogD),ResistivityData.nD-2));

  } else {
    V_D = D;
    V_B = B;
    V_T = T;

    Tidx = max(0,min(int(log10(T/ResistivityData.TempParam[0])/ResistivityData.dlogT),ResistivityData.nT-2));
    Bidx = max(0,min(int(log10(B/ResistivityData.MagParam[0])/ResistivityData.dlogB),ResistivityData.nB-2));
    Didx = max(0,min(int(log10(D/ResistivityData.DensParam[0])/ResistivityData.dlogD),ResistivityData.nD-2));

  }

  // Set the value so we're returning *something*
  eta_perp = 0.;


  

  for (p = 0;p <= 1;p++) {
    for (q = 0;q<=1;q++) {
      index = ResistivityData.nB*ResistivityData.nD*(Tidx+p) + ResistivityData.nB*(Didx+q) + Bidx;

      slope = (ResistivityData.ResistivityPerp[index+1]-ResistivityData.ResistivityPerp[index])/(ResistivityData.MagParam[Bidx+1]-ResistivityData.MagParam[Bidx]);

      value[q] = ResistivityData.ResistivityPerp[index] + slope*(V_B-ResistivityData.MagParam[Bidx]);
    }

    slope = (value[1]-value[0])/(ResistivityData.DensParam[Didx+1]-ResistivityData.DensParam[Didx]);

    value2[p] = value[0] + slope*(V_D-ResistivityData.DensParam[Didx]);
  }

  slope = (value2[1]-value2[0])/(ResistivityData.TempParam[Tidx+1]-ResistivityData.TempParam[Tidx]);

  value3 = value2[0] + slope*(V_T-ResistivityData.TempParam[Tidx]);

  eta_perp = pow(10.0,value3);

  
  return eta_perp;
}
