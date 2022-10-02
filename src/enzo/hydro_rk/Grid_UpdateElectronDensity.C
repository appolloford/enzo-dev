/***********************************************************************
/
/  GRID CLASS (CALCULATE ELECTRON DENSITY FROM SPECIES)
/
/  written by: John Wise
/  date:       July, 2009
/  modified1:
/
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#ifdef USE_NAUNET
#include "naunet_enzo.h"
#endif
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "TopGridData.h"
#include "Grid.h"

int FindField(int field, int farray[], int numfields);

int grid::UpdateElectronDensity(void)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  if (MultiSpecies == 0)
    return SUCCESS;

  int i, n, dim, size, nfield, n0;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;

  this->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, 
			      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);

#ifdef USE_NAUNET
  int GH2CNINum, GHNCINum, GNO2INum, GSiOINum, GCOINum, GHNCOINum, GMgINum,
      GNOINum, GO2INum, GO2HINum, GSiCINum, GSiC2INum, GSiC3INum, GCH3OHINum,
      GCO2INum, GH2SiOINum, GHNOINum, GN2INum, GH2COINum, GHCNINum, GH2OINum,
      GNH3INum, SiC3IINum, H2CNINum, GCH4INum, H2NOIINum, H2SiOINum, HeHIINum,
      HNCOINum, HOCIINum, SiC2IINum, GSiH4INum, SiC2INum, SiC3INum, SiH5IINum,
      SiH4IINum, SiCIINum, O2HINum, SiCINum, NO2INum, SiH3IINum, SiH2IINum,
      OCNINum, SiH2INum, SiOHIINum, SiHIINum, SiH4INum, SiHINum, SiH3INum,
      SiOIINum, HCO2IINum, HNOINum, CH3OHINum, MgINum, MgIINum, CH4IINum,
      SiOINum, CNIINum, HCNHIINum, N2HIINum, O2HIINum, SiIINum, SiINum, HNCINum,
      HNOIINum, N2IINum, H3COIINum, CH4INum, COIINum, NH3INum, CH3INum, CO2INum,
      NIINum, OIINum, HCNIINum, NH2IINum, NHIINum, O2IINum, CH3IINum, NH2INum,
      CH2IINum, H2OIINum, NH3IINum, NOIINum, H3OIINum, N2INum, CIINum, HCNINum,
      CHIINum, CH2INum, H2COIINum, NHINum, OHIINum, CNINum, H2COINum, HCOINum,
      CHINum, H3IINum, NOINum, NINum, OHINum, O2INum, CINum, HCOIINum, H2OINum,
      OINum, COINum;

  if (MultiSpecies == NAUNET_SPECIES) {
    if (IdentifyNaunetSpeciesFields(GH2CNINum, GHNCINum, GNO2INum, GSiOINum,
                                    GCOINum, GHNCOINum, GMgINum, GNOINum,
                                    GO2INum, GO2HINum, GSiCINum, GSiC2INum,
                                    GSiC3INum, GCH3OHINum, GCO2INum, GH2SiOINum,
                                    GHNOINum, GN2INum, GH2COINum, GHCNINum,
                                    GH2OINum, GNH3INum, SiC3IINum, H2CNINum,
                                    GCH4INum, H2NOIINum, H2SiOINum, HeHIINum,
                                    HNCOINum, HOCIINum, SiC2IINum, GSiH4INum,
                                    SiC2INum, SiC3INum, SiH5IINum, SiH4IINum,
                                    SiCIINum, O2HINum, SiCINum, NO2INum,
                                    SiH3IINum, SiH2IINum, OCNINum, SiH2INum,
                                    SiOHIINum, SiHIINum, SiH4INum, SiHINum,
                                    SiH3INum, SiOIINum, HCO2IINum, HNOINum,
                                    CH3OHINum, MgINum, MgIINum, CH4IINum,
                                    SiOINum, CNIINum, HCNHIINum, N2HIINum,
                                    O2HIINum, SiIINum, SiINum, HNCINum,
                                    HNOIINum, N2IINum, H3COIINum, CH4INum,
                                    COIINum, H2IINum, NH3INum, CH3INum, CO2INum,
                                    NIINum, OIINum, HCNIINum, NH2IINum, NHIINum,
                                    O2IINum, CH3IINum, NH2INum, CH2IINum,
                                    H2OIINum, NH3IINum, NOIINum, H3OIINum,
                                    N2INum, CIINum, HCNINum, CHIINum, CH2INum,
                                    H2COIINum, NHINum, OHIINum, CNINum,
                                    H2COINum, HCOINum, HeIINum, CHINum, H3IINum,
                                    HeINum, NOINum, NINum, OHINum, O2INum,
                                    CINum, HIINum, HCOIINum, H2OINum, OINum,
                                    DeNum, COINum, H2INum, HINum) == FAIL) {
      ENZO_FAIL("Error in grid->IdentifyNaunetSpeciesFields.");
    }
  }
#endif

  for (dim = 0, size = 1; dim < GridRank; dim++)
    size *= GridDimension[dim];
  
  for (i = 0; i < size; i++)
    BaryonField[DeNum][i] = BaryonField[HIINum][i] + 0.25*BaryonField[HeIINum][i] +
      0.5*BaryonField[HeIIINum][i]; 

  if (MultiSpecies > 1)
    for (i = 0; i < size; i++)
      BaryonField[DeNum][i] += 0.5*BaryonField[H2IINum][i] - BaryonField[HMNum][i];

  if (MultiSpecies > 2)
    for (i = 0; i < size; i++)
      BaryonField[DeNum][i] += 0.5*BaryonField[DIINum][i];

#ifdef USE_NAUNET
  if (MultiSpecies == NAUNET_SPECIES) {
    for (i = 0; i < size; i++) {
      BaryonField[DeNum][i] += 1.0 * BaryonField[SiC3IINum][i] / 64.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[H2NOIINum][i] / 32.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[HeHIINum][i] / 5.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[HOCIINum][i] / 29.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[SiC2IINum][i] / 52.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[SiH5IINum][i] / 33.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[SiH4IINum][i] / 32.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[SiCIINum][i] / 40.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[SiH3IINum][i] / 31.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[SiH2IINum][i] / 30.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[SiOHIINum][i] / 45.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[SiHIINum][i] / 29.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[SiOIINum][i] / 44.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[HCO2IINum][i] / 45.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[MgIINum][i] / 24.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[CH4IINum][i] / 16.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[CNIINum][i] / 26.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[HCNHIINum][i] / 28.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[N2HIINum][i] / 29.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[O2HIINum][i] / 33.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[SiIINum][i] / 28.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[HNOIINum][i] / 31.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[N2IINum][i] / 28.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[H3COIINum][i] / 31.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[COIINum][i] / 28.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[NIINum][i] / 14.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[OIINum][i] / 16.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[HCNIINum][i] / 27.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[NH2IINum][i] / 16.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[NHIINum][i] / 15.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[O2IINum][i] / 32.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[CH3IINum][i] / 15.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[CH2IINum][i] / 14.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[H2OIINum][i] / 18.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[NH3IINum][i] / 17.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[NOIINum][i] / 30.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[H3OIINum][i] / 19.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[CIINum][i] / 12.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[CHIINum][i] / 13.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[H2COIINum][i] / 30.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[OHIINum][i] / 17.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[H3IINum][i] / 3.0;
      BaryonField[DeNum][i] += 1.0 * BaryonField[HCOIINum][i] / 29.0;
      
    }
  }
#endif

  return SUCCESS;

}