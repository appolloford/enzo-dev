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
  int GCH3OHINum, GCH4INum, GCOINum, GCO2INum, GH2CNINum, GH2COINum, GH2OINum,
      GH2SiOINum, GHCNINum, GHNCINum, GHNCOINum, GHNOINum, GMgINum, GN2INum,
      GNH3INum, GNOINum, GNO2INum, GO2INum, GO2HINum, GSiCINum, GSiC2INum,
      GSiC3INum, GSiH4INum, GSiOINum, CINum, CIINum, CHINum, CHIINum, CH2INum,
      CH2IINum, CH3INum, CH3IINum, CH3OHINum, CH4INum, CH4IINum, CNINum,
      CNIINum, COINum, COIINum, CO2INum, H2CNINum, H2COINum, H2COIINum,
      H2NOIINum, H2OINum, H2OIINum, H2SiOINum, H3IINum, H3COIINum, H3OIINum,
      HCNINum, HCNIINum, HCNHIINum, HCOINum, HCOIINum, HCO2IINum, HeHIINum,
      HNCINum, HNCOINum, HNOINum, HNOIINum, HOCIINum, MgINum, MgIINum, NINum,
      NIINum, N2INum, N2IINum, N2HIINum, NHINum, NHIINum, NH2INum, NH2IINum,
      NH3INum, NH3IINum, NOINum, NOIINum, NO2INum, OINum, OIINum, O2INum,
      O2IINum, O2HINum, O2HIINum, OCNINum, OHINum, OHIINum, SiINum, SiIINum,
      SiCINum, SiCIINum, SiC2INum, SiC2IINum, SiC3INum, SiC3IINum, SiHINum,
      SiHIINum, SiH2INum, SiH2IINum, SiH3INum, SiH3IINum, SiH4INum, SiH4IINum,
      SiH5IINum, SiOINum, SiOIINum, SiOHIINum;

  if (MultiSpecies == NAUNET_SPECIES) {
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
      BaryonField[DeNum][i] += 1.0 * BaryonField[CIINum][i] / 12.0 + 1.0 *
                            BaryonField[CHIINum][i] / 13.0 + 1.0 *
                            BaryonField[CH2IINum][i] / 14.0 + 1.0 *
                            BaryonField[CH3IINum][i] / 15.0 + 1.0 *
                            BaryonField[CH4IINum][i] / 16.0 + 1.0 *
                            BaryonField[CNIINum][i] / 26.0 + 1.0 *
                            BaryonField[COIINum][i] / 28.0 + 1.0 *
                            BaryonField[H2COIINum][i] / 30.0 + 1.0 *
                            BaryonField[H2NOIINum][i] / 32.0 + 1.0 *
                            BaryonField[H2OIINum][i] / 18.0 + 1.0 *
                            BaryonField[H3IINum][i] / 3.0 + 1.0 *
                            BaryonField[H3COIINum][i] / 31.0 + 1.0 *
                            BaryonField[H3OIINum][i] / 19.0 + 1.0 *
                            BaryonField[HCNIINum][i] / 27.0 + 1.0 *
                            BaryonField[HCNHIINum][i] / 28.0 + 1.0 *
                            BaryonField[HCOIINum][i] / 29.0 + 1.0 *
                            BaryonField[HCO2IINum][i] / 45.0 + 1.0 *
                            BaryonField[HeHIINum][i] / 5.0 + 1.0 *
                            BaryonField[HNOIINum][i] / 31.0 + 1.0 *
                            BaryonField[HOCIINum][i] / 29.0 + 1.0 *
                            BaryonField[MgIINum][i] / 24.0 + 1.0 *
                            BaryonField[NIINum][i] / 14.0 + 1.0 *
                            BaryonField[N2IINum][i] / 28.0 + 1.0 *
                            BaryonField[N2HIINum][i] / 29.0 + 1.0 *
                            BaryonField[NHIINum][i] / 15.0 + 1.0 *
                            BaryonField[NH2IINum][i] / 16.0 + 1.0 *
                            BaryonField[NH3IINum][i] / 17.0 + 1.0 *
                            BaryonField[NOIINum][i] / 30.0 + 1.0 *
                            BaryonField[OIINum][i] / 16.0 + 1.0 *
                            BaryonField[O2IINum][i] / 32.0 + 1.0 *
                            BaryonField[O2HIINum][i] / 33.0 + 1.0 *
                            BaryonField[OHIINum][i] / 17.0 + 1.0 *
                            BaryonField[SiIINum][i] / 28.0 + 1.0 *
                            BaryonField[SiCIINum][i] / 40.0 + 1.0 *
                            BaryonField[SiC2IINum][i] / 52.0 + 1.0 *
                            BaryonField[SiC3IINum][i] / 64.0 + 1.0 *
                            BaryonField[SiHIINum][i] / 29.0 + 1.0 *
                            BaryonField[SiH2IINum][i] / 30.0 + 1.0 *
                            BaryonField[SiH3IINum][i] / 31.0 + 1.0 *
                            BaryonField[SiH4IINum][i] / 32.0 + 1.0 *
                            BaryonField[SiH5IINum][i] / 33.0 + 1.0 *
                            BaryonField[SiOIINum][i] / 44.0 + 1.0 *
                            BaryonField[SiOHIINum][i] / 45.0;
    }
  }
#endif

  return SUCCESS;

}
