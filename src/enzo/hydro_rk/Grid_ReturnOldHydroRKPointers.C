/***********************************************************************
/
/  GRID CLASS (RETURNS AN ARRAY OF POINTERS THAT ARE COMPATIBLE WITH
/              THE HYDRO_RK SOLVERS)
/
/  written by: John Wise
/  date:       July, 2009
/  modified1:
/
/
************************************************************************/
// clang-format off
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

int grid::ReturnOldHydroRKPointers(float **Prim, bool ReturnMassFractions)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  int i, n, dim, size, nfield, n0;
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int B1Num, B2Num, B3Num, PhiNum;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
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
#endif

  /* Add the physical quantities */

  if (HydroMethod == HD_RK) {
    this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				     Vel3Num, TENum);
    nfield = n0 = NEQ_HYDRO;
  }

  else if (HydroMethod == MHD_RK) {
    this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				     Vel3Num, TENum, B1Num, B2Num, B3Num, 
				     PhiNum);
    nfield = n0 = NEQ_MHD;
  }
  
  Prim[iden] = OldBaryonField[DensNum];
  Prim[ivx] = OldBaryonField[Vel1Num];
  if (GridRank > 1)
    Prim[ivy] = OldBaryonField[Vel2Num];
  if (GridRank > 2)
    Prim[ivz] = OldBaryonField[Vel3Num];
  Prim[ietot] = OldBaryonField[TENum];
  if (DualEnergyFormalism)
    Prim[ieint] = OldBaryonField[GENum];

  if (HydroMethod == MHD_RK) {
    Prim[iBx] = OldBaryonField[B1Num];
    Prim[iBy] = OldBaryonField[B2Num];
    Prim[iBz] = OldBaryonField[B3Num];
    Prim[iPhi] = OldBaryonField[PhiNum];
  }

  /* Add the species */

  if (MultiSpecies) {
    this->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, 
                                HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);
#ifdef USE_NAUNET
  if (MultiSpecies == NAUNET_SPECIES)
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
#endif

    //Prim[nfield++] = OldBaryonField[DeNum];
    Prim[nfield++] = OldBaryonField[HINum];
    Prim[nfield++] = OldBaryonField[HIINum];
    Prim[nfield++] = OldBaryonField[HeINum];
    Prim[nfield++] = OldBaryonField[HeIINum];
    Prim[nfield++] = OldBaryonField[HeIIINum];

    if (MultiSpecies > 1) {
      Prim[nfield++] = OldBaryonField[HMNum];
      Prim[nfield++] = OldBaryonField[H2INum];
      Prim[nfield++] = OldBaryonField[H2IINum];
    }

    if (MultiSpecies > 2) {
      Prim[nfield++] = OldBaryonField[DINum];
      Prim[nfield++] = OldBaryonField[DIINum];
      Prim[nfield++] = OldBaryonField[HDINum];
    }

#ifdef USE_NAUNET
    if (MultiSpecies == NAUNET_SPECIES) {
      Prim[nfield++] = OldBaryonField[GCH3OHINum];
      Prim[nfield++] = OldBaryonField[GCH4INum];
      Prim[nfield++] = OldBaryonField[GCOINum];
      Prim[nfield++] = OldBaryonField[GCO2INum];
      Prim[nfield++] = OldBaryonField[GH2CNINum];
      Prim[nfield++] = OldBaryonField[GH2COINum];
      Prim[nfield++] = OldBaryonField[GH2OINum];
      Prim[nfield++] = OldBaryonField[GH2SiOINum];
      Prim[nfield++] = OldBaryonField[GHCNINum];
      Prim[nfield++] = OldBaryonField[GHNCINum];
      Prim[nfield++] = OldBaryonField[GHNCOINum];
      Prim[nfield++] = OldBaryonField[GHNOINum];
      Prim[nfield++] = OldBaryonField[GMgINum];
      Prim[nfield++] = OldBaryonField[GN2INum];
      Prim[nfield++] = OldBaryonField[GNH3INum];
      Prim[nfield++] = OldBaryonField[GNOINum];
      Prim[nfield++] = OldBaryonField[GNO2INum];
      Prim[nfield++] = OldBaryonField[GO2INum];
      Prim[nfield++] = OldBaryonField[GO2HINum];
      Prim[nfield++] = OldBaryonField[GSiCINum];
      Prim[nfield++] = OldBaryonField[GSiC2INum];
      Prim[nfield++] = OldBaryonField[GSiC3INum];
      Prim[nfield++] = OldBaryonField[GSiH4INum];
      Prim[nfield++] = OldBaryonField[GSiOINum];
      Prim[nfield++] = OldBaryonField[CINum];
      Prim[nfield++] = OldBaryonField[CIINum];
      Prim[nfield++] = OldBaryonField[CHINum];
      Prim[nfield++] = OldBaryonField[CHIINum];
      Prim[nfield++] = OldBaryonField[CH2INum];
      Prim[nfield++] = OldBaryonField[CH2IINum];
      Prim[nfield++] = OldBaryonField[CH3INum];
      Prim[nfield++] = OldBaryonField[CH3IINum];
      Prim[nfield++] = OldBaryonField[CH3OHINum];
      Prim[nfield++] = OldBaryonField[CH4INum];
      Prim[nfield++] = OldBaryonField[CH4IINum];
      Prim[nfield++] = OldBaryonField[CNINum];
      Prim[nfield++] = OldBaryonField[CNIINum];
      Prim[nfield++] = OldBaryonField[COINum];
      Prim[nfield++] = OldBaryonField[COIINum];
      Prim[nfield++] = OldBaryonField[CO2INum];
      Prim[nfield++] = OldBaryonField[H2CNINum];
      Prim[nfield++] = OldBaryonField[H2COINum];
      Prim[nfield++] = OldBaryonField[H2COIINum];
      Prim[nfield++] = OldBaryonField[H2NOIINum];
      Prim[nfield++] = OldBaryonField[H2OINum];
      Prim[nfield++] = OldBaryonField[H2OIINum];
      Prim[nfield++] = OldBaryonField[H2SiOINum];
      Prim[nfield++] = OldBaryonField[H3IINum];
      Prim[nfield++] = OldBaryonField[H3COIINum];
      Prim[nfield++] = OldBaryonField[H3OIINum];
      Prim[nfield++] = OldBaryonField[HCNINum];
      Prim[nfield++] = OldBaryonField[HCNIINum];
      Prim[nfield++] = OldBaryonField[HCNHIINum];
      Prim[nfield++] = OldBaryonField[HCOINum];
      Prim[nfield++] = OldBaryonField[HCOIINum];
      Prim[nfield++] = OldBaryonField[HCO2IINum];
      Prim[nfield++] = OldBaryonField[HeHIINum];
      Prim[nfield++] = OldBaryonField[HNCINum];
      Prim[nfield++] = OldBaryonField[HNCOINum];
      Prim[nfield++] = OldBaryonField[HNOINum];
      Prim[nfield++] = OldBaryonField[HNOIINum];
      Prim[nfield++] = OldBaryonField[HOCIINum];
      Prim[nfield++] = OldBaryonField[MgINum];
      Prim[nfield++] = OldBaryonField[MgIINum];
      Prim[nfield++] = OldBaryonField[NINum];
      Prim[nfield++] = OldBaryonField[NIINum];
      Prim[nfield++] = OldBaryonField[N2INum];
      Prim[nfield++] = OldBaryonField[N2IINum];
      Prim[nfield++] = OldBaryonField[N2HIINum];
      Prim[nfield++] = OldBaryonField[NHINum];
      Prim[nfield++] = OldBaryonField[NHIINum];
      Prim[nfield++] = OldBaryonField[NH2INum];
      Prim[nfield++] = OldBaryonField[NH2IINum];
      Prim[nfield++] = OldBaryonField[NH3INum];
      Prim[nfield++] = OldBaryonField[NH3IINum];
      Prim[nfield++] = OldBaryonField[NOINum];
      Prim[nfield++] = OldBaryonField[NOIINum];
      Prim[nfield++] = OldBaryonField[NO2INum];
      Prim[nfield++] = OldBaryonField[OINum];
      Prim[nfield++] = OldBaryonField[OIINum];
      Prim[nfield++] = OldBaryonField[O2INum];
      Prim[nfield++] = OldBaryonField[O2IINum];
      Prim[nfield++] = OldBaryonField[O2HINum];
      Prim[nfield++] = OldBaryonField[O2HIINum];
      Prim[nfield++] = OldBaryonField[OCNINum];
      Prim[nfield++] = OldBaryonField[OHINum];
      Prim[nfield++] = OldBaryonField[OHIINum];
      Prim[nfield++] = OldBaryonField[SiINum];
      Prim[nfield++] = OldBaryonField[SiIINum];
      Prim[nfield++] = OldBaryonField[SiCINum];
      Prim[nfield++] = OldBaryonField[SiCIINum];
      Prim[nfield++] = OldBaryonField[SiC2INum];
      Prim[nfield++] = OldBaryonField[SiC2IINum];
      Prim[nfield++] = OldBaryonField[SiC3INum];
      Prim[nfield++] = OldBaryonField[SiC3IINum];
      Prim[nfield++] = OldBaryonField[SiHINum];
      Prim[nfield++] = OldBaryonField[SiHIINum];
      Prim[nfield++] = OldBaryonField[SiH2INum];
      Prim[nfield++] = OldBaryonField[SiH2IINum];
      Prim[nfield++] = OldBaryonField[SiH3INum];
      Prim[nfield++] = OldBaryonField[SiH3IINum];
      Prim[nfield++] = OldBaryonField[SiH4INum];
      Prim[nfield++] = OldBaryonField[SiH4IINum];
      Prim[nfield++] = OldBaryonField[SiH5IINum];
      Prim[nfield++] = OldBaryonField[SiOINum];
      Prim[nfield++] = OldBaryonField[SiOIINum];
      Prim[nfield++] = OldBaryonField[SiOHIINum];
          }
#endif

  } // ENDIF MultiSpecies

  /* Add the colours (treat them as species) */

  int SNColourNum, MetalNum, MetalIaNum, MetalIINum, MBHColourNum, Galaxy1ColourNum, 
    Galaxy2ColourNum; 

  if (this->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum, MetalIINum, MBHColourNum, 
				 Galaxy1ColourNum, Galaxy2ColourNum) == FAIL) {
    fprintf(stderr, "Error in grid->IdentifyColourFields.\n");
    return FAIL;
  }
  
  if (MetalNum != -1) {
    Prim[nfield++] = OldBaryonField[MetalNum];
    if (StarMakerTypeIaSNe)
      Prim[nfield++] = OldBaryonField[MetalIaNum];
    if (StarMakerTypeIISNeMetalField)
      Prim[nfield++] = OldBaryonField[MetalIINum];
    if (MultiMetals || TestProblemData.MultiMetals) {
      Prim[nfield++] = OldBaryonField[MetalNum+1];
      Prim[nfield++] = OldBaryonField[MetalNum+2];
    }
  }

  if (SNColourNum      != -1) Prim[nfield++] = OldBaryonField[SNColourNum];  
  /*   //##### These fields are currently not being used and only causing interpolation problems
  if (MBHColourNum     != -1) Prim[nfield++] = OldBaryonField[MBHColourNum];
  if (Galaxy1ColourNum != -1) Prim[nfield++] = OldBaryonField[Galaxy1ColourNum];
  if (Galaxy2ColourNum != -1) Prim[nfield++] = OldBaryonField[Galaxy2ColourNum];
  */

  /* Convert the species and color fields into mass fractions */

  for (dim = 0, size = 1; dim < GridRank; dim++)
    size *= GridDimension[dim];

  if (ReturnMassFractions)  
    for (n = n0; n < nfield; n++)
      for (i = 0; i < size; i++)
	Prim[n][i] /= Prim[iden][i];  

  return SUCCESS;

}