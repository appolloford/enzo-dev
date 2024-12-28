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

int grid::ReturnHydroRKPointers(float **Prim, bool ReturnMassFractions)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  int i, n, dim, size, nfield, n0;
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int B1Num, B2Num, B3Num, PhiNum, CRNum;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
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
#endif

  /* Add the physical quantities */

  if (HydroMethod == HD_RK) {
    this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				     Vel3Num, TENum);
    nfield = n0 = NEQ_HYDRO;
  }

  else if (HydroMethod == MHD_RK) {
    if (CRModel)
      this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum, B1Num, B2Num, B3Num, PhiNum, CRNum);
    else 
      this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				     Vel3Num, TENum, B1Num, B2Num, B3Num, 
				     PhiNum);
    nfield = n0 = NEQ_MHD;
  };

  Prim[iden] = BaryonField[DensNum];
  Prim[ivx]  = BaryonField[Vel1Num];
  Prim[ivy]  = BaryonField[Vel2Num];
  Prim[ivz]  = BaryonField[Vel3Num];
  Prim[ietot]= BaryonField[TENum];
  if (DualEnergyFormalism)
    Prim[ieint] = BaryonField[GENum];

  if (HydroMethod == MHD_RK) {
    Prim[iBx] = BaryonField[B1Num];
    Prim[iBy] = BaryonField[B2Num];
    Prim[iBz] = BaryonField[B3Num];
    Prim[iPhi]= BaryonField[PhiNum];
  }

  if (CRModel){
    Prim[iCR] = BaryonField[CRNum];
  }

  /*
  printf("Physical Quantities: %"ISYM" %"ISYM"  %"ISYM" %"ISYM" %"ISYM"  %"ISYM"  %"ISYM" %"ISYM" %"ISYM" %"ISYM"\n", 
	 DensNum, GENum, Vel1Num, Vel2Num, 
	 Vel3Num, TENum, B1Num, B2Num, B3Num, 
	 PhiNum);
  */
  /* Add the species */

  if (MultiSpecies) {
    this->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, 
                                HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);
#ifdef USE_NAUNET
  if (MultiSpecies == NAUNET_SPECIES)
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
#endif

    //Prim[nfield++] = BaryonField[DeNum];
    Prim[nfield++] = BaryonField[HINum];
    Prim[nfield++] = BaryonField[HIINum];
    Prim[nfield++] = BaryonField[HeINum];
    Prim[nfield++] = BaryonField[HeIINum];
    Prim[nfield++] = BaryonField[HeIIINum];

    if (MultiSpecies > 1) {
      Prim[nfield++] = BaryonField[HMNum];
      Prim[nfield++] = BaryonField[H2INum];
      Prim[nfield++] = BaryonField[H2IINum];
    }

    if (MultiSpecies > 2) {
      Prim[nfield++] = BaryonField[DINum];
      Prim[nfield++] = BaryonField[DIINum];
      Prim[nfield++] = BaryonField[HDINum];
    }

#ifdef USE_NAUNET
    if (MultiSpecies == NAUNET_SPECIES) {
      Prim[nfield++] = BaryonField[GH2CNINum];
      Prim[nfield++] = BaryonField[GHNCINum];
      Prim[nfield++] = BaryonField[GNO2INum];
      Prim[nfield++] = BaryonField[GSiOINum];
      Prim[nfield++] = BaryonField[GCOINum];
      Prim[nfield++] = BaryonField[GHNCOINum];
      Prim[nfield++] = BaryonField[GMgINum];
      Prim[nfield++] = BaryonField[GNOINum];
      Prim[nfield++] = BaryonField[GO2INum];
      Prim[nfield++] = BaryonField[GO2HINum];
      Prim[nfield++] = BaryonField[GSiCINum];
      Prim[nfield++] = BaryonField[GSiC2INum];
      Prim[nfield++] = BaryonField[GSiC3INum];
      Prim[nfield++] = BaryonField[GCH3OHINum];
      Prim[nfield++] = BaryonField[GCO2INum];
      Prim[nfield++] = BaryonField[GH2SiOINum];
      Prim[nfield++] = BaryonField[GHNOINum];
      Prim[nfield++] = BaryonField[GN2INum];
      Prim[nfield++] = BaryonField[GH2COINum];
      Prim[nfield++] = BaryonField[GHCNINum];
      Prim[nfield++] = BaryonField[GH2OINum];
      Prim[nfield++] = BaryonField[GNH3INum];
      Prim[nfield++] = BaryonField[SiC3IINum];
      Prim[nfield++] = BaryonField[H2CNINum];
      Prim[nfield++] = BaryonField[GCH4INum];
      Prim[nfield++] = BaryonField[H2NOIINum];
      Prim[nfield++] = BaryonField[H2SiOINum];
      Prim[nfield++] = BaryonField[HeHIINum];
      Prim[nfield++] = BaryonField[HNCOINum];
      Prim[nfield++] = BaryonField[HOCIINum];
      Prim[nfield++] = BaryonField[SiC2IINum];
      Prim[nfield++] = BaryonField[GSiH4INum];
      Prim[nfield++] = BaryonField[SiC2INum];
      Prim[nfield++] = BaryonField[SiC3INum];
      Prim[nfield++] = BaryonField[SiH5IINum];
      Prim[nfield++] = BaryonField[SiH4IINum];
      Prim[nfield++] = BaryonField[SiCIINum];
      Prim[nfield++] = BaryonField[O2HINum];
      Prim[nfield++] = BaryonField[SiCINum];
      Prim[nfield++] = BaryonField[NO2INum];
      Prim[nfield++] = BaryonField[SiH3IINum];
      Prim[nfield++] = BaryonField[SiH2IINum];
      Prim[nfield++] = BaryonField[OCNINum];
      Prim[nfield++] = BaryonField[SiH2INum];
      Prim[nfield++] = BaryonField[SiOHIINum];
      Prim[nfield++] = BaryonField[SiHIINum];
      Prim[nfield++] = BaryonField[SiH4INum];
      Prim[nfield++] = BaryonField[SiHINum];
      Prim[nfield++] = BaryonField[SiH3INum];
      Prim[nfield++] = BaryonField[SiOIINum];
      Prim[nfield++] = BaryonField[HCO2IINum];
      Prim[nfield++] = BaryonField[HNOINum];
      Prim[nfield++] = BaryonField[CH3OHINum];
      Prim[nfield++] = BaryonField[MgINum];
      Prim[nfield++] = BaryonField[MgIINum];
      Prim[nfield++] = BaryonField[CH4IINum];
      Prim[nfield++] = BaryonField[SiOINum];
      Prim[nfield++] = BaryonField[CNIINum];
      Prim[nfield++] = BaryonField[HCNHIINum];
      Prim[nfield++] = BaryonField[N2HIINum];
      Prim[nfield++] = BaryonField[O2HIINum];
      Prim[nfield++] = BaryonField[SiIINum];
      Prim[nfield++] = BaryonField[SiINum];
      Prim[nfield++] = BaryonField[HNCINum];
      Prim[nfield++] = BaryonField[HNOIINum];
      Prim[nfield++] = BaryonField[N2IINum];
      Prim[nfield++] = BaryonField[H3COIINum];
      Prim[nfield++] = BaryonField[CH4INum];
      Prim[nfield++] = BaryonField[COIINum];
      Prim[nfield++] = BaryonField[NH3INum];
      Prim[nfield++] = BaryonField[CH3INum];
      Prim[nfield++] = BaryonField[CO2INum];
      Prim[nfield++] = BaryonField[NIINum];
      Prim[nfield++] = BaryonField[OIINum];
      Prim[nfield++] = BaryonField[HCNIINum];
      Prim[nfield++] = BaryonField[NH2IINum];
      Prim[nfield++] = BaryonField[NHIINum];
      Prim[nfield++] = BaryonField[O2IINum];
      Prim[nfield++] = BaryonField[CH3IINum];
      Prim[nfield++] = BaryonField[NH2INum];
      Prim[nfield++] = BaryonField[CH2IINum];
      Prim[nfield++] = BaryonField[H2OIINum];
      Prim[nfield++] = BaryonField[NH3IINum];
      Prim[nfield++] = BaryonField[NOIINum];
      Prim[nfield++] = BaryonField[H3OIINum];
      Prim[nfield++] = BaryonField[N2INum];
      Prim[nfield++] = BaryonField[CIINum];
      Prim[nfield++] = BaryonField[HCNINum];
      Prim[nfield++] = BaryonField[CHIINum];
      Prim[nfield++] = BaryonField[CH2INum];
      Prim[nfield++] = BaryonField[H2COIINum];
      Prim[nfield++] = BaryonField[NHINum];
      Prim[nfield++] = BaryonField[OHIINum];
      Prim[nfield++] = BaryonField[CNINum];
      Prim[nfield++] = BaryonField[H2COINum];
      Prim[nfield++] = BaryonField[HCOINum];
      Prim[nfield++] = BaryonField[CHINum];
      Prim[nfield++] = BaryonField[H3IINum];
      Prim[nfield++] = BaryonField[NOINum];
      Prim[nfield++] = BaryonField[NINum];
      Prim[nfield++] = BaryonField[OHINum];
      Prim[nfield++] = BaryonField[O2INum];
      Prim[nfield++] = BaryonField[CINum];
      Prim[nfield++] = BaryonField[HCOIINum];
      Prim[nfield++] = BaryonField[H2OINum];
      Prim[nfield++] = BaryonField[OINum];
      Prim[nfield++] = BaryonField[COINum];
      
    }
#endif

  } // ENDIF MultiSpecies

  /* Add the colours (NColor is determined in EvolveLevel) */  

  int SNColourNum, MetalNum, MetalIaNum, MetalIINum, MBHColourNum, Galaxy1ColourNum, 
    Galaxy2ColourNum; 

  if (this->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum, MetalIINum, MBHColourNum, 
				 Galaxy1ColourNum, Galaxy2ColourNum) == FAIL) {
    fprintf(stderr, "Error in grid->IdentifyColourFields.\n");
    return FAIL;
  }
  
  if (MetalNum != -1) {
    Prim[nfield++] = BaryonField[MetalNum];
    if (StarMakerTypeIaSNe)
      Prim[nfield++] = BaryonField[MetalIaNum];
    if (StarMakerTypeIISNeMetalField)
      Prim[nfield++] = BaryonField[MetalIINum];
    if (MultiMetals || TestProblemData.MultiMetals) {
      Prim[nfield++] = BaryonField[MetalNum+1];
      Prim[nfield++] = BaryonField[MetalNum+2];
    }
  }

  if (SNColourNum      != -1) Prim[nfield++] = BaryonField[SNColourNum];  
  /*   //##### These fields are currently not being used and only causing interpolation problems
  if (MBHColourNum     != -1) Prim[nfield++] = BaryonField[MBHColourNum];
  if (Galaxy1ColourNum != -1) Prim[nfield++] = BaryonField[Galaxy1ColourNum];
  if (Galaxy2ColourNum != -1) Prim[nfield++] = BaryonField[Galaxy2ColourNum];
  */

  /* Convert the species and color fields into mass fractions */

  for (dim = 0, size = 1; dim < GridRank; dim++)
    size *= GridDimension[dim];

  if (ReturnMassFractions)  
    for (n = n0; n < nfield; n++)
      for (i = 0; i < size; i++) 
	Prim[n][i] /= Prim[iden][i];

  //  fprintf(stdout, "grid::ReturnHydroRKPointers: nfield = %"ISYM"\n", nfield);  

  return SUCCESS;

}