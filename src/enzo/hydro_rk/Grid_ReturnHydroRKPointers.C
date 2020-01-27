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

#include <stdio.h>
#include <math.h>
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
  int B1Num, B2Num, B3Num, PhiNum;
  int HeIIINum, HMNum, DINum, DIINum, HDINum;
      int DeNum, CHINum, OINum, HNCINum, HCNINum, H2INum,
 CINum, HINum, H2OINum, OHINum, O2INum, CH2INum,
 H2COINum, HCOINum, MGINum, NH3INum, NOINum,
 CNINum, COINum, N2INum, NH2INum, CH3INum,
 CH4INum, NINum, NHINum, HeINum, HNOINum,
 CH3OHINum, CO2INum, H2CNINum, HNCOINum, NO2INum,
 O2HINum, OCNINum, CH3OH_DUSTINum, HNCO_DUSTINum,
 H2CO_DUSTINum, CH4_DUSTINum, CO_DUSTINum,
 H2O_DUSTINum, NO_DUSTINum, CO2_DUSTINum,
 N2_DUSTINum, HCN_DUSTINum, NH3_DUSTINum,
 O2_DUSTINum, NO2_DUSTINum, HNO_DUSTINum,
 O2H_DUSTINum, H2CN_DUSTINum, MG_DUSTINum,
 HNC_DUSTINum, E_DUSTINum, HCOIINum, HIINum,
 HOCIINum, CIINum, CH2IINum, CHIINum, H2COIINum,
 MGIINum, NH3IINum, NOIINum, CNIINum, COIINum,
 N2IINum, O2IINum, H2OIINum, NH2IINum, OIINum,
 OHIINum, CH3IINum, CH4IINum, NIINum, HCNIINum,
 NHIINum, H2IINum, HeIINum, HNOIINum, H2NOIINum,
 H3IINum, H3COIINum, H3OIINum, HCNHIINum,
 HCO2IINum, HeHIINum, N2HIINum, O2HIINum;

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
  if (MultiSpecies)
    if (IdentifySpeciesFieldsKrome(
  DeNum, CHINum, OINum, HNCINum, HCNINum, H2INum,
 CINum, HINum, H2OINum, OHINum, O2INum, CH2INum,
 H2COINum, HCOINum, MGINum, NH3INum, NOINum,
 CNINum, COINum, N2INum, NH2INum, CH3INum,
 CH4INum, NINum, NHINum, HeINum, HNOINum,
 CH3OHINum, CO2INum, H2CNINum, HNCOINum, NO2INum,
 O2HINum, OCNINum, CH3OH_DUSTINum, HNCO_DUSTINum,
 H2CO_DUSTINum, CH4_DUSTINum, CO_DUSTINum,
 H2O_DUSTINum, NO_DUSTINum, CO2_DUSTINum,
 N2_DUSTINum, HCN_DUSTINum, NH3_DUSTINum,
 O2_DUSTINum, NO2_DUSTINum, HNO_DUSTINum,
 O2H_DUSTINum, H2CN_DUSTINum, MG_DUSTINum,
 HNC_DUSTINum, E_DUSTINum, HCOIINum, HIINum,
 HOCIINum, CIINum, CH2IINum, CHIINum, H2COIINum,
 MGIINum, NH3IINum, NOIINum, CNIINum, COIINum,
 N2IINum, O2IINum, H2OIINum, NH2IINum, OIINum,
 OHIINum, CH3IINum, CH4IINum, NIINum, HCNIINum,
 NHIINum, H2IINum, HeIINum, HNOIINum, H2NOIINum,
 H3IINum, H3COIINum, H3OIINum, HCNHIINum,
 HCO2IINum, HeHIINum, N2HIINum, O2HIINum
  ) == FAIL) {
            ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
    }


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
      Prim[nfield++] =  BaryonField[CHINum         ];
      Prim[nfield++] =  BaryonField[OINum          ];
      Prim[nfield++] =  BaryonField[HNCINum        ];
      Prim[nfield++] =  BaryonField[HCNINum        ];
      Prim[nfield++] =  BaryonField[CINum          ];
      Prim[nfield++] =  BaryonField[H2OINum        ];
      Prim[nfield++] =  BaryonField[OHINum         ];
      Prim[nfield++] =  BaryonField[O2INum         ];
      Prim[nfield++] =  BaryonField[CH2INum        ];
      Prim[nfield++] =  BaryonField[H2COINum       ];
      Prim[nfield++] =  BaryonField[HCOINum        ];
      Prim[nfield++] =  BaryonField[MGINum         ];
      Prim[nfield++] =  BaryonField[NH3INum        ];
      Prim[nfield++] =  BaryonField[NOINum         ];
      Prim[nfield++] =  BaryonField[CNINum         ];
      Prim[nfield++] =  BaryonField[COINum         ];
      Prim[nfield++] =  BaryonField[N2INum         ];
      Prim[nfield++] =  BaryonField[NH2INum        ];
      Prim[nfield++] =  BaryonField[CH3INum        ];
      Prim[nfield++] =  BaryonField[CH4INum        ];
      Prim[nfield++] =  BaryonField[NINum          ];
      Prim[nfield++] =  BaryonField[NHINum         ];
      Prim[nfield++] =  BaryonField[HNOINum        ];
      Prim[nfield++] =  BaryonField[CH3OHINum      ];
      Prim[nfield++] =  BaryonField[CO2INum        ];
      Prim[nfield++] =  BaryonField[H2CNINum       ];
      Prim[nfield++] =  BaryonField[HNCOINum       ];
      Prim[nfield++] =  BaryonField[NO2INum        ];
      Prim[nfield++] =  BaryonField[O2HINum        ];
      Prim[nfield++] =  BaryonField[OCNINum        ];
      Prim[nfield++] =  BaryonField[CH3OH_DUSTINum ];
      Prim[nfield++] =  BaryonField[HNCO_DUSTINum  ];
      Prim[nfield++] =  BaryonField[H2CO_DUSTINum  ];
      Prim[nfield++] =  BaryonField[CH4_DUSTINum   ];
      Prim[nfield++] =  BaryonField[CO_DUSTINum    ];
      Prim[nfield++] =  BaryonField[H2O_DUSTINum   ];
      Prim[nfield++] =  BaryonField[NO_DUSTINum    ];
      Prim[nfield++] =  BaryonField[CO2_DUSTINum   ];
      Prim[nfield++] =  BaryonField[N2_DUSTINum    ];
      Prim[nfield++] =  BaryonField[HCN_DUSTINum   ];
      Prim[nfield++] =  BaryonField[NH3_DUSTINum   ];
      Prim[nfield++] =  BaryonField[O2_DUSTINum    ];
      Prim[nfield++] =  BaryonField[NO2_DUSTINum   ];
      Prim[nfield++] =  BaryonField[HNO_DUSTINum   ];
      Prim[nfield++] =  BaryonField[O2H_DUSTINum   ];
      Prim[nfield++] =  BaryonField[H2CN_DUSTINum  ];
      Prim[nfield++] =  BaryonField[MG_DUSTINum    ];
      Prim[nfield++] =  BaryonField[HNC_DUSTINum   ];
      Prim[nfield++] =  BaryonField[E_DUSTINum     ];
      Prim[nfield++] =  BaryonField[HCOIINum       ];
      Prim[nfield++] =  BaryonField[HOCIINum       ];
      Prim[nfield++] =  BaryonField[CIINum         ];
      Prim[nfield++] =  BaryonField[CH2IINum       ];
      Prim[nfield++] =  BaryonField[CHIINum        ];
      Prim[nfield++] =  BaryonField[H2COIINum      ];
      Prim[nfield++] =  BaryonField[MGIINum        ];
      Prim[nfield++] =  BaryonField[NH3IINum       ];
      Prim[nfield++] =  BaryonField[NOIINum        ];
      Prim[nfield++] =  BaryonField[CNIINum        ];
      Prim[nfield++] =  BaryonField[COIINum        ];
      Prim[nfield++] =  BaryonField[N2IINum        ];
      Prim[nfield++] =  BaryonField[O2IINum        ];
      Prim[nfield++] =  BaryonField[H2OIINum       ];
      Prim[nfield++] =  BaryonField[NH2IINum       ];
      Prim[nfield++] =  BaryonField[OIINum         ];
      Prim[nfield++] =  BaryonField[OHIINum        ];
      Prim[nfield++] =  BaryonField[CH3IINum       ];
      Prim[nfield++] =  BaryonField[CH4IINum       ];
      Prim[nfield++] =  BaryonField[NIINum         ];
      Prim[nfield++] =  BaryonField[HCNIINum       ];
      Prim[nfield++] =  BaryonField[NHIINum        ];
      Prim[nfield++] =  BaryonField[HNOIINum       ];
      Prim[nfield++] =  BaryonField[H2NOIINum      ];
      Prim[nfield++] =  BaryonField[H3IINum        ];
      Prim[nfield++] =  BaryonField[H3COIINum      ];
      Prim[nfield++] =  BaryonField[H3OIINum       ];
      Prim[nfield++] =  BaryonField[HCNHIINum      ];
      Prim[nfield++] =  BaryonField[HCO2IINum      ];
      Prim[nfield++] =  BaryonField[HeHIINum       ];
      Prim[nfield++] =  BaryonField[N2HIINum       ];
      Prim[nfield++] =  BaryonField[O2HIINum       ];

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
