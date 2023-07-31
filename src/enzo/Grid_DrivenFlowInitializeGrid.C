/***********************************************************************
/
/  GRID CLASS: DrivenFlowInitializeGrid
/
/  written by: Wolfram Schmidt
/  date:       May, 2005
/  modified1:  Sep, 2014: modified to support enzo 2.4 // P. Grete
/
/  PURPOSE: Initializes grid for a driven flow simulation
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include "preincludes.h"
#ifdef USE_NAUNET
#include "naunet_enzo.h"
#endif
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
              float *TemperatureUnits, float *TimeUnits,
              float *VelocityUnits, FLOAT *MassUnits, FLOAT Time);

int grid::DrivenFlowInitializeGrid(float DrivenFlowDensity,
                   float DrivenFlowPressure, float DrivenFlowMagField, int SetBaryonFields)
{
  /* create fields */

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;

  int vel = NumberOfBaryonFields;

  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;

  if (EquationOfState == 0)
    FieldType[NumberOfBaryonFields++] = TotalEnergy;

  if (DualEnergyFormalism)
      FieldType[NumberOfBaryonFields++] = InternalEnergy;

  if (UseMHD) {
    FieldType[NumberOfBaryonFields++] = Bfield1;
    FieldType[NumberOfBaryonFields++] = Bfield2;
    FieldType[NumberOfBaryonFields++] = Bfield3;
  }
  if (HydroMethod == MHD_RK){
    FieldType[NumberOfBaryonFields++] = PhiField;
    
    if(UsePoissonDivergenceCleaning)
      FieldType[NumberOfBaryonFields++] = Phi_pField;
  }

  int accel = NumberOfBaryonFields;

  FieldType[NumberOfBaryonFields++] = DrivingField1;
  FieldType[NumberOfBaryonFields++] = DrivingField2;
  FieldType[NumberOfBaryonFields++] = DrivingField3;

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum,  kphHINum, gammaNum, kphHeINum,
      kphHeIINum, kdissH2INum, RPresNum1, RPresNum2, RPresNum3;
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
  int KspNum;
#endif

  if (MultiSpecies) {
    FieldType[DeNum    = NumberOfBaryonFields++] = ElectronDensity;
    FieldType[HINum    = NumberOfBaryonFields++] = HIDensity;
    FieldType[HIINum   = NumberOfBaryonFields++] = HIIDensity;
    FieldType[HeINum   = NumberOfBaryonFields++] = HeIDensity;
    FieldType[HeIINum  = NumberOfBaryonFields++] = HeIIDensity;
    FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
    if (MultiSpecies > 1) {
      FieldType[HMNum    = NumberOfBaryonFields++] = HMDensity;
      FieldType[H2INum   = NumberOfBaryonFields++] = H2IDensity;
      FieldType[H2IINum  = NumberOfBaryonFields++] = H2IIDensity;
    }
    if (MultiSpecies > 2) {
      FieldType[DINum   = NumberOfBaryonFields++] = DIDensity;
      FieldType[DIINum  = NumberOfBaryonFields++] = DIIDensity;
      FieldType[HDINum  = NumberOfBaryonFields++] = HDIDensity;
    }
#ifdef USE_NAUNET
    if (MultiSpecies == NAUNET_SPECIES) {
      FieldType[GH2CNINum    = NumberOfBaryonFields++] = GH2CNIDensity;
      FieldType[GHNCINum     = NumberOfBaryonFields++] = GHNCIDensity;
      FieldType[GNO2INum     = NumberOfBaryonFields++] = GNO2IDensity;
      FieldType[GSiOINum     = NumberOfBaryonFields++] = GSiOIDensity;
      FieldType[GCOINum      = NumberOfBaryonFields++] = GCOIDensity;
      FieldType[GHNCOINum    = NumberOfBaryonFields++] = GHNCOIDensity;
      FieldType[GMgINum      = NumberOfBaryonFields++] = GMgIDensity;
      FieldType[GNOINum      = NumberOfBaryonFields++] = GNOIDensity;
      FieldType[GO2INum      = NumberOfBaryonFields++] = GO2IDensity;
      FieldType[GO2HINum     = NumberOfBaryonFields++] = GO2HIDensity;
      FieldType[GSiCINum     = NumberOfBaryonFields++] = GSiCIDensity;
      FieldType[GSiC2INum    = NumberOfBaryonFields++] = GSiC2IDensity;
      FieldType[GSiC3INum    = NumberOfBaryonFields++] = GSiC3IDensity;
      FieldType[GCH3OHINum   = NumberOfBaryonFields++] = GCH3OHIDensity;
      FieldType[GCO2INum     = NumberOfBaryonFields++] = GCO2IDensity;
      FieldType[GH2SiOINum   = NumberOfBaryonFields++] = GH2SiOIDensity;
      FieldType[GHNOINum     = NumberOfBaryonFields++] = GHNOIDensity;
      FieldType[GN2INum      = NumberOfBaryonFields++] = GN2IDensity;
      FieldType[GH2COINum    = NumberOfBaryonFields++] = GH2COIDensity;
      FieldType[GHCNINum     = NumberOfBaryonFields++] = GHCNIDensity;
      FieldType[GH2OINum     = NumberOfBaryonFields++] = GH2OIDensity;
      FieldType[GNH3INum     = NumberOfBaryonFields++] = GNH3IDensity;
      FieldType[SiC3IINum    = NumberOfBaryonFields++] = SiC3IIDensity;
      FieldType[H2CNINum     = NumberOfBaryonFields++] = H2CNIDensity;
      FieldType[GCH4INum     = NumberOfBaryonFields++] = GCH4IDensity;
      FieldType[H2NOIINum    = NumberOfBaryonFields++] = H2NOIIDensity;
      FieldType[H2SiOINum    = NumberOfBaryonFields++] = H2SiOIDensity;
      FieldType[HeHIINum     = NumberOfBaryonFields++] = HeHIIDensity;
      FieldType[HNCOINum     = NumberOfBaryonFields++] = HNCOIDensity;
      FieldType[HOCIINum     = NumberOfBaryonFields++] = HOCIIDensity;
      FieldType[SiC2IINum    = NumberOfBaryonFields++] = SiC2IIDensity;
      FieldType[GSiH4INum    = NumberOfBaryonFields++] = GSiH4IDensity;
      FieldType[SiC2INum     = NumberOfBaryonFields++] = SiC2IDensity;
      FieldType[SiC3INum     = NumberOfBaryonFields++] = SiC3IDensity;
      FieldType[SiH5IINum    = NumberOfBaryonFields++] = SiH5IIDensity;
      FieldType[SiH4IINum    = NumberOfBaryonFields++] = SiH4IIDensity;
      FieldType[SiCIINum     = NumberOfBaryonFields++] = SiCIIDensity;
      FieldType[O2HINum      = NumberOfBaryonFields++] = O2HIDensity;
      FieldType[SiCINum      = NumberOfBaryonFields++] = SiCIDensity;
      FieldType[NO2INum      = NumberOfBaryonFields++] = NO2IDensity;
      FieldType[SiH3IINum    = NumberOfBaryonFields++] = SiH3IIDensity;
      FieldType[SiH2IINum    = NumberOfBaryonFields++] = SiH2IIDensity;
      FieldType[OCNINum      = NumberOfBaryonFields++] = OCNIDensity;
      FieldType[SiH2INum     = NumberOfBaryonFields++] = SiH2IDensity;
      FieldType[SiOHIINum    = NumberOfBaryonFields++] = SiOHIIDensity;
      FieldType[SiHIINum     = NumberOfBaryonFields++] = SiHIIDensity;
      FieldType[SiH4INum     = NumberOfBaryonFields++] = SiH4IDensity;
      FieldType[SiHINum      = NumberOfBaryonFields++] = SiHIDensity;
      FieldType[SiH3INum     = NumberOfBaryonFields++] = SiH3IDensity;
      FieldType[SiOIINum     = NumberOfBaryonFields++] = SiOIIDensity;
      FieldType[HCO2IINum    = NumberOfBaryonFields++] = HCO2IIDensity;
      FieldType[HNOINum      = NumberOfBaryonFields++] = HNOIDensity;
      FieldType[CH3OHINum    = NumberOfBaryonFields++] = CH3OHIDensity;
      FieldType[MgINum       = NumberOfBaryonFields++] = MgIDensity;
      FieldType[MgIINum      = NumberOfBaryonFields++] = MgIIDensity;
      FieldType[CH4IINum     = NumberOfBaryonFields++] = CH4IIDensity;
      FieldType[SiOINum      = NumberOfBaryonFields++] = SiOIDensity;
      FieldType[CNIINum      = NumberOfBaryonFields++] = CNIIDensity;
      FieldType[HCNHIINum    = NumberOfBaryonFields++] = HCNHIIDensity;
      FieldType[N2HIINum     = NumberOfBaryonFields++] = N2HIIDensity;
      FieldType[O2HIINum     = NumberOfBaryonFields++] = O2HIIDensity;
      FieldType[SiIINum      = NumberOfBaryonFields++] = SiIIDensity;
      FieldType[SiINum       = NumberOfBaryonFields++] = SiIDensity;
      FieldType[HNCINum      = NumberOfBaryonFields++] = HNCIDensity;
      FieldType[HNOIINum     = NumberOfBaryonFields++] = HNOIIDensity;
      FieldType[N2IINum      = NumberOfBaryonFields++] = N2IIDensity;
      FieldType[H3COIINum    = NumberOfBaryonFields++] = H3COIIDensity;
      FieldType[CH4INum      = NumberOfBaryonFields++] = CH4IDensity;
      FieldType[COIINum      = NumberOfBaryonFields++] = COIIDensity;
      FieldType[NH3INum      = NumberOfBaryonFields++] = NH3IDensity;
      FieldType[CH3INum      = NumberOfBaryonFields++] = CH3IDensity;
      FieldType[CO2INum      = NumberOfBaryonFields++] = CO2IDensity;
      FieldType[NIINum       = NumberOfBaryonFields++] = NIIDensity;
      FieldType[OIINum       = NumberOfBaryonFields++] = OIIDensity;
      FieldType[HCNIINum     = NumberOfBaryonFields++] = HCNIIDensity;
      FieldType[NH2IINum     = NumberOfBaryonFields++] = NH2IIDensity;
      FieldType[NHIINum      = NumberOfBaryonFields++] = NHIIDensity;
      FieldType[O2IINum      = NumberOfBaryonFields++] = O2IIDensity;
      FieldType[CH3IINum     = NumberOfBaryonFields++] = CH3IIDensity;
      FieldType[NH2INum      = NumberOfBaryonFields++] = NH2IDensity;
      FieldType[CH2IINum     = NumberOfBaryonFields++] = CH2IIDensity;
      FieldType[H2OIINum     = NumberOfBaryonFields++] = H2OIIDensity;
      FieldType[NH3IINum     = NumberOfBaryonFields++] = NH3IIDensity;
      FieldType[NOIINum      = NumberOfBaryonFields++] = NOIIDensity;
      FieldType[H3OIINum     = NumberOfBaryonFields++] = H3OIIDensity;
      FieldType[N2INum       = NumberOfBaryonFields++] = N2IDensity;
      FieldType[CIINum       = NumberOfBaryonFields++] = CIIDensity;
      FieldType[HCNINum      = NumberOfBaryonFields++] = HCNIDensity;
      FieldType[CHIINum      = NumberOfBaryonFields++] = CHIIDensity;
      FieldType[CH2INum      = NumberOfBaryonFields++] = CH2IDensity;
      FieldType[H2COIINum    = NumberOfBaryonFields++] = H2COIIDensity;
      FieldType[NHINum       = NumberOfBaryonFields++] = NHIDensity;
      FieldType[OHIINum      = NumberOfBaryonFields++] = OHIIDensity;
      FieldType[CNINum       = NumberOfBaryonFields++] = CNIDensity;
      FieldType[H2COINum     = NumberOfBaryonFields++] = H2COIDensity;
      FieldType[HCOINum      = NumberOfBaryonFields++] = HCOIDensity;
      FieldType[CHINum       = NumberOfBaryonFields++] = CHIDensity;
      FieldType[H3IINum      = NumberOfBaryonFields++] = H3IIDensity;
      FieldType[NOINum       = NumberOfBaryonFields++] = NOIDensity;
      FieldType[NINum        = NumberOfBaryonFields++] = NIDensity;
      FieldType[OHINum       = NumberOfBaryonFields++] = OHIDensity;
      FieldType[O2INum       = NumberOfBaryonFields++] = O2IDensity;
      FieldType[CINum        = NumberOfBaryonFields++] = CIDensity;
      FieldType[HCOIINum     = NumberOfBaryonFields++] = HCOIIDensity;
      FieldType[H2OINum      = NumberOfBaryonFields++] = H2OIDensity;
      FieldType[OINum        = NumberOfBaryonFields++] = OIDensity;
      FieldType[COINum       = NumberOfBaryonFields++] = COIDensity;

      FieldType[KspNum = NumberOfBaryonFields++] = Ksputtering;
    }
#endif
  }


  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (!SetBaryonFields)
    return SUCCESS;

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
        VelocityUnits;
  FLOAT MassUnits = 1;

  int size = 1;

  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];


  this->AllocateGrids();

  /* set density, total energy and velocity in problem dimension */

  float Energy = DrivenFlowPressure / ((Gamma-1.0) * DrivenFlowDensity);

  for (int i = 0; i < size; i++) {
    BaryonField[iden][i] = DrivenFlowDensity;
    BaryonField[ietot][i] = Energy;
  }

  if (DualEnergyFormalism) ///CF (internal energy = total energy, because v=0)
    for (int i = 0; i < size; i++) {
      BaryonField[ieint][i] = Energy;
    }
  printf("DrivenFlowInitializeGrid %"FSYM" %"FSYM" %"FSYM" %"FSYM"\n",
    DrivenFlowDensity,DrivenFlowPressure,Gamma,Energy);
  
  if (HydroMethod == MHD_RK) {
      for (int i = 0; i < size; i++) {
          BaryonField[iBx  ][i]  = DrivenFlowMagField;
      }
      Energy += 0.5 * pow(DrivenFlowMagField,2) / DrivenFlowDensity;
  }
  
  if ( UseMHDCT ){
    for ( int i = 0; i < MagneticSize[0]; i++){
       MagneticField[0][i] = DrivenFlowMagField;
    }
    Energy += 0.5 * pow(DrivenFlowMagField,2) / DrivenFlowDensity;
  this->CenterMagneticField();
  }
  
  if (EquationOfState == 0)
    for( int i = 0; i < size; i++)
      BaryonField[ietot][i] = Energy;

#ifdef USE_NAUNET
  if (MultiSpecies == NAUNET_SPECIES) {
    for( int i = 0; i < size; i++) {
      for (int speciesNum = DeNum; speciesNum <= COINum; speciesNum ++) {
        BaryonField[speciesNum][i] = 1e-40*BaryonField[0][i] / 1.4;
        // printf("speciesNum: %d\n", speciesNum);
      }
      if (1){
        /* set your preferable initial abundances */
        BaryonField[H2INum][i] = 1.00e+0 * BaryonField[0][i] / 1.4;
        BaryonField[HINum][i] = 5.00e-5 * BaryonField[0][i] / 1.4;
        BaryonField[HeINum][i] = 3.90e-1 * BaryonField[0][i] / 1.4;
        BaryonField[NINum][i] = 14.0 * 7.5e-5 * BaryonField[0][i] / 1.4;
        BaryonField[OINum][i] = 16.0 * 1.8e-4 * BaryonField[0][i] / 1.4;
        BaryonField[COINum][i] = 28.0 * 1.4e-4 * BaryonField[0][i] / 1.4;
        BaryonField[MgINum][i] = 24.0 * 7.0e-9 * BaryonField[0][i] / 1.4;
        BaryonField[SiINum][i] = 28.0 * 8.0e-9 * BaryonField[0][i] / 1.4;
        // BaryonField[GSiOINum][n] = 44.0 * 8.0e-9 * BaryonField[0][n] / 1.4;
      }
      else {
        ENZO_FAIL("Error in CollidingCloudInitialize[Sub]Grid: species table is not implemented.");
        // read in table
        // for (int abNum=0; abNum<NSpecies+1; abNum++) {
        //   int speciesNum = speciesMap[abNum];
        //   if (speciesNum != -1) {
        //     BaryonField[speciesNum][i] = PrestellarCoreInitAbundance[abNum] 
        //                                * PrestellarCoreMoleMass[abNum] 
        //                                * BaryonField[0][i] / 1.40045;
        //   }
        // }
      }
    }
  }
#endif


  return SUCCESS;
}
