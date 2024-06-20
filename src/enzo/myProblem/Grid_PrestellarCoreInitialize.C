/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR PRESTELLAR CORE) 
/
/  written by: Chia-Jung Hsu
/  date:       March, 2019
/  modified1:  
/
/  PURPOSE: Sets the density and turbulence in a prestellar core.
/
/  RETURNS: FAIL or SUCCESS
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
#include "Grid.h"

void TurbulenceGenerator_FFTW(const int field_type, const int gridsize, const int randomSeed,
        const float power_turb, float* vfield);

void Turbulence_Generator(float **vel, int dim0, int dim1, int dim2, int ind, 
        float kmin, float kmax, float dk,
        FLOAT **LeftEdge, FLOAT **CellWidth, int seed);

int grid::PrestellarCoreInitializeGrid(
            float PrestellarCoreRadius,
            float PrestellarCoreDensity,
            float PrestellarCoreSurfaceDensity,
            float PrestellarCoreDensityJump,
            float PrestellarCoreInternalEnergy,
            float PrestellarCoreAngularVelocity,
            float PrestellarCoreBzField,
            float PrestellarCoreAmbientBzField,
            float PrestellarCoreVelocityDispersion,
            float PrestellarCoreTurbulenceKStart,
            float PrestellarCoreTurbulenceKEnd,
            float PrestellarCoreOPR,
            float PrestellarCoreCODeplete,
            int PrestellarCoreRandomSeed,
            int level,
            int* baseDims,
            float* PrestellarCoreInitAbundance,
            float* PrestellarCoreMoleMass,
            float* Turbulence,
            bool SetBaryonField)
{

  /* declarations */

  int size = 1, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* set fields in the rotating core region: x^2+y^2+z^2 < dr^2. */

  // static int count = 0, count2 = 0;
  int index, jndex, kndex, i, j, k, n;
  float zonex, zoney, zonez, radius;
  float r_f = 0.15 * PrestellarCoreRadius,
        r_s = 0.05 * PrestellarCoreRadius,
        r_c = PrestellarCoreRadius;
  float r2 = PrestellarCoreRadius*PrestellarCoreRadius;
  float xcenter = (DomainRightEdge[0] + DomainLeftEdge[0])/2.0;
  float ycenter = (DomainRightEdge[1] + DomainLeftEdge[1])/2.0;
  float zcenter = (DomainRightEdge[2] + DomainLeftEdge[2])/2.0;
  float dSize   = (DomainRightEdge[0] - DomainLeftEdge[0]);
  float PrestellarCoreOuterDensity = PrestellarCoreDensityJump * PrestellarCoreSurfaceDensity;
  float ambientInternalEnergy = PrestellarCoreInternalEnergy / PrestellarCoreDensityJump;

  int DensNum, TENum, GENum, V1Num, V2Num, V3Num, B1Num, B2Num, B3Num, PhiNum;

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

  NumberOfBaryonFields = 0;
  FieldType[DensNum = NumberOfBaryonFields++] = Density;
  FieldType[TENum   = NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[GENum   = NumberOfBaryonFields++] = InternalEnergy;
  int vel = NumberOfBaryonFields;
  FieldType[V1Num   = NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1 || HydroMethod > 2)
    FieldType[V2Num   = NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2 || HydroMethod > 2)
    FieldType[V3Num   = NumberOfBaryonFields++] = Velocity3;
  if ( UseMHD ) {
    FieldType[B1Num = NumberOfBaryonFields++] = Bfield1;
    FieldType[B2Num = NumberOfBaryonFields++] = Bfield2;
    FieldType[B3Num = NumberOfBaryonFields++] = Bfield3;
    if( HydroMethod == MHD_RK ){
        FieldType[PhiNum = NumberOfBaryonFields++] = PhiField;
    }
    // if (UseDivergenceCleaning) {
    //   FieldType[NumberOfBaryonFields++] = Phi_pField;
    // }
  }

  if (WritePotential)
    FieldType[NumberOfBaryonFields++] = GravPotential;

  int colorfields = NumberOfBaryonFields;
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
      FieldType[GCH3OHINum  = NumberOfBaryonFields++] = GCH3OHIDensity;
      FieldType[GCH4INum  = NumberOfBaryonFields++] = GCH4IDensity;
      FieldType[GCOINum  = NumberOfBaryonFields++] = GCOIDensity;
      FieldType[GCO2INum  = NumberOfBaryonFields++] = GCO2IDensity;
      FieldType[GH2CNINum  = NumberOfBaryonFields++] = GH2CNIDensity;
      FieldType[GH2COINum  = NumberOfBaryonFields++] = GH2COIDensity;
      FieldType[GH2OINum  = NumberOfBaryonFields++] = GH2OIDensity;
      FieldType[GH2SiOINum  = NumberOfBaryonFields++] = GH2SiOIDensity;
      FieldType[GHCNINum  = NumberOfBaryonFields++] = GHCNIDensity;
      FieldType[GHNCINum  = NumberOfBaryonFields++] = GHNCIDensity;
      FieldType[GHNCOINum  = NumberOfBaryonFields++] = GHNCOIDensity;
      FieldType[GHNOINum  = NumberOfBaryonFields++] = GHNOIDensity;
      FieldType[GMgINum  = NumberOfBaryonFields++] = GMgIDensity;
      FieldType[GN2INum  = NumberOfBaryonFields++] = GN2IDensity;
      FieldType[GNH3INum  = NumberOfBaryonFields++] = GNH3IDensity;
      FieldType[GNOINum  = NumberOfBaryonFields++] = GNOIDensity;
      FieldType[GNO2INum  = NumberOfBaryonFields++] = GNO2IDensity;
      FieldType[GO2INum  = NumberOfBaryonFields++] = GO2IDensity;
      FieldType[GO2HINum  = NumberOfBaryonFields++] = GO2HIDensity;
      FieldType[GSiCINum  = NumberOfBaryonFields++] = GSiCIDensity;
      FieldType[GSiC2INum  = NumberOfBaryonFields++] = GSiC2IDensity;
      FieldType[GSiC3INum  = NumberOfBaryonFields++] = GSiC3IDensity;
      FieldType[GSiH4INum  = NumberOfBaryonFields++] = GSiH4IDensity;
      FieldType[GSiOINum  = NumberOfBaryonFields++] = GSiOIDensity;
      FieldType[CINum  = NumberOfBaryonFields++] = CIDensity;
      FieldType[CIINum  = NumberOfBaryonFields++] = CIIDensity;
      FieldType[CHINum  = NumberOfBaryonFields++] = CHIDensity;
      FieldType[CHIINum  = NumberOfBaryonFields++] = CHIIDensity;
      FieldType[CH2INum  = NumberOfBaryonFields++] = CH2IDensity;
      FieldType[CH2IINum  = NumberOfBaryonFields++] = CH2IIDensity;
      FieldType[CH3INum  = NumberOfBaryonFields++] = CH3IDensity;
      FieldType[CH3IINum  = NumberOfBaryonFields++] = CH3IIDensity;
      FieldType[CH3OHINum  = NumberOfBaryonFields++] = CH3OHIDensity;
      FieldType[CH4INum  = NumberOfBaryonFields++] = CH4IDensity;
      FieldType[CH4IINum  = NumberOfBaryonFields++] = CH4IIDensity;
      FieldType[CNINum  = NumberOfBaryonFields++] = CNIDensity;
      FieldType[CNIINum  = NumberOfBaryonFields++] = CNIIDensity;
      FieldType[COINum  = NumberOfBaryonFields++] = COIDensity;
      FieldType[COIINum  = NumberOfBaryonFields++] = COIIDensity;
      FieldType[CO2INum  = NumberOfBaryonFields++] = CO2IDensity;
      FieldType[H2CNINum  = NumberOfBaryonFields++] = H2CNIDensity;
      FieldType[H2COINum  = NumberOfBaryonFields++] = H2COIDensity;
      FieldType[H2COIINum  = NumberOfBaryonFields++] = H2COIIDensity;
      FieldType[H2NOIINum  = NumberOfBaryonFields++] = H2NOIIDensity;
      FieldType[H2OINum  = NumberOfBaryonFields++] = H2OIDensity;
      FieldType[H2OIINum  = NumberOfBaryonFields++] = H2OIIDensity;
      FieldType[H2SiOINum  = NumberOfBaryonFields++] = H2SiOIDensity;
      FieldType[H3IINum  = NumberOfBaryonFields++] = H3IIDensity;
      FieldType[H3COIINum  = NumberOfBaryonFields++] = H3COIIDensity;
      FieldType[H3OIINum  = NumberOfBaryonFields++] = H3OIIDensity;
      FieldType[HCNINum  = NumberOfBaryonFields++] = HCNIDensity;
      FieldType[HCNIINum  = NumberOfBaryonFields++] = HCNIIDensity;
      FieldType[HCNHIINum  = NumberOfBaryonFields++] = HCNHIIDensity;
      FieldType[HCOINum  = NumberOfBaryonFields++] = HCOIDensity;
      FieldType[HCOIINum  = NumberOfBaryonFields++] = HCOIIDensity;
      FieldType[HCO2IINum  = NumberOfBaryonFields++] = HCO2IIDensity;
      FieldType[HeHIINum  = NumberOfBaryonFields++] = HeHIIDensity;
      FieldType[HNCINum  = NumberOfBaryonFields++] = HNCIDensity;
      FieldType[HNCOINum  = NumberOfBaryonFields++] = HNCOIDensity;
      FieldType[HNOINum  = NumberOfBaryonFields++] = HNOIDensity;
      FieldType[HNOIINum  = NumberOfBaryonFields++] = HNOIIDensity;
      FieldType[HOCIINum  = NumberOfBaryonFields++] = HOCIIDensity;
      FieldType[MgINum  = NumberOfBaryonFields++] = MgIDensity;
      FieldType[MgIINum  = NumberOfBaryonFields++] = MgIIDensity;
      FieldType[NINum  = NumberOfBaryonFields++] = NIDensity;
      FieldType[NIINum  = NumberOfBaryonFields++] = NIIDensity;
      FieldType[N2INum  = NumberOfBaryonFields++] = N2IDensity;
      FieldType[N2IINum  = NumberOfBaryonFields++] = N2IIDensity;
      FieldType[N2HIINum  = NumberOfBaryonFields++] = N2HIIDensity;
      FieldType[NHINum  = NumberOfBaryonFields++] = NHIDensity;
      FieldType[NHIINum  = NumberOfBaryonFields++] = NHIIDensity;
      FieldType[NH2INum  = NumberOfBaryonFields++] = NH2IDensity;
      FieldType[NH2IINum  = NumberOfBaryonFields++] = NH2IIDensity;
      FieldType[NH3INum  = NumberOfBaryonFields++] = NH3IDensity;
      FieldType[NH3IINum  = NumberOfBaryonFields++] = NH3IIDensity;
      FieldType[NOINum  = NumberOfBaryonFields++] = NOIDensity;
      FieldType[NOIINum  = NumberOfBaryonFields++] = NOIIDensity;
      FieldType[NO2INum  = NumberOfBaryonFields++] = NO2IDensity;
      FieldType[OINum  = NumberOfBaryonFields++] = OIDensity;
      FieldType[OIINum  = NumberOfBaryonFields++] = OIIDensity;
      FieldType[O2INum  = NumberOfBaryonFields++] = O2IDensity;
      FieldType[O2IINum  = NumberOfBaryonFields++] = O2IIDensity;
      FieldType[O2HINum  = NumberOfBaryonFields++] = O2HIDensity;
      FieldType[O2HIINum  = NumberOfBaryonFields++] = O2HIIDensity;
      FieldType[OCNINum  = NumberOfBaryonFields++] = OCNIDensity;
      FieldType[OHINum  = NumberOfBaryonFields++] = OHIDensity;
      FieldType[OHIINum  = NumberOfBaryonFields++] = OHIIDensity;
      FieldType[SiINum  = NumberOfBaryonFields++] = SiIDensity;
      FieldType[SiIINum  = NumberOfBaryonFields++] = SiIIDensity;
      FieldType[SiCINum  = NumberOfBaryonFields++] = SiCIDensity;
      FieldType[SiCIINum  = NumberOfBaryonFields++] = SiCIIDensity;
      FieldType[SiC2INum  = NumberOfBaryonFields++] = SiC2IDensity;
      FieldType[SiC2IINum  = NumberOfBaryonFields++] = SiC2IIDensity;
      FieldType[SiC3INum  = NumberOfBaryonFields++] = SiC3IDensity;
      FieldType[SiC3IINum  = NumberOfBaryonFields++] = SiC3IIDensity;
      FieldType[SiHINum  = NumberOfBaryonFields++] = SiHIDensity;
      FieldType[SiHIINum  = NumberOfBaryonFields++] = SiHIIDensity;
      FieldType[SiH2INum  = NumberOfBaryonFields++] = SiH2IDensity;
      FieldType[SiH2IINum  = NumberOfBaryonFields++] = SiH2IIDensity;
      FieldType[SiH3INum  = NumberOfBaryonFields++] = SiH3IDensity;
      FieldType[SiH3IINum  = NumberOfBaryonFields++] = SiH3IIDensity;
      FieldType[SiH4INum  = NumberOfBaryonFields++] = SiH4IDensity;
      FieldType[SiH4IINum  = NumberOfBaryonFields++] = SiH4IIDensity;
      FieldType[SiH5IINum  = NumberOfBaryonFields++] = SiH5IIDensity;
      FieldType[SiOINum  = NumberOfBaryonFields++] = SiOIDensity;
      FieldType[SiOIINum  = NumberOfBaryonFields++] = SiOIIDensity;
      FieldType[SiOHIINum  = NumberOfBaryonFields++] = SiOHIIDensity;
    }
#endif

  }

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (!SetBaryonField)
    return SUCCESS;

  this->AllocateGrids(); 

  float coreMass1 = 0.0, coreMass2 = 0.0, coreMass3;
  float xlen, ylen, zlen, Vol;
  for (i = 0; i < size; i++) {

    index = i % GridDimension[0];
    jndex = ((i-index) % (GridDimension[0]*GridDimension[1]))/GridDimension[0];
    kndex = (i-index-jndex*GridDimension[0])/GridDimension[0]/GridDimension[1];

    zonex = *(CellLeftEdge[0] + index) + 0.5*(*(CellWidth[0] + index)) - xcenter;
    zoney = *(CellLeftEdge[1] + jndex) + 0.5*(*(CellWidth[1] + jndex)) - ycenter;
    zonez = *(CellLeftEdge[2] + kndex) + 0.5*(*(CellWidth[2] + kndex)) - zcenter;
    radius = sqrt(zonex*zonex + zoney*zoney + zonez*zonez);

    // if (count<1) 
    //     printf("coreDens: %13.7e, outerDens: %13.7e, radius: %13.7e, coreRad: %13.7e\n", 
    //         PrestellarCoreDensity, PrestellarCoreOuterDensity, 
    //         radius, r_c);

    BaryonField[DensNum][i] = PrestellarCoreOuterDensity
        + (PrestellarCoreDensity - PrestellarCoreOuterDensity)
        / (1.0 + pow( (radius/r_f), 1.5))
        * (0.5 - 0.5*tanh( (radius-r_c)/ r_s ));

    xlen = *(CellWidth[0] + index);
    ylen = *(CellWidth[1] + jndex);
    zlen = *(CellWidth[2] + kndex);
    Vol = xlen * ylen * zlen;
    if (!ParallelRootGridIO){
      if (radius < r_c) coreMass1 += BaryonField[DensNum][i] * Vol;
      if (BaryonField[DensNum][i] > PrestellarCoreOuterDensity*10.0) coreMass2 += BaryonField[DensNum][i] * Vol;
      if (BaryonField[DensNum][i] > PrestellarCoreOuterDensity*1.1) coreMass3 += BaryonField[DensNum][i] * Vol;
    }

    // ambient gas temperature is 10 times larger than core
    float transfer = 0.5*tanh( (radius-r_c)/ r_s );
    BaryonField[TENum][i] = ambientInternalEnergy * (0.5 + transfer)
                          + PrestellarCoreInternalEnergy * (0.5 - transfer);
    // BaryonField[TENum][i] = PrestellarCoreInternalEnergy * 10.0;

    if (DualEnergyFormalism) {
      BaryonField[GENum][i] = ambientInternalEnergy * (0.5 + transfer)
                            + PrestellarCoreInternalEnergy * (0.5 - transfer);
    }
    
    if (zonex*zonex + zoney*zoney + zonez*zonez < r2) {
      // BaryonField[TENum][i] /= 10.0;
      BaryonField[V1Num][i] =   PrestellarCoreAngularVelocity*zoney;
      BaryonField[V2Num][i] = - PrestellarCoreAngularVelocity*zonex;
      BaryonField[V3Num][i] = 0.0; // the rotation axis is || to z.

      // if there is no turbulence, calculate kinetic energy now
      if (PrestellarCoreTurbulenceKStart > PrestellarCoreTurbulenceKEnd)
        BaryonField[TENum][i] += (BaryonField[V1Num][i]*BaryonField[V1Num][i] 
                               +  BaryonField[V2Num][i]*BaryonField[V2Num][i]
                               +  BaryonField[V3Num][i]*BaryonField[V3Num][i])/2.0;
    }

    if (UseMHD) {
      radius = sqrt(zonex*zonex + zoney*zoney);
      BaryonField[B1Num][i] = 0.0;
      BaryonField[B2Num][i] = 0.0;
      BaryonField[B3Num][i] = PrestellarCoreAmbientBzField + 
        (PrestellarCoreBzField - PrestellarCoreAmbientBzField)
        / (1.0 + pow( (radius/r_f), 0.5))
        * (0.5 - 0.5*tanh( (radius-r_c)/ r_s ));
      BaryonField[TENum][i] += 0.5*BaryonField[B3Num][i]*BaryonField[B3Num][i]/BaryonField[DensNum][i];
      // printf("zonex: %13.7e, zoney: %13.7e, zonez: %13.7e, Bz: %13.7e \n", zonex, zoney, zonez, BaryonField[B3Num][i]);
    }

#ifdef USE_NAUNET
    if (MultiSpecies == NAUNET_SPECIES) {
      for (int speciesNum = DeNum; speciesNum <= SiOHIINum; speciesNum ++) {
        BaryonField[speciesNum][i] = 1e-40*BaryonField[0][i] / 1.4;
      }
      if (1){
        /* set your preferable initial abundances */
        BaryonField[H2INum][i] = 1.00e+0 * BaryonField[0][i] / 1.4;
        BaryonField[HINum][i] = 5.00e-5 * BaryonField[0][i] / 1.4;
        BaryonField[HeINum][i] = 3.90e-1 * BaryonField[0][i] / 1.4;
        BaryonField[NINum][i] = 14.0 * 7.5e-5 * BaryonField[0][i] / 1.4;
        BaryonField[OINum][i] = 16.0 * 1.8e-4 * BaryonField[0][i] / 1.4;
        BaryonField[COINum][i] = 28.0 * 1.4e-4 * BaryonField[0][i] / 1.4;
        BaryonField[SiINum][i] = 28.0 * 8.0e-9 * BaryonField[0][i] / 1.4;
        BaryonField[MgINum][i] = 24.0 * 7.0e-9 * BaryonField[0][i] / 1.4;
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
#else
    // initial fraction copied from Grid_TurbulenceInitializeGrid.C
    // feel free to change
    float InitialFractionHII = 1.2e-5;
    float InitialFractionHeII = 1.0e-14;
    float InitialFractionHeIII = 1.0e-17;

    if (MultiSpecies) {
      BaryonField[HIINum][n] = InitialFractionHII *
        CoolData.HydrogenFractionByMass * BaryonField[iden][n];
      BaryonField[HeIINum][n] = InitialFractionHeII*
        BaryonField[iden][n] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
      BaryonField[HeIIINum][n] = InitialFractionHeIII*
        BaryonField[iden][n] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
      BaryonField[HeINum][n] =
        (1.0 - CoolData.HydrogenFractionByMass)*BaryonField[iden][n] -
        BaryonField[HeIINum][n] - BaryonField[HeIIINum][n];

  
      BaryonField[HINum][n] =
        CoolData.HydrogenFractionByMass*BaryonField[iden][n]
        - BaryonField[HIINum][n];
      

      /* electron "density": n_e * m_p */
      
      BaryonField[DeNum][n] = BaryonField[HIINum][n] +
        0.25*BaryonField[HeIINum][n] + 0.5*BaryonField[HeIIINum][n];
      

      // for (int speciesNum = DeNum; speciesNum <= HeIIINum; speciesNum ++) {
      //   BaryonField[speciesNum][i] = 1e-20*BaryonField[0][i];
      // }
      // BaryonField[H2INum][i]          = 5.00e-1*BaryonField[0][i] / 1.412;
      // BaryonField[HINum][i]           = 5.00e-1*BaryonField[0][i] / 1.412;
      // BaryonField[HeINum][i]          = 4.00e-1*BaryonField[0][i] / 1.412;
    }
#endif


    // if (count < 1){
    //     printf("zonex: %13.7e, zoney: %13.7e\n", zonex, zoney);
    //     printf("Eint: %13.7e, xvel: %13.7e\n", BaryonField[1][i], BaryonField[2][i]);
    // }
    // count++;
  }

  if (!ParallelRootGridIO){
    printf("Total core mass inside Rc: %13.7e\n", coreMass1);
    printf("Total core mass where rho > 10.0 rho0: %13.7e\n", coreMass2);
    printf("Total core mass where rho > 1.1 rho0: %13.7e\n", coreMass3);
  }

  // float k1 = 65.0/dSize, 
  //       k2 = 73.0/dSize,
  //       //dk = 1.0/dSize;
  //       //k2 = float(int(GridDimension[0]/10.0)), 
  //       dk = max(1.0/dSize,(k2-k1)/13.0);
  float k1 = PrestellarCoreTurbulenceKStart, k2 = PrestellarCoreTurbulenceKEnd, 
        dk = 1.0;
        //dk = max(1.0, (k2-k1)/10.0);


  if (k1 < k2){

    // increase grid dimensional with level to initialize turbulence in higher resolution
    int activesize = 1, turbsize=1;
    for (dim = 0; dim < GridRank; dim++) {
      activesize *= (GridDimension[dim] - 2*NumberOfGhostZones);
      turbsize *= baseDims[dim] * pow(RefineBy, level);
    }
    int gridsize = baseDims[0] * pow(RefineBy, level);
    float *TurbulenceVelocity[3];
    for (dim = 0; dim < GridRank; dim++) {
      TurbulenceVelocity[dim] = new float[activesize];
    }

    // look for the sub-region from the grids of turbulence initialization
    int TurbGridStartIdx[3] = {0},
        TurbGridEndIdx[3] = {0};

    for (dim = 0; dim < GridRank; dim++) {
      TurbGridStartIdx[dim] = (GridLeftEdge[dim] - DomainLeftEdge[dim])/CellWidth[dim][0];
      TurbGridEndIdx[dim]   = TurbGridStartIdx[dim] + (GridDimension[dim]-2*NumberOfGhostZones) - 1;
      printf("StartIndex=%d, EndIndex=%d\n", TurbGridStartIdx[dim], TurbGridEndIdx[dim]);
    }
 
    for (dim = 0; dim < GridRank; dim++)
    for (i = TurbGridStartIdx[0], n=0; i <= TurbGridEndIdx[0]; i++)
    for (j = TurbGridStartIdx[1]; j <= TurbGridEndIdx[1]; j++)
    for (k = TurbGridStartIdx[2]; k <= TurbGridEndIdx[2]; k++, n++){
      int ref = dim + 3*i + 3*gridsize*j + 3*gridsize*gridsize*k;
      TurbulenceVelocity[dim][n] = Turbulence[ref];
    }

    // Old turbulence generator
    // Turbulence_Generator(TurbulenceVelocity, 
    //         GridDimension[0]-2*NumberOfGhostZones,
    //         GridDimension[1]-2*NumberOfGhostZones,
    //         GridDimension[2]-2*NumberOfGhostZones,
    //         4.0, k1, k2, dk,
    //         CellLeftEdge, CellWidth, 65536);    

    n = 0;
    for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++, n++) {
      int igrid = i + GridDimension[0]*(j+k*GridDimension[1]);
      // zonex = *(CellLeftEdge[0] + i) + 0.5*(*(CellWidth[0] + i)) - xcenter;
      // zoney = *(CellLeftEdge[1] + j) + 0.5*(*(CellWidth[1] + j)) - ycenter;
      // zonez = *(CellLeftEdge[2] + k) + 0.5*(*(CellWidth[2] + k)) - zcenter;
      // radius = sqrt(zonex*zonex + zoney*zoney + zonez*zonez);

      BaryonField[V1Num][igrid] += TurbulenceVelocity[0][n];
                                 //* (0.5 - 0.5*tanh( (radius-r_c)/ r_s ));
      BaryonField[V2Num][igrid] += TurbulenceVelocity[1][n];
                                 //* (0.5 - 0.5*tanh( (radius-r_c)/ r_s ));
      BaryonField[V3Num][igrid] += TurbulenceVelocity[2][n];
                                 //* (0.5 - 0.5*tanh( (radius-r_c)/ r_s ));
      // total energy
      BaryonField[TENum][igrid] += (BaryonField[V1Num][igrid]*BaryonField[V1Num][igrid] 
                                 +  BaryonField[V2Num][igrid]*BaryonField[V2Num][igrid]
                                 +  BaryonField[V3Num][igrid]*BaryonField[V3Num][igrid])/2.0;

      // float vel = sqrt(BaryonField[2][igrid]*BaryonField[2][igrid] 
      //           +      BaryonField[3][igrid]*BaryonField[3][igrid]
      //           +      BaryonField[4][igrid]*BaryonField[4][igrid]);
      // printf("xvel: %13.7e, yvel: %13.7e, zvel: %13.7e, vel: %13.7e \n",
      //         BaryonField[2][igrid], BaryonField[3][igrid], BaryonField[4][igrid], vel);
    }}}
  
    int idx = GridStartIndex[0] + GridDimension[0]*(GridStartIndex[1]+GridStartIndex[2]*GridDimension[1]);
    // if(count2<1)
    //   printf("Eint: %13.7e, xvel: %13.7e\n", BaryonField[1][idx], BaryonField[2][idx]);
    // count2 ++;
  
    for (dim = 0; dim < GridRank; dim++)
      delete [] TurbulenceVelocity[dim];

  } // if (k1 < k2), set turbulence
 
  return SUCCESS;
}

