/***********************************************************************
/
/  GRID CLASS (WRAP THE NAUNET CHEMISTRY SOLVER)
/
/  written by: Chia-Jung Hsu
/  date:       2021
/  modified1:
/
/  PURPOSE: Solve chemistry with naunet.
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
// clang-format off
#include "preincludes.h"
#ifdef USE_NAUNET
#include "naunet_enzo.h"
#endif
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "phys_constants.h"

// function prototypes

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);

#ifdef USE_NAUNET

int grid::NaunetWrapper()
{

  if (use_naunet == FALSE)
    return SUCCESS;

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  LCAPERF_START("grid_NaunetWrapper");

  int GCH3OHINum, GCH4INum, GCOINum, GCO2INum, GH2CNINum, GH2COINum, GH2OINum,
      GH2SiOINum, GHCNINum, GHNCINum, GHNCOINum, GHNOINum, GMgINum, GN2INum,
      GNH3INum, GNOINum, GNO2INum, GO2INum, GO2HINum, GSiCINum, GSiC2INum,
      GSiC3INum, GSiH4INum, GSiOINum, CINum, CIINum, CHINum, CHIINum, CH2INum,
      CH2IINum, CH3INum, CH3IINum, CH3OHINum, CH4INum, CH4IINum, CNINum,
      CNIINum, COINum, COIINum, CO2INum, DeNum, HINum, HIINum, H2INum, H2IINum,
      H2CNINum, H2COINum, H2COIINum, H2NOIINum, H2OINum, H2OIINum, H2SiOINum,
      H3IINum, H3COIINum, H3OIINum, HCNINum, HCNIINum, HCNHIINum, HCOINum,
      HCOIINum, HCO2IINum, HeINum, HeIINum, HeHIINum, HNCINum, HNCOINum,
      HNOINum, HNOIINum, HOCIINum, MgINum, MgIINum, NINum, NIINum, N2INum,
      N2IINum, N2HIINum, NHINum, NHIINum, NH2INum, NH2IINum, NH3INum, NH3IINum,
      NOINum, NOIINum, NO2INum, OINum, OIINum, O2INum, O2IINum, O2HINum,
      O2HIINum, OCNINum, OHINum, OHIINum, SiINum, SiIINum, SiCINum, SiCIINum,
      SiC2INum, SiC2IINum, SiC3INum, SiC3IINum, SiHINum, SiHIINum, SiH2INum,
      SiH2IINum, SiH3INum, SiH3IINum, SiH4INum, SiH4IINum, SiH5IINum, SiOINum,
      SiOIINum, SiOHIINum;

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;

  double dt_chem = dtFixed;
  
  // Compute the size of the fields.
 
  int i;
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  Eint32 *g_grid_dimension, *g_grid_start, *g_grid_end;
  g_grid_dimension = new Eint32[GridRank];
  g_grid_start = new Eint32[GridRank];
  g_grid_end = new Eint32[GridRank];
  for (i = 0; i < GridRank; i++) {
    g_grid_dimension[i] = (Eint32) GridDimension[i];
    g_grid_start[i] = (Eint32) GridStartIndex[i];
    g_grid_end[i] = (Eint32) GridEndIndex[i];
  }
 
  // Find fields: density, total energy, velocity1-3.
 
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  // Find Multi-species fields.

  GCH3OHINum = GCH4INum = GCOINum = GCO2INum = GH2CNINum = GH2COINum =
    GH2OINum = GH2SiOINum = GHCNINum = GHNCINum = GHNCOINum = GHNOINum = GMgINum
    = GN2INum = GNH3INum = GNOINum = GNO2INum = GO2INum = GO2HINum = GSiCINum =
    GSiC2INum = GSiC3INum = GSiH4INum = GSiOINum = CINum = CIINum = CHINum =
    CHIINum = CH2INum = CH2IINum = CH3INum = CH3IINum = CH3OHINum = CH4INum =
    CH4IINum = CNINum = CNIINum = COINum = COIINum = CO2INum = DeNum = HINum =
    HIINum = H2INum = H2IINum = H2CNINum = H2COINum = H2COIINum = H2NOIINum =
    H2OINum = H2OIINum = H2SiOINum = H3IINum = H3COIINum = H3OIINum = HCNINum =
    HCNIINum = HCNHIINum = HCOINum = HCOIINum = HCO2IINum = HeINum = HeIINum =
    HeHIINum = HNCINum = HNCOINum = HNOINum = HNOIINum = HOCIINum = MgINum =
    MgIINum = NINum = NIINum = N2INum = N2IINum = N2HIINum = NHINum = NHIINum =
    NH2INum = NH2IINum = NH3INum = NH3IINum = NOINum = NOIINum = NO2INum = OINum
    = OIINum = O2INum = O2IINum = O2HINum = O2HIINum = OCNINum = OHINum =
    OHIINum = SiINum = SiIINum = SiCINum = SiCIINum = SiC2INum = SiC2IINum =
    SiC3INum = SiC3IINum = SiHINum = SiHIINum = SiH2INum = SiH2IINum = SiH3INum
    = SiH3IINum = SiH4INum = SiH4IINum = SiH5IINum = SiOINum = SiOIINum =
    SiOHIINum = 0;
 
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
      ENZO_FAIL("Error in grid->IdentifyNaunetSpeciesFields.\n");
    }
 
  // Get easy to handle pointers for each variable.
 
  float *density     = BaryonField[DensNum];
  float *totalenergy = BaryonField[TENum];
  float *gasenergy   = BaryonField[GENum];
  float *velocity1   = BaryonField[Vel1Num];
  float *velocity2   = BaryonField[Vel2Num];
  float *velocity3   = BaryonField[Vel3Num];

  // Compute the cooling time.

  FLOAT a = 1.0, dadt;
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
        VelocityUnits = 1, TimeUnits = 1, aUnits = 1;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
           &TimeUnits, &VelocityUnits, Time);

  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(Time+0.5*dt_chem, &a, &dadt);
    aUnits = 1.0/(1.0 + InitialRedshift);
  } 
  else if (RadiationFieldRedshift > -1){
    a        = 1.0 / (1.0 + RadiationFieldRedshift);
    aUnits   = 1.0;
  }
  float afloat = float(a);

  /* Metal cooling codes. */
 
  int MetalNum = 0, SNColourNum = 0;
  int MetalFieldPresent = FALSE;

  // First see if there's a metal field (so we can conserve species in
  // the solver)
  MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields);
  SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields);
  MetalFieldPresent = (MetalNum != -1 || SNColourNum != -1);

  // Double check if there's a metal field when we have metal cooling
  if (MetalCooling && MetalFieldPresent == FALSE) {
    if (debug)
      fprintf(stderr, "Warning: No metal field found.  Turning OFF MetalCooling.\n");
    MetalCooling = FALSE;
    MetalNum = 0;
  }

  // If both metal fields (Pop I/II and III) exist, create a field
  // that contains their sum

  float *MetalPointer = NULL;
  float *TotalMetals = NULL;

  if (MetalNum != -1 && SNColourNum != -1) {
    TotalMetals = new float[size];
    for (i = 0; i < size; i++)
      TotalMetals[i] = BaryonField[MetalNum][i] + BaryonField[SNColourNum][i];
    MetalPointer = TotalMetals;
  } // ENDIF both metal types
  else {
    if (MetalNum != -1)
      MetalPointer = BaryonField[MetalNum];
    else if (SNColourNum != -1)
      MetalPointer = BaryonField[SNColourNum];
  } // ENDELSE both metal types

  int temp_thermal = FALSE;
  float *thermal_energy;
  if ( UseMHD ){
    iBx = FindField(Bfield1, FieldType, NumberOfBaryonFields);
    iBy = FindField(Bfield2, FieldType, NumberOfBaryonFields);
    iBz = FindField(Bfield3, FieldType, NumberOfBaryonFields);  
  }

  if (HydroMethod==Zeus_Hydro) {
    thermal_energy = BaryonField[TENum];
  }
  else if (DualEnergyFormalism) {
    thermal_energy = BaryonField[GENum];
  }
  else {
    temp_thermal = TRUE;
    thermal_energy = new float[size];
    for (i = 0; i < size; i++) {
      thermal_energy[i] = BaryonField[TENum][i] - 
        0.5 * POW(BaryonField[Vel1Num][i], 2.0);
      if(GridRank > 1)
        thermal_energy[i] -= 0.5 * POW(BaryonField[Vel2Num][i], 2.0);
      if(GridRank > 2)
        thermal_energy[i] -= 0.5 * POW(BaryonField[Vel3Num][i], 2.0);

      if( UseMHD ) {
        thermal_energy[i] -= 0.5 * (POW(BaryonField[iBx][i], 2.0) + 
                                    POW(BaryonField[iBy][i], 2.0) + 
                                    POW(BaryonField[iBz][i], 2.0)) / 
          BaryonField[DensNum][i];
      }
    } // for (int i = 0; i < size; i++)
  }

  float *temperature = new float[size]; 
  if (this->ComputeTemperatureField(temperature) == FAIL){
    ENZO_FAIL("Error in grid->ComputeTemperatureField.");
  }

  float NumberDensityUnits = DensityUnits / mh;

  Naunet naunet;
  naunet.Init();

  // TODO: comoving, heating/cooling
  
  // Set your parameters here
  NaunetData data;

  float y[NAUNET_NEQNS];

  for (i=0; i<size; i++) {
    data.nH       = BaryonField[iden][i] / (Mu * mh);
    data.Tgas     = temperature[i];
    data.zeta     = 1.3e-17;
    data.Av       = 10.0;
    data.omega    = 0.5;
    data.G0       = 4.0;
    data.rG       = 1.0e-5;
    data.gdens    = 7.6394373e-13 * data.nH;
    data.sites    = 1.5e15;
    data.fr       = 1.0;
    data.opt_thd  = 1.0;
    data.opt_crd  = 1.0;
    data.opt_uvd  = 1.0;
    data.opt_h2d  = 1.0;
    data.crdeseff = 1.0e5;
    data.h2deseff = 1.0e-2;
    data.uvcreff  = 1.0e-3;

    y[IDX_GCH3OHI] =  BaryonField[GCH3OHINum][i] / 32.0 * NumberDensityUnits;
    y[IDX_GCH4I] =  BaryonField[GCH4INum][i] / 16.0 * NumberDensityUnits;
    y[IDX_GCOI] =  BaryonField[GCOINum][i] / 28.0 * NumberDensityUnits;
    y[IDX_GCO2I] =  BaryonField[GCO2INum][i] / 44.0 * NumberDensityUnits;
    y[IDX_GH2CNI] =  BaryonField[GH2CNINum][i] / 28.0 * NumberDensityUnits;
    y[IDX_GH2COI] =  BaryonField[GH2COINum][i] / 30.0 * NumberDensityUnits;
    y[IDX_GH2OI] =  BaryonField[GH2OINum][i] / 18.0 * NumberDensityUnits;
    y[IDX_GH2SiOI] =  BaryonField[GH2SiOINum][i] / 46.0 * NumberDensityUnits;
    y[IDX_GHCNI] =  BaryonField[GHCNINum][i] / 27.0 * NumberDensityUnits;
    y[IDX_GHNCI] =  BaryonField[GHNCINum][i] / 27.0 * NumberDensityUnits;
    y[IDX_GHNCOI] =  BaryonField[GHNCOINum][i] / 43.0 * NumberDensityUnits;
    y[IDX_GHNOI] =  BaryonField[GHNOINum][i] / 31.0 * NumberDensityUnits;
    y[IDX_GMgI] =  BaryonField[GMgINum][i] / 24.0 * NumberDensityUnits;
    y[IDX_GN2I] =  BaryonField[GN2INum][i] / 28.0 * NumberDensityUnits;
    y[IDX_GNH3I] =  BaryonField[GNH3INum][i] / 17.0 * NumberDensityUnits;
    y[IDX_GNOI] =  BaryonField[GNOINum][i] / 30.0 * NumberDensityUnits;
    y[IDX_GNO2I] =  BaryonField[GNO2INum][i] / 46.0 * NumberDensityUnits;
    y[IDX_GO2I] =  BaryonField[GO2INum][i] / 32.0 * NumberDensityUnits;
    y[IDX_GO2HI] =  BaryonField[GO2HINum][i] / 33.0 * NumberDensityUnits;
    y[IDX_GSiCI] =  BaryonField[GSiCINum][i] / 40.0 * NumberDensityUnits;
    y[IDX_GSiC2I] =  BaryonField[GSiC2INum][i] / 52.0 * NumberDensityUnits;
    y[IDX_GSiC3I] =  BaryonField[GSiC3INum][i] / 64.0 * NumberDensityUnits;
    y[IDX_GSiH4I] =  BaryonField[GSiH4INum][i] / 32.0 * NumberDensityUnits;
    y[IDX_GSiOI] =  BaryonField[GSiOINum][i] / 44.0 * NumberDensityUnits;
    y[IDX_CI] =  BaryonField[CINum][i] / 12.0 * NumberDensityUnits;
    y[IDX_CII] =  BaryonField[CIINum][i] / 12.0 * NumberDensityUnits;
    y[IDX_CHI] =  BaryonField[CHINum][i] / 13.0 * NumberDensityUnits;
    y[IDX_CHII] =  BaryonField[CHIINum][i] / 13.0 * NumberDensityUnits;
    y[IDX_CH2I] =  BaryonField[CH2INum][i] / 14.0 * NumberDensityUnits;
    y[IDX_CH2II] =  BaryonField[CH2IINum][i] / 14.0 * NumberDensityUnits;
    y[IDX_CH3I] =  BaryonField[CH3INum][i] / 15.0 * NumberDensityUnits;
    y[IDX_CH3II] =  BaryonField[CH3IINum][i] / 15.0 * NumberDensityUnits;
    y[IDX_CH3OHI] =  BaryonField[CH3OHINum][i] / 32.0 * NumberDensityUnits;
    y[IDX_CH4I] =  BaryonField[CH4INum][i] / 16.0 * NumberDensityUnits;
    y[IDX_CH4II] =  BaryonField[CH4IINum][i] / 16.0 * NumberDensityUnits;
    y[IDX_CNI] =  BaryonField[CNINum][i] / 26.0 * NumberDensityUnits;
    y[IDX_CNII] =  BaryonField[CNIINum][i] / 26.0 * NumberDensityUnits;
    y[IDX_COI] =  BaryonField[COINum][i] / 28.0 * NumberDensityUnits;
    y[IDX_COII] =  BaryonField[COIINum][i] / 28.0 * NumberDensityUnits;
    y[IDX_CO2I] =  BaryonField[CO2INum][i] / 44.0 * NumberDensityUnits;
    y[IDX_EM] =  BaryonField[DeNum][i] / 1.0 * NumberDensityUnits;
    y[IDX_HI] =  BaryonField[HINum][i] / 1.0 * NumberDensityUnits;
    y[IDX_HII] =  BaryonField[HIINum][i] / 1.0 * NumberDensityUnits;
    y[IDX_H2I] =  BaryonField[H2INum][i] / 2.0 * NumberDensityUnits;
    y[IDX_H2II] =  BaryonField[H2IINum][i] / 2.0 * NumberDensityUnits;
    y[IDX_H2CNI] =  BaryonField[H2CNINum][i] / 28.0 * NumberDensityUnits;
    y[IDX_H2COI] =  BaryonField[H2COINum][i] / 30.0 * NumberDensityUnits;
    y[IDX_H2COII] =  BaryonField[H2COIINum][i] / 30.0 * NumberDensityUnits;
    y[IDX_H2NOII] =  BaryonField[H2NOIINum][i] / 32.0 * NumberDensityUnits;
    y[IDX_H2OI] =  BaryonField[H2OINum][i] / 18.0 * NumberDensityUnits;
    y[IDX_H2OII] =  BaryonField[H2OIINum][i] / 18.0 * NumberDensityUnits;
    y[IDX_H2SiOI] =  BaryonField[H2SiOINum][i] / 46.0 * NumberDensityUnits;
    y[IDX_H3II] =  BaryonField[H3IINum][i] / 3.0 * NumberDensityUnits;
    y[IDX_H3COII] =  BaryonField[H3COIINum][i] / 31.0 * NumberDensityUnits;
    y[IDX_H3OII] =  BaryonField[H3OIINum][i] / 19.0 * NumberDensityUnits;
    y[IDX_HCNI] =  BaryonField[HCNINum][i] / 27.0 * NumberDensityUnits;
    y[IDX_HCNII] =  BaryonField[HCNIINum][i] / 27.0 * NumberDensityUnits;
    y[IDX_HCNHII] =  BaryonField[HCNHIINum][i] / 28.0 * NumberDensityUnits;
    y[IDX_HCOI] =  BaryonField[HCOINum][i] / 29.0 * NumberDensityUnits;
    y[IDX_HCOII] =  BaryonField[HCOIINum][i] / 29.0 * NumberDensityUnits;
    y[IDX_HCO2II] =  BaryonField[HCO2IINum][i] / 45.0 * NumberDensityUnits;
    y[IDX_HeI] =  BaryonField[HeINum][i] / 4.0 * NumberDensityUnits;
    y[IDX_HeII] =  BaryonField[HeIINum][i] / 4.0 * NumberDensityUnits;
    y[IDX_HeHII] =  BaryonField[HeHIINum][i] / 5.0 * NumberDensityUnits;
    y[IDX_HNCI] =  BaryonField[HNCINum][i] / 27.0 * NumberDensityUnits;
    y[IDX_HNCOI] =  BaryonField[HNCOINum][i] / 43.0 * NumberDensityUnits;
    y[IDX_HNOI] =  BaryonField[HNOINum][i] / 31.0 * NumberDensityUnits;
    y[IDX_HNOII] =  BaryonField[HNOIINum][i] / 31.0 * NumberDensityUnits;
    y[IDX_HOCII] =  BaryonField[HOCIINum][i] / 29.0 * NumberDensityUnits;
    y[IDX_MgI] =  BaryonField[MgINum][i] / 24.0 * NumberDensityUnits;
    y[IDX_MgII] =  BaryonField[MgIINum][i] / 24.0 * NumberDensityUnits;
    y[IDX_NI] =  BaryonField[NINum][i] / 14.0 * NumberDensityUnits;
    y[IDX_NII] =  BaryonField[NIINum][i] / 14.0 * NumberDensityUnits;
    y[IDX_N2I] =  BaryonField[N2INum][i] / 28.0 * NumberDensityUnits;
    y[IDX_N2II] =  BaryonField[N2IINum][i] / 28.0 * NumberDensityUnits;
    y[IDX_N2HII] =  BaryonField[N2HIINum][i] / 29.0 * NumberDensityUnits;
    y[IDX_NHI] =  BaryonField[NHINum][i] / 15.0 * NumberDensityUnits;
    y[IDX_NHII] =  BaryonField[NHIINum][i] / 15.0 * NumberDensityUnits;
    y[IDX_NH2I] =  BaryonField[NH2INum][i] / 16.0 * NumberDensityUnits;
    y[IDX_NH2II] =  BaryonField[NH2IINum][i] / 16.0 * NumberDensityUnits;
    y[IDX_NH3I] =  BaryonField[NH3INum][i] / 17.0 * NumberDensityUnits;
    y[IDX_NH3II] =  BaryonField[NH3IINum][i] / 17.0 * NumberDensityUnits;
    y[IDX_NOI] =  BaryonField[NOINum][i] / 30.0 * NumberDensityUnits;
    y[IDX_NOII] =  BaryonField[NOIINum][i] / 30.0 * NumberDensityUnits;
    y[IDX_NO2I] =  BaryonField[NO2INum][i] / 46.0 * NumberDensityUnits;
    y[IDX_OI] =  BaryonField[OINum][i] / 16.0 * NumberDensityUnits;
    y[IDX_OII] =  BaryonField[OIINum][i] / 16.0 * NumberDensityUnits;
    y[IDX_O2I] =  BaryonField[O2INum][i] / 32.0 * NumberDensityUnits;
    y[IDX_O2II] =  BaryonField[O2IINum][i] / 32.0 * NumberDensityUnits;
    y[IDX_O2HI] =  BaryonField[O2HINum][i] / 33.0 * NumberDensityUnits;
    y[IDX_O2HII] =  BaryonField[O2HIINum][i] / 33.0 * NumberDensityUnits;
    y[IDX_OCNI] =  BaryonField[OCNINum][i] / 42.0 * NumberDensityUnits;
    y[IDX_OHI] =  BaryonField[OHINum][i] / 17.0 * NumberDensityUnits;
    y[IDX_OHII] =  BaryonField[OHIINum][i] / 17.0 * NumberDensityUnits;
    y[IDX_SiI] =  BaryonField[SiINum][i] / 28.0 * NumberDensityUnits;
    y[IDX_SiII] =  BaryonField[SiIINum][i] / 28.0 * NumberDensityUnits;
    y[IDX_SiCI] =  BaryonField[SiCINum][i] / 40.0 * NumberDensityUnits;
    y[IDX_SiCII] =  BaryonField[SiCIINum][i] / 40.0 * NumberDensityUnits;
    y[IDX_SiC2I] =  BaryonField[SiC2INum][i] / 52.0 * NumberDensityUnits;
    y[IDX_SiC2II] =  BaryonField[SiC2IINum][i] / 52.0 * NumberDensityUnits;
    y[IDX_SiC3I] =  BaryonField[SiC3INum][i] / 64.0 * NumberDensityUnits;
    y[IDX_SiC3II] =  BaryonField[SiC3IINum][i] / 64.0 * NumberDensityUnits;
    y[IDX_SiHI] =  BaryonField[SiHINum][i] / 29.0 * NumberDensityUnits;
    y[IDX_SiHII] =  BaryonField[SiHIINum][i] / 29.0 * NumberDensityUnits;
    y[IDX_SiH2I] =  BaryonField[SiH2INum][i] / 30.0 * NumberDensityUnits;
    y[IDX_SiH2II] =  BaryonField[SiH2IINum][i] / 30.0 * NumberDensityUnits;
    y[IDX_SiH3I] =  BaryonField[SiH3INum][i] / 31.0 * NumberDensityUnits;
    y[IDX_SiH3II] =  BaryonField[SiH3IINum][i] / 31.0 * NumberDensityUnits;
    y[IDX_SiH4I] =  BaryonField[SiH4INum][i] / 32.0 * NumberDensityUnits;
    y[IDX_SiH4II] =  BaryonField[SiH4IINum][i] / 32.0 * NumberDensityUnits;
    y[IDX_SiH5II] =  BaryonField[SiH5IINum][i] / 33.0 * NumberDensityUnits;
    y[IDX_SiOI] =  BaryonField[SiOINum][i] / 44.0 * NumberDensityUnits;
    y[IDX_SiOII] =  BaryonField[SiOIINum][i] / 44.0 * NumberDensityUnits;
    y[IDX_SiOHII] =  BaryonField[SiOHIINum][i] / 45.0 * NumberDensityUnits;
    naunet.Solve(y, dt_chem * TimeUnits, &data);

    BaryonField[GCH3OHINum][i] = y[IDX_GCH3OHI] * 32.0 / NumberDensityUnits;
    BaryonField[GCH4INum][i] = y[IDX_GCH4I] * 16.0 / NumberDensityUnits;
    BaryonField[GCOINum][i] = y[IDX_GCOI] * 28.0 / NumberDensityUnits;
    BaryonField[GCO2INum][i] = y[IDX_GCO2I] * 44.0 / NumberDensityUnits;
    BaryonField[GH2CNINum][i] = y[IDX_GH2CNI] * 28.0 / NumberDensityUnits;
    BaryonField[GH2COINum][i] = y[IDX_GH2COI] * 30.0 / NumberDensityUnits;
    BaryonField[GH2OINum][i] = y[IDX_GH2OI] * 18.0 / NumberDensityUnits;
    BaryonField[GH2SiOINum][i] = y[IDX_GH2SiOI] * 46.0 / NumberDensityUnits;
    BaryonField[GHCNINum][i] = y[IDX_GHCNI] * 27.0 / NumberDensityUnits;
    BaryonField[GHNCINum][i] = y[IDX_GHNCI] * 27.0 / NumberDensityUnits;
    BaryonField[GHNCOINum][i] = y[IDX_GHNCOI] * 43.0 / NumberDensityUnits;
    BaryonField[GHNOINum][i] = y[IDX_GHNOI] * 31.0 / NumberDensityUnits;
    BaryonField[GMgINum][i] = y[IDX_GMgI] * 24.0 / NumberDensityUnits;
    BaryonField[GN2INum][i] = y[IDX_GN2I] * 28.0 / NumberDensityUnits;
    BaryonField[GNH3INum][i] = y[IDX_GNH3I] * 17.0 / NumberDensityUnits;
    BaryonField[GNOINum][i] = y[IDX_GNOI] * 30.0 / NumberDensityUnits;
    BaryonField[GNO2INum][i] = y[IDX_GNO2I] * 46.0 / NumberDensityUnits;
    BaryonField[GO2INum][i] = y[IDX_GO2I] * 32.0 / NumberDensityUnits;
    BaryonField[GO2HINum][i] = y[IDX_GO2HI] * 33.0 / NumberDensityUnits;
    BaryonField[GSiCINum][i] = y[IDX_GSiCI] * 40.0 / NumberDensityUnits;
    BaryonField[GSiC2INum][i] = y[IDX_GSiC2I] * 52.0 / NumberDensityUnits;
    BaryonField[GSiC3INum][i] = y[IDX_GSiC3I] * 64.0 / NumberDensityUnits;
    BaryonField[GSiH4INum][i] = y[IDX_GSiH4I] * 32.0 / NumberDensityUnits;
    BaryonField[GSiOINum][i] = y[IDX_GSiOI] * 44.0 / NumberDensityUnits;
    BaryonField[CINum][i] = y[IDX_CI] * 12.0 / NumberDensityUnits;
    BaryonField[CIINum][i] = y[IDX_CII] * 12.0 / NumberDensityUnits;
    BaryonField[CHINum][i] = y[IDX_CHI] * 13.0 / NumberDensityUnits;
    BaryonField[CHIINum][i] = y[IDX_CHII] * 13.0 / NumberDensityUnits;
    BaryonField[CH2INum][i] = y[IDX_CH2I] * 14.0 / NumberDensityUnits;
    BaryonField[CH2IINum][i] = y[IDX_CH2II] * 14.0 / NumberDensityUnits;
    BaryonField[CH3INum][i] = y[IDX_CH3I] * 15.0 / NumberDensityUnits;
    BaryonField[CH3IINum][i] = y[IDX_CH3II] * 15.0 / NumberDensityUnits;
    BaryonField[CH3OHINum][i] = y[IDX_CH3OHI] * 32.0 / NumberDensityUnits;
    BaryonField[CH4INum][i] = y[IDX_CH4I] * 16.0 / NumberDensityUnits;
    BaryonField[CH4IINum][i] = y[IDX_CH4II] * 16.0 / NumberDensityUnits;
    BaryonField[CNINum][i] = y[IDX_CNI] * 26.0 / NumberDensityUnits;
    BaryonField[CNIINum][i] = y[IDX_CNII] * 26.0 / NumberDensityUnits;
    BaryonField[COINum][i] = y[IDX_COI] * 28.0 / NumberDensityUnits;
    BaryonField[COIINum][i] = y[IDX_COII] * 28.0 / NumberDensityUnits;
    BaryonField[CO2INum][i] = y[IDX_CO2I] * 44.0 / NumberDensityUnits;
    BaryonField[DeNum][i] = y[IDX_EM] * 1.0 / NumberDensityUnits;
    BaryonField[HINum][i] = y[IDX_HI] * 1.0 / NumberDensityUnits;
    BaryonField[HIINum][i] = y[IDX_HII] * 1.0 / NumberDensityUnits;
    BaryonField[H2INum][i] = y[IDX_H2I] * 2.0 / NumberDensityUnits;
    BaryonField[H2IINum][i] = y[IDX_H2II] * 2.0 / NumberDensityUnits;
    BaryonField[H2CNINum][i] = y[IDX_H2CNI] * 28.0 / NumberDensityUnits;
    BaryonField[H2COINum][i] = y[IDX_H2COI] * 30.0 / NumberDensityUnits;
    BaryonField[H2COIINum][i] = y[IDX_H2COII] * 30.0 / NumberDensityUnits;
    BaryonField[H2NOIINum][i] = y[IDX_H2NOII] * 32.0 / NumberDensityUnits;
    BaryonField[H2OINum][i] = y[IDX_H2OI] * 18.0 / NumberDensityUnits;
    BaryonField[H2OIINum][i] = y[IDX_H2OII] * 18.0 / NumberDensityUnits;
    BaryonField[H2SiOINum][i] = y[IDX_H2SiOI] * 46.0 / NumberDensityUnits;
    BaryonField[H3IINum][i] = y[IDX_H3II] * 3.0 / NumberDensityUnits;
    BaryonField[H3COIINum][i] = y[IDX_H3COII] * 31.0 / NumberDensityUnits;
    BaryonField[H3OIINum][i] = y[IDX_H3OII] * 19.0 / NumberDensityUnits;
    BaryonField[HCNINum][i] = y[IDX_HCNI] * 27.0 / NumberDensityUnits;
    BaryonField[HCNIINum][i] = y[IDX_HCNII] * 27.0 / NumberDensityUnits;
    BaryonField[HCNHIINum][i] = y[IDX_HCNHII] * 28.0 / NumberDensityUnits;
    BaryonField[HCOINum][i] = y[IDX_HCOI] * 29.0 / NumberDensityUnits;
    BaryonField[HCOIINum][i] = y[IDX_HCOII] * 29.0 / NumberDensityUnits;
    BaryonField[HCO2IINum][i] = y[IDX_HCO2II] * 45.0 / NumberDensityUnits;
    BaryonField[HeINum][i] = y[IDX_HeI] * 4.0 / NumberDensityUnits;
    BaryonField[HeIINum][i] = y[IDX_HeII] * 4.0 / NumberDensityUnits;
    BaryonField[HeHIINum][i] = y[IDX_HeHII] * 5.0 / NumberDensityUnits;
    BaryonField[HNCINum][i] = y[IDX_HNCI] * 27.0 / NumberDensityUnits;
    BaryonField[HNCOINum][i] = y[IDX_HNCOI] * 43.0 / NumberDensityUnits;
    BaryonField[HNOINum][i] = y[IDX_HNOI] * 31.0 / NumberDensityUnits;
    BaryonField[HNOIINum][i] = y[IDX_HNOII] * 31.0 / NumberDensityUnits;
    BaryonField[HOCIINum][i] = y[IDX_HOCII] * 29.0 / NumberDensityUnits;
    BaryonField[MgINum][i] = y[IDX_MgI] * 24.0 / NumberDensityUnits;
    BaryonField[MgIINum][i] = y[IDX_MgII] * 24.0 / NumberDensityUnits;
    BaryonField[NINum][i] = y[IDX_NI] * 14.0 / NumberDensityUnits;
    BaryonField[NIINum][i] = y[IDX_NII] * 14.0 / NumberDensityUnits;
    BaryonField[N2INum][i] = y[IDX_N2I] * 28.0 / NumberDensityUnits;
    BaryonField[N2IINum][i] = y[IDX_N2II] * 28.0 / NumberDensityUnits;
    BaryonField[N2HIINum][i] = y[IDX_N2HII] * 29.0 / NumberDensityUnits;
    BaryonField[NHINum][i] = y[IDX_NHI] * 15.0 / NumberDensityUnits;
    BaryonField[NHIINum][i] = y[IDX_NHII] * 15.0 / NumberDensityUnits;
    BaryonField[NH2INum][i] = y[IDX_NH2I] * 16.0 / NumberDensityUnits;
    BaryonField[NH2IINum][i] = y[IDX_NH2II] * 16.0 / NumberDensityUnits;
    BaryonField[NH3INum][i] = y[IDX_NH3I] * 17.0 / NumberDensityUnits;
    BaryonField[NH3IINum][i] = y[IDX_NH3II] * 17.0 / NumberDensityUnits;
    BaryonField[NOINum][i] = y[IDX_NOI] * 30.0 / NumberDensityUnits;
    BaryonField[NOIINum][i] = y[IDX_NOII] * 30.0 / NumberDensityUnits;
    BaryonField[NO2INum][i] = y[IDX_NO2I] * 46.0 / NumberDensityUnits;
    BaryonField[OINum][i] = y[IDX_OI] * 16.0 / NumberDensityUnits;
    BaryonField[OIINum][i] = y[IDX_OII] * 16.0 / NumberDensityUnits;
    BaryonField[O2INum][i] = y[IDX_O2I] * 32.0 / NumberDensityUnits;
    BaryonField[O2IINum][i] = y[IDX_O2II] * 32.0 / NumberDensityUnits;
    BaryonField[O2HINum][i] = y[IDX_O2HI] * 33.0 / NumberDensityUnits;
    BaryonField[O2HIINum][i] = y[IDX_O2HII] * 33.0 / NumberDensityUnits;
    BaryonField[OCNINum][i] = y[IDX_OCNI] * 42.0 / NumberDensityUnits;
    BaryonField[OHINum][i] = y[IDX_OHI] * 17.0 / NumberDensityUnits;
    BaryonField[OHIINum][i] = y[IDX_OHII] * 17.0 / NumberDensityUnits;
    BaryonField[SiINum][i] = y[IDX_SiI] * 28.0 / NumberDensityUnits;
    BaryonField[SiIINum][i] = y[IDX_SiII] * 28.0 / NumberDensityUnits;
    BaryonField[SiCINum][i] = y[IDX_SiCI] * 40.0 / NumberDensityUnits;
    BaryonField[SiCIINum][i] = y[IDX_SiCII] * 40.0 / NumberDensityUnits;
    BaryonField[SiC2INum][i] = y[IDX_SiC2I] * 52.0 / NumberDensityUnits;
    BaryonField[SiC2IINum][i] = y[IDX_SiC2II] * 52.0 / NumberDensityUnits;
    BaryonField[SiC3INum][i] = y[IDX_SiC3I] * 64.0 / NumberDensityUnits;
    BaryonField[SiC3IINum][i] = y[IDX_SiC3II] * 64.0 / NumberDensityUnits;
    BaryonField[SiHINum][i] = y[IDX_SiHI] * 29.0 / NumberDensityUnits;
    BaryonField[SiHIINum][i] = y[IDX_SiHII] * 29.0 / NumberDensityUnits;
    BaryonField[SiH2INum][i] = y[IDX_SiH2I] * 30.0 / NumberDensityUnits;
    BaryonField[SiH2IINum][i] = y[IDX_SiH2II] * 30.0 / NumberDensityUnits;
    BaryonField[SiH3INum][i] = y[IDX_SiH3I] * 31.0 / NumberDensityUnits;
    BaryonField[SiH3IINum][i] = y[IDX_SiH3II] * 31.0 / NumberDensityUnits;
    BaryonField[SiH4INum][i] = y[IDX_SiH4I] * 32.0 / NumberDensityUnits;
    BaryonField[SiH4IINum][i] = y[IDX_SiH4II] * 32.0 / NumberDensityUnits;
    BaryonField[SiH5IINum][i] = y[IDX_SiH5II] * 33.0 / NumberDensityUnits;
    BaryonField[SiOINum][i] = y[IDX_SiOI] * 44.0 / NumberDensityUnits;
    BaryonField[SiOIINum][i] = y[IDX_SiOII] * 44.0 / NumberDensityUnits;
    BaryonField[SiOHIINum][i] = y[IDX_SiOHII] * 45.0 / NumberDensityUnits;
    }
  
  
  if (HydroMethod != Zeus_Hydro) {
    for (i = 0; i < size; i++) {
      BaryonField[TENum][i] = thermal_energy[i] +
        0.5 * POW(BaryonField[Vel1Num][i], 2.0);
      if(GridRank > 1)
        BaryonField[TENum][i] += 0.5 * POW(BaryonField[Vel2Num][i], 2.0);
      if(GridRank > 2)
        BaryonField[TENum][i] += 0.5 * POW(BaryonField[Vel3Num][i], 2.0);

      if( UseMHD ) {
        BaryonField[TENum][i] += 0.5 * (POW(BaryonField[iBx][i], 2.0) + 
                                        POW(BaryonField[iBy][i], 2.0) + 
                                        POW(BaryonField[iBz][i], 2.0)) / 
          BaryonField[DensNum][i];
      }

    } // for (int i = 0; i < size; i++)
  } // if (HydroMethod != Zeus_Hydro)

  if (temp_thermal == TRUE) {
    delete [] thermal_energy;
  }
  delete [] temperature;

  delete [] TotalMetals;
  delete [] g_grid_dimension;
  delete [] g_grid_start;
  delete [] g_grid_end;

  LCAPERF_STOP("grid_GrackleWrapper");

  return SUCCESS;
}

#endif