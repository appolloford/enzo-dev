/***********************************************************************
/
/  INITIALIZE TURBULENT CLOUD
/
/  written by: Peng Wang
/  date:       September, 2007
/  modified1:
/
/
************************************************************************/
#ifdef USE_MPI
#include <mpi.h>
#endif /* USE_MPI */

#include <stdio.h>
#include <map>
#include <string>
#include <math.h>
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
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "CommunicationUtilities.h"


void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[], int level);
int GetUnits(float *DensityUnits, float *LengthUnits, float *TemperatureUnits,
             float *TimeUnits, float *VelocityUnits, FLOAT Time);
void MHDCTSetupFieldLabels();

int CollidingCloudInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid, 
        TopGridData &MetaData, int SetBaryonFields)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *BxName = "Bx";
  char *ByName = "By";
  char *BzName = "Bz";
  char *PhiName = "Phi";
  char *ColourName = "colour";
  char *ElectronName = "Electron_Density";
  char *HIName    = "HI_Density";
  char *HIIName   = "HII_Density";
  char *HeIName   = "HeI_Density";
  char *HeIIName  = "HeII_Density";
  char *HeIIIName = "HeIII_Density";
  char *HMName    = "HM_Density";
  char *H2IName   = "H2I_Density";
  char *H2IIName  = "H2II_Density";
  char *DIName    = "DI_Density";
  char *DIIName   = "DII_Density";
  char *HDIName   = "HDI_Density";
  char *kphHIName    = "HI_kph";
  char *gammaHIName  = "HI_gamma";
  char *kphHeIName   = "HeI_kph";
  char *gammaHeIName = "HeI_gamma";
  char *kphHeIIName  = "HeII_kph";
  char *gammaHeIIName= "HeII_gamma";
  char *kdissH2IName = "H2I_kdiss";
  char *RadAccel1Name = "RadAccel1";
  char *RadAccel2Name = "RadAccel2";
  char *RadAccel3Name = "RadAccel3";
  char *Drive1Name = "DrivingField1";
  char *Drive2Name = "DrivingField2";
  char *Drive3Name = "DrivingField3";
  char *GravPotenName = "PotentialField";
  char *Acce1Name = "AccelerationField1";
  char *Acce2Name = "AccelerationField2";
  char *Acce3Name = "AccelerationField3";
  char *MetalName = "Metal_Density";
  char *Phi_pName = "Phip";

  //char *GMC1Name = "GMC_1";
  //char *GMC2Name = "GMC_2";

#ifdef USE_NAUNET
  /* Additional species in deuterium network*/
  const char *GCH3OHIName            = "GCH3OHI_Density";
  const char *GCH4IName              = "GCH4I_Density";
  const char *GCOIName               = "GCOI_Density";
  const char *GCO2IName              = "GCO2I_Density";
  const char *GH2CNIName             = "GH2CNI_Density";
  const char *GH2COIName             = "GH2COI_Density";
  const char *GH2OIName              = "GH2OI_Density";
  const char *GH2SiOIName            = "GH2SiOI_Density";
  const char *GHCNIName              = "GHCNI_Density";
  const char *GHNCIName              = "GHNCI_Density";
  const char *GHNCOIName             = "GHNCOI_Density";
  const char *GHNOIName              = "GHNOI_Density";
  const char *GMgIName               = "GMgI_Density";
  const char *GN2IName               = "GN2I_Density";
  const char *GNH3IName              = "GNH3I_Density";
  const char *GNOIName               = "GNOI_Density";
  const char *GNO2IName              = "GNO2I_Density";
  const char *GO2IName               = "GO2I_Density";
  const char *GO2HIName              = "GO2HI_Density";
  const char *GSiCIName              = "GSiCI_Density";
  const char *GSiC2IName             = "GSiC2I_Density";
  const char *GSiC3IName             = "GSiC3I_Density";
  const char *GSiH4IName             = "GSiH4I_Density";
  const char *GSiOIName              = "GSiOI_Density";
  const char *CIName                 = "CI_Density";
  const char *CIIName                = "CII_Density";
  const char *CHIName                = "CHI_Density";
  const char *CHIIName               = "CHII_Density";
  const char *CH2IName               = "CH2I_Density";
  const char *CH2IIName              = "CH2II_Density";
  const char *CH3IName               = "CH3I_Density";
  const char *CH3IIName              = "CH3II_Density";
  const char *CH3OHIName             = "CH3OHI_Density";
  const char *CH4IName               = "CH4I_Density";
  const char *CH4IIName              = "CH4II_Density";
  const char *CNIName                = "CNI_Density";
  const char *CNIIName               = "CNII_Density";
  const char *COIName                = "COI_Density";
  const char *COIIName               = "COII_Density";
  const char *CO2IName               = "CO2I_Density";
  const char *H2CNIName              = "H2CNI_Density";
  const char *H2COIName              = "H2COI_Density";
  const char *H2COIIName             = "H2COII_Density";
  const char *H2NOIIName             = "H2NOII_Density";
  const char *H2OIName               = "H2OI_Density";
  const char *H2OIIName              = "H2OII_Density";
  const char *H2SiOIName             = "H2SiOI_Density";
  const char *H3IIName               = "H3II_Density";
  const char *H3COIIName             = "H3COII_Density";
  const char *H3OIIName              = "H3OII_Density";
  const char *HCNIName               = "HCNI_Density";
  const char *HCNIIName              = "HCNII_Density";
  const char *HCNHIIName             = "HCNHII_Density";
  const char *HCOIName               = "HCOI_Density";
  const char *HCOIIName              = "HCOII_Density";
  const char *HCO2IIName             = "HCO2II_Density";
  const char *HeHIIName              = "HeHII_Density";
  const char *HNCIName               = "HNCI_Density";
  const char *HNCOIName              = "HNCOI_Density";
  const char *HNOIName               = "HNOI_Density";
  const char *HNOIIName              = "HNOII_Density";
  const char *HOCIIName              = "HOCII_Density";
  const char *MgIName                = "MgI_Density";
  const char *MgIIName               = "MgII_Density";
  const char *NIName                 = "NI_Density";
  const char *NIIName                = "NII_Density";
  const char *N2IName                = "N2I_Density";
  const char *N2IIName               = "N2II_Density";
  const char *N2HIIName              = "N2HII_Density";
  const char *NHIName                = "NHI_Density";
  const char *NHIIName               = "NHII_Density";
  const char *NH2IName               = "NH2I_Density";
  const char *NH2IIName              = "NH2II_Density";
  const char *NH3IName               = "NH3I_Density";
  const char *NH3IIName              = "NH3II_Density";
  const char *NOIName                = "NOI_Density";
  const char *NOIIName               = "NOII_Density";
  const char *NO2IName               = "NO2I_Density";
  const char *OIName                 = "OI_Density";
  const char *OIIName                = "OII_Density";
  const char *O2IName                = "O2I_Density";
  const char *O2IIName               = "O2II_Density";
  const char *O2HIName               = "O2HI_Density";
  const char *O2HIIName              = "O2HII_Density";
  const char *OCNIName               = "OCNI_Density";
  const char *OHIName                = "OHI_Density";
  const char *OHIIName               = "OHII_Density";
  const char *SiIName                = "SiI_Density";
  const char *SiIIName               = "SiII_Density";
  const char *SiCIName               = "SiCI_Density";
  const char *SiCIIName              = "SiCII_Density";
  const char *SiC2IName              = "SiC2I_Density";
  const char *SiC2IIName             = "SiC2II_Density";
  const char *SiC3IName              = "SiC3I_Density";
  const char *SiC3IIName             = "SiC3II_Density";
  const char *SiHIName               = "SiHI_Density";
  const char *SiHIIName              = "SiHII_Density";
  const char *SiH2IName              = "SiH2I_Density";
  const char *SiH2IIName             = "SiH2II_Density";
  const char *SiH3IName              = "SiH3I_Density";
  const char *SiH3IIName             = "SiH3II_Density";
  const char *SiH4IName              = "SiH4I_Density";
  const char *SiH4IIName             = "SiH4II_Density";
  const char *SiH5IIName             = "SiH5II_Density";
  const char *SiOIName               = "SiOI_Density";
  const char *SiOIIName              = "SiOII_Density";
  const char *SiOHIIName             = "SiOHII_Density";
#endif


  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   ret, level, i;

  /* set default parameters */

  int RefineAtStart   = TRUE;
  int PutSink         = FALSE;
  int SetTurbulence = TRUE;
  int RandomSeed = 52761;
  float CloudDensity=1.0, CloudSoundSpeed=1.0, CloudMachNumber=1.0, 
    CloudAngularVelocity = 0.0, InitialBField = 0.0;
  FLOAT CloudRadius = 0.05;
  int CloudType = 1;
  float RelativeVelocity=0.0, Btheta=0.0, ImpactParameter=0.0;  // new params

  /* read input from parameter file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    ret += sscanf(line, "RefineAtStart = %"ISYM, &RefineAtStart);
    ret += sscanf(line, "PutSink = %"ISYM, &PutSink);
    ret += sscanf(line, "Density = %"FSYM, &CloudDensity);
    ret += sscanf(line, "SoundVelocity = %"FSYM, &CloudSoundSpeed);
    ret += sscanf(line, "MachNumber = %"FSYM, &CloudMachNumber);
    ret += sscanf(line, "AngularVelocity = %"FSYM, &CloudAngularVelocity);
    ret += sscanf(line, "CloudRadius = %"PSYM, &CloudRadius);
    ret += sscanf(line, "SetTurbulence = %"ISYM, &SetTurbulence);
    ret += sscanf(line, "RandomSeed = %"ISYM, &RandomSeed);
    ret += sscanf(line, "InitialBfield = %"FSYM, &InitialBField);
    ret += sscanf(line, "CloudType = %"ISYM, &CloudType);
    ret += sscanf(line, "RelativeVelocity = %"FSYM, &RelativeVelocity); // new params
    ret += sscanf(line, "Btheta = %"FSYM, &Btheta);                     // new params
    ret += sscanf(line, "ImpactParameter = %"FSYM, &ImpactParameter);   // new params
  }

  /* Convert to code units */
  
  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1.0, TimeUnits = 1.0, VelocityUnits = 1.0, 
    PressureUnits = 1.0, MagneticUnits = 1.0;
  if (UsePhysicalUnit) 
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, MetaData.Time);
  PressureUnits = DensityUnits * pow(VelocityUnits,2);
  MagneticUnits = sqrt(PressureUnits*4.0*M_PI);

  CloudDensity /= DensityUnits;
  CloudSoundSpeed /= VelocityUnits;
  InitialBField /= MagneticUnits;
  CloudAngularVelocity *= TimeUnits;
  RelativeVelocity /= VelocityUnits;  // new params

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    printf("Magnetic Units=%"GSYM"\n", MagneticUnits);  
    printf("B field=%"GSYM"\n", InitialBField);  
    printf("Plasma beta=%"GSYM"\n", CloudDensity*CloudSoundSpeed*CloudSoundSpeed/(InitialBField*InitialBField/2.0));
    printf("DensityUnits=%"GSYM",VelocityUnits=%"GSYM",LengthUnits=%"GSYM",TimeUnits=%"GSYM" (%"GSYM" yr),PressureUnits=%"GSYM"\n", 
      DensityUnits, VelocityUnits, LengthUnits, TimeUnits, TimeUnits/3.1558e7, PressureUnits);
    printf("CloudDensity=%"GSYM", CloudSoundSpeed=%"GSYM", CloudRadius=%"GSYM", CloudAngularVelocity=%"GSYM"\n", 
      CloudDensity, CloudSoundSpeed, CloudRadius, CloudAngularVelocity);
    printf("Temperature Unit: %10.3e\n", TemperatureUnits);
  }

  /* Begin grid initialization */

  HierarchyEntry *CurrentGrid; // all level 0 grids on this processor first
  CurrentGrid = &TopGrid;
  while (CurrentGrid != NULL) {
    if (CurrentGrid->GridData->CollidingCloudInitializeGrid(
        CloudDensity, CloudSoundSpeed, CloudRadius, CloudMachNumber, CloudAngularVelocity, InitialBField, 
        SetTurbulence, CloudType, RandomSeed, PutSink, 0, SetBaryonFields, RelativeVelocity, Btheta, 
        ImpactParameter) == FAIL) {   // 123456 New Parameters!!
      fprintf(stderr, "Error in CollidingCloudInitializeGrid.\n");
      return FAIL;
    }
    CurrentGrid = CurrentGrid->NextGridThisLevel;
  }

  if (SetBaryonFields) {
    // Compute Velocity Normalization
    double v_rms  = 0;
    double Volume = 0;
    Eflt fac = 1;
    
    if (SetTurbulence) {
      CurrentGrid = &TopGrid;
      while (CurrentGrid != NULL) {
        if (CurrentGrid->GridData->PrepareVelocityNormalization(&v_rms, &Volume) == FAIL) {
          fprintf(stderr, "Error in PrepareVelocityNormalization.\n");
          return FAIL;
        }
        CurrentGrid = CurrentGrid->NextGridThisLevel;
        fprintf(stderr, "v_rms, Volume: %"GSYM"  %"GSYM"\n", v_rms, Volume);
      }
      
#ifdef USE_MPI
      CommunicationAllReduceValues(&v_rms, 1, MPI_SUM);
      CommunicationAllReduceValues(&Volume, 1, MPI_SUM);
#endif
      fprintf(stderr, "v_rms, Volume: %"GSYM"  %"GSYM"\n", v_rms, Volume);
      // Carry out the Normalization
      v_rms = sqrt(v_rms/Volume); // actuall v_rms
      fac = CloudSoundSpeed*CloudMachNumber/v_rms;
      
      CurrentGrid = &TopGrid;
      while (CurrentGrid != NULL) {
        if (CurrentGrid->GridData->NormalizeVelocities(fac) == FAIL) {
          fprintf(stderr, "Error in grid::NormalizeVelocities.\n");
          return FAIL;
        }
        CurrentGrid = CurrentGrid->NextGridThisLevel;
      }
    }


  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set). */

  if (MinimumMassForRefinement[0] == FLOAT_UNDEFINED) {
    MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
    for (int dim = 0; dim < MetaData.TopGridRank; dim++){
      MinimumMassForRefinement[0] *=(DomainRightEdge[dim]-DomainLeftEdge[dim])/
      float(MetaData.TopGridDims[dim]);
    }
  }

  /*
  if (UsePhysicalUnit){
    MinimumMassForRefinement[0] /= DensityUnits * pow(LengthUnits,3);  //#####
  }
  */

  /* If requested, refine the grid to the desired level. */

  if (RefineAtStart) {

    /* Declare, initialize and fill out the LevelArray. */

    LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      LevelArray[level] = NULL;
    AddLevel(LevelArray, &TopGrid, 0);

    /* Add levels to the maximum depth or until no new levels are created,
       and re-initialize the level after it is created. */

    for (level = 0; level < MaximumRefinementLevel; level++) {
      if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
        fprintf(stderr, "Error in RebuildHierarchy.\n");
        return FAIL;
      }
      if (LevelArray[level+1] == NULL)
        break;

      LevelHierarchyEntry *Temp = LevelArray[level+1];
      while (Temp != NULL) {
        if (Temp->GridData->CollidingCloudInitializeGrid(
            CloudDensity, CloudSoundSpeed, CloudRadius, fac, CloudAngularVelocity, InitialBField,
            SetTurbulence, CloudType, RandomSeed, PutSink, level+1, SetBaryonFields,
            RelativeVelocity, Btheta, ImpactParameter) == FAIL) { // 123456 New Parameters !!
          fprintf(stderr, "Error in CollidingCloudInitializeGrid.\n");
          return FAIL;
        }
        Temp = Temp->NextGridThisLevel;
      }
    } // end: loop over levels

    /* Loop back from the bottom, restoring the consistency among levels. */

    for (level = MaximumRefinementLevel; level > 0; level--) {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
        // printf("Project solution back to parent:\n");
        if (Temp->GridData->ProjectSolutionToParentGrid(
            *LevelArray[level-1]->GridData) == FAIL) {
          fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
          return FAIL;
        }
        Temp = Temp->NextGridThisLevel;
      }
    }

  } // end: if (RefineAtStart)

  /* set up field names and units */

  int count = 0;
  DataLabel[count++] = DensName;
  DataLabel[count++] = Vel1Name;
  DataLabel[count++] = Vel2Name;
  DataLabel[count++] = Vel3Name;
  DataLabel[count++] = TEName;
  if (DualEnergyFormalism) {
    DataLabel[count++] = GEName;
  }
  if (UseMHD) {
    DataLabel[count++] = BxName;
    DataLabel[count++] = ByName;
    DataLabel[count++] = BzName;
  }
  if( HydroMethod == MHD_RK){
    DataLabel[count++] = PhiName;
    if(UsePoissonDivergenceCleaning){
      DataLabel[count++] = Phi_pName;
    }
  }
  if (MultiSpecies) {
    DataLabel[count++] = ElectronName;
    DataLabel[count++] = HIName;
    DataLabel[count++] = HIIName;
    DataLabel[count++] = HeIName;
    DataLabel[count++] = HeIIName;
    DataLabel[count++] = HeIIIName;
    if (MultiSpecies > 1) {
      DataLabel[count++] = HMName;
      DataLabel[count++] = H2IName;
      DataLabel[count++] = H2IIName;
    }
    if (MultiSpecies > 2) {
      DataLabel[count++] = DIName;
      DataLabel[count++] = DIIName;
      DataLabel[count++] = HDIName;
    }
#ifdef USE_NAUNET
    if (MultiSpecies == NAUNET_SPECIES) {
      DataLabel[count++] = (char*) GCH3OHIName;
      DataLabel[count++] = (char*) GCH4IName;
      DataLabel[count++] = (char*) GCOIName;
      DataLabel[count++] = (char*) GCO2IName;
      DataLabel[count++] = (char*) GH2CNIName;
      DataLabel[count++] = (char*) GH2COIName;
      DataLabel[count++] = (char*) GH2OIName;
      DataLabel[count++] = (char*) GH2SiOIName;
      DataLabel[count++] = (char*) GHCNIName;
      DataLabel[count++] = (char*) GHNCIName;
      DataLabel[count++] = (char*) GHNCOIName;
      DataLabel[count++] = (char*) GHNOIName;
      DataLabel[count++] = (char*) GMgIName;
      DataLabel[count++] = (char*) GN2IName;
      DataLabel[count++] = (char*) GNH3IName;
      DataLabel[count++] = (char*) GNOIName;
      DataLabel[count++] = (char*) GNO2IName;
      DataLabel[count++] = (char*) GO2IName;
      DataLabel[count++] = (char*) GO2HIName;
      DataLabel[count++] = (char*) GSiCIName;
      DataLabel[count++] = (char*) GSiC2IName;
      DataLabel[count++] = (char*) GSiC3IName;
      DataLabel[count++] = (char*) GSiH4IName;
      DataLabel[count++] = (char*) GSiOIName;
      DataLabel[count++] = (char*) CIName;
      DataLabel[count++] = (char*) CIIName;
      DataLabel[count++] = (char*) CHIName;
      DataLabel[count++] = (char*) CHIIName;
      DataLabel[count++] = (char*) CH2IName;
      DataLabel[count++] = (char*) CH2IIName;
      DataLabel[count++] = (char*) CH3IName;
      DataLabel[count++] = (char*) CH3IIName;
      DataLabel[count++] = (char*) CH3OHIName;
      DataLabel[count++] = (char*) CH4IName;
      DataLabel[count++] = (char*) CH4IIName;
      DataLabel[count++] = (char*) CNIName;
      DataLabel[count++] = (char*) CNIIName;
      DataLabel[count++] = (char*) COIName;
      DataLabel[count++] = (char*) COIIName;
      DataLabel[count++] = (char*) CO2IName;
      DataLabel[count++] = (char*) H2CNIName;
      DataLabel[count++] = (char*) H2COIName;
      DataLabel[count++] = (char*) H2COIIName;
      DataLabel[count++] = (char*) H2NOIIName;
      DataLabel[count++] = (char*) H2OIName;
      DataLabel[count++] = (char*) H2OIIName;
      DataLabel[count++] = (char*) H2SiOIName;
      DataLabel[count++] = (char*) H3IIName;
      DataLabel[count++] = (char*) H3COIIName;
      DataLabel[count++] = (char*) H3OIIName;
      DataLabel[count++] = (char*) HCNIName;
      DataLabel[count++] = (char*) HCNIIName;
      DataLabel[count++] = (char*) HCNHIIName;
      DataLabel[count++] = (char*) HCOIName;
      DataLabel[count++] = (char*) HCOIIName;
      DataLabel[count++] = (char*) HCO2IIName;
      DataLabel[count++] = (char*) HeHIIName;
      DataLabel[count++] = (char*) HNCIName;
      DataLabel[count++] = (char*) HNCOIName;
      DataLabel[count++] = (char*) HNOIName;
      DataLabel[count++] = (char*) HNOIIName;
      DataLabel[count++] = (char*) HOCIIName;
      DataLabel[count++] = (char*) MgIName;
      DataLabel[count++] = (char*) MgIIName;
      DataLabel[count++] = (char*) NIName;
      DataLabel[count++] = (char*) NIIName;
      DataLabel[count++] = (char*) N2IName;
      DataLabel[count++] = (char*) N2IIName;
      DataLabel[count++] = (char*) N2HIIName;
      DataLabel[count++] = (char*) NHIName;
      DataLabel[count++] = (char*) NHIIName;
      DataLabel[count++] = (char*) NH2IName;
      DataLabel[count++] = (char*) NH2IIName;
      DataLabel[count++] = (char*) NH3IName;
      DataLabel[count++] = (char*) NH3IIName;
      DataLabel[count++] = (char*) NOIName;
      DataLabel[count++] = (char*) NOIIName;
      DataLabel[count++] = (char*) NO2IName;
      DataLabel[count++] = (char*) OIName;
      DataLabel[count++] = (char*) OIIName;
      DataLabel[count++] = (char*) O2IName;
      DataLabel[count++] = (char*) O2IIName;
      DataLabel[count++] = (char*) O2HIName;
      DataLabel[count++] = (char*) O2HIIName;
      DataLabel[count++] = (char*) OCNIName;
      DataLabel[count++] = (char*) OHIName;
      DataLabel[count++] = (char*) OHIIName;
      DataLabel[count++] = (char*) SiIName;
      DataLabel[count++] = (char*) SiIIName;
      DataLabel[count++] = (char*) SiCIName;
      DataLabel[count++] = (char*) SiCIIName;
      DataLabel[count++] = (char*) SiC2IName;
      DataLabel[count++] = (char*) SiC2IIName;
      DataLabel[count++] = (char*) SiC3IName;
      DataLabel[count++] = (char*) SiC3IIName;
      DataLabel[count++] = (char*) SiHIName;
      DataLabel[count++] = (char*) SiHIIName;
      DataLabel[count++] = (char*) SiH2IName;
      DataLabel[count++] = (char*) SiH2IIName;
      DataLabel[count++] = (char*) SiH3IName;
      DataLabel[count++] = (char*) SiH3IIName;
      DataLabel[count++] = (char*) SiH4IName;
      DataLabel[count++] = (char*) SiH4IIName;
      DataLabel[count++] = (char*) SiH5IIName;
      DataLabel[count++] = (char*) SiOIName;
      DataLabel[count++] = (char*) SiOIIName;
      DataLabel[count++] = (char*) SiOHIIName;
      }
#endif

  }  // if Multispecies                                                                                                   
  //  DataLabel[count++] = MetalName;
  //if (PhotonTestUseColour)
  //DataLabel[count++] = ColourName;
#ifdef TRANSFER
  if (RadiativeTransfer)
    if (MultiSpecies) {
      DataLabel[count++]  = kphHIName;
      DataLabel[count++]  = gammaHIName;
      DataLabel[count++]  = kphHeIName;
      DataLabel[count++]  = gammaHeIName;
      DataLabel[count++]  = kphHeIIName;
      DataLabel[count++]  = gammaHeIIName;
      if (MultiSpecies > 1)
        DataLabel[count++]= kdissH2IName;
    } // if RadiativeTransfer                                                                                                

  if (RadiationPressure) {
    DataLabel[count++]  = RadAccel1Name;
    DataLabel[count++]  = RadAccel2Name;
    DataLabel[count++]  = RadAccel3Name;
  }
#endif
  if (UseDrivingField) {
    DataLabel[count++] = Drive1Name;
    DataLabel[count++] = Drive2Name;
    DataLabel[count++] = Drive3Name;
  }
  if (WritePotential) {
    DataLabel[count++] = GravPotenName;
    DataLabel[count++] = Acce1Name;
    DataLabel[count++] = Acce2Name;
    DataLabel[count++] = Acce3Name;
  }
  MHDCTSetupFieldLabels();

  //if (ProblemType == 106) {
  //  DataLabel[count++] = GMC1Name;
  //  DataLabel[count++] = GMC2Name;
  //}

  for (i = 0; i < count; i++) {
    DataUnits[i] = NULL;
  }

  /* If streaming movie output, write header file. */

//   FILE *header;
//   char *headerName = (char*) "movieHeader.dat";
//   int sizeOfFLOAT = sizeof(FLOAT);
//   int sizeOfRecord = (7+MAX_MOVIE_FIELDS)*sizeof(int) + sizeof(float) +
//     6*sizeof(FLOAT);
//   char *movieVersion = (char*) "1.3";
//   int nMovieFields = 0;
//   while (MovieDataField[nMovieFields] != INT_UNDEFINED &&
//          nMovieFields < MAX_MOVIE_FIELDS) nMovieFields++;

//   if (MovieSkipTimestep != INT_UNDEFINED) {
//     if ((header = fopen(headerName, "w")) == NULL) {
//       fprintf(stderr, "Error in opening movie header.\n");
//       return FAIL;
//     }
//     fprintf(header, "MovieVersion = %s\n", movieVersion);
//     fprintf(header, "RootReso = %"ISYM"\n", MetaData.TopGridDims[0]);
//     fprintf(header, "FLOATSize = %"ISYM"\n", sizeOfFLOAT);
//     fprintf(header, "RecordSize = %"ISYM"\n", sizeOfRecord);
//     fprintf(header, "NumFields = %"ISYM"\n", nMovieFields);
//     fprintf(header, "NumCPUs = %"ISYM"\n", NumberOfProcessors);
//     fprintf(header, "FileStem = %s\n", NewMovieName);
//     fclose(header);
//   } /* END: write movie header file */
//   /* Open Amira Data file, if requested */

//   if (MovieSkipTimestep != INT_UNDEFINED) {
//     char *AmiraFileName = new char[80];
//     char *TextureFileName = new char[80];
//     int *RefineByArray = new int[3];
//     bool error = FALSE;
//     //    hid_t DataType = H5T_NATIVE_FLOAT;
//     char fileID[4], pid[6];
//     float root_dx = 1.0 / MetaData.TopGridDims[0];

//     //    staggering stag = CELL_CENTERED;
//     //    fieldtype field_type = SCALAR;
//     for (dim = 0; dim < MAX_DIMENSION; dim++) RefineByArray[dim] = RefineBy;

//     sprintf(pid, "_P%3.3d", MyProcessorNumber);
//     sprintf(fileID, "%3.3d", NewMovieDumpNumber);

//     strcpy(AmiraFileName, "AmiraData");
//     strcat(AmiraFileName, fileID);
//     strcat(AmiraFileName, pid);
//     strcat(AmiraFileName, ".hdf5");

//     strcpy(TextureFileName, "TextureData");
//     strcat(TextureFileName, fileID);
//     strcat(TextureFileName, pid);
//     strcat(TextureFileName, ".hdf5");

//     int field, nBaryonFields;
//     int nFields = 0;
//     while (MovieDataField[nFields] != INT_UNDEFINED)
//       nFields++;
//     nBaryonFields = TopGrid.GridData->ReturnNumberOfBaryonFields();

//     char **FieldNames = new char*[nFields];
//     for (field = 0; field < nFields; field++) {
//       FieldNames[field] = new char[64];

//       if (MovieDataField[field] != TEMPERATURE_FIELD)
//         strcpy(FieldNames[field], DataLabel[MovieDataField[field]]);
//       else
//         strcpy(FieldNames[field], "Temperature");
//     }
//     if (Movie3DVolumes > 0)
//       MetaData.AmiraGrid.AMRHDF5Create(AmiraFileName, RefineByArray,
//                                        DataType, stag, field_type,
//                                        MetaData.CycleNumber, MetaData.Time,
//                                        0, root_dx, 1,
//                                        (MyProcessorNumber == ROOT_PROCESSOR),
//                                        nFields, (NewMovieParticleOn > 0),
//                                        NumberOfParticleAttributes, FieldNames,
//                                        error);

//     if (Movie2DTextures > 0)
//       MetaData.Textures.AMRHDF5Create(TextureFileName, RefineByArray,
//                                       DataType, stag, field_type,
//                                       MetaData.CycleNumber, MetaData.Time,
//                                       0, root_dx, 1,
//                                       (MyProcessorNumber == ROOT_PROCESSOR),
//                                       nFields, (NewMovieParticleOn > 0),
//                                       NumberOfParticleAttributes, FieldNames,
//                                       error);

//     if (error) {
//       fprintf(stderr, "Error in AMRHDF5Writer.\n");
//       return FAIL;
//     }

//     delete [] AmiraFileName;
//     delete [] TextureFileName;
//     delete [] RefineByArray;
//     for (field = 0; field < nFields; field++)
//       delete [] FieldNames[field];

//   } /* ENDIF Movie */

  } // endif SetBaryonFields

  return SUCCESS;

}
