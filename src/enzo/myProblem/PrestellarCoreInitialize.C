/***********************************************************************
/
/  INITIALIZE A PRESTELLAR CORE
/
/  written by: Chia-Jung Hsu
/  date:       March, 2019
/  modified1:  
/
/  PURPOSE:
/
/  REFERENCE: Goodson, M. D., Kong, S., Tan, J. C., Heitsch, F., 
/             & Caselli, P. (2016). 
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
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
#include "phys_constants.h"

void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,LevelHierarchyEntry *LevelArray[], 
        int level);
int GetUnits(float *DensityUnits, float *LengthUnits, float *TemperatureUnits,
        float *TimeUnits, float *VelocityUnits, FLOAT Time);
void TurbulenceGenerator_FFTW(const int field_type, const int gridsize, const int randomSeed,
        const int KStart, const int KEnd, const float power_turb, float* vfield);

int PrestellarCoreInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
        TopGridData &MetaData, bool SetBaryonField)
{
  char  *DensName  = "Density";
  char  *TEName    = "TotalEnergy";
  char  *GEName    = "GasEnergy";
  char  *Vel1Name  = "x-velocity";
  char  *Vel2Name  = "y-velocity";
  char  *Vel3Name  = "z-velocity";
  char  *B1Name    = "Bx";
  char  *B2Name    = "By";
  char  *B3Name    = "Bz";
  char  *PhiName   = "Phi";
  char  *GPotName  = "Grav_Potential";

  const char* DeName           = "Electron_Density";
  const char* CHIName          = "CHI_Density";
  const char* OIName           = "OI_Density";
  const char* HNCIName         = "HNCI_Density";
  const char* HCNIName         = "HCNI_Density";
  const char* H2IName          = "H2I_Density";
  const char* CIName           = "CI_Density";
  const char* HIName           = "HI_Density";
  const char* H2OIName         = "H2OI_Density";
  const char* OHIName          = "OHI_Density";
  const char* O2IName          = "O2I_Density";
  const char* CH2IName         = "CH2I_Density";
  const char* H2COIName        = "H2COI_Density";
  const char* HCOIName         = "HCOI_Density";
  const char* MGIName          = "MGI_Density";
  const char* NH3IName         = "NH3I_Density";
  const char* NOIName          = "NOI_Density";
  const char* CNIName          = "CNI_Density";
  const char* COIName          = "COI_Density";
  const char* N2IName          = "N2I_Density";
  const char* NH2IName         = "NH2I_Density";
  const char* CH3IName         = "CH3I_Density";
  const char* CH4IName         = "CH4I_Density";
  const char* NIName           = "NI_Density";
  const char* NHIName          = "NHI_Density";
  const char* HeIName          = "HeI_Density";
  const char* HNOIName         = "HNOI_Density";
  const char* CH3OHIName       = "CH3OHI_Density";
  const char* CO2IName         = "CO2I_Density";
  const char* H2CNIName        = "H2CNI_Density";
  const char* HNCOIName        = "HNCOI_Density";
  const char* NO2IName         = "NO2I_Density";
  const char* O2HIName         = "O2HI_Density";
  const char* OCNIName         = "OCNI_Density";
  const char* CH3OH_DUSTIName  = "CH3OH_DUSTI_Density";
  const char* HNCO_DUSTIName   = "HNCO_DUSTI_Density";
  const char* H2CO_DUSTIName   = "H2CO_DUSTI_Density";
  const char* CH4_DUSTIName    = "CH4_DUSTI_Density";
  const char* CO_DUSTIName     = "CO_DUSTI_Density";
  const char* H2O_DUSTIName    = "H2O_DUSTI_Density";
  const char* NO_DUSTIName     = "NO_DUSTI_Density";
  const char* CO2_DUSTIName    = "CO2_DUSTI_Density";
  const char* N2_DUSTIName     = "N2_DUSTI_Density";
  const char* HCN_DUSTIName    = "HCN_DUSTI_Density";
  const char* NH3_DUSTIName    = "NH3_DUSTI_Density";
  const char* O2_DUSTIName     = "O2_DUSTI_Density";
  const char* NO2_DUSTIName    = "NO2_DUSTI_Density";
  const char* HNO_DUSTIName    = "HNO_DUSTI_Density";
  const char* O2H_DUSTIName    = "O2H_DUSTI_Density";
  const char* H2CN_DUSTIName   = "H2CN_DUSTI_Density";
  const char* MG_DUSTIName     = "MG_DUSTI_Density";
  const char* HNC_DUSTIName    = "HNC_DUSTI_Density";
  const char* E_DUSTIName      = "E_DUSTI_Density";
  const char* HCOIIName        = "HCOII_Density";
  const char* HIIName          = "HII_Density";
  const char* HOCIIName        = "HOCII_Density";
  const char* CIIName          = "CII_Density";
  const char* CH2IIName        = "CH2II_Density";
  const char* CHIIName         = "CHII_Density";
  const char* H2COIIName       = "H2COII_Density";
  const char* MGIIName         = "MGII_Density";
  const char* NH3IIName        = "NH3II_Density";
  const char* NOIIName         = "NOII_Density";
  const char* CNIIName         = "CNII_Density";
  const char* COIIName         = "COII_Density";
  const char* N2IIName         = "N2II_Density";
  const char* O2IIName         = "O2II_Density";
  const char* H2OIIName        = "H2OII_Density";
  const char* NH2IIName        = "NH2II_Density";
  const char* OIIName          = "OII_Density";
  const char* OHIIName         = "OHII_Density";
  const char* CH3IIName        = "CH3II_Density";
  const char* CH4IIName        = "CH4II_Density";
  const char* NIIName          = "NII_Density";
  const char* HCNIIName        = "HCNII_Density";
  const char* NHIIName         = "NHII_Density";
  const char* H2IIName         = "H2II_Density";
  const char* HeIIName         = "HeII_Density";
  const char* HNOIIName        = "HNOII_Density";
  const char* H2NOIIName       = "H2NOII_Density";
  const char* H3IIName         = "H3II_Density";
  const char* H3COIIName       = "H3COII_Density";
  const char* H3OIIName        = "H3OII_Density";
  const char* HCNHIIName       = "HCNHII_Density";
  const char* HCO2IIName       = "HCO2II_Density";
  const char* HeHIIName        = "HeHII_Density";
  const char* N2HIIName        = "N2HII_Density";
  const char* O2HIIName        = "O2HII_Density";

  const char* HeIIIName        = "HeIII_Density";
  const char* HMName           = "HM_Density";
  const char* DIName           = "DI_Density";
  const char* DIIName          = "DII_Density";
  const char* HDIName          = "HDI_Density";

  /* parameter declarations */

  FLOAT PrestellarCoreSubgridLeft, PrestellarCoreSubgridRight;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
 
  /* local declarations */

  char line[MAX_LINE_LENGTH];
  int  dim, ret, NumberOfSubgridZones[MAX_DIMENSION],
       SubgridDims[MAX_DIMENSION];

  float PrestellarCoreDensity        = 10.0,
        PrestellarCoreSurfaceDensity = 1.0,
        PrestellarCoreDensityJump    = 0.1,
        PrestellarCoreTotalEnergy    = 1.0,
        PrestellarCoreGasEnergy      = 0.5,
        PrestellarCoreTemperature    = 10.0,
        PrestellarCoreAngularVelocity= 0.0,
        PrestellarCoreRadius         = 1.0,
        PrestellarCoreVelocity[3]    = {0.0, 0.0, 0.0},
        PrestellarCoreBField[3]      = {0.0, 0.0, 0.0},
        PrestellarCoreBzField        = 0.0,
        PrestellarCoreAmbientBzField = 0.0,
        PrestellarCoreVelocityDispersion = 0.0,
        PrestellarCoreTurbulenceKStart   = 1.0,
        PrestellarCoreTurbulenceKEnd     = -1.0,
        PrestellarCoreTurbulencePower    = -4.0,
        PrestellarCoreOPR = 3.0,
        PrestellarCoreCODeplete = 10.0;

  int PrestellarCoreRandomSeed = 65536,
      PrestellarCoreTurbulenceType = 2; // curl_free = 1, div_free = 2

  // Add one more index for electron
  float *PrestellarCoreInitAbundance = new float [NSpecies+1],
        *PrestellarCoreMoleMass = new float [NSpecies+1];

  const float Bzmconst = 1.38e-16,
              gasConst = 8.3144725e7;

  /* set no subgrids by default. */
  PrestellarCoreSubgridLeft         = 0.0;    // start of subgrid(s)
  PrestellarCoreSubgridRight        = 0.0;    // end of subgrid(s)

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */
    ret += sscanf(line, "PrestellarCoreRadius            = %"FSYM, 
          &PrestellarCoreRadius);
    ret += sscanf(line, "PrestellarCoreDensity           = %"FSYM, 
          &PrestellarCoreDensity);
    ret += sscanf(line, "PrestellarCoreSurfaceDensity    = %"FSYM, 
          &PrestellarCoreSurfaceDensity);
    ret += sscanf(line, "PrestellarCoreDensityJump       = %"FSYM, 
          &PrestellarCoreDensityJump);
    ret += sscanf(line, "PrestellarCoreTemperature       = %"FSYM, 
          &PrestellarCoreTemperature);
    ret += sscanf(line, "PrestellarCoreAngularVelocity   = %"FSYM, 
          &PrestellarCoreAngularVelocity);
    ret += sscanf(line, "PrestellarCoreBzField           = %"FSYM, 
          &PrestellarCoreBzField);
    ret += sscanf(line, "PrestellarCoreAmbientBzField    = %"FSYM, 
          &PrestellarCoreAmbientBzField);
    ret += sscanf(line, "PrestellarCoreVelocityDispersion= %"FSYM, 
          &PrestellarCoreVelocityDispersion);
    ret += sscanf(line, "PrestellarCoreSubgridLeft       = %"FSYM, 
          &PrestellarCoreSubgridLeft);
    ret += sscanf(line, "PrestellarCoreSubgridRight      = %"FSYM, 
          &PrestellarCoreSubgridRight);
    ret += sscanf(line, "PrestellarCoreTurbulenceKStart  = %"FSYM, 
          &PrestellarCoreTurbulenceKStart);
    ret += sscanf(line, "PrestellarCoreTurbulenceKEnd    = %"FSYM, 
          &PrestellarCoreTurbulenceKEnd);
    ret += sscanf(line, "PrestellarCoreTurbulencePower   = %"FSYM, 
          &PrestellarCoreTurbulencePower);
    ret += sscanf(line, "PrestellarCoreTurbulenceType    = %"ISYM, 
          &PrestellarCoreTurbulenceType);
    ret += sscanf(line, "PrestellarCoreOPR               = %"FSYM, 
          &PrestellarCoreOPR);
    ret += sscanf(line, "PrestellarCoreCODeplete         = %"FSYM, 
          &PrestellarCoreCODeplete);
    ret += sscanf(line, "PrestellarCoreRandomSeed        = %"ISYM, 
          &PrestellarCoreRandomSeed);

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "PrestellarCore") && 
        line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, 
          "warning: the following parameter line was not interpreted:\n%s\n", line);
  } // end input from parameter file

  for (int i=0; i<NSpecies+1; i++){
    PrestellarCoreInitAbundance[i] = 0.0;
    PrestellarCoreMoleMass[i] = 0.0;
  }
  float PrestellarCoreInternalEnergy;
  PrestellarCoreInternalEnergy = PrestellarCoreTemperature * gasConst
                            / (Gamma - 1.0) / Mu;

  /* get units */

  float DensityUnits     = 1.0, 
        LengthUnits      = 1.0, 
        TemperatureUnits = 1.0, 
        TimeUnits        = 1.0, 
        VelocityUnits    = 1.0;
  if (UsePhysicalUnit)
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, 
        &VelocityUnits, 0.0);

  PrestellarCoreRadius /= LengthUnits;
  PrestellarCoreDensity /= DensityUnits;
  PrestellarCoreSurfaceDensity /= DensityUnits;
  PrestellarCoreAngularVelocity *= TimeUnits;
  //PrestellarCoreTemperature = PrestellarCoreTemperature * Bzmconst * gamma
  //                          * pow(TimeUnits,2)/pow(LengthUnits,2)
  //                          / DensityUnits;
  PrestellarCoreInternalEnergy = PrestellarCoreInternalEnergy * pow(TimeUnits,2)/pow(LengthUnits, 2);

  float MassUnits = DensityUnits*pow(LengthUnits,3);
  GravitationalConstant = 4.0*M_PI*GravConst*MassUnits*pow(TimeUnits,2)/pow(LengthUnits,3);

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    printf ("Mass Unit: %13.7e\n", MassUnits);
    printf ("G: %13.7e\n", GravitationalConstant);
    printf ("Internal energy: %13.7e\n", PrestellarCoreInternalEnergy);
    printf ("Energy Unit %13.7e\n", pow(LengthUnits, 2)/pow(TimeUnits,2));
    printf ("velocity Unit %13.7e\n", VelocityUnits);
  }

  //float MagnetUnits = pow(MassUnits, 0.5) / pow(LengthUnits, 0.5) / TimeUnits;
  float MagnetUnits = sqrt(4.0*M_PI*DensityUnits)*VelocityUnits;
  PrestellarCoreBzField /= MagnetUnits;
  PrestellarCoreAmbientBzField /= MagnetUnits;

  PrestellarCoreVelocityDispersion /= VelocityUnits;

  /* Set up current problem time, ambient total energy. */

  MetaData.Time         = 0.0;

  /* set periodic boundaries,
     otherwise keep reflecting (the default) */

  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    MetaData.LeftFaceBoundaryCondition[dim]  = periodic;
    MetaData.RightFaceBoundaryCondition[dim] = periodic;
  }

  /* set up uniform grid without the core */

  int baseDims[3];
  for (dim = 0; dim < MetaData.TopGridRank; dim++){
    baseDims[dim] = MetaData.TopGridDims[dim];
  }

  float r_f = 0.15 * PrestellarCoreRadius,
        r_s = 0.05 * PrestellarCoreRadius,
        r_c = PrestellarCoreRadius;
  float xcenter = (DomainRightEdge[0] + DomainLeftEdge[0])/2.0;
  float ycenter = (DomainRightEdge[1] + DomainLeftEdge[1])/2.0;
  float zcenter = (DomainRightEdge[2] + DomainLeftEdge[2])/2.0;
  float PrestellarCoreOuterDensity = PrestellarCoreDensityJump * PrestellarCoreSurfaceDensity;

  bool initTurb = PrestellarCoreTurbulenceKStart < PrestellarCoreTurbulenceKEnd
                  && SetBaryonField;
  const int KStart = (int)PrestellarCoreTurbulenceKStart;
  const int KEnd   = (int)PrestellarCoreTurbulenceKEnd;
  const float Power = PrestellarCoreTurbulencePower;
  float* Turbulence;
  int GridRank = MetaData.TopGridRank;
  if (initTurb) {
    // initialize turbulence on root grid level
    int turbsize=1;
    for (dim = 0; dim < GridRank; dim++) {
      turbsize *= baseDims[dim];
    }

    int gridsize = baseDims[0];
    Turbulence = new float[GridRank*turbsize];
    if (MyProcessorNumber == ROOT_PROCESSOR)
      printf("turbsize: %d\n", turbsize);

    TurbulenceGenerator_FFTW( PrestellarCoreTurbulenceType, gridsize, PrestellarCoreRandomSeed,
                              KStart, KEnd, Power, Turbulence);
 
    int m, n, l;
    float x, y, z, rad, density, mass, coreMass;
    float cellWidth[3];
    for (int i=0; i<3; i++) 
      cellWidth[i] = (DomainRightEdge[i] - DomainLeftEdge[i]) / baseDims[i];
    // calculate the standard deviation used to normalize according to velocity dispersion
    float mean[3] = {0.0}, meanNW[3] = {0.0}, meanC[3] = {0.0}, stdDev[3] = {0.0};
    for (dim = 0; dim < GridRank; dim++){
      mass = 0.0;
      coreMass = 0.0;
      for (int i = 0; i < turbsize; i++){
        m = i % baseDims[0];
        n = ((i-m) % (baseDims[0]*baseDims[1]))/baseDims[0];
        l = (i-m-n*baseDims[0])/baseDims[0]/baseDims[1];
        x = DomainLeftEdge[0] + (m + 0.5)*(cellWidth[0]) - xcenter;
        y = DomainLeftEdge[1] + (n + 0.5)*(cellWidth[1]) - ycenter;
        z = DomainLeftEdge[2] + (l + 0.5)*(cellWidth[2]) - zcenter;
        rad = sqrt(x*x + y*y + z*z);
        density = PrestellarCoreOuterDensity
            + (PrestellarCoreDensity - PrestellarCoreOuterDensity)
            / (1.0 + pow( (rad/r_f), 1.5))
            * (0.5 - 0.5*tanh( (rad-r_c)/ r_s ));
        mass += density;

        mean[dim] += Turbulence[dim + i*GridRank]*density;
        meanNW[dim] += Turbulence[dim + i*GridRank];
        if ( rad <= PrestellarCoreRadius ) {
          meanC[dim] += Turbulence[dim + i*GridRank]*density;
          coreMass += density;
        }
      }
      mean[dim] /= mass;
      mean[dim] /= turbsize;
      meanNW[dim] /= turbsize;
      meanC[dim] /= coreMass;
      meanC[dim] /= turbsize;
      if (MyProcessorNumber == ROOT_PROCESSOR) {
        printf("vel mean value of dim %d: %13.7e \n", dim, mean[dim]);
        printf("no weighted vel mean value of dim %d: %13.7e \n", dim, meanNW[dim]);
        printf("core vel mean value of dim %d: %13.7e \n", dim, meanC[dim]);
      }
      for (int i = 0; i < turbsize; i++){
        stdDev[dim] += pow((Turbulence[dim + i*GridRank]-mean[dim]), 2.0);
      }
      stdDev[dim] = sqrt(stdDev[dim]/turbsize);
      if (MyProcessorNumber == ROOT_PROCESSOR)
        printf("vel standard deviation of dim %d: %13.7e \n", dim, stdDev[dim]);
    }
    for (dim = 0; dim < GridRank; dim++)
      for (int i=0; i<turbsize; i++) {
        Turbulence[dim + i*GridRank] -= meanC[dim];
        Turbulence[dim + i*GridRank] *= PrestellarCoreVelocityDispersion/stdDev[dim];
      }

  }

  HierarchyEntry *CurrentGrid; // all level 0 grids on this processor first
  CurrentGrid = &TopGrid;
  //int count = 0;
  while (CurrentGrid != NULL) {
    //printf("count %i %i\n", count++, MyProcessorNumber);
    if (CurrentGrid->GridData->PrestellarCoreInitializeGrid(
                                  PrestellarCoreRadius,
                                  PrestellarCoreDensity,
                                  PrestellarCoreSurfaceDensity,
                                  PrestellarCoreDensityJump,
                                  PrestellarCoreInternalEnergy,
                                  PrestellarCoreAngularVelocity,
                                  PrestellarCoreBzField,
                                  PrestellarCoreAmbientBzField,
                                  PrestellarCoreVelocityDispersion,
                                  PrestellarCoreTurbulenceKStart,
                                  PrestellarCoreTurbulenceKEnd,
                                  PrestellarCoreOPR,
                                  PrestellarCoreCODeplete,
                                  PrestellarCoreRandomSeed,
                                  0, baseDims,
                                  PrestellarCoreInitAbundance,
                                  PrestellarCoreMoleMass,
                                  Turbulence,
                                  SetBaryonField) == FAIL) {
      ENZO_FAIL("Error in PrestellarCoreInitialize[Sub]Grid.");
    } 
    CurrentGrid = CurrentGrid->NextGridThisLevel;
  }

  if (initTurb){
    delete [] Turbulence;
  }

  if (!SetBaryonField)
    return SUCCESS;
  /* Create as many subgrids as there are refinement levels 
     needed to resolve the initial explosion region upon the start-up. */

  LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];;

  for (int level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    LevelArray[level] = NULL;
  AddLevel(LevelArray, &TopGrid, 0);

  // Add levels to the maximum depth or until no new levels are created,
  // and re-initialize the level after it is created.

  for (int level = 0; level < MaximumRefinementLevel; level++) {

    if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
      ENZO_FAIL("Error in RebuildHierarchy.");
    }

    if (LevelArray[level+1] == NULL)
      break;

    if (initTurb) {
      // increase grid dimensional with level to initialize turbulence in higher resolution
      int turbsize=1;
      for (dim = 0; dim < GridRank; dim++) {
        turbsize *= baseDims[dim] * pow(RefineBy, level+1);
      }

      int gridsize = baseDims[0] * pow(RefineBy, level+1);
      Turbulence = new float[GridRank*turbsize];
      if (MyProcessorNumber == ROOT_PROCESSOR)
        printf("turbsize: %d\n", turbsize);
      TurbulenceGenerator_FFTW( PrestellarCoreTurbulenceType, gridsize, PrestellarCoreRandomSeed, 
                                KStart, KEnd, Power, Turbulence);
   
      int m, n, l;
      float x, y, z, rad, density, mass, coreMass;
      float cellWidth[3];
      for (int i=0; i<3; i++) 
        cellWidth[i] = (DomainRightEdge[i] - DomainLeftEdge[i]) / gridsize;
      // calculate the standard deviation used to normalize according to velocity dispersion
      float mean[3] = {0.0}, meanNW[3] = {0.0}, meanC[3] = {0.0}, stdDev[3] = {0.0};
      for (dim = 0; dim < GridRank; dim++){
        mass = 0.0;
        coreMass = 0.0;
        for (int i = 0; i < turbsize; i++){
          m = i % gridsize;
          n = ((i-m) % (gridsize*gridsize))/gridsize;
          l = (i-m-n*gridsize)/gridsize/gridsize;
          x = DomainLeftEdge[0] + (m + 0.5)*(cellWidth[0]) - xcenter;
          y = DomainLeftEdge[1] + (n + 0.5)*(cellWidth[1]) - ycenter;
          z = DomainLeftEdge[2] + (l + 0.5)*(cellWidth[2]) - zcenter;
          rad = sqrt(x*x + y*y + z*z);
          density = PrestellarCoreOuterDensity
              + (PrestellarCoreDensity - PrestellarCoreOuterDensity)
              / (1.0 + pow( (rad/r_f), 1.5))
              * (0.5 - 0.5*tanh( (rad-r_c)/ r_s ));
          mass += density;
  
          mean[dim] += Turbulence[dim + i*GridRank]*density;
          meanNW[dim] += Turbulence[dim + i*GridRank];
          if ( rad <= PrestellarCoreRadius ) {
            meanC[dim] += Turbulence[dim + i*GridRank]*density;
            coreMass += density;
          }
        }
        mean[dim] /= mass;
        mean[dim] /= turbsize;
        meanNW[dim] /= turbsize;
        meanC[dim] /= coreMass;
        meanC[dim] /= turbsize;
        if (MyProcessorNumber == ROOT_PROCESSOR) {
          printf("vel mean value of dim %d: %13.7e \n", dim, mean[dim]);
          printf("no weighted vel mean value of dim %d: %13.7e \n", dim, meanNW[dim]);
          printf("core vel mean value of dim %d: %13.7e \n", dim, meanC[dim]);
        }
        for (int i = 0; i < turbsize; i++){
          stdDev[dim] += pow((Turbulence[dim + i*GridRank]-mean[dim]), 2.0);
        }
        stdDev[dim] = sqrt(stdDev[dim]/turbsize);
        if (MyProcessorNumber == ROOT_PROCESSOR)
          printf("vel standard deviation of dim %d: %13.7e \n", dim, stdDev[dim]);
      }
      for (dim = 0; dim < GridRank; dim++)
        for (int i=0; i<turbsize; i++) {
          Turbulence[dim + i*GridRank] -= meanC[dim];
          Turbulence[dim + i*GridRank] *= PrestellarCoreVelocityDispersion/stdDev[dim];
        }
    }

    LevelHierarchyEntry *Temp = LevelArray[level+1];
    while (Temp) {
      if (Temp->GridData->PrestellarCoreInitializeGrid(
                                    PrestellarCoreRadius,
                                    PrestellarCoreDensity,
                                    PrestellarCoreSurfaceDensity,
                                    PrestellarCoreDensityJump,
                                    PrestellarCoreInternalEnergy,
                                    PrestellarCoreAngularVelocity,
                                    PrestellarCoreBzField,
                                    PrestellarCoreAmbientBzField,
                                    PrestellarCoreVelocityDispersion,
                                    PrestellarCoreTurbulenceKStart,
                                    PrestellarCoreTurbulenceKEnd,
                                    PrestellarCoreOPR,
                                    PrestellarCoreCODeplete,
                                    PrestellarCoreRandomSeed,
                                    level+1, baseDims,
                                    PrestellarCoreInitAbundance,
                                    PrestellarCoreMoleMass,
                                    Turbulence,
                                    SetBaryonField) == FAIL) {
        ENZO_FAIL("Error in PrestellarCoreInitialize[Sub]Grid.");
      } 
      Temp = Temp->NextGridThisLevel;
    }
    if (initTurb){
      delete [] Turbulence;
    }
  }

  // Loop back from the bottom, restoring the consistency among levels. */

  for (int level = MaximumRefinementLevel; level > 0; level--) {
    LevelHierarchyEntry *Temp = LevelArray[level];
    while (Temp != NULL) {
      if (Temp->GridData->ProjectSolutionToParentGrid(
                  *LevelArray[level-1]->GridData) == FAIL) {
        ENZO_FAIL("Error in grid->ProjectSolutionToParentGrid.");
      }
      Temp = Temp->NextGridThisLevel;
    }
  }

  /* set up field names and units */
  int count = 0;
  DataLabel[count++] = DensName;
  DataLabel[count++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[count++] = GEName;
  DataLabel[count++] = Vel1Name;
  DataLabel[count++] = Vel2Name;
  DataLabel[count++] = Vel3Name;
  if ( UseMHD ) {
    DataLabel[count++] = B1Name;
    DataLabel[count++] = B2Name;
    DataLabel[count++] = B3Name;
    DataLabel[count++] = (char*) PhiName;
  }
  if (WritePotential)
    DataLabel[count++] = GPotName;

  if (MultiSpecies == 4) {
    DataLabel[count++] = (char*) DeName;
    DataLabel[count++] = (char*) CHIName;
    DataLabel[count++] = (char*) OIName;
    DataLabel[count++] = (char*) HNCIName;
    DataLabel[count++] = (char*) HCNIName;
    DataLabel[count++] = (char*) H2IName;
    DataLabel[count++] = (char*) CIName;
    DataLabel[count++] = (char*) HIName;
    DataLabel[count++] = (char*) H2OIName;
    DataLabel[count++] = (char*) OHIName;
    DataLabel[count++] = (char*) O2IName;
    DataLabel[count++] = (char*) CH2IName;
    DataLabel[count++] = (char*) H2COIName;
    DataLabel[count++] = (char*) HCOIName;
    DataLabel[count++] = (char*) MGIName;
    DataLabel[count++] = (char*) NH3IName;
    DataLabel[count++] = (char*) NOIName;
    DataLabel[count++] = (char*) CNIName;
    DataLabel[count++] = (char*) COIName;
    DataLabel[count++] = (char*) N2IName;
    DataLabel[count++] = (char*) NH2IName;
    DataLabel[count++] = (char*) CH3IName;
    DataLabel[count++] = (char*) CH4IName;
    DataLabel[count++] = (char*) NIName;
    DataLabel[count++] = (char*) NHIName;
    DataLabel[count++] = (char*) HeIName;
    DataLabel[count++] = (char*) HNOIName;
    DataLabel[count++] = (char*) CH3OHIName;
    DataLabel[count++] = (char*) CO2IName;
    DataLabel[count++] = (char*) H2CNIName;
    DataLabel[count++] = (char*) HNCOIName;
    DataLabel[count++] = (char*) NO2IName;
    DataLabel[count++] = (char*) O2HIName;
    DataLabel[count++] = (char*) OCNIName;
    DataLabel[count++] = (char*) CH3OH_DUSTIName;
    DataLabel[count++] = (char*) HNCO_DUSTIName;
    DataLabel[count++] = (char*) H2CO_DUSTIName;
    DataLabel[count++] = (char*) CH4_DUSTIName;
    DataLabel[count++] = (char*) CO_DUSTIName;
    DataLabel[count++] = (char*) H2O_DUSTIName;
    DataLabel[count++] = (char*) NO_DUSTIName;
    DataLabel[count++] = (char*) CO2_DUSTIName;
    DataLabel[count++] = (char*) N2_DUSTIName;
    DataLabel[count++] = (char*) HCN_DUSTIName;
    DataLabel[count++] = (char*) NH3_DUSTIName;
    DataLabel[count++] = (char*) O2_DUSTIName;
    DataLabel[count++] = (char*) NO2_DUSTIName;
    DataLabel[count++] = (char*) HNO_DUSTIName;
    DataLabel[count++] = (char*) O2H_DUSTIName;
    DataLabel[count++] = (char*) H2CN_DUSTIName;
    DataLabel[count++] = (char*) MG_DUSTIName;
    DataLabel[count++] = (char*) HNC_DUSTIName;
    DataLabel[count++] = (char*) E_DUSTIName;
    DataLabel[count++] = (char*) HCOIIName;
    DataLabel[count++] = (char*) HIIName;
    DataLabel[count++] = (char*) HOCIIName;
    DataLabel[count++] = (char*) CIIName;
    DataLabel[count++] = (char*) CH2IIName;
    DataLabel[count++] = (char*) CHIIName;
    DataLabel[count++] = (char*) H2COIIName;
    DataLabel[count++] = (char*) MGIIName;
    DataLabel[count++] = (char*) NH3IIName;
    DataLabel[count++] = (char*) NOIIName;
    DataLabel[count++] = (char*) CNIIName;
    DataLabel[count++] = (char*) COIIName;
    DataLabel[count++] = (char*) N2IIName;
    DataLabel[count++] = (char*) O2IIName;
    DataLabel[count++] = (char*) H2OIIName;
    DataLabel[count++] = (char*) NH2IIName;
    DataLabel[count++] = (char*) OIIName;
    DataLabel[count++] = (char*) OHIIName;
    DataLabel[count++] = (char*) CH3IIName;
    DataLabel[count++] = (char*) CH4IIName;
    DataLabel[count++] = (char*) NIIName;
    DataLabel[count++] = (char*) HCNIIName;
    DataLabel[count++] = (char*) NHIIName;
    DataLabel[count++] = (char*) H2IIName;
    DataLabel[count++] = (char*) HeIIName;
    DataLabel[count++] = (char*) HNOIIName;
    DataLabel[count++] = (char*) H2NOIIName;
    DataLabel[count++] = (char*) H3IIName;
    DataLabel[count++] = (char*) H3COIIName;
    DataLabel[count++] = (char*) H3OIIName;
    DataLabel[count++] = (char*) HCNHIIName;
    DataLabel[count++] = (char*) HCO2IIName;
    DataLabel[count++] = (char*) HeHIIName;
    DataLabel[count++] = (char*) N2HIIName;
    DataLabel[count++] = (char*) O2HIIName;
    DataLabel[count++] = (char*) HeIIIName;
    DataLabel[count++] = (char*) HMName;
    DataLabel[count++] = (char*) DIName;
    DataLabel[count++] = (char*) DIIName;
    DataLabel[count++] = (char*) HDIName;
  }  // if Multispecies == 4


  for (int j=0; j< count; j++) 
    DataUnits[j] = NULL;

  delete [] PrestellarCoreInitAbundance;
  delete [] PrestellarCoreMoleMass;
  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "PrestellarCoreRadius          = %"FSYM"\n", PrestellarCoreRadius);
    fprintf(Outfptr, "PrestellarCoreDensity         = %"FSYM"\n", PrestellarCoreDensity);
    fprintf(Outfptr, "PrestellarCoreSurfaceDensity  = %"FSYM"\n", PrestellarCoreSurfaceDensity);
    fprintf(Outfptr, "PrestellarCoreDensityJump     = %"FSYM"\n", PrestellarCoreDensityJump);
    fprintf(Outfptr, "PrestellarCoreTemperature     = %"FSYM"\n", PrestellarCoreTemperature);
    fprintf(Outfptr, "PrestellarCoreAngularVelocity = %"FSYM"\n", PrestellarCoreAngularVelocity);
    fprintf(Outfptr, "PrestellarCoreBzField         = %"FSYM"\n", PrestellarCoreBzField);
    fprintf(Outfptr, "PrestellarCoreAmbientBzField  = %"FSYM"\n", PrestellarCoreAmbientBzField);
    fprintf(Outfptr, "PrestellarCoreVelocityDispersion  = %"FSYM"\n", PrestellarCoreVelocityDispersion);
    fprintf(Outfptr, "PrestellarCoreTurbulenceKStart= %"FSYM"\n", PrestellarCoreTurbulenceKStart);
    fprintf(Outfptr, "PrestellarCoreTurbulenceKEnd  = %"FSYM"\n", PrestellarCoreTurbulenceKEnd);
    fprintf(Outfptr, "PrestellarCoreTurbulenceType  = %"ISYM"\n", PrestellarCoreTurbulenceType);
    fprintf(Outfptr, "PrestellarCoreOPR             = %"FSYM"\n", PrestellarCoreOPR);
    fprintf(Outfptr, "PrestellarCoreCODeplete       = %"FSYM"\n", PrestellarCoreCODeplete);
    fprintf(Outfptr, "PrestellarCoreRandomSeed      = %"ISYM"\n", PrestellarCoreRandomSeed);
    fprintf(Outfptr, "Mass Unit                     = %"FSYM"\n", MassUnits);
    fprintf(Outfptr, "G                             = %"FSYM"\n", GravitationalConstant);
    fprintf(Outfptr, "Internal energy               = %"FSYM"\n", PrestellarCoreInternalEnergy);
    fprintf(Outfptr, "Energy Unit                   = %"FSYM"\n", pow(LengthUnits, 2)/pow(TimeUnits,2));
    fprintf(Outfptr, "PrestellarCoreInitAbundance, PrestellarCoreMoleMass\n");
    for (int i=0; i<NSpecies+1; i++){
      fprintf(Outfptr, "%13.7e                , %13.7e\n", PrestellarCoreInitAbundance[i], PrestellarCoreMoleMass[i]);
    }
  }

  return SUCCESS;

}
