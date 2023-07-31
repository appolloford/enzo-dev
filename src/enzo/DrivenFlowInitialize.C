/***********************************************************************
/
/  INITIALIZE DRIVEN FLOW SIMULATION
/
/  written by: Wolfram Schmidt
/  date:       May, 2005
/  modified1:  April, 2007
/  modified2: Sep, 2014: updated to support Enzo 2.4   // P. Grete
/
/  PURPOSE: Initializes simulation of flow driven by stochastic forcing
/
/  RETURNS: SUCCESS or FAIL
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
#include "Hierarchy.h"
#include "TopGridData.h"
#include "phys_constants.h"

StochasticForcing Forcing;
void MHDCTSetupFieldLabels();

int DrivenFlowInitialize(FILE *fptr, FILE *Outfptr, 
             HierarchyEntry &TopGrid, TopGridData &MetaData, 
             int SetBaryonFields)
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
  char *StochAccel1Name = "x-acceleration";
  char *StochAccel2Name = "y-acceleration";
  char *StochAccel3Name = "z-acceleration";

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

#ifdef USE_NAUNET
  const char *GH2CNIName             = "GH2CNI_Density";
  const char *GHNCIName              = "GHNCI_Density";
  const char *GNO2IName              = "GNO2I_Density";
  const char *GSiOIName              = "GSiOI_Density";
  const char *GCOIName               = "GCOI_Density";
  const char *GHNCOIName             = "GHNCOI_Density";
  const char *GMgIName               = "GMgI_Density";
  const char *GNOIName               = "GNOI_Density";
  const char *GO2IName               = "GO2I_Density";
  const char *GO2HIName              = "GO2HI_Density";
  const char *GSiCIName              = "GSiCI_Density";
  const char *GSiC2IName             = "GSiC2I_Density";
  const char *GSiC3IName             = "GSiC3I_Density";
  const char *GCH3OHIName            = "GCH3OHI_Density";
  const char *GCO2IName              = "GCO2I_Density";
  const char *GH2SiOIName            = "GH2SiOI_Density";
  const char *GHNOIName              = "GHNOI_Density";
  const char *GN2IName               = "GN2I_Density";
  const char *GH2COIName             = "GH2COI_Density";
  const char *GHCNIName              = "GHCNI_Density";
  const char *GH2OIName              = "GH2OI_Density";
  const char *GNH3IName              = "GNH3I_Density";
  const char *SiC3IIName             = "SiC3II_Density";
  const char *H2CNIName              = "H2CNI_Density";
  const char *GCH4IName              = "GCH4I_Density";
  const char *H2NOIIName             = "H2NOII_Density";
  const char *H2SiOIName             = "H2SiOI_Density";
  const char *HeHIIName              = "HeHII_Density";
  const char *HNCOIName              = "HNCOI_Density";
  const char *HOCIIName              = "HOCII_Density";
  const char *SiC2IIName             = "SiC2II_Density";
  const char *GSiH4IName             = "GSiH4I_Density";
  const char *SiC2IName              = "SiC2I_Density";
  const char *SiC3IName              = "SiC3I_Density";
  const char *SiH5IIName             = "SiH5II_Density";
  const char *SiH4IIName             = "SiH4II_Density";
  const char *SiCIIName              = "SiCII_Density";
  const char *O2HIName               = "O2HI_Density";
  const char *SiCIName               = "SiCI_Density";
  const char *NO2IName               = "NO2I_Density";
  const char *SiH3IIName             = "SiH3II_Density";
  const char *SiH2IIName             = "SiH2II_Density";
  const char *OCNIName               = "OCNI_Density";
  const char *SiH2IName              = "SiH2I_Density";
  const char *SiOHIIName             = "SiOHII_Density";
  const char *SiHIIName              = "SiHII_Density";
  const char *SiH4IName              = "SiH4I_Density";
  const char *SiHIName               = "SiHI_Density";
  const char *SiH3IName              = "SiH3I_Density";
  const char *SiOIIName              = "SiOII_Density";
  const char *HCO2IIName             = "HCO2II_Density";
  const char *HNOIName               = "HNOI_Density";
  const char *CH3OHIName             = "CH3OHI_Density";
  const char *MgIName                = "MgI_Density";
  const char *MgIIName               = "MgII_Density";
  const char *CH4IIName              = "CH4II_Density";
  const char *SiOIName               = "SiOI_Density";
  const char *CNIIName               = "CNII_Density";
  const char *HCNHIIName             = "HCNHII_Density";
  const char *N2HIIName              = "N2HII_Density";
  const char *O2HIIName              = "O2HII_Density";
  const char *SiIIName               = "SiII_Density";
  const char *SiIName                = "SiI_Density";
  const char *HNCIName               = "HNCI_Density";
  const char *HNOIIName              = "HNOII_Density";
  const char *N2IIName               = "N2II_Density";
  const char *H3COIIName             = "H3COII_Density";
  const char *CH4IName               = "CH4I_Density";
  const char *COIIName               = "COII_Density";
  const char *NH3IName               = "NH3I_Density";
  const char *CH3IName               = "CH3I_Density";
  const char *CO2IName               = "CO2I_Density";
  const char *NIIName                = "NII_Density";
  const char *OIIName                = "OII_Density";
  const char *HCNIIName              = "HCNII_Density";
  const char *NH2IIName              = "NH2II_Density";
  const char *NHIIName               = "NHII_Density";
  const char *O2IIName               = "O2II_Density";
  const char *CH3IIName              = "CH3II_Density";
  const char *NH2IName               = "NH2I_Density";
  const char *CH2IIName              = "CH2II_Density";
  const char *H2OIIName              = "H2OII_Density";
  const char *NH3IIName              = "NH3II_Density";
  const char *NOIIName               = "NOII_Density";
  const char *H3OIIName              = "H3OII_Density";
  const char *N2IName                = "N2I_Density";
  const char *CIIName                = "CII_Density";
  const char *HCNIName               = "HCNI_Density";
  const char *CHIIName               = "CHII_Density";
  const char *CH2IName               = "CH2I_Density";
  const char *H2COIIName             = "H2COII_Density";
  const char *NHIName                = "NHI_Density";
  const char *OHIIName               = "OHII_Density";
  const char *CNIName                = "CNI_Density";
  const char *H2COIName              = "H2COI_Density";
  const char *HCOIName               = "HCOI_Density";
  const char *CHIName                = "CHI_Density";
  const char *H3IIName               = "H3II_Density";
  const char *NOIName                = "NOI_Density";
  const char *NIName                 = "NI_Density";
  const char *OHIName                = "OHI_Density";
  const char *O2IName                = "O2I_Density";
  const char *CIName                 = "CI_Density";
  const char *HCOIIName              = "HCOII_Density";
  const char *H2OIName               = "H2OI_Density";
  const char *OIName                 = "OI_Density";
  const char *COIName                = "COI_Density";

  const char *KspName                = "Sputtering_Rate";
#endif


  /* declarations */

  char line[MAX_LINE_LENGTH];
  int ret;

  /* set default parameters specifying the random force field */

  int DrivenFlowAlpha[3]       = {1, 1, 1};       // ratio of domain size to characteristic length
  int DrivenFlowSeed = 20150418;                  // seed of random number generator
  float DrivenFlowBandWidth[3] = {1.0, 1.0, 1.0}; // band width (1.0 = maximal)
  float DrivenFlowAutoCorrl[3] = {1.0, 1.0, 1.0}; // ratio auto-correlation to large-eddy turn-over time scale
  float DrivenFlowMach[3]      = {1.0, 1.0, 1.0}; // Mach number
  float DrivenFlowWeight       = 1.0;             // weight of solenoidal components

  /* set other default parameters */

  float DrivenFlowDensity     = 1.0; // initial mass density
  float DrivenFlowPressure    = 1.0; // initial pressure

  float DrivenFlowMagField           = 0.0; // initial magnetic field

  forcing_type DrivenFlowProfile;            //defined in typedefs.h
  float DrivenFlowDomainLength[3];
  float DrivenFlowVelocity[3];
  float SoundSpeed;

  /* read input from file */
  if (debug) printf("DrivenFlowInitialize: reading problem-specific parameters.\n");

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "DrivenFlowProfile = %"ISYM, &DrivenFlowProfile);
    
    ret += sscanf(line, "DrivenFlowAlpha = %"ISYM" %"ISYM" %"ISYM, 
          DrivenFlowAlpha, DrivenFlowAlpha+1, DrivenFlowAlpha+2);
    ret += sscanf(line, "DrivenFlowSeed = %"ISYM, &DrivenFlowSeed);

    ret += sscanf(line, "DrivenFlowBandWidth = %"FSYM"%"FSYM"%"FSYM, 
          DrivenFlowBandWidth, DrivenFlowBandWidth+1, DrivenFlowBandWidth+2);

    ret += sscanf(line, "DrivenFlowAutoCorrl = %"FSYM"%"FSYM"%"FSYM, 
          DrivenFlowAutoCorrl, DrivenFlowAutoCorrl+1, DrivenFlowAutoCorrl+2);

    ret += sscanf(line, "DrivenFlowMach = %"FSYM"%"FSYM"%"FSYM, 
          DrivenFlowMach, DrivenFlowMach+1, DrivenFlowMach+2);

    ret += sscanf(line, "DrivenFlowWeight = %"FSYM, &DrivenFlowWeight);

    ret += sscanf(line, "DrivenFlowDensity = %"FSYM, &DrivenFlowDensity);
    ret += sscanf(line, "DrivenFlowPressure = %"FSYM, &DrivenFlowPressure);
    ret += sscanf(line, "DrivenFlowMagField = %"FSYM, &DrivenFlowMagField);

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "DrivenFlow"))
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  }


  /* thermodynamic initial values */

#if USE_NAUNET
  if (MultiSpecies == NAUNET_SPECIES) {
#else
  if (MultiSpecies == 0) {
#endif
    if (EquationOfState == 0)
      SoundSpeed = sqrt(Gamma * DrivenFlowPressure / DrivenFlowDensity);
    else
      SoundSpeed = IsothermalSoundSpeed;
    
  } else {
    fprintf(stderr,"DrivenFlowInitialize: Multispecies != 0 untested at this point.\n");
    return FALSE;
  }

  if (SelfGravity) {
      fprintf(stderr,"DrivenFlowInitialize: SelfGravity untested at this point.\n");
      return FALSE;
  }

  if ((HydroMethod != MHD_RK) && (HydroMethod != HD_RK) && (HydroMethod != MHD_Li)) {
      fprintf(stderr,"DrivenFlowInitialize: Only support for MUSCL framework and MHDCT at this point.\n");
      return FALSE;
  }

  if (MetaData.TopGridRank != 3) {
      fprintf(stderr,"DrivenFlowInitialize: Only 3D tested at this point.\n");
      return FALSE;
  }


  // set proper internal unit for magnetic field
  DrivenFlowMagField /= sqrt(4*pi);
  
  /* compute characteristic velocity from speed of sound and Mach numbers */

  for (int dim = 0; dim < MetaData.TopGridRank; dim++) {
      DrivenFlowVelocity[dim] = DrivenFlowMach[dim] * SoundSpeed;
      DrivenFlowDomainLength[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
      if (debug)
          printf("dim = %"ISYM" vel = %"FSYM" len = %"FSYM"\n",
            dim,DrivenFlowVelocity[dim],DrivenFlowDomainLength[dim]);
  }

  /* Begin grid initialization */

  HierarchyEntry *CurrentGrid; // all level 0 grids on this processor first
  CurrentGrid = &TopGrid;
  while (CurrentGrid != NULL) {
      if (CurrentGrid->GridData->DrivenFlowInitializeGrid(DrivenFlowDensity,
              DrivenFlowPressure,DrivenFlowMagField,SetBaryonFields) == FAIL) {
          fprintf(stderr, "Error in DrivenFlowInitializeGrid.\n");
          return FAIL;
      }
      CurrentGrid = CurrentGrid->NextGridThisLevel;
  }
  
  if (SetBaryonFields) {
  /* create a stochasitc forcing object with the specified parameters */
  Forcing.Init(MetaData.TopGridRank,
           DrivenFlowProfile,
           DrivenFlowAlpha,
           DrivenFlowDomainLength,
           DrivenFlowBandWidth,
           DrivenFlowVelocity,
           DrivenFlowAutoCorrl,
           DrivenFlowWeight,
           DrivenFlowSeed);

  /* set up field names and units */

  int count = 0;
  DataLabel[count++] = DensName;
  
  DataLabel[count++] = Vel1Name;
  DataLabel[count++] = Vel2Name;
  DataLabel[count++] = Vel3Name;


  if(EquationOfState == 0)
    DataLabel[count++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[count++] = GEName;

  if ( UseMHD ) {
    DataLabel[count++] = BxName;
    DataLabel[count++] = ByName;
    DataLabel[count++] = BzName;
  }
  if( HydroMethod == MHD_RK ){
    DataLabel[count++] = PhiName;
  }
  MHDCTSetupFieldLabels();

    DataLabel[count++] = StochAccel1Name;
    DataLabel[count++] = StochAccel2Name;
    DataLabel[count++] = StochAccel3Name;

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
      DataLabel[count++] = (char*) GH2CNIName;
      DataLabel[count++] = (char*) GHNCIName;
      DataLabel[count++] = (char*) GNO2IName;
      DataLabel[count++] = (char*) GSiOIName;
      DataLabel[count++] = (char*) GCOIName;
      DataLabel[count++] = (char*) GHNCOIName;
      DataLabel[count++] = (char*) GMgIName;
      DataLabel[count++] = (char*) GNOIName;
      DataLabel[count++] = (char*) GO2IName;
      DataLabel[count++] = (char*) GO2HIName;
      DataLabel[count++] = (char*) GSiCIName;
      DataLabel[count++] = (char*) GSiC2IName;
      DataLabel[count++] = (char*) GSiC3IName;
      DataLabel[count++] = (char*) GCH3OHIName;
      DataLabel[count++] = (char*) GCO2IName;
      DataLabel[count++] = (char*) GH2SiOIName;
      DataLabel[count++] = (char*) GHNOIName;
      DataLabel[count++] = (char*) GN2IName;
      DataLabel[count++] = (char*) GH2COIName;
      DataLabel[count++] = (char*) GHCNIName;
      DataLabel[count++] = (char*) GH2OIName;
      DataLabel[count++] = (char*) GNH3IName;
      DataLabel[count++] = (char*) SiC3IIName;
      DataLabel[count++] = (char*) H2CNIName;
      DataLabel[count++] = (char*) GCH4IName;
      DataLabel[count++] = (char*) H2NOIIName;
      DataLabel[count++] = (char*) H2SiOIName;
      DataLabel[count++] = (char*) HeHIIName;
      DataLabel[count++] = (char*) HNCOIName;
      DataLabel[count++] = (char*) HOCIIName;
      DataLabel[count++] = (char*) SiC2IIName;
      DataLabel[count++] = (char*) GSiH4IName;
      DataLabel[count++] = (char*) SiC2IName;
      DataLabel[count++] = (char*) SiC3IName;
      DataLabel[count++] = (char*) SiH5IIName;
      DataLabel[count++] = (char*) SiH4IIName;
      DataLabel[count++] = (char*) SiCIIName;
      DataLabel[count++] = (char*) O2HIName;
      DataLabel[count++] = (char*) SiCIName;
      DataLabel[count++] = (char*) NO2IName;
      DataLabel[count++] = (char*) SiH3IIName;
      DataLabel[count++] = (char*) SiH2IIName;
      DataLabel[count++] = (char*) OCNIName;
      DataLabel[count++] = (char*) SiH2IName;
      DataLabel[count++] = (char*) SiOHIIName;
      DataLabel[count++] = (char*) SiHIIName;
      DataLabel[count++] = (char*) SiH4IName;
      DataLabel[count++] = (char*) SiHIName;
      DataLabel[count++] = (char*) SiH3IName;
      DataLabel[count++] = (char*) SiOIIName;
      DataLabel[count++] = (char*) HCO2IIName;
      DataLabel[count++] = (char*) HNOIName;
      DataLabel[count++] = (char*) CH3OHIName;
      DataLabel[count++] = (char*) MgIName;
      DataLabel[count++] = (char*) MgIIName;
      DataLabel[count++] = (char*) CH4IIName;
      DataLabel[count++] = (char*) SiOIName;
      DataLabel[count++] = (char*) CNIIName;
      DataLabel[count++] = (char*) HCNHIIName;
      DataLabel[count++] = (char*) N2HIIName;
      DataLabel[count++] = (char*) O2HIIName;
      DataLabel[count++] = (char*) SiIIName;
      DataLabel[count++] = (char*) SiIName;
      DataLabel[count++] = (char*) HNCIName;
      DataLabel[count++] = (char*) HNOIIName;
      DataLabel[count++] = (char*) N2IIName;
      DataLabel[count++] = (char*) H3COIIName;
      DataLabel[count++] = (char*) CH4IName;
      DataLabel[count++] = (char*) COIIName;
      DataLabel[count++] = (char*) NH3IName;
      DataLabel[count++] = (char*) CH3IName;
      DataLabel[count++] = (char*) CO2IName;
      DataLabel[count++] = (char*) NIIName;
      DataLabel[count++] = (char*) OIIName;
      DataLabel[count++] = (char*) HCNIIName;
      DataLabel[count++] = (char*) NH2IIName;
      DataLabel[count++] = (char*) NHIIName;
      DataLabel[count++] = (char*) O2IIName;
      DataLabel[count++] = (char*) CH3IIName;
      DataLabel[count++] = (char*) NH2IName;
      DataLabel[count++] = (char*) CH2IIName;
      DataLabel[count++] = (char*) H2OIIName;
      DataLabel[count++] = (char*) NH3IIName;
      DataLabel[count++] = (char*) NOIIName;
      DataLabel[count++] = (char*) H3OIIName;
      DataLabel[count++] = (char*) N2IName;
      DataLabel[count++] = (char*) CIIName;
      DataLabel[count++] = (char*) HCNIName;
      DataLabel[count++] = (char*) CHIIName;
      DataLabel[count++] = (char*) CH2IName;
      DataLabel[count++] = (char*) H2COIIName;
      DataLabel[count++] = (char*) NHIName;
      DataLabel[count++] = (char*) OHIIName;
      DataLabel[count++] = (char*) CNIName;
      DataLabel[count++] = (char*) H2COIName;
      DataLabel[count++] = (char*) HCOIName;
      DataLabel[count++] = (char*) CHIName;
      DataLabel[count++] = (char*) H3IIName;
      DataLabel[count++] = (char*) NOIName;
      DataLabel[count++] = (char*) NIName;
      DataLabel[count++] = (char*) OHIName;
      DataLabel[count++] = (char*) O2IName;
      DataLabel[count++] = (char*) CIName;
      DataLabel[count++] = (char*) HCOIIName;
      DataLabel[count++] = (char*) H2OIName;
      DataLabel[count++] = (char*) OIName;
      DataLabel[count++] = (char*) COIName;

      DataLabel[count++] = (char*) KspName;
    }
#endif
  }

  /* Write parameters to parameter output file */

  if (debug) printf("DrivenFlowInitialize: writing parameters to output file.\n");

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "Profile    = %"ISYM"\n\n", DrivenFlowProfile);

    fprintf(Outfptr, "Alpha1     = %"ISYM"\n",   DrivenFlowAlpha[0]);
    fprintf(Outfptr, "Alpha2     = %"ISYM"\n",   DrivenFlowAlpha[1]);
    fprintf(Outfptr, "Alpha3     = %"ISYM"\n\n", DrivenFlowAlpha[2]);

    fprintf(Outfptr, "BandWidth1 = %"FSYM"\n",   DrivenFlowBandWidth[0]);
    fprintf(Outfptr, "BandWidth2 = %"FSYM"\n",   DrivenFlowBandWidth[1]);
    fprintf(Outfptr, "BandWidth3 = %"FSYM"\n\n", DrivenFlowBandWidth[2]);

    fprintf(Outfptr, "AutoCorrl1 = %"FSYM"\n",   DrivenFlowAutoCorrl[0]);
    fprintf(Outfptr, "AutoCorrl2 = %"FSYM"\n",   DrivenFlowAutoCorrl[1]);
    fprintf(Outfptr, "AutoCorrl3 = %"FSYM"\n\n", DrivenFlowAutoCorrl[2]);

    fprintf(Outfptr, "SolnWeight = %"FSYM"\n\n", DrivenFlowWeight);

    fprintf(Outfptr, "Density    = %"FSYM"\n",   DrivenFlowDensity);
    fprintf(Outfptr, "Pressure   = %"FSYM"\n\n", DrivenFlowPressure);

    fprintf(Outfptr, "Velocity1  = %"FSYM"\n",   DrivenFlowVelocity[0]);
    fprintf(Outfptr, "Velocity2  = %"FSYM"\n",   DrivenFlowVelocity[1]);
    fprintf(Outfptr, "Velocity3  = %"FSYM"\n\n", DrivenFlowVelocity[2]);
  }
  }
  return SUCCESS;
}
