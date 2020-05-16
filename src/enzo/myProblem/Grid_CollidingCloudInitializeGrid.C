/***********************************************************************
/
/  GRID CLASS (INITIALIZE HYDRO/MHD/RADHYDRO TURBULENCE SIMULATION)
/
/  written by: Peng Wang
/  date:       September, 2007
/  modified1: Tom Abel, allow for parallel generations of ICs (Sept 2009)
/
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "../hydro_rk/EOS.h"
#include "phys_constants.h"

//int FindField(int f, int farray[], int n);  // tried adding for GMC color using InitializeUniform

int GetUnits(float *DensityUnits, float *LengthUnits, float *TemperatureUnits, 
             float *TimeUnits, float *VelocityUnits, FLOAT Time);
void Turbulence_Generator(float **vel, int dim0, int dim1, int dim2, int ind, 
                          float kmin, float kmax, float dk,
                          FLOAT **LeftEdge, FLOAT **CellWidth, int seed);

int grid::CollidingCloudInitializeGrid(float CloudDensity, float CloudSoundSpeed, FLOAT CloudRadius, 
        float CloudMachNumber, float CloudAngularVelocity, float InitialBField,
        int SetTurbulence, int CloudType, int TurbulenceSeed, int PutSink, int level,
        int SetBaryonFields, float RelativeVelocity, float Btheta, float ImpactParameter)  // GMC collision parameters
{

  /* declarations */

  int dim, i, j, k,l, m, n, field, sphere, size, igrid, activesize;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum,  kphHINum, gammaNum, kphHeINum,
    kphHeIINum, kdissH2INum, RPresNum1, RPresNum2, RPresNum3;
  //int GMC1Num, GMC2Num;
  ////GMC1Num = FindField(GMC1, FieldType, NumberOfBaryonFields);
  ////GMC2Num = FindField(GMC2, FieldType, NumberOfBaryonFields);
  ////int ColourNum;
  ////ColourNum = 0; //
#ifdef USE_KROME
  int CHINum, OINum, HNCINum, HCNINum, CINum, H2OINum,
 OHINum, O2INum, CH2INum, H2COINum, HCOINum,
 MGINum, NH3INum, NOINum, CNINum, COINum,
 N2INum, NH2INum, CH3INum, CH4INum, NINum,
 NHINum, HNOINum, CH3OHINum, CO2INum, H2CNINum,
 HNCOINum, NO2INum, O2HINum, OCNINum, CH3OH_DUSTINum,
 HNCO_DUSTINum, H2CO_DUSTINum, CH4_DUSTINum,
 CO_DUSTINum, H2O_DUSTINum, NO_DUSTINum, CO2_DUSTINum,
 N2_DUSTINum, HCN_DUSTINum, NH3_DUSTINum,
 O2_DUSTINum, NO2_DUSTINum, HNO_DUSTINum,
 O2H_DUSTINum, H2CN_DUSTINum, MG_DUSTINum,
 HNC_DUSTINum, E_DUSTINum, HCOIINum, HOCIINum,
 CIINum, CH2IINum, CHIINum, H2COIINum, MGIINum,
 NH3IINum, NOIINum, CNIINum, COIINum, N2IINum,
 O2IINum, H2OIINum, NH2IINum, OIINum, OHIINum,
 CH3IINum, CH4IINum, NIINum, HCNIINum, NHIINum,
 HNOIINum, H2NOIINum, H3IINum, H3COIINum,
 H3OIINum, HCNHIINum, HCO2IINum, HeHIINum,
 N2HIINum, O2HIINum;
  int speciesMap[NKROMESPECIES] = {-1};
#endif


  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;
  //  if(HydroMethod == Zeus_Hydro)
  //  FieldType[NumberOfBaryonFields++] = InternalEnergy;
  //else  
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  

  if (DualEnergyFormalism) {
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  }
  if (UseMHD) {
    FieldType[NumberOfBaryonFields++] = Bfield1;
    FieldType[NumberOfBaryonFields++] = Bfield2;
    FieldType[NumberOfBaryonFields++] = Bfield3;
  }
  if ( HydroMethod == MHD_RK){
    FieldType[NumberOfBaryonFields++] = PhiField;
  }
  if(UsePoissonDivergenceCleaning)
    FieldType[NumberOfBaryonFields++] = Phi_pField;
  
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
#ifdef USE_KROME
    if (MultiSpecies == KROMESPECIES) {
      FieldType[CHINum          = NumberOfBaryonFields++] = CHIDensity;
      FieldType[OINum           = NumberOfBaryonFields++] = OIDensity;
      FieldType[HNCINum         = NumberOfBaryonFields++] = HNCIDensity;
      FieldType[HCNINum         = NumberOfBaryonFields++] = HCNIDensity;
      FieldType[CINum           = NumberOfBaryonFields++] = CIDensity;
      FieldType[H2OINum         = NumberOfBaryonFields++] = H2OIDensity;
      FieldType[OHINum          = NumberOfBaryonFields++] = OHIDensity;
      FieldType[O2INum          = NumberOfBaryonFields++] = O2IDensity;
      FieldType[CH2INum         = NumberOfBaryonFields++] = CH2IDensity;
      FieldType[H2COINum        = NumberOfBaryonFields++] = H2COIDensity;
      FieldType[HCOINum         = NumberOfBaryonFields++] = HCOIDensity;
      FieldType[MGINum          = NumberOfBaryonFields++] = MGIDensity;
      FieldType[NH3INum         = NumberOfBaryonFields++] = NH3IDensity;
      FieldType[NOINum          = NumberOfBaryonFields++] = NOIDensity;
      FieldType[CNINum          = NumberOfBaryonFields++] = CNIDensity;
      FieldType[COINum          = NumberOfBaryonFields++] = COIDensity;
      FieldType[N2INum          = NumberOfBaryonFields++] = N2IDensity;
      FieldType[NH2INum         = NumberOfBaryonFields++] = NH2IDensity;
      FieldType[CH3INum         = NumberOfBaryonFields++] = CH3IDensity;
      FieldType[CH4INum         = NumberOfBaryonFields++] = CH4IDensity;
      FieldType[NINum           = NumberOfBaryonFields++] = NIDensity;
      FieldType[NHINum          = NumberOfBaryonFields++] = NHIDensity;
      FieldType[HNOINum         = NumberOfBaryonFields++] = HNOIDensity;
      FieldType[CH3OHINum       = NumberOfBaryonFields++] = CH3OHIDensity;
      FieldType[CO2INum         = NumberOfBaryonFields++] = CO2IDensity;
      FieldType[H2CNINum        = NumberOfBaryonFields++] = H2CNIDensity;
      FieldType[HNCOINum        = NumberOfBaryonFields++] = HNCOIDensity;
      FieldType[NO2INum         = NumberOfBaryonFields++] = NO2IDensity;
      FieldType[O2HINum         = NumberOfBaryonFields++] = O2HIDensity;
      FieldType[OCNINum         = NumberOfBaryonFields++] = OCNIDensity;
      FieldType[CH3OH_DUSTINum  = NumberOfBaryonFields++] = CH3OH_DUSTIDensity;
      FieldType[HNCO_DUSTINum   = NumberOfBaryonFields++] = HNCO_DUSTIDensity;
      FieldType[H2CO_DUSTINum   = NumberOfBaryonFields++] = H2CO_DUSTIDensity;
      FieldType[CH4_DUSTINum    = NumberOfBaryonFields++] = CH4_DUSTIDensity;
      FieldType[CO_DUSTINum     = NumberOfBaryonFields++] = CO_DUSTIDensity;
      FieldType[H2O_DUSTINum    = NumberOfBaryonFields++] = H2O_DUSTIDensity;
      FieldType[NO_DUSTINum     = NumberOfBaryonFields++] = NO_DUSTIDensity;
      FieldType[CO2_DUSTINum    = NumberOfBaryonFields++] = CO2_DUSTIDensity;
      FieldType[N2_DUSTINum     = NumberOfBaryonFields++] = N2_DUSTIDensity;
      FieldType[HCN_DUSTINum    = NumberOfBaryonFields++] = HCN_DUSTIDensity;
      FieldType[NH3_DUSTINum    = NumberOfBaryonFields++] = NH3_DUSTIDensity;
      FieldType[O2_DUSTINum     = NumberOfBaryonFields++] = O2_DUSTIDensity;
      FieldType[NO2_DUSTINum    = NumberOfBaryonFields++] = NO2_DUSTIDensity;
      FieldType[HNO_DUSTINum    = NumberOfBaryonFields++] = HNO_DUSTIDensity;
      FieldType[O2H_DUSTINum    = NumberOfBaryonFields++] = O2H_DUSTIDensity;
      FieldType[H2CN_DUSTINum   = NumberOfBaryonFields++] = H2CN_DUSTIDensity;
      FieldType[MG_DUSTINum     = NumberOfBaryonFields++] = MG_DUSTIDensity;
      FieldType[HNC_DUSTINum    = NumberOfBaryonFields++] = HNC_DUSTIDensity;
      FieldType[E_DUSTINum      = NumberOfBaryonFields++] = E_DUSTIDensity;
      FieldType[HCOIINum        = NumberOfBaryonFields++] = HCOIIDensity;
      FieldType[HOCIINum        = NumberOfBaryonFields++] = HOCIIDensity;
      FieldType[CIINum          = NumberOfBaryonFields++] = CIIDensity;
      FieldType[CH2IINum        = NumberOfBaryonFields++] = CH2IIDensity;
      FieldType[CHIINum         = NumberOfBaryonFields++] = CHIIDensity;
      FieldType[H2COIINum       = NumberOfBaryonFields++] = H2COIIDensity;
      FieldType[MGIINum         = NumberOfBaryonFields++] = MGIIDensity;
      FieldType[NH3IINum        = NumberOfBaryonFields++] = NH3IIDensity;
      FieldType[NOIINum         = NumberOfBaryonFields++] = NOIIDensity;
      FieldType[CNIINum         = NumberOfBaryonFields++] = CNIIDensity;
      FieldType[COIINum         = NumberOfBaryonFields++] = COIIDensity;
      FieldType[N2IINum         = NumberOfBaryonFields++] = N2IIDensity;
      FieldType[O2IINum         = NumberOfBaryonFields++] = O2IIDensity;
      FieldType[H2OIINum        = NumberOfBaryonFields++] = H2OIIDensity;
      FieldType[NH2IINum        = NumberOfBaryonFields++] = NH2IIDensity;
      FieldType[OIINum          = NumberOfBaryonFields++] = OIIDensity;
      FieldType[OHIINum         = NumberOfBaryonFields++] = OHIIDensity;
      FieldType[CH3IINum        = NumberOfBaryonFields++] = CH3IIDensity;
      FieldType[CH4IINum        = NumberOfBaryonFields++] = CH4IIDensity;
      FieldType[NIINum          = NumberOfBaryonFields++] = NIIDensity;
      FieldType[HCNIINum        = NumberOfBaryonFields++] = HCNIIDensity;
      FieldType[NHIINum         = NumberOfBaryonFields++] = NHIIDensity;
      FieldType[HNOIINum        = NumberOfBaryonFields++] = HNOIIDensity;
      FieldType[H2NOIINum       = NumberOfBaryonFields++] = H2NOIIDensity;
      FieldType[H3IINum         = NumberOfBaryonFields++] = H3IIDensity;
      FieldType[H3COIINum       = NumberOfBaryonFields++] = H3COIIDensity;
      FieldType[H3OIINum        = NumberOfBaryonFields++] = H3OIIDensity;
      FieldType[HCNHIINum       = NumberOfBaryonFields++] = HCNHIIDensity;
      FieldType[HCO2IINum       = NumberOfBaryonFields++] = HCO2IIDensity;
      FieldType[HeHIINum        = NumberOfBaryonFields++] = HeHIIDensity;
      FieldType[N2HIINum        = NumberOfBaryonFields++] = N2HIIDensity;
      FieldType[O2HIINum        = NumberOfBaryonFields++] = O2HIIDensity;


      speciesMap[0] = DeNum;
      speciesMap[1] = CHINum;
      speciesMap[2] = OINum;
      speciesMap[3] = HNCINum;
      speciesMap[4] = HCNINum;
      speciesMap[5] = H2INum;
      speciesMap[6] = CINum;
      speciesMap[7] = HINum;
      speciesMap[8] = H2OINum;
      speciesMap[9] = OHINum;
      speciesMap[10] = O2INum;
      speciesMap[11] = CH2INum;
      speciesMap[12] = H2COINum;
      speciesMap[13] = HCOINum;
      speciesMap[14] = MGINum;
      speciesMap[15] = NH3INum;
      speciesMap[16] = NOINum;
      speciesMap[17] = CNINum;
      speciesMap[18] = COINum;
      speciesMap[19] = N2INum;
      speciesMap[20] = NH2INum;
      speciesMap[21] = CH3INum;
      speciesMap[22] = CH4INum;
      speciesMap[23] = NINum;
      speciesMap[24] = NHINum;
      speciesMap[25] = HeINum;
      speciesMap[26] = HNOINum;
      speciesMap[27] = CH3OHINum;
      speciesMap[28] = CO2INum;
      speciesMap[29] = H2CNINum;
      speciesMap[30] = HNCOINum;
      speciesMap[31] = NO2INum;
      speciesMap[32] = O2HINum;
      speciesMap[33] = OCNINum;
      speciesMap[34] = CH3OH_DUSTINum;
      speciesMap[35] = HNCO_DUSTINum;
      speciesMap[36] = H2CO_DUSTINum;
      speciesMap[37] = CH4_DUSTINum;
      speciesMap[38] = CO_DUSTINum;
      speciesMap[39] = H2O_DUSTINum;
      speciesMap[40] = NO_DUSTINum;
      speciesMap[41] = CO2_DUSTINum;
      speciesMap[42] = N2_DUSTINum;
      speciesMap[43] = HCN_DUSTINum;
      speciesMap[44] = NH3_DUSTINum;
      speciesMap[45] = O2_DUSTINum;
      speciesMap[46] = NO2_DUSTINum;
      speciesMap[47] = HNO_DUSTINum;
      speciesMap[48] = O2H_DUSTINum;
      speciesMap[49] = H2CN_DUSTINum;
      speciesMap[50] = MG_DUSTINum;
      speciesMap[51] = HNC_DUSTINum;
      speciesMap[52] = E_DUSTINum;
      speciesMap[53] = HCOIINum;
      speciesMap[54] = HIINum;
      speciesMap[55] = HOCIINum;
      speciesMap[56] = CIINum;
      speciesMap[57] = CH2IINum;
      speciesMap[58] = CHIINum;
      speciesMap[59] = H2COIINum;
      speciesMap[60] = MGIINum;
      speciesMap[61] = NH3IINum;
      speciesMap[62] = NOIINum;
      speciesMap[63] = CNIINum;
      speciesMap[64] = COIINum;
      speciesMap[65] = N2IINum;
      speciesMap[66] = O2IINum;
      speciesMap[67] = H2OIINum;
      speciesMap[68] = NH2IINum;
      speciesMap[69] = OIINum;
      speciesMap[70] = OHIINum;
      speciesMap[71] = CH3IINum;
      speciesMap[72] = CH4IINum;
      speciesMap[73] = NIINum;
      speciesMap[74] = HCNIINum;
      speciesMap[75] = NHIINum;
      speciesMap[76] = H2IINum;
      speciesMap[77] = HeIINum;
      speciesMap[78] = HNOIINum;
      speciesMap[79] = H2NOIINum;
      speciesMap[80] = H3IINum;
      speciesMap[81] = H3COIINum;
      speciesMap[82] = H3OIINum;
      speciesMap[83] = HCNHIINum;
      speciesMap[84] = HCO2IINum;
      speciesMap[85] = HeHIINum;
      speciesMap[86] = N2HIINum;
      speciesMap[87] = O2HIINum;

    }
#endif

  }
  //  FieldType[ColourNum = NumberOfBaryonFields++] = Metallicity;

  if (RadiativeTransfer && (MultiSpecies < 1)) {
    fprintf(stderr, "Grid_PhotonTestInitialize: Radiative Transfer but not MultiSpecies set");
    return FAIL;
  }

#ifdef TRANSFER
  if (RadiativeTransfer) {
    if (MultiSpecies) {
      FieldType[kphHINum    = NumberOfBaryonFields++] = kphHI;
      FieldType[gammaNum    = NumberOfBaryonFields++] = PhotoGamma;
      FieldType[kphHeINum   = NumberOfBaryonFields++] = kphHeI;
      FieldType[kphHeIINum  = NumberOfBaryonFields++] = kphHeII;
      if (MultiSpecies > 1) {
        FieldType[kdissH2INum    = NumberOfBaryonFields++] = kdissH2I;
      }
    }
    
    if (RadiationPressure) {
      FieldType[RPresNum1 = NumberOfBaryonFields++] = RadPressure0;
      FieldType[RPresNum2 = NumberOfBaryonFields++] = RadPressure1;
      FieldType[RPresNum3 = NumberOfBaryonFields++] = RadPressure2;
    }
    
    NumberOfPhotonPackages = 0;
    PhotonPackages-> NextPackage= NULL;
    
  }
#endif

  int idrivex, idrivey, idrivez;
  if (UseDrivingField) {
    idrivex = NumberOfBaryonFields;
    idrivey = idrivex + 1;
    idrivez = idrivex + 2;
    FieldType[NumberOfBaryonFields++] = DrivingField1;
    FieldType[NumberOfBaryonFields++] = DrivingField2;
    FieldType[NumberOfBaryonFields++] = DrivingField3;
  }

  if (WritePotential) {
    FieldType[NumberOfBaryonFields++] = GravPotential;
    FieldType[NumberOfBaryonFields++] = AccelerationField1;
    FieldType[NumberOfBaryonFields++] = AccelerationField2;
    FieldType[NumberOfBaryonFields++] = AccelerationField3;
  }

  //FieldType[GMC1Num=NumberOfBaryonFields++] = GMC1; 
  //FieldType[GMC2Num=NumberOfBaryonFields++] = GMC2;

  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1.0, TimeUnits = 1.0,
    VelocityUnits = 1.0;
  if (UsePhysicalUnit)
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, Time);
  double MassUnits = DensityUnits*pow(LengthUnits,3);
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    printf("Mass Units = %"GSYM" \n",MassUnits);
    printf("Time Units = %"GSYM" \n",TimeUnits);
    printf("Density Units = %"GSYM" \n",DensityUnits);
  }

  GravitationalConstant = 4.0*pi*GravConst*MassUnits*pow(TimeUnits,2)/pow(LengthUnits,3);

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  if (SetBaryonFields == 0) 
    return SUCCESS;


  size = 1;
  for (dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  AllocateGrids();

  /* Initialize radiation fields */
#ifdef TRANSFER
  if (this->InitializeRadiativeTransferFields() == FAIL) {
    fprintf(stderr, "\nError in InitializeRadiativeTransferFields.\n");
    return FAIL;
  }
#endif
  
  float *TurbulenceVelocity[3], *DrivingField[3];
  float CloudInternalEnergy, CloudPressure;
  activesize = 1;
  for (dim = 0; dim < GridRank; dim++) {
    activesize *= (GridDimension[dim] - 2*NumberOfGhostZones);
  }

  for (dim = 0; dim < GridRank; dim++) {
    TurbulenceVelocity[dim] = new float[size];
    DrivingField[dim] = new float[size];
    for (n = 0; n < activesize; n++) {
      TurbulenceVelocity[dim][n] = 0.0;
      DrivingField[dim][n] = 0.0;
    }
  }


  /* Set density, rotational speed, internal energy and chemical species. */

  /* assume isothermal initially */

  CloudPressure = CloudDensity * pow(CloudSoundSpeed,2);
  //  CloudInternalEnergy = pow(CloudSoundSpeed,2) / (Gamma - 1.0);
  float h, cs, dpdrho, dpde;
  EOS(CloudPressure, CloudDensity, CloudInternalEnergy, h, cs, dpdrho, dpde, EOSType, 1);

  float InitialFractionHII = 1.2e-5;
  float InitialFractionHeII = 1.0e-14;
  float InitialFractionHeIII = 1.0e-17;  

  /* Cloud center is box center. */
  float Cloud2Radius = CloudRadius*0.5; // r1 = 1/2*r2
  FLOAT xc = 0.5, yc = 0.5, zc = 0.5, xpos, ypos, zpos,
    x, y, z, r;
    //cosphi, sinphi, x, y, z, r;
  FLOAT r2d, r1, r2;
  FLOAT vcol, Bparam;
  FLOAT tcross = sqrt(1.0/(GravConst*(4.0/3.0)*pi*CloudDensity*DensityUnits))/TimeUnits;
  FLOAT Dsep, x1c, y1c, z1c, x2c, y2c, z2c;

  /* test: colliding, small, identical, no separation */
  if (CloudType >= 201 && CloudType <= 209) { 
      if (CloudType == 201) {
          Bparam=0.0;
      } else if (CloudType == 202) {
          Bparam=0.5;
      }
      //vcol = 5.75e5/VelocityUnits; //mach 25 collision test for small clouds
      vcol = 10e5/VelocityUnits; //10 km/s for large clouds (for filament_paper nonturbulent)
      Dsep = 0.0;
      x1c = 0.5-(CloudRadius), y1c = 0.5, z1c = 0.5, x2c = 0.5+(CloudRadius), y2c = 0.5+(Bparam*CloudRadius), z2c = 0.5;
      } else {     // with separation, different clouds
      vcol = 10e5/VelocityUnits; // v_rel=10km/s fiducial collision
      //Dsep = vcol*tcross;
      //x1c = 0.5-(CloudRadius+Dsep/2.0), y1c = 0.5, z1c = 0.5, x2c = 0.5+(Cloud2Radius+Dsep/2.0), y2c = 0.5, z2c = 0.5;
      }
  /* colliding, big, identical, no separation */
  if (CloudType >= 211 && CloudType <= 219) {
      Bparam = ImpactParameter;    //  GMC collision parameters
      vcol = RelativeVelocity;     //  GMC collision parameters
      Dsep = 0.0;
      x1c = 0.5-(CloudRadius), y1c = 0.5, z1c = 0.5, x2c = 0.5+(CloudRadius), y2c = 0.5+(Bparam*CloudRadius), z2c = 0.5;
      }
  /* same as before but cloud-only velocities */
  if (CloudType >= 221 && CloudType <= 229) { 
      Bparam = ImpactParameter;    //  GMC collision parameters
      vcol = RelativeVelocity;     //  GMC collision parameters
      Dsep = 0.0;
      x1c = 0.5-(CloudRadius), y1c = 0.5, z1c = 0.5, x2c = 0.5+(CloudRadius), y2c = 0.5+(Bparam*CloudRadius), z2c = 0.5;
      }


  float Density, eint, Velx, Vely, Velz;
  //float GMC1Val, GMC2Val;
  //float GMC1, GMC2;
  /* angled B-field */
  float theta = 90.0, phi = Btheta;   // GMC parameters
  float costheta = cos(theta*M_PI/180.0), sintheta = sin(theta*M_PI/180);
  float cosphi = cos(phi*M_PI/180.0), sinphi = sin(phi*M_PI/180);


  n = 0;
  for (k = 0; k < GridDimension[2]; k++) {
    for (j = 0; j < GridDimension[1]; j++) {
      for (i = 0; i < GridDimension[0]; i++, n++) {

        x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
        y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
        z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
        r = sqrt(pow(x-xc,2) + pow(y-yc,2) + pow(z-zc,2));
        r = max(r, CellWidth[0][0]);

        r2d = sqrt(pow(x-xc,2) + pow(y-yc,2));
        r2d = max(r2d, CellWidth[0][0]);

        r1 = sqrt(pow(x-x1c,2) + pow(y-y1c,2) + pow(z-z1c,2));
        r1 = max(r1, CellWidth[0][0]);
        r2 = sqrt(pow(x-x2c,2) + pow(y-y2c,2) + pow(z-z2c,2));
        r2 = max(r2, CellWidth[0][0]);

        xpos = x - xc;
        ypos = y - yc;
        zpos = z - zc;

        /* compute the azimuthal angle */
        //cosphi = xpos/sqrt(xpos*xpos+ypos*ypos);
        //sinphi = ypos/sqrt(xpos*xpos+ypos*ypos);

        Velx = Vely = Velz = 0.0;

        // default values:
        Density = CloudDensity;
        eint = CloudInternalEnergy;

        /* CloudType >=100 cylinder */

        //if (r < CloudRadius) {
        if (r < CloudRadius && CloudType < 100) {

          /* Type 0: uniform cloud, 7: uniform cloud (only turb k=1-2) */

          if (CloudType == 0 || CloudType == 7) {
            Density = CloudDensity;
            eint = CloudInternalEnergy;
          }

          /* Type 1: cloud density profile is the same as that of Nakamura & Li, ApJ, 662, 395 */

          if (CloudType == 1) {
            Density = CloudDensity / (1.0 + pow(3.0*r/CloudRadius,2));
            eint = CloudInternalEnergy;
            Velx = -CloudAngularVelocity * ypos;
            Vely =  CloudAngularVelocity * xpos;
            /*else {
              Velx = -CloudAngularVelocity * (CloudRadius-r) * 1.5 / CloudRadius * ypos;
              Vely =  CloudAngularVelocity * (CloudRadius-r) * 1.5 / CloudRadius * xpos;
              }*/
          }

          /* Type 2: cloud density is uniform sphere embeded in low density medium */

          if (CloudType == 2) {
            Density = CloudDensity;
            eint = CloudInternalEnergy;
            Velx = -CloudAngularVelocity * ypos;
            Vely =  CloudAngularVelocity * xpos;
          }

          /* Type 3: flattened 1/r^2 profile without ambient medium. */

          if (CloudType == 3) {
            Density = CloudDensity / (1.0 + pow(6.0*r/CloudRadius,2));
            eint = CloudInternalEnergy;
          }
 
          /* Type 4: 1/r^2 profile with a smaller core.
             This is a model for massive star formation with a seed
             protostar in the center */

          if (CloudType == 4) {
            Density = 4.2508525*CloudDensity / (1.0 + pow(9.0*r/CloudRadius,2));
            eint = CloudInternalEnergy;
            Velx = -CloudAngularVelocity * ypos;
            Vely =  CloudAngularVelocity * xpos;
          }

          if (CloudType == 5) {
            float drho = rand();
            drho = 0.0;//drho/RAND_MAX;
            if (r < 5) {
              Density = CloudDensity*(1.0+0.1*drho);
            } else {
              Density = CloudDensity/pow(r/0.05, 2)*(1.0+0.1*drho);
            }
            eint = CloudInternalEnergy/(1.0+0.1*drho);
          }

          /* Type 6: flattened 1/r^2 profile with large core and with ambient medium. */

          if (CloudType == 6) {
            Density = 1.0522054*CloudDensity / (1.0 + pow(4.0*r/CloudRadius,2));
            eint = CloudInternalEnergy;
          }

          /* Type 8: cooling test problem, sphere w high density range (not used) */

          if (CloudType == 8) {
            Density = CloudDensity / (0.001 + pow(99*r/CloudRadius,2));
            eint = CloudInternalEnergy;
          }


        } else if (r2d < CloudRadius && CloudType == 100) {
            Density = CloudDensity;
            eint = CloudInternalEnergy;

        } else if ((r1 < CloudRadius || r2 < Cloud2Radius) && CloudType == 200) {
            Density = CloudDensity;
            eint = CloudInternalEnergy;

        /* two identical clouds */
        } else if ((r1 < CloudRadius) && (CloudType >= 201 && CloudType <= 219)) { 
            Density = CloudDensity;
            eint = CloudInternalEnergy;
            //GMC1 = 1.0;  //
            //GMC2 = 0.01;  //

        } else if ((r2 < CloudRadius) && (CloudType >= 201 && CloudType <= 219)) { 
            Density = CloudDensity;
            eint = CloudInternalEnergy;
            //GMC1 = 0.01;  //
            //GMC2 = 1.0;  //

        } else if (r1 < CloudRadius && CloudType == 300 ) {
            Density = CloudDensity;
            eint = CloudInternalEnergy;

        /* centered isolated cloud */
        } else if (r < CloudRadius && (CloudType >= 301 && CloudType <=309)) { 
            Density = CloudDensity;
            eint = CloudInternalEnergy;

        } else {

          if (CloudType == 0 || CloudType == 7 ) {
            Density = CloudDensity/100.0;
            eint = CloudInternalEnergy*100.;
          }

          if (CloudType == 1) {
            Density = CloudDensity/20.0;
            eint = CloudInternalEnergy;
          }

          if (CloudType == 2) {
            Density = CloudDensity / 10.0;
            eint = CloudInternalEnergy * 10.0;
          }

          if (CloudType == 3) {
            Density = max(DensityUnits, CloudDensity/(1.0 + pow(6.0*r/CloudRadius,2)));
            eint = CloudInternalEnergy;
          }

          if (CloudType == 4) {
            Density = max(DensityUnits,0.5*4.25*CloudDensity/(1.0 + pow(9.0*r/CloudRadius,2)));
            eint = CloudInternalEnergy*200.0; //400.0;
          }


          if (CloudType == 6) {
            //Density = max(DensityUnits, 0.5*CloudDensity/(1.0 + pow(4.0*r/CloudRadius,2)));
            Density = 0.1*CloudDensity/(1.0 + pow(4.0,2));
            eint = CloudInternalEnergy*100.0; //400.0;
          }

          if (CloudType == 8) {//
            Density = CloudDensity/1e4;
            eint = CloudInternalEnergy*1e4;
          }


          /* cylinder */
          if (CloudType == 100) {
            Density = CloudDensity/10.0;
            eint = CloudInternalEnergy*10.;
          }

          /* two clouds */
          if (CloudType >= 200 && CloudType <= 299 ) {
            Density = CloudDensity/10.0;
            eint = CloudInternalEnergy*10.;
            //GMC1 = 0.01;  //
            //GMC2 = 0.01;  //
          }

          /* isolated cloud */
          if (CloudType >= 300) {
            Density = CloudDensity/10.0;
            eint = CloudInternalEnergy*10.;
          }

        }

        BaryonField[iden ][n] = Density;
        //        BaryonField[ColourNum][n] = Density*0.018477;
        BaryonField[ivx  ][n] = Velx;
        BaryonField[ivy  ][n] = Vely;
        BaryonField[ivz  ][n] = Velz;
        BaryonField[ietot][n] = eint;
        if (HydroMethod != Zeus_Hydro)
          BaryonField[ietot][n] += 0.5*(Velx*Velx + Vely*Vely + Velz*Velz);
        if (DualEnergyFormalism) {
          BaryonField[ieint][n] = eint;
        }
        //BaryonField[GMC1Num][n] = GMC1; //
        //BaryonField[GMC2Num][n] = GMC2; //


        if(Velx != 0.0) 
          printf("    PROBLEM!!!! eint = %g, Velx = %g, Vely = %g, Velz = %g \n", eint, Velx, Vely, Velz);

        if (UseMHD) {
          BaryonField[iBx  ][n]  = InitialBField*sintheta*cosphi;
          BaryonField[iBy  ][n]  = InitialBField*sintheta*sinphi;
          BaryonField[iBz  ][n]  = InitialBField*costheta;
          BaryonField[ietot][n] += 0.5 * pow(InitialBField,2) / Density;
        }
        if( HydroMethod == MHD_RK ){
          BaryonField[iPhi ][n]  = 0.0;
        }
        if ( UseMHDCT ){
          MagneticField[0][n] = 0.0;
          BaryonField[iBy][n] = 0.0;
          BaryonField[iBz][n] = InitialBField;
        }

        if (UseDrivingField) {
          BaryonField[idrivex][n] = 0.0;
          BaryonField[idrivey][n] = 0.0;
          BaryonField[idrivez][n] = 0.0;
        }
 #ifdef USE_KROME
        if (MultiSpecies == KROMESPECIES) {
          for (int abNum=0; abNum<NKROMESPECIES; abNum++) {
            int speciesNum = speciesMap[abNum];
            if (speciesNum != -1) {
              BaryonField[speciesNum][i] = 1.0e-20 //* PrestellarCoreMoleMass[abNum] 
                                         * BaryonField[0][i] / 1.40045;
            }
          }
          // for (int speciesNum = DeNum; speciesNum <= D3OIINum; speciesNum ++) {
          //   BaryonField[speciesNum][i] = 1e-20*BaryonField[0][i];
          // }
          if (1){
            /* set your preferable initial abundances */
            BaryonField[H2INum][i]          = 5.00e-1*BaryonField[0][i] / 1.412;
            BaryonField[HINum][i]           = 5.00e-1*BaryonField[0][i] / 1.412;
            BaryonField[HeINum][i]          = 4.00e-1*BaryonField[0][i] / 1.412;
            BaryonField[NINum][i]           = 14.0*6.1e-5*BaryonField[0][i] / 1.412;
            BaryonField[OINum][i]           = 16.0*4.6e-4*BaryonField[0][i] / 1.412;
            BaryonField[CIINum][i]          = 12.0*2.6e-4*BaryonField[0][i] / 1.412;
            BaryonField[MGINum][i]          = 24.0*3.981e-5*BaryonField[0][i] / 1.412;
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

      }
    }
  }

  /* Initialize turbulent velocity field */

  if (SetTurbulence) {

    float k1, k2, dk;
    if (CloudType == 0) {
      k1 = 2.0;
      k2 = 10.0;
      dk = 1.0;
    }
    if (CloudType == 1) {
      k1 = 2, k2 = 4, dk = 1;
      k1 = int(1.0/CloudRadius);
      dk = int(0.5/CloudRadius);
      k2 = k1 + 8.0*dk;
    }
    if (CloudType == 2) {
      k1 = 2.;
      k2 = 37;
      dk = 5;
    }
    if (CloudType == 3) {
      k1 = 2.0;
      k2 = 10.0;
      dk = 1.0;
    }
    if (CloudType == 4 || CloudType == 6) {
      k1 = 2.0;
      k2 = min(34.0, int(GridDimension[0]/10));
      printf("                GridDimension[0] = %"ISYM"\n",GridDimension[0] );
      dk = max(1.0,int((k2-k1)/10));
    }
    if (CloudType == 7) {
      k1 = 1.0;
      k2 = 2.0;
      dk = 0.5;
    }

    if (CloudType == 8) { //
      k1 = 0.0;
      k2 = 0.0;
      dk = 0.0;
    }


    if (CloudType == 100) {
      k1 = 2.0;
      k2 = 10.0;
      dk = 1.0;
    }

    if (CloudType == 200) {
      k1 = 10.0;
      k2 = 30.0;
      dk = 1.0;
    }

    if (CloudType == 201) {
      k1 = 0.0;
      k2 = 0.0;
      dk = 0.0;
    }

    if (CloudType == 300) {
      k1 = 10.0;
      k2 = 30.0;
      dk = 1.0;
    }

    if (CloudType == 301 || CloudType == 211) { // 
      k1 = 2.0;
      k2 = 20.0;
      dk = 1.0;
    }

    if (CloudType == 302 || CloudType == 212) {
      k1 = 20.0;
      k2 = 40.0;
      dk = 1.0;
    }

    if (CloudType == 303 || CloudType == 213) {
      k1 = 4.0;
      k2 = 20.0;
      dk = 1.0;
    }

    if (CloudType == 304 || CloudType == 214) {
      k1 = 4.0;
      k2 = 30.0;
      dk = 1.0;
    }



    printf("Begin generating turbulent velocity spectrum...\n");
    Turbulence_Generator(TurbulenceVelocity, GridDimension[0], 
                         GridDimension[1],
                         GridDimension[2],
                         4.0, k1, k2, dk,
                         CellLeftEdge, CellWidth, TurbulenceSeed);    
    /*Turbulence_Generator(TurbulenceVelocity, GridDimension[0]-2*NumberOfGhostZones, 
                         GridDimension[1]-2*NumberOfGhostZones,
                         GridDimension[2]-2*NumberOfGhostZones,
                         4.0, k1, k2, dk,
                         CellLeftEdge, CellWidth, TurbulenceSeed);    */

    printf("Turbulent spectrum generated\n");

    float VelocityNormalization = 1;
    // for level > 0 grids the CloudMachNumber passed in is actuall the Velocity normalization factor
    if (level > 0) VelocityNormalization = CloudMachNumber; 
    printf("Cloud Mach Number = %"GSYM" \n",CloudMachNumber);
    // for (i = 0; i < 3; i++) {
    //   for (n = 0; n < activesize; n++) {
    //     TurbulenceVelocity[i][n] *= VelocityNormalization;
    //   }
    // }


    /* Set turbulent velocity field */

    if (level == 0) VelocityNormalization = 1;

    float velfact = 0;

    n = 0;
    for (k = 0; k < GridDimension[2]; k++) {
      for (j = 0; j < GridDimension[1]; j++) {
        for (i = 0; i < GridDimension[0]; i++, n++) {
    /*for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
        for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {*/
          igrid = i + GridDimension[0]*(j+k*GridDimension[1]);
          x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
          y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
          z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
          
          r = sqrt(pow(fabs(x-xc),2)+pow(fabs(y-yc),2)+pow(fabs(z-zc),2));
          r = max(r, 0.1*CellWidth[0][0]);

          r2d = sqrt(pow(fabs(x-xc),2)+pow(fabs(y-yc),2));
          r2d = max(r2d, 0.1*CellWidth[0][0]);

          r1 = sqrt(pow(x-x1c,2) + pow(y-y1c,2) + pow(z-z1c,2));
          r1 = max(r1, CellWidth[0][0]);
          r2 = sqrt(pow(x-x2c,2) + pow(y-y2c,2) + pow(z-z2c,2));
          r2 = max(r2, CellWidth[0][0]);

          if (CloudType < 100) {
            if (r < CloudRadius) {

              TurbulenceVelocity[0][n] *= VelocityNormalization;
              TurbulenceVelocity[1][n] *= VelocityNormalization;
              TurbulenceVelocity[2][n] *= VelocityNormalization;

              BaryonField[ivx][igrid] += TurbulenceVelocity[0][n];
              BaryonField[ivy][igrid] += TurbulenceVelocity[1][n];
              BaryonField[ivz][igrid] += TurbulenceVelocity[2][n];
              BaryonField[ietot][igrid] +=
                0.5 * (pow(TurbulenceVelocity[0][n],2) +
                       pow(TurbulenceVelocity[1][n],2) +
                       pow(TurbulenceVelocity[2][n],2));
            }
          }

          if (CloudType == 100) {
            if (r2d < CloudRadius) {
              TurbulenceVelocity[0][n] *= VelocityNormalization;
              TurbulenceVelocity[1][n] *= VelocityNormalization;
              TurbulenceVelocity[2][n] *= VelocityNormalization;

              BaryonField[ivx][igrid] += TurbulenceVelocity[0][n];
              BaryonField[ivy][igrid] += TurbulenceVelocity[1][n];
              BaryonField[ivz][igrid] += TurbulenceVelocity[2][n];
              BaryonField[ietot][igrid] +=
                0.5 * (pow(TurbulenceVelocity[0][n],2) +
                       pow(TurbulenceVelocity[1][n],2) +
                       pow(TurbulenceVelocity[2][n],2));
            }
          }

          if (CloudType == 200) {
            if (r1 < CloudRadius) {
              TurbulenceVelocity[0][n] *= VelocityNormalization;
              TurbulenceVelocity[1][n] *= VelocityNormalization;
              TurbulenceVelocity[2][n] *= VelocityNormalization;

              BaryonField[ivx][igrid] += TurbulenceVelocity[0][n];
              BaryonField[ivy][igrid] += TurbulenceVelocity[1][n];
              BaryonField[ivz][igrid] += TurbulenceVelocity[2][n];
            }

            if (r2 < Cloud2Radius) {
              TurbulenceVelocity[0][n] *= VelocityNormalization*Cloud2Radius/CloudRadius;
              TurbulenceVelocity[1][n] *= VelocityNormalization*Cloud2Radius/CloudRadius;
              TurbulenceVelocity[2][n] *= VelocityNormalization*Cloud2Radius/CloudRadius;

              BaryonField[ivx][igrid] += TurbulenceVelocity[0][n];
              BaryonField[ivy][igrid] += TurbulenceVelocity[1][n];
              BaryonField[ivz][igrid] += TurbulenceVelocity[2][n];
            }
              /* velocity scale: 1=2km/s (prev:velfact=2.0) */
              if (level > 0 && x < 0.5) {
                velfact = 2.5;
              } else if (level > 0 && x >= 0.5) {
                velfact = -2.5;
              } else {
                velfact = 0;
              }
              BaryonField[ivx][igrid] += velfact;

              BaryonField[ietot][igrid] +=
                0.5 * (pow(BaryonField[ivx][igrid],2) +
                       pow(BaryonField[ivy][igrid],2) +
                       pow(BaryonField[ivz][igrid],2));
          }

          if (CloudType >= 211 && CloudType <= 219) {     // large colliding identical clouds
            if (r1 < CloudRadius || r2 < CloudRadius) {
              TurbulenceVelocity[0][n] *= VelocityNormalization;
              TurbulenceVelocity[1][n] *= VelocityNormalization;
              TurbulenceVelocity[2][n] *= VelocityNormalization;

              BaryonField[ivx][igrid] += TurbulenceVelocity[0][n];
              BaryonField[ivy][igrid] += TurbulenceVelocity[1][n];
              BaryonField[ivz][igrid] += TurbulenceVelocity[2][n];
            }
              /* velocity scale: 1=2km/s */
              //if (x < 0.5) {
              if (level > 0 && x < 0.5) {
                velfact = 0.5*vcol; 
              //} else if (x >= 0.5) {
              } else if (level > 0 && x >= 0.5) {
                velfact = -0.5*vcol;
              //} else {
                //velfact = 0;
              }
              BaryonField[ivx][igrid] += velfact; // (adds bulk flow; also changes in Grid_NormalizeVelocity.C)

              BaryonField[ietot][igrid] +=
                0.5 * (pow(BaryonField[ivx][igrid],2) +
                       pow(BaryonField[ivy][igrid],2) +
                       pow(BaryonField[ivz][igrid],2));
          }

          if (CloudType >= 221 && CloudType <= 229) {     // large colliding identical clouds, vel in cloud only
            if (r1 < CloudRadius || r2 < CloudRadius) {
              TurbulenceVelocity[0][n] *= VelocityNormalization;
              TurbulenceVelocity[1][n] *= VelocityNormalization;
              TurbulenceVelocity[2][n] *= VelocityNormalization;

              BaryonField[ivx][igrid] += TurbulenceVelocity[0][n];
              BaryonField[ivy][igrid] += TurbulenceVelocity[1][n];
              BaryonField[ivz][igrid] += TurbulenceVelocity[2][n];
            }
              /* velocity scale: 1=2km/s */
              if (level > 0 && r1 < CloudRadius ) {
                velfact = 0.5*vcol;
              } else if (level > 0 && r2 < CloudRadius) {
                velfact = -0.5*vcol;
              } else {
                velfact = 0;
              }
              BaryonField[ivx][igrid] += velfact;

              BaryonField[ietot][igrid] +=
                0.5 * (pow(BaryonField[ivx][igrid],2) +
                       pow(BaryonField[ivy][igrid],2) +
                       pow(BaryonField[ivz][igrid],2));
          }


          if (CloudType == 300) {
            if (r1 < CloudRadius) {
              TurbulenceVelocity[0][n] *= VelocityNormalization;
              TurbulenceVelocity[1][n] *= VelocityNormalization;
              TurbulenceVelocity[2][n] *= VelocityNormalization;

              BaryonField[ivx][igrid] += TurbulenceVelocity[0][n];
              BaryonField[ivy][igrid] += TurbulenceVelocity[1][n];
              BaryonField[ivz][igrid] += TurbulenceVelocity[2][n];
            }

              /* velocity scale: 1=2km/s (prev:velfact=2.0) */
              if (level > 0) {
                velfact = -2.5;
              } else {
                velfact = 0;
              }
              BaryonField[ivx][igrid] += velfact;

              BaryonField[ietot][igrid] +=
                0.5 * (pow(BaryonField[ivx][igrid],2) +
                       pow(BaryonField[ivy][igrid],2) +
                       pow(BaryonField[ivz][igrid],2));

          }

          if (CloudType >= 301 && CloudType <= 309) {
            if (r < CloudRadius) {
              TurbulenceVelocity[0][n] *= VelocityNormalization;
              TurbulenceVelocity[1][n] *= VelocityNormalization;
              TurbulenceVelocity[2][n] *= VelocityNormalization;

              BaryonField[ivx][igrid] += TurbulenceVelocity[0][n];
              BaryonField[ivy][igrid] += TurbulenceVelocity[1][n];
              BaryonField[ivz][igrid] += TurbulenceVelocity[2][n];
            }

              /* velocity scale: 1=2km/s (prev:velfact=2.0) */
              velfact = 0;
              BaryonField[ivx][igrid] += velfact;

              BaryonField[ietot][igrid] +=
                0.5 * (pow(BaryonField[ivx][igrid],2) +
                       pow(BaryonField[ivy][igrid],2) +
                       pow(BaryonField[ivz][igrid],2));

          }


          
          //if (r < CloudRadius) {
            //BaryonField[ivx][igrid] += TurbulenceVelocity[0][n];
            //BaryonField[ivy][igrid] += TurbulenceVelocity[1][n];
            //BaryonField[ivz][igrid] += TurbulenceVelocity[2][n];
            //BaryonField[ietot][igrid] += 
              //0.5 * (pow(TurbulenceVelocity[0][n],2) + 
                     //pow(TurbulenceVelocity[1][n],2) + 
                     //pow(TurbulenceVelocity[2][n],2));
          //}
        } 
      }
    }    

  }

  /* Initialize driving force field = efficiency * density * velocity / t_ff*/


  if (UseDrivingField) {
    float k1, k2, dk;
    k1 = 3.0;
    k2 = 4.0;
    dk = 1.0;
    printf("Begin generating driving force field ...\n");
    Turbulence_Generator(DrivingField, GridDimension[0]-2*NumberOfGhostZones, 
                         GridDimension[1]-2*NumberOfGhostZones,
                         GridDimension[2]-2*NumberOfGhostZones,
                         4.0, k1, k2, dk,
                         CellLeftEdge, CellWidth, TurbulenceSeed);
    printf("Driving force field generated\n");


    /* Renormalize to ensure <F>=0 */

    double Fx = 0.0, Fy = 0.0, Fz = 0.0;
    for (n = 0; n < activesize; n++) {
      Fx += DrivingField[0][n];
      Fy += DrivingField[1][n];
      Fz += DrivingField[2][n];
    }

    Fx /= activesize;
    Fy /= activesize;
    Fz /= activesize;
    
    for (n = 0; n < activesize; n++) {
      DrivingField[0][n] -= Fx;
      DrivingField[1][n] -= Fy;
      DrivingField[2][n] -= Fz;
    }
    
    /* Renormalize the mass-weighted 3D rms velocity inside the cloud */
    /*
    double VelRMS = 0.0, Mass = 0.0;
    n = 0;
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
        for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
          igrid = i + GridDimension[0]*(j+k*GridDimension[1]);
          x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
          y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
          z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
          
          r = sqrt(pow(fabs(x-xc),2)+pow(fabs(y-yc),2)+pow(fabs(z-zc),2));
          r = max(r, 0.1*CellWidth[0][0]);
          
          if (r <= CloudRadius) {
            VelRMS += BaryonField[iden][igrid] * 
              sqrt(pow(DrivingField[0][n],2) + pow(DrivingField[1][n],2) + pow(DrivingField[2][n],2));
            Mass += BaryonField[iden][igrid];
          }
          
        }
      }
    }
    */
    //printf("Grid_TubInit: Mass = %"FSYM"\n",Mass);
    /* VelRMS /= Mass;
    double t_ff = sqrt(32.0/(3.0*M_PI*CloudDensity));
    double NormFactor = CloudMachNumber * CloudSoundSpeed / VelRMS / t_ff;
    for (dim = 0; i < GridRank; dim++) {
      for (n = 0; n < activesize; n++) {
        DrivingField[dim][n] *= NormFactor;
      }
      }*/

    n = 0;
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
        for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
          igrid = i + GridDimension[0]*(j+k*GridDimension[1]);
          BaryonField[idrivex][igrid] = DrivingField[0][n];
          BaryonField[idrivey][igrid] = DrivingField[1][n];
          BaryonField[idrivez][igrid] = DrivingField[2][n];
        }
      }
    }

  }

  for (dim = 0; dim < GridRank; dim++) {
    delete [] TurbulenceVelocity[dim];
    delete [] DrivingField[dim];
  }    

  /* Put a sink particle if we are studying massive star formation */

  //  if (PutSink == 1 && level == MaximumRefinementLevel) {
  if (PutSink == 1 && level == 0) {  // set it up on level zero and make it mustrefine

    //    double mass_p = 20.0*1.989e33;
    double mass_p = 3.415*1.989e33;
    mass_p /= MassUnits;
    double dx = CellWidth[0][0];
    double den_p = mass_p / pow(dx,3);
    double t_dyn = sqrt(3*M_PI/(6.672e-8*den_p*DensityUnits));
    t_dyn /= TimeUnits;

    double dxm = dx / pow(RefineBy, MaximumRefinementLevel);

    printf("Adding a Sink Particle. \n");

    NumberOfParticles = 1;
    NumberOfStars = 1;    
    NumberOfParticleAttributes = 3;
    //    MaximumParticleNumber = 1;
    if (StellarWindFeedback) NumberOfParticleAttributes = 6;
    this->AllocateNewParticles(NumberOfParticles);

    ParticleMass[0] = den_p;
    ParticleNumber[0] = 0;
    ParticleType[0] = PARTICLE_TYPE_MUST_REFINE;
    ParticlePosition[0][0] = 0.35+0.5*dxm;
    ParticlePosition[1][0] = 0.45+0.5*dxm;
    ParticlePosition[2][0] = 0.70+0.5*dxm;
    ParticleVelocity[0][0] = 0.0;
    ParticleVelocity[1][0] = 0.0;
    ParticleVelocity[2][0] = 0.0;

    ParticleAttribute[0][0] = 0.0; // creation time             
    ParticleAttribute[1][0] = t_dyn;
    ParticleAttribute[2][0] = mass_p;

    if (StellarWindFeedback) {
      ParticleAttribute[3][0] = 1.0;  
      ParticleAttribute[4][0] = 0.0;
      ParticleAttribute[5][0] = 0.0;
    }

    for (i = 0; i< MAX_DIMENSION+1; i++){
      ParticleAcceleration[i] = NULL;
    }
    this->ClearParticleAccelerations();

  }


  if (PutSink == 2 && level == 0) {  // set it up on level zero and make it mustrefine

    NumberOfParticles = 64;
    NumberOfStars = 64;

    printf("Adding Sink Particles. \n");
    NumberOfParticleAttributes = 3;
    if (StellarWindFeedback) NumberOfParticleAttributes = 6;
    this->AllocateNewParticles(NumberOfParticles);

    //    double mass_p = 20.0*1.989e33;
    double mass_m = 3.415*1.989e33; //Mass of massive stars
    double mass_s = 0.01*1.989e33; //Mass of small stars
    mass_m /= MassUnits;
    mass_s /= MassUnits;
    double dx = CellWidth[0][0];
    double den_m = mass_m / pow(dx,3);
    double den_s = mass_s / pow(dx,3);
    double t_dyn_m = sqrt(3*M_PI/(6.672e-8*den_m*DensityUnits));
    double t_dyn_s = sqrt(3*M_PI/(6.672e-8*den_s*DensityUnits));
    t_dyn_m /= TimeUnits;
    t_dyn_s /= TimeUnits;
    double dxm = dx / pow(RefineBy, MaximumRefinementLevel);

    //    MaximumParticleNumber = 1;


    for (k=0; k<4; k++){
      for (j=0; j<4; j++){
        for (i=0; i<4; i++){
          l = i+4*j+16*k;
          printf("Creating particle %i \n",l);
          ParticleMass[l] = den_m;
          ParticleNumber[l] = l;
          ParticleType[l] = PARTICLE_TYPE_MUST_REFINE;
          ParticlePosition[0][l] = 0.125+0.25*i+0.5*dxm;
          ParticlePosition[1][l] = 0.125+0.25*j+0.5*dxm;
          ParticlePosition[2][l] = 0.125+0.25*k+0.5*dxm;
          ParticleVelocity[0][l] = 0.0;
          ParticleVelocity[1][l] = 0.0;
          ParticleVelocity[2][l] = 0.0;
          ParticleAcceleration[0] = NULL;
          ParticleAcceleration[1] = NULL;
          ParticleAcceleration[2] = NULL;

          ParticleAttribute[0][l] = 0.001; // creation time             
          ParticleAttribute[1][l] = t_dyn_m; // t_dyn
          ParticleAttribute[2][l] = mass_m;

          if (StellarWindFeedback) {
            ParticleAttribute[3][l] = 1.0;  
            ParticleAttribute[4][l] = 0.0;
            ParticleAttribute[5][l] = 0.0;
          }
          /*for (m = 0; m< MAX_DIMENSION+1; m++){
            ParticleAcceleration[m][l] = NULL;
            }*/

          this->ClearParticleAccelerations();
          printf("Completed particle %i, position %g,%g,%g \n",l,ParticlePosition[0][l],ParticlePosition[1][l],ParticlePosition[2][l]);
          //printf("Domain Right Edge = %g %g %g, dx = %g\n",DomainRightEdge[0],DomainRightEdge[1],DomainRightEdge[2],dx);
        }
      }
    }

  }

  if (PutSink == 3 && level == 0) {  // set it up on level zero and make it mustrefine

    NumberOfParticles = 64;
    NumberOfStars = 64;

    printf("Adding Dummy Sink Particles. \n");
    NumberOfParticleAttributes = 3;
    if (StellarWindFeedback) NumberOfParticleAttributes = 6;
    this->AllocateNewParticles(NumberOfParticles);

    //    double mass_p = 20.0*1.989e33;
    double mass_m = 1.989e13; //Mass of massive stars
    mass_m /= MassUnits;
    double dx = CellWidth[0][0];
    double den_m = mass_m / pow(dx,3);
    double t_dyn_m = sqrt(3*M_PI/(6.672e-8*den_m*DensityUnits));
    t_dyn_m /= TimeUnits;
    double dxm = dx / pow(RefineBy, MaximumRefinementLevel);

    for (k=0; k<4; k++){
      for (j=0; j<4; j++){
        for (i=0; i<4; i++){
          l = i+4*j+16*k;
          printf("Creating particle %i \n",l);
          ParticleMass[l] = den_m;
          ParticleNumber[l] = l;
          ParticleType[l] = PARTICLE_TYPE_MUST_REFINE;
          ParticlePosition[0][l] = 1.125+0.25*i;//+0.5*dx;
          ParticlePosition[1][l] = 1.125+0.25*j;//+0.5*dx;
          ParticlePosition[2][l] = 1.125+0.25*k;//+0.5*dx;
          ParticleVelocity[0][l] = 0.0;
          ParticleVelocity[1][l] = 0.0;
          ParticleVelocity[2][l] = 0.0;
          ParticleAcceleration[0] = NULL;
          ParticleAcceleration[1] = NULL;
          ParticleAcceleration[2] = NULL;

          ParticleAttribute[0][l] = 0.0; // creation time             
          ParticleAttribute[1][l] = 0.0; //t_dyn_m; // t_dyn
          ParticleAttribute[2][l] = mass_m;

          if (StellarWindFeedback) {
            ParticleAttribute[3][l] = 1.0;  
            ParticleAttribute[4][l] = 0.0;
            ParticleAttribute[5][l] = 0.0;
          }
          /*for (m = 0; m< MAX_DIMENSION+1; m++){
            ParticleAcceleration[m][l] = NULL;
            }*/

          this->ClearParticleAccelerations();
          printf("Completed particle %i, position %g,%g,%g \n",l,ParticlePosition[0][l],ParticlePosition[1][l],ParticlePosition[2][l]);
          //printf("Domain Right Edge = %g %g %g, dx = %g\n",DomainRightEdge[0],DomainRightEdge[1],DomainRightEdge[2],dx);
        }
      }
    }

  }




  /*  printf("XXX PutSinkParticle = %"ISYM"\n", PutSinkParticle);
  int PutSinkParticle = 0;
  printf("XXX PutSinkParticle = %"ISYM"\n", PutSinkParticle);
  if (PutSinkParticle == 1 && level == 0) {
    NumberOfParticleAttributes = 6;
    double mass_p = 1.1*1.989e33;
    mass_p /= MassUnits;
    double dx = CellWidth[0][0];
    double den_p = mass_p / pow(dx,3);
    double t_dyn = sqrt(3*Pi/(6.672e-8*den_p*DensityUnits));
    t_dyn /= TimeUnits;

    NumberOfParticles = 1;
    NumberOfStarParticles = 1;
    MaximumParticleNumber = 1;
    this->AllocateNewParticles(NumberOfParticles);
    ParticleMass[0] = den_p;
    ParticleNumber[0] = 0;
    ParticleType[0] = PARTICLE_TYPE_MUST_REFINE;
    ParticlePosition[0][0] = 0.50;                         
    ParticlePosition[1][0] = 0.50;
    ParticlePosition[2][0] = 0.50;

    ParticleVelocity[0][0] = 0.0;
    ParticleVelocity[1][0] = 0.0;
    ParticleVelocity[2][0] = 0.0;
    LibbyTouched = 1;
    for (i = 0; i< MAX_DIMENSION+1; i++){
      ParticleAcceleration[i] = NULL;
    }
    this->ClearParticleAccelerations();

    ParticleAttribute[0][0] = 0.0; // creation time    
    ParticleAttribute[1][0] = t_dyn; // dynamical time                                                                
    ParticleAttribute[2][0] = mass_p; //                                                                                 
    printf("XXX Sink Particle in, NumberOfParticles = %"ISYM" \n",NumberOfParticles);
    }*/


  int TestMerge = 0;
  if (TestMerge == 1 && level == 0) {

    double mass_p = 0.05*1.989e33;
    mass_p /= MassUnits;
    double dx = CellWidth[0][0];
    double den_p = mass_p / pow(dx,3);
    double t_dyn = sqrt(3*M_PI/(6.672e-8*den_p*DensityUnits));
    t_dyn /= TimeUnits;

    double dxm = dx / pow(2.0, MaximumRefinementLevel);

    NumberOfParticles = 4;
    NumberOfStars = 4;
    //    MaximumParticleNumber = 4;
    this->AllocateNewParticles(NumberOfParticles);
    for (int i = 0; i < 4; i++) {
      ParticleMass[i] = den_p;
      ParticleNumber[i] = i;
      ParticleType[i] = PARTICLE_TYPE_MUST_REFINE;
      for (int dim = 0; dim < 3; dim++)
        ParticleVelocity[dim][i] = 0.0;
      ParticleAttribute[0][i] = 0.0;
      ParticleAttribute[1][i] = 0.0;
      ParticleAttribute[2][0] = mass_p;
    }
    
    //ParticleMass[0] = 10.0*den_p;

    ParticlePosition[0][0] = 0.501;
    ParticlePosition[1][0] = 0.501;
    ParticlePosition[2][0] = 0.49;

    ParticlePosition[0][1] = 0.501;
    ParticlePosition[1][1] = 0.501;
    ParticlePosition[2][1] = 0.51;

    ParticlePosition[0][2] = 0.501;
    ParticlePosition[1][2] = 0.501;
    ParticlePosition[2][2] = 0.92;

    ParticlePosition[0][3] = 0.501;
    ParticlePosition[1][3] = 0.501;
    ParticlePosition[2][3] = 0.53;


  }

  return SUCCESS;
}

