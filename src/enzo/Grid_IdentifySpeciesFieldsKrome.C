/***********************************************************************
/
/  GRID CLASS (IDENTIFY THE MULTI-SPECIES FIELDS)
/
/  written by: KROME developers, based on the original ENZO routine
/  date:       2013
/
************************************************************************/
 
#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
 
 
int grid::IdentifySpeciesFieldsKrome(
 int &DeNum, int &CHINum, int &OINum, int &HNCINum,
 int &HCNINum, int &H2INum, int &CINum, int &HINum,
 int &H2OINum, int &OHINum, int &O2INum, int &CH2INum,
 int &H2COINum, int &HCOINum, int &MGINum,
 int &NH3INum, int &NOINum, int &CNINum, int &COINum,
 int &N2INum, int &NH2INum, int &CH3INum,
 int &CH4INum, int &NINum, int &NHINum, int &HeINum,
 int &HNOINum, int &CH3OHINum, int &CO2INum,
 int &H2CNINum, int &HNCOINum, int &NO2INum,
 int &O2HINum, int &OCNINum, int &CH3OH_DUSTINum,
 int &HNCO_DUSTINum, int &H2CO_DUSTINum, int &CH4_DUSTINum,
 int &CO_DUSTINum, int &H2O_DUSTINum, int &NO_DUSTINum,
 int &CO2_DUSTINum, int &N2_DUSTINum, int &HCN_DUSTINum,
 int &NH3_DUSTINum, int &O2_DUSTINum, int &NO2_DUSTINum,
 int &HNO_DUSTINum, int &O2H_DUSTINum, int &H2CN_DUSTINum,
 int &MG_DUSTINum, int &HNC_DUSTINum, int &E_DUSTINum,
 int &HCOIINum, int &HIINum, int &HOCIINum,
 int &CIINum, int &CH2IINum, int &CHIINum,
 int &H2COIINum, int &MGIINum, int &NH3IINum,
 int &NOIINum, int &CNIINum, int &COIINum,
 int &N2IINum, int &O2IINum, int &H2OIINum,
 int &NH2IINum, int &OIINum, int &OHIINum,
 int &CH3IINum, int &CH4IINum, int &NIINum,
 int &HCNIINum, int &NHIINum, int &H2IINum,
 int &HeIINum, int &HNOIINum, int &H2NOIINum,
 int &H3IINum, int &H3COIINum, int &H3OIINum,
 int &HCNHIINum, int &HCO2IINum, int &HeHIINum,
 int &N2HIINum, int &O2HIINum
){

 DeNum = CHINum = OINum = HNCINum = HCNINum =
 H2INum = CINum = HINum = H2OINum = OHINum =
 O2INum = CH2INum = H2COINum = HCOINum = MGINum =
 NH3INum = NOINum = CNINum = COINum = N2INum =
 NH2INum = CH3INum = CH4INum = NINum = NHINum =
 HeINum = HNOINum = CH3OHINum = CO2INum =
 H2CNINum = HNCOINum = NO2INum = O2HINum =
 OCNINum = CH3OH_DUSTINum = HNCO_DUSTINum =
 H2CO_DUSTINum = CH4_DUSTINum = CO_DUSTINum =
 H2O_DUSTINum = NO_DUSTINum = CO2_DUSTINum =
 N2_DUSTINum = HCN_DUSTINum = NH3_DUSTINum =
 O2_DUSTINum = NO2_DUSTINum = HNO_DUSTINum =
 O2H_DUSTINum = H2CN_DUSTINum = MG_DUSTINum =
 HNC_DUSTINum = E_DUSTINum = HCOIINum = HIINum =
 HOCIINum = CIINum = CH2IINum = CHIINum =
 H2COIINum = MGIINum = NH3IINum = NOIINum =
 CNIINum = COIINum = N2IINum = O2IINum = H2OIINum =
 NH2IINum = OIINum = OHIINum = CH3IINum =
 CH4IINum = NIINum = HCNIINum = NHIINum =
 H2IINum = HeIINum = HNOIINum = H2NOIINum =
 H3IINum = H3COIINum = H3OIINum = HCNHIINum =
 HCO2IINum = HeHIINum = N2HIINum = O2HIINum = 0;
 
  /* Find Fields for the KROME species */
 DeNum = FindField(ElectronDensity, FieldType, NumberOfBaryonFields);
 CHINum = FindField(CHIDensity, FieldType, NumberOfBaryonFields);
 OINum = FindField(OIDensity, FieldType, NumberOfBaryonFields);
 HNCINum = FindField(HNCIDensity, FieldType, NumberOfBaryonFields);
 HCNINum = FindField(HCNIDensity, FieldType, NumberOfBaryonFields);
 H2INum = FindField(H2IDensity, FieldType, NumberOfBaryonFields);
 CINum = FindField(CIDensity, FieldType, NumberOfBaryonFields);
 HINum = FindField(HIDensity, FieldType, NumberOfBaryonFields);
 H2OINum = FindField(H2OIDensity, FieldType, NumberOfBaryonFields);
 OHINum = FindField(OHIDensity, FieldType, NumberOfBaryonFields);
 O2INum = FindField(O2IDensity, FieldType, NumberOfBaryonFields);
 CH2INum = FindField(CH2IDensity, FieldType, NumberOfBaryonFields);
 H2COINum = FindField(H2COIDensity, FieldType, NumberOfBaryonFields);
 HCOINum = FindField(HCOIDensity, FieldType, NumberOfBaryonFields);
 MGINum = FindField(MGIDensity, FieldType, NumberOfBaryonFields);
 NH3INum = FindField(NH3IDensity, FieldType, NumberOfBaryonFields);
 NOINum = FindField(NOIDensity, FieldType, NumberOfBaryonFields);
 CNINum = FindField(CNIDensity, FieldType, NumberOfBaryonFields);
 COINum = FindField(COIDensity, FieldType, NumberOfBaryonFields);
 N2INum = FindField(N2IDensity, FieldType, NumberOfBaryonFields);
 NH2INum = FindField(NH2IDensity, FieldType, NumberOfBaryonFields);
 CH3INum = FindField(CH3IDensity, FieldType, NumberOfBaryonFields);
 CH4INum = FindField(CH4IDensity, FieldType, NumberOfBaryonFields);
 NINum = FindField(NIDensity, FieldType, NumberOfBaryonFields);
 NHINum = FindField(NHIDensity, FieldType, NumberOfBaryonFields);
 HeINum = FindField(HeIDensity, FieldType, NumberOfBaryonFields);
 HNOINum = FindField(HNOIDensity, FieldType, NumberOfBaryonFields);
 CH3OHINum = FindField(CH3OHIDensity, FieldType, NumberOfBaryonFields);
 CO2INum = FindField(CO2IDensity, FieldType, NumberOfBaryonFields);
 H2CNINum = FindField(H2CNIDensity, FieldType, NumberOfBaryonFields);
 HNCOINum = FindField(HNCOIDensity, FieldType, NumberOfBaryonFields);
 NO2INum = FindField(NO2IDensity, FieldType, NumberOfBaryonFields);
 O2HINum = FindField(O2HIDensity, FieldType, NumberOfBaryonFields);
 OCNINum = FindField(OCNIDensity, FieldType, NumberOfBaryonFields);
 CH3OH_DUSTINum = FindField(CH3OH_DUSTIDensity, FieldType, NumberOfBaryonFields);
 HNCO_DUSTINum = FindField(HNCO_DUSTIDensity, FieldType, NumberOfBaryonFields);
 H2CO_DUSTINum = FindField(H2CO_DUSTIDensity, FieldType, NumberOfBaryonFields);
 CH4_DUSTINum = FindField(CH4_DUSTIDensity, FieldType, NumberOfBaryonFields);
 CO_DUSTINum = FindField(CO_DUSTIDensity, FieldType, NumberOfBaryonFields);
 H2O_DUSTINum = FindField(H2O_DUSTIDensity, FieldType, NumberOfBaryonFields);
 NO_DUSTINum = FindField(NO_DUSTIDensity, FieldType, NumberOfBaryonFields);
 CO2_DUSTINum = FindField(CO2_DUSTIDensity, FieldType, NumberOfBaryonFields);
 N2_DUSTINum = FindField(N2_DUSTIDensity, FieldType, NumberOfBaryonFields);
 HCN_DUSTINum = FindField(HCN_DUSTIDensity, FieldType, NumberOfBaryonFields);
 NH3_DUSTINum = FindField(NH3_DUSTIDensity, FieldType, NumberOfBaryonFields);
 O2_DUSTINum = FindField(O2_DUSTIDensity, FieldType, NumberOfBaryonFields);
 NO2_DUSTINum = FindField(NO2_DUSTIDensity, FieldType, NumberOfBaryonFields);
 HNO_DUSTINum = FindField(HNO_DUSTIDensity, FieldType, NumberOfBaryonFields);
 O2H_DUSTINum = FindField(O2H_DUSTIDensity, FieldType, NumberOfBaryonFields);
 H2CN_DUSTINum = FindField(H2CN_DUSTIDensity, FieldType, NumberOfBaryonFields);
 MG_DUSTINum = FindField(MG_DUSTIDensity, FieldType, NumberOfBaryonFields);
 HNC_DUSTINum = FindField(HNC_DUSTIDensity, FieldType, NumberOfBaryonFields);
 E_DUSTINum = FindField(E_DUSTIDensity, FieldType, NumberOfBaryonFields);
 HCOIINum = FindField(HCOIIDensity, FieldType, NumberOfBaryonFields);
 HIINum = FindField(HIIDensity, FieldType, NumberOfBaryonFields);
 HOCIINum = FindField(HOCIIDensity, FieldType, NumberOfBaryonFields);
 CIINum = FindField(CIIDensity, FieldType, NumberOfBaryonFields);
 CH2IINum = FindField(CH2IIDensity, FieldType, NumberOfBaryonFields);
 CHIINum = FindField(CHIIDensity, FieldType, NumberOfBaryonFields);
 H2COIINum = FindField(H2COIIDensity, FieldType, NumberOfBaryonFields);
 MGIINum = FindField(MGIIDensity, FieldType, NumberOfBaryonFields);
 NH3IINum = FindField(NH3IIDensity, FieldType, NumberOfBaryonFields);
 NOIINum = FindField(NOIIDensity, FieldType, NumberOfBaryonFields);
 CNIINum = FindField(CNIIDensity, FieldType, NumberOfBaryonFields);
 COIINum = FindField(COIIDensity, FieldType, NumberOfBaryonFields);
 N2IINum = FindField(N2IIDensity, FieldType, NumberOfBaryonFields);
 O2IINum = FindField(O2IIDensity, FieldType, NumberOfBaryonFields);
 H2OIINum = FindField(H2OIIDensity, FieldType, NumberOfBaryonFields);
 NH2IINum = FindField(NH2IIDensity, FieldType, NumberOfBaryonFields);
 OIINum = FindField(OIIDensity, FieldType, NumberOfBaryonFields);
 OHIINum = FindField(OHIIDensity, FieldType, NumberOfBaryonFields);
 CH3IINum = FindField(CH3IIDensity, FieldType, NumberOfBaryonFields);
 CH4IINum = FindField(CH4IIDensity, FieldType, NumberOfBaryonFields);
 NIINum = FindField(NIIDensity, FieldType, NumberOfBaryonFields);
 HCNIINum = FindField(HCNIIDensity, FieldType, NumberOfBaryonFields);
 NHIINum = FindField(NHIIDensity, FieldType, NumberOfBaryonFields);
 H2IINum = FindField(H2IIDensity, FieldType, NumberOfBaryonFields);
 HeIINum = FindField(HeIIDensity, FieldType, NumberOfBaryonFields);
 HNOIINum = FindField(HNOIIDensity, FieldType, NumberOfBaryonFields);
 H2NOIINum = FindField(H2NOIIDensity, FieldType, NumberOfBaryonFields);
 H3IINum = FindField(H3IIDensity, FieldType, NumberOfBaryonFields);
 H3COIINum = FindField(H3COIIDensity, FieldType, NumberOfBaryonFields);
 H3OIINum = FindField(H3OIIDensity, FieldType, NumberOfBaryonFields);
 HCNHIINum = FindField(HCNHIIDensity, FieldType, NumberOfBaryonFields);
 HCO2IINum = FindField(HCO2IIDensity, FieldType, NumberOfBaryonFields);
 HeHIINum = FindField(HeHIIDensity, FieldType, NumberOfBaryonFields);
 N2HIINum = FindField(N2HIIDensity, FieldType, NumberOfBaryonFields);
 O2HIINum = FindField(O2HIIDensity, FieldType, NumberOfBaryonFields);


  /* Error if any not found. */
  if (DeNum<0 || CHINum<0 || OINum<0 || HNCINum<0 ||
 HCNINum<0 || H2INum<0 || CINum<0 || HINum<0 ||
 H2OINum<0 || OHINum<0 || O2INum<0 || CH2INum<0 ||
 H2COINum<0 || HCOINum<0 || MGINum<0 || NH3INum<0 ||
 NOINum<0 || CNINum<0 || COINum<0 || N2INum<0 ||
 NH2INum<0 || CH3INum<0 || CH4INum<0 || NINum<0 ||
 NHINum<0 || HeINum<0 || HNOINum<0 || CH3OHINum<0 ||
 CO2INum<0 || H2CNINum<0 || HNCOINum<0 ||
 NO2INum<0 || O2HINum<0 || OCNINum<0 || CH3OH_DUSTINum<0 ||
 HNCO_DUSTINum<0 || H2CO_DUSTINum<0 || CH4_DUSTINum<0 ||
 CO_DUSTINum<0 || H2O_DUSTINum<0 || NO_DUSTINum<0 ||
 CO2_DUSTINum<0 || N2_DUSTINum<0 || HCN_DUSTINum<0 ||
 NH3_DUSTINum<0 || O2_DUSTINum<0 || NO2_DUSTINum<0 ||
 HNO_DUSTINum<0 || O2H_DUSTINum<0 || H2CN_DUSTINum<0 ||
 MG_DUSTINum<0 || HNC_DUSTINum<0 || E_DUSTINum<0 ||
 HCOIINum<0 || HIINum<0 || HOCIINum<0 || CIINum<0 ||
 CH2IINum<0 || CHIINum<0 || H2COIINum<0 ||
 MGIINum<0 || NH3IINum<0 || NOIINum<0 || CNIINum<0 ||
 COIINum<0 || N2IINum<0 || O2IINum<0 || H2OIINum<0 ||
 NH2IINum<0 || OIINum<0 || OHIINum<0 || CH3IINum<0 ||
 CH4IINum<0 || NIINum<0 || HCNIINum<0 || NHIINum<0 ||
 H2IINum<0 || HeIINum<0 || HNOIINum<0 || H2NOIINum<0 ||
 H3IINum<0 || H3COIINum<0 || H3OIINum<0 ||
 HCNHIINum<0 || HCO2IINum<0 || HeHIINum<0 ||
 N2HIINum<0 || O2HIINum<0) {
    ENZO_VFAIL("De=%"ISYM", CHI=%"ISYM", OI=%"ISYM", HNCI=%"ISYM", HCNI=%"ISYM", H2I=%"ISYM", CI=%"ISYM", HI=%"ISYM", H2OI=%"ISYM", OHI=%"ISYM", O2I=%"ISYM", CH2I=%"ISYM", H2COI=%"ISYM", HCOI=%"ISYM", MGI=%"ISYM", NH3I=%"ISYM", NOI=%"ISYM", CNI=%"ISYM", COI=%"ISYM", N2I=%"ISYM", NH2I=%"ISYM", CH3I=%"ISYM", CH4I=%"ISYM", NI=%"ISYM", NHI=%"ISYM", HeI=%"ISYM", HNOI=%"ISYM", CH3OHI=%"ISYM", CO2I=%"ISYM", H2CNI=%"ISYM", HNCOI=%"ISYM", NO2I=%"ISYM", O2HI=%"ISYM", OCNI=%"ISYM", CH3OH_DUSTI=%"ISYM", HNCO_DUSTI=%"ISYM", H2CO_DUSTI=%"ISYM", CH4_DUSTI=%"ISYM", CO_DUSTI=%"ISYM", H2O_DUSTI=%"ISYM", NO_DUSTI=%"ISYM", CO2_DUSTI=%"ISYM", N2_DUSTI=%"ISYM", HCN_DUSTI=%"ISYM", NH3_DUSTI=%"ISYM", O2_DUSTI=%"ISYM", NO2_DUSTI=%"ISYM", HNO_DUSTI=%"ISYM", O2H_DUSTI=%"ISYM", H2CN_DUSTI=%"ISYM", MG_DUSTI=%"ISYM", HNC_DUSTI=%"ISYM", E_DUSTI=%"ISYM", HCOII=%"ISYM", HII=%"ISYM", HOCII=%"ISYM", CII=%"ISYM", CH2II=%"ISYM", CHII=%"ISYM", H2COII=%"ISYM", MGII=%"ISYM", NH3II=%"ISYM", NOII=%"ISYM", CNII=%"ISYM", COII=%"ISYM", N2II=%"ISYM", O2II=%"ISYM", H2OII=%"ISYM", NH2II=%"ISYM", OII=%"ISYM", OHII=%"ISYM", CH3II=%"ISYM", CH4II=%"ISYM", NII=%"ISYM", HCNII=%"ISYM", NHII=%"ISYM", H2II=%"ISYM", HeII=%"ISYM", HNOII=%"ISYM", H2NOII=%"ISYM", H3II=%"ISYM", H3COII=%"ISYM", H3OII=%"ISYM", HCNHII=%"ISYM", HCO2II=%"ISYM", HeHII=%"ISYM", N2HII=%"ISYM", O2HII=%"ISYM"\n",
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
 HCO2IINum, HeHIINum, N2HIINum, O2HIINum)
  }
 
  return SUCCESS;
}
