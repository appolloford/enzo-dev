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
float GetIceYield(int specidx, float velocity) {
  const float icemass = 18.0 * mh; // mass of water ice (target)
  const float icebinding = 0.53 * 1.6e-12;
  const float eff = 0.8;
  const float yieldconst = 8.3e-4;

  float projmass = mh * A_Table[specidx];

  float eta = 4.0 * eff * projmass * icemass / pow(projmass+icemass, 2.0);
  float eps = eta * 0.5 * projmass * velocity * velocity / icebinding;
  float eps0 = max(1.0, 4.0*eta);

  float iceyield = 0.0;
  if (eps > eps0) {
    iceyield = 2.0 * yieldconst * pow(eps-eps0, 2.0) / (1.0 + pow(eps/30.0, 4.0/3.0));
  }

  return iceyield;

}
#endif

int grid::ComputeSputteringRate(float *sputteringrate)
{

#ifdef USE_NAUNET

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

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
 
  int i, j, k, igrid;
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

  int KspNum = FindField(Ksputtering, FieldType, NumberOfBaryonFields);
 
  // Get easy to handle pointers for each variable.
 
  float *density     = BaryonField[DensNum];
  float *totalenergy = BaryonField[TENum];
  float *gasenergy   = BaryonField[GENum];
  float *velocity1   = BaryonField[Vel1Num];
  float *velocity2   = BaryonField[Vel2Num];
  float *velocity3   = BaryonField[Vel3Num];

  float *sputtering  = BaryonField[KspNum];

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

  float NumberDensityUnits = DensityUnits / mh;

  int ixm, ixp, iym, iyp, izm, izp;
  float vxl, vxu, vyl, vyu, vzl, vzu;
  int projidxlist[6] = {IDX_H2I, IDX_HeI, IDX_CI, IDX_OI, IDX_SiI, IDX_COI};
  int projnumlist[6] = {H2INum, HeINum, CINum, OINum, SiINum, COINum};
  float masslist[6] = {2.0, 4.0, 12.0, 16.0, 28.0, 28.0};
  float fluxes[6] = {0.0};
  float velin[6] = {0.0};

  // TODO: comoving
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      igrid = (k * GridDimension[1] + j) * GridDimension[0] + GridStartIndex[0];
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, igrid++) {
        sputteringrate[igrid] = 0.0;
        ixm = igrid - 1;
        ixp = igrid + 1;
        iym = igrid - GridDimension[0];
        iyp = igrid + GridDimension[0];
        izm = igrid - GridDimension[1] * GridDimension[0];
        izp = igrid + GridDimension[1] * GridDimension[0];
        velin[0] = min(velocity1[igrid] - velocity1[ixm], 0.0);
        velin[1] = max(velocity1[igrid] - velocity1[ixp], 0.0);
        velin[2] = min(velocity2[igrid] - velocity2[iym], 0.0);
        velin[3] = max(velocity2[igrid] - velocity2[iyp], 0.0);
        velin[4] = min(velocity3[igrid] - velocity3[izm], 0.0);
        velin[5] = max(velocity3[igrid] - velocity3[izp], 0.0);
        for (int idx = 0; idx < 6; idx ++) {
          int projidx = projidxlist[idx];
          int projnum = projnumlist[idx];
          float mass = masslist[idx];
          float *specdens = BaryonField[projnum];
          fluxes[0] = min(specdens[igrid] * velocity1[igrid] - specdens[ixm] * velocity1[ixm], 0.0);
          fluxes[1] = max(specdens[igrid] * velocity1[igrid] - specdens[ixp] * velocity1[ixp], 0.0);
          fluxes[2] = min(specdens[igrid] * velocity2[igrid] - specdens[iym] * velocity2[iym], 0.0);
          fluxes[3] = max(specdens[igrid] * velocity2[igrid] - specdens[iyp] * velocity2[iyp], 0.0);
          fluxes[4] = min(specdens[igrid] * velocity3[igrid] - specdens[izm] * velocity3[izm], 0.0);
          fluxes[5] = max(specdens[igrid] * velocity3[igrid] - specdens[izp] * velocity3[izp], 0.0);

          // printf("Species density=%13.7e\n", specdens[igrid] * NumberDensityUnits / mass);
          for (int idir = 0; idir < 6; idir ++) {
            float yields = GetIceYield(projidx, velin[idir] * VelocityUnits);
            // One cell count half of the inflow flux in one direction
            float flux = 0.5 * abs( fluxes[idir] * VelocityUnits * NumberDensityUnits / mass );
            sputteringrate[igrid] += yields * flux;
            // printf("projmass=%lf, velocity[%d]=%13.7e, yields[%d]=%13.7e, flux[%d]=%13.7e, sputtering=%13.7e\n",
            //         mass, idir, velin[idir]*VelocityUnits, idir, yields, idir, flux, sputteringrate[igrid]);
          }
        }
        sputtering[igrid] = sputteringrate[igrid];
      }
    }
  }

  delete [] g_grid_dimension;
  delete [] g_grid_start;
  delete [] g_grid_end;

  LCAPERF_STOP("grid_NaunetWrapper");

#endif

  return SUCCESS;
}

