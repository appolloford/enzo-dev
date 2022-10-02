#ifndef __NAUNET_ENZO_H__
#define __NAUNET_ENZO_H__
// clang-format off
#include "naunet.h"
#include "naunet_data.h"

// the MultiSpecies option of ENZO
#define NAUNET_SPECIES 4
// the number of species in ENZO, including the fields which are not present
// in the naunet header, not including electron as it is not renormed with 
// others
#define ENZO_NSPECIES 118


#define A_GH2CNI 28.0
#define A_GHNCI 27.0
#define A_GNO2I 46.0
#define A_GSiOI 44.0
#define A_GCOI 28.0
#define A_GHNCOI 43.0
#define A_GMgI 24.0
#define A_GNOI 30.0
#define A_GO2I 32.0
#define A_GO2HI 33.0
#define A_GSiCI 40.0
#define A_GSiC2I 52.0
#define A_GSiC3I 64.0
#define A_GCH3OHI 32.0
#define A_GCO2I 44.0
#define A_GH2SiOI 46.0
#define A_GHNOI 31.0
#define A_GN2I 28.0
#define A_GH2COI 30.0
#define A_GHCNI 27.0
#define A_GH2OI 18.0
#define A_GNH3I 17.0
#define A_SiC3II 64.0
#define A_H2CNI 28.0
#define A_GCH4I 16.0
#define A_H2NOII 32.0
#define A_H2SiOI 46.0
#define A_HeHII 5.0
#define A_HNCOI 43.0
#define A_HOCII 29.0
#define A_SiC2II 52.0
#define A_GSiH4I 32.0
#define A_SiC2I 52.0
#define A_SiC3I 64.0
#define A_SiH5II 33.0
#define A_SiH4II 32.0
#define A_SiCII 40.0
#define A_O2HI 33.0
#define A_SiCI 40.0
#define A_NO2I 46.0
#define A_SiH3II 31.0
#define A_SiH2II 30.0
#define A_OCNI 42.0
#define A_SiH2I 30.0
#define A_SiOHII 45.0
#define A_SiHII 29.0
#define A_SiH4I 32.0
#define A_SiHI 29.0
#define A_SiH3I 31.0
#define A_SiOII 44.0
#define A_HCO2II 45.0
#define A_HNOI 31.0
#define A_CH3OHI 32.0
#define A_MgI 24.0
#define A_MgII 24.0
#define A_CH4II 16.0
#define A_SiOI 44.0
#define A_CNII 26.0
#define A_HCNHII 28.0
#define A_N2HII 29.0
#define A_O2HII 33.0
#define A_SiII 28.0
#define A_SiI 28.0
#define A_HNCI 27.0
#define A_HNOII 31.0
#define A_N2II 28.0
#define A_H3COII 31.0
#define A_CH4I 16.0
#define A_COII 28.0
#define A_H2II 2.0
#define A_NH3I 17.0
#define A_CH3I 15.0
#define A_CO2I 44.0
#define A_NII 14.0
#define A_OII 16.0
#define A_HCNII 27.0
#define A_NH2II 16.0
#define A_NHII 15.0
#define A_O2II 32.0
#define A_CH3II 15.0
#define A_NH2I 16.0
#define A_CH2II 14.0
#define A_H2OII 18.0
#define A_NH3II 17.0
#define A_NOII 30.0
#define A_H3OII 19.0
#define A_N2I 28.0
#define A_CII 12.0
#define A_HCNI 27.0
#define A_CHII 13.0
#define A_CH2I 14.0
#define A_H2COII 30.0
#define A_NHI 15.0
#define A_OHII 17.0
#define A_CNI 26.0
#define A_H2COI 30.0
#define A_HCOI 29.0
#define A_HeII 4.0
#define A_CHI 13.0
#define A_H3II 3.0
#define A_HeI 4.0
#define A_NOI 30.0
#define A_NI 14.0
#define A_OHI 17.0
#define A_O2I 32.0
#define A_CI 12.0
#define A_HII 1.0
#define A_HCOII 29.0
#define A_H2OI 18.0
#define A_OI 16.0
#define A_EM 1.0
#define A_COI 28.0
#define A_H2I 2.0
#define A_HI 1.0


const float A_Table[NSPECIES] = {
    A_GH2CNI,
    A_GHNCI,
    A_GNO2I,
    A_GSiOI,
    A_GCOI,
    A_GHNCOI,
    A_GMgI,
    A_GNOI,
    A_GO2I,
    A_GO2HI,
    A_GSiCI,
    A_GSiC2I,
    A_GSiC3I,
    A_GCH3OHI,
    A_GCO2I,
    A_GH2SiOI,
    A_GHNOI,
    A_GN2I,
    A_GH2COI,
    A_GHCNI,
    A_GH2OI,
    A_GNH3I,
    A_SiC3II,
    A_H2CNI,
    A_GCH4I,
    A_H2NOII,
    A_H2SiOI,
    A_HeHII,
    A_HNCOI,
    A_HOCII,
    A_SiC2II,
    A_GSiH4I,
    A_SiC2I,
    A_SiC3I,
    A_SiH5II,
    A_SiH4II,
    A_SiCII,
    A_O2HI,
    A_SiCI,
    A_NO2I,
    A_SiH3II,
    A_SiH2II,
    A_OCNI,
    A_SiH2I,
    A_SiOHII,
    A_SiHII,
    A_SiH4I,
    A_SiHI,
    A_SiH3I,
    A_SiOII,
    A_HCO2II,
    A_HNOI,
    A_CH3OHI,
    A_MgI,
    A_MgII,
    A_CH4II,
    A_SiOI,
    A_CNII,
    A_HCNHII,
    A_N2HII,
    A_O2HII,
    A_SiII,
    A_SiI,
    A_HNCI,
    A_HNOII,
    A_N2II,
    A_H3COII,
    A_CH4I,
    A_COII,
    A_H2II,
    A_NH3I,
    A_CH3I,
    A_CO2I,
    A_NII,
    A_OII,
    A_HCNII,
    A_NH2II,
    A_NHII,
    A_O2II,
    A_CH3II,
    A_NH2I,
    A_CH2II,
    A_H2OII,
    A_NH3II,
    A_NOII,
    A_H3OII,
    A_N2I,
    A_CII,
    A_HCNI,
    A_CHII,
    A_CH2I,
    A_H2COII,
    A_NHI,
    A_OHII,
    A_CNI,
    A_H2COI,
    A_HCOI,
    A_HeII,
    A_CHI,
    A_H3II,
    A_HeI,
    A_NOI,
    A_NI,
    A_OHI,
    A_O2I,
    A_CI,
    A_HII,
    A_HCOII,
    A_H2OI,
    A_OI,
    A_EM,
    A_COI,
    A_H2I,
    A_HI
};

#endif