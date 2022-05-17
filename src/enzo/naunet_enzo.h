#ifndef __NAUNET_ENZO_H__
#define __NAUNET_ENZO_H__
// clang-format off
#include "naunet.h"
#include "naunet_data.h"

#define NAUNET_SPECIES 4
// number of equations to be solved
#define NAUNET_NSPECIES 118
#define NAUNET_NEQUATIONS 114

#define IDX_GCH3OHI 0
#define IDX_GCH4I 1
#define IDX_GCOI 2
#define IDX_GCO2I 3
#define IDX_GH2CNI 4
#define IDX_GH2COI 5
#define IDX_GH2OI 6
#define IDX_GH2SiOI 7
#define IDX_GHCNI 8
#define IDX_GHNCI 9
#define IDX_GHNCOI 10
#define IDX_GHNOI 11
#define IDX_GMgI 12
#define IDX_GN2I 13
#define IDX_GNH3I 14
#define IDX_GNOI 15
#define IDX_GNO2I 16
#define IDX_GO2I 17
#define IDX_GO2HI 18
#define IDX_GSiCI 19
#define IDX_GSiC2I 20
#define IDX_GSiC3I 21
#define IDX_GSiH4I 22
#define IDX_GSiOI 23
#define IDX_CI 24
#define IDX_CII 25
#define IDX_CHI 26
#define IDX_CHII 27
#define IDX_CH2I 28
#define IDX_CH2II 29
#define IDX_CH3I 30
#define IDX_CH3II 31
#define IDX_CH3OHI 32
#define IDX_CH4I 33
#define IDX_CH4II 34
#define IDX_CNI 35
#define IDX_CNII 36
#define IDX_COI 37
#define IDX_COII 38
#define IDX_CO2I 39
#define IDX_EM 40
#define IDX_HI 41
#define IDX_HII 42
#define IDX_H2I 43
#define IDX_H2II 44
#define IDX_H2CNI 45
#define IDX_H2COI 46
#define IDX_H2COII 47
#define IDX_H2NOII 48
#define IDX_H2OI 49
#define IDX_H2OII 50
#define IDX_H2SiOI 51
#define IDX_H3II 52
#define IDX_H3COII 53
#define IDX_H3OII 54
#define IDX_HCNI 55
#define IDX_HCNII 56
#define IDX_HCNHII 57
#define IDX_HCOI 58
#define IDX_HCOII 59
#define IDX_HCO2II 60
#define IDX_HeI 61
#define IDX_HeII 62
#define IDX_HeHII 63
#define IDX_HNCI 64
#define IDX_HNCOI 65
#define IDX_HNOI 66
#define IDX_HNOII 67
#define IDX_HOCII 68
#define IDX_MgI 69
#define IDX_MgII 70
#define IDX_NI 71
#define IDX_NII 72
#define IDX_N2I 73
#define IDX_N2II 74
#define IDX_N2HII 75
#define IDX_NHI 76
#define IDX_NHII 77
#define IDX_NH2I 78
#define IDX_NH2II 79
#define IDX_NH3I 80
#define IDX_NH3II 81
#define IDX_NOI 82
#define IDX_NOII 83
#define IDX_NO2I 84
#define IDX_OI 85
#define IDX_OII 86
#define IDX_O2I 87
#define IDX_O2II 88
#define IDX_O2HI 89
#define IDX_O2HII 90
#define IDX_OCNI 91
#define IDX_OHI 92
#define IDX_OHII 93
#define IDX_SiI 94
#define IDX_SiII 95
#define IDX_SiCI 96
#define IDX_SiCII 97
#define IDX_SiC2I 98
#define IDX_SiC2II 99
#define IDX_SiC3I 100
#define IDX_SiC3II 101
#define IDX_SiHI 102
#define IDX_SiHII 103
#define IDX_SiH2I 104
#define IDX_SiH2II 105
#define IDX_SiH3I 106
#define IDX_SiH3II 107
#define IDX_SiH4I 108
#define IDX_SiH4II 109
#define IDX_SiH5II 110
#define IDX_SiOI 111
#define IDX_SiOII 112
#define IDX_SiOHII 113


#define A_GCH3OHI 32.0
#define A_GCH4I 16.0
#define A_GCOI 28.0
#define A_GCO2I 44.0
#define A_GH2CNI 28.0
#define A_GH2COI 30.0
#define A_GH2OI 18.0
#define A_GH2SiOI 46.0
#define A_GHCNI 27.0
#define A_GHNCI 27.0
#define A_GHNCOI 43.0
#define A_GHNOI 31.0
#define A_GMgI 24.0
#define A_GN2I 28.0
#define A_GNH3I 17.0
#define A_GNOI 30.0
#define A_GNO2I 46.0
#define A_GO2I 32.0
#define A_GO2HI 33.0
#define A_GSiCI 40.0
#define A_GSiC2I 52.0
#define A_GSiC3I 64.0
#define A_GSiH4I 32.0
#define A_GSiOI 44.0
#define A_CI 12.0
#define A_CII 12.0
#define A_CHI 13.0
#define A_CHII 13.0
#define A_CH2I 14.0
#define A_CH2II 14.0
#define A_CH3I 15.0
#define A_CH3II 15.0
#define A_CH3OHI 32.0
#define A_CH4I 16.0
#define A_CH4II 16.0
#define A_CNI 26.0
#define A_CNII 26.0
#define A_COI 28.0
#define A_COII 28.0
#define A_CO2I 44.0
#define A_EM 1.0
#define A_HI 1.0
#define A_HII 1.0
#define A_H2I 2.0
#define A_H2II 2.0
#define A_H2CNI 28.0
#define A_H2COI 30.0
#define A_H2COII 30.0
#define A_H2NOII 32.0
#define A_H2OI 18.0
#define A_H2OII 18.0
#define A_H2SiOI 46.0
#define A_H3II 3.0
#define A_H3COII 31.0
#define A_H3OII 19.0
#define A_HCNI 27.0
#define A_HCNII 27.0
#define A_HCNHII 28.0
#define A_HCOI 29.0
#define A_HCOII 29.0
#define A_HCO2II 45.0
#define A_HeI 4.0
#define A_HeII 4.0
#define A_HeHII 5.0
#define A_HNCI 27.0
#define A_HNCOI 43.0
#define A_HNOI 31.0
#define A_HNOII 31.0
#define A_HOCII 29.0
#define A_MgI 24.0
#define A_MgII 24.0
#define A_NI 14.0
#define A_NII 14.0
#define A_N2I 28.0
#define A_N2II 28.0
#define A_N2HII 29.0
#define A_NHI 15.0
#define A_NHII 15.0
#define A_NH2I 16.0
#define A_NH2II 16.0
#define A_NH3I 17.0
#define A_NH3II 17.0
#define A_NOI 30.0
#define A_NOII 30.0
#define A_NO2I 46.0
#define A_OI 16.0
#define A_OII 16.0
#define A_O2I 32.0
#define A_O2II 32.0
#define A_O2HI 33.0
#define A_O2HII 33.0
#define A_OCNI 42.0
#define A_OHI 17.0
#define A_OHII 17.0
#define A_SiI 28.0
#define A_SiII 28.0
#define A_SiCI 40.0
#define A_SiCII 40.0
#define A_SiC2I 52.0
#define A_SiC2II 52.0
#define A_SiC3I 64.0
#define A_SiC3II 64.0
#define A_SiHI 29.0
#define A_SiHII 29.0
#define A_SiH2I 30.0
#define A_SiH2II 30.0
#define A_SiH3I 31.0
#define A_SiH3II 31.0
#define A_SiH4I 32.0
#define A_SiH4II 32.0
#define A_SiH5II 33.0
#define A_SiOI 44.0
#define A_SiOII 44.0
#define A_SiOHII 45.0


const float A_Table[114] = {
    A_GCH3OHI,
    A_GCH4I,
    A_GCOI,
    A_GCO2I,
    A_GH2CNI,
    A_GH2COI,
    A_GH2OI,
    A_GH2SiOI,
    A_GHCNI,
    A_GHNCI,
    A_GHNCOI,
    A_GHNOI,
    A_GMgI,
    A_GN2I,
    A_GNH3I,
    A_GNOI,
    A_GNO2I,
    A_GO2I,
    A_GO2HI,
    A_GSiCI,
    A_GSiC2I,
    A_GSiC3I,
    A_GSiH4I,
    A_GSiOI,
    A_CI,
    A_CII,
    A_CHI,
    A_CHII,
    A_CH2I,
    A_CH2II,
    A_CH3I,
    A_CH3II,
    A_CH3OHI,
    A_CH4I,
    A_CH4II,
    A_CNI,
    A_CNII,
    A_COI,
    A_COII,
    A_CO2I,
    A_EM,
    A_HI,
    A_HII,
    A_H2I,
    A_H2II,
    A_H2CNI,
    A_H2COI,
    A_H2COII,
    A_H2NOII,
    A_H2OI,
    A_H2OII,
    A_H2SiOI,
    A_H3II,
    A_H3COII,
    A_H3OII,
    A_HCNI,
    A_HCNII,
    A_HCNHII,
    A_HCOI,
    A_HCOII,
    A_HCO2II,
    A_HeI,
    A_HeII,
    A_HeHII,
    A_HNCI,
    A_HNCOI,
    A_HNOI,
    A_HNOII,
    A_HOCII,
    A_MgI,
    A_MgII,
    A_NI,
    A_NII,
    A_N2I,
    A_N2II,
    A_N2HII,
    A_NHI,
    A_NHII,
    A_NH2I,
    A_NH2II,
    A_NH3I,
    A_NH3II,
    A_NOI,
    A_NOII,
    A_NO2I,
    A_OI,
    A_OII,
    A_O2I,
    A_O2II,
    A_O2HI,
    A_O2HII,
    A_OCNI,
    A_OHI,
    A_OHII,
    A_SiI,
    A_SiII,
    A_SiCI,
    A_SiCII,
    A_SiC2I,
    A_SiC2II,
    A_SiC3I,
    A_SiC3II,
    A_SiHI,
    A_SiHII,
    A_SiH2I,
    A_SiH2II,
    A_SiH3I,
    A_SiH3II,
    A_SiH4I,
    A_SiH4II,
    A_SiH5II,
    A_SiOI,
    A_SiOII,
    A_SiOHII
};

#endif