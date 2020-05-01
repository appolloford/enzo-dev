/***********************************************************************
/
/  Tabular Resistivity Data
/
***********************************************************************/

#define RESISTIVITY_TABLE_MAX_DIMENSION 5

struct ResistivityDataType
{

  // Cooling grid file.
  char * ResistivityTableFile;

  // Rank of Cloudy dataset.
  int ResistivityTableRank;

  // Sizes of each dimension
  int nT,nD,nB;

  float dlogT,dlogD,dlogB;

  // Table parameters
  float64 *TempParam;
  float64 *DensParam;
  float64 *MagParam;

  // Perpendicular Resistivity
  float *ResistivityPerp;

  // Hall Resistivity
  float *ResistivityHall;

  // Parallel Resistivity
  float *ResistivityParallel;

  // Length of 1D flattened Cloudy data
  int ResistivityDataSize;
};
