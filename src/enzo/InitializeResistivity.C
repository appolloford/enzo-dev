/***********************************************************************
/
/  INITIALIZE TABULAR RESISTIVITY
/
/  written by: Duncan Christie
/  date:       Spring, 2016
/
/  PURPOSE:  Read in the resistivity from an HDF5 file and store it.
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <math.h>
#include "hdf5.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"


/**************************** Functions Prototypes ******************************/
float GetResistivity(float rho, float B, float T,float Time);
float ResistivityTableLookup(float rho,float T,float B);

int InitializeResistivity()
{
  float Time = 0;
  int q, w,i,j,k,n;
  float64 *temp_data;
  long_int temp_int;
  long_int *temp_int_arr;
  char parameter_name[MAX_LINE_LENGTH];

  if (debug) {
    fprintf(stderr,"Initializing Tabular Resistivity.\n");
    fprintf(stderr,"ResistivityTableFile: %s.\n",ResistivityData.ResistivityTableFile);
  }

  // Read cooling data in from hdf5 file.

  hid_t       file_id, dset_id, attr_id; 
  herr_t      status;
  herr_t      h5_error = -1;
  H5A_info_t  a_info;

  file_id = H5Fopen(ResistivityData.ResistivityTableFile, H5F_ACC_RDONLY, H5P_DEFAULT);

  // Resistivity Table

  dset_id =  H5Dopen(file_id, "/eta_AD");
  if (dset_id == h5_error) {
    fprintf(stderr,"Can't open 'eta_AD' in %s.\n",ResistivityData.ResistivityTableFile);
    return FAIL;
  }

  ResistivityData.ResistivityDataSize = 1;

  // Get the dimensions in each direction
  
  attr_id = H5Aopen_name(dset_id, "nT");
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to open 'nT' attribute in 'eta_AD' dataset.\n");
    return FAIL;
  }
  status = H5Aread(attr_id, HDF5_I8, &temp_int);
  if (status < 0) {
    fprintf(stderr,"Failed to read 'nT' attribute in 'eta_AD' dataset.\n");
    return FAIL;
  }
  
  ResistivityData.nT = temp_int;

  ResistivityData.ResistivityDataSize *= ResistivityData.nT;

  status = H5Aclose(attr_id);
  if (status < 0) {
    fprintf(stderr,"Failed to close 'eta_AD' attribute 'nT'.\n");
    return FAIL;
  }



  attr_id = H5Aopen_name(dset_id, "nD");
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to open 'nD' attribute in 'eta_AD' dataset.\n");
    return FAIL;
  }
  status = H5Aread(attr_id, HDF5_I8, &temp_int);
  if (status < 0) {
    fprintf(stderr,"Failed to read 'nD' attribute in 'eta_AD' dataset.\n");
    return FAIL;
  }
  
  ResistivityData.nD = temp_int;

  ResistivityData.ResistivityDataSize *= ResistivityData.nD;

  status = H5Aclose(attr_id);
  if (status < 0) {
    fprintf(stderr,"Failed to close 'eta_AD' attribute 'nD'.\n");
    return FAIL;
  }


  attr_id = H5Aopen_name(dset_id, "nB");
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to open 'nB' attribute in 'eta_AD' dataset.\n");
    return FAIL;
  }
  status = H5Aread(attr_id, HDF5_I8, &temp_int);
  if (status < 0) {
    fprintf(stderr,"Failed to read 'nB' attribute in 'eta_AD' dataset.\n");
    return FAIL;
  }
  ResistivityData.nB = temp_int;


  ResistivityData.ResistivityDataSize *= ResistivityData.nB;

  status = H5Aclose(attr_id);
  if (status < 0) {
    fprintf(stderr,"Failed to close 'eta_AD' attribute 'nB'.\n");
    return FAIL;
  }

  // Temperature Grid Dimension.

  temp_data = new float64[ResistivityData.nT];
  attr_id = H5Aopen_name(dset_id, "Temperature");
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to open 'Temperature' attribute in 'eta_AD' dataset.\n");
    return FAIL;
  }

  status = H5Aread(attr_id, HDF5_R8,temp_data);
  if (status < 0) {
    fprintf(stderr,"Failed to read 'Temperature' attribute in 'eta_AD' dataset.\n");
    return FAIL;
  }

  ResistivityData.TempParam = temp_data;

  ResistivityData.dlogT = log10(ResistivityData.TempParam[1]/ResistivityData.TempParam[0]);

  if (ADUseLogInterp) {
    for (n=0;n<ResistivityData.nT; n++) {
      ResistivityData.TempParam[n] = log10(ResistivityData.TempParam[n]);
    }
  } 


  status = H5Aclose(attr_id);
  if (status < 0) {
    fprintf(stderr,"Failed to close 'eta_AD' attribute 'Temperature'.\n");
    return FAIL;
  }


  // Density Grid Dimension.
  temp_data = new float64[ResistivityData.nD];
  attr_id = H5Aopen_name(dset_id, "Density");
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to open 'Density' attribute in 'eta_AD' dataset.\n");
    return FAIL;
  }

  status = H5Aread(attr_id, HDF5_R8,temp_data);
  if (status < 0) {
    fprintf(stderr,"Failed to read 'Density' attribute in 'eta_AD' dataset.\n");
    return FAIL;
  }
  
  ResistivityData.DensParam = temp_data;
  ResistivityData.dlogD = log10(ResistivityData.DensParam[1]/ResistivityData.DensParam[0]);

  if (ADUseLogInterp) {
    for (n=0;n<ResistivityData.nD; n++) {
      ResistivityData.DensParam[n] = log10(ResistivityData.DensParam[n]);
    }
  } 


  status = H5Aclose(attr_id);
  if (status < 0) {
    fprintf(stderr,"Failed to close 'eta_AD' attribute 'Density'.\n");
    return FAIL;
  }

  //delete [] temp_data;

  // Magnetic Field Grid Dimension.
  temp_data = new float64[ResistivityData.nB];
  attr_id = H5Aopen_name(dset_id, "BField");
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to open 'BField' attribute in 'eta_AD' dataset.\n");
    return FAIL;
  }


  status = H5Aread(attr_id, HDF5_R8,temp_data);
  if (status < 0) {
    fprintf(stderr,"Failed to read 'BField' attribute in 'eta_AD' dataset.\n");
    return FAIL;
  }

  ResistivityData.MagParam = temp_data;
  ResistivityData.dlogB = log10(ResistivityData.MagParam[1]/ResistivityData.MagParam[0]);

  if (ADUseLogInterp) {
    for (n=0;n<ResistivityData.nB; n++) {
      ResistivityData.MagParam[n] = log10(ResistivityData.MagParam[n]);
    }
  } 

  status = H5Aclose(attr_id);
  if (status < 0) {
    fprintf(stderr,"Failed to close 'eta_AD' attribute 'BField'.\n");
    return FAIL;
  }


  // Load the actual table
  temp_data = new float64[ResistivityData.ResistivityDataSize];
  ResistivityData.ResistivityPerp = new float64[ResistivityData.ResistivityDataSize];

  status = H5Dread(dset_id, HDF5_R8, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_data);
  
  if (status < 0) {
    fprintf(stderr,"Failed to read 'eta_AD' dataset.\n");
    return FAIL;
  }
  
  for (n=0;n<ResistivityData.ResistivityDataSize; n++) {
    ResistivityData.ResistivityPerp[n] = log10(temp_data[n]);
  }
  delete temp_data;



  status = H5Dclose(dset_id);
  if (status < 0) {
    fprintf(stderr,"Failed to close 'eta_AD' dataset.\n");
    return FAIL;
  }

  status = H5Fclose (file_id);

  if (debug) {
    fprintf(stdout,"Exiting Initializeresistivity()\n");
  }

  return SUCCESS;
}
