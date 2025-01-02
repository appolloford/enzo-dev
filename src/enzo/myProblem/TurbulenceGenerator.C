/***********************************************************************
/
/  FFTW Turbulence Generator (A copy from Gandalf)
/  ref : https://github.com/gandalfcode/gandalf
/
/  written by: Chia-Jung Hsu
/  date:       May, 2019
/  modified1:  
/
/  PURPOSE: Generate a turbulence by FFTW in 3d cube
/
/  RETURNS: NONE
/
************************************************************************/
#include <iostream>
#ifdef USE_FFTW
#include <fftw3.h>
#endif
#include <cstdlib>
#include <cmath>
#include "phys_constants.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

using namespace std;


//=================================================================================================
//  Class XorshiftRand
/// \brief   'xorshift' random number generator
/// \details 'xorshift' random number generator
///          (See Numerical Recipes 3rd Ed. Chap 3, wikipedia 'Xorshift')
/// \author  D. A. Hubber
/// \date    24/05/2014
//=================================================================================================
class XorshiftRand
{
public:
 //private:

  // Selected full-period triple (ID:A1 from Numerical Recipes)
  // and MLCG modulo 2^64 mapping (ID:D3 from Numerical Recipes)
  static const unsigned_long_int a1 = 21;
  static const unsigned_long_int a2 = 35;
  static const unsigned_long_int a3 = 4;
  static const unsigned_long_int amod = 4768777513237032717;
  static constexpr float invrandmax = 1.0/1.84467440737095e19;


 //public:

  // Internal variables for algorithm
  unsigned_long_int x;

  // Constructor and destructor
  XorshiftRand(unsigned_long_int _seed): x(_seed)
    {
      srand(_seed);
      for (int k=0; k<10; k++) xorshiftrand();
    };
  //~XorshiftRand() {};


  inline unsigned_long_int xorshiftrand(void)
  {
    x ^= x >> a1;
    x ^= x << a2;
    x ^= x >> a3;
    return x*amod;
  }

  inline int intrand(void) {return (int) xorshiftrand();}
  inline long_int longintrand(void) {return (long_int) xorshiftrand();}
  //inline float floatrand(void) {return invrandmax*(float) xorshiftrand();}
  inline float floatrand(void) {return (float) rand()/RAND_MAX;}

  inline float gaussrand(float mean, float sigma)
  {
    float U = 0.0;
    float V = 0.0;
    while (U == 0.0) {
      U = floatrand();
      V = floatrand();
    };
    return sqrt(-2.0*log(U))*cos(2*pi*V);
  }

  void PrintRandomNumberRange(void)
  {
    cout << "Integer range;         min : 1,    max : " << pow(2,64) << endl;
    cout << "Floating point range;  min : 0.0,  max : 1.0" << endl;
    return;
  }

};

// Declare invrandmax constant here (prevents warnings with some compilers)
//const float XorshiftRand::invrandmax = 1.0/1.84467440737095e19;



// XorshiftRand *randnumb = new XorshiftRand(time(NULL));

//=============================================================================
//  DotProduct
//  Calculates the dot product between two vectors, v1 and v2,
//  of given length 'ndim'
//=============================================================================
template <typename T>
static inline T DotProduct(T *v1, T *v2, int ndim)
{
  if (ndim == 1)
    return v1[0]*v2[0];
  else if (ndim == 2)
    return v1[0]*v2[0] + v1[1]*v2[1];
  else
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

void TurbulenceGenerator_FFTW
 (const int field_type,                ///< [in] Type of turbulent velocity field
  const int gridsize,                  ///< [in] Size of velocity grid
  const int randomSeed,
  const int KStart,
  const int KEnd,
  const float power_turb,             ///< [in] Power spectrum index
  float *vfield)                      ///< [out] Array containing velocity field grid
{
#ifdef USE_FFTW
  XorshiftRand *randnumb = new XorshiftRand(randomSeed);

  int ndim=3;
  // Only valid for 3 dimensions
  //-----------------------------------------------------------------------------------------------
  if (ndim == 3) {

    bool divfree;                      // Select div-free turbulence
    bool curlfree;                     // Select curl-free turbulence
    int kmax;                          // Max. extent of k (in 3D)
    int kmin;                          // Min. extent of k (in 3D)
    int krange;                        // Range of k values (kmax - kmin + 1)
    int i,j,k;                         // Grid counters
    int ii,jj,kk;                      // Aux. grid counters
    int kmod;                          // ..
    int kstart;
    int kend;
    int d;                             // Dimension counter
    float F[ndim];                    // Fourier vector component
    float unitk[3];                   // Unit k-vector
    float **power,**phase;            // ..
    float Rnd[3];                     // Random numbers, variable in Gaussian calculation
    //float k_rot[ndim];                // bulk rotation modes
    //float k_com[ndim];                // bulk compression modes
    fftw_plan plan;                    // ??
    fftw_complex *incomplexfield;      // ..
    fftw_complex *outcomplexfield;     // ..

    if(debug)
      printf("TurbulenceGenerator_FFTW: fieldType=%d, gridsize=%d, powerTurb=%lf\n", 
              field_type, gridsize, power_turb);


    divfree = false;
    curlfree = false;
    if (field_type == 1) curlfree = true;
    if (field_type == 2) divfree = true;

    kmin = -(gridsize/2 - 1);
    kmax = gridsize/2;
    krange = kmax - kmin + 1;
//    cout << "kmin : " << kmin << "    kmax : " << kmax << "    krange : "
//         << krange << "    gridsize : " << gridsize << endl;
    if (krange != gridsize) {
      ENZO_FAIL("Error : krange != gridsize in TurbulenceGenerator_FFTW");
    }
    kstart = max(KStart, 0);
    kend = min(KEnd, kmax);
    power = new float*[3];  phase = new float*[3];
    for (d=0; d<3; d++) power[d] = new float[krange*krange*krange];
    for (d=0; d<3; d++) phase[d] = new float[krange*krange*krange];
    incomplexfield = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*gridsize*gridsize*gridsize);
    outcomplexfield = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*gridsize*gridsize*gridsize);

    for (d=0; d<3; d++) {
      for (i=0; i<krange*krange*krange; i++) power[d][i] = 0.0;
      for (i=0; i<krange*krange*krange; i++) phase[d][i] = 0.0;
    }
    for (i=0; i<3*gridsize*gridsize*gridsize; i++) vfield[i] = 0.0;


    // Define wave vectors in Fourier space.  Each wave vector has coordinates in Fourier space,
    // random phases in three dimensions, and power in three dimensions, giving a power and a
    // phase vector field in Fourier space.  With a 'natural' field type there is no coordination
    // between x,y,z components.  With a div-free or curl-free field, there is a relation between
    // coordinates in Fourier space and power vector. For a div-free field, power vectors are
    // perpendicular to the Fourier coordinate vector, and for a curl-free field power vectors
    // are (anti-)parallel to the Fourier coordinate vector
    // (i,j,k) is the Fourier coordinate vector
    // power(1:3,i,j,k) is the power vector at Fourier coordinates i,j,k
    // phase(1:3,i,j,k) is the phase vector at Fourier coordinates i,j,k
    //---------------------------------------------------------------------------------------------
    for (kmod=kstart; kmod<=kend; kmod++) {

      for (i=kmin; i<=kmax; i++) {
        for (j=kmin; j<=kmax; j++) {
          for (k=kmin; k<=kmax; k++) {

            // Skip any k-vectors that have already been calculated
            if (abs(i) != kmod && abs(j) != kmod && abs(k) != kmod) continue;
            if (abs(i) > kmod || abs(j) > kmod || abs(k) > kmod) continue;

            //continue;
            ii = (i + krange)%krange;
            jj = (j + krange)%krange;
            kk = (k + krange)%krange;

            // cycle antiparallel k-vectors
            //if (k < 0) continue;
            //if (k == 0) {
            //  if (j < 0) continue;
            //  if (j == 0 && i < 0) continue;
            //}

            // Central power = 0
            if (i == 0 && j == 0 && k == 0) continue;
            if (i*i + j*j + k*k >= kmax*kmax) continue;

            // Power value, to be multipled by random power chosen from a Gaussian
            // This is what gives the slope of the eventual power spectrum
            for (d=0; d<3; d++) F[d] = sqrt(pow(sqrt((float)(i*i + j*j + k*k)),power_turb));

            for (d=0; d<3; d++) {
              Rnd[0] = randnumb->floatrand();

              // Random phase between 0 and 2*pi (actually -pi and +pi).
              phase[d][ii + krange*jj + krange*krange*kk] = (2.0*Rnd[0] - 1.0)*pi;

              // Create Gaussian distributed random numbers
              Rnd[1] = randnumb->gaussrand(0.0,1.0);
              Rnd[2] = randnumb->gaussrand(0.0,1.0);
              F[d] = Rnd[1]*F[d];
            }

            // Helmholtz decomposition!
            unitk[0] = (float) i;
            unitk[1] = (float) j;
            unitk[2] = (float) k;
            float ksqd = DotProduct(unitk,unitk,3);
            for (d=0; d<3; d++) unitk[d] /= sqrt(ksqd);

            // For curl free turbulence, vector F should be
            // parallel/anti-parallel to vector k
            if (curlfree) {
              for (d=0; d<3; d++) {
                power[d][ii + krange*jj + krange*krange*kk] = unitk[d]*DotProduct(F,unitk,3);
              }
            }

            // For divergence free turbulence, vector F should be perpendicular to vector k
            else if (divfree) {
              for (d=0; d<3; d++) {
                power[d][ii + krange*jj + krange*krange*kk] = F[d] - unitk[d]*DotProduct(F,unitk,3);
              }
            }
            else {
              for (d=0; d<3; d++) power[d][ii + krange*jj + krange*krange*kk] = F[d];
            }

          }
        }
      }

    }
    //---------------------------------------------------------------------------------------------


    float power_spectrum[kmax+1][3];
    for (i=0; i<kmax+1; i++) {
      for (d=0; d<3; d++) power_spectrum[i][d] = 0.0;
    }

    for (i=kmin; i<=kmax; i++) {
      for (j=kmin; j<=kmax; j++) {
        for (k=kmin; k<=kmax; k++) {
          ii = (i + krange)%krange;
          jj = (j + krange)%krange;
          kk = (k + krange)%krange;
          float kmag = sqrt((float)(i*i + j*j + k*k));
          if (kmag >= kmax) kmag = kmax;
          int ibin = (int) (1.0001*kmag);
          power_spectrum[ibin][0] += pow(power[0][ii + krange*jj + krange*krange*kk],2);
          power_spectrum[ibin][1] += pow(power[1][ii + krange*jj + krange*krange*kk],2);
          power_spectrum[ibin][2] += pow(power[2][ii + krange*jj + krange*krange*kk],2);
        }
      }
    }

    plan = fftw_plan_dft_3d(gridsize, gridsize, gridsize, incomplexfield,
                            outcomplexfield, FFTW_BACKWARD, FFTW_ESTIMATE);


    // reorder array: positive wavenumbers are placed in ascending order along
    // first half of dimension, i.e. 0 to k_max, negative wavenumbers are placed
    // along second half of dimension, i.e. -k_min to 1.
    for (d=0; d<3; d++) {

      for (i=kmin; i<=kmax; i++) {
        for (j=kmin; j<=kmax; j++) {
          for (k=kmin; k<=kmax; k++) {
            ii = (i + krange)%krange;
            jj = (j + krange)%krange;
            kk = (k + krange)%krange;
            //ii = i - kmin;
            //jj = j - kmin;
            //kk = k - kmin;
            incomplexfield[ii +krange*jj + krange*krange*kk][0] =
              power[d][ii + krange*jj + krange*krange*kk]*
              cos(phase[d][ii + krange*jj + krange*krange*kk]);
            incomplexfield[ii +krange*jj + krange*krange*kk][1] =
              power[d][ii + krange*jj + krange*krange*kk]*
              sin(phase[d][ii + krange*jj + krange*krange*kk]);

          }
        }
      }

      //fftw_execute_dft(plan, complexfield, complexfield);
      fftw_execute_dft(plan, incomplexfield, outcomplexfield);


      for (i=0; i<krange; i++) {
        for (j=0; j<krange; j++) {
          for (k=0; k<krange; k++) {
            vfield[d + 3*i + 3*krange*j + 3*krange*krange*k] =
              outcomplexfield[i + krange*j + krange*krange*k][0];

          }
        }
      }

    }


    fftw_destroy_plan(plan);

    for (d=0; d<3; d++) delete[] phase[d];
    for (d=0; d<3; d++) delete[] power[d];
    delete[] phase;
    delete[] power;
    fftw_free(outcomplexfield);
    fftw_free(incomplexfield);


  }
  //-----------------------------------------------------------------------------------------------
#endif
  return;
}

