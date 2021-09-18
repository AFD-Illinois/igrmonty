#define NDIM 4
#define NUMIN 1.e8
#define NUMAX 1.e24
#define LNUMIN log(NUMIN)
#define LNUMAX log(NUMAX)
#define DLNU ((LNUMAX-LNUMIN)/N_ESAMP)
#define THETAE_MAX 1000.
#define THETAE_MIN 0.3
#define WEIGHT_MIN (1.e28)

#define SYNCHROTRON (1)
#define BREMSSTRAHLUNG (1)
#define COMPTON (1)

#define NVAR (8)  // number of primitive variables
#define KRHO     0
#define UU       1
#define U1       2
#define U2       3
#define U3       4
#define B1       5
#define B2       6
#define B3       7
//#define KEL      8
//#define KTOT     9

#define SMALL (1.e-40)
#define MMW   (0.5)

#define N_ESAMP 200
#define N_EBINS 200
#define N_THBINS 6

#include "hdf5_utils.h"
#include "h5io.h"

#ifndef HAMR_READ
#define HAMR_READ (1)
#endif

/* Temperature ratio (T_p/T_e) */
#define Monika_Te (1)   /* electron temperature by plasma beta (Monika Moscibradzka) */
#if (Monika_Te)
#define Rhigh (20.)     // Rhigh & Rlow would be overwritten by trat_large & trat_small in the parameter file if exists. 
#define Rlow  (1.)
#else
#define TPoTE (3.33)   // Constant temperature ratio (Tp/Te): overwritten by TP_OVER_TE in the parmeter file if exists.
#endif

#define read_dscale (0) /* read the density scale from GRMHD data set */

#define MBH_in 4.1e6  /* Sgr A* is updated by Gravity Collaboration 2018,615,L15). This will be overloaded by the value in the parameter file if it exists. */
