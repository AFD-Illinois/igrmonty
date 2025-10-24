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
//Bremss (only supported for thermal electrons) options
//0 No bremsstrahlung
//1 Based on Rybicki and Lightman eq. 5.14a with eq. 5.25 corrective factor
//2 Straub+ 2012 section 3.2.1
//3 van Hoof 2015 + Nozawa 2009
#define BREMSSTRAHLUNG (1)
#define COMPTON (1)

#define KRHO     0
#define UU       1
#define U1       2
#define U2       3
#define U3       4
#define B1       5
#define B2       6
#define B3       7
#define KEL      8
#define KTOT     9

#define SMALL (1.e-40)
#define MMW   (0.5)

#define N_ESAMP 200
#define N_EBINS 200
#define N_THBINS 12

#include "hdf5_utils.h"
#include "h5io.h"

