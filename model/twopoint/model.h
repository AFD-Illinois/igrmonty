#define NDIM 4
#define NUMIN 1.e9
#define NUMAX 1.e24
#define LNUMIN log(NUMIN)
#define LNUMAX log(NUMAX)
#define DLNU ((LNUMAX-LNUMIN)/N_ESAMP)
#define THETAE_MAX (1000.)
#define THETAE_MIN (0.3)
#define WEIGHT_MIN (1.e28)

#define SYNCHROTRON (1)
#define COMPTON (1)
#define KAPPA (5.)
#define DIST_KAPPA (0)

// Bremss options
// 0 No bremsstrahlung
// 1 Rybicki and Lightman eq. 5.14a with eq. 5.25 corrective factor
// 2 Piecewise formula
// 3 van Hoof 2015 + Nozawa 2009
// Bremss is only supported for thermal electrons
#define BREMSSTRAHLUNG (3)

#define KRHO     0
#define UU       1
#define U1       2
#define U2       3
#define U3       4
#define B1       5
#define B2       6
#define B3       7

#define SMALL (1.e-40)
#define MMW   (0.5)

#define N_ESAMP 200
#define N_EBINS 200
#define N_THBINS 3

#include "h5io.h"

#define HDF5_OUTPUT (1)

#define MODEL_TRANSPARENT (1)

extern double lnumin, lnumax, alpha_spec;

