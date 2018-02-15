#define NDIM 4
#define NUMIN 1.e9
#define NUMAX 1.e16
#define LNUMIN log(NUMIN)
#define LNUMAX log(NUMAX)
#define DLNU ((LNUMAX-LNUMIN)/N_ESAMP)
#define THETAE_MAX (1000.)
#define THETAE_MIN (0.3)
#define WEIGHT_MIN (1.e28)

#define SYNCHROTRON (1)
#define BREMSSTRAHLUNG (0)
#define COMPTON (0)
#define KAPPA (5)
#define DIST_KAPPA (1)

#define KRHO     0
#define UU       1
#define U1       2
#define U2       3
#define U3       4
#define B1       5
#define B2       6
#define B3       7
#define NPRIM    8

#define SMALL (1.e-40)
#define MMW   (0.5)

#define N_ESAMP 200
#define N_EBINS 200
#define N_THBINS 3

