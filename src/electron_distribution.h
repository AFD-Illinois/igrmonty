#ifndef ELECTRON_DISTRIBUTION_H
#define ELECTRON_DISTRIBUTION_H
#include "constants.h"
#include "decs.h"
double dnd3p_bimaxwell(double tperp_over_tpar, double p, double ne, double thetae_perp, double p_par, double p_perp);
double dnd3p_bimaxwell_fast(double tperp_over_tpar, double p, double ne, double thetae_perp, double p_par, double p_perp);
double dnd3p_bimaxwell_prefactor(double tperp_over_tpar, double ne, double thetae_perp);


#endif // ELECTRON_DISTRIBUTION_H
