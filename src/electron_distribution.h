#ifndef ELECTRON_DISTRIBUTION_H
#define ELECTRON_DISTRIBUTION_H
#include "constants.h"
#include "decs.h"
double dnd3p_bimaxwell(double A, double p, double ne, double thetae_perp, double p_par, double p_perp);
double dnd3p_bimaxwell_fast(double A, double p, double ne, double thetae_perp, double p_par, double p_perp);
double dnd3p_bimaxwell_prefactor(double A, double ne, double thetae_perp);


#endif // ELECTRON_DISTRIBUTION_H
