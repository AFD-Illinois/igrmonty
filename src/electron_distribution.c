#include "electron_distribution.h"

/**
 * Compute normalized anisotropic distribution function value for a given 
 * electron momentum, anisotropy value, and perpendicular temperature.
 *
 * Parameters:
 *   tperp_over_tpar            - Anisotropy parameter, defined as T_perp/T_parallel
 *   p            - Momentum of electron (units of momentum/mc)
 *   theta        - Angle the momentum makes with respect to T_parallel (or B field)
 *   ne          - Number density, in cgs units
 *   thetae_perp  - Temperature of the perpendicular component, in dimensionless units
 *
 * Returns:
 *   dne_d3p - The normalized distribution function (dne/d^3p), see Treumann et al. 2016
 */
double dnd3p_bimaxwell(double tperp_over_tpar, double p, double ne, double thetae_perp, double p_par, double p_perp)
{
    double psq_perp = p_perp*p_perp;
    double psq_par = p_par*p_par;

    // multiply K2(1/Thetae) by e^(1/Thetae) for numerical purposes
    double K2f;
    if (thetae_perp > 1.e-2) {
        K2f = gsl_sf_bessel_Kn(2, 1. / thetae_perp) * exp(1. / thetae_perp);
    } else {
        K2f = sqrt(M_PI * thetae_perp / 2.);
    }

    double prefactor = ne * sqrt(tperp_over_tpar) / (4 * M_PI) / (thetae_perp*K2f);

    return prefactor * exp(-(sqrt(1 + psq_perp + tperp_over_tpar*psq_par)-1)/(thetae_perp));
}

// fast version of the bimaxwell distribution function dne_d3p, without the K2(1/Thetae) term or any prefactor.
// used for trapezoid integration routines
double dnd3p_bimaxwell_fast(double tperp_over_tpar, double p, double ne, double thetae_perp, double p_par, double p_perp)
{
    double psq_perp = p_perp*p_perp;
    double psq_par = p_par*p_par;

    return exp(-(sqrt(1.0 + psq_perp + tperp_over_tpar*psq_par)-1.0)/(thetae_perp));
}

// prefactor for the bimaxwell distribution function dne_d3p
double dnd3p_bimaxwell_prefactor(double tperp_over_tpar, double ne, double thetae_perp)
{
    // multiply K2(1/Thetae) by e^(1/Thetae) for numerical purposes
    double K2f;
    if (thetae_perp > 1.e-2) {
        K2f = gsl_sf_bessel_Kn(2, 1. / thetae_perp) * exp(1. / thetae_perp);
    } else {
        K2f = sqrt(M_PI * thetae_perp / 2.);
    }

    return ne * sqrt(tperp_over_tpar) / (4 * M_PI) / (thetae_perp*K2f);
}