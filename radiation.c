

/* 

model-independent radiation-related utilities.

*/

#include "decs.h"

double Bnu_inv(double nu, double Thetae)
{

	double x;

	x = HPL * nu / (ME * CL * CL * Thetae);

	if (x < 1.e-3)		/* Taylor expand */
		return ((2. * HPL / (CL * CL)) /
			(x / 24. * (24. + x * (12. + x * (4. + x)))));
	else
		return ((2. * HPL / (CL * CL)) / (exp(x) - 1.));
}

double jnu_inv(double nu, double Thetae, double Ne, double B, double theta)
{
	double j;

	j = jnu(nu, Ne, Thetae, B, theta);

	return (j / (nu * nu));
}

/* return Lorentz invariant scattering opacity */
double alpha_inv_scatt(double nu, double Thetae, double Ne)
{
  #if COMPTON
	double kappa;

	kappa = kappa_es(nu, Thetae);

	return (nu * kappa * Ne * MP);
  #else
  return 0.;
  #endif
}

/* return Lorentz invariant absorption opacity */
double alpha_inv_abs(double nu, double Thetae, double Ne, double B,
		     double theta)
{
	double j, bnu;

	j = jnu_inv(nu, Thetae, Ne, B, theta);
	bnu = Bnu_inv(nu, Thetae);

	return (j / (bnu + 1.e-100));
}


/* return electron scattering opacity, in cgs */
double kappa_es(double nu, double Thetae)
{
	double Eg;

	/* assume pure hydrogen gas to 
	   convert cross section to opacity */
	Eg = HPL * nu / (ME * CL * CL);
	return (total_compton_cross_lkup(Eg, Thetae) / MP);
}

/* get frequency in fluid frame, in Hz */
double get_fluid_nu(double X[4], double K[4], double Ucov[NDIM])
{
	double ener, nu;

	/* this is the energy in electron rest-mass units */
	ener = -(K[0] * Ucov[0] +
		 K[1] * Ucov[1] + K[2] * Ucov[2] + K[3] * Ucov[3]);

	nu = ener * ME * CL * CL / HPL;

	if (isnan(ener)) {
		fprintf(stderr, "isnan get_fluid_nu, K: %g %g %g %g\n",
			K[0], K[1], K[2], K[3]);
		fprintf(stderr, "isnan get_fluid_nu, X: %g %g %g %g\n",
			X[0], X[1], X[2], X[3]);
		fprintf(stderr, "isnan get_fluid_nu, U: %g %g %g %g\n",
			Ucov[0], Ucov[1], Ucov[2], Ucov[3]);
	}

	return nu;

}

/* return angle between magnetic field and wavevector */
double get_bk_angle(double X[NDIM], double K[NDIM], double Ucov[NDIM],
		    double Bcov[NDIM], double B)
{

	double k, mu;

	if (B == 0.)
		return (M_PI / 2.);

	k = fabs(K[0] * Ucov[0] + K[1] * Ucov[1] + K[2] * Ucov[2] +
		 K[3] * Ucov[3]);

	/* B is in cgs but Bcov is in code units */
	mu = (K[0] * Bcov[0] + K[1] * Bcov[1] + K[2] * Bcov[2] +
	      K[3] * Bcov[3]) / (k * B / B_unit);

	if (fabs(mu) > 1.)
		mu /= fabs(mu);

	return (acos(mu));
}
