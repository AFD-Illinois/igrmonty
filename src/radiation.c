

/* 

model-independent radiation-related utilities.

*/

#include "decs.h"

/* planck function */
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

/* return j_\nu/\nu^2, the invariant emissivity */
double jnu_inv(double nu, double Thetae, double Ne, double B, double theta)
{
	double j;

	j = jnu(nu, Ne, Thetae, B, theta);

	return (j / (nu * nu));
}

/* return invariant scattering opacity */
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

/* return invariant absorption opacity */
double alpha_inv_abs(double nu, double Thetae, double Ne, double B,
		     double theta)
{
  #if DIST_KAPPA && ( BREMSSTRAHLUNG != 0)
  printf("ERROR absorptivities not set up for bremss and kappa!\n");
  exit(-1);
  #endif
  
  #if DIST_KAPPA
  // Pandya+ 2016 absorptivity
  double Aslo, Ashi, As;

  double kap = KAPPA;
  double w = Thetae;
  double nuc = EE*B/(2.*M_PI*ME*CL);
  double nuk = nuc*pow(w*kap,2)*sin(theta);
  double Xk = nu/nuk;

  Aslo  = pow(Xk,-2./3.)*pow(3.,1./6.)*10./41.;
  Aslo *= 2.*M_PI/(pow(w*kap,10./3.-kap));
  Aslo *= (kap - 2.)*(kap - 1.)*kap/(3.*kap - 1.);
  Aslo *= gsl_sf_gamma(5./3.); 
  // Evaluate 2F1(a,b;c,z), using analytic continuation if |z| > 1
  double a = kap - 1./3.;
  double b = kap + 1.;
  double c = kap + 2./3.;
  double z = -kap*w;
  double hg2F1;
  if (fabs(z) == 1.) {
    hg2F1 = 0.;
  } else if (fabs(z) < 1.) {
    hg2F1 = gsl_sf_hyperg_2F1(a, b, c, z);
  } else {
    hg2F1  = pow(1.-z,-a)*gsl_sf_gamma(c)*gsl_sf_gamma(b-a)/(gsl_sf_gamma(b)*gsl_sf_gamma(c-a))*gsl_sf_hyperg_2F1(a,c-b,a-b+1,1./(1.-z));
    hg2F1 += pow(1.-z,-b)*gsl_sf_gamma(c)*gsl_sf_gamma(a-b)/(gsl_sf_gamma(a)*gsl_sf_gamma(c-b))*gsl_sf_hyperg_2F1(b,c-a,b-a+1,1./(1.-z));
  }
  Aslo *= hg2F1;

  Ashi  = pow(Xk,-(1. + kap)/2.)*pow(M_PI,3./2.)/3.;
  Ashi *= (kap - 2.)*(kap - 1.)*kap/pow(w*kap,3.);
  Ashi *= (2.*gsl_sf_gamma(2. + kap/2.)/(2. + kap) - 1.);
  Ashi *= (pow(3./kap,19./4.) + 3./5.);

  double xbr = pow(-7./4. + 8./5.*kap,-43./50.);

  As = pow(pow(Aslo,-xbr) + pow(Ashi,-xbr),-1./xbr);
  double alphas = Ne*EE*EE/(nu*ME*CL)*As;
  double cut = exp(-nu/NUCUT);

  return nu*alphas*cut;

  #else
	double j, bnu;

	j = jnu_inv(nu, Thetae, Ne, B, theta);
	bnu = Bnu_inv(nu, Thetae);

  //double alpha_kirch = j/(bnu + 1.e-100);

	return (j / (bnu + 1.e-100));
  #endif // DIST_KAPPA
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
