/*

model-independent radiation-related utilities.

*/

#include "decs.h"
#include "model_radiation.h"
#include "par.h"

// this file defines:
//
//   Bnu_inv
//   jnu_inv
//   alpha_inv_scatt
//   alpha_inv_abs
//   kappa_es
//   get_fluid_nu
//   get_bk_angle
//

double model_kappa = 4.;

double powerlaw_gamma_cut = 1.e10;
double powerlaw_gamma_min = 1.e2;
double powerlaw_gamma_max = 1.e5;
double powerlaw_p = 3.25;


void try_set_radiation_parameter(const char *line)
{
  read_param(line, "powerlaw_gamma_cut", &powerlaw_gamma_cut, TYPE_DBL);
  read_param(line, "powerlaw_gamma_min", &powerlaw_gamma_min, TYPE_DBL);
  read_param(line, "powerlaw_gamma_max", &powerlaw_gamma_max, TYPE_DBL);
  read_param(line, "powerlaw_p", &powerlaw_p, TYPE_DBL);
}

// determine w by finding effective w for total
// energy to match thermal (MJ) at Thetae
double kappa_w(double Thetae, double kappa)
{
  return (kappa - 3.)/kappa * Thetae;
}

// planck function
double Bnu_inv(double nu, double Thetae)
{
	double x = HPL * nu / (ME * CL * CL * Thetae);

	if (x < 1.e-3) { // Taylor expand if small
		return ((2. * HPL / (CL * CL)) /
			(x / 24. * (24. + x * (12. + x * (4. + x)))));
  }

	return (2. * HPL / (CL * CL)) / (exp(x) - 1.);
}

// return j_\nu/\nu^2, the invariant emissivity
double jnu_inv(double nu, double Thetae, double Ne, double B, double theta)
{
	double j = jnu(nu, Ne, Thetae, B, theta);

	return j / (nu * nu);
}

// return invariant scattering opacity if Compton scattering enabled
double alpha_inv_scatt(double nu, double Thetae, double Ne)
{
  #if COMPTON

	return nu * kappa_es(nu, Thetae) * Ne * MP;

  #else

  return 0.;

  #endif
}

// return invariant absorption opacity 
double alpha_inv_abs(double nu, double Thetae, double Ne, double B,
		     double theta)
{

#if ( ( BREMSSTRAHLUNG != 0 ) && (MODEL_EDF==EDF_KAPPA_FIXED) )
  fprintf(stderr, "ERROR absorptivities not set up for bremss and kappa!\n");
  exit(-1);
#endif

#if MODEL_EDF==EDF_KAPPA_FIXED

  // Pandya+ 2016 absorptivity
 
  double kap = model_kappa;
  double w = kappa_w(Thetae, kap);

  double nuc = EE*B/(2.*M_PI*ME*CL);
  double nuk = nuc*pow(w*kap,2)*sin(theta);
  double Xk = nu/nuk;

  double Aslo, Ashi, As;

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

#elif MODEL_EDF==EDF_POWER_LAW

  /*
  double p = powerlaw_p;
  double gmin = powerlaw_gamma_min;
  double gmax = powerlaw_gamma_max;

  double sth = sin(theta);
  double nuc = EE * B / (2.*M_PI*ME*CL);
  double factor = (Ne * EE*EE)/(nu * ME*CL);

  double X = nu/(nuc*sth);

  double As = pow(3.,(p+1)/2.)*(p-1)/(4*(pow(gmin,1-p)-pow(gmax,1-p)));
  As *= gsl_sf_gamma((3*p+2)/12.)*gsl_sf_gamma((3*p+22)/12.)*pow(1./3.*X,-(p+2)/2.);

  return As*factor;
   */

  double sth = sin(theta);
  double nu_c = EE * B / (2 * M_PI * ME * CL); 

  double prefactor = Ne * EE*EE / (nu * ME * CL);

  double t1 = pow(3., (powerlaw_p+1)/2.) * (powerlaw_p - 1.);
  double t2 = 4. * (pow(powerlaw_gamma_min, 1.-powerlaw_p) - 
                    pow(powerlaw_gamma_max, 1.-powerlaw_p));
  double t3 = tgamma((3*powerlaw_p+2)/12.) * tgamma((3.*powerlaw_p+22)/12.);
  double t4 = pow(nu/(nu_c * sth), -(powerlaw_p+2)/2.);

  return nu * prefactor * t1 / t2 * t3 * t4;

#elif MODEL_EDF==EDF_MAXWELL_JUTTNER

	double j = jnu_inv(nu, Thetae, Ne, B, theta);
	double bnu = Bnu_inv(nu, Thetae);

  if (j > 0) {
	  return j / (bnu + 1.e-100);
  }

  return 0;

#else

  fprintf(stderr, "must select valid MODEL_EDF\n");
  exit(3);

#endif 
}


// return electron scattering opacity in cgs
double kappa_es(double nu, double Thetae)
{

	// assume pure hydrogen gas to
	// convert cross section to opacity
	
	double Eg = HPL * nu / (ME * CL * CL);

  if (Eg > 1.e75) {
    fprintf(stderr, "out of bounds: %g %g %g\n", Eg, Thetae, nu);
  }

	return total_compton_cross_lkup(Eg, Thetae) / MP;
}

// get frequency in fluid frame, in Hz
double get_fluid_nu(const double X[NDIM], const double K[NDIM], const double Ucov[NDIM])
{
	// in electron rest-mass units 
	double energy = -(K[0]*Ucov[0] + K[1]*Ucov[1] + K[2]*Ucov[2] + K[3]*Ucov[3]);

  // in Hz
	double nu = energy * ME * CL * CL / HPL;

	if (isnan(energy)) {
		fprintf(stderr, "isnan get_fluid_nu, K: %g %g %g %g\n",
			K[0], K[1], K[2], K[3]);
		fprintf(stderr, "isnan get_fluid_nu, X: %g %g %g %g\n",
			X[0], X[1], X[2], X[3]);
		fprintf(stderr, "isnan get_fluid_nu, U: %g %g %g %g\n",
			Ucov[0], Ucov[1], Ucov[2], Ucov[3]);
	}

	return nu;
}

// return angle between magnetic field and wavevector
double get_bk_angle(double X[NDIM], double K[NDIM], double Ucov[NDIM],
		    double Bcov[NDIM], double B)
{
	double k, mu;

	if (B == 0.) {
		return M_PI / 2.;
  }

	k = fabs(K[0]*Ucov[0] + K[1]*Ucov[1] + K[2]*Ucov[2] + K[3]*Ucov[3]);

	// B is in cgs but Bcov is in code units
	mu = (K[0] * Bcov[0] + K[1] * Bcov[1] + K[2] * Bcov[2] + K[3] * Bcov[3]) / (k * B / B_unit);

	if (fabs(mu) > 1.) {
		mu /= fabs(mu);
  }

	return acos(mu);

	(void)X; // silence unused parameter warning 
}
