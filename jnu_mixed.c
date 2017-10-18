

#include "decs.h"
//#pragma omp threadprivate(r)
/* 

"mixed" emissivity formula 

interpolates between Petrosian limit and
classical thermal synchrotron limit

good for Thetae > 1

*/

double jnu_synch(double nu, double Ne, double Thetae, double B, double theta);
double jnu_kappa(double nu, double Ne, double Thetae, double B, double theta);
double jnu_bremss(double nu, double Ne, double Thetae);
double int_jnu_synch(double Ne, double Thetae, double Bmag, double nu);
double int_jnu_kappa(double Ne, double Thetae, double Bmag, double nu);
double int_jnu_bremss(double Ne, double Thetae, double nu);

double jnu(double nu, double Ne, double Thetae, double B, double theta)
{
  double j = 0.;
  
  #if SYNCHROTRON
  #if DIST_KAPPA
  j += jnu_kappa(nu, Ne, Thetae, B, theta);
  #else
  j += jnu_synch(nu, Ne, Thetae, B, theta);
  #endif
  #endif

  #if BREMSSTRAHLUNG
  j += jnu_bremss(nu, Ne, Thetae);
  #endif
  
  return j;
}

double int_jnu(double Ne, double Thetae, double B, double nu)
{
  double intj = 0.;
  
  #if SYNCHROTRON
  intj += int_jnu_synch(Ne, Thetae, B, nu);
  #endif

  #if BREMSSTRAHLUNG
  intj += int_jnu_bremss(Ne, Thetae, nu);
  #endif
  
  return intj;
}

double jnu_bremss(double nu, double Ne, double Thetae)
{
  double Te = Thetae*ME*CL*CL/KBOL;
  double rel = (1. + 4.4e-10*Te);
  double x, efac;
  double gff = 1.2;

  x = HPL*nu/(KBOL*Te);

  if (x < 1.e-3) {
    efac = (24 - 24*x + 12*x*x - 4.*x*x*x + x*x*x*x)/24.;
  } else {
    efac = exp(-x);
  }

  double jv = 1./(4.*M_PI)*pow(2,5)*M_PI*pow(EE,6)/(3.*ME*pow(CL,3));
  jv *= pow(2.*M_PI/(3.*KBOL*ME),1./2.);
  jv *= pow(Te,-1./2.)*Ne*Ne;
  jv *= efac*rel*gff;

  return jv;
}

#define CST 1.88774862536	/* 2^{11/12} */
double jnu_synch(double nu, double Ne, double Thetae, double B,
		 double theta)
{
	double K2, nuc, nus, x, f, j, sth, xp1, xx;
	double K2_eval(double Thetae);

	if (Thetae < THETAE_MIN)
		return 0.;

	K2 = K2_eval(Thetae);

	nuc = EE * B / (2. * M_PI * ME * CL);
	sth = sin(theta);
	nus = (2. / 9.) * nuc * Thetae * Thetae * sth;
	if (nu > 1.e12 * nus)
		return (0.);
	x = nu / nus;
	xp1 = pow(x, 1. / 3.);
	xx = sqrt(x) + CST * sqrt(xp1);
	f = xx * xx;
	j = (M_SQRT2 * M_PI * EE * EE * Ne * nus / (3. * CL * K2)) * f *
	    exp(-xp1);

	return (j);
}

#include <gsl/gsl_sf_gamma.h>
double jnu_kappa(double nu, double Ne, double Thetae, double B, double theta)
{
  double kap = KAPPA;
	double nuc = EE * B / (2. * M_PI * ME * CL);
  double js = Ne*pow(EE,2)*nuc/CL;
  double x = 3.*pow(kap,-3./2.);
  double Jslo, Jshi;

  double nuk = nuc*pow(Thetae*kap,2)*sin(theta);
  double Xk = nu/nuk;

  Jslo = pow(Xk,1./3.)*sin(theta)*4.*M_PI*gsl_sf_gamma(kap-4./3.)/(pow(3.,7./3.)*gsl_sf_gamma(kap-2.));
  Jshi = pow(Xk,-(kap-2.)/2.)*sin(theta)*pow(3.,(kap-1.)/2.);
  Jshi *= (kap-2.)*(kap-1.)/4.*gsl_sf_gamma(kap/4.-1./3.)*gsl_sf_gamma(kap/4.+4./3.);

  double Js = pow(pow(Jslo,-x) + pow(Jshi,-x),-1./x);
  
  if (isnan(js*Js) || js*Js < 0. || js*Js > 1.e200) {
    printf("BAD jkap! %e\n", js*Js);
  }

  if (isnan(Jslo) || isinf(Jslo) || Jslo < 0. || Jslo > 1.e100) {
    printf("JSLO ERROR! %e\n", Jslo);
  }

  //return js*Jslo;

  return jnu_synch(nu, Ne, Thetae, B, theta);

  return js*Js;
}

#undef CST

#define JCST	(M_SQRT2*EE*EE*EE/(27*ME*CL*CL))
double int_jnu_synch(double Ne, double Thetae, double Bmag, double nu)
{
/* Returns energy per unit time at							*
 * frequency nu in cgs										*/

	double j_fac, K2;
	double F_eval(double Thetae, double B, double nu);
	double K2_eval(double Thetae);


	if (Thetae < THETAE_MIN)
		return 0.;

	K2 = K2_eval(Thetae);
	if (K2 == 0.)
		return 0.;

	j_fac = Ne * Bmag * Thetae * Thetae / K2;

	return JCST * j_fac * F_eval(Thetae, Bmag, nu);
}

double int_jnu_kappa(double Ne, double Thetae, double Bmag, double nu)
{
  return 0.;
}

#undef JCST

double int_jnu_bremss(double Ne, double Thetae, double nu)
{
  return 4.*M_PI*jnu_bremss(nu, Ne, Thetae);
}

#define CST 1.88774862536	/* 2^{11/12} */
double jnu_integrand(double th, void *params)
{

	double K = *(double *) params;
	double sth = sin(th);
	double x = K / sth;

	if (sth < 1.e-150 || x > 2.e8)
		return 0.;

	return sth * sth * pow(sqrt(x) + CST * pow(x, 1. / 6.),
			       2.) * exp(-pow(x, 1. / 3.));
}

#undef CST

/* Tables */
double F[N_ESAMP + 1], K2[N_ESAMP + 1];
double lK_min, dlK;
double lT_min, dlT;

#define EPSABS 0.
#define EPSREL 1.e-6
#define KMIN (0.002)
#define KMAX (1.e7)
#define TMIN (THETAE_MIN)
#define TMAX (1.e2)
void init_emiss_tables(void)
{

	int k;
	double result, err, K, T;
	gsl_function func;
	gsl_integration_workspace *w;

	func.function = &jnu_integrand;
	func.params = &K;

	lK_min = log(KMIN);
	dlK = log(KMAX / KMIN) / (N_ESAMP);

	lT_min = log(TMIN);
	dlT = log(TMAX / TMIN) / (N_ESAMP);

	/*  build table for F(K) where F(K) is given by
	   \int_0^\pi ( (K/\sin\theta)^{1/2} + 2^{11/12}(K/\sin\theta)^{1/6})^2 \exp[-(K/\sin\theta)^{1/3}]
	   so that J_{\nu} = const.*F(K)
	 */
	w = gsl_integration_workspace_alloc(1000);
	for (k = 0; k <= N_ESAMP; k++) {
		K = exp(k * dlK + lK_min);
		gsl_integration_qag(&func, 0., M_PI / 2., EPSABS, EPSREL,
				    1000, GSL_INTEG_GAUSS61, w, &result,
				    &err);
		F[k] = log(4 * M_PI * result);
	}
	gsl_integration_workspace_free(w);

	/*  build table for quick evaluation of the bessel function K2 for emissivity */
	for (k = 0; k <= N_ESAMP; k++) {
		T = exp(k * dlT + lT_min);
		K2[k] = log(gsl_sf_bessel_Kn(2, 1. / T));

	}

	/* Avoid doing divisions later */
	dlK = 1. / dlK;
	dlT = 1. / dlT;

	fprintf(stderr, "done.\n\n");

	return;
}

/* rapid evaluation of K_2(1/\Thetae) */

double K2_eval(double Thetae)
{

	double linear_interp_K2(double);

	if (Thetae < THETAE_MIN)
		return 0.;
	if (Thetae > TMAX)
		return 2. * Thetae * Thetae;

	return linear_interp_K2(Thetae);
}

#define KFAC	(9*M_PI*ME*CL/EE)
double F_eval(double Thetae, double Bmag, double nu)
{

	double K, x;
	double linear_interp_F(double);

	K = KFAC * nu / (Bmag * Thetae * Thetae);

	if (K > KMAX) {
		return 0.;
	} else if (K < KMIN) {
		/* use a good approximation */
		x = pow(K, 0.333333333333333333);
		return (x * (37.67503800178 + 2.240274341836 * x));
	} else {
		return linear_interp_F(K);
	}
}

#undef KFAC
#undef KMIN
#undef KMAX
#undef EPSABS
#undef EPSREL

double linear_interp_K2(double Thetae)
{

	int i;
	double di, lT;

	lT = log(Thetae);

	di = (lT - lT_min) * dlT;
	i = (int) di;
	di = di - i;

	return exp((1. - di) * K2[i] + di * K2[i + 1]);
}

double linear_interp_F(double K)
{

	int i;
	double di, lK;

	lK = log(K);

	di = (lK - lK_min) * dlK;
	i = (int) di;
	di = di - i;

	return exp((1. - di) * F[i] + di * F[i + 1]);
}
