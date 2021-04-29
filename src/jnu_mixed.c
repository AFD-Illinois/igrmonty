#include "decs.h"
#include "model_radiation.h"

#include <assert.h>

/* 

"mixed" emissivity formula 

interpolates between Petrosian limit and
classical thermal synchrotron limit

good for Thetae > 1

*/

// exposed functions are:
//   jnu(double nu, double Ne, double Thetae, double B, double theta)
//      -> directly computes emissivity from formulae. used in radiation.c
//   jnu_ratio_brems(double nu, double Ne, double Thetae, double B, double theta)
//      -> directly computes by calling from above. used to track spectrum components
//   int_jnu(double Ne, double Thetae, double B, double nu)
//      -> precomputed numerically, saved in table, then interpolated. used to get nsph / zone.

// these functions are scoped here only and called from above "public" functions.
static double jnu_thermal(double nu, double Ne, double Thetae, double B, double theta);
static double jnu_kappa(double, double, double, double, double, radiation_params *);
static double jnu_powerlaw(double nu, double Ne, double Thetae, double B, double theta);
static double jnu_bremss(double nu, double Ne, double Thetae);
static double int_jnu_thermal(double Ne, double Thetae, double Bmag, double nu);
static double int_jnu_kappa(double, double, double, double, radiation_params *);
static double int_jnu_powerlaw(double Ne, double Thetae, double Bmag, double nu);
static double int_jnu_bremss(double Ne, double Thetae, double nu);

// these are helper functions
double K2_eval(double Thetae);
double F_eval(double Thetae, double B, double nu);
double G_eval(double Thetae, double B, double nu, double kappa);
double H_eval(double Thetae, double B, double nu);
double linear_interp_K2(double);
double linear_interp_F(double);
double linear_interp_G(double, double);
double linear_interp_H(double);
void populate_table_from_function(double table[], double (* function) (double x, void * params), radiation_params *rpars);

typedef struct int_rpars_struct {
  double value;
  radiation_params *rpars;
} int_rpars;

double jnu(double nu, double Ne, double Thetae, double B, double theta, radiation_params *rpars)
{
  double j = 0.;
  
#if SYNCHROTRON
 #if (MODEL_EDF==EDF_KAPPA_FIXED) || (MODEL_EDF==EDF_KAPPA_VARIABLE)
  j += jnu_kappa(nu, Ne, Thetae, B, theta, rpars);
 #elif MODEL_EDF==EDF_MAXWELL_JUTTNER
  j += jnu_thermal(nu, Ne, Thetae, B, theta);
 #elif MODEL_EDF==EDF_POWER_LAW
  j += jnu_powerlaw(nu, Ne, Thetae, B, theta);
 #else
  fprintf(stderr, "must choose valid MODEL_EDF\n");
  exit(3);
 #endif
#endif

  #if BREMSSTRAHLUNG
  j += jnu_bremss(nu, Ne, Thetae);
  #endif
  
  return j;
}

double jnu_ratio_brems(double nu, double Ne, double Thetae, double B, double theta, radiation_params *rpars)
{
  double synch = 0.;
  double brems = 0.;

#if SYNCHROTRON
 #if (MODEL_EDF==EDF_KAPPA_FIXED) || (MODEL_EDF==EDF_KAPPA_VARIABLE)
  synch = jnu_kappa(nu, Ne, Thetae, B, theta, rpars);
 #elif MODEL_EDF==EDF_MAXWELL_JUTTNER
  synch = jnu_thermal(nu, Ne, Thetae, B, theta);
 #elif MODEL_EDF==EDF_POWER_LAW
  synch = jnu_powerlaw(nu, Ne, Thetae, B, theta);
 #else
  fprintf(stderr, "must choose valid MODEL_EDF\n");
  exit(3);
 #endif  // MODEL_EDF
#endif  // SYNCHROTRON

  #if BREMSSTRAHLUNG
  brems = jnu_bremss(nu, Ne, Thetae);
  #endif // BREMSSTRAHLUNG

  if ( synch + brems == 0 ) return 0.;
  
  return brems / ( synch + brems );

  // silence unused warnings
  (void)jnu_bremss;
  (void)jnu_thermal;
  (void)jnu_kappa;
  (void)jnu_powerlaw;
}

double int_jnu(double Ne, double Thetae, double B, double nu, radiation_params *rpars)
{
  double intj = 0.;
  
#if SYNCHROTRON
 #if (MODEL_EDF==EDF_KAPPA_FIXED) || (MODEL_EDF==EDF_KAPPA_VARIABLE)
  intj += int_jnu_kappa(Ne, Thetae, B, nu, rpars);
 #elif MODEL_EDF==EDF_MAXWELL_JUTTNER
  intj += int_jnu_thermal(Ne, Thetae, B, nu);
 #elif MODEL_EDF==EDF_POWER_LAW
  intj += int_jnu_powerlaw(Ne, Thetae, B, nu);
 #else
  fprintf(stderr, "must choose valid MODEL_EDF\n");
  exit(3);
 #endif  // MODEL_EDF
#endif  // SYNCHROTRON

#if BREMSSTRAHLUNG
  intj += int_jnu_bremss(Ne, Thetae, nu);
#endif
  
  return intj;

  // silence unused warnings
  (void)int_jnu_bremss;
  (void)int_jnu_thermal;
  (void)int_jnu_kappa;
  (void)int_jnu_powerlaw;
}

static double jnu_bremss(double nu, double Ne, double Thetae)
{
  if (Thetae < THETAE_MIN) 
    return 0.;

  double Te = Thetae * ME * CL * CL / KBOL;
  double x = HPL*nu/(KBOL*Te);
  double efac = 0.;
  double gff = 1.2;
  
  if (x < 1.e-3) {
    efac = (24. - 24.*x + 12.*x*x - 4.*x*x*x + x*x*x*x) / 24.;
  } else {
    efac = exp(-x);
  }

#if 1   // following Straub+ 2012
  double Fei=0., Fee=0., fei=0., fee=0.;

  double SOMMERFELD_ALPHA = 1. / 137.036;
  double eta = 0.5616;
  double e_charge = 4.80e-10; // in esu
  double re = e_charge * e_charge / ME / CL / CL;
  double gammaE = 0.577; // = - Log[0.5616] 

  if (x > 1) {
    gff = sqrt(3. / M_PI / x);
  } else {
    gff = sqrt(3.) / M_PI * log(4 / gammaE / x);
  }

  if (Thetae < 1) {
    Fei = 4. * sqrt(2.*Thetae/M_PI/M_PI/M_PI) * (1. + 1.781*pow(Thetae,1.34));
    Fee = 20./9./sqrt(M_PI) * (44. - 3.*M_PI*M_PI) * pow(Thetae,1.5);
    Fee *= (1. + 1.1*Thetae + Thetae*Thetae - 1.25*pow(Thetae,2.5));
  } else {
    Fei = 9.*Thetae/(2.*M_PI) * ( log(1.123 * Thetae + 0.48) + 1.5 );
    Fee = 24. * Thetae * ( log(2.*eta*Thetae) + 1.28 );
  }
  fei = Ne * Ne * SIGMA_THOMSON * SOMMERFELD_ALPHA * ME * CL * CL * CL * Fei;
  fee = Ne * Ne * re * re * SOMMERFELD_ALPHA * ME * CL * CL * CL * Fee;

  return (fei+fee) / (4.*M_PI) * HPL/KBOL/Te * efac * gff;
 
#else 
  // Method from Rybicki & Lightman, ultimately from Novikov & Thorne

  double rel = (1. + 4.4e-10*Te);

  double jv = 1./(4.*M_PI)*pow(2,5)*M_PI*pow(EE,6)/(3.*ME*pow(CL,3));
  jv *= pow(2.*M_PI/(3.*KBOL*ME),1./2.);
  jv *= pow(Te,-1./2.)*Ne*Ne;
  jv *= efac*rel*gff;

  return jv;
#endif 

}

#define CST 1.88774862536	/* 2^{11/12} */
static double jnu_thermal(double nu, double Ne, double Thetae, double B,
			double theta)
{
	double K2, nuc, nus, x, f, j, sth, xp1, xx;
	double K2_eval(double Thetae);

	if (Thetae < THETAE_MIN) {
		return 0.;
  }

	K2 = K2_eval(Thetae);

	nuc = EE * B / (2. * M_PI * ME * CL);
	sth = sin(theta);
	nus = (2. / 9.) * nuc * Thetae * Thetae * sth;

	if (nu > 1.e12 * nus) {
		return 0.;
  }

	x = nu / nus;
	xp1 = pow(x, 1. / 3.);
	xx = sqrt(x) + CST * sqrt(xp1);
	f = xx * xx;
	j = (M_SQRT2 * M_PI * EE * EE * Ne * nus / (3. * CL * K2)) * f *
	    exp(-xp1);

	return j;
}

static double jnu_powerlaw(double nu, double Ne, double Thetae, double B, double theta)
{
  if (Thetae < THETAE_MIN) {
    return 0.;
  }
  if (theta < SMALL || theta > M_PI-SMALL) {
    return 0.;
  }

  double p = powerlaw_p;
  double gmin = powerlaw_gamma_min;
  double gmax = powerlaw_gamma_max;

  double sth = sin(theta);
  double nuc = EE * B / (2. * M_PI * ME * CL);
  double factor = (Ne * pow(EE,2.) * nuc)/CL;

  if (nu > 1.e8 * nuc) {
    return 0.;
  }

  double Xs = nu/(nuc*sth);

  double Js = pow(3.,p/2.)*(p-1)*sth/(2*(p+1)*(pow(gmin,1-p)-pow(gmax,1-p)));
  Js *= gsl_sf_gamma((3*p-1)/12.)*gsl_sf_gamma((3*p+19)/12.)*pow(Xs,-(p-1)/2.);

  return Js*factor;
}

#include <gsl/gsl_sf_gamma.h>
static double jnu_kappa(double nu, double Ne, double Thetae, double B, double theta, radiation_params *rpars)
{
	if (Thetae < THETAE_MIN) {
		return 0.;
  }
  if (theta < SMALL || theta > M_PI-SMALL) {
    return 0.;
  } 

  double kap = rpars->kappa;
	double nuc = EE * B / (2. * M_PI * ME * CL);
  double js = Ne*pow(EE,2)*nuc/CL;
  double x = 3.*pow(kap,-3./2.);
  double Jslo, Jshi;

  double w = kappa_w(Thetae, kap);
  double nuk = nuc * w*kap * w*kap * sin(theta);
  double Xk = nu/nuk;

  Jslo = pow(Xk,1./3.)*sin(theta)*4.*M_PI*gsl_sf_gamma(kap-4./3.)/(pow(3.,7./3.)*gsl_sf_gamma(kap-2.));
  Jshi = pow(Xk,-(kap-2.)/2.)*sin(theta)*pow(3.,(kap-1.)/2.);
  Jshi *= (kap-2.)*(kap-1.)/4.*gsl_sf_gamma(kap/4.-1./3.)*gsl_sf_gamma(kap/4.+4./3.);

  double Js = pow(pow(Jslo,-x) + pow(Jshi,-x),-1./x);

  if (isnan(js*Js) || js*Js < 0. || js*Js > 1.e200 || js*Js < 1.e-200) {
    printf("BAD jkap! %e\n", js*Js);
    printf("nu Ne Thetae B theta = %e %e %e %e %e\n", nu, Ne, Thetae, B, theta);
    exit(-1);
  }

  if (isnan(Jslo) || isinf(Jslo) || Jslo < 0. || Jslo > 1.e100) {
    printf("JSLO ERROR! %e\n", Jslo);
  }

  double cut = exp(-nu/NUCUT);

  return js * Js * cut;
}

#undef CST

#define JCST	(M_SQRT2*EE*EE*EE/(27*ME*CL*CL))
static double int_jnu_thermal(double Ne, double Thetae, double Bmag, double nu)
{
  // Returns energy per unit time at frequency nu, all in cgs

	double j_fac, K2;

	if (Thetae < THETAE_MIN) {
		return 0.;
  }

	K2 = K2_eval(Thetae);
	if (K2 == 0.) {
		return 0.;
  }

	j_fac = Ne * Bmag * Thetae * Thetae / K2;

	return JCST * j_fac * F_eval(Thetae, Bmag, nu);
}

static double jnu_kappa_integrand(double th, void *params)
{
	double K = ((int_rpars *)params)->value;
	double sth = sin(th);
	double Xk = K / sth;
  double kap = ((int_rpars *)params)->rpars->kappa;

	if (sth < 1.e-150 || Xk > 2.e8) {
		return 0.;
  }

  double GAM1 = gsl_sf_gamma(kap - 4./3.);
  double GAM2 = gsl_sf_gamma(kap - 2.);
  double GAM3 = gsl_sf_gamma(kap/4. - 1./3.);
  double GAM4 = gsl_sf_gamma(kap/4. + 4./3.);

  double Jslo = pow(Xk,1./3.)*sth*4.*M_PI*GAM1/(pow(3.,7./3.)*GAM2);
  double Jshi = pow(Xk,-(kap-2.)/2.)*sth*pow(3.,(kap-1.)/2.)*(kap-2.)*(kap-1.)/4.*GAM3*GAM4;

  double x = 3.*pow(kap,-3./2.);
  double Js = pow(pow(Jslo,-x) + pow(Jshi,-x),-1./x);

  return Js * sth;
}

static double jnu_powerlaw_integrand(double th, void *params)
{
  double K = ((int_rpars *)params)->value;
  double sth = sin(th);
  double x = K / sth;

  double p = powerlaw_p;
  double gmin = powerlaw_gamma_min;
  double gmax = powerlaw_gamma_max;

  double factor = sth;

  double Js = pow(3.,p/2.)*(p-1)*sth/(2*(p+1)*(pow(gmin,1-p)-pow(gmax,1-p)));
  Js *= gsl_sf_gamma((3*p-1)/12.)*gsl_sf_gamma((3*p+19)/12.)*pow(x,-(p-1)/2.);

  return Js*factor;
}

static double int_jnu_powerlaw(double Ne, double Thetae, double B, double nu)
{
  if (Thetae < THETAE_MIN) {
    return 0.;
  }

  double CONST = EE * EE * EE / (2. * M_PI * ME * CL * CL);
  return CONST * Ne * B * H_eval(Thetae, B, nu);
}

static double int_jnu_kappa(double Ne, double Thetae, double B, double nu, radiation_params *rpars)
{
	if (Thetae < THETAE_MIN) {
		return 0.;
  }

  double nuc = EE*B/(2.*M_PI*ME*CL);
	double js = Ne*EE*EE*nuc/CL;
  double cut = exp(-nu/NUCUT);

	return js*G_eval(Thetae, B, nu, rpars->kappa)*cut;
}
#undef JCST

static double int_jnu_bremss(double Ne, double Thetae, double nu)
{
  return 4.*M_PI*jnu_bremss(nu, Ne, Thetae);
}

#define CST 1.88774862536	/* 2^{11/12} */
static double jnu_thermal_integrand(double th, void *params)
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
static double _F[N_ESAMP + 1];
static double _G[KAPPA_NSAMP][N_ESAMP + 1];  // kappa edf
static double _H[N_ESAMP + 1];  // power law
static double _K2[N_ESAMP + 1];
static double lK_min, dlK;
static double lL_min, dlL;
static double lT_min, dlT;

#define EPSABS 0.
#define EPSREL 1.e-6
#define KMIN (0.002)
#define KMAX (1.e7)
#define LMIN (0.002)
#define LMAX (1.e7)
#define TMIN (THETAE_MIN)
#define TMAX (1.e2)
void init_emiss_tables(void)
{

	int k;
	double result, err, K, T;

  // Thermal synchrotron lookup table
  {
    gsl_function func;
    gsl_integration_workspace *w;

    func.function = &jnu_thermal_integrand;
    func.params = &K;

    lK_min = log(KMIN);
    dlK = log(KMAX / KMIN) / (N_ESAMP);

    //  build table for F(K) where F(K) is given by
    //   \int_0^\pi ( (K/\sin\theta)^{1/2} + 2^{11/12}(K/\sin\theta)^{1/6})^2 \exp[-(K/\sin\theta)^{1/3}]
    // so that J_{\nu} = const.*F(K)
    w = gsl_integration_workspace_alloc(1000);
    for (k = 0; k <= N_ESAMP; k++) {
      K = exp(k * dlK + lK_min);
      gsl_integration_qag(&func, 0., M_PI / 2., EPSABS, EPSREL,
              1000, GSL_INTEG_GAUSS61, w, &result,
              &err);
      _F[k] = log(4 * M_PI * result);
    }
    gsl_integration_workspace_free(w);
  }

  radiation_params rpars;

  // kappa lookup table
#if MODEL_EDF==EDF_KAPPA_FIXED
  rpars.kappa = model_kappa;
  populate_table_from_function(_G[0], &jnu_kappa_integrand, &rpars);
#elif MODEL_EDF==EDF_KAPPA_VARIABLE
  for (int i=0; i<KAPPA_NSAMP; ++i) {
    rpars.kappa = KAPPA_MIN + i*DKAPPA;
    populate_table_from_function(_G[i], &jnu_kappa_integrand, &rpars);
  }
#endif

  // powerlaw lookup table
  populate_table_from_function(_H, &jnu_powerlaw_integrand, &rpars);

	// Bessel K2 lookup table
  {
    lT_min = log(TMIN);
    dlT = log(TMAX / TMIN) / (N_ESAMP);
    for (k = 0; k <= N_ESAMP; k++) {
		  T = exp(k * dlT + lT_min);
		  _K2[k] = log(gsl_sf_bessel_Kn(2, 1. / T));
	  }
  }
}

void populate_table_from_function(double table[], double (* function) (double x, void * params), radiation_params *rpars)
{
  int_rpars irp;
  double result, err;
  gsl_function func;
  gsl_integration_workspace *w;

  func.function = function;
  func.params = &irp;
  irp.rpars = rpars;

  lL_min = log(LMIN);
  dlL = log(LMAX / LMIN) / (N_ESAMP);

  w = gsl_integration_workspace_alloc(1000);
  for (int k = 0; k <= N_ESAMP; k++) {
    irp.value = exp(k * dlL + lL_min);
    gsl_integration_qag(&func, 0., M_PI / 2., EPSABS, EPSREL, 1000,
      GSL_INTEG_GAUSS61, w, &result, &err);
    table[k] = log(4*M_PI*result);
  }
  gsl_integration_workspace_free(w);
}

// rapid evaluation of K_2(1/\Thetae) 

double K2_eval(double Thetae)
{
	if (Thetae < THETAE_MIN)
		return 0.;
	if (Thetae > TMAX)
		return 2. * Thetae * Thetae;

	return linear_interp_K2(Thetae);
}

#define KFAC (9*M_PI*ME*CL/EE)
double F_eval(double Thetae, double Bmag, double nu)
{
	double K, x;

	K = KFAC * nu / (Bmag * Thetae * Thetae);

	if (K > KMAX) {
		return 0.;
	} else if (K < KMIN) {
		// use a good approximation
		x = pow(K, 0.333333333333333333);
		return (x * (37.67503800178 + 2.240274341836 * x));
	} else {
		return linear_interp_F(K);
	}
}


double H_eval(double Thetae, double Bmag, double nu) {

  double K;
  double nuc = EE * Bmag / (2. * M_PI * ME * CL);

  K = nu / nuc;
  if (K > KMAX)
    return 0.;
  if (K < KMIN)
    return 0.;
  double H_value = linear_interp_H(K); 
  if (isnan(H_value))
    fprintf(stderr, " h_eval %e %e %e %e %e\n", nu, Thetae, nuc, Bmag, H_value);
  return H_value;
}


#define GFAC (2.*M_PI*ME*CL/EE)
double G_eval(double Thetae, double Bmag, double nu, double kappa)
{
  double w = kappa_w(Thetae, kappa);
	double L = GFAC*nu/(Bmag* w*kappa * w*kappa); 

	if (L > LMAX) {
		return 0.;
	} else if (L < LMIN) {
	  return 0.;
  } else {

		return linear_interp_G(L, kappa);
	}
}

#undef KFAC
#undef KMIN
#undef KMAX
#undef GFAC
#undef LMIN
#undef LMAX
#undef EPSABS
#undef EPSREL

double linear_interp_K2(double Thetae)
{
	int i;
	double di, lT;

	lT = log(Thetae);

	di = (lT - lT_min)/dlT;
	i = (int) di;
	di = di - i;

	return exp((1. - di) * _K2[i] + di * _K2[i + 1]);
}

double linear_interp_F(double K)
{
	int i;
	double di, lK;

	lK = log(K);

	di = (lK - lK_min)/dlK;
	i = (int) di;
	di = di - i;

	return exp((1. - di) * _F[i] + di * _F[i + 1]);
}

double linear_interp_G(double L, double kappa)
{
  // get energy bin
	double lL = log(L);
	double di = (lL - lL_min)/dlL;
	int i = (int) di;
	di = di - i;

#if MODEL_EDF==EDF_KAPPA_FIXED

  return exp((1. - di) * _G[0][i] + di * _G[0][i + 1]);

#elif MODEL_EDF==EDF_KAPPA_VARIABLE

  // get kappa bin
  double dj = (kappa - KAPPA_MIN)/DKAPPA;
  int j = (int) dj;
  dj = dj - j;

  double v_below = exp((1. - di) * _G[j][i] + di * _G[j][i + 1]);
  double v_above = exp((1. - di) * _G[j + 1][i] + di * _G[j + 1][i + 1]);

  return (1. - dj) * v_below + dj * v_above;

#endif

  assert(1==0);  // shouldn't call into this table if not using a kappa model
  return 0.;

}

double linear_interp_H(double L)
{
  int i;
  double di, lL;

  lL = log(L);

  di = (lL - lL_min)/dlL;
  i = (int) di;
  di = di - i;

  return exp((1. - di) * _H[i] + di * _H[i + 1]);
}

