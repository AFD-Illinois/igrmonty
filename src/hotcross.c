#include <assert.h>

#include "decs.h"
#include "hotcross.h"
#include "model_radiation.h"
#include "tablecache.h"

/*

   given energy of photon in fluid rest frame w, in units of electron rest mass
   energy, and temperature of plasma, again in electron rest-mass units, return hot
   cross section in cgs.

   This has been checked against Wienke's Table 1, with some disagreement at
   the one part in 10^{-3} level, see wienke_table_1 in the subdirectory hotcross.
   It is not clear what this is due to, but Table 1 does appear to have been evaluated
   using Monte Carlo integration (!).

   A better way to do this would be to make a table in w*thetae and w/thetae; most
   of the variation is accounted for by w*thetae.

*/

#define gamma_max (1000.)
#define MINW  1.e-12
#define MAXW  1.e15
#define MINT  0.0001
#define MAXT  1.e4
#define NW  220
#define NT  80

#if MODEL_EDF==EDF_KAPPA_VARIABLE
double table[KAPPA_NSAMP][NW + 1][NT + 1];
#else
double table[1][NW + 1][NT + 1];
#endif

double dlw, dlT, lminw, lmint;

double kappa_function_int(double beta, void *params);
double hc_klein_nishina(double we);
double boostcross(double w, double mue, double gammae);

typedef struct int_rpars_struct {
  double Thetae;
  radiation_params *rpars;
} int_rpars;

// always recompute this. slightly slower than if saved, but safer
// since we could switch eDF
void init_hotcross(void)
{
  dlw = log10(MAXW / MINW) / NW;
  dlT = log10(MAXT / MINT) / NT;
  lminw = log10(MINW);
  lmint = log10(MINT);

  fprintf(stderr, "table for compton cross section... ");

#if MODEL_EDF==EDF_KAPPA_VARIABLE

  size_t dims[3] = { KAPPA_NSAMP, NW+1, NT+1 };
  double start[3] = { KAPPA_MIN, lminw, lmint };
  double dx[3] = { DKAPPA, dlw, dlT };

  int version = 1;
  // if you change the FORMAT of the table, please increment the version number. otherwise,
  // the table loader automatically recognizes changes to the table's shape, start, and dx.
  if (load_table("hotcross_edf_kappa_var.h5", version, table, 3, dims, start, dx) != 0) {

    fprintf(stderr, "generating... ");

#pragma omp parallel for
    for (int k = 0; k < KAPPA_NSAMP; ++k) {
      for (int j = 0; j <= NT; j++) {
        radiation_params rpars;
        rpars.kappa = KAPPA_MIN + k * DKAPPA;
        double lT = lmint + j * dlT;
        double norm = getnorm_dNdg(pow(10., lT), &rpars);
        for (int i = 0; i <= NW; i++) {
          double lw = lminw + i * dlw;
          double value = total_compton_cross_num(pow(10.,lw),pow(10.,lT),norm,&rpars);
          // note: this table is in w and *thetae*
          table[k][i][j] = log10(value);
          if (isnan(table[k][i][j])) {
            fprintf(stderr, "%d %d %g %g\n", i, j, lw, lT);
            exit(0);
          }
        }
      }
    }
    write_table("hotcross_edf_kappa_var.h5", version, table, 3, dims, start, dx);

  } else {
    fprintf(stderr, "loading from file... ");
  }

#else

#pragma omp parallel for
  for (int j = 0; j <= NT; j++) {
    radiation_params rpars;
    rpars.kappa = model_kappa;
    double lT = lmint + j * dlT;
    double norm = getnorm_dNdg(pow(10., lT), &rpars);
    for (int i = 0; i <= NW; i++) {
      double lw = lminw + i * dlw;
      double value = total_compton_cross_num(pow(10.,lw),pow(10.,lT),norm,&rpars);
      // note: this table is in w and *thetae*
      table[0][i][j] = log10(value);
      if (isnan(table[0][i][j])) {
        fprintf(stderr, "%d %d %g %g\n", i, j, lw, lT);
        exit(0);
      }
    }
  }

#endif

  fprintf(stderr, "done.\n");
}



double total_compton_cross_lkup(double w, double thetae, radiation_params *rpars)
{
  int i, j;
  double lw, lT, di, dj;

  // cold/low-energy: just use thomson cross section
  if (w * thetae < 1.e-6) {
    return SIGMA_THOMSON;
  }

  // cold, but possible high energy photon: use klein-nishina
  if (thetae < MINT) {
    return hc_klein_nishina(w) * SIGMA_THOMSON;
  }

  // in-bounds for table ... do bilinear interpolation
  if ((w > MINW && w < MAXW) && (thetae > MINT && thetae < MAXT)) {
#if MODEL_EDF==EDF_KAPPA_VARIABLE
  if (rpars->kappa > KAPPA_MIN) {
#endif

    lw = log10(w);
    lT = log10(thetae);
    i = (int) ((lw - lminw) / dlw);
    j = (int) ((lT - lmint) / dlT);
    di = (lw - lminw) / dlw - i;
    dj = (lT - lmint) / dlT - j;

#if MODEL_EDF==EDF_KAPPA_VARIABLE

    // get kappa bin
    double dk = (rpars->kappa - KAPPA_MIN)/DKAPPA;
    int k = (int) fmin(dk, KAPPA_NSAMP-2);
    dk = fmin(dk - k, 1.0);

    double lc1 = (1.-di) * (1.-dj) * table[k][i][j]
           + di * (1.-dj) * table[k][i+1][j]
           + (1.-di) * dj * table[k][i][j+1]
           + di * dj * table[k][i+1][j+1];
    double lc2 = (1.-di) * (1.-dj) * table[k+1][i][j]
           + di * (1.-dj) * table[k+1][i+1][j]
           + (1.-di) * dj * table[k+1][i][j+1]
           + di * dj * table[k+1][i+1][j+1];

    double lcross = (1. - dk) * lc1 + dk * lc2;

#else

    double lcross = (1.-di) * (1.-dj) * table[0][i][j]
           + di * (1.-dj) * table[0][i+1][j] 
           + (1.-di) * dj * table[0][i][j+1] 
           + di * dj * table[0][i+1][j+1];

#endif

    if (isnan(lcross)) {
      fprintf(stderr, "%g %g %d %d %g %g\n", lw, lT, i, j, di, dj);
    }

    return pow(10., lcross);
#if MODEL_EDF==EDF_KAPPA_VARIABLE
  } else {fprintf(stderr, "kappa < MIN: %g\n", rpars->kappa); exit(-1);}
#endif
  }

  //fprintf(stderr, "out of bounds: %g %g\n", w, thetae);

  return total_compton_cross_num(w, thetae, getnorm_dNdg(thetae, rpars), rpars);

}

#define MAXGAMMA  12.
#define DMUE    0.05
#define DGAMMAE   0.05

double total_compton_cross_num(double w, double thetae, double norm, radiation_params *rpars)
{
  double dmue, dgammae, mue, gammae, f, cross;

  if (isnan(w)) {
    fprintf(stderr, "compton cross isnan: %g %g\n", w, thetae);
    return 0.;
  }

  // check for easy-to-do limits
  if (thetae < MINT && w < MINW) {
    return SIGMA_THOMSON;
  }

  if (thetae < MINT) {
    return hc_klein_nishina(w) * SIGMA_THOMSON;
  }

  dmue = DMUE;
  dgammae = thetae * DGAMMAE;

  // integrate over mu_e, gamma_e, where mu_e is the cosine of the
  // angle between k and u_e, and the angle k is assumed to lie,
  // wlog, along the z axis
  cross = 0.;
  for (mue = -1. + dmue/2.; mue < 1.; mue += dmue)
    for (gammae= 1. + dgammae/2; gammae < 1. + MAXGAMMA*thetae; gammae += dgammae) {

      f = 0.5 * norm*dNdgammae(thetae, gammae, rpars);
      cross += dmue * dgammae * boostcross(w, mue, gammae) * f;

      if (isnan(cross)) {
        fprintf(stderr, "cross is nan %g %g %g %g %g %g\n", w,
          thetae, mue, gammae,
          dNdgammae(thetae, gammae, rpars),
          boostcross(w, mue, gammae));
      }
    }

  return cross * SIGMA_THOMSON;
}

// use beta instead of gamma, because its nicely bounded between 0 and 1.
double dNdg_integrand(double beta, void *params) 
{
  double gammae = exp(beta); // integrating in log space

  int_rpars *irpars = (int_rpars *)params;

  return  gammae * dNdgammae(irpars->Thetae, gammae, irpars->rpars);
}

double dNdgammae_powerlaw(double thetae, double gammae)
{
   double p = powerlaw_p;
   double gmin = powerlaw_gamma_min;
   double gmax = powerlaw_gamma_max;

   double exp_cutoff = exp(-gammae / gamma_max);
   (void)exp_cutoff;

   if (gammae < gmin || gmax < gammae) return 0.;

  // note no exponential cutoff. this means we're not using powerlaw_gamma_cut
  // or gamma_max. this choice makes normalization easier and seems consistent
  // with the symphony emissivity formula
  return (p-1) * pow(gammae, -p) / ( pow(gmin, 1-p) - pow(gmax, 1-p) );
}

double dNdgammae_kappa(double thetae, double gammae, double kappa)
{
  double exp_cutoff = exp(-gammae / gamma_max);

  // determine w by finding effective w for total energy match
  // (thermal MJ) at Thetae
  double w = kappa_w(thetae, kappa);
  
  return gammae * 
         sqrt(gammae*gammae - 1.) * 
         pow(1. + (gammae - 1.)/(kappa*w), -(kappa+1.)) * 
         exp_cutoff;
}

double getnorm_dNdg(double thetae, radiation_params *rpars)
{
#if (MODEL_EDF==EDF_KAPPA_FIXED) || (MODEL_EDF==EDF_KAPPA_VARIABLE)
  
  int_rpars irpars;
  irpars.Thetae = thetae;
  irpars.rpars = rpars;

  double result, error;
  gsl_function F;
  F.function = *dNdg_integrand;
  F.params = &irpars;

  double absolute_error = 0.;
  double relative_error = 1.e-6;
  size_t limit = 5000;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(5000);
  gsl_integration_qag(&F, 
                      0, log(1. + 100 * thetae), 
                      absolute_error, relative_error, limit, 
                      GSL_INTEG_GAUSS61,
                      w, &result, &error);
  gsl_integration_workspace_free(w);

  return 1. / result;

#elif MODEL_EDF==EDF_POWER_LAW

  return 1.;
  (void)thetae;  // silence unused parameter warning

#elif MODEL_EDF==EDF_MAXWELL_JUTTNER

  return 1.;
  (void)thetae;  // silence unused parameter warning

#else

  fprintf(stderr, "must select valid MODEL_EDF\n");
  exit(3);

#endif
}

double dNdgammae(double thetae, double gammae, radiation_params *rpars)
{
#if (MODEL_EDF==EDF_KAPPA_FIXED) || (MODEL_EDF==EDF_KAPPA_VARIABLE)

  return dNdgammae_kappa(thetae, gammae, rpars->kappa);

#elif MODEL_EDF==EDF_POWER_LAW

  return dNdgammae_powerlaw(thetae, gammae);

#elif MODEL_EDF==EDF_MAXWELL_JUTTNER

  // multiply K2(1/Thetae) by e^(1/Thetae) for numerical purposes
  double K2f;
  if (thetae > 1.e-2) {
    K2f = gsl_sf_bessel_Kn(2, 1. / thetae) * exp(1. / thetae);
  } else {
    K2f = sqrt(M_PI * thetae / 2.);
  }

  return (gammae * sqrt(gammae * gammae - 1.) / (thetae * K2f)) *
    exp(-(gammae - 1.) / thetae);

#else

  fprintf(stderr, "must select valid MODEL_EDF\n");
  exit(3);

#endif
}

double boostcross(double w, double mue, double gammae)
{
  double we, boostcross, v;

  // energy in electron rest frame 
  v = sqrt(gammae * gammae - 1.) / gammae;
  we = w * gammae * (1. - mue * v);

  boostcross = hc_klein_nishina(we) * (1. - mue * v);

  if (boostcross > 2) {
    fprintf(stderr, "w,mue,gammae: %g %g %g\n", w, mue,
      gammae);
    fprintf(stderr, "v,we, boostcross: %g %g %g\n", v, we,
      boostcross);
    fprintf(stderr, "kn: %g %g %g\n", v, we, boostcross);
  }

  if (isnan(boostcross)) {
    fprintf(stderr, "isnan: %g %g %g\n", w, mue, gammae);
    exit(0);
  }

  return boostcross;
}

double hc_klein_nishina(double we)
{
  double sigma;

  if (we < 1.e-3)
    return 1. - 2. * we;

  sigma = (3. / 4.) * (2. / (we * we) +
           (1. / (2. * we) -
            (1. + we) / (we * we * we)) * log(1. + 2. * we) +
           (1. + we) / ((1. + 2. * we) * (1. + 2. * we))
      );

  return sigma;

}

