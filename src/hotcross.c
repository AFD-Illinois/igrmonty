#include <assert.h>

#include "decs.h"
#include "hotcross.h"
#include "model_radiation.h"
#include "tablecache.h"
#include "electron_distribution.h"
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
#define NW 220
#define NT 80
#define NA 2
#define NXI 2

#if MODEL_EDF==EDF_KAPPA_VARIABLE
double table[KAPPA_NSAMP][NW + 1][NT + 1];
#else
double table[1][NW + 1][NT + 1];
#endif

double ani_table[NW + 1][NT + 1][NA + 1][NXI + 1];

double dlw, dlT, lminw, lmint;
double lminA, lmaxA, dlA, minxi, maxxi, dxi;

double kappa_function_int(double beta, void *params);
double hc_klein_nishina(double we);
double boostcross(double w, double mue, double gammae);

typedef struct int_rpars_struct {
  double Thetae;
  radiation_params *rpars;
} int_rpars;

// function declarations for 3d quadrature integration
double integral_z(double z, void *params);
double integral_y(double y, void *params);
double integral_x(double x, void *params);
double integrate_3D(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, double p1, double p2, double p3, double p4, double p5);

// always recompute this. slightly slower than if saved, but safer
// since we could switch eDF
void init_hotcross(void)
{
  // do the entire table in 4 D for the anisotropic case if anisotropy is enabled. Only for thermal models for now.
  if(anisotropy){

    if (debug){
      double photon_energy = pow(10,-10);
      double thetae_perp = 1e4;
      double A = 1.0e0;
      double xi = 2.356119;
      double ne = 1.000000e+00;
      radiation_params *rpars;
      clock_t t1 = clock();
      double cross= compute_hotcross_anisotropic(photon_energy, A, xi, thetae_perp, ne);
      clock_t t2 = clock();
      fprintf(stderr,"hotcross for photon_energy: %e, A: %e, xi: %e, thetae_perp: %e, ne: %e is %e\n", photon_energy, A, xi, thetae_perp, ne, cross);
      printf("Time taken: %f seconds\n", (double)(t2 - t1) / CLOCKS_PER_SEC);
      double cross_iso = total_compton_cross_num(photon_energy, thetae_perp, 1.0, rpars);
      fprintf(stderr,"isotropic hotcross for photon_energy: %e, thetae_perp: %e, ne: %e is %e\n", photon_energy, thetae_perp, ne, cross_iso);
      fprintf(stderr, "relative difference is: %e\n", (cross - cross_iso)/cross_iso);
      exit(0);
    }

    size_t rank=4;
    size_t dims[4] = {NW + 1, NT + 1, NA + 1, NXI + 1};

    dlw = log10(MAXW / MINW) / NW;
    dlT = log10(MAXT / MINT) / NT;
    lminw = log10(MINW);
    lmint = log10(MINT);
    // anisotropy parameters A (ratio of temperatures) and xi (pitch angle of photon wrt local B field)
    lminA = -1.0;
    lmaxA = 1.0;
    dlA = (lmaxA - lminA) / NA;
    minxi = 0.;
    maxxi = M_PI-1e-4;
    dxi = (maxxi - minxi) / NXI;

    double start[4] = {lminw, lmint, lminA, minxi};
    double dx[4] = {dlw, dlT, dlA, dxi};
    
    fprintf(stderr, "table for anisotropic compton cross section... ");
    if(load_table("hotcross_anisotropic.h5", 1, ani_table, rank, dims, start, dx) != 0){
      fprintf(stderr, "generating table...\n");
      #pragma omp parallel for schedule(dynamic,2) collapse(2)
      for (int i = 0; i <= NT; i++) {
        for (int j = 0; j <= NW; j++) {
          if (j==0) fprintf(stderr, "%ld ", i);
          // int k=0;
          // int l=0;
          for (int k = 0; k <= NA; k++) {
            for (int l = 0; l <= NXI; l++) {
              double lT = lmint + i * dlT;
              double lw = lminw + j * dlw;
              double lA = lminA + k * dlA;
              double xi = minxi + l * dxi;
              double value = compute_hotcross_anisotropic(pow(10., lw), pow(10., lA), xi , pow(10., lT), 1.0);
              // fprintf(stderr,"value: %e\n", value);exit(0);
              // note: this table is in w, *thetae*, A and xi
              ani_table[j][i][k][l] = log10(value);
              // ani_table[j][i][k+1][l] = log10(value);
              // ani_table[j][i][k][l+1] = log10(value);
              // ani_table[j][i][k+1][l+1] = log10(value);
              if (isnan(ani_table[j][i][k][l]) || value==0) {
                fprintf(stderr, "lw%g lT%g lA%g xi%g\n", lw, lT, lA, xi);
                exit(0);
              }
            }
          }
        }
      }
      write_table("hotcross_anisotropic.h5", 1, ani_table, rank, dims, start, dx);
    }
    else{
      fprintf(stderr, "loading from file... ");
    }
    fprintf(stderr, "done.\n");
    return;
  }

  if (debug){
    double photon_energy = pow(10,-3);
    double thetae_perp = 1e4;
    double ne = 1.000000e+00;
    clock_t t1 = clock();
    radiation_params rpars;
    double cross= total_compton_cross_lkup(photon_energy, thetae_perp, &rpars);
    clock_t t2 = clock();
    fprintf(stderr,"hotcross for photon_energy: %e, thetae_perp: %e, ne: %e is %e\n", photon_energy, thetae_perp, ne, cross);
    printf("Time taken: %f seconds\n", (double)(t2 - t1) / CLOCKS_PER_SEC);
    exit(0);
  }
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
  write_table("hotcross_isotropic.h5",1,table,3,(size_t[]){1,NW+1,NT+1},(double[]){0.0,lminw,lmint},(double[]){0.0,dlw,dlT});

#endif

  fprintf(stderr, "done.\n");
}

// similar to isotropic hotcross lookup but now with 4 parameters, limits on w and thetae (here thetae_perp) are the same, might need to change for large anisotropy cases
double total_compton_cross_lkup_anisotropic(double w, double thetae, double A, double xi)
{
  int i, j, k, l;
  double lw, lT, lA, di, dj, dk, dl;

  double lc1;

  // cold/low-energy: just use thomson cross section
  if (w * thetae < 1.e-6) {
    return SIGMA_THOMSON;
  }

  // cold, but possible high energy photon: use klein-nishina
  if (thetae < MINT) {
    return hc_klein_nishina(w) * SIGMA_THOMSON;
  }

  // in-bounds for table ... do tetralinear interpolation
  // if ((w > MINW && w < MAXW) && (thetae > MINT && thetae < MAXT) && (A > pow(10., lminA) && A < pow(10., lmaxA)) && (xi > minxi && xi < maxxi)) {
  if ((w > MINW && w < MAXW) && (thetae > MINT && thetae < MAXT)) {
    lw = log10(w);
    lT = log10(thetae);
    lA = log10(A);
    i = (int) ((lw - lminw) / dlw);
    j = (int) ((lT - lmint) / dlT);
    k = (int) ((lA - lminA) / dlA);
    l = (int) ((xi - minxi) / dxi);
    // k=0;
    // l=0;
    di = (lw - lminw) / dlw - i;
    dj = (lT - lmint) / dlT - j;
    dk = (lA - lminA) / dlA - k;
    dl = (xi - minxi) / dxi - l;
    // dk=1;
    // dl=1;

    lc1 = (1.-di) * (1.-dj) * (1.-dk) * (1.-dl) * ani_table[i][j][k][l]
           + di * (1.-dj) * (1.-dk) * (1.-dl) * ani_table[i+1][j][k][l]
           + (1.-di) * dj * (1.-dk) * (1.-dl) * ani_table[i][j+1][k][l]
           + di * dj * (1.-dk) * (1.-dl) * ani_table[i+1][j+1][k][l]
           + (1.-di) * (1.-dj) * dk * (1.-dl) * ani_table[i][j][k+1][l]
           + di * (1.-dj) * dk * (1.-dl) * ani_table[i+1][j][k+1][l]
           + (1.-di) * dj * dk * (1.-dl) * ani_table[i][j+1][k+1][l]
           + di * dj * dk * (1.-dl) * ani_table[i+1][j+1][k+1][l]
           + (1.-di) * (1.-dj) * (1.-dk) * dl * ani_table[i][j][k][l+1]
           + di * (1.-dj) * (1.-dk) * dl * ani_table[i+1][j][k][l+1]
           + (1.-di) * dj * (1.-dk) * dl * ani_table[i][j+1][k][l+1]
           + di * dj * (1.-dk) * dl * ani_table[i+1][j+1][k][l+1]
           + (1.-di) * (1.-dj) * dk * dl * ani_table[i][j][k+1][l+1]
           + di * (1.-dj) * dk * dl * ani_table[i+1][j][k+1][l+1]
           + (1.-di) * dj * dk * dl * ani_table[i][j+1][k+1][l+1]
           + di * dj * dk * dl * ani_table[i+1][j+1][k+1][l+1];
    // lc1 = (1.-di) * (1.-dj) * ani_table[i][j][0][0]
    //     + di * (1.-dj) * ani_table[i+1][j][0][0] 
    //     + (1.-di) * dj * ani_table[i][j+1][0][0] 
    //     + di * dj * ani_table[i+1][j+1][0][0];

    if(isnan(lc1)){
      fprintf(stderr,"lc1 is nan\n");
      fprintf(stderr, "%g %g %g %g %d %d %d %d %g %g %g %g\n", lw, lT, lA, xi, i, j, k, l, di, dj, dk, dl);
    }
    return pow(10., lc1);
  }
  // otherwise just compute numerically
  return compute_hotcross_anisotropic(w, A, xi, thetae, 1.0);
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
  if (rpars->kappa >= KAPPA_MIN) {
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
#if MODEL_EDF==EDF_KAPPA_VARIA
          // }BLE
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

// anisotropic hotcross function definitions

// Function to check scattering limits
double check_scattering_limits(double photon_energy, double thetae) {
    if (thetae < MINT && photon_energy < MINW) {
        return SIGMA_THOMSON;
    } else if (thetae < MINT) {
        return hc_klein_nishina(photon_energy);
    } else {
        return -1.0;
    }
}

/**
 * Collect all terms of the integrand for hotcross integration of bi-maxwellian distributions.
 *
 * Parameters:
 *   p_par        - Parallel component of the electron momentum
 *   p_perp       - Perpendicular component of the electron momentum
 *   phi          - Azimuth of the electron momentum with respect to the photon k
 *   photon_energy - k^0 in the plasma rest frame (in units of electron rest mass energy)
 *   A            - Anisotropy factor
 *   xi           - Pitch angle of the photon
 *   thetae_perp  - Perpendicular temperature of the electron distribution
 *   ne           - Number density of electrons
 * Returns:
 *   The integrand for hotcross integration.
 */

double hotcross_integrand_bimaxwell(double p_par, double p_perp, double phi, double photon_energy, double A, double xi, double thetae_perp, double ne) {
    double psq = p_perp * p_perp + p_par * p_par;
    // if psq is zero no need to compute
    if (psq == 0) {
        return 0;
    }
    double gammae = sqrt(psq + 1);
    double beta = sqrt(1 - 1 / (gammae * gammae));
    double p = sqrt(psq);
    double k[3] = {cos(xi), sin(xi), 0}; // Direction of photon in plasma frame
    //choice of coordinates [z,r,phi]: the plane spanning k and B is phi=0, T_par and T_perp form remaining cylindrical axes. Assumed k is normalized to 1 here.
    double mu_photon = (p_perp * k[1] * cos(phi) + p_par * k[0]) / p;

    double boost_factor =  (1 - mu_photon * beta);
    double sigma_kn_eframe = hc_klein_nishina(photon_energy * gammae * boost_factor);
    double distr_fn_val = dnd3p_bimaxwell_fast(A, p, ne, thetae_perp, p_par, p_perp);
    if (isnan(boost_factor)){
      fprintf(stderr, "boost_factor is nan. variables p_par: %f, p_perp: %f, phi: %f\n", p_par, p_perp, phi);
      exit(0);
    }
    if (isnan(sigma_kn_eframe)){
      fprintf(stderr, "hc_klein_nishina is nan. variables p_par: %f, p_perp: %f, phi: %f\n", p_par, p_perp, phi);
      exit(0);
    }
    if (isnan(distr_fn_val)){
      fprintf(stderr, "dnd3p_bimaxwell is nan. variables p_par: %f, p_perp: %f, phi: %f\n", p_par, p_perp, phi);
      exit(0);
    }
    // factor of p_perp here comes from the Jacobian of the transformation to cylindrical coordinates
    return boost_factor * sigma_kn_eframe * distr_fn_val * p_perp;
}

// Function that does 3D trapezoidal integration of the cross section for an anisotropic edf
double tpltrap(double photon_energy, double A, double xi, double thetae_perp, double ne, double p_perp_max, double p_par_max)
{
  double result = 0.0;
  int num_steps = 200;
  double phi_step = 2 * M_PI / 100;
  // double p_perp_step = (sqrt(MAXGAMMA*thetae_perp) - 1)/(MAXGAMMA*thetae_perp-1) * DGAMMAE*thetae_perp;
  // double p_par_step = (sqrt(MAXGAMMA*thetae_perp/A) - 1)/(MAXGAMMA*thetae_perp/A-1) * DGAMMAE*thetae_perp/A;

  // double p_perp_step = thetae_perp*DGAMMAE;
  // double p_par_step = thetae_perp*DGAMMAE/A;
  double p_perp_step = p_perp_max / num_steps;
  double p_par_step = p_par_max / num_steps;
  // fprintf(stderr,"p_perp_step: %f, p_par_step: %f\n", p_perp_step, p_par_step);
  // fprintf(stderr,"number of steps in p_perp: %e, p_par: %e\n", p_perp_max/p_perp_step, p_par_max/p_par_step);
  for (double phi = phi_step/2; phi < 2 * M_PI; phi += phi_step)
  {
    for (double p_perp = p_perp_step/2; p_perp < p_perp_max; p_perp += p_perp_step)
    {
      for (double p_par = -p_par_max + p_perp_step/2; p_par < p_par_max; p_par += p_par_step)
      {
        result += hotcross_integrand_bimaxwell(p_par, p_perp, phi, photon_energy, A, xi, thetae_perp, ne) *
                  phi_step * p_perp_step * p_par_step;
        if (isnan(result)) {
          fprintf(stderr, "result is nan. variables p_par: %f, p_perp: %f, phi: %f\n", p_par, p_perp, phi);
          exit(0);
        }
      }
    }
  }
  return result;
}

/**
 * Compute the hot cross section by numerical integration over the momentum space 
 * of the distribution function for a given set of parameters.
 *
 * For the anisotropic case, the integral is performed in cylindrical coordinates 
 * (p_perp, p_par, and phi), as a 3D integral is necessary due to the lack of exploitable 
 * symmetries in the cosine of mu_photon.
 *
 * Parameters:
 *   photon_energy - k^0 in the plasma rest frame (in units of electron rest mass energy)
 *   A             - Anisotropy factor
 *   xi            - Pitch angle of the photon
 *
 * Returns:
 *   sigma_hot - The hot cross section.
 */
double compute_hotcross_anisotropic(double photon_energy, double A, double xi, double thetae_perp, double ne) {
  double sigma = check_scattering_limits(photon_energy, thetae_perp);
  
  if (sigma != -1) {
    return sigma * SIGMA_THOMSON;
  }
  
  double p_perp_max = sqrt((1+MAXGAMMA*thetae_perp)*(1.+MAXGAMMA*thetae_perp) - 1.)/sqrt(2);
  double p_par_max = sqrt((1.+MAXGAMMA*thetae_perp/A)*(1.+MAXGAMMA*thetae_perp/A) - 1.)/sqrt(2);
  // double p_par_max = p_perp_max;
  // double p_perp_max = sqrt(MAXGAMMA*thetae_perp);
  // double p_par_max = sqrt(MAXGAMMA*thetae_perp/A);
  // fprintf(stderr,"p_perp_max: %f, p_par_max: %f\n", p_perp_max, p_par_max);
  // fprintf(stderr,"photon_energy: %e, A: %e, xi: %e, thetae_perp: %e, ne: %e\n", photon_energy, A, xi, thetae_perp, ne);
  // clock_t t1 = clock();
  // double result = tpltrap(photon_energy, A, xi, thetae_perp, ne, p_perp_max, p_par_max) * SIGMA_THOMSON;
  double result = integrate_3D(0,2*M_PI,0,p_perp_max,-p_par_max,p_par_max,photon_energy,A,xi,thetae_perp,ne)*SIGMA_THOMSON;
  result *= dnd3p_bimaxwell_prefactor(A, ne, thetae_perp);
  // clock_t t2 = clock();
    // printf("Time taken: %f seconds\n", (double)(t2 - t1) / CLOCKS_PER_SEC);
  return result;
}

// functions for quadrature integration in 3d for anisotropic integral
// here outermost integral is over x (phi), then y (p_perp), then z (p_par)
// 1D integral over z
// the abuse of the variable data is unfortunate, take care to make sure the parameters are being passed through each subfunction correctly.

double dummy_fn(double x, double y, double z, double p1, double p2, double p3, double p4, double p5) {
    return 1;
}

// integrand of z (p_par) integral
double integral_z(double z, void *params) {
    double *data = (double *)params;
    double x = data[0], y = data[1];
    double p1 = data[2], p2 = data[3], p3 = data[4], p4 = data[5], p5=data[6];
    // fprintf(stderr,"p1: %f, p2: %f, p3: %f, p4: %f, p5: %f\n", p1, p2, p3, p4, p5);exit(0);
    // return dummy_fn(z,y,x, p1, p2, p3, p4, p5);
    return hotcross_integrand_bimaxwell(z, y, x, p1, p2, p3, p4, p5);
}

// returns integrand of y (p_perp) integral i.e., integral over z(p_par)
double integral_y(double y, void *params) {
    gsl_integration_romberg_workspace *w = gsl_integration_romberg_alloc(10);
    double result, error;
    size_t neval;
    double *params_y = (double *)params;
    double x = params_y[0];
    double data[] = {x, y, params_y[3], params_y[4], params_y[5], params_y[6], params_y[7]};
    gsl_function F;
    F.function = &integral_z;
    F.params = data;
    // fprintf(stderr,"p1: %f, p2: %f, p3: %f, p4: %f, p5: %f\n", params_y[3], params_y[4], params_y[5], params_y[6], params_y[7]);exit(0);
    // gsl_integration_qagi(&F, 0, 1e-3, 1000, w, &result, &error);
    gsl_integration_romberg(&F, params_y[1], params_y[2], 0, 1e-4, &result, &neval, w);
    gsl_integration_romberg_free(w);
    return result;
}

// returns integrand of x (phi) integral i.e., integral over y(p_perp) and z(p_par)
double integral_x(double x, void *params) {
    gsl_integration_romberg_workspace *w = gsl_integration_romberg_alloc(10);
    double result, error;
    size_t neval;
    double *params_x = (double *)params;
    double data[] = {x, params_x[0], params_x[1], params_x[4], params_x[5], params_x[6], params_x[7], params_x[8]};
    gsl_function F;
    F.function = &integral_y;
    F.params = data;
    // fprintf(stderr,"p1: %f, p2: %f, p3: %f, p4: %f, p5: %f\n", params_x[4], params_x[5], params_x[6], params_x[7], params_x[8]);exit(0);
    gsl_integration_romberg(&F, params_x[2], params_x[3], 0, 1e-4, &result, &neval, w);
    // gsl_integration_qagiu(&F, params_x[2], 0, 1e-3, 1000, w, &result, &error);
    gsl_integration_romberg_free(w);
    return result;
}

// Compute the full 3D integral \int_x \int_y \int_z f(x,y,z) dx dy dz
// where f(x,y,z) is the function hotcross_integrand_bimaxwell
// x is phi, y is p_perp, z is p_par
// p1, p2, p3, p4, p5 are the parameters of the function hotcross_integrand_bimaxwell
double integrate_3D(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, double p1, double p2, double p3, double p4, double p5) {
    gsl_integration_romberg_workspace *w = gsl_integration_romberg_alloc(10);
    double result, error;
    size_t neval;
    // gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    // double result, error;
    gsl_function F;
    double limits[9] = {z_min, z_max, y_min, y_max, p1, p2, p3, p4, p5};
    // fprintf(stderr,"p1: %f, p2: %f, p3: %f, p4: %f, p5: %f\n", p1, p2, p3, p4, p5);exit(0);

    F.function = &integral_x;
    F.params = limits;
    gsl_integration_romberg(&F, x_min, x_max, 0, 1e-4, &result, &neval, w);
    gsl_integration_romberg_free(w);

    // gsl_integration_qag(&F, x_min, x_max, 0, 1e-3, 1000, 6, w, &result, &error);
    // gsl_integration_workspace_free(w);
    return result;
}
