
#include "decs.h"

void get_fluid_zone(int i, int j, int k, double *Ne, double *Thetae, double *B,
        double Ucon[NDIM], double Bcon[NDIM]);

void init_model(int argc, char *argv[], Params *params)
{
  fprintf(stderr, "getting simulation data...\n");

  if (params->loaded) {

    Ns = (int) params->Ns;

  } else {

    report_bad_input(argc);

    double Ntot;
    sscanf(argv[1], "%lf", &Ntot);
	  Ns = (int) Ntot;

  }

  // Read dumpfile
  init_data(argc, argv, params);

  // make look-up table for hot cross sections
  init_hotcross();

  // make table for solid angle integrated emissivity and K2
  init_emiss_tables();

  // make table for superphoton weights
  init_weight_table();

  // Initialize random number generators
  init_monty_rand();

  // If using van Hoof 2015, initialize the electron-ion Gaunt factor spline
  #if  BREMSSTRAHLUNG == 3
  init_bremss_spline();
  #endif

}

int n2gen = -1;
double dnmax;
int zone_i, zone_j, zone_k;
void make_super_photon(struct of_photon *ph, int *quit_flag)
{
	#ifdef EMIT_ORIGIN
  if (n2gen < 0) {
    n2gen = Ns;
  }
  n2gen--;
  if (n2gen < 0) {
    *quit_flag = 1;
  }
  
  sample_origin_photon(ph);
  #else

  #pragma omp critical
  {
  if (zone_i != N1) {
    while (n2gen <= 0) {
      n2gen = get_zone(&zone_i, &zone_j, &zone_k, &dnmax);
    }
    n2gen--;
  }
  }

	if (zone_i == N1) {
		*quit_flag = 1;
  } else {
    sample_zone_photon(zone_i, zone_j, zone_k, dnmax, ph);
  }

  #endif // EMIT_ORIGIN
}

void init_weight_table(void)
{
  double sum[N_ESAMP+1], nu[N_ESAMP+1];
 
  fprintf(stderr, "Building table for superphoton weights\n");

  #pragma omp parallel for
  for (int i = 0; i <= N_ESAMP; i++) {
    sum[i] = 0.;
    nu[i] = exp(i * DLNU + LNUMIN);
  }

  double sfac = dx[1]*dx[2]*dx[3]*L_unit*L_unit*L_unit;

  #pragma omp parallel for shared(sum) collapse(3)
  ZLOOP {
    double Ne, Thetae, B, Ucon[NDIM], Bcon[NDIM];
    get_fluid_zone(i, j, k, &Ne, &Thetae, &B, Ucon, Bcon);

    if (Ne == 0.) continue;

    for (int l=0; l<N_ESAMP; ++l) {
     #pragma omp atomic
      sum[l] += int_jnu(Ne, Thetae, B, nu[l]) * sfac * geom[i][j].gzone;
    }
  }

  #pragma omp parallel for
  for (int i = 0; i <= N_ESAMP; i++)
    wgt[i] = log(sum[i]/(HPL*Ns) + WEIGHT_MIN);

  #pragma omp parallel for collapse(3) 
  ZLOOP {
    double Ne, Thetae, Bmag;
    double Ucon[NDIM], Bcon[NDIM];
    double ninterp = 0.;
    get_fluid_zone(i, j, k, &Ne, &Thetae, &Bmag, Ucon, Bcon);
    for (int m=0; m<=N_ESAMP; ++m) {
      ninterp += DLNU * int_jnu(Ne, Thetae, Bmag, exp(m*DLNU + LNUMIN)) / (HPL*exp(wgt[m]));
    }
    ninterp *= geom[i][j].g * dx[1]*dx[2]*dx[3] * L_unit*L_unit*L_unit;
    n2gens[i][j][k] = ninterp;
  }

  fprintf(stderr, "done.\n\n");
}

#define BTHSQMIN  (1.e-4)
#define BTHSQMAX  (1.e9)
#define NINT    (40000)
double lb_min, dlb;
double nint[NINT + 1];
double dndlnu_max[NINT + 1];

void init_zone(int i, int j, int k, double *nz, double *dnmax)
{
  //int l;
  double Ne, Thetae, Bmag;
  double dn, ninterp;
  double Ucon[NDIM], Bcon[NDIM];

  get_fluid_zone(i, j, k, &Ne, &Thetae, &Bmag, Ucon, Bcon);

  if (Ne == 0.) {// || Thetae < THETAE_MIN) {
    *nz = 0.;
    *dnmax = 0.;
    return;
  }

  *nz = n2gens[i][j][k];
  ninterp = *nz / geom[i][j].g / dx[1]/dx[2]/dx[3] / L_unit/L_unit/L_unit;

  if (*nz > Ns * log(NUMAX / NUMIN)) {
    fprintf(stderr,
      "Something very wrong in zone %d %d: \ng = %g B=%g  Thetae=%g  ninterp=%g nz = %e\n\n",
      i, j, geom[i][j].g, Bmag, Thetae, ninterp, *nz);
    exit(-1);
    *nz = 0.;
    *dnmax = 0.;
  }

  // *nz gives an idea of how many superphotons we want in this entire zone.
  // then *nz / N_ESAMP is the number of superphotons we want to create per
  // energy bin.
  for (int m=0; m<N_ESAMP; ++m) {
    double nu = exp(m*DLNU +LNUMIN);
    dn = int_jnu(Ne, Thetae, Bmag, nu)/(HPL*exp(wgt[m]));
    if (*nz > 0) {
      if (dn == 0) {
        zwgt[m] = 0.;
      } else {
        zwgt[m] = log( exp(wgt[m]) * dn * DLNU / ninterp * N_ESAMP );
      }
    } else {
      zwgt[m] = 0.;
    }
  }
}

int zone_flag;
static int zi = 0;
static int zj = 0;
static int zk = -1;

void reset_zones() 
{
  zi = 0;
  zj = 0;
  zk = -1;
  zone_i = 0;
  zone_j = 0;
  zone_k = -1;
}

int get_zone(int *i, int *j, int *k, double *dnmax)
{
  int in2gen;
  double n2gen;

  zone_flag = 1;
  zk++;
  if (zk >= N3) {
    zk = 0;
    zj++;
    if (zj >= N2) {
      zj = 0;
      zi++;
      if (zi >= N1) {
        in2gen = 1;
        *i = N1;
        return 1;
      }
    }
  }

  n2gen = n2gens[zi][zj][zk] * Ns_scale;
  if (fmod(n2gen, 1.) > monty_rand()) {
    in2gen = (int) n2gen + 1;
  } else {
    in2gen = (int) n2gen;
  }
 
  if (in2gen > 0) {
    init_zone(zi, zj, zk, &n2gen, dnmax);
  }

  *i = zi;
  *j = zj;
  *k = zk;

#ifdef MODEL_TRANSPARENT
  in2gen *= 10000;
  if (zi != 20) return 0;
  if (zj > 1 && zj < N2-2) return 0;
#endif

  return in2gen * Ns_scale;
}

#ifdef EMIT_ORIGIN
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
void sample_origin_photon(struct of_photon *ph)
{
  double K_tetrad[NDIM], tmpK[NDIM], E;//, Nln;
  double nu, /*th, */cth, sth, phi, sphi, cphi, /*jmax, */weight;
  //double Ne, Thetae, Bmag;//Ucon[NDIM], Bcon[NDIM], bhat[NDIM];
  static double Econ[NDIM][NDIM], Ecov[NDIM][NDIM];

  // Assume spherical coordinates
  ph->X[0] = 0.;
  ph->X[1] = 1.e-5;
  ph->X[2] = M_PI/2.;
  ph->X[3] = 0.;

  // Sample intensity uniformly in frequency
  nu = exp(monty_rand()*(LNUMAX - LNUMIN) + LNUMIN);
  weight = get_Inu(nu)/get_Imax();

  ph->w = weight;
  
  cth = 2.*monty_rand() - 1.;
  //th = acos(cth);
  sth = sqrt(1. - cth*cth);
  phi = 2.*M_PI*monty_rand();
  cphi = cos(phi);
  sphi = sin(phi);

  double th = 0.;
  phi = 0.;
  sth = sin(th);
  cth = cos(th);
  cphi = cos(phi);
  sphi = sin(phi);

  E = nu*HPL/(ME*CL*CL);
  K_tetrad[0] = E;
  K_tetrad[1] = E*cth;
  K_tetrad[2] = E*cphi*sth;
  K_tetrad[3] = E*sphi*sth;

  double Ucon[NDIM] = {1, 0, 0, 0};
  double ehat[NDIM] = {1, 0, 0, 1};
  double gcov[NDIM][NDIM];
  gcov_func(ph->X, gcov);
  make_tetrad(Ucon, ehat, gcov, Econ, Ecov);
 
  tetrad_to_coordinate(Econ, K_tetrad, ph->K);

  K_tetrad[0] *= -1.;
  tetrad_to_coordinate(Ecov, K_tetrad, tmpK);

  ph->E = ph->E0 = ph->E0s = -tmpK[0];
  ph->L = tmpK[3];
  ph->tau_scatt = 0.;
  ph->tau_abs = 0.;
  ph->X1i = ph->X[1];
  ph->X2i = ph->X[2];
  ph->nscatt = 0;
  ph->ne0 = 0.;
  ph->b0 = 0.;
  ph->thetae0 = 0.;
}
#endif // EMIT_ORIGIN

double zone_linear_interp_weight(double nu) {

  int i;
  double di, lnu;

  lnu = log(nu);
  di = (lnu - LNUMIN)/DLNU;
  i = (int)di;
  di = di - i;

  // intel compiler has issues if zwgt[i] = -inf 
  // and returns exp( EXPRESSION ) = -nan, so we
  // manually check here.
  if ( isinf(zwgt[i]) || isinf(zwgt[i+1]) ) return 0.;

  return exp( (1. - di)*zwgt[i] + di*zwgt[i + 1] ); 

  double wgt = exp( (1. - di)*zwgt[i] + di*zwgt[i + 1] );

  if ( isnan(wgt) ) return 0.;
  return wgt;
}

void sample_zone_photon(int i, int j, int k, double dnmax, struct of_photon *ph)
{
  double K_tetrad[NDIM], tmpK[NDIM], E, Nln;
  double nu, th, cth, sth, phi, sphi, cphi, jmax, weight;
  double Ne, Thetae, Bmag, Ucon[NDIM], Bcon[NDIM];

  ijktoX(i, j, k, ph->X);

  get_fluid_zone(i, j, k, &Ne, &Thetae, &Bmag, Ucon, Bcon);

#ifdef MODEL_TRANSPARENT

  // monochromatic
  if (lnumin == lnumax) {
    nu = pow(10., lnumin);
    ph->w = 1.e+40;
  }

  // power law spectrum
  else {
    double lnu = monty_rand() * (lnumax - lnumin) + lnumin;
    nu = pow(10., lnu);
    double numin = pow(10., lnumin);
    ph->w = 1.e+40 * pow(nu, alpha_spec) / pow(numin, alpha_spec); 
  }

  // isotropic emission
  cth = 2. * monty_rand() - 1.;

#else

  Nln = LNUMAX - LNUMIN;

  // Sample from superphoton distribution in current simulation zone
  nu = exp(monty_rand() * Nln + LNUMIN);
  weight = zone_linear_interp_weight(nu);

  ph->w = weight;
  jmax = jnu(nu, Ne, Thetae, Bmag, M_PI / 2.);
  do {
    cth = 2. * monty_rand() - 1.;
    th = acos(cth);
  } while (monty_rand() > jnu(nu, Ne, Thetae, Bmag, th) / jmax);

#endif

  sth = sqrt(1. - cth * cth);
  phi = 2. * M_PI * monty_rand();
  cphi = cos(phi);
  sphi = sin(phi);

  E = nu * HPL / (ME * CL * CL);
  K_tetrad[0] = E;
  K_tetrad[1] = E * cth;
  K_tetrad[2] = E * cphi * sth;
  K_tetrad[3] = E * sphi * sth;

  // This section only used if CUSTOM_AVG == 1 in custom.h. See that
  // file for more details.
  double blr, blh;
  bl_coord(ph->X, &blr, &blh);
  ph->QTY0 = blr;

  tetrad_to_coordinate(tetrads[i][j][k].Econ, K_tetrad, ph->K);

  K_tetrad[0] *= -1.;
  tetrad_to_coordinate(tetrads[i][j][k].Ecov, K_tetrad, tmpK);

  ph->E = ph->E0 = ph->E0s = -tmpK[0];
  ph->L = tmpK[3];
  ph->tau_scatt = 0.;
  ph->tau_abs = 0.;
  ph->X1i = ph->X[1];
  ph->X2i = ph->X[2];
  ph->nscatt = 0;
  ph->ne0 = Ne;
  ph->b0 = Bmag;
  ph->thetae0 = Thetae;
  ph->ratio_brems = jnu_ratio_brems(nu, Ne, Thetae, Bmag, th); // TODO uninitialized in transparent models
#ifdef TRACK_PH_CREATION
  ph->isrecorded = 0;
#endif // TRACK_PH_CREATION

}

void init_geometry()
{
  #pragma omp parallel for collapse(2)
  for (int i = 0; i < N1; i++) {
    for (int j = 0; j < N2; j++) {

      // save geometry for each zone (centered in zone). gcov/gcon/g should
      // be in geodesic coordinates, gzone should be in zone coordinates.
      double X[NDIM];
      geom[i][j].gzone = gdet_zone(i, j, 0);
      ijktoX(i, j, 0, X);
      gcov_func(X, geom[i][j].gcov);
      gcon_func(geom[i][j].gcov, geom[i][j].gcon);
      geom[i][j].g = gdet_func(geom[i][j].gcov);

    }
  }
}

void init_tetrads()
{
  #pragma omp parallel for collapse(3)
  for (int i=0; i<N1; ++i) {
    for (int j=0; j<N2; ++j) {
      for (int k=0; k<N3; ++k) {
        // precompute tetrads
        double Ne, Thetae, Bmag;
        double Ucon[NDIM], Bcon[NDIM], bhat[NDIM];

        get_fluid_zone(i, j, k, &Ne, &Thetae, &Bmag, Ucon, Bcon);

        if (Bmag > 0.) {
          for (int l = 0; l < NDIM; l++) {
            bhat[l] = Bcon[l] * B_unit / Bmag;
          }
        } else {
          for (int l = 1; l < NDIM; l++) {
            bhat[l] = 0.;
          }
          bhat[1] = 1.;
        }

        make_tetrad(Ucon, bhat, geom[i][j].gcov, tetrads[i][j][k].Econ, tetrads[i][j][k].Ecov);
      }
    }
  }
}




