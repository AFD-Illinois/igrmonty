
#include "decs.h"

void coord(int i, int j, int k, double *X);
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

  // Make table for quick evaluation of ns_zone
  init_nint_table();

  // Initialize random number generators
  init_monty_rand();
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
  while (n2gen <= 0) {
		n2gen = get_zone(&zone_i, &zone_j, &zone_k, &dnmax);
	}

	n2gen--;

	if (zone_i == N1)
		*quit_flag = 1;
	else
		*quit_flag = 0;

	if (*quit_flag != 1) {
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

  #pragma omp parallel for shared(sum) schedule(dynamic,1) collapse(3)
  ZLOOP {
    double Ne, Thetae, B, Ucon[NDIM], Bcon[NDIM];
    get_fluid_zone(i, j, k, &Ne, &Thetae, &B, Ucon, Bcon);

    if (Ne == 0.) continue;

    for (int l=0; l<N_ESAMP; ++l) {
     #pragma omp atomic
      sum[l] += int_jnu(Ne, Thetae, B, nu[l]) * sfac * geom[i][j].g;
    }
  }

  #pragma omp parallel for
  for (int i = 0; i <= N_ESAMP; i++)
    wgt[i] = log(sum[i]/(HPL*Ns) + WEIGHT_MIN);

  fprintf(stderr, "done.\n\n");
}

#define BTHSQMIN  (1.e-4)
#define BTHSQMAX  (1.e9)
#define NINT    (40000)
double lb_min, dlb;
double nint[NINT + 1];
double dndlnu_max[NINT + 1];

void init_nint_table(void)
{
  /*
  double Bmag, dn;
  static int firstc = 1;

  if (firstc) {
    lb_min = log(BTHSQMIN);
    dlb = log(BTHSQMAX / BTHSQMIN) / NINT;
    firstc = 0;
  }

  for (int i = 0; i <= NINT; i++) {
    nint[i] = 0.;
    Bmag = exp(i * dlb + lb_min);
    dndlnu_max[i] = 0.;
    for (int j = 0; j < N_ESAMP; j++) {
      dn = int_jnu(Ne_unit, 1., Bmag, exp(j*DLNU+LNUMIN))/(exp(wgt[j]) + 1.e-100);
      if (dn > dndlnu_max[i])
        dndlnu_max[i] = dn;
      nint[i] += DLNU*dn;
    }

    nint[i] *= dx[1]*dx[2]*dx[3]*L_unit*L_unit*L_unit;
    nint[i] = log(nint[i]);
    dndlnu_max[i] = log(dndlnu_max[i]);
  }
   */
}

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
  
  ninterp = 0.;
  *dnmax = 0.;
  for (int m = 0; m <= N_ESAMP; m++) {
    double nu = exp(m*DLNU +LNUMIN);
    dn = int_jnu(Ne, Thetae, Bmag, nu)/(HPL*exp(wgt[m]));
    if (dn > *dnmax) {
      *dnmax = dn;
    }
    ninterp += DLNU*dn;
  }

  *nz = geom[i][j].g * ninterp * dx[1]*dx[2]*dx[3] * L_unit*L_unit*L_unit;
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
        zwgt[m] = log(WEIGHT_MIN);
      } else {
        zwgt[m] = log( exp(wgt[m]) * dn * DLNU / ninterp * N_ESAMP + WEIGHT_MIN );
      }
    }
  }
}

int zone_flag;
int get_zone(int *i, int *j, int *k, double *dnmax)
{
  int in2gen;
  double n2gen;
  static int zi = 0;
  static int zj = 0;
  static int zk = -1;

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
  
  init_zone(zi, zj, zk, &n2gen, dnmax);
  if (fmod(n2gen, 1.) > monty_rand()) {
    in2gen = (int) n2gen + 1;
  } else {
    in2gen = (int) n2gen;
  }

  *i = zi;
  *j = zj;
  *k = zk;

#ifdef MODEL_TRANSPARENT
  in2gen *= 10000;
  if (zi != 20) return 0;
  if (zj > 1 && zj < N2-2) return 0;
#endif

  return in2gen;
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

  return exp((1. - di)*zwgt[i] + di*zwgt[i + 1]);
}

void sample_zone_photon(int i, int j, int k, double dnmax, struct of_photon *ph)
{
  double K_tetrad[NDIM], tmpK[NDIM], E, Nln;
  double nu, th, cth, sth, phi, sphi, cphi, jmax, weight;
  double Ne, Thetae, Bmag, Ucon[NDIM], Bcon[NDIM], bhat[NDIM];
  static double Econ[NDIM][NDIM], Ecov[NDIM][NDIM];

  coord(i, j, k, ph->X);

  

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

  /*
  if(E > 1.e-4) fprintf(stdout,"HOT: %d %d %g %g %g %g %g\n",
    i,j,E/(0.22*(EE*Bmag/(2.*M_PI*ME*CL))*(HPL/(ME*CL*CL))*Thetae*Thetae),
    ph->X[1],ph->X[2], Thetae,Bmag) ; 
  */

  if (zone_flag) { // First photon created in this zone, so make the tetrad
    if (Bmag > 0.) {
      for (int l = 0; l < NDIM; l++)
        bhat[l] = Bcon[l] * B_unit / Bmag;
    } else {
      for (int l = 1; l < NDIM; l++)
        bhat[l] = 0.;
      bhat[1] = 1.;
    }
    make_tetrad(Ucon, bhat, geom[i][j].gcov, Econ, Ecov);
    zone_flag = 0;
  }

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
  ph->ne0 = Ne;
  ph->b0 = Bmag;
  ph->thetae0 = Thetae;
  ph->ratio_brems = jnu_ratio_brems(nu, Ne, Thetae, Bmag, th);
#ifdef TRACK_PH_CREATION
  ph->isrecorded = 0;
#endif // TRACK_PH_CREATION
}

void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM])
{
  // Map X[3] into [0,stopx[3]), assuming startx[3] = 0
  double phi = fmod(X[3], stopx[3]);
  if (phi < 0.) phi = stopx[3] + phi;

  *i = (int)((X[1] - startx[1]) / dx[1] - 0.5 + 1000) - 1000;
  *j = (int)((X[2] - startx[2]) / dx[2] - 0.5 + 1000) - 1000;
  *k = (int)((X[3] - startx[3]) / dx[3] - 0.5 + 1000) - 1000;

  if (*i < 0) {
    *i = 0;
    del[1] = 0.;
  } else if (*i > N1 - 2) {
    *i = N1 - 2;
    del[1] = 1.;
  } else {
    del[1] = (X[1] - ((*i + 0.5) * dx[1] + startx[1])) / dx[1];
  }

  if (*j < 0) {
    *j = 0;
    del[2] = 0.;
  } else if (*j > N2 - 2) {
    *j = N2 - 2;
    del[2] = 1.;
  } else {
    del[2] = (X[2] - ((*j + 0.5) * dx[2] + startx[2])) / dx[2];
  }

  if (*k < 0) {
    *k = 0;
    del[3] = 0.;
  } else if (*k > N3 - 2) {
    *k = N3 - 2;
    del[3] = 1.;
  } else {
    del[3] = (X[3] - ((*k + 0.5)*dx[3] + startx[3]))/dx[3];
  }
}

/* return boyer-lindquist coordinate of point */
/*void bl_coord(double *X, double *r, double *th)
{

  *r = exp(X[1]) + R0;
  *th = M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);

  return;
}*/

void coord(int i, int j, int k, double *X)
{
  // Zone-centered coordinate values
  X[0] = startx[0];
  X[1] = startx[1] + (i + 0.5)*dx[1];
  X[2] = startx[2] + (j + 0.5)*dx[2];
  X[3] = startx[3] + (k + 0.5)*dx[3];
}

/*
void set_units(char *munitstr)
{
  double MBH;

  sscanf(munitstr, "%lf", &M_unit);

  // from this, calculate units of length, time, mass,
  //    and derivative units
  MBH = 4.6e6 * MSUN ;
  L_unit = GNEWT * MBH / (CL * CL);
  T_unit = L_unit / CL;

  fprintf(stderr, "\nUNITS\n");
  fprintf(stderr, "L,T,M: %g %g %g\n", L_unit, T_unit, M_unit);

  RHO_unit = M_unit / pow(L_unit, 3);
  U_unit = RHO_unit * CL * CL;
  B_unit = CL * sqrt(4. * M_PI * RHO_unit);

  fprintf(stderr, "rho,u,B: %g %g %g\n", RHO_unit, U_unit, B_unit);

  Ne_unit = RHO_unit / (MP + ME);

  max_tau_scatt = (6. * L_unit) * RHO_unit * 0.4;

  fprintf(stderr, "max_tau_scatt: %g\n", max_tau_scatt);

}*/

void init_geometry()
{
  int i, j;
  double X[NDIM];

  for (i = 0; i < N1; i++) {
    for (j = 0; j < N2; j++) {

      // Zone-centered, assume symmetry in X[3]
      coord(i, j, 0, X);

      gcov_func(X, geom[i][j].gcov);

      geom[i][j].g = gdet_func(geom[i][j].gcov);

      gcon_func(geom[i][j].gcov, geom[i][j].gcon);
    }
  }
}

