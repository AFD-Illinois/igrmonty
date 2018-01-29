
#include "decs.h"

void coord(int i, int j, int k, double *X);
void get_fluid_zone(int i, int j, int k, double *Ne, double *Thetae, double *B,
        double Ucon[NDIM], double Bcon[NDIM]);

#define OLD_WGT (0)

void init_model(int argc, char *argv[])
{
  fprintf(stderr, "getting simulation data...\n");

  double Ntot;
  sscanf(argv[1], "%lf", &Ntot);
	Ns = (int) Ntot;

  // Read dumpfile
  init_data(argc, argv);

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
// HOW CAN THIS BE THREADSAFE?!?!???!?!?!?!?!??!?!?!
void make_super_photon(struct of_photon *ph, int *quit_flag)
{
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
}

#if OLD_WGT
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
  double jcst = M_SQRT2*EE*EE*EE/(27*ME*CL*CL);

  // THIS IS BROKEN WITH OPENMP!!
  //#pragma omp parallel
  {
    int lstart, lend, myid, nthreads;
    double Ne, Thetae, K2, B, fac, Ucon[NDIM], Bcon[NDIM];
    nthreads = omp_get_num_threads();
    myid = omp_get_thread_num();
    lstart = myid * (N_ESAMP / nthreads);
    lend = (myid + 1) * (N_ESAMP / nthreads);
    if (myid == nthreads - 1)
      lend = N_ESAMP + 1;

    //#pragma omp for collapse(3)
    ZLOOP {
        get_fluid_zone(i, j, k, &Ne, &Thetae, &B, Ucon, Bcon);
        if (Ne == 0. || Thetae < THETAE_MIN)
          continue;

        K2 = K2_eval(Thetae);
        fac = (jcst*Ne*B*Thetae*Thetae/K2)*sfac*geom[i][j].g;
        for (int l = lstart; l < lend; l++)
          sum[l] += fac*F_eval(Thetae, B, nu[l]);
      }
  } // omp parallel


  #pragma omp parallel for
  for (int i = 0; i <= N_ESAMP; i++)
    wgt[i] = log(sum[i]/(HPL*Ns) + WEIGHT_MIN);

  fprintf(stderr, "done.\n\n");
}
#else
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
  //double jcst = M_SQRT2*EE*EE*EE/(27*ME*CL*CL);

  // THIS IS BROKEN WITH OPENMP!!
  //#pragma omp parallel
  {
    int lstart, lend, myid, nthreads;
    double Ne, Thetae, B, Ucon[NDIM], Bcon[NDIM];
    nthreads = omp_get_num_threads();
    myid = omp_get_thread_num();
    lstart = myid * (N_ESAMP / nthreads);
    lend = (myid + 1) * (N_ESAMP / nthreads);
    if (myid == nthreads - 1)
      lend = N_ESAMP + 1;

    //#pragma omp for collapse(3)
    ZLOOP {
      get_fluid_zone(i, j, k, &Ne, &Thetae, &B, Ucon, Bcon);
      if (Ne == 0. || Thetae < THETAE_MIN)
        continue;

      for (int l = lstart; l < lend; l++) {
        sum[l] += int_jnu(Ne, Thetae, B, nu[l])*sfac*geom[i][j].g;
      }
    }
  } // omp parallel


  #pragma omp parallel for
  for (int i = 0; i <= N_ESAMP; i++)
    wgt[i] = log(sum[i]/(HPL*Ns) + WEIGHT_MIN);

  fprintf(stderr, "done.\n\n");
}
#endif

#define BTHSQMIN  (1.e-4)
#define BTHSQMAX  (1.e9)
#define NINT    (40000)
double lb_min, dlb;
double nint[NINT + 1];
double dndlnu_max[NINT + 1];
#if OLD_WGT
void init_nint_table(void)
{
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
      dn = F_eval(1., Bmag,
            exp(j * DLNU +
          LNUMIN)) / (exp(wgt[j]) + 1.e-100);
      if (dn > dndlnu_max[i])
        dndlnu_max[i] = dn;
      nint[i] += DLNU * dn;
    }

    nint[i] *= dx[1] * dx[2] * dx[3] * L_unit * L_unit * L_unit
        * M_SQRT2 * EE * EE * EE / (27. * ME * CL * CL)
        * 1. / HPL;
    nint[i] = log(nint[i]);
    dndlnu_max[i] = log(dndlnu_max[i]);
  }
}
#else
void init_nint_table(void)
{
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
}
#endif

#if OLD_WGT
void init_zone(int i, int j, int k, double *nz, double *dnmax)
{
  int l;
  double Ne, Thetae, Bmag, lbth;
  double dl, dn, ninterp, K2;
  double Ucon[NDIM], Bcon[NDIM];

  get_fluid_zone(i, j, k, &Ne, &Thetae, &Bmag, Ucon, Bcon);

  if (Ne == 0. || Thetae < THETAE_MIN) {
    *nz = 0.;
    *dnmax = 0.;
    return;
  }

  lbth = log(Bmag * Thetae * Thetae);

  dl = (lbth - lb_min) / dlb;
  l = (int) dl;
  dl = dl - l;
  if (l < 0) {
    *dnmax = 0.;
    *nz = 0.;
    return;
  } else if (l >= NINT) {

    fprintf(stderr,
      "warning: outside of nint table range %g...change in harm_utils.c\n",
      Bmag * Thetae * Thetae);
    fprintf(stderr,"%g %g %g %g\n",Bmag,Thetae,lbth,(lbth - lb_min)/dlb) ;
    ninterp = 0.;
    *dnmax = 0.;
    for (l = 0; l <= N_ESAMP; l++) {
      dn = F_eval(Thetae, Bmag,
            exp(j * DLNU +
          LNUMIN)) / exp(wgt[l]);
      if (dn > *dnmax)
        *dnmax = dn;
      ninterp += DLNU * dn;
    }
    ninterp *= dx[1] * dx[2] * dx[3] * L_unit * L_unit * L_unit
        * M_SQRT2 * EE * EE * EE / (27. * ME * CL * CL)
        * 1. / HPL;
  } else {
    if (isinf(nint[l]) || isinf(nint[l + 1])) {
      ninterp = 0.;
      *dnmax = 0.;
    } else {
      ninterp =
          exp((1. - dl) * nint[l] + dl * nint[l + 1]);
      *dnmax =
          exp((1. - dl) * dndlnu_max[l] +
        dl * dndlnu_max[l + 1]);
    }
  }

  K2 = K2_eval(Thetae);
  if (K2 == 0.) {
    *nz = 0.;
    *dnmax = 0.;
    return;
  }

  *nz = geom[i][j].g * Ne * Bmag * Thetae * Thetae * ninterp / K2;
  if (*nz > Ns * log(NUMAX / NUMIN)) {
    fprintf(stderr,
      "Something very wrong in zone %d %d: \ng = %g B=%g  Thetae=%g  K2=%g  ninterp=%g nz = %e\n\n",
      i, j, geom[i][j].g, Bmag, Thetae, K2, ninterp, *nz);
    exit(-1);
    *nz = 0.;
    *dnmax = 0.;
  }

  printf("%i %i %i nz = %e\n", i,j,k,*nz);
}
#else
void init_zone(int i, int j, int k, double *nz, double *dnmax)
{
  //int l;
  double Ne, Thetae, Bmag;
  double dn, ninterp, K2;
  double Ucon[NDIM], Bcon[NDIM];

  get_fluid_zone(i, j, k, &Ne, &Thetae, &Bmag, Ucon, Bcon);

  if (Ne == 0. || Thetae < THETAE_MIN) {
    *nz = 0.;
    *dnmax = 0.;
    return;
  }


/*
  double lbth = log(Bmag * Thetae * Thetae);
  double dl = (lbth - lb_min) / dlb;
  int l = (int) dl;
  dl = dl - l;
  if (l < 0) {
    *dnmax = 0.;
    *nz = 0.;
    return;
  } else if (l >= NINT) {

    fprintf(stderr,
      "warning: outside of nint table range %g...change in harm_utils.c\n",
      Bmag * Thetae * Thetae);
    fprintf(stderr,"%g %g %g %g\n",Bmag,Thetae,lbth,(lbth - lb_min)/dlb) ;
    ninterp = 0.;
    *dnmax = 0.;
    for (int m = 0; m <= N_ESAMP; m++) {
      dn = int_jnu(Ne, Thetae, Bmag,  exp(m*DLNU +LNUMIN))/exp(wgt[m]);
      if (dn > *dnmax)
        *dnmax = dn;
      ninterp += DLNU*dn;
    }
    ninterp *= dx[1]*dx[2]*dx[3]*L_unit*L_unit*L_unit;
  } else {
    if (isinf(nint[l]) || isinf(nint[l + 1])) {
      ninterp = 0.;
      *dnmax = 0.;
    } else {
      ninterp =
          exp((1. - dl) * nint[l] + dl * nint[l + 1]);
      *dnmax =
          exp((1. - dl) * dndlnu_max[l] +
        dl * dndlnu_max[l + 1]);
    }
  }*/
  
  ninterp = 0.;
  *dnmax = 0.;
  for (int m = 0; m <= N_ESAMP; m++) {
    double nu = exp(m*DLNU +LNUMIN);
    dn = int_jnu(Ne, Thetae, Bmag, nu)/(HPL*exp(wgt[m]));
    if (dn > *dnmax)
      *dnmax = dn;
    ninterp += DLNU*dn;
  }
  ninterp *= dx[1]*dx[2]*dx[3]*L_unit*L_unit*L_unit;

  /*K2 = K2_eval(Thetae);
  if (K2 == 0.) {
    *nz = 0.;
    *dnmax = 0.;
    return;
  }*/

  //*nz = geom[i][j].g * Ne * Bmag * Thetae * Thetae * ninterp / K2;
  *nz = geom[i][j].g*ninterp/omp_get_num_threads();
  if (*nz > Ns * log(NUMAX / NUMIN)) {
    fprintf(stderr,
      "Something very wrong in zone %d %d: \ng = %g B=%g  Thetae=%g  K2=%g  ninterp=%g nz = %e\n\n",
      i, j, geom[i][j].g, Bmag, Thetae, K2, ninterp, *nz);
    exit(-1);
    *nz = 0.;
    *dnmax = 0.;
  }

  //printf("%i %i %i nz = %e dnmax = %e\n", i,j,k,*nz,*dnmax);
  //exit(-1);
}
#endif

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
  /*zj++;
  if (zj >= N2) {
    zj = 0;
    zi++;
    if (zi >= N1) {
      in2gen = 1;
      *i = N1;
      return 1;
    }
  }*/
  
  init_zone(zi, zj, zk, &n2gen, dnmax);
  if (fmod(n2gen, 1.) > monty_rand()) {
    in2gen = (int) n2gen + 1;
  } else {
    in2gen = (int) n2gen;
  }

  *i = zi;
  *j = zj;
  *k = zk;

  return in2gen;
}

#if OLD_WGT
void sample_zone_photon(int i, int j, int k, double dnmax, struct of_photon *ph)
{
  double K_tetrad[NDIM], tmpK[NDIM], E, Nln;
  double nu, th, cth, sth, phi, sphi, cphi, jmax, weight;
  double Ne, Thetae, Bmag, Ucon[NDIM], Bcon[NDIM], bhat[NDIM];
  static double Econ[NDIM][NDIM], Ecov[NDIM][NDIM];

  coord(i, j, k, ph->X);

  Nln = LNUMAX - LNUMIN;

  get_fluid_zone(i, j, k, &Ne, &Thetae, &Bmag, Ucon, Bcon);

  // Sample from superphoton distribution in current simulation zone
  do {
    nu = exp(monty_rand() * Nln + LNUMIN);
    weight = linear_interp_weight(nu);
  } while (monty_rand() > (F_eval(Thetae, Bmag, nu) / weight) / dnmax);

  ph->w = weight;
  jmax = jnu(nu, Ne, Thetae, Bmag, M_PI / 2.);
  do {
    cth = 2. * monty_rand() - 1.;
    th = acos(cth);

  } while (monty_rand() >
     jnu(nu, Ne, Thetae, Bmag, th) / jmax);

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
}
#else
void sample_zone_photon(int i, int j, int k, double dnmax, struct of_photon *ph)
{
  double K_tetrad[NDIM], tmpK[NDIM], E, Nln;
  double nu, th, cth, sth, phi, sphi, cphi, jmax, weight;
  double Ne, Thetae, Bmag, Ucon[NDIM], Bcon[NDIM], bhat[NDIM];
  static double Econ[NDIM][NDIM], Ecov[NDIM][NDIM];

  coord(i, j, k, ph->X);

  Nln = LNUMAX - LNUMIN;

  get_fluid_zone(i, j, k, &Ne, &Thetae, &Bmag, Ucon, Bcon);

  // Sample from superphoton distribution in current simulation zone
  do {
    nu = exp(monty_rand() * Nln + LNUMIN);
    weight = linear_interp_weight(nu);
  } while (monty_rand() > (int_jnu(Ne, Thetae, Bmag, nu)/(HPL*weight))/dnmax);

  ph->w = weight;
  jmax = jnu(nu, Ne, Thetae, Bmag, M_PI / 2.);
  do {
    cth = 2. * monty_rand() - 1.;
    th = acos(cth);

  } while (monty_rand() >
     jnu(nu, Ne, Thetae, Bmag, th) / jmax);

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
}
#endif

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

