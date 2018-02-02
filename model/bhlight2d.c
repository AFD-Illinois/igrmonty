#include "decs.h"

#define NVAR (NPRIM)

// Grid functions
double ****bcon;
double ****bcov;
double ****ucon;
double ****ucov;
double ****p;
double ***ne;
double ***thetae;
double ***b;

double TP_OVER_TE;

static double game;
static double MBH;
  
//static int with_radiation;
//static int with_derefine_poles;
static int with_electrons;

void safe_fscanf(FILE *stream, const char *format, ...)
{
  va_list args;
  va_start(args, format);
  int vfscanfReturn = vfscanf(stream, format, args);
  va_end(args);
  if (vfscanfReturn == -1) {
    fprintf(stderr, "fscanf() call failed! Exiting!\n");
    exit(-1);
  }
}

void report_bad_input(int argc) 
{
  if (argc < 3) {
    fprintf(stderr, "usage: \n");
    fprintf(stderr, "  HARM:    grmonty Ns fname M_unit[g] MBH[Msolar] Tp/Te\n");
    fprintf(stderr, "  bhlight: grmonty Ns fname\n");
    exit(0);
  }
}

///////////////////////////////// SUPERPHOTONS /////////////////////////////////

#define RMAX  1000.
#define ROULETTE  1.e4
int stop_criterion(struct of_photon *ph)
{
  double wmin, X1min, X1max;

  // Stop if weight below minimum weight
  wmin = WEIGHT_MIN;

  // Stop at event horizon
  X1min = log(Rh);

  // Stop at large distance
  X1max = log(RMAX);

  if (ph->X[1] < X1min)
    return 1;

  if (ph->X[1] > X1max) {
    if (ph->w < wmin) {
      if (monty_rand() <= 1. / ROULETTE) {
        ph->w *= ROULETTE;
      } else
        ph->w = 0.;
    }
    return 1;
  }

  if (ph->w < wmin) {
    if (monty_rand() <= 1. / ROULETTE) {
      ph->w *= ROULETTE;
    } else {
      ph->w = 0.;
      return 1;
    }
  }

  return (0);
}

int record_criterion(struct of_photon *ph)
{
  const double X1max = log(RMAX);

  if (ph->X[1] > X1max)
    return (1);

  else
    return (0);

}
#undef RMAX
#undef ROULETTE

#define EPS 0.04
double stepsize(double X[NDIM], double K[NDIM])
{
  double dl, dlx1, dlx2, dlx3;
  double idlx1, idlx2, idlx3;

  dlx1 = EPS * X[1] / (fabs(K[1]) + SMALL);
  dlx2 = EPS * GSL_MIN(X[2], stopx[2] - X[2]) / (fabs(K[2]) + SMALL);
  dlx3 = EPS / (fabs(K[3]) + SMALL);

  idlx1 = 1. / (fabs(dlx1) + SMALL);
  idlx2 = 1. / (fabs(dlx2) + SMALL);
  idlx3 = 1. / (fabs(dlx3) + SMALL);

  dl = 1. / (idlx1 + idlx2 + idlx3);

  return (dl);
}
#undef EPS

void record_super_photon(struct of_photon *ph)
{
  double lE, dx2;
  int iE, ix2;

  if (isnan(ph->w) || isnan(ph->E)) {
    fprintf(stderr, "record isnan: %g %g\n", ph->w, ph->E);
    return;
  }

  #pragma omp critical (MAXTAU)
  {
    if (ph->tau_scatt > max_tau_scatt)
      max_tau_scatt = ph->tau_scatt;
  }

  // Bin in X2 coord. Get theta bin, while folding around equator
  dx2 = (stopx[2] - startx[2]) / (2. * N_THBINS);
  if (ph->X[2] < 0.5 * (startx[2] + stopx[2]))
    ix2 = (int) (ph->X[2] / dx2);
  else
    ix2 = (int) ((stopx[2] - ph->X[2]) / dx2);

  // Check limits
  if (ix2 < 0 || ix2 >= N_THBINS)
    return;

  // Get energy bin (centered on iE*dlE + lE0)
  lE = log(ph->E);
  iE = (int) ((lE - lE0) / dlE + 2.5) - 2;

  // Check limits
  if (iE < 0 || iE >= N_EBINS)
    return;

  #pragma omp atomic
  N_superph_recorded++;
  #pragma omp atomic
  N_scatt += ph->nscatt;

  // Add superphoton to spectrum
  spect[ix2][iE].dNdlE += ph->w;
  spect[ix2][iE].dEdlE += ph->w * ph->E;
  spect[ix2][iE].tau_abs += ph->w * ph->tau_abs;
  spect[ix2][iE].tau_scatt += ph->w * ph->tau_scatt;
  spect[ix2][iE].X1iav += ph->w * ph->X1i;
  spect[ix2][iE].X2isq += ph->w * (ph->X2i * ph->X2i);
  spect[ix2][iE].X3fsq += ph->w * (ph->X[3] * ph->X[3]);
  spect[ix2][iE].ne0 += ph->w * (ph->ne0);
  spect[ix2][iE].b0 += ph->w * (ph->b0);
  spect[ix2][iE].thetae0 += ph->w * (ph->thetae0);
  spect[ix2][iE].nscatt += ph->w * ph->nscatt;
  spect[ix2][iE].nph += 1.;
}

struct of_spectrum shared_spect[N_THBINS][N_EBINS] = { };

void omp_reduce_spect()
{
  #pragma omp parallel
  {
    #pragma omp critical
    {
      for (int i = 0; i < N_THBINS; i++) {
        for (int j = 0; j < N_EBINS; j++) {
          shared_spect[i][j].dNdlE +=
              spect[i][j].dNdlE;
          shared_spect[i][j].dEdlE +=
              spect[i][j].dEdlE;
          shared_spect[i][j].tau_abs +=
              spect[i][j].tau_abs;
          shared_spect[i][j].tau_scatt +=
              spect[i][j].tau_scatt;
          shared_spect[i][j].X1iav +=
              spect[i][j].X1iav;
          shared_spect[i][j].X2isq +=
              spect[i][j].X2isq;
          shared_spect[i][j].X3fsq +=
              spect[i][j].X3fsq;
          shared_spect[i][j].ne0 += spect[i][j].ne0;
          shared_spect[i][j].b0 += spect[i][j].b0;
          shared_spect[i][j].thetae0 +=
              spect[i][j].thetae0;
          shared_spect[i][j].nscatt +=
              spect[i][j].nscatt;
          shared_spect[i][j].nph += spect[i][j].nph;
        }
      }
    } // omp critical

    #pragma omp barrier

    #pragma omp master
    {
      for (int i = 0; i < N_THBINS; i++) {
        for (int j = 0; j < N_EBINS; j++) {
          spect[i][j].dNdlE =
              shared_spect[i][j].dNdlE;
          spect[i][j].dEdlE =
              shared_spect[i][j].dEdlE;
          spect[i][j].tau_abs =
              shared_spect[i][j].tau_abs;
          spect[i][j].tau_scatt =
              shared_spect[i][j].tau_scatt;
          spect[i][j].X1iav =
              shared_spect[i][j].X1iav;
          spect[i][j].X2isq =
              shared_spect[i][j].X2isq;
          spect[i][j].X3fsq =
              shared_spect[i][j].X3fsq;
          spect[i][j].ne0 = shared_spect[i][j].ne0;
          spect[i][j].b0 = shared_spect[i][j].b0;
          spect[i][j].thetae0 =
              shared_spect[i][j].thetae0;
          spect[i][j].nscatt =
              shared_spect[i][j].nscatt;
          spect[i][j].nph = shared_spect[i][j].nph;
        }
      }
    } // omp master
  } // omp parallel
}

double bias_func(double Te, double w)
{
  double bias, max ;

  max = 0.5 * w / WEIGHT_MIN;

  //bias = Te*Te;
  bias = 10.*Te*Te/(5. * max_tau_scatt);
  //bias = 100. * Te * Te / (bias_norm * max_tau_scatt);

  //if (bias < TP_OVER_TE)
  //  bias = TP_OVER_TE;
  if (bias > max)
    bias = max;

  return bias;// / TP_OVER_TE;
}

void get_fluid_zone(int i, int j, int k, double *Ne, double *Thetae, double *B,
        double Ucon[NDIM], double Bcon[NDIM])
{
  double Ucov[NDIM], Bcov[NDIM];
  double Bp[NDIM], Vcon[NDIM], Vfac, VdotV, UdotBp;
  double sig ;

  *Ne = p[KRHO][i][j][k] * Ne_unit;
  if (with_electrons) {
    *Thetae = p[KELCOND][i][j][k]*pow(p[KRHO][i][j][k],game-1.)*Thetae_unit;
  } else {
    *Thetae = p[UU][i][j][k] / (*Ne) * Ne_unit * Thetae_unit;
  }

  Bp[1] = p[B1][i][j][k];
  Bp[2] = p[B2][i][j][k];
  Bp[3] = p[B3][i][j][k];

  Vcon[1] = p[U1][i][j][k];
  Vcon[2] = p[U2][i][j][k];
  Vcon[3] = p[U3][i][j][k];

  // Get Ucov
  VdotV = 0.;
  for (int l = 1; l < NDIM; l++)
    for (int m = 1; m < NDIM; m++)
      VdotV += geom[i][j].gcov[l][m] * Vcon[l] * Vcon[m];
  Vfac = sqrt(-1. / geom[i][j].gcon[0][0] * (1. + fabs(VdotV)));
  Ucon[0] = -Vfac * geom[i][j].gcon[0][0];
  for (int l = 1; l < NDIM; l++)
    Ucon[l] = Vcon[l] - Vfac * geom[i][j].gcon[0][l];
  lower(Ucon, geom[i][j].gcov, Ucov);

  // Get B and Bcov
  UdotBp = 0.;
  for (int l = 1; l < NDIM; l++)
    UdotBp += Ucov[l] * Bp[l];
  Bcon[0] = UdotBp;
  for (int l = 1; l < NDIM; l++)
    Bcon[l] = (Bp[l] + Ucon[l] * UdotBp) / Ucon[0];
  lower(Bcon, geom[i][j].gcov, Bcov);

  *B = sqrt(Bcon[0] * Bcov[0] + Bcon[1] * Bcov[1] +
      Bcon[2] * Bcov[2] + Bcon[3] * Bcov[3]) * B_unit;

  if (*Thetae > THETAE_MAX) *Thetae = THETAE_MAX;

  sig = pow(*B/B_unit,2)/(*Ne/Ne_unit);
  if(sig > 1.) *Thetae = SMALL;//*Ne = 1.e-10*Ne_unit;
}

void get_fluid_params(double X[NDIM], double gcov[NDIM][NDIM], double *Ne,
          double *Thetae, double *B, double Ucon[NDIM],
          double Ucov[NDIM], double Bcon[NDIM],
          double Bcov[NDIM])
{
  double rho, kel;
  double Bp[NDIM], Vcon[NDIM], Vfac, VdotV, UdotBp;
  double gcon[NDIM][NDIM];
  double interp_scalar(double X[NDIM], double ***var);
  double sig ;

  if (X[1] < startx[1] || X[1] > stopx[1] ||
      X[2] < startx[2] || X[2] > stopx[2]) {

    *Ne = 0.;

    return;
  }

  rho = interp_scalar(X, p[KRHO]);
  kel = interp_scalar(X, p[KELCOND]);

  *Ne = rho*Ne_unit;
  if (with_electrons) {
    *Thetae = kel*pow(rho,game-1.)*Thetae_unit;
  } else {
    double uu = interp_scalar(X, p[UU]);
    *Thetae = uu/rho*Thetae_unit;
  }

  Bp[1] = interp_scalar(X, p[B1]);
  Bp[2] = interp_scalar(X, p[B2]);
  Bp[3] = interp_scalar(X, p[B3]);

  Vcon[1] = interp_scalar(X, p[U1]);
  Vcon[2] = interp_scalar(X, p[U2]);
  Vcon[3] = interp_scalar(X, p[U3]);

  gcov_func(X, gcov);
  gcon_func(gcov, gcon);

  // Get Ucov
  VdotV = 0.;
  for (int i = 1; i < NDIM; i++)
    for (int j = 1; j < NDIM; j++)
      VdotV += gcov[i][j] * Vcon[i] * Vcon[j];
  Vfac = sqrt(-1. / gcon[0][0] * (1. + fabs(VdotV)));
  Ucon[0] = -Vfac * gcon[0][0];
  for (int i = 1; i < NDIM; i++)
    Ucon[i] = Vcon[i] - Vfac * gcon[0][i];
  lower(Ucon, gcov, Ucov);

  // Get B and Bcov
  UdotBp = 0.;
  for (int i = 1; i < NDIM; i++)
    UdotBp += Ucov[i] * Bp[i];
  Bcon[0] = UdotBp;
  for (int i = 1; i < NDIM; i++)
    Bcon[i] = (Bp[i] + Ucon[i] * UdotBp) / Ucon[0];
  lower(Bcon, gcov, Bcov);

  *B = sqrt(Bcon[0] * Bcov[0] + Bcon[1] * Bcov[1] +
      Bcon[2] * Bcov[2] + Bcon[3] * Bcov[3]) * B_unit;

  if (*Thetae > THETAE_MAX) *Thetae = THETAE_MAX ;

  sig = pow(*B/B_unit,2)/(*Ne/Ne_unit);
  if (sig > 1.) *Thetae = SMALL;//*Ne = 1.e-10*Ne_unit;
}

////////////////////////////////// COORDINATES /////////////////////////////////
void gcov_func(double X[NDIM], double gcov[NDIM][NDIM])
{
  MUNULOOP gcov[mu][nu] = 0.;

  double sth, cth, s2, rho2;
  double r, th;

  bl_coord(X, &r, &th);

  cth = cos(th);
  sth = sin(th);

  s2 = sth*sth;
  rho2 = r*r + a*a*cth*cth;

  // KS -> MKS transformation
  double tfac, rfac, hfac, pfac;
  tfac = 1.;
  rfac = r - R0;
  hfac = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2]);
  pfac = 1.;

  gcov[0][0] = (-1. + 2.*r/rho2)*tfac*tfac;
  gcov[0][1] = (2.*r/rho2)*tfac*rfac;
  gcov[0][3] = (-2.*a*r*s2/rho2)*tfac*pfac;

  gcov[1][0] = gcov[0][1];
  gcov[1][1] = (1. + 2.*r/rho2)*rfac*rfac;
  gcov[1][3] = (-a*s2*(1. + 2.*r/rho2))*rfac*pfac;

  gcov[2][2] = (rho2)*hfac*hfac;

  gcov[3][0] = gcov[0][3];
  gcov[3][1] = gcov[1][3];
  gcov[3][3] = (s2*(rho2 + a*a*s2*(1. + 2.*r/rho2)))*pfac*pfac;
}

void bl_coord(double *X, double *r, double *th)
{
  *r = exp(X[1]) + R0;
  *th = M_PI*X[2] + ((1. - hslope)/2.)*sin(2.*M_PI*X[2]);
}

//double dOmega_func(double Xi[NDIM], double Xf[NDIM])
double dOmega_func(int j)
{
  double dbin = (stopx[2]-startx[2])/(2.*N_THBINS);
  double Xi[NDIM] = {0., stopx[1], j*dbin, 0.};
  double Xf[NDIM] = {0., stopx[1], (j+1)*dbin, 0.};

  double ri, rf, thi, thf;
  bl_coord(Xi, &ri, &thi);
  bl_coord(Xf, &rf, &thf);

  return 2.*M_PI*(-cos(thf) + cos(thi));
}

//////////////////////////////// INITIALIZATION ////////////////////////////////
void init_data(int argc, char *argv[])
{
  FILE *fp;
  double x[NDIM];
  double rp, hp, V, dV;//, two_temp_gam;

  // Header variables not used except locally
  double t, tf, cour, DTd, DTl, DTi, dt;
  int nstep, DTr, dump_cnt, image_cnt, rdump_cnt, lim, failed;
  double r, h, divb, vmin, vmax, gdet;
  double Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
  //double J;

  char *fname = argv[2];
  fp = fopen(fname, "r");
  if (fp == NULL) {
    fprintf(stderr, "Can't open sim data file\n");
    exit(-1);
  } else {
    fprintf(stderr, "Successfully opened %s\n", fname);
  }

  int idum;
  double fdum;

  // Get standard HARM header
  safe_fscanf(fp, "%lf ", &t);
  safe_fscanf(fp, "%d ", &N1);
  safe_fscanf(fp, "%d ", &N2);
  N3 = 1;
  safe_fscanf(fp, "%d ", &idum); // N3
  safe_fscanf(fp, "%lf ", &startx[1]);
  safe_fscanf(fp, "%lf ", &startx[2]);
  safe_fscanf(fp, "%lf ", &fdum); // startx[3]
  safe_fscanf(fp, "%lf ", &dx[1]);
  safe_fscanf(fp, "%lf ", &dx[2]);
  safe_fscanf(fp, "%lf ", &fdum); // dx[3]
  safe_fscanf(fp, "%lf ", &tf);
  safe_fscanf(fp, "%d ", &nstep);
  safe_fscanf(fp, "%lf ", &MBH);
  safe_fscanf(fp, "%lf ", &a);
  safe_fscanf(fp, "%lf ", &L_unit);
  safe_fscanf(fp, "%lf ", &T_unit);
  safe_fscanf(fp, "%lf ", &M_unit);
  safe_fscanf(fp, "%lf ", &Thetae_unit);
  safe_fscanf(fp, "%lf ", &gam);
  safe_fscanf(fp, "%lf ", &game);
  safe_fscanf(fp, "%lf ", &fdum); // gamp
  safe_fscanf(fp, "%d ", &idum); // RADMODEL
  safe_fscanf(fp, "%lf ", &fdum); // tp_over_te
  safe_fscanf(fp, "%lf ", &cour); 
  safe_fscanf(fp, "%lf ", &DTd);
  safe_fscanf(fp, "%lf ", &DTl);
  safe_fscanf(fp, "%lf ", &DTi);
  safe_fscanf(fp, "%d ", &DTr);
  safe_fscanf(fp, "%d ", &dump_cnt);
  safe_fscanf(fp, "%d ", &image_cnt);
  safe_fscanf(fp, "%d ", &rdump_cnt);
  safe_fscanf(fp, "%lf ", &dt);
  safe_fscanf(fp, "%d ", &lim);
  safe_fscanf(fp, "%d ", &failed);
  safe_fscanf(fp, "%lf ", &Rin);
  safe_fscanf(fp, "%lf ", &Rout);
  safe_fscanf(fp, "%lf ", &hslope);
  safe_fscanf(fp, "%lf ", &R0);
  safe_fscanf(fp, "%d ", &idum); // WITH_ELECTRONS
  safe_fscanf(fp, "%d ", &idum); // SPEC_THETABINS
  safe_fscanf(fp, "%d ", &idum); // SPEC_FREQBINS
  safe_fscanf(fp, "%lf ", &fdum); // SPEC_NUMIN
  safe_fscanf(fp, "%lf ", &fdum); // SPEC_NUMAX
  safe_fscanf(fp, "%lf ", &idum); // MONIKA_TPTE

  Rh = 1. + sqrt(1. - a*a);

  with_electrons = 1;

  printf("L_unit = %e T_unit = %e M_unit = %e\n", L_unit, T_unit, M_unit);

  // Nominal values for axisymmetric simulations
  startx[0] = 0.;
  startx[3] = 0.;

  stopx[0] = 1.;
  stopx[1] = startx[1] + N1*dx[1];
  stopx[2] = startx[2] + N2*dx[2];
  stopx[3] = 2.*M_PI;

  dx[0] = 1.;
  dx[3] = 2.*M_PI;

  // Finish setting up unit conversions
  RHO_unit = M_unit / pow(L_unit, 3);
  U_unit = RHO_unit * CL * CL;
  B_unit = CL * sqrt(4. * M_PI * RHO_unit);
  Ne_unit = RHO_unit / (MP + ME);
  max_tau_scatt = (6. * L_unit) * RHO_unit * 0.4;

  // Allocate storage and set geometry
  double ****malloc_rank4_double(int n1, int n2, int n3, int n4);
  //p = (double****)malloc_rank4(NVAR, N1, N2, N3, sizeof(double));
  p = malloc_rank4_double(NVAR, N1, N2, N3);
  printf("NVAR N1 N2 N3 = %i %i %i %i\n", NVAR, N1, N2, N3);
  geom = (struct of_geom**)malloc_rank2(N1, N2, sizeof(struct of_geom));
  init_geometry();

  dMact = 0.;
  Ladv = 0.;
  bias_norm = 0.;
  V = 0.;
  dV = dx[1]*dx[2]*dx[3];
  for (int np = 0; np < N1*N2; np++) {
    int j = np % N2;
    int i = (np - j) / N2;
    int k = 0;
    safe_fscanf(fp, "%lf %lf %lf %lf", &x[1], &x[2], &r, &h);

    // Check that coordinates are right
    bl_coord(x, &rp, &hp);
    if (fabs(rp - r) > 1.e-5 * rp || fabs(hp - h) > 1.e-5) {
      printf("x[] = %e %e %e %e\n", x[0], x[1], x[2], x[3]);
      fprintf(stderr, "grid setup error [%i %i %i]\n", i, j, k);
      fprintf(stderr, "rp,r,hp,h: %g %g %g %g\n", rp, r, hp, h);
      printf("hslope = %e R0 = %e\n", hslope, R0);
      fprintf(stderr, "edit R0, hslope, compile, and continue\n");
      exit(-1);
    }

    safe_fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
      &p[KRHO][i][j][k],
      &p[UU][i][j][k],
      &p[U1][i][j][k],
      &p[U2][i][j][k],
      &p[U3][i][j][k],
      &p[B1][i][j][k],
      &p[B2][i][j][k],
      &p[B3][i][j][k],
      &p[KTOT][i][j][k],
      &p[KELCOND][i][j][k],
      &p[KELNOCOND][i][j][k],
      &p[PHI][i][j][k],
      &p[FLR][i][j][k]);

    safe_fscanf(fp, "%lf", &divb);

    safe_fscanf(fp, "%lf %lf %lf %lf",
      &Ucon[0], &Ucon[1], &Ucon[2], &Ucon[3]);
    safe_fscanf(fp, "%lf %lf %lf %lf", &Ucov[0],
      &Ucov[1], &Ucov[2], &Ucov[3]);
    safe_fscanf(fp, "%lf %lf %lf %lf", &Bcon[0],
      &Bcon[1], &Bcon[2], &Bcon[3]);
    safe_fscanf(fp, "%lf %lf %lf %lf", &Bcov[0],
      &Bcov[1], &Bcov[2], &Bcov[3]);

    safe_fscanf(fp, "%lf ", &vmin);
    safe_fscanf(fp, "%lf ", &vmax);
    safe_fscanf(fp, "%lf ", &vmin);
    safe_fscanf(fp, "%lf ", &vmax);
    safe_fscanf(fp, "%lf ", &gdet);

    safe_fscanf(fp, "%lf ", &fdum); // <G_0>
    safe_fscanf(fp, "%lf ", &fdum); // <G_1>
    safe_fscanf(fp, "%lf ", &fdum); // <G_2>
    safe_fscanf(fp, "%lf ", &fdum); // <G_3>

    safe_fscanf(fp, "%lf ", &fdum); // qud
    safe_fscanf(fp, "%lf ", &fdum); // qvisc
    safe_fscanf(fp, "%lf ", &fdum); // qcoul

    safe_fscanf(fp, "%lf ", &fdum); // N_esuper
    safe_fscanf(fp, "%lf ", &fdum); // N_esuper_electron

    safe_fscanf(fp, "%lf ", &fdum); // Thetae

    bias_norm += dV*gdet*pow(p[KELCOND][i][j][k]*pow(p[KRHO][i][j][k],game-1.)*Thetae_unit, 2.);
    V += dV*gdet;

    if (i <= 20) dMact += gdet*p[KRHO][i][j][k]*Ucon[1];
    if (i >= 20 && i < 40) Ladv += gdet*p[UU][i][j][k]*Ucon[1]*Ucov[0];
  }

  bias_norm /= V;
  dMact *= dx[3]*dx[2];
  dMact /= 21.;
  Ladv *= dx[3]*dx[2];
  Ladv /= 21.;
  fprintf(stderr, "dMact: %g, Ladv: %g\n", dMact, Ladv);
}

//////////////////////////////////// OUTPUT ////////////////////////////////////

#define SPECTRUM_FILE_NAME "spectrum.dat"
void report_spectrum(int N_superph_made)
{
  double dOmega, nuLnu, tau_scatt, L;//, Xi[NDIM], Xf[NDIM];
  FILE *fp;

  double nu0,nu1,nu,fnu ;
  double dsource = 8000*PC ;

  fp = fopen(SPECTRUM_FILE_NAME, "w");
  if (fp == NULL) {
    fprintf(stderr, "trouble opening spectrum file\n");
    exit(0);
  }

  /* output */
  max_tau_scatt = 0.;
  L = 0.;
  double dL = 0.;
  for (int i = 0; i < N_EBINS; i++) {
    // Output log_10(photon energy/(me c^2))
    fprintf(fp, "%10.5g ", (i * dlE + lE0) / M_LN10);

    for (int j = 0; j < N_THBINS; j++) {
      // Convert accumulated photon number to nuLnu, in units of Lsun
      //coord(N1-1, j, 0, Xi);
      //coord(N2-1, j+1, 0, Xf);
      dOmega = 2.*dOmega_func(j);
      //dOmega = 2.*dOmega_func(Xi, Xf);

      nuLnu = (ME*CL*CL)*(4.*M_PI/dOmega)*(1./dlE);

      nuLnu *= spect[j][i].dEdlE/LSUN;
      dL += ME*CL*CL*spect[j][i].dEdlE;

      tau_scatt = spect[j][i].tau_scatt/(spect[j][i].dNdlE + SMALL);

      fprintf(fp, "%10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g ",
        nuLnu,
        dOmega,
        spect[j][i].tau_abs/(spect[j][i].dNdlE + SMALL),
        tau_scatt,
        spect[j][i].X1iav/(spect[j][i].dNdlE + SMALL),
        sqrt(fabs(spect[j][i].X2isq/(spect[j][i].dNdlE + SMALL))),
        sqrt(fabs(spect[j][i].X3fsq/(spect[j][i].dNdlE + SMALL))));


      nu0 = ME * CL * CL * exp((i - 0.5) * dlE + lE0) / HPL ;
      nu1 = ME * CL * CL * exp((i + 0.5) * dlE + lE0) / HPL ;

      if(nu0 < 230.e9 && nu1 > 230.e9) {
        nu = ME * CL * CL * exp(i * dlE + lE0) / HPL ;
        fnu = nuLnu*LSUN/(4.*M_PI*dsource*dsource*nu*JY) ;
        fprintf(stderr,"fnu: %10.5g\n",fnu) ;
      }

      // Average # scatterings
      fprintf(fp,"%10.5g ",spect[j][i].nscatt/(spect[j][i].dNdlE + SMALL));

      if (tau_scatt > max_tau_scatt)
        max_tau_scatt = tau_scatt;

      L += nuLnu * dOmega * dlE / (4. * M_PI);
    }
    fprintf(fp, "\n");
  }
  printf("dL = %e\n", dL);
  fprintf(stderr,
    "luminosity %g, dMact %g, efficiency %g, L/Ladv %g, max_tau_scatt %g\n",
    L, dMact * M_unit / T_unit / (MSUN / YEAR),
    L * LSUN / (dMact * M_unit * CL * CL / T_unit),
    L * LSUN / (Ladv * M_unit * CL * CL / T_unit),
    max_tau_scatt);

  double LEdd = 4.*M_PI*GNEWT*MBH*MP*CL/(SIGMA_THOMSON);
  double MdotEdd = 4.*M_PI*GNEWT*MBH*MP/(SIGMA_THOMSON*CL*0.1);
  printf("MdotEdd = %e\n", MdotEdd);
  double Mdot = dMact*M_unit/T_unit;
  double mdot = Mdot/MdotEdd;
  printf("Mdot = %e mdot = %e\n", Mdot, mdot);
  double Lum = L*LSUN;
  double lum = Lum/LEdd;
  printf("L = %e lum = %e\n", Lum, lum);

  fprintf(stderr, "\n");
  fprintf(stderr, "N_superph_made: %d\n", N_superph_made);
  fprintf(stderr, "N_superph_scatt: %d\n", N_scatt);
  fprintf(stderr, "N_superph_recorded: %d\n", N_superph_recorded);

  fclose(fp);
}

