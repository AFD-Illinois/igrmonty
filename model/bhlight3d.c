#include "decs.h"

#define NVAR (10)

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

static double poly_norm, poly_xt, poly_alpha, mks_smooth, game, gamp;
static double MBH;
  
static int with_radiation;
static int with_derefine_poles;
static int with_electrons;

void report_bad_input(int argc) 
{
  if (argc < 6) {
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
  bias = Te*Te/(5. * max_tau_scatt);
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
    *Thetae = p[KEL][i][j][k]*pow(p[KRHO][i][j][k],game-1.)*Thetae_unit;
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
  kel = interp_scalar(X, p[KEL]);

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

void set_dxdX(double X[NDIM], double dxdX[NDIM][NDIM])
{
  MUNULOOP dxdX[mu][nu] = 0.;

  dxdX[0][0] = 1.;
  dxdX[1][1] = exp(X[1]);
  if (with_derefine_poles) {
  dxdX[2][1] = -exp(mks_smooth*(startx[1]-X[1]))*mks_smooth*(
    M_PI/2. -
    M_PI*X[2] +
    poly_norm*(2.*X[2]-1.)*(1+(pow((-1.+2*X[2])/poly_xt,poly_alpha))/(1 + poly_alpha)) -
    1./2.*(1. - hslope)*sin(2.*M_PI*X[2])
    );
  dxdX[2][2] = M_PI + (1. - hslope)*M_PI*cos(2.*M_PI*X[2]) +
    exp(mks_smooth*(startx[1]-X[1]))*(
      -M_PI +
      2.*poly_norm*(1. + pow((2.*X[2]-1.)/poly_xt,poly_alpha)/(poly_alpha+1.)) +
      (2.*poly_alpha*poly_norm*(2.*X[2]-1.)*pow((2.*X[2]-1.)/poly_xt,poly_alpha-1.))/((1.+poly_alpha)*poly_xt) -
      (1.-hslope)*M_PI*cos(2.*M_PI*X[2])
      );
  } else {
    dxdX[2][2] = M_PI - (hslope - 1.)*M_PI*cos(2.*M_PI*X[2]);
  }
  dxdX[3][3] = 1.;
}

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

  gcov[0][0] = -1. + 2.*r/rho2;
  gcov[0][1] = 2.*r/rho2;
  gcov[0][3] = -2.*a*r*s2/rho2;

  gcov[1][0] = gcov[0][1];
  gcov[1][1] = 1. + 2.*r/rho2;
  gcov[1][3] = -a*s2*(1. + 2.*r/rho2);

  gcov[2][2] = rho2;

  gcov[3][0] = gcov[0][3];
  gcov[3][1] = gcov[1][3];
  gcov[3][3] = s2*(rho2 + a*a*s2*(1. + 2.*r/rho2));

  // Apply coordinate transformation to code coordinates X
  double dxdX[NDIM][NDIM];
  set_dxdX(X, dxdX);

  double gcov_ks[NDIM][NDIM];
  MUNULOOP {
    gcov_ks[mu][nu] = gcov[mu][nu];
    gcov[mu][nu] = 0.;
  }

  MUNULOOP {
		for (int lam = 0; lam < NDIM; lam++) {
			for (int kap = 0; kap < NDIM; kap++) {
				gcov[mu][nu] += gcov_ks[lam][kap]*dxdX[lam][mu]*dxdX[kap][nu];
			}
		}
  }
}

void bl_coord(double *X, double *r, double *th)
{
  *r = exp(X[1]) + R0;
  double thG = M_PI*X[2] + ((1. - hslope)/2.)*sin(2.*M_PI*X[2]);
  double y = 2*X[2] - 1.;
  double thJ = poly_norm*y*(1. + pow(y/poly_xt,poly_alpha)/(poly_alpha+1.)) + 0.5*M_PI;
  *th = thG + exp(mks_smooth*(startx[1] - X[1]))*(thJ - thG);
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

#include <hdf5.h>
#include <hdf5_hl.h>

void init_data(int argc, char *argv[], Params *params)
{
  const char *fname = NULL;
  double dV, V;
  hid_t file_id;

  if (params->loaded && strlen(params->dump) > 0) {
    fname = params->dump;
  } else {
    fname = argv[2];
    strncpy((char *)params->dump, argv[2], 255);
  }

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) {
    fprintf(stderr, "File %s does not exist! Exiting...\n", fname);
    exit(-1);
  }

  // Read header
  H5LTread_dataset_int(file_id, "RADIATION", &with_radiation);
  H5LTread_dataset_int(file_id, "DEREFINE_POLES", &with_derefine_poles);
  H5LTread_dataset_int(file_id, "ELECTRONS", &with_electrons);
  H5LTread_dataset_int(file_id, "N1", &N1);
  H5LTread_dataset_int(file_id, "N2", &N2);
  H5LTread_dataset_int(file_id, "N3", &N3);
  H5LTread_dataset_double(file_id, "startx1",  &startx[1]);
  H5LTread_dataset_double(file_id, "startx2",  &startx[2]);
  H5LTread_dataset_double(file_id, "startx3",  &startx[3]);
  H5LTread_dataset_double(file_id, "dx1",  &dx[1]);
  H5LTread_dataset_double(file_id, "dx2",  &dx[2]);
  H5LTread_dataset_double(file_id, "dx3",  &dx[3]);
  H5LTread_dataset_double(file_id, "a", &a);
  H5LTread_dataset_double(file_id, "gam", &gam);
  if (with_electrons) {
    H5LTread_dataset_double(file_id, "game", &game);
    H5LTread_dataset_double(file_id, "gamp", &gamp);
  }
  H5LTread_dataset_double(file_id, "Rin", &Rin);
  H5LTread_dataset_double(file_id, "Rout", &Rout);
  H5LTread_dataset_double(file_id, "hslope", &hslope);
  H5LTread_dataset_double(file_id, "poly_xt", &poly_xt);
  H5LTread_dataset_double(file_id, "poly_alpha", &poly_alpha);
  H5LTread_dataset_double(file_id, "mks_smooth", &mks_smooth);

  printf("HDR!\n");

  // Set polylog grid normalization
  poly_norm = 0.5*M_PI*1./(1. + 1./(poly_alpha + 1.)*1./pow(poly_xt, poly_alpha));

  stopx[0] = 1.;
  stopx[1] = startx[1]+N1*dx[1];
  stopx[2] = startx[2]+N2*dx[2];
  stopx[3] = startx[3]+N3*dx[3];

  // SET UNITS IF NECESSARY
  // CHECK FOR NUMBER OF COMMAND LINE ARGS HERE!!!
  if (with_radiation) {
    H5LTread_dataset_double(file_id, "M_unit", &M_unit);
    H5LTread_dataset_double(file_id, "T_unit", &T_unit);
    H5LTread_dataset_double(file_id, "L_unit", &L_unit);
    H5LTread_dataset_double(file_id, "Thetae_unit", &Thetae_unit);
    H5LTread_dataset_double(file_id, "Mbh", &MBH);
    H5LTread_dataset_double(file_id, "tp_over_te", &TP_OVER_TE);
  } else {

    // Enough command line args?
    if (! params->loaded) { 
      report_bad_input(argc);
      sscanf(argv[3], "%lf", &M_unit);
      sscanf(argv[4], "%lf", &MBH);
      sscanf(argv[5], "%lf", &TP_OVER_TE);
    } else {
      M_unit = params->M_unit;
      MBH = params->MBH;
      TP_OVER_TE = params->TP_OVER_TE;
    }

    MBH *= MSUN;

    L_unit = GNEWT*MBH/(CL*CL);
    T_unit = L_unit/CL;

    if (with_electrons) {
      Thetae_unit = MP/ME;
    } else {
      Thetae_unit = MP/ME*(gam-1.)*1./(1. + TP_OVER_TE);
    }
  }

  // Set remaining units and constants
  RHO_unit = M_unit/pow(L_unit,3);
  U_unit = RHO_unit*CL*CL;
  B_unit = CL*sqrt(4.*M_PI*RHO_unit);
  Ne_unit = RHO_unit/(MP + ME);
  max_tau_scatt = (6.*L_unit)*RHO_unit*0.4;
  Rh = 1. + sqrt(1. - a * a);

  printf("M_unit = %e\n", M_unit);
  printf("Ne_unit = %e\n", Ne_unit);
  printf("RHO_unit = %e\n", RHO_unit);
  printf("L_unit = %e\n", L_unit);
  printf("T_unit = %e\n", T_unit);
  printf("B_unit = %e\n", B_unit);
  printf("Thetae_unit = %e\n", Thetae_unit);

  // Allocate storage and set geometry
  double ****malloc_rank4_double(int n1, int n2, int n3, int n4);
  //p = (double****)malloc_rank4(NVAR, N1, N2, N3, sizeof(double));
  p = malloc_rank4_double(NVAR, N1, N2, N3);
  printf("N1 N2 N3 = %i %i %i %i\n", NVAR, N1, N2, N3);
  geom = (struct of_geom**)malloc_rank2(N1, N2, sizeof(struct of_geom));
  init_geometry();

  // Read data
  H5LTread_dataset_double(file_id, "RHO",  &p[KRHO][0][0][0]);
  H5LTread_dataset_double(file_id, "UU",   &p[UU][0][0][0]);
  H5LTread_dataset_double(file_id, "U1",   &p[U1][0][0][0]);
  H5LTread_dataset_double(file_id, "U2",   &p[U2][0][0][0]);
  H5LTread_dataset_double(file_id, "U3",   &p[U3][0][0][0]);
  H5LTread_dataset_double(file_id, "B1",   &p[B1][0][0][0]);
  H5LTread_dataset_double(file_id, "B2",   &p[B2][0][0][0]);
  H5LTread_dataset_double(file_id, "B3",   &p[B3][0][0][0]);
  if (with_electrons) {
    H5LTread_dataset_double(file_id, "KEL",  &p[KEL][0][0][0]);
    H5LTread_dataset_double(file_id, "KTOT", &p[KTOT][0][0][0]);
  }

  H5Fclose(file_id);

  V = dMact = Ladv = 0.;
  dV = dx[1]*dx[2]*dx[3];
  ZLOOP {
    V += dV*geom[i][j].g;
    bias_norm += dV*geom[i][j].g*pow(p[UU][i][j][k]/p[KRHO][i][j][k]*Thetae_unit,2.);

    if (i <= 20) {
      double Ne, Thetae, Bmag, Ucon[NDIM], Ucov[NDIM], Bcon[NDIM];
      get_fluid_zone(i, j, k, &Ne, &Thetae, &Bmag, Ucon, Bcon);
      lower(Ucon, geom[i][j].gcov, Ucov);
      dMact += geom[i][j].g*dx[2]*dx[3]*p[KRHO][i][j][k]*Ucon[1];
      Ladv += geom[i][j].g*dx[2]*dx[3]*p[UU][i][j][k]*Ucon[1]*Ucov[0];
    }

    //double X[NDIM], r, th;
    //coord(i,j,k,X);
    //bl_coord(X, &r, &th);
    //if (j == N2/2) {
    //  printf("[%i] r = %e RHO = %e\n", i, r, p[KRHO][i][j][k]);
   // }
  }

  //exit(-1);

  dMact /= 21.;
  Ladv /= 21.;
  bias_norm /= V;
  fprintf(stderr, "dMact: %g, Ladv: %g\n", dMact, Ladv);
}

//////////////////////////////////// OUTPUT ////////////////////////////////////

#if HDF5_OUTPUT
void report_spectrum(int N_superph_made, Params *params)
{

  hid_t fid = -1;

  if (params->loaded && strlen(params->spectrum) > 0) {
    fid = H5Fcreate(params->spectrum, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  } else {
    fid = H5Fcreate("spectrum.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  }

  if (fid < 0) {
    fprintf(stderr, "! unable to open/create hdf5 file.\n");
    exit(-3);
  }

  h5io_add_attribute_str(fid, "/", "githash", xstr(VERSION));

  h5io_add_group(fid, "/params");

  h5io_add_data_dbl(fid, "/params/NUCUT", NUCUT);
  h5io_add_data_dbl(fid, "/params/GAMMACUT", GAMMACUT);
  h5io_add_data_dbl(fid, "/params/NUMAX", NUMAX);
  h5io_add_data_dbl(fid, "/params/NUMIN", NUMIN);
  h5io_add_data_dbl(fid, "/params/LNUMAX", LNUMAX);
  h5io_add_data_dbl(fid, "/params/LNUMIN", LNUMIN);
  h5io_add_data_dbl(fid, "/params/DLNU", DLNU);
  h5io_add_data_dbl(fid, "/params/THETAE_MAX", THETAE_MIN);
  h5io_add_data_dbl(fid, "/params/THETAE_MIN", THETAE_MAX);
  h5io_add_data_dbl(fid, "/params/TP_OVER_TE", TP_OVER_TE);
  h5io_add_data_dbl(fid, "/params/WEIGHT_MIN", WEIGHT_MIN);
  h5io_add_data_dbl(fid, "/params/KAPPA", KAPPA);
  h5io_add_data_dbl(fid, "/params/L_unit", L_unit);
  h5io_add_data_dbl(fid, "/params/T_unit", T_unit);
  h5io_add_data_dbl(fid, "/params/M_unit", M_unit);
  h5io_add_data_dbl(fid, "/params/B_unit", B_unit);
  h5io_add_data_dbl(fid, "/params/Ne_unit", Ne_unit);
  h5io_add_data_dbl(fid, "/params/RHO_unit", RHO_unit);
  h5io_add_data_dbl(fid, "/params/U_unit", U_unit);
  h5io_add_data_dbl(fid, "/params/Thetae_unit", Thetae_unit);
  h5io_add_data_dbl(fid, "/params/MBH", MBH);
  h5io_add_data_dbl(fid, "/params/a", a);
  h5io_add_data_dbl(fid, "/params/Rin", Rin);
  h5io_add_data_dbl(fid, "/params/Rout", Rout);
  h5io_add_data_dbl(fid, "/params/hslope", hslope);
  h5io_add_data_dbl(fid, "/params/t", t);

  h5io_add_data_int(fid, "/params/SYNCHROTRON", SYNCHROTRON);
  h5io_add_data_int(fid, "/params/BREMSSTRAHLUNG", BREMSSTRAHLUNG);
  h5io_add_data_int(fid, "/params/COMPTON", COMPTON);
  h5io_add_data_int(fid, "/params/DIST_KAPPA", DIST_KAPPA);
  h5io_add_data_int(fid, "/params/N_ESAMP", N_ESAMP);
  h5io_add_data_int(fid, "/params/N_EBINS", N_EBINS);
  h5io_add_data_int(fid, "/params/N_THBINS", N_THBINS);
  h5io_add_data_int(fid, "/params/N1", N1);
  h5io_add_data_int(fid, "/params/N2", N2);
  h5io_add_data_int(fid, "/params/N3", N3);
  h5io_add_data_int(fid, "/params/Ns", Ns);

  h5io_add_data_str(fid, "/params/dump", params->dump);

  // temporary data buffers
  double lnu_buf[N_EBINS];
  double dOmega_buf[N_THBINS];
  double nuLnu_buf[N_EBINS][N_THBINS];
  double tau_abs_buf[N_EBINS][N_THBINS];
  double tau_scatt_buf[N_EBINS][N_THBINS];
  double x1av_buf[N_EBINS][N_THBINS];
  double x2av_buf[N_EBINS][N_THBINS];
  double x3av_buf[N_EBINS][N_THBINS];
  double nscatt_buf[N_EBINS][N_THBINS];

  // normal output routine
  double dOmega, nuLnu, tau_scatt, L, nu0, nu1, nu, fnu, dL;
  double dsource = 8000 * PC;

  max_tau_scatt = 0.;
  L = 0.;
  dL = 0.;

  for (int j=0; j<N_THBINS; ++j) {
    dOmega_buf[j] = 2. * dOmega_func(j);
  }

  for (int i=0; i<N_EBINS; ++i) {
    lnu_buf[i] = (i * dlE + lE0) / M_LN10;
    for (int j=0; j<N_THBINS; ++j) {

      dOmega = 2. * dOmega_func(j);

      nuLnu = (ME * CL * CL) * (4. * M_PI / dOmega) * (1. / dlE);
      nuLnu *= spect[j][i].dEdlE/LSUN;

      tau_scatt = spect[j][i].tau_scatt/(spect[j][i].dNdlE + SMALL);

      nuLnu_buf[i][j] = nuLnu;
      tau_abs_buf[i][j] = spect[j][i].tau_abs/(spect[j][i].dNdlE + SMALL);
      tau_scatt_buf[i][j] = tau_scatt;
      x1av_buf[i][j] = spect[j][i].X1iav/(spect[j][i].dNdlE + SMALL);
      x2av_buf[i][j] = sqrt(fabs(spect[j][i].X2isq/(spect[j][i].dNdlE + SMALL)));
      x3av_buf[i][j] = sqrt(fabs(spect[j][i].X3fsq/(spect[j][i].dNdlE + SMALL)));
      nscatt_buf[i][j] = spect[j][i].nscatt / (spect[j][i].dNdlE + SMALL);

      if (tau_scatt > max_tau_scatt) max_tau_scatt = tau_scatt;

      dL += ME * CL * CL * spect[j][i].dEdlE;
      L += nuLnu * dOmega * dlE / (4. * M_PI);

      nu0 = ME * CL * CL * exp((i-0.5) * dlE + lE0) / HPL;
      nu1 = ME * CL * CL * exp((i+0.5) * dlE + lE0) / HPL;

      if (nu0 < 230.e9 && nu1 > 230.e9) {
        nu = ME * CL * CL * exp(i * dlE + lE0) / HPL;
        fnu = nuLnu * LSUN / (4. * M_PI * dsource * dsource * nu * JY);
        fprintf(stderr, "fnu: %10.5g\n", fnu);
      }

    }
  }

  h5io_add_group(fid, "/output");

  h5io_add_data_dbl_1d(fid, "/output/lnu", N_EBINS, lnu_buf);
  h5io_add_data_dbl_1d(fid, "/output/dOmega", N_THBINS, dOmega_buf);
  h5io_add_data_dbl_2d(fid, "/output/nuLnu", N_EBINS, N_THBINS, nuLnu_buf);
  h5io_add_data_dbl_2d(fid, "/output/tau_abs", N_EBINS, N_THBINS, tau_abs_buf);
  h5io_add_data_dbl_2d(fid, "/output/tau_scatt", N_EBINS, N_THBINS, tau_scatt_buf);
  h5io_add_data_dbl_2d(fid, "/output/x1av", N_EBINS, N_THBINS, x1av_buf);
  h5io_add_data_dbl_2d(fid, "/output/x2av", N_EBINS, N_THBINS, x2av_buf);
  h5io_add_data_dbl_2d(fid, "/output/x3av", N_EBINS, N_THBINS, x3av_buf);
  h5io_add_data_dbl_2d(fid, "/output/nscatt", N_EBINS, N_THBINS, nscatt_buf);

  h5io_add_data_int(fid, "/output/Nrecorded", N_superph_recorded);
  h5io_add_data_int(fid, "/output/Nmade", N_superph_made);
  h5io_add_data_int(fid, "/output/Nscattered", N_scatt);

  double LEdd = 4. * M_PI * GNEWT * MBH * MP * CL / SIGMA_THOMSON;
  double MdotEdd = 4. * M_PI * GNEWT * MBH * MP / ( SIGMA_THOMSON * CL * 0.1 );
  double Lum = L * LSUN;
  double lum = Lum / LEdd;
  double Mdot = dMact * M_unit / T_unit;
  double mdot = Mdot / MdotEdd;

  h5io_add_data_dbl(fid, "/output/L", Lum);
  h5io_add_data_dbl(fid, "/output/Mdot", Mdot);
  h5io_add_data_dbl(fid, "/output/LEdd", LEdd);
  h5io_add_data_dbl(fid, "/output/MdotEdd", MdotEdd);

  h5io_add_attribute_str(fid, "/output/L", "units", "erg/s");
  h5io_add_attribute_str(fid, "/output/LEdd", "units", "erg/s");
  h5io_add_attribute_str(fid, "/output/Mdot", "units", "g/s");
  h5io_add_attribute_str(fid, "/output/MdotEdd", "units", "g/s");

  // diagnostic output to screen
  fprintf(stderr, "\n");

  fprintf(stderr, "MBH = %g Msun\n", MBH/MSUN);
  fprintf(stderr, "a = %g\n", a);

  fprintf(stderr, "dL = %g\n", dL);
  fprintf(stderr, "dMact = %g\n", dMact * M_unit / T_unit / (MSUN / YEAR));
  fprintf(stderr, "efficiency = %g\n", L * LSUN / (dMact * M_unit * CL * CL / T_unit));
  fprintf(stderr, "L/Ladv = %g\n", L * LSUN / (Ladv * M_unit * CL * CL / T_unit));
  fprintf(stderr, "max_tau_scatt = %g\n", max_tau_scatt);
  fprintf(stderr, "Mdot = %g g/s, MdotEdd = %g g/s, mdot = %g MdotEdd\n", Mdot, MdotEdd, mdot);
  fprintf(stderr, "L = %g erg/s, LEdd = %g erg/s, lum = %g LEdd\n", Lum, LEdd, lum);

  fprintf(stderr, "\n");

  fprintf(stderr, "N_superph_made = %d\n", N_superph_made);
  fprintf(stderr, "N_superph_scatt = %d\n", N_scatt);
  fprintf(stderr, "N_superph_recorded = %d\n", N_superph_recorded);

  H5Fclose(fid);

}

#else

#define SPECTRUM_FILE_NAME "spectrum.dat"
void report_spectrum(int N_superph_made, Params *params)
{
  double dOmega, nuLnu, tau_scatt, L;//, Xi[NDIM], Xf[NDIM];
  FILE *fp;

  double nu0,nu1,nu,fnu ;
  double dsource = 8000*PC ;

  if (params->loaded && strlen(params->spectrum) > 0) {
    fp = fopen(params->spectrum, "w");
  } else {
    fp = fopen(SPECTRUM_FILE_NAME, "w");
  }
  
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

#endif // HDF5_OUTPUT

