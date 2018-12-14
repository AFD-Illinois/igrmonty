#include "decs.h"

#define NVAR (10)
#define USE_FIXED_TPTE (0)
#define USE_MIXED_TPTE (1)

// electron model. these values will be overwritten by anything found in par.c 
// or in the runtime parameter file.
// with_electrons ->
//     0 : constant TP_OVER_TE
//     1 : use dump file model (kawazura?)
//     2: use mixed TP_OVER_TE (moscibrodzka "beta" model)
static double tp_over_te = 3.;
static double trat_small = 2.;
static double trat_large = 70.;
static double beta_crit = 1.;
static int with_electrons;

double biasTuning = 1.;

// fluid data
double ****bcon;
double ****bcov;
double ****ucon;
double ****ucov;
double ****p;
double ***ne;
double ***thetae;
double ***b;

double TP_OVER_TE;

// coordinate functions
static inline __attribute__((always_inline)) void set_dxdX_metric(double X[NDIM], double dxdX[NDIM][NDIM], int metric);
static inline __attribute__((always_inline)) void gcov_ks(double r, double th, double gcov[NDIM][NDIM]);

// metric parameters
//  note: if METRIC_eKS, then the code will use "exponentialKS" coordinates
//        defined by x^a = { x^0, log(x^1), x^2, x^3 } where x^0, x^1, x^2,
//        x^3 are normal KS coordinates. in addition you must set METRIC_*
//        in order to specify how Xtoijk and gdet_zone should work.
int METRIC_eKS;
static int with_derefine_poles, METRIC_MKS3;
static double poly_norm, poly_xt, poly_alpha, mks_smooth; // mmks
static double mks3R0, mks3H0, mks3MY1, mks3MY2, mks3MP0; // mks3

static double MBH, game, gamp;

static hdf5_blob fluid_header = { 0 };
  
static int with_radiation;

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

#define ROULETTE  1.e4
int stop_criterion(struct of_photon *ph)
{
  // stop if weight below minimum weight
  if (ph->w < WEIGHT_MIN) {
    if (monty_rand() <= 1. / ROULETTE) {
      ph->w *= ROULETTE;
    } else {
      ph->w = 0.;
      return 1;
    }
  }

  // right now, only support X->KS exponential coordinates
  double X1min = log(Rh*1.05);
  double X1max = log(Rmax*1.1);

  if (ph->X[1] < X1min || ph->X[1] > X1max) {
    return 1;
  }

  return 0;
}

int record_criterion(struct of_photon *ph)
{
  double X1max = log(Rmax*1.1);

  if (ph->X[1] > X1max) {
    return 1;
  }

  return 0;
}
#undef ROULETTE

#define EPS 0.04
double stepsize(double X[NDIM], double K[NDIM])
{
  double dl, dlx1, dlx2, dlx3;
  double idlx1, idlx2, idlx3;

  /*
  dlx1 = EPS * X[1] / (fabs(K[1]) + SMALL);
  dlx2 = EPS * GSL_MIN(X[2], stopx[2] - X[2]) / (fabs(K[2]) + SMALL);
  dlx3 = EPS / (fabs(K[3]) + SMALL);
   */
 
  #define MIN(A,B) (A<B?A:B)

  dlx1 = EPS / (fabs(K[1]) + SMALL);
  dlx2 = EPS * MIN(X[2], 1. - X[2]) / (fabs(K[2]) + SMALL);
  dlx3 = EPS / (fabs(K[3]) + SMALL) ;

  #undef MIN

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
  int iE, ix2, ic;

  if (isnan(ph->w) || isnan(ph->E)) {
    fprintf(stderr, "record isnan: %g %g\n", ph->w, ph->E);
    return;
  }

  #pragma omp critical (MAXTAU)
  {
    if (ph->tau_scatt > max_tau_scatt)
      max_tau_scatt = ph->tau_scatt;
  }

  // bin in X[2] BL coord while folding around the equator
  double r, th;
  bl_coord(ph->X, &r, &th);
  dx2 = M_PI/2./N_THBINS;
  if (th > M_PI/2.) {
    ix2 = (int)( (M_PI - th) / dx2 );
  } else {
    ix2 = (int)( th / dx2 );
  }

  // Check limits
  if (ix2 < 0 || ix2 >= N_THBINS)
    return;

  // Get energy bin (centered on iE*dlE + lE0)
  lE = log(ph->E);
  iE = (int) ((lE - lE0) / dlE + 2.5) - 2;

  // Check limits
  if (iE < 0 || iE >= N_EBINS)
    return;

  // Get compton bin
  ic = ph->nscatt;
  if (ic > 3) ic = 3;

  #pragma omp atomic
  N_superph_recorded++;
  #pragma omp atomic
  N_scatt += ph->nscatt;

  double ratio_synch = 1. - ph->ratio_brems;
  
  // Add superphoton to synch spectrum
  spect[ic][ix2][iE].dNdlE += ph->w * ratio_synch;
  spect[ic][ix2][iE].dEdlE += ph->w * ph->E * ratio_synch;
  spect[ic][ix2][iE].tau_abs += ph->w * ph->tau_abs * ratio_synch;
  spect[ic][ix2][iE].tau_scatt += ph->w * ph->tau_scatt * ratio_synch;
  spect[ic][ix2][iE].X1iav += ph->w * ph->X1i * ratio_synch;
  spect[ic][ix2][iE].X2isq += ph->w * (ph->X2i * ph->X2i) * ratio_synch;
  spect[ic][ix2][iE].X3fsq += ph->w * (ph->X[3] * ph->X[3]) * ratio_synch;
  spect[ic][ix2][iE].ne0 += ph->w * (ph->ne0) * ratio_synch;
  spect[ic][ix2][iE].b0 += ph->w * (ph->b0) * ratio_synch;
  spect[ic][ix2][iE].thetae0 += ph->w * (ph->thetae0) * ratio_synch;
  spect[ic][ix2][iE].nscatt += ph->w * ph->nscatt * ratio_synch;
  spect[ic][ix2][iE].nph += 1.;
  // .. to brems spectrum
  spect[ic+(N_COMPTBINS+1)][ix2][iE].dNdlE += ph->w * ph->ratio_brems;
  spect[ic+(N_COMPTBINS+1)][ix2][iE].dEdlE += ph->w * ph->E * ph->ratio_brems;
  spect[ic+(N_COMPTBINS+1)][ix2][iE].tau_abs += ph->w * ph->tau_abs * ph->ratio_brems;
  spect[ic+(N_COMPTBINS+1)][ix2][iE].tau_scatt += ph->w * ph->tau_scatt * ph->ratio_brems;
  spect[ic+(N_COMPTBINS+1)][ix2][iE].X1iav += ph->w * ph->X1i * ph->ratio_brems;
  spect[ic+(N_COMPTBINS+1)][ix2][iE].X2isq += ph->w * (ph->X2i * ph->X2i) * ph->ratio_brems;
  spect[ic+(N_COMPTBINS+1)][ix2][iE].X3fsq += ph->w * (ph->X[3] * ph->X[3]) * ph->ratio_brems;
  spect[ic+(N_COMPTBINS+1)][ix2][iE].ne0 += ph->w * (ph->ne0) * ph->ratio_brems;
  spect[ic+(N_COMPTBINS+1)][ix2][iE].b0 += ph->w * (ph->b0) * ph->ratio_brems;
  spect[ic+(N_COMPTBINS+1)][ix2][iE].thetae0 += ph->w * (ph->thetae0) * ph->ratio_brems;
  spect[ic+(N_COMPTBINS+1)][ix2][iE].nscatt += ph->w * ph->nscatt * ph->ratio_brems;
  spect[ic+(N_COMPTBINS+1)][ix2][iE].nph += 1.;
}

struct of_spectrum shared_spect[N_TYPEBINS][N_THBINS][N_EBINS] = { };

void omp_reduce_spect()
{
  #pragma omp parallel
  {
    #pragma omp critical
    {
      for (int k = 0; k < N_TYPEBINS; k++) {
        for (int i = 0; i < N_THBINS; i++) {
          for (int j = 0; j < N_EBINS; j++) {
            shared_spect[k][i][j].dNdlE +=
                spect[k][i][j].dNdlE;
            shared_spect[k][i][j].dEdlE +=
                spect[k][i][j].dEdlE;
            shared_spect[k][i][j].tau_abs +=
                spect[k][i][j].tau_abs;
            shared_spect[k][i][j].tau_scatt +=
                spect[k][i][j].tau_scatt;
            shared_spect[k][i][j].X1iav +=
                spect[k][i][j].X1iav;
            shared_spect[k][i][j].X2isq +=
                spect[k][i][j].X2isq;
            shared_spect[k][i][j].X3fsq +=
                spect[k][i][j].X3fsq;
            shared_spect[k][i][j].ne0 += 
                spect[k][i][j].ne0;
            shared_spect[k][i][j].b0 += 
                spect[k][i][j].b0;
            shared_spect[k][i][j].thetae0 +=
                spect[k][i][j].thetae0;
            shared_spect[k][i][j].nscatt +=
                spect[k][i][j].nscatt;
            shared_spect[k][i][j].nph += 
                spect[k][i][j].nph;
          }
        }
      }
    } // omp critical

    #pragma omp barrier

    #pragma omp master
    {
      for (int k = 0; k < N_TYPEBINS; k++) {
        for (int i = 0; i < N_THBINS; i++) {
          for (int j = 0; j < N_EBINS; j++) {
            spect[k][i][j].dNdlE =
                shared_spect[k][i][j].dNdlE;
            spect[k][i][j].dEdlE =
                shared_spect[k][i][j].dEdlE;
            spect[k][i][j].tau_abs =
                shared_spect[k][i][j].tau_abs;
            spect[k][i][j].tau_scatt =
                shared_spect[k][i][j].tau_scatt;
            spect[k][i][j].X1iav =
                shared_spect[k][i][j].X1iav;
            spect[k][i][j].X2isq =
                shared_spect[k][i][j].X2isq;
            spect[k][i][j].X3fsq =
                shared_spect[k][i][j].X3fsq;
            spect[k][i][j].ne0 = 
                shared_spect[k][i][j].ne0;
            spect[k][i][j].b0 = 
                shared_spect[k][i][j].b0;
            spect[k][i][j].thetae0 =
                shared_spect[k][i][j].thetae0;
            spect[k][i][j].nscatt =
                shared_spect[k][i][j].nscatt;
            spect[k][i][j].nph = 
                shared_spect[k][i][j].nph;
          }
        }
      }
    } // omp master
  } // omp parallel
}

double bias_func(double Te, double w)
{
  /*
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
   */

  // use old method with bias tuning parameter
  double bias, max;

  max = 0.5 * w / WEIGHT_MIN;

  if (Te > SCATTERING_THETAE_MAX) Te = SCATTERING_THETAE_MAX;
  bias = 16. * Te * Te / (5. * max_tau_scatt);

  if (bias > max) bias = max;

  return bias * biasTuning;
}

double thetae_func(double uu, double rho, double B, double kel)
{
  // assumes uu, rho, B, kel in code units
  double thetae = 0.;

  if (with_electrons == 0) {
    // fixed tp/te ratio
    thetae = uu / rho * Thetae_unit; 
  } else if (with_electrons == 1) {
    // howes/kawazura model from IHARM electron thermodynamics
    thetae = kel * pow(rho, game-1.) * Thetae_unit;
  } else if (with_electrons == 2 ) {
    double beta = uu * (gam-1.) / 0.5 / B / B;
    double b2 = beta*beta / beta_crit/beta_crit;
    double trat = trat_large * b2/(1.+b2) + trat_small /(1.+b2);
    thetae = (MP/ME) * (game-1.) * (gamp-1.) / ( (gamp-1.) + (game-1.)*trat ) * uu / rho;
  }

  return thetae;
}

void get_fluid_zone(int i, int j, int k, double *Ne, double *Thetae, double *B,
        double Ucon[NDIM], double Bcon[NDIM])
{

  double Ucov[NDIM], Bcov[NDIM];
  double Bp[NDIM], Vcon[NDIM], Vfac, VdotV, UdotBp;
  double sig ;

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

  *Ne = p[KRHO][i][j][k] * Ne_unit;
  *Thetae = thetae_func(p[UU][i][j][k], p[KRHO][i][j][k], (*B)/B_unit, p[KEL][i][j][k]);

  if (*Thetae > THETAE_MAX) *Thetae = THETAE_MAX;

  sig = pow(*B/B_unit,2)/(*Ne/Ne_unit);
  if(sig > 1. || i < 9) {
    *Thetae = SMALL;
  }

}

void get_fluid_params(double X[NDIM], double gcov[NDIM][NDIM], double *Ne,
          double *Thetae, double *B, double Ucon[NDIM],
          double Ucov[NDIM], double Bcon[NDIM],
          double Bcov[NDIM])
{
  double rho, kel, uu;
  double Bp[NDIM], Vcon[NDIM], Vfac, VdotV, UdotBp;
  double gcon[NDIM][NDIM];
  double interp_scalar(double X[NDIM], double ***var);
  double sig ;

  if ( X_in_domain(X) == 0 ) {
    *Ne = 0.;
    return;
  }

  rho = interp_scalar(X, p[KRHO]);
  kel = interp_scalar(X, p[KEL]);
  uu = interp_scalar(X, p[UU]);

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

  *Ne = rho*Ne_unit;
  *Thetae = thetae_func(uu, rho, (*B)/B_unit, kel);
  if (*Thetae > THETAE_MAX) *Thetae = THETAE_MAX ;

  sig = pow(*B/B_unit,2)/(*Ne/Ne_unit);
  if (sig > 1.) *Thetae = SMALL;
}

////////////////////////////////// COORDINATES /////////////////////////////////

int X_in_domain(double X[NDIM]) {
  // returns 1 if X is within the computational grid.
  // checks different sets of coordinates depending on
  // specified grid coordinates

  if (METRIC_eKS) {
    double XG[4] = { 0 };
    double Xks[4] = { X[0], exp(X[1]), M_PI*X[2], X[3] };

    if (METRIC_MKS3) {
      // if METRIC_MKS3, ignore theta boundaries
      double H0 = mks3H0, MY1 = mks3MY1, MY2 = mks3MY2, MP0 = mks3MP0;
      double KSx1 = Xks[1], KSx2 = Xks[2];
      XG[0] = Xks[0];
      XG[1] = log(Xks[1] - mks3R0);
      XG[2] = (-(H0*pow(KSx1,MP0)*M_PI) - pow(2,1 + MP0)*H0*MY1*M_PI +
        2*H0*pow(KSx1,MP0)*MY1*M_PI + pow(2,1 + MP0)*H0*MY2*M_PI +
        2*pow(KSx1,MP0)*atan(((-2*KSx2 + M_PI)*tan((H0*M_PI)/2.))/M_PI))/(2.*
        H0*(-pow(KSx1,MP0) - pow(2,1 + MP0)*MY1 + 2*pow(KSx1,MP0)*MY1 +
          pow(2,1 + MP0)*MY2)*M_PI);
      XG[3] = Xks[3];

      if (XG[1] < startx[1] || XG[1] > stopx[1]) return 0;
    }

  } else {
    if(X[1] < startx[1] ||
       X[1] > stopx[1]  ||
       X[2] < startx[2] ||
       X[2] > stopx[2]) {
      return 0;
    }
  }

  return 1;
}


/*
 *  translates geodesic coordinates to a grid zone and returns offset
 *  for interpolation purposes. integer index corresponds to the zone
 *  center "below" the desired point and del[i] \in [0,1) returns the
 *  offset from that zone center. note that the indices and deltas at
 *  the edges are different for x1,x2,x3 and compared to ipole.
 *
 *  0    0.5    1
 *  [     |     ]
 *  A  B  C DE  F
 *
 *  startx = 0.
 *  dx = 0.5
 *
 *  A -> ( 0, 0.0)  [ x1, x2 ]
 *  A -> ( 1, 0.5)  [ x3 ]
 *  B -> ( 0, 0.0)
 *  C -> ( 0, 0.5)
 *  D -> ( 0, 0.9)
 *  E -> ( 1, 0.0)
 *  F -> ( 1, 0.0)  [ x1, x2 ]
 *  F -> ( 1, 0.5)  [ x3 ]
 *
 */
void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM])
{
  // unless we're reading from data, i,j,k are the normal expected thing
  double phi;
  double XG[4];

  if (METRIC_eKS) {
    // the geodesics are evolved in eKS so invert through KS -> zone coordinates
    double Xks[4] = { X[0], exp(X[1]), M_PI*X[2], X[3] };
    if (METRIC_MKS3) {
      double H0 = mks3H0, MY1 = mks3MY1, MY2 = mks3MY2, MP0 = mks3MP0;
      double KSx1 = Xks[1], KSx2 = Xks[2];
      XG[0] = Xks[0];
      XG[1] = log(Xks[1] - mks3R0);
      XG[2] = (-(H0*pow(KSx1,MP0)*M_PI) - pow(2.,1. + MP0)*H0*MY1*M_PI + 
        2.*H0*pow(KSx1,MP0)*MY1*M_PI + pow(2.,1. + MP0)*H0*MY2*M_PI + 
        2.*pow(KSx1,MP0)*atan(((-2.*KSx2 + M_PI)*tan((H0*M_PI)/2.))/M_PI))/(2.*
        H0*(-pow(KSx1,MP0) - pow(2.,1 + MP0)*MY1 + 2.*pow(KSx1,MP0)*MY1 + 
          pow(2.,1. + MP0)*MY2)*M_PI);
      XG[3] = Xks[3];
    }
  } else {
    MULOOP XG[mu] = X[mu];
  }

  // the X[3] coordinate is allowed to vary so first map it to [0, stopx[3])
  phi = fmod(XG[3], stopx[3]);
  if (phi < 0.0) phi = stopx[3]+phi;

  // get provisional zone index. see note above function for details. note we
  // shift to zone centers because that's where variables are most exact.
  *i = (int) ((XG[1] - startx[1]) / dx[1] - 0.5 + 1000) - 1000;
  *j = (int) ((XG[2] - startx[2]) / dx[2] - 0.5 + 1000) - 1000;
  *k = (int) ((phi  - startx[3]) / dx[3] - 0.5 + 1000) - 1000;  

  // don't allow "center zone" to be outside of [0,N*-1]. this will often fire
  // for exotic corodinate systems and occasionally for normal ones. wrap x3.
  if (*i < 0) *i = 0;
  if (*j < 0) *j = 0;
  if (*k < 0) *k = 0;
  if (*i > N1-2) *i = N1-2; 
  if (*j > N2-2) *j = N2-2; 
  if (*k > N3-1) *k = N3-1; 

  // now construct del
  del[1] = (XG[1] - ((*i + 0.5) * dx[1] + startx[1])) / dx[1];
  del[2] = (XG[2] - ((*j + 0.5) * dx[2] + startx[2])) / dx[2];
  del[3] = (phi - ((*k + 0.5) * dx[3] + startx[3])) / dx[3];

  // finally enforce limits on del
  if (del[1] > 1.) del[1] = 1.;
  if (del[1] < 0.) del[1] = 0.;
  if (del[2] > 1.) del[2] = 1.;
  if (del[2] < 0.) del[2] = 0.;
  if (del[3] > 1.) del[3] = 1.;
  if (del[3] < 0.) {
    int oldk = *k;
    *k = N3-1;
    del[3] += 1.;
    if (del[3] < 0) {
      fprintf(stderr, " ! unable to resolve X[3] coordinate to zone %d %d %g %g\n", oldk, *k, del[3], XG[3]);
      exit(-7);
    }
  }
}

// return geodesic coordinates associated with center of zone i,j,k
void ijktoX(int i, int j, int k, double *X)
{
  // first do the naive thing 
  X[1] = startx[1] + (i+0.5)*dx[1];
  X[2] = startx[2] + (j+0.5)*dx[2];
  X[3] = startx[3] + (k+0.5)*dx[3];

  // now transform to geodesic coordinates if necessary by first
  // converting to KS and then to destination coordinates (eKS).
  if (METRIC_eKS) {
      double xKS[4] = { 0 };
    if (METRIC_MKS3) {
      double x0 = X[0];
      double x1 = X[1];
      double x2 = X[2];
      double x3 = X[3];

      double H0 = mks3H0;
      double MY1 = mks3MY1;
      double MY2 = mks3MY2;
      double MP0 = mks3MP0;
      
      xKS[0] = x0;
      xKS[1] = exp(x1) + mks3R0;
      xKS[2] = (M_PI*(1+1./tan((H0*M_PI)/2.)*tan(H0*M_PI*(-0.5+(MY1+(pow(2,MP0)*(-MY1+MY2))/pow(exp(x1)+R0,MP0))*(1-2*x2)+x2))))/2.;
      xKS[3] = x3;
    }
    
    X[0] = xKS[0];
    X[1] = log(xKS[1]);
    X[2] = xKS[2] / M_PI;
    X[3] = xKS[3];
  }
}

void set_dxdX(double X[NDIM], double dxdX[NDIM][NDIM])
{
  set_dxdX_metric(X, dxdX, 0);
}

static inline __attribute__((always_inline)) void set_dxdX_metric(double X[NDIM], double dxdX[NDIM][NDIM], int metric)
{
  // Jacobian with respect to KS basis where X is given in
  // non-KS basis
  MUNULOOP dxdX[mu][nu] = 0.;

  if ( METRIC_eKS && metric==0 ) {

    MUNULOOP dxdX[mu][nu] = mu==nu ? 1 : 0;
    dxdX[1][1] = exp(X[1]);
    dxdX[2][2] = M_PI;

  } else if ( METRIC_MKS3 ) {
   
    // mks3 ..
    dxdX[0][0] = 1.;
    dxdX[1][1] = exp(X[1]);
    dxdX[2][1] = -(pow(2.,-1. + mks3MP0)*exp(X[1])*mks3H0*mks3MP0*(mks3MY1 - 
              mks3MY2)*pow(M_PI,2)*pow(exp(X[1]) + mks3R0,-1 - mks3MP0)*(-1 + 
              2*X[2])*1./tan((mks3H0*M_PI)/2.)*pow(1./cos(mks3H0*M_PI*(-0.5 + (mks3MY1 + 
              (pow(2,mks3MP0)*(-mks3MY1 + mks3MY2))/pow(exp(X[1]) + mks3R0,mks3MP0))*(1 - 
              2*X[2]) + X[2])),2));
    dxdX[2][2]= (mks3H0*pow(M_PI,2)*(1 - 2*(mks3MY1 + (pow(2,mks3MP0)*(-mks3MY1 +
             mks3MY2))/pow(exp(X[1]) + mks3R0,mks3MP0)))*1./tan((mks3H0*M_PI)/2.)*
             pow(1./cos(mks3H0*M_PI*(-0.5 + (mks3MY1 + (pow(2,mks3MP0)*(-mks3MY1 +
             mks3MY2))/pow(exp(X[1]) + mks3R0,mks3MP0))*(1 - 2*X[2]) + X[2])),2))/2.;
    dxdX[3][3] = 1.;

  } else if ( with_derefine_poles ) {

    // mmks
    dxdX[0][0] = 1.;
    dxdX[1][1] = exp(X[1]);
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
    dxdX[3][3] = 1.;

  } else {

    // mks
    dxdX[0][0] = 1.;
    dxdX[1][1] = exp(X[1]);
    dxdX[2][2] = M_PI - (hslope - 1.)*M_PI*cos(2.*M_PI*X[2]);
    dxdX[3][3] = 1.;

  }

}

void gcov_func(double X[NDIM], double gcov[NDIM][NDIM])
{
  MUNULOOP gcov[mu][nu] = 0.;

  double r, th;

  // despite the name, get equivalent values for
  // r, th for KS coordinates
  bl_coord(X, &r, &th);

  // compute ks metric
  gcov_ks(r, th, gcov);

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

// compute KS metric at point (r,th) in KS coordinates (cyclic in t, ph)
static inline __attribute__((always_inline)) void gcov_ks(double r, double th, double gcov[NDIM][NDIM])
{
  double cth = cos(th);
  double sth = sin(th);

  double s2 = sth*sth;
  double rho2 = r*r + a*a*cth*cth;

  // compute ks metric for ks coordinates (cyclic in t,phi)
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
}

// return the gdet associated with zone coordinates for the zone at
// i,j,k
double gdet_zone(int i, int j, int k) 
{
  // get the X for the zone (in geodesic coordinates for bl_coord) 
  // and in zone coordinates (for set_dxdX_metric)
  double X[NDIM], Xzone[NDIM];
  ijktoX(i,j,k, X);
  Xzone[0] = 0.;
  Xzone[1] = startx[1] + (i+0.5)*dx[1];
  Xzone[2] = startx[2] + (j+0.5)*dx[2];
  Xzone[3] = startx[3] + (k+0.5)*dx[3];

  // then get gcov for the zone (in zone coordinates)
  double gcovKS[NDIM][NDIM], gcov[NDIM][NDIM];
  double r, th;
  double dxdX[NDIM][NDIM];
  MUNULOOP gcovKS[mu][nu] = 0.;
  MUNULOOP gcov[mu][nu] = 0.;
  bl_coord(X, &r, &th);
  gcov_ks(r, th, gcovKS);
  set_dxdX_metric(Xzone, dxdX, 1);
  MUNULOOP {
    for (int lam=0; lam<NDIM; ++lam) {
      for (int kap=0; kap<NDIM; ++kap) {
        gcov[mu][nu] += gcovKS[lam][kap]*dxdX[lam][mu]*dxdX[kap][nu];
      }
    }
  }

  return gdet_func(gcov); 
}

// returns BL.{r,th} == KS.{r,th} of point with geodesic coordinates X
void bl_coord(double *X, double *r, double *th)
{
  *r = exp(X[1]);

  if (METRIC_eKS) {
    *r = exp(X[1]);
    *th = M_PI * X[2];
  } else if (METRIC_MKS3) {
    *r = exp(X[1]) + mks3R0;
    *th = (M_PI*(1. + 1./tan((mks3H0*M_PI)/2.)*tan(mks3H0*M_PI*(-0.5 + (mks3MY1 + (pow(2.,mks3MP0)*(-mks3MY1 + mks3MY2))/pow(exp(X[1])+mks3R0,mks3MP0))*(1. - 2.*X[2]) + X[2]))))/2.;
  } else if (with_derefine_poles) {
    double thG = M_PI*X[2] + ((1. - hslope)/2.)*sin(2.*M_PI*X[2]);
    double y = 2*X[2] - 1.;
    double thJ = poly_norm*y*(1. + pow(y/poly_xt,poly_alpha)/(poly_alpha+1.)) + 0.5*M_PI;
    *th = thG + exp(mks_smooth*(startx[1] - X[1]))*(thJ - thG);
  } else {
    *th = M_PI*X[2] + ((1. - hslope)/2.)*sin(2.*M_PI*X[2]);
  }
}

// warning: this function assumes startx = 0 and stopx = 1 (that we bin evenly in BL)
double dOmega_func(int j)
{
  double dx2 = M_PI/2./N_THBINS;
  double thi = j * dx2;
  double thf = (j+1) * dx2;

  return 2.*M_PI*(-cos(thf) + cos(thi));
}

//////////////////////////////// INITIALIZATION ////////////////////////////////

#include <hdf5.h>
#include <hdf5_hl.h>

void init_data(int argc, char *argv[], Params *params)
{
  const char *fname = NULL;
  double dV, V;
  int nprims = 0;

  NPRIM = 10;

  if (params->loaded && strlen(params->dump) > 0) {
    fname = params->dump;
    trat_small = params->trat_small;
    trat_large = params->trat_large;
    beta_crit = params->beta_crit;
    biasTuning = params->biasTuning;
  } else {
    fname = argv[2];
    strncpy((char *)params->dump, argv[2], 255);
  }

  if ( hdf5_open((char *)fname) < 0 ) {
    fprintf(stderr, "File %s does not exist! Exiting...\n", fname);
    exit(-1);
  }

  // get dump info to copy to grmonty output
  fluid_header = hdf5_get_blob("/header");

  // read header
  hdf5_set_directory("/header/");

  // flag reads
  with_electrons = 0;
  with_radiation = 0;
  with_derefine_poles = 0;
  if ( hdf5_exists("has_electrons") ) 
    hdf5_read_single_val(&with_electrons, "has_electrons", H5T_STD_I32LE);
  if ( hdf5_exists("has_radiation") ) 
    hdf5_read_single_val(&with_radiation, "has_radiation", H5T_STD_I32LE);

  // read geometry
  with_derefine_poles = 0;
  METRIC_MKS3 = 0;
  char metric_name[20];
  hid_t string_type = hdf5_make_str_type(20);
  hdf5_read_single_val(&metric_name, "metric", string_type);
  if ( strncmp(metric_name, "MMKS", 19) == 0 ) {
    with_derefine_poles = 1;
  } else if ( strncmp(metric_name, "MKS3", 19) == 0 ) {
    METRIC_eKS = 1;
    METRIC_MKS3 = 1;
    fprintf(stderr, "using eKS metric with exotic \"MKS3\" zones...\n");
  }

  hdf5_read_single_val(&nprims, "n_prim", H5T_STD_I32LE);
  hdf5_read_single_val(&N1, "n1", H5T_STD_I32LE);
  hdf5_read_single_val(&N2, "n2", H5T_STD_I32LE);
  hdf5_read_single_val(&N3, "n3", H5T_STD_I32LE);
  hdf5_read_single_val(&gam, "gam", H5T_IEEE_F64LE);

  // conditional reads
  game = 4./3;
  gamp = 5./3;
  if (with_electrons) {
    fprintf(stderr, "custom electron model loaded...\n");
    hdf5_read_single_val(&game, "gam_e", H5T_IEEE_F64LE);
    hdf5_read_single_val(&gamp, "gam_p", H5T_IEEE_F64LE);
  }

  if (!USE_FIXED_TPTE && !USE_MIXED_TPTE) {
    if (with_electrons != 1) {
      fprintf(stderr, "! no electron temperature model specified in model/iharm.c\n");
      exit(-3);
    }
    with_electrons = 1;
    Thetae_unit = MP/ME;
  } else if (USE_FIXED_TPTE && !USE_MIXED_TPTE) {
    with_electrons = 0; // force TP_OVER_TE to overwrite electron temperatures
    fprintf(stderr, "using fixed tp_over_te ratio = %g\n", tp_over_te);
    Thetae_unit = MP/ME * (gam-1.) / (1. + tp_over_te);
    Thetae_unit = 2./3. * MP/ME / (2. + tp_over_te);
  } else if (USE_MIXED_TPTE && !USE_FIXED_TPTE) {
    Thetae_unit = 2./3. * MP/ME / 5.; 
    with_electrons = 2;
    fprintf(stderr, "using mixed tp_over_te with trat_small = %g and trat_large = %g\n", trat_small, trat_large);
  } else {
    fprintf(stderr, "! please change electron model in model/iharm.c\n");
    exit(-3);
  }

  if (with_radiation) {
    fprintf(stderr, "custom radiation field tracking information loaded...\n");
    hdf5_set_directory("/header/units/");
    hdf5_read_single_val(&M_unit, "M_unit", H5T_IEEE_F64LE);
    hdf5_read_single_val(&T_unit, "T_unit", H5T_IEEE_F64LE);
    hdf5_read_single_val(&L_unit, "L_unit", H5T_IEEE_F64LE);
    if (!USE_FIXED_TPTE && !USE_MIXED_TPTE) {
      hdf5_read_single_val(&Thetae_unit, "Thetae_unit", H5T_IEEE_F64LE);
    }
    hdf5_read_single_val(&MBH, "Mbh", H5T_IEEE_F64LE);
    hdf5_read_single_val(&TP_OVER_TE, "tp_over_te", H5T_IEEE_F64LE);
  } else {
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
  }

  hdf5_set_directory("/header/geom/");
  hdf5_read_single_val(&startx[1], "startx1", H5T_IEEE_F64LE);
  hdf5_read_single_val(&startx[2], "startx2", H5T_IEEE_F64LE);
  hdf5_read_single_val(&startx[3], "startx3", H5T_IEEE_F64LE);
  hdf5_read_single_val(&dx[1], "dx1", H5T_IEEE_F64LE);
  hdf5_read_single_val(&dx[2], "dx2", H5T_IEEE_F64LE);
  hdf5_read_single_val(&dx[3], "dx3", H5T_IEEE_F64LE);

  hdf5_set_directory("/header/geom/mks/");
  if (with_derefine_poles) hdf5_set_directory("/header/geom/mmks/");
  if ( METRIC_MKS3 ) {
    hdf5_set_directory("/header/geom/mks3/");
    hdf5_read_single_val(&a, "a", H5T_IEEE_F64LE);
    hdf5_read_single_val(&mks3R0, "R0", H5T_IEEE_F64LE);
    hdf5_read_single_val(&mks3H0, "H0", H5T_IEEE_F64LE);
    hdf5_read_single_val(&mks3MY1, "MY1", H5T_IEEE_F64LE);
    hdf5_read_single_val(&mks3MY2, "MY2", H5T_IEEE_F64LE);
    hdf5_read_single_val(&mks3MP0, "MP0", H5T_IEEE_F64LE);
    Rout = 100.; 
  } else {
    hdf5_read_single_val(&a, "a", H5T_IEEE_F64LE);
    hdf5_read_single_val(&hslope, "hslope", H5T_IEEE_F64LE);
    hdf5_read_single_val(&Rin, "Rin", H5T_IEEE_F64LE);
    hdf5_read_single_val(&Rout, "Rout", H5T_IEEE_F64LE);
    if (with_derefine_poles) {
      fprintf(stderr, "custom refinement at poles loaded...\n");
      hdf5_read_single_val(&poly_xt, "poly_xt", H5T_IEEE_F64LE);
      hdf5_read_single_val(&poly_alpha, "poly_alpha", H5T_IEEE_F64LE);
      hdf5_read_single_val(&mks_smooth, "mks_smooth", H5T_IEEE_F64LE);
      poly_norm = 0.5*M_PI*1./(1. + 1./(poly_alpha + 1.)*1./pow(poly_xt, poly_alpha));
    }
  }

  // Set other geometry
  stopx[0] = 1.;
  stopx[1] = startx[1]+N1*dx[1];
  stopx[2] = startx[2]+N2*dx[2];
  stopx[3] = startx[3]+N3*dx[3];

  // Set remaining units and constants
  RHO_unit = M_unit/pow(L_unit,3);
  U_unit = RHO_unit*CL*CL;
  B_unit = CL*sqrt(4.*M_PI*RHO_unit);
  Ne_unit = RHO_unit/(MP + ME);
  max_tau_scatt = (6.*L_unit)*RHO_unit*0.4; // this doesn't make sense ... 
  max_tau_scatt = 0.0001; // TODO look at this in the future and figure out a smarter general value

  // Horizon and "max R for geodesic tracking" in KS coordinates
  Rh = 1. + sqrt(1. - a * a);
  Rmax = 1000;

  fprintf(stderr, "L_unit, T_unit, M_unit = %g %g %g\n", L_unit, T_unit, M_unit);
  fprintf(stderr, "B_unit, Ne_unit, RHO_unit = %g %g %g\n", B_unit, Ne_unit, RHO_unit);
  fprintf(stderr, "Thetae_unit = %g\n", Thetae_unit);

  // Allocate storage and set geometry
  double ****malloc_rank4_double(int n1, int n2, int n3, int n4);
  p = malloc_rank4_double(NVAR, N1, N2, N3);
  fprintf(stderr, "NVAR N1 N2 N3 = %i %i %i %i\n", NVAR, N1, N2, N3);
  n2gens = (double ***)malloc_rank3(N1, N2, N3, sizeof(double));
  geom = (struct of_geom**)malloc_rank2(N1, N2, sizeof(struct of_geom));
  tetrads = (struct of_tetrads***)malloc_rank3(N1, N2, N3, sizeof(struct of_tetrads));
  init_geometry();

  // Read prims. 
  // Assume standard ordering in iharm dump file, especially for
  // electron variables...
  hdf5_set_directory("/");

  hsize_t fdims[] = { N1, N2, N3, nprims };
  hsize_t fstart[] = { 0, 0, 0, 0 }; //{global_start[0], global_start[1], global_start[2], 0};
  hsize_t fcount[] = {N1, N2, N3, 1};
  hsize_t mstart[] = {0, 0, 0, 0};

  fstart[3] = 0;
  hdf5_read_array(p[KRHO][0][0], "prims", 4, fdims, fstart, fcount, fcount, mstart, H5T_IEEE_F64LE);
  fstart[3] = 1;
  hdf5_read_array(p[UU][0][0], "prims", 4, fdims, fstart, fcount, fcount, mstart, H5T_IEEE_F64LE);
  fstart[3] = 2;
  hdf5_read_array(p[U1][0][0], "prims", 4, fdims, fstart, fcount, fcount, mstart, H5T_IEEE_F64LE);
  fstart[3] = 3;
  hdf5_read_array(p[U2][0][0], "prims", 4, fdims, fstart, fcount, fcount, mstart, H5T_IEEE_F64LE);
  fstart[3] = 4;
  hdf5_read_array(p[U3][0][0], "prims", 4, fdims, fstart, fcount, fcount, mstart, H5T_IEEE_F64LE);
  fstart[3] = 5;
  hdf5_read_array(p[B1][0][0], "prims", 4, fdims, fstart, fcount, fcount, mstart, H5T_IEEE_F64LE);
  fstart[3] = 6;
  hdf5_read_array(p[B2][0][0], "prims", 4, fdims, fstart, fcount, fcount, mstart, H5T_IEEE_F64LE);
  fstart[3] = 7;
  hdf5_read_array(p[B3][0][0], "prims", 4, fdims, fstart, fcount, fcount, mstart, H5T_IEEE_F64LE);

  if (with_electrons == 1) {

    fstart[3] = 8;
    hdf5_read_array(p[KEL][0][0], "prims", 4, fdims, fstart, fcount, fcount, mstart, H5T_IEEE_F64LE);

    fstart[3] = 9;
    hdf5_read_array(p[KTOT][0][0], "prims", 4, fdims, fstart, fcount, fcount, mstart, H5T_IEEE_F64LE);

  }

  hdf5_close();

  V = dMact = Ladv = 0.;
  dV = dx[1]*dx[2]*dx[3];
  ZLOOP {
    V += dV*geom[i][j].gzone;
    bias_norm += dV*geom[i][j].gzone*pow(p[UU][i][j][k]/p[KRHO][i][j][k]*Thetae_unit,2.);

    if (10 <= i && i <= 20) {
      double Ne, Thetae, Bmag, Ucon[NDIM], Ucov[NDIM], Bcon[NDIM];
      get_fluid_zone(i, j, k, &Ne, &Thetae, &Bmag, Ucon, Bcon);
      lower(Ucon, geom[i][j].gcov, Ucov);
      dMact += geom[i][j].gzone*dx[2]*dx[3]*p[KRHO][i][j][k]*Ucon[1];
      Ladv += geom[i][j].gzone*dx[2]*dx[3]*p[UU][i][j][k]*Ucon[1]*Ucov[0];
    }

  }

  /*
  for (int i=0; i<100; ++i) {
    double Xp[NDIM];
    ijktoX(i,64,64,Xp);
    fprintf(stderr, "%d -> %g %g %g %g\n", i, Xp[0], exp(Xp[1]), Xp[2], Xp[3]);
  }
   */

  dMact /= 11.;
  Ladv /= 1.;
  bias_norm /= V;
  fprintf(stderr, "dMact: %g, Ladv: %g\n", dMact, Ladv);

  init_tetrads();
}

//////////////////////////////////// OUTPUT ////////////////////////////////////

void report_spectrum(int N_superph_made, Params *params)
{

  hid_t fid = -1;

  if (params->loaded && strlen(params->spectrum) > 0) {
    fid = H5Fcreate(params->spectrum, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  } else {
    fid = H5Fcreate("spectrum.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  }

  if (fid < 0) {
    fprintf(stderr, "! unable to open/create spectrum hdf5 file.\n");
    exit(-3);
  }

  h5io_add_attribute_str(fid, "/", "githash", xstr(VERSION));

  h5io_add_blob(fid, "/fluid_header", fluid_header);
  hdf5_close_blob(fluid_header);
   
  h5io_add_group(fid, "/params");

  h5io_add_data_dbl(fid, "/params/NUCUT", NUCUT);
  h5io_add_data_dbl(fid, "/params/GAMMACUT", GAMMACUT);
  h5io_add_data_dbl(fid, "/params/NUMAX", NUMAX);
  h5io_add_data_dbl(fid, "/params/NUMIN", NUMIN);
  h5io_add_data_dbl(fid, "/params/LNUMAX", LNUMAX);
  h5io_add_data_dbl(fid, "/params/LNUMIN", LNUMIN);
  h5io_add_data_dbl(fid, "/params/DLNU", DLNU);
  h5io_add_data_dbl(fid, "/params/THETAE_MIN", THETAE_MIN);
  h5io_add_data_dbl(fid, "/params/THETAE_MAX", THETAE_MAX);
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
  h5io_add_data_str(fid, "/params/model", xstr(MODEL));

  h5io_add_group(fid, "/params/electrons");
  if (with_electrons == 0) {
    h5io_add_data_dbl(fid, "/params/electrons/tp_over_te", tp_over_te);
  } else if (with_electrons == 2) {
    h5io_add_data_dbl(fid, "/params/electrons/rlow", trat_small);
    h5io_add_data_dbl(fid, "/params/electrons/rhigh", trat_large);
  }
  h5io_add_data_int(fid, "/params/electrons/type", with_electrons);

  // temporary data buffers
  double lnu_buf[N_EBINS];
  double dOmega_buf[N_THBINS];
  double nuLnu_buf[N_TYPEBINS][N_EBINS][N_THBINS];
  double tau_abs_buf[N_TYPEBINS][N_EBINS][N_THBINS];
  double tau_scatt_buf[N_TYPEBINS][N_EBINS][N_THBINS];
  double x1av_buf[N_TYPEBINS][N_EBINS][N_THBINS];
  double x2av_buf[N_TYPEBINS][N_EBINS][N_THBINS];
  double x3av_buf[N_TYPEBINS][N_EBINS][N_THBINS];
  double nscatt_buf[N_TYPEBINS][N_EBINS][N_THBINS];
  double Lcomponent_buf[N_TYPEBINS];

  // normal output routine
  double dOmega, nuLnu, tau_scatt, L, Lcomponent, dL;

  max_tau_scatt = 0.;
  L = 0.;
  dL = 0.;

  for (int j=0; j<N_THBINS; ++j) {
    // warning: this assumes geodesic X \in [0,1]
    dOmega_buf[j] = 2. * dOmega_func(j);
  }

  for (int k=0; k<N_TYPEBINS; ++k) {
    Lcomponent = 0.;
    for (int i=0; i<N_EBINS; ++i) {
      lnu_buf[i] = (i * dlE + lE0) / M_LN10;
      for (int j=0; j<N_THBINS; ++j) {

        dOmega = dOmega_buf[j];

        nuLnu = (ME * CL * CL) * (4. * M_PI / dOmega) * (1. / dlE);
        nuLnu *= spect[k][j][i].dEdlE/LSUN;

        tau_scatt = spect[k][j][i].tau_scatt/(spect[k][j][i].dNdlE + SMALL);

        nuLnu_buf[k][i][j] = nuLnu;
        tau_abs_buf[k][i][j] = spect[k][j][i].tau_abs/(spect[k][j][i].dNdlE + SMALL);
        tau_scatt_buf[k][i][j] = tau_scatt;
        x1av_buf[k][i][j] = spect[k][j][i].X1iav/(spect[k][j][i].dNdlE + SMALL);
        x2av_buf[k][i][j] = sqrt(fabs(spect[k][j][i].X2isq/(spect[k][j][i].dNdlE + SMALL)));
        x3av_buf[k][i][j] = sqrt(fabs(spect[k][j][i].X3fsq/(spect[k][j][i].dNdlE + SMALL)));
        nscatt_buf[k][i][j] = spect[k][j][i].nscatt / (spect[k][j][i].dNdlE + SMALL);

        if (tau_scatt > max_tau_scatt) max_tau_scatt = tau_scatt;

        dL += ME * CL * CL * spect[k][j][i].dEdlE;
        L += nuLnu * dOmega * dlE / (4. * M_PI);
        Lcomponent += nuLnu * dOmega * dlE / (4. * M_PI);
      }
    }
    Lcomponent_buf[k] = Lcomponent;
  }

  h5io_add_group(fid, "/output");

  h5io_add_data_dbl_1d(fid, "/output/lnu", N_EBINS, lnu_buf);
  h5io_add_data_dbl_1d(fid, "/output/dOmega", N_THBINS, dOmega_buf);
  h5io_add_data_dbl_3d(fid, "/output/nuLnu", N_TYPEBINS, N_EBINS, N_THBINS, nuLnu_buf);
  h5io_add_data_dbl_3d(fid, "/output/tau_abs", N_TYPEBINS, N_EBINS, N_THBINS, tau_abs_buf);
  h5io_add_data_dbl_3d(fid, "/output/tau_scatt", N_TYPEBINS, N_EBINS, N_THBINS, tau_scatt_buf);
  h5io_add_data_dbl_3d(fid, "/output/x1av", N_TYPEBINS, N_EBINS, N_THBINS, x1av_buf);
  h5io_add_data_dbl_3d(fid, "/output/x2av", N_TYPEBINS, N_EBINS, N_THBINS, x2av_buf);
  h5io_add_data_dbl_3d(fid, "/output/x3av", N_TYPEBINS, N_EBINS, N_THBINS, x3av_buf);
  h5io_add_data_dbl_3d(fid, "/output/nscatt", N_TYPEBINS, N_EBINS, N_THBINS, nscatt_buf);
  h5io_add_data_dbl_1d(fid, "/output/Lcomponent", N_TYPEBINS, Lcomponent_buf);

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
  h5io_add_data_dbl(fid, "/output/efficiency", L * LSUN / (dMact * M_unit * CL * CL / T_unit));

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

