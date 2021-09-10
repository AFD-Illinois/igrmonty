#include "decs.h"
#include "coordinates.h"
#include "model_radiation.h"

#define USE_FIXED_TPTE (1)   // don't use for HAMR dataset
#define USE_MIXED_TPTE (0)   // don't use for HAMR dataset

double interp_scalar(const double X[NDIM], double ***var);

// electron model. these values will be overwritten by anything found in par.c
// or in the runtime parameter file.
// with_electrons ->
//     0 : constant TP_OVER_TE
//     1 : use dump file model (kawazura?)  -> not compatible with HAMR dataset
//     2: use mixed TP_OVER_TE (moscibrodzka "beta" model)
static double tp_over_te;
static double trat_small;
static double trat_large;
static double beta_crit;
static double Thetae_max;
static int with_electrons;

// fluid data
double ****bcon;
double ****bcov;
double ****ucon;
double ****ucov;
double ****p;
double ***ne;
double ***thetae;
double ***b;

double ***sigma_array;
double ***beta_array;

double TP_OVER_TE;

static double MBH, game, gamp;

static int with_radiation;

void readattr(hid_t file_id, const char *attr_name, hid_t mem_type_id, void *buf);
void readdata(hid_t file_id, const char *attr_name, hid_t mem_type_id, hid_t memspace, void *buf);

void report_bad_input(int argc)
{
  if (argc < 6) {
    fprintf(stderr, "usage: \n");
    fprintf(stderr, "  HAMR (read_dscale):   grmonty Ns fname \n");
    fprintf(stderr, "  HAMR (!read_dscale):  grmonty Ns fname Rho_unit \n");
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

  // bin in X[2] BL coord while folding around the equator and check limit
  double r, th;
  bl_coord(ph->X, &r, &th);
  dx2 = M_PI/2./N_THBINS;
  if (th > M_PI/2.) {
    ix2 = (int)( (M_PI - th) / dx2 );
  } else {
    ix2 = (int)( th / dx2 );
  }
  if (ix2 < 0 || ix2 >= N_THBINS) return;

  #if CUSTOM_AVG==1
  double nu = ph->E * ME*CL*CL / HPL;
  if (nu < CA_MIN_FREQ || CA_MAX_FREQ < nu) return;
  // Get custom average bin
  dlE = (log(CA_MAX_QTY) - log(CA_MIN_QTY)) / CA_NBINS;
  lE = log(ph->QTY0);
  iE = (int) ((lE - log(CA_MIN_QTY)) / dlE + 2.5) - 2;
  if (iE < 0 || iE >= CA_NBINS) return;
  #else
  // Get energy bin (centered on iE*dlE + lE0)
  lE = log(ph->E);
  iE = (int) ((lE - lE0) / dlE + 2.5) - 2;
  if (iE < 0 || iE >= N_EBINS) return;
  #endif // CUSTOM_AVG

  // Get compton bin
  ic = ph->nscatt;
  if (ic > 3) ic = 3;

  #pragma omp atomic
  N_superph_recorded++;

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

  //bias = 0; // test: turn off compton scattering
  return bias * biasTuning;
}

double thetae_func(double uu, double rho, double B)
{
  // assumes uu, rho, B in code units
  double thetae = 0.;

  if (with_electrons == 0) {  
    // fixed tp/te ratio
    thetae = MP/ME * (gam-1.) * uu / rho / tp_over_te;
  } else if (with_electrons == 1) {
    // howes/kawazura model from IHARM electron thermodynamics:: not included in HAMR dataset
    fprintf(stderr, "hows/kawazura model is not compatible with h-amr dumps..\n");
	exit(-3);
    //thetae = kel * pow(rho, game-1.) * Thetae_unit;
  } else if (with_electrons == 2 ) {
	// Moscibrodzka beta model for eletron temperature. 
    double beta = uu * (gam-1.) / 0.5 / B / B;
    double b2 = beta*beta / beta_crit/beta_crit;
    double trat = trat_large * b2/(1.+b2) + trat_small /(1.+b2);
    if (B == 0) trat = trat_large;
    thetae = (MP/ME) * (gam-1.) * uu / rho / trat;
  }

  //return 1./(1./thetae + 1./Thetae_max);
  return thetae;
}

/*
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
    if (B == 0) trat = trat_large;
    thetae = (MP/ME) * (game-1.) * (gamp-1.) / ( (gamp-1.) + (game-1.)*trat ) * uu / rho;
  }

  return 1./(1./thetae + 1./Thetae_max);
}
*/

void get_fluid_zone(int i, int j, int k, double *Ne, double *Thetae, double *B,
        double Ucon[NDIM], double Bcon[NDIM])
{

  double Ucov[NDIM], Bcov[NDIM];
  double Bp[NDIM], Vcon[NDIM], Vfac, VdotV, UdotBp;
  double sig ;

  //Bp[1] = p[B1][i][j][k];
  //Bp[2] = p[B2][i][j][k];
  //Bp[3] = p[B3][i][j][k];
  Bp[1] = p[B1][i][j][k] * sqrt(-geom[i][j].gcon[0][0]);
  Bp[2] = p[B2][i][j][k] * sqrt(-geom[i][j].gcon[0][0]);
  Bp[3] = p[B3][i][j][k] * sqrt(-geom[i][j].gcon[0][0]);
  
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
  //*Thetae = thetae_func(p[UU][i][j][k], p[KRHO][i][j][k], (*B)/B_unit, p[KEL][i][j][k]);
  //*Thetae = MP/ME * (gam-1.) * p[UU][i][j][k] / p[KRHO][i][j][k] / TP_OVER_TE;
  *Thetae = thetae_func(p[UU][i][j][k], p[KRHO][i][j][k], (*B)/B_unit);

  if (*Thetae > THETAE_MAX) *Thetae = THETAE_MAX;

  sig = pow(*B/B_unit,2)/(*Ne/Ne_unit);
  //if(sig > 1. || i < 9) {
  if(sig > 1.) {
    *Thetae = SMALL;
  }

}

double get_model_sigma(const double X[NDIM])
{
  return interp_scalar(X, sigma_array);
}

double get_model_beta(const double X[NDIM])
{
  return interp_scalar(X, beta_array);
}

void get_fluid_params(const double X[NDIM], double gcov[NDIM][NDIM], double *Ne,
          double *Thetae, double *B, double Ucon[NDIM],
          double Ucov[NDIM], double Bcon[NDIM],
          double Bcov[NDIM])
{
  double rho, kel, uu;
  double Bp[NDIM], Vcon[NDIM], Vfac, VdotV, UdotBp;
  double gcon[NDIM][NDIM];
  double interp_scalar(const double X[NDIM], double ***var);
  double sig ;

  if ( X_in_domain(X) == 0 ) {
    *Ne = 0.;
    return;
  }

  rho = interp_scalar(X, p[KRHO]);
  //kel = interp_scalar(X, p[KEL]);
  uu = interp_scalar(X, p[UU]);

  //Bp[1] = interp_scalar(X, p[B1]);
  //Bp[2] = interp_scalar(X, p[B2]);
  //Bp[3] = interp_scalar(X, p[B3]);

  Vcon[1] = interp_scalar(X, p[U1]);
  Vcon[2] = interp_scalar(X, p[U2]);
  Vcon[3] = interp_scalar(X, p[U3]);

  gcov_func(X, gcov);
  gcon_func(gcov, gcon);

  Bp[1] = interp_scalar(X, p[B1]) * sqrt(-gcon[0][0]);
  Bp[2] = interp_scalar(X, p[B2]) * sqrt(-gcon[0][0]);
  Bp[3] = interp_scalar(X, p[B3]) * sqrt(-gcon[0][0]);

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
  //*Thetae = thetae_func(uu, rho, (*B)/B_unit, kel);
  //*Thetae = MP/ME * (gam-1.) * uu / rho / TP_OVER_TE;
  *Thetae = thetae_func(uu, rho, (*B)/B_unit);
  if (*Thetae > THETAE_MAX) *Thetae = THETAE_MAX ;

  sig = pow(*B/B_unit,2)/(*Ne/Ne_unit);
  if (sig > 1.) *Thetae = SMALL;
} 

////////////////////////////////// COORDINATES /////////////////////////////////

void gcov_func(const double X[NDIM], double gcov[NDIM][NDIM])
{
  // despite the name, get equivalent values for
  // r, th for KS coordinates
  double r, th;
  bl_coord(X, &r, &th);

  // compute ks metric
  double gcovKS[NDIM][NDIM];
  gcov_ks(r, th, gcovKS);

  // Apply coordinate transformation to code coordinates X
  double dxdX[NDIM][NDIM];
  set_dxdX(X, dxdX);

  MUNULOOP {
    gcov[mu][nu] = 0.;
    for (int lam = 0; lam < NDIM; lam++) {
      for (int kap = 0; kap < NDIM; kap++) {
        gcov[mu][nu] += gcovKS[lam][kap]*dxdX[lam][mu]*dxdX[kap][nu];
      }
    }
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

  hid_t    file_id;        /* File identifier */
  hid_t    memspace;       /* memory space identifier */
  hsize_t  dimsm[1];       /* memory space dimensions */
  herr_t   ret;            /* Return value */

  //double t,a,gam,Rin,Rout,hslope,R0;
  //int N1, N2, N3;
  int RANK_OUT=1;          /* dimension of data array for HDF5 dataset */
  int gridIndex,gridIndex2D;

  double *x1_in,*x2_in,*x3_in,*r_in,*h_in,*ph_in,*RHO_in,*UU_in, //*U0_in,
         *U1_in,*U2_in,*U3_in,*B1_in,*B2_in,*B3_in,*gdet_in,*Ucov0_in,*Ucon0_in;

  int i, j, z, ieh;
  double x[4], xp[4];
  double rin, hin, phin, gdet, Ucov0, Ucon0, dscale;
  double rp, hp, x2temp;

  NPRIM = NVAR;

  if (params->loaded && strlen(params->dump) > 0) {
    fname = params->dump;
    trat_small = params->trat_small;
    trat_large = params->trat_large;
    beta_crit = params->beta_crit;
    biasTuning = params->biasTuning;
    Thetae_max = params->Thetae_max;
  } else {      // preferrable condition for H-AMR
    fname = argv[2];
    strncpy((char *)params->dump, argv[2], 255);

	#if (Monika_Te)    // values from model.h (only for hamr)
	trat_small = Rlow;
	trat_large = Rhigh;
	params->trat_small = trat_small;
	params->trat_large = trat_large;
	#else
	tp_over_te = TPoTE;
	params->TP_OVER_TE = TPoTE;
	#endif
	beta_crit = params->beta_crit;
	biasTuning = params->biasTuning;
	Thetae_max = THETAE_MAX;
	params->Thetae_max = Thetae_max;
  }

  if ( hdf5_open((char *)fname) < 0 ) {
    fprintf(stderr, "File %s does not exist! Exiting...\n", fname);
    exit(-1);
  } else {
    file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
	fprintf(stderr, "successfully opened HAMR dataset: %s\n", fname);
  }

  // flag reads
  #if (Monika_Te)
  with_electrons = 2;   // Monika's Te beta model
  #else
  with_electrons = 0;   // constant Temperature ratio
  #endif
  // keep the flag 0 for HAMR data set.
  with_radiation = 0;
  // read geometry
  with_derefine_poles = 0;
  METRIC_MKS3 = 0;
  METRIC_eKS  = 0;

  /* read attributes */
  readattr(file_id, "t",      H5T_NATIVE_DOUBLE, &t);
  readattr(file_id, "N1",     H5T_NATIVE_INT,    &N1);
  readattr(file_id, "N2",     H5T_NATIVE_INT,    &N2);
  readattr(file_id, "N3",     H5T_NATIVE_INT,    &N3);
  readattr(file_id, "startx", H5T_NATIVE_DOUBLE, &startx[1]);
  readattr(file_id, "dx",     H5T_NATIVE_DOUBLE, &dx[1]);
  readattr(file_id, "a",      H5T_NATIVE_DOUBLE, &a);
  readattr(file_id, "gam",    H5T_NATIVE_DOUBLE, &gam);
  readattr(file_id, "Rin",    H5T_NATIVE_DOUBLE, &Rin);
  readattr(file_id, "Rout",   H5T_NATIVE_DOUBLE, &Rout);
  readattr(file_id, "hslope", H5T_NATIVE_DOUBLE, &hslope);
  readattr(file_id, "R0",     H5T_NATIVE_DOUBLE, &R0);

  /* read density scale unit RHO_unit from HAMR dataset */
  #if (read_dscale==1)
  readattr(file_id, "dscale", H5T_NATIVE_DOUBLE, &RHO_unit);
  if (dscale==0){
    fprintf(stderr,"density scale = 0; should be checked!!! \n");
    exit(-3);
  }
  #endif

  /* check the parameters */
  fprintf(stderr,"t: %g, N1: %d, N2: %d, N3: %d \n",t,N1,N2,N3);
  fprintf(stderr,"startx: %g %g %g \n", startx[1],startx[2],startx[3]);
  fprintf(stderr,"dx: %g %g %g \n", dx[1],dx[2],dx[3]);
  fprintf(stderr,"a: %g, gam: %g, Rin: %g, Rout: %g, hslope: %g, R0: %g \n",a,gam,Rin,Rout,hslope,R0);

  /* nominal non-zero values for axisymmetric simulations */
  startx[0] = 0.;
  startx[2] = 0.;
  startx[3] = 0.;

  dx[2]=dx[2]/2.0;
  stopx[0] = 1.;
  stopx[1] = startx[1] + N1 * dx[1];
  stopx[2] = startx[2] + N2 * dx[2];
  stopx[3] = startx[3] + N3 * dx[3];
  //stopx[3] = 2. * M_PI;

  fprintf(stderr, "Sim range x1, x2, x3:  %g %g, %g %g, %g %g\n", startx[1],
      stopx[1], startx[2], stopx[2], startx[3], stopx[3]);

  fprintf(stderr, "hslope: %f, Rin: %f, Rout:%f \n", hslope, Rin, Rout);

  dx[0] = 1.;
  //dx[3] = 2. * M_PI;

  // conditional reads
  game = gam;   // check gamma whether it has different values for electrons and ions
  gamp = gam;
  //game = 4./3;
  //gamp = 5./3;

  if (with_electrons == 0){
	Thetae_unit = MP/ME;
    fprintf(stderr, "using fixed Tp/Te ratio = %g\n", tp_over_te);
    //Thetae_unit = MP/ME * (gam-1.) / (1. + tp_over_te);
    //Thetae_unit = 2./3. * MP/ME / (2. + tp_over_te);
    Thetae_unit = MP/ME * (gam-1.) / tp_over_te;
  } else if (with_electrons == 2){
    Thetae_unit = 2./3. * MP/ME / 5.;
    with_electrons = 2;
    fprintf(stderr, "using beta model with Rlow = %g and Rhigh = %g\n", trat_small, trat_large);
  } else {
    fprintf(stderr, "! please change electron model. Check with_electrons in model.c  \n");
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
    if (! params->loaded) {    // perferrable condition for HAMR dataset
      //report_bad_input(argc);
      #if(read_dscale != 1)
	  if (argc < 4) report_bad_input(argc);
	  sscanf(argv[3], "%lf", &RHO_unit);
      #else
	  if (argc < 3) report_bad_input(argc);
	  #endif

	  M_unit = RHO_unit * pow(L_unit, 3);
	  MBH = MBH_in;
      params->MBH = MBH;
      TP_OVER_TE = params->TP_OVER_TE;
    } else {
      M_unit = params->M_unit;
      MBH = params->MBH;
      TP_OVER_TE = params->TP_OVER_TE;
    }
    MBH *= MSUN;
    L_unit = GNEWT*MBH/(CL*CL);
    T_unit = L_unit/CL;
  }
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
  double ***malloc_rank3_double(int n1, int n2, int n3);
  double ****malloc_rank4_double(int n1, int n2, int n3, int n4);
  p = malloc_rank4_double(NVAR, N1, N2, N3);
  sigma_array = malloc_rank3_double(N1, N2, N3);
  beta_array = malloc_rank3_double(N1, N2, N3);
  fprintf(stderr, "NVAR N1 N2 N3 = %i %i %i %i\n", NVAR, N1, N2, N3);
  n2gens = (double ***)malloc_rank3(N1, N2, N3, sizeof(double));
  geom = (struct of_geom**)malloc_rank2(N1, N2, sizeof(struct of_geom));
  tetrads = (struct of_tetrads***)malloc_rank3(N1, N2, N3, sizeof(struct of_tetrads));
  init_geometry();

  // Read prims.
  /* allocate the memory for dataset */
  x1_in    = (double *) malloc(N1*N2*N3 * sizeof(double));
  x2_in    = (double *) malloc(N1*N2*N3 * sizeof(double));
  x3_in    = (double *) malloc(N1*N2*N3 * sizeof(double));
  r_in     = (double *) malloc(N1*N2*N3 * sizeof(double));
  h_in     = (double *) malloc(N1*N2*N3 * sizeof(double));
  ph_in    = (double *) malloc(N1*N2*N3 * sizeof(double));
  RHO_in   = (double *) malloc(N1*N2*N3 * sizeof(double));
  UU_in    = (double *) malloc(N1*N2*N3 * sizeof(double));
  //U0_in    = (double *) malloc(N1*N2*N3 * sizeof(double)); 
  U1_in    = (double *) malloc(N1*N2*N3 * sizeof(double));
  U2_in    = (double *) malloc(N1*N2*N3 * sizeof(double));
  U3_in    = (double *) malloc(N1*N2*N3 * sizeof(double));
  B1_in    = (double *) malloc(N1*N2*N3 * sizeof(double));
  B2_in    = (double *) malloc(N1*N2*N3 * sizeof(double));
  B3_in    = (double *) malloc(N1*N2*N3 * sizeof(double));
  gdet_in  = (double *) malloc(N1*N2    * sizeof(double));
  Ucov0_in = (double *) malloc(N1*N2*N3 * sizeof(double));
  Ucon0_in = (double *) malloc(N1*N2*N3 * sizeof(double));

  /* memory size of the data */
  dimsm[0] = N1*N2*N3;
  memspace = H5Screate_simple(RANK_OUT,dimsm,NULL);

  /* read the datasets */
  readdata(file_id, "x1",    H5T_NATIVE_DOUBLE, memspace, &x1_in[0]);
  readdata(file_id, "x2",    H5T_NATIVE_DOUBLE, memspace, &x2_in[0]);
  readdata(file_id, "x3",    H5T_NATIVE_DOUBLE, memspace, &x3_in[0]);
  readdata(file_id, "r",     H5T_NATIVE_DOUBLE, memspace, &r_in[0]);
  readdata(file_id, "h",     H5T_NATIVE_DOUBLE, memspace, &h_in[0]);
  readdata(file_id, "ph",    H5T_NATIVE_DOUBLE, memspace, &ph_in[0]);
  readdata(file_id, "RHO",   H5T_NATIVE_DOUBLE, memspace, &RHO_in[0]);
  readdata(file_id, "UU",    H5T_NATIVE_DOUBLE, memspace, &UU_in[0]);
  //readdata(file_id, "U0",    H5T_NATIVE_DOUBLE, memspace, &U0_in[0]);
  readdata(file_id, "U1",    H5T_NATIVE_DOUBLE, memspace, &U1_in[0]);
  readdata(file_id, "U2",    H5T_NATIVE_DOUBLE, memspace, &U2_in[0]);
  readdata(file_id, "U3",    H5T_NATIVE_DOUBLE, memspace, &U3_in[0]);
  readdata(file_id, "B1",    H5T_NATIVE_DOUBLE, memspace, &B1_in[0]);
  readdata(file_id, "B2",    H5T_NATIVE_DOUBLE, memspace, &B2_in[0]);
  readdata(file_id, "B3",    H5T_NATIVE_DOUBLE, memspace, &B3_in[0]);
  readdata(file_id, "Ucov0", H5T_NATIVE_DOUBLE, memspace, &Ucov0_in[0]);
  readdata(file_id, "Ucon0", H5T_NATIVE_DOUBLE, memspace, &Ucon0_in[0]);

  /* the memory space of "gdet" is 2D */
  dimsm[0] = N1*N2;
  memspace = H5Screate_simple(RANK_OUT,dimsm,NULL);

  readdata(file_id, "gdet", H5T_NATIVE_DOUBLE, memspace, &gdet_in[0]);

  /* close HDF5 file */
  ret = H5Fclose(file_id);

  /* find the index for event horizon ridius */
  for (i=0;i<N1;i++){
    if (r_in[i*N2*N3] >= Rh){
      ieh = i;
      break;
    }
  }
  fprintf(stderr, "the radius of event horizon is %g and the index of X1 is %i \n", Rh, ieh);

  /* pass the 1D dataset to pointers */
  for (i=0;i<N1;i++) for(j=0;j<N2;j++) for(z=0;z<N3;z++){
      gridIndex   = i*N2*N3 + j*N3 + z;
      gridIndex2D = i*N2 + j;

      x[1] = x1_in[gridIndex];
      x[2] = x2_in[gridIndex];
      x[3] = x3_in[gridIndex];
      rin  = r_in[gridIndex];
      hin  = h_in[gridIndex];
      phin = ph_in[gridIndex];

      /* 
        H-AMR internal coordinates: x2c = (1+x2)/2 
        --> In grmonty, x2 is treated as x2c
      */
      x2temp=(1.0+x[2])/2.0;
      x[2]=x2temp;
      /* check that we've got the coordinate parameters right */
      coord(i,j,z,xp);
      bl_coord(x, &rp, &hp);
      if (fabs(x[1]-xp[1]) > 1.e5 * x[1] || fabs(x[2]-xp[2]) > 1.e5 || fabs(x[3]-xp[3]) > 1.e5){
          fprintf(stderr, "grid setup error\n");
          fprintf(stderr, "x[1],xp[1],x[2],xp[2],x[3],xp[3]: %g %g %g %g %g %g \n",
              x[1], xp[1], x[2], xp[2], x[3], xp[3]);
          fprintf(stderr,
              "check the internal coordinates, and continue\n");
          exit(1);
      } else if (fabs(rp - rin) > 1.e-3 * rp || fabs(hp - hin) > 1.e-5 || fabs(xp[3]-phin) > 1.e-5) {
          fprintf(stderr, "grid setup error\n");
          fprintf(stderr, "rp,r, (rp-r)/r, hp,h,php,ph: %g %g %g %g %g %g %g \n",
              rp, rin, (rp-rin)/rp, hp, hin, xp[3], phin);
          fprintf(stderr,
              "edit R0, hslope, compile, and continue\n");
          exit(1);
      }

      /* Since x2c = (1+x2)/2, the vectors in x2 direction should be corrected.
         Or, we need to correct the theta correction term in gcov_func  (pi -> pi/2)
      */
      p[KRHO][i][j][z] = RHO_in[gridIndex];
      p[UU][i][j][z]   = UU_in[gridIndex];
      //p[U0][i][j][z]   = U0_in[gridIndex];
      p[U1][i][j][z]   = U1_in[gridIndex];
      p[U2][i][j][z]   = U2_in[gridIndex]/2.;
      p[U3][i][j][z]   = U3_in[gridIndex];
      p[B1][i][j][z]   = B1_in[gridIndex];
      p[B2][i][j][z]   = B2_in[gridIndex]/2.;
      p[B3][i][j][z]   = B3_in[gridIndex];

	  /*
      if(i==20 && j==80 && z==20){
            fprintf(stderr, "B1= %g\n", p[B1][i][j][z]);
            fprintf(stderr, "B2= %g\n", p[B2][i][j][z]);
            fprintf(stderr, "B3= %g\n", p[B3][i][j][z]);
      }  
	  */

      gdet  = gdet_in[gridIndex2D];
      Ucov0 = Ucov0_in[gridIndex];
      Ucon0 = Ucon0_in[gridIndex];

      //bias_norm +=
      //   dV * gdet * pow(p[UU][i][j][z] / p[KRHO][i][j][z] *
      //         Thetae_unit, 2.);
      V += dV * gdet;

      /* check accretion rate */
      if (i == ieh)
          dMact += gdet * p[KRHO][i][j][z] * p[U1][i][j][z] *Ucon0;
      if (i >= 20 && i < 40)
          Ladv += gdet * p[UU][i][j][z] * p[U1][i][j][z] *Ucon0 * Ucov0;

      /* check the dataset */
      /*printf("%d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g \n",
              i,j,z,x[1],x[2],x[3],
              rr,h,ph,
              p[KRHO][i][j][z],p[UU][i][j][z],p[U1][i][j][z],p[U2][i][j][z],p[U3][i][j][z],
              p[B1][i][j][z],p[B2][i][j][z],p[B3][i][j][z],gdet,
              Ucov0,Ucon0);*/

      /*printf("%d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g \n",
              i,j,z,x1_in[gridIndex],x2_in[gridIndex],x3_in[gridIndex],
              r_in[gridIndex],h_in[gridIndex],ph_in[gridIndex],
              RHO_in[gridIndex],UU_in[gridIndex],U1_in[gridIndex],U2_in[gridIndex],U3_in[gridIndex],
              B1_in[gridIndex],B2_in[gridIndex],B3_in[gridIndex],gdet_in[gridIndex2D],
              Ucov0_in[gridIndex],Ucon0_in[gridIndex]);*/

	  double Ne, Thetae, Bmag, Ucon[NDIM], Ucov[NDIM], Bcon[NDIM];
      get_fluid_zone(i, j, z, &Ne, &Thetae, &Bmag, Ucon, Bcon);

      double bsq = Bmag*Bmag/B_unit/B_unit;  // in code units
      sigma_array[i][j][z] = bsq / (Ne / Ne_unit);
      beta_array[i][j][z] = p[UU][i][j][z] * (gam - 1.) * 2. / bsq;

  }

  /* deallocate memories */
  free(x1_in);
  free(x2_in);
  free(x3_in);
  free(r_in);
  free(h_in);
  free(ph_in);
  free(RHO_in);
  free(UU_in);
  //free(U0_in);
  free(U1_in);
  free(U2_in);
  free(U3_in);
  free(B1_in);
  free(B2_in);
  free(B3_in);
  free(gdet_in);
  free(Ucov0_in);
  free(Ucon0_in);

  //bias_norm /= V;
  //dMact *= dx[3] * dx[2];
  /* since dx[2] was rearranged by dx[2] = dx[2]/2 while using gdet from the data, 
    the accretion rate should be multiplied by 2 */
  dMact *= dx[3] * dx[2] * 2;
  //dMact /= 21.;
  Ladv *= dx[3] * dx[2];
  Ladv /= 21.;
  fprintf(stderr, "dMact: %g, Ladv: %g\n", dMact, Ladv);

  /* done! */

#if (MODEL_EDF==EDF_MAXWELL_JUTTNER)
  fprintf(stderr, "Model EDF: Thermal distribution (Maxwell Juttner). \n");
#elif (MODEL_EDF==EDF_KAPPA_FIXED)
  fprintf(stderr, "Model EDF: Kappa fixed.\n");
#elif (MODEL_EDF==EDF_POWER_LAW)
  fprintf(stderr, "Model EDF: POWER law.\n");
#endif

  init_tetrads();
}

void readattr(hid_t file_id, const char *attr_name, hid_t mem_type_id, void *buf)
{
    hid_t attr_id;        /* attribute identifier */
    herr_t ret;           /* Return value */

    attr_id = H5Aopen(file_id, attr_name, H5P_DEFAULT);
    ret     = H5Aread(attr_id, mem_type_id, buf);
    ret     = H5Aclose(attr_id);
}


void readdata(hid_t file_id, const char *attr_name, hid_t mem_type_id, hid_t memspace, void *buf)
{
    hid_t ds_id;          /* dataset identifier */
    herr_t ret;           /* Return value */
    hid_t dataspace;      /* data space identifier */

    ds_id     = H5Dopen(file_id, attr_name, H5P_DEFAULT);
    dataspace = H5Dget_space(ds_id);    /* dataspace handle */
    ret       = H5Dread(ds_id, mem_type_id, memspace, dataspace, H5P_DEFAULT, buf);
    ret       = H5Dclose(ds_id);
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

  h5io_add_group(fid, "/params");

  #if CUSTOM_AVG==1
  h5io_add_data_dbl(fid, "/params/CA_MIN_FREQ", CA_MIN_FREQ);
  h5io_add_data_dbl(fid, "/params/CA_MAX_FREQ", CA_MAX_FREQ);
  h5io_add_data_dbl(fid, "/params/CA_MIN_QTY", CA_MIN_QTY);
  h5io_add_data_dbl(fid, "/params/CA_MAX_QTY", CA_MAX_QTY);
  h5io_add_data_dbl(fid, "/params/CA_NBINS", CA_NBINS);
  #endif // CUSTOM_AVG

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
  h5io_add_data_dbl(fid, "/params/KAPPA", model_kappa);
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
  h5io_add_data_dbl(fid, "/params/bias", biasTuning);

  h5io_add_data_int(fid, "/params/SYNCHROTRON", SYNCHROTRON);
  h5io_add_data_int(fid, "/params/BREMSSTRAHLUNG", BREMSSTRAHLUNG);
  h5io_add_data_int(fid, "/params/COMPTON", COMPTON);
  h5io_add_data_int(fid, "/params/DIST_KAPPA", MODEL_EDF==EDF_KAPPA_FIXED?1:0);
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
