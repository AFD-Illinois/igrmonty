#include "decs.h"
#include "coordinates.h"
#include "model_radiation.h"


void report_bad_input(int argc)
{
  // this function is called very early on and can be 
  // used to set units / validate input.
  
  if (argc < 2) {
    fprintf(stderr, "usage: \n");
    fprintf(stderr, "  riaf:    Ns\n");
    exit(0);
  }
}

///////////////////////////////// SUPERPHOTONS /////////////////////////////////

#define ROULETTE  1.e4
int stop_criterion(struct of_photon *ph)
{
  // this function is called on each step and used to 
  // determine if we should stop evolving the photon.
  // generally triggered to prune the low-weight ones
  // and whenever the photon has left the domain

  // stop if weight below minimum weight
  if (ph->w < WEIGHT_MIN) {
    if (monty_rand() <= 1. / ROULETTE) {
      ph->w *= ROULETTE;
    } else {
      ph->w = 0.;
      return 1;
    }
  }

  double r, h;
  bl_coord(ph->X, &r, &h);

  // if ( out_of_bounds ) return 1;

  return 0;
}

int record_criterion(struct of_photon *ph)
{
  // when we stop evolving a superphoton, this function determines
  // whether or not we record it in the spectrum. in general, only
  // record photons that are "far enough" away so that any effects
  // due to geometry are small.

  double r, h;
  bl_coord(ph->X, &r, &h);

  // if ( should_record_to_spectrum ) return 1;
  
  return 0;
}
#undef ROULETTE

#define EPS 0.04
double stepsize(double X[NDIM], double K[NDIM])
{
  // determine the stepsize for each geodesic step. note especially
  // the way we set dlx2, which cares about the "upper" X2 limit.
  

  double dl, dlx1, dlx2, dlx3;
  double idlx1, idlx2, idlx3;

  #define MIN(A,B) (A<B?A:B)
  
  #define X2MAX 1

  dlx1 = EPS / (fabs(K[1]) + SMALL);
  dlx2 = EPS * MIN(X[2], X2MAX - X[2]) / (fabs(K[2]) + SMALL);
  dlx3 = EPS / (fabs(K[3]) + SMALL) ;

  #undef X2MAX

  #undef MIN

  idlx1 = 1. / (fabs(dlx1) + SMALL);
  idlx2 = 1. / (fabs(dlx2) + SMALL);
  idlx3 = 1. / (fabs(dlx3) + SMALL);

  dl = 1. / (idlx1 + idlx2 + idlx3);

  return dl;
}
#undef EPS

void record_super_photon(struct of_photon *ph)
{
  // record the superphoton. you probably don't need to modify this function
  // unless you want to change the binning scheme, which is currently set to
  // record in elevational bins with constant dtheta (not X2)

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
  // TODO: should shared_spect be explicitly set to zero (in all model files?)
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
  // bias tuning magic. you'll may need to play around with this to get
  // a nice spectrum for the Compton component(s)

  double bias, max;

  max = 0.5 * w / WEIGHT_MIN;

  bias = Te * Te / (5. * max_tau_scatt);

  if (bias > max)
    bias = max;

  return  bias * biasTuning;
}

void get_fluid_zone(int i, int j, int k, double *Ne, double *Thetae, double *B,
        double Ucon[NDIM], double Bcon[NDIM])
{
  // return Ne, Thetae, B (in cgs) as well as Ucon and Bcon at center of i,j,k-th 
  // zone. B must be equal to sqrt(Bcon.Bcov) * B_unit. vector components must be
  // in "raytracing" coordinates.

  // get grid zone center
  double X[NDIM] = { 0. };
  ijktoX(i, j, k, X);

  // get metric at grid center
  double gcov[NDIM][NDIM];
  gcov_func(X, gcov);

  // create dummy vectors
  double Ucov[NDIM] = { 0. };
  double Bcov[NDIM] = { 0. };

  // call into analytic model
  get_fluid_params(X, gcov, Ne, Thetae, B, Ucon, Ucov, Bcon, Bcov);
}

void get_fluid_params(const double X[NDIM], double gcov[NDIM][NDIM], double *Ne,
          double *Thetae, double *B, double Ucon[NDIM],
          double Ucov[NDIM], double Bcon[NDIM],
          double Bcov[NDIM])
{
  // return the values of Ne, Thetae, B (in cgs) as well as Ucon,Ucov,Bcon,
  // and Bcov at X (with covariant metric gcov) where the vector components
  // are in "raytracing" coordinates. note that B must be equal to 
  //    sqrt(Bcon.Bcov) * B_unit

  *Ne = 0.; 
  *Thetae = 0.;
  *B = 0.;

  double gcon[NDIM][NDIM];
  gcon_func(gcov, gcon);
  Ucov[0] = sqrt(-1 / gcon[0][0]);
  Ucov[1] = 0;
  Ucov[2] = 0;
  Ucov[3] = 0;

  Bcon[0] = 0;
  Bcon[1] = 0;
  Bcon[2] = 0;
  Bcon[3] = 0;

  lower(Ucov, gcon, Ucon);
  lower(Bcon, gcov, Bcov);

  // BOGUE RETURN
}

////////////////////////////////// COORDINATES /////////////////////////////////

void gcov_func(const double X[NDIM], double gcov[NDIM][NDIM])
{
  // set the covariant metric gcov based on location X. all components
  // are in "raytracing" coordinates. note that you can use any of the
  // provided functions like bl_coord, gcov_ks, set_dxdX

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

double dOmega_func(int j)
{
  // return solid angle subtended by j-th elevational record bin

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
  // set any general model parameters 
  biasTuning = params->biasTuning;
  model_kappa = 4.;

  // set metric information
  with_derefine_poles = 0; // if we don't set anything, code defaults to MKS
  hslope = 1.;

  // set units. note that B_unit is used elsewhere!
  double MBH_solar = 4.14e6;
  L_unit = GNEWT * MBH_solar * MSUN / (CL * CL);
  RHO_unit = 1. * (MP + ME);
  B_unit = CL * sqrt(4.*M_PI*RHO_unit); 
  T_unit = L_unit / CL;
  M_unit = RHO_unit * pow(L_unit, 3);
  max_tau_scatt = (6. * L_unit) * RHO_unit * 0.4;

  // set up grid. grmonty requires some notion of a grid in order to figure
  // out where to emit the original superphotons
  startx[0] = 0.0;
  startx[1] = 0.1;
  startx[2] = 0.0;
  startx[3] = 0.0;

  // make sure to think about objects with high optical depth when setting
  // radial zone extents!
  N1 = 512;
  N2 = 512;
  N3 = 1;

  dx[0] = 0.;
  dx[1] = 1. / N1;
  dx[2] = 1. / N2;
  dx[3] = 2. * M_PI / N3;

  stopx[0] = 1.;
  stopx[1] = startx[1]+N1*dx[1];
  stopx[2] = startx[2]+N2*dx[2];
  stopx[3] = startx[3]+N3*dx[3];

  fprintf(stderr, "using model TEMPLATE.\n");

  // call init_geometry to set up gcon/gcov/gdet/gdet_zone in each zone
  geom = (struct of_geom**)malloc_rank2(N1, N2, sizeof(struct of_geom));
  init_geometry();

  // call init_tetrads to set up Econ/Ecov in each zone. must be called after
  // we have a way to determine Ucon, Bcon in each zone
  tetrads = (struct of_tetrads***)malloc_rank3(N1, N2, N3, sizeof(struct of_tetrads));
  init_tetrads();

  // allocate memory for "number of superphotons to generate per zone"
  n2gens = (double ***)malloc_rank3(N1, N2, N3, sizeof(double));
}

//////////////////////////////////// OUTPUT ////////////////////////////////////

void report_spectrum(int N_superph_made, Params *params)
{
  // write spectrum file and print information to the screen. you'll probably
  // only want to modify this function to output model parameters.

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
  h5io_add_data_dbl(fid, "/params/WEIGHT_MIN", WEIGHT_MIN);
  h5io_add_data_dbl(fid, "/params/KAPPA", model_kappa);
  h5io_add_data_dbl(fid, "/params/L_unit", L_unit);
  h5io_add_data_dbl(fid, "/params/T_unit", T_unit);
  h5io_add_data_dbl(fid, "/params/Thetae_unit", Thetae_unit);
  h5io_add_data_dbl(fid, "/params/Rin", Rin);
  h5io_add_data_dbl(fid, "/params/Rout", Rmax);
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

  h5io_add_data_str(fid, "/params/model", xstr(MODEL));

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
  double nph_buf[N_TYPEBINS][N_EBINS][N_THBINS];

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

        nph_buf[k][i][j] = spect[k][j][i].nph;

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
  h5io_add_data_dbl_3d(fid, "/output/nph", N_TYPEBINS, N_EBINS, N_THBINS, nph_buf);
  h5io_add_data_dbl_1d(fid, "/output/Lcomponent", N_TYPEBINS, Lcomponent_buf);

  h5io_add_data_int(fid, "/output/Nrecorded", N_superph_recorded);
  h5io_add_data_int(fid, "/output/Nmade", N_superph_made);
  h5io_add_data_int(fid, "/output/Nscattered", N_scatt);

  double Lum = L * LSUN;
  h5io_add_data_dbl(fid, "/output/L", Lum);
  h5io_add_attribute_str(fid, "/output/L", "units", "erg/s");

  // diagnostic output to screen
  fprintf(stderr, "\n");

  fprintf(stderr, "max_tau_scatt = %g\n", max_tau_scatt);
  fprintf(stderr, "L = %g erg/s \n", Lum);

  fprintf(stderr, "\n");

  fprintf(stderr, "N_superph_made = %d\n", N_superph_made);
  fprintf(stderr, "N_superph_scatt = %d\n", N_scatt);
  fprintf(stderr, "N_superph_recorded = %d\n", N_superph_recorded);

  H5Fclose(fid);

}


