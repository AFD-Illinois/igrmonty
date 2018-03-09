
/* 

   grmonty Nph

   Using monte carlo method, estimate spectrum of an appropriately
   scaled axisymmetric GRMHD simulation as a function of 
   latitudinal viewing angle.

   Input simulation data is assumed to be in dump format provided by 
   HARM code.  Location of input file is, at present, hard coded
   (see init_sim_data.c).  

   Nph super-photons are generated in total and then allowed
   to propagate.  They are weighted according to the emissivity.
   The photons are pushed by the geodesic equation.
   Their weight decays according to the local absorption coefficient.
   The photons also scatter with probability related to the local
   scattering opacity.  

   The electrons are assumed to have a thermal distribution 
   function, and to be at the same temperature as the protons.

   CFG 8-17-06

   Implemented synchrotron sampling, 22 Jan 07

   fixed bugs in tetrad/compton scattering routines, 31 Jan 07

   Implemented normalization for output, 6 Feb 07

   Separated out different synchrotron sampling routines
   into separate files, 8 Mar 07

   fixed bug in energy recording; bug used Kcon[0] rather than 
   Kcov[0] as energy, 18 Mar 07

   major reorganization to encapsulate problem-dependent parts 5-6 Nov 07

 */

#include "decs.h"

/* defining declarations for global variables */
Params params = { 0 };
struct of_geom **geom;
int nthreads;
int N1, N2, N3, n_within_horizon;
double F[N_ESAMP + 1], wgt[N_ESAMP + 1];
int Ns, N_superph_recorded, N_scatt;
struct of_spectrum spect[N_THBINS][N_EBINS] = { };

double t;
double a;
double R0, Rin, Rh, Rout, Rms;
double hslope;
double startx[NDIM], stopx[NDIM], dx[NDIM];

double dlE, lE0;
double gam;
double dMsim;
double M_unit, L_unit, T_unit;
double RHO_unit, U_unit, B_unit, Ne_unit, Thetae_unit;
double max_tau_scatt, Ladv, dMact, bias_norm;

int main(int argc, char *argv[])
{
  //omp_set_num_threads(1);
  double N_superph_made;
  time_t currtime, starttime;

  // Spectral bin parameters
  dlE = 0.25;   // bin width
  lE0 = log(1.e-12);  // location of first bin, in electron rest-mass units

  // Load parameters
  for (int i=0; i<argc-1; ++i) {
    if ( strcmp(argv[i], "-par") == 0 ) {
      load_par(argv[i+1], &params);
    }
  }

  init_model(argc, argv, &params);

  N_superph_made = 0;
  N_superph_recorded = 0;
  N_scatt = 0;
  starttime = time(NULL);

  printf("SYNCH: %i\n", SYNCHROTRON);
  printf("BREMS: %i\n", BREMSSTRAHLUNG);
  printf("COMPT: %i\n", COMPTON);

  fprintf(stderr, "Entering main loop...\n");
 
  // Get number of superphotons for each zone
  /*void init_zone(int i, int j, int k, double *nz, double *dnmax);
  ZLOOP {
    double nz, dnmax;
    init_zone(i, j, k, &nz, &dnmax);
    if  (j == N2/2)
    printf("[%i %i %i] nz = %e dnmax = %e\n", i,j,k,nz,dnmax);
  }*/
  
  int quit_flag = 0;
  #pragma omp parallel private(spect)
  {
    struct of_photon ph;
    while (1) {

      // get pseudo-quanta 
#pragma omp critical (MAKE_SPHOT)
      {
        if (!quit_flag)
          make_super_photon(&ph, &quit_flag);
      }
      if (quit_flag)
        break;

      // push them around 
      track_super_photon(&ph);

      // step
#pragma omp atomic
      N_superph_made += 1;

      // give interim reports on rates
      if (((int) (N_superph_made)) % 100000 == 0
          && N_superph_made > 0) {
        currtime = time(NULL);
        fprintf(stderr, "time %g, rate %g ph/s\n",
          (double) (currtime - starttime),
          N_superph_made / (currtime -
                starttime));
      }
    }
  }

  currtime = time(NULL);
  fprintf(stderr, "Final time %g, rate %g ph/s\n",
    (double) (currtime - starttime),
    N_superph_made / (currtime - starttime));

  omp_reduce_spect();

  report_spectrum((int) N_superph_made, &params);

  return 0;
}

