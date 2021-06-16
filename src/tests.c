#include "decs.h"
#include "hotcross.h"
#include "compton.h"
#include "model_radiation.h"

// test dNdgammae function in hotcross.c this is the only 
// public "interface" for the sampling eDF, so it can act
// as a sort of regression test
void test_hotcross_dNdgammae(const char *ofname, double Thetae)
{
  fprintf(stderr, "testing src/hotcross.c for Thetae=%g\n", Thetae);

  init_hotcross();

  FILE *fp = fopen(ofname, "w");
  fprintf(fp, "# %g %g %g\n", Thetae, model_kappa, powerlaw_p);

  double norm = getnorm_dNdg(Thetae);
  for (double lge = 0.; lge < 5; lge += 0.001) {
    fprintf(fp, "%g %g\n", lge, dNdgammae(Thetae, pow(10., lge)) * norm);
  }

  fprintf(fp, "\n");
  fclose(fp);
}

// test sample_beta_distr in compton.c. as above, this is the
// only public "interface", so we use it as a regression test
void test_compton_sample_beta_dist(const char *ofname, double Thetae)
{
  fprintf(stderr, "testing src/compton.c for Thetae=%g\n", Thetae);

  init_monty_rand(64);
  int nsamp = 100000;

  FILE *fp;
  double ge, be;

  fp = fopen(ofname, "w");
  fprintf(fp, "# %g %g %g\n", Thetae, model_kappa, powerlaw_p);

  for (int i=0; i<nsamp; ++i) {
    sample_beta_distr(Thetae, &ge, &be);
    fprintf(fp, "%g ", ge);
  }

  fprintf(fp, "\n");
  fclose(fp);
}


double Thetae_from_kappa_w(double kappa, double w)
{
  return w * kappa / (kappa - 3.);
}

void run_all_tests() {

  // note: in the future, it might make sense to allow 
  // switching of the eDF at runtime
  
  // set eDF parameters
  model_kappa = 4.;
  powerlaw_p = 3.;
  powerlaw_gamma_min = 25.;
  powerlaw_gamma_max = 1.e7;
  powerlaw_gamma_cut = 1.e3;

  // test hotcross.c functionality
  test_hotcross_dNdgammae("test/dNdgammae_0.1.out", 0.1);
  test_hotcross_dNdgammae("test/dNdgammae_1.out", 1);
  test_hotcross_dNdgammae("test/dNdgammae_5.out", 5);
  test_hotcross_dNdgammae("test/dNdgammae_10.out", 10);

  // test compton.c functionality
  test_compton_sample_beta_dist("test/sample_beta_0.1.out", 0.1);
  test_compton_sample_beta_dist("test/sample_beta_1.out", 1);
  test_compton_sample_beta_dist("test/sample_beta_5.out", 5);
  test_compton_sample_beta_dist("test/sample_beta_10.out", 10);

  exit(42);
}

// functions below left in for legacy reasons
// primarily used to unit test the individual
// components of the above. not guaranteed to
// compile/work.

void test_compton_sampling_functions(double kappa)
{
  fprintf(stderr, "testing eDF in compton.c for kappa=%g\n", kappa);

  model_kappa = kappa;

  FILE *fp = fopen("test/compton_sampling_functions.out", "w");
  fprintf(fp, "%g\n", kappa);

  // test to make sure the distribution functions return reasonable values
  for (double Thetae = 0; Thetae<100; Thetae+=0.5) {
    double geofmin = -1;
    double vofmin = 1.e10;
    double dofmin = 0;
    for (double ge=1.; ge < 1001; ge+=0.01) {
      double dist = fdist(ge, Thetae);
      double ddist = fabs(dfdgam(ge, &Thetae));
      if ( ddist < vofmin ) {
        geofmin = ge;
        vofmin = ddist;
        dofmin = dist;
      }
    }
    fprintf(fp, "%g %g %g %g\n", Thetae, geofmin, vofmin, dofmin);
  }

  fprintf(fp, "\n");
  fclose(fp);
}

void test_compton_sampling(double Thetae, double kappa, const char *fname)
{
  fprintf(stderr, "testing Compton sampling for Thetae=%g and kappa=%g\n", Thetae, kappa);

  model_kappa = kappa;

  FILE *fp = fopen(fname, "w");
  fprintf(fp, "%g %g ", kappa, Thetae);
  
  // draw monte carlo samples for some value of Thetae
  // check that the proper distribution is recovered
  init_monty_rand(64);
  double gamma_e, beta_e;
  for (int i=0; i<10000; ++i) {
    sample_beta_distr(Thetae, &gamma_e, &beta_e); 
    fprintf(fp, "%g ", gamma_e);
  }

  fprintf(fp, "\n");
  fclose(fp);
}

void test_emiss_abs()
{
  init_emiss_tables();

  model_kappa = 4.;

  fprintf(stderr, "DIST_KAPPA %d\n", MODEL_EDF==EDF_KAPPA_FIXED?1:0);

  double B = 10;
  double Thetae = 10;
  double theta = M_PI/3.;

  for (double lnu=9; lnu < 15; lnu += 0.2) {
    double nu = pow(10., lnu);
    double iem = jnu_inv(nu, Thetae, 1., B, theta);
    double iabs = alpha_inv_abs(nu, Thetae, 1., B, theta);
    if (1==0) 
    fprintf(stderr, "%g %g %g\n", nu, iem, iabs);
    double ijnu = int_jnu(1., Thetae, B, nu);
    fprintf(stderr, "%g %g\n", nu, ijnu);
  }
}


void test_hotcross() 
{
  init_hotcross();

  double w = 1.2;
  double thetae = 3.4;
  double value = total_compton_cross_lkup(w, thetae);

  fprintf(stderr, "%g %g -> %g\n", w, thetae, value);
}


