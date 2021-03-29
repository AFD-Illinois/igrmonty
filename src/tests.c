
#include "decs.h"
#include "hotcross.h"
#include "compton.h"

void test_compton_sampling_functions();
void test_compton_sampling(double Thetae, double kappa, const char *fname);
void test_emiss_abs();
void test_hotcross();
void test_dNdgamma(double Thetae, double kappa);
void test_compare_beta_dist(double Thetae, double kappa);

double Thetae_from_kappa_w(double kappa, double w)
{
  return w * kappa / (kappa - 3.);
}

void run_all_tests() {

  double kappa = 4;
  double Thetae = 1;

  kappa = 20;
  Thetae = 5;

  test_compare_beta_dist(Thetae, kappa);

  // most recent debugging quits here
  exit(42);

  test_compton_sampling(Thetae, kappa, "test/compton_sampling.out");
  test_compton_sampling_functions(kappa);

  Thetae = Thetae_from_kappa_w(kappa, 2.5);
  test_dNdgamma(Thetae, kappa);

  // now reproduce figure 4
  kappa = 5;
  Thetae = Thetae_from_kappa_w(kappa, 1.);
  test_compton_sampling(Thetae, kappa, "test/sampling_5_1.out");
  Thetae = Thetae_from_kappa_w(kappa, 5.);
  test_compton_sampling(Thetae, kappa, "test/sampling_5_5.out");
  Thetae = Thetae_from_kappa_w(kappa, 50.);
  test_compton_sampling(Thetae, kappa, "test/sampling_5_50.out");

  test_emiss_abs();
  test_hotcross();

  exit(42);
}

void test_compare_beta_dist(double Thetae, double kappa)
{
  KAPPA = kappa;

  init_monty_rand(64);
  int nsamp = 100000;

  FILE *fp;
  double ge, be;

  fp = fopen("test/beta_dist_y.out", "w");
  fprintf(fp, "%g %g\n", kappa, Thetae);

  for (int i=0; i<nsamp; ++i) {
    sample_beta_distr_y(Thetae, &ge, &be);
    fprintf(fp, "%g ", ge);
  }

  fprintf(fp, "\n");
  fclose(fp);

  fp = fopen("test/beta_dist_num.out", "w");
  fprintf(fp, "%g %g\n", kappa, Thetae);

  for (int i=0; i<nsamp; ++i) {
    sample_beta_distr_num(Thetae, &ge, &be);
    fprintf(fp, "%g ", ge);
  }

  fprintf(fp, "\n");
  fclose(fp);
}


void test_compton_sampling_functions(double kappa)
{
  fprintf(stderr, "testing eDF in compton.c for kappa=%g\n", kappa);

  KAPPA = kappa;

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

  KAPPA = kappa;

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

void test_dNdgamma(double Thetae, double kappa) 
{
  fprintf(stderr, "testing eDF for Thetae=%g and kappa=%g\n", Thetae, kappa);

  // can set kappa parameter here. make sure to initialize tables first!
  KAPPA = kappa;
  init_hotcross();

  // print normalized dNdgamma for given Thetae
  double norm = getnorm_dNdg(Thetae);

  FILE *fp = fopen("test/dNdgamma.out", "w");
  fprintf(fp, "%g %g\n", Thetae, KAPPA);

  for (double lge = 0.; lge < 3; lge += 0.01) {
    fprintf(fp, "%g %g\n", lge, dNdgammae(Thetae, pow(10., lge)) * norm);
  }

  fclose(fp);
}

void test_emiss_abs()
{
  init_emiss_tables();

  KAPPA = 4.;

  fprintf(stderr, "DIST_KAPPA %d\n", DIST_KAPPA);

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


