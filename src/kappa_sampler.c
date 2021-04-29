//
//  main.c
//  kappasampler
//
//  Created by Jordy on 13/07/2017.
//  Copyright © 2017 Jordy. All rights reserved.
//

//  GNW included this file from jordydavelaar/kmonty for comparison purposes

#define gmin 25
#define gamma_max 1e3
#define gmax 1e7
#define pindex 3

#include "decs.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
// rewrite this

static double yhigh;


struct f_params {
  double u;
};



double hypergeom_eval(double X) {
  double kappa=KAPPA;
  double z,hyp2F1;
  double a = 1.;
  double b = -kappa - 1. / 2.;
  double c = 1. / 2.;

     if(fabs(X)<1e-6){
  z=-X;
        hyp2F1=1+a*b*z/c + a*(1+a)*b*(1+b)/(2*c*(1+c))*z*z;
      return hyp2F1;
     }
     else if(fabs(X)>1) {
          z = -X;
          hyp2F1 = pow(1.-z, -a) * gsl_sf_gamma(c) * gsl_sf_gamma(b-a)
                                     / (gsl_sf_gamma(b)*gsl_sf_gamma(c-a))
         * gsl_sf_hyperg_2F1(a, c-b, a-b+1., 1./(1.-z))
         + pow(1.-z, -b) * gsl_sf_gamma(c) * gsl_sf_gamma(a-b)
                                     / (gsl_sf_gamma(a) * gsl_sf_gamma(c-b))
         * gsl_sf_hyperg_2F1(b, c-a, b-a+1., 1./(1.-z));
                     return hyp2F1;
             }
     else if(fabs(X)<1) {
                     //fprintf(stderr,"too small %e %e\n", X,gsl_sf_hyperg_2F1(a,b,c,-X));
                     return gsl_sf_hyperg_2F1(a,b,c,-X);
             }

         /*
  di = 1000. * log(X / HYPMIN) / log(10);
  i = (int)di;
  double dX = (X - HYPMIN * pow(10, i / 1000.)) /
              (HYPMIN * pow(10, (i + 1) / 1000.) - HYPMIN * pow(10, i / 1000.));
//  if (kappa_synch >= kappa_max) {
//    printf("Kappa=%0.1lf not supported, maximum kappa is %0.1lf\n", kappa_synch,
//           kappa_max - dkappa);
//    exit(1);
//  }
  if (kappa_synch < kappa_min) {
    printf("Kappa=%0.1lf not supported, minimum kappa is %0.1lf\n", kappa_synch,
           kappa_min);
    exit(1);
  }
  int k = (int)((kappa_synch - kappa_min) / dkappa);
  double value = (1 - dX) * hypergeom[k][i] + dX * hypergeom[k][i + 1];

  return value;
*/
  return 0;
}


double find_y(double u, double (*df)(double), double (*f)(double, void *),
              double w) {
  // Get maximum for window
  int status, steps = 0, max_steps = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double solution;
  double low = 1e-5; // sqrt(0.001/w);
  double high = 1e5; // sqrt(1e4/w);
  gsl_function F;
  struct f_params params = {u};
  F.function = f;
  F.params = &params;
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(s, &F, low, high);
  // printf("trying to find %e\n",u);
  do {
    steps++;
    status = gsl_root_fsolver_iterate(s);
    solution = gsl_root_fsolver_root(s);
    low = gsl_root_fsolver_x_lower(s);
    high = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(low, high, 1e-3, 0);
  } while (status == GSL_CONTINUE && steps < max_steps);
  // fprintf(stderr,"steps %d solution %e\n",steps,solution);
  return solution;
}

double dF3(double y) {
  double value, denom, num;
  double y2 = y * y;
  double kappa = KAPPA;

  num =
      4. * y2 * pow((kappa + y2) / (kappa), -kappa - 1.) * gsl_sf_gamma(kappa);
  denom = sqrt(M_PI) * sqrt(kappa) * gsl_sf_gamma(kappa - 1. / 2.);
  value = num / denom;

  return value;
}
double dF4(double y) {
  double value, denom, num;
  double y2 = y * y;
  double kappa = KAPPA;

  num = 2. * (kappa - 1) * y2 * y * pow((kappa + y2) / (kappa), -kappa - 1.);
  denom = kappa;
  value = num / denom;

  return value;
}
double dF5(double y) {
  double value, denom, num;
  double y2 = y * y;
  double kappa = KAPPA;

  num = 8 * y2 * y2 * pow((kappa + y2) / (kappa), -kappa - 1) *
        gsl_sf_gamma(kappa);
  denom = 3 * sqrt(M_PI) * pow(kappa, 3. / 2.) * gsl_sf_gamma(kappa - 3. / 2.);
  value = num / denom;

  return value;
}
double dF6(double y) {
  double value, denom, num;
  double y2 = y * y;
  double kappa = KAPPA;

  num = (kappa * kappa - 3. * kappa + 2) * pow(y, 5.) *
        pow((kappa + y2) / (kappa), -kappa - 1);
  denom = kappa * kappa;
  value = num / denom;
  return value;
}

double F3(double y, void *params) {
  struct f_params *p = (struct f_params *)params;
  double u = p->u;
  double value, denom, num, hyp2F1;
  double kappa = KAPPA;

  double y2 = y * y;
  double z = -y2 / kappa;

  hyp2F1 = hypergeom_eval(-z);

  num = -sqrt(kappa) * pow(((y2 + kappa) / kappa), -kappa) *
        gsl_sf_gamma(kappa) * (-kappa * hyp2F1 + y2 * (2 * kappa + 1) + kappa);
  denom = y * sqrt(M_PI) * gsl_sf_gamma(3. / 2. + kappa);

  value = num / denom - u;

  return value;
}
double F4(double y, void *params) {
  struct f_params *p = (struct f_params *)params;
  double u = p->u;
  double value, denom, num;
  double y2 = y * y;
  double kappa = KAPPA;

  num = 1 - (1 + y2) * pow((kappa + y2) / (kappa), -kappa);
  denom = 1;

  value = num / denom - u;

  return value;
}
double F5(double y, void *params) {
  struct f_params *p = (struct f_params *)params;
  double u = p->u;

  double value, denom, num, hyp2F1;
  double kappa = KAPPA;

  double y2 = y * y;
  double z = -y2 / kappa;

  hyp2F1 = hypergeom_eval(-z);

  num =
      pow((y2 + kappa) / kappa, -kappa) * gsl_sf_gamma(kappa) *
      (3 * kappa * kappa * (hyp2F1 - 1) + (1. - 4. * kappa * kappa) * y2 * y2 -
       3. * kappa * (2. * kappa + 1.) * y2);
  denom =
      3. * pow(kappa, 1. / 2.) * y * sqrt(M_PI) * gsl_sf_gamma(3. / 2. + kappa);
  value = num / denom - u;

  return value;
}
double F6(double y, void *params) {
  struct f_params *p = (struct f_params *)params;
  double u = p->u;

  double value, denom, num;
  double y2 = y * y;
  double y4 = y2 * y2;
  double kappa = KAPPA;

  num =
      (y4 - (y4 + 2. * y2 + 2.) * kappa) * pow((kappa + y2) / (kappa), -kappa);
  denom = 2 * kappa;
  value = num / denom + 1 - u;

  return value;
}

double sample_y_distr_kappa(double Thetae, double kappa) {

  // this code is untested and seems to disagree with 
  // analytic formula for kappa eDF pdf
  assert(1==0);

  double w = (kappa - 3.)/kappa * Thetae;
  double S_3, pi_3, pi_4, pi_5, pi_6, y, x1, x2, prob;
  double num, den;
  yhigh = sqrt((10 * gamma_max - 1) / w);
  // gsl_set_error_handler_off();

  pi_3 = sqrt(kappa) * sqrt(M_PI) * gsl_sf_gamma(-1. / 2. + kappa) /
         (4. * gsl_sf_gamma(kappa));
  pi_4 = kappa / (2. * kappa - 2.) * sqrt(0.5 * w);
  pi_5 = 3. * pow(kappa, 3. / 2.) * sqrt(M_PI) *
         gsl_sf_gamma(-3. / 2. + kappa) / (8. * gsl_sf_gamma(kappa)) * w;
  pi_6 = kappa * kappa / (2. - 3. * kappa + kappa * kappa) * w *
         sqrt(0.5 * w);

  S_3 = pi_3 + pi_4 + pi_5 + pi_6;

  pi_3 /= S_3;
  pi_4 /= S_3;
  pi_5 /= S_3;
  pi_6 /= S_3;
  //    printf("%e\n",pi_3);
  //    printf("%e\n",pi_3+pi_4);
  //    printf("%e\n",pi_3+pi_4+pi_5);
  //    printf("%e\n",pi_3+pi_4+pi_5+pi_6);

  do {
    x1 = monty_rand();
    double u = monty_rand();
    if (x1 < pi_3) {
      y = find_y(u, dF3, F3, w);
    } else if (x1 < pi_3 + pi_4) {
      y = find_y(u, dF4, F4, w);
    } else if (x1 < pi_3 + pi_4 + pi_5) {
      y = find_y(u, dF5, F5, w);
    } else {
      y = find_y(u, dF6, F6, w);
    }

    x2 = monty_rand();
    num = sqrt(1. + 0.5 * w * y * y);
    den = (1. + y * sqrt(0.5 * w));
    prob = (num / den) *
           exp(-(w * y * y) / gamma_max); //* w*sqrt(2*w)/S_3;
    if (y != y)
      prob = 0;
    //     printf("nth %e %e prob %e
    //     %e\n",(1+w*y*y)/gamma_max,(1+w*y*y),prob,x2 );

  } while (x2 >= prob);

  return (y);
}

double sample_gamma_distr_pwl(){
double p = pindex;
double x1 = (double)rand()/(double)RAND_MAX;

return  pow(pow(gmin,1-p)*(1 - x1 ) + x1 * pow(gmax,1-p),1/(1-p));
}

#undef pindex
#undef gmin
#undef gamma_max
#undef gmax
