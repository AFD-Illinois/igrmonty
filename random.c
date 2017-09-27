#include "decs.h"

gsl_rng *r;
#pragma omp threadprivate(r)

void init_monty_rand()
{
  #pragma omp parallel
  {
    int seed = 139*omp_get_thread_num() + time(NULL);
    r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, seed);
  }
}

double monty_rand()
{
  return gsl_rng_uniform(r);
}

