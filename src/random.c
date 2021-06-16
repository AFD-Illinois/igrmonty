#include "decs.h"

gsl_rng *r;
#pragma omp threadprivate(r)

void init_monty_rand(int seed)
{
  if(seed < 0)
    seed = time(NULL);
  //gsl_rng = (gsl_rng*)malloc_rank1(nthreads, sizeof(gsl_rng*));
  #pragma omp parallel
  {
    r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, 139 * omp_get_thread_num() + seed);
  }
}

double monty_rand()
{
  return gsl_rng_uniform(r);
}

void monty_ran_dir_3d(double *n0x, double *n0y, double *n0z)
{
  gsl_ran_dir_3d(r, n0x, n0y, n0z);
}

double monty_ran_chisq(int n)
{
  return gsl_ran_chisq(r, n);
}
