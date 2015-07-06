

#include "decs.h"

/* 

   set up the metric (indicies both up and down), and the
   connection coefficients on a grid, based
   on a HARM simulation grid.  
   
   In principle the geometry could be sampled on a finer or
   coarser grid than the simulation grid. 

   These routines are taken directly out of HARM.

   They require the code be compiled against the 
   Gnu scientific library (GSL).

   CFG 21 July 06
   
*/

gsl_matrix *gsl_gcov, *gsl_gcon;
gsl_permutation *perm;
#pragma omp threadprivate (gsl_gcov, gsl_gcon, perm)

/* assumes gcov has been set first; returns determinant */
double gdet_func(double gcov[][NDIM])
{
	double d;
	int k, l, signum;

	if (gsl_gcov == NULL) {
		gsl_gcov = gsl_matrix_alloc(NDIM, NDIM);
		gsl_gcon = gsl_matrix_alloc(NDIM, NDIM);
		perm = gsl_permutation_alloc(NDIM);
	}

	DLOOP gsl_matrix_set(gsl_gcov, k, l, gcov[k][l]);

	gsl_linalg_LU_decomp(gsl_gcov, perm, &signum);

	d = gsl_linalg_LU_det(gsl_gcov, signum);

	return (sqrt(fabs(d)));
}

/* invert gcov to get gcon */

/*
void gcon_func(double gcov[][NDIM], double gcon[][NDIM])
{
	int k, l, signum;
	
	if (gsl_gcov  == NULL) {
		gsl_gcov = gsl_matrix_alloc(NDIM, NDIM);
		gsl_gcon = gsl_matrix_alloc(NDIM, NDIM);
		perm = gsl_permutation_alloc(NDIM);
	}
	

	DLOOP gsl_matrix_set(gsl_gcov, k, l, gcov[k][l]);

	gsl_linalg_LU_decomp(gsl_gcov, perm, &signum);

	gsl_linalg_LU_invert(gsl_gcov, perm, gsl_gcon);

	DLOOP gcon[k][l] = gsl_matrix_get(gsl_gcon, k, l);

	// done!
}
*/

#undef DELTA
