#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "constants.h"
#define BREMSSTRAHLUNG (3)
#include "bremss.c"

gsl_spline2d *bremss_spline;
gsl_interp_accel *bremss_xacc;
gsl_interp_accel *bremss_yacc;

int main(){

  init_bremss_spline();

// Make a grid in x = h nu / k T and tau = k T / me c^2  space, compute Gaunt factor, then test.

  gsl_vector *tau = gsl_vector_alloc(4);
  gsl_vector *x = gsl_vector_alloc(100);
  gsl_matrix *gff_ee = gsl_matrix_alloc (4, 100);

  for (int j = 0; j < 100; j++){
    gsl_vector_set (x, j, 1.e-4 + ( (10.-1.e-4)/99.) * j);
  }
  gsl_vector_set (tau, 0, 0.00195695118);
  gsl_vector_set (tau, 1, 0.00195695118 * 10);
  gsl_vector_set (tau, 2, 0.00195695118 * 1000);
  gsl_vector_set (tau, 3, 0.00195695118 * 7000);

  double te_val, nu_val;
  for (int i = 0; i < 4; i++){
    for (int j = 0; j < 100; j++){
      te_val = gsl_vector_get(tau, i) * ME * CL * CL / KBOL;
      nu_val = gsl_vector_get(x, j) * KBOL * te_val / HPL;
      gsl_matrix_set(gff_ee, i, j, gffee(te_val, nu_val));
    }
  }
  
//Verify the qualitative features of N09 Fig. 1

//Feature 1: 1 < gff < 10 at x = 1e-4 for all tau
  for (int i = 0; i < 4; i++){
    if ( ( gsl_matrix_get(gff_ee, i, 0) > 10. ) || ( gsl_matrix_get(gff_ee, i, 0) < 1. ) ) return -1;
  }

//Feature 2: 2 < gff < 2.6 at x = 10 for low tau
  for (int i = 0; i < 2; i++){
    if ( ( gsl_matrix_get(gff_ee, i, 99) < 2. ) || ( gsl_matrix_get(gff_ee, i, 99) > 2.6 ) ) return -2;
  }

//Feature 3: min(gff) \approx 1.8 for low tau
  for (int i = 0; i < 2; i++){
    gsl_vector_view row = gsl_matrix_row (gff_ee, i);
    if ( ( gsl_vector_min(&row.vector) > 2. ) || (gsl_vector_min(&row.vector) < 1.7) ) return -3;
  }

//Feature 4: 40 < gff < 50 at x = 10 for high tau
  for (int i = 2; i < 4; i++){
    if ( ( gsl_matrix_get(gff_ee, i, 99) < 40. ) || ( gsl_matrix_get(gff_ee, i, 99) > 50. ) ) return -4;
  }

//Save the matrix for plotting in Python if desired
//  FILE * fp = fopen("gff_ee.dat", "w");
//  gsl_matrix_fprintf(fp, gff_ee, "%.5e");

  gsl_vector_free (tau);
  gsl_vector_free (x);
  gsl_matrix_free (gff_ee);

  return 0;

}
