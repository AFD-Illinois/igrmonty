#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../../src/constants.h"
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "../../build_archive/bremss.c"

gsl_spline2d *bremss_spline;
gsl_interp_accel *bremss_xacc;
gsl_interp_accel *bremss_yacc;

int main(){

  init_bremss_spline();

// Make a grid in log_gamma_sq and log_nu space, compute Gaunt factor there, then test.
  double te_val, nu_val;

  gsl_vector *log_gamma_sq = gsl_vector_alloc(100);
  gsl_vector *log_u = gsl_vector_alloc(100);
  gsl_matrix *gff_ei = gsl_matrix_alloc (100, 100);

  for (int i = 0; i < 100; i++){
    gsl_vector_set (log_gamma_sq, i, -6. + (16./99.) * i);
    gsl_vector_set (log_u, i, -16. + (26./99.) * i);
  }

  for (int i = 0; i < 100; i++){
    for (int j = 0; j < 100; j++){
      te_val = 157900. / pow(10.,gsl_vector_get(log_gamma_sq, i));
      nu_val = pow(10., gsl_vector_get(log_u, j)) * KBOL * te_val / HPL;
      gsl_matrix_set(gff_ei, i, j, gffei(te_val, nu_val));
    }
  }
  
//Verify the qualitative features of VH15 Fig. 5 left panel (constant loggammasq, variable logu)

//Feature 1: 1 < log_10(gff) < 3 for log_10(u) = -16 at all temperatures.
  for (int i = 0; i < 100; i++){
    if ( ( gsl_matrix_get(gff_ei, i, 0) < 10. ) || ( gsl_matrix_get(gff_ei, i, 0) > 1000. ) ) return -1;
  }

//Feature 2: At 0 <= log_gamma_sq <= 5, min(gff(log_u)) \approx 0.1
  for (int i = 37; i < 69; i++){
    gsl_vector_view row = gsl_matrix_row (gff_ei, i);
    printf("%.5e\n", gsl_vector_min(&row.vector));
    if ( ( gsl_vector_min(&row.vector) > 0.1 ) || (gsl_vector_min(&row.vector) < 0.01) ) return -1;
  }


//Verify the qualitative features of VH15 Fig. 5 right panel (variable log_gamma_sq, constant log_u)

//Feature 1: 1.25 < log_10(gff) < 2.25 for log_gamma_sq = -6 and log_u <= 0
  for ( int j = 0; j < 61; j++ ){
    if ( ( gsl_matrix_get(gff_ei, 0, j) < 18. ) || ( gsl_matrix_get(gff_ei, 0, j) > 175. ) ) return -1;
    if ( ( gsl_matrix_get(gff_ei, 99, j) < 0.1 ) || ( gsl_matrix_get(gff_ei, 99, j) > 18. ) ) return -1;
  }

//Feature 2: min(gff(log_gamma_sq)) \approx 0.1 at log_u > 5
  for (int j = 80; j < 100; j++){
    gsl_vector_view column = gsl_matrix_column (gff_ei, j);
    printf("%.5e\n", gsl_vector_min(&column.vector));
    if ( ( gsl_vector_min(&column.vector) > 0.1 ) || (gsl_vector_min(&column.vector) < 0.01) ) return -1;
  }

//Save the matrix for plotting in Python if desired
//  FILE * fp = fopen("gff_ei.dat", "w");
//  gsl_matrix_fprintf(fp, gff_ei, "%e");

  gsl_vector_free (log_u);
  gsl_vector_free (log_gamma_sq);
  gsl_matrix_free (gff_ei);

  return 0;

}
