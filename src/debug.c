#include <stdio.h>
#include "decs.h"

// The world needed these
// Maybe not in this form
void print_matrix(char *name, double g[NDIM][NDIM])
{
  // Print a name and a matrix
  fprintf(stderr, "%s:\n", name);
  fprintf(stderr,
          "%g\t%g\t%g\t%g\n%g\t%g\t%g\t%g\n%g\t%g\t%g\t%g\n%g\t%g\t%g\t%g\n",
          g[0][0], g[0][1], g[0][2], g[0][3], g[1][0], g[1][1], g[1][2],
          g[1][3], g[2][0], g[2][1], g[2][2], g[2][3], g[3][0], g[3][1],
          g[3][2], g[3][3]);
  // Additionally kill things if/when we hit NaNs
  MUNULOOP
    if (isnan(g[mu][nu]))
      exit(-1);
}
void print_vector(char *name, double v[NDIM])
{
  fprintf(stderr, "%s: ", name);
  fprintf(stderr,
          "%.10g\t%.10g\t%.10g\t%.10g\n",
          v[0], v[1], v[2], v[3]);
  MUNULOOP
    if (isnan(v[nu]))
      exit(-1);
}

