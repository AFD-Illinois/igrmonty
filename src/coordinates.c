#include "decs.h"
#include "coordinates.h"

#include <assert.h>

////////////////////////////////// COORDINATES /////////////////////////////////

// metric parameters
//  note: if METRIC_eKS, then the code will use "exponentialKS" coordinates
//        defined by x^a = { x^0, log(x^1), x^2, x^3 } where x^0, x^1, x^2,
//        x^3 are normal KS coordinates. in addition you must set METRIC_*
//        in order to specify how Xtoijk and gdet_zone should work.
int METRIC_eKS = 0;
int with_derefine_poles, METRIC_MKS3, METRIC_sphMINK, METRIC_esphMINK = 0;
double poly_norm, poly_xt, poly_alpha, mks_smooth; // mmks
double mks3R0, mks3H0, mks3MY1, mks3MY2, mks3MP0; // mks3

int X_in_domain(const double X[NDIM]) {
  // returns 1 if X is within the computational grid.
  // checks different sets of coordinates depending on
  // specified grid coordinates

  if (METRIC_sphMINK) return 1;

  if (METRIC_eKS) {
    double XG[4] = { 0 };
    double Xks[4] = { X[0], exp(X[1]), M_PI*X[2], X[3] };

    if (METRIC_MKS3) {
      // if METRIC_MKS3, ignore theta boundaries
      double H0 = mks3H0, MY1 = mks3MY1, MY2 = mks3MY2, MP0 = mks3MP0;
      double KSx1 = Xks[1], KSx2 = Xks[2];
      XG[0] = Xks[0];
      XG[1] = log(Xks[1] - mks3R0);
      XG[2] = (-(H0*pow(KSx1,MP0)*M_PI) - pow(2,1 + MP0)*H0*MY1*M_PI +
        2*H0*pow(KSx1,MP0)*MY1*M_PI + pow(2,1 + MP0)*H0*MY2*M_PI +
        2*pow(KSx1,MP0)*atan(((-2*KSx2 + M_PI)*tan((H0*M_PI)/2.))/M_PI))/(2.*
        H0*(-pow(KSx1,MP0) - pow(2,1 + MP0)*MY1 + 2*pow(KSx1,MP0)*MY1 +
          pow(2,1 + MP0)*MY2)*M_PI);
      XG[3] = Xks[3];

      if (XG[1] < startx[1] || XG[1] > stopx[1]) return 0;
    }

  } else {
    if(X[1] < startx[1] ||
       X[1] > stopx[1]  ||
       X[2] < startx[2] ||
       X[2] > stopx[2]) {
      return 0;
    }
  }

  return 1;
}

void set_dxdX(const double X[NDIM], double dxdX[NDIM][NDIM])
{
  set_dxdX_metric(X, dxdX, 0);
}

// return the gdet associated with zone coordinates for the zone at
// i,j,k
double gdet_zone(int i, int j, int k)
{
  // TODO support extra zone coordinates (gdet versus gdetzone)

  // get the X for the zone (in geodesic coordinates for bl_coord)
  // and in zone coordinates (for set_dxdX_metric)
  double X[NDIM], Xzone[NDIM];
  ijktoX(i,j,k, X);
  Xzone[0] = 0.;
  Xzone[1] = startx[1] + (i+0.5)*dx[1];
  Xzone[2] = startx[2] + (j+0.5)*dx[2];
  Xzone[3] = startx[3] + (k+0.5)*dx[3];

  if (METRIC_sphMINK || METRIC_esphMINK) {
    double gcov[NDIM][NDIM];
    gcov_func(Xzone, gcov);
    return gdet_func(gcov);
  }

  // then get gcov for the zone (in zone coordinates)
  double gcovKS[NDIM][NDIM], gcov[NDIM][NDIM];
  double r, th;
  double dxdX[NDIM][NDIM];
  bl_coord(X, &r, &th);
  gcov_ks(r, th, gcovKS);
  set_dxdX_metric(Xzone, dxdX, 1);
  MUNULOOP {
    gcov[mu][nu] = 0.;
    for (int lam=0; lam<NDIM; ++lam) {
      for (int kap=0; kap<NDIM; ++kap) {
        gcov[mu][nu] += gcovKS[lam][kap]*dxdX[lam][mu]*dxdX[kap][nu];
      }
    }
  }

  return gdet_func(gcov);
}

void bl_to_ks(const double X[NDIM], double ucon_bl[NDIM], double ucon_ks[NDIM])
{
  double r, th;
  bl_coord(X, &r, &th);

  double trans[NDIM][NDIM];
  MUNULOOP
    trans[mu][nu] = delta(mu, nu);

  trans[0][1] = 2. * r / (r * r - 2. * r + a * a);
  trans[3][1] = a / (r * r - 2. * r + a * a);

  MULOOP
    ucon_ks[mu] = 0.;
  MUNULOOP
    ucon_ks[mu] += trans[mu][nu] * ucon_bl[nu];
}

void vec_from_ks(const double X[NDIM], double v_ks[NDIM], double v_nat[NDIM]) {
  double dXdx[NDIM][NDIM], dxdX[NDIM][NDIM];
  set_dxdX(X, dxdX);
  int sing = invert_matrix(dxdX, dXdx);
  if (sing != 0) {
    fprintf(stderr, "failure to transform from ks at\n");
    print_vector("X", X);
  }

  assert(sing == 0);

  MULOOP v_nat[mu] = 0.;
  MUNULOOP v_nat[mu] += dXdx[mu][nu] * v_ks[nu];
}

// returns BL.{r,th} == KS.{r,th} of point with geodesic coordinates X
void bl_coord(const double *X, double *r, double *th)
{
  if (METRIC_sphMINK) {
    *r = X[1];
    *th = X[2];
    return;
  } else if (METRIC_esphMINK) {
    *r = exp(X[1]);
    *th = X[2];
    return;
  }

  *r = exp(X[1]);

  if (METRIC_eKS) {
    *r = exp(X[1]);
    *th = M_PI * X[2];
  } else if (METRIC_MKS3) {
    *r = exp(X[1]) + mks3R0;
    *th = (M_PI*(1. + 1./tan((mks3H0*M_PI)/2.)*tan(mks3H0*M_PI*(-0.5 + (mks3MY1 + (pow(2.,mks3MP0)*(-mks3MY1 + mks3MY2))/pow(exp(X[1])+mks3R0,mks3MP0))*(1. - 2.*X[2]) + X[2]))))/2.;
  } else if (with_derefine_poles) {
    double thG = M_PI*X[2] + ((1. - hslope)/2.)*sin(2.*M_PI*X[2]);
    double y = 2*X[2] - 1.;
    double thJ = poly_norm*y*(1. + pow(y/poly_xt,poly_alpha)/(poly_alpha+1.)) + 0.5*M_PI;
    *th = thG + exp(mks_smooth*(startx[1] - X[1]))*(thJ - thG);
  } else {
    *th = M_PI*X[2] + ((1. - hslope)/2.)*sin(2.*M_PI*X[2]);
  }
}

// cribbed from ipole
/*
 invert_matrix():
 Uses LU decomposition and back substitution to invert a matrix
 A[][] and assigns the inverse to Ainv[][].  This routine does not
 destroy the original matrix A[][].
 Returns (1) if a singular matrix is found,  (0) otherwise.
 */
int invert_matrix(double Am[][NDIM], double Aminv[][NDIM])
{

  int i, j;
  int n = NDIM;
  int permute[NDIM];
  double dxm[NDIM], Amtmp[NDIM][NDIM];

  for (i = 0; i < NDIM; i++)
    for (j = 0; j < NDIM; j++) {
      Amtmp[i][j] = Am[i][j];
    }

  // Get the LU matrix:
  if (LU_decompose(Amtmp, permute) != 0) {
    //fprintf(stderr, "invert_matrix(): singular matrix encountered! \n");
    return (1);
  }

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      dxm[j] = 0.;
    }
    dxm[i] = 1.;

    /* Solve the linear system for the i^th column of the inverse matrix: :  */
    LU_substitution(Amtmp, dxm, permute);

    for (j = 0; j < n; j++) {
      Aminv[j][i] = dxm[j];
    }

  }

  return (0);
}

/*
 LU_decompose():
 Performs a LU decomposition of the matrix A using Crout's method
 with partial implicit pivoting.  The exact LU decomposition of the
 matrix can be reconstructed from the resultant row-permuted form via
 the integer array permute[]
 The algorithm closely follows ludcmp.c of "Numerical Recipes
 in C" by Press et al. 1992.
 This will be used to solve the linear system  A.x = B
 Returns (1) if a singular matrix is found,  (0) otherwise.
*/
int LU_decompose(double A[][NDIM], int permute[])
{

  const double absmin = 1.e-30; /* Value used instead of 0 for singular matrices */

  //static double row_norm[NDIM];
  double row_norm[NDIM];
  double absmax, maxtemp;

  int i, j, k, max_row;
  int n = NDIM;

  max_row = 0;

  /* Find the maximum elements per row so that we can pretend later
   we have unit-normalized each equation: */

  for (i = 0; i < n; i++) {
    absmax = 0.;

    for (j = 0; j < n; j++) {

      maxtemp = fabs(A[i][j]);

      if (maxtemp > absmax) {
        absmax = maxtemp;
      }
    }

    /* Make sure that there is at least one non-zero element in this row: */
    if (absmax == 0.) {
      //fprintf(stderr, "LU_decompose(): row-wise singular matrix!\n");
      return (1);
    }

    row_norm[i] = 1. / absmax; /* Set the row's normalization factor. */
  }

  /* The following the calculates the matrix composed of the sum
   of the lower (L) tridagonal matrix and the upper (U) tridagonal
   matrix that, when multiplied, form the original maxtrix.
   This is what we call the LU decomposition of the maxtrix.
   It does this by a recursive procedure, starting from the
   upper-left, proceding down the column, and then to the next
   column to the right.  The decomposition can be done in place
   since element {i,j} require only those elements with {<=i,<=j}
   which have already been computed.
   See pg. 43-46 of "Num. Rec." for a more thorough description.
   */

  /* For each of the columns, starting from the left ... */
  for (j = 0; j < n; j++) {

    /* For each of the rows starting from the top.... */

    /* Calculate the Upper part of the matrix:  i < j :   */
    for (i = 0; i < j; i++) {
      for (k = 0; k < i; k++) {
        A[i][j] -= A[i][k] * A[k][j];
      }
    }

    absmax = 0.0;

    /* Calculate the Lower part of the matrix:  i <= j :   */

    for (i = j; i < n; i++) {

      for (k = 0; k < j; k++) {
        A[i][j] -= A[i][k] * A[k][j];
      }

      /* Find the maximum element in the column given the implicit
       unit-normalization (represented by row_norm[i]) of each row:
       */
      maxtemp = fabs(A[i][j]) * row_norm[i];

      if (maxtemp >= absmax) {
        absmax = maxtemp;
        max_row = i;
      }

    }

    /* Swap the row with the largest element (of column j) with row_j.  absmax
     This is the partial pivoting procedure that ensures we don't divide
     by 0 (or a small number) when we solve the linear system.
     Also, since the procedure starts from left-right/top-bottom,
     the pivot values are chosen from a pool involving all the elements
     of column_j  in rows beneath row_j.  This ensures that
     a row  is not permuted twice, which would mess things up.
     */
    if (max_row != j) {

      /* Don't swap if it will send a 0 to the last diagonal position.
       Note that the last column cannot pivot with any other row,
       so this is the last chance to ensure that the last two
       columns have non-zero diagonal elements.
       */

      if ((j == (n - 2)) && (A[j][j + 1] == 0.)) {
        max_row = j;
      } else {
        for (k = 0; k < n; k++) {

          maxtemp = A[j][k];
          A[j][k] = A[max_row][k];
          A[max_row][k] = maxtemp;

        }

        /* Don't forget to swap the normalization factors, too...
         but we don't need the jth element any longer since we
         only look at rows beneath j from here on out.
         */
        row_norm[max_row] = row_norm[j];
      }
    }

    /* Set the permutation record s.t. the j^th element equals the
     index of the row swapped with the j^th row.  Note that since
     this is being done in successive columns, the permutation
     vector records the successive permutations and therefore
     index of permute[] also indexes the chronology of the
     permutations.  E.g. permute[2] = {2,1} is an identity
     permutation, which cannot happen here though.
     */

    permute[j] = max_row;

    if (A[j][j] == 0.) {
      A[j][j] = absmin;
    }

    /* Normalize the columns of the Lower tridiagonal part by their respective
     diagonal element.  This is not done in the Upper part because the
     Lower part's diagonal elements were set to 1, which can be done w/o
     any loss of generality.
     */
    if (j != (n - 1)) {
      maxtemp = 1. / A[j][j];

      for (i = (j + 1); i < n; i++) {
        A[i][j] *= maxtemp;
      }
    }

  }

  return (0);

  /* End of LU_decompose() */

}

/*
 LU_substitution():
 Performs the forward (w/ the Lower) and backward (w/ the Upper)
 substitutions using the LU-decomposed matrix A[][] of the original
 matrix A' of the linear equation:  A'.x = B.  Upon entry, A[][]
 is the LU matrix, B[] is the source vector, and permute[] is the
 array containing order of permutations taken to the rows of the LU
 matrix.  See LU_decompose() for further details.
 Upon exit, B[] contains the solution x[], A[][] is left unchanged.
*/
void LU_substitution(double A[][NDIM], double B[], int permute[])
{
  int i, j;
  int n = NDIM;
  double tmpvar;

  /* Perform the forward substitution using the LU matrix.
   */
  for (i = 0; i < n; i++) {

    /* Before doing the substitution, we must first permute the
     B vector to match the permutation of the LU matrix.
     Since only the rows above the currrent one matter for
     this row, we can permute one at a time.
     */
    tmpvar = B[permute[i]];
    B[permute[i]] = B[i];
    for (j = (i - 1); j >= 0; j--) {
      tmpvar -= A[i][j] * B[j];
    }
    B[i] = tmpvar;
  }

  /* Perform the backward substitution using the LU matrix.
   */
  for (i = (n - 1); i >= 0; i--) {
    for (j = (i + 1); j < n; j++) {
      B[i] -= A[i][j] * B[j];
    }
    B[i] /= A[i][i];
  }

  /* End of LU_substitution() */
}

