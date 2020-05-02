#include "decs.h"

double interp_scalar(const double X[NDIM], double ***var)
{
  double interp, del[NDIM], b1, b2;
  int i, j, k, ip1, jp1, kp1;

  Xtoijk(X, &i, &j, &k, del);

  ip1 = i + 1;
  jp1 = j + 1;
  kp1 = k + 1;

  // Hard-coded periodicity in X[3]
  if (k == N3 - 1) kp1 = 0;

  b1 = 1. - del[1];
  b2 = 1. - del[2];

  // Interpolate in x1, x2
  interp = var[i][j][k]*b1*b2 +
           var[i][jp1][k]*b1*del[2] +
           var[ip1][j][k]*del[1]*b2 +
           var[ip1][jp1][k]*del[1]*del[2];

  // Interpolate in x3
  interp = (1.-del[3])*interp +
           del[3]*(
           var[i  ][j  ][kp1]*b1*b2 +
           var[i  ][jp1][kp1]*del[2]*b1 +
           var[ip1][j  ][kp1]*del[1]*b2 +
           var[ip1][jp1][kp1]*del[1]*del[2]);

  return interp;
}

double linear_interp_weight(double nu)
{

  int i;
  double di, lnu;

  lnu = log(nu);

  di = (lnu - LNUMIN)/DLNU;
  i = (int)di;
  di = di - i;

  return exp((1. - di)*wgt[i] + di*wgt[i + 1]);
}

