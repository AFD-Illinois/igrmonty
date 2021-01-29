#include "decs.h"
#include "coordinates.h"

/*
 *  translates geodesic coordinates to a grid zone and returns offset
 *  for interpolation purposes. integer index corresponds to the zone
 *  center "below" the desired point and del[i] \in [0,1) returns the
 *  offset from that zone center. note that the indices and deltas at
 *  the edges are different for x1,x2,x3 and compared to ipole.
 *
 *  0    0.5    1
 *  [     |     ]
 *  A  B  C DE  F
 *
 *  startx = 0.
 *  dx = 0.5
 *
 *  A -> ( 0, 0.0)  [ x1, x2 ]
 *  A -> ( 1, 0.5)  [ x3 ]
 *  B -> ( 0, 0.0)
 *  C -> ( 0, 0.5)
 *  D -> ( 0, 0.9)
 *  E -> ( 1, 0.0)
 *  F -> ( 1, 0.0)  [ x1, x2 ]
 *  F -> ( 1, 0.5)  [ x3 ]
 *
 */
void Xtoijk(const double X[NDIM], int *i, int *j, int *k, double del[NDIM])
{
  // unless we're reading from data, i,j,k are the normal expected thing
  double phi;
  double XG[4];

  if (METRIC_eKS) {
    // the geodesics are evolved in eKS so invert through KS -> zone coordinates
    double Xks[4] = { X[0], exp(X[1]), M_PI*X[2], X[3] };
    if (METRIC_MKS3) {
      double H0 = mks3H0, MY1 = mks3MY1, MY2 = mks3MY2, MP0 = mks3MP0;
      double KSx1 = Xks[1], KSx2 = Xks[2];
      XG[0] = Xks[0];
      XG[1] = log(Xks[1] - mks3R0);
      XG[2] = (-(H0*pow(KSx1,MP0)*M_PI) - pow(2.,1. + MP0)*H0*MY1*M_PI + 
        2.*H0*pow(KSx1,MP0)*MY1*M_PI + pow(2.,1. + MP0)*H0*MY2*M_PI + 
        2.*pow(KSx1,MP0)*atan(((-2.*KSx2 + M_PI)*tan((H0*M_PI)/2.))/M_PI))/(2.*
        H0*(-pow(KSx1,MP0) - pow(2.,1 + MP0)*MY1 + 2.*pow(KSx1,MP0)*MY1 + 
          pow(2.,1. + MP0)*MY2)*M_PI);
      XG[3] = Xks[3];
    }
  } else {
    MULOOP XG[mu] = X[mu];
  }

  // the X[3] coordinate is allowed to vary so first map it to [0, stopx[3])
  phi = fmod(XG[3], stopx[3]);
  if (phi < 0.0) phi = stopx[3]+phi;

  // get provisional zone index. see note above function for details. note we
  // shift to zone centers because that's where variables are most exact.
  *i = (int) ((XG[1] - startx[1]) / dx[1] - 0.5 + 1000) - 1000;
  *j = (int) ((XG[2] - startx[2]) / dx[2] - 0.5 + 1000) - 1000;
  *k = (int) ((phi  - startx[3]) / dx[3] - 0.5 + 1000) - 1000;  

  // don't allow "center zone" to be outside of [0,N*-1]. this will often fire
  // for exotic corodinate systems and occasionally for normal ones. wrap x3.
  if (*i < 0) *i = 0;
  if (*j < 0) *j = 0;
  if (*k < 0) *k = 0;
  if (*i > N1-2) *i = N1-2; 
  if (*j > N2-2) *j = N2-2; 
  if (*k > N3-1) *k = N3-1; 

  // now construct del
  del[1] = (XG[1] - ((*i + 0.5) * dx[1] + startx[1])) / dx[1];
  del[2] = (XG[2] - ((*j + 0.5) * dx[2] + startx[2])) / dx[2];
  del[3] = (phi - ((*k + 0.5) * dx[3] + startx[3])) / dx[3];

  // finally enforce limits on del
  if (del[1] > 1.) del[1] = 1.;
  if (del[1] < 0.) del[1] = 0.;
  if (del[2] > 1.) del[2] = 1.;
  if (del[2] < 0.) del[2] = 0.;
  if (del[3] > 1.) del[3] = 1.;
  if (del[3] < 0.) {
    int oldk = *k;
    *k = N3-1;
    del[3] += 1.;
    if (del[3] < 0) {
      fprintf(stderr, " ! unable to resolve X[3] coordinate to zone %d %d %g %g\n", oldk, *k, del[3], XG[3]);
      exit(-7);
    }
  }
}

// return geodesic coordinates associated with center of zone i,j,k
void ijktoX(int i, int j, int k, double *X)
{
  // first do the naive thing 
  X[1] = startx[1] + (i+0.5)*dx[1];
  X[2] = startx[2] + (j+0.5)*dx[2];
  X[3] = startx[3] + (k+0.5)*dx[3];

  // now transform to geodesic coordinates if necessary by first
  // converting to KS and then to destination coordinates (eKS).
  if (METRIC_eKS) {
      double xKS[4] = { 0 };
    if (METRIC_MKS3) {
      double x0 = X[0];
      double x1 = X[1];
      double x2 = X[2];
      double x3 = X[3];

      double H0 = mks3H0;
      double MY1 = mks3MY1;
      double MY2 = mks3MY2;
      double MP0 = mks3MP0;
      
      xKS[0] = x0;
      xKS[1] = exp(x1) + mks3R0;
      xKS[2] = (M_PI*(1+1./tan((H0*M_PI)/2.)*tan(H0*M_PI*(-0.5+(MY1+(pow(2,MP0)*(-MY1+MY2))/pow(exp(x1)+R0,MP0))*(1-2*x2)+x2))))/2.;
      xKS[3] = x3;
    }
    
    X[0] = xKS[0];
    X[1] = log(xKS[1]);
    X[2] = xKS[2] / M_PI;
    X[3] = xKS[3];
  }
}

// Legacy naming
void coord(int i, int j, int k, double *X) { ijktoX(i, j, k, X); }
