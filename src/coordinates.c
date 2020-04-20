#include "decs.h"

////////////////////////////////// COORDINATES /////////////////////////////////

// metric parameters
//  note: if METRIC_eKS, then the code will use "exponentialKS" coordinates
//        defined by x^a = { x^0, log(x^1), x^2, x^3 } where x^0, x^1, x^2,
//        x^3 are normal KS coordinates. in addition you must set METRIC_*
//        in order to specify how Xtoijk and gdet_zone should work.
int METRIC_eKS;
int with_derefine_poles, METRIC_MKS3;
double poly_norm, poly_xt, poly_alpha, mks_smooth; // mmks
double mks3R0, mks3H0, mks3MY1, mks3MY2, mks3MP0; // mks3

// coordinate functions
static inline __attribute__((always_inline)) void set_dxdX_metric(double X[NDIM], double dxdX[NDIM][NDIM], int metric);
static inline __attribute__((always_inline)) void gcov_ks(double r, double th, double gcov[NDIM][NDIM]);

int X_in_domain(double X[NDIM]) {
  // returns 1 if X is within the computational grid.
  // checks different sets of coordinates depending on
  // specified grid coordinates

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

void set_dxdX(double X[NDIM], double dxdX[NDIM][NDIM])
{
  set_dxdX_metric(X, dxdX, 0);
}

static inline __attribute__((always_inline)) void set_dxdX_metric(double X[NDIM], double dxdX[NDIM][NDIM], int metric)
{
  // Jacobian with respect to KS basis where X is given in
  // non-KS basis
  MUNULOOP dxdX[mu][nu] = 0.;

  if ( METRIC_eKS && metric==0 ) {

    MUNULOOP dxdX[mu][nu] = mu==nu ? 1 : 0;
    dxdX[1][1] = exp(X[1]);
    dxdX[2][2] = M_PI;

  } else if ( METRIC_MKS3 ) {
   
    // mks3 ..
    dxdX[0][0] = 1.;
    dxdX[1][1] = exp(X[1]);
    dxdX[2][1] = -(pow(2.,-1. + mks3MP0)*exp(X[1])*mks3H0*mks3MP0*(mks3MY1 - 
              mks3MY2)*pow(M_PI,2)*pow(exp(X[1]) + mks3R0,-1 - mks3MP0)*(-1 + 
              2*X[2])*1./tan((mks3H0*M_PI)/2.)*pow(1./cos(mks3H0*M_PI*(-0.5 + (mks3MY1 + 
              (pow(2,mks3MP0)*(-mks3MY1 + mks3MY2))/pow(exp(X[1]) + mks3R0,mks3MP0))*(1 - 
              2*X[2]) + X[2])),2));
    dxdX[2][2]= (mks3H0*pow(M_PI,2)*(1 - 2*(mks3MY1 + (pow(2,mks3MP0)*(-mks3MY1 +
             mks3MY2))/pow(exp(X[1]) + mks3R0,mks3MP0)))*1./tan((mks3H0*M_PI)/2.)*
             pow(1./cos(mks3H0*M_PI*(-0.5 + (mks3MY1 + (pow(2,mks3MP0)*(-mks3MY1 +
             mks3MY2))/pow(exp(X[1]) + mks3R0,mks3MP0))*(1 - 2*X[2]) + X[2])),2))/2.;
    dxdX[3][3] = 1.;

  } else if ( with_derefine_poles ) {

    // mmks
    dxdX[0][0] = 1.;
    dxdX[1][1] = exp(X[1]);
    dxdX[2][1] = -exp(mks_smooth*(startx[1]-X[1]))*mks_smooth*(
      M_PI/2. -
      M_PI*X[2] +
      poly_norm*(2.*X[2]-1.)*(1+(pow((-1.+2*X[2])/poly_xt,poly_alpha))/(1 + poly_alpha)) -
      1./2.*(1. - hslope)*sin(2.*M_PI*X[2])
      );
    dxdX[2][2] = M_PI + (1. - hslope)*M_PI*cos(2.*M_PI*X[2]) +
      exp(mks_smooth*(startx[1]-X[1]))*(
        -M_PI +
        2.*poly_norm*(1. + pow((2.*X[2]-1.)/poly_xt,poly_alpha)/(poly_alpha+1.)) +
        (2.*poly_alpha*poly_norm*(2.*X[2]-1.)*pow((2.*X[2]-1.)/poly_xt,poly_alpha-1.))/((1.+poly_alpha)*poly_xt) -
        (1.-hslope)*M_PI*cos(2.*M_PI*X[2])
        );
    dxdX[3][3] = 1.;

  } else {

    // mks
    dxdX[0][0] = 1.;
    dxdX[1][1] = exp(X[1]);
    dxdX[2][2] = M_PI - (hslope - 1.)*M_PI*cos(2.*M_PI*X[2]);
    dxdX[3][3] = 1.;

  }

}

// compute KS metric at point (r,th) in KS coordinates (cyclic in t, ph)
static inline __attribute__((always_inline)) void gcov_ks(double r, double th, double gcov[NDIM][NDIM])
{
  double cth = cos(th);
  double sth = sin(th);

  double s2 = sth*sth;
  double rho2 = r*r + a*a*cth*cth;

  // compute ks metric for ks coordinates (cyclic in t,phi)
  gcov[0][0] = -1. + 2.*r/rho2;
  gcov[0][1] = 2.*r/rho2;
  gcov[0][3] = -2.*a*r*s2/rho2;

  gcov[1][0] = gcov[0][1];
  gcov[1][1] = 1. + 2.*r/rho2;
  gcov[1][3] = -a*s2*(1. + 2.*r/rho2);

  gcov[2][2] = rho2;

  gcov[3][0] = gcov[0][3];
  gcov[3][1] = gcov[1][3];
  gcov[3][3] = s2*(rho2 + a*a*s2*(1. + 2.*r/rho2));
}

// return the gdet associated with zone coordinates for the zone at
// i,j,k
double gdet_zone(int i, int j, int k) 
{
  // get the X for the zone (in geodesic coordinates for bl_coord) 
  // and in zone coordinates (for set_dxdX_metric)
  double X[NDIM], Xzone[NDIM];
  ijktoX(i,j,k, X);
  Xzone[0] = 0.;
  Xzone[1] = startx[1] + (i+0.5)*dx[1];
  Xzone[2] = startx[2] + (j+0.5)*dx[2];
  Xzone[3] = startx[3] + (k+0.5)*dx[3];

  // then get gcov for the zone (in zone coordinates)
  double gcovKS[NDIM][NDIM], gcov[NDIM][NDIM];
  double r, th;
  double dxdX[NDIM][NDIM];
  MUNULOOP gcovKS[mu][nu] = 0.;
  MUNULOOP gcov[mu][nu] = 0.;
  bl_coord(X, &r, &th);
  gcov_ks(r, th, gcovKS);
  set_dxdX_metric(Xzone, dxdX, 1);
  MUNULOOP {
    for (int lam=0; lam<NDIM; ++lam) {
      for (int kap=0; kap<NDIM; ++kap) {
        gcov[mu][nu] += gcovKS[lam][kap]*dxdX[lam][mu]*dxdX[kap][nu];
      }
    }
  }

  return gdet_func(gcov); 
}

// returns BL.{r,th} == KS.{r,th} of point with geodesic coordinates X
void bl_coord(double *X, double *r, double *th)
{
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