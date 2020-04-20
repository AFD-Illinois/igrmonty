#pragma once

// metric parameters
//  note: if METRIC_eKS, then the code will use "exponentialKS" coordinates
//        defined by x^a = { x^0, log(x^1), x^2, x^3 } where x^0, x^1, x^2,
//        x^3 are normal KS coordinates. in addition you must set METRIC_*
//        in order to specify how Xtoijk and gdet_zone should work.
extern int METRIC_eKS;
extern int with_derefine_poles, METRIC_MKS3;
extern double poly_norm, poly_xt, poly_alpha, mks_smooth; // mmks
extern double mks3R0, mks3H0, mks3MY1, mks3MY2, mks3MP0; // mks3