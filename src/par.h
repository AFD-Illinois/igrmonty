#ifndef PAR_H
#define PAR_H

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>

#include "model.h"

#define TYPE_INT (1)
#define TYPE_DBL (2)
#define TYPE_STR (3)

// feel free to change any part of this structure
typedef struct params_t {
  int seed;

  double Ns;
  double MBH;
  double M_unit;
  double mdot;  // in units of MdotEddington
  const char dump[256];
  const char spectrum[256];

  // bias
  double biasTuning;
  int    fitBias;
  double fitBiasNs;
  double targetRatio;

  // two point model
  double lnumin, lnumax, alpha_spec;

  // electron temperature models
  double TP_OVER_TE;
  double beta_crit;
  double trat_small;
  double trat_large;
  double Thetae_max;

  char loaded;
} Params;

// if you modify the 'parameters' struct above, you'll need to
// modify these functions as well to deal with the changes
void load_par_from_argv(int argc, char *argv[], Params *params);
void load_par(const char *, Params *);

// only modify if you add/modify types
void read_param(const char *, const char *, void *, int);

#endif // PAR_H
