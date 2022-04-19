
#include <math.h>
#include "par.h"
#include "model_radiation.h"

// call this to set defaults
void load_par_from_argv(int argc, char *argv[], Params *params) {

  // set default values here
  params->seed        = -1; // will use time() to randomize seed if set to -1

  params->biasTuning  = 1.;
  params->fitBias     = 0;
  params->fitBiasNs   = 10000.;
  params->targetRatio = M_SQRT2;

  params->TP_OVER_TE = 3.;
  params->beta_crit = 1.;
  params->trat_small = 1.;
  params->trat_large = 10.;
  params->Thetae_max = 1.e100;

  params->mdot = 0;  // zero default value means ignore

  // Load parameters
  for (int i=0; i<argc-1; ++i) {
    if ( strcmp(argv[i], "-par") == 0 ) {
      load_par(argv[i+1], params);
    }
  }

  // override parameters based on what has been loaded
  if (params->mdot > 0) {
    params->M_unit = 1.;
  }
}

// sets default values for elements of params (if desired) and loads
// from par file 'fname'
void load_par (const char *fname, Params *params) {
 
  char line[256];
  FILE *fp = fopen(fname, "r");

  if (fp == NULL) {
    fprintf(stderr, "! unable to load parameter file '%s'. (%d: %s)\n", fname, errno, strerror(errno));
    exit(-1);
  }

  // modify parameters/types below
  while (fgets(line, 255, fp) != NULL) {

    if (line[0] == '#') continue; 

    read_param(line, "seed", &(params->seed), TYPE_DBL);

    read_param(line, "Ns", &(params->Ns), TYPE_DBL);
    read_param(line, "MBH", &(params->MBH), TYPE_DBL);
    read_param(line, "M_unit", &(params->M_unit), TYPE_DBL);
    read_param(line, "mdot", &(params->mdot), TYPE_DBL);
    read_param(line, "dump", (void *)(params->dump), TYPE_STR);
    read_param(line, "spectrum", (void *)(params->spectrum), TYPE_STR);

    // bias
    read_param(line, "bias",        &(params->biasTuning),  TYPE_DBL);
    read_param(line, "fit_bias",    &(params->fitBias),     TYPE_INT);
    read_param(line, "fit_bias_ns", &(params->fitBiasNs),   TYPE_DBL);
    read_param(line, "ratio",       &(params->targetRatio), TYPE_DBL);

    // two point model
    read_param(line, "lnumin", &(params->lnumin), TYPE_DBL);
    read_param(line, "lnumax", &(params->lnumax), TYPE_DBL);
    read_param(line, "alpha_spec", &(params->alpha_spec), TYPE_DBL);    

    // electron temperature models
    read_param(line, "TP_OVER_TE", &(params->TP_OVER_TE), TYPE_DBL);
    read_param(line, "beta_crit", &(params->beta_crit), TYPE_DBL);
    read_param(line, "trat_small", &(params->trat_small), TYPE_DBL);
    read_param(line, "trat_large", &(params->trat_large), TYPE_DBL);
    read_param(line, "Thetae_max", &(params->Thetae_max), TYPE_DBL);

    // set model parameters
    try_set_radiation_parameter(line);
  }

  fclose(fp);

  params->loaded = 1;

}

// loads value -> (type *)*val if line corresponds to (key,value)
void read_param (const char *line, const char *key, void *val, int type) {

  // silence valgrind warnings
  char word[256] = {0};
  char value[256] = {0};

  sscanf(line, "%s %s", word, value);

  if (strcmp(word, key) == 0) {
    switch (type) {
    case TYPE_INT:
      sscanf(value, "%d", (int *)val);
      break;
    case TYPE_DBL:
      sscanf(value, "%lf", (double *)val);
      break;
    case TYPE_STR:
      sscanf(value, "%s", (char *)val);
      break;
    default:
      fprintf(stderr, "! attempt to load unknown type '%d'.\n", type);
    }
  }

}



