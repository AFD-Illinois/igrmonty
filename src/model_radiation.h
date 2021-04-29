#ifndef MODEL_RADIATION_H
#define MODEL_RADIATION_H

// these parameters control which radiation
// model is used during the calculation.

#define EDF_MAXWELL_JUTTNER (1)
#define EDF_KAPPA_FIXED (2)
#define EDF_KAPPA_VARIABLE (3)
#define EDF_POWER_LAW (4)

#define MODEL_EDF EDF_MAXWELL_JUTTNER

// for variable kappa. also set behavior for how to
// deal with out-of-bounds kappas in get_model_kappa(...)
#define KAPPA_MIN (2.5)
#define KAPPA_MAX (7.5)
#define KAPPA_NSAMP (50)
#define DKAPPA (((double)KAPPA_MAX - (double)KAPPA_MIN) / KAPPA_NSAMP)

// used during par loading
void try_set_radiation_parameter(const char *line);

// these parameters define properties for
// the different models

// fixed kappa model
extern double model_kappa;  // defined in radiation.c

// powerlaw model
extern double powerlaw_gamma_cut;  // all defined in radiation.c
extern double powerlaw_gamma_min;
extern double powerlaw_gamma_max;
extern double powerlaw_p;

// used to pass radiation model parameters
typedef struct radiation_param_struct {
  double kappa;
} radiation_params;
radiation_params get_model_radiation_params(const double X[NDIM]);

#endif // MODEL_RADIATION_H
