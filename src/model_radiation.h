#ifndef MODEL_RADIATION_H
#define MODEL_RADIATION_H

// these parameters control which radiation
// model is used during the calculation.

#define EDF_MAXWELL_JUTTNER (1)
#define EDF_KAPPA_FIXED (2)
#define EDF_POWER_LAW (4)

#define MODEL_EDF EDF_MAXWELL_JUTTNER

// these parameters define properties for
// the different models

// fixed kappa model
extern double model_kappa;  // defined in radiation.c

// powerlaw model
extern double powerlaw_gamma_cut;  // all defined in radiation.c
extern double powerlaw_gamma_min;
extern double powerlaw_gamma_max;
extern double powerlaw_p;

void try_set_radiation_parameter(const char *line);

#endif // MODEL_RADIATION_H
