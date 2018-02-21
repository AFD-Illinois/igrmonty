#ifndef PAR_H
#define PAR_H

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>

#define TYPE_INT (1)
#define TYPE_DBL (2)
#define TYPE_STR (3)

// feel free to change any part of this structure
typedef struct params_t {
  double Ns;
  double MBH;
  double M_unit;
  double TP_OVER_TE;
  const char dump[256];
  const char spectrum[256];
  char loaded;
} Params;

// if you modify the 'parameters' struct above, you'll need to
// modify this function as well to deal with the changes
void load_par(const char *, Params *);

// only modify if you add/modify types
void read_param(const char *, const char *, void *, int);

#endif // PAR_H
