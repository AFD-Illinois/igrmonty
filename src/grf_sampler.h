#ifndef GRF_SAMPLER_H
#define GRF_SAMPLER_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

void generate_grf_magnetic_field_cartesian(int N, double L, double**** grf_field, double L0, double amp_fac, double alpha);

void sample_grf_magnetic_field(int N, double L, double**** grf_field, double r, double theta, double phi, double* B);

void write_grf_field_to_file(const char* filename, int N, double L, double**** grf_field);

#endif