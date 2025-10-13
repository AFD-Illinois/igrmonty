#ifndef GRF_SAMPLER_H
#define GRF_SAMPLER_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

void generate_grf_magnetic_field_cartesian(int N, double L, double L0, double amp_fac, double alpha, double**** grf_field);

void sample_grf_magnetic_field(int N, double L, double r, double theta, double phi, double* B, double**** grf_field);

void write_grf_field_to_file(const char* filename, int N, double L, double**** grf_field);

double compute_amp_fac(int N, double L, double alpha, double L0, double B_rms_target);

#endif