#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <time.h>
#include <omp.h>

#define IDX(i,j,k) ((i)*N*N + (j)*N + (k))

static inline int wrap_k(int k, int N) {
    return (k < N/2) ? k : k - N;
}

static double gauss_rand() {
    static int has_spare = 0;
    static double spare;
    if (has_spare) {
        has_spare = 0;
        return spare;
    }
    has_spare = 1;
    double u, v, s;
    do {
        u = 2.0 * rand() / RAND_MAX - 1.0;
        v = 2.0 * rand() / RAND_MAX - 1.0;
        s = u*u + v*v;
    } while (s >= 1.0 || s == 0.0);
    s = sqrt(-2.0 * log(s) / s);
    spare = v * s;
    return u * s;
}

void generate_grf_magnetic_field_cartesian(int N, double L, double L0, double amp_fac, double alpha, double**** grf_field) {
/**
 * Generate a Gaussian random magnetic field in a periodic Cartesian box. Note the magnetic field is exactly divergence free in k-space, but not in real space.
 *
 * The power spectrum of B is
 *     P(k) = amp_fac / (k^2 + k0^2)^(alpha/2).
 *
 * Parameters
 * ----------
 * N : int
 *     Grid size (N x N x N).
 * L : double
 *     Domain length.
 * alpha : double
 *     Power spectrum slope parameter.
 * amp_fac : double
 *     Scaling factor for the power spectrum (of A).
 * L0 : double
 *     Correlation length cutoff for power spectrum. k0=2*pi/L0
 *
 * Returns
 * -------
 * grf_field : double**** (N x N x N x 3)
 *     Magnetic field components (Bx, By, Bz).
 */

    int i, j, k;
    double dk = 2.0 * M_PI / L;
    double k0 = 2.0 * M_PI / L0;
    // allocate Fourier arrays
    fftw_complex *B_kx = fftw_malloc(sizeof(fftw_complex) * N * N * N);
    fftw_complex *B_ky = fftw_malloc(sizeof(fftw_complex) * N * N * N);
    fftw_complex *B_kz = fftw_malloc(sizeof(fftw_complex) * N * N * N);

    fftw_plan plan_x, plan_y, plan_z;
    fftw_complex *Bx = fftw_malloc(sizeof(fftw_complex) * N * N * N);
    fftw_complex *By = fftw_malloc(sizeof(fftw_complex) * N * N * N);
    fftw_complex *Bz = fftw_malloc(sizeof(fftw_complex) * N * N * N);

    plan_x = fftw_plan_dft_3d(N, N, N, B_kx, Bx, FFTW_BACKWARD, FFTW_ESTIMATE);
    plan_y = fftw_plan_dft_3d(N, N, N, B_ky, By, FFTW_BACKWARD, FFTW_ESTIMATE);
    plan_z = fftw_plan_dft_3d(N, N, N, B_kz, Bz, FFTW_BACKWARD, FFTW_ESTIMATE);


    // compute Fourier coefficients
    for (i = 0; i < N; i++) {
        double kx = (i < N/2 ? i : i-N) * dk;
        for (j = 0; j < N; j++) {
            double ky = (j < N/2 ? j : j-N) * dk;
            for (k = 0; k < N; k++) {
                double kz = (k < N/2 ? k : k-N) * dk;
                int idx = (i*N + j)*N + k;

                double ksq = kx*kx + ky*ky + kz*kz;

                double amp = 0.0;
                if (ksq > 0) {
                    amp = sqrt(amp_fac) / pow(ksq + k0*k0, alpha/4.0);
                }


                fftw_complex tmpx = (gauss_rand() + I*gauss_rand());
                fftw_complex tmpy = (gauss_rand() + I*gauss_rand());
                fftw_complex tmpz = (gauss_rand() + I*gauss_rand());

                // project out k·B (divergence-free condition)
                if (ksq > 0) {
                    fftw_complex k_dot_B = kx*tmpx + ky*tmpy + kz*tmpz;
                    tmpx -= kx * k_dot_B / ksq;
                    tmpy -= ky * k_dot_B / ksq;
                    tmpz -= kz * k_dot_B / ksq;
                }

                // normalize magnitude to 1 (if not zero)
                // double norm = cabs(tmpx) + cabs(tmpy) + cabs(tmpz);
                double normx = cabs(tmpx);
                double normy = cabs(tmpy);
                double normz = cabs(tmpz);
                if (normx > 0) tmpx /= normx;
                if (normy > 0) tmpy /= normy;
                if (normz > 0) tmpz /= normz;

                // apply amplitude and sqrt(2/3) factor
                double scale = amp * sqrt(2.0/3.0);
                B_kx[idx] = tmpx * scale;
                B_ky[idx] = tmpy * scale;
                B_kz[idx] = tmpz * scale;
            }
        }
    }

    // inverse FFT
    fftw_execute(plan_x);
    fftw_execute(plan_y);
    fftw_execute(plan_z);

    // normalize (FFTW does not normalize by N^3)
    double norm = 1.0 / (N * N * N);

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            for (k = 0; k < N; k++) {
                int idx = (i*N + j)*N + k;
                (grf_field)[i][j][k][0] = creal(Bx[idx]) * norm;
                (grf_field)[i][j][k][1] = creal(By[idx]) * norm;
                (grf_field)[i][j][k][2] = creal(Bz[idx]) * norm;
            }
        }
    }

    // cleanup
    fftw_destroy_plan(plan_x);
    fftw_destroy_plan(plan_y);
    fftw_destroy_plan(plan_z);

    fftw_free(B_kx);
    fftw_free(B_ky);
    fftw_free(B_kz);
    fftw_free(Bx);
    fftw_free(By);
    fftw_free(Bz);
}

void sample_grf_magnetic_field(int N, double L,
                              double r, double theta, double phi,
                              double* B, double**** grf_field) {
/*
Given a location in r, theta, phi, return a B field interpolated from a cartesian grid grf_field of size L^3 and resolution L/N.
*/
    double Bx, By, Bz;
    double sth=sin(theta),cth=cos(theta);
    double sphi=sin(phi),cphi=cos(phi);
    // Convert spherical (r, theta, phi) to Cartesian (x, y, z)
    double x = r * sth * cphi;
    double y = r * sth * sphi;
    double z = r * cth;

    // Map (x, y, z) to grid indices
    double grid_x = (x + L/2) * (N-1) / L;
    double grid_y = (y + L/2) * (N-1) / L;
    double grid_z = (z + L/2) * (N-1) / L;

    int ix = (int)grid_x;
    int iy = (int)grid_y;
    int iz = (int)grid_z;

    // Clamp indices to valid range
    if (ix < 0) ix = 0; if (ix > N-2) ix = N-2;
    if (iy < 0) iy = 0; if (iy > N-2) iy = N-2;
    if (iz < 0) iz = 0; if (iz > N-2) iz = N-2;

    double dx = grid_x - ix;
    double dy = grid_y - iy;
    double dz = grid_z - iz;

    // Trilinear interpolation
    double bx = 0, by = 0, bz = 0;
    for (int i = 0; i <= 1; ++i)
    for (int j = 0; j <= 1; ++j)
    for (int k = 0; k <= 1; ++k) {
        double w = ((i ? dx : 1-dx) *
                    (j ? dy : 1-dy) *
                    (k ? dz : 1-dz));
        bx += w * grf_field[ix+i][iy+j][iz+k][0];
        by += w * grf_field[ix+i][iy+j][iz+k][1];
        bz += w * grf_field[ix+i][iy+j][iz+k][2];
    }

    Bx = bx;
    By = by;
    Bz = bz;
    // Convert (bx, by, bz) from Cartesian to spherical coordinates
    double Br = bx * sth * cphi + by * sth * sphi + bz * cth;
    double Btheta = bx * cth * cphi + by * cth * sphi - bz * sth;
    double Bphi = -bx * sphi + by * cphi;

    B[0]=Br;
    B[1]=Btheta;
    B[2]=Bphi;
}

void write_grf_field_to_file(const char* filename, int N, double L, double**** grf_field) {
    FILE* fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
        return;
    }
    double dx = L / N;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                double x = -L/2 + i * dx;
                double y = -L/2 + j * dx;
                double z = -L/2 + k * dx;
                double Bx = grf_field[i][j][k][0];
                double By = grf_field[i][j][k][1];
                double Bz = grf_field[i][j][k][2];
                fprintf(fp, "%.14g %.14g %.14g %.14g %.14g %.14g\n", x, y, z, Bx, By, Bz);
            }
        }
    }
    fclose(fp);
}

double compute_amp_fac(int N, double L, double alpha, double L0, double B_rms_target) {
    
    double k0 = 2.0*M_PI/L0;
    // Allocate k-vector
    double *kx = (double*)malloc(N * sizeof(double));
    if (kx == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    // --- Compute k-vector (like np.fft.fftfreq) ---
    for (int i = 0; i < N; i++) {
        kx[i] = (i < N/2 ? i : i - N) * (2.0 * M_PI / L);
    }

    double integral = 0.0;

    // --- Parallel loop over 3D grid ---
    #pragma omp parallel for collapse(3) reduction(+:integral) schedule(dynamic)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                double kxi = kx[i];
                double kyj = kx[j];
                double kzk = kx[k];

                double ksq = kxi*kxi + kyj*kyj + kzk*kzk;

                if (ksq > 0.0) {
                    double denom = pow(ksq + k0*k0, alpha/2.0);
                    double Pk = 1.0 / denom;
                    integral += Pk;
                }
            }
        }
    }

    free(kx);

    // --- Compute amp_fac ---
    double amp_fac = (B_rms_target * B_rms_target) * pow((double)N, 6) / integral;

    return amp_fac;
}
