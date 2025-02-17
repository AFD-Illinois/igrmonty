double fdist(double ge, double Thetae, double kappa);
double dfdgam(double ge, void *params);

//void sample_beta_distr_y(double, double *, double *, radiation_params *);
void sample_beta_distr_num(double, double *, double *, radiation_params *);
void sample_edf_distr_anisotropic(double Thetae_perp, double *gamma_e, double *mu, double A, double xi, radiation_params *rpars);
