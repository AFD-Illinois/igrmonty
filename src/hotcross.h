
double getnorm_dNdg(double thetae, radiation_params *rpars);
double total_compton_cross_num(double w, double thetae, double norm, radiation_params *rpars);
double dNdgammae(double thetae, double gammae, radiation_params *rpars);

// Function prototypes for anisotropic edf
double beta(double gammae);
double hc_klein_nishina(double we);
double dist_func_norm();
double check_scattering_limits(double photon_energy, double thetae_perp);
double hotcross_integrand_bimaxwell(double p_par, double p_perp, double phi, double photon_energy, double A, double xi, double thetae_perp, double ne);
double tpltrap(double photon_energy, double A, double xi, double thetae_perp, double ne, double p_perp_max, double p_par_max);
double compute_hotcross_anisotropic(double photon_energy, double A, double xi, double thetae_perp, double ne);
