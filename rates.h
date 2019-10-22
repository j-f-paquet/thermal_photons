#ifndef RATES_H
#define RATES_H

void init_rates(struct photonRate * currRate, enum rate_type id);
double eval_photon_rate(const struct photonRate * currRate, double kOverT, double T, double kOverTkOverTOver_e_P, double bulk_pressure, double eps_plus_P, double cs2);

// Rate functions
double rate_thermal_ideal(double,double,double);
double rate_qgp_ideal_born_AMYfit(double, double, double);
double rate_qgp_ideal_born_KLS(double, double, double);
double rate_qgp_ideal_born_JF_sqrtg(double kOverT, double T, double kkPiOver_e_P_k2);
double rate_qgp_viscous_only_born_JF_sqrtg(double kOverT, double T, double kkPiOver_e_P_k2);
double rate_hg_ideal_Turbide_fit(double kOverT, double T, double kkPiOver_e_P_k2);
double rate_qgp_ideal_LO_AMYfit(double, double, double);
double rate_qgp_ideal_born_AMYfit_with_cuts(double kOverT, double T, double kkPiOver_e_P_k2);
double rate_qgp_viscous_only_born_g2_sqrtg(double kOverT, double T, double kkPiOver_e_P_k2); 
double rate_hg_ideal_Turbide_fit_chem_pot_Boltz_Tch150(double kOverT, double T, double kkPiOver_e_P_k2);
double rate_hg_ideal_Turbide_fit_chem_pot_Boltz_Tch160(double kOverT, double T, double kkPiOver_e_P_k2);
double rate_hg_ideal_Turbide_fit_chem_pot_Boltz_Tch165(double kOverT, double T, double kkPiOver_e_P_k2);
double rate_qgp_viscous_only_born_g2_sqrtg_fit2(double kOverT, double T, double kkPiOver_e_P_k2); 
double rate_qgp_viscous_only_born_g2_sqrtg_new_deltaf(double kOverT, double T, double kkPiOver_e_P_k2); 
double rate_qgp_viscous_only_born_g2_sqrtg_new_deltaf_cuts(double kOverT, double T, double kkPiOver_e_P_k2); 
double rate_qgp_viscous_only_born_g2_sqrtg_new_deltaf_fit2(double kOverT, double T, double kkPiOver_e_P_k2); 
double rate_qgp_ideal_Dusling(double kOverT, double T, double kkPiOver_e_P_k2);
double rate_qgp_shear_viscous_Dusling(double kOverT, double T, double kkPiOver_e_P_k2);
double rate_qgp_bulk_viscous_Dusling(double kOverT, double T, double kkPiOver_e_P_k2);
double rate_hg_in_medium_rho_ideal_Rapp_fit(double kOverT, double T, double kkPiOver_e_P_k2); 
double rate_hg_ideal_Turbide_noPiPi_fit(double kOverT, double T, double kkPiOver_e_P_k2);
double rate_hg_pion_brem_ideal_Rapp_fit(double kOverT, double T, double kkPiOver_e_P_k2);
double rate_hg_ideal_Zahed_Dusling_2pi_fit(double kOverT, double T, double kkPiOver_e_P_k2);
double rate_hg_piRhoOmega_ideal_Rapp_fit(double kOverT, double T, double kkPiOver_e_P_k2);

// Rate-related functions
double T2_normalisation(double kOverT, double T, double kk);
double minus_sign_normalisation(double kOverT, double T, double kk);

// Helper functions
void tabulate_fit(struct photonRate * currRate); 
void load_rate_from_file(struct photonRate * currRate);

#endif
