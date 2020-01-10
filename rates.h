#ifndef RATES_H
#define RATES_H

// ###############################################################
// ### Define where partonic and hadronic rates should be used ###
// ###############################################################

//QGP fraction definition
const double CONST_pure_QGP_T=0.1801;	
const double CONST_pure_HG_T=0.1799;	

// ######################################################
// ### Decide how to handle large viscous corrections ###
// ######################################################

// Decide whether to regulate cases where 
// (viscous correction to photon rate)/(ideal photon rate) > viscous_rate_over_ideal_smaller_than
const bool regulate_negative_rate=true; 
const double viscous_rate_over_ideal_smaller_than=1.0;


// ##############################
// ### Available photon rates ###
// ##############################

enum chem_freezeout_temp {Tch150, Tch160, Tch165};

/****** Available rates ******/
enum rate_type {
        qgp_ideal_born_AMYfit,
        qgp_ideal_born_KLS,
        qgp_ideal_born_JF_sqrtg,
        qgp_viscous_only_born_JF_sqrtg, 
        hg_ideal_Turbide_fit,
        qgp_ideal_LO_AMYfit,
        qgp_ideal_LO_AMYfit_orig,
        qgp_ideal_born_AMY_table,
        qgp_ideal_born_AMYfit_with_cuts,
        qgp_ideal_born_AMYfit_tabulated,
        qgp_viscous_only_born_g2_sqrtg,
        qgp_viscous_only_born_g2_sqrtg_table,
        qgp_viscous_only_born_g2_sqrtg_fit_tabulated,
        hg_ideal_Turbide_fit_chem_pot_Boltz_Tch150,
        hg_ideal_Turbide_fit_chem_pot_Boltz_Tch160,
        hg_ideal_Turbide_fit_chem_pot_Boltz_Tch165,
        hg_ideal_Chun_table_CE,
        hg_viscous_Chun_table_CE,
        hg_ideal_Chun_table_PCE165,
        hg_viscous_Chun_table_PCE165,
        hg_ideal_Turbide_fit_tabulated,
        qgp_viscous_only_born_g2_sqrtg_new_deltaf,
        //hg_bulk_Chun_table_PCE165,
        hg_bulk_Chun_table_CE,
        qgp_ideal_Dusling,
        qgp_shear_viscous_Dusling,
        qgp_bulk_viscous_Dusling,
        hg_bulk_Chun_table_CE_eos_transport,
        hg_pion_brem_ideal_Rapp_fit_tabulated,
        hg_in_medium_rho_ideal_Rapp_fit_tabulated,
        hg_ideal_Turbide_fit_noPiPi_tabulated,
        hg_ideal_Zahed_Dusling_2pi_fit_tabulated,
        hg_piRhoOmega_ideal_Rapp_fit_tabulated,
        thermal_ideal,
        hadronic_ideal
};



// ###########################################
// ### struct used to define a photon rate ###
// ##########################################

struct photonRate {

	//Name used when saving the results in a file
	std::string name;

	//Set to 1 to read a table from a file, instead of using a hard-coded fit
	bool use_table_instead_of_fit;
	std::string filename_of_external_table;

	//Function returning the hard-coded fit
	double (*rate_fit_function)(double, double, double);

	//The rate can be further multiplied by this function, in case
	//additional normalisation is needed (e.g. for tables)
	double (*extra_normalisation_factor_function)(double, double, double);

	//Use fit, but tabulate it internally for speed
	bool tabulate_fit_for_speed;

	double ** tabulated_rate;
	double ** tabulated_rate_log;

	//Parameters to specify by what factor should the rate be multiplied by, if any
	bool is_qgp; //Rate is multiplied by QGP_fraction
	bool is_hg; //Rate is multiplied by (1-QGP_fraction)
	bool is_thermal; //Rate is multiplied by (1-QGP_fraction)

	bool is_shear_viscous; //Multiply rate by K_mu K_nu Pi^\mu\mu/k^2
	bool is_bulk_viscous; //Multiply rate by the bulk pressure
        // If this is a viscous rate, specify the ideal rate for the corresponding production channels
        // so that one can verify if the rate goes negative in certain parts of the fluid (which is unphysical, of course)
        enum rate_type ideal_rate_for_corresponding_production_channel;

	//Parameters used to specifiy either the table is...
	bool use_k_instead_of_kOverT_for_table;
	int number_of_points_in_kOverT;
	int number_of_points_in_temp;
	double min_temp, max_temp;
	double min_kOverT, max_kOverT;

	photonRate() {

		//Set some default values
		rate_fit_function=0;
		extra_normalisation_factor_function=0;
		is_qgp=false;
		is_hg=false;
		is_thermal=false;

		is_shear_viscous=false;
		is_bulk_viscous=false;

	}

};

// ###########################
// ### Forward declaration ###
// ###########################

void init_rates(std::map<enum rate_type, struct photonRate> * rate_list, enum rate_type rate_id);
void validate_rates(std::map<enum rate_type, struct photonRate> * rate_list);
double eval_photon_rate(std::map<enum rate_type, struct photonRate> * rate_list, enum rate_type rate_id, double kOverT, double T, double kOverTkOverTOver_e_P, double bulk_pressure, double eps_plus_P, double cs2);

// Rate functions
double rate_thermal_ideal(double,double,double);
double rate_hadronic_ideal(double kOverT, double T, double kkPiOver_e_P_k2);
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
