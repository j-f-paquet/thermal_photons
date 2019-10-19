#ifndef PHOTON_H
#define PHOTON_H

/*******************************
********** How to use **********
1) 

*********** Inputs *************
1) Spacetime grid: 
* text file: 5 floats per line (temperature, QGP fraction, ux/u0, uy/u0, uz/u0)
* binary file: blocks of 5 floats
* the spacetime position of the cell is given by the line
2) Shear grid:

*********************************/

//Libraries
#include <string>
#include <cmath>
#include <cstdio>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>

//Read parameters directly from a file generated by the hydro
#include "hydro_info_header_h"


/********* Inputs ********/
//Format of the input files
const bool CONST_binaryMode=MUSIC_outputBinaryEvolution; //0 for text, 1 for binary
//Location of the spacetime grid file
const std::string stGridFile="./evolution_xyeta.dat";

//Information about the spacetime grid
//Number of cells of the grid
const int cellNb_x=MUSIC_real_nx;
const int cellNb_y=MUSIC_real_ny;
const int cellNb_eta=MUSIC_real_neta;
//Size of the cells
const double CONST_cellsize_X=MUSIC_dx; //In fm
const double CONST_cellsize_Y=MUSIC_dy; //In fm
const double CONST_cellsize_Eta=MUSIC_deta; //In units of rapidity
//Initial time tau_0
const double CONST_tau0=MUSIC_tau0;
const double CONST_effective_dTau=MUSIC_dtau;

const bool CONST_boost_invariant=1;
const int CONST_nb_steps_eta_integration=30;
const double CONST_max_eta_integration=3.0;

//Run with viscosity or not
const bool CONST_with_shear_viscosity=MUSIC_with_shear_viscosity; //turn on and off shear viscosity
const bool CONST_with_bulk_viscosity=MUSIC_with_bulk_viscosity; //turn on and off bulk viscosity
const bool CONST_with_viscosity=CONST_with_shear_viscosity; //general flag for viscosity
//const double shear_to_s=0.08;
//Location of the viscous files
const std::string shearViscosityFile="evolution_Wmunu_over_epsilon_plus_P_xyeta.dat";
const std::string bulkViscosityFile="evolution_bulk_pressure_xyeta.dat";

//Discretization of photon spectrum
//kT
const double CONST_kTMin=0.2; //Minimum value for kT
const double CONST_kTMax=4.0; //Maximum value for kT
const int CONST_NkT=20;  //Warning: delta kT=(kTMax-kTmin)/(NkT-1) 
//const double kTdisc[3] = [0.2,4.0,0.2] //Discretization in kT: [kT min, kT max, delta kT]
const int CONST_Nphi=20;  //phi
//Rapidity
const double CONST_rapMin=-0.0; //Minimum value for rapidity
const double CONST_rapMax=0.0; //Maximum value for rapidity
const int CONST_Nrap=1;  //Warning: delta rap=(rapMax-rapmin)/(Nrap-1)

//Deltas used for the (uniform) discretization of the grid
const double CONST_delRap= CONST_Nrap > 1 ? (CONST_rapMax-CONST_rapMin)/(CONST_Nrap-1.0) : 0;
const double CONST_delPhi=(2*M_PI)/(CONST_Nphi-1.0);
const double CONST_delKt=(CONST_kTMax-CONST_kTMin)/(CONST_NkT-1.0);

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
};
/****************************/
//const enum rate_type CONST_rates_to_use[] = {qgp_ideal_born_AMYfit,qgp_ideal_born_KLS,hg_ideal_Turbide_fit,qgp_ideal_LO_AMYfit,qgp_ideal_born_AMYfit_with_cuts,qgp_ideal_born_AMY_table, qgp_ideal_born_AMYfit_tabulated};
const enum rate_type CONST_rates_to_use[] = {
qgp_ideal_born_AMYfit,
//qgp_ideal_born_KLS,
//qgp_ideal_born_JF_sqrtg,
//qgp_viscous_only_born_JF_sqrtg,
//hg_ideal_Turbide_fit,
qgp_ideal_LO_AMYfit,
//qgp_ideal_born_AMY_table,
//qgp_ideal_born_AMYfit_with_cuts,
//qgp_ideal_born_AMYfit_tabulated,
qgp_viscous_only_born_g2_sqrtg,
//qgp_viscous_only_born_g2_sqrtg_table,
//qgp_viscous_only_born_g2_sqrtg_fit_tabulated,
//hg_ideal_Turbide_fit_chem_pot_Boltz_Tch150,
//hg_ideal_Turbide_fit_chem_pot_Boltz_Tch160,
//hg_ideal_Turbide_fit_chem_pot_Boltz_Tch165,
//hg_ideal_Chun_table_CE,
hg_viscous_Chun_table_CE,
//hg_ideal_Chun_table_PCE165,
//hg_viscous_Chun_table_PCE165,
hg_ideal_Turbide_fit_tabulated,
//qgp_viscous_only_born_g2_sqrtg_new_deltaf,
//hg_bulk_Chun_table_PCE165,
hg_bulk_Chun_table_CE,
//qgp_ideal_Dusling,
//qgp_shear_viscous_Dusling,
qgp_bulk_viscous_Dusling,
//hg_bulk_Chun_table_CE_eos_transport,
hg_pion_brem_ideal_Rapp_fit_tabulated,
hg_in_medium_rho_ideal_Rapp_fit_tabulated,
hg_ideal_Turbide_fit_noPiPi_tabulated,
hg_ideal_Zahed_Dusling_2pi_fit_tabulated
};
//const int CONST_rates_to_use[] = {1,2,3,4,5,6};
const int CONST_N_rates = sizeof(CONST_rates_to_use)/sizeof(int);
//

//Mid-rapidity cut: the midrapidity result will be an average over approximatively -midRapCut to midRapCut
const double CONST_midRapCut = 0.5;
//Number of Fourier coefficient to compute
const int CONST_FourierNb = 6;

//QGP fraction definition
const double CONST_pure_QGP_T=0.1801;	
const double CONST_pure_HG_T=0.1799;	
//const double CONST_freezeout_T=MUSIC_kinetic_FO_temperature_in_GeV;
const double CONST_freezeout_T=0.145;

//Generate spectra sums from t0 to t_i with t_i \in CONST_tauList
//const double CONST_tauList[]={.6,1.0,2.0};
//const double CONST_tempList[]={.180,.500};

//General constants
const int CONST_Nc=3;
const int CONST_Nf=3;
const double CONST_CF=(CONST_Nc*CONST_Nc-1.0)/(2.0*CONST_Nc);
const double CONST_alphaEM=1/137.0;
const double CONST_alphaS=1.0/M_PI;
const double CONST_gs=sqrt(4*M_PI*CONST_alphaS);
const double CONST_mInfOverT=sqrt(CONST_CF)*CONST_gs/2.0;
const double CONST_twoPiCubed=pow(2*M_PI,3);
const double CONST_hbarc=0.1973;
const double CONST_GeV2_to_GeVm2_fmm4=1.0/(CONST_hbarc*CONST_hbarc*CONST_hbarc*CONST_hbarc);

/*************************/

//Use to store spacetime position and related informations
struct phaseSpace_pos {

	//
	double tau, x, y, eta;
	int itau, ix, iy, ieta;

	//
	int ikT, irap, iphi;

	//
	//int iTauList;
	//bool newiTau;

	//Spacetime integration weights
	double w_eta;

};

enum interp_type { linear, quadratic };

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

//Forward declaration
void photon_prod(const struct photonRate rate_list[]);
void init_rates(struct photonRate * currRate, enum rate_type id); 
void openFileRead(bool binary, std::string filename, void ** pointer);
bool spacetimeRead(bool binary, void * file, float T_and_boosts[]);
bool viscRead(bool binary, void * shearFile, void * bulkFile, float visc_info[]);
void update_position_info(int line, struct phaseSpace_pos *curr_pos);
void pre_computeDescretizedSpectrum(struct phaseSpace_pos *curr_pos, float T_and_boosts[], float visc_info[], const struct photonRate rate_list[], double discSpectra[CONST_N_rates][CONST_NkT][CONST_Nrap][CONST_Nphi][3]);
void compute_observables(const struct photonRate rate_list[], double discSpectra[CONST_N_rates][CONST_NkT][CONST_Nrap][CONST_Nphi][3]);
void computeDescretizedSpectrum(struct phaseSpace_pos *curr_pos, float T_and_boosts[], float visc_info[], const struct photonRate rate_list[], double discSpectra[CONST_N_rates][CONST_NkT][CONST_Nrap][CONST_Nphi][3]);
void fill_grid(struct phaseSpace_pos *curr_pos, double kR, double T, double Akk, double bulk_pressure, double eps_plus_P, double cs2, const struct photonRate * currRate, double discSpectra[CONST_NkT][CONST_Nrap][CONST_Nphi][3]);
void compute_midrapidity_yield_and_vn(const struct photonRate currRate[], double discSpectra[CONST_N_rates][CONST_NkT][CONST_Nrap][CONST_Nphi][3]);

#endif
