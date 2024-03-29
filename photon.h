#ifndef PHOTON_H
#define PHOTON_H

#include <string>
#include <cmath>
#include <cstdio>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include "rates.h"

//Read parameters directly from a file generated by the hydro
#include "hydro_info_header_h"


// #########################################################
// ### Location and format of the hydrodynamic evolution ###
// #########################################################
//
//const bool CONST_binaryMode=MUSIC_outputBinaryEvolution; //0 for text, 1 for binary
enum evolution_file_format { old_format, new_format };
// Old format: "T, dummy, betax, betay, betaz", plus separate files for shear tensor and bulk pressure
// New format: "volume=dx*dy*deta*dtau*tau, eta, T, ux, uy, ueta, Wxx, Wxy, Wxeta, Wyy,Wyeta, pi_b"
//Location of the spacetime grid file
const bool CONST_binaryMode=false; //0 for text, 1 for binary

//const enum evolution_file_format CONST_file_format=new_format;
//const std::string stGridFile="./evolution_xyeta_eos_qcd_new_format.dat";

const enum evolution_file_format CONST_file_format=old_format;
const std::string stGridFile="./evolution_xyeta_eos_qcd_old_format.dat";

const std::string shearViscosityFile="evolution_Wmunu_over_epsilon_plus_P_xyeta.dat";
const std::string bulkViscosityFile="evolution_bulk_pressure_xyeta.dat";

// Viscous or not? (that is, should we read the files with the shear tensor and the bulk pressure)
const bool CONST_with_shear_viscosity=MUSIC_with_shear_viscosity; //turn on and off shear viscosity
const bool CONST_with_bulk_viscosity=MUSIC_with_bulk_viscosity; //turn on and off bulk viscosity
const bool CONST_with_viscosity=CONST_with_shear_viscosity||CONST_with_bulk_viscosity; //general flag for viscosity


// Information about the spacetime grid, used only for the old evolution file format
//Number of cells of the grid
const int cellNb_x=MUSIC_real_nx;
const int cellNb_y=MUSIC_real_ny;
const int cellNb_eta=MUSIC_real_neta;
// Size of the cells
const double CONST_cellsize_X=MUSIC_dx; //In fm
const double CONST_cellsize_Y=MUSIC_dy; //In fm
const double CONST_cellsize_Eta=MUSIC_deta; //In units of rapidity
const double CONST_effective_dTau=MUSIC_dtau;
// Initial time tau_0
const double CONST_tau0=MUSIC_tau0;

// #################################
// ### 2+1D hydro or 3+1D hydro? ###
// #################################

// Is the evolution file a single slice in eta_s of a 2+1D hydro?
const bool CONST_boost_invariant=1;
const double CONST_eta_s_of_saved_slice=0.0;
const int CONST_nb_steps_eta_integration=30;
const double CONST_max_eta_integration=3.0;

// ##################################################
// ### Stop producing photons at this temperature ###
// ##################################################

//const double CONST_freezeout_T=MUSIC_kinetic_FO_temperature_in_GeV;
const double CONST_freezeout_T=0.145;


// #################################################################################
// ### Discretization of the photon "p d^3N/(p_T dp_T dphi dy)" will be computed ###
// #################################################################################

//Discretization
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

//If more than one point if used in rapidity: Mid-rapidity cut
// the midrapidity result will be an average over approximatively -midRapCut to midRapCut
const double CONST_midRapCut = 0.5;

//Number of Fourier coefficient to compute
const int CONST_FourierNb = 6;



// #########################################################################
// ### Calculate "p d^3N/(p_T dp_T dphi dy)" for all the following rates ###
// #########################################################################

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






// #########################
// ### General constants ###
// #########################

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




// ###################################################################
// ### struct to pass information about local hydrodynamic profile ###
// ###################################################################

struct hydro_info_t {

	float tau;
	float V4; //in fm^4
	float eta_s; // spatial rapidity
	float T, muB; //In GeV
	float ux, uy, tau_ueta; //u^x, u^y, u^eta

        // Additional thermodynamic information
        float epsilon_plus_P, cs2;

	// Viscous part
	float pitautau_over_eps_plus_p, pitaux_over_eps_plus_p, pitauy_over_eps_plus_p, tau_pitaueta_over_eps_plus_p, pixx_over_eps_plus_p, pixy_over_eps_plus_p, tau_pixeta_over_eps_plus_p, piyy_over_eps_plus_p, tau_piyeta_over_eps_plus_p, tau_tau_pietaeta_over_eps_plus_p;
	float Pi_b;

};




// ###########################
// ### Forward declaration ###
// ###########################

void photon_prod(std::map<enum rate_type, struct photonRate> * rate_list);
bool open_file_read(bool binary, std::string filename, std::FILE ** pointer);
bool init_hydro_field_files(std::FILE * hydro_fields_files[3]);
void close_hydro_field_files(std::FILE * hydro_fields_files[3]);
bool read_hydro_fields(std::FILE * hydro_fields_files[3], struct hydro_info_t & hydro_info);
bool read_hydro_fields_new_format(std::FILE * hydro_fields_files[3], struct hydro_info_t & hydro_info);
bool read_hydro_fields_old_format(std::FILE * hydro_fields_files[3], struct hydro_info_t & hydro_info);
void pre_computeDescretizedSpectrum(struct hydro_info_t & hydro_info, std::map<enum rate_type, struct photonRate> * rate_list, double discSpectra[CONST_N_rates][CONST_NkT][CONST_Nrap][CONST_Nphi][3]);
void computeDescretizedSpectrum(struct hydro_info_t & hydro_info, std::map<enum rate_type, struct photonRate> * rate_list, double discSpectra[CONST_N_rates][CONST_NkT][CONST_Nrap][CONST_Nphi][3]);
void fill_grid(int irap, int iphi, int ikT, double kR, double T, double V4, double kOverTkOverTOver_e_P, double bulk_pressure, double eps_plus_P, double cs2, std::map<enum rate_type, struct photonRate> * rate_list, enum rate_type, double discSpectra[CONST_NkT][CONST_Nrap][CONST_Nphi][3]);
void compute_observables(std::map<enum rate_type, struct photonRate> * rate_list, double discSpectra[CONST_N_rates][CONST_NkT][CONST_Nrap][CONST_Nphi][3]);
void compute_midrapidity_yield_and_vn(std::map<enum rate_type, struct photonRate> * rate_list, double discSpectra[CONST_N_rates][CONST_NkT][CONST_Nrap][CONST_Nphi][3]);

#endif
