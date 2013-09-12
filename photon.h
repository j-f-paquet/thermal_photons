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
const double CONST_effective_dTau=MUSIC_effective_dtau;

const bool CONST_boost_invariant=1;
const int CONST_nb_steps_eta_integration=30;
const double CONST_max_eta_integration=3.0;

//Run with viscosity or not
const bool CONST_with_shear_viscosity=MUSIC_with_shear_viscosity; //0 for thermal, 1 for anisotropic
const bool CONST_with_viscosity=CONST_with_shear_viscosity;
//const double shear_to_s=0.08;
//Location of the viscous files
const std::string viscosityFile="evolution_Wmunu_over_epsilon_plus_P_xyeta.dat";

//Discretization of photon spectrum
//kT
const double CONST_ktMin=0.2; //Minimum value for kT
const double CONST_ktMax=4.0; //Maximum value for kT
const int CONST_Nkt=20;  //Warning: delta kT=(ktMax-kTmin)/(Nkt-1) 
//const double kTdisc[3] = [0.2,4.0,0.2] //Discretization in kT: [kT min, kT max, delta kT]
const int CONST_Nphi=16;  //phi
//Rapidity
const double CONST_etaMin=-0.0; //Minimum value for eta
const double CONST_etaMax=0.0; //Maximum value for eta
const int CONST_Neta=1;  //Warning: delta eta=(etaMax-etamin)/(Nkt-1)

//Deltas used for the (uniform) discretization of the grid
const double CONST_delEta= CONST_Neta > 1 ? (CONST_etaMax-CONST_etaMin)/(CONST_Neta-1.0) : 0;
const double CONST_delPhi=(2*M_PI)/(CONST_Nphi-1.0);
const double CONST_delKt=(CONST_ktMax-CONST_ktMin)/(CONST_Nkt-1.0);

//Observables
const char CONST_available_rate[][100]={"rate_qgp_ideal_born_AMYfit","rate_qgp_ideal_born_KLS","rate_qgp_ideal_born_JF_sqrtg","rate_qgp_viscous_only_born_JF_sqrtg", "rate_hg_ideal_Turbide_fit","rate_qgp_ideal_LO_AMYfit"};
/*
Rates:
1: double rate_qgp_ideal_born_AMYfit(double kOverT, double T, double kkPiOver_e_P_k2);
2: double rate_qgp_ideal_born_KLS(double kOverT, double T, double kkPiOver_e_P_k2);
3: double rate_qgp_ideal_born_JF_sqrtg(double kOverT, double T, double kkPiOver_e_P_k2);
4: double rate_qgp_viscous_only_born_JF_sqrtg(double kOverT, double T, double kkPiOver_e_P_k2);
5: rate_hg_ideal_Turbide_fit
6: rate_qgp_ideal_LO_AMYfit
*/
const int CONST_rates_to_use[] = {1,2,5,6};
//const int CONST_rates_to_use[] = {1,2,3,4,5,6};
const int CONST_N_rates = sizeof(CONST_rates_to_use)/sizeof(int);
//

//Mid-rapidity cut: the midrapidity result will be an average over approximatively -midRapCut to midRapCut
const double CONST_midRapCut = 0.5;
//Number of Fourier coefficient to compute
const int CONST_FourierNb = 2;

//QGP fraction definition
const double CONST_pure_QGP_T=0.22;	
const double CONST_pure_HG_T=0.184;	
const double CONST_freezeout_T=MUSIC_kinetic_FO_temperature_in_GeV;

//Generate spectra sums from t0 to t_i with t_i \in CONST_tauList
//const double CONST_tauList[]={.6,1.0,2.0};
//const double CONST_tempList[]={.180,.500};

//General constants
const int CONST_Nc=3;
const int CONST_Nf=3;
const double CONST_CF=(CONST_Nc*CONST_Nc-1.0)/(2.0*CONST_Nc);
const double CONST_alphaEM=1/137.0;
const double CONST_alphaS=0.3;
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
	int ikt, irap, iphi;

	//
	//int iTauList;
	//bool newiTau;

	//Spacetime integration weights
	double w_eta;

};

struct photonRate {

	//Name used when saving the results in a file
	std::string name;

	//Set to 1 to read a table from a file, instead of using a hard-coded fit
	bool use_table_instead_of_fit;
	std::string filename_of_external_table;

	//Function returning the hard-coded fit
	double (*rate_fit_function)(double, double, double);

	//Use fit, but tabulate it internally for speed
	bool tabulate_fit_for_speed;

	double ** tabulated_rate;

	//Parameters to specify by what factor should the rate be multiplied by, if any
	bool is_qgp; //Rate is multiplied by QGP_fraction
	bool is_hg; //Rate is multiplied by (1-QGP_fraction)

	bool is_shear_viscous; //Multiply rate by K_mu K_nu Pi^\mu\mu/k^2

	//Parameters used to specifiy either the table is...
	bool use_k_instead_of_kOverT;
	int number_of_points_in_kOverT;
	int number_of_points_in_temp;
	double min_temp, max_temp;
	double min_kOverT, max_kOverT;

	//Function to use to find the nearest tabulated kOverT
	double (**index_from_kOverT)(double, double, double);
	double (**kOverT_from_index)(double, double, double);

	//Function to use to find the nearest tabulated temperature
	double (**index_from_temp)(double, double, double);
	double (**temp_from_index)(double, double, double);


	photonRate() {

		//Set some default values
		rate_fit_function=0;
		is_qgp=false;
		is_hg=false;

		is_shear_viscous=false;

	}

};

//const bool CONST_use_accel_rates[] = {0,0,0,0,1,0};
//const int accel_table_sample_x[] = {500,500,500,500,500,500};
//const int accel_table_sample_y[] = {500,500,500,500,250,500};
//const double accel_table_min_temperature[]={CONST_pure_HG_T,CONST_pure_HG_T,CONST_pure_HG_T,CONST_pure_HG_T,CONST_freezeout_T,CONST_pure_HG_T};
//const double accel_table_max_temperature[]={2.0,2.0,2.0,2.0,CONST_pure_QGP_T,2.0};
//struct rate_accel {
//
//	double *** tabulated_rates;
//
//	rate_accel() {
//
//		void get_photon_rate(int selector, double (**local_rate)(double, double, double));
//		double kOverT_from_index(int i, int size);
//		double temp_from_index(int i, int size, int rate_no);
//
//		tabulated_rates = new double ** [CONST_N_rates];
//
//		int tmp_sample_x,tmp_sample_y;
//		double (*local_rate)(double, double, double);
//
//		for(int rate_no=0;rate_no<CONST_N_rates;rate_no++) {
//			if (!CONST_use_accel_rates[CONST_rates_to_use[rate_no]-1]) continue;
//			get_photon_rate(CONST_rates_to_use[rate_no], &local_rate);
//			tmp_sample_x=accel_table_sample_x[CONST_rates_to_use[rate_no]-1];
//			tmp_sample_y=accel_table_sample_y[CONST_rates_to_use[rate_no]-1];
//			tabulated_rates[rate_no] = new double * [tmp_sample_y];	
//			for(int k=0; k<tmp_sample_y;k++) {
//				tabulated_rates[rate_no][k] = new double [tmp_sample_x];	
//				for(int j=0; j<tmp_sample_x;j++) {
//					tabulated_rates[rate_no][k][j]=(*local_rate)(kOverT_from_index(j,tmp_sample_x),temp_from_index(k,tmp_sample_y,rate_no),0.0);	
//				}
//			}
//		}
//
//	}
//
//};
//
//const struct rate_accel CONST_rate_tables;

#endif
