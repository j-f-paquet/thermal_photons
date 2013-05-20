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

/********* Inputs ********/
//Format of the input files
const bool CONST_binaryMode=0; //0 for text, 1 for binary
//Location of the spacetime grid file
const std::string stGridFile="./evolution_small.dat";

//Information about the spacetime grid
//Number of cells of the grid
const int cellNb_x=65;
const int cellNb_y=65;
const int cellNb_eta=64;
//Size of the cells
const double cellsize_X=10./(cellNb_x-1); //In fm
const double cellsize_Y=10./(cellNb_y-1); //In fm
const double cellsize_Eta=10./(cellNb_eta-1); //In units of rapidity
//Initial time tau_0
const double CONST_tau0=0.4;
const double deltaTau=0.001;

//Run with viscosity or not
const bool CONST_viscosity=0; //0 for thermal, 1 for anisotropic
const double shear_to_s=0.08;
//Location of the viscous files
const std::string viscosityFile="evolution_Wmunu_over_shear_xyeta.dat";

//Discretization of photon spectrum
//kT
const double CONST_ktMin=0.2; //Minimum value for kT
const double CONST_ktMax=4.0; //Maximum value for kT
const int CONST_Nkt=5;  //Warning: delta kT=(ktMax-kTmin)/(Nkt-1) 
//const double kTdisc[3] = [0.2,4.0,0.2] //Discretization in kT: [kT min, kT max, delta kT]
const int CONST_Nphi=16;  //phi
//Rapidity
const double CONST_etaMin=-1.0; //Minimum value for eta
const double CONST_etaMax=1.0; //Maximum value for eta
const int CONST_Neta=3;  //Warning: delta eta=(etaMax-etamin)/(Nkt-1)

//Deltas used for the (uniform) discretization of the grid
const double CONST_delEta=(CONST_etaMax-CONST_etaMin)/(CONST_Neta-1.0);
const double CONST_delPhi=(2*M_PI)/(CONST_Nphi-1.0);
const double CONST_delKt=(CONST_ktMax-CONST_ktMin)/(CONST_Nkt-1.0);

//Observables
//const std::vector<std::string> CONST_rateList = {"ideal","viscous","viscousDusling"};
//const int miaw[] = {1,2,3,4};
//const char char_rateList[4][100] = {"01", "02", "03", "04"};
//std::vector<std::string> v(char_rateList, char_rateList + 4);
const char char_rateList[][100] = {"rate_qgp_ideal_Born"};
//, "ideal_LL", "viscous_Dusling", "viscous_LL", "ideal_LL", "viscous_Dusling", "viscous_LL"};
const int CONST_N_rates=int(sizeof(char_rateList)/sizeof(char)/100.);
const std::vector<std::string> CONST_rateList(char_rateList, char_rateList + CONST_N_rates);

//Mid-rapidity cut: the midrapidity result will be an average over approximatively -midRapCut to midRapCut
const double CONST_midRapCut = 0.5;
//Number of Fourier coefficient to compute
const int CONST_FourierNb = 3;
	

//Generate spectra sums from t0 to t_i with t_i \in CONST_tauList
const double CONST_tauList[]={.6,1.0,2.0};
const double CONST_tempList[]={.180,.500};

//General constants
const int CONST_Nc=3;
const int CONST_Nf=3;
const double CONST_CF=(CONST_Nc*CONST_Nc-1.0)/(2.0*CONST_Nc);
const double CONST_alphaEM=1/137.0;
const double CONST_alphaS=0.3;
const double CONST_gs=sqrt(4*M_PI*CONST_alphaS);
const double CONST_mInfOverT=CONST_CF*CONST_gs/2.0;
const double CONST_twoPiCubed=pow(2*M_PI,3);

/*************************/

//Use to store spacetime position and related informations
struct phaseSpace_pos {

	//
	double tau, x, y, eta;

	//
	int ikt, ieta, iphi;

	//
	int iTauList;
	bool newiTau;

};

#endif
