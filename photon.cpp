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

/********* Inputs ********/
//Format of the input files
const bool binaryMode=0; //0 for binary, 1 for text
//Location of the spacetime grid file
const std::string stGridFile="stgrid.dat";

//Information about the spacetime grid
//Number of cells of the grid
const int cellNb_x=65;
const int cellNb_y=65;
const int cellNb_eta=65;
//Size of the cells
const double cellsize_X=10./(cellNb_x-1); //In fm
const double cellsize_Y=10./(cellNb_y-1); //In fm
const double cellsize_Eta=10./(cellNb_eta-1); //In units of rapidity
//Initial time tau_0
const double CONST_tau0=0.4;
const double deltaTau=0.001;

//Run with viscosity or not
const bool viscosity=0; //0 for thermal, 1 for anisotropic
const double shear_to_s=0.08;
//Location of the viscous files
const std::string viscosityFile="miaw";

//Discretization of photon spectrum
//kT
const double CONST_ktMin=0.2; //Minimum value for kT
const double CONST_ktMax=4.0; //Maximum value for kT
const int CONST_Nkt=10;  //Warning: delta kT=(ktMax-kTmin)/(Nkt-1) 
//const double kTdisc[3] = [0.2,4.0,0.2] //Discretization in kT: [kT min, kT max, delta kT]
const int CONST_Nphi=16;  //phi
//Rapidity
const double CONST_etaMin=-1.0; //Minimum value for eta
const double CONST_etaMax=1.0; //Maximum value for eta
const int CONST_Neta=3;  //Warning: delta eta=(etaMax-etamin)/(Nkt-1)

//Observables

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



//Main
int main() {


	//Forward declaration
	void photon_prod();

	//Grid information

	//Compute photon production
	photon_prod();


}


//Compute photon production
void photon_prod() {

	//Forward declaration of functions
	void openFileRead(bool binary, std::string filename, void ** pointer);
	bool spacetimeRead(bool binary, void * file, float T_and_boosts[]);
	void infer_position_info(int line, struct phaseSpace_pos *curr_pos);

	//Variables
	float T_and_boosts[5], shear_info[10];
	bool readRes; //Result of reading of the file
	std::FILE * stFile; //Spacetime file
	double tau;
	int line=1; //Spacetime position inferred from line number
	struct phaseSpace_pos curr_pos;
	//vector<double> discSpectra;
	//double discSpectra[neta][nphi][nkt][nSpec];

	//rate(kOverT)=rate_ideal(koverT)+\hat{k}_\alpha \hat{k}_\beta \Pi^{\alpha \beta}*factor_born_viscous(kOverT)
	//Frame of \Pi^{\alpha \beta}:

	//Initialise spacetime position
	curr_pos.tau=CONST_tau0;
	//Skip any time in CONST_tauList that is smaller than CONST_tau0
	for(curr_pos.iTauList=0;CONST_tauList[curr_pos.iTauList]<CONST_tau0;curr_pos.iTauList+=1);
	curr_pos.newiTau=0;

	//Open spacetime grid file
	openFileRead(binaryMode, stGridFile, (void **) &stFile);

	//Read the first line of the spacetime grid
	//readRes=spacetimeRead(binaryMode, stFile, &T, &qgp, &ux, &uy, &uz);
	readRes=spacetimeRead(binaryMode, stFile, &T_and_boosts[0]);
	//Loop over the rest of the file
	while (readRes) {

		//Compute (tau,x,y,eta) from line number
		infer_position_info(line,&curr_pos);

		//Do stuff
		//computeDescretizedSpectrum();
		printf("Temp=%f\n",T_and_boosts[0]);

		//Try to read the next line
		readRes=spacetimeRead(binaryMode, stFile, T_and_boosts);
		line+=1;
	//	readRes=spacetimeRead(binaryMode, stFile, &T, &qgp, &ux, &uy, &uz);

	}

	//Close spacetime grid file
	std::fclose(stFile);


	//Compute observables from the discretized photon spectra
	//compObservables();

}

/***** File reading stuff *****/
//Open file for reading
void openFileRead(bool binary, std::string filename, void ** pointer) {

	//If binary
	if (binary) {
		*pointer=std::fopen(filename.c_str(),"rb");
	}
	else {
		*pointer=std::fopen(filename.c_str(),"r");
	}

	//std::printf("Pointer=%p\n",*pointer);

	//Check if it opened correctly
	if (NULL==*pointer) {
		std::printf("Error: Could not open file \"%s\"",filename.c_str());
	}


}

//Read spacetime file 
//bool spacetimeRead(bool binary, void * file, float * T, float * qgp, float * ux, float * uy, float * uz) {
//bool spacetimeRead(bool binary, void * file, float * T_and_boosts) {
bool spacetimeRead(bool binary, void * file, float T_and_boosts[]) {

	int elemRead;

	//If binary
	if (binary) {
		elemRead=std::fread(&T_and_boosts,5*sizeof(float),5,(std::FILE *) file);
		//elemRead=std::fread(T,sizeof(float),1,(std::FILE *) file);
		//elemRead+=std::fread(qgp,sizeof(float),1,(std::FILE *) file);
		//elemRead+=std::fread(ux,sizeof(float),1,(std::FILE *) file);
		//elemRead+=std::fread(uy,sizeof(float),1,(std::FILE *) file);
		//elemRead+=std::fread(uz,sizeof(float),1,(std::FILE *) file);
	}
	else {
		elemRead=std::fscanf((std::FILE *) file, "%f %f %f %f %f", &T_and_boosts[0], &T_and_boosts[1], &T_and_boosts[2], &T_and_boosts[3], &T_and_boosts[4]);
	}

	//If fscanf couldn't read the five elements, it's the end of the file or there's a problem
	if (elemRead != 5) {
		return 0;
	}
	else {
		return 1;
	}

}

/***** Computation of the discretized spectrum *****/
void computeDescretizedSpectrum(bool viscosity, int * grid_pos, void * grid_info, float T_and_boosts[], float shear_info[], double * discSpectra) {

	//Local variables
	double T, qgpFrac, betax, betay, betaz, gamma;
	double eta, phi, kt;
	double kL, kLx, kLy, kLz;
	double kR, kkPiOverEta;
	double delEta, delPhi, delPt;
	double coshEta, sinhEta, cosPhi, sinPhi, invCoshEta;
	double stPosition[4];

	//Assign those to local variables for convenience
	T=T_and_boosts[0];
	qgpFrac=T_and_boosts[1];
	betax=T_and_boosts[2];
	betay=T_and_boosts[3];
	betaz=T_and_boosts[4];

	//Pre-compute gamma for efficiency
	gamma=1.0/sqrt(1+betax*betax+betay*betay+betaz*betaz);

	//Deltas used for the (uniform) discretization of the grid
	delEta=(CONST_etaMax-CONST_etaMin)/(CONST_Neta-1.0);
	delPhi=(2*M_PI)/(CONST_Nphi-1.0);
	delPt=(CONST_ktMax-CONST_ktMin)/(CONST_Nkt-1.0);

	//Compute (tau,x,y,eta) from line number
	//Check if


	//Loop over transverse momentum kT, azimuthal angle phi and rapidity eta
	//(note that there is no different here between the rapidity and the pseudorapidity, the photon being massless)
	//Loop over rapidity eta
	for(int ieta=0;ieta<CONST_Neta; ieta++) {

		eta=CONST_etaMin+ieta*delEta;	

		coshEta=cosh(eta);
		sinhEta=sinh(eta);
		invCoshEta=1.0/coshEta;

		//Loop over phi (uniform discretization - to be used with the trapezoidal method)
		for(int iphi=0;iphi<CONST_Nphi; iphi++) {

			phi=iphi*delPhi;	

			cosPhi=cos(phi);
			sinPhi=sin(phi);

			//Evaluate (A_L)_{alpha beta} \hat{k}_L^alpha \alpha{k}_L^beta
			if (viscosity) {
				//shear_info: Wtt,Wtx,Wty,Wtz,Wxx,Wxy,Wxz,Wyy,Wyz,Wzz
				//*shear_info+: 0   1   2   3   4   5   6   7   8   9
				//\hat{K}=(1,cos(phi)/cosh(eta),sin(phi)/cosh(eta),sinh(eta)/cosh(eta))
				//kkPiOverEta=A00 + 1/cosh(eta)*( 2*(A01*k1+A02*k2+A03*k3) + 1/cosheta*() )
				//kkPiOverEta=A00 + 1/cosh(eta)^2*(A11 cos(phi)^2+A22*sin(phi)^2+A33*sinh(eta)^2)+2/cosh(eta)*(A01*cos(phi)+A02*sin(phi)+A03*sinh(eta)+1/cosh(eta)*(A12*cos(phi)*sin(phi)+A13*cos(phi)*sinh(eta)+A23*sin(phi)*sinh(eta)))
				//kkPiOverEta=*(shear_info) + invCoshEta*invCoshEta*( *(shear_info+4)*cosPhi*cosPhi + *(shear_info+7)*sinPhi*sinPhi + *(shear_info+9)*sinhEta*sinhEta) + 2.0*invCoshEta*(*(shear_info+1)*cosPhi + *(shear_info+2)*sinPhi + *(shear_info+3)*sinhEta + invCoshEta*( *(shear_info+5)*cosPhi*sinPhi + *(shear_info+6)*cosPhi*sinhEta + *(shear_info+8)*sinPhi*sinhEta));
				kkPiOverEta=shear_info[0] + invCoshEta*invCoshEta*( shear_info[4]*cosPhi*cosPhi + shear_info[7]*sinPhi*sinPhi + shear_info[9]*sinhEta*sinhEta) + 2.0*invCoshEta*( shear_info[1]*cosPhi + shear_info[2]*sinPhi + shear_info[3]*sinhEta + invCoshEta*( shear_info[5]*cosPhi*sinPhi + shear_info[6]*cosPhi*sinhEta + shear_info[8]*sinPhi*sinhEta));
				
				//Akk=(*shear_info)+invCoshEta( (*shear_info+4)*cosPhi*cosPhi);
			}
			

			//Loop over kT
			for(int ikt=0;ikt<CONST_Nkt; ikt++) {

				kt=CONST_ktMin+ikt*delPt;	

				//Photon momentum in the lab frame
				//k=mT cosh(eta)=kT cosh(eta)
				kL=kt*coshEta;
				//kx=kT cos(phi)
				kLx=kt*cosPhi;
				//ky=kT sin(phi)
				kLy=kt*sinPhi;
				//kz=mT sinh(eta)=kT sinh(eta)
				kLz=kt*sinhEta;

				//kR.uR=kL.uL
				//k_rf=(k_L-\vec{u}/u0.\vec{k})/sqrt(1-u^2/u0^2)
				kR=gamma*(kL-betax*kLx-betay*kLy-betaz*kLz);

				//Our rate
				//dGamma(\vec{k}_L)=dGamma_0(k_rf)+(A_L)_{alpha beta} k_L^alpha k_L^beta Z(rf) 
				/*stPosition[0]=tau;
				stPosition[1]=ieta;
				stPosition[2]=iphi;
				stPosition[3]=ikt;*/
				//fill_grid(&stPostion, kR,Akk,discSpectra);	

				//

			}
		}	
	}

}

//Discretized spectra: array[times][Neta][Nphi][Npt][rates]
//Discretized spectra, version 2: array[times][Neta][Nphi][Npt][rates][value_and_remainder]
//stPos=[tau, ieta, iphi, ikt]
void fill_grid(struct phaseSpace_pos *curr_pos, double kR, double Akk, double * discSpectra) {

	//Compute rates
	//rate_dusling(double kR, double T, double kkPiOverEta double * res, double * remainder);

	//Store rates


}

void infer_position_info(int line, struct phaseSpace_pos *curr_pos) {

	//We don't need x/y/eta position, so we don't compute it
	//curr_pos->x=-1;
	//curr_pos->y=-1;
	//curr_pos->eta=-1;

	//There is (cellNb_x*cellNb_y*cellNb_eta) line per timestep
	if ((line-1)%cellNb_x*cellNb_y*cellNb_eta == 0) {
		//Update tau	
		curr_pos->tau+=deltaTau;

		printf("New tau=%d at line %i\n",curr_pos->tau,line);

		//Check if there's a change in iTauList
		if (curr_pos->tau > CONST_tauList[curr_pos->iTauList]) {
			curr_pos->iTauList+=1;
			curr_pos->newiTau=1;
		}
		else {
			curr_pos->newiTau=0;
		}
		
	}

}


/***** Rates *****/

//Pre-factor A(k)
double prefA(double kOverT, double T) {

	//Forward declaration of functions
	//inline double fermiDirac(double kOverT);
	double fermiDirac(double kOverT);

	//Nf, dF, 
	const double qCharge2[]={4.0/9.0,5.0/9.0,2.0/3.0,10.0/9.0,11.0/9.0,5.0/3.0};

	//Asymptotic thermal quark mass
	double mInf2=CONST_CF*CONST_gs*CONST_gs*T*T/4.0;

	//
	double res;

	res=2*CONST_alphaEM*CONST_Nc*qCharge2[CONST_Nf-1]*CONST_mInfOverT*CONST_mInfOverT*T*fermiDirac(kOverT)/kOverT;

}

//QGP ideal rate
double rate_qgp_born_ideal_val(double kOverT, double T) {

	return prefA(kOverT,T)/CONST_twoPiCubed*(log(1/CONST_mInfOverT)+C_hard(kOverT));

}

//QGP viscous rate
double factor_born_viscous(double kOverT) {



}

//AMY fits (from arXiv:1302.5970, see also arXiv:hep-ph/0111107)
//Valid for k/T>0.2
double C_hard(double kOverT) {

	return 1/2.0*log(2*kOverT)+0.041/kOverT-0.3615+1.01*exp(-1.35*kOverT);

}

//Valid for 0.2<k/T<50
double C_coll(double kOverT) {

	return sqrt(1+CONST_Nf/6.0)*(0.548*log(12.28+1/kOverT)/pow(kOverT,3./2.)+0.133*kOverT/sqrt(1+kOverT/16.27));

}


//
double rate_qgp_Born_viscous(double kROverT, double kkPiOverEta, double T, double * remainder) {

	//Limits of the fit
	const double kROverT_min=0.2;
	const double kROverT_max=50.0;

	double res=0.0, rem=0.0;

	//
	if (kROverT > kROverT_min) {
		//Inside validity range
		if (kROverT < kROverT_max) {
			//res=...
		}
		//Above validity range
		else {
			//rem=
		}
	}
	//Below validity range
	else {
		//rem=
	}

	//Return values
	*remainder=rem;
	return res;

}


//rate_dusling(double kR, double T, double kkPiOverEta double * res, double * remainder)
double rate_qgp_Dusling(double kR, double kkPiOverEta, double T) {

	//No limits of the fit, thus no remainder


}

//Kevin Dusling (arXiv:0903.1764)
double dusling_factor(double k, double T) {

	//E d^3 R/dk^3= k A(k)/(2 Pi)^3 1/2*log(3.7388 (kOverT)/g^2) (1+(1-nf(kOverT))*kOverT^2*eta/s*1/(2*T)*kkPOverEta)

}



//Fermi-Dirac
double fermiDirac(double kOverT) {

	return 1.0/(exp(-kOverT)+1.0);

}

//Bose-Einstien
double boseEin(double kOverT) {

	return 1.0/(exp(-kOverT)-1.0);

}

