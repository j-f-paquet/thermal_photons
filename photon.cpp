#include "photon.h"

//Main
int main() {


	//Forward declaration
	void photon_prod();

	//Grid information

	//Compute photon production
	photon_prod();

	std::cout << "Number of rates=" << CONST_rateList.size() << "\n";
	for(int i=0; i<CONST_rateList.size();i++) {
		std::cout << "rate " << i << ":" << CONST_rateList[i] << "\n";
	}


}


//Compute photon production
void photon_prod() {

	//Forward declaration of functions
	void openFileRead(bool binary, std::string filename, void ** pointer);
	bool spacetimeRead(bool binary, void * file, float T_and_boosts[]);
	void infer_position_info(int line, struct phaseSpace_pos *curr_pos);
	void computeDescretizedSpectrum(bool viscosity, struct phaseSpace_pos *curr_pos, float T_and_boosts[], float shear_info[], double discSpectra[CONST_Neta][CONST_Nphi][CONST_Nkt][3][CONST_N_rates]);

	//Variables
	float T_and_boosts[5], shear_info[10];
	bool readRes; //Result of reading of the file
	std::FILE * stFile; //Spacetime file
	double tau;
	int line=1; //Spacetime position inferred from line number
	struct phaseSpace_pos curr_pos;
	//The second to last dimension is meant for including an upper and a lower bound on the uncertainty, if possible
	//discSpectra[][][][0][] is for the lower bound
	//discSpectra[][][][1][] is for the value bound
	//discSpectra[][][][2][] is for the upper bound
	double discSpectra[CONST_Neta][CONST_Nphi][CONST_Nkt][3][CONST_N_rates] = {0.0};

	//rate(kOverT)=rate_ideal(koverT)+\hat{k}_\alpha \hat{k}_\beta \Pi^{\alpha \beta}*factor_born_viscous(kOverT)
	//Frame of \Pi^{\alpha \beta}:

	//Initialise spacetime position
	curr_pos.tau=CONST_tau0;
	//Skip any time in CONST_tauList that is smaller than CONST_tau0
	for(curr_pos.iTauList=0;CONST_tauList[curr_pos.iTauList]<CONST_tau0;curr_pos.iTauList+=1);
	curr_pos.newiTau=0;

	//Open spacetime grid file
	openFileRead(CONST_binaryMode, stGridFile, (void **) &stFile);


	//Read the first line of the spacetime grid
	//readRes=spacetimeRead(binaryMode, stFile, &T, &qgp, &ux, &uy, &uz);
	readRes=spacetimeRead(CONST_binaryMode, stFile, &T_and_boosts[0]);
	//Loop over the rest of the file
	while (readRes) {
		//Compute (tau,x,y,eta) from line number
		infer_position_info(line,&curr_pos);

		//Do stuff
		//printf("Temp=%f\n",T_and_boosts[0]);

		//Try to read the next line
		readRes=spacetimeRead(CONST_binaryMode, stFile, T_and_boosts);
		line+=1;
	//	readRes=spacetimeRead(binaryMode, stFile, &T, &qgp, &ux, &uy, &uz);
		computeDescretizedSpectrum(CONST_viscosity, &curr_pos, T_and_boosts, 0, discSpectra);

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
bool spacetimeRead(bool binary, void * file, float T_and_boosts[]) {

	int elemRead;

	//If binary
	if (binary) {
		elemRead=std::fread(&T_and_boosts,5*sizeof(float),5,(std::FILE *) file);
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

//Read shear file 
bool shearRead(bool binary, void * file, float shear_info[]) {

	int elemRead;

	//If binary
	if (binary) {
		elemRead=std::fread(&shear_info,10*sizeof(float),10,(std::FILE *) file);
	}
	else {
		elemRead=std::fscanf((std::FILE *) file, "%f %f %f %f %f %f %f %f %f %f", &shear_info[0], &shear_info[1], &shear_info[2], &shear_info[3], &shear_info[4], &shear_info[5], &shear_info[6], &shear_info[7], &shear_info[8], &shear_info[9]);
	}

	//If fscanf couldn't read the five elements, it's the end of the file or there's a problem
	if (elemRead != 10) {
		return 0;
	}
	else {
		return 1;
	}

}


/***** Computation of the discretized spectrum *****/
void computeDescretizedSpectrum(bool viscosity, struct phaseSpace_pos *curr_pos, float T_and_boosts[], float shear_info[], double discSpectra[CONST_Neta][CONST_Nphi][CONST_Nkt][3][CONST_N_rates]) {
	
	//Forward declaration
	//void fill_grid(struct phaseSpace_pos *curr_pos, double kR, double Akk, double * discSpectra);
	void fill_grid(struct phaseSpace_pos *curr_pos, double kR, double T, double Akk, double discSpectra[CONST_Neta][CONST_Nphi][CONST_Nkt][3][CONST_N_rates]);

	//Local variables
	double T, qgpFrac, betax, betay, betaz, gamma;
	double eta, phi, kt;
	double kL, kLx, kLy, kLz;
	double kR, kHatkHatPiOver_e_P;
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
	//Loop over kT
	for(int ikt=0;ikt<CONST_Nkt; ikt++) {

		kt=CONST_ktMin+ikt*delPt;	

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
					kHatkHatPiOver_e_P=shear_info[0] + invCoshEta*invCoshEta*( shear_info[4]*cosPhi*cosPhi + shear_info[7]*sinPhi*sinPhi + shear_info[9]*sinhEta*sinhEta) + 2.0*invCoshEta*( shear_info[1]*cosPhi + shear_info[2]*sinPhi + shear_info[3]*sinhEta + invCoshEta*( shear_info[5]*cosPhi*sinPhi + shear_info[6]*cosPhi*sinhEta + shear_info[8]*sinPhi*sinhEta));
					
					//Akk=(*shear_info)+invCoshEta( (*shear_info+4)*cosPhi*cosPhi);
				}
				else {
					kHatkHatPiOver_e_P=0.0;
				}
			
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
				curr_pos->ikt=ikt;
				curr_pos->iphi=iphi;
				curr_pos->ieta=ieta;
				fill_grid(curr_pos, kR, T, kHatkHatPiOver_e_P,discSpectra);	


			}
		}	
	}

}

//Discretized spectra: array[times][Neta][Nphi][Npt][rates]
//Discretized spectra, version 2: array[times][Neta][Nphi][Npt][rates][value_and_remainder]
//stPos=[tau, ieta, iphi, ikt]
void fill_grid(struct phaseSpace_pos *curr_pos, double kR, double T, double kHatkHatPiOver_e_P, double discSpectra[CONST_Neta][CONST_Nphi][CONST_Nkt][3][CONST_N_rates]) {

	//
	double rate_qgp_ideal_born(double, double, double);

	//
	int ieta=curr_pos->ieta;
	int iphi=curr_pos->iphi;
	int ikt=curr_pos->ikt;
	double tmpRate;

	//Loop over rates
	for(int iRate=0; iRate<CONST_rateList.size();iRate++) {
	//double discSpectra[CONST_Neta][CONST_Nphi][CONST_Nkt][3][CONST_rateList.size()];

		//double (*local_rate)(double, double, double) = CONST_rateList[iRate].c_str();
		double (*local_rate)(double, double, double) = rate_qgp_ideal_born;


		//tmpRate=CONST_rateList[iRate].c_str()(0.0,0.0,0.0);
		tmpRate=(*local_rate)(kR,T,kHatkHatPiOver_e_P);
		
		//Fill value
		discSpectra[ieta][iphi][ikt][1][iRate]+=tmpRate;

		//Fill lower bound uncertainty
		discSpectra[ieta][iphi][ikt][0][iRate]=0.0;

		//Fill upper bound uncertainty
		discSpectra[ieta][iphi][ikt][2][iRate]=0.0;

	}

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

		/*
		printf("New tau=%d at line %i\n",curr_pos->tau,line);

		//Check if there's a change in iTauList
		if (curr_pos->tau > CONST_tauList[curr_pos->iTauList]) {
			curr_pos->iTauList+=1;
			curr_pos->newiTau=1;
		}
		else {
			curr_pos->newiTau=0;
		}
		*/
	}

}
