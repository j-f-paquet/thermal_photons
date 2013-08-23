#include "photon.h"

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
	bool viscRead(bool binary, void * file, float visc_info[]);
	void update_position_info(int line, struct phaseSpace_pos *curr_pos);
	void pre_computeDescretizedSpectrum(bool viscosity, struct phaseSpace_pos *curr_pos, float T_and_boosts[], float visc_info[], double discSpectra[CONST_Neta][CONST_Nphi][CONST_Nkt][3][CONST_N_rates]);
	void compute_observables(double discSpectra[CONST_Neta][CONST_Nphi][CONST_Nkt][3][CONST_N_rates]);

	//Variables
	float T_and_boosts[5], visc_info[10];
	bool read_T_flag, read_visc_flag; //Result of reading of the file
	std::FILE * stFile; //Spacetime file
	std::FILE * viscFile; //Spacetime file
	//double tau;
	int line=1; //Spacetime position inferred from line number
	struct phaseSpace_pos curr_pos;
	//The second to last dimension is meant for including an upper and a lower bound on the uncertainty, if possible
	//discSpectra[][][][0/1/2][] is for the lower bound/value/upper bound
	double discSpectra[CONST_Neta][CONST_Nphi][CONST_Nkt][3][CONST_N_rates] = {0.0};

	//rate(kOverT)=rate_ideal(koverT)+\hat{k}_\alpha \hat{k}_\beta \Pi^{\alpha \beta}*factor_born_viscous(kOverT)
	//Frame of \Pi^{\alpha \beta}:

	//Initialise spacetime position
	curr_pos.tau=CONST_tau0;
	curr_pos.itau=1;
	curr_pos.ix=1;
	curr_pos.iy=1;
	curr_pos.ieta=1;

	//Open spacetime grid file
	openFileRead(CONST_binaryMode, stGridFile, (void **) &stFile);
	if (CONST_with_viscosity) openFileRead(CONST_binaryMode, viscosityFile, (void **) &viscFile);

	//Read the first line of the spacetime grid
	read_T_flag=spacetimeRead(CONST_binaryMode, stFile, T_and_boosts);
	if (CONST_with_viscosity) read_visc_flag=viscRead(CONST_binaryMode, viscFile, visc_info);
	//Loop over the rest of the file
	while ((read_T_flag)&&((!CONST_with_viscosity)||(read_visc_flag && CONST_with_viscosity))) {

		//std::cout << "Line " << line << "=Position (itau,ix,iy,ieta)=(" << curr_pos.itau << "," << curr_pos.ix << "," << curr_pos.iy << "," << curr_pos.ieta << ")\n";

		if (T_and_boosts[0]>=CONST_freezeout_T) {
			pre_computeDescretizedSpectrum(CONST_with_viscosity, &curr_pos, T_and_boosts, visc_info, discSpectra);
		}
		//Compute (tau,x,y,eta) from line number
		update_position_info(line,&curr_pos);

		//Try to read the next line
		read_T_flag=spacetimeRead(CONST_binaryMode, stFile, T_and_boosts);
		if (CONST_with_viscosity) read_visc_flag=viscRead(CONST_binaryMode, viscFile, visc_info);
		line+=1;

	}

	if ((!std::feof(stFile))||((CONST_with_viscosity)&&(!std::feof(viscFile)))) {
		std::cout << "!!!!!!!!!!!!!!! Warning !!!!!!!!!!!!!!!!!! Stopped reading the evolution files before the end of the file!\n";	
	}

	//Close spacetime grid file
	std::fclose(stFile);
	if (CONST_with_viscosity) std::fclose(viscFile);

	//Compute observables from the discretized photon spectra
	compute_observables(discSpectra);

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

//Read file with info about viscous hydro 
bool viscRead(bool binary, void * file, float visc_info[]) {

	int elemRead;

	//If binary
	if (binary) {
		elemRead=std::fread(&visc_info,10*sizeof(float),10,(std::FILE *) file);
	}
	else {
		elemRead=std::fscanf((std::FILE *) file, "%f %f %f %f %f %f %f %f %f %f", &visc_info[0], &visc_info[1], &visc_info[2], &visc_info[3], &visc_info[4], &visc_info[5], &visc_info[6], &visc_info[7], &visc_info[8], &visc_info[9]);
	}

	//If fscanf couldn't read the ten elements, it's the end of the file or there's a problem
	if (elemRead != 10) {
		return 0;
	}
	else {
		return 1;
	}

}


void pre_computeDescretizedSpectrum(bool viscosity, struct phaseSpace_pos *curr_pos, float T_and_boosts[], float visc_info[], double discSpectra[CONST_Neta][CONST_Nphi][CONST_Nkt][3][CONST_N_rates]) {

	//Forward declaration
	void computeDescretizedSpectrum(bool viscosity, struct phaseSpace_pos *curr_pos, float T_and_boosts[], float visc_info[], double discSpectra[CONST_Neta][CONST_Nphi][CONST_Nkt][3][CONST_N_rates]);

	if (!CONST_boost_invariant) {
		computeDescretizedSpectrum(viscosity, curr_pos, T_and_boosts, visc_info, discSpectra);
	}
	else {
		const int eta_slice_choice=cellNb_eta/2+1;
		//Integrate only over 1 slice in eta
		if (eta_slice_choice == curr_pos->ieta) {
			double eta, eta_slice;
			double old_u0, new_u0;
			double ux, uy;
			eta_slice=(-cellNb_eta/2+eta_slice_choice-1)*CONST_cellsize_Eta;
			old_u0=1.0/sqrt(1-T_and_boosts[2]*T_and_boosts[2]-T_and_boosts[3]*T_and_boosts[3]-T_and_boosts[4]*T_and_boosts[4]);
			ux=T_and_boosts[2]*old_u0;
			uy=T_and_boosts[3]*old_u0;

			//Integrate in eta
			for(int j=0;j<CONST_nb_steps_eta_integration;j++) {
				eta=-1*CONST_max_eta_integration+2*CONST_max_eta_integration*j/(CONST_nb_steps_eta_integration-1);
				new_u0=sqrt((1+ux*ux+uy*uy)/(1-pow(tanh(eta),2))); 
				T_and_boosts[2]=ux/new_u0;
				T_and_boosts[3]=uy/new_u0;
				T_and_boosts[4]=tanh(eta);

				curr_pos->eta=eta;

				computeDescretizedSpectrum(viscosity, curr_pos, T_and_boosts, visc_info, discSpectra);
			}
		}
	}

}


/***** Computation of the discretized spectrum *****/
void computeDescretizedSpectrum(bool viscosity, struct phaseSpace_pos *curr_pos, float T_and_boosts[], float visc_info[], double discSpectra[CONST_Neta][CONST_Nphi][CONST_Nkt][3][CONST_N_rates]) {
	
	//Forward declaration
	void fill_grid(struct phaseSpace_pos *curr_pos, double kR, double T, double Akk, double discSpectra[CONST_Neta][CONST_Nphi][CONST_Nkt][3][CONST_N_rates]);

	//Local variables
	double T, qgpFrac, betax, betay, betaz, gamma;
	double eta, phi, kt;
	double kL, kLx, kLy, kLz;
	double kR, kHatkHatPiOver_e_P;
	double coshEta, sinhEta, cosPhi, sinPhi, invCoshEta;
	//double stPosition[4];

	//Assign those to local variables for convenience
	T=T_and_boosts[0];
	qgpFrac=T_and_boosts[1];
	betax=T_and_boosts[2];
	betay=T_and_boosts[3];
	betaz=T_and_boosts[4];

	//Pre-compute gamma for efficiency
	gamma=1.0/sqrt(1-(betax*betax+betay*betay+betaz*betaz));

	//Compute (tau,x,y,eta) from line number

	//Loop over transverse momentum kT, azimuthal angle phi and rapidity eta
	//(note that there is no different here between the rapidity and the pseudorapidity, the photon being massless)
	//Loop over kT
	for(int ikt=0;ikt<CONST_Nkt; ikt++) {

		kt=CONST_ktMin+ikt*CONST_delKt;	

		//Loop over rapidity eta
		for(int ieta=0;ieta<CONST_Neta; ieta++) {

			eta=CONST_etaMin+ieta*CONST_delEta;	

			coshEta=cosh(eta);
			sinhEta=sinh(eta);
			invCoshEta=1.0/coshEta;

			//Loop over phi (uniform discretization - to be used with the trapezoidal method)
			for(int iphi=0;iphi<CONST_Nphi; iphi++) {

				phi=iphi*CONST_delPhi;	

				cosPhi=cos(phi);
				sinPhi=sin(phi);

				//Evaluate (A_L)_{alpha beta} \hat{k}_L^alpha \alpha{k}_L^beta
				if (viscosity) {
					//shear_info: Wtt,Wtx,Wty,Wtz,Wxx,Wxy,Wxz,Wyy,Wyz,Wzz
					//*shear_info+: 0   1   2   3   4   5   6   7   8   9
					//\hat{K}=(1,cos(phi)/cosh(eta),sin(phi)/cosh(eta),sinh(eta)/cosh(eta))
					if (!CONST_boost_invariant) {
						//kkPiOverEta=A00 + 1/cosh(eta)*( 2*(A01*k1+A02*k2+A03*k3) + 1/cosheta*() )
						//kkPiOverEta=A00 + 1/cosh(eta)^2*(A11 cos(phi)^2+A22*sin(phi)^2+A33*sinh(eta)^2)+2/cosh(eta)*(A01*cos(phi)+A02*sin(phi)+A03*sinh(eta)+1/cosh(eta)*(A12*cos(phi)*sin(phi)+A13*cos(phi)*sinh(eta)+A23*sin(phi)*sinh(eta)))
						//kkPiOverEta=*(shear_info) + invCoshEta*invCoshEta*( *(shear_info+4)*cosPhi*cosPhi + *(shear_info+7)*sinPhi*sinPhi + *(shear_info+9)*sinhEta*sinhEta) + 2.0*invCoshEta*(*(shear_info+1)*cosPhi + *(shear_info+2)*sinPhi + *(shear_info+3)*sinhEta + invCoshEta*( *(shear_info+5)*cosPhi*sinPhi + *(shear_info+6)*cosPhi*sinhEta + *(shear_info+8)*sinPhi*sinhEta));
						kHatkHatPiOver_e_P=visc_info[0] + invCoshEta*invCoshEta*( visc_info[4]*cosPhi*cosPhi + visc_info[7]*sinPhi*sinPhi + visc_info[9]*sinhEta*sinhEta) + 2.0*invCoshEta*( -1.0*visc_info[1]*cosPhi - visc_info[2]*sinPhi - visc_info[3]*sinhEta + invCoshEta*( visc_info[5]*cosPhi*sinPhi + visc_info[6]*cosPhi*sinhEta + visc_info[8]*sinPhi*sinhEta));
					}
					//In the boost-invariant case, \Pi^\mu\nu is the value a eta=0
					//k_\mu k\nu \Pi^\mu\nu must be calculed correctly
					else {
						//K=(kt cosh(y), kt cos(phi), kt sin(phi), kt sinh(y))
						//k=kt*cosh(eta)
						//K boosted from rapidity y to rapidity eta:
						//K=(kt cosh(y-eta), kt cos(phi), kt sin(phi), kt sinh(y-eta))
						//k boosted = kt*cosh(y-eta)
						//\hat{K}_boosted=(1,cos(phi)/cosh(y-eta),sin(phi)/cosh(y-eta),sinh(y-eta)/cosh(y-eta))
						//
						double boostedCoshEta=cosh(eta-curr_pos->eta);
						double boostedSinhEta=sinh(eta-curr_pos->eta);
						double boostedInvCoshEta=1.0/boostedCoshEta;
						kHatkHatPiOver_e_P=visc_info[0] + boostedInvCoshEta*boostedInvCoshEta*( visc_info[4]*cosPhi*cosPhi + visc_info[7]*sinPhi*sinPhi + visc_info[9]*boostedSinhEta*boostedSinhEta) + 2.0*boostedInvCoshEta*( -1.0*visc_info[1]*cosPhi - visc_info[2]*sinPhi - visc_info[3]*boostedSinhEta + boostedInvCoshEta*( visc_info[5]*cosPhi*sinPhi + visc_info[6]*cosPhi*boostedSinhEta + visc_info[8]*sinPhi*boostedSinhEta));
					}


					//double tr_check=(visc_info[0]-visc_info[4]-visc_info[7]-visc_info[9]);
					//if (tr_check > 1e-5) {
					//	std::cout << "Warning: Large deviation from tracelessness (" << tr_check << ")!\n";
					//}
					
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
				curr_pos->ikt=ikt;
				curr_pos->iphi=iphi;
				curr_pos->irap=ieta;
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
	double QGP_fraction(double T);
	void get_photon_rate(int selector, double (**local_rate)(double, double, double));
	double get_photon_rate_accel(double kOverT, double T, double kk, int rate_no);

	//
	int ieta=curr_pos->irap;
	int iphi=curr_pos->iphi;
	int ikt=curr_pos->ikt;
	double tmpRate;
	double (*local_rate)(double, double, double);
	double tmp_cellsize_eta;

	//Size of cell in eta different if integrating 2+1D or 3+1D
	if (CONST_boost_invariant) {
		tmp_cellsize_eta=2*CONST_max_eta_integration/CONST_nb_steps_eta_integration;
	}
	else {
		tmp_cellsize_eta=CONST_cellsize_Eta;
	}

	//Loop over rates
	for(int iRate=0; iRate<CONST_N_rates;iRate++) {
	//double discSpectra[CONST_Neta][CONST_Nphi][CONST_Nkt][3][CONST_rateList.size()];

		get_photon_rate(CONST_rates_to_use[iRate], &local_rate);

		//Either use the fits directly or a table made from the fits
		if (CONST_use_accel_rates[CONST_rates_to_use[iRate]]) {
			tmpRate=get_photon_rate_accel(kR/T, T, kHatkHatPiOver_e_P, iRate);
		}
		else {
			tmpRate=(*local_rate)(kR/T,T,kHatkHatPiOver_e_P);
		}
		//tmpRate=kR*(1.0+cos(kR*(2.0)*iphi*CONST_delPhi));

		//QGP fraction
		//tmpRate*=QGP_fraction(T);

		//Cell volume: dx*dy*dz*dt=dx*dy*dEta*dTau*tau
		tmpRate*=CONST_cellsize_X*CONST_cellsize_Y*tmp_cellsize_eta*CONST_effective_dTau*curr_pos->tau;

		//Fill value
		discSpectra[ieta][iphi][ikt][1][iRate]+=tmpRate;

		//Fill lower bound uncertainty
		discSpectra[ieta][iphi][ikt][0][iRate]=0.0;

		//Fill upper bound uncertainty
		discSpectra[ieta][iphi][ikt][2][iRate]=0.0;

	}


}

void update_position_info(int line, struct phaseSpace_pos *curr_pos) {

	//We don't need x/y/eta position, so we don't compute it
	//curr_pos->x=-1;
	//curr_pos->y=-1;
	//curr_pos->eta=-1;

	//There is (cellNb_x*cellNb_y*cellNb_eta) line per timestep
	if ((line)%(cellNb_x*cellNb_y*cellNb_eta) == 0) {

		std::cout << "Done with time slice " << curr_pos->tau << " fm^-1\n";

		//Update tau	
		curr_pos->tau+=CONST_effective_dTau;
		curr_pos->itau+=1;
		curr_pos->ix=1;
		curr_pos->iy=1;
		curr_pos->ieta=1;

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
	else if ((line)%(cellNb_x*cellNb_y)==0) {
		curr_pos->ieta+=1;
		curr_pos->ix=1;
		curr_pos->iy=1;
	}
	else if ((line)%(cellNb_x)==0) {
		curr_pos->iy+=1;
		curr_pos->ix=1;
	}
	else {
		curr_pos->ix+=1;
	}

}


void compute_observables(double discSpectra[CONST_Neta][CONST_Nphi][CONST_Nkt][3][CONST_N_rates]) {

	//
	void compute_midrapidity_yield_and_vn(int rate_id, double discSpectra[CONST_Neta][CONST_Nphi][CONST_Nkt][3][CONST_N_rates]);

	//One file per rate
	for(int i=0; i<CONST_N_rates; i++) {

		compute_midrapidity_yield_and_vn(i, discSpectra);

	}

}

//Output the phi-integrated, rapidity-averaged-around-0 yield as a function of pT
void compute_midrapidity_yield_and_vn(int rate_no, double discSpectra[CONST_Neta][CONST_Nphi][CONST_Nkt][3][CONST_N_rates]) {

	double kt, eta, phi, yFac, rap_interval;
	double yield, vn[CONST_FourierNb];
	int iEtamin=0, iEtamax;
	bool exact_midrap=false, bad_rap_discret = false;
	std::string rate_name = CONST_available_rate[CONST_rates_to_use[rate_no]-1];

	//Set output file name
	std::stringstream tmpStr;
	tmpStr << "vn_";
	tmpStr << rate_name.c_str(); //CONST_rate_names[rate_id];
	tmpStr << ".dat";

	//Open output file
	std::ofstream outfile;
	outfile.open(tmpStr.str().c_str());
	//Set the format of the output
	//outfile.width (10);
	outfile.precision(10);
	outfile.setf(std::ios::scientific);

	//Output result
	outfile << "#pt" << "\t" << "yield";
	for(int j=1;j<=CONST_FourierNb; j++) {
		outfile	<< "\tyield*vn[" << j << "]\tvn[" << j << "]";
	}
	outfile<< "\n";
		

	//Identify the cells in rapidity that should be averaged over
	for(int ieta=0;ieta<CONST_Neta; ieta++) {

		eta=CONST_etaMin+ieta*CONST_delEta;	

		if (fabs(eta) <= CONST_midRapCut) {
			iEtamin=ieta;
			iEtamax=ieta;
		}
		else if (fabs(eta) > CONST_midRapCut) {
			//iEtamax=ieta-1;
			continue;
		}
		else iEtamax=ieta;
	
	}

	//If there is a single point and it is eta=0.0, 
	if (iEtamin == iEtamax) {
		if (0.0 == CONST_etaMin+iEtamin*CONST_delEta) {
			exact_midrap=true;
			rap_interval=1.0;
		}
		else {
			bad_rap_discret=true;
		}
	}
	else {
		rap_interval=(iEtamax-iEtamin)*CONST_delEta;
	}

	if (bad_rap_discret) {
			outfile << "Can't evaluate midrapidity spectra with current rapidity discretization\n";
	}
	else {

		for(int ikt=0;ikt<CONST_Nkt; ikt++) {

			//Will contain the results of the phi integration and rapidity averaging
			yield=0.0;
			for(int i=1;i<=CONST_FourierNb; i++) {
				vn[i]=0;
			}

			kt=CONST_ktMin+ikt*CONST_delKt;	

			//Loop over rapidity eta
			for(int ieta=iEtamin;ieta<=iEtamax; ieta++) {

				eta=CONST_etaMin+ieta*CONST_delEta;	

				//Loop over phi (trapezoidal method)
				for(int iphi=0;iphi<CONST_Nphi-1; iphi++) {

					phi=iphi*CONST_delPhi;	
					
					if (exact_midrap) {
						yFac=1.0;
					}
					//Let's use a simple midpoint rule for now
					else {
						yFac=CONST_delEta;
					}

					//Finally, multiply by dNdydptdphi[NY][NPT][NPHI+1] 
					//tmpIntRes+=phiFac*yFac*particleList[j].dNdydptdphi[iy][ipt][iphi];
					yield+=discSpectra[ieta][iphi][ikt][1][rate_no]*yFac;
					for(int i=1;i<=CONST_FourierNb; i++) {
						vn[i]+=yFac*discSpectra[ieta][iphi][ikt][1][rate_no]*cos(i*phi);
					}
					
				}
				
			}

			//Multiply by delta_ph and divide by the rapidity integration range 
			//(to yield an average instead of an integral)
			yield*=CONST_delPhi/rap_interval/(2.0*M_PI);
			for(int j=1;j<=CONST_FourierNb; j++) {
				vn[j]*=CONST_delPhi/rap_interval/(2.0*M_PI);
			}

			//Output result
			outfile << kt << "\t" << yield;
			for(int j=1;j<=CONST_FourierNb; j++) {
				outfile	<< "\t" << vn[j] << "\t" << vn[j]/yield;
			}
			outfile<< "\n";


		}

	}

	//Close file
	outfile.close();

}
