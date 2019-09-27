#include "photon.h"

//Main
int main() {


	//Forward declaration
	void photon_prod(const struct photonRate rate_list[]);
	void init_rates(struct photonRate * currRate, enum rate_type id); 

	//
	struct photonRate rate_list[CONST_N_rates];

	//Grid information

	//Initialise rates
	for(int i=0;i<CONST_N_rates;i++) init_rates(&rate_list[i], CONST_rates_to_use[i]);


	//Compute photon production
	photon_prod(rate_list);
}


//Compute photon production
void photon_prod(const struct photonRate rate_list[]) {

	//Forward declaration of functions
	void openFileRead(bool binary, std::string filename, void ** pointer);
	bool spacetimeRead(bool binary, void * file, float T_and_boosts[]);
	bool viscRead(bool binary, void * shearFile, void * bulkFile, float visc_info[]);
	void update_position_info(int line, struct phaseSpace_pos *curr_pos);
	void pre_computeDescretizedSpectrum(struct phaseSpace_pos *curr_pos, float T_and_boosts[], float visc_info[], const struct photonRate rate_list[], double discSpectra[CONST_Nrap][CONST_Nphi][CONST_NkT][3][CONST_N_rates]);
	void compute_observables(const struct photonRate rate_list[], double discSpectra[CONST_Nrap][CONST_Nphi][CONST_NkT][3][CONST_N_rates]);

	//Variables
	float T_and_boosts[5], visc_info[13];
	bool read_T_flag, read_visc_flag=true; //Result of reading of the file
	std::FILE * stFile; //Spacetime file
	std::FILE * shearViscFile; //Shear viscosity file
	std::FILE * bulkViscFile; //Bulk viscosity file
	//double tau;
	int line=1; //Spacetime position inferred from line number
	struct phaseSpace_pos curr_pos;
	//The second to last dimension is meant for including an upper and a lower bound on the uncertainty, if possible
	//discSpectra[][][][0/1/2][] is for the lower bound/value/upper bound
	double discSpectra[CONST_Nrap][CONST_Nphi][CONST_NkT][3][CONST_N_rates] = {0.0};

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
	if (CONST_with_viscosity) {
		if (CONST_with_shear_viscosity) openFileRead(CONST_binaryMode, shearViscosityFile, (void **) &shearViscFile);
		if (CONST_with_bulk_viscosity) openFileRead(CONST_binaryMode, bulkViscosityFile, (void **) &bulkViscFile);
	}

	//Read the first line of the spacetime grid
	read_T_flag=spacetimeRead(CONST_binaryMode, stFile, T_and_boosts);
	if (CONST_with_viscosity) {
		read_visc_flag=viscRead(CONST_binaryMode, shearViscFile, bulkViscFile, visc_info);
	}
	//Loop over the rest of the file
	while ((read_T_flag)&&(read_visc_flag)) {

		//std::cout << "Line " << line << "=Position (itau,ix,iy,ieta)=(" << curr_pos.itau << "," << curr_pos.ix << "," << curr_pos.iy << "," << curr_pos.ieta << ")\n";
		
		if (T_and_boosts[0]>=CONST_freezeout_T) {
			pre_computeDescretizedSpectrum(&curr_pos, T_and_boosts, visc_info, rate_list, discSpectra);
		}
		//Compute (tau,x,y,eta) from line number
		update_position_info(line,&curr_pos);

		//Try to read the next line
		read_T_flag=spacetimeRead(CONST_binaryMode, stFile, T_and_boosts);
		if (CONST_with_viscosity) read_visc_flag=viscRead(CONST_binaryMode, shearViscFile, bulkViscFile, visc_info);
		line+=1;

	}

	if ((!std::feof(stFile))||((CONST_with_viscosity)&&(!std::feof(shearViscFile)))) {
		std::cout << "!!!!!!!!!!!!!!! Warning !!!!!!!!!!!!!!!!!! Stopped reading the evolution files before the end of the file!\n";	
	}

	//Close spacetime grid file
	std::fclose(stFile);
	if (CONST_with_viscosity) {
		if (CONST_with_shear_viscosity) std::fclose(shearViscFile);
		if (CONST_with_bulk_viscosity) std::fclose(bulkViscFile);
	}

	//Compute observables from the discretized photon spectra
	compute_observables(rate_list, discSpectra);

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
		float T,qgpFrac,vx,vy,vz;
		elemRead=std::fread(&T,sizeof(float),1,(std::FILE *) file);
		elemRead+=std::fread(&qgpFrac,sizeof(float),1,(std::FILE *) file);
		elemRead+=std::fread(&vx,sizeof(float),1,(std::FILE *) file);
		elemRead+=std::fread(&vy,sizeof(float),1,(std::FILE *) file);
		elemRead+=std::fread(&vz,sizeof(float),1,(std::FILE *) file);
		T_and_boosts[0]=T;
		T_and_boosts[1]=qgpFrac;
		T_and_boosts[2]=vx;
		T_and_boosts[3]=vy;
		T_and_boosts[4]=vz;
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
bool viscRead(bool binary, void * shearFile, void * bulkFile, float visc_info[]) {

	int shearElemRead, bulkElemRead;
	float Wtt=0.0, Wtx=0.0, Wty=0.0, Wtz=0.0, Wxx=0.0,Wxy=0.0,Wxz=0.0, Wyy=0.0, Wyz=0.0, Wzz=0.0;
	float bulk_pressure=0.0, eps_plus_P=0.0, cs2=0.0;
	bool return_value=true;

	//If binary
	if (binary) {
		if (CONST_with_shear_viscosity) {
			shearElemRead= std::fread(&Wtt,sizeof(float),1,(std::FILE *) shearFile);
			shearElemRead+=std::fread(&Wtx,sizeof(float),1,(std::FILE *) shearFile);
			shearElemRead+=std::fread(&Wty,sizeof(float),1,(std::FILE *) shearFile);
			shearElemRead+=std::fread(&Wtz,sizeof(float),1,(std::FILE *) shearFile);
			shearElemRead+=std::fread(&Wxx,sizeof(float),1,(std::FILE *) shearFile);
			shearElemRead+=std::fread(&Wxy,sizeof(float),1,(std::FILE *) shearFile);
			shearElemRead+=std::fread(&Wxz,sizeof(float),1,(std::FILE *) shearFile);
			shearElemRead+=std::fread(&Wyy,sizeof(float),1,(std::FILE *) shearFile);
			shearElemRead+=std::fread(&Wyz,sizeof(float),1,(std::FILE *) shearFile);
			shearElemRead+=std::fread(&Wzz,sizeof(float),1,(std::FILE *) shearFile);
			if (shearElemRead != 10) return_value=false;
		}
		if (CONST_with_bulk_viscosity) {
			bulkElemRead=std::fread(&bulk_pressure,sizeof(float),1,(std::FILE *) bulkFile);
			bulkElemRead+=std::fread(&eps_plus_P,sizeof(float),1,(std::FILE *) bulkFile);
			bulkElemRead+=std::fread(&cs2,sizeof(float),1,(std::FILE *) bulkFile);
			if (bulkElemRead != 3) return_value=false;
		}
	}
	else {
		if (CONST_with_shear_viscosity) {
			//shearElemRead=std::fscanf((std::FILE *) shearFile, "%f %f %f %f %f %f %f %f %f %f", &visc_info[0], &visc_info[1], &visc_info[2], &visc_info[3], &visc_info[4], &visc_info[5], &visc_info[6], &visc_info[7], &visc_info[8], &visc_info[9]);
			shearElemRead=std::fscanf((std::FILE *) shearFile, "%f %f %f %f %f %f %f %f %f %f", &Wtt, &Wtx, &Wty, &Wtz, &Wxx, &Wxy, &Wxz, &Wyy, &Wyz, &Wzz);
			if (shearElemRead != 10) return_value=false;
		}
		if (CONST_with_bulk_viscosity) {
			bulkElemRead=std::fscanf((std::FILE *) bulkFile, "%f %f %f", &bulk_pressure, &eps_plus_P, &cs2);
			if (bulkElemRead != 3) return_value=false;
		}
	}

	visc_info[0]=Wtt;
	visc_info[1]=Wtx;
	visc_info[2]=Wty;
	visc_info[3]=Wtz;
	visc_info[4]=Wxx;
	visc_info[5]=Wxy;
	visc_info[6]=Wxz;
	visc_info[7]=Wyy;
	visc_info[8]=Wyz;
	visc_info[9]=Wzz;
	visc_info[10]=bulk_pressure;
	visc_info[11]=eps_plus_P;
	visc_info[12]=cs2;

	//
	return return_value;
}


void pre_computeDescretizedSpectrum(struct phaseSpace_pos *curr_pos, float T_and_boosts[], float visc_info[], const struct photonRate rate_list[], double discSpectra[CONST_Nrap][CONST_Nphi][CONST_NkT][3][CONST_N_rates]) {

	//Forward declaration
	void computeDescretizedSpectrum(struct phaseSpace_pos *curr_pos, float T_and_boosts[], float visc_info[], const struct photonRate rate_list[], double discSpectra[CONST_Nrap][CONST_Nphi][CONST_NkT][3][CONST_N_rates]);

	if (!CONST_boost_invariant) {
		computeDescretizedSpectrum(curr_pos, T_and_boosts, visc_info, rate_list, discSpectra);
	}
	else {
		const int eta_slice_choice=cellNb_eta/2+1;
		//Integrate only over 1 slice in eta
		if (eta_slice_choice == curr_pos->ieta) {
			double eta, eta_slice;
			double old_u0, new_u0;
			double ux, uy;
			const int integration_step=2*int(CONST_nb_steps_eta_integration/2.0);
			const double delta_eta = 2.0*CONST_max_eta_integration/integration_step;
			eta_slice=(-cellNb_eta/2+eta_slice_choice-1)*CONST_cellsize_Eta;
			old_u0=1.0/sqrt(1-T_and_boosts[2]*T_and_boosts[2]-T_and_boosts[3]*T_and_boosts[3]-T_and_boosts[4]*T_and_boosts[4]);
			ux=T_and_boosts[2]*old_u0;
			uy=T_and_boosts[3]*old_u0;

			//Integrate in eta with trapezoidal method, using the symmetry around 0 to potentially speed-up the convergence
			//Integrate in eta
			for(int j=0;j<=integration_step;j++) {
				eta=-1*CONST_max_eta_integration+j*delta_eta;
				new_u0=sqrt((1+ux*ux+uy*uy)/(1-pow(tanh(eta),2))); 
				T_and_boosts[2]=ux/new_u0;
				T_and_boosts[3]=uy/new_u0;
				T_and_boosts[4]=tanh(eta);

				curr_pos->eta=eta;
				if ((0 == j)||(integration_step == j)) {
					curr_pos->w_eta=delta_eta/2.0;
				}
				else if (j%2 == 0) {
					curr_pos->w_eta=delta_eta;
				}
				else {
					curr_pos->w_eta=delta_eta;
				}

				computeDescretizedSpectrum(curr_pos, T_and_boosts, visc_info, rate_list, discSpectra);
			}
		}
	}

}


/***** Computation of the discretized spectrum *****/
void computeDescretizedSpectrum(struct phaseSpace_pos *curr_pos, float T_and_boosts[], float visc_info[], const struct photonRate rate_list[], double discSpectra[CONST_Nrap][CONST_Nphi][CONST_NkT][3][CONST_N_rates]) {
	
	//Forward declaration
	void fill_grid(struct phaseSpace_pos *curr_pos, double kR, double T, double Akk, double bulk_pressure, double eps_plus_P, double cs2, const struct photonRate rate_list[], double discSpectra[CONST_Nrap][CONST_Nphi][CONST_NkT][3][CONST_N_rates]);

	//Local variables
	double gamma;
	double rap, phi, kT;
	double kL, kLx, kLy, kLz;
	double kR, kOverTkOverTOver_e_P;
	double coshRap, sinhRap, cosPhi, sinPhi, invCoshRap;
	//double stPosition[4];

	//Assign those to local variables for convenience
	const double T=T_and_boosts[0];
	//const double qgpFrac=T_and_boosts[1];
	const double betax=T_and_boosts[2];
	const double betay=T_and_boosts[3];
	const double betaz=T_and_boosts[4];

	const double pitt=visc_info[0];
	const double pitx=visc_info[1];
	const double pity=visc_info[2];
	const double pitz=visc_info[3];
	const double pixx=visc_info[4];
	const double pixy=visc_info[5];
	const double pixz=visc_info[6];
	const double piyy=visc_info[7];
	const double piyz=visc_info[8];
	const double pizz=visc_info[9];

	const double bulk_pressure=visc_info[10];
	const double eps_plus_P=visc_info[11];
	const double cs2=visc_info[12];

	//pre-tabulate for speed
	double cosPhiArray[CONST_Nphi];
	for(int i=0; i<CONST_Nphi;i++) cosPhiArray[i]=cos(i*CONST_delPhi);

	double sinPhiArray[CONST_Nphi];
	for(int i=0; i<CONST_Nphi;i++) sinPhiArray[i]=sin(i*CONST_delPhi);

	//Pre-compute gamma for efficiency
	gamma=1.0/sqrt(1-(betax*betax+betay*betay+betaz*betaz));

	//Compute (tau,x,y,eta) from line number

	//Loop over transverse momentum kT, azimuthal angle phi and rapidity rap
	//(note that there is no different here between the rapidity and the pseudorapidity, the photon being massless)
	//Loop over kT
	for(int ikT=0;ikT<CONST_NkT; ikT++) {

		kT=CONST_kTMin+ikT*CONST_delKt;	

		//Loop over rapidity rap
		for(int irap=0;irap<CONST_Nrap; irap++) {

			rap=CONST_rapMin+irap*CONST_delRap;	

			coshRap=cosh(rap);
			sinhRap=sinh(rap);
			invCoshRap=1.0/coshRap;

			//Loop over phi (uniform discretization - to be used with the trapezoidal method)
			for(int iphi=0;iphi<CONST_Nphi; iphi++) {

				phi=iphi*CONST_delPhi;	

				//cosPhi=cos(phi);
				//sinPhi=sin(phi);
				cosPhi=cosPhiArray[iphi];
				sinPhi=sinPhiArray[iphi];


				//Evaluate (A_L)_{alpha beta} \hat{k}_L^alpha \alpha{k}_L^beta
				if (CONST_with_viscosity&&CONST_with_shear_viscosity) {
					//shear_info: Wtt,Wtx,Wty,Wtz,Wxx,Wxy,Wxz,Wyy,Wyz,Wzz
					//*shear_info+: 0   1   2   3   4   5   6   7   8   9
					//K=(kT cosh(y), kT cos(phi), kT sin(phi), kT sinh(y))
					//\hat{K}=(1,cos(phi)/cosh(rap),sin(phi)/cosh(rap),sinh(rap)/cosh(rap))
					if (!CONST_boost_invariant) {
						//kkPiOverEta=A00 + 1/cosh(rap)*( 2*(A01*k1+A02*k2+A03*k3) + 1/coshrap*() )
						//kkPiOverEta=A00 + 1/cosh(rap)^2*(A11 cos(phi)^2+A22*sin(phi)^2+A33*sinh(rap)^2)+2/cosh(rap)*(A01*cos(phi)+A02*sin(phi)+A03*sinh(rap)+1/cosh(rap)*(A12*cos(phi)*sin(phi)+A13*cos(phi)*sinh(rap)+A23*sin(phi)*sinh(rap)))
						//kkPiOverEta=*(shear_info) + invCoshEta*invCoshEta*( *(shear_info+4)*cosPhi*cosPhi + *(shear_info+7)*sinPhi*sinPhi + *(shear_info+9)*sinhEta*sinhEta) + 2.0*invCoshEta*(*(shear_info+1)*cosPhi + *(shear_info+2)*sinPhi + *(shear_info+3)*sinhEta + invCoshEta*( *(shear_info+5)*cosPhi*sinPhi + *(shear_info+6)*cosPhi*sinhEta + *(shear_info+8)*sinPhi*sinhEta));
						//kOverTkOverTOver_e_P=kT*kT*coshRap*coshRap*(visc_info[0] + invCoshRap*invCoshRap*( visc_info[4]*cosPhi*cosPhi + visc_info[7]*sinPhi*sinPhi + visc_info[9]*sinhRap*sinhRap) + 2.0*invCoshRap*( -1.0*visc_info[1]*cosPhi - visc_info[2]*sinPhi - visc_info[3]*sinhRap + invCoshRap*( visc_info[5]*cosPhi*sinPhi + visc_info[6]*cosPhi*sinhRap + visc_info[8]*sinPhi*sinhRap)));

						const double kt=kT*coshRap;
						const double kx=kT*cosPhi;
						const double ky=kT*sinPhi;
						const double kz=kT*sinhRap;

						kOverTkOverTOver_e_P=kt*(kt*pitt-2*kx*pitx-2*ky*pity-2*kz*pitz)+kx*(kx*pixx+2*ky*pixy+2*kz*pixz)+ky*(ky*piyy+2*kz*piyz)+kz*kz*pizz;

						kOverTkOverTOver_e_P/=T*T;





					}
					//In the boost-invariant case, \Pi^\mu\nu is the value a eta=0
					//k_\mu k\nu \Pi^\mu\nu must be calculed correctly
					else {
						//K=(kT cosh(y), kT cos(phi), kT sin(phi), kT sinh(y))
						//(k^tau,k^x,k^y,k^eta)=(kT cosh(y-eta),k^x,k^y,kT sinh(y-eta)/tau)
						//k=kT*cosh(rap)
						//shear_info: Wtt,Wtx,Wty,Wtz,Wxx,Wxy,Wxz,Wyy,Wyz,Wzz
						//*shear_info+: 0   1   2   3   4   5   6   7   8   9
						const double tau=curr_pos->tau;

						const double ktau=kT*cosh(rap-curr_pos->eta);
						const double kx=kT*cosPhi;
						const double ky=kT*sinPhi;
						const double keta=kT*sinh(rap-curr_pos->eta)/tau;

						const double eta_of_pimunu_slice=0.0;
						const double dtau_dt=cosh(eta_of_pimunu_slice);
						const double deta_dt=-1.0*sinh(eta_of_pimunu_slice)/tau;
						const double dtau_dz=-1.0*sinh(eta_of_pimunu_slice);
						const double deta_dz=cosh(eta_of_pimunu_slice)/tau;

						const double tau2=tau*tau;

						const double pitautau=dtau_dt*dtau_dt*pitt+2*dtau_dt*dtau_dz*pitz+dtau_dz*dtau_dz*pizz;
						const double pitaux=dtau_dt*pitx+dtau_dz*pixz;
						const double pitauy=dtau_dt*pity+dtau_dz*piyz;
						const double pitaueta=dtau_dt*deta_dt*pitt+(dtau_dt*deta_dz+dtau_dz*deta_dt)*pitz+dtau_dz*deta_dz*pizz;
						const double pixeta=deta_dt*pitx+deta_dz*pixz;
						const double piyeta=deta_dt*pity+deta_dz*piyz;
						const double pietaeta=deta_dt*deta_dt*pitt+2*deta_dt*deta_dz*pitz+deta_dz*deta_dz*pizz;

						kOverTkOverTOver_e_P=ktau*(ktau*pitautau-2*kx*pitaux-2*ky*pitauy-2*tau2*keta*pitaueta)+kx*(kx*pixx+2*ky*pixy+2*tau2*keta*pixeta)+ky*(ky*piyy+2*tau2*keta*piyeta)+tau2*tau2*keta*keta*pietaeta;

						kOverTkOverTOver_e_P/=T*T;

					}


					//double tr_check=(visc_info[0]-visc_info[4]-visc_info[7]-visc_info[9]);
					//if (tr_check > 1e-5) {
					//	std::cout << "Warning: Large deviation from tracelessness (" << tr_check << ")!\n";
					//}
					
					//Akk=(*shear_info)+invCoshRap( (*shear_info+4)*cosPhi*cosPhi);
				}
				else {
					kOverTkOverTOver_e_P=0.0;
				}
			
				//Photon momentum in the lab frame
				//k=mT cosh(rap)=kT cosh(rap)
				kL=kT*coshRap;
				//kx=kT cos(phi)
				kLx=kT*cosPhi;
				//ky=kT sin(phi)
				kLy=kT*sinPhi;
				//kz=mT sinh(rap)=kT sinh(rap)
				kLz=kT*sinhRap;

				//kR.uR=kL.uL
				//k_rf=(k_L-\vec{u}/u0.\vec{k})/sqrt(1-u^2/u0^2)
				kR=gamma*(kL-betax*kLx-betay*kLy-betaz*kLz);

				//Our rate
				//dGamma(\vec{k}_L)=dGamma_0(k_rf)+(A_L)_{alpha beta} k_L^alpha k_L^beta Z(rf) 
				curr_pos->ikT=ikT;
				curr_pos->iphi=iphi;
				curr_pos->irap=irap;
				fill_grid(curr_pos, kR, T, kOverTkOverTOver_e_P, bulk_pressure, eps_plus_P, cs2, rate_list,discSpectra);	


			}
		}	
	}

}

//Discretized spectra: array[times][Nrap][Nphi][Npt][rates]
//Discretized spectra, version 2: array[times][Nrap][Nphi][Npt][rates][value_and_remainder]
//stPos=[tau, irap, iphi, ikT]
void fill_grid(struct phaseSpace_pos *curr_pos, double kR, double T, double kOverTkOverTOver_e_P, double bulk_pressure, double eps_plus_P, double cs2, const struct photonRate rate_list[], double discSpectra[CONST_Nrap][CONST_Nphi][CONST_NkT][3][CONST_N_rates]) {

	//
	double QGP_fraction(double T);
	double eval_photon_rate(const struct photonRate * currRate, double kOverT, double T, double kkPiOver_e_P_k2, double bulk_pressure, double eps_plus_P, double cs2); 

	//
	int irap=curr_pos->irap;
	int iphi=curr_pos->iphi;
	int ikT=curr_pos->ikT;
	double tmpRate;
	//double (*local_rate)(double, double, double);
	double tmp_cellsize_eta;

	//Size of cell in eta different if integrating 2+1D or 3+1D
	if (CONST_boost_invariant) {
		tmp_cellsize_eta=curr_pos->w_eta;
	}
	else {
		tmp_cellsize_eta=CONST_cellsize_Eta;
	}

	//Loop over rates
	for(int iRate=0; iRate<CONST_N_rates;iRate++) {
	//double discSpectra[CONST_Neta][CONST_Nphi][CONST_NkT][3][CONST_rateList.size()];

		//get_photon_rate(CONST_rates_to_use[iRate], &local_rate);

		tmpRate=eval_photon_rate(&rate_list[iRate],kR/T,T,kOverTkOverTOver_e_P, bulk_pressure, eps_plus_P, cs2);

//		//Either use the fits directly or a table made from the fits
//		if (CONST_use_accel_rates[CONST_rates_to_use[iRate]-1]) {
//			tmpRate=get_photon_rate_accel(kR/T, T, kOverTkOverTOver_e_P, iRate);
//		}
//		else {
//			tmpRate=(*local_rate)(kR/T,T,kOverTkOverTOver_e_P);
//		}
		//tmpRate=kR*(1.0+cos(kR*(2.0)*iphi*CONST_delPhi));

		//QGP fraction
		//tmpRate*=QGP_fraction(T);

		//Cell volume: dx*dy*dz*dt=dx*dy*dEta*dTau*tau
		tmpRate*=CONST_cellsize_X*CONST_cellsize_Y*tmp_cellsize_eta*CONST_effective_dTau*curr_pos->tau;

		//Fill value
		discSpectra[irap][iphi][ikT][1][iRate]+=tmpRate;

		//Fill lower bound uncertainty
		discSpectra[irap][iphi][ikT][0][iRate]=0.0;

		//Fill upper bound uncertainty
		discSpectra[irap][iphi][ikT][2][iRate]=0.0;

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


void compute_observables(const struct photonRate rate_list[], double discSpectra[CONST_Nrap][CONST_Nphi][CONST_NkT][3][CONST_N_rates]) {

	//
	void compute_midrapidity_yield_and_vn(const struct photonRate currRate[], double discSpectra[CONST_Nrap][CONST_Nphi][CONST_NkT][3][CONST_N_rates]);

	compute_midrapidity_yield_and_vn(rate_list, discSpectra);

}

//Output the phi-integrated, rapidity-averaged-around-0 yield as a function of pT
void compute_midrapidity_yield_and_vn(const struct photonRate rate_list[], double discSpectra[CONST_Nrap][CONST_Nphi][CONST_NkT][3][CONST_N_rates]) {

	double kT, rap, phi, yFac, rap_interval;
	double yield, vn_sin[CONST_FourierNb], vn_cos[CONST_FourierNb];
	int iRapmin=0, iRapmax;
	bool exact_midrap=false, bad_rap_discret = false;

	//One file per rate
	for(int rate_no=0; rate_no<CONST_N_rates; rate_no++) {

		//Set output file name
		std::stringstream tmpStr;
		tmpStr << "vn_";
		tmpStr << rate_list[rate_no].name.c_str();
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
			//outfile	<< "\tyield*vn[" << j << "]\tvn[" << j << "]";
			outfile	<< "\tyield*vn_cos[" << j << "]\t yield*vn_sin[" << j << "]";
		}
		outfile<< "\n";
			

		//Identify the cells in rapidity that should be averaged over
		for(int irap=0;irap<CONST_Nrap; irap++) {

			rap=CONST_rapMin+irap*CONST_delRap;	

			if (fabs(rap) <= CONST_midRapCut) {
				iRapmin=irap;
				iRapmax=irap;
			}
			else if (fabs(rap) > CONST_midRapCut) {
				//iRapmax=irap-1;
				continue;
			}
			else iRapmax=irap;
		
		}

		//If there is a single point and it is rap=0.0, 
		if (iRapmin == iRapmax) {
			if (0.0 == CONST_rapMin+iRapmin*CONST_delRap) {
				exact_midrap=true;
				rap_interval=1.0;
			}
			else {
				bad_rap_discret=true;
			}
		}
		else {
			rap_interval=(iRapmax-iRapmin)*CONST_delRap;
		}

		if (bad_rap_discret) {
				outfile << "Can't evaluate midrapidity spectra with current rapidity discretization\n";
		}
		else {

			for(int ikT=0;ikT<CONST_NkT; ikT++) {

				//Will contain the results of the phi integration and rapidity averaging
				yield=0.0;
				for(int i=1;i<=CONST_FourierNb; i++) {
					vn_sin[i-1]=0;
					vn_cos[i-1]=0;
				}

				kT=CONST_kTMin+ikT*CONST_delKt;	

				//Loop over rapidity rap
				for(int irap=iRapmin;irap<=iRapmax; irap++) {

					rap=CONST_rapMin+irap*CONST_delRap;	

					//Loop over phi (trapezoidal method)
					for(int iphi=0;iphi<CONST_Nphi-1; iphi++) {

						phi=iphi*CONST_delPhi;	
						
						if (exact_midrap) {
							yFac=1.0;
						}
						//Let's use a simple midpoint rule for now
						else {
							yFac=CONST_delRap;
						}

						//Finally, multiply by dNdydptdphi[NY][NPT][NPHI+1] 
						//tmpIntRes+=phiFac*yFac*particleList[j].dNdydptdphi[iy][ipt][iphi];
						yield+=discSpectra[irap][iphi][ikT][1][rate_no]*yFac;
						for(int i=1;i<=CONST_FourierNb; i++) {
							vn_cos[i-1]+=yFac*discSpectra[irap][iphi][ikT][1][rate_no]*cos(i*phi);
							vn_sin[i-1]+=yFac*discSpectra[irap][iphi][ikT][1][rate_no]*sin(i*phi);
						}
						
					}
					
				}

				//Multiply by delta_ph and divide by the rapidity integration range 
				//(to yield an average instead of an integral)
				yield*=CONST_delPhi/rap_interval/(2.0*M_PI);
				for(int j=1;j<=CONST_FourierNb; j++) {
					vn_sin[j-1]*=CONST_delPhi/rap_interval/(2.0*M_PI);
					vn_cos[j-1]*=CONST_delPhi/rap_interval/(2.0*M_PI);
				}

				//Output result
				outfile << kT << "\t" << yield;
				for(int j=1;j<=CONST_FourierNb; j++) {
					//outfile	<< "\t" << vn[j] << "\t" << vn[j]/yield;
					outfile	<< "\t" << vn_cos[j-1] << "\t" << vn_sin[j-1];
				}
				outfile<< "\n";


			}

		}

		//Close file
		outfile.close();
	}

}
