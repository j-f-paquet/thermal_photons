#include "photon.h"
#include "rates.h"

// #####################################################################################################################
// ### This code convolute a near-thermal photon emission rate (photon emission per spacetime volume)
// ### with a hydrodynamic spacetime profile.
// ###
// ### Rates are assumed to take the form
// ### k d3Gamma/d3k = k d3Gamma_{ideal}/d3k 
// ###                  + \pi^{\mu\nu} K^\mu K^\nu (shear_correction)
// ###                  + \Pi (bulk_correction)
// ### where the thermal rate "k d3Gamma_{ideal}/d3k" only depends on the temperature T and the flow velocity u^\mu
// ### while the shear/bulk_corrections can also depend on other thermodynamic quantities
// ###
// ### All photon rates are defined in the "rates.cpp" file
// #####################################################################################################################


unsigned long int GLOBAL_line_number=0;

//Main
int main() {

	//
	struct photonRate rate_list[CONST_N_rates];

	//Initialise rates
	for(int i=0;i<CONST_N_rates;i++) init_rates(&rate_list[i], CONST_rates_to_use[i]);


	//Compute photon production
	photon_prod(rate_list);
}


//Compute photon production
void photon_prod(const struct photonRate rate_list[]) {

	//Variables
	struct hydro_info_t hydro_info;

	bool read_T_flag; //Result of reading of the file
	std::FILE * hydro_fields_files[3]; 

	//The second to last dimension is meant for including an upper and a lower bound on the uncertainty, if possible
	//discSpectra[][][][0/1/2][] is for the lower bound/value/upper bound
	double discSpectra[CONST_N_rates][CONST_NkT][CONST_Nrap][CONST_Nphi][3] = {0.0};

	// Code under construction...
	if (CONST_with_viscosity) {
			std::cout << "Many parts of the code still have to be finished so that viscous corrections can be calculated correctly. Not working now.\n;";
			exit(1);
	}

        std::cout << "Computing thermal photons...\n";

	// Loop over file containing hydro fields
	init_hydro_field_files(hydro_fields_files);
	read_T_flag=read_hydro_fields(hydro_fields_files, hydro_info);
	while (read_T_flag) {

		if (hydro_info.T >= CONST_freezeout_T) {
			pre_computeDescretizedSpectrum(hydro_info, rate_list, discSpectra);
		}

		//Try to read the next line
		read_T_flag=read_hydro_fields(hydro_fields_files, hydro_info);

	}
	//if ((!std::feof(stFile))||((CONST_with_viscosity)&&(!std::feof(shearViscFile)))) {
	//	std::cout << "!!!!!!!!!!!!!!! Warning !!!!!!!!!!!!!!!!!! Stopped reading the evolution files before the end of the file!\n";	
	//}
	close_hydro_field_files(hydro_fields_files);

	//Compute observables from the discretized photon spectra
	compute_observables(rate_list, discSpectra);

}

/***** File reading stuff *****/

//Open file for reading
bool open_file_read(bool binary, std::string filename, std::FILE ** pointer) {

	bool return_value=true;

	//If binary
	if (binary) {
		*pointer=std::fopen(filename.c_str(),"rb");
	}
	else {
		*pointer=std::fopen(filename.c_str(),"r");
		//std::cout << "pointer=" << pointer << "\n";
		//float test;
		//int elem_read=std::fscanf(pointer, "%f", &test); 
		//std::cout << "elemread" << elem_read << " & float=" << test << "\n";
		//elem_read=std::fscanf(pointer, "%f", &test); 
		//std::cout << "elemread" << elem_read << " & float=" << test << "\n";
		//elem_read=std::fscanf(pointer, "%f", &test); 
		//std::cout << "elemread" << elem_read << " & float=" << test << "\n";
	//	exit(1);
	}

	//Check if it opened correctly
	if (NULL == *pointer) {
	       //	std::ferror(pointer)) {
		std::printf("Error: Could not open file \"%s\"",filename.c_str());
		return_value=false;
	}

	return return_value;

}

bool init_hydro_field_files(std::FILE * hydro_fields_files[3]) {

	bool return_value;

	if (CONST_file_format == new_format) {
		return_value=open_file_read(CONST_binaryMode,stGridFile,&hydro_fields_files[0]);
	}
	else {
		return_value=open_file_read(CONST_binaryMode,stGridFile,&hydro_fields_files[0]);
		bool shear_return_value=true;
		bool bulk_return_value=true;
		if (CONST_with_viscosity) {
			if (CONST_with_shear_viscosity) shear_return_value=open_file_read(CONST_binaryMode,shearViscosityFile,&hydro_fields_files[1]);
			if (CONST_with_bulk_viscosity) bulk_return_value=open_file_read(CONST_binaryMode,bulkViscosityFile,&hydro_fields_files[2]);
		}
		return_value = return_value && shear_return_value && bulk_return_value;
	}

	return return_value;

}

void close_hydro_field_files(std::FILE * hydro_fields_files[3]) {

	std::fclose(hydro_fields_files[0]);
	if ((CONST_with_viscosity)&&(CONST_file_format == old_format)) {
		if (CONST_with_shear_viscosity) std::fclose(hydro_fields_files[1]);
		if (CONST_with_bulk_viscosity) std::fclose(hydro_fields_files[2]);
	}
}

//Read spacetime file 
bool read_hydro_fields(std::FILE * hydro_fields_files[3], struct hydro_info_t & hydro_info) {

	bool return_value;

	if (CONST_file_format == new_format) {
		return_value=read_hydro_fields_new_format(hydro_fields_files,hydro_info);
	}
	else {
		return_value=read_hydro_fields_old_format(hydro_fields_files,hydro_info);
	}


	return return_value;

}

bool read_hydro_fields_new_format(std::FILE * hydro_fields_files[3], struct hydro_info_t & hydro_info) {

	const int elem_to_read=12;

	int elem_read;

	const bool binary=CONST_binaryMode;

	std::FILE * tmp_file=hydro_fields_files[0];

	//float ideal[] = {static_cast<float>(volume),
	//	static_cast<float>(eta_local),
	//	static_cast<float>(T_local*hbarc),
	//	static_cast<float>(ux),
	//	static_cast<float>(uy),
	//	static_cast<float>(ueta)};
	//fwrite(ideal, sizeof(float), 6, out_file_xyeta);
	//if (DATA.turn_on_shear == 1) {
	//	float shear_pi[] = {static_cast<float>(Wxx),
	//		static_cast<float>(Wxy),
	//		static_cast<float>(Wxeta),
	//		static_cast<float>(Wyy),
	//		static_cast<float>(Wyeta)};
	//	fwrite(shear_pi, sizeof(float), 5, out_file_xyeta);
	//}
	//if (DATA.turn_on_bulk == 1) {
	//	float bulk_pi[] = {static_cast<float>(pi_b)};
	//	fwrite(bulk_pi, sizeof(float), 1, out_file_xyeta);
	//}
	float volume, T, eta_s, ux, uy, tau_ueta;
	float Wxx, Wxy, Wxeta, Wyy, Wyeta;
	float Pi_b;

	//If binary
	if (binary) {
		//float muB;
		elem_read=std::fread(&volume,sizeof(float),1,tmp_file);
		elem_read+=std::fread(&eta_s,sizeof(float),1,tmp_file);
		elem_read+=std::fread(&T,sizeof(float),1,tmp_file);
		elem_read+=std::fread(&ux,sizeof(float),1,tmp_file);
		elem_read+=std::fread(&uy,sizeof(float),1,tmp_file);
		elem_read+=std::fread(&tau_ueta,sizeof(float),1,tmp_file);
		//elem_read+=std::fread(&muB,sizeof(float),1,tmp_file);
		
		elem_read+=std::fread(&Wxx,sizeof(float),1,tmp_file);
		elem_read+=std::fread(&Wxy,sizeof(float),1,tmp_file);
		elem_read+=std::fread(&Wxeta,sizeof(float),1,tmp_file);
		elem_read+=std::fread(&Wyy,sizeof(float),1,tmp_file);
		elem_read+=std::fread(&Wyeta,sizeof(float),1,tmp_file);

		elem_read+=std::fread(&Pi_b,sizeof(float),1,tmp_file);


	}
	else {
		elem_read=std::fscanf(tmp_file, "%f %f %f %f %f %f %f %f %f %f %f %f", &volume, &eta_s, &T, &ux, &uy, &tau_ueta, &Wxx, &Wxy, &Wxeta, &Wyy, &Wyeta, &Pi_b);
	}

	// For the boost-invariant case, don't include "delta_eta" in volume for now --- will be added later
	// Since the volume is already pre-computed in this file format, divide it of its fake "delta_eta"
	const double fake_deta=0.1;
	if (CONST_boost_invariant) volume/=fake_deta;

	hydro_info.V4=volume;
	hydro_info.eta_s=eta_s;
	hydro_info.T=T;
	hydro_info.ux=ux;
	hydro_info.uy=uy;
	hydro_info.tau_ueta=tau_ueta;

	hydro_info.pixx=Wxx;
	hydro_info.pixy=Wxy;
	hydro_info.pixeta=Wxeta;
	hydro_info.piyy=Wyy;
	hydro_info.piyeta=Wyeta;

	hydro_info.Pi_b=Pi_b;

	// The current format doesn't have enough information for my code to compute the viscous correction to photons
	// Need to: reconstruct other components of pimunu from the five available, which is straightforward
	// Plus need: EOS information to get epsilon+P and cs2 from T for bulk correction, which isn't how this code was designed... but it can always be fixed. I don't like this because there's always the risk that the wrong EOS information is passed. But there's no way around that without saving a lot of redundant information in the hydro files.
	if (CONST_with_viscosity) {
			std::cout << "This part of the code isn't currently working...\n;";
			exit(1);
	}

	//If fscanf couldn't read exactly the right number of elements, it's the end of the file or there's a problem
	if (elem_read != elem_to_read) {
		return false;
	}
	else {
		return true;
	}

}

//Read spacetime file 
bool read_hydro_fields_old_format(std::FILE * hydro_fields_files[3], struct hydro_info_t & hydro_info) {

	bool return_value=true;

	const bool binary=CONST_binaryMode;

	int elem_to_read=5;
	//if (viscosity) element_to_read+=...

	int elem_read;

	std::FILE * tmp_ideal_file=hydro_fields_files[0];

	float T,qgpFrac,vx,vy,vz;

	// Read first file with ideal hydro field
	if (binary) {
		elem_read=std::fread(&T,sizeof(float),1,tmp_ideal_file);
		elem_read+=std::fread(&qgpFrac,sizeof(float),1,tmp_ideal_file);
		elem_read+=std::fread(&vx,sizeof(float),1,tmp_ideal_file);
		elem_read+=std::fread(&vy,sizeof(float),1,tmp_ideal_file);
		elem_read+=std::fread(&vz,sizeof(float),1,tmp_ideal_file);
	}
	else {
		elem_read=std::fscanf(tmp_ideal_file, "%f %f %f %f %f", &T, &qgpFrac, &vx, &vy, &vz);
	}

	if (elem_read != elem_to_read) return_value=false;

	float bulk_pressure=0.0, eps_plus_P=0.0, cs2=0.0;
	float Wtt=0.0, Wtx=0.0, Wty=0.0, Wtz=0.0, Wxx=0.0,Wxy=0.0,Wxz=0.0, Wyy=0.0, Wyz=0.0, Wzz=0.0;

	if (CONST_with_viscosity) {


		if (CONST_with_shear_viscosity) { 

			int shear_elem_read=0;

			const int shear_elem_to_read=10;

			std::FILE * tmp_shear_file=hydro_fields_files[1];
			//If binary
			if (binary) {
				shear_elem_read= std::fread(&Wtt,sizeof(float),1,tmp_shear_file);
				shear_elem_read+=std::fread(&Wtx,sizeof(float),1,tmp_shear_file);
				shear_elem_read+=std::fread(&Wty,sizeof(float),1,tmp_shear_file);
				shear_elem_read+=std::fread(&Wtz,sizeof(float),1,tmp_shear_file);
				shear_elem_read+=std::fread(&Wxx,sizeof(float),1,tmp_shear_file);
				shear_elem_read+=std::fread(&Wxy,sizeof(float),1,tmp_shear_file);
				shear_elem_read+=std::fread(&Wxz,sizeof(float),1,tmp_shear_file);
				shear_elem_read+=std::fread(&Wyy,sizeof(float),1,tmp_shear_file);
				shear_elem_read+=std::fread(&Wyz,sizeof(float),1,tmp_shear_file);
				shear_elem_read+=std::fread(&Wzz,sizeof(float),1,tmp_shear_file);
			}
			else {
				shear_elem_read=std::fscanf(tmp_shear_file, "%f %f %f %f %f %f %f %f %f %f", &Wtt, &Wtx, &Wty, &Wtz, &Wxx, &Wxy, &Wxz, &Wyy, &Wyz, &Wzz);
			}

			if (shear_elem_read != shear_elem_to_read) return_value=false;

		}
		if (CONST_with_bulk_viscosity) {

			int bulk_elem_read=0;

			const int bulk_elem_to_read=3;

			std::FILE * tmp_bulk_file=hydro_fields_files[2];

			//If binary
			if (binary) {
				bulk_elem_read=std::fread(&bulk_pressure,sizeof(float),1,tmp_bulk_file);
				bulk_elem_read+=std::fread(&eps_plus_P,sizeof(float),1,tmp_bulk_file);
				bulk_elem_read+=std::fread(&cs2,sizeof(float),1,tmp_bulk_file);
			}
			else {
				bulk_elem_read=std::fscanf(tmp_bulk_file, "%f %f %f", &bulk_pressure, &eps_plus_P, &cs2);
			}
			if (bulk_elem_read != bulk_elem_to_read) return_value=false;
		}

	}

	//If fscanf couldn't read exactly the right number of elements, it's the end of the file or there's a problem
	if (return_value) {
		
		//float ux, uy, ueta, tau, volume, eta_s;
		// Determine tau and then volume
		const int itau=int((GLOBAL_line_number/(cellNb_x*cellNb_y*cellNb_eta)));
		const double tau=MUSIC_tau0+CONST_effective_dTau*itau; //get_tau_from_linenumber();

		// For the boost-invariant case, don't include deta in volume for now --- will be added later
		double volume=CONST_cellsize_X*CONST_cellsize_Y*CONST_effective_dTau*tau;
		if (!CONST_boost_invariant) volume*=CONST_cellsize_Eta;

//		std::cout << "Tau is" << tau << "\n";

		// For boost-invariant hydro,  eta_s will be integrated over
		// but we need to know which slice in eta was saved: CONST_eta_s_of_saved_slice
		double eta_s;
		if (CONST_boost_invariant) {
			eta_s=CONST_eta_s_of_saved_slice;
		}
		else {
			// If the hydro fields are not boost-invariant, we don't need to save eta_s,
			// but we need it in this function
			const int ieta=int((GLOBAL_line_number % (cellNb_x*cellNb_y*cellNb_eta) )/(cellNb_x*cellNb_y));
			eta_s=(ieta-cellNb_eta/2);
		}

		// Get ux, uy, ueta from vx, vy, vz
		const double ut=1.0/sqrt(1-vx*vx-vy*vy-vz*vz);
		const double ux=vx*ut;
		const double uy=vy*ut;
		const double uz=vz*ut;
		const double tau_ueta=-sinh(eta_s)*ut+cosh(eta_s)*uz;


		hydro_info.V4=volume;
		hydro_info.eta_s=eta_s;
		hydro_info.T=T;
		hydro_info.ux=ux;
		hydro_info.uy=uy;
		hydro_info.tau_ueta=tau_ueta;

		// Shear viscosity related

		const double dtau_dt=cosh(eta_s);
		const double deta_dt=-1.0*sinh(eta_s)/tau;
		const double dtau_dz=-1.0*sinh(eta_s);
		const double deta_dz=cosh(eta_s)/tau;

		hydro_info.pitautau=dtau_dt*dtau_dt*Wtt+2*dtau_dt*dtau_dz*Wtz+dtau_dz*dtau_dz*Wzz;
		hydro_info.pitaux=dtau_dt*Wtx+dtau_dz*Wxz;
		hydro_info.pitauy=dtau_dt*Wty+dtau_dz*Wyz;
		hydro_info.pitaueta=dtau_dt*deta_dt*Wtt+(dtau_dt*deta_dz+dtau_dz*deta_dt)*Wtz+dtau_dz*deta_dz*Wzz;
		hydro_info.pixeta=deta_dt*Wtx+deta_dz*Wxz;
		hydro_info.piyeta=deta_dt*Wty+deta_dz*Wyz;
		hydro_info.pietaeta=deta_dt*deta_dt*Wtt+2*deta_dt*deta_dz*Wtz+deta_dz*deta_dz*Wzz;
		hydro_info.pixx=Wxx;
		hydro_info.pixy=Wxy;
		hydro_info.piyy=Wyy;

		// Bulk viscosity plus other information needed
		hydro_info.Pi_b=bulk_pressure;
		hydro_info.epsilon_plus_P=eps_plus_P;
		hydro_info.cs2=cs2;

		// Remember that one line was read
		GLOBAL_line_number++;

	}

	return return_value;

}

void pre_computeDescretizedSpectrum(struct hydro_info_t & hydro_info, const struct photonRate rate_list[], double discSpectra[CONST_N_rates][CONST_NkT][CONST_Nrap][CONST_Nphi][3]) {

	if (!CONST_boost_invariant) {
		computeDescretizedSpectrum(hydro_info, rate_list, discSpectra);
	}
	else {
		const int integration_steps=2*int(CONST_nb_steps_eta_integration/2.0);
		const double delta_eta = 2.0*CONST_max_eta_integration/integration_steps;

		// Value of V4 without "delta_eta"
		const double pre_V4=hydro_info.V4;

		//Integrate in eta with trapezoidal method, using the symmetry around 0 to potentially speed-up the convergence
		for(int j=0;j<=integration_steps;j++) {
			const double eta_s=-1*CONST_max_eta_integration+j*delta_eta;

			hydro_info.eta_s=eta_s;

			double deta;
			if ((0 == j)||(integration_steps == j)) {
				deta=delta_eta/2.0;
			}
			else if (j%2 == 0) {
				deta=delta_eta;
			}
			else {
				deta=delta_eta;
			}

			hydro_info.V4=pre_V4*deta;

			computeDescretizedSpectrum(hydro_info, rate_list, discSpectra);
		}

		//const int eta_slice_choice=cellNb_eta/2+1;
		////Integrate only over 1 slice in eta
		//if (eta_slice_choice == curr_pos->ieta) {
		//	double eta, eta_slice;
		//	double old_u0, new_u0;
		//	double ux, uy;
		//	const int integration_step=2*int(CONST_nb_steps_eta_integration/2.0);
		//	const double delta_eta = 2.0*CONST_max_eta_integration/integration_step;
		//	eta_slice=(-cellNb_eta/2+eta_slice_choice-1)*CONST_cellsize_Eta;
		//	old_u0=1.0/sqrt(1-T_and_boosts[2]*T_and_boosts[2]-T_and_boosts[3]*T_and_boosts[3]-T_and_boosts[4]*T_and_boosts[4]);
		//	ux=T_and_boosts[2]*old_u0;
		//	uy=T_and_boosts[3]*old_u0;

		//	//Integrate in eta with trapezoidal method, using the symmetry around 0 to potentially speed-up the convergence
		//	//Integrate in eta
		//	for(int j=0;j<=integration_step;j++) {
		//		eta=-1*CONST_max_eta_integration+j*delta_eta;
		//		new_u0=sqrt((1+ux*ux+uy*uy)/(1-pow(tanh(eta),2))); 
		//		T_and_boosts[2]=ux/new_u0;
		//		T_and_boosts[3]=uy/new_u0;
		//		T_and_boosts[4]=tanh(eta);

		//		curr_pos->eta=eta;
		//		if ((0 == j)||(integration_step == j)) {
		//			curr_pos->w_eta=delta_eta/2.0;
		//		}
		//		else if (j%2 == 0) {
		//			curr_pos->w_eta=delta_eta;
		//		}
		//		else {
		//			curr_pos->w_eta=delta_eta;
		//		}

		//		computeDescretizedSpectrum(curr_pos, T_and_boosts, visc_info, rate_list, discSpectra);
		//	}
		//}
	}

}


/***** Computation of the discretized spectrum *****/
void computeDescretizedSpectrum(struct hydro_info_t & hydro_info, const struct photonRate rate_list[], double discSpectra[CONST_N_rates][CONST_NkT][CONST_Nrap][CONST_Nphi][3]) {
	
	//Assign those to local variables for convenience
	const double V4=hydro_info.V4;
	const double T=hydro_info.T;
//	const double muB=hydro_info.muB;
	//const double qgpFrac=V_T_and_boosts[1];
	const double ux=hydro_info.ux;
	const double uy=hydro_info.uy;
	const double tau_ueta=hydro_info.tau_ueta;
	const double utau=sqrt(1.+ux*ux+uy*uy+tau_ueta*tau_ueta);
	//
	const double coshEtaS=cosh(hydro_info.eta_s);
	const double sinhEtaS=sinh(hydro_info.eta_s);

	//const double pitt=visc_info[0];
	//const double pitx=visc_info[1];
	//const double pity=visc_info[2];
	//const double pitz=visc_info[3];
	//const double pixx=visc_info[4];
	//const double pixy=visc_info[5];
	//const double pixz=visc_info[6];
	//const double piyy=visc_info[7];
	//const double piyz=visc_info[8];
	//const double pizz=visc_info[9];

	//const double bulk_pressure=visc_info[10];
	//const double eps_plus_P=visc_info[11];
	//const double cs2=visc_info[12];

	//pre-tabulate for speed
	double cosPhiArray[CONST_Nphi];
	for(int i=0; i<CONST_Nphi;i++) cosPhiArray[i]=cos(i*CONST_delPhi);

	double sinPhiArray[CONST_Nphi];
	for(int i=0; i<CONST_Nphi;i++) sinPhiArray[i]=sin(i*CONST_delPhi);

        //Loop over rates
        for(int iRate=0; iRate<CONST_N_rates;iRate++) {

		const struct photonRate current_rate=rate_list[iRate]; 

                //Loop over transverse momentum kT, azimuthal angle phi and rapidity rap
                //(note that there is no different here between the rapidity and the pseudorapidity, the photon being massless)
                //Loop over kT
		#pragma omp parallel for collapse(3) shared(discSpectra)
                for(int ikT=0;ikT<CONST_NkT; ikT++) {

                        //Loop over rapidity rap
                        for(int irap=0;irap<CONST_Nrap; irap++) {

                                //Loop over phi (uniform discretization - to be used with the trapezoidal method)
                                for(int iphi=0;iphi<CONST_Nphi; iphi++) {

					const double kT=CONST_kTMin+ikT*CONST_delKt;	

					const double rap=CONST_rapMin+irap*CONST_delRap;	

					const double coshRap=cosh(rap);
					const double sinhRap=sinh(rap);
					//invCoshRap=1.0/coshRap;

                                        //const double phi=iphi*CONST_delPhi;	

                                        //cosPhi=cos(phi);
                                        //sinPhi=sin(phi);
                                        const double cosPhi=cosPhiArray[iphi];
                                        const double sinPhi=sinPhiArray[iphi];

					double kOverTkOverTOver_e_P=0.0;
					double bulk_pressure=0.0;
					double eps_plus_P=0.0;
					double cs2=0.0;


                                        //Evaluate (A_L)_{alpha beta} \hat{k}_L^alpha \alpha{k}_L^beta
//                                        if (CONST_with_viscosity&&CONST_with_shear_viscosity) {
//                                                //shear_info: Wtt,Wtx,Wty,Wtz,Wxx,Wxy,Wxz,Wyy,Wyz,Wzz
//                                                //*shear_info+: 0   1   2   3   4   5   6   7   8   9
//                                                //K=(kT cosh(y), kT cos(phi), kT sin(phi), kT sinh(y))
//                                                //\hat{K}=(1,cos(phi)/cosh(rap),sin(phi)/cosh(rap),sinh(rap)/cosh(rap))
//                                                if (!CONST_boost_invariant) {
//                                                        //kkPiOverEta=A00 + 1/cosh(rap)*( 2*(A01*k1+A02*k2+A03*k3) + 1/coshrap*() )
//                                                        //kkPiOverEta=A00 + 1/cosh(rap)^2*(A11 cos(phi)^2+A22*sin(phi)^2+A33*sinh(rap)^2)+2/cosh(rap)*(A01*cos(phi)+A02*sin(phi)+A03*sinh(rap)+1/cosh(rap)*(A12*cos(phi)*sin(phi)+A13*cos(phi)*sinh(rap)+A23*sin(phi)*sinh(rap)))
//                                                        //kkPiOverEta=*(shear_info) + invCoshEta*invCoshEta*( *(shear_info+4)*cosPhi*cosPhi + *(shear_info+7)*sinPhi*sinPhi + *(shear_info+9)*sinhEta*sinhEta) + 2.0*invCoshEta*(*(shear_info+1)*cosPhi + *(shear_info+2)*sinPhi + *(shear_info+3)*sinhEta + invCoshEta*( *(shear_info+5)*cosPhi*sinPhi + *(shear_info+6)*cosPhi*sinhEta + *(shear_info+8)*sinPhi*sinhEta));
//                                                        //kOverTkOverTOver_e_P=kT*kT*coshRap*coshRap*(visc_info[0] + invCoshRap*invCoshRap*( visc_info[4]*cosPhi*cosPhi + visc_info[7]*sinPhi*sinPhi + visc_info[9]*sinhRap*sinhRap) + 2.0*invCoshRap*( -1.0*visc_info[1]*cosPhi - visc_info[2]*sinPhi - visc_info[3]*sinhRap + invCoshRap*( visc_info[5]*cosPhi*sinPhi + visc_info[6]*cosPhi*sinhRap + visc_info[8]*sinPhi*sinhRap)));
//
//                                                        const double kt=kT*coshRap;
//                                                        const double kx=kT*cosPhi;
//                                                        const double ky=kT*sinPhi;
//                                                        const double kz=kT*sinhRap;
//
//                                                        kOverTkOverTOver_e_P=kt*(kt*pitt-2*kx*pitx-2*ky*pity-2*kz*pitz)+kx*(kx*pixx+2*ky*pixy+2*kz*pixz)+ky*(ky*piyy+2*kz*piyz)+kz*kz*pizz;
//
//                                                        kOverTkOverTOver_e_P/=T*T;
//
//
//
//
//
//                                                }
//                                                //In the boost-invariant case, \Pi^\mu\nu is the value a eta=0
//                                                //k_\mu k\nu \Pi^\mu\nu must be calculed correctly
//                                                else {
//                                                        //K=(kT cosh(y), kT cos(phi), kT sin(phi), kT sinh(y))
//                                                        //(k^tau,k^x,k^y,k^eta)=(kT cosh(y-eta),k^x,k^y,kT sinh(y-eta)/tau)
//                                                        //k=kT*cosh(rap)
//                                                        //shear_info: Wtt,Wtx,Wty,Wtz,Wxx,Wxy,Wxz,Wyy,Wyz,Wzz
//                                                        //*shear_info+: 0   1   2   3   4   5   6   7   8   9
//                                                        const double tau=curr_pos_copy.tau;
//
//                                                        const double ktau=kT*cosh(rap-curr_pos_copy.eta);
//                                                        const double kx=kT*cosPhi;
//                                                        const double ky=kT*sinPhi;
//                                                        const double keta=kT*sinh(rap-curr_pos_copy.eta)/tau;
//
//                                                        const double eta_of_pimunu_slice=0.0;
//                                                        const double dtau_dt=cosh(eta_of_pimunu_slice);
//                                                        const double deta_dt=-1.0*sinh(eta_of_pimunu_slice)/tau;
//                                                        const double dtau_dz=-1.0*sinh(eta_of_pimunu_slice);
//                                                        const double deta_dz=cosh(eta_of_pimunu_slice)/tau;
//
//                                                        const double tau2=tau*tau;
//
//                                                        const double pitautau=dtau_dt*dtau_dt*pitt+2*dtau_dt*dtau_dz*pitz+dtau_dz*dtau_dz*pizz;
//                                                        const double pitaux=dtau_dt*pitx+dtau_dz*pixz;
//                                                        const double pitauy=dtau_dt*pity+dtau_dz*piyz;
//                                                        const double pitaueta=dtau_dt*deta_dt*pitt+(dtau_dt*deta_dz+dtau_dz*deta_dt)*pitz+dtau_dz*deta_dz*pizz;
//                                                        const double pixeta=deta_dt*pitx+deta_dz*pixz;
//                                                        const double piyeta=deta_dt*pity+deta_dz*piyz;
//                                                        const double pietaeta=deta_dt*deta_dt*pitt+2*deta_dt*deta_dz*pitz+deta_dz*deta_dz*pizz;
//
//                                                        kOverTkOverTOver_e_P=ktau*(ktau*pitautau-2*kx*pitaux-2*ky*pitauy-2*tau2*keta*pitaueta)+kx*(kx*pixx+2*ky*pixy+2*tau2*keta*pixeta)+ky*(ky*piyy+2*tau2*keta*piyeta)+tau2*tau2*keta*keta*pietaeta;
//
//                                                        kOverTkOverTOver_e_P/=T*T;
//
//                                                }
//
//
//                                                //double tr_check=(visc_info[0]-visc_info[4]-visc_info[7]-visc_info[9]);
//                                                //if (tr_check > 1e-5) {
//                                                //	std::cout << "Warning: Large deviation from tracelessness (" << tr_check << ")!\n";
//                                                //}
//
//                                                //Akk=(*shear_info)+invCoshRap( (*shear_info+4)*cosPhi*cosPhi);
//                                        }
//                                        else {
//                                                kOverTkOverTOver_e_P=0.0;
//                                        }

					
					//Photon momentum in the lab frame
					//k=mT cosh(rap)=kT cosh(rap)
					const double kLt=kT*coshRap;
					//kx=kT cos(phi)
					const double kLx=kT*cosPhi;
					//ky=kT sin(phi)
					const double kLy=kT*sinPhi;
					//kz=mT sinh(rap)=kT sinh(rap)
					const double kLz=kT*sinhRap;

					const double kLtau=kLt*coshEtaS-kLz*sinhEtaS;
					const double tau_kLeta=-1*kLt*sinhEtaS+kLz*coshEtaS;

					//kR.uR=kL.uL
					//k_rf=(k_L-\vec{u}/u0.\vec{k})/sqrt(1-u^2/u0^2)
	//				kR=gamma*(kL-betax*kLx-betay*kLy-betaz*kLz);
					const double kR=utau*kLtau-ux*kLx-uy*kLy-tau_ueta*tau_kLeta;

                                        //Our rate
                                        //dGamma(\vec{k}_L)=dGamma_0(k_rf)+(A_L)_{alpha beta} k_L^alpha k_L^beta Z(rf) 
                                        fill_grid(irap, iphi, ikT, kR, T, V4, kOverTkOverTOver_e_P, bulk_pressure, eps_plus_P, cs2, &current_rate,discSpectra[iRate]);	

                                }

			}
		}	
	}

}

//Discretized spectra: array[times][Nrap][Nphi][Npt][rates]
//Discretized spectra, version 2: array[times][Nrap][Nphi][Npt][rates][value_and_remainder]
//stPos=[tau, irap, iphi, ikT]
void fill_grid(int irap, int iphi, int ikT, double kR, double T, double V4, double kOverTkOverTOver_e_P, double bulk_pressure, double eps_plus_P, double cs2, const struct photonRate * currRate, double discSpectra[CONST_NkT][CONST_Nrap][CONST_Nphi][3]) {

	//
	double tmpRate;
	//double (*local_rate)(double, double, double);

        tmpRate=eval_photon_rate(currRate,kR/T,T,kOverTkOverTOver_e_P, bulk_pressure, eps_plus_P, cs2);

        //Cell volume: dx*dy*dz*dt=dx*dy*dEta*dTau*tau
        tmpRate*=V4;

        //Fill value
        discSpectra[ikT][irap][iphi][1]+=tmpRate;

        //Fill lower bound uncertainty
        discSpectra[ikT][irap][iphi][0]=0.0;

        //Fill upper bound uncertainty
        discSpectra[ikT][irap][iphi][2]=0.0;



}

void compute_observables(const struct photonRate rate_list[], double discSpectra[CONST_N_rates][CONST_NkT][CONST_Nrap][CONST_Nphi][3]) {

	compute_midrapidity_yield_and_vn(rate_list, discSpectra);

}

//Output the phi-integrated, rapidity-averaged-around-0 yield as a function of pT
void compute_midrapidity_yield_and_vn(const struct photonRate rate_list[], double discSpectra[CONST_N_rates][CONST_NkT][CONST_Nrap][CONST_Nphi][3]) {

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
                                                const double tmp=discSpectra[rate_no][ikT][irap][iphi][1];
						yield+=tmp*yFac;
						for(int i=1;i<=CONST_FourierNb; i++) {
							vn_cos[i-1]+=yFac*tmp*cos(i*phi);
							vn_sin[i-1]+=yFac*tmp*sin(i*phi);
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
