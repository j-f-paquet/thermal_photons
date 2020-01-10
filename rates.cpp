#include "photon.h"
#include "rates.h"

/***** Rates *****/
void init_rates(std::map<enum rate_type, struct photonRate> * rate_list, enum rate_type rate_id) {

        struct photonRate tmp_rate;
        struct photonRate * currRate=&tmp_rate;

	switch(rate_id) {

		case thermal_ideal:
			
			currRate->name="rate_thermal_ideal";
			
			currRate->is_thermal=true;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=true;
			currRate->rate_fit_function=rate_thermal_ideal;

			currRate->use_k_instead_of_kOverT_for_table=false;
			currRate->number_of_points_in_kOverT=1000;
			currRate->number_of_points_in_temp=500;
			currRate->min_temp=.1;
			currRate->max_temp=1.;
			currRate->min_kOverT=0.0;
			currRate->max_kOverT=80.0;

			tabulate_fit(currRate);

			break;

		case hadronic_ideal:
			
			currRate->name="rate_hadronic_ideal";
			
			currRate->is_thermal=true;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=true;
			currRate->rate_fit_function=rate_hadronic_ideal;

			currRate->use_k_instead_of_kOverT_for_table=false;
			currRate->number_of_points_in_kOverT=1000;
			currRate->number_of_points_in_temp=500;
			currRate->min_temp=.1;
			currRate->max_temp=1.;
			currRate->min_kOverT=0.0;
			currRate->max_kOverT=80.0;


			tabulate_fit(currRate);

			break;


		case qgp_ideal_born_AMYfit_tabulated:
			
			currRate->name="rate_qgp_ideal_born_AMYfit_tabulated";
			
			currRate->is_qgp=true;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=true;
			currRate->rate_fit_function=rate_qgp_ideal_born_AMYfit;

			currRate->use_k_instead_of_kOverT_for_table=true;
			currRate->number_of_points_in_kOverT=80;
			currRate->number_of_points_in_temp=351;
			currRate->min_temp=0.1;
			currRate->max_temp=0.8;
			currRate->min_kOverT=0.05;
			currRate->max_kOverT=4.0;

			tabulate_fit(currRate);

			break;

		case qgp_ideal_born_AMY_table:
			
			currRate->name="rate_qgp_ideal_born_AMY_table";
			
			currRate->is_qgp=true;

			currRate->use_table_instead_of_fit=true;
			currRate->tabulate_fit_for_speed=false;
			currRate->rate_fit_function=rate_qgp_ideal_born_AMYfit;


			currRate->filename_of_external_table="./partonic_rates/rate_QGP_2to2_total_eqrate.dat";
			currRate->use_k_instead_of_kOverT_for_table=true;
			currRate->number_of_points_in_kOverT=80;
			currRate->number_of_points_in_temp=351;
			currRate->min_temp=0.1;
			currRate->max_temp=0.8;
			currRate->min_kOverT=0.05;
			currRate->max_kOverT=4.0;
			load_rate_from_file(currRate);

			break;


		case qgp_ideal_born_AMYfit_with_cuts:

			currRate->name="rate_qgp_ideal_born_AMYfit_cuts";
			
			currRate->is_qgp=true;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=false;
			currRate->rate_fit_function=rate_qgp_ideal_born_AMYfit_with_cuts;

			break;
			

		//rate_qgp_ideal_born_AMYfit
		case qgp_ideal_born_AMYfit:
			
			currRate->name="rate_qgp_ideal_born_AMYfit";
			
			currRate->is_qgp=true;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=true;
			currRate->rate_fit_function=rate_qgp_ideal_born_AMYfit;

			currRate->use_k_instead_of_kOverT_for_table=false;
			currRate->number_of_points_in_kOverT=500;
			currRate->number_of_points_in_temp=300;
			currRate->min_temp=CONST_pure_HG_T;
			currRate->max_temp=1.;
			currRate->min_kOverT=0.0;
			currRate->max_kOverT=80.0;

			tabulate_fit(currRate);

			break;

		//
		case qgp_ideal_Dusling:
			
			currRate->name="rate_qgp_ideal_Dusling";
			
			currRate->is_qgp=true;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=true;
			currRate->rate_fit_function=rate_qgp_ideal_Dusling;

			currRate->use_k_instead_of_kOverT_for_table=false;
			currRate->number_of_points_in_kOverT=500;
			currRate->number_of_points_in_temp=300;
			currRate->min_temp=CONST_pure_HG_T;
			currRate->max_temp=1.;
			currRate->min_kOverT=0.0;
			currRate->max_kOverT=80.0;

			tabulate_fit(currRate);

			break;

		//
		case qgp_shear_viscous_Dusling:
			
			currRate->name="rate_qgp_shear_viscous_Dusling";
			
			currRate->is_qgp=true;
			currRate->is_shear_viscous=true;
                        currRate->ideal_rate_for_corresponding_production_channel=qgp_ideal_Dusling;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=true;
			currRate->rate_fit_function=rate_qgp_shear_viscous_Dusling;

			currRate->use_k_instead_of_kOverT_for_table=false;
			currRate->number_of_points_in_kOverT=500;
			currRate->number_of_points_in_temp=300;
			currRate->min_temp=CONST_pure_HG_T;
			currRate->max_temp=1.;
			currRate->min_kOverT=0.0;
			currRate->max_kOverT=80.0;

			tabulate_fit(currRate);

			break;

		//
		case qgp_bulk_viscous_Dusling:
			
			currRate->name="rate_qgp_bulk_viscous_Dusling";
			
			currRate->is_qgp=true;
			currRate->is_bulk_viscous=true;
                        currRate->ideal_rate_for_corresponding_production_channel=qgp_ideal_Dusling;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=true;
			currRate->rate_fit_function=rate_qgp_bulk_viscous_Dusling;

			currRate->use_k_instead_of_kOverT_for_table=false;
			currRate->number_of_points_in_kOverT=500;
			currRate->number_of_points_in_temp=300;
			currRate->min_temp=CONST_pure_HG_T;
			currRate->max_temp=1.;
			currRate->min_kOverT=0.0;
			currRate->max_kOverT=80.0;

			tabulate_fit(currRate);

			break;



		//rate_qgp_ideal_born_KLS
		case qgp_ideal_born_KLS:
			
			currRate->name="rate_qgp_ideal_born_KLS";
			
			currRate->is_qgp=true;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=false;
			currRate->rate_fit_function=rate_qgp_ideal_born_KLS;
			break;

		//rate_qgp_ideal_born_JF_sqrtg
		case qgp_ideal_born_JF_sqrtg:
			
			currRate->name="rate_qgp_ideal_born_JF_sqrtg";
			
			currRate->is_qgp=true;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=false;
			currRate->rate_fit_function=rate_qgp_ideal_born_JF_sqrtg;
			break;

		//rate_qgp_viscous_only_born_JF_sqrtg
		case qgp_viscous_only_born_JF_sqrtg:
			
			currRate->name="rate_qgp_viscous_only_born_JF_sqrtg";
			
			currRate->is_qgp=true;
			currRate->is_shear_viscous=true;
                        currRate->ideal_rate_for_corresponding_production_channel=qgp_ideal_born_JF_sqrtg;
                        //currRate->ideal_rate_for_corresponding_production_channel=rate_qgp_ideal_born_AMYfit;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=false;
			currRate->rate_fit_function=rate_qgp_viscous_only_born_JF_sqrtg;
			break;

		//rate_hg_ideal_Turbide_fit
		case hg_ideal_Turbide_fit:

			currRate->name="rate_hg_ideal_Turbide_fit";
			
			currRate->is_hg=true;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=true;
			currRate->rate_fit_function=rate_hg_ideal_Turbide_fit;

			currRate->use_k_instead_of_kOverT_for_table=false;
			currRate->number_of_points_in_kOverT=500;
			currRate->number_of_points_in_temp=250;
			currRate->min_temp=CONST_freezeout_T;
			currRate->max_temp=CONST_pure_QGP_T;
			currRate->min_kOverT=0.0;
			currRate->max_kOverT=80.0;

			tabulate_fit(currRate);
			break;

		//rate_qgp_ideal_LO_AMYfit
		case qgp_ideal_LO_AMYfit_orig:
			
			currRate->name="qgp_ideal_LO_AMYfit_orig";
			
			currRate->is_qgp=true;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=false;
			currRate->rate_fit_function=rate_qgp_ideal_LO_AMYfit;
			break;

		//rate_qgp_ideal_LO_AMYfit
		case qgp_ideal_LO_AMYfit:
			
			currRate->name="rate_qgp_ideal_LO_AMYfit";
			
			currRate->is_qgp=true;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=true;
			currRate->rate_fit_function=rate_qgp_ideal_LO_AMYfit;

			currRate->use_k_instead_of_kOverT_for_table=false;
			currRate->number_of_points_in_kOverT=1000;
			currRate->number_of_points_in_temp=500;
			currRate->min_temp=CONST_pure_HG_T;
			currRate->max_temp=1.;
			currRate->min_kOverT=0.0;
			currRate->max_kOverT=80.0;

			tabulate_fit(currRate);
			break;

		case qgp_viscous_only_born_g2_sqrtg:
			
			currRate->name="rate_qgp_viscous_only_born_g2_sqrtg";
			
			currRate->is_qgp=true;
			currRate->is_shear_viscous=true;
                        currRate->ideal_rate_for_corresponding_production_channel=qgp_ideal_born_AMYfit;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=true;
			currRate->rate_fit_function=rate_qgp_viscous_only_born_g2_sqrtg;

			currRate->use_k_instead_of_kOverT_for_table=false;
			currRate->number_of_points_in_kOverT=500;
			currRate->number_of_points_in_temp=300;
			currRate->min_temp=CONST_pure_HG_T;
			currRate->max_temp=1.;
			currRate->min_kOverT=0.0;
			currRate->max_kOverT=80.0;

			tabulate_fit(currRate);

			break;

		case qgp_viscous_only_born_g2_sqrtg_fit_tabulated:
			
			currRate->name="rate_qgp_viscous_only_born_g2_sqrtg_fit_tabulated";
			
			currRate->is_qgp=true;
			currRate->is_shear_viscous=true;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=true;
			currRate->rate_fit_function=rate_qgp_viscous_only_born_g2_sqrtg;
                        currRate->ideal_rate_for_corresponding_production_channel=qgp_ideal_born_AMYfit;

			currRate->use_k_instead_of_kOverT_for_table=true;
			currRate->number_of_points_in_kOverT=80;
			currRate->number_of_points_in_temp=351;
			currRate->min_temp=0.1;
			currRate->max_temp=0.8;
			currRate->min_kOverT=0.05;
			currRate->max_kOverT=4.0;

			tabulate_fit(currRate);

			break;

		case qgp_viscous_only_born_g2_sqrtg_table:
			
			currRate->name="rate_qgp_viscous_only_born_g2_sqrtg_table";
			
			currRate->is_qgp=true;
			currRate->is_shear_viscous=true;
                        currRate->ideal_rate_for_corresponding_production_channel=qgp_ideal_born_AMYfit;

			currRate->use_table_instead_of_fit=true;
			currRate->tabulate_fit_for_speed=false;

			currRate->filename_of_external_table="./partonic_rates/rate_QGP_2to2_total_viscous.dat";
			currRate->use_k_instead_of_kOverT_for_table=true;
			currRate->number_of_points_in_kOverT=80;
			currRate->number_of_points_in_temp=351;
			currRate->min_temp=0.1;
			currRate->max_temp=0.8;
			currRate->min_kOverT=0.05;
			currRate->max_kOverT=4.0;

			currRate->extra_normalisation_factor_function=T2_normalisation;

			load_rate_from_file(currRate);

			break;

		//rate_hg_ideal_Turbide_fit
		case hg_ideal_Turbide_fit_chem_pot_Boltz_Tch150:

			currRate->name="rate_hg_ideal_Turbide_fit_chem_pot_Boltz_Tch150";
			
			currRate->is_hg=true;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=true;
			currRate->rate_fit_function=rate_hg_ideal_Turbide_fit_chem_pot_Boltz_Tch150;

			currRate->use_k_instead_of_kOverT_for_table=false;
			currRate->number_of_points_in_kOverT=500;
			currRate->number_of_points_in_temp=250;
			currRate->min_temp=CONST_freezeout_T;
			currRate->max_temp=CONST_pure_QGP_T;
			currRate->min_kOverT=0.0;
			currRate->max_kOverT=80.0;

			tabulate_fit(currRate);
			break;

		//rate_hg_ideal_Turbide_fit
		case hg_ideal_Turbide_fit_chem_pot_Boltz_Tch160:

			currRate->name="rate_hg_ideal_Turbide_fit_chem_pot_Boltz_Tch160";
			
			currRate->is_hg=true;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=true;
			currRate->rate_fit_function=rate_hg_ideal_Turbide_fit_chem_pot_Boltz_Tch160;

			currRate->use_k_instead_of_kOverT_for_table=false;
			currRate->number_of_points_in_kOverT=500;
			currRate->number_of_points_in_temp=250;
			currRate->min_temp=CONST_freezeout_T;
			currRate->max_temp=CONST_pure_QGP_T;
			currRate->min_kOverT=0.0;
			currRate->max_kOverT=80.0;

			tabulate_fit(currRate);
			break;

		//rate_hg_ideal_Turbide_fit
		case hg_ideal_Turbide_fit_chem_pot_Boltz_Tch165:

			currRate->name="rate_hg_ideal_Turbide_fit_chem_pot_Boltz_Tch165";
			
			currRate->is_hg=true;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=true;
			currRate->rate_fit_function=rate_hg_ideal_Turbide_fit_chem_pot_Boltz_Tch165;

			currRate->use_k_instead_of_kOverT_for_table=false;
			currRate->number_of_points_in_kOverT=500;
			currRate->number_of_points_in_temp=250;
			currRate->min_temp=CONST_freezeout_T;
			currRate->max_temp=CONST_pure_QGP_T;
			currRate->min_kOverT=0.0;
			currRate->max_kOverT=80.0;

			tabulate_fit(currRate);

			break;

		case hg_ideal_Chun_table_CE:
			
			currRate->name="rate_hg_ideal_Chun_table_CE";
			
			currRate->is_hg=true;
			currRate->is_shear_viscous=false;

			currRate->use_table_instead_of_fit=true;
			currRate->tabulate_fit_for_speed=false;

			currRate->filename_of_external_table="./hadronic_rates/chun_eqrate_photon_rate_CE.dat";
			currRate->use_k_instead_of_kOverT_for_table=true;
			currRate->number_of_points_in_kOverT=81;
			currRate->number_of_points_in_temp=351;
			currRate->min_temp=0.1;
			currRate->max_temp=0.8;
			currRate->min_kOverT=0.05;
			currRate->max_kOverT=4.05;

			load_rate_from_file(currRate);

			break;

		case hg_viscous_Chun_table_CE:
			
			currRate->name="rate_hg_viscous_Chun_table_CE";
			
			currRate->is_hg=true;
			currRate->is_shear_viscous=true;
                        currRate->ideal_rate_for_corresponding_production_channel=hg_ideal_Chun_table_CE;

			currRate->use_table_instead_of_fit=true;
			currRate->tabulate_fit_for_speed=false;

			currRate->filename_of_external_table="./hadronic_rates/chun_shear_photon_rate_CE.dat";
			currRate->use_k_instead_of_kOverT_for_table=true;
			currRate->number_of_points_in_kOverT=81;
			currRate->number_of_points_in_temp=351;
			currRate->min_temp=0.1;
			currRate->max_temp=0.8;
			currRate->min_kOverT=0.05;
			currRate->max_kOverT=4.05;

			currRate->extra_normalisation_factor_function=T2_normalisation;

			load_rate_from_file(currRate);

			break;

		case hg_ideal_Chun_table_PCE165:
			
			currRate->name="rate_hg_ideal_Chun_table_PCE165";
			
			currRate->is_hg=true;
			currRate->is_shear_viscous=false;

			currRate->use_table_instead_of_fit=true;
			currRate->tabulate_fit_for_speed=false;

			currRate->filename_of_external_table="./hadronic_rates/rate_HG_2to2_total_eqrate_PCE165.dat";
			currRate->use_k_instead_of_kOverT_for_table=true;
			currRate->number_of_points_in_kOverT=80;
			currRate->number_of_points_in_temp=351;
			currRate->min_temp=0.1;
			currRate->max_temp=0.8;
			currRate->min_kOverT=0.05;
			currRate->max_kOverT=4.0;

			load_rate_from_file(currRate);

			break;

		case hg_viscous_Chun_table_PCE165:
			
			currRate->name="rate_hg_viscous_Chun_table_PCE165";
			
			currRate->is_hg=true;
			currRate->is_shear_viscous=true;
                        currRate->ideal_rate_for_corresponding_production_channel=hg_ideal_Chun_table_PCE165;

			currRate->use_table_instead_of_fit=true;
			currRate->tabulate_fit_for_speed=false;

			currRate->filename_of_external_table="./partonic_rates/rate_HG_2to2_total_viscous_PCE165.dat";
			currRate->use_k_instead_of_kOverT_for_table=true;
			currRate->number_of_points_in_kOverT=80;
			currRate->number_of_points_in_temp=351;
			currRate->min_temp=0.1;
			currRate->max_temp=0.8;
			currRate->min_kOverT=0.05;
			currRate->max_kOverT=4.0;

			currRate->extra_normalisation_factor_function=T2_normalisation;

			load_rate_from_file(currRate);

			break;

		case hg_bulk_Chun_table_CE:
			
			currRate->name="rate_hg_bulk_Chun_table_CE";
			
			currRate->is_hg=true;
			currRate->is_shear_viscous=false;
			currRate->is_bulk_viscous=true;
                        currRate->ideal_rate_for_corresponding_production_channel=hg_ideal_Chun_table_CE;

			currRate->use_table_instead_of_fit=true;
			currRate->tabulate_fit_for_speed=false;

			currRate->filename_of_external_table="./hadronic_rates/chun_bulk_photon_rate_CE.dat";
			currRate->use_k_instead_of_kOverT_for_table=true;
			currRate->number_of_points_in_kOverT=81;
			currRate->number_of_points_in_temp=75;
			currRate->min_temp=0.1;
			currRate->max_temp=0.248;
			currRate->min_kOverT=0.05;
			currRate->max_kOverT=4.05;

			//currRate->extra_normalisation_factor_function=minus_sign_normalisation;

			load_rate_from_file(currRate);

			break;

		case hg_bulk_Chun_table_CE_eos_transport:
			
			currRate->name="rate_hg_bulk_Chun_table_CE_eos_transport";
			
			currRate->is_hg=true;
			currRate->is_shear_viscous=false;
			currRate->is_bulk_viscous=true;
                        currRate->ideal_rate_for_corresponding_production_channel=hg_ideal_Chun_table_CE;

			currRate->use_table_instead_of_fit=true;
			currRate->tabulate_fit_for_speed=false;

			currRate->filename_of_external_table="./hadronic_rates/chun_bulk_photon_rate_CE_eos_transport.dat";
			currRate->use_k_instead_of_kOverT_for_table=true;
			currRate->number_of_points_in_kOverT=81;
			currRate->number_of_points_in_temp=75;
			currRate->min_temp=0.1;
			currRate->max_temp=0.248;
			currRate->min_kOverT=0.05;
			currRate->max_kOverT=4.05;

			//currRate->extra_normalisation_factor_function=minus_sign_normalisation;

			load_rate_from_file(currRate);

			break;

//		case hg_bulk_Chun_table_PCE165:
//			
//			currRate->name="rate_hg_bulk_Chun_table_PCE165";
//			
//			currRate->is_hg=true;
//			currRate->is_shear_viscous=false;
//			currRate->is_bulk_viscous=true;
//
//			currRate->use_table_instead_of_fit=true;
//			currRate->tabulate_fit_for_speed=false;
//
//			currRate->filename_of_external_table="./partonic_rates/chun_bulk_photon_rate_PCE165.dat";
//			currRate->use_k_instead_of_kOverT_for_table=true;
//			currRate->number_of_points_in_kOverT=81;
//			currRate->number_of_points_in_temp=76;
//			currRate->min_temp=0.1;
//			currRate->max_temp=0.25;
//			currRate->min_kOverT=0.05;
//			currRate->max_kOverT=4.05;
//
//			//currRate->extra_normalisation_factor_function=T2_normalisation;
//
//			load_rate_from_file(currRate);
//
//			break;

		//rate_hg_ideal_Turbide_fit
		case hg_ideal_Turbide_fit_tabulated:

			currRate->name="rate_hg_ideal_Turbide_fit_tabulated";

			currRate->is_hg=true;
			currRate->is_shear_viscous=false;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=true;
			currRate->rate_fit_function=rate_hg_ideal_Turbide_fit;

			currRate->use_k_instead_of_kOverT_for_table=true;
			currRate->number_of_points_in_kOverT=80;
			currRate->number_of_points_in_temp=351;
			currRate->min_temp=0.1;
			currRate->max_temp=0.8;
			currRate->min_kOverT=0.05;
			currRate->max_kOverT=4.0;

			tabulate_fit(currRate);
			break;


		//rate_hg_ideal_Turbide_fit
		case hg_ideal_Turbide_fit_noPiPi_tabulated:

			currRate->name="rate_hg_ideal_Turbide_fit_noPiPi_tabulated";

			currRate->is_hg=true;
			currRate->is_shear_viscous=false;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=true;
			currRate->rate_fit_function=rate_hg_ideal_Turbide_noPiPi_fit;

			currRate->use_k_instead_of_kOverT_for_table=true;
			currRate->number_of_points_in_kOverT=80;
			currRate->number_of_points_in_temp=351;
			currRate->min_temp=0.1;
			currRate->max_temp=0.8;
			currRate->min_kOverT=0.05;
			currRate->max_kOverT=4.0;

			tabulate_fit(currRate);
			break;

		//rate_hg_ideal_Turbide_fit
		case hg_in_medium_rho_ideal_Rapp_fit_tabulated:

			currRate->name="rate_hg_in_medium_rho_ideal_Rapp_fit_tabulated";

			currRate->is_hg=true;
			currRate->is_shear_viscous=false;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=true;
			currRate->rate_fit_function=rate_hg_in_medium_rho_ideal_Rapp_fit;

			currRate->use_k_instead_of_kOverT_for_table=true;
			currRate->number_of_points_in_kOverT=150;
			currRate->number_of_points_in_temp=351;
			currRate->min_temp=0.1;
			currRate->max_temp=0.22;
			currRate->min_kOverT=0.05;
			currRate->max_kOverT=4.0;

			tabulate_fit(currRate);
			break;




		//rate_hg_ideal_Turbide_fit
		case hg_piRhoOmega_ideal_Rapp_fit_tabulated:

			currRate->name="rate_hg_in_piRhoOmega_ideal_Rapp_fit_tabulated";

			currRate->is_hg=true;
			currRate->is_shear_viscous=false;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=true;
			currRate->rate_fit_function=rate_hg_piRhoOmega_ideal_Rapp_fit;

			currRate->use_k_instead_of_kOverT_for_table=true;
			currRate->number_of_points_in_kOverT=150;
			currRate->number_of_points_in_temp=351;
			currRate->min_temp=0.1;
			currRate->max_temp=0.22;
			currRate->min_kOverT=0.05;
			currRate->max_kOverT=5.0;

			tabulate_fit(currRate);
			break;



		//rate pions bremstralung
		case hg_pion_brem_ideal_Rapp_fit_tabulated:

			currRate->name="rate_hg_pion_brem_ideal_Rapp_fit_tabulated";

			currRate->is_hg=true;
			currRate->is_shear_viscous=false;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=true;
			currRate->rate_fit_function=rate_hg_pion_brem_ideal_Rapp_fit;

			currRate->use_k_instead_of_kOverT_for_table=true;
			currRate->number_of_points_in_kOverT=150;
			currRate->number_of_points_in_temp=351;
			currRate->min_temp=0.1;
			currRate->max_temp=0.22;
			currRate->min_kOverT=0.05;
			currRate->max_kOverT=4.0;

			tabulate_fit(currRate);
			break;


		//rate pions bremstralung
		case hg_ideal_Zahed_Dusling_2pi_fit_tabulated:

			currRate->name="rate_hg_ideal_Zahed_Dusling_2pi_fit_tabulated";

			currRate->is_hg=true;
			currRate->is_shear_viscous=false;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=true;
			currRate->rate_fit_function=rate_hg_ideal_Zahed_Dusling_2pi_fit;

			currRate->use_k_instead_of_kOverT_for_table=false;
			currRate->number_of_points_in_kOverT=500;
			currRate->number_of_points_in_temp=250;
			currRate->min_temp=CONST_freezeout_T;
			currRate->max_temp=CONST_pure_QGP_T;
			currRate->min_kOverT=0.0;
			currRate->max_kOverT=80.0;

			tabulate_fit(currRate);


			break;

		case qgp_viscous_only_born_g2_sqrtg_new_deltaf:
			
			currRate->name="rate_qgp_viscous_only_born_g2_sqrtg_new_deltaf";
			
			currRate->is_qgp=true;
			currRate->is_shear_viscous=true;
                        currRate->ideal_rate_for_corresponding_production_channel=qgp_ideal_born_AMYfit;

			currRate->use_table_instead_of_fit=false;
			currRate->tabulate_fit_for_speed=false;
			currRate->rate_fit_function=rate_qgp_viscous_only_born_g2_sqrtg_new_deltaf;
			break;

	}

        (*rate_list)[rate_id] = tmp_rate;

}

void tabulate_fit(struct photonRate * currRate) {

	void get_photon_rate(int selector, double (**local_rate)(double, double, double));
	double interp_x_from_index(double xmin, double xmax, int i, int size); 

	if (currRate->tabulate_fit_for_speed) {

		if (currRate->rate_fit_function == 0) {
			std::cout << "Fit function undefined, can't tabulate it! Aborting...\n";
			exit(1);
		}
		else {

                        std::cout << "Tabulating rate " << currRate->name << "\n";

			const int size_x=currRate->number_of_points_in_kOverT;
			const int size_y=currRate->number_of_points_in_temp;

			currRate->tabulated_rate = new double * [size_y];	
			currRate->tabulated_rate_log = new double * [size_y];	

			for(int k=0; k<size_y;k++) {
				const double temp=interp_x_from_index(currRate->min_temp,currRate->max_temp,k,currRate->number_of_points_in_temp);
				currRate->tabulated_rate[k] = new double [size_x];	
				currRate->tabulated_rate_log[k] = new double [size_x];	
				for(int j=0; j<size_x;j++) {
					const double ku=interp_x_from_index(currRate->min_kOverT,currRate->max_kOverT,j,currRate->number_of_points_in_kOverT);
					//const double temp=temp_from_index(currRate,k);
					//const double ku=kOverT_from_index(currRate,j);
					double kOverT;
					//If the table is tabulated with k, ku is k, not k/T
					if (currRate->use_k_instead_of_kOverT_for_table) {
						kOverT=ku/temp;
					}
					else {
						kOverT=ku;
					}
					const double tmp_res=(*(currRate->rate_fit_function))(kOverT,temp,0.0);
					currRate->tabulated_rate[k][j]=tmp_res;
					currRate->tabulated_rate_log[k][j]=log(tmp_res);
				}
			}

		}

	}
}

//Assume that each line correspond to a temperature and each column to a photon energy (or photon energy over T)
void load_rate_from_file(struct photonRate * currRate) {

	float tmpRate;
	std::ifstream rateFile;
	rateFile.open(currRate->filename_of_external_table.c_str(),std::ios::in);

	//Check if the file exists
	if (!rateFile.good()) {
		std::cout << "Can't read \"" << currRate->filename_of_external_table.c_str() << "\"... Aborting.\n";
		exit(1);
	}
	else  {

		const int size_x=currRate->number_of_points_in_kOverT;
		const int size_y=currRate->number_of_points_in_temp;

		currRate->tabulated_rate = new double * [size_y];	
		currRate->tabulated_rate_log = new double * [size_y];	

		for(int k=0; k<size_y;k++) {
			currRate->tabulated_rate[k] = new double [size_x];	
			currRate->tabulated_rate_log[k] = new double [size_x];	
			for(int j=0; j<size_x;j++) {
				rateFile >> tmpRate;
				currRate->tabulated_rate[k][j]=tmpRate;
				currRate->tabulated_rate_log[k][j]=log(tmpRate);
			}
		}

		rateFile.close();
	}
}

double eval_photon_rate(std::map<enum rate_type, struct photonRate> * rate_list, enum rate_type rate_id, double kOverT, double T, double kOverTkOverTOver_e_P, double bulk_pressure, double eps_plus_P, double cs2) {

	//void get_photon_rate(int selector, double (**local_rate)(double, double, double));
        double get_photon_rate_accel(const struct photonRate * currRate, double kOverT, double T, double kk);
	double QGP_fraction(double T); 

	//
	double res=1.0;

        struct photonRate * currRate=&((*rate_list)[rate_id]);

	if (!currRate->is_thermal) {
	//For speed, check if a QGP fraction is used first
		if (currRate->is_qgp) {
			res=QGP_fraction(T); 

		}
		else if (currRate->is_hg) {
			res=1-QGP_fraction(T);
		}
	}

	//Multiply by shear viscosity factor?
	if (currRate->is_shear_viscous) {
		res*=kOverTkOverTOver_e_P/2.0;
	}
	//Multiply by bulk pressure?
	else if (currRate->is_bulk_viscous) {
		//Multiply QGP rate by 1/(1/3-cs^2)*(bulk pressure)/(epsilon+P) 
		if (currRate->is_qgp) {
			res*=bulk_pressure/eps_plus_P/(1./3.0-cs2);
		}
		//and hadron gas rate by just the bulk pressure
		else if (currRate->is_hg) {
			res*=bulk_pressure;
		}
	}

	//Only compute the rate, which is the slowest part of the calculation, if all the above factors are non-zero
	if (res != 0.0) {

		//For tabulated rate fit or for external rate tables 
		if ((currRate->tabulate_fit_for_speed)||(currRate->use_table_instead_of_fit)) {
			//if the table is w.r.t. k instead of k/T, we have to retrieve the correct value
			double ku;
			if (currRate->use_k_instead_of_kOverT_for_table) {
				ku=kOverT*T;
			}
			else {
				ku=kOverT;
			}
			res*=get_photon_rate_accel(currRate, ku, T, kOverTkOverTOver_e_P);
		}
		//Use plain fit
		else {
			res*=(*(currRate->rate_fit_function))(kOverT,T,kOverTkOverTOver_e_P);
		}


		if (currRate->extra_normalisation_factor_function != 0) {
			res*=(*(currRate->extra_normalisation_factor_function))(kOverT,T,kOverTkOverTOver_e_P);
		}


                // Regulation of viscous corrections
                if (((currRate->is_shear_viscous)||(currRate->is_bulk_viscous))&&(regulate_negative_rate)) {
                
                        // Get the corresponding ideal rate...
                        double ideal_rate=eval_photon_rate(rate_list, currRate->ideal_rate_for_corresponding_production_channel, kOverT, T, kOverTkOverTOver_e_P, bulk_pressure, eps_plus_P, cs2); 

                        // Cap the viscous rate to (viscous_rate_over_ideal_smaller_than*ideal_rate)
                        if (fabs(res)>viscous_rate_over_ideal_smaller_than*fabs(ideal_rate)) {
                                res=(res > 0 ? 1 : -1)*viscous_rate_over_ideal_smaller_than*fabs(ideal_rate);
                        }

                }

	}

	return res;
}

double interp_x_from_index(double xmin, double xmax, int i, int size) {

	double res=xmin+(xmax-xmin)*i*1.0/(size-1);

	return res;

}

int interp_index_from_x(double xmin, double xmax, double x, int size) {

	const int index=floor((x-xmin)*(size-1)*1.0/(xmax-xmin)+0.5);

	return index;

}

//double kOverT_from_index(const struct photonRate * currRate, int i) {
//	const double k_min=currRate->min_kOverT;
//	const double k_max=currRate->max_kOverT;
//	const int size=currRate->number_of_points_in_kOverT;
//	//return 80*i*i*1.0/(size*size);
//	return k_min+(k_max-k_min)*i*i*1.0/(size*size);
//}
//
//int index_from_kOverT(const struct photonRate * currRate, double kOverT) {
//	const double k_min=currRate->min_kOverT;
//	const double k_max=currRate->max_kOverT;
//	const int size=currRate->number_of_points_in_kOverT;
//	return floor(sqrt((kOverT-k_min)*size*size*1.0/(k_max-k_min)));
//	//return floor(sqrt(kOverT*size*size*1.0/80.0));
//}
//
//double temp_from_index(const struct photonRate * currRate, int i) {
//	const double t_min=currRate->min_temp;
//	const double t_max=currRate->max_temp;
//	const int size=currRate->number_of_points_in_temp;
//	return t_min+(t_max-t_min)*i*i*1.0/(size*size);
//}
//
//int index_from_temp(const struct photonRate * currRate, double temp) {
//	const double t_min=currRate->min_temp;
//	const double t_max=currRate->max_temp;
//	const int size=currRate->number_of_points_in_temp;
//	return floor(sqrt((temp-t_min)*size*size*1.0/(t_max-t_min)));
//}

double get_photon_rate_accel(const struct photonRate * currRate, double kOverT, double T, double kk) {

	double interp_x_from_index(double xmin, double xmax, int i, int size);
	int interp_index_from_x(double xmin, double xmax, double x, int size); 

	double res;

	const int size_x=currRate->number_of_points_in_kOverT;
	const int size_y=currRate->number_of_points_in_temp;

	//int a1=index_from_kOverT(currRate,kOverT);
	//int b1=index_from_temp(currRate,T);		
	int a1=interp_index_from_x(currRate->min_kOverT,currRate->max_kOverT,kOverT,currRate->number_of_points_in_kOverT);
	int b1=interp_index_from_x(currRate->min_temp,currRate->max_temp,T,currRate->number_of_points_in_temp);

	//Must decide what to do when the requested value is outside the table's range
	if ((a1+1>=size_x)||(b1+1>=size_y)||(a1<0)||(b1<0)) {
		res=0.0;
	}
	else {

		double fx1y1=currRate->tabulated_rate[b1][a1];
		double fx2y1=currRate->tabulated_rate[b1][a1+1];
		double fx1y2=currRate->tabulated_rate[b1+1][a1];
		double fx2y2=currRate->tabulated_rate[b1+1][a1+1];

		//double x1=kOverT_from_index(currRate,a1);
		//double x2=kOverT_from_index(currRate,a1+1);
		//double y1=temp_from_index(currRate,b1);
		//double y2=temp_from_index(currRate,b1+1);
		double x1=interp_x_from_index(currRate->min_kOverT,currRate->max_kOverT,a1,currRate->number_of_points_in_kOverT);
		double x2=interp_x_from_index(currRate->min_kOverT,currRate->max_kOverT,a1+1,currRate->number_of_points_in_kOverT);
		double y1=interp_x_from_index(currRate->min_temp,currRate->max_temp,b1,currRate->number_of_points_in_temp);
		double y2=interp_x_from_index(currRate->min_temp,currRate->max_temp,b1+1,currRate->number_of_points_in_temp);

		//Robust but inefficient
		if (fx1y1>0&&fx2y1>0&&fx1y2>0&&fx2y2>0) {

			//Exact version
			//double log11=log(fx1y1);
			//double log12=log(fx1y2);
			//double log21=log(fx2y1);
			//double log22=log(fx2y2);
			double log11=currRate->tabulated_rate_log[b1][a1];
			double log12=currRate->tabulated_rate_log[b1+1][a1];
			double log21=currRate->tabulated_rate_log[b1][a1+1];
			double log22=currRate->tabulated_rate_log[b1+1][a1+1];

			//Slighly faster version
			//Log[a + x] = Log[a] + x/a - x^2/(2 a^2) + x^3/(3 a^3) + ...
			//PadeApproximant[Log[ap + R ap], {R, 0, 2}]=(Log[ap]+R (1+Log[ap])+1/6 R^2 (3+Log[ap]))/(1+R+R^2/6)
			//double log11=log(fx1y1);
			//double ratio=(fx1y2-fx1y1)/fx1y1;
			//double log12=log11+ratio-ratio*ratio*0.5+1./3.*ratio*ratio*ratio;
			//ratio=(fx2y1-fx1y1)/fx1y1;
			//double log21=log11+ratio-ratio*ratio*0.5+1./3.*ratio*ratio*ratio;
			//ratio=(fx2y2-fx1y1)/fx1y1;
			//double log22=log11+ratio-ratio*ratio*0.5+1./3.*ratio*ratio*ratio;
		
			res=exp((log11*(x2-kOverT)*(y2-T)+log21*(kOverT-x1)*(y2-T)+log12*(x2-kOverT)*(T-y1)+log22*(kOverT-x1)*(T-y1))/((x2-x1)*(y2-y1)));
			//res=exp((fx1y1*(x2-kOverT)*(y2-T)+fx2y1*(kOverT-x1)*(y2-T)+fx1y2*(x2-kOverT)*(T-y1)+fx2y2*(kOverT-x1)*(T-y1))/((x2-x1)*(y2-y1)));
			//res=exp((log(fx1y1)*(x2-kOverT)*(y2-T)+log(fx2y1)*(kOverT-x1)*(y2-T)+log(fx1y2)*(x2-kOverT)*(T-y1)+log(fx2y2)*(kOverT-x1)*(T-y1))/((x2-x1)*(y2-y1)));

		}
		else {
			res=(fx1y1*(x2-kOverT)*(y2-T)+fx2y1*(kOverT-x1)*(y2-T)+fx1y2*(x2-kOverT)*(T-y1)+fx2y2*(kOverT-x1)*(T-y1))/((x2-x1)*(y2-y1));
		}
	}
	
	return res;

}

//Template for rate E d^3 Gamma/d k^3
//last argument can be dummy in ideal case
//double rate_template(double kOverT, double T, double kkPiOver_e_P_k2) {
//
//	//Compute the whole rate
//
//	//Return it
//
//}

void validate_rates(std::map<enum rate_type, struct photonRate> * rate_list) {

        // Loop over all the rate objects
        for(std::map<enum rate_type, struct photonRate>::iterator it = rate_list->begin(); it != rate_list->end(); it++ ) {

                // Get one rate object
                struct photonRate currRate=it->second;

                // Only if the rate is viscous and regulation is on...
                if (((currRate.is_shear_viscous)||(currRate.is_bulk_viscous))&&(regulate_negative_rate)) {
                        enum rate_type ideal_rate_id=currRate.ideal_rate_for_corresponding_production_channel;
                        // Check if the corresponding ideal rate is initialized...
                        if (rate_list->count(ideal_rate_id)==0) {
                                // If not, initialize it
                                init_rates(rate_list, ideal_rate_id);
                        }
                }
        }

}


double QGP_fraction(double T) {

	double res;

	//
	if (T>=CONST_pure_QGP_T) {
		res=1.0;
	}
	else if (T>=CONST_pure_HG_T) {
		//Linear interpolation?
		res=(T-CONST_pure_HG_T)/(CONST_pure_QGP_T-CONST_pure_HG_T);
	}
	else {
		res=0.0;
	}

	return res;

}


//Pre-factor A(k)
double prefA(double kOverT, double T) {

	//Forward declaration of functions
	//inline double fermiDirac(double kOverT);
	double fermiDirac(double kOverT);

	//Nf, dF, 
	const double qCharge2[]={4.0/9.0,5.0/9.0,2.0/3.0,10.0/9.0,11.0/9.0,5.0/3.0};

	//
	double res;

	res=2*CONST_alphaEM*CONST_Nc*qCharge2[CONST_Nf-1]*CONST_mInfOverT*CONST_mInfOverT*T*T*fermiDirac(kOverT)/kOverT;

	return res;

}

//QGP ideal rate - Born - AMY fit [arXiv:hep-ph/0111107, section I]
double rate_qgp_ideal_born_AMYfit(double kOverT, double T, double kkPiOver_e_P_k2) {

	//Forward declaration
	double C_hard(double kOverT);

	double res= CONST_GeV2_to_GeVm2_fmm4*kOverT*prefA(kOverT,T)/CONST_twoPiCubed*(log(1/CONST_mInfOverT)+C_hard(kOverT));

	return res;

}

//QGP ideal rate - Born - AMY fit [arXiv:hep-ph/0111107, section I]
double rate_qgp_ideal_born_AMYfit_with_cuts(double kOverT, double T, double kkPiOver_e_P_k2) {

	//Forward declaration
	double C_hard(double kOverT);
	double res;

	if ((kOverT*T < 0.05)||(kOverT*T>4.0)||(T<0.1)||(T>0.8)) {
		res=0.0;
	} 
	else {
		res=CONST_GeV2_to_GeVm2_fmm4*kOverT*prefA(kOverT,T)/CONST_twoPiCubed*(log(1/CONST_mInfOverT)+C_hard(kOverT));
	}

	return res;

}

//QGP ideal rate - LO - AMY fit [arXiv:hep-ph/0111107, section I]
double rate_qgp_ideal_LO_AMYfit(double kOverT, double T, double kkPiOver_e_P_k2) {

	//Forward declaration
	double C_hard(double kOverT);
	double C_LPM(double kOverT);

	double res= CONST_GeV2_to_GeVm2_fmm4*kOverT*prefA(kOverT,T)/CONST_twoPiCubed*(log(1/CONST_mInfOverT)+C_hard(kOverT)+C_LPM(kOverT));

	return res;

}

//QGP ideal rate - Dusling [arXiv:0903.1764, eq.5]
double rate_qgp_ideal_Dusling(double kOverT, double T, double kkPiOver_e_P_k2) {

	//Forward declaration
	double fermiDirac(double kOverT);

	//
	const double qCharge2[]={4.0/9.0,5.0/9.0,2.0/3.0,10.0/9.0,11.0/9.0,5.0/3.0};

	double res=CONST_GeV2_to_GeVm2_fmm4*1./(2.0*M_PI*M_PI)*qCharge2[CONST_Nf-1]*CONST_alphaEM*CONST_alphaS*fermiDirac(kOverT)*T*T*log(3.7388*kOverT*1.0/(CONST_gs*CONST_gs));

	return res;

}

//Dusling's QGP shear rate: df_shear/f*(Dusling's ideal rate)
//df_shear=K_mu K_nu Pi^munu/(2*(eps+P)*T^2)*(1-f(k))*f(k) 
//A factor K_mu K_nu Pi^munu/(T^2*(eps+P))/2.0 is included in function "eval_photon_rate()" when shear viscosity flag is on
//This function should thus be (1-f(k))*(ideal rate)
double rate_qgp_shear_viscous_Dusling(double kOverT, double T, double kkPiOver_e_P_k2) {

	//Forward declaration
	double fermiDirac(double kOverT);
	double rate_qgp_ideal_Dusling(double kOverT, double T, double kkPiOver_e_P_k2);


	//
	double res=(1.0-fermiDirac(kOverT))*rate_qgp_ideal_Dusling(kOverT,T,0.0);

	return res;


}

//Dusling's QGP bulk rate: df_bulk/f*(Dusling's ideal rate)
//df_bulk=-f(k)*(1-f(k))*(m^2/T^2/(k/T)-k/T)*1/(1/3-c_s^2)*bulk_pressure/15./(eps+P)
//A factor bulk_pressure/eps_plus_P/(1./3.0-cs2) is included in function "eval_photon_rate()" when bulk viscosity flag is on
//This function should thus be -(1-f(k))*(m^2/T^2/(k/T)-k/T)/15.0*(ideal rate)
double rate_qgp_bulk_viscous_Dusling(double kOverT, double T, double kkPiOver_e_P_k2) {

	//Forward declaration
	double fermiDirac(double kOverT);
	double rate_qgp_ideal_Dusling(double kOverT, double T, double kkPiOver_e_P_k2);

	const double EOverT=sqrt(kOverT*kOverT+CONST_mInfOverT*CONST_mInfOverT);


	//
	//double res=-1.0*(1.0-fermiDirac(kOverT))*(CONST_mInfOverT*CONST_mInfOverT/kOverT-kOverT)/15.0*rate_qgp_ideal_Dusling(kOverT,T,0.0);
	double res=-1.0*(1.0-fermiDirac(kOverT))*(CONST_mInfOverT*CONST_mInfOverT/EOverT-EOverT)/15.0*rate_qgp_ideal_Dusling(kOverT,T,0.0);

	return res;


}

//QGP ideal rate - KLS high k/T, low g formula [Kapusta et al, PRD44, 9 (1991), eq.41]
double rate_qgp_ideal_born_KLS(double kOverT, double T, double kkPiOver_e_P_k2) {

	const double qCharge2[]={4.0/9.0,5.0/9.0,2.0/3.0,10.0/9.0,11.0/9.0,5.0/3.0};

	double res= CONST_GeV2_to_GeVm2_fmm4*qCharge2[CONST_Nf-1]*CONST_alphaEM*CONST_alphaS/(2.0*M_PI*M_PI)*T*T*exp(-kOverT)*log(2.912*kOverT/(CONST_gs*CONST_gs));

	return res;

}

//QGP ideal rate - JF fit - q^*=sqrt(g)
double rate_qgp_ideal_born_JF_sqrtg(double kOverT, double T, double kkPiOver_e_P_k2) {

	double res=CONST_GeV2_to_GeVm2_fmm4*kOverT*prefA(kOverT,T)/CONST_twoPiCubed*(0.8452052719374467 + 0.06345436545672481*kOverT + 0.20266453593373313*kOverT*kOverT + 0.007103855524696941*kOverT*kOverT*kOverT)/(1 + 0.3137709585719375*kOverT + 0.12623968017081683*kOverT*kOverT + 0.0021744062978126125*kOverT*kOverT*kOverT);

	return res;

}

//viscous correction to rate: A_\alpha\beta K^\alpha K^\beta/k^2 k A(k)/(2 pi)^3 * viscous_correction_born_JF_sqrtg()
double rate_qgp_viscous_only_born_JF_sqrtg(double kOverT, double T, double kkPiOver_e_P_k2) {

	double res = 1.0/(kOverT*kOverT)*CONST_GeV2_to_GeVm2_fmm4*kOverT*prefA(kOverT,T)/CONST_twoPiCubed*exp(-0.5041041126181884 + (-0.5335015716121183 + 1.9967643068761307*kOverT - 0.5616138941792664*kOverT*kOverT - 0.0009120108228910325*kOverT*kOverT*kOverT)/(1 - 2.607918425474197*kOverT - 0.8369709712322181*kOverT*kOverT))*pow(kOverT,2.1309931380115588);
	
	return res;

}

//viscous correction to rate
//q^*=sqrt(g_s), g_s=2 and 0.23<k/T<50
//Set to 0 outside its validity range
double rate_qgp_viscous_only_born_g2_sqrtg(double kOverT, double T, double kkPiOver_e_P_k2) {
	
	const double kOverT_min=.23, kOverT_max=55.;

	double res = 0.0;

	if ((kOverT>kOverT_min)&&(kOverT<kOverT_max)) res=1.0/(kOverT*kOverT)*CONST_GeV2_to_GeVm2_fmm4*kOverT*prefA(kOverT,T)/CONST_twoPiCubed*exp(0.295523 +(122.469 -713.728*kOverT+1034.67*kOverT*kOverT-1437.68*kOverT*kOverT*kOverT-650.401*kOverT*kOverT*kOverT*kOverT+562.913*kOverT*kOverT*kOverT*kOverT*kOverT+1.86032*kOverT*kOverT*kOverT*kOverT*kOverT*kOverT)/(1-155.427*kOverT+794.21*kOverT*kOverT-743.499*kOverT*kOverT*kOverT+910.294*kOverT*kOverT*kOverT*kOverT+109.401*kOverT*kOverT*kOverT*kOverT*kOverT))*pow(kOverT,0.797818);
	
	return res;

}

////viscous correction to rate
////q^*=sqrt(g_s), g_s=2
////Set to 0 outside its validity range
//double rate_qgp_viscous_only_born_g2_sqrtg_fit2(double kOverT, double T, double kkPiOver_e_P_k2) {
//	
//	const double kOverT_min=.23, kOverT_max=56.23;
//
//	double res = 0.0;
//
//	if ((kOverT>kOverT_min)&&(kOverT<kOverT_max)) res=1.0/(kOverT*kOverT)*CONST_GeV2_to_GeVm2_fmm4*kOverT*prefA(kOverT,T)/CONST_twoPiCubed*exp(-0.6931751991916247 + (7.87763445794781 - 42.411575069925824*kOverT + 48.51030946174481*pow(kOverT,2) - 65.86346604276692*pow(kOverT,3) - 100.25264753957238*pow(kOverT,4) + 58.71186065436014*pow(kOverT,5) + 61.16364489051317*pow(kOverT,6) + 0.1738872251791384*pow(kOverT,7))/(1 - 16.38152525372209*kOverT + 41.463359273491264*pow(kOverT,2) + 46.0077555440831*pow(kOverT,3) + 5.826535120681859*pow(kOverT,4) + 83.13189143375328*pow(kOverT,5) + 9.484958438301062*pow(kOverT,6)))*pow(kOverT,0.7128199956159046);
//	
//	return res;
//
//}



//viscous correction to rate
//Assume \delta f = f_{eq}*(1+a*f_{eq})*g(u_{\mu}k^{\mu}/T)*pi_{\mu\nu}/(e+P)/2*k^{\mu}*k^{\nu} /(T**2)
//with g(x)=122.2408622*(0.06892336-0.02450075*x+0.05561067*x*x+0.0016648148*x*x*x)/(1+x)**2/(0.1+x)
//q^*=sqrt(g_s), g_s=2 and 0.1<k/T<54.95
//Set to 0 outside its validity range
double rate_qgp_viscous_only_born_g2_sqrtg_new_deltaf(double kOverT, double T, double kkPiOver_e_P_k2) {
	
	const double kOverT_min=.1, kOverT_max=54.95;

	double res = 0.0;

	if ((kOverT>kOverT_min)&&(kOverT<kOverT_max)) res=1.0/(kOverT*kOverT)*CONST_GeV2_to_GeVm2_fmm4*kOverT*prefA(kOverT,T)/CONST_twoPiCubed*exp(1.4274971976983575 + (-0.8991812638895219 - 1.5871546798742364*kOverT - 0.8076355682623856*pow(kOverT,2) + 0.16857175093784463*pow(kOverT,3) + 0.0014374161291589363*pow(kOverT,4))/(1 - 0.03663738034947868*kOverT + 0.8338180260100825*pow(kOverT,2) + 0.12919671811881156*pow(kOverT,3)))*pow(kOverT,1.1322592873421042);
	
	return res;

}

////QGP viscous rate
//double factor_born_viscous(double kOverT) {
//
//
//
//}

//AMY fits (from arXiv:1302.5970, see also arXiv:hep-ph/0111107)
//Valid for k/T>0.2
double C_hard(double kOverT) {

	return 1/2.0*log(2*kOverT)+0.041/kOverT-0.3615+1.01*exp(-1.35*kOverT);

}

//AMY fits (from arXiv:1302.5970, see also arXiv:hep-ph/0111107)
//Valid for 0.2<k/T<50
double C_LPM(double kOverT) {

	return sqrt(1+CONST_Nf/6.0)*(0.548*log(12.28+1.0/kOverT)*pow(kOverT,-3.0/2.0)+0.133*kOverT/sqrt(1+kOverT/16.27));

}

////Valid for 0.2<k/T<50
//double C_coll(double kOverT) {
//
//	return sqrt(1+CONST_Nf/6.0)*(0.548*log(12.28+1/kOverT)/pow(kOverT,3./2.)+0.133*kOverT/sqrt(1+kOverT/16.27));
//
//}
//
//
////
//double rate_qgp_Born_viscous(double kROverT, double kkPiOverEta, double T, double * remainder) {
//
//	//Limits of the fit
//	const double kROverT_min=0.2;
//	const double kROverT_max=50.0;
//
//	double res=0.0, rem=0.0;
//
//	//
//	if (kROverT > kROverT_min) {
//		//Inside validity range
//		if (kROverT < kROverT_max) {
//			//res=...
//		}
//		//Above validity range
//		else {
//			//rem=
//		}
//	}
//	//Below validity range
//	else {
//		//rem=
//	}
//
//	//Return values
//	*remainder=rem;
//	return res;
//
//}
//
//
////rate_dusling(double kR, double T, double kkPiOverEta double * res, double * remainder)
//double rate_qgp_Dusling(double kR, double kkPiOverEta, double T) {
//
//	//No limits of the fit, thus no remainder
//
//
//}
//
////Kevin Dusling (arXiv:0903.1764)
//double dusling_factor(double k, double T) {
//
//	//E d^3 R/dk^3= k A(k)/(2 Pi)^3 1/2*log(3.7388 (kOverT)/g^2) (1+(1-nf(kOverT))*kOverT^2*eta/s*1/(2*T)*kkPOverEta)
//
//}



//Fermi-Dirac
double fermiDirac(double kOverT) {

	return 1.0/(exp(kOverT)+1.0);

}

//Bose-Einstien
double boseEin(double kOverT) {

	return 1.0/(exp(kOverT)-1.0);

}

//Rapp et al's photon rate from in-medium rho mesons, as parametrized in ...
double rate_hg_in_medium_rho_ideal_Rapp_fit(double kOverT, double T, double kkPiOver_e_P_k2) {

	//
	const double a=-31.21+353.61*T-1739.4*T*T+3105*T*T*T;
	const double b=-5.513-42.2*T+333*T*T-570*T*T*T;
	const double c=-6.153+57*T-134.61*T*T+8.31*T*T*T;

	const double k=kOverT*T;

	double res=exp(a*k+b+c/(k+0.2));

	return res;

}

//Rapp et al's photon rate from pi-rho-omega
double rate_hg_piRhoOmega_ideal_Rapp_fit(double kOverT, double T, double kkPiOver_e_P_k2) {

	//
	const double T2=T*T;
	const double T3=T2*T;
	const double k=kOverT*T;

//c--------------------------------------

	const double q0=k;
//	const double q0q=q0*q0;
 
// c********** pi+rho->gamma+omega
       const double a1=-35.8991+460.425*T-2592.04*T2+5342.32*T3;
       const double a2=-41.9725+601.952*T-3587.8 *T2+7604.97*T3;
       const double a3=0.740436-16.7159*T+133.526*T2-347.589*T3;
       const double a4= 2.00611-3.79343*T+29.3101*T2-72.8725*T3;
       const double a5=-8.33046+121.091*T-801.676*T2+1712.16*T3;
       const double a6=17.9029 -  388.5*T+2779.03*T2- 6448.4*T3;
       const double a7=-15.622 +340.651*T-2483.18*T2+5870.61*T3;
       const double FFpiro=exp(a1*q0+a2+a3*pow(q0,a4)+a5*pow((q0+a6),a7));
// c********** rho+omega->gamma+pi
      const double b1=-29.6866+331.769*T-1618.66*T2+2918.53*T3;
      const double b2=-15.3332+90.2225*T-300.185*T2+428.386*T3;
      const double b3=-7.35061+109.288*T-630.396*T2+1227.69*T3;
      const double b4=-10.6044+  109.1*T-500.718*T2+872.951*T3;
       const double FFomro=exp(b1*q0+b2+b3/(q0+0.2)+b4/pow((q0+0.2),2));
// c********** pi+omega->gamma+rho
       const double d1=-29.4663 +291.356*T-1301.27*T2+2102.12*T3;
       const double d2=-45.081  +688.929*T-4150.15*T2+8890.76*T3;
       const double d3=-0.260076+8.92875*T- 60.868*T2+ 136.57*T3;
       const double d4= 2.2663  -8.30596*T+49.3342*T2-90.8501*T3;
       const double d5= 10.2955 -317.077*T+2412.15*T2- 6020.9*T3;
       const double d6= 3.12251 -47.5277*T+ 222.61*T2-  241.9*T3;
       const double d7=-3.39045 +56.5927*T- 336.97*T2+622.756*T3;
       const double FFompi=exp(d1*q0+d2+d3*pow(q0,d4)+d5*pow((q0+d6),d7));

//      const double dR1d3q= fugaro*fugapi* (FFpiro + fugapi*FFompi + fugaro*FFomro);

	double res=FFpiro + FFompi + FFomro;

	return res;

}


//Rapp et al's photon rate from in-medium rho mesons, as parametrized in ...
double rate_hg_pion_brem_ideal_Rapp_fit(double kOverT, double T, double kkPiOver_e_P_k2) {

	//
	const double T2=T*T;
	const double T3=T2*T;
	const double a=-16.28+62.45*T-93.4*T2-7.5*T3;
	const double b=-35.54+414.8*T-2054*T2+3718.8*T3;
	const double c=0.7364-10.72*T+56.32*T2-103.5*T3;
	const double d=-2.51+58.152*T-318.24*T2+610.7*T3;

	const double k=kOverT*T;

	double res=exp(a+b*k+c*k*k+d/(k+0.2));

	return res;

}

//Fit of Zahed-Dusling's low T photon rate. 100 < T < 200 MeV, 0.2 < k < 5 GeV
double rate_hg_ideal_Zahed_Dusling_2pi_fit(double kOverT, double T, double kkPiOver_e_P_k2) {



	//
	const double lT=log(T);
	const double lk=log(kOverT*T);
	const double lT2=lT*lT;
	const double lT3=lT2*lT;
	const double lk2=lk*lk;
	const double lk3=lk2*lk;
	const double lk4=lk3*lk;

	double res=exp(3.9242043019455433
	               + 4.785547499773419*(1 + 0.41363881867911095*lk - 0.04811524892833928*lk2 + 0.3808205573095817*lk3 +0.11256358970445153*lk4)*lT 
	               - 1.8846748688191692*(1 - 0.06623529535737044*lk + 0.389375957931405*lk2 - 0.821714627504325*lk3 - 0.21057052677851257*lk4)*lT2 
		       + 0.06875606519946258*(1 + 5.228449112565715*lk + 1.9302736990409124*lk2 + 7.541355952332698*lk3 + 1.5150073214876982*lk4)*lT3);

	return res;

}

//Rapp et al's photon rate from in-medium rho mesons, as parametrized in ...
double rate_hg_in_medium_rho_ideal_Rapp_fit_PCE160(double kOverT, double T, double kkPiOver_e_P_k2) {

	//
	//const double r=(1+exp(-2*mu
	return 1.0/0.0;

}

//HG ideal rate - Simon Turbide's fit [???]
double rate_hg_ideal_Turbide_fit(double kOverT, double T, double kkPiOver_e_P_k2) {

	//
	double rate_hg_ideal_Turbide_fit_piRho_a1PiRho_piGamma(double E, double T);
	double rate_hg_ideal_Turbide_fit_piRho_omegatChan_piGamma(double E, double T);
	double rate_hg_ideal_Turbide_fit_piPi_a1PiRho_rhoGamma(double E, double T);
	double rate_hg_ideal_Turbide_fit_rho_a1PiRho_piPiGamma(double E, double T);
	double rate_hg_ideal_Turbide_fit_piKstar_KGamma(double E, double T);
	double rate_hg_ideal_Turbide_fit_piK_KstarGamma(double E, double T);
	double rate_hg_ideal_Turbide_fit_rhoK_KGamma(double E, double T);
	double rate_hg_ideal_Turbide_fit_KKstar_piGamma(double E, double T);

	const double E=kOverT*T;

	double res=0.0;

	res=rate_hg_ideal_Turbide_fit_piRho_a1PiRho_piGamma(E,T)+ \
	    rate_hg_ideal_Turbide_fit_piRho_omegatChan_piGamma(E,T)+ \
	    rate_hg_ideal_Turbide_fit_piPi_a1PiRho_rhoGamma(E,T)+ \
	    rate_hg_ideal_Turbide_fit_rho_a1PiRho_piPiGamma(E,T)+ \
	    rate_hg_ideal_Turbide_fit_piKstar_KGamma(E,T)+ \
	    rate_hg_ideal_Turbide_fit_piK_KstarGamma(E,T)+ \
	    rate_hg_ideal_Turbide_fit_rhoK_KGamma(E,T)+ \
	    rate_hg_ideal_Turbide_fit_KKstar_piGamma(E,T);

	return res;

}



//HG ideal rate - Simon Turbide's fit without pion bremstr. doublecontains
double rate_hg_ideal_Turbide_noPiPi_fit(double kOverT, double T, double kkPiOver_e_P_k2) {

	//
	double rate_hg_ideal_Turbide_fit_piRho_a1PiRho_piGamma(double E, double T);
	double rate_hg_ideal_Turbide_fit_piRho_omegatChan_piGamma(double E, double T);
	double rate_hg_ideal_Turbide_fit_piKstar_KGamma(double E, double T);
	double rate_hg_ideal_Turbide_fit_piK_KstarGamma(double E, double T);
	double rate_hg_ideal_Turbide_fit_rhoK_KGamma(double E, double T);
	double rate_hg_ideal_Turbide_fit_KKstar_piGamma(double E, double T);

	const double E=kOverT*T;

	double res=0.0;

	res=rate_hg_ideal_Turbide_fit_piRho_a1PiRho_piGamma(E,T)+ \
	    rate_hg_ideal_Turbide_fit_piRho_omegatChan_piGamma(E,T)+ \
	    rate_hg_ideal_Turbide_fit_piKstar_KGamma(E,T)+ \
	    rate_hg_ideal_Turbide_fit_piK_KstarGamma(E,T)+ \
	    rate_hg_ideal_Turbide_fit_rhoK_KGamma(E,T)+ \
	    rate_hg_ideal_Turbide_fit_KKstar_piGamma(E,T);

	return res;

}

//From Maxime's program
double FFpion(double E)
{
	double tbar = 34.5096*pow(E,0.737) - 67.556*pow(E,0.7584) + 32.858*pow(E,0.7806); 
	return pow(2.0/(2.0-tbar),2);
}

double FFw(double E)
{
	double tbar = -61.595*pow(E,0.9979) + 28.592*pow(E,1.1579) + 37.738*pow(E,0.9317) - 5.282*pow(E,1.3686);
	return pow(2.0/(2.0-tbar),2);
}

double FFk(double E)
{
	double tbar = -76.424*pow(E,0.6236) + 36.944*pow(E,0.6604) + 39.0448*pow(E,0.5873);
	return pow(2.0/(2.0-tbar),2);
}

// Hadronic rates whitout viscosity
double rate_hg_ideal_Turbide_fit_piRho_a1PiRho_piGamma(double E, double T) {
	return pow(FFpion(E),4)*pow(T,2.8)*exp(-(1.461*pow(T,2.3094)+0.727)/pow(2*T*E,0.86) + (0.566*pow(T,1.4094)-0.9957)*E/T);
}

double rate_hg_ideal_Turbide_fit_piRho_omegatChan_piGamma(double E, double T) {
	return pow(FFw(E),4)*pow(T,1.0)*exp((1.865*pow(T,1.02)-2.6)/pow(2*E*T,0.62) + (3.053*pow(T,1.8)-1.038)*E/T );

}

double rate_hg_ideal_Turbide_fit_piPi_a1PiRho_rhoGamma(double E, double T) {
	return pow(FFpion(E),4)*pow(T,-5)*exp(-(9.314*pow(T,-0.584)-5.328)*pow(2*T*E,0.088) + (0.3189*pow(T,0.721)-0.8998)*E/T);
}

double rate_hg_ideal_Turbide_fit_rho_a1PiRho_piPiGamma(double E, double T) {
	return pow(FFpion(E),4)*pow(T,-2)*exp(-(-35.459*pow(T,1.126)+18.827)/pow(2*E*T,-1.44*pow(T,0.142)+0.9996) -1.21*E/T);

}

double rate_hg_ideal_Turbide_fit_piKstar_KGamma(double E, double T) {
	return pow(FFpion(E),4)*pow(T,3.75)*exp(-0.35/pow(2*T*E,1.05)+(2.3894*pow(T,0.03435)-3.222)*E/T);
}

double rate_hg_ideal_Turbide_fit_piK_KstarGamma(double E, double T) {
	return pow(FFpion(E),4)*pow(T,-3)*exp(-(5.4018*pow(T,-0.6864) - 1.51)*pow(2*T*E,0.07) - 0.91*E/T);
}

double rate_hg_ideal_Turbide_fit_rhoK_KGamma(double E, double T) {
	return pow(FFk(E),4)*pow(T,3.5)*exp(-(0.9386*pow(T,1.551)+0.634)/pow(2*T*E,1.01) + (0.568*pow(T,0.5397)-1.164)*E/T);
}

double rate_hg_ideal_Turbide_fit_KKstar_piGamma(double E, double T) {
	return pow(FFk(E),4)*pow(T,3.7)*exp(-(6.096*pow(T,1.889)+1.0299)/pow(2*E*T,-1.613*pow(T,2.162)+0.975) -0.96*E/T);
}

double rate_hg_ideal_Turbide_fit_chem_pot_Boltz_Tch150(double kOverT, double T, double kkPiOver_e_P_k2) {

	double rate_hg_ideal_Turbide_fit_chem_pot_Boltz(double kOverT, double T, double kkPiOver_e_P_k2, enum chem_freezeout_temp Tch); 

	return rate_hg_ideal_Turbide_fit_chem_pot_Boltz(kOverT, T, kkPiOver_e_P_k2, Tch150);

}

double rate_hg_ideal_Turbide_fit_chem_pot_Boltz_Tch160(double kOverT, double T, double kkPiOver_e_P_k2) {

	double rate_hg_ideal_Turbide_fit_chem_pot_Boltz(double kOverT, double T, double kkPiOver_e_P_k2, enum chem_freezeout_temp Tch); 

	return rate_hg_ideal_Turbide_fit_chem_pot_Boltz(kOverT, T, kkPiOver_e_P_k2, Tch160);

}

double rate_hg_ideal_Turbide_fit_chem_pot_Boltz_Tch165(double kOverT, double T, double kkPiOver_e_P_k2) {

	double rate_hg_ideal_Turbide_fit_chem_pot_Boltz(double kOverT, double T, double kkPiOver_e_P_k2, enum chem_freezeout_temp Tch); 

	return rate_hg_ideal_Turbide_fit_chem_pot_Boltz(kOverT, T, kkPiOver_e_P_k2, Tch165);

}

double rate_hg_ideal_Turbide_fit_chem_pot_Boltz(double kOverT, double T, double kkPiOver_e_P_k2, enum chem_freezeout_temp Tch) {

	//
	double rate_hg_ideal_Turbide_fit_piRho_a1PiRho_piGamma(double E, double T);
	double rate_hg_ideal_Turbide_fit_piRho_omegatChan_piGamma(double E, double T);
	double rate_hg_ideal_Turbide_fit_piPi_a1PiRho_rhoGamma(double E, double T);
	double rate_hg_ideal_Turbide_fit_rho_a1PiRho_piPiGamma(double E, double T);
	double rate_hg_ideal_Turbide_fit_piKstar_KGamma(double E, double T);
	double rate_hg_ideal_Turbide_fit_piK_KstarGamma(double E, double T);
	double rate_hg_ideal_Turbide_fit_rhoK_KGamma(double E, double T);
	double rate_hg_ideal_Turbide_fit_KKstar_piGamma(double E, double T);
	double get_PH_chem_pot_pions(double T, enum chem_freezeout_temp Tch); 
	double get_PH_chem_pot_kaons(double T, enum chem_freezeout_temp Tch);

	const double E=kOverT*T;
	
	double res=0.0;
	double muPion=get_PH_chem_pot_pions(T, Tch);
	double muKaon=get_PH_chem_pot_kaons(T, Tch);
	double muRho=2*muPion;
	double muKstar=muPion+muKaon;
	double expMuPionOverT=exp(muPion/T);
	double expMuKaonOverT=exp(muKaon/T);
	double expMuRhoOverT=exp(muRho/T);
	double expMuKstarOverT=exp(muKstar/T);

	res=rate_hg_ideal_Turbide_fit_piRho_a1PiRho_piGamma(E,T)*expMuPionOverT*expMuRhoOverT+ \
	    rate_hg_ideal_Turbide_fit_piRho_omegatChan_piGamma(E,T)*expMuPionOverT*expMuRhoOverT+ \
	    rate_hg_ideal_Turbide_fit_piPi_a1PiRho_rhoGamma(E,T)*expMuPionOverT*expMuPionOverT+ \
	    rate_hg_ideal_Turbide_fit_rho_a1PiRho_piPiGamma(E,T)*expMuRhoOverT+ \
	    rate_hg_ideal_Turbide_fit_piKstar_KGamma(E,T)*expMuPionOverT*expMuKstarOverT+ \
	    rate_hg_ideal_Turbide_fit_piK_KstarGamma(E,T)*expMuPionOverT*expMuKaonOverT+ \
	    rate_hg_ideal_Turbide_fit_rhoK_KGamma(E,T)*expMuRhoOverT*expMuKaonOverT+ \
	    rate_hg_ideal_Turbide_fit_KKstar_piGamma(E,T)*expMuKaonOverT*expMuKstarOverT;

	return res;

}


//pi, rho, K, K^*, but rho and K^* are actually unstable
//Take T in GeV, return mu in GeV
double get_PH_chem_pot_pions(double T, enum chem_freezeout_temp Tch) {

	//
	double res;

	switch(Tch) {

		case Tch150:
			if (T>0.15) {
				res=0.0;
			}
			else {
				res=(0.1391802540030629 - 2.5394267515127775*T + 20.879628260096943*T*T - 67.56361472739776*T*T*T)/(1 - 11.157083798775501*T + 51.749504569561154*T*T);
			}
			break;
			
		case Tch160:
			if (T>0.16) {
				res=0.0;
			}
			else {
				res=(0.13989273701516192 - 2.3597467055955885*T + 19.326406907950684*T*T - 62.7506772821009*T*T*T)/(1 - 10.758531025652175*T + 50.30650344070788*T*T); 
			}
			break;

		case Tch165:
			if (T>0.165) {
				res=0.0;
			}
			else {
				res=(0.14147237963580053 - 2.3888305421239244*T + 23.09944275257463*T*T - 83.70949379912233*T*T*T)/(1 - 10.87028053710845*T + 64.05564313498687*T*T);
			}
			break;

	}
	
	return res;

}

double get_PH_chem_pot_kaons(double T, enum chem_freezeout_temp Tch) {

	//
	double res;

	switch(Tch) {

		case Tch150:
			if (T>0.15) {
				res=0.0;
			}
			else {
				res=(0.4907686964365895 - 5.8950237384113695*T + 48.53455196905417*T*T - 207.0194972599132*T*T*T)/(1 - 3.640490397093469*T + 55.61735320330468*T*T);
			}
			break;
			
		case Tch160:
			if (T>0.16) {
				res=0.0;
			}
			else {
				res=(0.49039646786244145 - 6.445326572927662*T + 63.72787484481105*T*T - 266.31358950869026*T*T*T)/(1 - 5.345400242205774*T + 79.774948475777*T*T);
			}
			break;

		case Tch165:
			if (T>0.165) {
				res=0.0;
			}
			else {
				res=(0.48692564193025784 - 10.07234004754924*T + 115.61157740249656*T*T - 439.19569740354456*T*T*T)/(1 - 13.489702407741524*T + 143.5787763163101*T*T);
			}
			break;

	}
	

	return res;

}


double T2_normalisation(double kOverT, double T, double kk) {
	const double res=T;
	return res*res;
}


double minus_sign_normalisation(double kOverT, double T, double kk) {
	return -1;
}

double rate_thermal_ideal(double kOverT, double T, double kkPiOver_e_P_k2) {

	//
	double qgpFrac=QGP_fraction(T);
	double hgFrac=1.-qgpFrac;

	//
	double res=0.0;
	if (hgFrac > 0.0) {
		res+=hgFrac*rate_hg_ideal_Turbide_noPiPi_fit(kOverT,T,kkPiOver_e_P_k2);
		res+=hgFrac*rate_hg_in_medium_rho_ideal_Rapp_fit(kOverT,T,kkPiOver_e_P_k2);
		res+=hgFrac*rate_hg_pion_brem_ideal_Rapp_fit(kOverT,T,kkPiOver_e_P_k2);
		res+=hgFrac*rate_hg_piRhoOmega_ideal_Rapp_fit(kOverT,T,kkPiOver_e_P_k2);
	}


	if (qgpFrac > 0.0) {
		res+=qgpFrac*rate_qgp_ideal_LO_AMYfit(kOverT,T,kkPiOver_e_P_k2);
	}

	return res;

}

double rate_hadronic_ideal(double kOverT, double T, double kkPiOver_e_P_k2) {

	//
	double qgpFrac=QGP_fraction(T);
	double hgFrac=1.-qgpFrac;

	//
	double res=0.0;
	if (hgFrac > 0.0) {
		res+=hgFrac*rate_hg_ideal_Turbide_noPiPi_fit(kOverT,T,kkPiOver_e_P_k2);
		res+=hgFrac*rate_hg_in_medium_rho_ideal_Rapp_fit(kOverT,T,kkPiOver_e_P_k2);
		res+=hgFrac*rate_hg_pion_brem_ideal_Rapp_fit(kOverT,T,kkPiOver_e_P_k2);
		res+=hgFrac*rate_hg_piRhoOmega_ideal_Rapp_fit(kOverT,T,kkPiOver_e_P_k2);
	}

	return res;

}
