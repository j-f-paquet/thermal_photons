#include "photon.h"

/***** Rates *****/

/*
Rates:
1: double rate_qgp_ideal_born_AMYfit(double kOverT, double T, double kkPiOver_e_P_k2);
2: double rate_qgp_ideal_born_KLS(double kOverT, double T, double kkPiOver_e_P_k2);
3: double rate_qgp_ideal_born_JF_sqrtg(double kOverT, double T, double kkPiOver_e_P_k2);
4: double rate_qgp_viscous_only_born_JF_sqrtg(double kOverT, double T, double kkPiOver_e_P_k2);
5: rate_hg_ideal_Turbide_fit
6: rate_qgp_ideal_LO_AMYfit
*/


double get_photon_rate(int selector, double (**local_rate)(double, double, double)) {

	double rate_qgp_ideal_born_AMYfit(double, double, double);
	double rate_qgp_ideal_born_KLS(double, double, double);
	double rate_qgp_ideal_born_JF_sqrtg(double kOverT, double T, double kkPiOver_e_P_k2);
	double rate_qgp_viscous_only_born_JF_sqrtg(double kOverT, double T, double kkPiOver_e_P_k2);
	double rate_hg_ideal_Turbide_fit(double kOverT, double T, double kkPiOver_e_P_k2);
	double rate_qgp_ideal_LO_AMYfit(double, double, double);

	//double (*local_rate)(double, double, double) = CONST_rateList[iRate].c_str();
	switch(selector) {
		case 1:
			//double (*local_rate)(double, double, double) = rate_qgp_ideal_born_AMYfit;
			*local_rate = rate_qgp_ideal_born_AMYfit;
//			*rate_name= "rate_qgp_ideal_born_AMYfit";
			break;
		case 2: 
			*local_rate = rate_qgp_ideal_born_KLS;
//			*rate_name="rate_qgp_ideal_born_KLS";
			break;
		case 3: 
			*local_rate = rate_qgp_ideal_born_JF_sqrtg;
//			*rate_name="rate_qgp_ideal_born_JF_sqrtg";
			break;
		case 4:
			*local_rate = rate_qgp_viscous_only_born_JF_sqrtg;
//			*rate_name="rate_qgp_viscous_only_born_JF_sqrtg";
			break;
		case 5:
			*local_rate = rate_hg_ideal_Turbide_fit;
//			*rate_name="rate_hg_ideal_Turbide_fit";
			break;
		case 6:
			*local_rate = rate_qgp_ideal_LO_AMYfit;
//			*rate_name="rate_qgp_ideal_LO_AMYfit";
			break;
	}



}

double kOverT_from_index(int i, int size) {
	return 80*i*i*1.0/(size*size);
}

int index_from_kOverT(double kOverT, int size) {
	return floor(sqrt(kOverT*size*size*1.0/80.0));
}

double temp_from_index(int i, int size) {
	return 1.0*i*i*1.0/(size*size);
}

int index_from_temp(double temp, int size) {
	return floor(sqrt(temp*size*size*1.0/1.0));
}

double get_photon_rate_accel(double kOverT, double T, double kk, int rate_no) {

	double res;

	const int size_x=accel_table_sample_x[CONST_rates_to_use[rate_no]-1];
	const int size_y=accel_table_sample_y[CONST_rates_to_use[rate_no]-1];

	int a1=index_from_kOverT(kOverT,size_x);
	int b1=index_from_temp(T,size_y);		

	//Must decide what to do when the requested value is outside the table's range
	if ((a1+1>=size_x)||(b1+1>=size_y)||(a1<0)||(b1<0)) {
		res=0.0;
	}
	else {

		double fx1y1=CONST_rate_tables.tabulated_rates[rate_no][b1][a1];
		double fx2y1=CONST_rate_tables.tabulated_rates[rate_no][b1][a1+1];
		double fx1y2=CONST_rate_tables.tabulated_rates[rate_no][b1+1][a1];
		double fx2y2=CONST_rate_tables.tabulated_rates[rate_no][b1+1][a1+1];

		double x1=kOverT_from_index(a1,size_x);
		double x2=kOverT_from_index(a1+1,size_x);
		double y1=temp_from_index(b1,size_y);
		double y2=temp_from_index(b1+1,size_y);

		res=(fx1y1*(x2-kOverT)*(y2-T)+fx2y1*(kOverT-x1)*(y2-T)+fx1y2*(x2-kOverT)*(T-y1)+fx2y2*(kOverT-x1)*(T-y1))/((x2-x1)*(y2-y1));

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

	double res=QGP_fraction(T);

	if (res > 0) res*=CONST_GeV2_to_GeVm2_fmm4*kOverT*prefA(kOverT,T)/CONST_twoPiCubed*(log(1/CONST_mInfOverT)+C_hard(kOverT));

	return res;

}

//QGP ideal rate - LO - AMY fit [arXiv:hep-ph/0111107, section I]
double rate_qgp_ideal_LO_AMYfit(double kOverT, double T, double kkPiOver_e_P_k2) {

	//Forward declaration
	double C_hard(double kOverT);
	double C_LPM(double kOverT);

	double res=QGP_fraction(T);

	if (res > 0) res*=CONST_GeV2_to_GeVm2_fmm4*kOverT*prefA(kOverT,T)/CONST_twoPiCubed*(log(1/CONST_mInfOverT)+C_hard(kOverT)+C_LPM(kOverT));

	return res;

}


//QGP ideal rate - KLS high k/T, low g formula [Kapusta et al, PRD44, 9 (1991), eq.41]
double rate_qgp_ideal_born_KLS(double kOverT, double T, double kkPiOver_e_P_k2) {

	const double qCharge2[]={4.0/9.0,5.0/9.0,2.0/3.0,10.0/9.0,11.0/9.0,5.0/3.0};

	double res=QGP_fraction(T);

	if (res > 0) res*=CONST_GeV2_to_GeVm2_fmm4*qCharge2[CONST_Nf-1]*CONST_alphaEM*CONST_alphaS/(2.0*M_PI*M_PI)*T*T*exp(-kOverT)*log(2.912*kOverT/(CONST_gs*CONST_gs));

	return res;

}

//QGP ideal rate - JF fit - q^*=sqrt(g)
double rate_qgp_ideal_born_JF_sqrtg(double kOverT, double T, double kkPiOver_e_P_k2) {

	double res=QGP_fraction(T);

	if (res > 0) res*=CONST_GeV2_to_GeVm2_fmm4*kOverT*prefA(kOverT,T)/CONST_twoPiCubed*(0.8452052719374467 + 0.06345436545672481*kOverT + 0.20266453593373313*kOverT*kOverT + 0.007103855524696941*kOverT*kOverT*kOverT)/(1 + 0.3137709585719375*kOverT + 0.12623968017081683*kOverT*kOverT + 0.0021744062978126125*kOverT*kOverT*kOverT);

	return res;

}

//viscous correction to rate: A_\alpha\beta K^\alpha K^\beta/k^2 k A(k)/(2 pi)^3 * viscous_correction_born_JF_sqrtg()
double rate_qgp_viscous_only_born_JF_sqrtg(double kOverT, double T, double kkPiOver_e_P_k2) {

	double res = QGP_fraction(T);
	//
	if (res > 0) res*=CONST_GeV2_to_GeVm2_fmm4*kOverT*prefA(kOverT,T)/CONST_twoPiCubed*kkPiOver_e_P_k2*exp(-0.5041041126181884 + (-0.5335015716121183 + 1.9967643068761307*kOverT - 0.5616138941792664*kOverT*kOverT - 0.0009120108228910325*kOverT*kOverT*kOverT)/(1 - 2.607918425474197*kOverT - 0.8369709712322181*kOverT*kOverT))*pow(kOverT,2.1309931380115588);
	
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


//HG ideal rate - Simon Turbide's fit [???]
double rate_hg_ideal_Turbide_fit(double kOverT, double T, double kkPiOver_e_P_k2) {

	//
	double HadronicPhase(double E, double T, int process);

	double res=(1-QGP_fraction(T));
	double tmp=0.0;

	if (res > 0) {
		for(int i=1; i<=8;i++) tmp+=HadronicPhase(kOverT*T, T, i);
		res*=tmp;
	}

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

// Hadronic rate whitout viscosity

double HadronicPhase(double E, double T, int process)
{

	double A2=0, A3=0, A4=0, A5=0, A6=0, A7=0, A8=0, A9=0;
	switch (process) {
		case 1:
			A2 = pow(FFpion(E),4)*pow(T,2.8)*exp(-(1.461*pow(T,2.3094)+0.727)/pow(2*T*E,0.86) + (0.566*pow(T,1.4094)-0.9957)*E/T);
			break;
		case 2: 
			A3 = pow(FFpion(E),4)*pow(T,-5)*exp(-(9.314*pow(T,-0.584)-5.328)*pow(2*T*E,0.088) + (0.3189*pow(T,0.721)-0.8998)*E/T);
			break;
		case 3:
			A4 = pow(FFpion(E),4)*pow(T,-2)*exp(-(-35.459*pow(T,1.126)+18.827)/pow(2*E*T,-1.44*pow(T,0.142)+0.9996) -1.21*E/T);
			break;
		case 4:
			A5 = pow(FFpion(E),4)*pow(T,3.75)*exp(-0.35/pow(2*T*E,1.05)+(2.3894*pow(T,0.03435)-3.222)*E/T);
			break;
		case 5:
			A6 = pow(FFpion(E),4)*pow(T,-3)*exp(-(5.4018*pow(T,-0.6864) - 1.51)*pow(2*T*E,0.07) - 0.91*E/T);
			break;
		case 6:
			A7 = pow(FFk(E),4)*pow(T,3.5)*exp(-(0.9386*pow(T,1.551)+0.634)/pow(2*T*E,1.01) + (0.568*pow(T,0.5397)-1.164)*E/T);
			break;
		case 7:
			A8 = pow(FFk(E),4)*pow(T,3.7)*exp(-(6.096*pow(T,1.889)+1.0299)/pow(2*E*T,-1.613*pow(T,2.162)+0.975) -0.96*E/T);
			break;
		case 8: 
			A9 = pow(FFw(E),4)*pow(T,1.0)*exp((1.865*pow(T,1.02)-2.6)/pow(2*E*T,0.62) + (3.053*pow(T,1.8)-1.038)*E/T );
			break;
	}
	return (A2 + A3 + A4 + A5 + A6 + A7 + A8 + A9);
}



