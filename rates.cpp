#include "photon.h"

/***** Rates *****/

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

//QGP ideal rate - AMY fit [arXiv:hep-ph/0111107, section I]
double rate_qgp_ideal_born_AMYfit(double kOverT, double T, double kkPiOver_e_P_k2) {

	//Forward declaration
	double C_hard(double kOverT);

	return kOverT*prefA(kOverT,T)/CONST_twoPiCubed*(log(1/CONST_mInfOverT)+C_hard(kOverT));

}


//QGP ideal rate - KLS high k/T, low g formula [Kapusta et al, PRD44, 9 (1991), eq.41]
double rate_qgp_ideal_born_KLS(double kOverT, double T, double kkPiOver_e_P_k2) {

	const double qCharge2[]={4.0/9.0,5.0/9.0,2.0/3.0,10.0/9.0,11.0/9.0,5.0/3.0};

	return QGP_fraction(T)*qCharge2[CONST_Nf-1]*CONST_alphaEM*CONST_alphaS/(2.0*M_PI*M_PI)*T*T*exp(-kOverT)*log(2.912*kOverT/(CONST_gs*CONST_gs));

}

//QGP ideal rate - JF fit - q^*=sqrt(g)
double rate_qgp_ideal_born_JF_sqrtg(double kOverT, double T, double kkPiOver_e_P_k2) {

	return kOverT*prefA(kOverT,T)/CONST_twoPiCubed*(0.8452052719374467 + 0.06345436545672481*kOverT + 0.20266453593373313*kOverT*kOverT + 0.007103855524696941*kOverT*kOverT*kOverT)/(1 + 0.3137709585719375*kOverT + 0.12623968017081683*kOverT*kOverT + 0.0021744062978126125*kOverT*kOverT*kOverT);

}

//viscous correction to rate: A_\alpha\beta K^\alpha K^\beta/k^2 k A(k)/(2 pi)^3 * viscous_correction_born_JF_sqrtg()
double rate_qgp_viscous_only_born_JF_sqrtg(double kOverT, double T, double kkPiOver_e_P_k2) {

	//
	return  kOverT*prefA(kOverT,T)/CONST_twoPiCubed*kkPiOver_e_P_k2*exp(-0.5041041126181884 + (-0.5335015716121183 + 1.9967643068761307*kOverT - 0.5616138941792664*kOverT*kOverT - 0.0009120108228910325*kOverT*kOverT*kOverT)/(1 - 2.607918425474197*kOverT - 0.8369709712322181*kOverT*kOverT))*pow(kOverT,2.1309931380115588);

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

