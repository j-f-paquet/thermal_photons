#include "photon.h"

/***** Rates *****/

//Template for rate E d^3 Gamma/d k^3
//last argument can be dummy in ideal case
double rate_template(double kOverT, double T, double kkPiOver_e_P_k2) {

	//Compute the whole rate

	//Return it

}





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

	//Forward declaration
	double C_hard(double kOverT);

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

