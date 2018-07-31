#include "../interface/pz_calculator.h"


//Pz calculator of neutrino candidate from missing transverse energy
double pz_calculator:: getPz()const{

	double M_W=80.4;

	double emu = lepton.E();
	double pxmu = lepton.Px();
	double pymu = lepton.Py();
 	double pzmu = lepton.Pz();
  	double pxnu = met.Px();
  	double pynu = met.Py();
	double pznu = 0.;

	double a = M_W*M_W - M_lepton*M_lepton + 2.0*(pxmu*pxnu + pymu*pynu);
  	double A = 4.0*(emu*emu - pzmu*pzmu);
  	double B = -4.0*a*pzmu;
	double C = 4.0*emu*emu*(pxnu*pxnu + pynu*pynu) - a*a;
			
	double tmproot = B*B - 4.0*A*C;

	if (tmproot<0) {
		pznu = - B/(2*A); //take real part of complex roots
	}
	else{
    		double tmpsol1 = (-B + TMath::Sqrt(tmproot))/(2.0*A);
		double tmpsol2 = (-B - TMath::Sqrt(tmproot))/(2.0*A);

		if (fabs(tmpsol2-pzmu) < fabs(tmpsol1-pzmu)) { pznu = tmpsol2;}
		else pznu = tmpsol1;
	}
	return pznu;
}
