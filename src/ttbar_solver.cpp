#include "../interface/ttbar_solver.h"

double ttbar_solver:: getChi2()const{
//calculate chi2

//masses
double M_W=80.4;
double M_t=173;

//resolutions
double sigmabjet=5.0;
double sigmaljet=5.0;
double sigmalepton=1.0;
double sigmaneutrino=10.0;

double Chi2;
if(!FH){
	TLorentzVector w1=neutrino+lepton;
	TLorentzVector top1=w1+bjeta;
	TLorentzVector w2=ljeta+ljetb;
	TLorentzVector top2=w2+bjetb;

	Chi2=(M_W-w1.M())*(M_W-w1.M())/((sigmalepton+sigmaneutrino)*(sigmalepton+sigmaneutrino)) + (M_W-w2.M())*(M_W-w2.M())/((2*sigmaljet)*(2*sigmaljet)) + (M_t-top1.M())*(M_t-top1.M())/((sigmalepton+sigmaneutrino+sigmabjet)*(sigmalepton+sigmaneutrino+sigmabjet)) + (M_t-top2.M())*(M_t-top2.M())/((2*sigmaljet+sigmabjet)*(2*sigmaljet+sigmabjet));

}else{

	TLorentzVector w1=ljeta+ljetb;
        TLorentzVector top1=w1+bjeta;
        TLorentzVector w2=ljetc+ljetd;
        TLorentzVector top2=w2+bjetb;

	Chi2=(M_W-w1.M())*(M_W-w1.M())/((2*sigmaljet)*(2*sigmaljet)) + (M_W-w2.M())*(M_W-w2.M())/((2*sigmaljet)*(2*sigmaljet)) + (M_t-top1.M())*(M_t-top1.M())/((2*sigmaljet+sigmabjet)*(2*sigmaljet+sigmabjet)) + (M_t-top2.M())*(M_t-top2.M())/((2*sigmaljet+sigmabjet)*(2*sigmaljet+sigmabjet));

}
return Chi2;
}

