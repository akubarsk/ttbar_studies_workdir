#include "TLorentzVector.h"
#include "TMath.h"

class pz_calculator{
public:
    
    void setLepton(const TLorentzVector lep){lepton=lep;}
    void setMET(TLorentzVector metP4){met=metP4;}
    void setLeptonMass(double Mlep){M_lepton=Mlep;}

    double getPz()const;

private:
    TLorentzVector getNeutrinoVector(TLorentzVector lepton, TLorentzVector met);

    TLorentzVector lepton, met;
    double M_lepton;

};
