#include "TLorentzVector.h"

class ttbar_solver{
public:
    void setLightJetA(const TLorentzVector& jet){ljeta=jet;}
    void setLightJetB(TLorentzVector jet){ljetb=jet;}

    //for full hadronic events
    void setLightJetC(TLorentzVector jet){
    ljetc=jet;
    FH=true;}
    void setLightJetD(TLorentzVector jet){ljetd=jet;}

    void setBJetA(TLorentzVector jet){bjeta=jet;}
    void setBJetB(TLorentzVector jet){bjetb=jet;}

    void setLepton(TLorentzVector lep){
    lepton=lep;
    FH=false;}
    void setNeutrino(TLorentzVector neu){neutrino=neu;}

    double getChi2()const;

private:
    TLorentzVector getNeutrinoVector(TLorentzVector lepton, TLorentzVector met);

    TLorentzVector ljeta, ljetb,ljetc,ljetd, bjeta, bjetb, lepton, neutrino;

    bool FH=false;

};

