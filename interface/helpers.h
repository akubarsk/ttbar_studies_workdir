#include "TLorentzVector.h"
#include "TMath.h"

double deltaPhi(const double & a, const double & b);
double dR(const double & a, const double & b, const double & c, const double & d);

//makes Lorentz vector from jets
template<class T>
TLorentzVector makeTLorentzVector(T* c){
    TLorentzVector out;
    out.SetPtEtaPhiM(c->PT,c->Eta,c->Phi,c->Mass);
    return out;
}

//calculeting angular distance between two objects
template<class T>
double deltaR(T* a, T* b){
double dphi= deltaPhi(a->Phi, b->Phi);
double deta=a->Eta - b->Eta;
return TMath::Sqrt(dphi*dphi+deta*deta);
}
//calculating angular distance if we cannot use leptons class
template<class T>
double deltaR(T* a, double lep_phi, double lep_eta){
double dphi= deltaPhi(a->Phi, lep_phi);
double deta=a->Eta - lep_eta;
return TMath::Sqrt(dphi*dphi+deta*deta);
}

//calculating angular width of jet
template<class T>
double JetDeltaR(T* a){
double deta = a->DeltaEta;
double dphi = a->DeltaPhi;
return TMath::Sqrt(dphi*dphi+deta*deta);
}

//checking minimal conditions for B-jets
template<class T>
bool eligibleBJet(T* a){
bool eligible=true;
if(fabs(a->Eta)>2.4){eligible=false;}
if(a->PT < 50){eligible=false;}

return eligible;
}

//checking minimal conditions for next-to leading jet
template<class T>
bool eligibleLeadingJet(T* a){
bool eligible=true;
if(fabs(a->Eta)>2.5){eligible=false;}
if(a->PT < 150){eligible=false;}

return eligible;
}

//checking conditions for light jets
template<class T>
bool eligibleLightJet(T* a){
bool eligible=true;
if(a->BTag){eligible=false;}
if(fabs(a->Eta)>2.5){eligible=false;}

return eligible;
}

//checking conditions for leptons(electrons and muons)
template<class T>
bool eligibleLepton(T* a){
bool eligible=true;
if(a->PT <45){eligible=false;}
if(fabs(a->Eta)>2.1){eligible=false;}

return eligible;
}
