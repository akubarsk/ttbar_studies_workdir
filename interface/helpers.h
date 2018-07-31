#include "TLorentzVector.h"
#include "TMath.h"

double deltaPhi(const double & a, const double & b);

template<class T>
TLorentzVector makeTLorentzVector(T* c){
    TLorentzVector out;
    out.SetPtEtaPhiM(c->PT,c->Eta,c->Phi,c->Mass);
    return out;
}

template<class T>
double deltaR(T* a, T* b){
double dphi= deltaPhi(a->Phi, b->Phi);
double deta=a->Eta - b->Eta;
return TMath::Sqrt(dphi*dphi+deta*deta);
}

template<class T>
double deltaR(T* a, double lep_phi, double lep_eta){
double dphi= deltaPhi(a->Phi, lep_phi);
double deta=a->Eta - lep_eta;
return TMath::Sqrt(dphi*dphi+deta*deta);
}


template<class T>
double JetDeltaR(T* a){
double deta = a->DeltaEta;
double dphi = a->DeltaPhi;
return TMath::Sqrt(dphi*dphi+deta*deta);
}
