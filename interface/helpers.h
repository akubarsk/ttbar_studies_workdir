#include "TLorentzVector.h"
#include "TMath.h"

double deltaPhi(const double & a, const double & b);

template<class T>
TLorentzVector makeTLorentzVector(T* c){
    TLorentzVector out;
    out.SetPtEtaPhiM(c->PT,c->Eta,c->Phi,c->Mass);
    return out;
}

template<class K, class U>
double deltaR(K* a, U* b){
double dphi= deltaPhi(a.Phi(), b.Phi());
double deta=a.PseudoRapidity() - b.PseudoRapidity();
return TMath::Sqrt(dphi*dphi+deta*deta);
}

template<class T>
double JetDeltaR(T* a){
double deta = a->DeltaEta;
double dphi = a->DeltaPhi;
return TMath::Sqrt(dphi*dphi+deta*deta);
}
