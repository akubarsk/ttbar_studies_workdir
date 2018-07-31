/*
 * testanalyser.cpp
 *
 *  Created on: 24 Aug 2016
 *      Author: jkiesele
 */

#include "interface/testanalyser.h"
#include "../interface/helpers.h"
#include "../interface/ttbar_solver.h"
#include "../interface/pz_calculator.h"

void testanalyser::analyze(size_t childid /* this info can be used for printouts */){

	/*
	 * This skeleton analyser runs directly on the Delphes output.
	 * It can be used to create histograms directly or a skim.
	 * If a skim is created, a new input configuration will be written automatically
	 * and stored in the output directory together with the ntuples.
	 * The skim can contain delphes objects again or can be flat. This is up
	 * to the user.
	 * Examples for both are given here.
	 *
	 * The same skeleton can be used to read the skim. Please refer to the comments
	 * marked with "==SKIM=="
	 *
	 * These parts are commented, since the code is supposed to work as an example without
	 * modifications on Delphes output directly.
	 */



	/*
	 * Define the branches that are to be considered for the analysis
	 * This branch handler (notice the "d")
	 * is used to run directly in Delphes output.
	 * For skimmed ntuples, see below
	 */
	d_ana::dBranchHandler<Electron> elecs(tree(),"Electron");
	/*
	 * Other branches might be the following
	 * (for a full list, please inspect the Delphes sample root file with root)
	 * For the Delphes class description, see $DELPHES_PATH/classes/DelphesClasses.h
	 */
	//
	d_ana::dBranchHandler<HepMCEvent>  event(tree(),"Event");
	d_ana::dBranchHandler<Jet>         genjet(tree(),"GenJet");
	d_ana::dBranchHandler<Jet>         jet(tree(),"Jet");
	d_ana::dBranchHandler<Muon>        muontight(tree(),"Muon");
	d_ana::dBranchHandler<MissingET>   met(tree(),"MissingET");


	/* ==SKIM==
	 *
	 * If a skim of the Delphes outout was created in a way indicated
	 * further below, use the tBranchHandler (please notive the "t")
	 * to access vectors of objects...
	 *
	 */
	// d_ana::tBranchHandler<std::vector<Electron> > electrons(tree(),"Electrons");

	/*==SKIM==
	 *
	 * Or an object directly
	 *
	 */
	//d_ana::tBranchHandler<MissingET> met(tree(),"MET");



	/*
	 * Always use this function to add a new histogram (can also be 2D)!
	 * Histograms created this way are automatically added to the output file
	 */
	//Playing with values
	TH1* histoEr=addPlot(new TH1D("Errors","number of double counting",10,0,10),"n","N_{event}");

	//sorted invariant masses
	TH1* histoW=addPlot(new TH1D("W boson","jets from W decay",100,10,150),"M [GeV]","N_{jets}");
	TH1* histoW3C=addPlot(new TH1D("Wplusjet","m_{inv} of 3 jets(W sorted)",100,0,500),"M [GeV]","N_{jets}");
	TH1* histoT=addPlot(new TH1D("tquark","m_{inv} of jets from t-quark",100,10,500),"M [GeV]","N");
	/*TH1* histoSaT=addPlot(new TH1D("tquark_Sa","m_{inv} of 3 jets, where 1 is W",100,0,500),"M [GeV]","N");*/

	//boosted top quark
	TH1* histoSL_Boosted=addPlot(new TH1D("tquark_SL_boosted","t-quark m_{inv} mass in semileptonic decay",100,0,300),"M [GeV]","N");
	TH1* histoFromLep=addPlot(new TH1D("tquark_FromLep","t-quark m_{inv} mass from blv in semileptonic decay",100,0,500),"M [GeV]","N");	
	TH1* histoFH=addPlot(new TH1D("tquark_FH","t-quark m_{inv} mass in full hadronic decay",100,0,500),"M [GeV]","N");	

	//Chi minimazing method
	TH1* histoSL=addPlot(new TH1D("tquark_SL","t-quark m_{inv} mass in semileptonic decay",50,100,400),"M [GeV]","N");
	
	/*
	 * If (optionally) a skim or a flat ntuple is to be created, please use the following function to initialize
	 * the tree.
	 * The output files will be written automatically, and a config file will be created.
	 */
	TTree* myskim=addTree();
	/*
	 * Add a simple branch to the skim
	 */
	Double_t elecPt=0;
	myskim->Branch("elecPt", &elecPt);
	/*
	 * Or store a vector of objects (also possible to store only one object)
	 */
	std::vector<Electron> skimmedelecs;
	myskim->Branch("Electrons",&skimmedelecs);



	size_t nevents=tree()->entries();
	if(isTestMode())
		nevents/=100;
	for(size_t eventno=0;eventno<nevents;eventno++){
		/*
		 * The following two lines report the status and set the event link
		 * Do not remove!
		 */
		reportStatus(eventno,nevents);
		tree()->setEntry(eventno);

		//Basic "constants"
		double M_W  = 80.4;
  		double M_mu =  0.10566;
  		double M_e = 0.511e-3;
		double M_t = 173.1;

		/*
		 * Begin the event-by-event analysis
		 */		

		
		//Event analysis
		//Limits for event:
		//|eta|< 2.5 for jets
		//At least 2 b-jets
		//At least 2 light jets(4 for FH)

		//Counting eligible jets
		size_t bjets=0;
		size_t ljets=0;
		for(size_t i=0;i<jet.size();i++){
			if(fabs(jet.at(i)->Eta)<2.5){
				if(jet.at(i)->BTag){
					bjets++;
				}else{ljets++;}
			}
		}
		
		if(bjets<2){continue;}
		if(ljets<2){continue;}

		//Semieptonic
		//Limits for events:
		//|Eta|<2.1 for leptons
		//P_T>30GeV for leptons
		//Only one lepton
		
		//Class for calculating neutrino Pz
		pz_calculator neutrPZ;

		//Class for finding Chi^2
		ttbar_solver solve;

		bool interesting=true;
		TLorentzVector lepton;
		double M_lepton=M_e; //setting lepton mass to elektron, changed further if that not the case
		double lepton_eta;
		double lepton_phi;
		bool second=false;
		
		//Searching events with only one lepton
		
		//Searching eligible electron
		for(size_t i=0;i<elecs.size();i++){
			if((elecs.at(i)->PT)>30 && fabs(elecs.at(i)->Eta)<2.1){
				if(second){
					interesting=false;
				}else{
					second=true;
					lepton=elecs.at(i)->P4();
					lepton_eta=elecs.at(i)->Eta;
					lepton_phi=elecs.at(i)->Phi;
				}
			}
		 }
		//Searching eligible muon
		for(size_t i=0;i<muontight.size();i++){
			if((muontight.at(i)->PT)>30 && fabs(muontight.at(i)->Eta)<2.1){
                        	if(second){
                                	interesting=false;
                                }else{
                                       	second=true;
                                        lepton=muontight.at(i)->P4();
					M_lepton=M_mu;
					lepton_eta=muontight.at(i)->Eta;
					lepton_phi=muontight.at(i)->Phi;
				}
			}
		}

		if(interesting){			
			//missing transversal energy
			TLorentzVector metP4=(met.at(0)->P4());
			
			//p_z calculating
			neutrPZ.setLepton(lepton);
			neutrPZ.setMET(metP4);
			neutrPZ.setLeptonMass(M_lepton);			

			//neutrino (candidate) Lorentz vector
			TLorentzVector neutrino(metP4.Px(),metP4.Py(),neutrPZ.getPz(),TMath::Sqrt(TMath::Power(metP4.Pt(),2)+TMath::Power(neutrPZ.getPz(),2)));

			//fitting leptonic decay(boosted)
			bool check = true; //testing is there multiple bjets in range(for myself)

			for(size_t i=0;i<jet.size();i++){
				if(jet.at(i)->BTag && (jet.at(i)->Eta)<2.5){
					TLorentzVector v1 = makeTLorentzVector(jet.at(i));
					//double DeltaR = deltaR(v1,lepton);
					double DEta=jet.at(i)->Eta - lepton.PseudoRapidity();
					double DPhi=jet.at(i)->Phi - lepton.Phi();
					double DeltaR=TMath::Sqrt(TMath::Power(DEta,2)+TMath::Power(DPhi,2));

					if(DeltaR<1.2){ // mass_t from leptonic decay
						v1=v1+lepton+neutrino;	
						//if lepton in jet, then substracking its vector
						double JetR=JetDeltaR(jet.at(i));
						if(DeltaR<JetR){v1=v1-lepton;}
						
						//histogram from blv decay
						histoFromLep->Fill(v1.M());
						
						//checking how much I double read lepton decay with this algorythm
						if(!check){histoEr->Fill(1);}
						check = false;

					}else{//mass_t from hadronic decay
							for(size_t j=0;j<jet.size();j++){
								if(!(jet.at(j)->BTag) && (jet.at(j)->Eta)<2.5){
									DEta=jet.at(i)->Eta - jet.at(j)->Eta;
									DPhi=jet.at(i)->Phi - jet.at(j)->Phi;
									DeltaR=TMath::Sqrt(TMath::Power(DEta,2)+TMath::Power(DPhi,2));
									if(DeltaR<1.2){
										for(size_t k=j;k<jet.size();k++){
											if(!(jet.at(k)->BTag) && (jet.at(k)->Eta)<2.5 && k!=j){
												int leading; //leading jet
												if(jet.at(i)->PT > jet.at(j)->PT){leading=i;}
												else{leading=j;}
												if(jet.at(leading)->PT < jet.at(k)->PT){leading=k;}
															
												if(leading!=k){
													DEta=jet.at(leading)->Eta - jet.at(k)->Eta;
													DPhi=jet.at(leading)->Phi - jet.at(k)->Phi;
													DeltaR=TMath::Sqrt(TMath::Power(DEta,2)+TMath::Power(DPhi,2));
												}else{
													DEta=jet.at(i)->Eta - jet.at(k)->Eta;
													DPhi=jet.at(i)->Phi - jet.at(k)->Phi;
													DeltaR=TMath::Sqrt(TMath::Power(DEta,2)+TMath::Power(DPhi,2));

													DEta=jet.at(j)->Eta - jet.at(k)->Eta;
													DPhi=jet.at(j)->Phi - jet.at(k)->Phi;
													double DeltaR_l=TMath::Sqrt(TMath::Power(DEta,2)+TMath::Power(DPhi,2));
													if(DeltaR<DeltaR_l){DeltaR=DeltaR_l;}//Taking into account highest value
													
												}
												if(DeltaR<1.2){
													TLorentzVector v1;
													TLorentzVector v2;
													TLorentzVector v3;
													v1.SetPtEtaPhiM(jet.at(i)->PT,jet.at(i)->Eta,jet.at(i)->Phi,jet.at(i)->Mass);
													v2.SetPtEtaPhiM(jet.at(j)->PT,jet.at(j)->Eta,jet.at(j)->Phi,jet.at(j)->Mass);
													v3.SetPtEtaPhiM(jet.at(k)->PT,jet.at(k)->Eta,jet.at(k)->Phi,jet.at(k)->Mass);
													v1=v1+v2+v3;
													histoSL_Boosted->Fill(v1.M());
												}
											}
										}
										
									}
								}
							}
					}
				}
			}
			//Minimum of chi method:

			double bestchi=1000.0;
			double topmass1,topmass2;
			solve.setLepton(lepton);
			solve.setNeutrino(neutrino);

			for(size_t i=0;i<jet.size();i++){
				//bjet from blv
                                if(jet.at(i)->BTag && (jet.at(i)->Eta)<2.5){
					TLorentzVector bjet1,bjet2,ljet1,ljet2;
					bjet1=makeTLorentzVector(jet.at(i));
					solve.setBJetA(bjet1);

					for(size_t j=0;j<jet.size();j++){
						//bjet from bqq'
                               			if(jet.at(j)->BTag && (jet.at(j)->Eta)<2.5 && j!=i){
							bjet2=makeTLorentzVector(jet.at(j));
							solve.setBJetB(bjet2);

							for(size_t k=0;k<jet.size();k++){
								//first lightjet
								if(!(jet.at(k)->BTag) && (jet.at(k)->Eta)<2.5){
									ljet1=makeTLorentzVector(jet.at(k));
									solve.setLightJetA(ljet1);
				
									for(size_t l=k;l<jet.size();l++){
										//second lightjet
										if(!(jet.at(l)->BTag) && (jet.at(l)->Eta)<2.5 && l!=k){

											ljet2=makeTLorentzVector(jet.at(l));
											solve.setLightJetB(ljet2);
											
											//compering chi2 with previus combination
											double chi2=solve.getChi2();
											if(chi2<bestchi){
												bestchi=chi2;
												topmass1=(bjet1+neutrino+lepton).M();
												topmass2=(bjet2+ljet1+ljet2).M();
											}
										}
									}
								}
							}
						}
					}
				}
			}
			//histogram of semileptonic decay
			histoSL->Fill(topmass1); //from blv
			histoSL->Fill(topmass2); //from bqq'
		}

		//Full hadronic - more sophisticated way:
		
		//Sorting out interesting full hadronic events:
		if(ljets>3){	//minimum for reconstruction ttbar
			//Minimum of chi method:
			double sigmabjet=10.0;
                        double sigmaljet=10.0;

                        double bestchi=1000.0;
                        double topmass1,topmass2,chi2;
			//First bjet
                        for(size_t i=0;i<jet.size();i++){
                                if(jet.at(i)->BTag && (jet.at(i)->Eta)<2.5){
                                        TLorentzVector top1,top2,bjet1,bjet2,ljet1,ljet2,ljet3,ljet4,w1,w2;
                                        bjet1.SetPtEtaPhiM(jet.at(i)->PT,jet.at(i)->Eta,jet.at(i)->Phi,jet.at(i)->Mass);
					//Second bjet
                                        for(size_t j=0;j<jet.size();j++){
                                                if(jet.at(j)->BTag && (jet.at(j)->Eta)<2.5 && j!=i){
                                                        bjet2.SetPtEtaPhiM(jet.at(j)->PT,jet.at(j)->Eta,jet.at(j)->Phi,jet.at(j)->Mass);
							//first lightjet
                                                        for(size_t k=0;k<jet.size();k++){
                                                                if(!(jet.at(k)->BTag) && (jet.at(k)->Eta)<2.5){
                                                                        ljet1.SetPtEtaPhiM(jet.at(k)->PT,jet.at(k)->Eta,jet.at(k)->Phi,jet.at(k)->Mass);
									//Second lightjet
                                                                        for(size_t l=k;l<jet.size();l++){
                                                                                if(!(jet.at(l)->BTag) && (jet.at(l)->Eta)<2.5 && l!=k){
                                                                                        ljet2.SetPtEtaPhiM(jet.at(l)->PT,jet.at(l)->Eta,jet.at(l)->Phi,jet.at(l)->Mass);
                                                                                        top1=bjet1+ljet1+ljet2;
                                                                                        w1=ljet1+ljet2;
											//Third lightjet
											for(size_t m=0;m<jet.size();m++){
												if(!(jet.at(m)->BTag) && (jet.at(m)->Eta)<2.5 && m!=k && m!=l){
													ljet3.SetPtEtaPhiM(jet.at(m)->PT,jet.at(m)->Eta,jet.at(m)->Phi,jet.at(m)->Mass);
													//fourth lightjet
													for(size_t n=m;n<jet.size();n++){
														if(!(jet.at(n)->BTag) && (jet.at(n)->Eta)<2.5 && n!=k && n!=l && n!=m){
															ljet4.SetPtEtaPhiM(jet.at(n)->PT,jet.at(n)->Eta,jet.at(n)->Phi,jet.at(n)->Mass);
															top2=bjet2+ljet3+ljet4;
                                                                                        				w2=ljet3+ljet4;
															chi2=(M_W-w1.M())*(M_W-w1.M())/((2*sigmaljet)*(2*sigmaljet)) + (M_W-w2.M())*(M_W-w2.M())/((2*sigmaljet)*(2*sigmaljet)) + (M_t-top1.M())*(M_t-top1.M())/((2*sigmaljet+sigmabjet)*(2*sigmaljet+sigmabjet)) + (M_t-top2.M())*(M_t-top2.M())/((2*sigmaljet+sigmabjet)*(2*sigmaljet+sigmabjet));
															if(chi2<bestchi){
                                                                                                				bestchi=chi2;
                                                                				                                topmass1=top1.M();
                                				                                                                topmass2=top2.M();
				                                                                                        }

														}
													}
												}
											}
                                                                                }
                                                                        }
                                                                }
                                                        }
                                                }
                                        }
                                }
                        }
			histoFH->Fill(topmass1);
			histoFH->Fill(topmass2);
		}


		//Full hadronic
		bool used[jet.size()+1]={false};	
		int bestint=0;
		int wjets[jet.size()+1];
		for(size_t i=0;i<jet.size();i++){
		
		TLorentzVector v1;
                v1.SetPtEtaPhiM(jet.at(i)->PT,jet.at(i)->Eta,jet.at(i)->Phi,jet.at(i)->Mass);
		if(!used[i]){
			wjets[i]=-1;
		}

		TLorentzVector bestfit;
		if(!(jet.at(i)->BTag) && (i!=jet.size()) && !used[i]){
			for(size_t j=i;j<jet.size();j++){
				if(j!=i && !used[j]){
					TLorentzVector v2;
					v2.SetPtEtaPhiM(jet.at(j)->PT,jet.at(j)->Eta,jet.at(j)->Phi,jet.at(j)->Mass);
					v2=v2+v1;
					if(fabs(v2.M()-80.39)<fabs(bestfit.M()-80.39)){ //selecting 2nd jet to get inv mass equal W as close as possiable
						bestfit=v2;
						bestint=j;
					}
				}

			}
			if(bestfit.M()>10){ //Not considering casses where invariant mass is extremaly small compere to W
				histoW->Fill(bestfit.M());
				used[bestint]=true;
				wjets[i]=bestint;
				wjets[bestint]=i;
			}
			
		}		

		}
		/*3 jet masses if W jets sorted*/
		int wjets2[jet.size()+1];
		for(size_t i=0;i<jet.size()+1;i++){
			wjets2[i]=wjets[i];
		}
		int bjetcount=0;
		for(size_t i=0;i<jet.size();i++){
		
		TLorentzVector v1;
		TLorentzVector v2;
		TLorentzVector v3;
		v1.SetPtEtaPhiM(jet.at(i)->PT,jet.at(i)->Eta,jet.at(i)->Phi,jet.at(i)->Mass);
		
		if(wjets[i]==-1){
			/*for(size_t j=i;j<jet.size();j++){
				if(j!=i && wjets[j]==-1){
					v2.SetPtEtaPhiM(jet.at(j)->PT,jet.at(j)->Eta,jet.at(j)->Phi,jet.at(j)->Mass);
					for(size_t k=j;k<jet.size();k++){
						if(k!=j && wjets[k]==-1){
							v3.SetPtEtaPhiM(jet.at(k)->PT,jet.at(k)->Eta,jet.at(k)->Phi,jet.at(k)->Mass);
							v3=v2+v1;
							histoW3C->Fill(v3.M());
						}
					}
				}
			}*/
		}else{
			size_t j=wjets[i];
			v2.SetPtEtaPhiM(jet.at(j)->PT,jet.at(j)->Eta,jet.at(j)->Phi,jet.at(j)->Mass);
			for(size_t k=i;k<jet.size();k++){
				if(k!=j && jet.at(k)->BTag){
					v3.SetPtEtaPhiM(jet.at(k)->PT,jet.at(k)->Eta,jet.at(k)->Phi,jet.at(k)->Mass);
					v3=v3+v2+v1;
					histoW3C->Fill(v3.M());
                                }
			}
		}
		/*t-quark mass(sorting b+W jets)*/
		TLorentzVector bestfit;
		if(jet.at(i)->BTag){
			bjetcount++;
			for(size_t j=i;j<jet.size();j++){
				if(wjets2[j]!=-1){
					size_t k=wjets2[j];
					v3.SetPtEtaPhiM(jet.at(k)->PT,jet.at(k)->Eta,jet.at(k)->Phi,jet.at(k)->Mass);
					v2.SetPtEtaPhiM(jet.at(j)->PT,jet.at(j)->Eta,jet.at(j)->Phi,jet.at(j)->Mass);
					v3=v3+v2+v1;
					if(fabs(v3.M()-173)<fabs(bestfit.M()-173)){
						bestfit=v3;
						bestint=j;
					}
				}
			}
			if(bestfit.M()>10){
				histoT->Fill(bestfit.M());
				wjets2[bestint]=-1;
			}
		}
		}
				
		 /* Or to fill the skim
		 */
		skimmedelecs.clear();
		for(size_t i=0;i<elecs.size();i++){
			//flat info
			elecPt=elecs.at(i)->PT;
			if(elecs.at(i)->PT < 20) continue;
			//or objects
			skimmedelecs.push_back(*elecs.at(i));
		}

		myskim->Fill();


		/*==SKIM==
		 * Access the branches of the skim
		 */
		//std::vector<Electron> * skimelecs=electrons.content();
		//for(size_t i=0;i<skimelecs->size();i++){
		//	histo->Fill(skimelecs->at(i).PT);
		//}
	}


	/*
	 * Must be called in the end, takes care of thread-safe writeout and
	 * call-back to the parent process
	 */
	processEndFunction();
}



void testanalyser::postProcess(){
	/*
	 * This function can be used to analyse the output histograms, e.g. extract a signal contribution etc.
	 * The function can also be called directly on an output file with the histograms, if
	 * RunOnOutputOnly = true
	 * is set in the analyser's config file
	 *
	 * This function also represents an example of how the output of the analyser can be
	 * read-back in an external program.
	 * Just include the sampleCollection.h header and follow the procedure below
	 *
	 */

	/*
	 * Here, the input file to the extraction of parameters from the histograms is the output file
	 * of the parallelised analysis.
	 * The sampleCollection class can also be used externally for accessing the output consistently
	 */
	d_ana::sampleCollection samplecoll;
	samplecoll.readFromFile(getOutPath());

	std::vector<TString> alllegends = samplecoll.listAllLegends();

	/*
	 * Example how to process the output.
	 * Usually, one would define the legendname of the histogram to be used here
	 * by hand, e.g. "signal" or "background".
	 * To make this example run in any case, I am using alllegends.at(0), which
	 * could e.g. be the signal legend.
	 *
	 * So in practise, the following would more look like
	 * samplecoll.getHistos("signal");
	 */
	if(alllegends.size()>0){
		d_ana::histoCollection histos=samplecoll.getHistos(alllegends.at(0));

		/*
		 * here, the histogram created in the analyze() function is selected and evaluated
		 * The histoCollection maintains ownership (you don't need to delete the histogram)
		 */
		const TH1* myplot=histos.getHisto("histoname1");
		/*const TH1* myplot3=histos.getHisto("histoname2");*/
		/*const TH1* myplot4=histos.getHisto("histoname3");*/

		std::cout << "(example output): the integral is " << myplot->Integral() <<std::endl;

		/*
		 * If the histogram is subject to changes, please clone it and take ownership
		 */

		TH1* myplot2=histos.cloneHisto("histoname1");

		/*
		 * do something with the histogram
		 */

		delete myplot2;
	}

	/*
	 * do the extraction here.
	 */



}



