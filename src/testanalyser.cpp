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

	d_ana::dBranchHandler<Jet>	   ParticleFlowJet04(tree(),"ParticleFlowJet04");
	d_ana::dBranchHandler<Jet>         ParticleFlowJet05(tree(),"ParticleFlowJet05");
	d_ana::dBranchHandler<Jet>	   ParticleFlowJet08(tree(),"ParticleFlowJet08");
	d_ana::dBranchHandler<Jet>         ParticleFlowJet12(tree(),"ParticleFlowJet12");
	d_ana::dBranchHandler<Jet>	   ParticleFlowJet15(tree(),"ParticleFlowJet15");
	
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
	TH1* histoname1=addPlot(new TH1D("histoname1","histotitle1",100,0,100)," "," ");
	//Masses of different jets
	TH1* histoAK4M=addPlot(new TH1D("AK4_mass","Mass of AK4 jets",100,0,500),"GeV","N_jet");
	TH1* histoAK5M=addPlot(new TH1D("AK5_mass","Mass of AK5 jets",100,0,500),"GeV","N_jet");
	TH1* histoAK8M=addPlot(new TH1D("AK8_mass","Mass of AK8 jets",100,0,500),"GeV","N_jet");
	TH1* histoAK12M=addPlot(new TH1D("AK12_mass","Mass of AK12 jets",100,0,500),"GeV","N_jet");
	TH1* histoAK15M=addPlot(new TH1D("AK15_mass","Mass of AK15 jets",100,0,500),"GeV","N_jet");
	//Pt of different jets
	TH1* histoAK5PT=addPlot(new TH1D("AK5_PT","Pt of AK5 jets",100,0,500),"GeV","N_jet");
	TH1* histoAK12PT=addPlot(new TH1D("AK12_PT","Pt of AK12 jets",100,0,500),"GeV","N_jet");
	TH1* histoLead=addPlot(new TH1D("Leading_PT","Pt of Leading jet",100,400,1000),"GeV","N_jet");
	//Boosted top quark analysis
	TH1* histoSL=addPlot(new TH1D("tquark_SL","t-quark m_{inv} mass in semileptonic decay",40,50,400),"M [GeV]","N");
        TH1* histoSL_H=addPlot(new TH1D("tquark_SL_H","t-quark m_{inv} mass in semileptonic decay(only bqq')",40,50,400),"M [GeV]","N");
        TH1* histoSL_L=addPlot(new TH1D("tquark_SL_L","t-quark m_{inv} mass in semileptonic decay(only blv)",40,50,400),"M [GeV]","N");
	
	
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

		//checking AK jets properties
		for(size_t i=0;i<ParticleFlowJet04.size();i++){
			histoAK4M->Fill(ParticleFlowJet04.at(i)->Mass);
		}
		for(size_t i=0;i<ParticleFlowJet05.size();i++){
                        histoAK5M->Fill(ParticleFlowJet05.at(i)->Mass);
			histoAK5PT->Fill(ParticleFlowJet05.at(i)->PT);
                }
		for(size_t i=0;i<ParticleFlowJet08.size();i++){
                        histoAK8M->Fill(ParticleFlowJet08.at(i)->Mass);
                }
		for(size_t i=0;i<ParticleFlowJet15.size();i++){
                        histoAK15M->Fill(ParticleFlowJet15.at(i)->Mass);
                }
		for(size_t i=0;i<ParticleFlowJet12.size();i++){
                        histoAK12M->Fill(ParticleFlowJet12.at(i)->Mass);
			histoAK12PT->Fill(ParticleFlowJet12.at(i)->PT);
                }



		/*
		 * Begin the event-by-event analysis
		 */		

		
		//Event analysis
		//Limits for event:
		//|eta|< 2.5 for jets(2.4 for b-jets)
		//At least 2 AK5 b-jets with pt>50GeV(at least one over 150)
		//Exactly 1 leading and 1 next-to leading jets
		//leading jets Pt>400GeV and next-to leading pt>150Gev
		//
		//|Eta|<2.1 for leptons
		//P_T>45GeV for leptons
		//Only one lepton(muon or electron)
		//
		//MET P_T>20GeV
		//Scalar sum of PT of MET and PT of lepton>150GeV
		//deltaPhi between lepton/jet and met has to be small enough


		//Checking is there enough fast b-jets
		size_t bjets=0;
		bool interesting=false;

		for(size_t i=0;i<ParticleFlowJet05.size();i++){
			if(eligibleBJet(ParticleFlowJet05.at(i))){
				if(ParticleFlowJet05.at(i)->PT >150){ //at least one b-jet has to be with Pt over 150GeV
					interesting=true;
				}
				bjets++;
			}
		}
		if(bjets<2){continue;}

		//Semieptonic decay
		
		//Class for calculating neutrino Pz
		pz_calculator neutrPZ;

		TLorentzVector lepton;
		double M_e = 0.511e-3;
		double M_mu = 0.10566;

		double M_lepton=M_e; //setting lepton mass to electron, changed further if that not the case
		double lepton_eta;
		double lepton_phi;
		bool second=false;

		//Searching events with only one lepton		
		//
		//Searching eligible electron
		for(size_t i=0;i<elecs.size();i++){
			if(eligibleLepton(elecs.at(i))){
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
			if(eligibleLepton(muontight.at(i))){
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
		if(!interesting){continue;}
		//Searching for leading and next-to leading jets
		TLorentzVector LeadingJet,LepLeadingJet;
		int leadingjets=0;
		for(size_t i=0;i<ParticleFlowJet12.size();i++){
                        if(eligibleLeadingJet(ParticleFlowJet12.at(i))){
				leadingjets++;
				if(ParticleFlowJet12.at(i)->PT > LeadingJet.Pt()){
					LepLeadingJet=LeadingJet;
					LeadingJet=makeTLorentzVector(ParticleFlowJet12.at(i));
	
				}else{
					LepLeadingJet=makeTLorentzVector(ParticleFlowJet12.at(i));
				}
			}
		}
		//There need to be exactly 1 leading and one next-to leading. Leading jet Pt has to be >400GeV
		if(leadingjets!=2 || LeadingJet.Pt()<400 || LeadingJet.M()<LepLeadingJet.M()){continue;}
		//There should be less than 1.2 angular distance between lepton and next-to leading jet
		if(dR(LepLeadingJet.Phi(),LepLeadingJet.Eta(),lepton_phi,lepton_eta)>1.2){continue;}
		//missing transversal energy
              	TLorentzVector metP4=(met.at(0)->P4());
		//limits for missing transversal momentum
		double metPT=metP4.Pt();
		if(metPT<20){continue;}
		if((metPT+lepton.Pt())<150){continue;}
		if(M_lepton==M_e && (deltaPhi(LepLeadingJet.Phi(),metP4.Phi())>(metPT/50) || deltaPhi(lepton_phi,metP4.Phi())>(metPT/50))){continue;}		
			
		//p_z calculating(for neutrino)
		neutrPZ.setLepton(lepton);
		neutrPZ.setMET(metP4);
		neutrPZ.setLeptonMass(M_lepton);			

		//neutrino (candidate)
		TLorentzVector neutrino(metP4.Px(),metP4.Py(),neutrPZ.getPz(),TMath::Sqrt(TMath::Power(metP4.Pt(),2)+TMath::Power(neutrPZ.getPz(),2)));
			
		LepLeadingJet+=neutrino; //adding neutrino vector to next-to leading jet as we assume this jet is blv product			

		//histograms			
		//leading jet Pt(only from bqq)
		histoLead->Fill(LeadingJet.Pt());
		//mass of leading jets
		histoSL->Fill(LeadingJet.M());
                histoSL->Fill(LepLeadingJet.M());
		//mass of main leading jet
		histoSL_H->Fill(LeadingJet.M());			
		//mass of subleading jet
		histoSL_L->Fill(LepLeadingJet.M());



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



