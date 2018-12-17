#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <TH1D.h>
#include <TH2F.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include "TRandom3.h"

#include "puweicalc.h"
#include "untuplizer.h"
#include "Utilities.h"

#include "Module/ElectronMatch.h"
#include "Module/PhotonMatch.h"
#include "Module/MuonMatch.h"

#include "Module/MuonSelection.h"
#include "ElectronSelections.h"
#include "ElectronSelections_HEEP.h"

#include "Roch_Moriond17/RoccoR.cc"


void xAna(const char* inpaths, string fileName = "minitree.root", Long64_t ev1 = 0, Long64_t ev2 = -1, 
	  float xs = 1, int isEleChannel = 1, int RunFile = 0, string RunProcess = "nominal") {

  /* The analysis entry point.
   *
   * inpaths = array of paths to files with ggNtuplizer's TTrees;
   * outpath = path to output root file;
   * [ev1, ev2) = event region to process (ev1 is included, ev2 is excluded).
   */
  char outpath[400];
  sprintf(outpath, "%s_%s.root", fileName.c_str(), RunProcess.c_str());

  // Get Trigger SFs histogram
  TFile *f_DoubleMu_Trig_Mu17LegSFs[4];
  f_DoubleMu_Trig_Mu17LegSFs[0] = new TFile("SFs/DoubleMuonTriggerSF/sf_Mu17Leg_Eta0to09.root", "READ");
  f_DoubleMu_Trig_Mu17LegSFs[1] = new TFile("SFs/DoubleMuonTriggerSF/sf_Mu17Leg_Eta09to12.root", "READ");
  f_DoubleMu_Trig_Mu17LegSFs[2] = new TFile("SFs/DoubleMuonTriggerSF/sf_Mu17Leg_Eta12to21.root", "READ");
  f_DoubleMu_Trig_Mu17LegSFs[3] = new TFile("SFs/DoubleMuonTriggerSF/sf_Mu17Leg_Eta21to24.root", "READ");
  TFile *f_DoubleMu_Trig_Mu8LegSFs[4];
  f_DoubleMu_Trig_Mu8LegSFs[0] = new TFile("SFs/DoubleMuonTriggerSF/sf_Mu8Leg_Eta0to09.root", "READ");
  f_DoubleMu_Trig_Mu8LegSFs[1] = new TFile("SFs/DoubleMuonTriggerSF/sf_Mu8Leg_Eta09to12.root", "READ");
  f_DoubleMu_Trig_Mu8LegSFs[2] = new TFile("SFs/DoubleMuonTriggerSF/sf_Mu8Leg_Eta12to21.root", "READ");
  f_DoubleMu_Trig_Mu8LegSFs[3] = new TFile("SFs/DoubleMuonTriggerSF/sf_Mu8Leg_Eta21to24.root", "READ");

  TH1F *h_DoubleMu_Trig_Mu17LegSFs[4];
  TH1F *h_DoubleMu_Trig_Mu8LegSFs[4];
  for(int i=0 ; i<4 ; i++){
    h_DoubleMu_Trig_Mu17LegSFs[i] = (TH1F*) f_DoubleMu_Trig_Mu17LegSFs[i]->Get("scale_factor");
    h_DoubleMu_Trig_Mu8LegSFs[i] = (TH1F*) f_DoubleMu_Trig_Mu8LegSFs[i]->Get("scale_factor");
  }

  //Get Muon ID SFs histogram
  TFile *f_MuonIDSFs = new TFile("SFs/HZZMuonIDSFs/final_HZZ_Moriond17Preliminary_v4.root", "READ");
  TH2F *h_MuonIDSFs = (TH2F*) f_MuonIDSFs->Get("FINAL");

  //Get Photon SFs histogram
  TFile *f_phoSFs = new TFile("SFs/PhotonSFs_Moriond16/egammaEffi.txt_EGM2DMVA90wp.root", "READ");
  TH2F *h_MVA_phoSFs = (TH2F*) f_phoSFs->Get("EGamma_SF2D");

  TRandom3 random; // gaus random number for energy smearing
  
  // open tree(s) with events to be processed
  TreeReader data(inpaths);
  
  // print out which type of variables should be associated with tree branches
  //data.Print();

  // prepare output tree
  TFile* fo = TFile::Open(outpath, "RECREATE");
  TTree* outtree_CountNegWei = new TTree("t_CountNegWei", "mini tree");
  TTree* outtree = new TTree("t", "mini tree");
  TTree* outtree_SetLimit = new TTree("t_Limit", "mini tree");
  TTree* outtree_BkgModeling = new TTree("t_BkgM", "mini tree");

  // do whathever preparations are necessary for, if MC information is present
  // pileup reweighting for MC
  PUWeightCalculator puCalc; // pu nominal
  PUWeightCalculator puCalc_Unup; // pu uncertainty (+4.6%)
  PUWeightCalculator puCalc_Undn; // pu uncertainty (-4.6%)
  if (data.HasMC()) {
    puCalc.Init("puweights_Moriond17/spring16/PU_histo_13TeV_GoldenJSON_69200nb.root"); 
    puCalc_Unup.Init("puweights_Moriond17/spring16/PU_histo_13TeV_GoldenJSON_72383nb.root"); 
    puCalc_Undn.Init("puweights_Moriond17/spring16/PU_histo_13TeV_GoldenJSON_66016nb.root"); 
  }
  Float_t puwei_Unup, puwei_Undn;

  // define the tool of rochaster muon pT correction
  RoccoR rc("Roch_Moriond17/rcdata.2016.v3"); // contains all variations

  // histogram for work flow
  TH1D * WorkFlow = new TH1D("WorkFlow", "WorkFlow", 10, 1, 11);

  // variables to be associated with the output tree branches
  Int_t LumiSection, RunNumber;
  Float_t puwei, generatorWeight;
  Float_t Ele1Pt, Ele1Eta, Ele1Phi;
  Float_t Ele2Pt, Ele2Eta, Ele2Phi;
  Float_t Muon1Pt, Muon1Eta, Muon1Phi;
  Float_t Muon2Pt, Muon2Eta, Muon2Phi;
  Float_t Photon1Pt, Photon1Eta, Photon1Phi, Photon1P;
  Float_t diEleM, diElePhotonM;
  Float_t diMuonM, diMuonPhotonM;
  Float_t diEleDr, diMuonDr, diMuonPt, diElePt, diMuonP, diEleP;
  Float_t diElePhoDr, diMuonPhoDr;
  Int_t MuonPassAccep, PhotonPassAccep;
  Int_t hasGoodGenProcess;
  Int_t MCPassHLT;
  Int_t Ele1Counter, Ele2Counter, Muon1Counter, Muon2Counter, Photon1Counter;
  Double_t Entries;
  Float_t PhoSFs; 
  Float_t Mu1IDSFs, Mu2IDSFs;
  Float_t DoubleMu_Trig_Mu17LegSFs, DoubleMu_Trig_Mu8LegSFs;
  Float_t DoubleMu_Trig_SF; //leg1  leg2
  Float_t totWeiForFit;
  Float_t totWeiForFit_UnMuSF_up, totWeiForFit_UnMuSF_dn, totWeiForFit_UnPhoSF_up, totWeiForFit_UnPhoSF_dn, totWeiForFit_UnPU_up, totWeiForFit_UnPU_dn;


  Int_t nVtx;

  if (data.HasMC()){
    outtree->Branch("hasGoodGenProcess", &hasGoodGenProcess);
    outtree->Branch("MCPassHLT", &MCPassHLT);
    outtree->Branch("puwei", &puwei);
    outtree->Branch("xs", &xs);
    outtree->Branch("generatorWeight", &generatorWeight);
    outtree->Branch("Entries", &Entries);
    outtree->Branch("PhoSFs", &PhoSFs);
    outtree->Branch("Mu1IDSFs", &Mu1IDSFs);
    outtree->Branch("Mu2IDSFs", &Mu2IDSFs);
    outtree->Branch("DoubleMu_Trig_Mu17LegSFs", &DoubleMu_Trig_Mu17LegSFs);
    outtree->Branch("DoubleMu_Trig_Mu8LegSFs", &DoubleMu_Trig_Mu8LegSFs);
    outtree->Branch("DoubleMu_Trig_SF", &DoubleMu_Trig_SF);
    outtree->Branch("totWeiForFit", &totWeiForFit);
    outtree->Branch("totWeiForFit_UnMuSF_up", &totWeiForFit_UnMuSF_up);
    outtree->Branch("totWeiForFit_UnMuSF_dn", &totWeiForFit_UnMuSF_dn);
    outtree->Branch("totWeiForFit_UnPhoSF_up", &totWeiForFit_UnPhoSF_up);
    outtree->Branch("totWeiForFit_UnPhoSF_dn", &totWeiForFit_UnPhoSF_dn);
    outtree->Branch("totWeiForFit_UnPU_up", &totWeiForFit_UnPU_up);
    outtree->Branch("totWeiForFit_UnPU_dn", &totWeiForFit_UnPU_dn);

    outtree_SetLimit->Branch("totWeiForFit", &totWeiForFit);
    outtree_CountNegWei->Branch("generatorWeight", &generatorWeight);
  }
  outtree->Branch("LumiSection", &LumiSection);
  outtree->Branch("RunNumber", &RunNumber);
  outtree->Branch("Ele1Pt", &Ele1Pt);
  outtree->Branch("Ele1Eta", &Ele1Eta);
  outtree->Branch("Ele1Phi", &Ele1Phi);
  outtree->Branch("Ele2Pt", &Ele2Pt);
  outtree->Branch("Ele2Eta", &Ele2Eta);
  outtree->Branch("Ele2Phi", &Ele2Phi);
  outtree->Branch("Muon1Pt", &Muon1Pt);
  outtree->Branch("Muon1Eta", &Muon1Eta);
  outtree->Branch("Muon1Phi", &Muon1Phi);
  outtree->Branch("Muon2Pt", &Muon2Pt);
  outtree->Branch("Muon2Eta", &Muon2Eta);
  outtree->Branch("Muon2Phi", &Muon2Phi);
  outtree->Branch("Photon1Pt", &Photon1Pt);
  outtree->Branch("Photon1Eta", &Photon1Eta);
  outtree->Branch("Photon1Phi", &Photon1Phi);
  outtree->Branch("Photon1P", &Photon1P);
  outtree->Branch("diEleDr", &diEleDr);
  outtree->Branch("diElePt", &diElePt);
  outtree->Branch("diEleP", &diEleP);
  outtree->Branch("diEleM", &diEleM); // two Electron invariant mass
  outtree->Branch("diElePhotonM", &diElePhotonM); // two Electron + one photon invariant mass
  outtree->Branch("diElePhoDr", &diElePhoDr); // two Electron + one photon invariant mass
  outtree->Branch("diMuonDr", &diMuonDr);
  outtree->Branch("diMuonPt", &diMuonPt);
  outtree->Branch("diMuonP", &diMuonP);
  outtree->Branch("diMuonM", &diMuonM); // two Muon invariant mass
  outtree->Branch("diMuonPhotonM", &diMuonPhotonM); // two Muon + one photon invariant mass
  outtree->Branch("diMuonPhoDr", &diMuonPhoDr);
  outtree->Branch("MuonPassAccep", &MuonPassAccep);
  outtree->Branch("PhotonPassAccep", &PhotonPassAccep);
  outtree->Branch("Ele1Counter", &Ele1Counter);
  outtree->Branch("Ele2Counter", &Ele2Counter);
  outtree->Branch("Muon1Counter", &Muon1Counter);
  outtree->Branch("Muon2Counter", &Muon2Counter);
  outtree->Branch("Photon1Counter", &Photon1Counter);
  outtree->Branch("nVtx", &nVtx);

  outtree_SetLimit->Branch("diMuonM", &diMuonM);
  outtree_BkgModeling->Branch("diMuonM", &diMuonM);

  Entries=data.GetEntriesFast();

  // random number for rochaster
  TRandom *Ran_Rochcor = new TRandom3();

  // event loop
  // update upper boundary of the region of events to process, if necessary
  if (ev2 < 0) ev2 = data.GetEntriesFast();
  if (ev2 > data.GetEntriesFast()) ev2 = data.GetEntriesFast();


  for (Long64_t ev = ev1; ev < ev2; ++ev) {

    /*
      Event Loop processing steps
                                                                                     ////////////////
      ///////////////////////////////      ///////////////     ////////////////  --> // Election   //     ///////////////     /////////////////
      // Initialized the variables //      // Apply HLT //     // in-time PU //      // Selection  // --> // Photon    //     // Calculating //
      // which is going to store   //  --> //           // --> //            //      ////////////////     // Selection // --> // Mll, Mllg,  //
      // into the minitree         //      //           //     //            //                           ///////////////     // dR, .. etc  //
      ///////////////////////////////      ///////////////     ////////////////  --> //////////////// -->                     /////////////////
                                                                                     // Muon       //
                                                                                     // Selection  //
                                                                                     ////////////////

     */

    // initialized the variable for store those information
    LumiSection=RunNumber=0;
    MuonPassAccep=PhotonPassAccep=0;
    Ele1Counter=Ele2Counter=Muon1Counter=Muon2Counter=Photon1Counter=0;
    Ele1Pt=Ele1Eta=Ele1Phi=0;
    Ele2Pt=Ele2Eta=Ele2Phi=0;
    Muon1Pt=Muon1Eta=Muon1Phi=0;
    Muon2Pt=Muon2Eta=Muon2Phi=0;
    Photon1Pt=Photon1Eta=Photon1Phi=Photon1P=0;
    diEleM=diElePhotonM=0;
    diEleDr=diElePt=0;
    diElePhoDr=diMuonPhoDr=0;
    diMuonM=diMuonPhotonM=0;
    diMuonDr=diMuonPt=0;
    diMuonP=diEleP=0;
    DoubleMu_Trig_Mu17LegSFs=DoubleMu_Trig_Mu8LegSFs=1;
    DoubleMu_Trig_SF=1;
    PhoSFs=1;
    Mu1IDSFs=Mu2IDSFs=1;
    totWeiForFit=1;
    totWeiForFit_UnMuSF_up=totWeiForFit_UnMuSF_dn=totWeiForFit_UnPhoSF_up=totWeiForFit_UnPhoSF_dn=totWeiForFit_UnPU_up=totWeiForFit_UnPU_dn=1;

    // SFs uncertaintes
    Float_t DoubleMu_Trig_Mu17LegSFs_Un=0, DoubleMu_Trig_Mu8LegSFs_Un=0;
    Float_t Mu1IDSFs_Un=0, Mu2IDSFs_Un=0;
    Float_t PhoSFs_Un=0;

    // print progress
    if (ev % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());
    
    data.GetEntry(ev);

    Int_t run = data.GetInt("run");
    RunNumber = run;
    LumiSection = data.GetInt("lumis");
    nVtx = data.GetInt("nVtx"); // nVtx has been defiined before

    if (data.HasMC()) { 
      Float_t genWeight = data.GetFloat("genWeight"); 
      if (genWeight < 0) {generatorWeight = -1;}
      if (genWeight >= 0) {generatorWeight = 1;}
    }

    if(data.HasMC() == 0){
      // HLT
      Long64_t HLTEleMuX = (Long64_t) data.GetLong64("HLTEleMuX");
      if (isEleChannel == 0){ // muon channel
	if (ev % 50000 == 0) {cout<<"applying HLT"<<endl;}
	if ( ( abs(HLTEleMuX>>14&1) != 1  ) && ( abs(HLTEleMuX>>15&1) != 1  ) ) continue; //Mu17_Mu8_TrkIsoVVL_DZ or TkMu8
      }
    }

    MCPassHLT=0;
    if(data.HasMC()){
      // HLT
      Long64_t HLTEleMuX = (Long64_t) data.GetLong64("HLTEleMuX");
      if (isEleChannel == 0){ // muon channel
	if (ev % 50000 == 0) {cout<<"applying HLT"<<endl;}
	if ( ( abs(HLTEleMuX>>14&1) == 1  ) || ( abs(HLTEleMuX>>15&1) == 1  ) ) MCPassHLT=1; //Mu17_Mu8_TrkIsoVVL_DZ or TkMu8
      }
    }

    // PU reweighting
    if (data.HasMC()) {
      float* puTrue = data.GetPtrFloat("puTrue");
      puwei = (float) puCalc.GetWeight(run, puTrue[1]); // in-time PU
      puwei_Unup = (float) puCalc_Unup.GetWeight(run, puTrue[1]); // PU uncertainty (+4.6%)
      puwei_Undn = (float) puCalc_Undn.GetWeight(run, puTrue[1]); // PU uncertainty (-4.6%)
    }

    Int_t    nMC = 0;
    Int_t*   mcPID = NULL;
    Int_t*   mcMomPID = NULL;
    Int_t*   mcGMomPID = NULL;
    UShort_t*   mcStatusFlag = NULL;
    Float_t* mcPt = NULL;
    Float_t* mcEta = NULL;
    Float_t* mcPhi = NULL;

    // genLevel
    if (data.HasMC()){

      nMC     = data.GetInt("nMC");
      mcPID   = data.GetPtrInt("mcPID");
      mcMomPID = data.GetPtrInt("mcMomPID");
      mcGMomPID = data.GetPtrInt("mcGMomPID");
      mcStatusFlag   = (UShort_t*) data.GetPtrShort("mcStatusFlag");
      mcPt = data.GetPtrFloat("mcPt");
      mcEta = data.GetPtrFloat("mcEta");
      mcPhi = data.GetPtrFloat("mcPhi");

      Int_t goodGenEventCounter=0;
      Int_t goodGenPho=0;
      hasGoodGenProcess=0;

      for (Int_t i=0; i<nMC; ++i) {

	if (((mcStatusFlag[i]>>0)&1) == 0 ) continue;
        if (((mcStatusFlag[i]>>1)&1) == 0 ) continue;

	// has mu mu photon
	if (abs(mcPID[i]) == 13 && abs(mcMomPID[i]) == 4900023 && abs(mcGMomPID[i]) == 25) goodGenEventCounter++;
	if (abs(mcPID[i]) == 22 && abs(mcMomPID[i]) == 25) goodGenPho=1;
      }// mc loop
      if (goodGenEventCounter == 2 && goodGenPho==1) hasGoodGenProcess=1;
    } // has MC

    // electron channel reco level study
    // electron
    Int_t nEle          = NULL;
    Float_t* elePt      = NULL;
    Float_t* eleEta     = NULL;
    Float_t* elePhi       = NULL;
    Float_t* eleSCEta     = NULL;
    Float_t* eleSCPhi     = NULL;
    UShort_t* eleIDbit    = NULL;
    Float_t* eleSCEn      = NULL;
    Float_t* eleHoverE   = NULL;
    Int_t* eleMissHits   = NULL;
    Float_t* eleD0       = NULL;
    
    // electron
    nEle          = data.GetInt("nEle");
    elePt      = data.GetPtrFloat("eleCalibPt");
    //elePt      = data.GetPtrFloat("elePt");
    eleEta     = data.GetPtrFloat("eleEta");
    elePhi       = data.GetPtrFloat("elePhi");
    eleSCEta     = data.GetPtrFloat("eleSCEta");
    eleSCPhi     = data.GetPtrFloat("eleSCPhi");
    eleIDbit    = (UShort_t*) data.GetPtrShort("eleIDbit");
    eleSCEn      = data.GetPtrFloat("eleSCEn");
    eleHoverE   = data.GetPtrFloat("eleHoverE");
    eleMissHits   = data.GetPtrInt("eleMissHits");
    eleD0       = data.GetPtrFloat("eleD0");
    
    Int_t ele_index[2]={-999, -999};
    Int_t eleAcceptance = 0;
    
    Int_t LeadingEleFound = 0; // whether Leading muon is found or not (0=not found , 1=found)
    for (int i=0 ; i<nEle ; ++i) {
      
      if (fabs(eleSCEta[i]) > 1.4442 && fabs(eleSCEta[i]) < 1.566) continue; // Eta cut
      if (fabs(eleSCEta[i]) > 2.5) continue; // Eta cut
      
      if (LeadingEleFound == 1){
	if (elePt[i]<5.) continue; // leading electron pt cut is 5 GeV
      }
      
      if (LeadingEleFound == 0){
	if (elePt[i]<30.) continue; // leading electron pt cut is 25 GeV
	LeadingEleFound=1;
      }
      
      // count acceptance efficiency      
      eleAcceptance++;
      
      if (eleAcceptance == 2){
	WorkFlow->Fill(3); // pass two electron preselection
      }
      
      //if (abs(eleIDbit[i]>>4&1) == 0) continue; // HEEP cut
      //if (isHEEPele(data, i) == 0) continue;
      if (ElectronIDCutBasedRunII(1, data, i) == 0) continue; // mini isolation from Shilpi # 08/18/16  
      
      // sort
      if (ele_index[0] == -999 && ele_index[1] == -999) { ele_index[0]=i; }
      else if(ele_index[0] != -999 && ele_index[1] == -999) {
	if (elePt[i] >  elePt[ele_index[0]]) { ele_index[1] = ele_index[0]; ele_index[0] = i; }
	else if (elePt[i] <= elePt[ele_index[0]]) { ele_index[1] = i; }
      }
      else {
	if (elePt[i] >  elePt[ele_index[0]]    )                                  { ele_index[1] = ele_index[0]; ele_index[0] = i; }
	else if (elePt[i] <= elePt[ele_index[0]] &&  elePt[i] > elePt[ele_index[1]])  { ele_index[1] = i;                                }
	else if (elePt[i] <= elePt[ele_index[0]] &&  elePt[i] < elePt[ele_index[1]]) continue;
      } // end of sort
    } // end of electron loop
    
    if (ele_index[0] != -999){
      Ele1Pt=elePt[ele_index[0]]; Ele1Eta=eleEta[ele_index[0]]; Ele1Phi=elePhi[ele_index[0]];
      Ele1Counter=1;
    }
    if (ele_index[1] != -999){
      Ele2Pt=elePt[ele_index[1]]; Ele2Eta=eleEta[ele_index[1]]; Ele2Phi=elePhi[ele_index[1]];
      Ele2Counter=1;
    } // the section pick up two electrons after sorting     
    
      
    // Muon selection
    
    Int_t    nMu             = data.GetInt("nMu");
    Float_t* muPt            = data.GetPtrFloat("muPt");
    Float_t* muEta           = data.GetPtrFloat("muEta");
    Float_t* muPhi           = data.GetPtrFloat("muPhi");
    Int_t*   muCharge        = data.GetPtrInt("muCharge");
    Int_t*   muType          = data.GetPtrInt("muType");   
    Float_t* muChi2NDF       = data.GetPtrFloat("muChi2NDF"); // chi square / Number Of Degree
    Int_t*   muStations      = data.GetPtrInt("muStations");
    Int_t*   muMatches       = data.GetPtrInt("muMatches");
    Float_t* muPFChIso       = data.GetPtrFloat("muPFChIso"); // Sum of Charged Hadron Pt
    Float_t* muPFPhoIso      = data.GetPtrFloat("muPFPhoIso"); // Sum of Photon Et
    Float_t* muPFNeuIso      = data.GetPtrFloat("muPFNeuIso"); // Sum of Neutral Hadron Et
    Float_t* muPFPUIso       = data.GetPtrFloat("muPFPUIso"); // Sum of PU(pile up) Pt
    Float_t* muIsoTrk        = data.GetPtrFloat("muIsoTrk"); // Tracker-base isolation cut
    Float_t* muD0            = data.GetPtrFloat("muD0"); // Impact parameter about xy plane
    Float_t* muDz            = data.GetPtrFloat("muDz"); // Impact parameter about z axis
    Int_t*   muTrkLayers     = data.GetPtrInt("muTrkLayers");
    Float_t* muSIP           = data.GetPtrFloat("muSIP"); // SIP
    Int_t*   muBestTrkType    = data.GetPtrInt("muBestTrkType");
    
    Int_t mu_index[2]={-999, -999};
    
    Int_t LeadingMuFound = 0; // whether Leading muon is found or not (0=not found , 1=found)
    // muon loop, selection up the event within acceptance area only. Do mc true muon matching

    double u1, u2; // Random number
    double muPt_SF;

    for (int i=0 ; i<nMu ; ++i) {

      
      /////////////////// start to rochaster muon pT correction
      muPt_SF=1;
      u1 = Ran_Rochcor->Rndm();
      u2 = Ran_Rochcor->Rndm();
      if(data.HasMC()){
	  muPt_SF = rc.kScaleAndSmearMC(muCharge[i], muPt[i], muEta[i], muPhi[i], muTrkLayers[i], u1, u2, 0, 0);
	}
      else{
	  muPt_SF = rc.kScaleDT(muCharge[i], muPt[i], muEta[i], muPhi[i], 0, 0);
	}
      muPt[i] = muPt[i] * muPt_SF; //muon pT correction
      //////////////////// finish rochaster muon pT correction
      

      // muon acceptance
      if (fabs(muEta[i]) > 2.4) continue; // Eta cut      
      if (LeadingMuFound == 1){
	if (muPt[i]<10.) continue; // trailing muon in EB pt cut is 10 GeV
      }
      if (LeadingMuFound == 0){
	if (muPt[i]<20.) continue; // leading muon in EB pt cut is 20 GeV
	LeadingMuFound=1;
      }
      
      MuonPassAccep ++;

      /*
      // Loose muon ID
      if ( ((muType[i] >> 5) &1) != 1) continue; // isPF Muon?
      if ( ( ((muType[i] >> 1) &1) != 1 ) && ( ((muType[i] >> 2) &1) != 1 ) ) continue; // isTracker Muon or isGlobal Muon
      */

      // HZZ analysis Muon Tight ID
      if ( ((muType[i] >> 1) &1) == 0 && ( ((muType[i] >> 2) &1) == 0 || muMatches[i] <= 0 ) ) continue; // HZZ analysis Loose Muon ID
      if ( muBestTrkType[i] == 2 ) continue; // HZZ analysis Loose Muon ID
      if ( ((muType[i] >> 5) &1) == 0 && muPt[i] < 200) continue;
      if ( ((muType[i] >> 5) &1) == 0 && !TrkHighPtID(data, i) && muPt[i] > 200 ) continue;


      // muon Dz and D0 cut
      if (muDz[i] > 1 || muD0[i]> 0.5) continue;

      // muon SIP cut at 4 or 10
      if (muSIP[i] > 10) continue;

      // PF Isolation cut
      Float_t MuPFIso = 0;
      MuPFIso = (muPFChIso[i] +
		 TMath::Max(0.0, (muPFNeuIso[i] + muPFPhoIso[i] - 0.5*muPFPUIso[i]) ) ) / muPt[i] ;
      if (MuPFIso > 0.35) continue;

      /////////////
      // Trk Iso //
      /////////////
      Float_t MuTrkIso = 0;
      // The trk Isolation is more tricky, it need to do the in-cone muon pT remove,
      // loop all the muon and check every muon's Dr. If the others muon is within 0.3 cone, substract the pT.
      MuTrkIso = muIsoTrk[i];
      for (int j=0 ; j<nMu ; j++){
	if (j == i) continue; // the muon is not itself
	if (deltaR(muEta[i], muPhi[i], muEta[j], muPhi[j]) < 0.3){
	  MuTrkIso = MuTrkIso - muPt[j]; // substract the muon pT if that muon within Trk Iso cone (0.3)
	}
      } // end of second muon loop
      MuTrkIso = MuTrkIso/muPt[i] ; // divide by the itself's pT
      //if (MuTrkIso > 0.1) continue; // tracker Isolation

      // sort
      if (mu_index[0] == -999 && mu_index[1] == -999) { mu_index[0]=i; }
      else if(mu_index[0] != -999 && mu_index[1] == -999) {
	if (muPt[i] >  muPt[mu_index[0]]) { mu_index[1] = mu_index[0]; mu_index[0] = i; }
	else if (muPt[i] <= muPt[mu_index[0]]) { mu_index[1] = i; }
      }
      else {
	if (muPt[i] >  muPt[mu_index[0]]    )                                  { mu_index[1] = mu_index[0]; mu_index[0] = i; }
	else if (muPt[i] <= muPt[mu_index[0]] &&  muPt[i] > muPt[mu_index[1]])  { mu_index[1] = i;                                }
	else if (muPt[i] <= muPt[mu_index[0]] &&  muPt[i] < muPt[mu_index[1]]) continue;
      } // end of sort
    } // end of muon loop
    
    if (mu_index[0] != -999){
      Muon1Pt=muPt[mu_index[0]]; Muon1Eta=muEta[mu_index[0]]; Muon1Phi=muPhi[mu_index[0]];
      Muon1Counter=1;
    }
    if (mu_index[1] != -999){
      Muon2Pt=muPt[mu_index[1]]; Muon2Eta=muEta[mu_index[1]]; Muon2Phi=muPhi[mu_index[1]];
      Muon2Counter=1;
    } // the section for pick up mcTrue muon after sorting
    
    int Mu1_TrigSFs_binPt=-999, Mu1_TrigSFs_binEta=-999;
    int Mu2_TrigSFs_binPt=-999, Mu2_TrigSFs_binEta=-999;

    if (fabs(Muon1Eta) < 0.9) { Mu1_TrigSFs_binEta=0; }
    else if (fabs(Muon1Eta) > 0.9 && fabs(Muon1Eta) < 1.2) { Mu1_TrigSFs_binEta=1; }
    else if (fabs(Muon1Eta) > 1.2 && fabs(Muon1Eta) < 2.1) { Mu1_TrigSFs_binEta=2; }
    else if (fabs(Muon1Eta) > 2.1 && fabs(Muon1Eta) < 2.4) { Mu1_TrigSFs_binEta=3; }

    if (fabs(Muon2Eta) < 0.9) { Mu2_TrigSFs_binEta=0; }
    else if (fabs(Muon2Eta) > 0.9 && fabs(Muon2Eta) < 1.2) { Mu2_TrigSFs_binEta=1; }
    else if (fabs(Muon2Eta) > 1.2 && fabs(Muon2Eta) < 2.1) { Mu2_TrigSFs_binEta=2; }
    else if (fabs(Muon2Eta) > 2.1 && fabs(Muon2Eta) < 2.4) { Mu2_TrigSFs_binEta=3; }

    Mu1_TrigSFs_binPt = h_DoubleMu_Trig_Mu17LegSFs[Mu1_TrigSFs_binEta]->FindBin(Muon1Pt);
    Mu2_TrigSFs_binPt = h_DoubleMu_Trig_Mu8LegSFs[Mu2_TrigSFs_binEta]->FindBin(Muon2Pt);

    DoubleMu_Trig_Mu17LegSFs    = h_DoubleMu_Trig_Mu17LegSFs[Mu1_TrigSFs_binEta]->GetBinContent(Mu1_TrigSFs_binPt);
    DoubleMu_Trig_Mu8LegSFs     = h_DoubleMu_Trig_Mu8LegSFs[Mu2_TrigSFs_binEta]->GetBinContent(Mu2_TrigSFs_binPt);
    DoubleMu_Trig_Mu17LegSFs_Un = h_DoubleMu_Trig_Mu17LegSFs[Mu1_TrigSFs_binEta]->GetBinError(Mu1_TrigSFs_binPt);
    DoubleMu_Trig_Mu8LegSFs_Un  = h_DoubleMu_Trig_Mu8LegSFs[Mu2_TrigSFs_binEta]->GetBinError(Mu2_TrigSFs_binPt);

    int Mu1_binPtEta=-999, Mu2_binPtEta=-999;
    Mu1_binPtEta=h_MuonIDSFs->FindBin(Muon1Eta, Muon1Pt); // FindBin(bin x, bin y)
    Mu2_binPtEta=h_MuonIDSFs->FindBin(Muon2Eta, Muon2Pt); // get pt eta bin
    
    Mu1IDSFs    = h_MuonIDSFs->GetBinContent(Mu1_binPtEta);
    Mu2IDSFs    = h_MuonIDSFs->GetBinContent(Mu2_binPtEta);
    Mu1IDSFs_Un = h_MuonIDSFs->GetBinError(Mu1_binPtEta);
    Mu2IDSFs_Un = h_MuonIDSFs->GetBinError(Mu2_binPtEta);

    if (DoubleMu_Trig_Mu17LegSFs == 0) {DoubleMu_Trig_Mu17LegSFs=1; DoubleMu_Trig_Mu17LegSFs_Un=0;}
    if (DoubleMu_Trig_Mu8LegSFs == 0)  {DoubleMu_Trig_Mu8LegSFs=1; DoubleMu_Trig_Mu8LegSFs_Un=0;}

    if (Mu1IDSFs == 0) {Mu1IDSFs=1; Mu1IDSFs_Un=0;}
    if (Mu2IDSFs == 0) {Mu2IDSFs=1; Mu2IDSFs_Un=0;}
  

    // photon selection    
    //photon

    Float_t* phoScale_stat_up = NULL;
    Float_t* phoScale_stat_dn = NULL;
    Float_t* phoScale_syst_up = NULL;
    Float_t* phoScale_syst_dn = NULL;
    Float_t* phoScale_gain_up = NULL;
    Float_t* phoScale_gain_dn = NULL;
    Float_t* phoResol_rho_up = NULL;   
    Float_t* phoResol_rho_dn = NULL;
    Float_t* phoResol_phi_up = NULL;
    Float_t* phoResol_phi_dn = NULL;

    Int_t nPho        = NULL;
    Int_t* phoEleVeto = NULL;
    Float_t* phoEt    = NULL;
    Float_t* phoEta   = NULL;
    Float_t* phoPhi   = NULL;
    Float_t* phoSCEta = NULL;
    Float_t* phoSCPhi = NULL;
    Float_t* phoIDMVA = NULL;
    UShort_t* phoIDbit = (UShort_t*) NULL;
    
    Float_t* phoPFChIso = NULL;
    Float_t* phoPFPhoIso = NULL;
    Float_t* phoPFNeuIso = NULL;
    Float_t* phoHoverE = NULL;
    Float_t* phoPFChWorstIso = NULL;
    Float_t* phoR9 = NULL;
   
    // photon
    ///////////////////////////////////////
    // Photon Energy (Nominal and SysUn) //
    ///////////////////////////////////////
    phoEt    = data.GetPtrFloat("phoCalibEt");
    phoScale_stat_up = data.GetPtrFloat("phoScale_stat_up");
    phoScale_stat_dn = data.GetPtrFloat("phoScale_stat_dn");
    phoScale_syst_up = data.GetPtrFloat("phoScale_syst_up");
    phoScale_syst_dn = data.GetPtrFloat("phoScale_syst_dn");
    phoScale_gain_up = data.GetPtrFloat("phoScale_gain_up");
    phoScale_gain_dn = data.GetPtrFloat("phoScale_gain_dn");
    phoResol_rho_up = data.GetPtrFloat("phoResol_rho_up");
    phoResol_rho_dn = data.GetPtrFloat("phoResol_rho_dn");
    phoResol_phi_up = data.GetPtrFloat("phoResol_phi_up");
    phoResol_phi_dn = data.GetPtrFloat("phoResol_phi_dn");
  
    nPho        = data.GetInt("nPho");
    phoEleVeto = data.GetPtrInt("phoEleVeto");
    phoEta   = data.GetPtrFloat("phoEta");
    phoPhi   = data.GetPtrFloat("phoPhi");
    phoSCEta = data.GetPtrFloat("phoSCEta");
    phoSCPhi = data.GetPtrFloat("phoSCPhi");
    phoIDMVA = data.GetPtrFloat("phoIDMVA");
    phoIDbit = (UShort_t*) data.GetPtrShort("phoIDbit");
      
    phoHoverE = data.GetPtrFloat("phoHoverE");
    phoPFChIso = data.GetPtrFloat("phoPFChIso");
    phoPFPhoIso = data.GetPtrFloat("phoPFPhoIso");
    phoPFNeuIso = data.GetPtrFloat("phoPFNeuIso");
    phoPFChWorstIso = data.GetPtrFloat("phoPFChWorstIso");
    phoR9 = data.GetPtrFloat("phoR9");
    
    Int_t phoPreSele=0, phoDr=0, phoVeto=0;
    for (Int_t i=0 ; i<nPho ; ++i){
      
      if (Photon1Counter == 1) continue;
      
      if (RunProcess == "phoScal_stat_up"){ phoEt[i] = phoEt[i] * phoScale_stat_up[i]; }
      if (RunProcess == "phoScal_stat_dn"){ phoEt[i] = phoEt[i] * phoScale_stat_dn[i]; }
      if (RunProcess == "phoScal_syst_up"){ phoEt[i] = phoEt[i] * phoScale_syst_up[i]; }
      if (RunProcess == "phoScal_syst_dn"){ phoEt[i] = phoEt[i] * phoScale_syst_dn[i]; }
      if (RunProcess == "phoScal_gain_up"){ phoEt[i] = phoEt[i] * phoScale_gain_up[i]; }
      if (RunProcess == "phoScal_gain_dn"){ phoEt[i] = phoEt[i] * phoScale_gain_dn[i]; }
      if (RunProcess == "phoSmear_rho_up"){ phoEt[i] = phoEt[i] * random.Gaus(1, phoResol_rho_up[i]); }
      if (RunProcess == "phoSmear_rho_dn"){ phoEt[i] = phoEt[i] * random.Gaus(1, phoResol_rho_dn[i]); }
      if (RunProcess == "phoSmear_phi_up"){ phoEt[i] = phoEt[i] * random.Gaus(1, phoResol_phi_up[i]); }
      if (RunProcess == "phoSmear_phi_dn"){ phoEt[i] = phoEt[i] * random.Gaus(1, phoResol_phi_dn[i]); }

      // Photon acceptance cut
      if (phoEt[i]<15.) continue;
      if (fabs(phoSCEta[i]) > 1.4442 && fabs(phoSCEta[i]) < 1.56) continue;
      if (fabs(phoSCEta[i]) > 2.5) continue;
      PhotonPassAccep++;

      phoPreSele=1;
      
      if (Ele1Counter == 1 || Ele2Counter == 1){
	if (deltaR(Ele1Eta, Ele1Phi, phoEta[i], phoPhi[i]) < 0.8 ||
	    deltaR(Ele2Eta, Ele2Phi, phoEta[i], phoPhi[i]) < 0.8) continue;
      } // photon , electron dR cut

      if (Muon1Counter == 1 || Muon2Counter == 1){
	if (deltaR(Muon1Eta, Muon1Phi, phoEta[i], phoPhi[i]) < 0.8 ||
            deltaR(Muon2Eta, Muon2Phi, phoEta[i], phoPhi[i]) < 0.8) continue;
      } // photon , muon dR cut

    
      phoDr=1;
      if (phoEleVeto[i] == 0) continue;
      phoVeto=1;
      
      {
	if (fabs(phoSCEta[i]) < 1.4442 && phoIDMVA[i] < 0.2) continue;
	if (fabs(phoSCEta[i]) > 1.566 && phoIDMVA[i] < 0.2) continue;
      }
      
      if (RunFile == 1 && PhoMatch(data, phoEta[i], phoPhi[i]) == 1 ) continue;

      Photon1Eta=phoEta[i]; Photon1Phi=phoPhi[i];
      Photon1Pt=phoEt[i];
      Photon1Counter=1;

      int binEtaEt=-999;
      binEtaEt  = h_MVA_phoSFs->FindBin(phoSCEta[i], Photon1Pt);
      PhoSFs    = h_MVA_phoSFs->GetBinContent(binEtaEt);
      PhoSFs_Un = h_MVA_phoSFs->GetBinError(binEtaEt);

      break;
      
    } // the loop for pick up photon
    
    if(data.HasMC()){

      DoubleMu_Trig_SF = DoubleMu_Trig_Mu17LegSFs * DoubleMu_Trig_Mu8LegSFs;

      totWeiForFit=DoubleMu_Trig_SF*Mu1IDSFs*Mu2IDSFs*PhoSFs*puwei; // this is for signal sample, that is why it don't have MC weight

      // Uncertainties sorce from Mu and Pho SFs
      totWeiForFit_UnMuSF_up  = DoubleMu_Trig_SF*(Mu1IDSFs+Mu1IDSFs_Un)*(Mu2IDSFs+Mu2IDSFs_Un)*PhoSFs*puwei;
      totWeiForFit_UnMuSF_dn  = DoubleMu_Trig_SF*(Mu1IDSFs-Mu1IDSFs_Un)*(Mu2IDSFs-Mu2IDSFs_Un)*PhoSFs*puwei;
      totWeiForFit_UnPhoSF_up = DoubleMu_Trig_SF*Mu1IDSFs*Mu2IDSFs*(PhoSFs+PhoSFs_Un)*puwei;
      totWeiForFit_UnPhoSF_dn = DoubleMu_Trig_SF*Mu1IDSFs*Mu2IDSFs*(PhoSFs-PhoSFs_Un)*puwei;
      totWeiForFit_UnPU_up    = DoubleMu_Trig_SF*Mu1IDSFs*Mu2IDSFs*PhoSFs*puwei_Unup;
      totWeiForFit_UnPU_dn    = DoubleMu_Trig_SF*Mu1IDSFs*Mu2IDSFs*PhoSFs*puwei_Undn;
    }
  
    // diLepton Dr
    diEleDr=deltaR(Ele1Eta, Ele1Phi, Ele2Eta, Ele2Phi);  // electron pair dR
    diMuonDr=deltaR(Muon1Eta, Muon1Phi, Muon2Eta, Muon2Phi);  // muon pair dR

    // diLepton Mass and Three body invariant mass
    TLorentzVector fourVec_Ele1, fourVec_Ele2, fourVec_Muon1, fourVec_Muon2,  fourVec_Mee, fourVec_Mmumu, fourVec_Photon1, 
      fourVec_Meeg, fourVec_Mmumug;

    fourVec_Ele1.SetPtEtaPhiM(Ele1Pt, Ele1Eta, Ele1Phi, 0.511*0.001);
    fourVec_Ele2.SetPtEtaPhiM(Ele2Pt, Ele2Eta, Ele2Phi, 0.511*0.001);
    fourVec_Muon1.SetPtEtaPhiM(Muon1Pt, Muon1Eta, Muon1Phi, 0.105);
    fourVec_Muon2.SetPtEtaPhiM(Muon2Pt, Muon2Eta, Muon2Phi, 0.105);
    fourVec_Photon1.SetPtEtaPhiM(Photon1Pt, Photon1Eta, Photon1Phi, 0);

    // FSR recovery
    //--------------------
    TLorentzVector fourVec_FSRpho;
    if (Muon1Counter == 1 && Muon2Counter == 1){
      for(int i=0 ; i<nMC ; i++){
	///////////////////////////////
	// PF photon ID need to wait //
	///////////////////////////////
	if (abs(mcPID[i]) != 22) continue;
	//if ( ((phoIDbit[i]>>0) &1) == 0) continue; // if photon didn;t pass loose photon ID
	//if (phoEt[i]<2) continue; // pt cut
	//if (fabs(phoEta[i]) < 2.4) continue; // eta cut
	// photon Isolation
	//Float_t PhoPFIso = 0;
	//PhoPFIso = (phoPFChIso[i] + (phoPFNeuIso[i] + phoPFPhoIso[i]) ) / phoEt[i] ;
	//if (PhoPFIso > 1.8) continue;

	if (deltaR(Muon1Eta, Muon1Phi, mcEta[i], mcPhi[i]) < 0.5 ) {
	  fourVec_FSRpho.SetPtEtaPhiM(mcPt[i], mcEta[i], mcPhi[i], 0);
	  fourVec_Muon1=fourVec_Muon1+fourVec_FSRpho;
	  continue; // if this photon has been add into first photon, continue;
	}
	
	if (deltaR(Muon2Eta, Muon2Phi, mcEta[i], mcPhi[i]) < 0.5 ) {
	  fourVec_FSRpho.SetPtEtaPhiM(mcPt[i], mcEta[i], mcPhi[i], 0);
	  fourVec_Muon2=fourVec_Muon2+fourVec_FSRpho;
	  continue; // this photon has been add into second photon, continue;
	}
      }
    } // end of FSR photon recovery

    // Reset Muon pT eta phi after photon FSR recovery
    Muon1Pt = fourVec_Muon1.Pt(); Muon1Eta = fourVec_Muon1.Eta(); Muon1Phi = fourVec_Muon1.Phi();
    Muon2Pt = fourVec_Muon2.Pt(); Muon2Eta = fourVec_Muon2.Eta(); Muon2Phi = fourVec_Muon2.Phi();

    fourVec_Mee=fourVec_Ele1+fourVec_Ele2;
    fourVec_Mmumu=fourVec_Muon1+fourVec_Muon2;
    fourVec_Meeg=fourVec_Ele1+fourVec_Ele2+fourVec_Photon1;
    fourVec_Mmumug=fourVec_Muon1+fourVec_Muon2+fourVec_Photon1;
    
    Photon1P=fourVec_Photon1.P();

    diEleM=fourVec_Mee.M();
    diMuonM=fourVec_Mmumu.M();

    diElePt =fourVec_Mee.Pt();
    diMuonPt=fourVec_Mmumu.Pt();

    diEleP = fourVec_Mee.P();
    diMuonP= fourVec_Mmumu.P();

    diElePhotonM=fourVec_Meeg.M();
    diMuonPhotonM=fourVec_Mmumug.M();

    diElePhoDr=deltaR(fourVec_Mee.Eta(), fourVec_Mee.Phi(), Photon1Eta, Photon1Phi);
    diMuonPhoDr=deltaR(fourVec_Mmumu.Eta(), fourVec_Mmumu.Phi(), Photon1Eta, Photon1Phi);

    //if (Muon1Counter == 1 && Muon2Counter == 1 && Photon1Counter == 1) {outtree->Fill();}
    outtree->Fill();

    if(data.HasMC()){
      outtree_CountNegWei->Fill();
    }

    if (data.HasMC()){
      if (hasGoodGenProcess == 1 &&Muon1Counter == 1 && Muon2Counter == 1 && Photon1Counter == 1 && diMuonPhotonM > 123 && diMuonPhotonM < 127 && diMuonM < 70)
	{outtree_SetLimit->Fill();}
    }
    else{
      if (Muon1Counter == 1 && Muon2Counter == 1 && Photon1Counter == 1 && diMuonPhotonM > 123 && diMuonPhotonM < 127 && diMuonM < 70) {outtree_SetLimit->Fill();}
    }

    if (Muon1Counter == 1 && Muon2Counter == 1 && Photon1Counter == 1 && ((diMuonPhotonM > 100 && diMuonPhotonM < 123) || (diMuonPhotonM > 127 && diMuonPhotonM < 180)) && 
	( (diMuonM > 0.21 && diMuonM < 2.9) || (diMuonM > 3.1 && diMuonM < 9.5) || (diMuonM > 10.5 && diMuonM < 70) )){outtree_BkgModeling->Fill(); }

  } // event loop
  
  fprintf(stderr, "Processed all events\n");

  // flush caches, close the output
  outtree->Write("", TObject::kOverwrite);
  outtree_SetLimit->Write("", TObject::kOverwrite);
  outtree_CountNegWei->Write("", TObject::kOverwrite);
  outtree_BkgModeling->Write("", TObject::kOverwrite);
  delete outtree;
  delete outtree_SetLimit;
  delete outtree_CountNegWei;
  delete outtree_BkgModeling;
  delete fo;
  delete h_MVA_phoSFs;
  //delete f_MuonIDSFs_RunBCDEF;
  delete h_MuonIDSFs;
}
