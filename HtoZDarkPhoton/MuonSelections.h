#include <TLorentzVector.h>
#include "external/rochcor2012jan22.h"  // latest
#include "TMath.h"

// #include "rochcor2012v2.h"  // previous

using namespace std;

//______________________________________________________________________________
void select_muons(TreeReader &data, vector<int> &accepted,
                  vector<float> &idLoose, vector<float> &idTight, 
		  vector<float> &idminiTight)
{
  /* Identification, isolation and energy scale correction for muons (for 2012).
   *
   * NOTE: the energy scale correction is applied to muPt and muInnerPtErr
   * branches only. Other energy-related branches (e.g muPz) are not corrected!
   *
   * NOTE: for multi-muon analysis, the loose muon ID should be complemented
   * with DeltaR > 0.02 cut between the muon pairs in order to suppress
   * contribution from split tracks.
   *
   * Documentation on the muon ID and isolation:
   * https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId?rev=49
   *
   * Documentation on the muon energy correction:
   * https://twiki.cern.ch/twiki/bin/viewauth/CMS/RochcorMuon?rev=2
   *
   * data = handle providing access to an input event;
   * accepted = array to fill with indices of accepted muons;
   * idLoose, idTight = arrays to fill synchronously with results of cut-based
   *   loose and tight identification criteria (0=not passed, 1=passed).
   */

  // load necessary tree branches
  Int_t  nMu            = data.GetInt("nMu");
  float* muPt           = data.GetPtrFloat("muPt");
  float* muEta          = data.GetPtrFloat("muEta");
  Int_t* muType         = data.GetPtrInt("muType");
  float* muChi2NDF      = data.GetPtrFloat("muChi2NDF");
  Int_t* muStations     = data.GetPtrInt("muStations");
  //float* muInnerD0GV    = data.GetPtrFloat("muInnerD0GV");
  //float* muInnerDzGV    = data.GetPtrFloat("muInnerDzGV");
  float* muD0    = data.GetPtrFloat("muD0");
  float* muDz    = data.GetPtrFloat("muDz");
  
  //Int_t* muNumberOfValidMuonHits = data.GetPtrInt("muNumberOfValidMuonHits");
  Int_t* muNumberOfValidMuonHits = data.GetPtrInt("muMuonHits");
  //Int_t* muNumberOfValidPixelHits = data.GetPtrInt("muNumberOfValidPixelHits");
  Int_t* muNumberOfValidPixelHits = data.GetPtrInt("muPixelHits");
  //Int_t* muNumberOfValidTrkLayers = data.GetPtrInt("muNumberOfValidTrkLayers");
  Int_t* muNumberOfValidTrkLayers = data.GetPtrInt("muTrkLayers");

  float* muPFIsoR04_CH  = data.GetPtrFloat("muPFChIso");
  float* muPFIsoR04_NH  = data.GetPtrFloat("muPFNeuIso");
  float* muPFIsoR04_Pho = data.GetPtrFloat("muPFPhoIso");
  float* muPFIsoR04_PU  = data.GetPtrFloat("muPFPUIso");
  
  
  // for the energy scale correction
  Int_t  run      = data.GetInt("run");
  float* muPhi    = data.GetPtrFloat("muPhi");
  Int_t* muCharge = data.GetPtrInt("muCharge");
  //float* muInnerPt    = data.GetPtrFloat("muInnerPt");
  //float* muInnerPtErr = data.GetPtrFloat("muInnerPtErr");

  float* muInnerPt    = data.GetPtrFloat("muBestTrkPt");
  float* muInnerPtErr = data.GetPtrFloat("muBestTrkPtError");
  float *muPFMiniIso      = data.GetPtrFloat("muPFMiniIso"); 

  // loop over muon candidates
  for (int i = 0; i < nMu; ++i) {
    if (muPt[i] < 8) continue;

    // NOTE: applied at a primary vertex, not at the location of muon stations
    //also make sure that it lies in the eta region of interest
    if ( fabs(muEta[i])>1.4442 && fabs(muEta[i])<1.566 ) continue;

    if( fabs(muEta[i]) > 2.5 ) continue;
    
    //if (fabs(muEta[i]) > 2.4) continue;

    float isLoose = 0;

    float isTight = 0;

    float isTightIDonly = 0;
    float isTightIsoonly = 0;

    // cut-based loose ID for 2012
    //
    // muon types are defined here: CMSSW/DataFormats/MuonReco/interface/Muon.h
    // namely: GlobalMuon     = 1<<1
    //         TrackerMuon    = 1<<2
    //         StandAloneMuon = 1<<3
    //         CaloMuon       = 1<<4
    //         PFMuon         = 1<<5
    //         RPCMuon        = 1<<6
    if ((muType[i] & 1<<5) != 0 && // PFMuon && (GlobalMuon || TrackerMuon)
        ((muType[i] & 1<<1) != 0 || (muType[i] & 1<<2) != 0))
      isLoose = 1;


    
    // cut-based tight ID for 2012
    if ((muType[i] & 1<<1) != 0 && // GlobalMuon && PFMuon
        (muType[i] & 1<<5) != 0 &&
        muChi2NDF[i] < 10 &&
        muNumberOfValidMuonHits[i] > 0 &&
        muStations[i] > 1 &&
	fabs(muD0[i]) < 0.2 && // FIXME: GV is used historically,
        fabs(muDz[i]) < 0.5 && // inconsistent with vtx used for photons
        muNumberOfValidPixelHits[i] > 0 &&
        muNumberOfValidTrkLayers[i] > 5) isTightIDonly = 1;

    
    double reliso = (muPFIsoR04_CH[i] + TMath::Max(0.,muPFIsoR04_NH[i]+muPFIsoR04_Pho[i]-0.5*muPFIsoR04_PU[i]) )/muPt[i];

    //if(reliso<0.25 && isTightIDonly) isTight = 1; 

    double minireliso = muPFMiniIso[i]/muPt[i];

    if(reliso<0.25 && isTightIDonly) isTight = 1; 

    int isminiTight = 0;
    if(minireliso<0.2 && isTightIDonly) isminiTight = 1; 


    if (isLoose < 0.5 && isTight < 0.5 && isminiTight < 0.5)
      continue;

    // isolation for 2012
    /*if ((muPFIsoR04_CH[i] + max(0., muPFIsoR04_NH[i] + muPFIsoR04_Pho[i]
         - 0.5 * muPFIsoR04_PU[i])) / muPt[i] > 0.4)
      continue;
    */

    TLorentzVector mu;
    mu.SetPtEtaPhiM(muPt[i], muEta[i], muPhi[i], 0.105658);

    // instance of the muon energy corrector; random seed=229
    /*static rochcor2012* muEneCorr = new rochcor2012(229);

    // Correct the energy scale.
    //
    // NOTE: the corrector simply multiplies the Lorentz vector by a constant
    // factor, assuming that the muon candidate is massless
    float pt_error; // unused
    */

    
    /*if (data.HasMC())
      muEneCorr->momcor_mc(mu, muCharge[i], 0, pt_error);
    else {
      // NOTE: no separate treatment of 2012ABC vs 2012D runs is performed with
      // the newest rochcor2012jan22 corrector; the differentiating code is
      // kept, however, in order to do comparisons with the old rochcor2012v2
      // corrector
      if (run < 203768)
         muEneCorr->momcor_data(mu, muCharge[i], 0, pt_error); // 2012ABC
      else
         muEneCorr->momcor_data(mu, muCharge[i], 1, pt_error); // 2012D
    }
    */

    muPt[i] = mu.Pt();

    // for MC: correct muInnerPtErr to be in agreement with additional smearing
    /*if (data.HasMC()) {
      // smearing factors are taken from rochcor2012jan22.C
      double sf[8] = {0.0153729, 0.0103115, 0.00701322, 0.00472529, 0.00460413,
                      0.00700919, 0.00967325, 0.0147967};
      double err1 = muInnerPtErr[i]/muInnerPt[i];
      double err2 = sf[muEneCorr->etabin(muEta[i])];
      muInnerPtErr[i] = muInnerPt[i] * sqrt(err1*err1 + err2*err2);
    }
    */
    
    // fill results
    accepted.push_back(i);
    idLoose.push_back(isLoose);
    idTight.push_back(isTight);
    idminiTight.push_back(isminiTight);

  }
}
