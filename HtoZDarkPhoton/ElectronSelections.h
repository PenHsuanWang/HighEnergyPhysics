#include <TRandom3.h>
#include "LeptonScaleCorrections.hh"

using namespace std;

//______________________________________________________________________________
float ElectronIDCutBased2012(int iWP, TreeReader &data, int i)
{
  /* Standard cut-based loose/medium/tight electron (and positron)
   * identification and isolation for 2012.
   *
   * Returns 1 if candidate is accepted or 0 otherwise.
   *
   * Documentation on the cut-based electron ID:
   * https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification?rev=30
   *
   * Documentation on the electron isolation:
   * https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaEARhoCorrection?rev=14
   *
   * iWP = working point to take: 0/loose, 1/medium, 2/tight;
   * data = handle providing access to an input event;
   * i = index of an electron/positron candidate to consider.
   */

  // load necessary tree branches
  float* elePt            = data.GetPtrFloat("elePt");
  float* eleSCEta         = data.GetPtrFloat("eleSCEta");
  float* eledEtaAtVtx     = data.GetPtrFloat("eledEtaAtVtx");
  float* eledPhiAtVtx     = data.GetPtrFloat("eledPhiAtVtx");
  float* eleSigmaIEtaIEta = data.GetPtrFloat("eleSigmaIEtaIEta");
  float* eleHoverE        = data.GetPtrFloat("eleHoverE");
  float* eleEcalEn        = data.GetPtrFloat("eleEcalEn");
  float* elePin           = data.GetPtrFloat("elePin");
  Int_t* eleConvVtxFit    = data.GetPtrInt("eleConvVtxFit");
  Int_t* eleMissHits      = data.GetPtrInt("eleMissHits");
  float* elePFChIso04     = data.GetPtrFloat("elePFChIso04");
  float* elePFPhoIso04    = data.GetPtrFloat("elePFPhoIso04");
  float* elePFNeuIso04    = data.GetPtrFloat("elePFNeuIso04");
  float  rho2012          = data.GetFloat("rho2012");
  vector<float>* eleD0Vtx = data.GetPtrVectorFloat("eleD0Vtx");
  vector<float>* eleDzVtx = data.GetPtrVectorFloat("eleDzVtx");

  // to improve the code readability
  float absEta = fabs(eleSCEta[i]);

  // NOTE: no |eleSCeta| < 2.5 and 1.4442 < |eleSCeta| < 1.566 cuts

  // 0=ECAL barrel or 1=ECAL endcaps
  int iBE = (absEta < 1.479) ? 0 : 1;

  // NOTE: the two-dimensional arrays below are [loose,medium,tight][EB,EE]

  // dEtaIn
  float dEtaInCut[3][2] = {{0.007, 0.009}, {0.004, 0.007}, {0.004, 0.005}};
  if (fabs(eledEtaAtVtx[i]) >= dEtaInCut[iWP][iBE]) return 0;

  // dPhiIn
  float dPhiInCut[3][2] = {{0.15, 0.10}, {0.06, 0.03}, {0.03, 0.02}};
  if (fabs(eledPhiAtVtx[i]) >= dPhiInCut[iWP][iBE]) return 0;

  // sigmaIEtaIEta
  float sigmaIEIECut[2] = {0.01, 0.03}; // [EB,EE]
  if (eleSigmaIEtaIEta[i] >= sigmaIEIECut[iBE]) return 0;

  // HCAL/ECAL
  float hoverECut[2] = {0.12, 0.10}; // [EB,EE]
  if (eleHoverE[i] >= hoverECut[iBE]) return 0;

  // |d0(vtx)|, vtx = primary (number zero) vertex from beam spot
  if (fabs(eleD0Vtx[i][0]) >= 0.02) return 0;

  // |dZ(vtx)|, vtx = primary (number zero) vertex from beam spot
  float dZCut[3] = {0.2, 0.1, 0.1}; // [loose,medium,tight]
  if (fabs(eleDzVtx[i][0]) >= dZCut[iWP]) return 0;

  // |1/E - 1/p|; uncorrected ECAL energy is used
  if (fabs(1/eleEcalEn[i] - 1/elePin[i]) >= 0.05) return 0;

  // if has matched conversion
  if (eleConvVtxFit[i] == 1) return 0;

  // missing hits
  Int_t missHitsCut[3] = {1, 1, 0};
  if (eleMissHits[i] > missHitsCut[iWP]) return 0;

  // find eta bin which determines effective area for the rho correction to the
  // isolation
  int bin;
  if (absEta < 1.0) bin = 0;
  else if (absEta >= 1.0 && absEta < 1.479) bin = 1;
  else if (absEta >= 1.479 && absEta < 2.0) bin = 2;
  else if (absEta >= 2.0 && absEta < 2.2) bin = 3;
  else if (absEta >= 2.2 && absEta < 2.3) bin = 4;
  else if (absEta >= 2.3 && absEta < 2.4) bin = 5;
  else bin = 6; // absEta >= 2.4

  // rho-corrected PF isolation
  double effArea[] = {0.208, 0.209, 0.115, 0.143, 0.183, 0.194, 0.261};
  double pfIsoCut[3] = {0.15, 0.15, 0.10}; // [loose,medium,tight]
  double neu = elePFPhoIso04[i] + elePFNeuIso04[i] - effArea[bin] * rho2012;
  if ((elePFChIso04[i] + max(neu, 0.)) / elePt[i] >= pfIsoCut[iWP]) return 0;

  // passed all ID and isolation cuts
  return 1;
}

//______________________________________________________________________________
float ElectronIDMVA(TreeReader &data, int i)
{
  /* Preselection, identification with MVA and isolation for electron (and
   * positron) candidates.
   *
   * Returns 1 if candidate is accepted or 0 otherwise.
   *
   * The preselection is necessary for matching with the online trigger
   * conditions.
   *
   * Documentation on the electron ID with MVA:
   * https://twiki.cern.ch/twiki/bin/viewauth/CMS/MultivariateElectronIdentification?rev=53
   *
   * Documentation on the electron isolation:
   * https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaEARhoCorrection?rev=14
   *
   * data = handle providing access to an input event;
   * i = index of an electron/positron candidate to consider.
   */

  // load necessary tree branches
  Int_t* eleMissHits      = data.GetPtrInt("eleMissHits");
  float* elePt            = data.GetPtrFloat("elePt");
  float* eleSCEta         = data.GetPtrFloat("eleSCEta");
  float* eleSigmaIEtaIEta = data.GetPtrFloat("eleSigmaIEtaIEta");
  float* eleHoverE        = data.GetPtrFloat("eleHoverE");
  float* eleIsoEcalDR03   = data.GetPtrFloat("eleIsoEcalDR03");
  float* eleIsoHcalDR03   = data.GetPtrFloat("eleIsoHcalDR03");
  float* eleIsoTrkDR03    = data.GetPtrFloat("eleIsoTrkDR03");
  float* eleIDMVATrig     = data.GetPtrFloat("eleIDMVATrig");
  float* elePFChIso04     = data.GetPtrFloat("elePFChIso04");
  float* elePFPhoIso04    = data.GetPtrFloat("elePFPhoIso04");
  float* elePFNeuIso04    = data.GetPtrFloat("elePFNeuIso04");
  float  rho2012          = data.GetFloat("rho2012");

  if (eleMissHits[i] != 0) return 0;

  // detector-driven pre-isolation
  if (eleIsoTrkDR03[i]/elePt[i] >= 0.2) return 0;
  if (eleIsoEcalDR03[i]/elePt[i] >= 0.2) return 0;
  if (eleIsoHcalDR03[i]/elePt[i] >= 0.2) return 0;

  // to improve the code readability
  float absEta = fabs(eleSCEta[i]);

  // ECAL barrel vs ECAL endcaps
  if (absEta < 1.479) {
    if (eleSigmaIEtaIEta[i] >= 0.014) return 0;
    if (eleHoverE[i] >= 0.15) return 0;
  } else {
    if (eleSigmaIEtaIEta[i] >= 0.035) return 0;
    if (eleHoverE[i] >= 0.10) return 0;
  }

  // MVA ID
  if (elePt[i] <= 20 && eleIDMVATrig[i] < -0.9) return 0;
  if (elePt[i] > 20 && eleIDMVATrig[i] < -0.5) return 0;

  // find eta bin which determines effective area for the rho correction to the
  // isolation
  int bin;
  if (absEta < 1.0) bin = 0;
  else if (absEta >= 1.0 && absEta < 1.479) bin = 1;
  else if (absEta >= 1.479 && absEta < 2.0) bin = 2;
  else if (absEta >= 2.0 && absEta < 2.2) bin = 3;
  else if (absEta >= 2.2 && absEta < 2.3) bin = 4;
  else if (absEta >= 2.3 && absEta < 2.4) bin = 5;
  else bin = 6; // absEta >= 2.4

  // rho-corrected PF isolation
  double effArea[] = {0.208, 0.209, 0.115, 0.143, 0.183, 0.194, 0.261};
  double neu = elePFPhoIso04[i] + elePFNeuIso04[i] - effArea[bin] * rho2012;
  if ((elePFChIso04[i] + max(neu, 0.)) / elePt[i] > 0.4) return 0;

  // passed MVA ID and isolation cuts
  return 1;
}

//______________________________________________________________________________
void select_electrons(TreeReader &data, vector<int> &accepted,
                      vector<float> &idLoose, vector<float> &idMedium,
                      vector<float> &idTight, vector<float> &idMVATrig)
{
  /* Identification, isolation and energy scale correction for electrons (and
   * positrons) for 2012.
   *
   * NOTE: the energy scale correction is applied to the elePt branch only.
   * Other energy-related branches (e.g eleEn) are not corrected!
   *
   * data = handle providing access to an input event;
   * accepted = array to fill with indices of accepted electrons/positrons;
   *
   * idLoose, idMedium, idTight = arrays to fill with results of cut-based
   *   loose, medium and tight identification criteria (0=not passed, 1=passed);
   *
   * idMVATrig = array to fill with results of triggering MVA identification.
   */

  // load necessary tree branches
  Int_t  nEle       = data.GetInt("nEle");
  float* elePt      = data.GetPtrFloat("elePt");

  // for the energy scale correction
  Int_t  run        = data.GetInt("run");
  float* eleSCEta   = data.GetPtrFloat("eleSCEta");

  ///SJ
  //float* eleRegrE   = data.GetPtrFloat("eleRegrE");

  float* eleE3x3    = data.GetPtrFloat("eleE3x3");
  float* eleSCRawEn = data.GetPtrFloat("eleSCRawEn");
  float* eleEta     = data.GetPtrFloat("eleEta");
  float* eleEn      = data.GetPtrFloat("eleEn");

  // loop over electron/positron candidates
  for (int i = 0; i < nEle; ++i) {
    if (elePt[i] < 8) continue;

    float isMVATrig = ElectronIDMVA(data, i);
    float isLoose = ElectronIDCutBased2012(0, data, i);
    float isMedium = ElectronIDCutBased2012(1, data, i);
    float isTight = ElectronIDCutBased2012(2, data, i);

    // skip candidate if none of the IDs accepted it
    if (isLoose < 0.5 && isMedium < 0.5 && isTight < 0.5 && isMVATrig < 0.5)
       continue;

    // correct the energy scale
    // NOTE: eleE3x3[i]/eleSCRawEn[i] somehow differs from eleR9[i]
    static TRandom3* rand = new TRandom3(230);
    
    ///SJ
    /*if (eleRegrE[i] != 0) {
      double eneCorr = correctedElectronEnergy(eleRegrE[i], eleSCEta[i],
         eleE3x3[i]/eleSCRawEn[i], run, 1, "Moriond2013", data.HasMC(), rand);
      elePt[i] = eneCorr/cosh(eleEta[i]);
    }
    else {
    double eneCorr = correctedElectronEnergy(eleEn[i], eleSCEta[i],
    eleE3x3[i]/eleSCRawEn[i], run, 0, "HCP2012", data.HasMC(), rand);
    elePt[i] = eneCorr/cosh(eleEta[i]);
    }
    */

    double eneCorr = correctedElectronEnergy(eleEn[i], eleSCEta[i],
					     eleE3x3[i]/eleSCRawEn[i], run, 0, "HCP2012", data.HasMC(), rand);
    elePt[i] = eneCorr/cosh(eleEta[i]);

    // fill results
    accepted.push_back(i);
    idLoose.push_back(isLoose);
    idMedium.push_back(isMedium);
    idTight.push_back(isTight);
    idMVATrig.push_back(isMVATrig);
  }
}




///// For run 2
float ElectronIDCutBasedRunII(int iWP, TreeReader &data, int i)
{
  /* Standard cut-based loose/medium/tight electron (and positron)
   * identification and isolation for 2012.
   *
   * Returns 1 if candidate is accepted or 0 otherwise.
   *
   * Documentation on the cut-based electron ID:
   * https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification?rev=30
   *
   * Documentation on the electron isolation:
   * https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaEARhoCorrection?rev=14
   *
   **** NOT THIS* iWP = working point to take: 0/loose, 1/medium, 2/tight;***
   **** THIS* iWP = working point to take: 0/veto, 1/loose, 2/medium, 3/tight;***
   * data = handle providing access to an input event;
   * i = index of an electron/positron candidate to consider.
   */


  
  
  
  // load necessary tree branches
  float* elePt            = data.GetPtrFloat("eleCalibPt");
  float* eleSCEta         = data.GetPtrFloat("eleSCEta");
  float* eledEtaAtVtx     = data.GetPtrFloat("eledEtaAtVtx");
  float* eledPhiAtVtx     = data.GetPtrFloat("eledPhiAtVtx");
  float* eleSigmaIEtaIEtaFull5x5 = data.GetPtrFloat("eleSigmaIEtaIEtaFull5x5");

  float* eleHoverE        = data.GetPtrFloat("eleHoverE");
  //float* eleRelIsoWithDBeta = data.GetPtrFloat("eleRelIsoWithDBeta");  

  float* eleD0        = data.GetPtrFloat("eleD0");
  float* eleDz        = data.GetPtrFloat("eleDz");
  float* eleEoverPInv        = data.GetPtrFloat("eleEoverPInv");

  Int_t* eleConvVeto    = data.GetPtrInt("eleConvVeto");
  Int_t* eleMissHits      = data.GetPtrInt("eleMissHits");

  float *elePFMiniIso      = data.GetPtrFloat("elePFMiniIso");

  // to improve the code readability
  float absEta = fabs(eleSCEta[i]);
  
  // NOTE: no |eleSCeta| < 2.5 and 1.4442 < |eleSCeta| < 1.566 cuts

  // 0=ECAL barrel or 1=ECAL endcaps
  int iBE = (absEta < 1.479) ? 0 : 1;

  // NOTE: the two-dimensional arrays below are [loose,medium,tight][EB,EE]

  // dEtaIn
  float dEtaInCut[4][2] = {{0.016315, 0.010671}, {0.012442, 0.010654}, {0.007641, 0.009285}, {0.006574,0.005681}};
  if (fabs(eledEtaAtVtx[i]) >= dEtaInCut[iWP][iBE]) return 0;

  // dPhiIn
  float dPhiInCut[4][2] = {{0.252044, 0.245263}, {0.072624, 0.145129}, {0.032643, 0.042447}, {0.022868, 0.032046}};
  if (fabs(eledPhiAtVtx[i]) >= dPhiInCut[iWP][iBE]) return 0;

  // sigmaIEtaIEta
  float sigmaIEIECut[4][2] = {{0.011100, 0.033987}, {0.010557, 0.032602}, {0.010399, 0.029524}, {0.010181, 0.028766}}; // [EB,EE]
  if (eleSigmaIEtaIEtaFull5x5[i] >= sigmaIEIECut[iWP][iBE]) return 0;

  // HCAL/ECAL
  float hoverECut[4][2] = {{0.345843, 0.134691}, {0.121476,0.131862}, {0.060662, 0.104263}, {0.037553, 0.081902}}; // [EB,EE]
  if (eleHoverE[i] >= hoverECut[iWP][iBE]) return 0;

  // |d0(vtx)|, vtx = primary (number zero) vertex from beam spot
  float d0Cut[4][2] = {{0.060279, 0.273097}, {0.022664, 0.097358}, {0.011811, 0.051682}, {0.009924, 0.027261}};
  if (fabs(eleD0[i]) >= d0Cut[iWP][iBE]) return 0;

  // |dZ(vtx)|, vtx = primary (number zero) vertex from beam spot
  float dZCut[4][2] = {{0.800538, 0.885860}, {0.173670, 0.198444}, {0.070775, 0.180720}, {0.15310, 0.147154}};
  if (fabs(eleDz[i]) >= dZCut[iWP][iBE]) return 0;

  // |1/E - 1/p|; uncorrected ECAL energy is used
  float ooEoopCut[4][2] = {{0.248070, 0.157160}, {0.221803, 0.142283}, {0.153987, 0.137468}, {0.131191, 0.106055}};
  if (fabs(eleEoverPInv[i]) >= ooEoopCut[iWP][iBE]) return 0;

  //PF isolation with beta
  //float relisoCut[4][2] = {{0.164369, 0.212604}, {0.120026, 0.162914}, {0.097213, 0.116708}, {0.074355, 0.090185}};
  //if(eleRelIsoWithDBeta[i] >= relisoCut[iWP][iBE]) return 0;

  // mini isolation
  if( elePFMiniIso[i]/elePt[i] > 0.1 ) return 0;

  // missing hits
  Int_t missHitsCut[4][2] = {{2, 3}, {1, 1}, {1,1}, {1,1}};
  if (eleMissHits[i] > missHitsCut[iWP][iBE]) return 0;

  
  if (!eleConvVeto[i])  return 0;
  
  
  // passed all ID and isolation cuts
  return 1;
}


////HEEP ID
float ElectronHEEPIDRunII(TreeReader &data, int i)
{
  /* Standard HEEPID loose/medium/tight electron (and positron)
   * identification and isolation for 2012.
   *
   * Returns 1 if candidate is accepted or 0 otherwise.
   *
   * Documentation on the cut-based electron ID:
   * https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification?rev=30
   *
   * Documentation on the electron isolation:
   * https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaEARhoCorrection?rev=14
   *
   **** NOT THIS* iWP = working point to take: 0/loose, 1/medium, 2/tight;***
   **** THIS* iWP = working point to take: 0/veto, 1/loose, 2/medium, 3/tight;***
   * data = handle providing access to an input event;
   * i = index of an electron/positron candidate to consider.
   */

  // load necessary tree branches
  //float* elePt            = data.GetPtrFloat("eleCalibPt");

  //bool debug = true;
  bool debug = false;
  
  float* elecaloEnergy    = data.GetPtrFloat("elecaloEnergy");
  float* eleSCEta         = data.GetPtrFloat("eleSCEta");

  //float* eleTheta         = data.GetPtrFloat("eleTheta");

  float* eleSCEn          = data.GetPtrFloat("eleSCEn");
  float* eledEtaseedAtVtx     = data.GetPtrFloat("eledEtaseedAtVtx");
  float* eledPhiAtVtx     = data.GetPtrFloat("eledPhiAtVtx");
  float* eleSigmaIEtaIEtaFull5x5 = data.GetPtrFloat("eleSigmaIEtaIEtaFull5x5");

  float* eleHoverE        = data.GetPtrFloat("eleHoverE");

  float* eleTrkdxy        = data.GetPtrFloat("eleTrkdxy");
  float* eleE1x5Full5x5     = data.GetPtrFloat("eleE1x5Full5x5");
  float* eleE2x5Full5x5     = data.GetPtrFloat("eleE2x5Full5x5");
  float* eleE5x5Full5x5     = data.GetPtrFloat("eleE5x5Full5x5");
  
  float* eleDr03EcalRecHitSumEt = data.GetPtrFloat("eleDr03EcalRecHitSumEt");
  float* eleDr03HcalDepth1TowerSumEt = data.GetPtrFloat("eleDr03HcalDepth1TowerSumEt");
  float* eleDr03TkSumPt = data.GetPtrFloat("eleDr03TkSumPt");

  float  rho2012          = data.GetFloat("rho");

  Int_t* eleMissHits      = data.GetPtrInt("eleMissHits");
  Int_t* eleEcalDrivenSeed= data.GetPtrInt("eleEcalDrivenSeed");
  
  

  // to improve the code readability
  float absEta = fabs(eleSCEta[i]);
  //float pt     = elecaloEnergy[i]*sin(eleTheta[i]);
  float pt     = elecaloEnergy[i]/cosh(eleSCEta[i]);
  // NOTE: no |eleSCeta| < 2.5 and 1.4442 < |eleSCeta| < 1.566 cuts

  // 0=ECAL barrel or 1=ECAL endcaps
  int iBE = (absEta < 1.4442) ? 0 : 1;

  // NOTE: the two-dimensional arrays below are [loose,medium,tight][EB,EE]

  // isEcalDriven?
  if(!eleEcalDrivenSeed[i]) return 0;

  // dEtaIn
  float dEtaInCut[2] = {0.004, 0.006};
  if (fabs(eledEtaseedAtVtx[i]) >= dEtaInCut[iBE]) return 0;

  // dPhiIn
  float dPhiInCut[2] = {0.06, 0.06};
  if (fabs(eledPhiAtVtx[i]) >= dPhiInCut[iBE]) return 0;

  // HCAL/ECAL
  float hoverECut[2] = { (float)(2./eleSCEn[i]+0.05), (float)(12.5/eleSCEn[i]+0.05)};
  if (eleHoverE[i] >= hoverECut[iBE]) return 0;

  // sigmaIEtaIEta
  float sigmaIEIECut[2] = {99999, 0.03};
  if (eleSigmaIEtaIEtaFull5x5[i] >= sigmaIEIECut[iBE]) return 0;

  // shower shape
  float e1x5Cut[2] = {0.83, -1.0};
  float e2x5Cut[2] = {0.94, -1.0};
  if(eleE1x5Full5x5[i]<e1x5Cut[iBE]  && eleE2x5Full5x5[i]/eleE5x5Full5x5[i]<e2x5Cut[iBE]) return 0;
    

  ///isolation
  //ECAL+HCAL
  float ehcalCut[2] = {(float)(2.0+0.03*pt+0.28*rho2012)
		       , (pt<50) ? (float)(2.5+0.28*rho2012) : (float)(2.5+0.03*(pt-50) + 0.28*rho2012) };


  if(debug) cout<<"isEE? pt : rho : ehcal cut : "<<iBE<<":"<<pt<<":"<<rho2012<<":"<<ehcalCut[iBE]<<endl;

  if( eleDr03EcalRecHitSumEt[i]+eleDr03HcalDepth1TowerSumEt[i] > ehcalCut[iBE] ) return 0;

  //track
  float trkCut[2] = {5, 5};
  if(eleDr03TkSumPt[iBE] > trkCut[iBE]) return 0;

  // missing hits
  Int_t missHitsCut[2] = {1, 1};
  if (eleMissHits[i] > missHitsCut[iBE]) return 0;

  //dxy cut
  float dxyCut[2] = {0.02, 0.05};
  if (fabs(eleTrkdxy[i]) >= dxyCut[iBE]) return 0;

  
  // passed all ID and isolation cuts
  return 1;
}




//_____________________________________MVA______________________________________
float ElectronIDMVARunII(TreeReader &data, int i)
{
  /* Preselection, identification with MVA and isolation for electron (and
   * positron) candidates.
   *
   * Returns 1 if candidate is accepted or 0 otherwise.
   *
   * The preselection is necessary for matching with the online trigger
   * conditions.
   *
   * Documentation on the electron ID with MVA:
   * https://twiki.cern.ch/twiki/bin/viewauth/CMS/MultivariateElectronIdentification?rev=53
   *
   * Documentation on the electron isolation:
   * https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaEARhoCorrection?rev=14
   *
   * data = handle providing access to an input event;
   * i = index of an electron/positron candidate to consider.
   */

  // load necessary tree branches
  Int_t* eleMissHits      = data.GetPtrInt("eleMissHits");
  float* elePt            = data.GetPtrFloat("eleCalibPt");
  float* eleSCEta         = data.GetPtrFloat("eleSCEta");
  float* eleSigmaIEtaIEta = data.GetPtrFloat("eleSigmaIEtaIEta");
  float* eleHoverE        = data.GetPtrFloat("eleHoverE");
  float* eleIsoEcalDR03   = data.GetPtrFloat("eleIsoEcalDR03");
  float* eleIsoHcalDR03   = data.GetPtrFloat("eleIsoHcalDR03");
  float* eleIsoTrkDR03    = data.GetPtrFloat("eleIsoTrkDR03");
  float* eleIDMVATrig     = data.GetPtrFloat("eleIDMVATrig");
  float* elePFChIso04     = data.GetPtrFloat("elePFChIso04");
  float* elePFPhoIso04    = data.GetPtrFloat("elePFPhoIso04");
  float* elePFNeuIso04    = data.GetPtrFloat("elePFNeuIso04");
  float  rho2012          = data.GetFloat("rho");

  if (eleMissHits[i] != 0) return 0;

  // detector-driven pre-isolation
  if (eleIsoTrkDR03[i]/elePt[i] >= 0.2) return 0;
  if (eleIsoEcalDR03[i]/elePt[i] >= 0.2) return 0;
  if (eleIsoHcalDR03[i]/elePt[i] >= 0.2) return 0;

  // to improve the code readability
  float absEta = fabs(eleSCEta[i]);

  // ECAL barrel vs ECAL endcaps
  if (absEta < 1.479) {
    if (eleSigmaIEtaIEta[i] >= 0.014) return 0;
    if (eleHoverE[i] >= 0.15) return 0;
  } else {
    if (eleSigmaIEtaIEta[i] >= 0.035) return 0;
    if (eleHoverE[i] >= 0.10) return 0;
  }

  // MVA ID
  if (elePt[i] <= 20 && eleIDMVATrig[i] < -0.9) return 0;
  if (elePt[i] > 20 && eleIDMVATrig[i] < -0.5) return 0;

  // find eta bin which determines effective area for the rho correction to the
  // isolation
  int bin;
  if (absEta < 1.0) bin = 0;
  else if (absEta >= 1.0 && absEta < 1.479) bin = 1;
  else if (absEta >= 1.479 && absEta < 2.0) bin = 2;
  else if (absEta >= 2.0 && absEta < 2.2) bin = 3;
  else if (absEta >= 2.2 && absEta < 2.3) bin = 4;
  else if (absEta >= 2.3 && absEta < 2.4) bin = 5;
  else bin = 6; // absEta >= 2.4

  // rho-corrected PF isolation
  double effArea[] = {0.208, 0.209, 0.115, 0.143, 0.183, 0.194, 0.261};
  double neu = elePFPhoIso04[i] + elePFNeuIso04[i] - effArea[bin] * rho2012;
  if ((elePFChIso04[i] + max(neu, 0.)) / elePt[i] > 0.4) return 0;

  // passed MVA ID and isolation cuts
  return 1;
}



//______________________________________________________________________________
void select_electronsRunII(TreeReader &data, vector<int> &accepted,
                      vector<float> &idLoose, vector<float> &idMedium,
                      vector<float> &idTight, vector<float> &idMVATrig)
{
  /* Identification, isolation and energy scale correction for electrons (and
   * positrons) for 2012.
   *
   * NOTE: the energy scale correction is applied to the elePt branch only.
   * Other energy-related branches (e.g eleEn) are not corrected!
   *
   * data = handle providing access to an input event;
   * accepted = array to fill with indices of accepted electrons/positrons;
   *
   * idLoose, idMedium, idTight = arrays to fill with results of cut-based
   *   loose, medium and tight identification criteria (0=not passed, 1=passed);
   *
   * idMVATrig = array to fill with results of triggering MVA identification.
   */

  // load necessary tree branches
  Int_t  nEle       = data.GetInt("nEle");
  float* elePt      = data.GetPtrFloat("eleCalibPt");

  // for the energy scale correction
  Int_t  run        = data.GetInt("run");
  float* eleSCEta   = data.GetPtrFloat("eleSCEta");

  ///SJ
  //float* eleRegrE   = data.GetPtrFloat("eleRegrE");

  //float* eleE3x3    = data.GetPtrFloat("eleE3x3");
  float* eleSCRawEn = data.GetPtrFloat("eleSCRawEn");
  float* eleEta     = data.GetPtrFloat("eleEta");
  float* eleEn      = data.GetPtrFloat("eleCalibEn");
  Short_t* eleid = data.GetPtrShort("eleIDbit");
  float*  eleIDMVATrg = data.GetPtrFloat("eleIDMVATrg");
  // loop over electron/positron candidates



  vector<int> accepted_uns;
  vector<float> idLoose_uns, idMedium_uns, idTight_uns, idMVATrig_uns;
  vector<float> elePt_uns; // array to be sorted

  for (int i = 0; i < nEle; ++i) {
    if (elePt[i] < 8) continue;
    
    //also make sure that it lies in the eta region of interest
    //if ( fabs(eleSCEta[i])>1.4442 && fabs(eleSCEta[i])<1.566 ) continue;

    if( fabs(eleSCEta[i]) > 2.5 ) continue;

    //float isMVATrig = ElectronIDMVARunII(data, i);
    int isMVATrig = 1;
    /*float isVeto    = ElectronIDCutBasedRunII(0, data, i);
    float isLoose   = ElectronIDCutBasedRunII(1, data, i);
    float isMedium  = ElectronIDCutBasedRunII(2, data, i);
    float isTight   = ElectronIDCutBasedRunII(3, data, i);
    float isHEEP    = ElectronHEEPIDRunII(data,i);
    */


    //Long64_t* eleid = data.GetPtrLong64("eleIDbit");
    float isVeto    = eleid[i]>>0&1; 
    float isLoose   = eleid[i]>>1&1; 
    float isMedium  = eleid[i]>>2&1; 
    float isTight   = eleid[i]>>3&1; 
    float isHEEP    = eleid[i]>>4&1; 
    
    double mvacut = 0.972153;
    if(fabs(eleSCEta[i]) > 0.8 && fabs(eleSCEta[i]) <1.566) mvacut = 0.922126;

    if(fabs(eleSCEta[i]) > 1.566 && fabs(eleSCEta[i]) <2.5) mvacut = 0.610764;

    isMVATrig = eleIDMVATrg[i] > mvacut;
    
    // skip candidate if none of the IDs accepted it
    if (isVeto < 0.5 && isLoose < 0.5 && isMedium < 0.5 && isTight < 0.5 && !isMVATrig)
       continue;

    // correct the energy scale
    // NOTE: eleE3x3[i]/eleSCRawEn[i] somehow differs from eleR9[i]
    static TRandom3* rand = new TRandom3(230);
    


    // fill results
    accepted_uns.push_back(i);
    idLoose_uns.push_back(isLoose);
    idMedium_uns.push_back(isMedium);
    idTight_uns.push_back(isTight);
    idMVATrig_uns.push_back(isMVATrig);
    elePt_uns.push_back(elePt[i]);

    ///SJ
    /*if (eleRegrE[i] != 0) {
      double eneCorr = correctedElectronEnergy(eleRegrE[i], eleSCEta[i],
         eleE3x3[i]/eleSCRawEn[i], run, 1, "Moriond2013", data.HasMC(), rand);
      elePt[i] = eneCorr/cosh(eleEta[i]);
    }
    else {
    double eneCorr = correctedElectronEnergy(eleEn[i], eleSCEta[i],
    eleE3x3[i]/eleSCRawEn[i], run, 0, "HCP2012", data.HasMC(), rand);
    elePt[i] = eneCorr/cosh(eleEta[i]);
    }
    

    double eneCorr = correctedElectronEnergy(eleEn[i], eleSCEta[i],
					     eleE3x3[i]/eleSCRawEn[i], run, 0, "HCP2012", data.HasMC(), rand);
    elePt[i] = eneCorr/cosh(eleEta[i]);
    */

    //SJ for now
    //elePt[i] = eleSCRawEn[i]/cosh(eleEta[i]);
  }



  int siz = (int) accepted_uns.size();
  if (siz < 1) return;
  
  // sort accepted photons in descending order of their pt.
  int ind[siz];
  TMath::Sort(siz, &elePt_uns.front(), ind);

  // fill to-be-filled arrays with sorted results
  for (int i = 0; i < siz; ++i) {

    // fill results
    accepted.push_back(accepted_uns[ind[i]]);
    idLoose.push_back(idLoose_uns[ind[i]]);
    idMedium.push_back(idMedium_uns[ind[i]]);
    idTight.push_back(idTight_uns[ind[i]]);
    idMVATrig.push_back(idMVATrig_uns[ind[i]]);
  }
}


///preselection for triggering electrons
float ElectronIDPreselRunII(TreeReader &data, int i)
{
  
  // load necessary tree branches
  float* elePt            = data.GetPtrFloat("eleCalibPt");
  float* eleSCEta         = data.GetPtrFloat("eleSCEta");
  float* eledEtaAtVtx     = data.GetPtrFloat("eledEtaAtVtx");
  float* eledPhiAtVtx     = data.GetPtrFloat("eledPhiAtVtx");
  float* eleSigmaIEtaIEtaFull5x5 = data.GetPtrFloat("eleSigmaIEtaIEtaFull5x5");

  float* eleHoverE        = data.GetPtrFloat("eleHoverE");
  float* elePFClusEcalIso = data.GetPtrFloat("elePFClusEcalIso");  
  float* elePFClusHcalIso = data.GetPtrFloat("elePFClusHcalIso");  
  float* eleDr03TkSumPt   = data.GetPtrFloat("eleDr03TkSumPt");  
  float *elePFMiniIso      = data.GetPtrFloat("elePFMiniIso"); 
  

  // to improve the code readability
  float absEta = fabs(eleSCEta[i]);
  
  // NOTE: no |eleSCeta| < 2.5 and 1.4442 < |eleSCeta| < 1.566 cuts

  // 0=ECAL barrel or 1=ECAL endcaps
  int iBE = (absEta < 1.479) ? 0 : 1;

  // NOTE: the two-dimensional arrays below are [loose,medium,tight][EB,EE]

  //pT 
  float ptCut[2] = {15., 15.};
  //float ptCut[2] = {5., 5.}; // for Zd sample, lower pT cut is needed # added by PH in 18/08/16
  if( elePt[i] < ptCut[iBE] ) return 0;
  
  // dEtaIn
  float dEtaInCut[2] = {0.0095,999};
  if (fabs(eledEtaAtVtx[i]) >= dEtaInCut[iBE]) return 0;

  // dPhiIn
  float dPhiInCut[2] = {0.065, 999};
  if (fabs(eledPhiAtVtx[i]) >= dPhiInCut[iBE]) return 0;

  // sigmaIEtaIEta
  float sigmaIEIECut[2] = {0.012, 0.033};
  if (eleSigmaIEtaIEtaFull5x5[i] >= sigmaIEIECut[iBE]) return 0;

  // HCAL/ECAL
  float hoverECut[2] = {0.09, 0.09};
  if (eleHoverE[i] >= hoverECut[iBE]) return 0;

  /*
  //ECAL PF isolation with 
  float pfEcalIsoCutOverPt[2] = {0.37, 0.45};
  if(elePFClusEcalIso[i]/elePt[i] >= pfEcalIsoCutOverPt[iBE]) return 0;

  //HCAL PF isolation with 
  float pfHcalIsoCutOverPt[2] = {0.25, 0.28};
  if(elePFClusHcalIso[i]/elePt[i] >= pfHcalIsoCutOverPt[iBE]) return 0;
  
  //tracker isolation with 
  float trk03IsoCutOverPt[2] = {0.18, 0.18};
  if(eleDr03TkSumPt[i]/elePt[i] >= trk03IsoCutOverPt[iBE]) return 0;
  */


  ///just put a cut on miniIsolation


  if( elePFMiniIso[i]/elePt[i] > 0.1 ) return 0;

  // passed all ID and isolation cuts
  return 1;
}


///COMMENTS: Write ElectronIDMVARunII
///Complete HEEP ID
///COmment out energy correction for now
