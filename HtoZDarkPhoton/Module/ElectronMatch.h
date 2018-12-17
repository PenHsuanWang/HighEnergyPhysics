
int EleMatch(TreeReader &data, Float_t eleEta, Float_t elePhi){

  int matched;
  matched=0;

  Int_t    nMC = 0;
  Int_t*   mcPID = NULL;
  Int_t*   mcMomPID = NULL;
  Int_t*   mcGMomPID = NULL;
  Float_t* mcPt = NULL;
  Float_t* mcEta = NULL;
  Float_t* mcPhi = NULL;
  UShort_t*   mcStatusFlag = NULL;

  if (data.HasMC()) {
    nMC            = data.GetInt("nMC");
    mcPID          = data.GetPtrInt("mcPID");
    mcMomPID       = data.GetPtrInt("mcMomPID");
    mcGMomPID      = data.GetPtrInt("mcGMomPID");
    mcPt           = data.GetPtrFloat("mcPt");
    mcEta          = data.GetPtrFloat("mcEta");
    mcPhi          = data.GetPtrFloat("mcPhi");
    mcStatusFlag   = (UShort_t*) data.GetPtrShort("mcStatusFlag");
  }

  for (int i=0 ; i<nMC ; ++i){
    if (deltaR(eleEta, elePhi, mcEta[i], mcPhi[i]) < 0.1 ){
      if (abs(mcPID[i]) == 11){
	matched = 1; // matched
	break;
      }
    }
  }

  if (matched == 1) return 1;
  if (matched == 0) return 0;

} // end of function


int EleMatchZdSample(TreeReader &data, Float_t eleEta, Float_t elePhi){

  int matched;
  matched=0;

  Int_t    nMC = 0;
  Int_t*   mcPID = NULL;
  Int_t*   mcMomPID = NULL;
  Int_t*   mcGMomPID = NULL;
  Float_t* mcPt = NULL;
  Float_t* mcEta = NULL;
  Float_t* mcPhi = NULL;
  UShort_t*   mcStatusFlag = NULL;

  if (data.HasMC()) {
    nMC            = data.GetInt("nMC");
    mcPID          = data.GetPtrInt("mcPID");
    mcMomPID       = data.GetPtrInt("mcMomPID");
    mcGMomPID      = data.GetPtrInt("mcGMomPID");
    mcPt           = data.GetPtrFloat("mcPt");
    mcEta          = data.GetPtrFloat("mcEta");
    mcPhi          = data.GetPtrFloat("mcPhi");
    mcStatusFlag   = (UShort_t*) data.GetPtrShort("mcStatusFlag");
  }

  for (int i=0 ; i<nMC ; ++i){
    if (deltaR(eleEta, elePhi, mcEta[i], mcPhi[i]) < 0.1 ){
      if (abs(mcPID[i]) == 11 && abs(mcMomPID[i]) == 4900023 && abs(mcGMomPID[i]) == 25){
        matched = 1; // matched
        break;
      }
    }
  }

  if (matched == 1) return 1;
  if (matched == 0) return 0;

} // end of function


int EleMatchXZGSample(TreeReader &data, Float_t eleEta, Float_t elePhi){

  int matched;
  matched=0;

  Int_t    nMC = 0;
  Int_t*   mcPID = NULL;
  Int_t*   mcMomPID = NULL;
  Int_t*   mcGMomPID = NULL;
  Float_t* mcPt = NULL;
  Float_t* mcEta = NULL;
  Float_t* mcPhi = NULL;
  UShort_t*   mcStatusFlag = NULL;

  if (data.HasMC()) {
    nMC            = data.GetInt("nMC");
    mcPID          = data.GetPtrInt("mcPID");
    mcMomPID       = data.GetPtrInt("mcMomPID");
    mcGMomPID      = data.GetPtrInt("mcGMomPID");
    mcPt           = data.GetPtrFloat("mcPt");
    mcEta          = data.GetPtrFloat("mcEta");
    mcPhi          = data.GetPtrFloat("mcPhi");
    mcStatusFlag   = (UShort_t*) data.GetPtrShort("mcStatusFlag");
  }

  for (int i=0 ; i<nMC ; ++i){
    if (deltaR(eleEta, elePhi, mcEta[i], mcPhi[i]) < 0.1 ){
      if (abs(mcPID[i]) == 11 && abs(mcMomPID[i]) == 23 && abs(mcGMomPID[i]) == 25){
        matched = 1; // matched
        break;
      }
    }
  }

  if (matched == 1) return 1;
  if (matched == 0) return 0;

} // end of function


/*

  This function is using to find which particle are matched. if return is 0 meaning it doesn't match with any particle

 */
int WhatEleMatch(TreeReader &data, Float_t eleEta, Float_t elePhi){

  int matched;
  matched=0;

  Int_t    nMC = 0;
  Int_t*   mcPID = NULL;
  Int_t*   mcMomPID = NULL;
  Int_t*   mcGMomPID = NULL;
  Float_t* mcPt = NULL;
  Float_t* mcEta = NULL;
  Float_t* mcPhi = NULL;
  UShort_t*   mcStatusFlag = NULL;

  if (data.HasMC()) {
    nMC            = data.GetInt("nMC");
    mcPID          = data.GetPtrInt("mcPID");
    mcMomPID       = data.GetPtrInt("mcMomPID");
    mcGMomPID      = data.GetPtrInt("mcGMomPID");
    mcPt           = data.GetPtrFloat("mcPt");
    mcEta          = data.GetPtrFloat("mcEta");
    mcPhi          = data.GetPtrFloat("mcPhi");
    mcStatusFlag   = (UShort_t*) data.GetPtrShort("mcStatusFlag");
  }

  for (int i=0 ; i<nMC ; ++i){
    if (deltaR(eleEta, elePhi, mcEta[i], mcPhi[i]) < 0.1 ){
      if (abs(mcPID[i]) == 11){
        matched = 1; // matched
        break;
      }
    }
  }

  if (matched == 1) return 11; 

  if (matched == 0){
    Int_t WhichParticle;
    WhichParticle=0;
    for (int i=0 ; i<nMC ; ++i){
      if (deltaR(eleEta, elePhi, mcEta[i], mcPhi[i]) < 0.1 ){
	WhichParticle = abs(mcPID[i]);
	break;
      }
    }
    return WhichParticle;
  }

} // end of function
