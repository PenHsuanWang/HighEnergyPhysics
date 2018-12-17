
int PhoMatch(TreeReader &data, Float_t phoEta, Float_t phoPhi){

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
    //if (((mcStatusFlag[i]>>0)&1) == 0 ) continue;
    if (((mcStatusFlag[i]>>1)&1) == 0 ) continue;

    if (deltaR(phoEta, phoPhi, mcEta[i], mcPhi[i]) < 0.1 ){
      if (abs(mcPID[i]) == 22){
	matched = 1; // matched
	break;
      }
    }
  }

  if (matched == 1) return 1;
  if (matched == 0) return 0;

} // end of function



int PhoMatchXZGSample(TreeReader &data, Float_t phoEta, Float_t phoPhi){

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
    if (deltaR(phoEta, phoPhi, mcEta[i], mcPhi[i]) < 0.1 ){
      if (abs(mcPID[i]) == 22 && abs(mcMomPID[i]) == 25){
        matched = 1; // matched
        break;
      }
    }
  }

  if (matched == 1) return 1;
  if (matched == 0) return 0;

} // end of function
