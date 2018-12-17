

bool TrkHighPtID(TreeReader &data, int i ){
  Float_t* muD0             = data.GetPtrFloat("muD0"); // Impact parameter about xy plane                                                                
  Float_t* muDz             = data.GetPtrFloat("muDz"); // Impact parameter about z axis                                                                  
  Int_t*   muMatches        = data.GetPtrInt("muMatches"); // muon number of match                                                                        
  Int_t*   muBestTrkType    = data.GetPtrInt("muBestTrkType");
  Int_t*   muStations       = data.GetPtrInt("muStations");
  Float_t* muBestTrkPtError = data.GetPtrFloat("muBestTrkPtError"); // ptError                                                                            
  Float_t* muBestTrkPt      = data.GetPtrFloat("muBestTrkPt"); // best track pt                                                                           
  Int_t*   muPixelHits      = data.GetPtrInt("muPixelHits");
  Int_t*   muTrkLayers      = data.GetPtrInt("muTrkLayers");

  if (muMatches[i] <= 1) return false;
  if (muBestTrkPtError[i]/muBestTrkPt[i] >= 0.3) return false;
  if (muD0[i] >= 0.2 ) return false;
  if (muDz[i] >= 0.5 ) return false;
  if (muPixelHits[i] == 0) return false;
  if (muTrkLayers[i] <= 5) return false;
  return true;
}

