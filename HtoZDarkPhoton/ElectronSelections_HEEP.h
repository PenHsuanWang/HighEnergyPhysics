float GetEleEt(float elePt){
  //Float_t Et = sqrt(elePt*elePt+(0.000511*0.000511));
  return elePt;
}

int isHEEPele(TreeReader &data, int i){

  Int_t    nEle               = data.GetInt("nEle");
  Float_t* elePt              = data.GetPtrFloat("elePt");
  Float_t* eleSCEn              = data.GetPtrFloat("eleSCEn");
  Float_t* eleSCEta           = data.GetPtrFloat("eleSCEta");
  Float_t* eleSigmaIEtaIEtaFull5x5   = data.GetPtrFloat("eleSigmaIEtaIEtaFull5x5");
  Float_t* eleE1x5Full5x5       = data.GetPtrFloat("eleE1x5Full5x5");
  Float_t* eleE2x5Full5x5       = data.GetPtrFloat("eleE2x5Full5x5");
  Float_t* eleE5x5Full5x5       = data.GetPtrFloat("eleE5x5Full5x5");
  Float_t* eleHoverE          = data.GetPtrFloat("eleHoverE");
  Float_t* eleDr03TkSumPt     = data.GetPtrFloat("eleDr03TkSumPt");
  Float_t* eleDr03HcalDepth1TowerSumEt = data.GetPtrFloat("eleDr03HcalDepth1TowerSumEt");
  Float_t* eleDr03EcalRecHitSumEt = data.GetPtrFloat("eleDr03EcalRecHitSumEt");
  Float_t* eledEtaseedAtVtx   = data.GetPtrFloat("eledEtaseedAtVtx");
  Int_t* eleEcalDrivenSeed  = data.GetPtrInt("eleEcalDrivenSeed");
  Int_t*   eleMissHits        = data.GetPtrInt("eleMissHits");
  Float_t* eleD0              = data.GetPtrFloat("eleD0");
  Float_t* eleDz              = data.GetPtrFloat("eleDz");
  Float_t* eleEoverPInv       = data.GetPtrFloat("eleEoverPInv");
  Float_t* eledEtaAtVtx       = data.GetPtrFloat("eledEtaAtVtx");
  Float_t* eledPhiAtVtx       = data.GetPtrFloat("eledPhiAtVtx");
  Int_t*   eleConvVeto        = data.GetPtrInt("eleConvVeto");
  Float_t* elePFChIso         = data.GetPtrFloat("elePFChIso");
  Float_t* elePFPhoIso        = data.GetPtrFloat("elePFPhoIso");
  Float_t* elePFNeuIso        = data.GetPtrFloat("elePFNeuIso");
  Float_t* elePFPUIso         = data.GetPtrFloat("elePFPUIso");
  Float_t  rho                = data.GetFloat("rho");



  UShort_t* eleIDbit = (UShort_t*) data.GetPtrShort("eleIDbit");

    if ( elePt[i] < 15.) return 0; // pt cut
    if (fabs(eleSCEta[i]) > 1.4442 && fabs(eleSCEta[i]) < 1.566) return 0; // EE and EB gap
    if (fabs(eleSCEta[i]) > 2.5) return 0; //eta range

    if ( fabs(eleSCEta[i]) <= 1.4442 && eleEcalDrivenSeed[i] != 1 ) return 0;
    if ( fabs(eleSCEta[i]) <= 1.4442 && fabs(eledEtaseedAtVtx[i]) >= 0.004 ) return 0;
    if ( fabs(eleSCEta[i]) <= 1.4442 && fabs(eledPhiAtVtx[i]) >= 0.06 ) return 0;
    if ( fabs(eleSCEta[i]) <= 1.4442 && eleHoverE[i] >= ((1/eleSCEn[i])+0.05) ) return 0;
    // if ( fabs(eleSCEta[i]) <= 1.4442 && eeleSigmaIEtaIEta[i] > HEEPcut[0][4] ) return 0;   n/a
    if ( fabs(eleSCEta[i]) <= 1.4442 && (eleE2x5Full5x5[i]/eleE5x5Full5x5[i]) <= 0.94 && (eleE1x5Full5x5[i]/eleE5x5Full5x5[i]) <= 0.83 ) return 0;
    //if ( fabs(eleSCEta[i]) <= 1.4442 && (eleDr03EcalRecHitSumEt[i]+eleDr03HcalDepth1TowerSumEt[i]) >= ((2+(0.03*GetEleEt(elePt[i])) + (0.28*rho))) ) return 0;
    if ( fabs(eleSCEta[i]) <= 1.4442 && GetEleEt(elePt[i]) <95 && eleDr03TkSumPt[i] >= 5 ) return 0;
    if ( fabs(eleSCEta[i]) <= 1.4442 && GetEleEt(elePt[i]) >=95 && eleDr03TkSumPt[i] >= 5 + 1.5*rho ) return 0;

    if ( fabs(eleSCEta[i]) <= 1.4442 && eleMissHits[i] > 1 ) return 0;
    if ( fabs(eleSCEta[i]) <= 1.4442 && fabs(eleD0[i]) >= 0.02 ) return 0;


    if ( fabs(eleSCEta[i]) >= 1.566 && eleEcalDrivenSeed[i] != 1 ) return 0;
    if ( fabs(eleSCEta[i]) >= 1.566 && fabs(eledEtaseedAtVtx[i]) >= 0.006 ) return 0;
    if ( fabs(eleSCEta[i]) >= 1.566 && fabs(eledPhiAtVtx[i]) >= 0.06 ) return 0;
    if ( fabs(eleSCEta[i]) >= 1.566 && eleHoverE[i] >= ((5/eleSCEn[i])+0.05) ) return 0;
    if ( fabs(eleSCEta[i]) >= 1.566 && eleSigmaIEtaIEtaFull5x5[i] >= 0.03 ) return 0;
    //if ( fabs(eleSCEta[i]) >= 1.566 && (eleE2x5Full5x5[i]/eleE5x5Full5x5[i]) >= HEEPcut[1][5] ) return 0;   n/a
    //if ( fabs(eleSCEta[i]) >= 1.566 && GetEleEt(elePt[i])<50 && (eleDr03EcalRecHitSumEt[i]+eleDr03HcalDepth1TowerSumEt[i]) >= ((2.5+(0.28*rho))) ) return 0;

    if ( fabs(eleSCEta[i]) >= 1.566 && GetEleEt(elePt[i])>=50 
	 && eleDr03EcalRecHitSumEt[i]+eleDr03HcalDepth1TowerSumEt[i] >= (2.5+(0.03*(GetEleEt(elePt[i])-50))+(0.28*rho)) ) return 0;
    
    if ( fabs(eleSCEta[i]) >= 1.566 && GetEleEt(elePt[i]) <100 && eleDr03TkSumPt[i] >= 5 ) return 0;
    if ( fabs(eleSCEta[i]) >= 1.566 && GetEleEt(elePt[i]) >=100 && eleDr03TkSumPt[i] >= 5 + 1.5*rho ) return 0;

    if ( fabs(eleSCEta[i]) >= 1.566 && eleMissHits[i] > 1 ) return 0;
    if ( fabs(eleSCEta[i]) >= 1.566 && fabs(eleD0[i]) >= 0.05 ) return 0;

    return 1;

}
