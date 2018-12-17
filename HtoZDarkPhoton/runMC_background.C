{
  gSystem->SetBuildDir("tmpdir", kTRUE);

  gROOT->LoadMacro("xAna.C+");

  gStyle->SetOptStat(0);


  for (Int_t i=1 ; i<=2 ; i++){


    if (i == 1){
      const char* inpaths[] =
        {"/data1/pwang/MC/V08_00_11_01/job_spring16_DYJetsToLL_m50_aMCatNLO.root" };
      xAna(inpaths, 1,"minitree_spring16_DYJetsToLL_m50_aMCatNLO_Mu.root", 0, -1, 6025.2, 0, 1);
    }

    if (i == 2){
      const char* inpaths[] =
        {"/data1/pwang/MC/V08_00_11_01/job_spring16_Zg_aMCatNLO.root",
	 "/data1/pwang/MC/V08_00_11_01/job_spring16_Zg_aMCatNLO_ext1.root"};
      xAna(inpaths, 2,"minitree_spring16_Zg_aMCatNLO_Mu.root", 0, -1, 117.864, 0, 2);
    }

    if (i == 3){
      const char* inpaths[] =
        {"/data1/pwang/MC/V08_00_11_01/job_spring16_WWTo2L2Nu.root",
	 "/data1/pwang/MC/V08_00_11_01/job_spring16_WWToLNuQQ_ext1.root",
	 "/data1/pwang/MC/V08_00_11_01/job_spring16_WWToLNuQQ.root" };
      xAna(inpaths, 3,"minitree_spring16_WW_Mu.root", 0, -1, 63.21, 0, 3);
    }

    if (i == 4){
      const char* inpaths[] =
        {"/data1/pwang/MC/V08_00_11_01/job_spring16_WZTo2L2Q.root",
	 "/data1/pwang/MC/V08_00_11_01/job_spring16_WZTo3LNu.root" };
      xAna(inpaths, 2,"minitree_spring16_WZ_Mu.root", 0, -1, 22.82, 0, 4);
    }

    if (i == 5){
      const char* inpaths[] =
        {"/data1/pwang/MC/V08_00_11_01/job_spring16_ZZTo2L2Nu.root",
	 "/data1/pwang/MC/V08_00_11_01/job_spring16_ZZTo2L2Q.root",
	 "/data1/pwang/MC/V08_00_11_01/job_spring16_ZZTo4L.root"};
      xAna(inpaths, 3,"minitree_spring16_ZZ_Mu.root", 0, -1, 10.32, 0, 5);
    }

    if (i == 6){
      const char* inpaths[] =
        {"/data1/pwang/MC/V08_00_11_01/job_spring16_TT_powheg_ext3.root" };
      xAna(inpaths, 1,"minitree_spring16_TT_powheg_ext3_Mu.root", 0, -1, 87.31, 0, 6);
    }

    if (i == 7){
      const char* inpaths[] =
        {"/data1/pwang/MC/V08_00_11_01/job_spring16_WJetsToLNu_aMCatNLO_ext1.root",
	 "/data1/pwang/MC/V08_00_11_01/job_spring16_WJetsToLNu_aMCatNLO.root"};
      xAna(inpaths, 1,"minitree_spring16_WJetsToLNu_aMCatNLO_Mu.root", 0, -1, 60290, 0, 7);
    }

    if (i == 8){
      const char* inpaths[] =
        {"/data1/pwang/MC/V08_00_11_01/job_spring16_Wg_MG.root" };
      xAna(inpaths, 1,"minitree_spring16_Wg_MG_Mu.root", 0, -1, 489, 0, 8);
    }

    if (i == 9){
      const char* inpaths[] =
	{"/data1/pwang/MC/V08_00_11_01/job_spring16_Dalitz_mmg_m125.root",
	 "/data1/pwang/MC/V08_00_11_01/job_spring16_Dalitz_mmg_m125_ext1.root" };
      xAna(inpaths, 2,"minitree_spring16_Dalitz_mmg_m125_Mu.root", 0, -1, 0.00196, 0, 9);
    }

    if (i == 11){
      const char* inpaths[] =
        {"/data1/pwang/MC/V08_00_11_01/job_QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf.root"};
      xAna(inpaths, 1,"minitree_QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_Mu.root", 0, -1, 108000000*0.000225, 0, 11);
    }

    if (i == 12){
      const char* inpaths[] =
        {"/data1/pwang/MC/V08_00_11_01/job_QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80.root"};
      xAna(inpaths, 1,"minitree_QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_Mu.root", 0, -1, 162060000*0.0016, 0, 12);
    }

    if (i == 13){
      const char* inpaths[] =
        {"/data1/pwang/MC/V08_00_11_01/job_QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf.root"};
      xAna(inpaths, 1,"minitree_QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_Mu.root", 0, -1, 54120000*0.002, 0, 13);
    }

  }
}
