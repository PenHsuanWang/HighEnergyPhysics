{
  gSystem->SetBuildDir("tmpdir", kTRUE);

  gROOT->LoadMacro("xAna.C+");

  gStyle->SetOptStat(0);


  for (Int_t i=13 ; i<=24 ; i++){


    if (i == 1){
      char inpaths[300] =
        {"/data6/ggNtuples/V08_00_26_05/job_spring16_Zdgamma_ZToEE_m_1/" };
      xAna(inpaths, "minitree_Zdgamma_ZToEE_m_1", 0, -1, -1, 1);
    }

    if (i == 2){
      char inpaths[300] =
        {"/data6/ggNtuples/V08_00_26_05/job_spring16_Zdgamma_ZToEE_m_3/" };
      xAna(inpaths, "minitree_Zdgamma_ZToEE_m_3", 0, -1, -1, 1);
    }

    if (i == 3){
      char inpaths[300] =
        {"/data6/ggNtuples/V08_00_26_05/job_spring16_Zdgamma_ZToEE_m_5/" };
      xAna(inpaths, "minitree_Zdgamma_ZToEE_m_5", 0, -1, -1, 1);
    }

    if (i == 4){
      char inpaths[300] =
        {"/data6/ggNtuples/V08_00_26_05/job_spring16_Zdgamma_ZToEE_m_10/" };
      xAna(inpaths, "minitree_Zdgamma_ZToEE_m_10", 0, -1, -1, 1);
    }

    if (i == 5){
      char inpaths[300] =
        {"/data6/ggNtuples/V08_00_26_05/job_spring16_Zdgamma_ZToEE_m_15/" };
      xAna(inpaths, "minitree_Zdgamma_ZToEE_m_15", 0, -1, -1, 1);
    }

    if (i == 6){
      char inpaths[300] =
        {"/data6/ggNtuples/V08_00_26_05/job_spring16_Zdgamma_ZToEE_m_20/" };
      xAna(inpaths, "minitree_Zdgamma_ZToEE_m_20", 0, -1, -1, 1);
    }

    if (i == 7){
      char inpaths[300] =
        {"/data6/ggNtuples/V08_00_26_05/job_spring16_Zdgamma_ZToEE_m_25/" };
      xAna(inpaths, "minitree_Zdgamma_ZToEE_m_25", 0, -1, -1, 1);
    }

    if (i == 8){
      char inpaths[300] =
        {"/data6/ggNtuples/V08_00_26_05/job_spring16_Zdgamma_ZToEE_m_30/" };
      xAna(inpaths, "minitree_Zdgamma_ZToEE_m_30", 0, -1, -1, 1);
    }

    if (i == 9){
      char inpaths[300] =
        {"/data6/ggNtuples/V08_00_26_05/job_spring16_Zdgamma_ZToEE_m_35/" };
      xAna(inpaths, "minitree_Zdgamma_ZToEE_m_35", 0, -1, -1, 1);
    }

    if (i == 10){
      char inpaths[300] =
        {"/data6/ggNtuples/V08_00_26_05/job_spring16_Zdgamma_ZToEE_m_40/" };
      xAna(inpaths, "minitree_Zdgamma_ZToEE_m_40", 0, -1, -1, 1);
    }

    if (i == 11){
      char inpaths[300] =
        {"/data6/ggNtuples/V08_00_26_05/job_spring16_Zdgamma_ZToEE_m_45/" };
      xAna(inpaths, "minitree_Zdgamma_ZToEE_m_45", 0, -1, -1, 1);
    }

    if (i == 12){
      char inpaths[300] =
        {"/data6/ggNtuples/V08_00_26_05/job_spring16_Zdgamma_ZToEE_m_50/" };
      xAna(inpaths, "minitree_Zdgamma_ZToEE_m_50", 0, -1, -1, 1);
    }
    

    if (i == 13){
      char inpaths[300] =
        {"/data1/pwang/MC/V08_00_26_05/job_spring16_Zdgamma_ZToMuMu_m_1/" };
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_1", 0, -1, -1, 0, 0);
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_1", 0, -1, -1, 0, 0, "phoScal_stat_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_1", 0, -1, -1, 0, 0, "phoScal_stat_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_1", 0, -1, -1, 0, 0, "phoScal_syst_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_1", 0, -1, -1, 0, 0, "phoScal_syst_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_1", 0, -1, -1, 0, 0, "phoScal_gain_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_1", 0, -1, -1, 0, 0, "phoScal_gain_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_1", 0, -1, -1, 0, 0, "phoSmear_rho_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_1", 0, -1, -1, 0, 0, "phoSmear_rho_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_1", 0, -1, -1, 0, 0, "phoSmear_phi_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_1", 0, -1, -1, 0, 0, "phoSmear_phi_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_1", 0, -1, -1, 0, 0, "PU_un_up"        );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_1", 0, -1, -1, 0, 0, "PU_un_dn"        );
    }

    if (i == 14){
      char inpaths[300] =
        {"/data1/pwang/MC/V08_00_26_05/job_spring16_Zdgamma_ZToMuMu_m_3/" };
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_3", 0, -1, -1, 0, 0);
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_3", 0, -1, -1, 0, 0, "phoScal_stat_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_3", 0, -1, -1, 0, 0, "phoScal_stat_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_3", 0, -1, -1, 0, 0, "phoScal_syst_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_3", 0, -1, -1, 0, 0, "phoScal_syst_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_3", 0, -1, -1, 0, 0, "phoScal_gain_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_3", 0, -1, -1, 0, 0, "phoScal_gain_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_3", 0, -1, -1, 0, 0, "phoSmear_rho_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_3", 0, -1, -1, 0, 0, "phoSmear_rho_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_3", 0, -1, -1, 0, 0, "phoSmear_phi_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_3", 0, -1, -1, 0, 0, "phoSmear_phi_dn" );
    }

    if (i == 15){
      char inpaths[300] =
        {"/data1/pwang/MC/V08_00_26_05/job_spring16_Zdgamma_ZToMuMu_m_5/" };
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_5", 0, -1, -1, 0, 0);
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_5", 0, -1, -1, 0, 0, "phoScal_stat_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_5", 0, -1, -1, 0, 0, "phoScal_stat_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_5", 0, -1, -1, 0, 0, "phoScal_syst_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_5", 0, -1, -1, 0, 0, "phoScal_syst_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_5", 0, -1, -1, 0, 0, "phoScal_gain_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_5", 0, -1, -1, 0, 0, "phoScal_gain_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_5", 0, -1, -1, 0, 0, "phoSmear_rho_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_5", 0, -1, -1, 0, 0, "phoSmear_rho_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_5", 0, -1, -1, 0, 0, "phoSmear_phi_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_5", 0, -1, -1, 0, 0, "phoSmear_phi_dn" );
    }

    if (i == 16){
      char inpaths[300] =
        {"/data1/pwang/MC/V08_00_26_05/job_spring16_Zdgamma_ZToMuMu_m_10/" };
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_10", 0, -1, -1, 0, 0);
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_10", 0, -1, -1, 0, 0, "phoScal_stat_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_10", 0, -1, -1, 0, 0, "phoScal_stat_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_10", 0, -1, -1, 0, 0, "phoScal_syst_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_10", 0, -1, -1, 0, 0, "phoScal_syst_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_10", 0, -1, -1, 0, 0, "phoScal_gain_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_10", 0, -1, -1, 0, 0, "phoScal_gain_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_10", 0, -1, -1, 0, 0, "phoSmear_rho_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_10", 0, -1, -1, 0, 0, "phoSmear_rho_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_10", 0, -1, -1, 0, 0, "phoSmear_phi_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_10", 0, -1, -1, 0, 0, "phoSmear_phi_dn" );
    }

    if (i == 17){
      char inpaths[300] =
        {"/data1/pwang/MC/V08_00_26_05/job_spring16_Zdgamma_ZToMuMu_m_15/" };
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_15", 0, -1, -1, 0, 0);
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_15", 0, -1, -1, 0, 0, "phoScal_stat_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_15", 0, -1, -1, 0, 0, "phoScal_stat_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_15", 0, -1, -1, 0, 0, "phoScal_syst_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_15", 0, -1, -1, 0, 0, "phoScal_syst_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_15", 0, -1, -1, 0, 0, "phoScal_gain_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_15", 0, -1, -1, 0, 0, "phoScal_gain_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_15", 0, -1, -1, 0, 0, "phoSmear_rho_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_15", 0, -1, -1, 0, 0, "phoSmear_rho_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_15", 0, -1, -1, 0, 0, "phoSmear_phi_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_15", 0, -1, -1, 0, 0, "phoSmear_phi_dn" );
    }

    if (i == 18){
      char inpaths[300] =
        {"/data1/pwang/MC/V08_00_26_05/job_spring16_Zdgamma_ZToMuMu_m_20/" };
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_20", 0, -1, -1, 0, 0);
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_20", 0, -1, -1, 0, 0, "phoScal_stat_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_20", 0, -1, -1, 0, 0, "phoScal_stat_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_20", 0, -1, -1, 0, 0, "phoScal_syst_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_20", 0, -1, -1, 0, 0, "phoScal_syst_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_20", 0, -1, -1, 0, 0, "phoScal_gain_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_20", 0, -1, -1, 0, 0, "phoScal_gain_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_20", 0, -1, -1, 0, 0, "phoSmear_rho_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_20", 0, -1, -1, 0, 0, "phoSmear_rho_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_20", 0, -1, -1, 0, 0, "phoSmear_phi_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_20", 0, -1, -1, 0, 0, "phoSmear_phi_dn" );
    }

    if (i == 19){
      char inpaths[300] =
        {"/data1/pwang/MC/V08_00_26_05/job_spring16_Zdgamma_ZToMuMu_m_25/" };
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_25", 0, -1, -1, 0, 0);
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_25", 0, -1, -1, 0, 0, "phoScal_stat_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_25", 0, -1, -1, 0, 0, "phoScal_stat_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_25", 0, -1, -1, 0, 0, "phoScal_syst_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_25", 0, -1, -1, 0, 0, "phoScal_syst_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_25", 0, -1, -1, 0, 0, "phoScal_gain_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_25", 0, -1, -1, 0, 0, "phoScal_gain_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_25", 0, -1, -1, 0, 0, "phoSmear_rho_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_25", 0, -1, -1, 0, 0, "phoSmear_rho_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_25", 0, -1, -1, 0, 0, "phoSmear_phi_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_25", 0, -1, -1, 0, 0, "phoSmear_phi_dn" );
    }

    if (i == 20){
      char inpaths[300] =
        {"/data1/pwang/MC/V08_00_26_05/job_spring16_Zdgamma_ZToMuMu_m_30/" };
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_30", 0, -1, -1, 0, 0);
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_30", 0, -1, -1, 0, 0, "phoScal_stat_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_30", 0, -1, -1, 0, 0, "phoScal_stat_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_30", 0, -1, -1, 0, 0, "phoScal_syst_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_30", 0, -1, -1, 0, 0, "phoScal_syst_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_30", 0, -1, -1, 0, 0, "phoScal_gain_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_30", 0, -1, -1, 0, 0, "phoScal_gain_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_30", 0, -1, -1, 0, 0, "phoSmear_rho_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_30", 0, -1, -1, 0, 0, "phoSmear_rho_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_30", 0, -1, -1, 0, 0, "phoSmear_phi_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_30", 0, -1, -1, 0, 0, "phoSmear_phi_dn" );
    }

    if (i == 21){
      char inpaths[300] =
        {"/data1/pwang/MC/V08_00_26_05/job_spring16_Zdgamma_ZToMuMu_m_35/" };
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_35", 0, -1, -1, 0, 0);
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_35", 0, -1, -1, 0, 0, "phoScal_stat_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_35", 0, -1, -1, 0, 0, "phoScal_stat_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_35", 0, -1, -1, 0, 0, "phoScal_syst_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_35", 0, -1, -1, 0, 0, "phoScal_syst_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_35", 0, -1, -1, 0, 0, "phoScal_gain_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_35", 0, -1, -1, 0, 0, "phoScal_gain_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_35", 0, -1, -1, 0, 0, "phoSmear_rho_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_35", 0, -1, -1, 0, 0, "phoSmear_rho_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_35", 0, -1, -1, 0, 0, "phoSmear_phi_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_35", 0, -1, -1, 0, 0, "phoSmear_phi_dn" );
    }

    if (i == 22){
      char inpaths[300] =
        {"/data1/pwang/MC/V08_00_26_05/job_spring16_Zdgamma_ZToMuMu_m_40/" };
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_40", 0, -1, -1, 0, 0);
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_40", 0, -1, -1, 0, 0, "phoScal_stat_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_40", 0, -1, -1, 0, 0, "phoScal_stat_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_40", 0, -1, -1, 0, 0, "phoScal_syst_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_40", 0, -1, -1, 0, 0, "phoScal_syst_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_40", 0, -1, -1, 0, 0, "phoScal_gain_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_40", 0, -1, -1, 0, 0, "phoScal_gain_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_40", 0, -1, -1, 0, 0, "phoSmear_rho_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_40", 0, -1, -1, 0, 0, "phoSmear_rho_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_40", 0, -1, -1, 0, 0, "phoSmear_phi_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_40", 0, -1, -1, 0, 0, "phoSmear_phi_dn" );
    }

    if (i == 23){
      char inpaths[300] =
        {"/data1/pwang/MC/V08_00_26_05/job_spring16_Zdgamma_ZToMuMu_m_45/" };
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_45", 0, -1, -1, 0, 0);
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_45", 0, -1, -1, 0, 0, "phoScal_stat_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_45", 0, -1, -1, 0, 0, "phoScal_stat_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_45", 0, -1, -1, 0, 0, "phoScal_syst_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_45", 0, -1, -1, 0, 0, "phoScal_syst_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_45", 0, -1, -1, 0, 0, "phoScal_gain_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_45", 0, -1, -1, 0, 0, "phoScal_gain_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_45", 0, -1, -1, 0, 0, "phoSmear_rho_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_45", 0, -1, -1, 0, 0, "phoSmear_rho_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_45", 0, -1, -1, 0, 0, "phoSmear_phi_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_45", 0, -1, -1, 0, 0, "phoSmear_phi_dn" );
    }

    if (i == 24){
      char inpaths[300] =
        {"/data1/pwang/MC/V08_00_26_05/job_spring16_Zdgamma_ZToMuMu_m_50/" };
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_50", 0, -1, -1, 0, 0);
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_50", 0, -1, -1, 0, 0, "phoScal_stat_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_50", 0, -1, -1, 0, 0, "phoScal_stat_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_50", 0, -1, -1, 0, 0, "phoScal_syst_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_50", 0, -1, -1, 0, 0, "phoScal_syst_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_50", 0, -1, -1, 0, 0, "phoScal_gain_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_50", 0, -1, -1, 0, 0, "phoScal_gain_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_50", 0, -1, -1, 0, 0, "phoSmear_rho_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_50", 0, -1, -1, 0, 0, "phoSmear_rho_dn" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_50", 0, -1, -1, 0, 0, "phoSmear_phi_up" );
      xAna(inpaths, "minitree_Zdgamma_ZToMuMu_m_50", 0, -1, -1, 0, 0, "phoSmear_phi_dn" );
    }

  }
}
