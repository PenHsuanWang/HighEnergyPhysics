{
  gSystem->SetBuildDir("tmpdir", kTRUE);
  
  gROOT->LoadMacro("xAna.C+");
  
  // data 2012ABCD, double electron trigger
  const char* inpaths[] = {
    "/data3/ggNtuples/V08_00_11_01/job_SingleMu_Run2016B_PRv2.root",
    "/data3/ggNtuples/V08_00_11_01/job_SingleMu_Run2016C_PRv2.root",
    "/data3/ggNtuples/V08_00_11_01/job_SingleMu_Run2016D_PRv2.root"
  };
  
  xAna(inpaths, 3, "minitree_2016RunBCD_SingleMu.root", 0, -1, -1, 0);
}
