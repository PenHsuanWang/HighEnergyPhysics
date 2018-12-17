#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"

void DrawEfficiency(const TH1F* D, const TH1F* N){


TGraphAsymmErrors *gEff = new TGraphAsymmErrors();
gEff->BayesDivide(D, N);


TLegend *L2 = new TLegend(0.35, 0.25, 0.75, 0.45);
TCanvas *c2 = new TCanvas("c2", "c2", 700, 700);
c2->SetTitle("  ");
c2->cd(0);
gEff->GetXaxis()->SetTitle("p_{T}^{#gamma} GeV");
gEff->GetXaxis()->SetTitleSize(0.05);
gEff->GetXaxis()->SetLabelSize(0.03);
gEff->GetYaxis()->SetTitle("eleVeto Eff");
gEff->SetMarkerStyle(20);
gEff->SetMarkerSize(1);
gEff->SetMarkerColor(4);
gEff->SetLineColor(4);
gEff->Draw("ap");
c2->Print("Test.png");

}
