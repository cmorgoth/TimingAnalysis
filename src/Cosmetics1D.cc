#include "Cosmetics1D.hh"

bool Makeup(TH1F* h, TString fname = "dummy", TString xTitle = "xTitle", TString yTitle = "yTitle"){
  TCanvas* C1 = new TCanvas("C1", "C1", 400, 500);
  C1->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,1);
  pad2->SetTopMargin(0.1);
  pad2->SetBottomMargin(0.15);
  pad2->SetLeftMargin(0.15);
  pad2->SetRightMargin(0.05);
  pad2->Draw();
  pad2->cd();
  h->SetTitleOffset(1.3, "X");
  h->SetTitleOffset(1.4, "Y");
  h->SetTitleSize(0.04, "X");
  h->SetTitleSize(0.04, "Y");
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->SetTitle("");
  h->SetTitle("");
  h->SetYTitle(yTitle);
  h->SetXTitle(xTitle);
  h->Draw();
  C1->Modified();
  C1->Update();
  TPaveStats* ps = (TPaveStats*)h->GetListOfFunctions()->FindObject("stats");
  ps->SetX1NDC(0.7);
  ps->SetX2NDC(0.95);
  ps->SetY1NDC(0.65);
  ps->SetY2NDC(0.89);
  ps->SetTextSize(.025);
  C1->SaveAs(fname+".png");
  C1->SaveAs(fname+".pdf");
  C1->SaveAs(fname+".C");
  C1->SaveAs(fname+".root");
  
  return true;
};
