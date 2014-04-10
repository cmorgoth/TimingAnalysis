#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TString.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"

#include "Cosmetics1D.hh"

int main(int argc,  char** argv){
    
  std::cout<< "argc: " << argc << std::endl;
  TString fname(argv[1]);
  std::cout << fname << std::endl;
  gStyle->SetOptStat(1100);
  gStyle->SetOptFit(11);
  
  gROOT->Reset();
  
  TFile* f = new TFile(fname);
  TTree* tree =  (TTree*)f->Get("tree");
  
  tree->Draw("t2gausroot - t1gausroot>>dT(50,-0.6,0.6)",
	     "ch1QualityBit==0 && ch2QualityBit==0 && ch4QualityBit>=0", "goff");
  TH1F* dT = (TH1F*)gDirectory->Get("dT");
  
  tree->Draw("t1gausroot - t4gausroot>>tof1(60,10.3,10.9)",
	     "ch1QualityBit==0 && ch2QualityBit>=0 && ch4QualityBit==0", "goff");
  TH1F* tof1 = (TH1F*)gDirectory->Get("tof1");

  tree->Draw("t3gausroot - t4gausroot>>tof3(45,9.5,9.95)",
	     "ch3QualityBit==0 && ch2QualityBit>=0 && ch4QualityBit==0", "goff");
  TH1F* tof3 = (TH1F*)gDirectory->Get("tof3");
  
  tree->Draw("t2gausroot - t4gausroot>>tof2(20,10.2,10.6)",
	     "ch1QualityBit>=0 && ch2QualityBit==0 && ch4QualityBit==0", "goff");
  TH1F* tof2 = (TH1F*)gDirectory->Get("tof2");
  
  tree->Draw("(t4gausroot - t2gausroot):(t2gausroot - t1gausroot)>>tof2_dT(20, 0.35, 0.75, 17, 4.4, 4.75)"
	     ,"ch1QualityBit==0 && ch2QualityBit==0 && ch4QualityBit==0", "goff");
  TH2F* tof2_dT = (TH2F*)gDirectory->Get("tof2_dT");
  
  TProfile* profX = tof2_dT->ProfileX("profx", 1, -1, "o");
  profX->SetAxisRange(4.4, 4.75, "Y");
  profX->Fit("pol1");
  TF1* g = new TF1("g","gaus(0)",0.49,0.6);
  TF1* g1 = new TF1("g1","gaus(0)",10.3,10.9);
  TF1* g2 = new TF1("g2","gaus(0)",4.475,4.7);
  TF1* g3 = new TF1("g3","gaus(0)",9.57,9.85);
  g->SetParNames("Norm_{fit}","#mu_{fit}", "#sigma_{fit}");
  g1->SetParNames("Norm_{fit}","#mu_{fit}", "#sigma_{fit}");
  g2->SetParNames("Norm_{fit}","#mu_{fit}", "#sigma_{fit}");
  g3->SetParNames("Norm_{fit}","#mu_{fit}", "#sigma_{fit}");

  //TString sn("Perp_8GeV_Ele_run123/");
  TString sn("Perp_8GeV_Ele_run122/");
  
  TCanvas* C = new TCanvas("C", "C	", 400, 500);
  C->cd();
  
  g->SetParameter(0,dT->GetMaximum());
  g->SetParameter(1,dT->GetMean());
  g->SetParameter(2,dT->GetRMS());
  dT->Fit(g,"MWLR");
  Makeup(dT, sn+"dT", "t_{2} - t_{1} (nsec)", "Events/10ps");
  
  g1->SetParameter(0,tof1->GetMaximum());
  g1->SetParameter(1,tof1->GetMean());
  g1->SetParameter(2,tof1->GetRMS());
  tof1->Fit(g1,"MWLR");
  Makeup(tof1, sn+"TOF1", "t_{1} - t_{0} (nsec)", "Events/10ps");

  g3->SetParameter(0,tof3->GetMaximum());
  g3->SetParameter(1,tof3->GetMean());
  g3->SetParameter(2,tof3->GetRMS());
  tof3->Fit(g3,"MWLR");
  Makeup(tof3, sn+"TOF3", "t_{3} - t_{0} (nsec)", "Events/10ps");
  
  g2->SetParameter(0,tof2->GetMaximum());
  g2->SetParameter(1,tof2->GetMean());
  g2->SetParameter(2,tof2->GetRMS());
  tof2->Fit(g2,"MWLR");
  Makeup(tof1, sn+"TOF2", "t_{2} - t_{0} (nsec)", "Events/10ps");  
  
  tof2_dT->Draw("colz");
  C->SaveAs(sn+"tof2_dT.png");
  C->SaveAs(sn+"tof2_dT.pdf");
  C->SaveAs(sn+"tof2_dT.C");
  C->SaveAs(sn+"tof2_dT.root");

  profX->Draw();
  C->SaveAs(sn+"profX.png");
  C->SaveAs(sn+"profX.pdf");
  C->SaveAs(sn+"profX.C");
  C->SaveAs(sn+"profX.root");
  return 0;

}
