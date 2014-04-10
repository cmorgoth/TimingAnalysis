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
#include "TPaveStats.h"
#include "TPad.h"

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
  
  tree->Draw("t2gausroot - t1gausroot>>dT(35,0.3,1.1)",
	     "ch1QualityBit==0 && ch2QualityBit==0 && ch3QualityBit>=0", "goff");
  TH1F* dT = (TH1F*)gDirectory->Get("dT");
  //tof1
  tree->Draw("t4gausroot - t1gausroot>>tof1(60,-4.7,15.3)",
	     "ch1QualityBit==0 && ch4QualityBit>=0 && ch3QualityBit>=0", "goff");
  TH1F* tof1 = (TH1F*)gDirectory->Get("tof1");
  //tof1_amp0.4
  tree->Draw("t3gausroot - t1gausroot>>tof1_0p04(30,4.7,5.3)",
	     "ch1QualityBit==0 && ch2QualityBit>=0 && ch3QualityBit==0 && ch1Amp>0.01 && ch1Amp<0.038", "goff");
  TH1F* tof1_0p04 = (TH1F*)gDirectory->Get("tof1_0p04");
  //tof1_amp_0.6
  tree->Draw("t3gausroot - t1gausroot>>tof1_0p06(30,4.7,5.3)",
	     "ch1QualityBit==0 && ch2QualityBit>=0 && ch3QualityBit==0 && ch1Amp > 0.038 && ch1Amp<0.063", "goff");
  TH1F* tof1_0p06 = (TH1F*)gDirectory->Get("tof1_0p06");
  //tof1_amp1.0
  tree->Draw("t3gausroot - t1gausroot>>tof1_0p1(30,4.7,5.3)",
	     "ch1QualityBit==0 && ch2QualityBit>=0 && ch3QualityBit==0 && ch1Amp > 0.063 && ch1Amp<0.105", "goff");
  TH1F* tof1_0p1 = (TH1F*)gDirectory->Get("tof1_0p1");
   //tof1_amp2.0
  tree->Draw("t3gausroot - t1gausroot>>tof1_0p2(30,4.7,5.3)",
	     "ch1QualityBit==0 && ch2QualityBit>=0 && ch3QualityBit==0 && ch1Amp>0.105 && ch1Amp<0.160", "goff");
  TH1F* tof1_0p2 = (TH1F*)gDirectory->Get("tof1_0p2");
  
  //tof1_amp5.0
  tree->Draw("t3gausroot - t1gausroot>>tof1_0p5(30,4.7,5.3)",
	     "ch1QualityBit==0 && ch2QualityBit>=0 && ch3QualityBit==0 && ch1Amp>0.16 && ch1Amp<0.25", "goff");
  TH1F* tof1_0p5 = (TH1F*)gDirectory->Get("tof1_0p5");

  //tof1_ampf
  tree->Draw("t3gausroot - t1gausroot>>tof1_0pf(30,4.7,5.3)",
	     "ch1QualityBit==0 && ch2QualityBit>=0 && ch3QualityBit==0 && ch1Amp>0.25 && ch1Amp<0.5", "goff");
  TH1F* tof1_0pf = (TH1F*)gDirectory->Get("tof1_0pf");
  
  //tof2
  tree->Draw("t4gausroot - t2gausroot>>tof2(75,-4.2,14.95)",
	     "ch1QualityBit==0 && ch4QualityBit==0 && ch3QualityBit>=0", "goff");
  TH1F* tof2 = (TH1F*)gDirectory->Get("tof2");

  //ch1Amp
  tree->Draw("ch1Amp>>ch1Amp(100,0.01,0.5)",
	     "ch1QualityBit==0 && ch2QualityBit==0 && ch3QualityBit==0", "goff");
  TH1F* ch1Amp = (TH1F*)gDirectory->Get("ch1Amp");
  //ch2Amp
  tree->Draw("ch2Amp>>ch2Amp(100,0.01,0.5)",
	     "ch1QualityBit==0 && ch2QualityBit==0 && ch3QualityBit==0", "goff");
  TH1F* ch2Amp = (TH1F*)gDirectory->Get("ch2Amp");
  //ch3Amp
  tree->Draw("ch3Amp>>ch3Amp(100,0.01,0.5)",
	     "ch1QualityBit==0 && ch2QualityBit==0 && ch3QualityBit==0", "goff");
  TH1F* ch3Amp = (TH1F*)gDirectory->Get("ch3Amp");
  
  tree->Draw("(t3gausroot - t2gausroot):(t2gausroot - t1gausroot)>>tof2_dT(20, 0.35, 0.75, 17, 4.4, 4.75)"
	     ,"ch1QualityBit==0 && ch2QualityBit==0 && ch3QualityBit==0", "goff");
  TH2F* tof2_dT = (TH2F*)gDirectory->Get("tof2_dT");
  
  TProfile* profX = tof2_dT->ProfileX("profx", 1, -1, "o");
  profX->SetAxisRange(4.4, 4.75, "Y");
  profX->Fit("pol1");
  
  float nsigma = 2.80;

  TF1* g = new TF1("g","gaus(0)",dT->GetMean()-nsigma*dT->GetRMS(),dT->GetMean()+nsigma*dT->GetRMS());
  TF1* g1 = new TF1("g1","gaus(0)",tof1->GetMean()-nsigma*tof1->GetRMS(),tof1->GetMean()+nsigma*tof1->GetRMS());
  TF1* g2 = new TF1("g2","gaus(0)",tof2->GetMean()-nsigma*tof2->GetRMS(),tof2->GetMean()+nsigma*tof2->GetRMS());
  g->SetParNames("Norm_{fit}","#mu_{fit}", "#sigma_{fit}");
  g1->SetParNames("Norm_{fit}","#mu_{fit}", "#sigma_{fit}");
  g2->SetParNames("Norm_{fit}","#mu_{fit}", "#sigma_{fit}");
  TString sn("Long_120GeV_Proton_TOF_Dist/");
  
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
  Makeup(tof1, sn+"TOF1", "t_{0} - t_{1} (nsec)", "Events/10ps");
  
  tof1_0p04->Fit(g1,"MWLR");
  Makeup(tof1_0p04, sn+"TOF1_10to38mv", "t_{0} - t_{1} (nsec)", "Events/20ps");
 
  tof1_0p06->Fit(g1,"MWLR");
  Makeup(tof1_0p06, sn+"TOF1_38to63mv", "t_{0} - t_{1} (nsec)", "Events/20ps");
 
  tof1_0p1->Fit(g1,"MWLR");
  Makeup(tof1_0p1, sn+"TOF1_63to105mv", "t_{0} - t_{1} (nsec)", "Events/20ps");
  
  tof1_0p2->Fit(g1,"MWLR");
  Makeup(tof1_0p2, sn+"TOF1_105to160mv", "t_{0} - t_{1} (nsec)", "Events/20ps");
 
  tof1_0p5->Fit(g1,"MWLR");
  Makeup(tof1_0p5, sn+"TOF1_160to250mv", "t_{0} - t_{1} (nsec)", "Events/20ps");
  
  tof1_0pf->Fit(g1,"MWLR");
  Makeup(tof1_0pf, sn+"TOF1_250to500mv", "t_{0} - t_{1} (nsec)", "Events/20ps"); 

  g2->SetParameter(0,tof2->GetMaximum());
  g2->SetParameter(1,tof2->GetMean());
  g2->SetParameter(2,tof2->GetRMS());
  tof2->Fit(g2,"MWLR");
  Makeup(tof2, sn+"TOF2", "t_{0} - t_{2} (nsec)", "Events/10ps");
  
  tof2_dT->Draw("colz");
  C->SaveAs(sn+"tof2_dT.png");
  C->SaveAs(sn+"tof2_dT.pdf");
  C->SaveAs(sn+"tof2_dT.C");
  C->SaveAs(sn+"tof2_dT.root");
  
  profX->Draw();
  //Makeup(profX, sn+"_profX", "t_{0} - t_{1} (nsec)", "Events/10ps");
  C->SaveAs(sn+"profX.png");
  C->SaveAs(sn+"profX.pdf");
  C->SaveAs(sn+"profX.C");
  C->SaveAs(sn+"profX.root");
  
  ch1Amp->Draw();
  Makeup(ch1Amp, sn+"ch1Amp", "Amp [Volts]", "Events/5mV");
 
  ch2Amp->Draw();
  Makeup(ch1Amp, sn+"ch2Amp", "Amp [Volts]", "Events/5mV");
  
  ch3Amp->Draw();
  Makeup(ch1Amp, sn+"ch3Amp", "Amp [Volts]", "Events/5mV");
  
  return 0;

}
