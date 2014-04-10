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


  float t1, t2, t3, t4, Amp1, Amp2, Amp3, Amp4;
  unsigned int bit1, bit2, bit3, bit4;
  
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("t1gausroot",1);
  tree->SetBranchStatus("t2gausroot",1);
  tree->SetBranchStatus("t3gausroot",1);
  tree->SetBranchStatus("t4gausroot",1);
  tree->SetBranchStatus("ch1Amp",1);
  tree->SetBranchStatus("ch2Amp",1);
  tree->SetBranchStatus("ch3Amp",1);
  tree->SetBranchStatus("ch4Amp",1);
  tree->SetBranchStatus("ch1QualityBit",1);
  tree->SetBranchStatus("ch2QualityBit",1);
  tree->SetBranchStatus("ch3QualityBit",1);
  tree->SetBranchStatus("ch4QualityBit",1);

  tree->SetBranchAddress("t1gausroot", &t1);
  tree->SetBranchAddress("t2gausroot", &t2);
  tree->SetBranchAddress("t3gausroot", &t3);
  tree->SetBranchAddress("t4gausroot", &t4);
  tree->SetBranchAddress("ch1Amp", &Amp1);
  tree->SetBranchAddress("ch2Amp", &Amp2);
  tree->SetBranchAddress("ch3Amp", &Amp3);
  tree->SetBranchAddress("ch4Amp", &Amp4);
  tree->SetBranchAddress("ch1QualityBit", &bit1);
  tree->SetBranchAddress("ch2QualityBit", &bit2);
  tree->SetBranchAddress("ch3QualityBit", &bit3);
  tree->SetBranchAddress("ch4QualityBit", &bit4);
  
  /*
    h_tof[0]=ch1-ch2
    h_tof[1]=ch1-ch3
    h_tof[2]=ch1-ch4
    h_tof[3]=ch2-ch3
    h_tof[4]=ch2-ch4
    h_tof[5]=ch3-ch4
  */
  TH1F* h_tof[6];
  for(int j = 0; j < 6; j++){
    TString s(Form("tof_%d",j));
    h_tof[j] = new TH1F(s,s, 5000, -25, 25);
  }
  
  //Filling Histograms for the first time
  int nentries = tree->GetEntries();
  for(int i = 0; i < nentries; i++){
    tree->GetEntry(i);
    if(bit1==0 && bit2==0)h_tof[0]->Fill(t1-t2);
    if(bit1==0 && bit3==0)h_tof[1]->Fill(t1-t3);
    if(bit1==0 && bit4==0)h_tof[2]->Fill(t1-t4);
    if(bit2==0 && bit3==0)h_tof[3]->Fill(t2-t3);
    if(bit2==0 && bit4==0)h_tof[4]->Fill(t2-t4);
    if(bit3==0 && bit4==0)h_tof[5]->Fill(t3-t4);
  }
  for(int j = 0; j < 6; j++){
    std::cout << "Mean: " << h_tof[j]->GetMean() << 
      " Max: " << h_tof[j]->GetBinCenter(h_tof[j]->GetMaximumBin()) <<
      "  RMS: " << h_tof[j]->GetRMS() << std::endl;
  }
  
  for(int j = 0; j < 6; j++){
    std::cout << "Mean: " << h_tof[j]->GetMean() << "  RMS: " << h_tof[j]->GetRMS() << std::endl;
    TString s(Form("tof_%d",j));
    int nb = -1;
    
    //rebinning and calculating range
    if(h_tof[j]->Integral()>500){
      nb = (int)(h_tof[j]->Integral())/20;
    }else{
      nb = (int)(h_tof[j]->Integral());
    }
    //float hi = h_tof[j]->GetMean()+2.0*h_tof[j]->GetRMS();
    //float lo = h_tof[j]->GetMean()-2.0*h_tof[j]->GetRMS();
    float hi = h_tof[j]->GetBinCenter(h_tof[j]->GetMaximumBin()) + 1.0;
    float lo = h_tof[j]->GetBinCenter(h_tof[j]->GetMaximumBin()) - 1.0;
    std::cout << "nb: " << nb << std::endl;
    delete h_tof[j];
    h_tof[j] = new TH1F(s,s, nb, lo, hi);
  }
  
  //Re-Filling
  for(int i = 0; i < nentries; i++){
    tree->GetEntry(i);
    if(bit1==0 && bit2==0)h_tof[0]->Fill(t1-t2);
    if(bit1==0 && bit3==0)h_tof[1]->Fill(t1-t3);
    if(bit1==0 && bit4==0)h_tof[2]->Fill(t1-t4);
    if(bit2==0 && bit3==0)h_tof[3]->Fill(t2-t3);
    if(bit2==0 && bit4==0)h_tof[4]->Fill(t2-t4);
    if(bit3==0 && bit4==0)h_tof[5]->Fill(t3-t4);
  }
  
  //Re-iterate if RMS < 80ps
  int map_ix[6] = {-1,-1,-1,-1,-1,-1};//Map Index
  int mi_ctr = 0;//Map index counter
  for(int j = 0; j < 6; j++){
    std::cout << j <<  " Mean: " << h_tof[j]->GetMean() << "  RMS: " << h_tof[j]->GetRMS() << std::endl;
    TString s(Form("tof_%d",j));
    int nb = -1;
    if(h_tof[j]->GetRMS() < 0.20){//is smaller than 90 ps?
      map_ix[mi_ctr] = j;
      mi_ctr++;
      if(h_tof[j]->Integral()>300){
	//nb = (int)(h_tof[j]->Integral())/40;
	nb = 80;
      }else{
	nb = 80;
	//nb = (int)(h_tof[j]->Integral());
      }
      float nsigma = 4.0;
      float hi = h_tof[j]->GetMean()+nsigma*h_tof[j]->GetRMS();
      float lo = h_tof[j]->GetMean()-nsigma*h_tof[j]->GetRMS();
      //float hi = h_tof[j]->GetMean()+.5;
      //float lo = h_tof[j]->GetMean()-0.5;
      std::cout << "nb: " << nb << std::endl;
      delete h_tof[j];
      h_tof[j] = new TH1F(s,s, nb, lo, hi);
    }
  }
  
  std::cout << "Reiterate: " << mi_ctr << std::endl;
  
  for(int i = 0; i < nentries; i++){
    tree->GetEntry(i);
    //std::cout << "===========================" << std::endl;
    for(int j = 0; j < mi_ctr; j++){
      //std::cout << map_ix[j] << std::endl;
      switch(map_ix[j]){
      case 0:
	//std::cout << map_ix[j] << std::endl;
	if(bit1==0 && bit2==0)h_tof[0]->Fill(t1-t2);
	break;
      case 1:
	//std::cout << map_ix[j] << std::endl;
	if(bit1==0 && bit3==0)h_tof[1]->Fill(t1-t3);
	break;
      case 2:
	//std::cout << map_ix[j] << std::endl;
	if(bit1==0 && bit4==0)h_tof[2]->Fill(t1-t4);
	break;
      case 3:
	//std::cout << map_ix[j] << std::endl;
	if(bit2==0 && bit3==0)h_tof[3]->Fill(t2-t3);
	break;
      case 4:
	//std::cout << map_ix[j] << std::endl;
	if(bit2==0 && bit4==0)h_tof[4]->Fill(t2-t4);
	break;
      case 5:
	//std::cout << map_ix[j] << std::endl;
	if(bit3==0 && bit4==0)h_tof[5]->Fill(t3-t4);
	break;
      default:
	std::cout << "DEFAULT!!" << std::endl;
	break;
      }
    }
  }
  
  TF1* gf;
  TString sn("Long_120GeV_Proton_TOF_Dist/");
  for(int i = 0; i < mi_ctr; i++){
    float nsigma = 2.0;
    float t_low = h_tof[map_ix[i]]->GetMean() - nsigma*h_tof[map_ix[i]]->GetRMS();
    float t_high = h_tof[map_ix[i]]->GetMean() + nsigma*h_tof[map_ix[i]]->GetRMS();
    gf = new TF1("gf","gaus(0)",t_low, t_high);
    gf->SetParNames("Norm_{fit}","#mu_{fit}", "#sigma_{fit}");
    gf->SetParameter(0,h_tof[map_ix[i]]->GetMaximum());
    gf->SetParameter(1,h_tof[map_ix[i]]->GetMean());
    gf->SetParameter(2,h_tof[map_ix[i]]->GetRMS());
    h_tof[map_ix[i]]->Fit(gf,"MWLR");
    TString s(Form("tof_%d",map_ix[i]));
    Makeup(h_tof[map_ix[i]], sn+s, "t_{2} - t_{1} (nsec)", "Events/10ps");
  }


  TFile* tmp = new TFile("tmp.root", "RECREATE");
  for(int j = 0; j < 6; j++){
    h_tof[j]->Write();
  }
  
  return 0;
  
}
