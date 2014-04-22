#include <iostream>
#include <string>
#include <stdlib.h>
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
#include "TFitResult.h"

#include "Cosmetics1D.hh"

int main(int argc,  char** argv){
    
  std::cout<< "argc: " << argc << std::endl;
  if(argc == 1){
    std::cout << "Error in usage!" << std::endl;
    std::cout << "Please do:" << std::endl;
    std::cout << "./GenAna <filename.root>" << std::endl;
    std::cout << "The program will store the results in Plots/filename" << std::endl;
    return 0;
  }
  TString fname(argv[1]);
  std::string fsname(argv[1]);
  std::string ls("ls ");
  std::string mkdir("mkdir ");
  int pi = fsname.rfind("/")+1;
  int pf = fsname.find("_ana.root");
  std::cout << fsname.substr(pi,pf-pi) << std::endl;
  TString sn(fsname.substr(pi,pf-pi).c_str());
  sn = "Plots/"+sn;
  std::cout << "sn: " << sn << std::endl;
  ls += "Plots/"+fsname.substr(pi, pf-pi);
  mkdir += "Plots/"+fsname.substr(pi,pf-pi);
  std::cout << ls << std::endl;
  std::cout << mkdir << std::endl;
  //return 0;
  const char* dir = ls.c_str();
  if(system(dir) > 0){//check if directory exist if it doesn't it creates one
    std::cout << mkdir << std::endl;
    system(mkdir.c_str());
  }
  gStyle->SetOptStat(1100);
  gStyle->SetOptFit(11);
  
  gROOT->Reset();
  
  //Reads tree and variables
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
  
  //Histogram Mapping
  /*
    h_tof[0]=ch1-ch2
    h_tof[1]=ch1-ch3
    h_tof[2]=ch1-ch4
    h_tof[3]=ch2-ch3
    h_tof[4]=ch2-ch4
    h_tof[5]=ch3-ch4
  */
  
  //Xtitles for Histograms
  TString xtit[6];
  xtit[0] = "t_{1} - t_{2} (nsec)";//tof0
  xtit[1] = "t_{1} - t_{3} (nsec)";//tof1
  xtit[2] = "t_{1} - t_{4} (nsec)";//tof2
  xtit[3] = "t_{2} - t_{3} (nsec)";//tof3
  xtit[4] = "t_{2} - t_{4} (nsec)";//tof4
  xtit[5] = "t_{3} - t_{4} (nsec)";//tof5
  
  /*
    h_tof_2d[0] = tof0-tof1
    h_tof_2d[1] = tof0-tof2
    h_tof_2d[2] = tof0-tof3
    h_tof_2d[3] = tof0-tof4
    h_tof_2d[4] = tof0-tof5
    h_tof_2d[5] = tof1-tof2
    h_tof_2d[6] = tof1-tof3
    h_tof_2d[7] = tof1-tof4
    h_tof_2d[8] = tof1-tof5
    h_tof_2d[9] = tof2-tof3
    h_tof_2d[10] = tof2-tof4
    h_tof_2d[11] = tof2-tof5
    h_tof_2d[12] = tof3-tof4
    h_tof_2d[13] = tof3-tof5
    h_tof_2d[14] = tof4-tof5
  */
  
  TH1F* h_tof[6];
  TH2F* h_tof_2d[15];
  TH1F* h_amp[4];
  //Creating TOF histograms
  for(int j = 0; j < 6; j++){
    TString s(Form("tof_%d",j));
    h_tof[j] = new TH1F(s,s, 5000, -25, 25);
  }
  //Creating Amp histograms
  for(int j = 0; j < 4; j++){
    TString s(Form("amp_%d",j));
    h_amp[j] = new TH1F(s,s, 100, 1, 500);
  }
  
  //Filling Histograms for the first time(Very wide range)
  int nentries = tree->GetEntries();
  for(int i = 0; i < nentries; i++){
    tree->GetEntry(i);
    if(bit1==0 && bit2==0)h_tof[0]->Fill(t1-t2);
    if(bit1==0 && bit3==0)h_tof[1]->Fill(t1-t3);
    if(bit1==0 && bit4==0)h_tof[2]->Fill(t1-t4);
    if(bit2==0 && bit3==0)h_tof[3]->Fill(t2-t3);
    if(bit2==0 && bit4==0)h_tof[4]->Fill(t2-t4);
    if(bit3==0 && bit4==0)h_tof[5]->Fill(t3-t4);
    if(bit1==0)h_amp[0]->Fill(Amp1*1000.0);//Converting amplitude to mVolts
    if(bit2==0)h_amp[1]->Fill(Amp2*1000.0);
    if(bit3==0)h_amp[2]->Fill(Amp3*1000.0);
    if(bit4==0)h_amp[3]->Fill(Amp4*1000.0);
  }
  
  //Saving Amplitude Histograms(no further iteration needed for them)
  for(int j = 0; j < 4; j++){
    TString s(Form("/Amp_%d",j+1));
    Makeup(h_amp[j], sn+s, "Amp [mV]", "Events/1mV");
  }
  
  //Iterate the range of the TOF histograms according the the MAX and the RMS
  for(int j = 0; j < 6; j++){
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
  
  //Re-Filling TOF Histograms re-binned
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
  float b_low[6];
  float b_high[6];
  int n_bins[6];
  for(int j = 0; j < 6; j++){
    TString s(Form("tof_%d",j));
    int nb = -1;
    if(h_tof[j]->GetRMS() < 0.30){//is smaller than 90 ps?
      map_ix[mi_ctr] = j;
      mi_ctr++;
      if(h_tof[j]->Integral()>300){
	//nb = (int)(h_tof[j]->Integral())/40;
	nb = 80;
      }else{
	nb = 80;
	//nb = (int)(h_tof[j]->Integral());
      }
      float nsigma = -1.0;
      if(h_tof[j]->GetRMS() > 0.20){
	nsigma = 1.0;
      }else{
	nsigma = 4.0;
      }
      float hi = h_tof[j]->GetMean()+nsigma*h_tof[j]->GetRMS();
      float lo = h_tof[j]->GetMean()-nsigma*h_tof[j]->GetRMS();
      b_low[j] = lo;
      b_high[j] = hi;
      n_bins[j] = nb;
      //float hi = h_tof[j]->GetMean()+.5;
      //float lo = h_tof[j]->GetMean()-0.5;
      std::cout << "nb: " << nb << std::endl;
      delete h_tof[j];
      h_tof[j] = new TH1F(s,s, nb, lo, hi);
    }
  }
  
  //Create 2d TOF
  //tof0_tofX
  std::cout << "debug 0" << std::endl;
  for(int i = 1; i < 6; i++){
    TString h_2d_n(Form("t0_t%d",i));
    h_tof_2d[i-1] = new TH2F(h_2d_n, h_2d_n, (int)n_bins[0]/2, b_low[0], b_high[0], 
			     (int)n_bins[i]/2, b_low[i], b_high[i]);
  }
  //tof1_tofX
  std::cout << "debug 1" << std::endl;
  for(int i = 2; i < 6; i++){
    TString h_2d_n(Form("t1_t%d",i));
    h_tof_2d[5+(i-2)] = new TH2F(h_2d_n, h_2d_n, (int)n_bins[1]/2, b_low[1], b_high[1], 
				 (int)n_bins[i]/2, b_low[i], b_high[i]);
  }
  //tof2_tofX
  std::cout << "debug 2" << std::endl;
  for(int i = 3; i < 6; i++){
    TString h_2d_n(Form("t2_t%d",i));
    h_tof_2d[9+(i-3)] = new TH2F(h_2d_n, h_2d_n, (int)n_bins[2]/2, b_low[2], b_high[2], 
				 (int)n_bins[i]/2, b_low[i], b_high[i]);
  }
  //tof3_tofX
  std::cout << "debug 3" << std::endl;
  for(int i = 4; i < 6; i++){
    TString h_2d_n(Form("t3_t%d",i));
    h_tof_2d[12+(i-4)] = new TH2F(h_2d_n, h_2d_n, (int)n_bins[3]/2, b_low[3], b_high[3], 
				  (int)n_bins[i]/2, b_low[i], b_high[i]);
  }
  //tof4_tofX
  std::cout << "debug 4" << std::endl;
  for(int i = 5; i < 6; i++){
    TString h_2d_n(Form("t4_t%d",i));
    h_tof_2d[14] = new TH2F(h_2d_n, h_2d_n, (int)n_bins[4]/2, b_low[4], b_high[4], 
			    (int)n_bins[i]/2, b_low[i], b_high[i]);
  }
  
  std::cout << "Reiterate: " << mi_ctr << std::endl;
  
  for(int i = 0; i < nentries; i++){
    tree->GetEntry(i);
    float tof0 = t1-t2;
    float tof1 = t1-t3;
    float tof2 = t1-t4;
    float tof3 = t2-t3;
    float tof4 = t2-t4;
    float tof5 = t3-t4;
    for(int j = 0; j < mi_ctr; j++){
      //std::cout << map_ix[j] << std::endl;
      switch(map_ix[j]){
      case 0:
	//std::cout << map_ix[j] << std::endl;
	if(bit1==0 && bit2==0)h_tof[0]->Fill(tof0);
	break;
      case 1:
	//std::cout << map_ix[j] << std::endl;
	if(bit1==0 && bit3==0)h_tof[1]->Fill(tof1);
	break;
      case 2:
	//std::cout << map_ix[j] << std::endl;
	if(bit1==0 && bit4==0)h_tof[2]->Fill(tof2);
	break;
      case 3:
	//std::cout << map_ix[j] << std::endl;
	if(bit2==0 && bit3==0)h_tof[3]->Fill(tof3);
	break;
      case 4:
	//std::cout << map_ix[j] << std::endl;
	if(bit2==0 && bit4==0)h_tof[4]->Fill(tof4);
	break;
      case 5:
	//std::cout << map_ix[j] << std::endl;
	if(bit3==0 && bit4==0)h_tof[5]->Fill(tof5);
	break;
      default:
	std::cout << "DEFAULT!!" << std::endl;
	break;
      }
    }
    
    if(bit1 == 0 && bit2 == 0 && bit3==0){
      h_tof_2d[0]->Fill(tof0, tof1);
      h_tof_2d[2]->Fill(tof0, tof3);
      h_tof_2d[6]->Fill(tof1, tof3);
    }
    if(bit1 == 0 && bit2 == 0 && bit4 == 0){
      h_tof_2d[1]->Fill(tof0, tof2);
      h_tof_2d[3]->Fill(tof0, tof4);
      h_tof_2d[10]->Fill(tof2, tof4);
    }
    if(bit1 == 0 && bit2 == 0 && bit3 == 0 && bit4 == 0){
      h_tof_2d[4]->Fill(tof0, tof5);
      h_tof_2d[7]->Fill(tof1, tof4);
      h_tof_2d[9]->Fill(tof2, tof3);
    }
    if(bit1 == 0 && bit3 == 0 && bit4 ==0){
       h_tof_2d[5]->Fill(tof1, tof2);
       h_tof_2d[8]->Fill(tof1, tof5);
       h_tof_2d[11]->Fill(tof2, tof5);
    }
    if(bit2 == 0 && bit3 == 0 && bit4==0){
      h_tof_2d[12]->Fill(tof3, tof4);
      h_tof_2d[13]->Fill(tof3, tof5);
      h_tof_2d[14]->Fill(tof4, tof5);
    }
    
  }
  
  TF1* gf;
  
  //Creating FIT function and Saving TOF Histograms
  for(int i = 0; i < 6/*mi_ctr*/; i++){
    float nsigma = 2.0;
    float t_low = h_tof[i]->GetMean() - nsigma*h_tof[i]->GetRMS();
    float t_high = h_tof[i]->GetMean() + nsigma*h_tof[i]->GetRMS();
    gf = new TF1("gf","gaus(0)",t_low, t_high);
    gf->SetParNames("Norm_{fit}","#mu_{fit}", "#sigma_{fit}");
    gf->SetParameter(0,h_tof[i]->GetMaximum());
    gf->SetParameter(1,h_tof[i]->GetMean());
    gf->SetParameter(2,h_tof[i]->GetRMS());
    h_tof[i]->Fit(gf,"MWLR");
    TString s(Form("/tof_%d",i));
    Makeup(h_tof[i], sn+s, xtit[i], "Events/10ps");
  }
  
  TProfile* pro_x[15];
  TProfile* pro_y[15];
  double par0_x[15];
  double par1_x[15];
  double mean_x[15];
  double par0_y[15];
  double par1_y[15];
  double mean_y[15];
  for(int j = 0; j < 15; j++){
    pro_x[j] = h_tof_2d[j]->ProfileX("_pfx", 1, -1, "o");
    pro_y[j] = h_tof_2d[j]->ProfileY("_pfy", 1, -1, "o");
    TFitResultPtr rx = pro_x[j]->Fit("pol1", "S", "", 
				     pro_x[j]->GetMean(1)-2.0*pro_x[j]->GetRMS(1),
				     pro_x[j]->GetMean(1)+2.0*pro_x[j]->GetRMS(1));
    pro_x[j]->SetAxisRange(pro_x[j]->GetMean(2)-2.0*pro_x[j]->GetRMS(2),
			   pro_x[j]->GetMean(2)+2.0*pro_x[j]->GetRMS(2),"Y");
    
    TFitResultPtr ry = pro_y[j]->Fit("pol1", "S", "", 
				     pro_y[j]->GetMean(1)-2.0*pro_y[j]->GetRMS(1),
				     pro_y[j]->GetMean(1)+2.0*pro_y[j]->GetRMS(1));
    pro_y[j]->SetAxisRange(pro_y[j]->GetMean(2)-2.0*pro_y[j]->GetRMS(2),
			   pro_y[j]->GetMean(2)+2.0*pro_y[j]->GetRMS(2),"Y");
    mean_x[j] = pro_x[j]->GetMean(1);
    par0_x[j] = rx->Parameter(0);//retrieve the value for the parameter 0
    par1_x[j] = rx->Parameter(1);
    mean_y[j] = pro_y[j]->GetMean(1);
    par0_y[j] = ry->Parameter(0);//retrieve the value for the parameter 0
    par1_y[j] = ry->Parameter(1);
  }
  
  TH1F* h_tof_c[30];
  for(int j = 0; j < 30; j++){
    TString sn(Form("tof_%d_c_v%d",j%6,j/6));
    h_tof_c[j] = new TH1F(sn, sn, n_bins[j%6], b_low[j%6], b_high[j%6]);
  }
  for(int i = 0; i < nentries; i++){
    tree->GetEntry(i);
    float tof0 = t1-t2;
    float tof1 = t1-t3;
    float tof2 = t1-t4;
    float tof3 = t2-t3;
    float tof4 = t2-t4;
    float tof5 = t3-t4;
    if(bit1==0 && bit2==0 && bit3==0){
      h_tof_c[0]->Fill(tof0-(tof1-mean_y[0])*par1_y[0]);
      h_tof_c[1]->Fill(tof1-(tof0-mean_x[0])*par1_x[0]);
      
      h_tof_c[3]->Fill(tof3-(tof0-mean_x[2])*par1_x[2]);
      h_tof_c[12]->Fill(tof0-(tof3-mean_y[2])*par1_x[2]);
      
      h_tof_c[13]->Fill(tof1-(tof3-mean_y[6])*par1_y[6]);
      h_tof_c[9]->Fill(tof3-(tof1-mean_x[6])*par1_x[6]);
     
    }
    if(bit1==0 && bit2==0 && bit4==0){
      h_tof_c[6]->Fill(tof0-(tof2-mean_y[1])*par1_y[1]);
      h_tof_c[2]->Fill(tof2-(tof0-mean_x[1])*par1_x[1]);
      
      h_tof_c[18]->Fill(tof0-(tof4-mean_y[3])*par1_y[3]);
      h_tof_c[4]->Fill(tof4-(tof0-mean_x[3])*par1_x[3]);

      h_tof_c[20]->Fill(tof2-(tof4-mean_y[10])*par1_y[10]);
      h_tof_c[16]->Fill(tof4-(tof2-mean_x[10])*par1_x[10]);
    }
    
    if(bit1==0 && bit2==0 && bit4==0){
      h_tof_c[24]->Fill(tof0-(tof5-mean_y[4])*par1_y[4]);
      h_tof_c[5]->Fill(tof5-(tof0-mean_x[4])*par1_x[4]);
      
      h_tof_c[19]->Fill(tof1-(tof4-mean_y[7])*par1_y[7]);
      h_tof_c[10]->Fill(tof4-(tof1-mean_x[7])*par1_x[7]);
      
      h_tof_c[14]->Fill(tof2-(tof3-mean_y[9])*par1_y[9]);
      h_tof_c[15]->Fill(tof3-(tof2-mean_x[9])*par1_x[9]);
    }
    if(bit1==0 && bit3==0 && bit4==0){
      h_tof_c[7]->Fill(tof1-(tof2-mean_y[5])*par1_y[5]);
      h_tof_c[8]->Fill(tof2-(tof1-mean_x[5])*par1_x[5]);
      
      h_tof_c[25]->Fill(tof1-(tof5-mean_y[8])*par1_y[8]);
      h_tof_c[11]->Fill(tof5-(tof1-mean_x[8])*par1_x[8]);
      
      h_tof_c[26]->Fill(tof2-(tof5-mean_y[11])*par1_y[11]);
      h_tof_c[17]->Fill(tof5-(tof2-mean_x[11])*par1_x[11]);
    }
    if(bit1==0 && bit3==0 && bit4==0){
      h_tof_c[21]->Fill(tof3-(tof4-mean_y[12])*par1_y[12]);
      h_tof_c[22]->Fill(tof4-(tof3-mean_x[12])*par1_x[12]);
      
      h_tof_c[27]->Fill(tof3-(tof5-mean_y[13])*par1_y[13]);
      h_tof_c[23]->Fill(tof5-(tof3-mean_x[13])*par1_x[13]);

      h_tof_c[28]->Fill(tof4-(tof5-mean_y[14])*par1_y[14]);
      h_tof_c[29]->Fill(tof5-(tof4-mean_x[14])*par1_x[14]);
    }
  }
  //Save figures for corrected TOF
  for(int i = 0; i < 30; i++){
    float nsigma = 2.0;
    float t_low = h_tof_c[i]->GetMean() - nsigma*h_tof_c[i]->GetRMS();
    float t_high = h_tof_c[i]->GetMean() + nsigma*h_tof_c[i]->GetRMS();
    gf = new TF1("gf","gaus(0)",t_low, t_high);
    gf->SetParNames("Norm_{fit}","#mu_{fit}", "#sigma_{fit}");
    gf->SetParameter(0,h_tof_c[i]->GetMaximum());
    gf->SetParameter(1,h_tof_c[i]->GetMean());
    gf->SetParameter(2,h_tof_c[i]->GetRMS());
    h_tof_c[i]->Fit(gf,"MWLR");
    TString s(Form("/tof_%d_c_v%d",i%6,i/6));
    Makeup(h_tof_c[i], sn+s, xtit[i%6], "Events/10ps");
  }
  
  TFile* tmp = new TFile(sn+"/tmp.root", "RECREATE");
 
  for(int j = 0; j < 6; j++){
    h_tof[j]->Write();
  }
  for(int j = 0; j < 15; j++){
    h_tof_2d[j]->Write();
    pro_x[j]->Write();
    pro_y[j]->Write();
  }
  for(int j = 0; j < 30; j++){
    h_tof_c[j]->Write();
  }
  return 0;
  
}
