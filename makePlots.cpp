#ifdef __MAKECINT__
#pragma link C++ class vector<vector<float> >+;
#endif

#include <iostream>
#include <fstream> 
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TColor.h"

float GetLastXBin(TH2F * histo);


int main (int argc, char **argv)
{
  // Get the tree

  TFile *f;

  if (argc >= 3)
    {
      f = new TFile(argv[1]);

      std::cout << ">> Opening file " << argv[1] << " ......" << std::endl;
      if (!f->IsOpen())
	{			// terminate if the file can't be opened
	  std::cerr << "!! File open error:" << argv[1] << std::endl;
	  return 1;
	}
    }
  else
    {				// terminate if there is no input file or more than 1 input file
      std::cerr << "!! No input file" << std::endl;
      return 1;
    }

  TTree *tree = (TTree*)f->Get("tree");

  // get the variables from the ntuple
  float t1gausroot = 0;
  float t2gausroot = 0;
  float t3gausroot = 0;
  float t4gausroot = 0;
  float t5gausroot = 0;
  float t6gausroot = 0;
  float t7gausroot = 0;
  float t8gausroot = 0;
  float ch1Amp = 0;
  float ch2Amp = 0;
  float ch3Amp = 0;
  float ch4Amp = 0;
  float ch5Amp = 0;
  float ch6Amp = 0;
  float ch7Amp = 0;
  float ch8Amp = 0;
  unsigned int ch1QualityBit = 0;
  unsigned int ch2QualityBit = 0;
  unsigned int ch3QualityBit = 0;
  unsigned int ch4QualityBit = 0;
  unsigned int ch5QualityBit = 0;
  unsigned int ch6QualityBit = 0;
  unsigned int ch7QualityBit = 0;
  unsigned int ch8QualityBit = 0;

  tree->SetBranchAddress("t1gausroot",&t1gausroot);
  tree->SetBranchAddress("t2gausroot",&t2gausroot);
  tree->SetBranchAddress("t3gausroot",&t3gausroot);
  tree->SetBranchAddress("t4gausroot",&t4gausroot);
  tree->SetBranchAddress("t5gausroot",&t5gausroot);
  tree->SetBranchAddress("t6gausroot",&t6gausroot);
  tree->SetBranchAddress("t7gausroot",&t7gausroot);
  tree->SetBranchAddress("t8gausroot",&t8gausroot);
  tree->SetBranchAddress("ch1Amp",&ch1Amp);
  tree->SetBranchAddress("ch2Amp",&ch2Amp);
  tree->SetBranchAddress("ch3Amp",&ch3Amp);
  tree->SetBranchAddress("ch4Amp",&ch4Amp);
  tree->SetBranchAddress("ch5Amp",&ch5Amp);
  tree->SetBranchAddress("ch6Amp",&ch6Amp);
  tree->SetBranchAddress("ch7Amp",&ch7Amp);
  tree->SetBranchAddress("ch8Amp",&ch8Amp);
  tree->SetBranchAddress("ch1QualityBit",&ch1QualityBit);
  tree->SetBranchAddress("ch2QualityBit",&ch2QualityBit);
  tree->SetBranchAddress("ch3QualityBit",&ch3QualityBit);
  tree->SetBranchAddress("ch4QualityBit",&ch4QualityBit);
  tree->SetBranchAddress("ch5QualityBit",&ch5QualityBit);
  tree->SetBranchAddress("ch6QualityBit",&ch6QualityBit);
  tree->SetBranchAddress("ch7QualityBit",&ch7QualityBit);
  tree->SetBranchAddress("ch8QualityBit",&ch8QualityBit);

  //create histograms
  TH1F *dt[6];
  TH1F *resolution[6];
  TH1F *dt_slopeCorrected[6];
  TH1F *dt_shower[6];
  TH1F *dt_shower_AmpCut[6];

  TH1F *amp_all[6];
  TH1F *amp_good[6];
  TProfile *hProfCorr[6];
  TH2F *h2DCorr[6];
  TH2F *h2DCorr_Fit[6];
  TH2F *h2DCorr_slopeCorrected[6];

  for(Int_t i = 0; i < 6; i++)
    {
      dt[i] = new TH1F(Form("dt_%d",i), Form("dt_%d",i), 1500,-3.0,3.0);
      dt_slopeCorrected[i] = new TH1F(Form("dt_slopeCorrected_%d",i), Form("dt_slopeCorrected_%d",i), 700,-3.0,3.0);
      dt_shower[i] = new TH1F(Form("dt_shower_%d",i), Form("dt_shower_%d",i), 3000,-3.0,3.0);
      dt_shower_AmpCut[i] = new TH1F(Form("dt_shower_AmpCut%d",i), Form("dt_shower_AmpCut%d",i), 3000,-3.0,3.0);
      amp_all[i] = new TH1F(Form("amp_all_%d",i), Form("amp_all_%d",i), 200,0,0.55);
      amp_good[i] = new TH1F(Form("amp_good_%d",i), Form("amp_good_%d",i), 200,0,0.55);

      h2DCorr_Fit[i] = new TH2F(Form("dt_vs_amp_2D_fit_%d",i), Form("dt_vs_amp_2D_fit_%d",i), 10, 0,0.5, 300, -3.0, 3.0);
      resolution[i] = new TH1F(Form("resolution_%d",i), Form("resolution_%d",i), 10, 0, 0.05);

      hProfCorr[i] = new TProfile(Form("dt_vs_amp_%d",i), Form("dt_vs_amp_%d",i), 100, 0,0.5, -3.0, 3.0);
      h2DCorr[i] = new TH2F(Form("dt_vs_amp_2D_%d",i), Form("dt_vs_amp_2D_%d",i), 100, 0,0.5, 300, -3.0, 3.0);
      h2DCorr_slopeCorrected[i] = new TH2F(Form("dt_vs_amp_corrected_%d",i), Form("dt_vs_amp_corrected_%d",i), 100, 0,0.5, 300, -3.0, 3.0);
    }

  std::string names[6] = {"CH1-CH2", "CH1-CH3", "CH2-CH3", "CH5-CH6", "CH5-CH7", "CH6-CH7"};
  std::string Amps[6] = {"ch2Amp", "ch3Amp", "ch3Amp", "ch6Amp", "ch7Amp", "ch7Amp"};

  std::string amps_xaxis[6] = {"CH1", "CH2", "CH3", "CH5", "CH6", "CH7"};

  //read all entries and fill the histograms
  Long64_t nentries = tree->GetEntries();

  std::cout<<"Number of events in Physics Sample: "<<nentries<<std::endl;
  
  for (Long64_t iEntry=0;iEntry<nentries;iEntry++) 
    {
      tree->GetEntry(iEntry);
      
      // any good pulses
      if( ch1QualityBit==0 && ch2QualityBit==0 )
	{
	  dt[0]->Fill(t1gausroot-t2gausroot);
	  hProfCorr[0]->Fill(ch2Amp,t1gausroot-t2gausroot);
	  h2DCorr[0]->Fill(ch2Amp,t1gausroot-t2gausroot);
	  h2DCorr_Fit[0]->Fill(ch2Amp,t1gausroot-t2gausroot);
	}

      if( ch1QualityBit==0 && ch3QualityBit==0 )
	{
	  dt[1]->Fill(t1gausroot-t3gausroot);
	  hProfCorr[1]->Fill(ch3Amp,t1gausroot-t3gausroot);
	  h2DCorr[1]->Fill(ch3Amp,t1gausroot-t3gausroot);
	  h2DCorr_Fit[1]->Fill(ch3Amp,t1gausroot-t3gausroot);
	}

      if( ch2QualityBit==0 && ch3QualityBit==0 )
	{
	  dt[2]->Fill(t2gausroot-t3gausroot);
	  hProfCorr[2]->Fill(ch3Amp,t2gausroot-t3gausroot);
	  h2DCorr[2]->Fill(ch3Amp,t2gausroot-t3gausroot);
	  h2DCorr_Fit[2]->Fill(ch3Amp,t2gausroot-t3gausroot);
	}

      if( ch5QualityBit==0 && ch6QualityBit==0 )
	{
	  dt[3]->Fill(t5gausroot-t6gausroot);
	  hProfCorr[3]->Fill(ch6Amp,t5gausroot-t6gausroot);
	  h2DCorr[3]->Fill(ch6Amp,t5gausroot-t6gausroot);
	  h2DCorr_Fit[3]->Fill(ch6Amp,t5gausroot-t6gausroot);
	}

      if( ch5QualityBit==0 && ch7QualityBit==0 )
	{
	  dt[4]->Fill(t5gausroot-t7gausroot);
	  hProfCorr[4]->Fill(ch7Amp,t5gausroot-t7gausroot);
	  h2DCorr[4]->Fill(ch7Amp,t5gausroot-t7gausroot);
	  h2DCorr_Fit[4]->Fill(ch7Amp,t5gausroot-t7gausroot);
	}

      if( ch6QualityBit==0 && ch7QualityBit==0 )
	{
	  dt[5]->Fill(t6gausroot-t7gausroot);
	  hProfCorr[5]->Fill(ch7Amp,t6gausroot-t7gausroot);
	  h2DCorr[5]->Fill(ch7Amp,t6gausroot-t7gausroot);
	  h2DCorr_Fit[5]->Fill(ch7Amp,t6gausroot-t7gausroot);
	}
      
      // amplitude plots
      amp_all[0]->Fill(ch1Amp);
      amp_all[1]->Fill(ch2Amp);
      amp_all[2]->Fill(ch3Amp);
      amp_all[3]->Fill(ch5Amp);
      amp_all[4]->Fill(ch6Amp);
      amp_all[5]->Fill(ch7Amp);

      if( ch1QualityBit==0 )
	amp_good[0]->Fill(ch1Amp);
      if( ch2QualityBit==0 )
      	amp_good[1]->Fill(ch2Amp);
      if( ch3QualityBit==0 )
      	amp_good[2]->Fill(ch3Amp);
      if( ch5QualityBit==0 )
      	amp_good[3]->Fill(ch5Amp);
      if( ch6QualityBit==0 )
      	amp_good[4]->Fill(ch6Amp);
      if( ch7QualityBit==0 )
      	amp_good[5]->Fill(ch7Amp);

      // good pulses with shower requirement
      if( ch1QualityBit==0 && ch3QualityBit==0 && ch6Amp>0.49)
	dt_shower[1]->Fill(t1gausroot-t3gausroot);

      if( ch1QualityBit==0 && ch2QualityBit==0 && ch6Amp>0.49) //cherenkov selection
	dt_shower[0]->Fill(t1gausroot-t2gausroot);

      if( ch5QualityBit==0 && ch6QualityBit==0 && ch6Amp>0.49)  //cherenkov selection
	dt_shower[3]->Fill(t5gausroot-t6gausroot);

      if( ch1QualityBit==0 && ch3QualityBit==0 && ch6Amp>0.49)
	if(ch3Amp>0.15)
	  dt_shower_AmpCut[1]->Fill(t1gausroot-t3gausroot);

    }
  
  // Touch up and save
  TCanvas * c = new TCanvas("c","c",600,600);
  for(Int_t i = 0; i < 6; i++)
    {
      c->cd();
      
      // TOFs
      gPad->SetLogy(0);
      dt[i]->SetAxisRange(dt[i]->GetMean()-0.5*fabs(dt[i]->GetMean()),dt[i]->GetMean()+0.5*fabs(dt[i]->GetMean()),"X");
      dt[i]->GetXaxis()->SetTitle(Form("%s [ns]",names[i].c_str()));
      dt[i]->SetTitle(names[i].c_str());
      dt[i]->GetYaxis()->SetTitle("Number of Events");
      dt[i]->Draw();
      dt[i]->Fit("gaus");
      TVirtualFitter * fitter = TVirtualFitter::GetFitter();
      TPaveText *pt = new TPaveText(.15,.75,.35,.85,"NDC");
      pt->AddText(Form("#sigma=%.2f psec",1000*fitter->GetParameter(2)));
      pt->SetFillColor(0);
      pt->Draw();

      c->SaveAs( Form("TOF_%s_Run%s.pdf", names[i].c_str(), argv[2]) );

      // amplitudes
      amp_good[i]->SetLineWidth(2);
      amp_good[i]->SetLineColor(1);
      amp_all[i]->SetLineWidth(3);
      amp_all[i]->SetLineColor(1);
      amp_all[i]->SetLineStyle(2);
      amp_all[i]->SetStats(0);
      amp_all[i]->Draw();
      amp_all[i]->GetXaxis()->SetTitle(Form("%s  Amp [mV]", amps_xaxis[i].c_str()));
      amp_all[i]->GetYaxis()->SetTitle("Number of Events");
      amp_all[i]->Draw();
      amp_good[i]->SetLineColor(1);
      amp_good[i]->Draw("same");
      TLegend *leg = new TLegend(0.5,0.8,0.8,0.9);
      leg->AddEntry(amp_all[i],"All pulses");
      leg->AddEntry(amp_good[i],"Pass all quality bits");
      leg->SetFillColor(0);
      leg->Draw();
      gPad->SetLogy();
  

      c->SaveAs( Form("Amplitudes_%s_Run%s.root", amps_xaxis[i].c_str(), argv[2]) );
    }
  
  gPad->SetLogy(0);
  // shower requirement
  for(Int_t i = 0; i < 4; i++)
    {
      c->cd();
      
      if(i==2) continue;

      dt_shower[i]->SetAxisRange(dt_shower[i]->GetMean()-0.2*fabs(dt_shower[i]->GetMean()),
				 dt_shower[i]->GetMean()+0.2*fabs(dt_shower[i]->GetMean()),"X");
      dt_shower[i]->GetXaxis()->SetTitle(Form("%s [ns]",names[i].c_str()));
      dt_shower[i]->SetTitle(names[i].c_str());
      dt_shower[i]->GetYaxis()->SetTitle("Number of Events");
      dt_shower[i]->Draw();
      dt_shower[i]->SetLineColor(1);
      dt_shower[i]->SetStats(0);
      dt_shower[i]->SetTitle("");
      dt_shower[i]->Draw();
      dt_shower[i]->Fit("gaus");
      dt_shower[i]->GetFunction("gaus")->SetLineColor(1);
      dt_shower[i]->GetFunction("gaus")->SetLineWidth(4);


      TVirtualFitter * fitter = TVirtualFitter::GetFitter();
      TPaveText *pt = new TPaveText(.15,.75,.35,.85,"NDC");
      pt->AddText(Form("#sigma=%.2f psec",1000*fitter->GetParameter(2)));
      pt->SetFillColor(0);
      pt->Draw();

      c->SaveAs( Form("TOF_Shower_%s_Run%s.root", names[i].c_str(), argv[2]) );
    }

  c->cd();
  dt_shower_AmpCut[1]->SetAxisRange(dt_shower_AmpCut[1]->GetMean()-0.5*fabs(dt_shower_AmpCut[1]->GetMean()),
				    dt_shower_AmpCut[1]->GetMean()+0.5*fabs(dt_shower_AmpCut[1]->GetMean()),"X");
  dt_shower_AmpCut[1]->GetXaxis()->SetTitle(Form("%s [ns]",names[1].c_str()));
  dt_shower_AmpCut[1]->SetTitle(names[1].c_str());
  dt_shower_AmpCut[1]->GetYaxis()->SetTitle("Number of Events");
  dt_shower_AmpCut[1]->Draw();
  dt_shower_AmpCut[1]->Fit("gaus");
  TVirtualFitter * fitter3 = TVirtualFitter::GetFitter();
  TPaveText *pt2 = new TPaveText(.15,.75,.35,.85,"NDC");
  pt2->AddText(Form("#sigma=%.2f psec",1000*fitter3->GetParameter(2)));
  pt2->SetFillColor(0);
  pt2->Draw();

  c->SaveAs( Form("TOF_Shower_AmpCut_%s_Run%s.pdf", names[1].c_str(), argv[2]) );

  // slopes
  TF1 * slopes[6];
  for(Int_t i = 0; i < 6; i++)
    {
      slopes[i] = new TF1(Form("slopes_%d",i),"pol1", 0, h2DCorr[i]->GetNbinsX());
      slopes[i]->SetLineColor(1);
      slopes[i]->SetLineWidth(3);
    }

  for(Int_t i = 0; i < 6; i++)
    {
      c->cd();

      hProfCorr[i]->Fit(slopes[i],"Q","", 0, h2DCorr[i]->GetNbinsX());

      UInt_t Number = 3;
      Double_t Red[3]   = { 0.15, 0.6, 0.9};
      Double_t Green[3] = { 0.15, 0.6, 0.9};
      Double_t Blue[3]  = { 0.15, 0.6, 0.9};
      Double_t Stops[3] = { 0.00, 0.25, 1.0};
      
      Int_t nb=50;
      TColor::CreateGradientColorTable(Number,Stops,Red,Green,Blue,nb);
      
      h2DCorr[i]->Draw("colz");
      h2DCorr[i]->SetLineColor(1);
      h2DCorr[i]->SetAxisRange(h2DCorr[i]->GetMean(2) - 0.3*fabs(h2DCorr[i]->GetMean(2)),
			       h2DCorr[i]->GetMean(2) + 0.3*fabs(h2DCorr[i]->GetMean(2)),"Y");
      h2DCorr[i]->SetAxisRange(0, GetLastXBin(h2DCorr[i]) + 0.05,"X");
      h2DCorr[i]->GetYaxis()->SetTitle(Form("Time Difference %s", names[i].c_str()));
      h2DCorr[i]->GetXaxis()->SetTitle(Form("%s", Amps[i].c_str()));
      h2DCorr[i]->SetTitle("");
      h2DCorr[i]->SetStats(0);
      hProfCorr[i]->SetLineColor(1);
      hProfCorr[i]->SetMarkerStyle(8);

      hProfCorr[i]->Draw("same");
      
      c->SaveAs( Form("TOF_vs_Amp_%s_%s_Run%s.pdf", names[i].c_str(), Amps[i].c_str(), argv[2]) );

      
      // resolution vs amplitude plots
      for(int j = 0; j < h2DCorr_Fit[i]->GetNbinsX(); j++)
	{
	  TH1D *py = h2DCorr_Fit[i]->ProjectionY("py", j, j + 1 ); 
	  py->Fit("gaus");
	  TVirtualFitter * fitter22 = TVirtualFitter::GetFitter();
	  resolution[i] -> SetBinContent(j, 1000*fitter22->GetParameter(2));
	  if (1000*fitter22->GetParameter(2) > 100. ) 
	    resolution[i] -> SetBinContent(j, 0.0);
	}
      resolution[i]->GetXaxis()->SetTitle(Form("%s", Amps[i].c_str()));
      resolution[i]->GetYaxis()->SetTitle(Form("Resolution in %s [psec]", names[i].c_str()));
      resolution[i]->Draw();
	    
      c->SaveAs( Form("Resolution_vs_Amp_%s_%s_Run%s.pdf", names[i].c_str(), Amps[i].c_str(), argv[2]) );
    }

  // correction for slopes
  std::cout<<"Second round: slope correction "<<nentries<<std::endl;
  
  for (Long64_t iEntry=0;iEntry<nentries;iEntry++) 
    {
      tree->GetEntry(iEntry);
      
      // any good pulses
      if( ch1QualityBit==0 && ch3QualityBit==0 && ch6Amp>0.49)
	{
	  float a = slopes[1]->GetParameter(0);
	  float b = slopes[1]->GetParameter(1);
	  dt_slopeCorrected[1]->Fill(t1gausroot-t3gausroot + (a+b*0.02) - (a+b*ch3Amp));
	  h2DCorr_slopeCorrected[1]->Fill(ch3Amp,t1gausroot-t3gausroot + (a+b*0.02) - (a+b*ch3Amp));
	}
    }
  
  c->cd();
  dt_slopeCorrected[1]->SetAxisRange(dt_slopeCorrected[1]->GetMean()-0.2*fabs(dt_slopeCorrected[1]->GetMean()),
				     dt_slopeCorrected[1]->GetMean()+0.2*fabs(dt_slopeCorrected[1]->GetMean()),"X");
  dt_slopeCorrected[1]->GetXaxis()->SetTitle(Form("%s [ns]",names[1].c_str()));
  dt_slopeCorrected[1]->SetTitle(names[1].c_str());
  dt_slopeCorrected[1]->GetYaxis()->SetTitle("Number of Events");
  dt_slopeCorrected[1]->SetLineColor(1);
  // dt_slopeCorrected[1]->SetStats(0);
  dt_slopeCorrected[1]->SetTitle("");
  dt_slopeCorrected[1]->Draw();
  dt_slopeCorrected[1]->Fit("gaus");
  dt_slopeCorrected[1]->GetFunction("gaus")->SetLineColor(1);
  dt_slopeCorrected[1]->GetFunction("gaus")->SetLineWidth(4);
  TVirtualFitter * fitter2 = TVirtualFitter::GetFitter();
  TPaveText *pt1 = new TPaveText(.15,.75,.35,.85,"NDC");
  pt1->AddText(Form("#sigma=%.2f psec",1000*fitter2->GetParameter(2)));
  pt1->SetFillColor(0);
  pt1->Draw();
  
  c->SaveAs( Form("TOF_SlopeCorrected_%s_Run%s.pdf", names[1].c_str(), argv[2]) );
  
  // 2D
  c->cd();
  
  h2DCorr_slopeCorrected[1]->Draw("colz");
  h2DCorr_slopeCorrected[1]->SetAxisRange(h2DCorr_slopeCorrected[1]->GetMean(2) - 0.7*fabs(h2DCorr_slopeCorrected[1]->GetMean(2)),
  			   h2DCorr_slopeCorrected[1]->GetMean(2) + 0.7*fabs(h2DCorr_slopeCorrected[1]->GetMean(2)),"Y");
  h2DCorr_slopeCorrected[1]->SetAxisRange(0, GetLastXBin(h2DCorr_slopeCorrected[1]) + 0.05,"X");
  h2DCorr_slopeCorrected[1]->GetYaxis()->SetTitle(Form("Time Difference %s", names[1].c_str()));
  h2DCorr_slopeCorrected[1]->GetXaxis()->SetTitle(Form("%s", Amps[1].c_str()));
  
  c->SaveAs( Form("TOF_vs_Amp_Corrected_%s_%s_Run%s.pdf", names[1].c_str(), Amps[1].c_str(), argv[2]) );
 
}

float GetLastXBin(TH2F * histo)
{
  float entry = 0.;
  int lastbin = histo->GetNbinsX();

  for (int i=histo->GetNbinsX(); i>0; i--)
    for (int j=0; j<histo->GetNbinsY(); j++)
      if(histo->GetBinContent(i,j) > 0) {lastbin = i; goto breakpoint;}

 breakpoint:
  return histo->GetXaxis()->GetBinLowEdge(lastbin);
}



