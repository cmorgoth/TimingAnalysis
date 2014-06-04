#include <iostream>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"

#define binsize 1

enum PulseQuality {
  kNegativePolarity       = 0x000001, // bit 0
  kSuddenJump             = 0x000002, // bit 1
  kFlatTop                = 0x000004, // bit 2
  kSecondPulse            = 0x000008, // bit 3
  kNoPulse                = 0x000010, // bit 4
  kLargeNegativeAmplitude = 0x000020, // bit 5
  kSaturated              = 0x000040  // bit 6
};

float find_CFD(TH1F * pulse, int first_bin, int last_bin, int minBin);
int FindMin( int n, float *a);
int FindMax( int n, float *a);
int FindRisingEdge( int n, int binMax, float *a);
int FindFirstPulsePeak( int n, float *a);
unsigned int CheckPulseQuality( int binMin, int binMax, float *a, float minPulse);
float ChannelIntegral(float *a, int peak);

float LinearFit_Baseline(TH1F * pulse, const int index_min, const int range);
float LinearFit_Intercept(TH1F * pulse, const float base, const int index_first, const int index_last);
float GausFit_MeanTime(TH1F * pulse, const int index_first, const int index_last);
float FitRisingEdge(TH1F* pulse, int nbinsL, int nbinsH);
float FitRisingEdge(TH1F* pulse, float Max);
void FitRisingEdge(TH1F* pulse, float Max, float &THM, float &t0, float baseline);
void FitFullPulse(TH1F* pulse, float &par0, float &par1, float &par2);

const int Nsamples = 1024;

int main (int argc, char **argv)
{
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

  float Channel1Voltages_[Nsamples];
  float Channel2Voltages_[Nsamples];
  float Channel3Voltages_[Nsamples];
  float Channel4Voltages_[Nsamples];
  float Channel5Voltages_[Nsamples];
  float Channel6Voltages_[Nsamples];
  float Channel7Voltages_[Nsamples];
  float Channel8Voltages_[Nsamples];

  bool convert2Volts = true;

  TTree* t1 = (TTree*)f->Get("p");      // andriy's converter
  if (t1) 
    {
      t1->SetBranchAddress("c1",Channel1Voltages_);
      t1->SetBranchAddress("c2",Channel2Voltages_);
      t1->SetBranchAddress("c3",Channel3Voltages_);
      t1->SetBranchAddress("c4",Channel4Voltages_);
      t1->SetBranchAddress("c5",Channel5Voltages_);
      t1->SetBranchAddress("c6",Channel6Voltages_);
      t1->SetBranchAddress("c7",Channel7Voltages_);
      t1->SetBranchAddress("c8",Channel8Voltages_);
    }
  
  if (!t1) 
    {
      t1 = (TTree*)f->Get("T");   // artur's converter
      convert2Volts = false;

      t1->SetBranchAddress("c1",Channel1Voltages_);
      t1->SetBranchAddress("c2",Channel2Voltages_);
      t1->SetBranchAddress("c3",Channel3Voltages_);
      t1->SetBranchAddress("c4",Channel4Voltages_);
    }

  //create two histograms
  TH1F *CH1pulse   = new TH1F("CH1pulse","CH1pulse",Nsamples,0,Nsamples);
  TH1F *CH2pulse   = new TH1F("CH2pulse","CH2pulse",Nsamples,0,Nsamples);
  TH1F *CH3pulse   = new TH1F("CH3pulse","CH3pulse",Nsamples,0,Nsamples);
  TH1F *CH4pulse   = new TH1F("CH4pulse","CH4pulse",Nsamples,0,Nsamples);
  
  TH1F *CH1Amp     = new TH1F("CH1Amp","CH1Amp",40,-0.3,-0.6);
  TH1F *CH2Amp     = new TH1F("CH2Amp","CH2Amp",40,-0.3,-0.6);

  TH1F *GausPeak_CH12_dt   = new TH1F("GausPeak_CH12_dt","GausPeak_CH12_dt; t1-t2 [ns]; Events",4000,-4,4);
  TH1F *GausPeak_CH34_dt   = new TH1F("GausPeak_CH34_dt","GausPeak_CH34_dt; t3-t4 [ns]; Events",4000,-4,4);
  TH1F *GausPeak_TOF_CH13   = new TH1F("GausPeak_TOF_CH13","GausPeak_TOF_CH13; t1-t3 [ns]; Events",4000,-4,4);
  TH1F *GausPeak_TOF_CH14   = new TH1F("GausPeak_TOF_CH14","GausPeak_TOF_CH14; t1-t4 [ns]; Events",4000,-40,4);

  // Create the output file with a TTree
  TFile* fout;
  if(strncmp(argv[2], "same", 100) == 0){
    std::string fn(argv[1]);
    int pf = fn.find(".root");
    int pi = fn.rfind("/")+1;
    fn = "AnaFiles/" + fn.substr(pi, pf-pi) + "_ana.root";
    std::cout << "fname: " << fn << std::endl;
    //return 0;
    fout = new TFile(fn.c_str(),"recreate");
  }else{
    fout = new TFile(argv[2],"recreate");
  }
  TTree* treeOut = new TTree("tree","tree");
  
  unsigned int eventNumber = 0;
  float ch1Time_gausfitroot = 0;
  float ch2Time_gausfitroot = 0;
  float ch3Time_gausfitroot = 0;
  float ch4Time_gausfitroot = 0;
  
  float ch1Amp = 0;
  float ch2Amp = 0;
  float ch3Amp = 0;
  float ch4Amp = 0;

  float ch1THM = 0;//Time at half the Maximum
  float ch2THM = 0;//Time at half the Maximum
  float ch3THM = 0;//Time at half the Maximum
  float ch4THM = 0;//Time at half the Maximum

  float ch1_TFF = 0.0;
  float ch2_TFF = 0.0;
  float ch3_TFF = 0.0;
  float ch4_TFF = 0.0;

  float ch1_TFF_v2 = 0.0;
  float ch2_TFF_v2 = 0.0;
  float ch3_TFF_v2 = 0.0;
  float ch4_TFF_v2 = 0.0;
  
  float ch1BL = 0.0;
  float ch2BL = 0.0;
  float ch3BL = 0.0;
  float ch4BL = 0.0;

  float ch1_AFF = 0.0;
  float ch2_AFF = 0.0;
  float ch3_AFF = 0.0;
  float ch4_AFF = 0.0;

  float ch1Int = 0;
  float ch2Int = 0;
  float ch3Int = 0;
  float ch4Int = 0;
  
  unsigned int ch1QualityBit = 0;
  unsigned int ch2QualityBit = 0;
  unsigned int ch3QualityBit = 0;
  unsigned int ch4QualityBit = 0;
  
  float ch1chisq = -1;
  float ch2chisq = -1;
  float ch3chisq = -1;
  float ch4chisq = -1;
  
  treeOut->Branch("event",&eventNumber,"event/i");
  treeOut->Branch("t1gausroot",&ch1Time_gausfitroot,"t1gausroot/F");
  treeOut->Branch("t2gausroot",&ch2Time_gausfitroot,"t2gausroot/F");
  treeOut->Branch("t3gausroot",&ch3Time_gausfitroot,"t3gausroot/F");
  treeOut->Branch("t4gausroot",&ch4Time_gausfitroot,"t4gausroot/F");
  
  treeOut->Branch("ch1Amp",&ch1Amp,"ch1Amp/F");
  treeOut->Branch("ch2Amp",&ch2Amp,"ch2Amp/F");
  treeOut->Branch("ch3Amp",&ch3Amp,"ch3Amp/F");
  treeOut->Branch("ch4Amp",&ch4Amp,"ch4Amp/F");

  treeOut->Branch("ch1THM",&ch1THM,"ch1THM/F");
  treeOut->Branch("ch2THM",&ch2THM,"ch2THM/F");
  treeOut->Branch("ch3THM",&ch3THM,"ch3THM/F");
  treeOut->Branch("ch4THM",&ch4THM,"ch4THM/F");

  treeOut->Branch("ch1BL",&ch1BL,"ch1BL/F");
  treeOut->Branch("ch2BL",&ch2BL,"ch2BL/F");
  treeOut->Branch("ch3BL",&ch3BL,"ch3BL/F");
  treeOut->Branch("ch4BL",&ch4BL,"ch4BL/F");
  
  treeOut->Branch("ch1_TFF", &ch1_TFF, "ch1_TFF/F");
  treeOut->Branch("ch2_TFF", &ch2_TFF, "ch2_TFF/F");
  treeOut->Branch("ch3_TFF", &ch3_TFF, "ch3_TFF/F");
  treeOut->Branch("ch4_TFF", &ch4_TFF, "ch4_TFF/F");
  
  treeOut->Branch("ch1_TFF_v2", &ch1_TFF_v2, "ch1_TFF_v2/F");
  treeOut->Branch("ch2_TFF_v2", &ch2_TFF_v2, "ch2_TFF_v2/F");
  treeOut->Branch("ch3_TFF_v2", &ch3_TFF_v2, "ch3_TFF_v2/F");
  treeOut->Branch("ch4_TFF_v2", &ch4_TFF_v2, "ch4_TFF_v2/F");

  
  treeOut->Branch("ch1_AFF", &ch1_AFF, "ch1_AFF/F");
  treeOut->Branch("ch2_AFF", &ch2_AFF, "ch2_AFF/F");
  treeOut->Branch("ch3_AFF", &ch3_AFF, "ch3_AFF/F");
  treeOut->Branch("ch4_AFF", &ch4_AFF, "ch4_AFF/F");
  
  treeOut->Branch("ch1QualityBit",&ch1QualityBit,"ch1QualityBit/i");
  treeOut->Branch("ch2QualityBit",&ch2QualityBit,"ch2QualityBit/i");
  treeOut->Branch("ch3QualityBit",&ch3QualityBit,"ch3QualityBit/i");
  treeOut->Branch("ch4QualityBit",&ch4QualityBit,"ch4QualityBit/i");
  
  treeOut->Branch("ch1Int",&ch1Int,"ch1Int/F");
  treeOut->Branch("ch2Int",&ch2Int,"ch2Int/F");
  treeOut->Branch("ch3Int",&ch3Int,"ch3Int/F");
  treeOut->Branch("ch4Int",&ch4Int,"ch4Int/F");
  
  treeOut->Branch("ch1chisq",&ch1chisq,"ch1chisq/F");
  treeOut->Branch("ch2chisq",&ch2chisq,"ch2chisq/F");
  treeOut->Branch("ch3chisq",&ch3chisq,"ch3chisq/F");
  treeOut->Branch("ch4chisq",&ch4chisq,"ch4chisq/F");
  

  //read all entries and fill the histograms
  Long64_t nentries = t1->GetEntries();
  //Long64_t nentries = 10;

  for (Long64_t iEntry=0;iEntry<nentries;iEntry++) 
    {
      if(iEntry%100==0) std::cout<<"Processing Event: "<<iEntry<<" out of: "<<nentries<<std::endl;
      
      t1->GetEntry(iEntry);
      eventNumber = iEntry+1;

      //////////////////////
      // convert to Volts
      //////////////////////
      if(convert2Volts){
	for (int ii=0;ii<Nsamples;ii++){
	  Channel1Voltages_[ii] = 0.001*Channel1Voltages_[ii];
	  Channel2Voltages_[ii] = 0.001*Channel2Voltages_[ii];
	  Channel3Voltages_[ii] = 0.001*Channel3Voltages_[ii];
	  Channel4Voltages_[ii] = 0.001*Channel4Voltages_[ii];
	}
      }
      //////////////////////
      // end convert to Volts
      //////////////////////
      
      // Find Min of the Channel data (Voltage)
      int index_min1 = FindMin (Nsamples, Channel1Voltages_);	// return index of the min
      int index_min2 = FindMin (Nsamples, Channel2Voltages_);	// return index of the min
      int index_min3 = FindMin (Nsamples, Channel3Voltages_);	// return index of the min
      int index_min4 = FindMin (Nsamples, Channel4Voltages_);	// return index of the min
      // Set histograms bins and errors
      for (int ii=0;ii<Nsamples;ii++)
	{
	  CH1pulse->SetBinContent(ii+1,Channel1Voltages_[ii]);
	  CH2pulse->SetBinContent(ii+1,Channel2Voltages_[ii]);
	  CH3pulse->SetBinContent(ii+1,Channel3Voltages_[ii]);
	  CH4pulse->SetBinContent(ii+1,Channel4Voltages_[ii]);
	  
	  CH1pulse->SetBinError(ii+1, 0.05*Channel1Voltages_[index_min1]);
	  CH2pulse->SetBinError(ii+1, 0.05*Channel2Voltages_[index_min2]);
	  CH3pulse->SetBinError(ii+1, 0.05*Channel3Voltages_[index_min3]);
	  CH4pulse->SetBinError(ii+1, 0.05*Channel4Voltages_[index_min4]);
	}     

      // Find Max of the Channel data (Voltage)
      int index_max1 = FindMax (Nsamples, Channel1Voltages_);	// return index of the max
      int index_max2 = FindMax (Nsamples, Channel2Voltages_);	// return index of the max
      int index_max3 = FindMax (Nsamples, Channel3Voltages_);	// return index of the max
      int index_max4 = FindMax (Nsamples, Channel4Voltages_);	// return index of the max
      
      
      //////////////////
      // Find the rising edge on CH1
      int fbin1 = FindRisingEdge(Nsamples, index_min1, Channel1Voltages_);
      int fbin2 = FindRisingEdge(Nsamples, index_min2, Channel2Voltages_);
      int fbin3 = FindRisingEdge(Nsamples, index_min3, Channel3Voltages_);
      int fbin4 = FindRisingEdge(Nsamples, index_min4, Channel4Voltages_);
      
      // Find the quality of the pulse
      ch1QualityBit = CheckPulseQuality ( index_min1, index_max1, Channel1Voltages_, 0.02);
      ch2QualityBit = CheckPulseQuality ( index_min2, index_max2, Channel2Voltages_, 0.02);
      ch3QualityBit = CheckPulseQuality ( index_min3, index_max3, Channel3Voltages_, 0.02);
      ch4QualityBit = CheckPulseQuality ( index_min4, index_max4, Channel4Voltages_, 0.02);
      
      // For the first version of Photonis data the pulse has many peaks --> find the first one
      int index_firstPulse1 = FindFirstPulsePeak(Nsamples, Channel1Voltages_); // this is useful ONLY for Photonis MCP
      int index_firstPulse2 = FindFirstPulsePeak(Nsamples, Channel2Voltages_); // this is useful ONLY for Photonis MCP
      int index_firstPulse3 = FindFirstPulsePeak(Nsamples, Channel3Voltages_); // this is useful ONLY for Photonis MCP
      int index_firstPulse4 = FindFirstPulsePeak(Nsamples, Channel4Voltages_); // this is useful ONLY for Photonis MCP
       
      //////////////////////////////////////////////
      ///////// Done with setup, start the fits
      //////////////////////////////////////////////

      // invert the pulses for easy fitting
      CH1pulse->Scale(-1.);
      CH2pulse->Scale(-1.);
      CH3pulse->Scale(-1.);
      CH4pulse->Scale(-1.);
      
      // fit the baseline
      float base1 = LinearFit_Baseline( CH1pulse, index_min1, 10 );
      float base2 = LinearFit_Baseline( CH2pulse, index_min2, 10 );
      float base3 = LinearFit_Baseline( CH3pulse, index_min3, 10 );
      float base4 = LinearFit_Baseline( CH4pulse, index_min4, 10 );
            
      ///////////////////
      // Gaussian fit
      // fit the gaussian peak
      ///////////////////
      float timepeak1 =  GausFit_MeanTime(CH1pulse, index_min1 - 3, index_min1+4);
      float timepeak2 =  GausFit_MeanTime(CH2pulse, index_min2 - 3, index_min2+4);
      float timepeak3 =  GausFit_MeanTime(CH3pulse, index_min3 - 3, index_min3+4);
      float timepeak4 =  GausFit_MeanTime(CH4pulse, index_min4 - 3, index_min4+4);
      
      
      ch1Time_gausfitroot = timepeak1*0.2;
      ch2Time_gausfitroot = timepeak2*0.2;
      ch3Time_gausfitroot = timepeak3*0.2;
      ch4Time_gausfitroot = timepeak4*0.2;
            
      /////////////////////////
      // Find the amplitudes
      /////////////////////////
      ch1Amp = -1 * Channel1Voltages_[index_min1] - base1;
      ch2Amp = -1 * Channel2Voltages_[index_min2] - base2;
      ch3Amp = -1 * Channel3Voltages_[index_min3] - base3;
      ch4Amp = -1 * Channel4Voltages_[index_min4] - base4;
      
      ch1Int = -1 * ChannelIntegral(Channel1Voltages_, index_min1) - 7 * base1;
      ch2Int = -1 * ChannelIntegral(Channel2Voltages_, index_min2) - 7 * base2;
      ch3Int = -1 * ChannelIntegral(Channel3Voltages_, index_min3) - 7 * base3;
      ch4Int = -1 * ChannelIntegral(Channel4Voltages_, index_min4) - 7 * base4;

      //FitFullPulse(CH2pulse, ch1_TFF, ch1_AFF);//ch1 Time Full Fit
      //FitFullPulse(CH3pulse, ch2_TFF, ch2_AFF);//ch2 Time Full Fit
      //FitFullPulse(CH4pulse, ch4_TFF, ch4_AFF, ch4_TFF_v2);//ch4 Time Full Fit
      
      //std::cout << "ch4_AFF_v2: " << ch4_TFF_v2 << std::endl;
      //std::cout << "ch2_AFF: " << ch2_AFF << std::endl;
      
      //Fit Rising Edge
      ch1THM = FitRisingEdge(CH1pulse, -1, 0);
      ch2THM = FitRisingEdge(CH2pulse, 0, 0);
      //FitRisingEdge(CH1pulse, ch1_AFF, ch1THM, ch1BL, base1);
      //FitRisingEdge(CH2pulse, ch1_AFF, ch1THM, ch2BL, base2);
      //FitRisingEdge(CH3pulse, ch2_AFF, ch2THM, ch2BL, base2);
      ch3THM = FitRisingEdge(CH3pulse, 0, 0);
      ch4THM = FitRisingEdge(CH4pulse, 0, 0);
      //FitRisingEdge(CH4pulse, ch4_AFF, ch4THM, ch4BL, base4);
      //std::cout << "ch4THM: " << ch4THM << "  ch4BL: " << ch4BL << std::endl;
      
      if(ch1Amp < 0.05 || ch2Amp < 0.05){
        //continue;
      }
      
      //Fill the tree
      treeOut->Fill();
    }
  
  CH1pulse->Write();
  CH2pulse->Write();
  CH3pulse->Write();
  CH4pulse->Write();
  
  GausPeak_CH12_dt->Write();
  GausPeak_CH34_dt->Write();
  GausPeak_TOF_CH13->Write();
  GausPeak_TOF_CH14->Write();
  
  CH1Amp->Write();
  CH2Amp->Write();
  
  treeOut->Write();

  fout->Close();
}

////////////////////////////////////////////
// find minimum of the pulse
// aa added protection against pulses with single high bin
////////////////////////////////////////////
int FindMin( int n, float *a) {
  
  if (n <= 0 || !a) return -1;
  float xmin = a[5];
  int loc = 0;
  for  (int i = 5; i < n-5; i++) {
    if (xmin > a[i] && a[i+1] < 0.5*a[i])  {
      xmin = a[i];
      loc = i;
    }
  }
  
  return loc;
}

////////////////////////////////////////////
// find maximum of the pulse 
////////////////////////////////////////////
int FindMax( int n, float *a) {

  if (n <= 0 || !a) return -1;
  float xmax = a[0];
  int loc = 0;
  for  (int i = 0; i < n; i++) {
    if (xmax < a[i]) {
      xmax = a[i];
      loc = i;
    }
  }
  return loc;

}

////////////////////////////////////////////
// find rising edge of the pulse
////////////////////////////////////////////
int FindRisingEdge( int n, int binMax, float *a) {

  if (n <= 0 || !a) return -1;
  float xmin = a[0];
  int loc = -99;
  for (int i = binMax-10; i <binMax; i++)
  { // sometimes there is noise, check that it is rising in three bins
    if ( a[i] < -0.01 && a[i+1] < a[i] && a[i+2] < a[i+1] )
      // if ( Channel1Voltages_[i+2] < 0.3*Channel1Voltages_[index_min1])
    {
      loc = i; 
      break;
    }
  }  
  return loc;  
}

////////////////////////////////////////////
// for photonis: find the first pulse 
////////////////////////////////////////////
int FindFirstPulsePeak( int n, float *a) {

  if (n <= 0 || !a) return -1;

  int loc = 0;
  for  (int i = 20; i < 1000; i++) {
    if ( a[i] < -0.015
         && a[i] <= a[i-1] && a[i] <= a[i+1]
         && a[i-1] <= a[i-2]
         && a[i-2] <= a[i-3]
      ) {

      loc = i;
      break;
    }
  }
  return loc;

}

////////////////////////////////////////////
// assign pulse quality bit
////////////////////////////////////////////
unsigned int CheckPulseQuality( int binMin, int binMax, float *a, float minPulse) {
  unsigned int answer = 0;

  //*******************************************
  //Check if there is real pulse in the event
  //*******************************************
  if (!(a[binMin] < -minPulse)) {
    answer |= kNoPulse;
    // std::cout << "No Pulse: " << answer << "\n";
  }

  //*******************************************
  //Check pulses with large opposite amplitude 
  //*******************************************
  if (a[binMax] > 0.1) {
    answer |= kLargeNegativeAmplitude;
    //cout << "kLargeNegativeAmplitude: " << answer << "\n";
  }

  //*******************************************
  //Check pulses with saturation
  //*******************************************
  if ( (a[binMin] == a[binMin+1] || a[binMin] == a[binMin-1])
    ) {
    answer |= kSaturated;
    //cout << "kSaturated: " << answer << "\n";
  }


  //************************************************************************************
  //check that no points near the peak has negative polarity
  //************************************************************************************
  for (int ii=binMin-3;ii<binMin+4;ii++)
  {
    //cout << "bin " << ii << " : " << a[ii] << "\n";
    if (a[ii] > 0) {
      answer |= kNegativePolarity;
      //cout << "Negative Polarity\n";
    }
  }


  //************************************************************************************
  //check that no points near the peak has sudden jump
  //************************************************************************************
  bool hasSuddenJump = false;
  if (!(a[binMin] < a[binMin+1] && a[binMin+1] < a[binMin+2] && a[binMin+2] < a[binMin+3]  && a[binMin+3] < a[binMin+4])) hasSuddenJump = true;
  if (!(a[binMin] < a[binMin-1] && a[binMin-1] < a[binMin-2] && a[binMin-2] < a[binMin-3]  && a[binMin-3] < a[binMin-4])) hasSuddenJump = true;

  if (hasSuddenJump == true) {
//     cout << "Sudden Jump\n";
//     cout << a[binMin] << " " 
//          <<  a[binMin+1] << " "
//          <<  a[binMin+2] << " "
//          <<  a[binMin+3] << " "
//            <<  a[binMin+4] << " "
//          << "\n";
//     cout << a[binMin] << " " 
//          <<  a[binMin-1] << " "
//          <<  a[binMin-2] << " "
//          <<  a[binMin-3] << " "
//          <<  a[binMin-4] << " "
//          << "\n";
    answer |= kSuddenJump;
    //cout << "kSuddenJump: " << answer << "\n";
  }
  
  //************************************************************************************
  //check that second point away from peak is at least 10% lower - prevents strange pulses where there is flattish (but not completely flat) top
  //************************************************************************************
  bool hasFlatTop = false;
  if (!(fabs(a[binMin] - a[binMin+2])/fabs(a[binMin]) > 0.05)) hasFlatTop = true;
  if (!(fabs(a[binMin] - a[binMin-2])/fabs(a[binMin]) > 0.05)) hasFlatTop = true;

  if (hasFlatTop == true) {
    // std::cout << "flat pulse\n";
    // std::cout << a[binMin] << " " <<binMin<< " " 
    //      <<  a[binMin+2] << " "
    //      <<  fabs(a[binMin] - a[binMin+2])/fabs(a[binMin]) << " "
    //      << "\n";
    // std::cout << a[binMin] << " " 
    //      <<  a[binMin-2] << " "
    //      <<  fabs(a[binMin] - a[binMin-2])/fabs(a[binMin]) << " "

    //      << "\n";
    answer |= kFlatTop;
    //cout << "kFlatTop: " << answer << "\n";
  }

  //************************************************************************************
  //check for presence of 2nd pulse
  //************************************************************************************
  bool hasSecondPulse = false;
  int secondPulseIndex = -99;
  float secondPulseHeight = 0;
  for  (int i = 20; i < 1000; i++) {
    if (secondPulseHeight > a[i] 
        && fabs(a[i] - a[i-1])/fabs(a[i]) < 0.5
        && fabs(a[i] - a[i+1])/fabs(a[i]) < 0.5
        && abs(binMin - i) > 20
        && fabs(fabs(a[i]) - fabs(a[binMin]))/fabs(a[binMin]) < 0.15
      ) {
      secondPulseHeight = a[i];
      secondPulseIndex = i;
      hasSecondPulse = true;
    }
  }

//     cout << "Second Pulse\n";
//     cout << secondPulseIndex << " : " << secondPulseHeight << "\n";
  if (hasSecondPulse) {
    // std::cout << "First Pulse\n";
    // std::cout << binMin << " : " << a[binMin] << "\n";
    // std::cout << "Second Pulse\n";
    // std::cout << secondPulseIndex << " : " << secondPulseHeight << "\n";
    answer |= kSecondPulse;
  }
  
  //cout << "final answer : " << answer << "\n";
  return answer; 
}

// find the baseline
float LinearFit_Baseline(TH1F * pulse, const int index_min, const int range)
{
  TF1 *fBaseline = new TF1("fBaseline","pol0",0, index_min-range);
  pulse->Fit("fBaseline","Q","", 0, index_min-range);
  float base = fBaseline->GetParameter(0);
  delete fBaseline;
  
  return base;
}

// find the intercept of the linear fit on the rising edge
float LinearFit_Intercept(TH1F* pulse, const float base, const int index_first, const int index_last)
{
  TF1* fRise = new TF1("fRise","pol1", index_first, index_last);
  pulse->Fit("fRise","Q","", index_first, index_last);
  float timeIntercept = (base - fRise->GetParameter(0))/fRise->GetParameter(1);
  delete fRise;

  return timeIntercept;
}

// find the mean time from gaus fit
float GausFit_MeanTime(TH1F* pulse, const int index_first, const int index_last)
{
  TF1* fpeak = new TF1("fpeak","gaus", index_first, index_last);
  pulse->Fit("fpeak","Q","", index_first, index_last);
  float timepeak = fpeak->GetParameter(1);
  delete fpeak;
  
  return timepeak;
}

float ChannelIntegral(float *a, int peak) 
{
  float integral = 0.;

  //for sharp gaussian peak type pulses
  //integral  = a[peak - 3] + a[peak - 2] + a[peak - 1] + a[peak] + a[peak + 1] + a[peak + 2] + a[peak + 3];

  //for scintillation type pulses
  for (int i= std::max(peak - 100,2); i < std::min (peak + 800, 1023); ++i) integral += a[i];

  return integral;
}

float FitRisingEdge(TH1F* pulse, int nbinsL, int nbinsH){
  int bM = pulse->FindFirstBinAbove(0.6*pulse->GetMaximum());
  int bL = pulse->FindFirstBinAbove(0.1*pulse->GetMaximum());
  TF1* f = new TF1("f", "[0]+x*[1]", pulse->GetBinCenter(bL+nbinsL), pulse->GetBinCenter(bM+nbinsH));
  //pulse->Fit(f,"MWLR");
  pulse->Fit(f,"RQ");
  float m = f->GetParameter(1);
  float b = f->GetParameter(0);
  delete f;
  //std::cout << "HTM: " << 0.2*(0.2*pulse->GetMaximum()-b)/m << std::endl;
  return  0.2*(0.2*pulse->GetMaximum()-b)/m;//converted to picoseconds
  
}

float FitRisingEdge(TH1F* pulse, float Max){
  int bM = pulse->FindFirstBinAbove(0.8*Max);
  int bL = pulse->FindFirstBinAbove(0.125*Max);
  TF1* f = new TF1("f", "[0]+x*[1]", pulse->GetBinCenter(bL), pulse->GetBinCenter(bM));
  //pulse->Fit(f,"MWLR");
  pulse->Fit(f,"RQ");
  float m = f->GetParameter(1);
  float b = f->GetParameter(0);
  delete f;
  //std::cout << "HTM: " << 0.2*(0.3*pulse->GetMaximum()-b)/m << std::endl;
  return  0.2*(0.5*pulse->GetMaximum()-b)/m;//converted to picoseconds
}

void FitRisingEdge(TH1F* pulse, float Max, float &THM, float &t0, float baseline){
  int bM = pulse->FindFirstBinAbove(0.8*Max);
  int bL = pulse->FindFirstBinAbove(0.125*Max);
  TF1* f = new TF1("f", "[0]+x*[1]", pulse->GetBinCenter(bL), pulse->GetBinCenter(bM));
  //pulse->Fit(f,"MWLR");
  pulse->Fit(f,"RQ");
  float m = f->GetParameter(1);
  float b = f->GetParameter(0);
  THM = 0.20*(0.15*pulse->GetMaximum()-b)/m;//converted to picoseconds
  t0 = 0.20*(baseline-b)/m;//converted to picoseconds
  delete f;
}

void FitFullPulse(TH1F* pulse, float &par0, float &par1, float &par2){
  TF1* f = new TF1("f","[2]+[1]*0.004217/2*exp(0.004217/2*(2.0*357.0+1.0*2.06**2.0-2.0*(x-[0])))*ROOT::Math::erfc((357.0+0.004217*2.06**2.0-(x-[0]))/(1.41*2.06))",10,1000);
  f->SetParameter(0, 40.0);
  f->SetParameter(1, 60.0);
  f->SetParameter(0, 0.0);
  //pulse->Fit(f,"MWLR");
  pulse->Fit(f,"RQ");
  par0 = f->GetParameter(0);
  par1 = f->GetMaximum(10,1000);
  par2 = f->GetX(10,0.3*par1);
  delete f;
}
