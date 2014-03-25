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

  float Channel1Voltages_[1024];
  float Channel2Voltages_[1024];
  float Channel3Voltages_[1024];
  float Channel4Voltages_[1024];
  float Channel5Voltages_[1024];
  float Channel6Voltages_[1024];
  float Channel7Voltages_[1024];
  float Channel8Voltages_[1024];

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

      t1->SetBranchAddress("Channel1Voltages",Channel1Voltages_);
      t1->SetBranchAddress("Channel2Voltages",Channel2Voltages_);
      t1->SetBranchAddress("Channel3Voltages",Channel3Voltages_);
      t1->SetBranchAddress("Channel4Voltages",Channel4Voltages_);
      t1->SetBranchAddress("Channel5Voltages",Channel5Voltages_);
      t1->SetBranchAddress("Channel6Voltages",Channel6Voltages_);
      t1->SetBranchAddress("Channel7Voltages",Channel7Voltages_);
      t1->SetBranchAddress("Channel8Voltages",Channel8Voltages_);
    }

  //create two histograms
  TH1F *CH1pulse   = new TH1F("CH1pulse","CH1pulse",1024,0,1024);
  TH1F *CH2pulse   = new TH1F("CH2pulse","CH2pulse",1024,0,1024);
  TH1F *CH3pulse   = new TH1F("CH3pulse","CH3pulse",1024,0,1024);
  TH1F *CH4pulse   = new TH1F("CH4pulse","CH4pulse",1024,0,1024);
  TH1F *CH5pulse   = new TH1F("CH5pulse","CH5pulse",1024,0,1024);
  TH1F *CH6pulse   = new TH1F("CH6pulse","CH6pulse",1024,0,1024);
  TH1F *CH7pulse   = new TH1F("CH7pulse","CH7pulse",1024,0,1024);
  TH1F *CH8pulse   = new TH1F("CH8pulse","CH8pulse",1024,0,1024);

  TH1F *CH1Amp     = new TH1F("CH1Amp","CH1Amp",40,-0.3,-0.6);
  TH1F *CH2Amp     = new TH1F("CH2Amp","CH2Amp",40,-0.3,-0.6);

  TH1F *GausPeak_CH12_dt   = new TH1F("GausPeak_CH12_dt","GausPeak_CH12_dt; t1-t2 [ns]; Events",4000,-4,4);
  TH1F *GausPeak_CH34_dt   = new TH1F("GausPeak_CH34_dt","GausPeak_CH34_dt; t3-t4 [ns]; Events",4000,-4,4);
  TH1F *GausPeak_TOF_CH13   = new TH1F("GausPeak_TOF_CH13","GausPeak_TOF_CH13; t1-t3 [ns]; Events",4000,-4,4);
  TH1F *GausPeak_TOF_CH14   = new TH1F("GausPeak_TOF_CH14","GausPeak_TOF_CH14; t1-t4 [ns]; Events",4000,-40,4);

  // Create the output file with a TTree
  TFile fout(argv[2],"recreate");
  TTree *treeOut = new TTree("tree","tree");
  
  unsigned int eventNumber = 0;
  float ch1Time_gausfitroot = 0;
  float ch2Time_gausfitroot = 0;
  float ch3Time_gausfitroot = 0;
  float ch4Time_gausfitroot = 0;
  float ch5Time_gausfitroot = 0;
  float ch6Time_gausfitroot = 0;
  float ch7Time_gausfitroot = 0;
  float ch8Time_gausfitroot = 0;
  float ch1Amp = 0;
  float ch2Amp = 0;
  float ch3Amp = 0;
  float ch4Amp = 0;
  float ch5Amp = 0;
  float ch6Amp = 0;
  float ch7Amp = 0;
  float ch8Amp = 0;
  float ch1Int = 0;
  float ch2Int = 0;
  float ch3Int = 0;
  float ch4Int = 0;
  float ch5Int = 0;
  float ch6Int = 0;
  float ch7Int = 0;
  float ch8Int = 0;
  unsigned int ch1QualityBit = 0;
  unsigned int ch2QualityBit = 0;
  unsigned int ch3QualityBit = 0;
  unsigned int ch4QualityBit = 0;
  unsigned int ch5QualityBit = 0;
  unsigned int ch6QualityBit = 0;
  unsigned int ch7QualityBit = 0;
  unsigned int ch8QualityBit = 0;
  float ch1chisq = -1;
  float ch2chisq = -1;
  float ch3chisq = -1;
  float ch4chisq = -1;
  float ch5chisq = -1;
  float ch6chisq = -1;
  float ch7chisq = -1;
  float ch8chisq = -1;

  treeOut->Branch("event",&eventNumber,"event/i");
  treeOut->Branch("t1gausroot",&ch1Time_gausfitroot,"t1gausroot/F");
  treeOut->Branch("t2gausroot",&ch2Time_gausfitroot,"t2gausroot/F");
  treeOut->Branch("t3gausroot",&ch3Time_gausfitroot,"t3gausroot/F");
  treeOut->Branch("t4gausroot",&ch4Time_gausfitroot,"t4gausroot/F");
  treeOut->Branch("t5gausroot",&ch5Time_gausfitroot,"t5gausroot/F");
  treeOut->Branch("t6gausroot",&ch6Time_gausfitroot,"t6gausroot/F");
  treeOut->Branch("t7gausroot",&ch7Time_gausfitroot,"t7gausroot/F");
  treeOut->Branch("t8gausroot",&ch8Time_gausfitroot,"t8gausroot/F");
  treeOut->Branch("ch1Amp",&ch1Amp,"ch1Amp/F");
  treeOut->Branch("ch2Amp",&ch2Amp,"ch2Amp/F");
  treeOut->Branch("ch3Amp",&ch3Amp,"ch3Amp/F");
  treeOut->Branch("ch4Amp",&ch4Amp,"ch4Amp/F");
  treeOut->Branch("ch5Amp",&ch5Amp,"ch5Amp/F");
  treeOut->Branch("ch6Amp",&ch6Amp,"ch6Amp/F");
  treeOut->Branch("ch7Amp",&ch7Amp,"ch7Amp/F");
  treeOut->Branch("ch8Amp",&ch8Amp,"ch8Amp/F");
  treeOut->Branch("ch1QualityBit",&ch1QualityBit,"ch1QualityBit/i");
  treeOut->Branch("ch2QualityBit",&ch2QualityBit,"ch2QualityBit/i");
  treeOut->Branch("ch3QualityBit",&ch3QualityBit,"ch3QualityBit/i");
  treeOut->Branch("ch4QualityBit",&ch4QualityBit,"ch4QualityBit/i");
  treeOut->Branch("ch5QualityBit",&ch5QualityBit,"ch5QualityBit/i");
  treeOut->Branch("ch6QualityBit",&ch6QualityBit,"ch6QualityBit/i");
  treeOut->Branch("ch7QualityBit",&ch7QualityBit,"ch7QualityBit/i");
  treeOut->Branch("ch8QualityBit",&ch8QualityBit,"ch8QualityBit/i");
  treeOut->Branch("ch1Int",&ch1Int,"ch1Int/F");
  treeOut->Branch("ch2Int",&ch2Int,"ch2Int/F");
  treeOut->Branch("ch3Int",&ch3Int,"ch3Int/F");
  treeOut->Branch("ch4Int",&ch4Int,"ch4Int/F");
  treeOut->Branch("ch5Int",&ch5Int,"ch5Int/F");
  treeOut->Branch("ch6Int",&ch6Int,"ch6Int/F");
  treeOut->Branch("ch7Int",&ch7Int,"ch7Int/F");
  treeOut->Branch("ch8Int",&ch8Int,"ch8Int/F");
  treeOut->Branch("ch1chisq",&ch1chisq,"ch1chisq/F");
  treeOut->Branch("ch2chisq",&ch2chisq,"ch2chisq/F");
  treeOut->Branch("ch3chisq",&ch3chisq,"ch3chisq/F");
  treeOut->Branch("ch4chisq",&ch4chisq,"ch4chisq/F");
  treeOut->Branch("ch5chisq",&ch5chisq,"ch5chisq/F");
  treeOut->Branch("ch6chisq",&ch6chisq,"ch6chisq/F");
  treeOut->Branch("ch7chisq",&ch7chisq,"ch7chisq/F");
  treeOut->Branch("ch8chisq",&ch8chisq,"ch8chisq/F");

  //read all entries and fill the histograms
  Long64_t nentries = t1->GetEntries();

  for (Long64_t iEntry=0;iEntry<nentries;iEntry++) 
    {
      if(iEntry%100==0) std::cout<<"Processing Event: "<<iEntry<<" out of: "<<nentries<<std::endl;
      
      t1->GetEntry(iEntry);
      eventNumber = iEntry;

      //////////////////////
      // convert to Volts
      //////////////////////
      if(convert2Volts)
	for (int ii=0;ii<1024;ii++)
	  {
	    Channel1Voltages_[ii] = 0.001*Channel1Voltages_[ii];
	    Channel2Voltages_[ii] = 0.001*Channel2Voltages_[ii];
	    Channel3Voltages_[ii] = 0.001*Channel3Voltages_[ii];
	    Channel4Voltages_[ii] = 0.001*Channel4Voltages_[ii];
	    Channel5Voltages_[ii] = 0.001*Channel5Voltages_[ii];
	    Channel6Voltages_[ii] = 0.001*Channel6Voltages_[ii];
	    Channel7Voltages_[ii] = 0.001*Channel7Voltages_[ii];
	    Channel8Voltages_[ii] = 0.001*Channel8Voltages_[ii];
	  }
      //////////////////////
      // end convert to Volts
      //////////////////////
      
      // Find Min of the Channel data (Voltage)
      int index_min1 = FindMin (1024, Channel1Voltages_);	// return index of the min
      int index_min2 = FindMin (1024, Channel2Voltages_);	// return index of the min
      int index_min3 = FindMin (1024, Channel3Voltages_);	// return index of the min
      int index_min4 = FindMin (1024, Channel4Voltages_);	// return index of the min
      int index_min5 = FindMin (1024, Channel5Voltages_);	// return index of the min
      int index_min6 = FindMin (1024, Channel6Voltages_);	// return index of the min
      int index_min7 = FindMin (1024, Channel7Voltages_);	// return index of the min
      int index_min8 = FindMin (1024, Channel8Voltages_);	// return index of the min

      // Set histograms bins and errors
      for (int ii=0;ii<1024;ii++)
	{
	  CH1pulse->SetBinContent(ii+1,Channel1Voltages_[ii]);
	  CH2pulse->SetBinContent(ii+1,Channel2Voltages_[ii]);
	  CH3pulse->SetBinContent(ii+1,Channel3Voltages_[ii]);
	  CH4pulse->SetBinContent(ii+1,Channel4Voltages_[ii]);
	  CH5pulse->SetBinContent(ii+1,Channel5Voltages_[ii]);
	  CH6pulse->SetBinContent(ii+1,Channel6Voltages_[ii]);
	  CH7pulse->SetBinContent(ii+1,Channel7Voltages_[ii]);
	  CH8pulse->SetBinContent(ii+1,Channel8Voltages_[ii]);

	  CH1pulse->SetBinError(ii+1, 0.05*Channel1Voltages_[index_min1]);
	  CH2pulse->SetBinError(ii+1, 0.05*Channel2Voltages_[index_min2]);
	  CH3pulse->SetBinError(ii+1, 0.05*Channel3Voltages_[index_min3]);
	  CH4pulse->SetBinError(ii+1, 0.05*Channel4Voltages_[index_min4]);
	  CH5pulse->SetBinError(ii+1, 0.05*Channel5Voltages_[index_min1]);
	  CH6pulse->SetBinError(ii+1, 0.05*Channel6Voltages_[index_min2]);
	  CH7pulse->SetBinError(ii+1, 0.05*Channel7Voltages_[index_min3]);
	  CH8pulse->SetBinError(ii+1, 0.05*Channel8Voltages_[index_min4]);
	}     

      // Find Max of the Channel data (Voltage)
      int index_max1 = FindMax (1024, Channel1Voltages_);	// return index of the max
      int index_max2 = FindMax (1024, Channel2Voltages_);	// return index of the max
      int index_max3 = FindMax (1024, Channel3Voltages_);	// return index of the max
      int index_max4 = FindMax (1024, Channel4Voltages_);	// return index of the max
      int index_max5 = FindMax (1024, Channel5Voltages_);	// return index of the max
      int index_max6 = FindMax (1024, Channel6Voltages_);	// return index of the max
      int index_max7 = FindMax (1024, Channel7Voltages_);	// return index of the max
      int index_max8 = FindMax (1024, Channel8Voltages_);	// return index of the max

      
      //////////////////
      // Find the rising edge on CH1
      int fbin1 = FindRisingEdge(1024, index_min1, Channel1Voltages_);
      int fbin2 = FindRisingEdge(1024, index_min2, Channel2Voltages_);
      int fbin3 = FindRisingEdge(1024, index_min3, Channel3Voltages_);
      int fbin4 = FindRisingEdge(1024, index_min4, Channel4Voltages_);
      int fbin5 = FindRisingEdge(1024, index_min5, Channel5Voltages_);
      int fbin6 = FindRisingEdge(1024, index_min6, Channel6Voltages_);
      int fbin7 = FindRisingEdge(1024, index_min7, Channel7Voltages_);
      int fbin8 = FindRisingEdge(1024, index_min8, Channel8Voltages_);

      // Find the quality of the pulse
      ch1QualityBit = CheckPulseQuality ( index_min1, index_max1, Channel1Voltages_, 0.02);
      ch2QualityBit = CheckPulseQuality ( index_min2, index_max2, Channel2Voltages_, 0.02);
      ch3QualityBit = CheckPulseQuality ( index_min3, index_max3, Channel3Voltages_, 0.02);
      ch4QualityBit = CheckPulseQuality ( index_min4, index_max4, Channel4Voltages_, 0.02);
      ch5QualityBit = CheckPulseQuality ( index_min5, index_max5, Channel5Voltages_, 0.02);
      ch6QualityBit = CheckPulseQuality ( index_min6, index_max6, Channel6Voltages_, 0.02);
      ch7QualityBit = CheckPulseQuality ( index_min7, index_max7, Channel7Voltages_, 0.02);
      ch8QualityBit = CheckPulseQuality ( index_min8, index_max8, Channel8Voltages_, 0.02);

      // For the first version of Photonis data the pulse has many peaks --> find the first one
      int index_firstPulse1 = FindFirstPulsePeak(1024, Channel1Voltages_); // this is useful ONLY for Photonis MCP
      int index_firstPulse2 = FindFirstPulsePeak(1024, Channel2Voltages_); // this is useful ONLY for Photonis MCP
      int index_firstPulse3 = FindFirstPulsePeak(1024, Channel3Voltages_); // this is useful ONLY for Photonis MCP
      int index_firstPulse4 = FindFirstPulsePeak(1024, Channel4Voltages_); // this is useful ONLY for Photonis MCP
      int index_firstPulse5 = FindFirstPulsePeak(1024, Channel5Voltages_); // this is useful ONLY for Photonis MCP
      int index_firstPulse6 = FindFirstPulsePeak(1024, Channel6Voltages_); // this is useful ONLY for Photonis MCP
      int index_firstPulse7 = FindFirstPulsePeak(1024, Channel7Voltages_); // this is useful ONLY for Photonis MCP
      int index_firstPulse8 = FindFirstPulsePeak(1024, Channel8Voltages_); // this is useful ONLY for Photonis MCP
 
      //////////////////////////////////////////////
      ///////// Done with setup, start the fits
      //////////////////////////////////////////////

      // invert the pulses for easy fitting
      CH1pulse->Scale(-1.);
      CH2pulse->Scale(-1.);
      CH3pulse->Scale(-1.);
      CH4pulse->Scale(-1.);
      CH5pulse->Scale(-1.);
      CH6pulse->Scale(-1.);
      CH7pulse->Scale(-1.);
      CH8pulse->Scale(-1.);

      // fit the baseline
      float base1 = LinearFit_Baseline( CH1pulse, index_min1, 10 );
      float base2 = LinearFit_Baseline( CH2pulse, index_min2, 10 );
      float base3 = LinearFit_Baseline( CH3pulse, index_min3, 10 );
      float base4 = LinearFit_Baseline( CH4pulse, index_min4, 10 );
      float base5 = LinearFit_Baseline( CH5pulse, index_min5, 10 );
      float base6 = LinearFit_Baseline( CH6pulse, index_min6, 10 );
      float base7 = LinearFit_Baseline( CH7pulse, index_min7, 10 );
      float base8 = LinearFit_Baseline( CH8pulse, index_min8, 10 );
      
      ///////////////////
      // Gaussian fit
      // fit the gaussian peak
      ///////////////////
      float timepeak1 =  GausFit_MeanTime(CH1pulse, index_min1 - 3, index_min1+4);
      float timepeak2 =  GausFit_MeanTime(CH2pulse, index_min2 - 3, index_min2+4);
      float timepeak3 =  GausFit_MeanTime(CH3pulse, index_min3 - 3, index_min3+4);
      float timepeak4 =  GausFit_MeanTime(CH4pulse, index_min4 - 3, index_min4+4);
      float timepeak5 =  GausFit_MeanTime(CH5pulse, index_min5 - 3, index_min5+4);
      float timepeak6 =  GausFit_MeanTime(CH6pulse, index_min6 - 3, index_min6+4);
      float timepeak7 =  GausFit_MeanTime(CH7pulse, index_min7 - 3, index_min7+4);
      float timepeak8 =  GausFit_MeanTime(CH8pulse, index_min8 - 3, index_min8+4);

      ch1Time_gausfitroot = timepeak1*0.2;
      ch2Time_gausfitroot = timepeak2*0.2;
      ch3Time_gausfitroot = timepeak3*0.2;
      ch4Time_gausfitroot = timepeak4*0.2;
      ch5Time_gausfitroot = timepeak5*0.2;
      ch6Time_gausfitroot = timepeak6*0.2;
      ch7Time_gausfitroot = timepeak7*0.2;
      ch8Time_gausfitroot = timepeak8*0.2;
      
      /////////////////////////
      // Find the amplitudes
      /////////////////////////
      ch1Amp = -1 * Channel1Voltages_[index_min1] - base1;
      ch2Amp = -1 * Channel2Voltages_[index_min2] - base2;
      ch3Amp = -1 * Channel3Voltages_[index_min3] - base3;
      ch4Amp = -1 * Channel4Voltages_[index_min4] - base4;
      ch5Amp = -1 * Channel5Voltages_[index_min5] - base5;
      ch6Amp = -1 * Channel6Voltages_[index_min6] - base6;
      ch7Amp = -1 * Channel7Voltages_[index_min7] - base7;
      ch8Amp = -1 * Channel8Voltages_[index_min8] - base8;

      ch1Int = -1 * ChannelIntegral(Channel1Voltages_, index_min1) - 7 * base1;
      ch2Int = -1 * ChannelIntegral(Channel2Voltages_, index_min2) - 7 * base2;
      ch3Int = -1 * ChannelIntegral(Channel3Voltages_, index_min3) - 7 * base3;
      ch4Int = -1 * ChannelIntegral(Channel4Voltages_, index_min4) - 7 * base4;
      ch5Int = -1 * ChannelIntegral(Channel5Voltages_, index_min5) - 7 * base5;
      ch6Int = -1 * ChannelIntegral(Channel6Voltages_, index_min6) - 7 * base6;
      ch7Int = -1 * ChannelIntegral(Channel7Voltages_, index_min7) - 7 * base7;
      ch8Int = -1 * ChannelIntegral(Channel8Voltages_, index_min8) - 7 * base8;
      
      //Fill the tree
      treeOut->Fill();
    }
  
  CH1pulse->Write();
  CH2pulse->Write();
  CH3pulse->Write();
  CH4pulse->Write();
  CH5pulse->Write();
  CH6pulse->Write();
  CH7pulse->Write();
  CH8pulse->Write();

  GausPeak_CH12_dt->Write();
  GausPeak_CH34_dt->Write();
  GausPeak_TOF_CH13->Write();
  GausPeak_TOF_CH14->Write();
  
  CH1Amp->Write();
  CH2Amp->Write();
  
  treeOut->Write();

  fout.Close();
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
  
  return base;
}

// find the intercept of the linear fit on the rising edge
float LinearFit_Intercept(TH1F * pulse, const float base, const int index_first, const int index_last)
{
  TF1 *fRise = new TF1("fRise","pol1", index_first, index_last);
  pulse->Fit("fRise","Q","", index_first, index_last);
  float timeIntercept = (base - fRise->GetParameter(0))/fRise->GetParameter(1);

  return timeIntercept;
}

// find the mean time from gaus fit
float GausFit_MeanTime(TH1F * pulse, const int index_first, const int index_last)
{
  TF1 *fpeak = new TF1("fpeak","gaus", index_first, index_last);
  pulse->Fit("fpeak","Q","", index_first, index_last);
  float timepeak = fpeak->GetParameter(1);
  
  return timepeak;
}

float ChannelIntegral(float *a, int peak) 
{
  float integral = 0.;
  integral  = a[peak - 3] + a[peak - 2] + a[peak - 1] + a[peak] + a[peak + 1] + a[peak + 2] + a[peak + 3];

  return integral;
}
