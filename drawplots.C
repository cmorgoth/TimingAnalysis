// run 012
// CH12
TCanvas * c = new TCanvas("c","c",600,600)
h1 = new TH1F("h1","h1",1000,-1.5,0.5)
tree->Draw("t1gausroot-t2gausroot>>h1","ch1QualityBit==0&&ch2QualityBit==0")
h1->SetAxisRange(h1->GetMean()-0.5*fabs(h1->GetMean()),h1->GetMean()+0.5*fabs(h1->GetMean()),"X")
h1->GetXaxis()->SetTitle("CH1-CH2 [ns]")
h1->GetYaxis()->SetTitle("Number of Events")
h1->Draw()
h1->Fit("gaus")
TVirtualFitter * fitter = TVirtualFitter::GetFitter();
TPaveText *pt = new TPaveText(.15,.75,.35,.85,"NDC");
pt->AddText(Form("#sigma=%.2f psec",1000*fitter->GetParameter(2)));
pt->SetFillColor(0);
pt->Draw()

// CH56
h1 = new TH1F("h1","h1",1000,-1.5,0.5)
tree->Draw("t5gausroot-t6gausroot>>h1","ch5QualityBit==0&&ch6QualityBit==0")
h1->SetAxisRange(h1->GetMean()-0.5*fabs(h1->GetMean()),h1->GetMean()+0.5*fabs(h1->GetMean()),"X")
h1->GetXaxis()->SetTitle("CH5-CH6 [ns]")
h1->GetYaxis()->SetTitle("Number of Events")
h1->Draw()
h1->Fit("gaus")
TVirtualFitter * fitter = TVirtualFitter::GetFitter();
TPaveText *pt = new TPaveText(.15,.75,.35,.85,"NDC");
pt->AddText(Form("#sigma=%.2f psec",1000*fitter->GetParameter(2)));
pt->SetFillColor(0);
pt->Draw()

// CH13
h1 = new TH1F("h1","h1",1000,0.5,2.5)
tree->Draw("t1gausroot-t3gausroot>>h1","ch1QualityBit==0&&ch3QualityBit==0")
h1->SetAxisRange(h1->GetMean()-0.5*fabs(h1->GetMean()),h1->GetMean()+0.5*fabs(h1->GetMean()),"X")
h1->GetXaxis()->SetTitle("CH1-CH3 [ns]")
h1->GetYaxis()->SetTitle("Number of Events")
h1->Draw()
h1->Fit("gaus")
TVirtualFitter * fitter = TVirtualFitter::GetFitter();
TPaveText *pt = new TPaveText(.15,.75,.35,.85,"NDC");
pt->AddText(Form("#sigma=%.2f psec",1000*fitter->GetParameter(2)));
pt->SetFillColor(0);
pt->Draw()

// CH13 (ch2Amp>0.49)
h1 = new TH1F("h1","h1",1000,0.5,2.5)
tree->Draw("t1gausroot-t3gausroot>>h1","ch1QualityBit==0&&ch3QualityBit==0&&ch6Amp>0.49&&ch3Amp>0.15")
h1->SetAxisRange(h1->GetMean()-0.5*fabs(h1->GetMean()),h1->GetMean()+0.5*fabs(h1->GetMean()),"X")
h1->GetXaxis()->SetTitle("CH1-CH3 [ns]")
h1->GetYaxis()->SetTitle("Number of Events")
h1->Draw()
h1->Fit("gaus")
TVirtualFitter * fitter = TVirtualFitter::GetFitter();
TPaveText *pt = new TPaveText(.15,.75,.35,.85,"NDC");
pt->AddText(Form("#sigma=%.2f psec",1000*fitter->GetParameter(2)));
pt->SetFillColor(0);
pt->Draw()

// CH23
h1 = new TH1F("h1","h1",1000,0.5,2.5)
tree->Draw("t2gausroot-t3gausroot>>h1","ch2QualityBit==0&&ch3QualityBit==0")
h1->SetAxisRange(h1->GetMean()-0.5*fabs(h1->GetMean()),h1->GetMean()+0.5*fabs(h1->GetMean()),"X")
h1->GetXaxis()->SetTitle("CH2-CH3 [ns]")
h1->GetYaxis()->SetTitle("Number of Events")
h1->Draw()
h1->Fit("gaus")
TVirtualFitter * fitter = TVirtualFitter::GetFitter();
TPaveText *pt = new TPaveText(.15,.75,.35,.85,"NDC");
pt->AddText(Form("#sigma=%.2f psec",1000*fitter->GetParameter(2)));
pt->SetFillColor(0);
pt->Draw()

// CH57
h1 = new TH1F("h1","h1",1000,0.5,2.5)
tree->Draw("t5gausroot-t7gausroot>>h1","ch5QualityBit==0&&ch7QualityBit==0")
h1->SetAxisRange(h1->GetMean()-0.5*fabs(h1->GetMean()),h1->GetMean()+0.5*fabs(h1->GetMean()),"X")
h1->GetXaxis()->SetTitle("CH5-CH7 [ns]")
h1->GetYaxis()->SetTitle("Number of Events")
h1->Draw()
h1->Fit("gaus")
TVirtualFitter * fitter = TVirtualFitter::GetFitter();
TPaveText *pt = new TPaveText(.15,.75,.35,.85,"NDC");
pt->AddText(Form("#sigma=%.2f psec",1000*fitter->GetParameter(2)));
pt->SetFillColor(0);
pt->Draw()

// CH57 (ch2Amp>0.49)
h1 = new TH1F("h1","h1",1000,0.5,2.5)
tree->Draw("t5gausroot-t7gausroot>>h1","ch5QualityBit==0&&ch7QualityBit==0&&ch6Amp>0.49")
h1->SetAxisRange(h1->GetMean()-0.5*fabs(h1->GetMean()),h1->GetMean()+0.5*fabs(h1->GetMean()),"X")
h1->GetXaxis()->SetTitle("CH5-CH7 [ns]")
h1->GetYaxis()->SetTitle("Number of Events")
h1->Draw()
h1->Fit("gaus")
TVirtualFitter * fitter = TVirtualFitter::GetFitter();
TPaveText *pt = new TPaveText(.15,.75,.35,.85,"NDC");
pt->AddText(Form("#sigma=%.2f psec",1000*fitter->GetParameter(2)));
pt->SetFillColor(0);
pt->Draw()

// CH67
h1 = new TH1F("h1","h1",1000,1.2,3.2)
tree->Draw("t6gausroot-t7gausroot>>h1","ch6QualityBit==0&&ch7QualityBit==0")
h1->SetAxisRange(h1->GetMean()-0.5*fabs(h1->GetMean()),h1->GetMean()+0.5*fabs(h1->GetMean()),"X")
h1->GetXaxis()->SetTitle("CH6-CH7 [ns]")
h1->GetYaxis()->SetTitle("Number of Events")
h1->Draw()
h1->Fit("gaus")
TVirtualFitter * fitter = TVirtualFitter::GetFitter();
TPaveText *pt = new TPaveText(.15,.75,.35,.85,"NDC");
pt->AddText(Form("#sigma=%.2f psec",1000*fitter->GetParameter(2)));
pt->SetFillColor(0);
pt->Draw()

/////////////////////////////////
/////// Amplitudes
/////////////////////////////////
gStyle->SetOptStat(0)
h1 = new TH1F("h1","h1",200,0,0.55)
h2 = new TH1F("h2","h2",200,0,0.55)
tree->Draw("ch1Amp>>h1","");
tree->Draw("ch1Amp>>h2","ch1QualityBit==0");
h1->GetXaxis()->SetTitle("CH1 Amp [mV]")
h1->GetYaxis()->SetTitle("Number of Events")
h1->Draw()
h2->SetLineColor(2)
h2->Draw("same")
leg = new TLegend(0.5,0.8,0.8,0.9);
leg->AddEntry(h1,"All pulses")
leg->AddEntry(h2,"Pass all quality bits")
leg->SetFillColor(0)
leg->Draw();
gPad->SetLogy()

h1 = new TH1F("h1","h1",200,0,0.55)
h2 = new TH1F("h2","h2",200,0,0.55)
tree->Draw("ch2Amp>>h1","");
tree->Draw("ch2Amp>>h2","ch2QualityBit==0");
h1->GetXaxis()->SetTitle("CH2 Amp [mV]")
h1->GetYaxis()->SetTitle("Number of Events")
h1->Draw()
h2->SetLineColor(2)
h2->Draw("same")
leg->Draw();
gPad->SetLogy()

h1 = new TH1F("h1","h1",200,0,0.55)
h2 = new TH1F("h2","h2",200,0,0.55)
tree->Draw("ch3Amp>>h1","");
tree->Draw("ch3Amp>>h2","ch3QualityBit==0");
h1->GetXaxis()->SetTitle("CH3 Amp [mV]")
h1->GetYaxis()->SetTitle("Number of Events")
h1->Draw()
h2->SetLineColor(2)
h2->Draw("same")
leg->Draw();
gPad->SetLogy()

h1 = new TH1F("h1","h1",200,0,0.55)
h2 = new TH1F("h2","h2",200,0,0.55)
tree->Draw("ch5Amp>>h1","");
tree->Draw("ch5Amp>>h2","ch5QualityBit==0");
h1->GetXaxis()->SetTitle("CH5 Amp [mV]")
h1->GetYaxis()->SetTitle("Number of Events")
h1->Draw()
h2->SetLineColor(2)
h2->Draw("same")
leg->Draw();
gPad->SetLogy()

h1 = new TH1F("h1","h1",200,0,0.55)
h2 = new TH1F("h2","h2",200,0,0.55)
tree->Draw("ch6Amp>>h1","");
tree->Draw("ch6Amp>>h2","ch6QualityBit==0");
h1->GetXaxis()->SetTitle("CH6 Amp [mV]")
h1->GetYaxis()->SetTitle("Number of Events")
h1->Draw()
h2->SetLineColor(2)
h2->Draw("same")
leg->Draw();
gPad->SetLogy()

h1 = new TH1F("h1","h1",200,0,0.55)
h2 = new TH1F("h2","h2",200,0,0.55)
tree->Draw("ch7Amp>>h1","");
tree->Draw("ch7Amp>>h2","ch7QualityBit==0");
h1->GetXaxis()->SetTitle("CH7 Amp [mV]")
h1->GetYaxis()->SetTitle("Number of Events")
h1->Draw()
h2->SetLineColor(2)
h2->Draw("same")
leg->Draw();
gPad->SetLogy()

h1 = new TH1F("h1","h1",200,0,0.55)
h2 = new TH1F("h2","h2",200,0,0.55)
tree->Draw("ch7Amp>>h1","");
tree->Draw("ch7Amp>>h2","ch7QualityBit==0");
h1->GetXaxis()->SetTitle("CH7 Amp [mV]")
h1->GetYaxis()->SetTitle("Number of Events")
h1->Draw()
h2->SetLineColor(2)
h2->Draw("same")
leg->Draw();
gPad->SetLogy()
