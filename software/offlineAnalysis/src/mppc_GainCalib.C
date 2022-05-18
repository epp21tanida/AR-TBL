#include <iostream>
#include <fstream>

#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"

TString root_file = "../../../../artbl_data/20220512_2_Sr_1000000_af.root";
TString gain_file = "./output/peakFind_out_X.dat";
TString ofilename = "./output/mppc_gain/20220512_2_Sr_1000000_af.pdf";
TString mip_file  = "./output/mppc_gain/20220512_2_Sr_1000000_af.dat";

Int_t mppc_GainCalib() {

  gStyle->SetOptFit(1);
  
  std::cout << "mppc_GainCalib ..." << std::endl;
  //  mkdir ./output/GainCalib
  
  const Int_t N = 64;
  

  // INPUT FILE
  TFile* ifile = TFile::Open(root_file.Data(), "READ");
  if (!ifile->IsOpen()) {
    std::cerr << "[error] file is not opened" << std::endl;
    return 1;
  }
  
  TTree* itree = dynamic_cast<TTree*>(ifile->Get("eventtree"));
  if (!itree) {
    std::cerr << "[error] tree is not found" << std::endl;
    return 1;
  }
  
  // itree->Print();
  
  // LOAD PEDESTAL AND GAIN FILE
  Double_t pedestal[N] = {}, gain[N] = {};
  
  std::ifstream gainfile(gain_file);
  if (!gainfile) {
    std::cerr << "[error] x gain file is not found" << std::endl;
    return 1;
  }
  
  Int_t ch, ipeak;
  Double_t adc;
  while (gainfile >> ch >> ipeak >> adc) {
    
    if (ipeak == 0) {pedestal[ch] = adc;} 
    else            {gain[ch]     = adc - pedestal[ch];}
  }
  
  // SET TEMPORAL VALUES FOR PEDESTAL AND GAIN
  // TH1F* temp_array[128]={0x0};  
  for (Int_t ch = 0; ch < N; ++ch) {
    itree->Draw(Form("adc_ch[%d] >> htemp(4000,0,4000)", ch));
    TH1* htemp = ((TH1F*)gROOT->FindObject("htemp"));
    pedestal[ch] = htemp->GetBinCenter(htemp->GetMaximumBin());
    //    temp_array[ch] = ((TH1F*)gROOT->FindObject("htemp"))->Clone(Form("htemp_%d", ch));
    //    gain[ch] = 50;
  }
  
  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1000);
  
  // TF1 GAUS FOR FITTING
  TF1* fGauss = new TF1("fGauss", "gaus(0)");
  
  // open output mip file
  std::ofstream mipfile(mip_file);
  if (!mipfile) {
    std::cerr << "[error] mip file is not found" << std::endl;
    return 1;
  }
  
  // draw and fit for each channel
  TCanvas::MakeDefCanvas();
  gPad->SetLogy();
  gPad->Print(ofilename + "[");
  
  for (Int_t ch = 0; ch < N; ++ch) {
    
    // draw raw hist
    itree->Draw(Form("adc_ch[%d] >> htemp_%d(4000,0,4000)", ch, ch));
    gPad->Print(ofilename);
    
    // draw normalized hist
    itree->Draw(Form("(adc_ch[%d]-%lf)/%lf >> htemp_%d_nm(52,-1.5,50.5)", ch, pedestal[ch], gain[ch], ch));
    gPad->Print(ofilename);
    
    // get normalized hist
    TH1* htemp_nm = dynamic_cast<TH1*>(gROOT->FindObject(Form("htemp_%d_nm", ch)));
    htemp_nm->SetTitle(Form("%s;[p.e.];Entries", htemp_nm->GetTitle()));
    
    // (tmp) fit mip using gaussian
    htemp_nm->Fit(fGauss, "", "", 10, 35);
    gPad->Print(ofilename);
    
    // write mip value
    mipfile << ch << "\t" << fGauss->GetParameter(1) << std::endl;
  }
  
  gPad->Print(ofilename + "]");
  
  return 0;
}
