#include "myAnalysis.h"

void myAnalysis(std::string str_inf_X="", std::string str_inf_Y=""){

  //BGN--input data file name
  //str_inf_X="/Users/yamaguchinaoki/cernbox2/BeamLine/software/output/20220129_LED_noDacTune_F.root";
  //str_inf_Y="/Users/yamaguchinaoki/cernbox2/BeamLine/software/output/20220129_LED_noDacTune_F.root";
  std::string str_ouf  = "output.root";
  //END--input data file name

  std::cout << "inputFile(X):: " << str_inf_X << std::endl;
  std::cout << "inputFile(Y):: " << str_inf_Y << std::endl;

  //Main
  foutput = new TFile(Form("./output/%s",str_ouf.c_str()),"RECREATE");
  
  ///load file
  TFile *fin_X = new TFile(Form("%s",str_inf_X.c_str()), "read");
  if (!fin_X || !fin_X->IsOpen()){
    std::cout << "!!!ERROR!!! The input file(X) is not set correctly" << std::endl;
    return;
  }
  TFile *fin_Y = new TFile(Form("%s",str_inf_Y.c_str()), "read");
  if (!fin_Y || !fin_Y->IsOpen()){
    std::cout << "!!!ERROR!!! The input file(Y) is not set correctly" << std::endl;
    return;
  }
  TTree *tree_X = (TTree*) fin_X->Get("eventtree");
  TTree *tree_Y = (TTree*) fin_Y->Get("eventtree");
  TDirectoryFile *folder_X = (TDirectoryFile*)fin_X->Get("hists");
  TDirectoryFile *folder_Y = (TDirectoryFile*)fin_Y->Get("hists");
  
  ///peak Find
  ///init.
  nMaxPeak_X = 0, nMaxPeak_Y = 0;
  x_peak_X = 0.,  x_peak_Y = 0.;
  x_mean_X = 0.,  x_mean_Y = 0.;

  std::ofstream fout_X("./output/peakFind_out_X.dat");
  std::ofstream fout_Y("./output/peakFind_out_Y.dat");
  h_X = new TH1F();
  h_Y = new TH1F();
  
  TSpectrum *s_X = new TSpectrum(maxpeaks);
  TSpectrum *s_Y = new TSpectrum(maxpeaks);

  for(int i_count=0; i_count < maxpeaks; i_count++){
    f_X[i_count] = new TF1(Form("f_X%d",i_count), "gaus");
    f_Y[i_count] = new TF1(Form("f_Y%d",i_count), "gaus");
  }
  ///main
  c1 = new TCanvas("c1","c1");
  //BEGIN--Channel loop
  for(int i_Ch=0; i_Ch < max_Ch; i_Ch++){

    std::cout << "##################" << std::endl;
    h_X = (TH1F*)folder_X->Get(Form("ch%i", i_Ch));
    h_Y = (TH1F*)folder_Y->Get(Form("ch%i", i_Ch));
    h_X->SetTitle(Form("ch%i X;ADC;Counts", i_Ch));
    h_Y->SetTitle(Form("ch%i Y;ADC;Counts", i_Ch));
    
    std::cout << Form("ch%i",i_Ch) << std::endl;
    s_X->Search(h_X, tS_sigma, "", tS_threshold);
    s_Y->Search(h_Y, tS_sigma, "", tS_threshold);
    nMaxPeak_X = s_X->GetNPeaks();
    nMaxPeak_Y = s_Y->GetNPeaks();

    x_MaxPeak_X = s_X->GetPositionX()[0];
    x_MaxPeak_Y = s_Y->GetPositionX()[0];
    
    //BEGIN--Peak Loop
    for(int i_peak=0; i_peak < nMaxPeak_X; i_peak++){//for X
      x_peak_X = s_X->GetPositionX()[i_peak];
      if(x_peak_X < x_MaxPeak_X) continue;
      f_X[i_peak]->SetParameter(1, x_peak_X);
      h_X->Fit(Form("f_X%i", i_peak),"+","", x_peak_X-10, x_peak_X+10);
      x_mean_X = f_X[i_peak]->GetParameter(1);
      std::cout << "X: "<< i_Ch << "\t" << i_peak << "\t" << x_mean_X << std::endl;
      fout_X << i_Ch << "\t" << i_peak << "\t" << x_mean_X << std::endl;
    }
    h_X->GetXaxis()->SetRangeUser(700,2000);
    if(savePDF==true) c1->Print(Form("./output/x_ch%i.pdf",i_Ch));
    
    for(int i_peak=0; i_peak < nMaxPeak_Y; i_peak++){//for Y
      x_peak_Y = s_Y->GetPositionX()[i_peak];
      if(x_peak_Y < x_MaxPeak_Y) continue;
      f_Y[i_peak]->SetParameter(1, x_peak_Y);
      h_Y->Fit(Form("f_Y%i", i_peak),"+","", x_peak_Y-10, x_peak_Y+10);
      x_mean_Y = f_Y[i_peak]->GetParameter(1);
      std::cout << "Y: " << i_Ch << "\t" << i_peak << "\t" << x_mean_Y << std::endl;
      fout_Y << i_Ch << "\t" << i_peak << "\t" << x_mean_Y << std::endl;
    }    
    //END--Peak Loop    
   
    c1->SetLogy();

    h_Y->GetXaxis()->SetRangeUser(700,2000);
    if(savePDF==true) c1->Print(Form("./output/y_ch%i.pdf",i_Ch));

    foutput->cd();
    h_X->Write();
    h_Y->Write();
  }
   fout_X.close();    
   fout_Y.close();    
  //END--Channel loop
   
   
   c1->SetLogy(false);
   std::ifstream fin2_X("././output/peakFind_out_X.dat");
   std::ifstream fin2_Y("././output/peakFind_out_Y.dat");
   std::ofstream fout2_X("./output/ch_1pe_X.dat");
   std::ofstream fout2_Y("./output/ch_1pe_Y.dat");

  pre_Ch = -1;  
  while (fin2_X >> i_Ch >> i_peak >> mean){
    sum_mean = 0.;
    if(pre_Ch==i_Ch){//same Ch(MPPC)
      v_mean.push_back(mean);
    }
    else if (pre_Ch==-1){//1st Ch(MPPC)
      v_mean.clear();
      v_mean.push_back(mean);
      pre_Ch = i_Ch;
    }
    else{//next Ch(MPPC) (2 or more)
      diff_mean = (*std::max_element(v_mean.begin(), v_mean.end())) - (*std::min_element(v_mean.begin(), v_mean.end()));

      //std::cout << "Max:: " <<*std::max_element(v_mean.begin(), v_mean.end()) << " Min:: " <<*std::min_element(v_mean.begin(), v_mean.end()) << std::endl;
      if(v_mean.size()==1) avg_mean = (diff_mean)/(v_mean.size());
      else                 avg_mean = (diff_mean)/(v_mean.size()-1);
      
      std::cout << diff_mean << "/" << v_mean.size()-1 << " Average:: " << avg_mean << std::endl;
      h_1pe_mppc_X->Fill(avg_mean);
      std::cout << "XX" << std::endl;
      h_1pe_Ch_X->Fill(pre_Ch,avg_mean);
      fout2_X << pre_Ch << "\t" << avg_mean << "\t" << *std::min_element(v_mean.begin(), v_mean.end())  << std::endl;
      std::cout << "***************" << std::endl;
      v_mean.clear();
      v_mean.push_back(mean);
      pre_Ch = i_Ch;
    }
    std::cout << "X: " <<  i_Ch << "\t" << i_peak << "\t" << mean << std::endl;
  }
  //BGN-- for ch.63
  diff_mean = (*std::max_element(v_mean.begin(), v_mean.end())) - (*std::min_element(v_mean.begin(), v_mean.end()));
  
  if(v_mean.size()==1) avg_mean = (diff_mean)/(v_mean.size());
  else                 avg_mean = (diff_mean)/(v_mean.size()-1);
  std::cout << diff_mean << "/" << v_mean.size()-1 << " Average:: " << avg_mean << std::endl;
  h_1pe_mppc_X->Fill(avg_mean);
  h_1pe_Ch_X->Fill(pre_Ch,avg_mean);
  std::cout << "XX" << std::endl;
  fout2_X << pre_Ch << "\t" << avg_mean << "\t" << *std::min_element(v_mean.begin(), v_mean.end())  << std::endl;
  std::cout << "***************" << std::endl;
  v_mean.clear();
  std::cout << "X: " <<  i_Ch << "\t" << i_peak << "\t" << mean << std::endl;
  //END-- for ch.63
  
  fout2_X.close();
  std::cout << "#########################" << std::endl;
  
  pre_Ch = -1;  
  while (fin2_Y >> i_Ch >> i_peak >> mean){
    sum_mean = 0.;
    if(pre_Ch==i_Ch){//same Ch(MPPC)
      v_mean.push_back(mean);
    }
    else if (pre_Ch==-1){//1st Ch(MPPC)
      v_mean.clear();
      v_mean.push_back(mean);
      pre_Ch = i_Ch;
    }
    else{//next Ch(MPPC) (2 or more)
      diff_mean = (*std::max_element(v_mean.begin(), v_mean.end())) - (*std::min_element(v_mean.begin(), v_mean.end()));

      if(v_mean.size()==1) avg_mean = (diff_mean)/(v_mean.size());
      else                 avg_mean = (diff_mean)/(v_mean.size()-1);
      
      std::cout << diff_mean << "/" << v_mean.size()-1 << " Average:: " << avg_mean << std::endl;
      h_1pe_mppc_Y->Fill(avg_mean);
      h_1pe_Ch_Y->Fill(pre_Ch,avg_mean);
      
      fout2_Y << pre_Ch << "\t" << avg_mean << "\t" << *std::min_element(v_mean.begin(), v_mean.end())  << std::endl;
      
      std::cout << "***************" << std::endl;
      v_mean.clear();
      v_mean.push_back(mean);
      pre_Ch = i_Ch;
    }
    std::cout << "Y: " << i_Ch << "\t" << i_peak << "\t" << mean << std::endl;
  }
  //BGN-- for ch.63-Y
  diff_mean = (*std::max_element(v_mean.begin(), v_mean.end())) - (*std::min_element(v_mean.begin(), v_mean.end()));
  
  if(v_mean.size()==1) avg_mean = (diff_mean)/(v_mean.size());
  else                 avg_mean = (diff_mean)/(v_mean.size()-1);
  std::cout << diff_mean << "/" << v_mean.size()-1 << " Average:: " << avg_mean << std::endl;
  h_1pe_mppc_Y->Fill(avg_mean);
  h_1pe_Ch_Y->Fill(pre_Ch,avg_mean);
  fout2_Y << pre_Ch << "\t" << avg_mean << "\t" << *std::min_element(v_mean.begin(), v_mean.end())  << std::endl;
  std::cout << "***************" << std::endl;
  v_mean.clear();
  std::cout << "Y: " <<  i_Ch << "\t" << i_peak << "\t" << mean << std::endl;
  //END-- for ch.63-Y
  
  fout2_Y.close();
  std::cout << "#########################" << std::endl;
  //h_1pe_mppc->Draw();

  foutput->cd();
  h_1pe_mppc_X->Write();
  h_1pe_mppc_Y->Write();
  h_1pe_Ch_X->Write();
  h_1pe_Ch_Y->Write();

  //Calibrator //wp
 
  std::ifstream fin3_X("./output/ch_1pe_X.dat");
  std::ifstream fin3_Y("./output/ch_1pe_Y.dat");
  std::ifstream fin4_X("./output/peakFind_X.dat");
  std::ifstream fin4_Y("./output/peakFind_Y.dat");
  
  int maxBin =0;//ADC max bin

  tmp_i_Ch=0, tmp_i_peak=0, nBin=0;
  tmp_mean =0.0, pe_ADC=0., pe_Energy=0., x_adc_0pe=0.;

  for (int count_i=0; count_i<max_Ch; count_i++){//init.
    adc_1pe_X[count_i] = -1.0;
    adc_1pe_Y[count_i] = -1.0;
    adc_0pe_X[count_i] = -1.0;
    adc_0pe_Y[count_i] = -1.0;
  }
  
  while (fin3_X >> i_Ch >> pe_ADC >> x_adc_0pe){
    //calibFactor = pe_Energy/pe_ADC;
    std::cout << "X: " << i_Ch << "\t" << pe_ADC << "\t" << x_adc_0pe << std::endl;
    adc_1pe_X[i_Ch] = pe_ADC;
    adc_0pe_X[i_Ch] = x_adc_0pe;
  }
  while (fin3_Y >> i_Ch >> pe_ADC >> x_adc_0pe){
    //calibFactor = pe_Energy/pe_ADC;
    std::cout << "Y: " << i_Ch << "\t" << pe_ADC << "\t" << x_adc_0pe << std::endl;
    adc_1pe_Y[i_Ch] = pe_ADC;
    adc_0pe_Y[i_Ch] = x_adc_0pe;
  }
  
  /*    
	for (int count_i=0; count_i<64; count_i++){
	std::cout <<  adc_1pe[count_i]  << " " <<  adc_0pe[count_i] << std::endl;
	}
  */
  /*
    if (i_Ch == 0){//test region      
    
    nBin = 350/calibFactor;
    h_eV  = new TH1F("h_ev","chX; photons;Counts",nBin,0,350);
    h = (TH1F*)folder->Get("ch0");
    
    while(fin4 >> tmp_i_Ch >> tmp_i_peak >> tmp_mean){
    if (tmp_i_Ch == 0 && tmp_i_peak == 0){
    maxBin = tmp_mean;
    break;
    }
    }      
    std::cout << maxBin << "\t" << h->GetBinContent(maxBin) << std::endl;
    for(int j_count=0; j_count < 3200; j_count++){
    h_eV->Fill(calibFactor*0.5+calibFactor*j_count,h->GetBinContent(maxBin+j_count));
    }      
    }
  */
  //  h_eV->Draw();
  
  ///hit Find
  ///init.
  
  max_adc_X=-1, max_fiberCh_X=-1;
  max_adc_Y=-1, max_fiberCh_Y=-1;
  nEvent=0;
  
  for (int cnt_i=0; cnt_i < max_Ch; cnt_i++){
    adc_X[cnt_i] = -1;
    adc_Y[cnt_i] = -1;
    photon_X[cnt_i] = -1;
    photon_Y[cnt_i] = -1;
  }
  hitmap = new TH2F("hitmap","hitmap;X;Y",64,0,64,64,0,64);
  hitX   = new TH1F("hitx", "hit X",64,0,64);
  hitY   = new TH1F("hity", "hit Y",64,0,64);
  
  ///main
  tree_X->SetBranchAddress("adc_ch", &adc_X);
  tree_Y->SetBranchAddress("adc_ch", &adc_Y);
  
  if(tree_X->GetEntries() != tree_Y->GetEntries()) {
    std::cout << "X and Y file are incorrelct, number of entry is different" << std::endl;
  }
  const int Nentry = tree_X->GetEntries();
  std::cout << "Entry:: " <<  Nentry << std::endl;
  
  for(int ientry=0; ientry < Nentry; ientry++){//event loop
    tree_X->GetEntry(ientry);
    tree_Y->GetEntry(ientry);
    
    //convert[adc->Nphoton]
    for (int cnt_i=0; cnt_i<max_Ch; cnt_i++){
      photon_X[cnt_i] = (adc_X[cnt_i] - adc_0pe_X[cnt_i])/(adc_1pe_X[cnt_i]);
      photon_Y[cnt_i] = (adc_Y[cnt_i] - adc_0pe_Y[cnt_i])/(adc_1pe_Y[cnt_i]);
    }
    
    max_adc_X = *std::max_element(adc_X, adc_X+64);
    max_adc_Y = *std::max_element(adc_Y, adc_Y+64);
    
    max_photon_X = *std::max_element(photon_X, photon_X+64);
    max_photon_Y = *std::max_element(photon_Y, photon_Y+64);
    
    max_fiberCh_X = std::max_element(photon_X, photon_X+64) - photon_X;
    max_fiberCh_Y = std::max_element(photon_Y, photon_Y+64) - photon_Y;

    hitmap->Fill(max_fiberCh_X+0.5, max_fiberCh_Y+0.5);
    hitX->Fill(max_fiberCh_X+0.5);
    hitY->Fill(max_fiberCh_Y+0.5);
  }

  hitmap->Draw("colz");
  foutput->cd();
  hitmap->Write();
  hitX->Write();
  hitY->Write(); 
  //foutput->Close();   
}
