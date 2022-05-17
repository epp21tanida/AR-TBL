#include "myAnalysis.h"

void myAnalysis(std::string str_inf_X="", std::string str_inf_Y=""){

  /////////////////////////////////
  // FILE I/O
  /////////////////////////////////

  // OUTPUT FILE
  foutput = new TFile("./output/output.root","RECREATE");
  
  // INPUT FILE
  std::cout << "inputFile(X):: " << str_inf_X << std::endl;
  std::cout << "inputFile(Y):: " << str_inf_Y << std::endl;

  TFile *fin_X = new TFile(Form("%s",str_inf_X.c_str()), "READ");
  if (!fin_X || !fin_X->IsOpen()){
    std::cout << "ERROR! The input file(X) is not set correctly" << std::endl;
    return;
  }

  TFile *fin_Y = new TFile(Form("%s",str_inf_Y.c_str()), "READ");
  if (!fin_Y || !fin_Y->IsOpen()){
    std::cout << "ERROR! The input file(Y) is not set correctly" << std::endl;
    return;
  }

  // RETRIEVE TTREE
  TTree *tree_X = (TTree*) fin_X->Get("eventtree");
  TTree *tree_Y = (TTree*) fin_Y->Get("eventtree");
  
  TDirectoryFile *folder_X = (TDirectoryFile*)fin_X->Get("hists");
  TDirectoryFile *folder_Y = (TDirectoryFile*)fin_Y->Get("hists");
  
  /////////////////////////////////
  // PEAK ANALYSIS
  // We scan the histograms using the TSpectrum::Search method and find the first two
  // peaks in the distributions. The position of the peaks along the x-axis is then used
  // to set the initial mean value for a TF1 Gaussian to be later fitted within a range
  // of the peaks. The subsequent fitted mean value is then saved, to be later used.
  /////////////////////////////////
  
  // INITIALIZATION VALUES
  nMaxPeak_X = 0,   nMaxPeak_Y = 0;
  x_peak_X   = 0.,  x_peak_Y   = 0.;
  x_mean_X   = 0.,  x_mean_Y   = 0.;

  std::ofstream fout_X("./output/peakFind_out_X.dat");
  std::ofstream fout_Y("./output/peakFind_out_Y.dat");
  
  // DECLARE HISTO
  h_X = new TH1F();
  h_Y = new TH1F();
  
  // DECLARE SPECTRUM (maxpeaks = 2)
  TSpectrum *s_X = new TSpectrum(maxpeaks);
  TSpectrum *s_Y = new TSpectrum(maxpeaks);

  // DECLARE TF1 GAUSSIANS FOR FITTING
  for(int i_count=0; i_count < maxpeaks; i_count++){
    f_X[i_count] = new TF1(Form("f_X%d",i_count), "gaus");
    f_Y[i_count] = new TF1(Form("f_Y%d",i_count), "gaus");
  }

  c1 = new TCanvas("c1","c1");
  
  // LOOP BY CHANNEL
  for(int i_Ch=0; i_Ch < max_Ch; i_Ch++){

    std::cout << "##################" << std::endl;

    // RETRIEVE INPUT HISTOS
    h_X = (TH1F*)folder_X->Get(Form("ch%i", i_Ch));
    h_Y = (TH1F*)folder_Y->Get(Form("ch%i", i_Ch));

    h_X->SetTitle(Form("ch%i X;ADC;Counts", i_Ch));
    h_Y->SetTitle(Form("ch%i Y;ADC;Counts", i_Ch));
    
    std::cout << Form("ch%i",i_Ch) << std::endl;

    // TSPECTRUM::SEARCH (sigma = 2, threshold = 0.0001)
    s_X->Search(h_X, tS_sigma, "", tS_threshold);
    s_Y->Search(h_Y, tS_sigma, "", tS_threshold);
    
    nMaxPeak_X = s_X->GetNPeaks();
    nMaxPeak_Y = s_Y->GetNPeaks();

    x_MaxPeak_X = s_X->GetPositionX()[0];
    x_MaxPeak_Y = s_Y->GetPositionX()[0];
    
    // X-AXIS FIBERS
    for(int i_peak=0; i_peak < nMaxPeak_X; i_peak++){//for X
      
      x_peak_X = s_X->GetPositionX()[i_peak];
      
      if(x_peak_X < x_MaxPeak_X) continue;
      
      // TF1 GAUSSIAN FIT
      f_X[i_peak]->SetParameter(1, x_peak_X);
      h_X->Fit(Form("f_X%i", i_peak),"+","", x_peak_X-10, x_peak_X+10);
      
      // GAUSSIAN POST-FIT MEAN VALUE
      x_mean_X = f_X[i_peak]->GetParameter(1);
      
      std::cout << "X: "<< i_Ch << "\t" << i_peak << "\t" << x_mean_X << std::endl;
      
      fout_X << i_Ch << "\t" << i_peak << "\t" << x_mean_X << std::endl;
    }
    
    // PLOT X
    h_X->GetXaxis()->SetRangeUser(700,2000);
    if(savePDF==true) c1->Print("./output/Form(x_ch%i.pdf,i_Ch)");
    
    // Y-AXIS FIBERS
    for(int i_peak=0; i_peak < nMaxPeak_Y; i_peak++){//for Y
      
      x_peak_Y = s_Y->GetPositionX()[i_peak];
      
      if(x_peak_Y < x_MaxPeak_Y) continue;
      
      // TF1 GAUSSIAN FIT
      f_Y[i_peak]->SetParameter(1, x_peak_Y);
      h_Y->Fit(Form("f_Y%i", i_peak),"+","", x_peak_Y-10, x_peak_Y+10);
      
      // GAUSSIAN POST-FIT MEAN VALUE
      x_mean_Y = f_Y[i_peak]->GetParameter(1);
      
      std::cout << "Y: " << i_Ch << "\t" << i_peak << "\t" << x_mean_Y << std::endl;
      
      fout_Y << i_Ch << "\t" << i_peak << "\t" << x_mean_Y << std::endl;
    }    
        
    // PLOT Y
    c1->SetLogy();
    h_Y->GetXaxis()->SetRangeUser(700,2000);
    if(savePDF==true) c1->Print("./output/Form(y_ch%i.pdf,i_Ch)");

    foutput->cd();
    h_X->Write();
    h_Y->Write();
  }

  fout_X.close();    
  fout_Y.close();    
  
  c1->SetLogy(false);
   
  /////////////////////////////////
  // MEAN ANALYSIS
  // With the previously obtained fitted mean values for the first two peaks,
  // we calculate the amplitude difference between the peaks and save the lowest
  // one of the two for the next step.
  /////////////////////////////////

  std::ifstream fin2_X("./output/peakFind_out_X.dat");
  std::ifstream fin2_Y("./output/peakFind_out_Y.dat");
  std::ofstream fout2_X("./ch_1pe_X.dat");
  std::ofstream fout2_Y("./ch_1pe_Y.dat");

  // X-AXIS FIBERS
  pre_Ch = -1;  
  while (fin2_X >> i_Ch >> i_peak >> mean){
    
    sum_mean = 0.;
    
    // same Ch(MPPC)
    if(pre_Ch==i_Ch){ v_mean.push_back(mean);}
    
    // 1st Ch(MPPC)
    else if (pre_Ch==-1){
      
      v_mean.clear();
      v_mean.push_back(mean);
      
      pre_Ch = i_Ch;
    }

    // next Ch(MPPC) (2 or more)    
    else{

      diff_mean = (*std::max_element(v_mean.begin(), v_mean.end())) - (*std::min_element(v_mean.begin(), v_mean.end()));
      
      if(v_mean.size()==1) avg_mean = (diff_mean)/(v_mean.size());
      else                 avg_mean = (diff_mean)/(v_mean.size()-1);
      
      std::cout << diff_mean << "/" << v_mean.size()-1 << " Average:: " << avg_mean << std::endl;      
      std::cout << "XX" << std::endl;
      
      h_1pe_mppc_X->Fill(avg_mean);
      h_1pe_Ch_X->Fill(pre_Ch,avg_mean);
      
      fout2_X << pre_Ch << "\t" << avg_mean << "\t" << *std::min_element(v_mean.begin(), v_mean.end())  << std::endl;
      
      std::cout << "***************" << std::endl;
      
      v_mean.clear();
      v_mean.push_back(mean);
      
      pre_Ch = i_Ch;
    }
    
    std::cout << "X: " <<  i_Ch << "\t" << i_peak << "\t" << mean << std::endl;
  }

  // BGN X-- for ch.63
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
  
  fout2_X.close();
  
  std::cout << "#########################" << std::endl;
  

  // Y-AXIS FIBERS
  pre_Ch = -1;  
  while (fin2_Y >> i_Ch >> i_peak >> mean){
    
    sum_mean = 0.;
    
    // SAME Ch(MPPC)
    if(pre_Ch==i_Ch){
      v_mean.push_back(mean);
    }
    
    // 1st Ch(MPPC)
    else if (pre_Ch==-1){

      v_mean.clear();
      v_mean.push_back(mean);
      
      pre_Ch = i_Ch;
    }
    
    // NEXT Ch(MPPC) (2 or more)
    else{

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
  
  fout2_Y.close();
  std::cout << "#########################" << std::endl;
  
  //h_1pe_mppc->Draw();

  foutput->cd();

  h_1pe_mppc_X->Write();
  h_1pe_mppc_Y->Write();
  h_1pe_Ch_X->Write();
  h_1pe_Ch_Y->Write();

  // Calibrator // wp
 
  std::ifstream fin3_X("./ch_1pe_X.dat");
  std::ifstream fin3_Y("./ch_1pe_Y.dat");

  // std::ifstream fin4_X("./output/peakFind_X.dat");
  // std::ifstream fin4_Y("./output/peakFind_Y.dat");
  
  int maxBin = 0;//ADC max bin

  nBin = 0;
  tmp_i_Ch = 0,  tmp_i_peak = 0,  tmp_mean  = 0.;
  pe_ADC   = 0., pe_Energy  = 0., x_adc_0pe = 0.;

  // INITIALIZE ADC ARRAYS
  for (int count_i=0; count_i<max_Ch; count_i++){//init.
    
    adc_1pe_X[count_i] = -1.0;
    adc_1pe_Y[count_i] = -1.0;
    adc_0pe_X[count_i] = -1.0;
    adc_0pe_Y[count_i] = -1.0;
  }
  
  //calibFactor = pe_Energy/pe_ADC;
  
  while (fin3_X >> i_Ch >> pe_ADC >> x_adc_0pe){
    
    std::cout << "X: " << i_Ch << "\t" << pe_ADC << "\t" << x_adc_0pe << std::endl;
    
    adc_1pe_X[i_Ch] = pe_ADC;
    adc_0pe_X[i_Ch] = x_adc_0pe;
  }

  while (fin3_Y >> i_Ch >> pe_ADC >> x_adc_0pe){
    
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
  
  //////////////////////////
  // HITS ANALYSIS 
  //////////////////////////

  nEvent = 0;
  max_adc_X = -1, max_fiberCh_X = -1;
  max_adc_Y = -1, max_fiberCh_Y = -1;
  
  // INITIALIZE ARRAYS & HISTOS
  for (int cnt_i=0; cnt_i < max_Ch; cnt_i++){
    
    adc_X[cnt_i] = -1;
    adc_Y[cnt_i] = -1;
    
    photon_X[cnt_i] = -1;
    photon_Y[cnt_i] = -1;
  }

  // INITIALIZE HIT HISTOS
  hitmap = new TH2F("hitmap", "hitmap;X;Y", 64,0,64,64,0,64);
  hitX   = new TH1F("hitx",   "hit X",      64,0,64);
  hitY   = new TH1F("hity",   "hit Y",      64,0,64);
  
  // RETRIEVE TTREES
  tree_X->SetBranchAddress("adc_ch", &adc_X);
  tree_Y->SetBranchAddress("adc_ch", &adc_Y);
  
  if(tree_X->GetEntries() != tree_Y->GetEntries()) {
    std::cout << "X and Y file are incorrect, number of entry is different" << std::endl;
  }
  
  const int Nentry = tree_X->GetEntries();
  std::cout << "Entry:: " <<  Nentry << std::endl;
  
  // EVENT LOOP
  for(int ientry=0; ientry < Nentry; ientry++){
  
    tree_X->GetEntry(ientry);
    tree_Y->GetEntry(ientry);
    
    // CONVERT [adc->Nphoton]
    for (int cnt_i=0; cnt_i<max_Ch; cnt_i++){
  
      photon_X[cnt_i] = (adc_X[cnt_i] - adc_0pe_X[cnt_i])/(adc_1pe_X[cnt_i]);
      photon_Y[cnt_i] = (adc_Y[cnt_i] - adc_0pe_Y[cnt_i])/(adc_1pe_Y[cnt_i]);
    }
    
    // GET THE ADC WITH THE HIGHEST VALUE
    max_adc_X = *std::max_element(adc_X, adc_X+64);
    max_adc_Y = *std::max_element(adc_Y, adc_Y+64);
    
    max_photon_X = *std::max_element(photon_X, photon_X+64);
    max_photon_Y = *std::max_element(photon_Y, photon_Y+64);
    
    max_fiberCh_X = std::max_element(photon_X, photon_X+64) - photon_X;
    max_fiberCh_Y = std::max_element(photon_Y, photon_Y+64) - photon_Y;

    // FILL HIT MAP ACCORDING TO FIBER POSITION
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
