#include <fstream>
#include <iostream>
#include <string>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <cmath>
#include <TF1.h>
#include <stdlib.h>

#include <TROOT.h>
#include <TFile.h>
#include <TString.h>
#include <vector>
#include <TH2.h>
#include <map>
#include <math.h>
#include <TTree.h>
#include <algorithm>
#include <cassert>
#include <utility>
#include <TSpectrum.h>

//Const. val.
//Config
//BGN--Configuration
bool savePDF = true;

//END--Configuration

///peak find
//const int   maxpeaks = 10;
const int   maxpeaks = 2;
const float tS_sigma = 2;//about TSpectrum
const float tS_threshold = 0.0001;//about TSpectrum

//int  adc_pedestal = 770;//tmp
int  testMode = 1;//special for debug

//hit find
const int   max_Ch = 64;


//variables
TFile *foutput;

///peak find
int   nMaxPeak_X, nMaxPeak_Y;
float x_peak_X, x_mean_X, x_MaxPeak_X;
float x_peak_Y, x_mean_Y, x_MaxPeak_Y;
TH1F  *h_X;
TH1F  *h_Y;
TF1   *f_X[maxpeaks];
TF1   *f_Y[maxpeaks];
TCanvas *c1;


///photo electron convertor
int i_Ch=0, i_peak=0, pre_Ch=-1;
float mean=0., avg_mean=-1., sum_mean=0., diff_mean=0.;
std::vector<float> v_mean;
TH1F *h_1pe_mppc_X = new TH1F("h_1pe_mppc_X", "1pe distribution X",200,20,40);
TH1F *h_1pe_mppc_Y = new TH1F("h_1pe_mppc_Y", "1pe distribution Y",200,20,40);
TH1F *h_1pe_Ch_X = new TH1F("h_1pe_ch_X", "1pe distribution X;Ch",64,0,64);
TH1F *h_1pe_Ch_Y = new TH1F("h_1pe_ch_Y", "1pe distribution Y;Ch",64,0,64);

//TH1F *h_eV;

int tmp_i_Ch, tmp_i_peak, nBin;
double tmp_mean, pe_ADC, pe_Energy, x_adc_0pe;

double  adc_1pe_X[max_Ch],adc_1pe_Y[max_Ch];
double  adc_0pe_X[max_Ch],adc_0pe_Y[max_Ch];

///hit Find
int adc_X[max_Ch], adc_Y[max_Ch];
int photon_X[max_Ch], photon_Y[max_Ch];

int max_adc_X,max_adc_Y; 
int max_photon_X,max_photon_Y; 
int max_fiberCh_X, max_fiberCh_Y; 
int nEvent;

TH2F* hitmap;
TH1F* hitX;
TH1F* hitY;
