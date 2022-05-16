ch39->Draw();
gPad->SetLogy();
ch39->GetXaxis()->UnZoom();
ch39->Rebin(10);
ch39->Fit("gaus","","",1300,3000);
gStyle->SetOptFit();
