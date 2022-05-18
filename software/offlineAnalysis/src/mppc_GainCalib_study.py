import ROOT
import pandas as pd

ROOT.gStyle.SetOptFit(1)

outDir = "./output/"
gainPath   = outDir + "datfile/peakFind_out_X."
outputFile = outDir + "mppc_gain/20220512_2_Sr_1000000_af."

dataDir  = "../../../dataset/"
dataVer  = "20220512_2_Sr_1000000"
dataCh   = "_af"
rootFile = dataDir + dataVer + dataCh + ".root"

def convert_DAT_CSV(path):

	# Read .dat file and change into lists
	datContent = [[float(x) for x in i.strip().split()] for i in open(path + "dat").readlines()]

	# col = ["ch", "ipeak", "adc"]
	col = ["a", "b"]

	# save the contents as a csv file via pandas
	pd.DataFrame(datContent, columns = col).to_csv(path + "csv", index=False)

def main():

	nCh = 64
	mip = {}
	pedgain = {"ped":[], "gain":[]}


	c = ROOT.TCanvas("c", "", 800, 600)

	# INPUT FILE
	inFile = ROOT.TFile(rootFile, "OPEN")
	if not inFile.IsOpen():
		print("Error: File is not opened")
		return 1

	# INPUT TREE
	inTree = inFile.Get("eventtree")
	if not inTree:
		print("Error: Tree is not found")
		return 1

	# PEDESTAL & GAIN FILE
	gain_df = pd.read_csv(gainPath+"csv").astype({"ch":"int","ipeak":"int"})
	if gain_df.empty:
		print("Error: Empty dataframe")
		return 1

	i = 0
	for idx, row in gain_df.iterrows():
		if row["ipeak"]	== 0: pedgain["ped"].append(row["adc"])
		else:
			pedgain["gain"].append(round(row["adc"] - pedgain["ped"][i], 3))
			i += 1


	for i in range(nCh):
		htemp = ROOT.TH1F("htemp"+str(i), "", 4000, 0, 4000)
		inTree.Draw("adc_ch["+str(i)+"]>>htemp"+str(i))
		pedgain["ped"][i] = htemp.GetBinCenter(htemp.GetMaximumBin())

	
	
	canv = ROOT.TCanvas("canv", "", 1600, 600)
	canv.Divide(2,1)
	
	p1 = canv.cd(1)
	p2 = canv.cd(2)
	p1.SetLogy()
	p2.SetLogy()

	# c.Print(outputFile + ".pdf[")
	canv.Print("test.pdf[")

	hbare = {}
	hnm   = {}

	fitHistoFile = ROOT.TFile("./output/FitHistos.root", "RECREATE")
	
	for j in range(nCh):

		hbare[j] = ROOT.TH1F("hbare"+str(j), "", 4000, 0,    4000)
		hnm[j]   = ROOT.TH1F("hnm"+str(j),   "", 82,   -1.5, 80.5)

		canv.cd(1)
		inTree.Draw("adc_ch["+str(j)+"]>>hbare"+str(j))
			
		canv.cd(2)
		inTree.Draw("(adc_ch["+str(j)+"]-"+str(pedgain["ped"][j])+")/"+str(pedgain["gain"][j])+">>hnm"+str(j))
		canv.Print("test.pdf")
		hnm[j].Write()

	canv.Print("test.pdf]")

	fitHistoFile.Close()

def FitGausExp():

	nCh = 64

	rdhDict  = {}
	histDict = {}
	frmDict  = {}
	pdfDict  = {}
	pdfDict["gaus"] = {}
	pdfDict["exp"]  = {}
	pdfDict["comp"] = {}

	inFile = ROOT.TFile("./output/FitHistos.root", "OPEN")

	xvar  = ROOT.RooRealVar("x", "p.e", -2, 82)
	xvar.setRange("fit", 5, 50)

	c = ROOT.TCanvas("c", "", 800, 600)
	c.SetLogy()

	c.Print("test2.pdf[")

	for i in range(nCh):

		si = str(i)

		pdfDict[i] = {}
		pdfDict[i]["g1"] = {}
		pdfDict[i]["g2"] = {}
		pdfDict[i]["g3"] = {}
		pdfDict[i]["exp"]  = {}
		pdfDict[i]["comp"] = {}

		frmDict[i]  = xvar.frame()
		histDict[i] = inFile.Get("hnm"+si).Clone()
		rdhDict[i]  = ROOT.RooDataHist("rdh"+str(i), "rdh"+si, ROOT.RooArgList(xvar), ROOT.RooFit.Import(histDict[i]))
		
		pdfDict[i]["g1"]["mean" +si] = ROOT.RooRealVar("mean1" +si, "mean1" +si, 20, 10, 35)
		pdfDict[i]["g1"]["sigma"+si] = ROOT.RooRealVar("sigma1"+si, "sigma1"+si, 5, 3, 20)
		pdfDict[i]["g1"]["pdf"  +si] = ROOT.RooGaussian("g1"   +si, "gaus1" +si, xvar, pdfDict[i]["g1"]["mean"+si], pdfDict[i]["g1"]["sigma"+si] )

		pdfDict[i]["g2"]["mean" +si] = ROOT.RooRealVar("mean2" +si, "mean2" +si, 30, 20, 50)
		pdfDict[i]["g2"]["sigma"+si] = ROOT.RooRealVar("sigma2"+si, "sigma2"+si, 15, 3, 40)
		pdfDict[i]["g2"]["pdf"  +si] = ROOT.RooGaussian("g2"   +si, "gaus2" +si, xvar, pdfDict[i]["g2"]["mean"+si], pdfDict[i]["g2"]["sigma"+si] )

		# pdfDict[i]["g3"]["mean" +si] = ROOT.RooRealVar("mean3" +si, "mean3" +si, 40, 20, 60)
		# pdfDict[i]["g3"]["sigma"+si] = ROOT.RooRealVar("sigma3"+si, "sigma3"+si, 10, 3, 20)
		# pdfDict[i]["g3"]["pdf"  +si] = ROOT.RooGaussian("g3"   +si, "gaus3" +si, xvar, pdfDict[i]["g3"]["mean"+si], pdfDict[i]["g3"]["sigma"+si] )

		pdfDict[i]["exp"]["const"+si] = ROOT.RooRealVar(  "const"+si, "const"+si, -1, -5, -0.001)
		pdfDict[i]["exp"]["pdf"  +si] = ROOT.RooExponential("exp"+si, "exp"  +si, xvar, pdfDict[i]["exp"]["const"+si])

		pdfDict[i]["comp"]["1w"+si] = ROOT.RooRealVar("1wgt"+si, "1wgt"+si, 0.5, 0.01,  0.7)
		pdfDict[i]["comp"]["2w"+si] = ROOT.RooRealVar("2wgt"+si, "2wgt"+si, 0.2, 0.01,  0.4)
		# pdfDict[i]["comp"]["3w"+si] = ROOT.RooRealVar("3wgt"+si, "3wgt"+si, 0.1, 0.05, 0.4)
		pdfDict[i]["comp"]["wl"]    = ROOT.RooArgList(pdfDict[i]["comp"]["1w"+si], pdfDict[i]["comp"]["2w"+si])

		pdfDict[i]["comp"]["pl"]    = ROOT.RooArgList(pdfDict[i]["g1"]["pdf"+si], pdfDict[i]["g2"]["pdf"+si], pdfDict[i]["exp"]["pdf"+si])
		pdfDict[i]["comp"]["cp"+si] = ROOT.RooAddPdf("comp" +si, "comp"+si, pdfDict[i]["comp"]["pl"], pdfDict[i]["comp"]["wl"])

		fr = pdfDict[i]["comp"]["cp"+str(i)].fitTo(rdhDict[i], ROOT.RooFit.Save(1), ROOT.RooFit.PrintLevel(-1), ROOT.RooFit.Range("fit"))
		fr.Print()	

		rdhDict[i].plotOn(frmDict[i])
		pdfDict[i]["comp"]["cp"+si].plotOn(frmDict[i])

		frmDict[i].SetTitle("adc_ch"+si)
		frmDict[i].SetMinimum(1e-1)
		frmDict[i].SetMaximum(1e7)
		frmDict[i].Draw()
		c.Print("test2.pdf")

	c.Print("test2.pdf]")



if __name__ == "__main__":

	timer = ROOT.TStopwatch()
	timer.Start()
	
	# main()

	FitGausExp()

	# convert_DAT_CSV(gainPath)

	timer.Print()	
	print("****** DONE ******")
