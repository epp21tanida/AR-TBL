import ROOT
import pandas as pd

import argparse

# input argument parsing
usage = "usage: %prog [options]"
parser = argparse.ArgumentParser(usage)
parser.add_argument("-a", "--ana", action="store_true",  dest="ana", default=False )
parser.add_argument("-f", "--fit", action="store_true",  dest="fit", default=False )
parser.add_argument("-n", "--nor", action="store_true",  dest="nor", default=False )
options = parser.parse_args()

ANA = options.ana
FIT = options.fit
NOR = options.nor

# global ROOT settings
ROOT.gStyle.SetOptFit(1)

# global paths
outDir     = "./output/"
gainPath   = outDir + "datfile/peakFind_out_X."
outputFile = outDir + "mppc_gain/20220512_2_Sr_1000000_af."

dataDir  = "../../../dataset/"
dataVer  = "20220512_2_Sr_1000000"
dataCh   = "_af"
rootFile = dataDir + dataVer + dataCh + ".root"

def convert_DAT_CSV(path):
	"""
	To convert .dat files into .csv files for ease of use later when using pandas
	args:
	- path : path towards target file
	"""


	# Read .dat file and change into lists
	datContent = [[float(x) for x in i.strip().split()] for i in open(path + "dat").readlines()]

	# col = ["ch", "ipeak", "adc"]
	col = ["a", "b"]

	# save the contents as a csv file via pandas
	pd.DataFrame(datContent, columns = col).to_csv(path + "csv", index=False)

def Normalize():
	"""
	To normalize the ADC channel histograms using extracted pedestal and gain values
	"""

	# num of channels to normalize
	nCh = 64

	# dict declarations
	mip   = {}
	hbare = {}
	hnm   = {}	
	pedgain = {"ped":[], "gain":[]}

	# dummy canvas for htemp plotting
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

	# Extract gain values
	i = 0
	for idx, row in gain_df.iterrows():
		if row["ipeak"]	== 0: pedgain["ped"].append(row["adc"])
		else:
			pedgain["gain"].append(round(row["adc"] - pedgain["ped"][i], 3))
			i += 1

	# Update pedestal value with highest peak
	for i in range(nCh):
		htemp = ROOT.TH1F("htemp"+str(i), "", 4000, 0, 4000)
		inTree.Draw("adc_ch["+str(i)+"]>>htemp"+str(i))
		pedgain["ped"][i] = htemp.GetBinCenter(htemp.GetMaximumBin())

	# canvas to plot pre-normalize histos and post-normalize	
	canv = ROOT.TCanvas("canv", "", 1600, 600)
	canv.Divide(2,1)
	
	# set logY plots
	p1 = canv.cd(1)
	p2 = canv.cd(2)
	p1.SetLogy()
	p2.SetLogy()

	outFileName = outDir+"GainCalib_norm.pdf"
	canv.Print(outFileName+"[")

	# save normalized histos for separate fitting
	fitHistoFile = ROOT.TFile("./output/GainCalib_fitHistos.root", "RECREATE")
	
	for j in range(nCh):

		hbare[j] = ROOT.TH1F("hbare"+str(j), "", 4000, 0,    4000)
		hnm[j]   = ROOT.TH1F("hnm"  +str(j), "", 82,   -1.5, 80.5)

		# pre-normalize
		canv.cd(1)
		inTree.Draw("adc_ch["+str(j)+"]>>hbare"+str(j))
			
		# post-normalize	
		canv.cd(2)
		inTree.Draw("(adc_ch["+str(j)+"]-"+str(pedgain["ped"][j])+")/"+str(pedgain["gain"][j])+">>hnm"+str(j))
		canv.Print(outFileName)
		hnm[j].Write()

	canv.Print(outFileName+"]")
	fitHistoFile.Close()

def Fit():
	"""
	To fit the normalized histos with a gaussian (+ gaussian) function to see the mean value
	Ideally the mean values for all channels are within a range
	"""

	# silence ROOT
	ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

	nCh = 64

	# dict declarations
	rdhDict  = {}
	histDict = {}
	frmDict  = {}
	pdfDict  = {}

	# input file
	inFile = ROOT.TFile("./output/GainCalib_fitHistos.root", "OPEN")

	# base variable for roodatahist import
	xvar  = ROOT.RooRealVar("x", "p.e", -6, 82)
	xvar.setRange("fit", 10, 50)

	c = ROOT.TCanvas("c", "", 800, 600)
	c.SetLogy()

	# latex for in-plot text
	t = ROOT.TLatex()
	t.SetNDC()
	t.SetTextSize(0.03)
	t.SetTextColor(1)
	t.SetTextFont(42)
	t.SetTextAlign(33)

	outFileName = outDir+"GainCalib_fit.pdf"
	c.Print(outFileName+"[")

	for i in range(nCh):

		si = str(i)

		pdfDict[i] = {}
		pdfDict[i]["g1"]   = {}
		pdfDict[i]["g2"]   = {}
		pdfDict[i]["exp"]  = {}
		pdfDict[i]["comp"] = {}

		# setup frame
		frmDict[i]  = xvar.frame()
		
		# import histo
		histDict[i] = inFile.Get("hnm"+si).Clone()
		rdhDict[i]  = ROOT.RooDataHist("rdh"+str(i), "rdh"+si, ROOT.RooArgList(xvar), ROOT.RooFit.Import(histDict[i]))
		
		# first gaussian
		pdfDict[i]["g1"]["mean" +si] = ROOT.RooRealVar("mean1" +si, "mean1" +si, 20, 15, 40)
		pdfDict[i]["g1"]["sigma"+si] = ROOT.RooRealVar("sigma1"+si, "sigma1"+si, 10, 3, 40)
		pdfDict[i]["g1"]["pdf"  +si] = ROOT.RooGaussian("g1"   +si, "gaus1" +si, xvar, pdfDict[i]["g1"]["mean"+si], pdfDict[i]["g1"]["sigma"+si] )

		# second gaussian
		pdfDict[i]["g2"]["mean" +si] = ROOT.RooRealVar("mean2" +si, "mean2" +si, 40, 30, 50)
		pdfDict[i]["g2"]["sigma"+si] = ROOT.RooRealVar("sigma2"+si, "sigma2"+si, 15, 3, 50)
		pdfDict[i]["g2"]["pdf"  +si] = ROOT.RooGaussian("g2"   +si, "gaus2" +si, xvar, pdfDict[i]["g2"]["mean"+si], pdfDict[i]["g2"]["sigma"+si] )

		# WIP: complementary exponential
		pdfDict[i]["exp"]["const"+si] = ROOT.RooRealVar(  "const"+si, "const"+si, -1, -5, -0.001)
		pdfDict[i]["exp"]["pdf"  +si] = ROOT.RooExponential("exp"+si, "exp"  +si, xvar, pdfDict[i]["exp"]["const"+si])

		# PDF weightage
		pdfDict[i]["comp"]["1w"+si] = ROOT.RooRealVar("1wgt"+si, "1wgt"+si, 0.5, 0.01,  1.0)
		pdfDict[i]["comp"]["2w"+si] = ROOT.RooRealVar("2wgt"+si, "2wgt"+si, 0.2, 0.01,  0.4)		
		pdfDict[i]["comp"]["wl"]    = ROOT.RooArgList(pdfDict[i]["comp"]["1w"+si])

		# composite PDF
		pdfDict[i]["comp"]["pl"]    = ROOT.RooArgList(pdfDict[i]["g1"]["pdf"+si], pdfDict[i]["g2"]["pdf"+si])
		pdfDict[i]["comp"]["cp"+si] = ROOT.RooAddPdf("comp" +si, "comp"+si, pdfDict[i]["comp"]["pl"], pdfDict[i]["comp"]["wl"])

		pdf = pdfDict[i]["comp"]["cp"+si]

		# fitting & save fit result
		fr = pdf.fitTo(rdhDict[i], ROOT.RooFit.Save(1), ROOT.RooFit.PrintLevel(-1), ROOT.RooFit.Range("fit"))
		fr.Print()	

		# plot histo and fitted pdf
		rdhDict[i].plotOn(frmDict[i])
		pdf.plotOn(frmDict[i])

		print("chi2"+si+" = ", frmDict[i].chiSquare(5))

		frmDict[i].SetTitle("adc_ch"+si)
		frmDict[i].SetMinimum(1e-1)
		frmDict[i].SetMaximum(1e7)
		frmDict[i].Draw()

		t.DrawLatex(0.85, 0.85, "chi2_"+si+" = "+str(round(frmDict[i].chiSquare(5), 6)))
		t.DrawLatex(0.85, 0.80, "mean_"+si+" = "+str(round(pdfDict[i]["g1"]["mean" +si].getVal(), 6)))

		c.Print(outFileName)

	c.Print(outFileName+"]")

def AnalyzeFit():
	"""
	To analyze the fit results to see the distribution of the variables
	"""

	def cleanup(line):

		words = ["RooFit", "covariance", "Status", "Floating", "---------"]

		return True if not (any(word in line for word in words) and line) else False


	import inspect
	def retrieve_name(var):
	    callers_local_vars = inspect.currentframe().f_back.f_locals.items()
	    return [var_name for var_name, var_val in callers_local_vars if var_val is var]

	# load fit results 
	inFile = outDir + "GainCalib_fitResult.txt"
	lines  = [x for x in open(inFile, "r").readlines() if cleanup(x)]

	# separate fit results according to variable
	wgt1   = [[float(i) for i in x.replace("+/-", "").replace("1wgt",   "").split()][1] for x in lines if "1wgt"   in x]
	mean1  = [[float(i) for i in x.replace("+/-", "").replace("mean1",  "").split()][1] for x in lines if "mean1"  in x]
	mean2  = [[float(i) for i in x.replace("+/-", "").replace("mean2",  "").split()][1] for x in lines if "mean2"  in x]
	sigma1 = [[float(i) for i in x.replace("+/-", "").replace("sigma1", "").split()][1] for x in lines if "sigma1" in x]
	sigma2 = [[float(i) for i in x.replace("+/-", "").replace("sigma2", "").split()][1] for x in lines if "sigma2" in x]
	chi2   = [[float(i) for i in x.replace("=",   "").replace("chi2",   "").split()][1] for x in lines if "chi2"   in x]
	
	canv = ROOT.TCanvas("canv", "", 800, 600)

	outFileName = outDir + "GainCalib_varAna.pdf"
	canv.Print(outFileName + "[")

	# plot variables on a graph
	for varList in [wgt1, mean1, sigma1, mean2, sigma2, chi2]:

		gr = ROOT.TGraph(len(varList))

		for idx, var in enumerate(varList):
			gr.SetPoint(idx, idx, var)

		gr.SetTitle(retrieve_name(varList)[0])
		gr.SetMarkerStyle(20)
		gr.SetMarkerSize(1)
		gr.SetMinimum(min(varList)-0.5)
		gr.SetMaximum(max(varList)*1.5)
		gr.Draw("ap")
		canv.Print(outFileName)

	canv.Print(outFileName + "]")

if __name__ == "__main__":

	timer = ROOT.TStopwatch()
	timer.Start()
	
	if NOR: Normalize()
	if FIT: Fit()
	if ANA: AnalyzeFit()

	# convert_DAT_CSV(gainPath)

	timer.Print()	
	print("****** DONE ******")
