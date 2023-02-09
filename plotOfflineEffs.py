import ROOT
from setTDRStyle import setTDRStyle

ROOT.gStyle.SetTitleXSize(0.05) 
ROOT.gStyle.SetTitleYOffset(1.25) 
ROOT.gStyle.SetTitleXOffset(1.0) 
ROOT.gStyle.SetTitleYSize(0.05)
ROOT.gStyle.SetLabelColor(1, "XYZ")
ROOT.gStyle.SetLabelFont(42, "XYZ")
ROOT.gStyle.SetLabelOffset(0.005, "XYZ")
ROOT.gStyle.SetLabelSize(0.035, "XYZ")



f= ROOT.TFile("offlineEffs.root","READ")

plots = ["Pt", "Eta", "Phi"]

rangesX = {"Pt":[0,10],"Eta":[-2.4,2.4],"Phi":[-3.14,3.14]}
rangesY = {"Pt":[0,1.4],"Eta":[0.,1.4],"Phi":[0.,1.4]}

labelsX = {"Pt": "p_{T}^{gen #mu} [GeV]", "Eta": "#eta^{gen #mu}", "Phi": "#phi^{gen #mu}"}

for plot in plots:

	effMu = f.Get("muon%s"%plot)
	effSoftMu = f.Get("softMuon%s"%plot)
	effLooseMu = f.Get("looseMuon%s"%plot)
	effTrack = f.Get("trackEff%s"%plot)

	effMu.SetMarkerColor(ROOT.kRed)
	effSoftMu.SetMarkerColor(ROOT.kBlue)
	effLooseMu.SetMarkerColor(ROOT.kGreen+2)
	effTrack.SetMarkerColor(ROOT.kBlack)
	effMu.SetLineColor(ROOT.kRed)
	effSoftMu.SetLineColor(ROOT.kBlue)
	effLooseMu.SetLineColor(ROOT.kGreen+2)
	effTrack.SetLineColor(ROOT.kBlack)

	effMu.SetMarkerStyle(20)
	effSoftMu.SetMarkerStyle(21)
	effLooseMu.SetMarkerStyle(23)
	effTrack.SetMarkerStyle(22)


	latex = ROOT.TLatex()
	latex.SetTextFont(42)
	latex.SetTextAlign(31)
	latex.SetTextSize(0.04)
	latex.SetNDC(True)
	latexCMS = ROOT.TLatex()
	latexCMS.SetTextFont(61)
	latexCMS.SetTextSize(0.055)
	latexCMS.SetNDC(True)
	latexCMSExtra = ROOT.TLatex()
	latexCMSExtra.SetTextFont(52)
	latexCMSExtra.SetTextSize(0.03)
	latexCMSExtra.SetNDC(True) 

	latexInfo = ROOT.TLatex()
	latexInfo.SetTextFont(42)
	latexInfo.SetTextSize(0.05)
	latexInfo.SetNDC(True) 


	c1= ROOT.TCanvas("c1","c1",800,800)
	c1.cd()
	plotPad = ROOT.TPad("plotPad","plotPad",0.01,0.01,0.99,0.99)

	setTDRStyle()		
	plotPad.UseCurrentStyle()
	plotPad.Draw()	
	plotPad.cd()

	plotPad.DrawFrame(rangesX[plot][0],rangesY[plot][0],rangesX[plot][1],rangesY[plot][1],"; %s; Efficiency"%labelsX[plot])

	leg = ROOT.TLegend(0.2586957,0.754878,0.7341137,0.9320557,"","brNDC")
	leg.SetTextFont(62)
	leg.SetLineColor(1)
	leg.SetLineStyle(1)
	leg.SetLineWidth(3)
	leg.SetFillColor(0)
	leg.SetFillStyle(1001)
	leg.SetShadowColor(0)
	# ~ leg.SetDrawOption(0)
	leg.SetBorderSize(0)
	leg.SetTextSize(0.04)
	 
	leg.AddEntry(effMu,"reco muon",'l')
	leg.AddEntry(effSoftMu,"reco muon + soft ID",'l')
	leg.AddEntry(effLooseMu,"reco muon + loose ID",'l')
	leg.AddEntry(effTrack,"reco track",'l')


	latex.DrawLatex(0.95, 0.96, "(14 TeV)")
	latexCMS.DrawLatex(0.15,0.955,"CMS")
	latexCMSExtra.DrawLatex(0.275,0.955,"Simulation Preliminary")

	leg.Draw()

	effMu.Draw("same")
	effSoftMu.Draw("same")
	effLooseMu.Draw("same")
	effTrack.Draw("same")

	ROOT.gPad.RedrawAxis()
	 
	c1.SaveAs("eff_%s.pdf"%plot)
	c1.SaveAs("eff_%s.png"%plot)



mMu = f.Get("m3Mu")
mSoftMu = f.Get("m3MuSoft")
mLooseMu = f.Get("m3MuLoose")
mTrack = f.Get("m3Track")

mMu.SetMarkerColor(ROOT.kRed)
mSoftMu.SetMarkerColor(ROOT.kBlue)
mLooseMu.SetMarkerColor(ROOT.kGreen+2)
mTrack.SetMarkerColor(ROOT.kBlack)
mMu.SetLineColor(ROOT.kRed)
mSoftMu.SetLineColor(ROOT.kBlue)
mLooseMu.SetLineColor(ROOT.kGreen+2)
mTrack.SetLineColor(ROOT.kBlack)

mMu.SetMarkerStyle(20)
mSoftMu.SetMarkerStyle(21)
mLooseMu.SetMarkerStyle(23)
mTrack.SetMarkerStyle(22)


latex = ROOT.TLatex()
latex.SetTextFont(42)
latex.SetTextAlign(31)
latex.SetTextSize(0.04)
latex.SetNDC(True)
latexCMS = ROOT.TLatex()
latexCMS.SetTextFont(61)
latexCMS.SetTextSize(0.055)
latexCMS.SetNDC(True)
latexCMSExtra = ROOT.TLatex()
latexCMSExtra.SetTextFont(52)
latexCMSExtra.SetTextSize(0.03)
latexCMSExtra.SetNDC(True) 

latexInfo = ROOT.TLatex()
latexInfo.SetTextFont(42)
latexInfo.SetTextSize(0.05)
latexInfo.SetNDC(True) 


c1= ROOT.TCanvas("c1","c1",800,800)
c1.cd()
plotPad = ROOT.TPad("plotPad","plotPad",0.01,0.01,0.99,0.99)

setTDRStyle()		
plotPad.UseCurrentStyle()
plotPad.Draw()	
plotPad.cd()

plotPad.DrawFrame(0,0,3,40000,"; M_{3#mu} [GeV]; Events")

leg = ROOT.TLegend(0.2586957,0.754878,0.7341137,0.9320557,"","brNDC")
leg.SetTextFont(62)
leg.SetLineColor(1)
leg.SetLineStyle(1)
leg.SetLineWidth(3)
leg.SetFillColor(0)
leg.SetFillStyle(1001)
leg.SetShadowColor(0)
# ~ leg.SetDrawOption(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)
	
leg.AddEntry(mMu,"reco muon",'l')
leg.AddEntry(mSoftMu,"reco muon + soft ID",'l')
leg.AddEntry(mLooseMu,"reco muon + loose ID",'l')
leg.AddEntry(mTrack,"reco track",'l')


latex.DrawLatex(0.95, 0.96, "(14 TeV)")
latexCMS.DrawLatex(0.15,0.955,"CMS")
latexCMSExtra.DrawLatex(0.275,0.955,"Simulation Preliminary")

leg.Draw()

mMu.Draw("same")
mSoftMu.Draw("same")
mLooseMu.Draw("same")
mTrack.Draw("same")

ROOT.gPad.RedrawAxis()
	
c1.SaveAs("m3Mu.pdf")
c1.SaveAs("m3Mu.png")



mTrack = f.Get("m3Track")

mTrack.SetMarkerColor(ROOT.kBlack)
mTrack.SetLineColor(ROOT.kBlack)
mTrack.SetMarkerStyle(22)


latex = ROOT.TLatex()
latex.SetTextFont(42)
latex.SetTextAlign(31)
latex.SetTextSize(0.04)
latex.SetNDC(True)
latexCMS = ROOT.TLatex()
latexCMS.SetTextFont(61)
latexCMS.SetTextSize(0.055)
latexCMS.SetNDC(True)
latexCMSExtra = ROOT.TLatex()
latexCMSExtra.SetTextFont(52)
latexCMSExtra.SetTextSize(0.03)
latexCMSExtra.SetNDC(True) 

latexInfo = ROOT.TLatex()
latexInfo.SetTextFont(42)
latexInfo.SetTextSize(0.05)
latexInfo.SetNDC(True) 


c1= ROOT.TCanvas("c1","c1",800,800)
c1.cd()
plotPad = ROOT.TPad("plotPad","plotPad",0.01,0.01,0.99,0.99)

setTDRStyle()		
plotPad.UseCurrentStyle()
plotPad.Draw()	
plotPad.cd()

plotPad.DrawFrame(0,0,3,5e5,"; M_{3#mu} [GeV]; Events")

leg = ROOT.TLegend(0.2586957,0.754878,0.7341137,0.9320557,"","brNDC")
leg.SetTextFont(62)
leg.SetLineColor(1)
leg.SetLineStyle(1)
leg.SetLineWidth(3)
leg.SetFillColor(0)
leg.SetFillStyle(1001)
leg.SetShadowColor(0)
# ~ leg.SetDrawOption(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)
	
leg.AddEntry(mTrack,"reco track",'l')


latex.DrawLatex(0.95, 0.96, "(14 TeV)")
latexCMS.DrawLatex(0.15,0.955,"CMS")
latexCMSExtra.DrawLatex(0.275,0.955,"Simulation Preliminary")

leg.Draw()

mTrack.Draw("same")

ROOT.gPad.RedrawAxis()
	
c1.SaveAs("m3Tracks.pdf")
c1.SaveAs("m3Tracks.png")




mMu = f.Get("nGenMuon")
mSoftMu = f.Get("nRecoMuon")
mLooseMu = f.Get("nRecoMuonMatched")

mMu.SetMarkerColor(ROOT.kRed)
mSoftMu.SetMarkerColor(ROOT.kBlue)
mLooseMu.SetMarkerColor(ROOT.kGreen+2)
mMu.SetLineColor(ROOT.kRed)
mSoftMu.SetLineColor(ROOT.kBlue)
mLooseMu.SetLineColor(ROOT.kGreen+2)

mMu.SetMarkerStyle(20)
mSoftMu.SetMarkerStyle(21)
mLooseMu.SetMarkerStyle(23)


latex = ROOT.TLatex()
latex.SetTextFont(42)
latex.SetTextAlign(31)
latex.SetTextSize(0.04)
latex.SetNDC(True)
latexCMS = ROOT.TLatex()
latexCMS.SetTextFont(61)
latexCMS.SetTextSize(0.055)
latexCMS.SetNDC(True)
latexCMSExtra = ROOT.TLatex()
latexCMSExtra.SetTextFont(52)
latexCMSExtra.SetTextSize(0.03)
latexCMSExtra.SetNDC(True) 

latexInfo = ROOT.TLatex()
latexInfo.SetTextFont(42)
latexInfo.SetTextSize(0.05)
latexInfo.SetNDC(True) 


c1= ROOT.TCanvas("c1","c1",800,800)
c1.cd()
plotPad = ROOT.TPad("plotPad","plotPad",0.01,0.01,0.99,0.99)

setTDRStyle()		
plotPad.UseCurrentStyle()
plotPad.Draw()	
plotPad.cd()

plotPad.DrawFrame(-0.5,0.5,500.5,5000000,"; N^{#mu}; Events")
plotPad.SetLogy()
leg = ROOT.TLegend(0.2586957,0.754878,0.7341137,0.9320557,"","brNDC")
leg.SetTextFont(62)
leg.SetLineColor(1)
leg.SetLineStyle(1)
leg.SetLineWidth(3)
leg.SetFillColor(0)
leg.SetFillStyle(1001)
leg.SetShadowColor(0)
# ~ leg.SetDrawOption(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)
	
leg.AddEntry(mMu,"gen muon",'l')
leg.AddEntry(mSoftMu,"reco muon",'l')
leg.AddEntry(mLooseMu,"reco muon matched to gen",'l')


latex.DrawLatex(0.95, 0.96, "(14 TeV)")
latexCMS.DrawLatex(0.15,0.955,"CMS")
latexCMSExtra.DrawLatex(0.275,0.955,"Simulation Preliminary")

leg.Draw()

mMu.Draw("same")
mSoftMu.Draw("same")
mLooseMu.Draw("same")

ROOT.gPad.RedrawAxis()
	
c1.SaveAs("nMu.pdf")
c1.SaveAs("nMu.png")






mMu = f.Get("nGenMuon")
mSoftMu = f.Get("nRecoMuon")
mLooseMu = f.Get("nRecoMuonMatched")

mMu.SetMarkerColor(ROOT.kRed)
mSoftMu.SetMarkerColor(ROOT.kBlue)
mLooseMu.SetMarkerColor(ROOT.kGreen+2)
mMu.SetLineColor(ROOT.kRed)
mSoftMu.SetLineColor(ROOT.kBlue)
mLooseMu.SetLineColor(ROOT.kGreen+2)

mMu.SetMarkerStyle(20)
mSoftMu.SetMarkerStyle(21)
mLooseMu.SetMarkerStyle(23)


latex = ROOT.TLatex()
latex.SetTextFont(42)
latex.SetTextAlign(31)
latex.SetTextSize(0.04)
latex.SetNDC(True)
latexCMS = ROOT.TLatex()
latexCMS.SetTextFont(61)
latexCMS.SetTextSize(0.055)
latexCMS.SetNDC(True)
latexCMSExtra = ROOT.TLatex()
latexCMSExtra.SetTextFont(52)
latexCMSExtra.SetTextSize(0.03)
latexCMSExtra.SetNDC(True) 

latexInfo = ROOT.TLatex()
latexInfo.SetTextFont(42)
latexInfo.SetTextSize(0.05)
latexInfo.SetNDC(True) 


c1= ROOT.TCanvas("c1","c1",800,800)
c1.cd()
plotPad = ROOT.TPad("plotPad","plotPad",0.01,0.01,0.99,0.99)

setTDRStyle()		
plotPad.UseCurrentStyle()
plotPad.Draw()	
plotPad.cd()

plotPad.DrawFrame(-0.5,0.5,50.5,5000000,"; N^{#mu}; Events")
plotPad.SetLogy()
leg = ROOT.TLegend(0.2586957,0.754878,0.7341137,0.9320557,"","brNDC")
leg.SetTextFont(62)
leg.SetLineColor(1)
leg.SetLineStyle(1)
leg.SetLineWidth(3)
leg.SetFillColor(0)
leg.SetFillStyle(1001)
leg.SetShadowColor(0)
# ~ leg.SetDrawOption(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)
	
leg.AddEntry(mMu,"gen muon",'l')
leg.AddEntry(mSoftMu,"reco muon",'l')
leg.AddEntry(mLooseMu,"reco muon matched to gen",'l')


latex.DrawLatex(0.95, 0.96, "(14 TeV)")
latexCMS.DrawLatex(0.15,0.955,"CMS")
latexCMSExtra.DrawLatex(0.275,0.955,"Simulation Preliminary")

leg.Draw()

mMu.Draw("same")
mSoftMu.Draw("same")
mLooseMu.Draw("same")

ROOT.gPad.RedrawAxis()
	
c1.SaveAs("nMu.pdf")
c1.SaveAs("nMu.png")





pTGen1 = f.Get("etaGen1")
pTGen2 = f.Get("etaGen2")
pTGen3 = f.Get("etaGen3")

pTGen1Matched = f.Get("etaGenMatched1")
pTGen2Matched = f.Get("etaGenMatched2")
pTGen3Matched = f.Get("etaGenMatched3")


pTGen1.SetMarkerColor(ROOT.kRed)
pTGen1Matched.SetMarkerColor(ROOT.kRed)
pTGen1.SetLineColor(ROOT.kRed)
pTGen1Matched.SetLineColor(ROOT.kRed)
pTGen2.SetMarkerColor(ROOT.kBlue)
pTGen2Matched.SetMarkerColor(ROOT.kBlue)
pTGen2.SetLineColor(ROOT.kBlue)
pTGen2Matched.SetLineColor(ROOT.kBlue)
pTGen3.SetMarkerColor(ROOT.kGreen+2)
pTGen3Matched.SetMarkerColor(ROOT.kGreen+2)
pTGen3.SetLineColor(ROOT.kGreen+2)
pTGen3Matched.SetLineColor(ROOT.kGreen+2)

pTGen1.SetMarkerStyle(20)
pTGen2.SetMarkerStyle(21)
pTGen3.SetMarkerStyle(22)
pTGen1Matched.SetMarkerStyle(23)
pTGen2Matched.SetMarkerStyle(24)
pTGen3Matched.SetMarkerStyle(25)

pTGen1Matched.SetLineStyle(ROOT.kDashed)
pTGen2Matched.SetLineStyle(ROOT.kDashed)
pTGen3Matched.SetLineStyle(ROOT.kDashed)


latex = ROOT.TLatex()
latex.SetTextFont(42)
latex.SetTextAlign(31)
latex.SetTextSize(0.04)
latex.SetNDC(True)
latexCMS = ROOT.TLatex()
latexCMS.SetTextFont(61)
latexCMS.SetTextSize(0.055)
latexCMS.SetNDC(True)
latexCMSExtra = ROOT.TLatex()
latexCMSExtra.SetTextFont(52)
latexCMSExtra.SetTextSize(0.03)
latexCMSExtra.SetNDC(True) 

latexInfo = ROOT.TLatex()
latexInfo.SetTextFont(42)
latexInfo.SetTextSize(0.05)
latexInfo.SetNDC(True) 


c1= ROOT.TCanvas("c1","c1",800,800)
c1.cd()
plotPad = ROOT.TPad("plotPad","plotPad",0.01,0.01,0.99,0.99)

setTDRStyle()		
plotPad.UseCurrentStyle()
plotPad.Draw()	
plotPad.cd()

plotPad.DrawFrame(-2.4,0,2.4,100000,"; #eta^{gen #mu}; Events")
# plotPad.SetLogy()
leg = ROOT.TLegend(0.2586957,0.754878,0.7341137,0.9320557,"","brNDC")
leg.SetTextFont(62)
leg.SetLineColor(1)
leg.SetLineStyle(1)
leg.SetLineWidth(3)
leg.SetFillColor(0)
leg.SetFillStyle(1001)
leg.SetShadowColor(0)
# ~ leg.SetDrawOption(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)
	
leg.AddEntry(pTGen1,"leading muon",'l')
leg.AddEntry(pTGen1Matched,"leading muon reco matched",'l')
leg.AddEntry(pTGen2,"subleading muon",'l')
leg.AddEntry(pTGen2Matched,"subleading muon reco matched",'l')
leg.AddEntry(pTGen3,"trailing muon",'l')
leg.AddEntry(pTGen3Matched,"trailing muon reco matched",'l')


latex.DrawLatex(0.95, 0.96, "(14 TeV)")
latexCMS.DrawLatex(0.15,0.955,"CMS")
latexCMSExtra.DrawLatex(0.275,0.955,"Simulation Preliminary")

leg.Draw()

pTGen1.Draw("samehist")
pTGen1Matched.Draw("samehist")
pTGen2.Draw("samehist")
pTGen2Matched.Draw("samehist")
pTGen3.Draw("samehist")
pTGen3Matched.Draw("samehist")

ROOT.gPad.RedrawAxis()
	
c1.SaveAs("genEta.pdf")
c1.SaveAs("genEta.png")


pTGen1 = f.Get("pt1Gen")
pTGen2 = f.Get("pt2Gen")
pTGen3 = f.Get("pt3Gen")

pTGen1Matched = f.Get("ptGenMatched1")
pTGen2Matched = f.Get("ptGen2Matched")
pTGen3Matched = f.Get("ptGenMatched3")


pTGen1.SetMarkerColor(ROOT.kRed)
pTGen1Matched.SetMarkerColor(ROOT.kRed)
pTGen1.SetLineColor(ROOT.kRed)
pTGen1Matched.SetLineColor(ROOT.kRed)
pTGen2.SetMarkerColor(ROOT.kBlue)
pTGen2Matched.SetMarkerColor(ROOT.kBlue)
pTGen2.SetLineColor(ROOT.kBlue)
pTGen2Matched.SetLineColor(ROOT.kBlue)
pTGen3.SetMarkerColor(ROOT.kGreen+2)
pTGen3Matched.SetMarkerColor(ROOT.kGreen+2)
pTGen3.SetLineColor(ROOT.kGreen+2)
pTGen3Matched.SetLineColor(ROOT.kGreen+2)

pTGen1.SetMarkerStyle(20)
pTGen2.SetMarkerStyle(21)
pTGen3.SetMarkerStyle(22)
pTGen1Matched.SetMarkerStyle(23)
pTGen2Matched.SetMarkerStyle(24)
pTGen3Matched.SetMarkerStyle(25)

pTGen1Matched.SetLineStyle(ROOT.kDashed)
pTGen2Matched.SetLineStyle(ROOT.kDashed)
pTGen3Matched.SetLineStyle(ROOT.kDashed)


latex = ROOT.TLatex()
latex.SetTextFont(42)
latex.SetTextAlign(31)
latex.SetTextSize(0.04)
latex.SetNDC(True)
latexCMS = ROOT.TLatex()
latexCMS.SetTextFont(61)
latexCMS.SetTextSize(0.055)
latexCMS.SetNDC(True)
latexCMSExtra = ROOT.TLatex()
latexCMSExtra.SetTextFont(52)
latexCMSExtra.SetTextSize(0.03)
latexCMSExtra.SetNDC(True) 

latexInfo = ROOT.TLatex()
latexInfo.SetTextFont(42)
latexInfo.SetTextSize(0.05)
latexInfo.SetNDC(True) 


c1= ROOT.TCanvas("c1","c1",800,800)
c1.cd()
plotPad = ROOT.TPad("plotPad","plotPad",0.01,0.01,0.99,0.99)

setTDRStyle()		
plotPad.UseCurrentStyle()
plotPad.Draw()	
plotPad.cd()

plotPad.DrawFrame(0,0.5,20,500000,"; p_{T}^{gen #mu}; Events")
plotPad.SetLogy()
leg = ROOT.TLegend(0.2586957,0.754878,0.7341137,0.9320557,"","brNDC")
leg.SetTextFont(62)
leg.SetLineColor(1)
leg.SetLineStyle(1)
leg.SetLineWidth(3)
leg.SetFillColor(0)
leg.SetFillStyle(1001)
leg.SetShadowColor(0)
# ~ leg.SetDrawOption(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)
	
leg.AddEntry(pTGen1,"leading muon",'l')
leg.AddEntry(pTGen1Matched,"leading muon reco matched",'l')
leg.AddEntry(pTGen2,"subleading muon",'l')
leg.AddEntry(pTGen2Matched,"subleading muon reco matched",'l')
leg.AddEntry(pTGen3,"trailing muon",'l')
leg.AddEntry(pTGen3Matched,"trailing muon reco matched",'l')


latex.DrawLatex(0.95, 0.96, "(14 TeV)")
latexCMS.DrawLatex(0.15,0.955,"CMS")
latexCMSExtra.DrawLatex(0.275,0.955,"Simulation Preliminary")

leg.Draw()

pTGen1.Draw("samehist")
pTGen1Matched.Draw("samehist")
pTGen2.Draw("samehist")
pTGen2Matched.Draw("samehist")
pTGen3.Draw("samehist")
pTGen3Matched.Draw("samehist")

ROOT.gPad.RedrawAxis()
	
c1.SaveAs("genPt.pdf")
c1.SaveAs("genPt.png")


pTGen1 = f.Get("pt1")
pTGen2 = f.Get("pt2")
pTGen3 = f.Get("pt3")

pTGen1Matched = f.Get("ptGenMatched1")
pTGen2Matched = f.Get("ptGen2Matched")
pTGen3Matched = f.Get("ptGenMatched3")


pTGen1.SetMarkerColor(ROOT.kRed)
pTGen1Matched.SetMarkerColor(ROOT.kRed)
pTGen1.SetLineColor(ROOT.kRed)
pTGen1Matched.SetLineColor(ROOT.kRed)
pTGen2.SetMarkerColor(ROOT.kBlue)
pTGen2Matched.SetMarkerColor(ROOT.kBlue)
pTGen2.SetLineColor(ROOT.kBlue)
pTGen2Matched.SetLineColor(ROOT.kBlue)
pTGen3.SetMarkerColor(ROOT.kGreen+2)
pTGen3Matched.SetMarkerColor(ROOT.kGreen+2)
pTGen3.SetLineColor(ROOT.kGreen+2)
pTGen3Matched.SetLineColor(ROOT.kGreen+2)

pTGen1.SetMarkerStyle(20)
pTGen2.SetMarkerStyle(21)
pTGen3.SetMarkerStyle(22)
pTGen1Matched.SetMarkerStyle(23)
pTGen2Matched.SetMarkerStyle(24)
pTGen3Matched.SetMarkerStyle(25)

pTGen1Matched.SetLineStyle(ROOT.kDashed)
pTGen2Matched.SetLineStyle(ROOT.kDashed)
pTGen3Matched.SetLineStyle(ROOT.kDashed)


latex = ROOT.TLatex()
latex.SetTextFont(42)
latex.SetTextAlign(31)
latex.SetTextSize(0.04)
latex.SetNDC(True)
latexCMS = ROOT.TLatex()
latexCMS.SetTextFont(61)
latexCMS.SetTextSize(0.055)
latexCMS.SetNDC(True)
latexCMSExtra = ROOT.TLatex()
latexCMSExtra.SetTextFont(52)
latexCMSExtra.SetTextSize(0.03)
latexCMSExtra.SetNDC(True) 

latexInfo = ROOT.TLatex()
latexInfo.SetTextFont(42)
latexInfo.SetTextSize(0.05)
latexInfo.SetNDC(True) 


c1= ROOT.TCanvas("c1","c1",800,800)
c1.cd()
plotPad = ROOT.TPad("plotPad","plotPad",0.01,0.01,0.99,0.99)

setTDRStyle()		
plotPad.UseCurrentStyle()
plotPad.Draw()	
plotPad.cd()

plotPad.DrawFrame(0,0.5,20,500000,"; p_{T}^{reco #mu}; Events")
plotPad.SetLogy()
leg = ROOT.TLegend(0.2586957,0.754878,0.7341137,0.9320557,"","brNDC")
leg.SetTextFont(62)
leg.SetLineColor(1)
leg.SetLineStyle(1)
leg.SetLineWidth(3)
leg.SetFillColor(0)
leg.SetFillStyle(1001)
leg.SetShadowColor(0)
# ~ leg.SetDrawOption(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)
	
leg.AddEntry(pTGen1,"leading muon",'l')
# leg.AddEntry(pTGen1Matched,"leading muon reco matched",'l')
leg.AddEntry(pTGen2,"subleading muon",'l')
# leg.AddEntry(pTGen2Matched,"subleading muon reco matched",'l')
leg.AddEntry(pTGen3,"trailing muon",'l')
# leg.AddEntry(pTGen3Matched,"trailing muon reco matched",'l')


latex.DrawLatex(0.95, 0.96, "(14 TeV)")
latexCMS.DrawLatex(0.15,0.955,"CMS")
latexCMSExtra.DrawLatex(0.275,0.955,"Simulation Preliminary")

leg.Draw()

pTGen1.Draw("samehist")
# pTGen1Matched.Draw("samehist")
pTGen2.Draw("samehist")
# pTGen2Matched.Draw("samehist")
pTGen3.Draw("samehist")
# pTGen3Matched.Draw("samehist")

ROOT.gPad.RedrawAxis()
	
c1.SaveAs("recoPt.pdf")
c1.SaveAs("recoPt.png")