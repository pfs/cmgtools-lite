import ROOT
import os
import sys

dcDir='/afs/cern.ch/work/m/mdunser/public/cmssw/heavyIons/CMSSW_9_4_6_patch1/src/CMGTools/TTHAnalysis/python/plotter/hin-ttbar/datacards_2019-10-21_forApproval_jetAnalysis/bdtcomb'
pfix='btagged'

if len(sys.argv)>1:
    dcDir=sys.argv[1]
    pfix=sys.argv[2]

ptitle='leptonic' if pfix=='leptonic' else 'leptonic+b-tagged'

inF=ROOT.TFile.Open(os.path.join(dcDir,'combinedCard/fitDiagnosticsobs.root'))

allBins=[]
for ch in inF.Get('shapes_fit_s').GetListOfKeys():
    chDir=ch.ReadObj()
    
    tot_sig=None
    tot_bkg=None
    data=None
    for proc in chDir.GetListOfKeys():
        pname=proc.GetName()
        plot=proc.ReadObj()
        if 'total_signal' in pname:     tot_sig=proc.ReadObj()
        if 'total_background' in pname: tot_bkg=proc.ReadObj()
        if pname=='data':               data=proc.ReadObj()

    x,y=ROOT.Double(0),ROOT.Double(0)
    for xbin in range(1,tot_sig.GetNbinsX()+1):

        try:
            y=data.GetBinContent(xbin)
            print xbin,y
        except:
            data.GetPoint(xbin-1,x,y)

        nbkg=tot_bkg.GetBinContent(xbin)
        sob=tot_sig.GetBinContent(xbin)/nbkg if nbkg>0 else -1

        allBins.append( [sob,
                         tot_sig.GetBinContent(xbin),
                         tot_bkg.GetBinContent(xbin),
                         tot_bkg.GetBinError(xbin),
                         float(y),
                         (ch.GetName(),xbin) ] )

allsob=[ ROOT.TMath.Log10(x[0]) for x in allBins if x[0]>=0. ]
sobTempl=ROOT.TH1F('sob',';log_{10}(S/B);Events',6,min(allsob),max(allsob))
sobH={}
for key in ['sig','bkg','bkgsumw2','data']:
    sobH[key]=sobTempl.Clone('sob'+key)
sobH['data'].SetBinErrorOption(ROOT.TH1.kPoisson)

for sob,sig,bkg,bkgunc,data,_ in allBins:
    if sob<0 : continue
    sob=ROOT.TMath.Log10(sob)
    sobH['sig'].Fill(sob,sig)
    xbin=sobH['bkg'].Fill(sob,bkg)
    sobH['bkgsumw2'].SetBinContent(xbin,sobH['bkgsumw2'].GetBinContent(xbin)+bkgunc**2)
    for i in range(int(data)):
        sobH['data'].Fill(sob)

#finalize background uncertainty
for xbin in range(1,sobH['bkg'].GetNbinsX()+1):
    finalUnc=ROOT.TMath.Sqrt(sobH['bkgsumw2'].GetBinContent(xbin))
    sobH['bkg'].SetBinError(xbin,finalUnc)


ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)
c = ROOT.TCanvas("c","c",600,600)
c.SetTopMargin(0)
c.SetLeftMargin(0)
c.SetRightMargin(0)
c.SetBottomMargin(0)

#data/MC
p1 = ROOT.TPad("p1", "p1", 0., 0.35, 1., 1.)
p1.SetTopMargin(0.1)
p1.SetRightMargin(0.03)
p1.SetLeftMargin(0.12)
p1.SetBottomMargin(0.03)
p1.Draw()
p1.SetLogy()

p1.cd()
sobTempl.GetYaxis().SetTitle('Events')
sobTempl.GetYaxis().SetTitleOffset(0.85)
sobTempl.GetYaxis().SetNdivisions(5)
sobTempl.GetYaxis().SetTitleSize(0.07)
sobTempl.GetYaxis().SetLabelSize(0.07)
sobTempl.GetXaxis().SetTitleSize(0)
sobTempl.GetXaxis().SetLabelSize(0)
sobTempl.GetYaxis().SetRangeUser(1,1e4) #5*sobH['bkg'].GetMaximum())
sobTempl.Draw()

stack = ROOT.THStack("stack","stack")
sobH['sig'].SetFillStyle(1001)
sobH['sig'].SetFillColor(ROOT.kRed+2)
sobH['sig'].SetLineColor(1)
sobH['bkg'].SetFillStyle(1001)
sobH['bkg'].SetFillColor(ROOT.kGray)
sobH['bkg'].SetLineColor(1)
stack.Add(sobH['bkg'])
stack.Add(sobH['sig'])
stack.Draw('histsame')

buncBand=sobH['bkg'].Clone('bunc')
buncBand.SetFillStyle(3244)
buncBand.SetFillColor(1)
buncBand.SetLineColor(1)
buncBand.SetMarkerStyle(1)
buncBand.Draw('e2same')

#totalUnc=ROOT.TGraphErrors(plotsPostfit['total'])
#totalUnc.SetFillStyle(3444)
#totalUnc.SetFillColor(1)
#totalUnc.SetMarkerStyle(1)
#totalUnc.Draw('e2')

sobH['data'].SetLineColor(1)
sobH['data'].SetLineWidth(2)
sobH['data'].SetMarkerColor(1)
sobH['data'].SetMarkerStyle(20)
sobH['data'].Draw('X0e1same')

leg = ROOT.TLegend(0.62,0.55,0.97,0.88)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.08)
leg.AddEntry(sobH['data'],'Data','ep')
leg.AddEntry(sobH['sig'],'t#bar{t}','f')
leg.AddEntry(sobH['bkg'],'Bkg.','f')
leg.AddEntry(buncBand,'Bkg. unc.','f')
leg.Draw('same')

tex=ROOT.TLatex()
tex.SetTextFont(42)
tex.SetTextSize(0.08)
tex.SetNDC()
tex.SetTextAlign(ROOT.kHAlignLeft+ROOT.kVAlignCenter)
tex.DrawLatex(0.12,0.95,'#bf{CMS} #it{Preliminary}')
tex.DrawLatex(0.15,0.85,'#scale[0.8]{#it{%s}}'%ptitle)
tex.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
tex.DrawLatex(0.95,0.95,'#scale[0.8]{1.7 nb^{-1} (#sqrt{s_{NN}}=5.02 TeV)}')

p1.RedrawAxis()

#prefit/postfit
c.cd()

p2 = ROOT.TPad("p2", "p1", 0., 0., 1., 0.35)
p2.SetTopMargin(0.01)
p2.SetRightMargin(0.03)
p2.SetLeftMargin(0.12)
p2.SetBottomMargin(0.32)
p2.SetGridy()
p2.Draw()
p2.cd()

frame=sobTempl.Clone()
frame.GetYaxis().SetRangeUser(0.5,3.2)
frame.GetYaxis().CenterTitle()
frame.GetYaxis().SetTitle('Data/#scale[0.7]{#bf{#sum}}Bkg')
frame.GetYaxis().SetTitleOffset(0.4)
frame.GetYaxis().SetTitleSize(0.13)
frame.GetYaxis().SetLabelSize(0.13)
frame.GetXaxis().SetTitleSize(0.13)
frame.GetXaxis().SetLabelSize(0.13)
frame.GetYaxis().SetNdivisions(4)
frame.Draw()

breluncBand=buncBand.Clone('brelunc')
for xbin in range(1,buncBand.GetNbinsX()+1):
    breluncBand.SetBinContent(xbin,1)
    breluncBand.SetBinError(xbin,buncBand.GetBinError(xbin)/buncBand.GetBinContent(xbin))
breluncBand.Draw('e2same')

spbob=sobH['sig'].Clone('spbob')
spbob.Add(sobH['bkg'])
spbob.Divide(sobH['bkg'])
spbob.SetFillStyle(0)
spbob.SetLineColor(ROOT.kRed+2)
spbob.SetLineWidth(3)
spbob.Draw('histsame')

dataob=sobH['data'].Clone('dob')
dataob.Divide(sobH['bkg'])
dataob.Draw('X0e1same')

p2.RedrawAxis()
c.cd()
c.Modified()
c.Update()
for ext in ['png','pdf']: c.SaveAs('sob_%s.%s'%(pfix,ext))

