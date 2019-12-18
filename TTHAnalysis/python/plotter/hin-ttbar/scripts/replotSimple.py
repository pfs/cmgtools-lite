import ROOT
import os
import sys

def getNorms(dcDir):

    """reads prefit and post-fit normalizations and returns in form of scale factors"""
    
    inF=ROOT.TFile.Open(os.path.join(dcDir,'combinedCard/fitDiagnosticsobs.root'))

    norms={}
    for tag,rasName in [('prefit','norm_prefit'),
                        ('postfit','norm_fit_s')]:
        ras=inF.Get(rasName)
        iter = ras.createIterator()
        var = iter.Next()
        while var :
            proc=var.GetName().split('/')[-1]
            if not proc in norms:
                norms[proc]={'prefit':[0,0],
                             'postfit':[0,0]}
            norms[proc][tag][0] += var.getVal()
            norms[proc][tag][1] += (var.getError())**2
            var = iter.Next()

    inF.Close()

    for proc in norms:
        for tag in norms[proc]:
            norms[proc][tag][1]=ROOT.TMath.Sqrt(norms[proc][tag][1])

        #convert to scale factors instead
        norms[proc]['postfit'][1] /= norms[proc]['prefit'][0]
        norms[proc]['postfit'][0] /= norms[proc]['prefit'][0]
        norms[proc]['prefit'][1]  /= norms[proc]['prefit'][0]
        norms[proc]['prefit'][0]   = 1.

    return norms

def getSpecifics(plotter,pname):
    """a quite fancy customization function :p """

    lmargin=0.12
    if '-leponZ' in plotter or '-sf' in plotter:
        lmargin=0.14

    maxSF=1.2
    if pname=='mllz' or pname=='sphericity' or (pname=='acoplan' and '-anyflavor' in plotter):
        maxSF=1.5

    yoffset=1.0
    if '-leponZ' in plotter or '-sf' in plotter:
        yoffset=1.2

    p3xndc=0.5
    if '-anyflavor' in plotter and 'sphericity' in pname:
        p3xndc=0.15
    if (pname=='acoplan' and '-anyflavor' in plotter):
        p3xndc=0.55


    return lmargin,maxSF,yoffset,p3xndc

def redrawWithNorms(norms,tag,plotter,tagTitle):

    baseDir=os.path.dirname(plotter)

    inF=ROOT.TFile.Open(plotter)
    for key in inF.GetListOfKeys():
        kname=key.GetName()
        if not 'canvas' in kname: continue
        pname=kname.replace('_canvas','')

        lmargin,maxSF,yoffset,p3xndc=getSpecifics(plotter,pname)

        c=key.ReadObj()

        c.Draw()

        #scale individual processes and add a normalization error band
        p1=c.GetListOfPrimitives().At(0)
        p1.SetTopMargin(0.07)
        p1.SetRightMargin(0.03)
        p1.SetLeftMargin(lmargin)
        p1.cd()

        #update labels with categories
        p1.GetListOfPrimitives().At(8).SetTitle('')
        p1.GetListOfPrimitives().At(6).SetTitle('')
        p1.GetListOfPrimitives().At(7).SetTitle('')      
        
        tex=ROOT.TLatex()
        tex.SetTextFont(42)
        tex.SetTextSize(0.055)
        tex.SetNDC()
        tex.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
        tex.DrawLatex(0.95,0.96,'#scale[0.8]{1.7 nb^{-1} (#sqrt{s_{NN}}=5.02 TeV)}')
        tex.SetTextAlign(ROOT.kHAlignLeft+ROOT.kVAlignCenter)
        tex.DrawLatex(lmargin,0.96,'#bf{CMS} #it{Preliminary}')        
        tex.SetTextSize(0.04)
        extraTitle='e#mu'
        if '-leponZee'  in plotter : extraTitle='Z#rightarrow ee'
        if '-leponZmm'  in plotter : extraTitle='Z#rightarrow #mu#mu'
        if '-sf'        in plotter : extraTitle='ee,#mu#mu'
        if '-anyflavor' in plotter : extraTitle='ee,#mu#mu,e#mu (BDT>0.5)'
        tex.DrawLatex(lmargin+0.04,0.87,'#it{%s}'%extraTitle)
        tex.DrawLatex(lmargin+0.04,0.82,'#it{%s}'%(tag))
        if tag=='postfit':
           tex.DrawLatex(lmargin+0.04,0.77,'#it{#scale[0.8]{(%s)}}'%tagTitle)

        data=p1.GetListOfPrimitives().At(4)

        #compute total before and after and sum up errors in each bin
        totalPre,totalPost=0.,0.
        totalH=p1.GetListOfPrimitives().At(1)
        maxY=totalH.GetMaximum()
        newtotalH=totalH.Clone('newtotal')
        newtotalH.Reset('ICE')
        totalUncH=totalH.Clone('totalunch')
        totalUncH.Reset('ICE')
        bkgH=totalH.Clone('newbkg')
        bkgH.Reset('ICE')
        bkgUncH=totalH.Clone('bkgunc')
        bkgUncH.Reset('ICE')
        stack=p1.GetListOfPrimitives().At(2)
        ttbarH=None
        for i in range(stack.GetNhists()):
            h=stack.GetHists().At(i)
            proc=h.GetName()
            proc=proc.replace(pname+'_','')
            
            totalPre+=h.Integral()
            h.Scale( norms[proc][tag][0] )
            totalPost+=h.Integral()

            newtotalH.Add(h)
            if proc!='ttbar':
                bkgH.Add(h)
            else:
                ttbarH=h.Clone('ttbar')

            for xbin in range(1,totalUncH.GetNbinsX()+1):
                curSum=totalUncH.GetBinContent(xbin)
                totalUncH.SetBinContent(xbin,
                                        curSum+(h.GetBinContent(xbin)*norms[proc][tag][1])**2)
                if proc!='ttbar':
                    curSum=bkgUncH.GetBinContent(xbin)
                    bkgUncH.SetBinContent(xbin,
                                          curSum+(h.GetBinContent(xbin)*norms[proc][tag][1])**2)

        #update contents and errors
        for xbin in range(1,totalUncH.GetNbinsX()+1):
            totalH.SetBinContent(xbin,
                                 newtotalH.GetBinContent(xbin))
            totalH.SetBinError(xbin,
                               ROOT.TMath.Sqrt(totalUncH.GetBinContent(xbin)))
            bkgH.SetBinError(xbin,
                             ROOT.TMath.Sqrt(bkgUncH.GetBinContent(xbin)))
        
        bkgUncH.Delete()
        totalUncH.Delete()
        newtotalH.Delete()
        totalH.SetMarkerStyle(1)
        totalH.SetFillStyle(3244)
        totalH.SetFillColor(1)
        totalH.SetLineColor(1)
        totalH.GetXaxis().SetLabelSize(0)
        totalH.GetYaxis().SetRangeUser(0,maxSF*maxY)
        totalH.GetYaxis().SetTitleOffset(yoffset)
        totalH.GetYaxis().SetTitleSize(0.06)
        totalH.Draw('e2same')

        #update legend
        leg=p1.GetListOfPrimitives().At(5)
        leg.AddEntry(totalH,'Norm. unc.','f')
        leg.SetY2NDC(0.9)
        leg.SetY1NDC(0.65)

        
        p1.RedrawAxis()

        #redo the ratio
        c.cd()
        p2=c.GetListOfPrimitives().At(1)
        p2.SetRightMargin(0.03)
        p2.SetLeftMargin(0.12)
        p2.cd()

        relUnc=totalH.Clone('relUnc')
        for xbin in range(1,relUnc.GetNbinsX()+1):
            val=relUnc.GetBinContent(xbin)
            relUnc.SetBinError(xbin,relUnc.GetBinError(xbin)/val if val!=0 else 0.)
            relUnc.SetBinContent(xbin,1)
        
        relUnc.Draw('e2')
        relUnc.GetYaxis().SetTitle('Data/Pred.')
        relUnc.GetXaxis().SetTitleSize(0.12)
        relUnc.GetXaxis().SetLabelSize(0.12)
        relUnc.GetYaxis().SetTitleOffset(0.4)
        relUnc.GetYaxis().SetTitleSize(0.15)
        relUnc.GetYaxis().SetLabelSize(0.12)
        #relUnc.GetYaxis().SetRangeUser(0.44,1.56)
        relUnc.GetYaxis().SetRangeUser(0.,2.)
        relUnc.GetYaxis().SetNdivisions(4)

        ratio=data.Clone('ratio')
        bsdata=data.Clone('bsdata')        
        x,y=ROOT.Double(0),ROOT.Double(0)
        minSub,maxSub=999,-999
        for i in range(ratio.GetN()):
            ratio.GetPoint(i,x,y)
            den=totalH.GetBinContent(i+1)
            
            eylo=abs(ratio.GetErrorYlow(i)/den if den!=0 else 0.)
            eyhi=abs(ratio.GetErrorYhigh(i)/den if den !=0 else 0.)
            r=float(y)/den if den!=0. else 0.
            ratio.SetPoint(i,x,r)
            ratio.SetPointError(i,0,0,eylo,eyhi)

            subVal=float(y)-bkgH.GetBinContent(i+1)
            minSub=min(subVal,minSub)
            maxSub=max(subVal,maxSub)
            bsdata.SetPoint(i,x,subVal)

        ratio.Draw('p')

        p2.RedrawAxis()

        c.cd()
        c.Modified()
        c.Update()
        for ext in ['png','pdf']:
            c.SaveAs(os.path.join(baseDir,'%s_%s.%s'%(pname,tag,ext)))

        #add background subtracted inset
        p1.cd()
        p3=ROOT.TPad('p3','p3',p3xndc,0.66,p3xndc+0.4,0.35)
        p3.SetFillStyle(0)
        p3.SetBorderSize(0)
        p3.Draw()
        p3.cd()
        ttbarH.Draw('hist')
        ttbarH.GetYaxis().SetRangeUser(minSub*1.5 if minSub<0 else minSub*0.5,maxSub*1.5)
        ttbarH.GetYaxis().SetNdivisions(4)
        ttbarH.GetXaxis().SetNdivisions(5)
        ttbarH.GetYaxis().SetLabelSize(0.1)
        ttbarH.GetXaxis().SetLabelSize(0.1)
        ttbarH.GetYaxis().SetTitleSize(0.)
        ttbarH.GetXaxis().SetTitleSize(0.)


        bsdata.Draw('p')
        p3.RedrawAxis()

        c.cd()
        c.Modified()
        c.Update()
        for ext in ['png','pdf']:
            c.SaveAs(os.path.join(baseDir,'%s_%s_bkgsubinset.%s'%(pname,tag,ext)))




ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)

plotter='/eos/user/p/psilva/www/HIN-19-001/simple_plots/2019-10-29-em/mll_AND_mllz_AND_llpt_AND_sphericity_AND_acoplan_AND_minmlb_AND_minmlb2_AND_centrality.root'
if len(sys.argv)>1:
    plotter=sys.argv[1]

dcDir='/afs/cern.ch/work/m/mdunser/public/cmssw/heavyIons/CMSSW_9_4_6_patch1/src/CMGTools/TTHAnalysis/python/plotter/hin-ttbar/datacards_2019-10-21_forApproval_jetAnalysis/bdtcomb'
if len(sys.argv)>2:
    dcDir=sys.argv[2]

postfitTitle='leptonic+b-tagged' if '_jetAnalysis' in dcDir else 'leptonic'

norms=getNorms(dcDir)
for key in ['prefit','postfit']:
    redrawWithNorms(norms,key,plotter,postfitTitle)
