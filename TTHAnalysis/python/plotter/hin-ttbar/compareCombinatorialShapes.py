import ROOT
import sys
import os

def getShapeVars(target,kfact,name):

    khisto=kfact[0].Clone('khisto')
    khisto.Divide(kfact[1])
    
    vars=[]
    for v in ['Up','Down']:
        vars.append( target.Clone(target.GetName()+name+'Up' ) )
        for xbin in range(vars[-1].GetNbinsX()):
            val=vars[-1].GetBinContent(xbin+1)
            delta=val*(khisto.GetBinContent(xbin+1)-1.)
            if v=='Down': delta *=-1
            vars[-1].SetBinContent(xbin+1,max(val+delta,1e-6))
            vars[-1].SetBinError(xbin+1,0)

    khisto.Delete()
    
    gr=ROOT.TGraphAsymmErrors()
    gr.SetName(target.GetName()+name)
    gr.SetTitle(name)
    for xbin in range(target.GetNbinsX()):       
        xcen=target.GetXaxis().GetBinCenter(xbin+1)
        xwid=target.GetXaxis().GetBinWidth(xbin+1)
        yval=target.GetBinContent(xbin+1)
        delta=[vars[i].GetBinContent(xbin+1)-yval for i in range(2)]
        gr.SetPoint(xbin,xcen,yval)
        gr.SetPointError(xbin,0.5*xwid,0.5*xwid,abs(delta[1]),abs(delta[0]))
                    
    return vars,gr


dir='/eos/user/p/psilva/www/HIN-19-001/combinatorialBackground/2019-05-22'
if len(sys.argv)>1 : dir=sys.argv[1]

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
c=ROOT.TCanvas('c','c',500,500)
c.SetTopMargin(0.05)
c.SetLeftMargin(0.12)
c.SetRightMargin(0.03)
c.SetBottomMargin(0.1)
for v in ['sphericity','bdt','bdtrarity']:
    for f in ['em','flavsf']:

        #read distributions from file
        dists={}
        for i in ['','-noiso']:
            fIn=ROOT.TFile.Open(dir+'-{flav}{iso}/{var}/{flav}.input.root'.format(flav=f,iso=i,var=v))
            for p in ['data_mix','data','W','data_mix_syst']:
                dists[(p,i)]=fIn.Get('x_'+p).Clone(p+i)
                dists[(p,i)].SetDirectory(0)
                dists[(p,i)].SetFillStyle(0)
            fIn.Close()

            nss=dists[('data',i)].Integral()
            dists[('data_mix',i)].Scale(nss/dists[('data_mix',i)].Integral())
            dists[('data_mix_syst',i)].Scale(nss/dists[('data_mix_syst',i)].Integral())


        #compute the variations
        shape_vars_mix, unc_mix = getShapeVars( target=dists[('data_mix','')],
                                                kfact=(dists[('data_mix_syst','')],dists[('data_mix','')]),
                                                name='mixshape' )
        shape_vars_ss, unc_ss   = getShapeVars( target=dists[('data_mix','')],
                                                kfact=(dists[('data','-noiso')],dists[('data_mix','-noiso')]),
                                                name='ssshape' )

        #show the variations
        c.Clear()
        frame=dists[('data_mix','')].Clone('frame')
        frame.Reset('ICE')
        frame.Draw()

        unc_ss.SetFillStyle(1001)
        unc_ss.SetFillColor(ROOT.kGray)
        unc_mix.SetTitle('ss shape')
        #unc_ss.Draw('e2')

        unc_mix.SetFillStyle(3353)
        unc_mix.SetFillColor(1)
        unc_mix.SetTitle('mix shape')
        unc_mix.Draw('e2')

        unc_stat=ROOT.TGraphErrors(dists[('data_mix','')])
        unc_stat.SetFillStyle(3325)
        unc_stat.SetFillColor(ROOT.kCyan)
        unc_stat.Draw('e2')

        dists[('data_mix','')].SetTitle('comb. background')
        dists[('data_mix','')].SetLineWidth(2)
        dists[('data_mix','')].Draw('histsame')

        dists[('W','')].SetLineColor(ROOT.kRed)
        dists[('W','')].SetLineWidth(2)
        dists[('W','')].SetLineStyle(9)
        dists[('W','')].SetTitle('W (MC)')
        dists[('W','')].Draw('histsame')

        dists[('data','')].SetMarkerStyle(20)
        dists[('data','')].Draw('e1same')

        xtit='Sphericity'
        if v=='bdt' : xtit='BDT'
        if v=='bdtrarity' : xtit='BDT (rarity)'
        frame.GetXaxis().SetTitle(xtit)

        leg=ROOT.TLegend(0.15,0.92,0.7,0.7)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.035)
        leg.AddEntry(dists[('data','')],       'Data (SS)',   'ep')
        leg.AddEntry(dists[('W','')],          'W+jets (MC)', 'l')
        leg.AddEntry(dists[('data_mix','')],   'Comb. model', 'l')  
        leg.AddEntry(unc_stat,                 'stat unc.', 'f')      
        leg.AddEntry(unc_ss,                   'SS shape unc.', 'f')
        leg.AddEntry(unc_mix,                  'mix shape unc.', 'f')
        leg.SetNColumns(2)
        leg.Draw()

        txt=ROOT.TLatex()
        txt.SetNDC(True)
        txt.SetTextFont(42)
        txt.SetTextSize(0.045)
        txt.SetTextAlign(12)
        txt.DrawLatex(0.12,0.97,'#bf{CMS} #it{preliminary}')
        txt.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
        txt.DrawLatex(0.95,0.97,'#scale[0.8]{1.6 #mub^{-1} (#sqrt{s_{NN}}=5.02 TeV)}')

        c.Modified()
        c.Update()

        c.SetLogy(False)
        frame.GetYaxis().SetRangeUser(0,1.5*max([dists[(x,'')].GetMaximum() for x in ['data','W','data_mix']]))
        for ext in ['png','pdf']:
            c.SaveAs('{var}_{flav}_combmodel.{ext}'.format(var=v,flav=f,ext=ext))

        c.SetLogy(True)
        frame.GetYaxis().SetRangeUser(1e-1,frame.GetMaximum()*5)
        for ext in ['png','pdf']:
            c.SaveAs('{var}_{flav}_combmodel_log.{ext}'.format(var=v,flav=f,ext=ext))
