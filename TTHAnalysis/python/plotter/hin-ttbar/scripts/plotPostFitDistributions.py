import sys, os, optparse, copy, math
import ROOT
from array import array

## usage:
## python plotPostFitDistributions.py /eos/cms/store/cmst3/group/hintt/fits/datacards_2019-09-19_elPt25muPt20_ttbarBRFix/sphericity/allFlavors.card/fitDiagnosticsobs.root --outdir <directoryForPlots> (--inclusive)

doPrel=False

mcStack=[('zg'        , 'Z/#gamma*'         , ROOT.kAzure+6) ,
         ('data_comb' , 'Nonprompt'            , ROOT.kOrange ) ,
         ('tW'        , 'tW'                   , ROOT.kTeal+9 ) ,
         ('VV'        , 'VV'                   , ROOT.kAzure-5) ,
         ('ttbar'     , 't#bar{t}' , 633)           , ]

#yranges={'ee0b': [0., 450.], 'mm0b': [0., 1300.], 'em0b': [0., 55.],
#         'ee1b': [0.,  45.], 'mm1b': [0.,  170.], 'em1b': [0., 21.],
#         'ee2b': [0.,  12.], 'mm2b': [0.,   32.], 'em2b': [0., 12.]}

yrangesInclusive           = {'ee': [0., 300.], 'mm': [0., 1000.], 'em': [0., 40.]}
yranges                    = {'ee': [1., 5000.], 'mm': [1., 10000.], 'em': [0., 50.]}


yrangesInclusiveSphericity = {'ee': [0., 250.], 'mm': [0.,  500.], 'em': [0., 27.]}
yrangesNb = {'ee': [0.5, 10000.], 'mm': [0.5,  50000.], 'em': [0.5, 1000.]}
yrangesSmall = {'ee': [0., 40.], 'mm': [0.,  55.], 'em': [0., 20.]}

valsAndErrors = {}


def saveTable(tex):
    f = open(options.outdir+'/yieldTable{a}.txt'.format(a='Inclusive' if options.inclusive else ''), 'w')
    f.write(tex)
    f.close()

def printTable(vae):

    for k,v in vae.items():
        if not 'obs' in k:
            vae[k] = '{a:.1f}'.format(a=v)
        else:
            vae[k] = '{a:.0f}'.format(a=v)

    vae['when' ] = 'after the fit' if not options.prefit else 'before the fit'
    vae['when2'] = 'Postfit' if not options.prefit else ''

    tabletex = '''
\\begin{{table}}[!htb]
\\centering
\\small
\\topcaption{{
    The number of expected background and signal events and the observed event yields in the different
    channels of $ee$, $\\mu\\mu$, and $e\\mu$, \\textbf{{{when}}}.
    \\label{{tab:yieldsJet{when2}}}}}
    \\resizebox{{\\textwidth}}{{!}}{{
    \\begin{{tabular}}{{lccc|ccc|ccc}}
Process        & $ee$ 0b    & $ee$ 1b   & $ee$ 2b   & $\\mu\\mu$ 0b   & $\\mu\\mu$ 1b  & $\\mu\\mu$ 2b  & $e\\mu$ 0b    & $e\\mu$  1b    & $e\\mu$  2b   \\\\ \\hline \\hline

Z/$\\gamma^{{*}}$ & {zg_ee0b} $\\pm$ {zg_ee0b_e} & {zg_ee1b} $\\pm$ {zg_ee1b_e} & {zg_ee2b} $\\pm$ {zg_ee2b_e} & {zg_mm0b} $\\pm$ {zg_mm0b_e} & {zg_mm1b} $\\pm$ {zg_mm1b_e} & {zg_mm2b} $\\pm$ {zg_mm2b_e} & {zg_em0b} $\\pm$ {zg_em0b_e} & {zg_em1b} $\\pm$ {zg_em1b_e} & {zg_em2b} $\\pm$ {zg_em2b_e}  \\\\
Nonprompt & {data_comb_ee0b} $\\pm$ {data_comb_ee0b_e} & {data_comb_ee1b} $\\pm$ {data_comb_ee1b_e} & {data_comb_ee2b} $\\pm$ {data_comb_ee2b_e} & {data_comb_mm0b} $\\pm$ {data_comb_mm0b_e} & {data_comb_mm1b} $\\pm$ {data_comb_mm1b_e} & {data_comb_mm2b} $\\pm$ {data_comb_mm2b_e} & {data_comb_em0b} $\\pm$ {data_comb_em0b_e} & {data_comb_em1b} $\\pm$ {data_comb_em1b_e} & {data_comb_em2b} $\\pm$ {data_comb_em2b_e}  \\\\
tW & {tW_ee0b} $\\pm$ {tW_ee0b_e} & {tW_ee1b} $\\pm$ {tW_ee1b_e} & {tW_ee2b} $\\pm$ {tW_ee2b_e} & {tW_mm0b} $\\pm$ {tW_mm0b_e} & {tW_mm1b} $\\pm$ {tW_mm1b_e} & {tW_mm2b} $\\pm$ {tW_mm2b_e} & {tW_em0b} $\\pm$ {tW_em0b_e} & {tW_em1b} $\\pm$ {tW_em1b_e} & {tW_em2b} $\\pm$ {tW_em2b_e}  \\\\
VV & {VV_ee0b} $\\pm$ {VV_ee0b_e} & {VV_ee1b} $\\pm$ {VV_ee1b_e} & {VV_ee2b} $\\pm$ {VV_ee2b_e} & {VV_mm0b} $\\pm$ {VV_mm0b_e} & {VV_mm1b} $\\pm$ {VV_mm1b_e} & {VV_mm2b} $\\pm$ {VV_mm2b_e} & {VV_em0b} $\\pm$ {VV_em0b_e} & {VV_em1b} $\\pm$ {VV_em1b_e} & {VV_em2b} $\\pm$ {VV_em2b_e}  \\\\ \\hline
Total background & {total_background_ee0b} $\\pm$ {total_background_ee0b_e} & {total_background_ee1b} $\\pm$ {total_background_ee1b_e} & {total_background_ee2b} $\\pm$ {total_background_ee2b_e} & {total_background_mm0b} $\\pm$ {total_background_mm0b_e} & {total_background_mm1b} $\\pm$ {total_background_mm1b_e} & {total_background_mm2b} $\\pm$ {total_background_mm2b_e} & {total_background_em0b} $\\pm$ {total_background_em0b_e} & {total_background_em1b} $\\pm$ {total_background_em1b_e} & {total_background_em2b} $\\pm$ {total_background_em2b_e}  \\\\
\\ttbar signal & {ttbar_ee0b} $\\pm$ {ttbar_ee0b_e} & {ttbar_ee1b} $\\pm$ {ttbar_ee1b_e} & {ttbar_ee2b} $\\pm$ {ttbar_ee2b_e} & {ttbar_mm0b} $\\pm$ {ttbar_mm0b_e} & {ttbar_mm1b} $\\pm$ {ttbar_mm1b_e} & {ttbar_mm2b} $\\pm$ {ttbar_mm2b_e} & {ttbar_em0b} $\\pm$ {ttbar_em0b_e} & {ttbar_em1b} $\\pm$ {ttbar_em1b_e} & {ttbar_em2b} $\\pm$ {ttbar_em2b_e}  \\\\ \\hline
Observed (data) & \\textbf{{{obs_ee0b}}} & \\textbf{{{obs_ee1b}}} & \\textbf{{{obs_ee2b}}} & \\textbf{{{obs_mm0b}}} & \\textbf{{{obs_mm1b}}} & \\textbf{{{obs_mm2b}}} & \\textbf{{{obs_em0b}}} & \\textbf{{{obs_em1b}}} & \\textbf{{{obs_em2b}}}   \\\\

\\end{{tabular}} }}

\\end{{table}}'''.format(**vae)
    print tabletex

    saveTable(tabletex)

def printTableInclusive(vae):

    for k,v in vae.items():
        if not 'obs' in k:
            vae[k] = '{a:.1f}'.format(a=v)
        else:
            vae[k] = '{a:.0f}'.format(a=v)

    vae['when' ] = 'after the fit' if not options.prefit else 'before the fit'
    vae['when2'] = 'Postfit' if not options.prefit else ''

    tabletex = '''\\begin{table}[!htb]
    \\centering
    \\topcaption{ The number of expected background and signal events and the observed event yields in the different
    channels of $ee$, $\\mu\\mu$, and $e\\mu$, \\textbf{{ {when} }}.
    \\label{tab:yields{when2}}}
        \\begin{tabular}{lccc}
            Process        & $ee$  & $\\mu\\mu$ & $e\\mu$ \\\\ \\hline \\hline\n'''
        
    tabletex += '''
    Nonprompt         & {data_comb_ee} $\\pm$  {data_comb_ee_e} & {data_comb_mm} $\\pm$ {data_comb_mm_e} & {data_comb_em} $\\pm$ {data_comb_em_e}       \\\\
    Z/$\\gamma^{{*}}$ & {zg_ee} $\\pm$  {zg_ee_e} & {zg_mm} $\\pm$  {zg_mm_e} & {zg_em} $\\pm$  {zg_em_e} \\\\
    $\\cPqt\\PW$      & {tW_ee} $\\pm$  {tW_ee_e} & {tW_mm} $\\pm$  {tW_mm_e} & {tW_em} $\\pm$  {tW_em_e} \\\\
    VV                & {VV_ee} $\\pm$  {VV_ee_e} & {VV_mm} $\\pm$  {VV_mm_e} & {VV_em} $\\pm$  {VV_em_e} \\\\ \\hline \\hline
    Total background  & {total_background_ee} $\\pm$ {total_background_ee_e} & {total_background_mm} $\\pm$ {total_background_mm_e} & {total_background_em} $\\pm$ {total_background_em_e}\\\\
    \\ttbar signal    & {ttbar_ee} $\\pm$  {ttbar_ee_e} & {ttbar_mm} $\\pm$  {ttbar_mm_e} & {ttbar_em} $\\pm$  {ttbar_em_e} \\\\ \\hline \\hline
    Observed (data)   & {obs_ee} & {obs_mm} & {obs_em} \\\\ \\hline \\hline'''.format(**vae)
    
    tabletex += '''
    
    \\end{tabular}
\\end{table}'''

    print tabletex

    saveTable(tabletex)


def convertGraph(inGraph, xtitle):
    if 'BDT' not in xtitle:
        return inGraph

    outGraph = inGraph.Clone(inGraph.GetName()+'_converted')
    npts=inGraph.GetN()
    for ip in range(npts):
        tmp_x, tmp_y = ROOT.Double(), ROOT.Double()
        inGraph.GetPoint(ip, tmp_x, tmp_y)
        #print 'for point {i} found x {x} and y {y}'.format(i=ip,x=tmp_x,y=tmp_y)
        outGraph.SetPoint(ip, (ip+0.5)/float(npts), tmp_y)

    return outGraph

def convertHisto(inHisto, plotName, xtitle):

    if 'BDT' not in xtitle:
        return inHisto

    plotName = plotName.split('/')[-1]

    nbins=inHisto.GetNbinsX()

    #nbins = 3 
    #if '2b' in plotName:
    #    nbins = 1

    if options.inclusive:
        nbins = 10
    
    outHisto = ROOT.TH1F(inHisto.GetName()+'_converted', '', nbins, 0, 1) #0.5, nbins+0.5)

    for ip in range(inHisto.GetXaxis().GetNbins()):
        outHisto.SetBinContent(ip+1, inHisto.GetBinContent(ip+1))
        outHisto.SetBinError  (ip+1, inHisto.GetBinError(ip+1)  )
        if not options.inclusive:
            binLabel='[0,1/3['
            if ip%3==1:
                binLabel='[1/3,2/3['
            if ip%3==2:
                binLabel='[2/3,1]'
            if ip==inHisto.GetXaxis().GetNbins()-1:
                binLabel='[0,1]'
            outHisto.GetXaxis().SetBinLabel(ip+1,binLabel)

        #print 'set bincontent', outHisto.
    
    
    return outHisto

def convertValsAndErrorsToHistos(ch, vae):

    returnDict = {}

    tmp_hist_total = ROOT.TH1F('total', 'total', 3, -0.5, 2.5)

    for proc in ['ttbar', 'total_background', 'VV', 'obs', 'tW', 'zg', 'data_comb']:
        if not proc == 'obs':
            tmp_hist = ROOT.TH1F(proc, proc, 3, -0.5, 2.5)
            for ib in range(3):
                tmp_hist.SetBinContent(ib+1, float(vae[proc+'_'+ch+str(ib)+'b']))
                tmp_hist.SetBinError  (ib+1, float(vae[proc+'_'+ch+str(ib)+'b_e']))
            if not proc == 'total_background':
                tmp_hist_total.Add(tmp_hist)
        else:
            tmp_hist = ROOT.TGraphAsymmErrors(3)
            tmp_hist.SetName('data')
            tmp_hist_data = ROOT.TH1F(proc, proc, 3, -0.5, 2.5)
            tmp_hist_data.Sumw2(0)
            tmp_hist_data.SetBinErrorOption(ROOT.TH1.kPoisson)
            for ib in range(3):
                val = int(vae[proc+'_'+ch+str(ib)+'b'])
                tmp_hist_data.SetBinContent(ib+1, val)
            for ib in range(3):
                val = int(vae[proc+'_'+ch+str(ib)+'b'])
                tmp_hist.SetPoint(ib, ib, val)
                tmp_hist.SetPointError(ib, 0., 0., tmp_hist_data.GetBinErrorLow(ib+1), tmp_hist_data.GetBinErrorUp(ib+1))

        returnDict[proc if not proc == 'obs' else 'data'] = copy.deepcopy(tmp_hist)

    returnDict['total'] = copy.deepcopy(tmp_hist_total)

    return returnDict
    

def doPostFitPlot(url,inclusive):

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)

    plotNameAdd = ''

    inF = ROOT.TFile(url)
    xtitle='BDT' if 'bdtcomb' in url else 'Sphericity'
    if 'ptll' in url: 
        xtitle = 'p_{T}^{ll}'
        plotNameAdd = 'ptll'
    if 'sphericity' in url:
        plotNameAdd = 'sphericity'

    print 'this is xtitle', xtitle

    megaPlotsPrefit={}
    megaPlotsPostfit={}
    for ich,ch in enumerate(inF.Get('shapes_fit_s' if not options.prefit else 'shapes_prefit').GetListOfKeys()):
        chName=ch.GetName()


        flav=chName.split('_')[0]

        chTitle=flav.replace('mm','#mu#mu')
        chTitle=chTitle.replace('em','e#mu')
        chTitle=chTitle.replace('_',' ')
        #if url.find('_jetAnalysis')>0 : chTitle += '+b-tags'

        chDir=ch.ReadObj()
        plotsPrefit={}
        plotsPostfit={}
        for proc in chDir.GetListOfKeys():

            pname=proc.GetName()
            if pname=='total_covar' : continue

            plotsPostfit[pname]=proc.ReadObj()
            try:
                plotsPostfit[pname].SetDirectory(0)
            except:
                pass
            plotsPrefit[pname]=inF.Get('shapes_prefit/%s/%s'%(chName,pname))
            try:
                plotsPrefit[pname].SetDirectory(0)
            except:
                pass
        
                
        prePost = 'postfit_' if not options.prefit else 'prefit_'
        compareFitResult(plotsPostfit=plotsPostfit,plotsPrefit=plotsPrefit,
                         plotName=options.outdir+prePost+plotNameAdd+chName,
                         xtitle=xtitle,extraTxt=[chTitle])

        xoffset=0
        if '1b' in chName:
            xoffset=3
        if '2b' in chName:
            xoffset=6

        if not flav in megaPlotsPrefit:
            megaPlotsPrefit[flav]={}
            megaPlotsPostfit[flav]={}

        for proc in plotsPrefit:

            cname=plotsPrefit[proc].ClassName()
            if cname=='TH1F':
                
                if not proc in megaPlotsPrefit[flav]:
                    megaPlotsPrefit[flav][proc]=ROOT.TH1F(flav+proc,';%s;Events'%xtitle,7,0,7) 
                    megaPlotsPrefit[flav][proc].SetDirectory(0)
                    megaPlotsPostfit[flav][proc]=megaPlotsPrefit[flav][proc].Clone(flav+proc+'postfit')
                    megaPlotsPostfit[flav][proc].SetDirectory(0)

                for xbin in range(1,plotsPrefit[proc].GetNbinsX()+1):
                    if xbin+xoffset>7 : break
                    megaPlotsPrefit[flav][proc].SetBinContent(xbin+xoffset,plotsPrefit[proc].GetBinContent(xbin))
                    megaPlotsPrefit[flav][proc].SetBinError(xbin+xoffset,plotsPrefit[proc].GetBinError(xbin))
                    megaPlotsPostfit[flav][proc].SetBinContent(xbin+xoffset,plotsPostfit[proc].GetBinContent(xbin))
                    megaPlotsPostfit[flav][proc].SetBinError(xbin+xoffset,plotsPostfit[proc].GetBinError(xbin))

            else:
                
                if not proc in megaPlotsPrefit[flav]:
                    megaPlotsPostfit[flav][proc]=ROOT.TGraphAsymmErrors()
                    megaPlotsPrefit[flav][proc]=ROOT.TGraphAsymmErrors()

                x,y=ROOT.Double(0),ROOT.Double(0)
                for ix in range(plotsPrefit[proc].GetN()):
                    if ix+xoffset>6: break
                    plotsPrefit[proc].GetPoint(ix,x,y)
                    megaPlotsPrefit[flav][proc].SetPoint(ix+xoffset,x,y)
                    megaPlotsPrefit[flav][proc].SetPointError(ix+xoffset,
                                                              0,
                                                              0,
                                                              plotsPrefit[proc].GetErrorYlow(ix),
                                                              plotsPrefit[proc].GetErrorYhigh(ix))
                    plotsPostfit[proc].GetPoint(ix,x,y)
                    megaPlotsPostfit[flav][proc].SetPoint(ix+xoffset,x,y)
                    megaPlotsPostfit[flav][proc].SetPointError(ix+xoffset,
                                                               0,
                                                               0,
                                                               plotsPostfit[proc].GetErrorYlow(ix),
                                                               plotsPostfit[proc].GetErrorYhigh(ix))


    if not inclusive:
        for flav in megaPlotsPostfit: 
            chTitle=flav.replace('mm','#mu#mu')
            chTitle=chTitle.replace('em','e#mu')
            chTitle=chTitle.replace('_',' ')
            compareFitResult(plotsPostfit=megaPlotsPostfit[flav],
                             plotsPrefit=megaPlotsPrefit[flav],
                             plotName=options.outdir+'mega'+prePost+plotNameAdd+'_'+flav,
                             xtitle=xtitle,extraTxt=[chTitle])

        # print 'this is plotsPostfit', plotsPostfit
     
def compareFitResult(plotsPrefit,plotsPostfit,plotName,xtitle,extraTxt=[]):

    if not options.inclusive:
        channel = ''.join(plotName.split('_')[-2:])
        channel = channel.replace('_','')
    else:
        channel = plotName.split('_')[-1]

    if 'nbjets' in channel:
        channel = channel[-2:]
    if 'ptll' in channel:
        channel = channel[-2:]
    if 'sphericity' in channel:
        channel = channel[-2:]


    print '=================='
    print '== at plotName', plotName, '==='
    print '== in channel ', channel, '==='
    print '=================='

    marginL = 0.14
    marginR = 0.03
    
    c = ROOT.TCanvas("c","c",800,800)
    c.SetTopMargin(0)
    c.SetLeftMargin(0)
    c.SetRightMargin(0)
    c.SetBottomMargin(0)

    #data/MC
    p1 = ROOT.TPad("p1", "p1", 0., 0.3, 1., 1.)
    p1.SetTopMargin(0.08)
    p1.SetRightMargin(marginR)
    p1.SetLeftMargin(marginL)
    p1.SetBottomMargin(0.03)
    p1.Draw()
    p1.cd()
    if 'b-tag' in xtitle:
        p1.SetLogy()
    if 'mega' in plotName and channel!='em':
        p1.SetLogy()


    frame=plotsPostfit['total'].Clone('frame')
    frame=convertHisto(frame,plotName,xtitle)
    frame.Reset('ICE')
    frame.GetYaxis().SetTitle('Events')
    frame.GetYaxis().SetTitleOffset(1.05)
    frame.GetYaxis().SetNdivisions(505)
    frame.GetYaxis().SetTitleSize(0.07)
    frame.GetYaxis().SetLabelSize(0.07)
    frame.GetXaxis().SetTitleSize(0)
    frame.GetXaxis().SetLabelSize(0)

    myyranges = yranges
    if options.inclusive and 'BDT' in xtitle:
        myyranges = yrangesInclusive
    if options.inclusive and 'hericity' in xtitle:
        myyranges = yrangesInclusiveSphericity
    if 'nbjets' in plotName:
        myyranges = yrangesNb
    if 'ptll' in plotName:
        myyranges = yrangesSmall 
    if 'phericity' in plotName:
        myyranges = yrangesSmall 

    print 'this is myranges', myyranges
    print 'this is channel', channel

    try:
        frame.GetYaxis().SetRangeUser(myyranges[channel][0], myyranges[channel][1]) 
    except:
        pass#frame.GetYaxis().SetRangeUser(0.,1.3*frame.GetMaximum())
    frame.Draw()

    if 'mega' in plotName:
        l=ROOT.TLine()
        l.SetNDC(True)
        l.SetLineColor(ROOT.kGray)
        l.SetLineStyle(9)
        lymax=frame.GetMaximum()
        l.DrawLine(3./7.,0.,3./7.,lymax)
        l.DrawLine(6./7.,0.,6./7.,lymax)


    leg = ROOT.TLegend(0.55,0.64,0.95,0.91)
    leg.SetNColumns(2)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.06)
    leg.AddEntry(plotsPostfit['data'],'Data','ep')

    #print plotsPostfit

    stack = ROOT.THStack("stack","stack")
    ratios=[]
    for proc,title,ci in mcStack:
        try:
            plotsPostfit[proc].SetTitle(title)
            newHisto = convertHisto(plotsPostfit[proc], plotName, xtitle)                        
            newHisto.SetLineColor(1)
            newHisto.SetFillStyle(1001)
            newHisto.SetFillColor(ci)
            newHisto.SetTitle(title)
            stack.Add(newHisto)

            error = ROOT.Double()
            integral = newHisto.IntegralAndError(1,newHisto.GetNbinsX()+1,error)
            valsAndErrors[proc+'_'+channel]      = integral
            valsAndErrors[proc+'_'+channel+'_e'] = error
            

            ratios.append( newHisto.Clone(proc+'_2prefit') )
            ratios[-1].SetDirectory(0)
            ratios[-1].Divide( newHisto )
            ratios[-1].SetFillStyle(0)
            ratios[-1].SetFillColor(0)
            ratios[-1].SetLineColor(ci)
            ratios[-1].SetLineWidth(3)
        except:
            valsAndErrors[proc+'_'+channel]      = 0.
            valsAndErrors[proc+'_'+channel+'_e'] = 0.
            print 'did not find ', proc, 'in channel', plotName.split('/')[-1]

    for ih in range(stack.GetStack().GetEntriesFast()):
        idx=stack.GetStack().GetEntriesFast()-1-ih
        leg.AddEntry(stack.GetStack().At(idx),stack.GetStack().At(idx).GetTitle(),'f')

    stack.Draw('histsame')

    error = ROOT.Double()
    integral = plotsPostfit['total_background'].IntegralAndError(1, plotsPostfit['total_background'].GetNbinsX()+1, error)
    valsAndErrors['total_background_'+channel]      = integral
    valsAndErrors['total_background_'+channel+'_e'] = error

    newTotal = convertHisto(plotsPostfit['total'], plotName, xtitle)
    plotsPostfit['newTotal']=newTotal
    totalUnc=ROOT.TGraphErrors(newTotal)
    totalUnc.SetFillStyle(3244)
    totalUnc.SetFillColor(1)
    totalUnc.SetMarkerStyle(1)
    totalUnc.Draw('e2same')
    #totalUncNew = convertGraph(totalUnc,xtitle)
    #totalUncNew.Draw('e2')
    leg.AddEntry(totalUnc,'Total unc.','f')

    plotsPostfit['data'].SetMarkerStyle(20)
    plotsPostfit['data'].SetLineWidth(2)

    integral = 0.
    for i in range(plotsPostfit['data'].GetN()):
        plotsPostfit['data'].SetPointEXhigh(i,0)
        plotsPostfit['data'].SetPointEXlow(i,0)
        tmp_x, tmp_y = ROOT.Double(), ROOT.Double()
        plotsPostfit['data'].GetPoint(i, tmp_x, tmp_y)
        integral += tmp_y
    valsAndErrors['obs_'+channel] = integral
    #plotsPostfit['data'].Draw('PZ')
    dataGraphNew = convertGraph(plotsPostfit['data'],xtitle)
    dataGraphNew.Draw('PZ')

    leg.Draw('same')

    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.07)
    tex.SetNDC()
    tex.SetTextAlign(ROOT.kHAlignLeft+ROOT.kVAlignCenter)
    if doPrel:
        tex.DrawLatex(marginL,0.96,'#bf{CMS} #it{Preliminary}')
    else:
        tex.DrawLatex(marginL,0.96,'#bf{CMS}')
    tex.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
    tex.DrawLatex(0.95,0.96,'#scale[0.8]{1.7 nb^{-1} (#sqrt{s_{NN}}=5.02 TeV)}')
    tex.SetTextAlign(ROOT.kHAlignLeft+ROOT.kVAlignCenter)
    tex.SetTextSize(0.07)
    for it in range(len(extraTxt)):
        tex.DrawLatex(marginL+0.04,0.87+it*0.05,'#it{%s}'%extraTxt[it])

    if 'mega' in plotName:

        tex.SetTextSize(0.06)
        tex.DrawLatex(0.3,0.58,'#it{0 b}')
        tex.DrawLatex(0.63,0.58,'#it{1 b}')
        tex.DrawLatex(0.89,0.58,'#it{2 b}')

    p1.RedrawAxis()

    #data/MC
    c.cd()
    p3 = ROOT.TPad("p3", "p3", 0., 0.0, 1., 0.3)
    p3.SetTopMargin(0.04)
    p3.SetRightMargin(marginR)
    p3.SetLeftMargin(marginL)
    p3.SetBottomMargin(0.4)
    #p3.SetGridy()
    p3.Draw()
    p3.cd()

    frame3=frame.Clone()
    frame3.GetYaxis().SetTitle('Data/Pred.')
    frame3.GetYaxis().SetTitleOffset(0.4)
    frame3.GetYaxis().SetTitleSize(0.16)
    frame3.GetYaxis().SetLabelSize(0.16)
    frame3.GetYaxis().SetNdivisions(5)
    frame3.GetXaxis().SetTitleSize(0.16)
    frame3.GetXaxis().SetLabelSize(0.16)
    frame3.GetXaxis().SetTitle(xtitle)
    if frame3.GetNbinsX() > 1:
        frame3.GetXaxis().SetNdivisions(505)
    else:
        frame3.GetXaxis().SetNdivisions(3)
    frame3.GetYaxis().SetRangeUser(0., 2.04)
    frame3.Draw()

    if 'mega' in plotName:
        l3=ROOT.TLine()
        l3.SetNDC(True)
        l3.SetLineColor(ROOT.kGray)
        l3.SetLineStyle(9)
        l3.DrawLine(3./7.,0.,3./7.,2.04)
        l3.DrawLine(6./7.,0.,6./7.,2.04)


    def getRelUnc(plotColl,name,ci,fill):
        totalNoUnc=plotColl[name].Clone(name+'_nounc')
        for i in range(totalNoUnc.GetNbinsX()):
            totalNoUnc.SetBinError(i+1,0)
        relUnc=plotColl[name].Clone(name+'_relUnc')
        relUnc.Divide(totalNoUnc)
        relUncGr=ROOT.TGraphErrors(relUnc)
        relUncGr.SetFillStyle(fill)
        relUncGr.SetFillColor(ci)
        relUncGr.SetMarkerStyle(1)
        relUnc.Delete()
        totalNoUnc.Delete()
        return relUncGr

    relPreUncGr=getRelUnc(plotsPrefit,'total',ROOT.kGray,1001)
    #relPreUncGr.Draw('e2')
    relFitUncGr=getRelUnc(plotsPostfit, 'newTotal', 1, 3244) # ROOT.kAzure-3, 3001)
    #relFitUncGr.Draw('e2')
    relFitUncGr2 = convertGraph(relFitUncGr, xtitle)
    relFitUncGr2.Draw('e2')
    #print 'this is relFitUncGr', relFitUncGr
    #relFitUncGr.GetYaxis().SetRangeUser(0., 2.)

    data2fitGr=plotsPostfit['data'].Clone('data2fit')
    x,y=ROOT.Double(0),ROOT.Double(0)
    for i in range(data2fitGr.GetN()):
        den=plotsPostfit['total'].GetBinContent(i+1)
        if float(den)==0 : 
            data2fitGr.SetPointEYhigh(i,0)
            data2fitGr.SetPointEYlow(i,0)
        else:
            denUnc=plotsPostfit['total'].GetBinError(i+1)
            data2fitGr.GetPoint(i,x,y)
            data2fitGr.SetPoint(i,x,y/den)
            eyhi=data2fitGr.GetErrorYhigh(i)
            eylo=data2fitGr.GetErrorYlow(i)
            data2fitGr.SetPointEYhigh(i,eyhi/den)
            data2fitGr.SetPointEYlow(i,eylo/den)

    #data2fitGr.Draw('PZ')
    data2fitGrNew = convertGraph(data2fitGr, xtitle)
    data2fitGrNew.Draw('PZ')


    leg3 = ROOT.TLegend(0.15,0.85,0.8,0.95)
    leg3.SetNColumns(3)
    leg3.SetFillStyle(0)
    leg3.SetBorderSize(0)
    leg3.SetTextSize(0.13)
    leg3.AddEntry(data2fitGr,'Data','ep')
    #leg3.AddEntry(relPreUncGr,'Prefit unc.','f')
    leg3.AddEntry(relFitUncGr,'Postfit unc.' if not options.prefit else 'Prefit unc.','f')
    #leg3.Draw()

    #p3.RedrawAxis()
    
    c.cd()
    c.Modified()
    c.Update()
    prepost = '_postfit' if not options.prefit else '_prefit'
    finalxtitle = xtitle.replace(' ','').replace('{','').replace('}','').replace('-','').replace('^','')
    for ext in ['png','pdf']:
        c.SaveAs('{a}{b}_{c}.{e}'.format(a=plotName,b=prepost,e=ext,c=finalxtitle))

    p1.Delete()
    #p2.Delete()
    p3.Delete()


if __name__ == "__main__":
    parser = optparse.OptionParser(usage='usage: %prog [opts] ', version='%prog 1.0')
    parser.add_option('--outdir'   , type='string'       , default=''    , help='output directory where the postfit plots are saved.')
    parser.add_option('--inclusive', action='store_true' , default=False , help='do it for inclusive plots (i.e. non b-tagged)')
    parser.add_option('--prefit', action='store_true' , default=False , help='make prefit distributions')
    (options, args) = parser.parse_args()

    url=sys.argv[1]
    if not os.path.isfile(url):
        print url,'is not a file'
        sys.exit()

    if not os.path.isdir(options.outdir):
        os.system('mkdir -p {od}'.format(od=options.outdir))
        os.system('cp ~mdunser/public/index.php {od}/'.format(od=options.outdir))


    doPostFitPlot(url,options.inclusive)

    if options.inclusive:
        printTableInclusive(valsAndErrors)
    else:
        for flav in ['em', 'mm', 'ee']:
            tmp_histos = convertValsAndErrorsToHistos(flav, valsAndErrors)
            compareFitResult(plotsPostfit=tmp_histos,plotsPrefit=tmp_histos,
                             plotName=options.outdir+'/nbjets_'+flav,
                             xtitle='b-tagged multiplicity',extraTxt=[flav.replace('m','#mu')])
        printTable(valsAndErrors)

    #sys.exit()
