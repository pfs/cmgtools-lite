import ROOT
from math import sqrt

ROOT.gROOT.SetBatch()

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

c=ROOT.TCanvas('c','c',900,700)
c.SetBottomMargin(0.12)
c.SetLeftMargin(0.05)
c.SetRightMargin(0.05)
c.SetTopMargin(0.03)

A=208
P2PBPBSCALE=1.#(A*1e-3)**2
C90=1.645

datasetleg={'PbPb': '#scale[1.1]{PbPb, XX.X #mub^{-1}, (#sqrt{s_{NN}}=5.02 TeV)}',
            'pp'  : '#scale[1.1]{pp, 27.4 pb^{-1}, (#sqrt{s}=5.02 TeV)}'}

theo={
    'PbPb' : [
        ('#left(#splitline{CT14+EPPS16}{NLO MCFM}#right) #upoint K_{NNLO+NNLL}(Top++)',     59.0, 5.3, 5.3, 2.0, 1.6),
        #('#left(#splitline{CT10+EPS09}{NLO MCFM}#right) #upoint K_{NNLO+NNLL}(Top++)',      57.5, 3.3, 4.3, 2.0, 1.5)
        
        ],
    'pp'  : [
        ('#splitline{CT14}{NNLO+NNLL Top++}',     75., 15.3/C90, 17.2/C90, 9.5, 7.0),
        ('#splitline{CT10}{NNLO+NNLL Top++}',     75., 14.3/C90, 17.3/C90, 9.3, 6.9),
        ]
    }

exp={'PbPb':
         [ 
        ('2l'     ,  59.,  3.,  3.,  10., 10.),
        ('2l+jets',  59.,  3.,  3.,  10., 10.),
        ],
     'pp':
         [
        ('1l/2l+jets #scale[0.8]{JHEP 03 (2016) 115}', 69.5, 6.1, 6.1, sqrt(5.82**2+6.1**2), sqrt(5.82**2+6.1**2)), ## summing up syst and lumi for the total
        #('l+jets #scale[0.8]{EPJC 77 (2017) 15}',  228.5,  3.8, 3.8, 15.43, 15.43)
        ]
     }

colors  = {'PbPb':1            , 'pp':1} #ROOT.TColor.GetColor('#fc8d59')}
fill    = {'PbPb':ROOT.kAzure+7, 'pp':ROOT.kGreen-5} 
markers = {'PbPb':20           , 'pp':20}

dy=3.5
for key in exp:
    dy+=len(exp[key])
    if key in theo:
        dy+=len(theo[key])-1

frame=ROOT.TH2F('frame',';#sigma [nb]',1,0,115,1,-1.25,dy)
frame.GetXaxis().SetTitleOffset(0.85)
frame.GetXaxis().SetTitleSize(0.06)
frame.GetXaxis().SetLabelSize(0.04)
frame.GetYaxis().SetNdivisions(0)
frame.Draw()

cms=ROOT.TLatex()
cms.SetTextFont(42)
cms.SetTextSize(0.06)
cms.SetNDC()
cms.DrawLatex(0.82,0.9,'#bf{CMS}') # #it{Preliminary}')

labels=ROOT.TLatex()
labels.SetTextFont(42)
labels.SetTextSize(0.05)

theoPDFGr,theoTotGr=[],[]
expStatGr,expTotalGr=[],[]
imeas=0

dy=0
for key in exp:

    if imeas>0: dy+=1
    expdy=len(exp[key])
    theody=0

    #check if there are theory references acompanying
    if key in theo:

        theody=len(theo[key])

        for i in xrange(0,len(theo[key])):
            title,mu,pdfDn,pdfUp,scaleDn,scaleUp=theo[key][i]
            if key=='pp':
                mu      *= P2PBPBSCALE
                pdfUp   *= P2PBPBSCALE
                pdfDn   *= P2PBPBSCALE
                scaleUp *= P2PBPBSCALE
                scaleDn *= P2PBPBSCALE
            totDn=ROOT.TMath.Sqrt(pdfDn**2+scaleDn**2)
            totUp=ROOT.TMath.Sqrt(pdfUp**2+scaleUp**2)
            print key,mu,totDn,totUp
            theoPDFGr.append( ROOT.TGraph() )
            theoPDFGr[-1].SetTitle(title)
            theoPDFGr[-1].SetName('theopdf_%s_%d'%(key,i))
            #theoPDFGr[-1].SetFillColor(fill[key]-i)
            theoPDFGr[-1].SetFillColorAlpha(fill[key]-i, 0.5);
            theoPDFGr[-1].SetFillStyle(1001)

            theoTotGr.append( theoPDFGr[-1].Clone('theotot_%s_%d'%(key,i)) )
            #theoTotGr[-1].SetFillColor(fill[key]-i)
            theoTotGr[-1].SetFillColorAlpha(fill[key]-i, 0.5);
            theoTotGr[-1].SetFillStyle(3001)

            ymin,ymax=dy+0.5,dy+expdy+0.5
            if i>0: ymin,ymax=dy+expdy+i,dy+expdy+1+i

            theoPDFGr[-1].SetPoint(0,mu-pdfDn,ymin)
            theoPDFGr[-1].SetPoint(1,mu-pdfDn,ymax)
            theoPDFGr[-1].SetPoint(2,mu+pdfUp,ymax)
            theoPDFGr[-1].SetPoint(3,mu+pdfUp,ymin)
            theoPDFGr[-1].SetPoint(4,mu-pdfDn,ymin)

            theoTotGr[-1].SetPoint(0,mu-totDn,ymin)
            theoTotGr[-1].SetPoint(1,mu-totDn,ymax)
            theoTotGr[-1].SetPoint(2,mu+totUp,ymax)
            theoTotGr[-1].SetPoint(3,mu+totUp,ymin)
            theoTotGr[-1].SetPoint(4,mu-totDn,ymin)
                
            theoTotGr[-1].Draw('f')
            theoPDFGr[-1].Draw('f')

            labels.DrawLatex(68,ymax-1.07,'#scale[0.6]{%s}'%title)

    dy+=expdy+theody
    labels.DrawLatex(5,dy-0.0,'#scale[0.7]{#bf{%s}}'%datasetleg[key])

    ROOT.gStyle.SetEndErrorSize(5)
    for i in xrange(0,len(exp[key])):
        expStatGr.append( ROOT.TGraphAsymmErrors() )
        expStatGr[-1].SetMarkerStyle(0)
        expStatGr[-1].SetLineWidth(3)
        expStatGr[-1].SetLineColor(colors[key])
        expStatGr[-1].SetFillStyle(0)

        expTotalGr.append( ROOT.TGraphAsymmErrors() )
        expTotalGr[-1].SetLineWidth(3)
        expTotalGr[-1].SetMarkerStyle(markers[key])
        expTotalGr[-1].SetMarkerColor(colors[key])
        expTotalGr[-1].SetLineColor(colors[key])
        expTotalGr[-1].SetFillStyle(0)

        title,xsec,statUp,statDn,totalUp,totalDn=exp[key][i]
        if key=='pp':
            xsec    *= P2PBPBSCALE
            statUp  *= P2PBPBSCALE
            statDn  *= P2PBPBSCALE
            totalUp *= P2PBPBSCALE
            totalDn *= P2PBPBSCALE

        yanchor=imeas+i+1 if imeas==0 else imeas+i

        expStatGr[-1].SetPoint(0,xsec,yanchor)
        expStatGr[-1].SetPointError(0,statUp,statDn,0,0)
        expStatGr[-1].Draw('p[]')

        expTotalGr[-1].SetPoint(0,xsec,yanchor)
        expTotalGr[-1].SetPointError(0,totalUp,totalDn,0,0)
        expTotalGr[-1].Draw('p')

        labels.DrawLatex(5,yanchor,'#scale[0.8]{%s}'%title)

    imeas += dy+2        

grlegstat=ROOT.TGraphErrors()
grlegstat.SetMarkerStyle(0)
grlegstat.SetLineWidth(2)
grlegstat.SetLineColor(1) #ROOT.kGray)
grlegstat.SetFillStyle(0)
grlegstat.SetPoint(0,97.5,0.5)
grlegstat.SetPointError(0,2.5,0)

grlegtot=ROOT.TGraphErrors()
grlegtot.SetLineWidth(2)
grlegtot.SetLineColor(1) #ROOT.kGray)
grlegtot.SetFillStyle(0)
grlegtot.SetPoint(0,97.5,0.5)
grlegtot.SetPointError(0,8,0)
grlegtot.SetMarkerStyle(20)
grlegtot.SetMarkerColor(1)

thlegpdf=ROOT.TGraphErrors()
thlegpdf.SetMarkerStyle(0)
#thlegpdf.SetFillColor(ROOT.kGray)
#thlegpdf.SetFillStyle(1001)
thlegpdf.SetFillColorAlpha(ROOT.kGray,0.5)
thlegpdf.SetFillStyle(1001)
thlegpdf.SetPoint(0,97.5,-0.5)
thlegpdf.SetPointError(0,2.5,0.25)

thlegtot=ROOT.TGraphErrors()
thlegtot.SetMarkerStyle(0)
thlegtot.SetFillColorAlpha(ROOT.kGray,0.5)
thlegtot.SetFillStyle(3001)
#thlegtot.SetFillColor(ROOT.kGray)
#thlegtot.SetFillStyle(3001)
thlegtot.SetPoint(0,97.5,-0.5)
thlegtot.SetPointError(0,8,0.25)


grlegstat.Draw('p[]')
grlegtot.Draw('p')
labels.SetTextSize(0.045)
#labels.DrawLatex(84,0.75,'#scale[0.6]{#color[15]{#it{Exp. unc.: stat  stat#oplussyst}}}')
labels.DrawLatex(84,0.75,'#scale[0.6]{#it{Exp. unc.: stat  stat#oplussyst}}')

thlegpdf.Draw('e2')
thlegtot.Draw('e2')
labels.SetTextSize(0.045)
#labels.DrawLatex(85.5,-0.15,'#scale[0.6]{#color[15]{#it{Th. unc.: pdf  pdf#oplusscales}}}')
labels.DrawLatex(85.5,-0.15,'#scale[0.6]{#it{Th. unc.: pdf  pdf#oplusscales}}')

#labels.DrawLatex(5,-0.6,'#scale[0.7]{#it{All theory predictions use: #mu_{F}=#mu_{R}=m_{t}=172.5 GeV, #alpha_{s} = 0.1180}}')

c.RedrawAxis()
c.Modified()
c.Update()
c.SaveAs('~/www/private/heavyIons/plots/card_inputs/2019-05-08-test/xsecsummary.pdf')
c.SaveAs('~/www/private/heavyIons/plots/card_inputs/2019-05-08-test/xsecsummary.png')
## raw_input()

cms.DrawLatex(0.755,0.85,'#scale[0.75]{#it{Preliminary}}')
c.Modified()
c.Update()
c.SaveAs('~/www/private/heavyIons/plots/card_inputs/2019-05-08-test/xsecsummary_prel.pdf')
c.SaveAs('~/www/private/heavyIons/plots/card_inputs/2019-05-08-test/xsecsummary_prel.png')
## raw_input()
