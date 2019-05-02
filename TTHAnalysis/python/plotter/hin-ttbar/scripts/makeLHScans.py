import ROOT, os, subprocess, sys, optparse, datetime
from array import array


def graphStyle(graph, mode):
    markerstyle = 20    if mode == 'bdt' else 21             if mode == 'sph' else 22
    color = ROOT.kBlack if mode == 'obs' else ROOT.kOrange-3 if mode == 'sph' else ROOT.kAzure+2
    style = 1
    graph.SetMarkerStyle(markerstyle)
    graph.SetMarkerColor(color)
    graph.SetLineColor  (color)
    graph.SetLineWidth  (3)
    graph.SetLineStyle  (style)
    graph.SetMarkerSize(0.8)
    graph.GetXaxis().SetTitle('#hat{#mu}_{t#bar{t}}')
    graph.GetYaxis().SetTitle('-2 #Delta ln L')
    ##graph.GetYaxis().SetRangeUser(0.0, 1.5)
    ##graph.GetXaxis().SetRangeUser(0.0, 2.5)
    return graph


def getGraph(tree, mode = 'obs'):
    xsec = 1.
    n = tree.Draw('2*dnllval:(ttbar_mu*{xs})'.format(xs=xsec), '', 'goff')
    vals = []
    minXsec = -999.
    for ev in tree:
        vals.append( [ev.ttbar_mu*xsec, (2.*ev.dnllval)] )
        if 2.*ev.dnllval == 0.:
            minXsec = ev.ttbar_mu*xsec
    vals = sorted(vals)
    graph = ROOT.TGraph(len(vals), array('d', [x[0] for x in vals]), array('d', [y[1] for y in vals]) )
    graph = graphStyle(graph, mode)
    graph_alt = ROOT.TGraph(len(vals), array('d', [y[1] for y in vals]), array('d', [x[0] for x in vals]) )
    err  = graph_alt.Eval(1.)
    return graph, minXsec


def makeScans(n = 100):

    f_ex1 = ROOT.TFile('../datacards_2019-04-29_withIsoApplied_scaledToEff/bdt/ttbar/fitresults_allFlavors_scan.root'       , 'READ'); l_ex1 = f_ex1.Get('fitresults')
    f_ex2 = ROOT.TFile('../datacards_2019-04-29_withIsoApplied_scaledToEff/sphericity/ttbar/fitresults_allFlavors_scan.root', 'READ'); l_ex2 = f_ex2.Get('fitresults')
    
    (g_ex1, xsec_ex1) = getGraph(l_ex1, 'bdt')
    (g_ex2, xsec_ex2) = getGraph(l_ex2, 'sph')

    canv = ROOT.TCanvas('foo', 'bar', 600, 600)
    canv.cd()
    canv.GetPad(0).SetLeftMargin(0.15)
    canv.GetPad(0).SetRightMargin(0.05)
    canv.GetPad(0).SetTopMargin(0.07)
    canv.GetPad(0).SetBottomMargin(0.15)

    xmax = 2.0
    dummy = ROOT.TH1F('foo', '', 100, 0., xmax)
    dummy.GetYaxis().SetRangeUser(0., 26.)
    dummy.GetXaxis().SetTitle(g_ex1.GetXaxis().GetTitle())
    dummy.GetYaxis().SetTitle(g_ex1.GetYaxis().GetTitle())

    dummy.GetXaxis().SetTitleSize(0.06)
    dummy.GetYaxis().SetTitleSize(0.06)

    dummy.GetXaxis().SetTitleOffset(1.00)
    dummy.GetYaxis().SetTitleOffset(1.00)

    dummy.GetXaxis().SetLabelSize(0.04)
    dummy.GetYaxis().SetLabelSize(0.04)

    dummy.GetXaxis().SetRangeUser(0., xmax)
    dummy.GetYaxis().SetNdivisions(505)
    dummy.GetXaxis().SetNdivisions(510)

    ROOT.gStyle.SetOptStat(0)
    dummy.Draw('AXIS')

    mg = ROOT.TMultiGraph()
    mg.Add(g_ex1)
    mg.Add(g_ex2)
    mg.Draw('l')

    leg = ROOT.TLegend(0.55, 0.64, 0.80, 0.85)
    #leg.SetFillStyle(0)
    leg.SetTextSize(0.04)
    leg.AddEntry(g_ex1, 'leptonic BDT', 'l')
    leg.AddEntry(g_ex2, 'sphericity'  , 'l')
    leg.SetLineColor(ROOT.kWhite)
    leg.Draw('same')

    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextAlign(31)
    lat.SetTextSize(0.025)
    lat.SetTextColor(ROOT.kGray+1)

    line = ROOT.TLine()
    line.SetLineStyle(7)
    line.SetLineWidth(2)
    line.SetLineColor(ROOT.kGray+1)

    yvals_sig = [0.135, 0.235, 0.395, 0.62, 0.905]
    for isig in [1,2,3,4,5]:
        line.DrawLine(0., float(isig)**2, xmax, float(isig)**2)
        lat.DrawLatex(0.15, yvals_sig[isig-1], '{n} #sigma'.format(n=isig))

    lat1 = ROOT.TLatex()
    lat1.SetNDC()
    lat1.SetTextColor(ROOT.kBlack)
    lat1.SetTextFont(42)
    lat1.SetTextSize(0.035)
    lat1.DrawLatex(0.10, 0.94, '#bf{CMS} #it{Preliminary}')
    lat1.DrawLatex(0.68, 0.94, 'XX #mub^{-1} (5.02 TeV)')


    tarr = ROOT.TArrow(0.43, 0.25, 0.43, 0.15, 0.03, ">")
    #tarr.SetNDC()
    tarr.SetLineWidth(2)
    lat.SetTextAlign(22)
    lat.SetTextFont(42)
    lat.SetTextSize(0.031)
    ## arry1, arry2 = -0.35, -0.05 #0.12, 0.30
    arry1, arry2 = 4.05, 0.12
    laty = 0.25

    ## no observed yet tarr.SetLineColor(ROOT.kBlack)
    ## no observed yet tarr.DrawArrow(xsec_obs, arry1, xsec_obs, arry2)
    ## no observed yet lat .DrawLatex(0.45, laty, '{obs:.2f}'.format(obs=xsec_obs))

    ## tarr.SetLineColor(ROOT.kRed-4)
    ## lat .SetTextColor(ROOT.kRed-4)
    ## tarr.DrawArrow(xsec_ex1, arry1, xsec_ex1, arry2)
    ## lat .DrawLatex(0.75, laty, '{obs:.2f}'.format(obs=xsec_ex1))

    ## tarr.SetLineColor(ROOT.kGreen+2)
    ## lat .SetTextColor(ROOT.kGreen+2)
    ## tarr.DrawArrow(xsec_ex2, arry1, xsec_ex2, arry2)
    ## lat .DrawLatex(0.40, laty, '{obs:.2f}'.format(obs=xsec_ex2))

    ## tarr.SetLineColor(ROOT.kBlack)
    ## lat .SetTextColor(ROOT.kBlack)
    ## tarr.DrawArrow(xsec_obs, arry1, xsec_obs, arry2)
    ## lat .DrawLatex(0.57, laty, '{obs:.2f}'.format(obs=xsec_obs))
    

    date = datetime.date.today().isoformat()
    canv.SaveAs('/afs/cern.ch/user/m/mdunser/www/private/heavyIons/plots/lhScans/{date}-lhScan.pdf'.format(date=date))
    canv.SaveAs('/afs/cern.ch/user/m/mdunser/www/private/heavyIons/plots/lhScans/{date}-lhScan.png'.format(date=date))


if __name__ == '__main__':
    parser = optparse.OptionParser(usage='usage: %prog [opts] ', version='%prog 1.0')
    parser.add_option('--makeScans', '--scans',  action='store_true', dest='makeScans'        , default=False, help='make likelihood scans')
    global opts
    (opts, args) = parser.parse_args()


    makeScans()
