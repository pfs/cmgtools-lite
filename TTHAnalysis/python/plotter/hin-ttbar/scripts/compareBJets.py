import ROOT, sys, optparse, datetime

## usage:
##     python compareBJets.py -i /eos/cms/store/cmst3/group/hintt/PbPb2018_skim13August/TTJets_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8.root  --postfix finalJets


date = datetime.date.today().isoformat()

parser = optparse.OptionParser(usage='usage: %prog [opts] ', version='%prog 1.0')
parser.add_option('-i', '--infile' , type='string', default='', help='infile name of the ttbar tree')
parser.add_option('-p', '--postfix', type='string', default='', help='postfix for the output file name')
(options, args) = parser.parse_args()

if not options.infile:
    print 'give an input ntuple file of the ttbar process!!'
    exit()

ROOT.gROOT.SetBatch()

ROOT.gStyle.SetOptStat(0)
infile = ROOT.TFile(options.infile,'read')

nom    = ROOT.TH1F('bjets' , ''      , 4, -0.5, 3.5); nom    .SetLineWidth(2); nom    .SetLineColor(ROOT.kBlack)
jecup  = ROOT.TH1F('jecup' , 'jecup' , 4, -0.5, 3.5); jecup  .SetLineWidth(2); jecup  .SetLineColor(ROOT.kAzure+8)
jecdn  = ROOT.TH1F('jecdn' , 'jecdn' , 4, -0.5, 3.5); jecdn  .SetLineWidth(2); jecdn  .SetLineColor(ROOT.kAzure-4)
jerup  = ROOT.TH1F('jerup' , 'jerup' , 4, -0.5, 3.5); jerup  .SetLineWidth(2); jerup  .SetLineColor(ROOT.kRed+1)
jerdn  = ROOT.TH1F('jerdn' , 'jerdn' , 4, -0.5, 3.5); jerdn  .SetLineWidth(2); jerdn  .SetLineColor(ROOT.kRed-4)
bup    = ROOT.TH1F('bup'   , 'bup'   , 4, -0.5, 3.5); bup    .SetLineWidth(2); bup    .SetLineColor(ROOT.kGreen-4)
bdn    = ROOT.TH1F('bdn'   , 'bdn'   , 4, -0.5, 3.5); bdn    .SetLineWidth(2); bdn    .SetLineColor(ROOT.kGreen-2)
udsgup = ROOT.TH1F('udsgup', 'udsgup', 4, -0.5, 3.5); udsgup .SetLineWidth(2); udsgup .SetLineColor(ROOT.kOrange)
udsgdn = ROOT.TH1F('udsgdn', 'udsgdn', 4, -0.5, 3.5); udsgdn .SetLineWidth(2); udsgdn .SetLineColor(ROOT.kOrange-3)
quendn = ROOT.TH1F('quendn', 'quendn', 4, -0.5, 3.5); quendn .SetLineWidth(2); quendn .SetLineColor(ROOT.kGray)


tree = infile.Get('tree')

tree.Draw('nbjet_sel>>bjets')
tree.Draw('nbjet_sel_jecup>>jecup')
tree.Draw('nbjet_sel_jecdn>>jecdn')
tree.Draw('nbjet_sel_jerup>>jerup')
tree.Draw('nbjet_sel_jerdn>>jerdn')
tree.Draw('nbjet_sel_bup>>bup'  )
tree.Draw('nbjet_sel_bdn>>bdn'  )
tree.Draw('nbjet_sel_udsgup>>udsgup')
tree.Draw('nbjet_sel_udsgdn>>udsgdn')
tree.Draw('nbjet_sel_quenchdn>>quendn')

c1 = ROOT.TCanvas('foo','n_{b-tags} with systematic variations', 800, 1200)
c1.GetPad(0).SetTopMargin(0.09)
c1.GetPad(0).SetBottomMargin(0.35)
c1.GetPad(0).SetLeftMargin(0.12)
c1.GetPad(0).SetRightMargin(0.04)
c1.GetPad(0).SetTickx(1)
c1.GetPad(0).SetTicky(1)

scaling = {}

nominal_nb = [ nom.GetBinContent(1), nom.GetBinContent(2), nom.GetBinContent(3) ]

for h in [jecup, jecdn, jerup, jerdn, bup, bdn, udsgup, udsgdn, quendn]:
    for nb in range(3):
        scaling[h.GetName()+'_'+str(nb)+'b'] = float('{a:.3f}'.format(a=h.GetBinContent(nb+1)/nominal_nb[nb]))

for i,j in sorted(sorted(scaling.items()), key=lambda x: 0 if 'dn' in x else -1):
    print i,j

nom    .Scale(1./nom    .Integral())
jecup  .Scale(1./jecup  .Integral())
jecdn  .Scale(1./jecdn  .Integral())
jerup  .Scale(1./jerup  .Integral())
jerdn  .Scale(1./jerdn  .Integral())
bup    .Scale(1./bup    .Integral())
bdn    .Scale(1./bdn    .Integral())
udsgup .Scale(1./udsgup .Integral())
udsgdn .Scale(1./udsgdn .Integral())
quendn .Scale(1./quendn .Integral())


nom.GetYaxis().SetRangeUser(0, 0.75)
nom.GetXaxis().SetTitle('n_{b-tags}')
nom.GetYaxis().SetTitle('pdf')

leg = ROOT.TLegend(0.6, 0.50, 0.87, 0.87)
leg.SetLineStyle(0)
leg.SetLineWidth(0)
leg.AddEntry(nom  , 'nominal t#bar{t} (mixed)', 'l')
leg.AddEntry(jecup, 'JEC up'  , 'l')
leg.AddEntry(jecdn, 'JEC down', 'l')
leg.AddEntry(jerup, 'JER up'  , 'l')
leg.AddEntry(jerdn, 'JER down', 'l')
leg.AddEntry(bup, 'b-tagging (b) up'  , 'l')
leg.AddEntry(bdn, 'b-tagging (b) down', 'l')
leg.AddEntry(udsgup, 'b-tagging (udsg) up'  , 'l')
leg.AddEntry(udsgdn, 'b-tagging (udsg) down', 'l')
leg.AddEntry(quendn, 'quenching down', 'l')

nom    .Draw('hist  ')
jecup  .Draw('hist same ')
jecdn  .Draw('hist same ')
jerup  .Draw('hist same ')
jerdn  .Draw('hist same ')
bup    .Draw('hist same ')
bdn    .Draw('hist same ')
udsgup .Draw('hist same ')
udsgdn .Draw('hist same ')
quendn .Draw('hist same ')
nom    .Draw('hist same ')

leg.Draw('same')

lat = ROOT.TLatex()
lat.SetNDC()
lat.SetTextSize(0.05)
lat.SetTextFont(42)
lat.DrawLatex(0.2, 0.95, 'n_{b-tags} with systematic variations')
            

pad2 = ROOT.TPad("pad2","pad2",0,0.,1,0.9)

pad2.SetTopMargin(0.70)
pad2.SetRightMargin(0.04)
pad2.SetLeftMargin(0.12)
pad2.SetBottomMargin(0.10)
pad2.SetFillColor(0)
pad2.SetGridy(0)
pad2.SetFillStyle(0)
pad2.SetTicky(1)
pad2.SetTickx(1)           

pad2.Draw()
pad2.cd()

ratio_nom    = nom    .Clone(nom    .GetName()+'_ratio')
ratio_jecup  = jecup  .Clone(jecup  .GetName()+'_ratio')
ratio_jecdn  = jecdn  .Clone(jecdn  .GetName()+'_ratio')
ratio_jerup  = jerup  .Clone(jerup  .GetName()+'_ratio')
ratio_jerdn  = jerdn  .Clone(jerdn  .GetName()+'_ratio')
ratio_bup    = bup    .Clone(bup    .GetName()+'_ratio')
ratio_bdn    = bdn    .Clone(bdn    .GetName()+'_ratio')
ratio_udsgup = udsgup .Clone(udsgup .GetName()+'_ratio')
ratio_udsgdn = udsgdn .Clone(udsgdn .GetName()+'_ratio')
ratio_quendn = quendn .Clone(quendn .GetName()+'_ratio')

ratio_nom    .Divide(nom)
ratio_jecup  .Divide(nom)
ratio_jecdn  .Divide(nom)
ratio_jerup  .Divide(nom)
ratio_jerdn  .Divide(nom)
ratio_bup    .Divide(nom)
ratio_bdn    .Divide(nom)
ratio_udsgup .Divide(nom)
ratio_udsgdn .Divide(nom)
ratio_quendn .Divide(nom)

ratio_nom    .GetYaxis().SetRangeUser(0.90,1.10)
ratio_nom    .GetYaxis().SetTitle('ratio')
ratio_nom    .GetYaxis().SetNdivisions(505)

ratio_nom    .Draw('pl')
ratio_jecup  .Draw('pl same')
ratio_jecdn  .Draw('pl same')
ratio_jerup  .Draw('pl same')
ratio_jerdn  .Draw('pl same')
ratio_bup    .Draw('pl same')
ratio_bdn    .Draw('pl same')
ratio_udsgup .Draw('pl same')
ratio_udsgdn .Draw('pl same')
ratio_quendn .Draw('pl same')


c1.SaveAs('~/www/private/heavyIons/plots/jetPlots/{d}_nbSystematics{pf}.png'.format(d=date,pf='_'+options.postfix if options.postfix else ''))
c1.SaveAs('~/www/private/heavyIons/plots/jetPlots/{d}_nbSystematics{pf}.pdf'.format(d=date,pf='_'+options.postfix if options.postfix else ''))
            
            
            
            
            
