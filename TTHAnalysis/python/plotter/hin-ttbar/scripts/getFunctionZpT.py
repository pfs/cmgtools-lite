import ROOT, os, math

ROOT.gStyle.SetOptStat(0)

file_inc= ROOT.TFile('~/www/private/heavyIons/plots/dy_plots/2019-06-25-mixedMC_newIso_withSF-mm/dyllpt_AND_dyleppt_AND_dyl1pt_AND_dyl2pt_AND_dysphericity_AND_dydphi_AND_dyllm.root')
file_chi= ROOT.TFile('~/www/private/heavyIons/plots/dy_plots/2019-06-25-mixedMC_newIso_withSF_centralityHi-mm/nleptons_AND_l1iso02_AND_l2iso02_AND_dyllpt_AND_dyleppt_AND_dyl1pt_AND_dyl2pt_AND_dysphericity_AND_dydphi_AND_dyllm_AND_lepiso02_AND_l1d0_AND_l2d0.root')
file_clo= ROOT.TFile('~/www/private/heavyIons/plots/dy_plots/2019-06-25-mixedMC_newIso_withSF_centralityLo-mm/nleptons_AND_l1iso02_AND_l2iso02_AND_dyllpt_AND_dyleppt_AND_dyl1pt_AND_dyl2pt_AND_dysphericity_AND_dydphi_AND_dyllm_AND_lepiso02_AND_l1d0_AND_l2d0.root')


plotname = 'dyllpt_zg'

h_mc_inc = file_inc.Get('dyllpt_zg'  ); h_mc_inc.SetLineWidth(2); h_mc_inc.SetMarkerColor(ROOT.kOrange+7); h_mc_inc.SetLineColor(ROOT.kOrange+6); h_mc_inc.SetFillStyle(0)
h_da_inc = file_inc.Get('dyllpt_data'); h_da_inc.SetLineWidth(2); h_da_inc.SetMarkerColor(ROOT.kBlack   ); h_da_inc.SetLineColor(ROOT.kBlack   )
h_da_chi = file_chi.Get('dyllpt_data'); h_da_chi.SetLineWidth(2); h_da_chi.SetMarkerColor(ROOT.kAzure+2 ); h_da_chi.SetLineColor(ROOT.kAzure+1 )
h_da_clo = file_clo.Get('dyllpt_data'); h_da_clo.SetLineWidth(2); h_da_clo.SetMarkerColor(ROOT.kSpring-6); h_da_clo.SetLineColor(ROOT.kSpring-5)

h_mc_inc.Scale(1./h_mc_inc.Integral()); h_mc_inc.SetTitle('')
h_da_inc.Scale(1./h_da_inc.Integral()); h_da_inc.SetTitle('')
h_da_chi.Scale(1./h_da_chi.Integral()); h_da_chi.SetTitle('')
h_da_clo.Scale(1./h_da_clo.Integral()); h_da_clo.SetTitle('')


tmp_canv= ROOT.TCanvas('whatever','',600,750)
tmp_canv.GetPad(0).SetTopMargin   (0.04);
tmp_canv.GetPad(0).SetBottomMargin(0.32);
tmp_canv.GetPad(0).SetLeftMargin  (0.18);
tmp_canv.GetPad(0).SetRightMargin (0.04);
tmp_canv.GetPad(0).SetTickx(1);
tmp_canv.GetPad(0).SetTicky(1);

tmp_canv.cd(1)

h_mc_inc.GetXaxis().SetLabelSize(0.)
h_mc_inc.GetYaxis().SetTitle('pdf')
h_mc_inc.Draw('hist')
h_da_inc.Draw('hist same')
h_da_chi.Draw('hist same')
h_da_clo.Draw('hist same')

leg=ROOT.TLegend(0.40, 0.65, 0.60, 0.85)
leg.SetTextSize(0.04)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.AddEntry(h_mc_inc, 'MC (inc.)', 'l')
leg.AddEntry(h_da_inc, 'data (inc.)', 'l')
leg.AddEntry(h_da_chi, 'data, low centrality', 'l')
leg.AddEntry(h_da_clo, 'data, high centrality', 'l')
leg.Draw('same')


ratio_mc_inc = h_mc_inc.Clone('ratio_mc_inc'); ratio_mc_inc.Divide(h_mc_inc)
ratio_da_inc = h_da_inc.Clone('ratio_da_inc'); ratio_da_inc.Divide(h_mc_inc)
ratio_da_chi = h_da_chi.Clone('ratio_da_chi'); ratio_da_chi.Divide(h_mc_inc)
ratio_da_clo = h_da_clo.Clone('ratio_da_clo'); ratio_da_clo.Divide(h_mc_inc)

ratiopad = ROOT.TPad('padratio', '', 0.,0.,1.,0.28)
ratiopad.SetTopMargin   (0.00);
ratiopad.SetBottomMargin(0.32);
ratiopad.SetLeftMargin  (0.18);
ratiopad.SetRightMargin (0.04);
ratiopad.SetTickx(1);
ratiopad.SetTicky(1);

ratiopad.Draw()
ratiopad.cd()
ratio_da_inc.Draw('axis')
ratio_da_inc.GetXaxis().SetLabelSize(0.12)
ratio_da_inc.GetXaxis().SetTitleOffset(1.20)
ratio_da_inc.GetXaxis().SetTitleSize(0.12)
ratio_da_inc.GetYaxis().SetTitle('ratio w/r/t MC')
ratio_da_inc.GetYaxis().SetTitleSize(0.09)
ratio_da_inc.GetYaxis().SetTitleOffset(0.8)
ratio_da_inc.GetYaxis().SetRangeUser(0., 2.)

line = ROOT.TLine()
line.SetLineWidth(2)
line.DrawLine(ratio_da_inc.GetXaxis().GetXmin(),1.,ratio_da_inc.GetXaxis().GetXmax(),1.)
ratio_mc_inc .Draw('pe same')
ratio_da_inc .Draw('pe same')
ratio_da_chi .Draw('pe same')
ratio_da_clo .Draw('pe same')

erf = ROOT.TF1('erf', '[0]*TMath::Erf((x-[1])/[2])', 0., 150.)
erf.SetLineColor(ROOT.kGray+2)
erf.SetLineWidth(2)
erf.SetParameters(1., -5., 15.)

ff = ROOT.TF1('myfunc', '[0]*x*x*x + [1]*x*x + [2]*x + [3]', -50., 150.)
ff.SetParameter(0, -0.00050785)
ff.SetParameter(1, -0.010089)
ff.SetParameter(2, 1.2293)
ff.SetParameter(3, 0.8)
ff.SetLineColor(ROOT.kAzure+2)
ff.SetLineWidth(2)

lan = ROOT.TF1('lan', 'TMath::Landau(x,[0],[1],0)',0., 50.)
lan.SetParameters(1.3, 10.)
lan.SetLineColor(ROOT.kPink+2)
lan.SetLineWidth(2)

wat = ROOT.TF1('wat', 'pol1(0)*expo(2)', 0., 100.)
wat.SetLineColor(ROOT.kPink+2)
wat.SetLineWidth(2)

ratio_da_inc .Fit('pol3', '', '', 0, 50.)
##ff.Draw('same')

#tmp_canv.cd(0)
lat = ROOT.TLatex()
lat.SetNDC(); lat.SetTextSize(0.08); lat.SetTextFont(42)
#x1=erf.GetParameter(1)
#lat.DrawLatex(0.4, 0.35, 'fit: {x0:.2f}*Erf( (p_{{T}}(ll) {s} {x1:.2f}) / {x2:.2f} ) '.format(x0=erf.GetParameter(0), x1=x1 if x1  >0 else abs(x1), x2=erf.GetParameter(2), s='-'if x1 >0. else '+' ))

tmp_canv.SaveAs('~/www/private/heavyIons/plots/dy_plots/ptllFunction_mixedWithSFs_2019-06-25.png')
tmp_canv.SaveAs('~/www/private/heavyIons/plots/dy_plots/ptllFunction_mixedWithSFs_2019-06-25.pdf')
