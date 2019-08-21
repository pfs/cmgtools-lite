import ROOT
import numpy as np
import sys

fIn=ROOT.TFile.Open(sys.argv[1])
pol='pol6'

origmva={'S':fIn.Get('bdtfine_ttbar'),
         'B':fIn.Get('bdtfine_data_comb')}

p=np.array([x/100. for x in np.arange(5, 100, 10)])
q=np.array([0. for i in range(len(p))])
origmva['B'].GetQuantiles(len(p),q,p)

cdfInv=ROOT.TGraph()
cdfInv.SetMarkerStyle(20)
cdfInv.SetPoint(0,-1,0)
for i in xrange(0,len(p)):
    cdfInv.SetPoint(i+1,q[i],p[i])
cdfInv.SetPoint(len(p)+1,1,1)
cdfInv.Fit(pol,'CFW')
f=cdfInv.GetFunction(pol)


newmva={}
for t in ['S','B']:
    newmva[t]=ROOT.TH1F(t+'_new',t+'_new',10,0,1)
    newmva[t].Reset('ICE')
    for xbin in range(origmva[t].GetNbinsX()):
        val=origmva[t].GetXaxis().GetBinCenter(xbin+1)
        newval=f.Eval(val)
        newmva[t].Fill(newval,origmva[t].GetBinContent(xbin+1))

c=ROOT.TCanvas('c','c',1000,500)
c.Divide(2,1)
c.cd(1)
cdfInv.Draw('ap')
c.cd(2)
newmva['S'].SetLineWidth(2)
newmva['S'].SetLineColor(2)
newmva['S'].Draw('hist')
newmva['B'].SetLineColor(ROOT.kOrange)
newmva['S'].SetLineWidth(2)
newmva['B'].Draw('histsame')
c.cd()
c.Modified()
c.Update()

print '-'*50
print 'The following formula should be added to functions.cc'
print '-'*50
print 'TF1 *f=new TF1("bdtcdfinv",%s",-1,1)'%f.GetFormula().GetExpFormula();
for i in range(f.GetNpar()):
    print 'f->SetParameter(%d,%f);'%(i,f.GetParameter(i))
print '-'*50
raw_input()
    
