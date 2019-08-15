import ROOT
import sys

def getBfindEff(h):

    """fit the b-matched multiplicity to extract eb and eb*"""

    ntot=sum([h.GetBinContent(xbin+1) for xbin in range(3)])

    form =  '('
    form += 'x<0.5 ? (1-[0])*(1-[1])'
    form += ' : (x<1.5 ? (1-[0])*[1]+[1]*(1-[0])'
    form += ' : [0]*[1])'
    form += ')*%f'%ntot
    func=ROOT.TF1('effb',form,-0.5,2.5)
    func.SetParLimits(0,0.,1)
    func.SetParLimits(1,0.,1)
    h.Fit(func,'MRQ+','',-0.5,2.5)

    #there is ambiguity in which one is which so sort the final values obtained
    effb=[(func.GetParameter(i),func.GetParError(i)) for i in range(2)]
    effb.sort(key=lambda x: x[0], reverse=True)

    return effb

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
    
c=ROOT.TCanvas('c','c',500,500)
c.SetTopMargin(0.05)
c.SetLeftMargin(0.12)
c.SetRightMargin(0.03)
c.SetBottomMargin(0.1)

print '%25s & %15s & %15s\\\\'%('Sample','$\\varepsilon_{b}$','$\\varepsilon_{b}*$')
fIn=ROOT.TFile.Open(sys.argv[1])
grfind={}
beffVals={}
for key in fIn.GetListOfKeys():
    name=key.GetName()
    if not 'nbjets' in name: continue
    if 'background' in name : continue
    obj=key.ReadObj()
    if not obj.InheritsFrom('TH1') : continue

    beff=getBfindEff(obj)

    tag=''.join(name.split('_')[2:])
    tag=tag.replace('ttbar','')
    if not tag in grfind:
        grfind[tag]=ROOT.TGraphErrors()
        grfind[tag].SetName(tag)
        title=tag
        title=title.replace('cen',' (0-30)')
        title=title.replace('periph',' (30-100)')
        grfind[tag].SetTitle(title)
        grfind[tag]=grfind[tag].Clone(tag+'_find')
        beffVals[tag]={}

    beffVals[tag][name]=beff
    np=grfind[tag].GetN()
    grfind[tag].SetPoint(np,beff[0][0],beff[1][0])
    grfind[tag].SetPointError(np,beff[0][1],beff[1][1])

fIn.Close()

mg=ROOT.TMultiGraph()
c.Clear()
frame=ROOT.TH1F('frame','frame',1,0.3,0.6)
frame.Draw()
frame.GetXaxis().SetRangeUser(0.3,0.6)
frame.GetYaxis().SetRangeUser(0.3,0.6)
frame.GetXaxis().SetTitle('#varepsilon_{b}')
frame.GetYaxis().SetTitle('#varepsilon_{b}*')

colors={'PbPb':1,'PbPbperiph':ROOT.kCyan+1,'PbPbcen':ROOT.kViolet-1}
markers={'PbPb':20,'PbPbperiph':32,'PbPbcen':26}

leg=ROOT.TLegend(0.15,0.7,0.5,0.5)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.035)
for tag in ['PbPb','PbPbperiph','PbPbcen']:
    grfind[tag].SetLineColor(colors[tag])
    grfind[tag].SetMarkerColor(colors[tag])
    grfind[tag].SetMarkerStyle(markers[tag])
    mg.Add(grfind[tag],'p')
    leg.AddEntry(grfind[tag],grfind[tag].GetTitle(),'p')
mg.Draw('p')
leg.Draw()

l=ROOT.TLine()
l.SetLineColor(ROOT.kGray)
l.DrawLine(0.3,0.3,0.6,0.6)

txt=ROOT.TLatex()
txt.SetTextFont(42)

txt.SetNDC(False)
txt.SetTextSize(0.025)

txt.SetNDC(True)
txt.SetTextSize(0.045)
txt.SetTextAlign(12)
txt.DrawLatex(0.12,0.97,'#bf{CMS} #it{preliminary}')
txt.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
txt.DrawLatex(0.95,0.97,'#scale[0.8]{1.6 #mub^{-1} (#sqrt{s_{NN}}=5.02 TeV)}')


c.Modified()
c.Update()
for ext in ['png','pdf']:
    c.SaveAs('bfindeff.%s'%ext)


#print variations
for tag in ['PbPb','PbPbperiph','PbPbcen']:
    beff=beffVals[tag]['nbjets_sell_ttbar_'+tag]
    print '%25s & %3.4f \\pm %3.4f & %3.4f \\pm %3.4f \\\\'%(tag,
                                                             beff[0][0],beff[0][1],
                                                             beff[1][0],beff[1][1])
    for name in beffVals[tag]:
        syst=name.split('_')[1]
        if 'sell' in syst : continue
        beffp=beffVals[tag][name]
        print '\t %25s & %3.4f %3.4f \\\\'%(syst, beffp[0][0]-beff[0][0],beffp[1][0]-beff[1][0])
