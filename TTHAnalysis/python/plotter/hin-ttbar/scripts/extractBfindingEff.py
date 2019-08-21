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
    h.Fit(func,'MRQ0+','',-0.5,2.5)

    #there is ambiguity in which one is which so sort the final values obtained
    effb=[(func.GetParameter(i),func.GetParError(i)) for i in range(2)]
    effb.sort(key=lambda x: x[0], reverse=True)

    return effb

def showBfindingSummary(grfind):

    """makes a summary plot for the AN """

    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.03)
    c.SetBottomMargin(0.1)

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


def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    
    #compute the b-finding efficiencies for the nominal expectations and variations
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

        varName=name.split('_')[1]
        beffVals[tag][varName]=beff
        np=grfind[tag].GetN()
        grfind[tag].SetPoint(np,beff[0][0],beff[1][0])
        grfind[tag].SetPointError(np,beff[0][1],beff[1][1])
    fIn.Close()

    #show a summary
    showBfindingSummary(grfind)
    
    #dump a table to the screen
    print '%20s %49s %78s'%('Variation','$\\varepsilon_{b}$','$\\varepsilon_{b}*$')

    print '%20s %25s %25s %25s %25s %25s %25s'%('Centrality','inclusive','0-30','$>$30','inclusive','0-30','$>$30')


    systList=['sell','jer','jec','b','udsg','quench']
    systTitles={'sell':'Nominal','jer':'JER','jec':'JES','b':'b-tag','udsg':'light/unmatched tag','quench':'quenching'}
    for syst in systList:

        print '%20s '%systTitles[syst],

        #loop over first/second jet
        for idx in range(2):
        
            #loop over categories
            for tag in ['PbPb','PbPbcen','PbPbperiph']:        
                beff=beffVals[tag]['sell']
        
                if syst=='sell':
                    print '%25s'%('%3.4f $\\pm$ %3.4f '%(beff[idx][0],beff[idx][1])),
                else:
                    beff_up=beffVals[tag][syst+'up']
                    beff_dn=beffVals[tag][syst+'dn']
                    var_up=abs(beff_up[idx][0]-beff[idx][0])
                    var_dn=abs(beff_dn[idx][0]-beff[idx][0])
                    print '%25s'%(' ~~$\\pm$ %3.4f '%max(var_up,var_dn)),
            
        print ''

if __name__ == "__main__":
    main()
