import ROOT
import sys

url=sys.argv[1]

fIn=ROOT.TFile(url)
t=fIn.Get('tree')

for hdef,cut in [('abs(bjet_flavorB[0])==0 ? 0 : 1+(abs(bjet_flavorB[0])==5)',
                  'nbjet>0 && bjet_pt[0]>30.'),
                 ('abs(bjet_flavorB[0])==0 ? 0 : 1+(abs(bjet_flavorB[0])==5)',
                  'nbjet>0 && bjet_pt[0]>30. && bjet_csvv2[0]>0.81'),
                 ('abs(bjet_flavorB[1])==0 ? 0 : 1+(abs(bjet_flavorB[1])==5)',
                  'nbjet>1 && bjet_pt[0]>30. && bjet_pt[1]>30.'),
                 ('abs(bjet_flavorB[1])==0 ? 0 : 1+(abs(bjet_flavorB[1])==5)',
                  'nbjet>1 && bjet_pt[0]>30. && bjet_pt[1]>30. && bjet_csvv2[1]>0.81'),
             ]:

    t.Draw("%s>>histo(3,0,3)"%hdef,"ncollWgt*(%s)"%cut,'goff')
    histo=t.GetHistogram()
    histo.Scale(1./histo.Integral())
    histo.Draw('hist')
    print [histo.GetBinContent(i) for i in [1,34,67]]
    histo.Reset('ICE')
