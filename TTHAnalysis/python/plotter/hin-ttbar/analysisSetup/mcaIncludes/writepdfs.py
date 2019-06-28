import os

f = open('mca_signal_theory.txt','r')

lines = f.readlines()

f.close()

outf = open('mca_signal_theory_withpdfs.txt','w')

for l in lines:
    outf.write(l)

temp = 'ttbar_pdfXXX_ZZZ : TTJets_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8 : 69*(208**2) : lep_charge[0]*lep_charge[1]<0 ; AddWeight="AAAmeWeights[YYY]", Label="t\#bar{t} pdfXXX ZZZ", SkipMe=True\n'

for i in range(0,48):
    pdfind = 10+1+2*i
    outf.write(temp.replace('XXX',str(i+1)).replace('YYY',str(pdfind  )).replace('ZZZ','Up').replace('AAA', ''))
    outf.write(temp.replace('XXX',str(i+1)).replace('YYY',str(pdfind+1)).replace('ZZZ','Dn').replace('AAA', ''))
    outf.write('\n')


outf.close()

