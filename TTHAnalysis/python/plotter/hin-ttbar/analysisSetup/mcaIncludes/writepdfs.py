import os

f = open('mca_signal_theory.txt','r')

lines = f.readlines()

f.close()

outf = open('mca_signal_theory_withpdfs.txt','w')

for l in lines:
    outf.write(l)

temp = 'ttbar_pdfXXX_ZZZ : TT_TuneCP5_5p02TeV-powheg-pythia8 : 69*(208**2) : lep_charge[0]*lep_charge[1]<0 ; AddWeight="AAAmeWeights[YYY]", Label="t\#bar{t} pdfXXX ZZZ", SkipMe=True\n'

for i in range(1,101):
    pdfind = 11+i
    outf.write(temp.replace('XXX',str(i)).replace('YYY',str(pdfind)).replace('ZZZ','Up').replace('AAA', ''   ))
    outf.write(temp.replace('XXX',str(i)).replace('YYY',str(pdfind)).replace('ZZZ','Dn').replace('AAA', '1./'))
    outf.write('\n')


outf.close()

