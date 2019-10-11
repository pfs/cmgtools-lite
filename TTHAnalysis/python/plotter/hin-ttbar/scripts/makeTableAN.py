import ROOT, os, sys


## usage: 
##        python makeTableAN.py ~/www/private/heavyIons/plots/card_inputs/2019-06-26-mixedWithSFs/ sphericity

def translate(name):
    if 'nonpro' in name: return 'comb'
    if 'gamma' in name: return 'zg'
    else: return name

indir = sys.argv[1]
var   = sys.argv[2]

tabletex = '''\\begin{table}[!htb]
\\centering
\\topcaption{ The number of expected background and signal events and the observed event yields in the different
channels of $ee$, $\\mu\\mu$, and $e\\mu$, prior to the fit.
\\label{tab:yields}}
    \\begin{tabular}{lccc}
        Process        & $ee$  & $\\mu\\mu$ & $e\\mu$ \\\\ \\hline \\hline\n'''
        

val = {}

for chan in ['ee', 'mm', 'em']:

    tmp_f = open(indir+'/'+chan+var+'.txt', 'r')
    lines = tmp_f.readlines()
    for line in lines:
        if '-------------' in line: continue
        if not 'DATA' in line:
            name = ''.join(line.split()[:-4])
            goodname = translate(name)
            val[chan+goodname ] = '{v:10s} $\\pm$ {e:.2f}'.format(v=line.split()[-4], e=float(line.split()[-4])*0.3) #line.split()[-2])
        else:
            val[chan+ 'DATA'] = line.split()[-1]

for k,v in  val.items():
    print k, v
    
tabletex += '''
Nonprompt    & {eecomb} & {mmcomb} & {emcomb} \\\\
Z/$\\gamma^{{*}}$   & {eezg} & {mmzg} & {emzg} \\\\
$\\cPqt\\PW$       & {eetW} & {mmtW} & {emtW} \\\\
VV               & {eeVV} & {mmVV} & {emVV} \\\\ \\hline \\hline
Total bkg        & {eeBACKGROUND} & {mmBACKGROUND} & {emBACKGROUND} \\\\
\\ttbar signal    & {eeSIGNAL} & {mmSIGNAL} & {emSIGNAL} \\\\ \\hline \\hline
Observed (data)    & {eeDATA}                     & {mmDATA}                     & {emDATA}                      \\\\ \\hline \\hline'''.format(**val)

tabletex += '''

\\end{tabular}
\\end{table}'''

print tabletex
