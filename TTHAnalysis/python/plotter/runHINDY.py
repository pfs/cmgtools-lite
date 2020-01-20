import optparse, subprocess, ROOT, datetime, math, array, copy, os, itertools
import numpy as np


def runplots(trees, friends, targetdir, fmca, fcut, fplots, enabledcuts, disabledcuts, processes, scaleprocesses, fitdataprocess, plotlist, showratio, extraopts = '', invertedcuts = [], name = '', newlumi=0.):
    
    if not type(trees)==list: trees = [trees]
    treestring = ' '.join(' -P '+ t for t in list(trees))
    cmd  = ' mcPlots.py --s2v -f -j 6 -l {lumi} --pdir {td} {trees} {fmca} {fcut} {fplots}'.format(lumi=lumi if not newlumi else newlumi, td=targetdir, trees=treestring, fmca=fmca, fcut=fcut, fplots=fplots)
    if friends:
        if not type(friends) == list: friends = [friends]
        #cmd += ' -F Friends {friends}/tree_Friend_{{cname}}.root'.format(friends=friends)
        for fff in friends:
            cmd += ' -F Friends {friends}/tree_Friend_{{cname}}.root '.format(friends=fff)
            
    cmd += ''.join(' -E ^'+cut for cut in enabledcuts )
    cmd += ''.join(' -X ^'+cut for cut in disabledcuts)
    if len(plotlist):
        cmd += ' --sP '+','.join(plot for plot in plotlist)
        cmd += ' -o '+targetdir+'/'+'_AND_'.join(plot for plot in plotlist)+'.root'
    else:
        cmd += ' -o '+targetdir+'/ALL.root'
    if len(processes):
        cmd += ' -p '+','.join(processes)
    if invertedcuts:
        cmd += ''.join(' -I ^'+cut for cut in invertedcuts )
    #if doPUandSF and not '-W ' in extraopts: cmd += ' -W puw2016_nTrueInt_36fb(nTrueInt)*LepGood_effSF[0] '
    if fitdataprocess:
        cmd+= ' --fitData '
        cmd+= ''.join(' --flp '+proc for proc in fitdataprocess)
    if scaleprocesses:
        for proc,scale in scaleprocesses.items():
            cmd += ' --scale-process {proc} {scale} '.format(proc=proc, scale=scale)
    showrat   = ''
    if showratio:
        showrat = ' --showRatio '
    cmd += showrat
    if extraopts:
        cmd += ' '+extraopts

    print '\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    print 'running: python', cmd
    subprocess.call(['python']+cmd.split())#+['/dev/null'],stderr=subprocess.PIPE)
    print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'


def simplePlot():
    print '=========================================='
    print 'running simple plots'
    print '=========================================='
    trees     = treedir
    friends   = ''

    fmca          = 'hin-dy/analysisSetup/mca.txt'
    fcut          = 'hin-dy/analysisSetup/cuts.txt'
    fplots        = 'hin-dy/analysisSetup/plots.txt'

    makeplots  = ['mll','llpt','acoplan','sphericity','st','maxleta','minleta','centrality']
    makeplots += ['l1pt','l1eta','l1sip2d','l1d0','l1dz','l1iso02']
    makeplots += ['l2pt','l2eta','l2sip2d','l2d0','l2dz','l2iso02']

    sf = 'ncollWgt*trigSF[0]*lepSF[0]*lepSF[1]'

    for iflav,flav in enumerate(['em', 'ee', 'mm']):
        targetdir = basedir+'/simple_plots/{date}{pf}-{flav}/'.format(date=date, pf=('-'+postfix if postfix else ''), flav=flav )

        enable    = [flav]
        disable   = []
        processes = []
        fittodata = []
        scalethem = {}

        extraopts = ' --maxRatioRange 0. 2. --fixRatioRange --legendColumns 2 --showIndivSigs '
        extraopts += ' -W {eff} '.format(eff=sf)
        showratio = True
        runplots(trees, friends, targetdir, fmca, fcut, fplots, enable, disable, processes, scalethem, fittodata, makeplots, showratio, extraopts)


if __name__ == '__main__':
    parser = optparse.OptionParser(usage='usage: %prog [opts] ', version='%prog 1.0')
    parser.add_option('--pf'        , '--postfix'    , dest='postfix'      , type='string'       , default=''    , help='postfix for running each module')
    parser.add_option('-d'          , '--date'       , dest='date'         , type='string'       , default=''    , help='run with specified date instead of today')
    parser.add_option('-l'          , '--lumi'       , dest='lumi'         , type='float'        , default=0.    , help='change lumi by hand')
    parser.add_option('--simple'    ,                  dest='simple'       , action='store_true' , default=False , help='make simple plot')
    parser.add_option('--combinatorial'    , action='store_true' , default=False , help='compare data comb with wjets')
    parser.add_option('--jetvars', action='store_true' , default=False , help='plot jet variables')
    parser.add_option('--replot' , action='store_true' , default=False , help='replot all the input jet plots')
    parser.add_option('--fit'    , action='store_true' , default=False , help='perform the fits to data with combine')
    parser.add_option('--compareSignals'    , action='store_true' , default=False , help='compareSignals')
    parser.add_option('--compareData'    , action='store_true' , default=False , help='compare data periods')
    parser.add_option('--dyPlots'    , action='store_true' , default=False , help='make plots for onZ ee/mm')
    parser.add_option('--jetAnalysis'    , action='store_true' , default=False , help='make the jet analysis datacards')
    parser.add_option('--dyReweighting'    , action='store_true' , default=False , help='make Z-pT reweighting test plots')
    parser.add_option('--checkSphericity'    , action='store_true' , default=False , help='check low sphericity events')
    parser.add_option('--jetFlavor'    , action='store_true' , default=False , help='plots for jet flavor studies')
    (opts, args) = parser.parse_args()


    global date, postfix, lumi, date, basedir, treedir
    postfix = opts.postfix
    fulllumi = 1751.83*1e-9
    lumidiffmu = 32.5*1e-9
    lumi = fulllumi if not opts.lumi else float(opts.lumi)
    date = datetime.date.today().isoformat()
    treedir='/eos/cms/store/cmst3/group/hintt/PbPb2018_200120_loose'

    user = os.environ['USER']
    if user == 'mdunser':
        basedir = '/afs/cern.ch/user/m/mdunser/www/private/heavyIons/plots/'
    elif user == 'psilva':
        basedir = '/eos/user/p/psilva/www/HIN-DY'

    if opts.date:
        date = opts.date

    if opts.simple:
        print 'making simple plots'
        simplePlot()
