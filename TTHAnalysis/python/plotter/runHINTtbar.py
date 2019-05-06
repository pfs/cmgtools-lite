import optparse, subprocess, ROOT, datetime, math, array, copy, os, itertools
import numpy as np



def getFitresult(url):
    res=[]
    fIn=ROOT.TFile.Open(url)
    t=fIn.Get('limit')
    for i in range(3):
        t.GetEntry(i)
        res.append( t.r )
        if i==0 : continue
        res[-1]=res[-1]-res[0]
    fIn.Close()
    return res

#doPUreweighting = True
doPUandSF = False

def submitIt(ODIR, name, cmd, dryRun=True):
    outdir=ODIR+"/jobs/"
    if not os.path.isdir(outdir): 
        os.system('mkdir -p '+outdir)
    srcfile = outdir+name+".sh"
    srcfile_op = open(srcfile,"w")
    srcfile_op.write("#! /bin/sh\n")
    srcfile_op.write("ulimit -c 0\n")
    srcfile_op.write("cd {cmssw};\neval $(scramv1 runtime -sh);\ncd {d};\n".format( 
            d = os.getcwd(), cmssw = os.environ['CMSSW_BASE']))
    srcfile_op.write(cmd+'\n')
    srcfile_op.close()
    os.system("chmod a+x "+srcfile)

    condorf_name = srcfile.replace('.sh','.condor')
    condorf_op = open(condorf_name,'w')
    condorf_op.write('''Universe = vanilla
Executable = {scriptName}
use_x509userproxy = $ENV(X509_USER_PROXY)
Log        = {od}{name}.log
Output     = {od}{name}.out
Error      = {od}{name}.error
getenv      = True

request_memory = 4000
+MaxRuntime = 14400
+AccountingGroup = "group_u_CMST3.all"
queue 1\n
 '''.format(scriptName=os.path.abspath(srcfile), name=name, od=outdir ) )
    condorf_op.close()
    subcmd = 'condor_submit {cf} '.format(cf=condorf_name)
    if dryRun: 
        print "[DRY-RUN]: ", subcmd
    else: os.system(subcmd)

def printAggressive(s):
    print '='.join('' for i in range(len(s)+1))
    print s
    print '='.join('' for i in range(len(s)+1))

def readScaleFactor(path, process, reterr = False):
    infile = open(path,'r')
    lines = infile.readlines()
    
    for line in lines:
        if 'Process {proc} scaled by'.format(proc=process) in line:
            scale = float(line.split()[4])
            scaleerr = float(line.split()[-1])
    if not reterr:
        return scale
    else:
        return scale, scaleerr

def runefficiencies(trees, friends, targetdir, fmca, fcut, ftight, fxvar, enabledcuts, disabledcuts, scaleprocesses, compareprocesses, showratio, extraopts = ''):
    
    if not type(trees)==list: trees = [trees]
    treestring = ' '.join(' -P '+ t for t in list(trees))
    cmd  = ' mcEfficiencies.py --s2v -f -j 6 -l {lumi} -o {td} {trees} {fmca} {fcut} {ftight} {fxvar}'.format(lumi=lumi, td=targetdir, trees=treestring, fmca=fmca, fcut=fcut, ftight=ftight, fxvar=fxvar)
    if friends:
        #cmd += ' --Fs {friends}'.format(friends=friends)
        #cmd += ' -F mjvars/t {friends}/friends_evVarFriend_{{cname}}.root --FMC sf/t {friends}/friends_sfFriend_{{cname}}.root  '.format(friends=friends)
        cmd += ' -F Friends {friends}/tree_Friend_{{cname}}.root'.format(friends=friends)
    # not needed here cmd += ' --mcc ttH-multilepton/mcc-eleIdEmu2.txt --mcc dps-ww/mcc-tauveto.txt '
    ## cmd += ' --obj treeProducerWMassEle ' ## the tree is called 'treeProducerWMassEle' not 'tree'
    cmd += ' --groupBy cut '
    #if doPUandSF and not '-W ' in extraopts: cmd += ' -W puw2016_nTrueInt_36fb(nTrueInt)*LepGood_effSF[0] '
    cmd += ''.join(' -E ^'+cut for cut in enabledcuts )
    cmd += ''.join(' -X ^'+cut for cut in disabledcuts)
    cmd += ' --compare {procs}'.format(procs=(','.join(compareprocesses)  ))
    if scaleprocesses:
        for proc,scale in scaleprocesses.items():
            cmd += ' --scale-process {proc} {scale} '.format(proc=proc, scale=scale)
    showrat   = ''
    if showratio:
        showrat = ' --showRatio '
    cmd += showrat
    if extraopts:
        cmd += ' '+extraopts

    print 'running: python', cmd
    subprocess.call(['python']+cmd.split())#+['/dev/null'],stderr=subprocess.PIPE)


def runplots(trees, friends, targetdir, fmca, fcut, fplots, enabledcuts, disabledcuts, processes, scaleprocesses, fitdataprocess, plotlist, showratio, extraopts = '', invertedcuts = [], submitit = False, name = ''):
    
    if not type(trees)==list: trees = [trees]
    treestring = ' '.join(' -P '+ t for t in list(trees))
    cmd  = ' mcPlots.py --s2v -f -j 6 -l {lumi} --pdir {td} {trees} {fmca} {fcut} {fplots}'.format(lumi=lumi, td=targetdir, trees=treestring, fmca=fmca, fcut=fcut, fplots=fplots)
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

    if not submitit:
        print 'running: python', cmd
        subprocess.call(['python']+cmd.split())#+['/dev/null'],stderr=subprocess.PIPE)
    elif submitit == 'return':
        return 'python '+cmd
    else:
        submitIt(targetdir, name if name else ''.join(plotlist), 'python '+cmd, False)

def compareCombBackgrounds():
    print '=========================================='
    print 'comparing wjets with data combinatorial'
    print '=========================================='
    trees     = treedir
    friends   = ''

    fmca          = 'hin-ttbar/checkCombinatorial/mca.txt'
    fcut          = 'hin-ttbar/analysisSetup/cuts.txt'
    fplots        = 'hin-ttbar/analysisSetup/plots.txt'

    for flav in ['ee']:
        targetdir = basedir+'/combinatorialBackground/{date}{pf}-{flav}/'.format(date=date, pf=('-'+postfix if postfix else ''), flav=flav )

        enable    = [flav]
        disable   = []
        processes = []
        fittodata = []
        scalethem = {}

        extraopts = ' --maxRatioRange 0. 2. --fixRatioRange --plotmode=norm --showRatio --ratioNums W,data_comb,ttbar --ratioDen data_comb '#--preFitData bdt '
        makeplots = []
        showratio = True
        runplots(trees, friends, targetdir, fmca, fcut, fplots, enable, disable, processes, scalethem, fittodata, makeplots, showratio, extraopts)

def plotJetVariables(replot):
    print '=========================================='
    print 'plotting jet variables with DY from data'
    print '=========================================='
    trees     = treedir
    friends   = ''

    fmca          = 'hin-ttbar/analysisSetup/mca.txt'
    fcut          = 'hin-ttbar/analysisSetup/cuts.txt'
    fplots        = 'hin-ttbar/analysisSetup/plots.txt'
    fsysts        = 'hin-ttbar/analysisSetup/systs.txt'

    disable   = []
    processes = []
    fittodata = []
    scalethem = {}

    jetvars = ['nbjets']#, 'njets']

    eff_e = 0.75
    eff_m = 0.90

    for flav,mass in itertools.product(['flavmm', 'flavee', 'flavem'],['onZ','offZ']):
        targetdir = basedir+'/jetPlots/{date}{pf}/{flav}_{mass}/'.format(date=date, pf=('-'+postfix if postfix else ''), flav=flav, mass=mass )
        enable    = [flav, mass]

    
        extraopts = ' --maxRatioRange 0. 2. --fixRatioRange --legendColumns 2 --showIndivSigs ' #--plotmode=norm '#--preFitData bdt '

        ## reweight by hand to expected efficiency for isolation!
        effscale  = eff_m**2 if 'mm' in flav else eff_e*eff_m if 'em' in flav else eff_e**2
        extraopts += ' -W {eff} '.format(eff=effscale)

        makeplots = jetvars
        showratio = True
        if (replot): runplots(trees, friends, targetdir, fmca, fcut, fplots, enable, disable, processes, scalethem, fittodata, makeplots, showratio, extraopts)

    yields = {}

    for flav,var in itertools.product(['mm', 'ee', 'em'],jetvars):
        targetdir = basedir+'/jetPlots/{date}{pf}/'.format(date=date, pf=('-'+postfix if postfix else ''))

        fname_base = '_AND_'.join(jetvars)

        ## first get the easy ones from the offZ file (this is enough for ee and mm)
        ## ====================

        f_nominal = ROOT.TFile(targetdir+'/flav{flav}_offZ/{fname}.root'.format(flav=flav,fname=fname_base),'read')
    
        tmp_data  = f_nominal.Get(var+'_data'     ); tmp_data.SetTitle('data'       )
        tmp_sig   = f_nominal.Get(var+'_ttbar'    ); tmp_sig .SetTitle('t#bar{t}'   )
        tmp_tw    = f_nominal.Get(var+'_tW'       ); tmp_tw  .SetTitle('tW'         )
        tmp_comb  = f_nominal.Get(var+'_data_comb'); tmp_comb.SetTitle('comb (data)')
        tmp_vv    = f_nominal.Get(var+'_VV'       ); tmp_vv  .SetTitle('VV'         )
        tmp_zg    = f_nominal.Get(var+'_zg'       )

        ## now construct the Z histogram
        ## first get onZ data
        f_onZ = ROOT.TFile(targetdir+'/flav{flav}_onZ/{fname}.root'.format(flav=flav,fname=fname_base),'read')
        tmp_zconstructed    = f_onZ.Get(var+'_data'     ); tmp_zconstructed.SetName('zg_constructed'); tmp_zconstructed.SetTitle('Z/#gamma (data)')
        tmp_nonZ_tt         = f_onZ.Get(var+'_ttbar'    )
        tmp_nonZ_tw         = f_onZ.Get(var+'_tW'       )
        tmp_nonZ_comb       = f_onZ.Get(var+'_data_comb')
        tmp_nonZ_vv         = f_onZ.Get(var+'_VV'       )
        tmp_onZ_zg          = f_onZ.Get(var+'_zg'       )

        ## for em add all the non-Z backgrounds back to the off-Z ones
        if (flav == 'em'):
            tmp_sig   .Add(tmp_nonZ_tt)
            tmp_tw    .Add(tmp_nonZ_tw)
            tmp_comb  .Add(tmp_nonZ_comb)
            tmp_vv    .Add(tmp_nonZ_vv)
            tmp_data  .Add(tmp_zconstructed) ## can do, because not yet scaled!!

        ## ATTENTION!!! FIXME WHAT TO DO WITH THE Z HERE...


        if not flav == 'em':
            ## subtract everything that isn't Z (mostly the SS i guess)
            tmp_zconstructed.Add(tmp_nonZ_tt  , -1.)
            tmp_zconstructed.Add(tmp_nonZ_tw  , -1.)
            tmp_zconstructed.Add(tmp_nonZ_comb, -1.)
            tmp_zconstructed.Add(tmp_nonZ_vv  , -1.)

            ## get the numerator and denominator for the scaling
            tmp_denominator = tmp_zconstructed.Integral()
            tmp_numerator   = tmp_zg          .Integral()

            ## scale it
            tmp_zconstructed.Scale(tmp_numerator/tmp_denominator)

            ## save the mm spectrum of the on-Z data (subtracted)
            if flav == 'mm':
                em_zconstructed = copy.deepcopy(tmp_zconstructed.Clone('em_zconstructed'))

        else:
            tmp_zconstructed = em_zconstructed.Clone('zg_constructed'); tmp_zconstructed.SetTitle('Z/#gamma (data #mu#mu)')

            ## get the numerator and denominator for the scaling
            tmp_denominator = tmp_zconstructed.Integral()
            tmp_numerator   = tmp_zg.Integral() + tmp_onZ_zg.Integral()

            ## scale it
            tmp_zconstructed.Scale(tmp_numerator/tmp_denominator)

        ## putting all backgrounds into a list for sorting
        backgrounds = []
        backgrounds.append(tmp_comb)
        backgrounds.append(tmp_vv)
        backgrounds.append(tmp_tw)
        backgrounds.append(tmp_zconstructed)
        ## don't sort... backgrounds = sorted(backgrounds, key = lambda x: x.Integral())#, reverse = True if flav == 'em' else False)
        backgrounds.append(tmp_sig)


        ## now on to the plotting
        ## have to reset the fill colors?
        tmp_zconstructed.SetFillColor(ROOT.kOrange+6)
        tmp_zconstructed.SetFillStyle(1001)
        tmp_sig         .SetFillColor(633)
        tmp_sig         .SetLineColor(633)
        tmp_tw          .SetFillColor(ROOT.kTeal+9)
        tmp_vv          .SetFillColor(ROOT.kAzure-5)
        tmp_comb        .SetFillColor(17)

        ## more plotting. this is really boring to code...
        leg = ROOT.TLegend(0.5,0.80,0.95,0.95)
        leg.SetNColumns(2); leg.SetLineWidth(0); leg.SetFillStyle(0)

        tmp_total = backgrounds[0].Clone('total_bkg'); tmp_total.SetMarkerSize(0.); tmp_total.SetTitle('')
        tmp_stack = ROOT.THStack()
        for ib,bkg in enumerate(backgrounds):
            tmp_stack.Add(bkg)
            if ib: tmp_total.Add(bkg)
            leg.AddEntry(bkg, bkg.GetTitle(), 'f')
        leg.AddEntry(tmp_data, 'Data', 'pe')

        if var == 'nbjets':
            for bkg in backgrounds+[tmp_data, tmp_total]:
                err = ROOT.Double(-1.)
                yields[bkg.GetTitle()+' '+flav+' 0b'  ] = bkg.GetBinContent(1)
                yields[bkg.GetTitle()+' '+flav+' 0b_e'] = bkg.GetBinError  (1)
                yields[bkg.GetTitle()+' '+flav+' 1b'  ] = bkg.IntegralAndError(2, bkg.GetXaxis().GetNbins()+1, err)
                yields[bkg.GetTitle()+' '+flav+' 1b_e'] = err
            
        print yields

        tmp_sigcopy = tmp_sig.Clone('sigcopy')
        tmp_sigcopy.SetFillStyle(0)
        tmp_sigcopy.SetLineWidth(4)

        ROOT.gROOT.SetBatch(); ROOT.gStyle.SetOptStat(0)
        ##tmp_canv = f_nominal.Get(var+'_canvas')
        tmp_canv= ROOT.TCanvas('whatever'+var+flav,'',600,750)
        tmp_canv.GetPad(0).SetTopMargin   (0.04);
        tmp_canv.GetPad(0).SetBottomMargin(0.32);
        tmp_canv.GetPad(0).SetLeftMargin  (0.18);
        tmp_canv.GetPad(0).SetRightMargin (0.04);
        tmp_canv.GetPad(0).SetTickx(1);
        tmp_canv.GetPad(0).SetTicky(1);

        ## now to the main PAD
        
        tmp_canv.cd(1)
        tmp_canv.SetLogy()
        tmp_total.SetMaximum(2.*max(tmp_total.GetMaximum(),tmp_data.GetMaximum()))
        tmp_total.SetMinimum(0.9)
        tmp_total.GetXaxis().SetLabelSize(0)
        tmp_total.Draw()
        
        tmp_stack.Draw('hist same')
        tmp_sigcopy.Draw('hist same')
        tmp_data.Draw('same')
        leg.Draw('same')
        tmp_total.Draw('AXIS same')

        lat = ROOT.TLatex()
        lat.SetNDC(); lat.SetTextFont(42); lat.SetTextSize(0.03)
        lat.DrawLatex(0.18, 0.97, '#font[61]{CMS} #font[52]{Preliminary}')
        lat.DrawLatex(0.68, 0.97, '446.9 #mub^{-1} (5.02 TeV)')

        ## FIRST!!! to the ratio:
        tmp_ratio = tmp_data.Clone('ratio'); tmp_ratio.SetTitle('')
        tmp_ratio.Divide(tmp_total)

        tmp_ratio.GetYaxis().SetRangeUser(0.,2.)
        tmp_ratio.GetYaxis().SetTitle('Data/pred.')
        tmp_ratio.GetYaxis().SetTitleSize(0.16)
        tmp_ratio.GetYaxis().SetTitleOffset(0.5)
        tmp_ratio.GetYaxis().SetLabelSize(0.10)
        tmp_ratio.GetYaxis().SetNdivisions(505)
        tmp_ratio.GetXaxis().SetLabelSize(0.10)
        tmp_ratio.GetXaxis().SetTitleSize(0.15)
        tmp_ratio.GetXaxis().SetTitleOffset(1.1)

        ratiopad = ROOT.TPad('padratio', '', 0.,0.,1.,0.28)
        ratiopad.SetTopMargin   (0.00);
        ratiopad.SetBottomMargin(0.32);
        ratiopad.SetLeftMargin  (0.18);
        ratiopad.SetRightMargin (0.04);
        ratiopad.SetTickx(1);
        ratiopad.SetTicky(1);

        ratiopad.Draw()
        ratiopad.cd()
        tmp_ratio.Draw('pe')
        line = ROOT.TLine()
        line.SetLineWidth(2)
        line.DrawLine(tmp_ratio.GetXaxis().GetXmin(),1.,tmp_ratio.GetXaxis().GetXmax(),1.)
        tmp_ratio.Draw('pe same')

        for k,v in sorted(yields.items()):
            if '_e' in k: continue
            print '{k} {c:.2f} $\pm$ {e:.2f}'.format(k=k,c=v, e=yields[k+'_e'])


        for ext in ['pdf','png','root']:
            tmp_canv.SaveAs(targetdir+'/{f}_{v}.{e}'.format(v=var,e=ext,f=flav))

        with open(targetdir+'/{f}_{v}.txt'.format(v=var,f=flav), 'w') as f:
            for b in backgrounds:
                f.write('{n}: {integral:.2f} \n'.format(n=b.GetName(),integral=b.Integral()))
            f.write('\ndata: {integral:.0f} \n'.format(integral=tmp_data.Integral()))

        os.system(' cp ~mdunser/public/index.php '+targetdir)


def compareSignals():
    print '=========================================='
    print 'running comparison of pp vs embedded signals'
    print '=========================================='
    trees     = treedir
    friends   = ''

    fmca          = 'hin-ttbar/analysisSetup/mca_signals.txt'
    fcut          = 'hin-ttbar/analysisSetup/cuts.txt'
    fplots        = 'hin-ttbar/analysisSetup/plots.txt'

    for iflav,flav in enumerate(['em', 'ee', 'mm']):
        targetdir = basedir+'/signalComparison/{date}{pf}-{flav}/'.format(date=date, pf=('-'+postfix if postfix else ''), flav=flav )
        enable    = [flav]
        disable   = []
        processes = ['ttbar', 'ttEmbedded', 'bjetQCD']
        fittodata = []
        scalethem = {}

        extraopts = ' --maxRatioRange 0. 2. --fixRatioRange --legendColumns 2 --showIndivSigs --plotmode=norm --ratioDen ttbar --ratioNums ttbar,ttEmbedded '
        makeplots = []
        showratio = True
        runplots(trees, friends, targetdir, fmca, fcut, fplots, enable, disable, processes, scalethem, fittodata, makeplots, showratio, extraopts)

def makeZplots():
    print '=========================================='
    print 'running DY plots'
    print '=========================================='
    trees     = treedir
    friends   = ''

    fmca          = 'hin-ttbar/analysisSetup/mca.txt'
    fcut          = 'hin-ttbar/analysisSetup/cuts.txt'
    fplots        = 'hin-ttbar/analysisSetup/plots.txt'
    fsysts        = 'hin-ttbar/analysisSetup/systs.txt'

    eff_e = 0.75
    eff_m = 0.90

    for iflav,flav in enumerate(['ee', 'mm']):
        targetdir = basedir+'/dy_plots/{date}{pf}-{flav}/'.format(date=date, pf=('-'+postfix if postfix else ''), flav=flav )

        enable    = ['flav'+flav, 'onZ']
        disable   = [flav]
        processes = []
        fittodata = []
        scalethem = {}

        extraopts = ' --maxRatioRange 0. 2. --fixRatioRange --legendColumns 2 --showIndivSigs ' #--plotmode=norm '#--preFitData bdt '

        effscale  = eff_m**2 if flav == 'mm' else eff_e*eff_m if flav == 'em' else eff_e**2
        extraopts += ' -W {eff} '.format(eff=effscale)

        makeplots = ['dyllpt', 'dyleppt', 'dyl1pt', 'dyl2pt', 'dysphericity', 'dydphi', 'dyllm']
        showratio = True
        runplots(trees, friends, targetdir, fmca, fcut, fplots, enable, disable, processes, scalethem, fittodata, makeplots, showratio, extraopts)


def simplePlot(makeCards):
    print '=========================================='
    print 'running simple plots'
    print '=========================================='
    trees     = treedir
    friends   = ''

    fmca          = 'hin-ttbar/analysisSetup/mca.txt'
    fmca_forCards = 'hin-ttbar/analysisSetup/mca_forCards.txt'
    fcut          = 'hin-ttbar/analysisSetup/cuts.txt'
    fplots        = 'hin-ttbar/analysisSetup/plots.txt'
    fsysts        = 'hin-ttbar/analysisSetup/systs.txt'

    nbinsForFit = 10
    fitVars = [('bdt'         , 'bdtrarity                      {n},0.,1.'.format(n=nbinsForFit)), 
               ('acoplanarity', '\'(1.-abs(dphi)/TMath::Pi())\' {n},0.,1.'.format(n=nbinsForFit)),
               ('sphericity'  , '\'llpt/(lep_pt[0]+lep_pt[1])\' {n},0.,1.'.format(n=nbinsForFit)),
               ('distance5'   , '\'getAvgDistance(llpt,abs(dphi),abs(lleta),abs(lep_eta[0]+lep_eta[1]),lep_pt[0],5)\' {n},1.,4.'.format(n=nbinsForFit)),
              ]

    eff_e = 0.75
    eff_m = 0.90

    for iflav,flav in enumerate(['em']):#, 'ee', 'mm']):
        targetdir = basedir+'/simple_plots/{date}{pf}-{flav}/'.format(date=date, pf=('-'+postfix if postfix else ''), flav=flav )

        enable    = [flav]
        disable   = []
        processes = []
        fittodata = []
        scalethem = {}


        extraopts = ' --maxRatioRange 0. 2. --fixRatioRange --legendColumns 2 --showIndivSigs ' #--plotmode=norm '#--preFitData bdt '
        effscale  = eff_m**2 if flav == 'mm' else eff_e*eff_m if flav == 'em' else eff_e**2
        extraopts += ' -W {eff} '.format(eff=effscale)
        makeplots = []
        showratio = True
        runplots(trees, friends, targetdir, fmca, fcut, fplots, enable, disable, processes, scalethem, fittodata, makeplots, showratio, extraopts)

        if makeCards:

            for fitVar in fitVars:
                outdirCards = 'hin-ttbar/datacards_{date}_{pf}/{fitVarName}/'.format(date=date, pf=postfix,fitVarName=fitVar[0])
                ## os.system('python makeShapeCardsSusy.py --s2v -f -j 6 -l {lumi} --od {od} -P {tdir} {mca} {cuts} -E ^{flav} {fitVar} {systs} -v 3 -o {flav} -b {flav} -W {scale} '.format(lumi=lumi,tdir=trees,mca=fmca_forCards,cuts=fcut,systs=fsysts,flav=flav,fitVar=fitVar[1],od=outdirCards, scale=effscale))

                if iflav == 2:
                    print 'running combine cards'
                    os.system('combineCards.py em={p}/ttbar/em.card.txt mm={p}/ttbar/mm.card.txt ee={p}/ttbar/ee.card.txt > {p}/ttbar/allFlavors.card.txt'           .format(p=outdirCards))
                    os.system('text2hdf5.py {p}/ttbar/allFlavors.card.txt --out {p}/ttbar/allFlavors.hdf5 '.format(p=outdirCards))
                    os.system('combinetf.py --binByBinStat --computeHistErrors --saveHists --doImpacts -t -1 {p}/ttbar/allFlavors.hdf5 '.format(p=outdirCards))
                    os.system('mv fitresults_123456789.root {p}/ttbar/fitresults_allFlavors.root'.format(p=outdirCards))
                
                    f_res = ROOT.TFile(' {p}/ttbar/fitresults_allFlavors.root'.format(p=outdirCards), 'read')
                    t_res = f_res.Get('fitresults')
                    for ev in t_res:
                        print '================================================='
                        print 'for fitvar', fitVar[0]
                        print 'RESULT: mu(ttbar) = {mu:.2f} +- {err:.2f}'.format(mu=ev.ttbar_mu, err=ev.ttbar_mu_err)
                        print '================================================='

                ## print 'running combine with systs'
                ## os.system('combine -M MultiDimFit {p}/ttbar/allFlavors.card.txt -t -1 --expectSignal=1 --saveFitResult --robustFit=1 --algo=cross --cl=0.68'     .format(p=outdirCards))
                ## resTot=getFitresult('higgsCombineTest.MultiDimFit.mH120.root')

                ## print 'running combine without systs'
                ## os.system('combine -M MultiDimFit {p}/ttbar/allFlavors.card.txt -t -1 --expectSignal=1 --saveFitResult --robustFit=1 --algo=cross --cl=0.68 -S 0'.format(p=outdirCards))
                ## resStat=getFitresult('higgsCombineTest.MultiDimFit.mH120.root')

                ## systHi=math.sqrt(resTot[1]**2-resStat[1]**2)
                ## systLo=math.sqrt(resTot[2]**2-resStat[2]**2)

                ## print '%3.3f +%3.3f-%3.3f (syst) +%3.3f-%3.3f (stat)'%(resTot[0],systHi,systLo,resStat[1],resStat[2])
                ## print resTot
                ## print resStat

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
    parser.add_option('--dyPlots'    , action='store_true' , default=False , help='make plots for onZ ee/mm')
    (opts, args) = parser.parse_args()

## LUMI=1618.466*(1e-6)
## CHLUMI={'mm':1587.941*(1e-6),
##         'ee':1664.148*(1e-6),
##         'blind':446.931*(1e-6)}


    global date, postfix, lumi, date, basedir, treedir
    postfix = opts.postfix
    #lumi = 1618.466*1e-9 if not opts.lumi else opts.lumi
    lumi = 446.931*1e-9 if not opts.lumi else float(opts.lumi)
    date = datetime.date.today().isoformat()

    user = os.environ['USER']
    if user == 'mdunser':
        basedir = '/afs/cern.ch/user/m/mdunser/www/private/heavyIons/plots/'
    elif user == 'psilva':
        basedir = 'foobar'
    
    treedir = '/eos/cms/store/cmst3/group/hintt/PbPb2018_skim27Apr/'

    if opts.date:
        date = opts.date

    if opts.simple:
        print 'making simple plots'
        simplePlot(opts.fit)
    if opts.combinatorial:
        print 'making comb comparison plots'
        compareCombBackgrounds()
    if opts.jetvars:
        print 'plotting jet related variables'
        plotJetVariables(opts.replot)
    if opts.compareSignals:
        print 'plotting jet related variables'
        compareSignals()
    if opts.dyPlots:
        print 'plotting jet related variables'
        makeZplots()
