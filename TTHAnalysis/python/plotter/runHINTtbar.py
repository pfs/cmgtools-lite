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


def runplots(trees, friends, targetdir, fmca, fcut, fplots, enabledcuts, disabledcuts, processes, scaleprocesses, fitdataprocess, plotlist, showratio, extraopts = '', invertedcuts = [], submitit = False, name = '', newlumi=0.):
    
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

    if not submitit:
        print '\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        print 'running: python', cmd
        subprocess.call(['python']+cmd.split())#+['/dev/null'],stderr=subprocess.PIPE)
        print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
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

    sf = 'ncollWgt*trigSF[0]*lepSF[0]*lepSF[1]*lepIsoSF[0]*lepIsoSF[1]'

    for flav in ['em', 'ee', 'mm', 'sf']:
        targetdir = basedir+'/combinatorialBackground/{date}{pf}-{flav}/'.format(date=date, pf=('-'+postfix if postfix else ''), flav=flav )

        enable    = [flav]
        disable   = []
        processes = []
        fittodata = []
        scalethem = {}

        extraopts = ' --maxRatioRange 0. 2. --fixRatioRange --plotmode=nostack --showRatio --ratioNums W,data_comb,ttbarSS,data_mixed_comb --ratioDen data_comb '#--preFitData bdt '
        extraopts += ' -W {sfs} '.format(sfs=sf)

        makeplots = ['bdtrarity', 'calsphericity', 'l1pt', 'l2pt', 'calllpt', 'dphi']
        showratio = True
        runplots(trees, friends, targetdir, fmca, fcut, fplots, enable, disable, processes, scalethem, fittodata, makeplots, showratio, extraopts)

def doBFindingControlPlots():
    print '=========================================='
    print 'plotting b-finding control variables'
    print '=========================================='
    trees     = treedir
    friends   = ''

    fmca   = 'hin-ttbar/analysisSetup/mca_bfinding.txt'
    fplots = 'hin-ttbar/analysisSetup/plots_bfinding.txt'
    fcut   = 'hin-ttbar/analysisSetup/cuts.txt'
    targetdir = basedir+'/bfinding/{date}/'.format(date=date)

    enable    = []
    disable   = []
    processes = ['ttbar_PbPb','ttbar_PbPb_cen','ttbar_PbPb_periph']
    fittodata = []
    scalethem = {}

    extraopts = '--maxRatioRange 0. 2. --fixRatioRange --plotmode=norm --showRatio --ratioNums %s --ratioDen %s'%(','.join(processes[1:]),processes[0])
    makeplots = []
    showratio = True
    runplots(trees, friends, targetdir, fmca, fcut, fplots, enable, disable, processes, scalethem, fittodata, makeplots, showratio, extraopts, newlumi=fulllumi)


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

    sf = 'ncollWgt*trigSF[0]*lepSF[0]*lepSF[1]*lepIsoSF[0]*lepIsoSF[1]'

    for flav,mass in itertools.product(['flavmm', 'flavee', 'flavem'],['onZ','offZ']):
        targetdir = basedir+'/jetPlots/{date}{pf}/{flav}_{mass}/'.format(date=date, pf=('-'+postfix if postfix else ''), flav=flav, mass=mass )
        enable    = [flav, mass]

    
        extraopts = ' --maxRatioRange 0. 2. --fixRatioRange --legendColumns 2 --showIndivSigs ' #--plotmode=norm '#--preFitData bdt '

        extraopts += ' -W {sfs} '.format(sfs=sf)

        makeplots = jetvars
        showratio = True
        thislumi  = fulllumi if 'onZ' in mass and not 'em' in flav else 0.
        if 'mm' in flav and 'onZ' in mass:
            thislumi = thislumi-lumidiffmu
        if (replot): runplots(trees, friends, targetdir, fmca, fcut, fplots, enable, disable, processes, scalethem, fittodata, makeplots, showratio, extraopts, newlumi=thislumi)

    yields = {}

    zscaling = {} 

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
            if tmp_denominator:
                tmp_scale_ratio = tmp_numerator/tmp_denominator
            else:
                tmp_scale_ratio = -1.
            tmp_zconstructed.Scale(tmp_scale_ratio)

            ## save the mm spectrum of the on-Z data (subtracted)
            if flav == 'mm':
                em_zconstructed = copy.deepcopy(tmp_zconstructed.Clone('em_zconstructed'))

            #marc
            tmp_zconstructed.Scale(1./tmp_zconstructed.Integral())
            tmp_zg          .Scale(1./tmp_zg          .Integral())

            print '========================================'
            print ' the overall scaling is', tmp_scale_ratio
            print 'THIS IS THE SCALING PER BIN IN', var, 'FOR ', flav
            for ibin in range(1,tmp_zconstructed.GetXaxis().GetNbins()+1):
                t_binc = int(tmp_zconstructed.GetXaxis().GetBinCenter(ibin))
                t_scale = tmp_zconstructed.GetBinContent(ibin)/tmp_zg.GetBinContent(ibin) if tmp_zg.GetBinContent(ibin) else -1.
                print 'bincenter {a:.1f}: scaling is {b:.3f}/{c:.3f} = {d:.3f}'.format(a=t_binc,
                                                                                       b=tmp_zconstructed.GetBinContent(ibin), 
                                                                                       c=tmp_zg.GetBinContent(ibin), 
                                                                                       d=t_scale)
                if ibin <= 3:
                    zscaling[flav+'_'+str(t_binc)+'b'] = float('{a:.3f}'.format(a=t_scale))

        else:
            tmp_zconstructed = em_zconstructed.Clone('zg_constructed'); tmp_zconstructed.SetTitle('Z/#gamma (data #mu#mu)')

            ## get the numerator and denominator for the scaling
            tmp_denominator = tmp_zconstructed.Integral()
            tmp_numerator   = tmp_zg.Integral() + tmp_onZ_zg.Integral()

            ## scale it
            if tmp_denominator:
                tmp_scale_ratio = tmp_numerator/tmp_denominator
            else:
                tmp_scale_ratio = -1.
            tmp_zconstructed.Scale(tmp_scale_ratio)

            #marc
            tmp_zconstructed.Scale(1./tmp_zconstructed.Integral())
            tmp_zg          .Scale(1./tmp_zg          .Integral())


            print '========================================'
            print ' the overall scaling is', tmp_scale_ratio
            print 'THIS IS THE SCALING PER BIN IN', var, 'FOR ', flav
            for ibin in range(1,tmp_zconstructed.GetXaxis().GetNbins()+1):
                t_binc = int(tmp_zconstructed.GetXaxis().GetBinCenter(ibin))
                t_scale = tmp_zconstructed.GetBinContent(ibin)/tmp_zg.GetBinContent(ibin) if tmp_zg.GetBinContent(ibin) else -1.
                print 'bincenter {a:.1f}: scaling is {b:.3f}/{c:.3f} = {d:.3f}'.format(a=t_binc, 
                                                                                       b=tmp_zconstructed.GetBinContent(ibin), 
                                                                                       c=tmp_zg.GetBinContent(ibin), 
                                                                                       d=t_scale)
                if ibin <= 3:
                    zscaling[flav+'_'+str(t_binc)+'b'] = float('{a:.3f}'.format(a=t_scale))

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
        tmp_zconstructed.SetFillColor(ROOT.kAzure+6)
        tmp_zconstructed.SetFillStyle(1001)
        tmp_sig         .SetFillColor(633)
        tmp_sig         .SetLineColor(633)
        tmp_tw          .SetFillColor(ROOT.kTeal+9)
        tmp_vv          .SetFillColor(ROOT.kAzure-5)
        tmp_comb        .SetFillColor(ROOT.kOrange)

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
        lat.DrawLatex(0.68, 0.97, '{l:.1f} #mub^{{-1}} (5.02 TeV)'.format(l=lumi))

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

    print 'this is zscaling', zscaling, '. writing it to a file'
    t_f = open('zjetscaling.data', 'w')
    for k,v in zscaling.items():
        t_f.write('{k} = {v} \n'.format(k=k, v=v)) #zscaling = '+ str(zscaling)+' \n')
    t_f.close()

    print 'done'
    


def compareSignals():
    print '=========================================='
    print 'running comparison of pp vs embedded signals'
    print '=========================================='
    trees     = treedir
    friends   = ''

    fmca          = 'hin-ttbar/analysisSetup/mca_signals.txt'
    fcut          = 'hin-ttbar/analysisSetup/cuts.txt'
    fplots        = 'hin-ttbar/analysisSetup/plots.txt'

    sf = 'ncollWgt*trigSF[0]*lepSF[0]*lepSF[1]*lepIsoSF[0]*lepIsoSF[1]'

    flavs=['mm']
    for iflav,flav in enumerate(flavs):
        targetdir = basedir+'/signalComparison/{date}{pf}-{flav}/'.format(date=date, pf=('-'+postfix if postfix else ''), flav=flav )
        enable    = [flav]
        disable   = []
        processes = ['ttbar', 'ttPowheg', 'ttnCTEQ', 'ttpdf1']
        fittodata = []
        scalethem = {}

        extraopts = ' --maxRatioRange 0. 2. --fixRatioRange --legendColumns 2  --plotmode=nostack --ratioDen ttbar --ratioNums ttbar,ttPowheg,ttnCTEQ,ttpdf1 '
        extraopts += ' -W {sfs} '.format(sfs=sf)
        
        makeplots = ['bdtrarity', 'sphericity']
        showratio = True
        runplots(trees, friends, targetdir, fmca, fcut, fplots, enable, disable, processes, scalethem, fittodata, makeplots, showratio, extraopts)


    #a comparison of the nominal shape with different systematics
    systList=[('pttop','ttptup,ttptdn'),
              ('mtop','ttmup,ttmdn'),
              ('qcdscale','ttmufup,ttmufdn,ttmurup,ttmufdn,ttmurufup,ttmurufdn'),
              #('alphaS', 'ttbar_asup,ttbar_asdn')
              ('pdf','ttpdf1,ttpdf42,ttnCTEQ'),
              ]
    for flav,syst in itertools.product(flavs,systList):
        sname,num=syst
        targetdir = basedir+'/tt_systs/{date}{pf}-{flav}-{sname}/'.format(date=date, pf=('-'+postfix if postfix else ''), flav=flav,sname=sname )
        enable    = [flav]
        disable   = []
        processes = ['ttbar']+num.split(',')
        fittodata = []
        scalethem = {}

        extraopts = ' --maxRatioRange 0.9 1.1 --fixRatioRange --legendColumns 2 --showIndivSigs --plotmode=nostack --ratioDen ttbar --ratioNums %s'%num
        extraopts += ' -W {sfs} '.format(sfs=sf)

        makeplots = ['llpt','sphericity','bdt','bdtrarity']
        showratio = True
        runplots(trees, friends, targetdir, fmca, fcut, fplots, enable, disable, processes, scalethem, fittodata, makeplots, showratio, extraopts)


def compareDataPeriods():
    print '=========================================='
    print 'comparing data periods'
    print '=========================================='
    trees     = treedir
    friends   = ''

    fmca          = 'hin-ttbar/analysisSetup/mca_data.txt'
    fcut          = 'hin-ttbar/analysisSetup/cuts.txt'
    fplots        = 'hin-ttbar/analysisSetup/plots.txt'

    flavs=['ee','mm']
    for iflav,flav in enumerate(flavs):
        targetdir = basedir+'/data/{date}{pf}-{flav}/'.format(date=date, pf=('-'+postfix if postfix else ''), flav=flav )
        enable    = ['flav'+flav,'onZ']
        disable   = [flav]
        processes = ['data','dt_hicen','dt_locen']
        fittodata = []
        scalethem = {}

        extraopts = ' --maxRatioRange 0.5 1.5 --fixRatioRange --legendColumns 2 --showIndivSigs --plotmode=norm --ratioDen data --ratioNums dt_hicen,dt_locen'
        makeplots = []
        showratio = True
        runplots(trees, friends, targetdir, fmca, fcut, fplots, enable, disable, processes, scalethem, fittodata, makeplots, showratio, extraopts)


def dyReweighting():
    print '=========================================='
    print 'running DY Z-pT plots'
    print '=========================================='
    trees     = treedir
    friends   = ''

    fmca          = 'hin-ttbar/analysisSetup/mca_Z.txt'
    fcut          = 'hin-ttbar/analysisSetup/cuts.txt'
    fplots        = 'hin-ttbar/analysisSetup/plots.txt'
    fsysts        = 'hin-ttbar/analysisSetup/systs.txt'

    sf = 'ncollWgt*trigSF[0]*lepSF[0]*lepSF[1]*lepIsoSF[0]*lepIsoSF[1]'


    for (iflav,flav),centrality in itertools.product(enumerate(['mm', 'ee']), ['inclusive', 'centralityLo', 'centralityHi']):
        targetdir = basedir+'/zptReweighting/{date}{pf}-{c}-{flav}/'.format(c=centrality,date=date, pf=('-'+postfix if postfix else ''), flav=flav )

        enable    = ['flav'+flav, 'onZ']
        if not 'inclusive' in centrality:
            enable.append(centrality)

        disable   = []
        processes = ['data', 'zg_raw', 'zg_ptZ_On', 'zg_ptZ_Up', 'zg_ptZ_Dn']
        fittodata = []
        scalethem = {}

        extraopts = ' --maxRatioRange 0.5 1.5 --fixRatioRange --legendColumns 2 --plotmode=norm --ratioDen data --ratioNums data,zg_raw,zg_ptZ_On,zg_ptZ_Up,zg_ptZ_Dn '

        extraopts += ' -W {sf} '.format(sf=sf)

        makeplots = ['dyllpt']#, 'dyllptraw', 'dyllptnocal', 'dyleppt', 'dyl1pt', 'dyl2pt', 'dysphericity', 'dydphi', 'dyllm']
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

    sf = 'ncollWgt*trigSF[0]*lepSF[0]*lepSF[1]*lepIsoSF[0]*lepIsoSF[1]'

    for iflav,flav in enumerate(['ee', 'mm']):
        targetdir = basedir+'/dy_plots/{date}{pf}-{flav}/'.format(date=date, pf=('-'+postfix if postfix else ''), flav=flav )

        enable    = ['flav'+flav, 'onZ']
        disable   = []
        processes = []
        fittodata = []
        scalethem = {}

        extraopts = ' --maxRatioRange 0. 2. --fixRatioRange --legendColumns 2 --showIndivSigs ' #--plotmode=norm '#--preFitData bdt '

        extraopts += ' -W {sf} '.format(sf=sf)

        makeplots = ['dyacoplanarity', 'dyllpt', 'dyleppt', 'dyl1pt', 'dyl2pt', 'dysphericity', 'dydphi', 'dyllm', 'dyllmcal', 'dyl1ptcal', 'dyl2ptcal', 'dyllmnocal', 'dyl1eta', 'dyl2eta', 'dylepeta']
        showratio = True

        unblinded = fulllumi - (0. if 'ee' in flav else lumidiffmu)
        #unblinded = unblinded*1e-9

        runplots(trees, friends, targetdir, fmca, fcut, fplots, enable, disable, processes, scalethem, fittodata, makeplots, showratio, extraopts, newlumi=unblinded)

def makeJetTable(inputdir):
    

    tableDir = {}

    for _file in os.listdir(inputdir):
        if not 'card.txt' in _file: continue

        region = _file.split('.')[0]
        tableDir[region] = {}

        f = open(inputdir+'/'+_file,'r')
        for line in f.readlines():
            if 'observation' in line: 
                tableDir[region]['obs'] = int(float(line.split()[1]))
            if 'process' in line and 'ttbar' in line:
                tableDir[region]['proc'] = line.split()[1:]
            if line.startswith('rate'):
                tableDir[region]['vals'] = line.split()[1:]

    line_obs = 'data ' 
    line_tt  = '\\ttbar ' 
    line_zg  = 'DY ' 
    line_tw  = 'tW ' 
    line_vv  = 'VV ' 
    line_co  = 'comb. ' 

    for ch in ['ee' ,'mm', 'em']:
        for nb in ['0b' ,'1b', '2b']:
            line_obs  += ' & \\textbf{{{o}}} '.format(o=tableDir[ch+'_'+nb]['obs'])
            
            ## ttbar
            tt_ind = tableDir[ch+'_'+nb]['proc'].index('ttbar') if 'ttbar' in tableDir[ch+'_'+nb]['proc'] else -1
            tt_val = float(tableDir[ch+'_'+nb]['vals'][tt_ind]) if tt_ind > -1 else 0.
            line_tt   += ' & {v:.1f} $\\pm$ {e:.1f} '.format(v=tt_val,e=0.3*tt_val)

            ## dy
            zg_ind = tableDir[ch+'_'+nb]['proc'].index('zg') if 'zg' in tableDir[ch+'_'+nb]['proc'] else -1
            zg_val = float(tableDir[ch+'_'+nb]['vals'][zg_ind]) if zg_ind > -1 else 0.
            line_zg   += ' & {v:.1f} $\\pm$ {e:.1f} '.format(v=zg_val,e=0.3*zg_val)

            ## vv
            vv_ind = tableDir[ch+'_'+nb]['proc'].index('VV') if 'VV' in tableDir[ch+'_'+nb]['proc'] else -1
            vv_val = float(tableDir[ch+'_'+nb]['vals'][vv_ind]) if vv_ind > -1 else 0.
            line_vv   += ' & {v:.1f} $\\pm$ {e:.1f} '.format(v=vv_val,e=0.3*vv_val)

            ## tw
            tw_ind = tableDir[ch+'_'+nb]['proc'].index('tW') if 'tW' in tableDir[ch+'_'+nb]['proc'] else -1
            tw_val = float(tableDir[ch+'_'+nb]['vals'][tw_ind]) if tw_ind > -1 else 0.
            line_tw   += ' & {v:.1f} $\\pm$ {e:.1f} '.format(v=tw_val,e=0.3*tw_val)

            ## comb
            co_ind = tableDir[ch+'_'+nb]['proc'].index('data_comb') if 'data_comb' in tableDir[ch+'_'+nb]['proc'] else -1
            co_val = float(tableDir[ch+'_'+nb]['vals'][co_ind]) if co_ind > -1 else 0.
            line_co   += ' & {v:.1f} $\\pm$ {e:.1f} '.format(v=co_val,e=0.3*co_val)

    line_obs += '\\\\ \\hline \\hline'
    line_tt  += '\\\\ \\hline'
    line_zg  += '\\\\'
    line_tw  += '\\\\' 
    line_vv  += '\\\\' 
    line_co  += '\\\\' 

    tabletex = '''
\\begin{{table}}[!htb]
\\centering
\\small
\\topcaption{{
    The number of expected background and signal events and the observed event yields in the different
    channels of $ee$, $\\mu\\mu$, and $e\\mu$, prior to the fit.
    \\label{{tab:yieldsJet}}}}
    \\begin{{tabular}}{{lccc|ccc|ccc}}
Process        & $ee$ 0b    & $ee$ 1b   & $ee$ 2b   & $\\mu\\mu$ 0b   & $\\mu\\mu$ 1b  & $\\mu\\mu$ 2b  & $e\\mu$ 0b    & $e\\mu$  1b    & $e\\mu$  2b   \\\\ \\hline \\hline
                                                                                                                                                                         
{line_co}
{line_zg}
{line_tw}
{line_vv}
{line_tt}
{line_obs}
\\end{{tabular}}

\\end{{table}}'''.format(line_obs=line_obs,line_tt=line_tt,line_zg=line_zg,line_co=line_co,line_tw=line_tw,line_vv=line_vv)

    print tabletex

##Combinatorial  & 6.89 $\pm$ 1.38   & 0.83 $\pm$ 0.48  & 0.83 $\pm$ 0.48  & 5.51 $\pm$ 1.23   & 0.28 $\pm$ 0.28  & 0.28 $\pm$ 0.28  & 16.81 $\pm$ 2.15 & 0.00 $\pm$ 0.00   & 0.00 $\pm$ 0.00  \\
##Z/$\gamma^{*}$ & 199.04 $\pm$ 2.59 & 6.70 $\pm$ 0.48  & 6.70 $\pm$ 0.48  & 443.12 $\pm$ 3.78 & 14.26 $\pm$ 0.68 & 14.26 $\pm$ 0.68 & 21.91 $\pm$ 0.19 & 0.71 $\pm$ 0.03   & 0.71 $\pm$ 0.03  \\
##$\cPqt\PW$     & 0.35 $\pm$ 0.00   & 0.33 $\pm$ 0.00  & 0.33 $\pm$ 0.00  & 0.83 $\pm$ 0.01   & 0.64 $\pm$ 0.01  & 0.64 $\pm$ 0.01  & 1.47 $\pm$ 0.01  & 1.25 $\pm$ 0.01   & 1.25 $\pm$ 0.01  \\
##VV             & 0.75 $\pm$ 0.00   & 0.00 $\pm$ 0.00  & 0.00 $\pm$ 0.00  & 1.53 $\pm$ 0.01   & 0.00 $\pm$ 0.00  & 0.00 $\pm$ 0.00  & 2.33 $\pm$ 0.01  & 0.01 $\pm$ 0.00   & 0.01 $\pm$ 0.00  \\
##\ttbar signal  & 1.42 $\pm$ 0.03   & 2.78 $\pm$ 0.04  & 2.78 $\pm$ 0.04  & 4.27 $\pm$ 0.06   & 6.37 $\pm$ 0.08  & 6.37 $\pm$ 0.08  & 6.81 $\pm$ 0.07  & 11.55 $\pm$ 0.10  & 11.55 $\pm$ 0.10 \\ \hline \hline
##Total          & 208.45 $\pm$ 2.94 & 10.64 $\pm$ 0.68 & 10.64 $\pm$ 0.68 & 455.26 $\pm$ 3.98 & 21.56 $\pm$ 0.74 & 21.56 $\pm$ 0.74 & 49.33 $\pm$ 2.16 & 13.51 $\pm$ 0.10  & 13.51 $\pm$ 0.10 \\

    ##print tableDir

        

def makeJetAnalysis():
    print '=========================================='
    print 'making datacards for jet analysis'
    print '=========================================='
    trees     = treedir
    friends   = ''

    fmca   = 'hin-ttbar/analysisSetup/mca_forCards.txt'
    fcut   = 'hin-ttbar/analysisSetup/cuts.txt'
    fplots = 'hin-ttbar/analysisSetup/plots_jetanalysis.txt'
    fsysts = 'hin-ttbar/analysisSetup/systs.txt'

    nbinsForFit = 3
    fitVars = []
    fitVars.append( ('sphericity'  , '\'llpt/(lep_pt[0]+lep_pt[1])\' {n},0.,1.'.format(n=nbinsForFit)) )
    #fitVars.append( ('bdt'         , '\'bdtrarity\'                      {n},0.,1.'.format(n=nbinsForFit)) )
    fitVars.append( ('bdtcomb'       , '\'bdtcdfinv(bdt)\'                 {n},0.,1.'.format(n=nbinsForFit)) )

    sf = 'ncollWgt*trigSF[0]*lepSF[0]*lepSF[1]*lepIsoSF[0]*lepIsoSF[1]'

    regions = ['ee_0b', 'ee_1b', 'ee_2b',
               'mm_0b', 'mm_1b', 'mm_2b',
               'em_0b', 'em_1b', 'em_2b']

    ## with elpt 25 zscaling ={'ee_0b': 0.983, 'ee_1b': 1.160, 'ee_2b': 1.207,
    ## with elpt 25            'mm_0b': 0.979, 'mm_1b': 1.147, 'mm_2b': 1.380,
    ## with elpt 25            'em_0b': 1.101, 'em_1b': 1.548, 'em_2b': 2.196 }
    ## with both pT 20 zscaling ={'ee_0b': 0.975, 'ee_1b': 1.268, 'ee_2b': 1.495,
    ## with both pT 20            'mm_0b': 0.980, 'mm_1b': 1.145, 'mm_2b': 1.312,
    ## with both pT 20            'em_0b': 1.089, 'em_1b': 1.331, 'em_2b': 2.339 }
    ## replaced by the read from the file zscaling ={'ee_0b': 0.967, 'ee_1b': 1.433, 'ee_2b': 1.118,
    ## replaced by the read from the file            'mm_0b': 0.982, 'mm_1b': 1.129, 'mm_2b': 1.375,
    ## replaced by the read from the file            'em_0b': 1.108, 'em_1b': 1.597, 'em_2b': 1.237 }

    ## get the scaling for the mixed combinatorial from the nbjets file
    combscaling = {}
    combscale_f = ROOT.TFile('combinatorialrescaling.root','read')
    h_ssda = combscale_f.Get('nbjets_data_comb')
    h_mix1 = combscale_f.Get('nbjets_data_mixed_comb1')
    h_mix5 = combscale_f.Get('nbjets_data_mixed_comb5')
    h_sstt = combscale_f.Get('nbjets_ttbarSS')
    for ibin in range(1,4):
        combscaling[str(ibin-1)+'b'] = max(0., (h_ssda.GetBinContent(ibin)-h_sstt.GetBinContent(ibin)) / h_mix1.GetBinContent(ibin) )

    print 'this is combscaling', combscaling
    combscale_f.Close()
    

    zscaling = {}
    t_f = open('zjetscaling.data', 'r')
    for line in t_f.readlines():
        zscaling[line.split()[0]] = float(line.split()[-1])


    for iflav,flav in enumerate(regions):

        targetdir = basedir+'/card_inputs_jetanalysis/{date}{pf}/'.format(date=date, pf=('-'+postfix if postfix else ''), flav=flav )

        enable    = flav.split('_')
        disable   = []
        processes = []
        fittodata = []
        scalethem = {}#'zg': zscaling[flav]}#, 'data_comb': combscaling[flav.split('_')[1]]}

        extraopts = ' --maxRatioRange 0. 2. --fixRatioRange --legendColumns 2 --showIndivSigs ' #--plotmode=norm '#--preFitData bdt '

        extraopts += ' -W {sfs} '.format(sfs=sf)

        makeplots = [flav+'_bdtcomb', flav+'_sphericity'] #[flav+'_bdt', flav+'_sphericity', flav+'_bdtcomb']
        showratio = True

        runplots(trees, friends, targetdir, fmca, fcut, fplots, enable, disable, processes, scalethem, fittodata, makeplots, showratio, extraopts, newlumi=fulllumi if 'onZ' in flav else 0.)

        cmd_cards_base = 'python makeShapeCards.py --s2v -f -j 6 -v 3 -l {lumi} -P {tdir} {mca} {cuts} '.format(lumi=lumi if not 'onZ' in flav else fulllumi,tdir=trees,mca=fmca,cuts=fcut)
        for fitVar in fitVars:
            outdirCards = 'hin-ttbar/datacards_{date}_{pf}_jetAnalysis/{fitVarName}/'.format(date=date, pf=postfix,fitVarName=fitVar[0])

            cmd_cards  = cmd_cards_base + ' -W {sfs} '.format(sfs=sf)
            cmd_cards += ' --od {od} '                  .format(od=outdirCards)
            cmd_cards += ' -o {flav} '                  .format(flav=flav)
            cmd_cards += ' '.join([' -E ^'+i+' ' for i in flav.split('_')])
            ## cmd_cards += ' --scale-process zg {f:.3f} '       .format(f=zscaling[flav])

            systfile = fsysts if not '2b' in flav else 'hin-ttbar/analysisSetup/systs2b.txt'

            cmd_cards += ' {fitVar} {systs} '.format(fitVar=fitVar[1].replace(str(nbinsForFit),'1') if '2b' in flav else fitVar[1], systs=systfile)

            print 'running the cards with command'
            print '==========================================='
            print cmd_cards
            os.system(cmd_cards)

    makeJetTable(outdirCards)

            ## if flav == regions[-1]:
            ##     print 'running combine cards'
            ##     cmd_combinecards = 'combineCards.py '
            ##     for r in regions:
            ##         cmd_combinecards += ' {f}={p}/{f}.card.txt '.format(p=outdirCards, f=r)
            ##     f_card_allFlavors = '{p}/allFlavors.card.txt'.format(p=outdirCards)
            ##     cmd_combinecards += ' > '+f_card_allFlavors
            ##     print 'combining cards with this command:'
            ##     print '===================================='
            ##     print cmd_combinecards
            ##     os.system(cmd_combinecards)
            ##     
            ##     print 'now putting the groups and stuff in the datacards'
            ##     os.system('cp {f} {f}.tmp'.format(f=f_card_allFlavors))

            ##     with open(f_card_allFlavors, 'a') as tmp_file:
            ##         tmp_file.write('theory      group = alphaS muR muF muRmuF ' + ' '.join('pdf'+str(x) for x in range(1,101)) +'\n')
            ##         tmp_file.write('ptModel     group = ptTop ptZ \n')  
            ##         tmp_file.write('topMass     group = mTop \n')  
            ##         tmp_file.write('luminosity  group = lumi \n')  
            ##         tmp_file.write('bkg         group = VV_lnN data_comb_lnN tW_lnN zg_lnN \n')  

            ##     tmp_file.close()

            ##     print '=========================================='
            ##     print '====== DONE SO FAR. NOW RUN COMBINE ======'
            ##     print ' ... load the right version, and then run:'
            ##     
            ##     print 'text2hdf5.py {p}/allFlavors.card.txt --out {p}/allFlavors.hdf5 '.format(p=outdirCards)
            ##     print 'combinetf.py --binByBinStat --computeHistErrors --saveHists --doImpacts -t -1 {p}/allFlavors.hdf5 '.format(p=outdirCards)
            ##     f_fitresults = '{p}/fitresults_allFlavors.root'.format(p=outdirCards)
            ##     print 'mv fitresults_123456789.root {f}'.format(f=f_fitresults)
            ##     
            ##     print '=========================================='
            ##     print '====== ONCE DONE WITH FITTING, MAKE ======'
            ##     print ' ... some plots for the postfit stuff'
            ##     print 'python hin-ttbar/scripts/diffNuisances.py --infile {inf} --pois "pdf.*,muR,muF,muRmuF,alphaS,ptTop,ptZ,mTop" --outdir {pd} -a --format html > nuisances_theory.html'.format(inf=f_fitresults, pd=targetdir)
            ##     print 'python hin-ttbar/scripts/diffNuisances.py --infile {inf} --pois ".*lnN.*,lumi" --outdir {pd} -a --format html > nuisances_experimental.html'.format(inf=f_fitresults, pd=targetdir)
            ##     print 'python hin-ttbar/scripts/subMatrix.py {inf} --params "ttbar_mu,pdf.*" --outdir {pd} '.format(inf=f_fitresults,pd=targetdir)
            ##     print 'python hin-ttbar/scripts/subMatrix.py {inf} --params "ttbar_mu,muR,muF,muRmuF,alphaS,ptTop,ptZ,mTop,.*lnN.*,lumi" --outdir {pd} '.format(inf=f_fitresults,pd=targetdir)

            ##     print '=========================================='
            ##     print 'that should be all for now... need to find a way to also do the numbers automatically...'

def checkSphericity():
    print '=========================================='
    print 'checking sphericity < 0.1 events'
    print '=========================================='
    trees     = treedir
    friends   = ''

    fmca   = 'hin-ttbar/analysisSetup/mca_forCards.txt'
    fcut   = 'hin-ttbar/analysisSetup/cuts.txt'
    fplots = 'hin-ttbar/analysisSetup/plots.txt'
    fsysts = 'hin-ttbar/analysisSetup/systs.txt'

    flavors = ['mm']#'ee', 'mm', 'em', ]

    sf = 'ncollWgt*trigSF[0]*lepSF[0]*lepSF[1]*lepIsoSF[0]*lepIsoSF[1]'


    makeplots = ['l1pt' 'l2pt', 'l1eta', 'l2eta', 'lepeta', 'dphi', 'dphifine', 'acoplanfine']

    for iflav,flav in enumerate(flavors):
        for cen in ['centralityLo', 'centralityHi']:
            targetdir = basedir+'/low_sphericity_events/{date}{pf}-{flav}-{cen}/'.format(date=date, pf=('-'+postfix if postfix else ''), flav=flav, cen=cen )

            enable    = [flav, cen]
            disable   = []
            processes = []
            fittodata = []
            scalethem = {}

            extraopts = ' --maxRatioRange 0. 2. --fixRatioRange --legendColumns 2 --showIndivSigs ' #--plotmode=norm '#--preFitData bdt '

            extraopts += ' -W {sfs} -I coplan0p01 '.format(sfs=sf)

            ## if 'onZ' in flav:
            ##     makeplots = [flav+'mll']
            ## else:
            ##     makeplots = [flav+x[0] for x in fitVars]

            showratio = True

            fs_lumi = lumi
            if 'onZ' in flav and not 'em' in flav:
                fs_lumi = fulllumi
            if 'mm' in flav:
                fs_lumi = fs_lumi - lumidiffmu

            runplots(trees, friends, targetdir, fmca, fcut, fplots, enable, disable, processes, scalethem, fittodata, makeplots, showratio, extraopts, newlumi=fs_lumi)

def makeCards():
    print '=========================================='
    print 'making datacards, including the Z for mm and ee'
    print '=========================================='
    trees     = treedir
    friends   = ''

    fmca   = 'hin-ttbar/analysisSetup/mca_forCards.txt'
    fcut   = 'hin-ttbar/analysisSetup/cuts.txt'
    fplots = 'hin-ttbar/analysisSetup/plots_results.txt'
    fsysts = 'hin-ttbar/analysisSetup/systs.txt'

    nbinsForFit = 10
    fitVars = [## ('bdt'         , 'bdtrarity                      {n},0.,1.'.format(n=nbinsForFit)), 
               ## ('minmlb' , '\'min(mass_2(lep_pt[0],lep_eta[0],lep_phi[0],0.,bjet_pt[0],bjet_eta[0],bjet_phi[0],bjet_mass[0]),min(mass_2(lep_pt[1],lep_eta[1],lep_phi[1],0.,bjet_pt[1],bjet_eta[1],bjet_phi[1],bjet_mass[1]),min(mass_2(lep_pt[0],lep_eta[0],lep_phi[0],0.,bjet_pt[1],bjet_eta[1],bjet_phi[1],bjet_mass[1]),mass_2(lep_pt[1],lep_eta[1],lep_phi[1],0.,bjet_pt[0],bjet_eta[0],bjet_phi[0],bjet_mass[0]))))\' 15,0.,150.'),
               #('sphericity'  , '\'llpt/(lep_pt[0]+lep_pt[1])\' {n},0.,1.'.format(n=nbinsForFit)),
               ('sphericity'  , '\'pt_2(lep_calpt[0],lep_phi[0],lep_calpt[1],lep_phi[1])/(lep_calpt[0]+lep_calpt[1])\' {n},0.,1.'.format(n=nbinsForFit)),
               #('bdtcomb'  , '\'bdtcdfinv(bdt)\' {n},0.,1.'.format(n=nbinsForFit)),
               ('ptll'  , '\'llpt\' 15,0.,150.'),
              ]

    regions = ['ee', 'mm', 'em']#, 'leponZee', 'leponZmm']

    sf = 'ncollWgt*trigSF[0]*lepSF[0]*lepSF[1]*lepIsoSF[0]*lepIsoSF[1]'

    for iflav,flav in enumerate(regions):
        targetdir = basedir+'/card_inputs/{date}{pf}/'.format(date=date, pf=('-'+postfix if postfix else ''), flav=flav )

        enable    = [flav, 'highBDT']
        disable   = []
        processes = []
        fittodata = []
        scalethem = {}

        extraopts = ' --maxRatioRange 0. 2. --fixRatioRange --legendColumns 2 --showIndivSigs ' #--plotmode=norm '#--preFitData bdt '

        extraopts += ' -W {sfs} '.format(sfs=sf)

        if 'onZ' in flav:
            makeplots = [flav+'mll']
        else:
            makeplots = [flav+x[0] for x in fitVars]

        showratio = True

        fs_lumi = lumi
        if 'onZ' in flav:
            fs_lumi = fulllumi
        if 'mm' in flav:
            fs_lumi = lumi - lumidiffmu

        runplots(trees, friends, targetdir, fmca, fcut, fplots, enable, disable, processes, scalethem, fittodata, makeplots, showratio, extraopts, newlumi=fs_lumi)

        cmd_cards_base = 'python makeShapeCards.py --s2v -f -j 6 -v 3 -l {lumi} -P {tdir} {mca} {cuts} '.format(lumi=fs_lumi,tdir=trees,mca=fmca,cuts=fcut)
        for fitVar in fitVars:
            outdirCards = 'hin-ttbar/datacards_{date}_{pf}/{fitVarName}/'.format(date=date, pf=postfix,fitVarName=fitVar[0])

            cmd_cards  = cmd_cards_base + ' -W {sfs} '.format(sfs=sf)
            cmd_cards += ' --od {od} '                  .format(od=outdirCards)
            cmd_cards += ' -E ^{flav} -o {flav} '       .format(flav=flav)
            cmd_cards += ' -E ^highBDT '

            ## use the fitvar for offZ and em, use mll for onZ
            if not 'onZ' in flav:
                cmd_cards += ' {fitVar} {systs} '.format(fitVar=fitVar[1], systs=fsysts)
            else:
                cmd_cards += ' llm  30,76.,106. '

            
            print 'running the cards with command'
            print '==========================================='
            print cmd_cards
            os.system(cmd_cards)

            if flav == regions[-1]:
                print 'running combine cards'
                cmd_combinecards = 'combineCards.py '
                for r in regions:
                    cmd_combinecards += ' {f}={p}/{f}.card.txt '.format(p=outdirCards, f=r)
                f_card_allFlavors = '{p}/combinedCard.txt'.format(p=outdirCards)
                cmd_combinecards += ' > '+f_card_allFlavors
                print 'combining cards with this command:'
                print '===================================='
                print cmd_combinecards
                os.system(cmd_combinecards)
                
                print 'now putting the groups and stuff in the datacards'
                os.system('cp {f} {f}.tmp'.format(f=f_card_allFlavors))

                with open(f_card_allFlavors, 'a') as tmp_file:
                    tmp_file.write('* autoMCStats 0 1 1 \n') # alphaS + ' '.join('pdf'+str(x) for x in range(1,101)) +'\n')
                    tmp_file.write('theory      group = muR muF muRmuF \n') # alphaS + ' '.join('pdf'+str(x) for x in range(1,101)) +'\n')
                    tmp_file.write('ptModel     group = ptTop ptZ \n')  
                    tmp_file.write('pdfs        group = {pdf} \n'.format(pdf=' '.join('pdf'+str(x) for x in range(1,49)) )  )
                    tmp_file.write('topMass     group = mTop \n')  
                    tmp_file.write('luminosity  group = lumi \n')  
                    tmp_file.write('bkg         group = VV_lnN data_comb_lnN tW_lnN zg_lnN \n')  

                tmp_file.close()

                print '=========================================='
                print '====== DONE SO FAR. NOW RUN COMBINE ======'
                print ' ... load the right version, and then run:'
                
                print 'text2hdf5.py {p}/allFlavors.card.txt --out {p}/allFlavors.hdf5 '.format(p=outdirCards)
                print 'combinetf.py --binByBinStat --computeHistErrors --saveHists --doImpacts -t -1 {p}/allFlavors.hdf5 '.format(p=outdirCards)
                f_fitresults = '{p}/fitresults_allFlavors.root'.format(p=outdirCards)
                print 'mv fitresults_123456789.root {f}'.format(f=f_fitresults)
                
                print '=========================================='
                print '====== ONCE DONE WITH FITTING, MAKE ======'
                print ' ... some plots for the postfit stuff'
                print 'python hin-ttbar/scripts/diffNuisances.py --infile {inf} --pois "pdf.*,muR,muF,muRmuF,alphaS,ptTop,ptZ,mTop" --outdir {pd} -a --format html > nuisances_theory.html'.format(inf=f_fitresults, pd=targetdir)
                print 'python hin-ttbar/scripts/diffNuisances.py --infile {inf} --pois ".*lnN.*,lumi" --outdir {pd} -a --format html > nuisances_experimental.html'.format(inf=f_fitresults, pd=targetdir)
                print 'python hin-ttbar/scripts/subMatrix.py {inf} --params "ttbar_mu,pdf.*" --outdir {pd} '.format(inf=f_fitresults,pd=targetdir)
                print 'python hin-ttbar/scripts/subMatrix.py {inf} --params "ttbar_mu,muR,muF,muRmuF,alphaS,ptTop,ptZ,mTop,.*lnN.*,lumi" --outdir {pd} '.format(inf=f_fitresults,pd=targetdir)

                print '=========================================='
                print 'that should be all for now... need to find a way to also do the numbers automatically...'
                
            
                ##f_res = ROOT.TFile(' {p}/fitresults_allFlavors.root'.format(p=outdirCards), 'read')
                ##t_res = f_res.Get('fitresults')
                ##for ev in t_res:
                ##    print '================================================='
                ##    print 'for fitvar', fitVar[0]
                ##    print 'RESULT: mu(ttbar) = {mu:.2f} +- {err:.2f}'.format(mu=ev.ttbar_mu, err=ev.ttbar_mu_err)
                ##    print '================================================='

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

def jetFlavor():
    print '=========================================='
    print 'running some jet flavor plots'
    print '=========================================='
    trees     = treedir
    friends   = ''

    fmca          = 'hin-ttbar/analysisSetup/mca.txt'
    fcut          = 'hin-ttbar/analysisSetup/cuts.txt'
    fplots        = 'hin-ttbar/analysisSetup/plots_results.txt'
    fsysts        = 'hin-ttbar/analysisSetup/systs.txt'

    makeplots = ['j1genpt','j2genpt','j1csv','j2csv'] #'jetptflavor999', 'jetcsvflavor999']#'jet1flavorpt30', 'jet2flavorpt30']#, 'jetpt', 'jetpt1', 'jetpt2', 'jetflavorpt30', 'jetflavor1', 'jetflavor2', 'jeteta', 'jetflavor', 'centrality']#'ptjet1vsptjet2', 'ptvsjetflavor1', 'ptvsjetflavor2', 'ptvscsvjet1', 'ptvscsvjet2']

    sf = 'ncollWgt*trigSF[0]*lepSF[0]*lepSF[1]*lepIsoSF[0]*lepIsoSF[1]'

    for (iflav,flav),centrality in itertools.product(enumerate(['anyflavor']), ['inclusive', 'centralityLo', 'centralityHi']):
        targetdir = basedir+'/jet_studies/{date}{pf}-{c}-{flav}/'.format(date=date, pf=('-'+postfix if postfix else ''), flav=flav, c=centrality )

        enable    = [flav]

        if not 'inclusive' in centrality:
            enable.append(centrality)

        disable   = []
        processes = ['ttbar', 'zg']#, 'tW', 'VV']
        fittodata = []
        scalethem = {}

        #extraopts = ' --legendColumns 2 --plotmode=norm '#--preFitData bdt '
        extraopts = ' --plotmode=norm '#--preFitData bdt '
        extraopts += ' -W {eff} '.format(eff=sf)
        extraopts += ' --lspam \"#bf{CMS}  #it{Simulation Preliminary}\"'
        extraopts += ' --rspam \"(#sqrt{s_{NN}}=5.02 TeV)\"'
        showratio = False
        runplots(trees, friends, targetdir, fmca, fcut, fplots, enable, disable, processes, scalethem, fittodata, makeplots, showratio, extraopts)


def simplePlot():
    print '=========================================='
    print 'running simple plots'
    print '=========================================='
    trees     = treedir
    friends   = ''

    fmca          = 'hin-ttbar/analysisSetup/mca.txt'
    fcut          = 'hin-ttbar/analysisSetup/cuts.txt'
    fplots        = 'hin-ttbar/analysisSetup/plots_results.txt'
    fsysts        = 'hin-ttbar/analysisSetup/systs.txt'

    makeplots = ['mll', 'mllz', 'llpt', 'sphericity', 'acoplan', 'minmlb', 'minmlb2', 'centrality'] #'centrality']#'llpt', 'l1pt', 'l2pt', 'l1eta', 'l2eta', 'acoplan', 'mll', #'mtll',  'l1d0', 'l2d0', 
                 #'l1dz', 'l2dz', 'l1sip2d', 'l2sip2d', 'njets', #'jetpt', 'jeteta', 'jetbtag', 
                 #'jetmass', 'centrality', 'bdt', 'bdtrarity', 'j1btag', 'j2btag']

    sf = 'ncollWgt*trigSF[0]*lepSF[0]*lepSF[1]*lepIsoSF[0]*lepIsoSF[1]'

    for iflav,flav in enumerate(['em','sf','leponZee','leponZmm','anyflavor']) : #em', 'sf','leponZll']): # 'ee', 'mm']):
        targetdir = basedir+'/simple_plots/{date}{pf}-{flav}/'.format(date=date, pf=('-'+postfix if postfix else ''), flav=flav )

        enable    = [flav]
        if flav=='anyflavor' :
            enable += 'highBDT'
        disable   = []
        processes = []
        fittodata = []
        scalethem = {}

        extraopts = '--maxRatioRange 0. 2. --fixRatioRange --legendColumns 2 --showIndivSigs --noStatTotLegendOnRatio' #--plotmode=norm '#--preFitData bdt '
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

## LUMI=1618.466*(1e-6)
## CHLUMI={'mm':1587.941*(1e-6),
##         'ee':1664.148*(1e-6),
##         'blind':446.931*(1e-6)}


    global date, postfix, lumi, date, basedir, treedir
    postfix = opts.postfix
    #lumi = 1618.466*1e-9 if not opts.lumi else opts.lumi
    blindedlumi = 457.99*1e-9
    fulllumi = 1751.83*1e-9
    lumidiffmu = 32.5*1e-9
    lumi = fulllumi if not opts.lumi else float(opts.lumi)
    date = datetime.date.today().isoformat()

    user = os.environ['USER']
    if user == 'mdunser':
        basedir = '/afs/cern.ch/user/m/mdunser/www/private/heavyIons/plots/'
    elif user == 'psilva':
        basedir = '/eos/user/p/psilva/www/HIN-19-001'
    
    ## this is pp MC treedir = '/eos/cms/store/cmst3/group/hintt/PbPb2018_skim27Apr/'
    ## this is with old eleID and stuff treedir = '/eos/cms/store/cmst3/group/hintt/PbPb2018_skim21June/' ## rereco data and mixed, official MC
    #treedir = '/eos/cms/store/cmst3/group/hintt/PbPb2018_skim13August/' ## newest on 07/08/2019
    #treedir = '/eos/cms/store/cmst3/group/hintt/PbPb2018_skim02Sep_loose/' ## newest on 03/09/2019
    #treedir = '/eos/cms/store/cmst3/group/hintt/PbPb2018_skimsFinal/' ## newest. unblinded
    treedir = '/eos/cms/store/cmst3/group/hintt/PbPb2018_skim09Oct_loose/' ## newest. unblinded

    if opts.date:
        date = opts.date

    if opts.simple:
        print 'making simple plots'
        simplePlot()
    if opts.combinatorial:
        print 'making comb comparison plots'
        compareCombBackgrounds()
    if opts.jetvars:
        print 'plotting jet related variables'
        plotJetVariables(opts.replot)
        #doBFindingControlPlots()
    if opts.compareSignals:
        print 'plotting jet related variables'
        compareSignals()
    if opts.dyPlots:
        print 'plotting jet related variables'
        makeZplots()
    if opts.dyReweighting:
        print 'plotting some control plots for dy pT reweighting'
        dyReweighting()
    if opts.compareData:
        print 'Comparing data periods'
        compareDataPeriods()
    if opts.fit:
        print 'Making input for the datacards'
        makeCards()
    if opts.jetAnalysis:
        print 'making the cards for the jet based analysis'
        makeJetAnalysis()
    if opts.checkSphericity:
        print 'chekcing low sphericity events'
        checkSphericity()
    if opts.jetFlavor:
        print 'checking low sphericity events'
        jetFlavor()
