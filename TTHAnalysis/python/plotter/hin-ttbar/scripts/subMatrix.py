import ROOT, os, datetime, re, operator, math
from array import array
ROOT.gROOT.SetBatch(True)

import utilities
utilities = utilities.util()

## ===================================================================
## USAGE:
## needs as infile a toys.root with limit tree from toys
## takes a comma separated list of regular expressions as input via --params
## if no output directory is given, it will just plot the smaller correlation matrix
## if output directory is given, it will save it there as pdf and png

## example:
## python w-helicity-13TeV/subMatrix.py toys.root --params alph,muR,muF,.*Ybin.*2,pdf12,pdf56,pdf42 --outdir <output_directory> --type toys/hessian
## examples of common regexps:
## NORMs: 'norm_.*'
## PDFs: range [1,20]: '^pdf([1-9]|1[0-9]|20)$'
## ===================================================================
#def SetCorrMatrixPalette():
#ROOT.gStyle.SetPalette()


if __name__ == "__main__":

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPaintTextFormat('.3f')

    date = datetime.date.today().isoformat()

    from optparse import OptionParser
    parser = OptionParser(usage='%prog toys.root [options] ')
    parser.add_option('-o','--outdir', dest='outdir',    default='', type='string', help='outdput directory to save the matrix')
    parser.add_option('-p','--params', dest='params',    default='', type='string', help='parameters for which you want to show the correlation matrix. comma separated list of regexps')
    parser.add_option('-t','--type'  , dest='type'  ,    default='hessian', type='string', help='which type of input file: toys(default),scans, or hessian')
    parser.add_option(     '--suffix', dest='suffix',    default='', type='string', help='suffix for the correlation matrix')
    parser.add_option(     '--parNameCanvas', dest='parNameCanvas',    default='', type='string', help='The canvas name is built using the parameters selected with --params. If they are many, better to pass a name, like QCDscales or PDF for example')
    parser.add_option(     '--nContours', dest='nContours',    default=51, type=int, help='Number of contours in palette. Default is 51 (keep it odd: no correlation is white with our palette)')
    parser.add_option(     '--palette'  , dest='palette',      default=0, type=int, help='Set palette: default is a built-in one, 55 is kRainbow')
    parser.add_option(     '--vertical-labels-X', dest='verticalLabelsX',    default=False, action='store_true', help='Set labels on X axis vertically (sometimes they overlap if rotated)')
    parser.add_option(     '--title'  , dest='title',    default='', type='string', help='Title for matrix ("small correlation matrix" is used as default). Use 0 to remove title')
    parser.add_option(     '--show-more-correlated' , dest='showMoreCorrelated',    default=0, type=int, help='Show the N nuisances more correlated (in absolute value) with the parameters given with --params. If 0, do not do this part')
    (options, args) = parser.parse_args()

    ROOT.TColor.CreateGradientColorTable(3,
                                      array ("d", [0.00, 0.50, 1.00]),
                                      ##array ("d", [1.00, 1.00, 0.00]),
                                      ##array ("d", [0.70, 1.00, 0.34]),
                                      ##array ("d", [0.00, 1.00, 0.82]),
                                      array ("d", [0.00, 1.00, 1.00]),
                                      array ("d", [0.34, 1.00, 0.65]),
                                      array ("d", [0.82, 1.00, 0.00]),
                                      255,  0.95)


    if not options.type in ['toys', 'scans', 'hessian']:
        print 'the given type needs to be either "toys", "scans", or "hessian"!!'
        sys.exit()

    if len(args) < 1:
        print 'You have to pass a root file with the fit result, either for toys or hessian'
        sys.exit()

    if options.outdir:
        ROOT.gROOT.SetBatch()
        if not os.path.isdir(options.outdir):
            os.system('mkdir -p {od}'.format(od=options.outdir))
        os.system('cp {pf} {od}'.format(pf='/afs/cern.ch/user/g/gpetrucc/php/index.php',od=options.outdir))

    pois_regexps = list(options.params.split(','))
    print "Filtering POIs with the following regex: ",pois_regexps

    params = []; indices = []

    ## directly store the mean and RMS into a dictionary
    fitvals = {}; fiterrs = {}

    cov = {}; corr = {}

    ### GET LIST OF PARAMETERS THAT MATCH THE SPECIFIED OPTION IN THE TOYFILE
    if options.type == 'toys':
        toyfile = ROOT.TFile(args[0], 'read')
        _tree = toyfile.Get('fitresults')
        lol = _tree.GetListOfLeaves()

        for l in lol:
            ## skip a bunch of those we don't want
            if '_err'   in l.GetName(): continue
            if '_minos' in l.GetName(): continue
            if '_gen'   in l.GetName(): continue
            if '_In'    in l.GetName(): continue
            for poi in pois_regexps:
                if re.match(poi, l.GetName()):
                    ## draw the parameter into a histogram
                    _tree.Draw(l.GetName()+'>>h_'+l.GetName())
                    ## find that histogram and clone it
                    h = ROOT.gROOT.FindObject('h_'+l.GetName()).Clone()
                    ## store mean and rms into the dictionaries from before
                    fitvals[l.GetName()] = h.GetMean()
                    fiterrs[l.GetName()] = h.GetRMS()
                    ## also keep a list of the parameter names, for sorting
                    params.append(l.GetName())

    elif options.type == 'hessian':
        hessfile = ROOT.TFile(args[0],'read')
        suffix = 'channelmu'
        corrmatrix = hessfile.Get('correlation_matrix_'+suffix)
        covmatrix  = hessfile.Get('covariance_matrix_'+suffix)
        for ib in range(1+corrmatrix.GetNbinsX()+1):
            for poi in pois_regexps:
                if re.match(poi, corrmatrix.GetXaxis().GetBinLabel(ib)):
                    ## store mean and rms into the dictionaries from before
                    ## also keep a list of the parameter names, for sorting
                    params .append(corrmatrix.GetXaxis().GetBinLabel(ib))
                    indices.append(ib)
            
    ## construct the covariances and the correlations in one go.
    for ip1, p1 in enumerate(params):
        for ip2, p2 in enumerate(params):
            if options.type == 'toys':
                var = '({x}-{x0})*({y}-{y0})'.format(x=p1,x0=fitvals[p1],y=p2,y0=fitvals[p2])
                _tree.Draw('{var}>>h_{x}_{y}'.format(var=var,x=p1,y=p2))
                h = ROOT.gROOT.FindObject('h_{x}_{y}'.format(x=p1,y=p2)).Clone()
                cov [(p1,p2)] = h.GetMean()
                corr[(p1,p2)] = cov[(p1,p2)]/(fiterrs[p1]*fiterrs[p2])
            elif options.type == 'hessian':
                cov [(p1,p2)] = covmatrix .GetBinContent(indices[ip1],indices[ip2])
                corr[(p1,p2)] = corrmatrix.GetBinContent(indices[ip1],indices[ip2])
        

    print "===> Build covariance matrix from this set of params: ", params

    p_tmp = set(params)
    params = list(p_tmp)

    # to help sorting with helicity
    # if using more helicity and Y bins, sort by hel,Ybin
    helSorted = { "left" : 1, "right" : 2, "long" : 3}
    chargeSorted = { "Wplus" : 1, "Wminus" : 2}

    ## sort the floatParams. alphabetically, except for pdfs, which are sorted by number
    ## for mu* QCD scales, distinguish among muR and muRXX with XX in 1-10
    
    # why is this commented? Isn't it nice that different charge and polarizations are grouped together?
    # see old example here: 
    # http://mciprian.web.cern.ch/mciprian/wmass/13TeV/helicityAnalysis/electron/fromEmanuele/13-12-18/fitresults_poim1_exp1_bbb1/subMatrix/smallCorrelation_Wplusminus_leftrightYbin_2-5.png
    #params = sorted(params, key= lambda x: (int(chargeSorted[x.split('_')[0]]),int(helSorted[x.split('_')[1]]),int(x.split('_')[-1])) if '_Ybin_' in x else 0)
    # one might want to invert the order of charge and helicity for the sorting

    params = sorted(params, key= lambda x: utilities.getNFromString(x) if '_Ybin_' in x else 0)
    params = sorted(params, key= lambda x: int(x.replace('pdf','')) if 'pdf' in x else 0)
    params = sorted(params, key= lambda x: int(x.replace('muRmuF','')) if ('muRmuF' in x and x != "muRmuF")  else 0)
    params = sorted(params, key= lambda x: int(x.replace('muR','')) if (''.join([j for j in x if not j.isdigit()]) == 'muR' and x != "muR") else 0)
    params = sorted(params, key= lambda x: int(x.replace('muF','')) if (''.join([j for j in x if not j.isdigit()]) == 'muF' and x != "muF") else 0)
    params = sorted(params, key= lambda x: utilities.getNFromString(x) if 'EffStat' in x else 0)            
    params = sorted(params, key= lambda x: utilities.getNFromString(x) if 'FakesEtaUncorrelated' in x else 0)            
    params = sorted(params, key= lambda x: utilities.getNFromString(x) if 'FakesPtUncorrelated' in x else 0)            

    # sort by charge if needed     
    # I think it is useful that different charges are separated in the matrix, it is easier to read it 
    params = sorted(params, key = lambda x: 0 if "plus" in x else 1 if "minus" in x else 2)
    print "sorted params = ", params

    c = ROOT.TCanvas("c","",1200,800)
    c.SetGridx()
    c.SetGridy()
    #ROOT.gStyle.SetPalette(55)
    #ROOT.gStyle.SetNumberContours(200); # default is 20 (values on palette go from -1 to 1)
    if options.nContours: ROOT.gStyle.SetNumberContours(options.nContours)
    if options.palette:   ROOT.gStyle.SetPalette(options.palette)


    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.11)
    c.SetBottomMargin(0.15)

    ## make the new, smaller TH2F correlation matrix
    nbins = len(params)
    th2_sub = ROOT.TH2F('sub_corr_matrix', 'small correlation matrix', nbins, 0., nbins, nbins, 0., nbins)
    th2_sub.GetXaxis().SetTickLength(0.)
    th2_sub.GetYaxis().SetTickLength(0.)
    
    ## pretty nested loop. enumerate the tuples
    for i,x in enumerate(params):
        for j,y in enumerate(params):
            ## set it into the new sub-matrix
            th2_sub.SetBinContent(i+1, j+1, corr[(x,y)])
            ## set the labels correctly
            th2_sub.GetXaxis().SetBinLabel(i+1, x)
            th2_sub.GetYaxis().SetBinLabel(j+1, y)

    th2_sub.GetZaxis().SetRangeUser(-1, 1)
    ROOT.gStyle.SetPaintTextFormat('1.2f')
    if len(params)<30: th2_sub.Draw('colz text45')
    else: th2_sub.Draw('colz')

    if options.verticalLabelsX: th2_sub.LabelsOption("v","X")
    if nbins >= 20: th2_sub.LabelsOption("v","X")

    if options.title: 
        if options.title == "0":
            th2_sub.SetTitle("")
        else:
            th2_sub.SetTitle(options.title)


    if options.parNameCanvas: 
        paramsName = options.parNameCanvas
    else : 
        paramsName = options.params.replace(',','AND')
        for x in ['.', '*', '$', '^', '|', '[', ']', '(', ')']:
            paramsName = paramsName.replace(x,'')
        
    if options.outdir:
        for i in ['pdf', 'png']:
            suff = '' if not options.suffix else '_'+options.suffix
            c.SaveAs(options.outdir+'/smallCorrelation{suff}_{pn}.{i}'.format(suff=suff,i=i,pn=paramsName))
        os.system('cp {pf} {od}'.format(pf='/afs/cern.ch/user/g/gpetrucc/php/index.php',od=options.outdir))


    if options.showMoreCorrelated:
        pass
