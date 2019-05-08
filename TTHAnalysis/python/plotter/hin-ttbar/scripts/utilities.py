
import ROOT, os, sys, re, array, math, json
import numpy as np

class util:

    def __init__(self):
        self.cmssw_version = os.environ['CMSSW_VERSION']
        self.isRecentRelease = (len(self.cmssw_version) and int(self.cmssw_version.split('_')[1]) > 8)

    def solvePol2(self, a,b,c):
    
        # calculate the discriminant
        d = (b**2) - (4*a*c)
    
        if not a or d < 0:
            return (0,0,0)
        
        # find two solutions
        sol1 = (-b-math.sqrt(d))/(2*a)
        sol2 = (-b+math.sqrt(d))/(2*a)
    
        bestfit = -1.*b/(2.*a)
    
        return (bestfit, sol1, sol2)
    
    def graphStyle(self, graph, style=20, color=ROOT.kOrange+7, size=1.0, titleY='-2 #Delta ln L', rangeY=(-0.01,4.0) ):
        graph.SetMarkerStyle(style)
        graph.SetMarkerColor(color)
        graph.SetLineWidth  (2)
        graph.SetMarkerSize(size)
        graph.GetYaxis().SetTitle(titleY)
        if rangeY:
            graph.GetYaxis().SetRangeUser(rangeY[0], rangeY[1])
    
    
    def getGraph(self, infile, par, norm, treename='limit'):
        f = ROOT.TFile(infile,'read')
        tree = f.Get(treename)
        vals = []
        normval = norm if norm else 1.
        for ev in tree:
            vals.append( [getattr(ev, par)/normval, 2.*ev.deltaNLL] )
        vals = sorted(vals)
        graph = ROOT.TGraph(len(vals), array.array('d', [x[0] for x in vals]), array.array('d', [y[1] for y in vals]) )
        self.graphStyle(graph)
        graph.GetXaxis().SetTitle(par)
        graph.SetTitle('scan for '+par)
        return graph

    def getErrorFromGraph(self, graph):
        graph.Fit('pol2')
        tmp_fit = graph.GetFunction('pol2')
        (best, sol1, sol2) = self.solvePol2(tmp_fit.GetParameter(2), tmp_fit.GetParameter(1), tmp_fit.GetParameter(0)-1)
        return (best, sol1, sol2)

    def getRebinned(self, ybins, charge, infile, ip):
        histo_file = ROOT.TFile(infile, 'READ')
    
        pstr = 'central' if not ip else 'pdf{ip}'.format(ip=ip)
    
        histos = {}
        for pol in ['left','right','long']:
            cp = '{ch}_{pol}'.format(ch=charge,pol=pol if not pol == 'long' else 'right')
    
            keys = histo_file.GetListOfKeys()
            for k in keys:
                if 'w{ch}'.format(ch=charge) in k.GetName() and pol in k.GetName() and pstr in k.GetName():
                    name = k.GetName()
            histo = histo_file.Get(name)# 'w{ch}_wy_W{ch}_{pol}'.format(ch=charge, pol=pol))
            conts = []
            epsilon = 0.000001 
            # val is a bin boundary, add epsilon to be sure to catch the correct bin (float can be truncated and findBin might return incorrect bin)
            for iv, val in enumerate(ybins[cp][:-1]):
                err = ROOT.Double()
                istart = histo.FindBin(val+epsilon)
                iend   = histo.FindBin(ybins[cp][iv+1]+epsilon)
                val = histo.IntegralAndError(istart, iend-1, err)/36000. ## do not include next bin
                conts.append(float(val))
            histos[pol] = conts
        histo_file.Close()
        return histos

    def getXSecFromShapes(self, ybins, charge, infile, ip, nchannels=1, polarizations = ['left','right','long'], excludeYbins = []):
        values = {}
        if not infile:
            for pol in polarizations: 
                cp = '{ch}_{pol}'.format(ch=charge,pol=pol)
                xsecs = []
                for iv in xrange(len(ybins[cp][:-1])):
                    if any(iv == x for x in excludeYbins): continue
                    xsecs.append(0.)
                values[pol] = xsecs
            return values

        histo_file = ROOT.TFile(infile, 'READ')
    
        pstr = '' if not ip else '_pdf{ip}Up'.format(ip=ip)
    
        for pol in polarizations:
            cp = '{ch}_{pol}'.format(ch=charge,pol=pol)
            xsecs = []
            for iv, val in enumerate(ybins[cp][:-1]):                
                if any(iv == x for x in excludeYbins): continue
                name = 'x_W{ch}_{pol}_Ybin_{iy}{suffix}'.format(ch=charge,pol=pol,iy=iv,ip=ip,suffix=pstr)                
                #print name
                histo = histo_file.Get(name)
                val = histo.Integral()/36000./float(nchannels) # xsec file yields normalized to 36 fb-1
                xsecs.append(float(val))
            values[pol] = xsecs
        histo_file.Close()
        return values

    def getPDFbandFromXsec(self, histoPDF, charge, infile, netabins, nptbins, firstPtBin=0):

        print "Inside getPDFbandFromXsec() ..."
        histo_file = ROOT.TFile(infile, 'READ')            

        for ieta in range(netabins):
            for ipt in range(firstPtBin,nptbins):

                nomi = "x_W{ch}_lep_ieta_{ie}_ipt_{ip}".format(ch=charge, ie=ieta, ip=ipt)
                hnomi = histo_file.Get(nomi)
                if not hnomi:
                    print "Error in getPDFbandFromXsec(): I couldn't find histogram " + nomi
                    quit()

                hpdftmp = None            
                pdfquadrsum = 0.0 
                xsecnomi = hnomi.Integral()                

                for ipdf in range(1, 61):
                    pdfvar = nomi + '_pdf{ip}Up'.format(ip=ipdf)
                    hpdftmp = histo_file.Get(pdfvar)
                    if not hpdftmp:
                        print "Error in getPDFbandFromXsec(): I couldn't find histogram " + pdfvar
                        quit()                    
                    tmpval = hpdftmp.Integral() - xsecnomi                     
                    pdfquadrsum += tmpval * tmpval

                histoPDF.SetBinError(ieta+1, ipt+1, math.sqrt(pdfquadrsum))
                histoPDF.SetBinContent(ieta+1, ipt+1, xsecnomi)
                        
        histo_file.Close()
        return 0


    def getPDFbandFromXsecEta(self, histoPDF, charge, infile, netabins, nptbins, firstPtBin=0):

        print "Inside getPDFbandFromXsecEta() ..."
        histo_file = ROOT.TFile(infile, 'READ')            

        for ieta in range(netabins):
            pdfquadrsum = 0.0
            xsecnomi = 0.0
            xsecpdf = [0.0 for i in range(60)]
            for ipt in range(firstPtBin,nptbins):

                nomi = "x_W{ch}_lep_ieta_{ie}_ipt_{ip}".format(ch=charge, ie=ieta, ip=ipt)
                hnomi = histo_file.Get(nomi)
                if not hnomi:
                    print "Error in getPDFbandFromXsecEta(): I couldn't find histogram " + nomi
                    quit()

                hpdftmp = None            
                xsecnomi += hnomi.Integral()                

                for ipdf in range(1, 61):
                    pdfvar = nomi + '_pdf{ip}Up'.format(ip=ipdf)
                    hpdftmp = histo_file.Get(pdfvar)
                    if not hpdftmp:
                        print "Error in getPDFbandFromXsecEta(): I couldn't find histogram " + pdfvar
                        quit()
                    
                    xsecpdf[ipdf-1] += hpdftmp.Integral()

            pdfquadrsum = 0.0
            for ipdf in range(60):
                tmpval = xsecpdf[ipdf] - xsecnomi 
                pdfquadrsum += tmpval*tmpval
            histoPDF.SetBinError(ieta+1, math.sqrt(pdfquadrsum))
            histoPDF.SetBinContent(ieta+1, xsecnomi)                        

        histo_file.Close()
        return 0



    def getParametersFromWS(self, ws, regexp):

        ## get all the nuisance parameters from the workspace
        pars = ws.allVars()
        pars = ROOT.RooArgList(pars)
        ## this has to be a loop over a range... doesn't work otherwise
        parameters = []
        all_parameters = []
        pois_regexps = list(regexp.split(','))
        ## get the parameters to scan from the list of allVars and match them
        ## to the given regexp
        for i in range(len(pars)):
            tmp_name = pars[i].GetName()
            if '_In' in tmp_name: ## those are the input parameters
                continue
            if tmp_name in ['CMS_th1x', 'r']: ## don't want those
                continue
            all_parameters.append(tmp_name)
            for poi in pois_regexps:
                if re.match(poi, tmp_name):
                    parameters.append(pars[i].GetName())

        return parameters

    def translateWStoTF(self, pname):
        if not 'r_' in pname:
            return pname
     
        if 'Wplus' in pname or 'Wminus' in pname:
            pnew = pname.replace('r_','')
            pnew = pnew.split('_')
            pnew.insert(-2, 'mu' )#if 'Wmu' in options.tensorflow else 'el')
            pnew = '_'.join(pnew)
            return pnew
     
        return -1

    def getFromHessian(self, infile, keepGen=False):
        _dict = {}
        
        f = ROOT.TFile(infile, 'read')
        tree = f.Get('fitresults')
        lok  = tree.GetListOfLeaves()
        for p in lok:
            if '_err'   in p.GetName(): continue
            if '_minos' in p.GetName(): continue
            if '_gen'   in p.GetName() and not keepGen: continue
            if '_In'    in p.GetName(): continue

            if not p.GetName()+'_err' in lok and not keepGen: continue

            if tree.GetEntries() > 1:
                print 'YOUR INPUT FILE HAS MORE THAN ONE FIT INSIDE. THIS IS PROBABLY NOT A HESSIAN FILE!!!'
                sys.exit()
            for ev in tree:
                mean = getattr(ev, p.GetName())
                err  = getattr(ev, p.GetName()+'_err') if hasattr(ev, p.GetName()+'_err') else 0

            _dict[p.GetName()] = (mean, mean+err, mean-err)
     
        return _dict


    def getFromToys(self, infile, keepGen=False, params=[]):
        _dict = {}
        
        f = ROOT.TFile(infile, 'read')
        tree = f.Get('fitresults')
        lok  = tree.GetListOfLeaves()
        
        for p in lok:
            if '_err'   in p.GetName(): continue
            if '_minos' in p.GetName(): continue
            if '_gen'   in p.GetName() and not keepGen: continue
            if '_In'    in p.GetName(): continue

            if len(params):
                match = [re.match(param,p.GetName()) for param in params]
                if not any(match): continue

            print 'gettin parameter ', p.GetName(), 'from toys file'
            
            tmp_hist = ROOT.TH1F(p.GetName(),p.GetName(), 100000, -5000., 5000.)
            tree.Draw(p.GetName()+'>>'+p.GetName())
            #tmp_hist = ROOT.gPad.GetPrimitive('foob')
            mean = tmp_hist.GetMean()
            err  = tmp_hist.GetRMS()
            _dict[p.GetName()] = (mean, mean+err, mean-err)
            del tmp_hist

        return _dict

    def getHistosFromToys(self, infile, nbins=100, xlow=-3.0, xup=3.0, getPull=False, matchBranch=None,excludeBranch=None, selection="", 
                          setStatOverflow=False, getMedian=False):

        # getPull = True will return a histogram centered at 0 and with expected rms=1, obtained as (x-x_gen)/x_err

        _dict = {}
        
        f = ROOT.TFile(infile, 'read')
        tree = f.Get('fitresults')
        lok  = tree.GetListOfLeaves()

        #nMaxBranch = 10
        #np = 0
        im = 1  # for median

        tree.SetBranchStatus("*",0)  # disabling and enabling branches makes the loop on events faster, at least if the median is used

        for p in lok:

            tree.SetBranchStatus(p.GetName(),1)
            
            if '_err'   in p.GetName(): continue
            if '_minos' in p.GetName(): continue
            if '_gen'   in p.GetName(): continue
            if '_In'    in p.GetName(): continue
            if matchBranch and not any(re.match(poi,p.GetName()) for poi in matchBranch.split(',')): continue
            if excludeBranch and any(re.match(excl,p.GetName()) for excl in excludeBranch.split(',')): continue

            # mainly for tests
            #if np == nMaxBranch: break
            #np += 1
            
            #print "Loading parameter --> %s " % p.GetName()
            #self.cmssw_version = os.environ['CMSSW_VERSION']
            #self.isRecentRelease = (len(self.cmssw_version) and int(self.cmssw_version.split('_')[1]) > 8)

            if getPull and (p.GetName()+"_gen") in lok and (p.GetName()+"_err") in lok:                
                tree.SetBranchStatus(p.GetName()+"_gen",1)
                tree.SetBranchStatus(p.GetName()+"_err",1)
                #print " Making pull --> (x-x_gen)/x_err for parameter %s" % p.GetName()
                tmp_hist_tmp = ROOT.TH1F(p.GetName()+"_tmp",p.GetName()+"_tmp", nbins, xlow, xup)
                if self.isRecentRelease:
                    if setStatOverflow: tmp_hist_tmp.SetStatOverflows(1)
                    else              : tmp_hist_tmp.SetStatOverflows(0)
                tmp_hist = ROOT.TH1F(p.GetName(),p.GetName(), 100, -3, 3)
                expression = "({p}-{pgen})/{perr}".format(p=p.GetName(),pgen=p.GetName()+"_gen",perr=p.GetName()+"_err")
                tree.Draw(expression+'>>'+p.GetName(),selection)
                tree.Draw(p.GetName()+'>>'+p.GetName()+"_tmp",selection)
                mean = tmp_hist_tmp.GetMean()
                err  = tmp_hist_tmp.GetRMS()
            else:
                tmp_hist = ROOT.TH1F(p.GetName(),p.GetName(), nbins, xlow, xup)
                if self.isRecentRelease:
                    if setStatOverflow: tmp_hist.SetStatOverflows(1)
                    else              : tmp_hist.SetStatOverflows(0)
                tree.Draw(p.GetName()+'>>'+p.GetName(),selection)
                mean = tmp_hist.GetMean()
                err  = tmp_hist.GetRMS()
            tmp_hist.SetDirectory(None)

            if getMedian:
                print "{n}) Computing median for {pn}".format(n=im,pn=p.GetName())
                im += 1
                vals = []
                #binCount = 1
                tot = tree.GetEntries()
                for ev in tree:
                    #sys.stdout.write('Bin {num}/{tot}   \r'.format(num=binCount,tot=tot))
                    #sys.stdout.flush()
                    #binCount += 1
                    vals.append(getattr(ev, p.GetName()))
                    
                vals.sort()
                nElements = len(vals)            
                if nElements%2: median = vals[(nElements-1)/2]
                else:           median = 0.5 * (vals[nElements/2] +  vals[nElements/2 + 1])
                _dict[p.GetName()] = (median, median+err, median-err, tmp_hist)
            else:
                _dict[p.GetName()] = (mean, mean+err, mean-err, tmp_hist)

            tree.SetBranchStatus(p.GetName(),0)
            tree.SetBranchStatus(p.GetName()+"_gen",0)
            tree.SetBranchStatus(p.GetName()+"_err",0)
     
        return _dict



    def getExprFromToys(self, name, expression, infile):
        f = ROOT.TFile(infile, 'read')
        tree = f.Get('fitresults')        
        tmp_hist = ROOT.TH1F(name,name, 100000, -100., 5000.)
        tree.Draw(expression+'>>'+name)
        mean = tmp_hist.GetMean()
        err  = tmp_hist.GetRMS()
        return (mean, mean+err, mean-err)


    def getNormalizedXsecFromToys(self, ybins, charge, pol, channel, iy, infile, absYmax=6.0):
        cp = '{ch}_{pol}'.format(ch=charge,pol=pol)
        ybins_expr = []
        for allpol in ['left','right', 'long']:
            for iv, val in enumerate(ybins[cp][:-1]):
                if abs(val)<absYmax:
                    ybins_expr.append('W{charge}_{pol}_W{charge}_{pol}_{ch}_Ybin_{iy}_pmaskedexp'.format(charge=charge,pol=allpol,ch=channel,iy=iv))
        num = 'W{charge}_{pol}_W{charge}_{pol}_{ch}_Ybin_{iy}_pmaskedexp'.format(charge=charge,pol=pol,ch=channel,iy=iy)
        den = '('+'+'.join(ybins_expr)+')'        
        ret = self.getExprFromToys(charge+pol+channel+str(iy),'{num}/{den}'.format(num=num,den=den),infile)
        return ret

    def getAsymmetryFromToys(self, pol, channel, iy, infile):
        expr = '(Wplus_{pol}_Wplus_{pol}_{ch}_Ybin_{iy}_pmaskedexp - Wminus_{pol}_Wminus_{pol}_{ch}_Ybin_{iy}_pmaskedexp)/(Wplus_{pol}_Wplus_{pol}_{ch}_Ybin_{iy}_pmaskedexp + Wminus_{pol}_Wminus_{pol}_{ch}_Ybin_{iy}_pmaskedexp)'.format(pol=pol,ch=channel,iy=iy)
        ret = self.getExprFromToys('chargeAsym',expr,infile)
        return ret

    def getDiffXsecAsymmetryFromToys(self, ieta, ipt, infile):
        xplus  = "Wplus_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(ieta=ieta,ipt=ipt)
        xminus = "Wminus_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(ieta=ieta,ipt=ipt)
        expr = '({pl}-{mn})/({pl}+{mn})'.format(pl=xplus,mn=xminus)
        ret = self.getExprFromToys('chargeAsym',expr,infile)
        return ret

    def getDenExpressionForNormDiffXsec(self, charge, netabins, nptbins):
        binsToNormalize = []
        for ieta in range(netabins):
            for ipt in range(nptbins):
                denChunk = "W{c}_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(c=charge,ieta=ieta,ipt=ipt)
                binsToNormalize.append(denChunk)
        den = "+".join(x for x in binsToNormalize)
        den = "(" + den + ")"
        return den

    def getNormalizedDiffXsecFromToys(self, charge, ieta, ipt, infile, den, friendTree=""):
        num = "W{c}_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(c=charge,ieta=ieta,ipt=ipt)
        expr = '{num}/{den}'.format(num=num,den=den)
        ret = self.getExprFromToys('normDiffXsec',expr,infile)
        return ret


    def getDiffXsecFromToys(self, charge, ieta, ipt, infile):
        expr = "W{c}_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(c=charge,ieta=ieta,ipt=ipt)
        ret = self.getExprFromToys('diffXsec',expr,infile)
        return ret

    def getExprFromToysFast(self, name, expression, nHistBins=100000, minHist=-100., maxHist=5000., tree=None):
        tmp_hist = ROOT.TH1F(name,name, nHistBins, minHist, maxHist)
        #self.cmssw_version = os.environ['CMSSW_VERSION']
        #self.isRecentRelease = (len(self.cmssw_version) and int(self.cmssw_version.split('_')[1]) > 8)
        if self.isRecentRelease: tmp_hist.SetStatOverflows(1)
        tree.Draw(expression+'>>'+name)
        mean = tmp_hist.GetMean()
        err  = tmp_hist.GetRMS()
        return (mean, mean+err, mean-err)

    def getNormalizedDiffXsecFromToysFast(self, charge, ieta, ipt, den, nHistBins=1000, minHist=0., maxHist=0.1, tree=None):
        num = "W{c}_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(c=charge,ieta=ieta,ipt=ipt)
        expr = '{num}/{den}'.format(num=num,den=den)
        ret = self.getExprFromToysFast('normDiffXsec',expr, nHistBins, minHist, maxHist, tree=tree)
        return ret

    def getDiffXsecFromToysFast(self, charge, ieta, ipt, nHistBins=2000, minHist=0., maxHist=200., tree=None):
        expr = "W{c}_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(c=charge,ieta=ieta,ipt=ipt)
        ret = self.getExprFromToysFast('diffXsec',expr,nHistBins, minHist, maxHist, tree=tree)
        return ret

    def getSignalStrengthFromToysFast(self, channel, charge, ieta, ipt, nHistBins=400, minHist=0., maxHist=2., tree=None):
        expr = "W{c}_lep_ieta_{ieta}_ipt_{ipt}_mu".format(c=charge,ieta=ieta,ipt=ipt)
        ret = self.getExprFromToysFast('diffXsec',expr, nHistBins, minHist, maxHist, tree=tree)
        return ret


    def getDiffXsecAsymmetryFromToysFast(self, ieta, ipt, nHistBins=2000, minHist=0., maxHist=1.0, tree=None):
        xplus  = "Wplus_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(ieta=ieta,ipt=ipt)
        xminus = "Wminus_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(ieta=ieta,ipt=ipt)
        expr = '({pl}-{mn})/({pl}+{mn})'.format(pl=xplus,mn=xminus)
        ret = self.getExprFromToysFast('chargeAsym',expr, nHistBins, minHist, maxHist, tree=tree)
        return ret

    def getExprFromHessian(self, name, expression, infile):
        f = ROOT.TFile(infile, 'read')
        tree = f.Get('fitresults')        
        tmp_hist = ROOT.TH1F(name,name, 100000, -100., 5000.)
        tree.Draw(expression+'>>'+name)
        mean = tmp_hist.GetMean()  # if this is hessian and not toys, there is just one entry, so the mean is the entry
        #err  = tmp_hist.GetRMS()  # not used in this context, (we are going to use this expression mainly for charge asymmetry, the uncertainty must be taken from toys)
        return mean

    def getDiffXsecAsymmetryFromHessian(self, ieta, ipt, infile):
        xplus  = "Wplus_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(ieta=ieta,ipt=ipt)
        xminus = "Wminus_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(ieta=ieta,ipt=ipt)
        expr = '({pl}-{mn})/({pl}+{mn})'.format(pl=xplus,mn=xminus)
        ret = self.getExprFromHessian('chargeAsym',expr,infile)
        return ret


    def getNormalizedDiffXsecFromHessian(self, charge, ieta, ipt, infile,den):
        num = "W{c}_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(c=charge,ieta=ieta,ipt=ipt)
        expr = '{num}/{den}'.format(num=num,den=den)
        ret = self.getExprFromHessian('normDiffXsec',expr,infile)
        return ret

    def getDiffXsecFromHessian(self, charge, ieta, ipt, infile):
        expr = "W{c}_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(c=charge,ieta=ieta,ipt=ipt)
        ret = self.getExprFromHessian('diffXsec',expr,infile)
        return ret

    ######## HESSIAN FAST

    def getExprFromHessianFast(self, name, expression, nHistBins=100000, minHist=-100., maxHist=5000., tree=None):
        tmp_hist = ROOT.TH1F(name,name, nHistBins, minHist, maxHist)
        if self.isRecentRelease: tmp_hist.SetStatOverflows(1)
        tree.Draw(expression+'>>'+name)
        mean = tmp_hist.GetMean()  # if this is hessian and not toys, there is just one entry, so the mean is the entry
        return mean

    def getDiffXsecAsymmetryFromHessianFast(self, ieta, ipt, nHistBins=2000, minHist=0., maxHist=1.0, tree=None, getErr=False, getGen=False):
        expr = "W_lep_ieta_{ieta}_ipt_{ipt}_chargeasym".format(ieta=ieta,ipt=ipt)
        if getErr: expr += "_err"
        elif getGen: expr += "_gen"
        ret = self.getExprFromHessianFast('chargeAsym',expr, nHistBins, minHist, maxHist, tree=tree)
        return ret


    def getNormalizedDiffXsecFromHessianFast(self, charge, ieta, ipt, nHistBins=1000, minHist=0., maxHist=0.1, tree=None, getErr=False, getGen=False):
        expr = "W{c}_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexpnorm".format(c=charge,ieta=ieta,ipt=ipt)
        if getErr: expr += "_err"
        elif getGen: expr += "_gen"
        ret = self.getExprFromHessianFast('normDiffXsec',expr, nHistBins, minHist, maxHist, tree=tree)
        return ret

    def getDiffXsecFromHessianFast(self, charge, ieta, ipt,  nHistBins=2000, minHist=0., maxHist=200., tree=None, getErr=False, getGen=False):
        expr = "W{c}_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(c=charge,ieta=ieta,ipt=ipt)
        if getErr: expr += "_err"
        elif getGen: expr += "_gen"
        ret = self.getExprFromHessianFast('diffXsec',expr, nHistBins, minHist, maxHist, tree=tree)
        return ret


    # this is for single charge
    def getDiffXsec1DFromHessianFast(self, charge, ietaORipt, isIeta=True, nHistBins=5000, minHist=0., maxHist=5000., tree=None, getErr=False, getGen=False):
        expr = "W{c}_lep_i{var}_{ivar}_sumxsec".format(c=charge,var="eta" if isIeta else "pt", ivar=ietaORipt)
        if getErr: expr += "_err"
        elif getGen: expr += "_gen"
        ret = self.getExprFromHessianFast('diffXsec1D_{var}'.format(var="eta" if isIeta else "pt"),expr, nHistBins, minHist, maxHist, tree=tree)
        return ret

    def getNormalizedDiffXsec1DFromHessianFast(self, charge, ietaORipt, isIeta=True, nHistBins=1000, minHist=0., maxHist=1., tree=None, getErr=False, getGen=False):
        expr = "W{c}_lep_i{var}_{ivar}_sumxsecnorm".format(c=charge,var="eta" if isIeta else "pt", ivar=ietaORipt)
        if getErr: expr += "_err"
        elif getGen: expr += "_gen"
        ret = self.getExprFromHessianFast('diffXsec1D_{var}'.format(var="eta" if isIeta else "pt"),expr, nHistBins, minHist, maxHist, tree=tree)
        return ret


    def getDiffXsecAsymmetry1DFromHessianFast(self, ietaORipt, isIeta=True, 
                                              nHistBins=1000, minHist=0., maxHist=1., tree=None, getErr=False, getGen=False):
        expr = "W_lep_i{var}_{ivar}_chargemetaasym".format(var="eta" if isIeta else "pt", ivar=ietaORipt)
        if getErr: expr += "_err"
        elif getGen: expr += "_gen"
        ret = self.getExprFromHessianFast('chargeAsym1D_{var}'.format(var="eta" if isIeta else "pt"),expr, nHistBins, minHist, maxHist, tree=tree)
        return ret


    def getSignalStrengthFromHessianFast(self, charge, ieta, ipt, nHistBins=200, minHist=0.5, maxHist=1.5, tree=None, getErr=False, getGen=False):
        expr = "W{c}_lep_ieta_{ieta}_ipt_{ipt}_mu".format(c=charge,ieta=ieta,ipt=ipt)
        if getErr: expr += "_err"
        elif getGen: expr += "_gen"
        ret = self.getExprFromHessianFast('normDiffXsec',expr, nHistBins, minHist, maxHist, tree=tree)
        return ret


    #################################### 

    def getFromScans(self, indir):
        _dict = {}
        
        for sd in os.listdir(indir):
     
            if 'jobs' in sd: continue
     
            par = self.translateWStoTF(sd) ## parameter name different than in TF
            f = ROOT.TFile(indir+'/'+sd+'/scan_'+sd+'.root', 'read')
            tree = f.Get('fitresults')
     
            vals = []
            for ev in tree:
                vals.append( [getattr(ev, par), 2.*ev.nllval  ] )
            vals = sorted(vals)
            lvals = vals[:len(vals)/2]
            rvals = vals[len(vals)/2:]
     
            graph = ROOT.TGraph(len(vals), array.array('d', [x[1] for x in vals]), array.array('d', [y[0] for y in vals]) )
     
            best = graph.Eval(0.)
            lgraph = ROOT.TGraph(len(lvals), array.array('d', [x[1] for x in lvals]), array.array('d', [y[0] for y in lvals]) )
            rgraph = ROOT.TGraph(len(rvals), array.array('d', [x[1] for x in rvals]), array.array('d', [y[0] for y in rvals]) )
            sol1  = lgraph.Eval(1.)
            sol2  = rgraph.Eval(1.)
     
            _dict[par] = (best, sol1, sol2)
     
        return _dict


    def effSigma(self, histo):
        xaxis = histo.GetXaxis()
        nb = xaxis.GetNbins()
        xmin = xaxis.GetXmin()
        ave = histo.GetMean()
        rms = histo.GetRMS()
        total=histo.Integral()
        if total < 100: 
            print "effsigma: Too few entries to compute it: ", total
            return 0.
        ierr=0
        ismin=999
        rlim=0.683*total
        bwid = xaxis.GetBinWidth(1)
        nrms=int(rms/bwid)
        if nrms > nb/10: nrms=int(nb/10) # Could be tuned...
        widmin=9999999.
        for iscan in xrange(-nrms,nrms+1): # // Scan window centre 
            ibm=int((ave-xmin)/bwid)+1+iscan
            x=(ibm-0.5)*bwid+xmin
            xj=x; xk=x;
            jbm=ibm; kbm=ibm;
            bin=histo.GetBinContent(ibm)
            total=bin
            for j in xrange(1,nb):
                if jbm < nb:
                    jbm += 1
                    xj += bwid
                    bin=histo.GetBinContent(jbm)
                    total += bin
                    if total > rlim: break
                else: ierr=1
                if kbm > 0:
                    kbm -= 1
                    xk -= bwid
                    bin=histo.GetBinContent(kbm)
                    total+=bin
                if total > rlim: break
                else: ierr=1
            dxf=(total-rlim)*bwid/bin
            wid=(xj-xk+bwid-dxf)*0.5
            if wid < widmin:
                widmin=wid
                ismin=iscan
        if ismin == nrms or ismin == -nrms: ierr=3
        if ierr != 0: print "effsigma: Error of type ", ierr
        return widmin

    def doShadedUncertainty(self,h):
        xaxis = h.GetXaxis()
        points = []; errors = []
        for i in xrange(h.GetNbinsX()):
            N = h.GetBinContent(i+1); dN = h.GetBinError(i+1);
            if N == 0 and dN == 0: continue
            x = xaxis.GetBinCenter(i+1);
            points.append( (x,N) )
            EYlow, EYhigh  = dN, min(dN,N);
            EXhigh, EXlow = (xaxis.GetBinUpEdge(i+1)-x, x-xaxis.GetBinLowEdge(i+1))
            errors.append( (EXlow,EXhigh,EYlow,EYhigh) )
        ret = ROOT.TGraphAsymmErrors(len(points))
        ret.SetName(h.GetName()+"_errors")
        for i,((x,y),(EXlow,EXhigh,EYlow,EYhigh)) in enumerate(zip(points,errors)):
            ret.SetPoint(i, x, y)
            ret.SetPointError(i, EXlow,EXhigh,EYlow,EYhigh)
        ret.SetFillStyle(3244);
        ret.SetFillColor(ROOT.kGray+2)
        ret.SetMarkerStyle(0)
        ret.Draw("PE2 SAME")
        return ret

    def safecolor(self, index):
        SAFE_COLOR_LIST=[ROOT.kBlack, ROOT.kRed, ROOT.kGreen+2, ROOT.kBlue, ROOT.kMagenta+1, ROOT.kOrange+7, ROOT.kCyan+1, ROOT.kGray+2, ROOT.kViolet+5, ROOT.kSpring+5, ROOT.kAzure+1, ROOT.kPink+7, ROOT.kOrange+3, ROOT.kBlue+3, ROOT.kMagenta+3, ROOT.kRed+2]+range(11,40)
        if index<len(SAFE_COLOR_LIST): return SAFE_COLOR_LIST[index]
        else: return index

    def getCoeffs(self,xL,xR,x0,err_xL,err_xR,err_x0, toyEvents=10000):        
        histos = { 'a0': ROOT.TH1D('a0','',100,-0.4,0.6),
                   'a4': ROOT.TH1D('a4','',100,-0.4,3.0)
                   }
        # given that we use mean and RMS, for security reason use Under/Overflow values to compute them
        # so we don't have to tune the range (even though the one about is already very sensible for a0 and a4)
        histos['a0'].SetStatOverflows(1)
        histos['a4'].SetStatOverflows(1)
        #print "getCoeffs: ",toyEvents," toyMC running..."
        for i in xrange(toyEvents):
            ixL = np.random.normal(xL,err_xL)
            ixR = np.random.normal(xR,err_xR)
            ix0 = np.random.normal(x0,err_x0)
            sumPol = ixL+ixR+ix0
            histos['a0'].Fill(2*ix0/sumPol)
            #histos['a0'].Fill(2*(1-ixL-ixR)/sumPol)
            histos['a4'].Fill(2*(ixL-ixR)/sumPol)
        #print "toyMC done"
        ret = {}
        for k,h in histos.iteritems():
            ret[k] = (h.GetMean(),h.GetRMS())
        return ret
            
    def getChargeAsy(self,xplus,xminus,err_xplus,err_xminus,toyEvents=10000):
        histo = ROOT.TH1D('asy','',100,-0.05,0.5)
        #print "getChargeAsy: ",toyEvents," toyMC running..."
        for i in xrange(toyEvents):
            ixplus  = np.random.normal(xplus,err_xplus)
            ixminus = np.random.normal(xminus,err_xminus)
            histo.Fill((ixplus-ixminus)/(ixplus+ixminus))
        #print "toyMC done"
        ret = {'asy': (histo.GetMean(),histo.GetRMS())}
        return ret

    def getChargeAsyFromTH1pair(self,h1, h2, toyEvents=10000, name='asy'):
        # assuming h1 and h2 have same binning
        print "Inside getChargeAsyFromTH1pair()"
        binsx = [h1.GetXaxis().GetBinLowEdge(i) for i in range(1,2+h1.GetNbinsX())]
        histo = ROOT.TH1D(name,'',len(binsx)-1, array.array('d',binsx))
        histotmp = ROOT.TH1D(name+'tmp','',200,-1.0,1.0)
        for i in range(1,1+h1.GetNbinsX()):
            histotmp.Reset("ICESM")
            for j in xrange(toyEvents):
                ixplus  = np.random.normal(h1.GetBinContent(i),h1.GetBinError(i))
                ixminus = np.random.normal(h2.GetBinContent(i),h2.GetBinError(i))
                histotmp.Fill((ixplus-ixminus)/(ixplus+ixminus))
            histo.SetBinContent(i, histotmp.GetMean()  )
            histo.SetBinError(  i, histotmp.GetStdDev())
        #print "toyMC done"
        return histo


    def getNEffStat(self, s):
        a = s.split('EffStat')[1]
        a = a.replace('minus','').replace('plus','')
        a = a.replace('mu','').replace('el','')
        return int(a)

    def getNFromString(self, s):
        los = [ int(i) for i in re.findall(r'\d+', s) ]
        if len(los) == 0: return 0
        if len(los) == 1: return los[0]
        if len(los)  > 1: return los[1] if 'EffStat' in s else los[0]
        return 0
        

