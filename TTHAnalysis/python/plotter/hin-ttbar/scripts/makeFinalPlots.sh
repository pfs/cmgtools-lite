#!/bin/bash

outDir=finalPlots
mkdir -p ${outDir}


function genlevel {
    python runHINTtbar.py --jetFlavor
    cp /eos/user/p/psilva/www/HIN-19-001/jet_studies/2019-10-30-inclusive-anyflavor/*.{png,pdf} ${outDir}
}

function impacts {
cat <<EOF

This is just a summary of instructions. You need to use a CMSSW version with combineHarvester to run this

plotImpacts.py -i /afs/cern.ch/work/m/mdunser/public/cmssw/heavyIons/CMSSW_9_4_6_patch1/src/CMGTools/TTHAnalysis/python/plotter/hin-ttbar/datacards_2019-10-21_forApproval/bdtcomb/combinedCard/impacts_obs.json  -o impacts_inc --cms-label-l "#it{leptonic}" --cms-label "Preliminary          #scale[0.8]{#it{1.7 nb^{-1} (#sqrt{s_{NN}}=5.02 TeV)}}" --transparent --blind --translate hin-ttbar/scripts/translate.json --per-page 15

plotImpacts.py -i /afs/cern.ch/work/m/mdunser/public/cmssw/heavyIons/CMSSW_9_4_6_patch1/src/CMGTools/TTHAnalysis/python/plotter/hin-ttbar/datacards_2019-10-21_forApproval_jetAnalysis/bdtcomb/combinedCard/impacts_obs.json  -o impacts_jet --cms-label-l "#it{leptonic+b-tagged}" --cms-label "Preliminary                 #scale[0.8]{#it{1.7 nb^{-1} (#sqrt{s_{NN}}=5.02 TeV)}}" --transparent --blind --translate hin-ttbar/scripts/translate.json --per-page 15

cp impacts_*pdf paper/

EOF

}
function postFit {

    python hin-ttbar/scripts/plotPostFitDistributions.py ~/work/HIN-19-001/fits/datacards_2019-10-21_forApproval/bdtcomb/combinedCard/fitDiagnosticsobs.root --outdir paper/ --inclusive
    python hin-ttbar/scripts/plotPostFitDistributions.py ~/work/HIN-19-001/fits/datacards_2019-10-21_forApproval/bdtcomb/combinedCard/fitDiagnosticsobs.root --outdir paper/ --inclusive --prefit

    python hin-ttbar/scripts/plotPostFitDistributions.py ~/work/HIN-19-001/fits/datacards_2019-10-21_forApproval_jetAnalysis/bdtcomb/combinedCard/fitDiagnosticsobs.root --outdir paper/
    python hin-ttbar/scripts/plotPostFitDistributions.py ~/work/HIN-19-001/fits/datacards_2019-10-21_forApproval_jetAnalysis/bdtcomb/combinedCard/fitDiagnosticsobs.root --outdir paper/ --prefit
    
    for ch in ee mm em; do
        cp paper/*${ch}_p*fit_BDT* ${outDir}
    done
    cp  paper/*btaggedmultiplicity* ${outDir}
}

function statAna {
    python hin-ttbar/scripts/plotStatAnalysisComparison.py
    for i in nllscan sig_toys r_toys; do
        cp StatAnaSummary/allFlavors*_${i}.png ${outDir}
    done
}


function simplePlots {
 
    #python runHINTtbar.py --simple

    #replot with prefit/postfit uncertainties from combine
    dcDir="/afs/cern.ch/work/m/mdunser/public/cmssw/heavyIons/CMSSW_9_4_6_patch1/src/CMGTools/TTHAnalysis/python/plotter/hin-ttbar/datacards_2019-10-21_forApproval_jetAnalysis/bdtcomb"
    for i in em sf leponZee leponZmm anyflavor; do
        plotterDir="/eos/user/p/psilva/www/HIN-19-001/simple_plots/2019-10-29-${i}"
        plotter="${plotterDir}/mll_AND_mllz_AND_llpt_AND_sphericity_AND_acoplan_AND_minmlb_AND_minmlb2_AND_centrality.root"
      
        #python hin-ttbar/scripts/replotSimple.py ${plotter} ${dcDir}; 


        for ext in pdf png; do
            if [[ ${i} == *"leponZ"* ]]; then
                cp ${plotterDir}/mllz_prefit.${ext} ${outDir}/mllz_${i}_prefit.${ext};
            else
                for d in llpt mll sphericity acoplan; do
                    cp ${plotterDir}/${d}_prefit.${ext} ${outDir}/${d}_${i}_prefit.${ext};
                    cp ${plotterDir}/${d}_postfit_bkgsubinset.${ext} ${outDir}/${d}_${i}_postfit_bkgsubinset.${ext};
                done
            fi
        done
    done
}

function xsecSummary {
    python hin-ttbar/scripts/xsecsummary.py
    mv xsecsummary.png ${outDir}
    mv xsecsummary.pdf ${outDir}
}

function sbplots {
    python hin-ttbar/scripts/plotSB.py /afs/cern.ch/work/m/mdunser/public/cmssw/heavyIons/CMSSW_9_4_6_patch1/src/CMGTools/TTHAnalysis/python/plotter/hin-ttbar/datacards_2019-10-21_forApproval/bdtcomb leptonic
    python hin-ttbar/scripts/plotSB.py /afs/cern.ch/work/m/mdunser/public/cmssw/heavyIons/CMSSW_9_4_6_patch1/src/CMGTools/TTHAnalysis/python/plotter/hin-ttbar/datacards_2019-10-21_forApproval_jetAnalysis/bdtcomb btagged

    mv sob*.p* ${outDir}
}

#genlevel
#impacts
#statAna
#simplePlots
#sbplots

#paper plots
postFit
xsecSummary

