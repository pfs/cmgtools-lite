#!/bin/bash

PO="--PO sigma_th=69e-6 --PO A=208 --PO eb_MC=0.6"
model="HiggsAnalysis.CombinedLimit.TopBtagInMediumModel:topBtags"

floatOtherPOIs="--floatOtherPOIs=1"
addOpts="--setParameters sigma=2.98,delta=0,eb=0.6 -t -1" #set to empty string to run observed


baseFitOpts="--robustFit=1 --setRobustFitStrategy=1"
baseFitOpts="${baseFitOpts} -m 172.5"
baseFitOpts="${baseFitOpts} --setParameterRanges sigma=0,5:eb=0,1:delta=0,1"
baseFitOpts="${baseFitOpts} ${floatOtherPOIs}"


echo "Make sure ${model} is in its place before running"

echo "Creating workspace"
text2workspace.py datacard.txt -P ${model} ${PO}


#single POI likelihoods
for p in sigma eb delta; do
    echo "[ $p ]"
    combine datacard.root -M MultiDimFit --algo singles --cl=0.68 ${baseFitOpts} ${addOpts}\
        -P ${p} -n ${p}_singles
    combine datacard.root -M MultiDimFit --algo grid --points 100 ${baseFitOpts} ${addOpts}\
        -P ${p} -n ${p}_grid
done

#2-parameter scan
echo "[ sigma vs eb ]"
combine datacard.root -M MultiDimFit --algo grid --points 2000 ${baseFitOpts} ${addOpts}\
     -P sigma -P eb -n sigma_eb_grid     