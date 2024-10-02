#!/usr/bin/env bash
date

# rm comb*
# rm higgs*
# rm Ac.pdf
# rm r_neg.pdf
# rm impacts.json

declare -a POIS=(
  "r_neg"
  "Ac"
)

export WORKSPACE=Ac_UL18_ele_0_500.root
export VERBOSITY=0

export SetParameters="rgx{r.+}=1,Ac=0.7"
export SetParametersExplicit="r_neg=1,Ac=0.7"
export SetParameterRanges="rgx{r.+}=0.5,2:Ac=-5,5"
export redefineSignalPOIs="Ac,r_neg"



export ASIMOV="-t -1"

# combine -M MultiDimFit -d datacard_ele_UL18_0_500_SR.txt -m 125 --freezeParameters MH -n .robustHesse.part6_multiSignalModel_v3 --cminDefaultMinimizerStrategy 0 --robustHesse 1 --robustHesseSave 1 --saveFitResult
# combine -M FitDiagnostics -d datacard_ele_UL18_0_500_SR.txt


###IMPACTS
echo
echo
echo "DO INITIAL FIT - Performing a fit using the MultiDimFit method. It scans each parameter (defined by --algo singles) and saves the workspace after the fit."
echo
echo
combine -M MultiDimFit --algo singles -d $WORKSPACE -n .part3E.snapshot -v $VERBOSITY --redefineSignalPOIs $redefineSignalPOIs --setParameterRanges $SetParameterRanges --setParameters $SetParametersExplicit --robustFit 1 --cminDefaultMinimizerStrategy 0 -m 125 --saveWorkspace $ASIMOV
mv higgsCombine.part3E.snapshot.MultiDimFit.mH125.root higgsCombine_initialFit_Test.MultiDimFit.mH125.root

# combine -M MultiDimFit --algo grid --points 100 -d $WORKSPACE -n .AcScan -v $VERBOSITY --redefineSignalPOIs $redefineSignalPOIs --setParameterRanges Ac=-5,5 --setParameters $SetParameters --robustFit 1 --cminDefaultMinimizerStrategy 0 -m 125



echo
echo
echo "DO FITS - calculate the impact of each nuisance parameter on parameters of interest"
echo
echo
combineTool.py -M Impacts -d $WORKSPACE --robustFit 1 --doFits -m 125 --redefineSignalPOIs $redefineSignalPOIs  --setParameters $SetParametersExplicit --setParameterRanges $SetParameterRanges --cminDefaultMinimizerStrategy 0 $ASIMOV --parallel 20



echo
echo
echo "STAT ONLY UNCERTAINTY (ALL NUISANCES FROZEN) - performs another MultiDimFit, but this time with all constrained nuisance parameters frozen (--freezeParameters allConstrainedNuisances). This can provide insight into how the fit behaves when the nuisances are not allowed to float."
echo
echo
combine -M MultiDimFit --algo singles -d $WORKSPACE -v $VERBOSITY --redefineSignalPOIs $redefineSignalPOIs --setParameterRanges $SetParameterRanges --setParameters $SetParameters --robustFit 1 --cminDefaultMinimizerStrategy 0 -m 125 --saveWorkspace -n _paramFit_Test_allConstrainedNuisancesFrozen --freezeParameters allConstrainedNuisances $ASIMOV



echo
echo
echo "CREATE IMPACTS JSON"
echo
echo
combineTool.py -M Impacts -d $WORKSPACE -o impacts.json -m 125 --redefineSignalPOIs $redefineSignalPOIs



echo
echo
echo "CREATE PLOTS"
echo
echo
for POI in ${POIS[@]}; do
  plotImpacts.py -i impacts.json -o $POI --POI $POI
done
echo
echo

mkdir output_combine
mv higgs* impacts.json Ac_UL18_ele_0_500.root combine_logger.out output_combine