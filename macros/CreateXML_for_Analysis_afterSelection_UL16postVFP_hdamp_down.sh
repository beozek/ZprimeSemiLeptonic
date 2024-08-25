#!/bin/bash

#where UHH2 code installed
pathGL_code=/nfs/dust/cms/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/
#where (NOT MERGED) trees after postselection stored
path_data=/nfs/dust/cms/group/zprime-uhh/Analysis_UL16postVFP/hdamp_down/electron/workdir_Analysis_UL16postVFP_electron_hdamp_down/uhh2.AnalysisModuleRunner.

mkdir $pathGL_code/ZprimeSemiLeptonic/data/Skimming_datasets_UL16postVFP_electron_hdamp_down
cd $pathGL_code/ZprimeSemiLeptonic/data/Skimming_datasets_UL16postVFP_electron_hdamp_down


# MC
for sample_name in TTToSemiLeptonic_UL16postVFP TTToHadronic_UL16postVFP TTTo2L2Nu_UL16postVFP
do
    echo $sample_name
    $pathGL_code/scripts/create-dataset-xmlfile ${path_data}"MC."${sample_name}"*.root" MC_$sample_name.xml
    python $pathGL_code/scripts/crab/readaMCatNloEntries.py 10 MC_$sample_name.xml True
done

pwd
cd $pathGL_code/ZprimeSemiLeptonic/macros
