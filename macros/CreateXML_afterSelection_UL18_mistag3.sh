#!/bin/bash
#where UHH2 code installed
pathGL_code=/nfs/dust/cms/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/
#where (NOT MERGED) trees after preselection stored
path_data=/nfs/dust/cms/group/zprime-uhh/Presel_UL18/workdir_Preselection_UL18/uhh2.AnalysisModuleRunner.
#path_data=/nfs/dust/cms/group/zprime-uhh/Analysis_UL18/electron/workdir_Zprime_Analysis_UL18_electron/uhh2.AnalysisModuleRunner.

# mkdir $pathGL_code/ZprimeSemiLeptonic/data/Skimming_datasets_UL18_muon_mistag
cd $pathGL_code/ZprimeSemiLeptonic/data/Skimming_datasets_UL18_muon_mistag
#mkdir $pathGL_code/ZprimeSemiLeptonic/data/Skimming_datasets_UL18_electron
#cd $pathGL_code/ZprimeSemiLeptonic/data/Skimming_datasets_UL18_electron


# #MC

for sample_name in RSGluonToTT_M-500_UL18 RSGluonToTT_M-1000_UL18 RSGluonToTT_M-1500_UL18 RSGluonToTT_M-2000_UL18 RSGluonToTT_M-2500_UL18 RSGluonToTT_M-3000_UL18 RSGluonToTT_M-3500_UL18 RSGluonToTT_M-4000_UL18 RSGluonToTT_M-4500_UL18 RSGluonToTT_M-5000_UL18 RSGluonToTT_M-5500_UL18 RSGluonToTT_M-6000_UL18

do
    echo $sample_name

       $pathGL_code/scripts/create-dataset-xmlfile ${path_data}"MC."${sample_name}"*.root" MC_$sample_name.xml
       python $pathGL_code/scripts/crab/readaMCatNloEntries.py 10 $sample_name.xml True
done

# # #DATA
#for sample_name in DATA_EGamma_RunA_UL18_blinded DATA_EGamma_RunB_UL18_blinded DATA_EGamma_RunC_UL18_blinded DATA_EGamma_RunD_UL18_blinded 
# for sample_name in DATA_EGamma_RunA_UL18 DATA_EGamma_RunB_UL18 DATA_EGamma_RunC_UL18 DATA_EGamma_RunD_UL18 DATA_SingleMuon_RunA_UL18 DATA_SingleMuon_RunB_UL18 DATA_SingleMuon_RunC_UL18 DATA_SingleMuon_RunD_UL18 

# do
#     echo $sample_name 
#     $pathGL_code/scripts/create-dataset-xmlfile ${path_data}"DATA."${sample_name}"*.root" $sample_name.xml
#     python $pathGL_code/scripts/crab/readaMCatNloEntries.py 10 $sample_name.xml True

# done

# pwd
# cd $pathGL_code/ZprimeSemiLeptonic/macros
