#!/bin/bash

text2workspace.py datacard_muon_UL16preVFP_750_1000_SR.txt -o Ac_UL16preVFP_750_1000.root  -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel -m 125 --PO map='.*/TTbar_1:r_neg[1,0,20]' --PO map='.*/TTbar_2:r_pos=expr;;r_pos("1190640500/1200535900*@0*(100+@1)/(100-@1)",r_neg,Ac[-2,-5,0])' --PO verbose 
