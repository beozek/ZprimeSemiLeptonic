#!/bin/bash

text2workspace.py datacard_muon_UL18_750_1000.txt -o Ac_UL18_muon_750_1000.root  -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel -m 125 --PO map='.*/TTbar_1:r_neg[1,0,20]' --PO map='.*/TTbar_2:r_pos=expr;;r_pos("3570287600/3601730000*@0*(100+@1)/(100-@1)",r_neg,Ac[-2,-5,0])' --PO verbose 
