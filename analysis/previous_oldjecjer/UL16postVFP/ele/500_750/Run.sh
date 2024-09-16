#!/bin/bash

text2workspace.py datacard_ele_UL16postVFP_500_750_SR.txt -o Ac_UL16postVFP_500_750.root  -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel -m 125 --PO map='.*/TTbar_1:r_neg[1,0,20]' --PO map='.*/TTbar_2:r_pos=expr;;r_pos("6027718100/6082951200*@0*(100+@1)/(100-@1)",r_neg,Ac[-2,-5,0])' --PO verbose