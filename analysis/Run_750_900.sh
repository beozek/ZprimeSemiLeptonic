#!/bin/bash

text2workspace.py datacard_test.txt -o Ac_muon_750_900.root  -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel -m 125 --PO map='.*/TTbar_1:r_neg[1,0,20]' --PO map='.*/TTbar_2:r_pos=expr;;r_pos("53446/54239*@0*(100+@1)/(100-@1)",r_neg,Ac[-2,-5,0])' --PO verbose
