#!/bin/bash

# rm -f *.root
# rm -f *.json
#rm -f *.txt

#text2workspace.py Ac_750.txt -o Ac_750.root  -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel -m 125 --PO map='.*/Ttbar_1:r_neg[1,0,20]' --PO map='.*/Ttbar_2:r_pos=expr;;r_pos("68351/69410*@0*(100+@1)/(100-@1)",r_neg,Ac[-2,-5,0])' --PO verbose 
 
#text2workspace.py Ac_900.txt -o Ac_900.root  -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel -m 125 --PO map='.*/Ttbar_1:r_neg[1,0,20]' --PO map='.*/Ttbar_2:r_pos=expr;;r_pos("36099/36764*@0*(100+@1)/(100-@1)",r_neg,Ac[-2,-5,0])' --PO verbose


#text2workspace.py Ac_750_900.txt -o Ac_750_900.root  -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel -m 125 --PO map='.*/Ttbar_1:r_neg[1,0,20]' --PO map='.*/Ttbar_2:r_pos=expr;;r_pos("32252/32646*@0*(100+@1)/(100-@1)",r_neg,Ac[-2,-5,0])' --PO verbose

text2workspace.py Ac_750_900.txt -o Ac_750_900.root  -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel -m 125 --PO map='.*/Ttbar_1:r_neg[1,0,20]' --PO map='.*/Ttbar_2:r_pos=expr;;r_pos("53446/54239*@0*(100+@1)/(100-@1)",r_neg,Ac[-2,-5,0])' --PO verbose

#text2workspace.py Ac_ele.txt -o Ac_ele.root  -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel -m 125 --PO map='.*/Ttbar_1:r_neg[1,0,20]' --PO map='.*/Ttbar_2:r_pos=expr;;r_pos("14905/15171*@0*(100+@1)/(100-@1)",r_neg,Ac[-2,-5,0])' --PO verbose


#bash impacts_750.sh
#mv Ac.pdf Ac_750.pdf 
#mv r_neg.pdf r_neg_750.pdf 

#bash impacts_750_900.sh
#mv Ac.pdf Ac_750_900.pdf
#mv r_neg.pdf r_neg_750_900.pdf

#bash impacts_900.sh
#mv Ac.pdf Ac_900.pdf
#mv r_neg.pdf r_neg_900.pdf

bash impacts_750_900.sh
# mv Ac.pdf Ac_muon.pdf
# mv r_neg.pdf r_neg_muon.pdf

#bash impacts_ele.sh
#mv Ac.pdf Ac_ele.pdf
#mv r_neg.pdf r_neg_ele.pdf








 







