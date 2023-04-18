#! /usr/bin/env python
from ROOT import *
import ROOT
import sys
import numpy

# As an input, Higgs Combine takes a txt based file containing the observed and expected yields.


# list of DeltaY ROOT files 
input_files = {
    "DeltaY_UL18_muon_750_900_Semileptonic.root": "Semileptonic", 
    "DeltaY_UL18_muon_750_900_ttbar.root": "TTbar",
    } 
input_files2 = {
    "DeltaY_UL18_muon_750_900_WJets.root" : "WJets",
    "DeltaY_UL18_muon_750_900_ST.root": "ST", 
    "DeltaY_UL18_muon_750_900_DATA.root" : "data_obs", 
    "DeltaY_UL18_muon_750_900_Diboson.root": "Diboson", 
    "DeltaY_UL18_muon_750_900_DY.root": "DY", 
    "DeltaY_UL18_muon_750_900_QCD.root": "QCD", 
    }


# Output ROOT file
output_file = ROOT.TFile("DeltaY_UL18_muon_750_900.root", "RECREATE")


# Loop over the.q input root files

for input_file_name,histogram_name in input_files.items():
    # Opening root file
    root_file = ROOT.TFile(input_file_name, "READ")

    # Getting DeltaY histogram from the file
    hist = root_file.Get("px1")
    if hist:
        # The Clone() method creates a copy of the histogram. One can save multiple histograms with the same name from different input files in the same output file.
        # Clone the histogram and give it a different name
        output_file.cd()
        hist_clone = hist.Clone(histogram_name)
        # Writing the histogram to the output file
        hist_clone.Write()
    root_file.Close()
    
for input_file_name2,histogram_name in input_files2.items():
    # Opening root file
    root_file = ROOT.TFile(input_file_name2, "READ")

    # Getting DeltaY histogram from the file
    hist = root_file.Get("DeltaY_reco_mass")
    if hist:
        # The Clone() method creates a copy of the histogram. One can save multiple histograms with the same name from different input files in the same output file.
        # Clone the histogram and give it a different name
        output_file.cd()
        hist_clone = hist.Clone(histogram_name)
        # Writing the histogram to the output file
        hist_clone.Write()
    root_file.Close()

# #Close the output file
output_file.Close()


#  # Loop over the histograms in the root file
#     for key in root_file.GetListOfKeys():
#         #Getting the histogram
#         hist = key.ReadObj()
        
#         # Writing the histogram to the output root file
#         hist.Write()

# output_file.Close()


