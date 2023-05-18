import ROOT


# Open the ROOT file
file = ROOT.TFile.Open("uhh2.AnalysisModuleRunner.Semileptonic.root")

# Get the TTree from the file
tree1 = file.Get("AnalysisTree")
tree2 = file.Get("AnalysisTree")
tree3 = file.Get("AnalysisTree")


# Create a histogram
hist1 = ROOT.TH1F("hist1", "Histogram", 100, -10, 10)
hist2 = ROOT.TH1F("hist2", "Histogram", 100, -10, 10)
hist3 = ROOT.TH1F("hist3", "Histogram", 100, -10, 10)


# Fill the histogram with the values from a branch in the tree
tree1.Draw("weight_pu>>hist1")
tree2.Draw("weight_pu_up>>hist2")
tree3.Draw("weight_pu_down>>hist3")

# Get the mean of the histogram
mean1 = hist1.GetMean()
mean2 = hist2.GetMean()
mean3 = hist3.GetMean()


print("weight_pu: ", mean1)
print("weight_pu_up: ", mean2)
print("weight_pu_down: ", mean3)


# # Don't forget to close the file!
# file.Close()

# def calculate_histogram_mean(file_path, tree_name, histogram_branch_name):
#     # Open the ROOT file
#     file = ROOT.TFile.Open(file_path)

#     # Access the TTree
#     tree = file.Get(tree_name)

#     # Set the branch address to the histogram object
#     histogram = ROOT.TH1F()
#     tree.SetBranchAddress(histogram_branch_name, histogram)

#     mean = histogram.GetMean()

#     return mean

# Usage example
# file_path = "uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_0.root"
# tree_name = "AnalysisTree"
# histogram_branch_name1 = "weight_pu"
# histogram_branch_name2 = "weight_pu_up"
# histogram_branch_name3 = "weight_pu_down"

# mean1 = calculate_histogram_mean(file_path, tree_name, histogram_branch_name1)
# mean2 = calculate_histogram_mean(file_path, tree_name, histogram_branch_name2)
# mean3 = calculate_histogram_mean(file_path, tree_name, histogram_branch_name3)
# print("weight_pu:", mean1)
# print("weight_pu_up:", mean2)
# print("weight_pu_down:", mean3)
