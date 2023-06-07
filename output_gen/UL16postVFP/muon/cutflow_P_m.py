import ROOT

# Open the ROOT file and access the TTree
file = ROOT.TFile.Open("/nfs/dust/cms/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/output_gen/UL16postVFP/muon/workdir_Analysis_UL16postVFP_muon_gen/nominal/uhh2.AnalysisModuleRunner.Semileptonic1.root")
tree = file.Get("AnalysisTree")

# Create histograms for each branch
hist1 = ROOT.TH1F("hist1", "DeltaY_P_gen_muon", 1, -2, 2)
hist2 = ROOT.TH1F("hist2", "DeltaY_P_gen_pt_muon", 1, -2, 2)
hist3 = ROOT.TH1F("hist3", "DeltaY_P_gen_eta_muon", 1, -2, 2)
hist4 = ROOT.TH1F("hist4", "DeltaY_P_gen_jet_muon", 1, -2, 2)

# Fill the histograms
for i, event in enumerate(tree):
    # Load the current tree entry
    hist1.Fill(event.DeltaY_P_gen_muon)
    hist2.Fill(event.DeltaY_P_gen_pt_muon)
    hist3.Fill(event.DeltaY_P_gen_eta_muon)
    hist4.Fill(event.DeltaY_P_gen_jet_muon)

        
    if i % 10000 == 0 and i > 0: 
        print("Processing entry {}/{}".format(i, tree.GetEntries()))

print("UL16postVFP muon -> Number of events without cuts: {}".format(hist1.Integral()))
print("UL16postVFP muon -> Number of events with lepton_pt cut: {}".format(hist2.Integral()))
print("UL16postVFP muon -> Number of events with lepton_pt & eta cut : {}".format(hist3.Integral()))
print("UL16postVFP muon -> Number of events with lepton_pt & eta & jet pteta cut : {}".format(hist4.Integral()))

# Create a new histogram with 3 bins
summary_hist = ROOT.TH1F("summary_hist", "UL16postVFP muon DeltaY > 0", 4, 0, 4)
summary_hist.GetXaxis().SetBinLabel(1, "All")
summary_hist.GetXaxis().SetBinLabel(2, "Leptopn Pt")
summary_hist.GetXaxis().SetBinLabel(3, "Lepton Pt&Eta")
summary_hist.GetXaxis().SetBinLabel(4, "Lepton Pt&Eta and Jet Pt&Eta")

# Fill the new histogram with the integral of each branch histogram
summary_hist.SetBinContent(1, hist1.Integral())
summary_hist.SetBinContent(2, hist2.Integral())
summary_hist.SetBinContent(3, hist3.Integral())
summary_hist.SetBinContent(4, hist4.Integral())

# Draw the new histogram
canvas = ROOT.TCanvas("canvas", "UL16postVFP muon DeltaY > 0", 800, 600)
summary_hist.Draw()

canvas.Draw()
canvas.SaveAs("DeltaY_P_gen_muon_new.png")



