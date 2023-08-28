#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

void ttree() {

    // Open the ROOT file
    TFile* file = TFile::Open("semi.root");
    if (!file) {
        std::cerr << "Failed to open the file." << std::endl;
        return;
    }

    // Retrieve the TTree
    TTree* tree = (TTree*)file->Get("AnalysisTree");
    if (!tree) {
        std::cerr << "Failed to retrieve the TTree." << std::endl;
        return;
    }

    // Histogram pointer
    Float_t not_reconstructed_900Inf;
    tree->SetBranchAddress("not_reconstructed_900Inf", &not_reconstructed_900Inf);

    int count = 0;
    for(int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        if (not_reconstructed_900Inf > 0) {
            count++;
        }
    }

    std::cout << "Number of events where not_reconstructed_900Inf is greater than 0: " << count << std::endl;

    file->Close();
}
