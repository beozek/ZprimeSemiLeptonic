
//  importing libraries

#include <algorithm>
#include <iterator>
#include <TROOT.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLatex.h>
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH1D.h"
#include "THStack.h"
#include "TRandom.h"
#include "TUnfoldDensity.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFrame.h"
#include "TPaveLabel.h"
#include "TPad.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TMath.h"

using namespace RooFit ;

void nom_sys_Diboson()

{

    gStyle->SetOptStat(0);

    // A chain is a collection of files containing TTree objects. 
    // TChain(const char *name, const char *title="", Mode mode=kWithGlobalRegistration or kWithoutGlobalReg)
    // TTree tree(name, title)

    TChain *reco = new TChain("AnalysisTree","");

    // reco-> Add("/nfs/dust/cms/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/output_combine/UL18/muon/workdir_Zprime_Analysis_UL18_muon_combine_ttbar/nominal/uhh2.AnalysisModuleRunner.ttbar1.root");
    // reco-> Add("/nfs/dust/cms/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/output_combine/UL18/muon/workdir_Zprime_Analysis_UL18_muon_combine_ttbar/nominal/uhh2.AnalysisModuleRunner.ttbar2.root");
    // reco-> Add("/nfs/dust/cms/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/output_combine/UL18/muon/workdir_Zprime_Analysis_UL18_muon_combine_ttbar/nominal/uhh2.AnalysisModuleRunner.ttbar3.root");
    // reco-> Add("/nfs/dust/cms/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/output_combine/UL18/muon/workdir_Zprime_Analysis_UL18_muon_combine_ttbar/nominal/uhh2.AnalysisModuleRunner.ttbar4.root");

    // reco-> Add("/nfs/dust/cms/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/output_combine/UL18/muon/workdir_Zprime_Analysis_UL18_muon_combine_latest/DATA/nominal/uhh2.AnalysisModuleRunner.DATA.root");
    // reco-> Add("/nfs/dust/cms/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/output_combine/UL18/muon/workdir_Zprime_Analysis_UL18_muon_combine_latest/WJets/nominal/uhh2.AnalysisModuleRunner.WJets.root");
    // reco-> Add("/nfs/dust/cms/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/output_combine/UL18/muon/workdir_Zprime_Analysis_UL18_muon_combine_latest/DY/nominal/uhh2.AnalysisModuleRunner.DY.root");
    reco-> Add("/nfs/dust/cms/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/output_combine/UL18/muon/workdir_Zprime_Analysis_UL18_muon_combine_latest/Diboson/nominal/uhh2.AnalysisModuleRunner.Diboson.root");
    // reco-> Add("/nfs/dust/cms/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/output_combine/UL18/muon/workdir_Zprime_Analysis_UL18_muon_combine_latest/ST/nominal/uhh2.AnalysisModuleRunner.ST.root");
    // reco-> Add("/nfs/dust/cms/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/output_combine/UL18/muon/workdir_Zprime_Analysis_UL18_muon_combine_latest/QCD/nominal/uhh2.AnalysisModuleRunner.QCD.root");

    TTree *treereco = (TTree*) reco;

    cout << "Number of Events:"<< treereco-> GetEntries()<<endl;

    
    // TH1D DeltaY Plots

    //DeltaY gen without mass cut
    TH1D *h_DeltaY_gen = new TH1D("DeltaY_gen","#Delta_Y_{gen}",10,-2.5,2.5);
    //DeltaY gen with mass cut
    TH1D *h_DeltaY_gen_mass = new TH1D("DeltaY_gen_mass","(#Delta_Y)_{gen}, M > 750",10,-2.5,2.5);
    //POSITIVE gen without mass
    TH1D *h_DeltaY_P_gen_nomass = new TH1D("DeltaY_P_gen_nomass","#Delta_Y_{gen}>0",1,0,2.5);
    //POSITIVE gen with mass
    TH1D *h_DeltaY_P_gen = new TH1D("DeltaY_P_gen","#Delta_Y_{gen}>0, M > 750",1,0,2.5);
    //NEGATIVE gen without mass
    TH1D *h_DeltaY_N_gen_nomass = new TH1D("DeltaY_N_gen_nomass","#Delta_Y_{gen} < 0",1,-2.5,0);
    //NEGATIVE gen with mass
    TH1D *h_DeltaY_N_gen = new TH1D("DeltaY_N_gen","(#Delta_Y_{gen} < 0, M > 750",1,-2.5,0);

    //DeltaY reco without mass cut
    TH1D *h_DeltaY_reco = new TH1D("DeltaY_reco","#Delta_Y_{reco}",2,-2.5,2.5);
    //DeltaY with mass cut
    TH1D *h_DeltaY_reco_mass = new TH1D("DeltaY_reco_mass","(#Delta_Y)_{reco}, M > 750",2,-2.5,2.5);
    //POSITIVE reco without mass
    TH1D *h_DeltaY_P_reco_nomass = new TH1D("DeltaY_P_reco_nomass","#Delta_Y_{reco}>0",1,0,2.5);
    //POSITIVE reco with mass
    TH1D *h_DeltaY_P_reco = new TH1D("DeltaY_P_reco","#Delta_Y_{reco}>0, M>750",1,0,2.5);
    //NEGATIVE reco without mass
    TH1D *h_DeltaY_N_reco_nomass = new TH1D("DeltaY_N_reco_nomass","#Delta_Y_{reco}<0",1,-2.5,0);
    //NEGATIVE reco with mass
    TH1D *h_DeltaY_N_reco = new TH1D("DeltaY_N_reco","#Delta_Y_{reco}<0, M>750",1,-2.5,0);


    // POSITIVE gen, POSITIVE reco, without mass cut
    TH1D *h_DeltaY_P_P_nomass = new TH1D("DeltaY_P_P_nomass","#Delta_Y_{gen} > 0, #Delta_Y_{reco} > 0 ",1,0,2.5);
    // POSITIVE gen, POSITIVE reco, with mass cut
    TH1D *h_DeltaY_P_P = new TH1D("DeltaY_P_P","#Delta_Y_{gen} > 0, #Delta_Y_{reco} > 0, M >750",1,0,2.5);
    // POSITIVE gen, NEGATIVE reco, without mass cut
    TH1D *h_DeltaY_P_N_nomass = new TH1D("DeltaY_P_N_nomass","#Delta_Y_{gen} > 0, #Delta_Y_{reco} < 0",1,-2.5,2.5);
    // POSITIVE gen, NEGATIVE reco, with mass cut
    TH1D *h_DeltaY_P_N = new TH1D("DeltaY_P_N","#Delta_Y_{gen} > 0, #Delta_Y_{reco} < 0, M >750",1,-2.5,2.5);
    // NEAGATIVE gen, POSITIVE reco, without mass cut
    TH1D *h_DeltaY_N_P_nomass = new TH1D("DeltaY_N_P_nomass","#Delta_Y_{gen} < 0, #Delta_Y_{reco} > 0",1,-2.5,2.5);
    // NEGATIVE gen, POSITIVE reco, with mass cut
    TH1D *h_DeltaY_N_P = new TH1D("DeltaY_N_P","#Delta_Y_{gen} < 0, #Delta_Y_{reco} > 0, M >750",1,-2.5,2.5);
    // NEGATIVE gen, NEGATIVE reco, without mass cut
    TH1D *h_DeltaY_N_N_nomass = new TH1D("DeltaY_N_N_nomass","#Delta_Y_{gen} < 0, #Delta_Y_{reco} < 0",1,-2.5,0);
    // NEGATIVE gen, NEGATIVE reco, with mass cut
    TH1D *h_DeltaY_N_N = new TH1D("DeltaY_N_N","#Delta_Y_{gen} < 0, #Delta_Y_{reco} < 0, M >750",1,-2.5,0);
    

    // TH1D Projection Plots

    TH1D *ProjY_1 = new TH1D("ProjY_1","Project along Y , #Delta_Y_{reco} < 0 ",2,-2.5,2.5);
    TH1D *ProjY_2 = new TH1D("ProjY_2","Project along Y , #Delta_Y_{reco} > 0 ",2,-2.5,2.5);

    TH1D *ProjX_1 = new TH1D("ProjX_1","Project along X , #Delta_Y_{gen} < 0 ",2,-2.5,2.5);
    TH1D *ProjX_2 = new TH1D("ProjX_2","Project along X ,#Delta_Y_{gen} > 0 ",2,-2.5,2.5);

   // TH1F Pileup 
    
    TH1D *h_weight_pu = new TH1D("weight_pu","weight_pu",1,-1,1);
    TH1D *h_weight_pu_down = new TH1D("weight_pu_down","weight_pu_down",1,-1,1);
    TH1D *h_weight_pu_up = new TH1D("weight_pu_up","weight_pu_up",1,-1,1);

    // TH1F Muon ID 
    
    TH1D *h_weight_sfmu_id = new TH1D("weight_sfmu_id","weight_sfmu_id",1,-1,1);
    TH1D *h_weight_sfmu_id_down = new TH1D("weight_sfmu_id_down","weight_sfmu_id_down",1,-1,1);
    TH1D *h_weight_sfmu_id_up = new TH1D("weight_sfmu_id_up","weight_sfmu_id_up",1,-1,1);

    // TH1F Muon Reconstruction

    TH1D *h_muonrecSF_nominal = new TH1D("muonrecSF_nominal","muonrecSF_nominal",1,-1,1);
    TH1D *h_muonrecSF_down = new TH1D("muonrecSF_down","muonrecSF_down",1,-1,1);
    TH1D *h_muonrecSF_up = new TH1D("muonrecSF_up","muonrecSF_up",1,-1,1);

    
    // TH1F Muon Trigger
    TH1D *h_weight_sfmu_trigger = new TH1D("weight_sfmu_trigger","weight_sfmu_trigger",1,-1,1);
    TH1D *h_weight_sfmu_trigger_down = new TH1D("weight_sfmu_trigger_down","weight_sfmu_trigger_down",1,-1,1);
    TH1D *h_weight_sfmu_trigger_up = new TH1D("weight_sfmu_trigger_up","weight_sfmu_trigger_up",1,-1,1);

    
    // TH1F Muon Isolation
    TH1D *h_weight_sfmu_iso = new TH1D("weight_sfmu_iso","weight_sfmu_iso",1,-1,1);
    TH1D *h_weight_sfmu_iso_down = new TH1D("weight_sfmu_iso_down","weight_sfmu_iso_down",1,-1,1);
    TH1D *h_weight_sfmu_iso_up = new TH1D("weight_sfmu_iso_up","weight_sfmu_iso_up",1,-1,1);

    
    // TH1F Electron Reconstruction

    TH1D *h_weight_sfelec_reco = new TH1D("weight_sfelec_reco","weight_sfelec_reco",1,-1,1);
    TH1D *h_weight_sfelec_reco_down = new TH1D("weight_sfelec_reco_down","weight_sfelec_reco_down",1,-1,1);
    TH1D *h_weight_sfelec_reco_up = new TH1D("weight_sfelec_reco_up","weight_sfelec_reco_up",1,-1,1);

   
    // TH1F Electron ID+isolation
    TH1D *h_weight_sfelec_id = new TH1D("weight_sfelec_id","weight_sfelec_id",1,-1,1);
    TH1D *h_weight_sfelec_id_down = new TH1D("weight_sfelec_id_down","weight_sfelec_id_down",1,-1,1);
    TH1D *h_weight_sfelec_id_up = new TH1D("weight_sfelec_id_up","weight_sfelec_id_up",1,-1,1);

    //Trigger Prefiring

    TH1D *h_prefiringWeight = new TH1D("prefiringWeight","prefiringWeight",1,-1,1);
    TH1D *h_prefiringWeightDown = new TH1D("prefiringWeightDown","prefiringWeightDown",1,-1,1);
    TH1D *h_prefiringWeightUp = new TH1D("prefiringWeightUp","prefiringWeightUp",1,-1,1);


    // TH2D Matrix 
    TH2D *Matrix = new TH2D("Matrix","", 2,-2.5,2.5,2,-2.5,2.5);

    float DeltaY_gen;
    float DeltaY_gen_mass;
    float DeltaY_P_gen;
    float DeltaY_P_gen_nomass;
    float DeltaY_N_gen;
    float DeltaY_N_gen_nomass;
    
    float DeltaY_reco;
    float DeltaY_reco_mass;
    float DeltaY_P_reco;
    float DeltaY_P_reco_nomass;
    float DeltaY_N_reco;
    float DeltaY_N_reco_nomass;

    float DeltaY_N_N;
    float DeltaY_N_P;
    float DeltaY_P_N;
    float DeltaY_P_P;
    float DeltaY_N_N_nomass;
    float DeltaY_N_P_nomass;
    float DeltaY_P_N_nomass;
    float DeltaY_P_P_nomass;

    float weight_pu;
    float weight_pu_down;
    float weight_pu_up;
    
    float weight_sfmu_id;
    float weight_sfmu_id_down;
    float weight_sfmu_id_up;

    float muonrecSF_nominal;
    float muonrecSF_down;
    float muonrecSF_up;
    
    float weight_sfmu_trigger;
    float weight_sfmu_trigger_down;
    float weight_sfmu_trigger_up;

    float weight_sfmu_iso;
    float weight_sfmu_iso_down;
    float weight_sfmu_iso_up;

    float weight_sfelec_reco;
    float weight_sfelec_reco_down;
    float weight_sfelec_reco_up;

    float weight_sfelec_id;
    float weight_sfelec_id_down;
    float weight_sfelec_id_up;

    float prefiringWeight;
    float prefiringWeightDown;
    float prefiringWeightUp;
    
    treereco->SetBranchAddress("DeltaY_gen", &DeltaY_gen);
    treereco->SetBranchAddress("DeltaY_gen_mass", &DeltaY_gen_mass);
    treereco->SetBranchAddress("DeltaY_P_gen", &DeltaY_P_gen);
    treereco->SetBranchAddress("DeltaY_P_gen_nomass", &DeltaY_P_gen_nomass);
    treereco->SetBranchAddress("DeltaY_N_gen", &DeltaY_N_gen);
    treereco->SetBranchAddress("DeltaY_N_gen_nomass", &DeltaY_N_gen_nomass);

    treereco->SetBranchAddress("DeltaY_reco", &DeltaY_reco);
    treereco->SetBranchAddress("DeltaY_reco_mass", &DeltaY_reco_mass);
    treereco->SetBranchAddress("DeltaY_P_reco", &DeltaY_P_reco);
    treereco->SetBranchAddress("DeltaY_P_reco_nomass", &DeltaY_P_reco_nomass);
    treereco->SetBranchAddress("DeltaY_N_reco", &DeltaY_N_reco);
    treereco->SetBranchAddress("DeltaY_N_reco_nomass", &DeltaY_N_reco_nomass);

    treereco->SetBranchAddress("DeltaY_N_N", &DeltaY_N_N);
    treereco->SetBranchAddress("DeltaY_N_P", &DeltaY_N_P);
    treereco->SetBranchAddress("DeltaY_P_P", &DeltaY_P_P);
    treereco->SetBranchAddress("DeltaY_P_N", &DeltaY_P_N);
    treereco->SetBranchAddress("DeltaY_N_N_nomass", &DeltaY_N_N_nomass);
    treereco->SetBranchAddress("DeltaY_N_P_nomass", &DeltaY_N_P_nomass);
    treereco->SetBranchAddress("DeltaY_P_N_nomass", &DeltaY_P_N_nomass);
    treereco->SetBranchAddress("DeltaY_P_P_nomass", &DeltaY_P_P_nomass);

    treereco->SetBranchAddress("weight_pu", &weight_pu);
    treereco->SetBranchAddress("weight_pu_down", &weight_pu_down);
    treereco->SetBranchAddress("weight_pu_up", &weight_pu_up);

    treereco->SetBranchAddress("weight_sfmu_id", &weight_sfmu_id);
    treereco->SetBranchAddress("weight_sfmu_id_down", &weight_sfmu_id_down);
    treereco->SetBranchAddress("weight_sfmu_id_up", &weight_sfmu_id_up);

    treereco->SetBranchAddress("muonrecSF_nominal", &muonrecSF_nominal);
    treereco->SetBranchAddress("muonrecSF_down", &muonrecSF_down);
    treereco->SetBranchAddress("muonrecSF_up", &muonrecSF_up);

    treereco->SetBranchAddress("weight_sfmu_trigger", &weight_sfmu_trigger);
    treereco->SetBranchAddress("weight_sfmu_trigger_down", &weight_sfmu_trigger_down);
    treereco->SetBranchAddress("weight_sfmu_trigger_up", &weight_sfmu_trigger_up);

    treereco->SetBranchAddress("weight_sfmu_iso", &weight_sfmu_iso);
    treereco->SetBranchAddress("weight_sfmu_iso_down", &weight_sfmu_iso_down);
    treereco->SetBranchAddress("weight_sfmu_iso_up", &weight_sfmu_iso_up);

    treereco->SetBranchAddress("weight_sfelec_reco", &weight_sfelec_reco);
    treereco->SetBranchAddress("weight_sfelec_reco_down", &weight_sfelec_reco_down);
    treereco->SetBranchAddress("weight_sfelec_reco_up", &weight_sfelec_reco_up);

    treereco->SetBranchAddress("weight_sfelec_id", &weight_sfelec_id);
    treereco->SetBranchAddress("weight_sfelec_id_down", &weight_sfelec_id_down);
    treereco->SetBranchAddress("weight_sfelec_id_up", &weight_sfelec_id_up);

    treereco->SetBranchAddress("prefiringWeight", &prefiringWeight);
    treereco->SetBranchAddress("prefiringWeightDown", &prefiringWeightDown);
    treereco->SetBranchAddress("prefiringWeightUp", &prefiringWeightUp);

    for (Int_t i = 0; i < treereco->GetEntries(); i++){
    // for (Int_t i = 0; i < 10000; i++){

        treereco->GetEntry(i);
        if (i%1000000 == 0) std::cout << "--- ... Processing event: " << i <<std::endl;
       
        h_DeltaY_gen->Fill(DeltaY_gen);
        h_DeltaY_gen_mass->Fill(DeltaY_gen_mass);
        h_DeltaY_P_gen_nomass->Fill(DeltaY_P_gen_nomass);
        h_DeltaY_P_gen->Fill(DeltaY_P_gen);
        h_DeltaY_N_gen_nomass->Fill(DeltaY_N_gen_nomass);
        h_DeltaY_N_gen->Fill(DeltaY_N_gen);

        h_DeltaY_reco->Fill(DeltaY_reco);
        h_DeltaY_reco_mass->Fill(DeltaY_reco_mass);
        h_DeltaY_P_reco_nomass->Fill(DeltaY_P_reco_nomass);
        h_DeltaY_P_reco->Fill(DeltaY_P_reco);
        h_DeltaY_N_reco_nomass->Fill(DeltaY_N_reco_nomass);
        h_DeltaY_N_reco->Fill(DeltaY_N_reco);

        h_DeltaY_P_P->Fill(DeltaY_P_P);
        h_DeltaY_P_N->Fill(DeltaY_P_N);
        h_DeltaY_N_N->Fill(DeltaY_N_N);
        h_DeltaY_N_P->Fill(DeltaY_N_P);
        h_DeltaY_P_P_nomass->Fill(DeltaY_P_P_nomass);
        h_DeltaY_P_N_nomass->Fill(DeltaY_P_N_nomass);
        h_DeltaY_N_P_nomass->Fill(DeltaY_N_P_nomass);
        h_DeltaY_N_N_nomass->Fill(DeltaY_N_N_nomass);

        h_weight_pu->Fill(weight_pu);
        h_weight_pu_down->Fill(weight_pu_down);
        h_weight_pu_up->Fill(weight_pu_up);
        
        h_weight_sfmu_id->Fill(weight_sfmu_id);
        h_weight_sfmu_id_down->Fill(weight_sfmu_id_down);
        h_weight_sfmu_id_up->Fill(weight_sfmu_id_up);

        h_muonrecSF_nominal->Fill(muonrecSF_nominal);
        h_muonrecSF_down->Fill(muonrecSF_down); 
        h_muonrecSF_up->Fill(muonrecSF_up);

        h_weight_sfmu_trigger->Fill(weight_sfmu_trigger); 
        h_weight_sfmu_trigger_down->Fill(weight_sfmu_trigger_down); 
        h_weight_sfmu_trigger_up->Fill(weight_sfmu_trigger_up); 

        h_weight_sfmu_iso->Fill(weight_sfmu_iso); 
        h_weight_sfmu_iso_down->Fill(weight_sfmu_iso_down); 
        h_weight_sfmu_iso_up->Fill(weight_sfmu_iso_up);

        h_weight_sfelec_reco->Fill(weight_sfelec_reco);
        h_weight_sfelec_reco_down->Fill(weight_sfelec_reco_down);
        h_weight_sfelec_reco_up->Fill(weight_sfelec_reco_up); 

        h_weight_sfelec_id->Fill(weight_sfelec_id);
        h_weight_sfelec_id_down->Fill(weight_sfelec_id_down); 
        h_weight_sfelec_id_up->Fill(weight_sfelec_id_up); 

        h_prefiringWeight->Fill(prefiringWeight);
        h_prefiringWeightDown->Fill(prefiringWeightDown); 
        h_prefiringWeightUp->Fill(prefiringWeightUp); 
    }

    double integral [2][2] = {{h_DeltaY_N_N->Integral(),h_DeltaY_P_N->Integral()},{h_DeltaY_N_P->Integral(),h_DeltaY_P_P->Integral()}};

     for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
              Matrix->SetBinContent(i+1,j+1,integral[i][j]);
       }
    }

    
  
    ProjY_1->GetXaxis()->SetTitle("#Delta_Y_{gen}");
    ProjY_2->GetXaxis()->SetTitle("#Delta_Y_{gen}");
    ProjX_1->GetXaxis()->SetTitle("#Delta_Y_{reco}");
    ProjX_2->GetXaxis()->SetTitle("#Delta_Y_{reco}");

    ProjY_1 = Matrix->ProjectionY("py1",1,1);
    ProjY_2 = Matrix->ProjectionY("py2",2,2);

    ProjX_1 = Matrix->ProjectionX("px1",1,1);
    ProjX_2 = Matrix->ProjectionX("px2", 2,2);


  // --------------- Output File ------------------

    TFile* myFile = new TFile("dY_UL18_muon_750_900_Diboson.root", "RECREATE");
    
    h_DeltaY_gen->Write();
    h_DeltaY_gen_mass->Write();
    h_DeltaY_P_gen_nomass->Write();
    h_DeltaY_P_gen->Write();
    h_DeltaY_N_gen_nomass->Write();
    h_DeltaY_N_gen->Write();

    h_DeltaY_reco->Write();
    h_DeltaY_reco_mass->Write();
    h_DeltaY_P_reco_nomass->Write();
    h_DeltaY_P_reco->Write();
    h_DeltaY_N_reco_nomass->Write();
    h_DeltaY_N_reco->Write();  
    
    h_DeltaY_P_P->Write();
    h_DeltaY_P_N->Write();
    h_DeltaY_N_P->Write();
    h_DeltaY_N_N->Write();
    h_DeltaY_P_P_nomass->Write();
    h_DeltaY_P_N_nomass->Write();
    h_DeltaY_N_P_nomass->Write();
    h_DeltaY_N_N_nomass->Write();

    Matrix->Write();

    ProjY_1->Write();
    ProjY_2->Write();
    ProjX_1->Write();
    ProjX_2->Write();

    h_weight_pu->Write();
    h_weight_pu_down->Write();
    h_weight_pu_up->Write();

    h_weight_sfmu_id->Write();
    h_weight_sfmu_id_down->Write();
    h_weight_sfmu_id_up->Write();

    h_muonrecSF_nominal->Write();
    h_muonrecSF_down->Write();
    h_muonrecSF_up->Write();

    h_weight_sfmu_trigger->Write();
    h_weight_sfmu_trigger_down->Write(); 
    h_weight_sfmu_trigger_up->Write();

    h_weight_sfmu_iso->Write();
    h_weight_sfmu_iso_down->Write();
    h_weight_sfmu_iso_up->Write();

    h_weight_sfelec_reco->Write();
    h_weight_sfelec_reco_down->Write();
    h_weight_sfelec_reco_up->Write();

    h_weight_sfelec_id->Write();
    h_weight_sfelec_id_down->Write();
    h_weight_sfelec_id_up->Write();

    h_prefiringWeight->Write();
    h_prefiringWeightDown->Write();
    h_prefiringWeightUp->Write();


myFile->Close();
}