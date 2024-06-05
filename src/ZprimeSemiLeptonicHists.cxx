#include "UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicHists.h"
#include "UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicModules.h"
#include "UHH2/core/include/Event.h"
#include <UHH2/core/include/Utils.h>
#include <UHH2/common/include/Utils.h>
#include "UHH2/common/include/JetIds.h"
#include <math.h>

#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/TTbarReconstruction.h>
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>

#include <UHH2/core/include/LorentzVector.h>

#include "TH1F.h"
#include "TH2D.h"
#include <iostream>



using namespace std;
using namespace uhh2;

ZprimeSemiLeptonicHists::ZprimeSemiLeptonicHists(uhh2::Context& ctx, const std::string& dirname):
Hists(ctx, dirname) {
  is_mc = ctx.get("dataset_type") == "MC";
  ishotvr = (ctx.get("is_hotvr") == "true");
  isdeepAK8 = (ctx.get("is_deepAK8") == "true");
 
  isUL16preVFP=false; isUL16postVFP=false; isUL17=false; isUL18 =false;

  isUL16preVFP  = (ctx.get("dataset_version").find("UL16preVFP")  != std::string::npos);
  isUL16postVFP = (ctx.get("dataset_version").find("UL16postVFP") != std::string::npos);
  isUL17        = (ctx.get("dataset_version").find("UL17")        != std::string::npos);
  isUL18        = (ctx.get("dataset_version").find("UL18")        != std::string::npos);

  debug = false;
  NN = false;
  isMuon = false; isElectron = false;
  if(ctx.get("channel") == "muon") isMuon = true;
  if(ctx.get("channel") == "electron") isElectron = true;
  is_tt = ctx.get("dataset_version").find("TTTo") == 0;
  if(isdeepAK8){
    h_AK8TopTags = ctx.get_handle<std::vector<TopJet>>("DeepAK8TopTags");
  }else if(ishotvr){
    h_AK8TopTags = ctx.get_handle<std::vector<TopJet>>("HOTVRTopTags");
  }
  h_CHSjets_matched = ctx.get_handle<std::vector<Jet>>("CHS_matched");

  h_BestZprimeCandidateChi2 = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestChi2");
  h_BestZprimeCandidateCorrectMatch = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestCorrectMatch");
  h_is_zprime_reconstructed_chi2 = ctx.get_handle<bool>("is_zprime_reconstructed_chi2");
  h_is_zprime_reconstructed_correctmatch = ctx.get_handle<bool>("is_zprime_reconstructed_correctmatch");
  //  h_chi2 = ctx.get_handle<float>("chi2");
  init();
}

void ZprimeSemiLeptonicHists::init(){

  // jets
  N_jets   = book<TH1F>("N_jets", "N_{jets}", 21, -0.5, 20.5);
  pt_jet   = book<TH1F>("pt_jet", "p_{T}^{jets} [GeV]", 45, 0, 900);
  pt_jet1  = book<TH1F>("pt_jet1", "p_{T}^{jet 1} [GeV]", 45, 0, 900);
  pt_jet2  = book<TH1F>("pt_jet2", "p_{T}^{jet 2} [GeV]", 45, 0, 900);
  pt_jet3  = book<TH1F>("pt_jet3", "p_{T}^{jet 3} [GeV]", 45, 0, 900);
  eta_jet  = book<TH1F>("eta_jet", "#eta^{jets}", 50, -2.5, 2.5);
  eta_jet1 = book<TH1F>("eta_jet1", "#eta^{jet 1}", 50, -2.5, 2.5);
  eta_jet2 = book<TH1F>("eta_jet2", "#eta^{jet 2}", 50, -2.5, 2.5);
  eta_jet3 = book<TH1F>("eta_jet3", "#eta^{jet 3}", 50, -2.5, 2.5);
  phi_jet  = book<TH1F>("phi_jet", "#phi^{jets}", 35, -3.5, 3.5);
  phi_jet1 = book<TH1F>("phi_jet1", "#phi^{jet 1}", 35, -3.5, 3.5);
  phi_jet2 = book<TH1F>("phi_jet2", "#phi^{jet 2}", 35, -3.5, 3.5);
  phi_jet3 = book<TH1F>("phi_jet3", "#phi^{jet 3}", 35, -3.5, 3.5);
  m_jet    = book<TH1F>("m_jet", "m^{jets}", 50, 0, 500);
  m_jet1   = book<TH1F>("m_jet1", "m^{jet 1}", 50, 0, 500);
  m_jet2   = book<TH1F>("m_jet2", "m^{jet 2}", 50, 0, 500);
  m_jet3   = book<TH1F>("m_jet3", "m^{jet 3}", 50, 0, 500);

  deepjetbscore_jet    = book<TH1F>("deepjetbscore_jet", "DeepJet b-tag score all AK4 jets", 20, 0, 1);
  deepjetbscore_jet1   = book<TH1F>("deepjetbscore_jet1", "DeepJet b-tag score AK4 jet 1", 20, 0, 1);
  deepjetbscore_jet2   = book<TH1F>("deepjetbscore_jet2", "DeepJet b-tag score AK4 jet 2", 20, 0, 1);
  deepjetbscore_jet3   = book<TH1F>("deepjetbscore_jet3", "DeepJet b-tag score AK4 jet 3", 20, 0, 1);
  N_bJetsDeepJet_loose = book<TH1F>("N_bJetsDeepJet_loose", "N_{jets}^{DeepJet loose}", 11, -0.5, 10.5);
  N_bJetsDeepJet_med   = book<TH1F>("N_bJetsDeepJet_med", "N_{jets}^{DeepJet medium}", 11, -0.5, 10.5);
  N_bJetsDeepJet_tight = book<TH1F>("N_bJetsDeepJet_tight", "N_{jets}^{DeepJet tight}", 11, -0.5, 10.5);

  // leptons
  N_mu             = book<TH1F>("N_mu", "N^{#mu}", 11, -0.5, 10.5);
  pt_mu            = book<TH1F>("pt_mu", "p_{T}^{#mu} [GeV]", 90, 0, 900);
  pt_mu1           = book<TH1F>("pt_mu1", "p_{T}^{#mu 1} [GeV]", 90, 0, 900);
  pt_mu2           = book<TH1F>("pt_mu2", "p_{T}^{#mu 2} [GeV]", 90, 0, 900);
  eta_mu           = book<TH1F>("eta_mu", "#eta^{#mu}", 50, -2.5, 2.5);
  eta_mu1          = book<TH1F>("eta_mu1", "#eta^{#mu 1}", 50, -2.5, 2.5);
  eta_mu2          = book<TH1F>("eta_mu2", "#eta^{#mu 2}", 50, -2.5, 2.5);
  phi_mu           = book<TH1F>("phi_mu", "#phi^{#mu}", 35, -3.5, 3.5);
  phi_mu1          = book<TH1F>("phi_mu1", "#phi^{#mu 1}", 35, 3.5, 3.5);
  phi_mu2          = book<TH1F>("phi_mu2", "#phi^{#mu 2}", 35, 3.5, 3.5);
  reliso_mu        = book<TH1F>("reliso_mu", "#mu rel. Iso", 40, 0, 0.5);
  reliso_mu1       = book<TH1F>("reliso_mu1", "#mu 1 rel. Iso", 40, 0, 0.5);
  reliso_mu2       = book<TH1F>("reliso_mu2", "#mu 2 rel. Iso", 40, 0, 0.5);
  reliso_mu_rebin  = book<TH1F>("reliso_mu_rebin", "#mu rel. Iso ", 400, 0, 5);
  reliso_mu1_rebin = book<TH1F>("reliso_mu1_rebin", "#mu 1 rel. Iso ", 400, 0, 5);
  reliso_mu2_rebin = book<TH1F>("reliso_mu2_rebin", "#mu 2 rel. Iso ", 400, 0, 5);

  N_mu_charge      = book<TH1F>("N_mu_charge", "N^{#mu} charge", 2, -1., 1.);
  N_lep_charge     = book<TH1F>("N_lep_charge", "lepton charge", 2, -2., 2.);
  N_ele_charge     = book<TH1F>("N_ele_charge", "N^{e} charge", 2, -1., 1.);
  


  N_ele             = book<TH1F>("N_ele", "N^{e}", 11, -0.5, 10.5);
  pt_ele            = book<TH1F>("pt_ele", "p_{T}^{e} [GeV]", 90, 0, 900);
  pt_ele1           = book<TH1F>("pt_ele1", "p_{T}^{e 1} [GeV]", 90, 0, 900);
  pt_ele2           = book<TH1F>("pt_ele2", "p_{T}^{e 2} [GeV]", 90, 0, 900);
  eta_ele           = book<TH1F>("eta_ele", "#eta^{e}", 50, -2.5, 2.5);
  eta_ele1          = book<TH1F>("eta_ele1", "#eta^{ele 1}", 50, -2.5, 2.5);
  eta_ele2          = book<TH1F>("eta_ele2", "#eta^{ele 2}", 50, -2.5, 2.5);
  phi_ele           = book<TH1F>("phi_ele", "#phi^{e}", 35, -3.5, 3.5);
  phi_ele1          = book<TH1F>("phi_ele1", "#phi^{e 1}", 35, -3.5, 3.5);
  phi_ele2          = book<TH1F>("phi_ele2", "#phi^{e 2}", 35, -3.5, 3.5);
  reliso_ele        = book<TH1F>("reliso_ele", "e rel. Iso", 40, 0, 0.5);
  reliso_ele1       = book<TH1F>("reliso_ele1", "e 1 rel. Iso", 40, 0, 0.5);
  reliso_ele2       = book<TH1F>("reliso_ele2", "e 2 rel. Iso", 40, 0, 0.5);
  reliso_ele_rebin  = book<TH1F>("reliso_ele_rebin", "e rel. Iso", 400, 0, 5);
  reliso_ele1_rebin = book<TH1F>("reliso_ele1_rebin", "e 1 rel. Iso", 400, 0, 5);
  reliso_ele2_rebin = book<TH1F>("reliso_ele2_rebin", "e 2 rel. Iso", 400, 0, 5);
  pt_ele_lowpt      = book<TH1F>("pt_ele_lowpt", "p_{T}^{e} [GeV] p_{T}<120 GeV", 180, 0, 900);
  pt_ele_midpt      = book<TH1F>("pt_ele_midpt", "p_{T}^{e} [GeV] 120<p_{T}<200 GeV", 180, 0, 900);
  pt_ele_highpt     = book<TH1F>("pt_ele_highpt", "p_{T}^{e} [GeV] p_{T}>200 GeV", 180, 0, 900);
  eta_ele_lowpt     = book<TH1F>("eta_ele_lowpt", "#eta^{e} p_{T}<120 GeV", 50, -2.5, 2.5);
  eta_ele_midpt     = book<TH1F>("eta_ele_midpt", "#eta^{e} 120<p_{T}<200 GeV", 50, -2.5, 2.5);
  eta_ele_highpt    = book<TH1F>("eta_ele_highpt", "#eta^{e} p_{T}>200 GeV", 50, -2.5, 2.5);
  phi_ele_lowpt     = book<TH1F>("phi_ele_lowpt", "#phi^{e} p_{T}<120 GeV", 35, -3.5, 3.5);
  phi_ele_midpt     = book<TH1F>("phi_ele_midpt", "#phi^{e} 120<p_{T}<200 GeV", 35, -3.5, 3.5);
  phi_ele_highpt    = book<TH1F>("phi_ele_highpt", "#phi^{e} p_{T}>200 GeV", 35, -3.5, 3.5);
  M_mumu            = book<TH1F>("M_mumu", "M_{#mu#mu} [GeV]",75 , 0, 500);
  M_ee              = book<TH1F>("M_ee", "M_{ee} [GeV]",75 , 0, 500);

  dRmin_mu_jet          = book<TH1F>("dRmin_mu_jet", "#DeltaR_{min}(#mu, jet)", 60, 0, 3);
  dRmin_mu_jet_scaled   = book<TH1F>("dRmin_mu_jet_scaled", "#DeltaR_{min}(#mu, jet) #times p_{T}^{jet 1}", 60, 0, 1000);
  dRmin_ele_jet         = book<TH1F>("dRmin_ele_jet", "#DeltaR_{min}(e, jet)", 60, 0, 3);
  dRmin_ele_jet_scaled  = book<TH1F>("dRmin_ele_jet_scaled", "#DeltaR_{min}(e, jet) #times p_{T}^{jet 1}", 60, 0, 1000);
  ptrel_mu_jet          = book<TH1F>("ptrel_mu_jet", "p_{T}^{rel}(#mu, jet)", 50, 0, 500);
  ptrel_ele_jet         = book<TH1F>("ptrel_ele_jet", "p_{T}^{rel}(e, jet)", 50, 0, 500);
  dRmin_mu1_jet         = book<TH1F>("dRmin_mu1_jet", "#DeltaR_{min}(#mu1, jet)", 60, 0, 3);
  dRmin_mu1_jet_scaled  = book<TH1F>("dRmin_mu1_jet_scaled", "#DeltaR_{min}(#mu1, jet) #times p_{T}^{jet 1}", 60, 0, 1000);
  dRmin_ele1_jet        = book<TH1F>("dRmin_ele1_jet", "#DeltaR_{min}(e1, jet)", 60, 0, 3);
  dRmin_ele1_jet_scaled = book<TH1F>("dRmin_ele1_jet_scaled", "#DeltaR_{min}(e1, jet) #times p_{T}^{jet 1}", 60, 0, 1000);
  ptrel_mu1_jet         = book<TH1F>("ptrel_mu1_jet", "p_{T}^{rel}(#mu1, jet)", 50, 0, 500);
  ptrel_ele1_jet        = book<TH1F>("ptrel_ele1_jet", "p_{T}^{rel}(e1, jet)", 50, 0, 500);

  dR_mu_nearjet  = book<TH1F>("dR_mu_nearjet", "#DeltaR(#mu, nearest jet)", 60, 0, 3);
  dR_ele_nearjet = book<TH1F>("dR_ele_nearjet", "#DeltaR(e, nearest jet)", 60, 0, 3);
  pt_mu_nearjet  = book<TH1F>("pt_mu_nearjet", "p_{T} (nearest jet to #mu)", 50, 0, 500);
  pt_ele_nearjet = book<TH1F>("pt_ele_nearjet", "p_{T} (nearest jet to e)", 50, 0, 500);

  dRmin_ptrel_mu   = book<TH2F>("dRmin_ptrel_mu", "#DeltaR_{min}(#mu, jet);p_{T}^{rel}(#mu, jet);p_{T}^{rel}(#mu, jet) vs. #DeltaR_{min}(#mu, jet)", 60, 0, 3, 50, 0, 500);
  dRmin_ptrel_mu1  = book<TH2F>("dRmin_ptrel_mu1", "#DeltaR_{min}(#mu1, jet);p_{T}^{rel}(#mu1, jet);p_{T}^{rel}(#mu1, jet) vs. #DeltaR_{min}(#mu1, jet)", 60, 0, 3, 50, 0, 500);
  dRmin_ptrel_ele  = book<TH2F>("dRmin_ptrel_ele", "#DeltaR_{min}(e, jet);p_{T}^{rel}(e, jet);p_{T}^{rel}(e, jet) vs. #DeltaR_{min}(e, jet)", 60, 0, 3, 50, 0, 500);
  dRmin_ptrel_ele1 = book<TH2F>("dRmin_ptrel_ele1", "#DeltaR_{min}(e1, jet);p_{T}^{rel}(e1, jet);p_{T}^{rel}(e1, jet) vs. #DeltaR_{min}(e1, jet)", 60, 0, 3, 50, 0, 500);

  // HOTVR jets
  N_HOTVRjets              = book<TH1F>("N_HOTVRjets", "N_{HOTVR jets}", 6, -0.5, 5.5);
  pt_HOTVRjet              = book<TH1F>("pt_HOTVRjet", "p_{T}^{HOTVR jets} [GeV]", 45, 0, 900);
  pt_HOTVRjet1             = book<TH1F>("pt_HOTVRjet1", "p_{T}^{HOTVR jet 1} [GeV]", 45, 0, 900);
  pt_HOTVRjet2             = book<TH1F>("pt_HOTVRjet2", "p_{T}^{HOTVR jet 2} [GeV]", 45, 0, 900);
  pt_HOTVRjet3             = book<TH1F>("pt_HOTVRjet3", "p_{T}^{HOTVR jet 3} [GeV]", 45, 0, 900);
  eta_HOTVRjet             = book<TH1F>("eta_HOTVRjet", "#eta^{HOTVR jets}", 50, -2.5, 2.5);
  eta_HOTVRjet1            = book<TH1F>("eta_HOTVRjet1", "#eta^{HOTVR jet 1}", 50, -2.5, 2.5);
  eta_HOTVRjet2            = book<TH1F>("eta_HOTVRjet2", "#eta^{HOTVR jet 2}", 50, -2.5, 2.5);
  eta_HOTVRjet3            = book<TH1F>("eta_HOTVRjet3", "#eta^{HOTVR jet 3}", 50, -2.5, 2.5);
  phi_HOTVRjet             = book<TH1F>("phi_HOTVRjet", "#phi^{HOTVR jets}", 35, -3.5, 3.5);
  phi_HOTVRjet1            = book<TH1F>("phi_HOTVRjet1", "#phi^{HOTVR jet 1}", 35, -3.5, 3.5);
  phi_HOTVRjet2            = book<TH1F>("phi_HOTVRjet2", "#phi^{HOTVR jet 2}", 35, -3.5, 3.5);
  phi_HOTVRjet3            = book<TH1F>("phi_HOTVRjet3", "#phi^{HOTVR jet 3}", 35, -3.5, 3.5);
  m_HOTVRjet               = book<TH1F>("m_HOTVRjet", "m^{HOTVR jets}", 50, 0, 500);
  m_HOTVRjet1              = book<TH1F>("m_HOTVRjet1", "m^{HOTVR jet 1}", 50, 0, 500);
  m_HOTVRjet2              = book<TH1F>("m_HOTVRjet2", "m^{HOTVR jet 2}", 50, 0, 500);
  m_HOTVRjet3              = book<TH1F>("m_HOTVRjet3", "m^{HOTVR jet 3}", 50, 0, 500);
  N_subjets_HOTVRjet       = book<TH1F>("N_subjets_HOTVRjet", "N_{subjets}^{HOTVR jets}", 6, -0.5, 5.5);
  N_subjets_HOTVRjet1      = book<TH1F>("N_subjets_HOTVRjet1", "N_{subjets}^{HOTVR jet 1}", 6, -0.5, 5.5);
  N_subjets_HOTVRjet2      = book<TH1F>("N_subjets_HOTVRjet2", "N_{subjets}^{HOTVR jet 2}", 6, -0.5, 5.5);
  N_subjets_HOTVRjet3      = book<TH1F>("N_subjets_HOTVRjet3", "N_{subjets}^{HOTVR jet 3}", 6, -0.5, 5.5);
  N_daughters_HOTVRjet     = book<TH1F>("N_daughters_HOTVRjet", "N_{daughters}^{HOTVR jets}", 51, -0.5, 50.5);
  N_daughters_HOTVRjet1    = book<TH1F>("N_daughters_HOTVRjet1", "N_{daughters}^{HOTVR jet 1}", 51, -0.5, 50.5);
  N_daughters_HOTVRjet2    = book<TH1F>("N_daughters_HOTVRjet2", "N_{daughters}^{HOTVR jet 2}", 51, -0.5, 50.5);
  N_daughters_HOTVRjet3    = book<TH1F>("N_daughters_HOTVRjet3", "N_{daughters}^{HOTVR jet 3}", 51, -0.5, 50.5);
  dRmin_AK8Puppi_HOTVRjet  = book<TH1F>("dRmin_AK8Puppi_HOTVRjet", "#DeltaR_{min}(HOTVR jet, AK8Puppi jet)", 60, 0, 3);
  dRmin_AK8Puppi_HOTVRjet1 = book<TH1F>("dRmin_AK8Puppi_HOTVRjet1", "#DeltaR_{min}(HOTVR jet 1, AK8Puppi jet)", 60, 0, 3);
  dRmin_AK8Puppi_HOTVRjet2 = book<TH1F>("dRmin_AK8Puppi_HOTVRjet2", "#DeltaR_{min}(HOTVR jet 2, AK8Puppi jet)", 60, 0, 3);
  dRmin_AK8Puppi_HOTVRjet3 = book<TH1F>("dRmin_AK8Puppi_HOTVRjet3", "#DeltaR_{min}(HOTVR jet 3, AK8Puppi jet)", 60, 0, 3);
  dRmin_mu_HOTVRjet        = book<TH1F>("dRmin_mu_HOTVRjet", "#DeltaR_{min}(HOTVR jet, #mu)", 60, 0, 3);
  dRmin_mu_HOTVRjet1       = book<TH1F>("dRmin_mu_HOTVRjet1", "#DeltaR_{min}(HOTVR jet 1, #mu)", 60, 0, 3);
  dRmin_mu_HOTVRjet2       = book<TH1F>("dRmin_mu_HOTVRjet2", "#DeltaR_{min}(HOTVR jet 2, #mu)", 60, 0, 3);
  dRmin_mu_HOTVRjet3       = book<TH1F>("dRmin_mu_HOTVRjet3", "#DeltaR_{min}(HOTVR jet 3, #mu)", 60, 0, 3);
  tau1_HOTVRjet            = book<TH1F>("tau1_HOTVRjet", "#tau_{1}^{HOTVR jets}", 24, 0, 1.2);
  tau1_HOTVRjet1           = book<TH1F>("tau1_HOTVRjet1", "#tau_{1}^{HOTVR jet 1}", 24, 0, 1.2);
  tau1_HOTVRjet2           = book<TH1F>("tau1_HOTVRjet2", "#tau_{1}^{HOTVR jet 2}", 24, 0, 1.2);
  tau1_HOTVRjet3           = book<TH1F>("tau1_HOTVRjet3", "#tau_{1}^{HOTVR jet 3}", 24, 0, 1.2);
  tau2_HOTVRjet            = book<TH1F>("tau2_HOTVRjet", "#tau_{2}^{HOTVR jets}", 24, 0, 1.2);
  tau2_HOTVRjet1           = book<TH1F>("tau2_HOTVRjet1", "#tau_{2}^{HOTVR jet 1}", 24, 0, 1.2);
  tau2_HOTVRjet2           = book<TH1F>("tau2_HOTVRjet2", "#tau_{2}^{HOTVR jet 2}", 24, 0, 1.2);
  tau2_HOTVRjet3           = book<TH1F>("tau2_HOTVRjet3", "#tau_{2}^{HOTVR jet 3}", 24, 0, 1.2);
  tau3_HOTVRjet            = book<TH1F>("tau3_HOTVRjet", "#tau_{3}^{HOTVR jets}", 24, 0, 1.2);
  tau3_HOTVRjet1           = book<TH1F>("tau3_HOTVRjet1", "#tau_{3}^{HOTVR jet 1}", 24, 0, 1.2);
  tau3_HOTVRjet2           = book<TH1F>("tau3_HOTVRjet2", "#tau_{3}^{HOTVR jet 2}", 24, 0, 1.2);
  tau3_HOTVRjet3           = book<TH1F>("tau3_HOTVRjet3", "#tau_{3}^{HOTVR jet 3}", 24, 0, 1.2);
  tau21_HOTVRjet           = book<TH1F>("tau21_HOTVRjet", "#tau_{2/1}^{HOTVR jets}", 24, 0, 1.2);
  tau21_HOTVRjet1          = book<TH1F>("tau21_HOTVRjet1", "#tau_{2/1}^{HOTVR jet 1}", 24, 0, 1.2);
  tau21_HOTVRjet2          = book<TH1F>("tau21_HOTVRjet2", "#tau_{2/1}^{HOTVR jet 2}", 24, 0, 1.2);
  tau21_HOTVRjet3          = book<TH1F>("tau21_HOTVRjet3", "#tau_{2/1}^{HOTVR jet 3}", 24, 0, 1.2);
  tau32_HOTVRjet           = book<TH1F>("tau32_HOTVRjet", "#tau_{3/2}^{HOTVR jets}", 24, 0, 1.2);
  tau32_HOTVRjet1          = book<TH1F>("tau32_HOTVRjet1", "#tau_{3/2}^{HOTVR jet 1}", 24, 0, 1.2);
  tau32_HOTVRjet2          = book<TH1F>("tau32_HOTVRjet2", "#tau_{3/2}^{HOTVR jet 2}", 24, 0, 1.2);
  tau32_HOTVRjet3          = book<TH1F>("tau32_HOTVRjet3", "#tau_{3/2}^{HOTVR jet 3}", 24, 0, 1.2);

  // AK8Puppi jets
  N_AK8Puppijets           = book<TH1F>("N_AK8Puppijets", "N_{AK8Puppi jets}", 6, -0.5, 5.5);
  pt_AK8Puppijet           = book<TH1F>("pt_AK8Puppijet", "p_{T}^{AK8Puppi jets} [GeV]", 45, 0, 900);
  pt_AK8Puppijet1          = book<TH1F>("pt_AK8Puppijet1", "p_{T}^{AK8Puppi jet 1} [GeV]", 45, 0, 900);
  pt_AK8Puppijet2          = book<TH1F>("pt_AK8Puppijet2", "p_{T}^{AK8Puppi jet 2} [GeV]", 45, 0, 900);
  pt_AK8Puppijet3          = book<TH1F>("pt_AK8Puppijet3", "p_{T}^{AK8Puppi jet 3} [GeV]", 45, 0, 900);
  eta_AK8Puppijet          = book<TH1F>("eta_AK8Puppijet", "#eta^{AK8Puppi jets}", 50, -2.5, 2.5);
  eta_AK8Puppijet1         = book<TH1F>("eta_AK8Puppijet1", "#eta^{AK8Puppi jet 1}", 50, -2.5, 2.5);
  eta_AK8Puppijet2         = book<TH1F>("eta_AK8Puppijet2", "#eta^{AK8Puppi jet 2}", 50, -2.5, 2.5);
  eta_AK8Puppijet3         = book<TH1F>("eta_AK8Puppijet3", "#eta^{AK8Puppi jet 3}", 50, -2.5, 2.5);
  phi_AK8Puppijet          = book<TH1F>("phi_AK8Puppijet", "#phi^{AK8Puppi jets}", 35, -3.5, 3.5);
  phi_AK8Puppijet1         = book<TH1F>("phi_AK8Puppijet1", "#phi^{AK8Puppi jet 1}", 35, -3.5, 3.5);
  phi_AK8Puppijet2         = book<TH1F>("phi_AK8Puppijet2", "#phi^{AK8Puppi jet 2}", 35, -3.5, 3.5);
  phi_AK8Puppijet3         = book<TH1F>("phi_AK8Puppijet3", "#phi^{AK8Puppi jet 3}", 35, -3.5, 3.5);
  mSD_AK8Puppijet          = book<TH1F>("mSD_AK8Puppijet", "m_{SD}^{AK8Puppi jets}", 50, 0, 500);
  mSD_AK8Puppijet1         = book<TH1F>("mSD_AK8Puppijet1", "m_{SD}^{AK8Puppi jet 1}", 50, 0, 500);
  mSD_AK8Puppijet2         = book<TH1F>("mSD_AK8Puppijet2", "m_{SD}^{AK8Puppi jet 2}", 50, 0, 500);
  mSD_AK8Puppijet3         = book<TH1F>("mSD_AK8Puppijet3", "m_{SD}^{AK8Puppi jet 3}", 50, 0, 500);
  N_subjets_AK8Puppijet    = book<TH1F>("N_subjets_AK8Puppijet", "N_{subjets}^{AK8Puppi jets}", 6, -0.5, 5.5);
  N_subjets_AK8Puppijet1   = book<TH1F>("N_subjets_AK8Puppijet1", "N_{subjets}^{AK8Puppi jet 1}", 6, -0.5, 5.5);
  N_subjets_AK8Puppijet2   = book<TH1F>("N_subjets_AK8Puppijet2", "N_{subjets}^{AK8Puppi jet 2}", 6, -0.5, 5.5);
  N_subjets_AK8Puppijet3   = book<TH1F>("N_subjets_AK8Puppijet3", "N_{subjets}^{AK8Puppi jet 3}", 6, -0.5, 5.5);
  N_daughters_AK8Puppijet  = book<TH1F>("N_daughters_AK8Puppijet", "N_{daughters}^{AK8Puppi jets}", 51, -0.5, 50.5);
  N_daughters_AK8Puppijet1 = book<TH1F>("N_daughters_AK8Puppijet1", "N_{daughters}^{AK8Puppi jet 1}", 51, -0.5, 50.5);
  N_daughters_AK8Puppijet2 = book<TH1F>("N_daughters_AK8Puppijet2", "N_{daughters}^{AK8Puppi jet 2}", 51, -0.5, 50.5);
  N_daughters_AK8Puppijet3 = book<TH1F>("N_daughters_AK8Puppijet3", "N_{daughters}^{AK8Puppi jet 3}", 51, -0.5, 50.5);
  dRmin_AK8_AK8Puppijet    = book<TH1F>("dRmin_AK8_AK8Puppijet", "#DeltaR_{min}(AK8Puppi jet, AK8 jet)", 60, 0, 3);
  dRmin_AK8_AK8Puppijet1   = book<TH1F>("dRmin_AK8_AK8Puppijet1", "#DeltaR_{min}(AK8Puppi jet 1, AK8 jet)", 60, 0, 3);
  dRmin_AK8_AK8Puppijet2   = book<TH1F>("dRmin_AK8_AK8Puppijet2", "#DeltaR_{min}(AK8Puppi jet 2, AK8 jet)", 60, 0, 3);
  dRmin_AK8_AK8Puppijet3   = book<TH1F>("dRmin_AK8_AK8Puppijet3", "#DeltaR_{min}(AK8Puppi jet 3, AK8 jet)", 60, 0, 3);
  dRmin_mu_AK8Puppijet     = book<TH1F>("dRmin_mu_AK8Puppijet", "#DeltaR_{min}(AK8Puppi jet, #mu)", 60, 0, 3);
  dRmin_mu_AK8Puppijet1    = book<TH1F>("dRmin_mu_AK8Puppijet1", "#DeltaR_{min}(AK8Puppi jet 1, #mu)", 60, 0, 3);
  dRmin_mu_AK8Puppijet2    = book<TH1F>("dRmin_mu_AK8Puppijet2", "#DeltaR_{min}(AK8Puppi jet 2, #mu)", 60, 0, 3);
  dRmin_mu_AK8Puppijet3    = book<TH1F>("dRmin_mu_AK8Puppijet3", "#DeltaR_{min}(AK8Puppi jet 3, #mu)", 60, 0, 3);
  tau1_AK8Puppijet         = book<TH1F>("tau1_AK8Puppijet", "#tau_{1}^{AK8Puppi jets}", 24, 0, 1.2);
  tau1_AK8Puppijet1        = book<TH1F>("tau1_AK8Puppijet1", "#tau_{1}^{AK8Puppi jet 1}", 24, 0, 1.2);
  tau1_AK8Puppijet2        = book<TH1F>("tau1_AK8Puppijet2", "#tau_{1}^{AK8Puppi jet 2}", 24, 0, 1.2);
  tau1_AK8Puppijet3        = book<TH1F>("tau1_AK8Puppijet3", "#tau_{1}^{AK8Puppi jet 3}", 24, 0, 1.2);
  tau2_AK8Puppijet         = book<TH1F>("tau2_AK8Puppijet", "#tau_{2}^{AK8Puppi jets}", 24, 0, 1.2);
  tau2_AK8Puppijet1        = book<TH1F>("tau2_AK8Puppijet1", "#tau_{2}^{AK8Puppi jet 1}", 24, 0, 1.2);
  tau2_AK8Puppijet2        = book<TH1F>("tau2_AK8Puppijet2", "#tau_{2}^{AK8Puppi jet 2}", 24, 0, 1.2);
  tau2_AK8Puppijet3        = book<TH1F>("tau2_AK8Puppijet3", "#tau_{2}^{AK8Puppi jet 3}", 24, 0, 1.2);
  tau3_AK8Puppijet         = book<TH1F>("tau3_AK8Puppijet", "#tau_{3}^{AK8Puppi jets}", 24, 0, 1.2);
  tau3_AK8Puppijet1        = book<TH1F>("tau3_AK8Puppijet1", "#tau_{3}^{AK8Puppi jet 1}", 24, 0, 1.2);
  tau3_AK8Puppijet2        = book<TH1F>("tau3_AK8Puppijet2", "#tau_{3}^{AK8Puppi jet 2}", 24, 0, 1.2);
  tau3_AK8Puppijet3        = book<TH1F>("tau3_AK8Puppijet3", "#tau_{3}^{AK8Puppi jet 3}", 24, 0, 1.2);
  tau21_AK8Puppijet        = book<TH1F>("tau21_AK8Puppijet", "#tau_{2/1}^{AK8Puppi jets}", 24, 0, 1.2);
  tau21_AK8Puppijet1       = book<TH1F>("tau21_AK8Puppijet1", "#tau_{2/1}^{AK8Puppi jet 1}", 24, 0, 1.2);
  tau21_AK8Puppijet2       = book<TH1F>("tau21_AK8Puppijet2", "#tau_{2/1}^{AK8Puppi jet 2}", 24, 0, 1.2);
  tau21_AK8Puppijet3       = book<TH1F>("tau21_AK8Puppijet3", "#tau_{2/1}^{AK8Puppi jet 3}", 24, 0, 1.2);
  tau32_AK8Puppijet        = book<TH1F>("tau32_AK8Puppijet", "#tau_{3/2}^{AK8Puppi jets}", 24, 0, 1.2);
  tau32_AK8Puppijet1       = book<TH1F>("tau32_AK8Puppijet1", "#tau_{3/2}^{AK8Puppi jet 1}", 24, 0, 1.2);
  tau32_AK8Puppijet2       = book<TH1F>("tau32_AK8Puppijet2", "#tau_{3/2}^{AK8Puppi jet 2}", 24, 0, 1.2);
  tau32_AK8Puppijet3       = book<TH1F>("tau32_AK8Puppijet3", "#tau_{3/2}^{AK8Puppi jet 3}", 24, 0, 1.2);

  deepjet_topscore_jet    = book<TH1F>("deepjet_topscore_jet", "DeepJet top score all AK8 jets", 20, 0, 1);
  deepjet_topscore_jet1   = book<TH1F>("deepjet_topscore_jet1", "DeepJet top score AK8 jet 1", 20, 0, 1);
  deepjet_topscore_jet2   = book<TH1F>("deepjet_topscore_jet2", "DeepJet top score AK8 jet 2", 20, 0, 1);
  deepjet_topscore_jet3   = book<TH1F>("deepjet_topscore_jet3", "DeepJet top score AK8 jet 3", 20, 0, 1);
  deepjet_wscore_jet      = book<TH1F>("deepjet_wscore_jet", "DeepJet W score all AK8 jets", 20, 0, 1);
  deepjet_wscore_jet1     = book<TH1F>("deepjet_wscore_jet1", "DeepJet W score AK8 jet 1", 20, 0, 1);
  deepjet_wscore_jet2     = book<TH1F>("deepjet_wscore_jet2", "DeepJet W score AK8 jet 2", 20, 0, 1);
  deepjet_zscore_jet      = book<TH1F>("deepjet_zscore_jet", "DeepJet Z score all AK8 jets", 20, 0, 1);
  deepjet_zscore_jet1     = book<TH1F>("deepjet_zscore_jet1", "DeepJet Z score AK8 jet 1", 20, 0, 1);
  deepjet_zscore_jet2     = book<TH1F>("deepjet_zscore_jet2", "DeepJet Z score AK8 jet 2", 20, 0, 1);
  deepjet_higgsscore_jet  = book<TH1F>("deepjet_higgsscore_jet", "DeepJet Higgs score all AK8 jets", 20, 0, 1);
  deepjet_higgsscore_jet1 = book<TH1F>("deepjet_higgsscore_jet1", "DeepJet Higgs score AK8 jet 1", 20, 0, 1);
  deepjet_higgsscore_jet2 = book<TH1F>("deepjet_higgsscore_jet2", "DeepJet Higgs score AK8 jet 2", 20, 0, 1);
  deepjet_qcdscore_jet    = book<TH1F>("deepjet_qcdscore_jet", "DeepJet qcd score all AK8 jets", 20, 0, 1);
  deepjet_qcdscore_jet1   = book<TH1F>("deepjet_qcdscore_jet1", "DeepJet qcd score AK8 jet 1", 20, 0, 1);
  deepjet_qcdscore_jet2   = book<TH1F>("deepjet_qcdscore_jet2", "DeepJet qcd score AK8 jet 2}", 20, 0, 1);

  deepjet_TvsQCD_jet  = book<TH1F>("deepjet_topProb_jet", "DeepJet  TvsQCD all AK8 jets", 20, 0, 1);
  deepjet_TvsQCD_jet1 = book<TH1F>("deepjet_topProb_jet1", "DeepJet TvsQCD AK8 jet 1", 20, 0, 1);
  deepjet_TvsQCD_jet2 = book<TH1F>("deepjet_topProb_jet2", "DeepJet TvsQCD AK8 jet 2", 20, 0, 1);
  deepjet_TvsQCD_jet3 = book<TH1F>("deepjet_topProb_jet3", "DeepJet TvsQCD AK8 jet 3", 20, 0, 1);

  massdecordeepjet_TvsQCD_jet  = book<TH1F>("massdecordeepjet_topProb_jet",  "Mass decor deepJet TvsQCD all AK8 jets", 20, 0, 1);
  massdecordeepjet_TvsQCD_jet1 = book<TH1F>("massdecordeepjet_topProb_jet1", "Mass decor DeepJet TvsQCD AK8 jet 1", 20, 0, 1);
  massdecordeepjet_TvsQCD_jet2 = book<TH1F>("massdecordeepjet_topProb_jet2", "Mass decor DeepJet TvsQCD AK8 jet 2", 20, 0, 1);
  massdecordeepjet_TvsQCD_jet3 = book<TH1F>("massdecordeepjet_topProb_jet3", "Mass decor DeepJet TvsQCD AK8 jet 3", 20, 0, 1);

  massdecordeepjet_topscore_jet    = book<TH1F>("massdecordeepjet_topscore_jet", "Mass decor deepjet top score all AK8 jets", 20, 0, 1);
  massdecordeepjet_topscore_jet1   = book<TH1F>("massdecordeepjet_topscore_jet1", "Mass decor deepjet top score AK8 jet 1", 20, 0, 1);
  massdecordeepjet_topscore_jet2   = book<TH1F>("massdecordeepjet_topscore_jet2", "Mass decor deepjet top score AK8 jet 2", 20, 0, 1);
  massdecordeepjet_topscore_jet3   = book<TH1F>("massdecordeepjet_topscore_jet3", "Mass decor deepjet top score AK8 jet 3", 20, 0, 1);
  massdecordeepjet_wscore_jet      = book<TH1F>("massdecordeepjet_wscore_jet", "Mass decor deepjet W score all AK8 jets", 20, 0, 1);
  massdecordeepjet_wscore_jet1     = book<TH1F>("massdecordeepjet_wscore_jet1", "Mass decor deepjet  W score AK8 jet 1", 20, 0, 1);
  massdecordeepjet_wscore_jet2     = book<TH1F>("massdecordeepjet_wscore_jet2", "Mass decor deepjet W score AK8 jet 2", 20, 0, 1);
  massdecordeepjet_zscore_jet      = book<TH1F>("massdecordeepjet_zscore_jet", "Mass decor deepjet Z score all AK8 jets", 20, 0, 1);
  massdecordeepjet_zscore_jet1     = book<TH1F>("massdecordeepjet_zscore_jet1", "Mass decor deepjet  Z score AK8 jet 1", 20, 0, 1);
  massdecordeepjet_zscore_jet2     = book<TH1F>("massdecordeepjet_zscore_jet2", "Mass decor deepjet Z score AK8 jet 2", 20, 0, 1);
  massdecordeepjet_higgsscore_jet  = book<TH1F>("massdecordeepjet_higgsscore_jet", "Mass decor deepjet Higgs score all AK8 jets", 20, 0, 1);
  massdecordeepjet_higgsscore_jet1 = book<TH1F>("massdecordeepjet_higgsscore_jet1", "Mass decor deepjet  Higgs score AK8 jet 1", 20, 0, 1);
  massdecordeepjet_higgsscore_jet2 = book<TH1F>("massdecordeepjet_higgsscore_jet2", "Mass decor deepjet Higgs score AK8 jet 2", 20, 0, 1);
  massdecordeepjet_qcdscore_jet    = book<TH1F>("massdecordeepjet_qcdscore_jet", "Mass decor deepjet qcd score all AK8 jets", 20, 0, 1);
  massdecordeepjet_qcdscore_jet1   = book<TH1F>("massdecordeepjet_qcdscore_jet1", "Mass decor deepjet qcd score AK8 jet 1", 20, 0, 1);
  massdecordeepjet_qcdscore_jet2   = book<TH1F>("massdecordeepjet_qcdscore_jet2", "Mass decor deepjet qcd score AK8 jet 2", 20, 0, 1);


  // HOTVR t-tagged
  N_HOTVRTaggedjets           = book<TH1F>("N_HOTVRTaggedjets", "N_{HOTVR Tagged jets}", 6, -0.5, 5.5);
  pt_HOTVRTaggedjet           = book<TH1F>("pt_HOTVRTaggedjet", "p_{T}^{HOTVR Tagged jets} [GeV]", 45, 0, 900);
  pt_HOTVRTaggedjet1          = book<TH1F>("pt_HOTVRTaggedjet1", "p_{T}^{HOTVR Tagged jet 1} [GeV]", 45, 0, 900);
  pt_HOTVRTaggedjet2          = book<TH1F>("pt_HOTVRTaggedjet2", "p_{T}^{HOTVR Tagged jet 2} [GeV]", 45, 0, 900);
  pt_HOTVRTaggedjet3          = book<TH1F>("pt_HOTVRTaggedjet3", "p_{T}^{HOTVR Tagged jet 3} [GeV]", 45, 0, 900);
  eta_HOTVRTaggedjet          = book<TH1F>("eta_HOTVRTaggedjet", "#eta^{HOTVR Tagged jets}", 50, -2.5, 2.5);
  eta_HOTVRTaggedjet1         = book<TH1F>("eta_HOTVRTaggedjet1", "#eta^{HOTVR Tagged jet 1}", 50, -2.5, 2.5);
  eta_HOTVRTaggedjet2         = book<TH1F>("eta_HOTVRTaggedjet2", "#eta^{HOTVR Tagged jet 2}", 50, -2.5, 2.5);
  eta_HOTVRTaggedjet3         = book<TH1F>("eta_HOTVRTaggedjet3", "#eta^{HOTVR Tagged jet 3}", 50, -2.5, 2.5);
  phi_HOTVRTaggedjet          = book<TH1F>("phi_HOTVRTaggedjet", "#phi^{HOTVR Tagged jets}", 35, -3.5, 3.5);
  phi_HOTVRTaggedjet1         = book<TH1F>("phi_HOTVRTaggedjet1", "#phi^{HOTVR Tagged jet 1}", 35, -3.5, 3.5);
  phi_HOTVRTaggedjet2         = book<TH1F>("phi_HOTVRTaggedjet2", "#phi^{HOTVR Tagged jet 2}", 35, -3.5, 3.5);
  phi_HOTVRTaggedjet3         = book<TH1F>("phi_HOTVRTaggedjet3", "#phi^{HOTVR Tagged jet 3}", 35, -3.5, 3.5);
  m_HOTVRTaggedjet            = book<TH1F>("m_HOTVRTaggedjet", "m^{HOTVR Tagged jets}", 50, 0, 500);
  m_HOTVRTaggedjet1           = book<TH1F>("m_HOTVRTaggedjet1", "m^{HOTVR Tagged jet 1}", 50, 0, 500);
  m_HOTVRTaggedjet2           = book<TH1F>("m_HOTVRTaggedjet2", "m^{HOTVR Tagged jet 2}", 50, 0, 500);
  m_HOTVRTaggedjet3           = book<TH1F>("m_HOTVRTaggedjet3", "m^{HOTVR Tagged jet 3}", 50, 0, 500);
  N_subjets_HOTVRTaggedjet    = book<TH1F>("N_subjets_HOTVRTaggedjet", "N_{subjets}^{HOTVR Tagged jets}", 6, -0.5, 5.5);
  N_subjets_HOTVRTaggedjet1   = book<TH1F>("N_subjets_HOTVRTaggedjet1", "N_{subjets}^{HOTVR Tagged jet 1}", 6, -0.5, 5.5);
  N_subjets_HOTVRTaggedjet2   = book<TH1F>("N_subjets_HOTVRTaggedjet2", "N_{subjets}^{HOTVR Tagged jet 2}", 6, -0.5, 5.5);
  N_subjets_HOTVRTaggedjet3   = book<TH1F>("N_subjets_HOTVRTaggedjet3", "N_{subjets}^{HOTVR Tagged jet 3}", 6, -0.5, 5.5);
  N_daughters_HOTVRTaggedjet  = book<TH1F>("N_daughters_HOTVRTaggedjet", "N_{daughters}^{HOTVR Tagged jets}", 51, -0.5, 50.5);
  N_daughters_HOTVRTaggedjet1 = book<TH1F>("N_daughters_HOTVRTaggedjet1", "N_{daughters}^{HOTVR Tagged jet 1}", 51, -0.5, 50.5);
  N_daughters_HOTVRTaggedjet2 = book<TH1F>("N_daughters_HOTVRTaggedjet2", "N_{daughters}^{HOTVR Tagged jet 2}", 51, -0.5, 50.5);
  N_daughters_HOTVRTaggedjet3 = book<TH1F>("N_daughters_HOTVRTaggedjet3", "N_{daughters}^{HOTVR Tagged jet 3}", 51, -0.5, 50.5);
  dRmin_AK8_HOTVRTaggedjet    = book<TH1F>("dRmin_AK8_HOTVRTaggedjet", "#DeltaR_{min}(HOTVR Tagged jet, AK8 jet)", 60, 0, 3);
  dRmin_AK8_HOTVRTaggedjet1   = book<TH1F>("dRmin_AK8_HOTVRTaggedjet1", "#DeltaR_{min}(HOTVR Tagged jet 1, AK8 jet)", 60, 0, 3);
  dRmin_AK8_HOTVRTaggedjet2   = book<TH1F>("dRmin_AK8_HOTVRTaggedjet2", "#DeltaR_{min}(HOTVR Tagged jet 2, AK8 jet)", 60, 0, 3);
  dRmin_AK8_HOTVRTaggedjet3   = book<TH1F>("dRmin_AK8_HOTVRTaggedjet3", "#DeltaR_{min}(HOTVR Tagged jet 3, AK8 jet)", 60, 0, 3);
  dRmin_mu_HOTVRTaggedjet     = book<TH1F>("dRmin_mu_HOTVRTaggedjet", "#DeltaR_{min}(HOTVR Tagged jet, #mu)", 60, 0, 3);
  dRmin_mu_HOTVRTaggedjet1    = book<TH1F>("dRmin_mu_HOTVRTaggedjet1", "#DeltaR_{min}(HOTVR Tagged jet 1, #mu)", 60, 0, 3);
  dRmin_mu_HOTVRTaggedjet2    = book<TH1F>("dRmin_mu_HOTVRTaggedjet2", "#DeltaR_{min}(HOTVR Tagged jet 2, #mu)", 60, 0, 3);
  dRmin_mu_HOTVRTaggedjet3    = book<TH1F>("dRmin_mu_HOTVRTaggedjet3", "#DeltaR_{min}(HOTVR Tagged jet 3, #mu)", 60, 0, 3);
  tau1_HOTVRTaggedjet         = book<TH1F>("tau1_HOTVRTaggedjet", "#tau_{1}^{HOTVR Tagged jets}", 24, 0, 1.2);
  tau1_HOTVRTaggedjet1        = book<TH1F>("tau1_HOTVRTaggedjet1", "#tau_{1}^{HOTVR Tagged jet 1}", 24, 0, 1.2);
  tau1_HOTVRTaggedjet2        = book<TH1F>("tau1_HOTVRTaggedjet2", "#tau_{1}^{HOTVR Tagged jet 2}", 24, 0, 1.2);
  tau1_HOTVRTaggedjet3        = book<TH1F>("tau1_HOTVRTaggedjet3", "#tau_{1}^{HOTVR Tagged jet 3}", 24, 0, 1.2);
  tau2_HOTVRTaggedjet         = book<TH1F>("tau2_HOTVRTaggedjet", "#tau_{2}^{HOTVR Tagged jets}", 24, 0, 1.2);
  tau2_HOTVRTaggedjet1        = book<TH1F>("tau2_HOTVRTaggedjet1", "#tau_{2}^{HOTVR Tagged jet 1}", 24, 0, 1.2);
  tau2_HOTVRTaggedjet2        = book<TH1F>("tau2_HOTVRTaggedjet2", "#tau_{2}^{HOTVR Tagged jet 2}", 24, 0, 1.2);
  tau2_HOTVRTaggedjet3        = book<TH1F>("tau2_HOTVRTaggedjet3", "#tau_{2}^{HOTVR Tagged jet 3}", 24, 0, 1.2);
  tau3_HOTVRTaggedjet         = book<TH1F>("tau3_HOTVRTaggedjet", "#tau_{3}^{HOTVR Tagged jets}", 24, 0, 1.2);
  tau3_HOTVRTaggedjet1        = book<TH1F>("tau3_HOTVRTaggedjet1", "#tau_{3}^{HOTVR Tagged jet 1}", 24, 0, 1.2);
  tau3_HOTVRTaggedjet2        = book<TH1F>("tau3_HOTVRTaggedjet2", "#tau_{3}^{HOTVR Tagged jet 2}", 24, 0, 1.2);
  tau3_HOTVRTaggedjet3        = book<TH1F>("tau3_HOTVRTaggedjet3", "#tau_{3}^{HOTVR Tagged jet 3}", 24, 0, 1.2);
  tau21_HOTVRTaggedjet        = book<TH1F>("tau21_HOTVRTaggedjet", "#tau_{2/1}^{HOTVR Tagged jets}", 24, 0, 1.2);
  tau21_HOTVRTaggedjet1       = book<TH1F>("tau21_HOTVRTaggedjet1", "#tau_{2/1}^{HOTVR Tagged jet 1}", 24, 0, 1.2);
  tau21_HOTVRTaggedjet2       = book<TH1F>("tau21_HOTVRTaggedjet2", "#tau_{2/1}^{HOTVR Tagged jet 2}", 24, 0, 1.2);
  tau21_HOTVRTaggedjet3       = book<TH1F>("tau21_HOTVRTaggedjet3", "#tau_{2/1}^{HOTVR Tagged jet 3}", 24, 0, 1.2);
  tau32_HOTVRTaggedjet        = book<TH1F>("tau32_HOTVRTaggedjet", "#tau_{3/2}^{HOTVR Tagged jets}", 24, 0, 1.2);
  tau32_HOTVRTaggedjet1       = book<TH1F>("tau32_HOTVRTaggedjet1", "#tau_{3/2}^{HOTVR Tagged jet 1}", 24, 0, 1.2);
  tau32_HOTVRTaggedjet2       = book<TH1F>("tau32_HOTVRTaggedjet2", "#tau_{3/2}^{HOTVR Tagged jet 2}", 24, 0, 1.2);
  tau32_HOTVRTaggedjet3       = book<TH1F>("tau32_HOTVRTaggedjet3", "#tau_{3/2}^{HOTVR  Tagged jet 3}", 24, 0, 1.2);

  // AK8Puppi t-tagged
  N_AK8PuppiTaggedjets           = book<TH1F>("N_AK8PuppiTaggedjets", "N_{AK8Puppi Tagged jets}", 6, -0.5, 5.5);
  pt_AK8PuppiTaggedjet           = book<TH1F>("pt_AK8PuppiTaggedjet", "p_{T}^{AK8Puppi Tagged jets} [GeV]", 45, 0, 900);
  pt_AK8PuppiTaggedjet1          = book<TH1F>("pt_AK8PuppiTaggedjet1", "p_{T}^{AK8Puppi Tagged jet 1} [GeV]", 45, 0, 900);
  pt_AK8PuppiTaggedjet2          = book<TH1F>("pt_AK8PuppiTaggedjet2", "p_{T}^{AK8Puppi Tagged jet 2} [GeV]", 45, 0, 900);
  pt_AK8PuppiTaggedjet3          = book<TH1F>("pt_AK8PuppiTaggedjet3", "p_{T}^{AK8Puppi Tagged jet 3} [GeV]", 45, 0, 900);
  eta_AK8PuppiTaggedjet          = book<TH1F>("eta_AK8PuppiTaggedjet", "#eta^{AK8Puppi Tagged jets}", 50, -2.5, 2.5);
  eta_AK8PuppiTaggedjet1         = book<TH1F>("eta_AK8PuppiTaggedjet1", "#eta^{AK8Puppi Tagged jet 1}", 50, -2.5, 2.5);
  eta_AK8PuppiTaggedjet2         = book<TH1F>("eta_AK8PuppiTaggedjet2", "#eta^{AK8Puppi Tagged jet 2}", 50, -2.5, 2.5);
  eta_AK8PuppiTaggedjet3         = book<TH1F>("eta_AK8PuppiTaggedjet3", "#eta^{AK8Puppi Tagged jet 3}", 50, -2.5, 2.5);
  phi_AK8PuppiTaggedjet          = book<TH1F>("phi_AK8PuppiTaggedjet", "#phi^{AK8Puppi Tagged jets}", 35, -3.5, 3.5);
  phi_AK8PuppiTaggedjet1         = book<TH1F>("phi_AK8PuppiTaggedjet1", "#phi^{AK8Puppi Tagged jet 1}", 35, -3.5, 3.5);
  phi_AK8PuppiTaggedjet2         = book<TH1F>("phi_AK8PuppiTaggedjet2", "#phi^{AK8Puppi Tagged jet 2}", 35, -3.5, 3.5);
  phi_AK8PuppiTaggedjet3         = book<TH1F>("phi_AK8PuppiTaggedjet3", "#phi^{AK8Puppi Tagged jet 3}", 35, -3.5, 3.5);
  mSD_AK8PuppiTaggedjet          = book<TH1F>("mSD_AK8PuppiTaggedjet", "m_{SD}^{AK8Puppi Tagged jets}", 50, 0, 500);
  mSD_AK8PuppiTaggedjet1         = book<TH1F>("mSD_AK8PuppiTaggedjet1", "m_{SD}^{AK8Puppi Tagged jet 1}", 50, 0, 500);
  mSD_AK8PuppiTaggedjet2         = book<TH1F>("mSD_AK8PuppiTaggedjet2", "m_{SD}^{AK8Puppi Tagged jet 2}", 50, 0, 500);
  mSD_AK8PuppiTaggedjet3         = book<TH1F>("mSD_AK8PuppiTaggedjet3", "m_{SD}^{AK8Puppi Tagged jet 3}", 50, 0, 500);
  N_subjets_AK8PuppiTaggedjet    = book<TH1F>("N_subjets_AK8PuppiTaggedjet", "N_{subjets}^{AK8Puppi Tagged jets}", 6, -0.5, 5.5);
  N_subjets_AK8PuppiTaggedjet1   = book<TH1F>("N_subjets_AK8PuppiTaggedjet1", "N_{subjets}^{AK8Puppi Tagged jet 1}", 6, -0.5, 5.5);
  N_subjets_AK8PuppiTaggedjet2   = book<TH1F>("N_subjets_AK8PuppiTaggedjet2", "N_{subjets}^{AK8Puppi Tagged jet 2}", 6, -0.5, 5.5);
  N_subjets_AK8PuppiTaggedjet3   = book<TH1F>("N_subjets_AK8PuppiTaggedjet3", "N_{subjets}^{AK8Puppi Tagged jet 3}", 6, -0.5, 5.5);
  N_daughters_AK8PuppiTaggedjet  = book<TH1F>("N_daughters_AK8PuppiTaggedjet", "N_{daughters}^{AK8Puppi Tagged jets}", 51, -0.5, 50.5);
  N_daughters_AK8PuppiTaggedjet1 = book<TH1F>("N_daughters_AK8PuppiTaggedjet1", "N_{daughters}^{AK8Puppi Tagged jet 1}", 51, -0.5, 50.5);
  N_daughters_AK8PuppiTaggedjet2 = book<TH1F>("N_daughters_AK8PuppiTaggedjet2", "N_{daughters}^{AK8Puppi Tagged jet 2}", 51, -0.5, 50.5);
  N_daughters_AK8PuppiTaggedjet3 = book<TH1F>("N_daughters_AK8PuppiTaggedjet3", "N_{daughters}^{AK8Puppi Tagged jet 3}", 51, -0.5, 50.5);
  dRmin_HOTVR_AK8PuppiTaggedjet  = book<TH1F>("dRmin_HOTVR_AK8PuppiTaggedjet", "#DeltaR_{min}(AK8Puppi Tagged jet, HOTVR jet)", 60, 0, 3);
  dRmin_HOTVR_AK8PuppiTaggedjet1 = book<TH1F>("dRmin_HOTVR_AK8PuppiTaggedjet1", "#DeltaR_{min}(AK8Puppi Tagged jet 1, HOTVR jet)", 60, 0, 3);
  dRmin_HOTVR_AK8PuppiTaggedjet2 = book<TH1F>("dRmin_HOTVR_AK8PuppiTaggedjet2", "#DeltaR_{min}(AK8Puppi Tagged jet 2, HOTVR jet)", 60, 0, 3);
  dRmin_HOTVR_AK8PuppiTaggedjet3 = book<TH1F>("dRmin_HOTVR_AK8PuppiTaggedjet3", "#DeltaR_{min}(AK8Puppi Tagged jet 3, HOTVR jet)", 60, 0, 3);
  dRmin_mu_AK8PuppiTaggedjet     = book<TH1F>("dRmin_mu_AK8PuppiTaggedjet", "#DeltaR_{min}(AK8Puppi Tagged jet, #mu)", 60, 0, 3);
  dRmin_mu_AK8PuppiTaggedjet1    = book<TH1F>("dRmin_mu_AK8PuppiTaggedjet1", "#DeltaR_{min}(AK8Puppi Tagged jet 1, #mu)", 60, 0, 3);
  dRmin_mu_AK8PuppiTaggedjet2    = book<TH1F>("dRmin_mu_AK8PuppiTaggedjet2", "#DeltaR_{min}(AK8Puppi Tagged jet 2, #mu)", 60, 0, 3);
  dRmin_mu_AK8PuppiTaggedjet3    = book<TH1F>("dRmin_mu_AK8PuppiTaggedjet3", "#DeltaR_{min}(AK8Puppi Tagged jet 3, #mu)", 60, 0, 3);
  tau1_AK8PuppiTaggedjet         = book<TH1F>("tau1_AK8PuppiTaggedjet", "#tau_{1}^{AK8Puppi Tagged jets}", 24, 0, 1.2);
  tau1_AK8PuppiTaggedjet1        = book<TH1F>("tau1_AK8PuppiTaggedjet1", "#tau_{1}^{AK8Puppi Tagged jet 1}", 24, 0, 1.2);
  tau1_AK8PuppiTaggedjet2        = book<TH1F>("tau1_AK8PuppiTaggedjet2", "#tau_{1}^{AK8Puppi Tagged jet 2}", 24, 0, 1.2);
  tau1_AK8PuppiTaggedjet3        = book<TH1F>("tau1_AK8PuppiTaggedjet3", "#tau_{1}^{AK8Puppi Tagged jet 3}", 24, 0, 1.2);
  tau2_AK8PuppiTaggedjet         = book<TH1F>("tau2_AK8PuppiTaggedjet", "#tau_{2}^{AK8Puppi Tagged jets}", 24, 0, 1.2);
  tau2_AK8PuppiTaggedjet1        = book<TH1F>("tau2_AK8PuppiTaggedjet1", "#tau_{2}^{AK8Puppi Tagged jet 1}", 24, 0, 1.2);
  tau2_AK8PuppiTaggedjet2        = book<TH1F>("tau2_AK8PuppiTaggedjet2", "#tau_{2}^{AK8Puppi Tagged jet 2}", 24, 0, 1.2);
  tau2_AK8PuppiTaggedjet3        = book<TH1F>("tau2_AK8PuppiTaggedjet3", "#tau_{2}^{AK8Puppi Tagged jet 3}", 24, 0, 1.2);
  tau3_AK8PuppiTaggedjet         = book<TH1F>("tau3_AK8PuppiTaggedjet", "#tau_{3}^{AK8Puppi Tagged jets}", 24, 0, 1.2);
  tau3_AK8PuppiTaggedjet1        = book<TH1F>("tau3_AK8PuppiTaggedjet1", "#tau_{3}^{AK8Puppi Tagged jet 1}", 24, 0, 1.2);
  tau3_AK8PuppiTaggedjet2        = book<TH1F>("tau3_AK8PuppiTaggedjet2", "#tau_{3}^{AK8Puppi Tagged jet 2}", 24, 0, 1.2);
  tau3_AK8PuppiTaggedjet3        = book<TH1F>("tau3_AK8PuppiTaggedjet3", "#tau_{3}^{AK8Puppi Tagged jet 3}", 24, 0, 1.2);
  tau21_AK8PuppiTaggedjet        = book<TH1F>("tau21_AK8PuppiTaggedjet", "#tau_{2/1}^{AK8Puppi Tagged jets}", 24, 0, 1.2);
  tau21_AK8PuppiTaggedjet1       = book<TH1F>("tau21_AK8PuppiTaggedjet1", "#tau_{2/1}^{AK8Puppi Tagged jet 1}", 24, 0, 1.2);
  tau21_AK8PuppiTaggedjet2       = book<TH1F>("tau21_AK8PuppiTaggedjet2", "#tau_{2/1}^{AK8Puppi Tagged jet 2}", 24, 0, 1.2);
  tau21_AK8PuppiTaggedjet3       = book<TH1F>("tau21_AK8PuppiTaggedjet3", "#tau_{2/1}^{AK8Puppi Tagged jet 3}", 24, 0, 1.2);
  tau32_AK8PuppiTaggedjet        = book<TH1F>("tau32_AK8PuppiTaggedjet", "#tau_{3/2}^{AK8Puppi Tagged jets}", 24, 0, 1.2);
  tau32_AK8PuppiTaggedjet1       = book<TH1F>("tau32_AK8PuppiTaggedjet1", "#tau_{3/2}^{AK8Puppi Tagged jet 1}", 24, 0, 1.2);
  tau32_AK8PuppiTaggedjet2       = book<TH1F>("tau32_AK8PuppiTaggedjet2", "#tau_{3/2}^{AK8Puppi Tagged jet 2}", 24, 0, 1.2);
  tau32_AK8PuppiTaggedjet3       = book<TH1F>("tau32_AK8PuppiTaggedjet3", "#tau_{3/2}^{AK8Puppi  Tagged jet 3}", 24, 0, 1.2);

  deepjet_topscore_Taggedjet  = book<TH1F>("deepjet_topscore_Taggedjet", "DeepJet top score all TaggedAK8 jets", 20, 0, 1);
  deepjet_topscore_Taggedjet1 = book<TH1F>("deepjet_topscore_Taggedjet1", "DeepJet top score TaggedAK8 jet 1", 20, 0, 1);
  deepjet_topscore_Taggedjet2 = book<TH1F>("deepjet_topscore_Taggedjet2", "DeepJet top score TaggedAK8 jet 2", 20, 0, 1);
  deepjet_topscore_Taggedjet3 = book<TH1F>("deepjet_topscore_Taggedjet3", "DeepJet top score TaggedAK8 jet 3", 20, 0, 1);

  massdecordeepjet_topscore_Taggedjet  = book<TH1F>("massdecordeepjet_topscore_Taggedjet", "Mass decor deepjet top score all TaggedAK8 jets", 20, 0, 1);
  massdecordeepjet_topscore_Taggedjet1 = book<TH1F>("massdecordeepjet_topscore_Taggedjet1", "Mass decor deepjet top score TaggedAK8 jet 1", 20, 0, 1);
  massdecordeepjet_topscore_Taggedjet2 = book<TH1F>("massdecordeepjet_topscore_Taggedjet2", "Mass decor deepjet top score TaggedAK8 jet 2", 20, 0, 1);
  massdecordeepjet_topscore_Taggedjet3 = book<TH1F>("massdecordeepjet_topscore_Taggedjet3", "Mass decor deepjet top score TaggedAK8 jet 3", 20, 0, 1);

  deepjet_TvsQCD_Taggedjet  = book<TH1F>("deepjet_TvsQCD_Taggedjet", "DeepJet  TvsQCD all TaggedAK8 jets", 20, 0, 1);
  deepjet_TvsQCD_Taggedjet1 = book<TH1F>("deepjet_TvsQCD_Taggedjet1", "DeepJet TvsQCD TaggedAK8 jet 1", 20, 0, 1);
  deepjet_TvsQCD_Taggedjet2 = book<TH1F>("deepjet_TvsQCD_Taggedjet2", "DeepJet TvsQCD TaggedAK8 jet 2", 20, 0, 1);
  deepjet_TvsQCD_Taggedjet3 = book<TH1F>("deepjet_TvsQCD_Taggedjet3", "DeepJet TvsQCD TaggedAK8 jet 3", 20, 0, 1);

  massdecordeepjet_TvsQCD_Taggedjet  = book<TH1F>("massdecordeepjet_TvsQCD_Taggedjet", "Mass decor deepjet  TvsQCD all TaggedAK8 jets", 20, 0, 1);
  massdecordeepjet_TvsQCD_Taggedjet1 = book<TH1F>("massdecordeepjet_TvsQCD_Taggedjet1", "Mass decor deepjet TvsQCD TaggedAK8 jet 1", 20, 0, 1);
  massdecordeepjet_TvsQCD_Taggedjet2 = book<TH1F>("massdecordeepjet_TvsQCD_Taggedjet2", "Mass decor deepjet TvsQCD TaggedAK8 jet 2", 20, 0, 1);
  massdecordeepjet_TvsQCD_Taggedjet3 = book<TH1F>("massdecordeepjet_TvsQCD_Taggedjet3", "Mass decor deepjet TvsQCD TaggedAK8 jet 3", 20, 0, 1);

  // general
  NPV           = book<TH1F>("NPV", "number of primary vertices", 91, -0.50, 90.5);
  MET           = book<TH1F>("MET", "missing E_{T} [GeV]", 50, 0, 7000);
  MET_rebin     = book<TH1F>("MET_rebin", "missing E_{T} [GeV]", 45, 0, 900);
  MET_rebin2    = book<TH1F>("MET_rebin2", "missing E_{T} [GeV]", 30, 0, 1500);
  MET_rebin3    = book<TH1F>("MET_rebin3", "missing E_{T} [GeV]", 15, 0, 1500);
  ST            = book<TH1F>("ST", "S_{T} [GeV]", 50, 0, 7000);
  ST_rebin      = book<TH1F>("ST_rebin", "S_{T} [GeV]", 200, 0, 5000);
  ST_rebin2     = book<TH1F>("ST_rebin2", "S_{T} [GeV]", 100, 0, 5000);
  ST_rebin3     = book<TH1F>("ST_rebin3", "S_{T} [GeV]", 50, 0, 5000);
  STjets        = book<TH1F>("STjets", "S_{T}^{jets} [GeV]", 50, 0, 7000);
  STjets_rebin  = book<TH1F>("STjets_rebin", "S_{T}^{jets} [GeV]", 200, 0, 5000);
  STjets_rebin2 = book<TH1F>("STjets_rebin2", "S_{T}^{jets} [GeV]", 100, 0, 5000);
  STjets_rebin3 = book<TH1F>("STjets_rebin3", "S_{T}^{jets} [GeV]", 50, 0, 5000);
  STlep         = book<TH1F>("STlep", "S_{T}^{lep} [GeV]", 50, 0, 7000);
  STlep_rebin   = book<TH1F>("STlep_rebin", "S_{T}^{lep} [GeV]", 45, 0, 900);
  STlep_rebin2  = book<TH1F>("STlep_rebin2", "S_{T}^{lep} [GeV]", 30, 0, 1500);
  STlep_rebin3  = book<TH1F>("STlep_rebin3", "S_{T}^{lep} [GeV]", 15, 0, 1500);

  // Zprime reconstruction
  toplep_pt         = book<TH1F>("toplep_pt", "p_{T}^{t,lep} [GeV]", 70, 0, 7000);
  toplep_eta        = book<TH1F>("toplep_eta", "#eta^{t,lep}", 60, -3.0, 3.0);
  toplep_phi        = book<TH1F>("toplep_phi", "#phi^{t,lep}", 70, -3.5, 3.5);
  toplep_m          = book<TH1F>("toplep_m", "m^{t,lep} [GeV]", 70, 0, 7000);
  tophad_pt         = book<TH1F>("tophad_pt", "p_{T}^{t,had} [GeV]", 70, 0, 7000);
  tophad_eta        = book<TH1F>("tophad_eta", "#eta^{t,had}", 60, -3.0, 3.0);
  tophad_phi        = book<TH1F>("tophad_phi", "#phi^{t,had}", 70, -3.5, 3.5);
  tophad_m          = book<TH1F>("tophad_m", "m^{t,had} [GeV]", 70, 0, 7000);
  ditop_mass        = book<TH1F>("ditop_mass", "m_{t#bar{t}} [GeV]", 70, 0, 7000);
  ditop_absDeltaPhi = book<TH1F>("ditop_absDeltaPhi", "|#Delta#phi(t,#bar{t})|", 35, 0, 3.5);
  ditop_deltaEta    = book<TH1F>("ditop_deltaEta", "#Delta#eta(t,#bar{t})", 120, -6.0, 6.0);
  ditop_absDeltaEta = book<TH1F>("ditop_absDeltaEta", "|#Delta#eta(t,#bar{t})|", 60, 0, 6.0);
  ditop_deltaR      = book<TH1F>("ditop_deltaR", "#DeltaR(t,#bar{t})", 100, 0, 10.0);
  
  // DeltaY
  DeltaY_reco            = book<TH1F>("DeltaY_reco", "#Delta Y_{(t,#bar{t})}",2,-2.5,2.5);
  DeltaY_reco_high       = book<TH1F>("DeltaY_reco", "#Delta Y_{(t,#bar{t})}",2,-2.5,2.5);
  DeltaY_reco_low        = book<TH1F>("DeltaY_reco", "#Delta Y_{(t,#bar{t})}",2,-2.5,2.5);
  DeltaY_reco_d1         = book<TH1F>("DeltaY_reco_d1", "#Delta Y_{(t,#bar{t})}",2,-2.5,2.5);
  DeltaY_reco_d2         = book<TH1F>("DeltaY_reco_d2", "#Delta Y_{(t,#bar{t})}",2,-2.5,2.5);
  DeltaY_reco_s1         = book<TH1F>("DeltaY_reco_s1", "#Delta Y_{(t,#bar{t})}",2,-2.5,2.5);
  DeltaY_reco_s2         = book<TH1F>("DeltaY_reco_s2", "#Delta Y_{(t,#bar{t})}",2,-2.5,2.5);
   //SpinCorr
  Sigma_phi             = book<TH1F>("Sigma_phi", "#Sigma #phi ",16,-3.2,3.2);
  Sigma_phi_high        = book<TH1F>("Sigma_phi_high", "#Sigma #phi ",16,-3.2,3.2);
  Sigma_phi_low         = book<TH1F>("Sigma_phi_low", "#Sigma #phi ",16,-3.2,3.2);
  Sigma_phi_1           = book<TH1F>("Sigma_phi_1", "#Sigma #phi ",16,-3.2,3.2);
  Sigma_phi_2           = book<TH1F>("Sigma_phi_2", "#Sigma #phi ",16,-3.2,3.2);

  Delta_phi             = book<TH1F>("Delta_phi", "#Delta #phi ",16,-3.2,3.2);
  Delta_phi_high        = book<TH1F>("Delta_phi_high", "#Delta #phi ",16,-3.2,3.2);
  Delta_phi_low         = book<TH1F>("Delta_phi_low", "#Delta #phi ",16,-3.2,3.2);
  Delta_phi_1           = book<TH1F>("Delta_phi_1", "#Delta #phi ",16,-3.2,3.2);
  Delta_phi_2           = book<TH1F>("Delta_phi_2", "#Delta #phi ",16,-3.2,3.2);


  DeltaY_gen            = book<TH1F>("DeltaY_gen", "#Delta Y_{(t,#bar{t})}",2,-2.5,2.5);
  // DeltaY_reco_0_500     = book<TH1F>("DeltaY_reco_0_500", "#Delta Y_{(t,#bar{t})} 0<Mtt<500",2,-2.5,2.5);
  // DeltaY_reco_500_750   = book<TH1F>("DeltaY_reco_500_750", "#Delta Y_{(t,#bar{t})} 500<Mtt<750",2,-2.5,2.5);
  // DeltaY_reco_750_1000  = book<TH1F>("DeltaY_reco_750_1000", "#Delta Y_{(t,#bar{t})} 750<Mtt<1000",2,-2.5,2.5);
  // DeltaY_reco_1000_1500 = book<TH1F>("DeltaY_reco_1000_1500", "#Delta Y_{(t,#bar{t})} 1000<Mtt<1500",2,-2.5,2.5);
  // DeltaY_reco_1500Inf   = book<TH1F>("DeltaY_reco_1500Inf", "#Delta Y_{(t,#bar{t})} Mtt>1500",2,-2.5,2.5);
  DeltaY_reco_best_plot = book<TH1F>("DeltaY_reco_best_plot", "#Delta Y_{(t,#bar{t}) Reco best}",2,-2.5,2.5);
  DeltaY_gen_best_plot  = book<TH1F>("DeltaY_gen_best_plot", "#Delta Y_{(t,#bar{t}) Gen Best}",2,-2.5,2.5);
  DeltaY_notMatched     = book<TH1F>("DeltaY_notMatched", "Reco Events Nor MATCHED",1,0,2);
  
 
  
  
  vector<float> bins_Zprime4 = {0,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200,3400,3600,3800,4000,4400,4800,5200,5600,6000,6100};
  vector<float> bins_Zprime5 = {0,200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3300,3600,3900,4200,4500,5000,5100};
  vector<float> bins_Zprime6 = {0,200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3300,3600,3900,4200,4500,5000};
  vector<float> bins_Zprime7 = {0,400,600,800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,5000};
  M_Zprime                 = book<TH1F>("M_Zprime", "M_{t#bar{t}} [GeV]", 280, 0, 7000);
  M_Zprime_rebin           = book<TH1F>("M_Zprime_rebin", "M_{t#bar{t}} [GeV]", 140, 0, 7000);
  M_Zprime_rebin2          = book<TH1F>("M_Zprime_rebin2", "M_{t#bar{t}} [GeV]", 70, 0, 7000);
  M_Zprime_rebin3          = book<TH1F>("M_Zprime_rebin3", "M_{t#bar{t}} [GeV]", 35, 0, 7000);
  M_Zprime_rebin4          = book<TH1F>("M_Zprime_rebin4", "M_{t#bar{t}} [GeV]", bins_Zprime4.size()-1, &bins_Zprime4[0]);
  M_Zprime_rebin5          = book<TH1F>("M_Zprime_rebin5", "M_{t#bar{t}} [GeV]", bins_Zprime5.size()-1, &bins_Zprime5[0]);
  M_Zprime_rebin6          = book<TH1F>("M_Zprime_rebin6", "M_{t#bar{t}} [GeV]", bins_Zprime6.size()-1, &bins_Zprime6[0]);
  M_Zprime_rebin7          = book<TH1F>("M_Zprime_rebin7", "M_{t#bar{t}} [GeV]", bins_Zprime7.size()-1, &bins_Zprime7[0]);
  M_Zprime_rebin8          = book<TH1F>("M_Zprime_rebin8", "M_{t#bar{t}} [GeV]", 400, 0, 10000);
  M_tophad                 = book<TH1F>("M_tophad", "M_{t}^{had} [GeV]", 70, 0, 700);
  M_toplep                 = book<TH1F>("M_toplep", "M_{t}^{lep} [GeV]", 70, 0, 700);
  chi2_Zprime              = book<TH1F>("chi2_Zprime", "#chi^{2}", 300, 0, 600);
  chi2_Zprime_rebin        = book<TH1F>("chi2_Zprime_rebin", "#chi^{2}", 20, 0, 100);
  chi2_Zprime_rebin2       = book<TH1F>("chi2_Zprime_rebin2", "#chi^{2}", 15, 0, 30);
  M_Zprime_ak4             = book<TH1F>("M_Zprime_ak4", "M_{t#bar{t}} (AK4 reconstruction) [GeV]", 280, 0, 7000);
  M_Zprime_ak4_rebin       = book<TH1F>("M_Zprime_ak4_rebin", "M_{t#bar{t}} (AK4 reconstruction) [GeV]", 140, 0, 7000);
  M_Zprime_ak4_rebin2      = book<TH1F>("M_Zprime_ak4_rebin2", "M_{t#bar{t}} (AK4 reconstruction) [GeV]", 70, 0, 7000);
  M_Zprime_ak4_rebin3      = book<TH1F>("M_Zprime_ak4_rebin3", "M_{t#bar{t}} (AK4 reconstruction) [GeV]", 35, 0, 7000);
  M_Zprime_ak4_rebin4      = book<TH1F>("M_Zprime_ak4_rebin4", "M_{t#bar{t}} (AK4 reconstruction [GeV]", bins_Zprime4.size()-1, &bins_Zprime4[0]);
  M_Zprime_ak4_rebin5      = book<TH1F>("M_Zprime_ak4_rebin5", "M_{t#bar{t}} (AK4 reconstruction [GeV]", bins_Zprime5.size()-1, &bins_Zprime5[0]);
  M_Zprime_ak4_rebin6      = book<TH1F>("M_Zprime_ak4_rebin6", "M_{t#bar{t}} (AK4 reconstruction [GeV]", bins_Zprime6.size()-1, &bins_Zprime6[0]);
  M_Zprime_ak4_rebin7      = book<TH1F>("M_Zprime_ak4_rebin7", "M_{t#bar{t}} (AK4 reconstruction [GeV]", bins_Zprime7.size()-1, &bins_Zprime7[0]);
  M_tophad_ak4             = book<TH1F>("M_tophad_ak4", "M_{t}^{had, AK4} (AK4 reconstruction) [GeV]", 70, 0, 700);
  M_toplep_ak4             = book<TH1F>("M_toplep_ak4", "M_{t}^{lep, AK4} (AK4 reconstruction) [GeV]", 70, 0, 700);
  chi2_Zprime_ak4          = book<TH1F>("chi2_Zprime_ak4", "#chi^{2} (AK4 reconstruction)", 50, 0, 500);
  chi2_Zprime_ak4_rebin    = book<TH1F>("chi2_Zprime_ak4_rebin", "#chi^{2} (AK4 reconstruction)", 20, 0, 100);
  chi2_Zprime_ak4_rebin2   = book<TH1F>("chi2_Zprime_ak4_rebin2", "#chi^{2} (AK4 reconstruction)", 15, 0, 30);
  M_Zprime_ttag            = book<TH1F>("M_Zprime_ttag", "M_{t#bar{t}} (t-tag reconstruction) [GeV]", 280, 0, 7000);
  M_Zprime_ttag_rebin      = book<TH1F>("M_Zprime_ttag_rebin", "M_{t#bar{t}} (t-tag reconstruction) [GeV]", 140, 0, 7000);
  M_Zprime_ttag_rebin2     = book<TH1F>("M_Zprime_ttag_rebin2", "M_{t#bar{t}} (t-tag reconstruction) [GeV]", 70, 0, 7000);
  M_Zprime_ttag_rebin3     = book<TH1F>("M_Zprime_ttag_rebin3", "M_{t#bar{t}} (t-tag reconstruction) [GeV]", 35, 0, 7000);
  M_Zprime_ttag_rebin4     = book<TH1F>("M_Zprime_ttag_rebin4", "M_{t#bar{t}} (t-tag reconstruction [GeV]", bins_Zprime4.size()-1, &bins_Zprime4[0]);
  M_Zprime_ttag_rebin5     = book<TH1F>("M_Zprime_ttag_rebin5", "M_{t#bar{t}} (t-tag reconstruction [GeV]", bins_Zprime5.size()-1, &bins_Zprime5[0]);
  M_Zprime_ttag_rebin6     = book<TH1F>("M_Zprime_ttag_rebin6", "M_{t#bar{t}} (t-tag reconstruction [GeV]", bins_Zprime6.size()-1, &bins_Zprime6[0]);
  M_Zprime_ttag_rebin7     = book<TH1F>("M_Zprime_ttag_rebin7", "M_{t#bar{t}} (t-tag reconstruction [GeV]", bins_Zprime7.size()-1, &bins_Zprime7[0]);
  M_tophad_ttag            = book<TH1F>("M_tophad_ttag", "M_{t}^{had, top-tag} (t-tag reconstruction) [GeV]", 70, 0, 700);
  M_tophad_ttag_v4         = book<TH1F>("M_tophad_ttag_v4", "M_{t}^{had, top-tag} (t-tag reconstruction, v4) [GeV]", 70, 0, 700);
  M_tophad_ttag_M          = book<TH1F>("M_tophad_ttag_M", "M_{t}^{had, top-tag} (t-tag reconstruction, M) [GeV]", 70, 0, 700);
  M_toplep_ttag            = book<TH1F>("M_toplep_ttag", "M_{t}^{lep, top-tag} (t-tag reconstruction) [GeV]", 70, 0, 700);
  chi2_Zprime_ttag         = book<TH1F>("chi2_Zprime_ttag", "#chi^{2} (t-tag reconstruction)", 50, 0, 500);
  chi2_Zprime_ttag_rebin   = book<TH1F>("chi2_Zprime_ttag_rebin", "#chi^{2} (t-tag reconstruction)", 20, 0, 100);
  chi2_Zprime_ttag_rebin2  = book<TH1F>("chi2_Zprime_ttag_rebin2", "#chi^{2} (t-tag reconstruction)", 15, 0, 30);
  M_tophad_dr_ak4          = book<TH1F>("M_tophad_dr_ak4", "M_{t}^{had, AK4} (correctly matched) [GeV]", 700, 0, 700);
  M_toplep_dr_ak4          = book<TH1F>("M_toplep_dr_ak4", "M_{t}^{lep, AK4} (correctly matched) [GeV]", 700, 0, 700);
  M_tophad_dr_ttag         = book<TH1F>("M_tophad_dr_ttag", "M_{t}^{had, top-tag} (correctly matched) [GeV]", 700, 0, 700);
  M_toplep_dr_ttag         = book<TH1F>("M_toplep_dr_ttag", "M_{t}^{lep, top-tag} (correctly matched) [GeV]", 700, 0, 700);
  dr_discr_Zprime          = book<TH1F>("dr_discr_Zprime", "dR best hypothesis (correctly matched)", 30, 0, 3);
  M_Zprime_dr              = book<TH1F>("M_Zprime_dr", "M_{t#bar{t}} (correctly matched) [GeV]", 280, 0, 7000);
  M_Zprime_dr_rebin        = book<TH1F>("M_Zprime_dr_rebin", "M_{t#bar{t}} (correctly matched) [GeV]", 140, 0, 7000);
  M_Zprime_dr_rebin2       = book<TH1F>("M_Zprime_dr_rebin2", "M_{t#bar{t}} (correctly matched) [GeV]", 70, 0, 7000);
  M_Zprime_dr_rebin3       = book<TH1F>("M_Zprime_dr_rebin3", "M_{t#bar{t}} (correctly matched) [GeV]", 35, 0, 7000);

  // 2D histogram of mttbar x cos(theta*) to optimize binning
  mttbar_vs_costhetastar = book<TH2F>("mttbar_vs_costhetastar", "m_{t#bar{t}} vs cos(#theta*)", 20, -1, 1, 1000, 0, 10000);
  costhetastar_vs_mttbar = book<TH2F>("costhetastar_vs_mttbar", "cos(#theta*) vs m_{t#bar{t}}", 1000, 0, 10000, 20, -1, 1);

  response_matrix = book<TH2F>("response_matrix", "#Delta Y_{(t,#bar{t})}_reco ;#Delta Y_{(t,#bar{t})}_gen",  2, -2.5, 2.5, 2, -2.5, 2.5);
  // response_matrix->GetXaxis()->SetBinLabel(1, "Negative");
  // response_matrix->GetXaxis()->SetBinLabel(2, "Positive");
  // response_matrix->GetYaxis()->SetBinLabel(1, "Negative");
  // response_matrix->GetYaxis()->SetBinLabel(2, "Positive");

  // Sphericity tensor
  S11 = book<TH1F>("S11", "S_{11}", 50, 0, 1);
  S12 = book<TH1F>("S12", "S_{12}", 50, 0, 1);
  S13 = book<TH1F>("S13", "S_{13}", 50, 0, 1);
  S22 = book<TH1F>("S22", "S_{22}", 50, 0, 1);
  S23 = book<TH1F>("S23", "S_{23}", 50, 0, 1);
  S33 = book<TH1F>("S33", "S_{33}", 50, 0, 1);

  sum_event_weights = book<TH1F>("sum_event_weights", "counting experiment", 1, 0.5, 1.5);

  // define theta star angles wrt ttbar rest frame
  // hadtop_thetastar     = book<TH1F>("hadtop_thetastar", "hadtop #theta^*", 70, -3.5, 3.5);
  cos_hadtop_thetastar = book<TH1F>("cos_hadtop_thetastar", "cos(hadtop #theta^*)", 100, -1.0, 1.0);
  // leptop_thetastar     = book<TH1F>("leptop_thetastar", "leptop #theta^*", 70, -3.5, 3.5);
  cos_leptop_thetastar = book<TH1F>("cos_leptop_thetastar", "cos(leptop #theta^*)", 100, -1.0, 1.0);
  abs_cos_leptop_thetastar = book<TH1F>("abs_cos_leptop_thetastar", "|cos(leptop #theta^*)|", 50, 0., 1.0);

  // mttbar histograms from TOP-20-001: https://www.hepdata.net/record/ins1901295
  vector<float> bins_ditopmass_Fig19 = {250,400,480,560,640,720,800,900,1000,1150,1300,1500,1700,2000,2300,3500};
  vector<float> bins_ditopmass_Fig25 = {250,420,520,620,800,1000,3500};
  TOP_20_001_ditopmass_Fig19 = book<TH1F>("TOP_20_001_ditopmass_Fig19","m_{t#bar{t} [GeV]}",bins_ditopmass_Fig19.size()-1,&bins_ditopmass_Fig19[0]);
  TOP_20_001_ditopmass_Fig25 = book<TH1F>("TOP_20_001_ditopmass_Fig25","m_{t#bar{t} [GeV]}",bins_ditopmass_Fig25.size()-1,&bins_ditopmass_Fig25[0]);

  // 2D ditribution NJets/HT to extract custom btag SF
  N_Jets_vs_HT  = book<TH2F>("N_Jets_vs_HT", "N_Jets_vs_HT", 21, 0., 21., 50, 0., 7000.);

  if(NN){
    // NN Hists
    NN_Mu_pt            = book<TH1F>("NN_Mu_pt", "NN_Mu_pt", 50, 0, 1000);
    NN_Mu_eta           = book<TH1F>("NN_Mu_eta", "NN_Mu_eta", 50, -2.5, 2.5);
    NN_Mu_phi           = book<TH1F>("NN_Mu_phi", "NN_Mu_phi", 35, -3.5, 3.5);
    NN_Mu_E             = book<TH1F>("NN_Mu_E", "NN_Mu_E", 150, 0, 3000);
    NN_Ele_pt           = book<TH1F>("NN_Ele_pt", "NN_Ele_pt", 50, 0, 1000);
    NN_Ele_eta          = book<TH1F>("NN_Ele_eta", "NN_Ele_eta", 50, -2.5, 2.5);
    NN_Ele_phi          = book<TH1F>("NN_Ele_phi", "NN_Ele_phi", 35, -3.5, 3.5);
    NN_Ele_E            = book<TH1F>("NN_Ele_E", "NN_Ele_E", 150, 0, 3000);
    NN_MET_pt           = book<TH1F>("NN_MET_pt", "NN_MET_pt", 150, 0, 1500);
    NN_MET_phi          = book<TH1F>("NN_MET_phi", "NN_MET_phi", 35, -3.5, 3.5);
    NN_N_Ak4            = book<TH1F>("NN_N_Ak4", "NN_N_Ak4", 20, 0, 20);
    NN_Ak4_j1_pt        = book<TH1F>("NN_Ak4_j1_pt", "NN_Ak4_j1_pt", 150, 0, 3000);
    NN_Ak4_j1_eta       = book<TH1F>("NN_Ak4_j1_eta", "NN_Ak4_j1_eta", 50, -2.5, 2.5);
    NN_Ak4_j1_phi       = book<TH1F>("NN_Ak4_j1_phi", "NN_Ak4_j1_phi", 35, -3.5, 3.5);
    NN_Ak4_j1_E         = book<TH1F>("NN_Ak4_j1_E", "NN_Ak4_j1_E", 100, 0, 5000);
    NN_Ak4_j1_m         = book<TH1F>("NN_Ak4_j1_m", "NN_Ak4_j1_m", 50, 0, 300);
    NN_Ak4_j1_btag      = book<TH1F>("NN_Ak4_j1_btag", "NN_Ak4_j1_btag", 50, 0, 1);
    NN_Ak4_j2_pt        = book<TH1F>("NN_Ak4_j2_pt", "NN_Ak4_j2_pt", 150, 0, 3000);
    NN_Ak4_j2_eta       = book<TH1F>("NN_Ak4_j2_eta", "NN_Ak4_j2_eta", 50, -2.5, 2.5);
    NN_Ak4_j2_phi       = book<TH1F>("NN_Ak4_j2_phi", "NN_Ak4_j2_phi", 35, -3.5, 3.5);
    NN_Ak4_j2_E         = book<TH1F>("NN_Ak4_j2_E", "NN_Ak4_j2_E", 100, 0, 5000);
    NN_Ak4_j2_m         = book<TH1F>("NN_Ak4_j2_m", "NN_Ak4_j2_m", 50, 0, 300);
    NN_Ak4_j2_btag      = book<TH1F>("NN_Ak4_j2_btag", "NN_Ak4_j2_btag", 50, 0, 1);
    NN_Ak4_j3_pt        = book<TH1F>("NN_Ak4_j3_pt", "NN_Ak4_j3_pt", 150, 0, 3000);
    NN_Ak4_j3_eta       = book<TH1F>("NN_Ak4_j3_eta", "NN_Ak4_j3_eta", 50, -2.5, 2.5);
    NN_Ak4_j3_phi       = book<TH1F>("NN_Ak4_j3_phi", "NN_Ak4_j3_phi", 35, -3.5, 3.5);
    NN_Ak4_j3_E         = book<TH1F>("NN_Ak4_j3_E", "NN_Ak4_j3_E", 100, 0, 5000);
    NN_Ak4_j3_m         = book<TH1F>("NN_Ak4_j3_m", "NN_Ak4_j3_m", 50, 0, 300);
    NN_Ak4_j3_btag      = book<TH1F>("NN_Ak4_j3_btag", "NN_Ak4_j3_btag", 50, 0, 1);
    NN_Ak4_j4_pt        = book<TH1F>("NN_Ak4_j4_pt", "NN_Ak4_j4_pt", 150, 0, 3000);
    NN_Ak4_j4_eta       = book<TH1F>("NN_Ak4_j4_eta", "NN_Ak4_j4_eta", 50, -2.5, 2.5);
    NN_Ak4_j4_phi       = book<TH1F>("NN_Ak4_j4_phi", "NN_Ak4_j4_phi", 35, -3.5, 3.5);
    NN_Ak4_j4_E         = book<TH1F>("NN_Ak4_j4_E", "NN_Ak4_j4_E", 100, 0, 5000);
    NN_Ak4_j4_m         = book<TH1F>("NN_Ak4_j4_m", "NN_Ak4_j4_m", 50, 0, 300);
    NN_Ak4_j4_btag      = book<TH1F>("NN_Ak4_j4_btag", "NN_Ak4_j4_btag", 50, 0, 1);
    NN_Ak4_j5_pt        = book<TH1F>("NN_Ak4_j5_pt", "NN_Ak4_j5_pt", 150, 0, 3000);
    NN_Ak4_j5_eta       = book<TH1F>("NN_Ak4_j5_eta", "NN_Ak4_j5_eta", 50, -2.5, 2.5);
    NN_Ak4_j5_phi       = book<TH1F>("NN_Ak4_j5_phi", "NN_Ak4_j5_phi", 35, -3.5, 3.5);
    NN_Ak4_j5_E         = book<TH1F>("NN_Ak4_j5_E", "NN_Ak4_j5_E", 100, 0, 5000);
    NN_Ak4_j5_m         = book<TH1F>("NN_Ak4_j5_m", "NN_Ak4_j5_m", 50, 0, 300);
    NN_Ak4_j5_btag      = book<TH1F>("NN_Ak4_j5_btag", "NN_Ak4_j5_btag", 50, 0, 1);
    NN_Ak4_j6_pt        = book<TH1F>("NN_Ak4_j6_pt", "NN_Ak4_j6_pt", 150, 0, 3000);
    NN_Ak4_j6_eta       = book<TH1F>("NN_Ak4_j6_eta", "NN_Ak4_j6_eta", 50, -2.5, 2.5);
    NN_Ak4_j6_phi       = book<TH1F>("NN_Ak4_j6_phi", "NN_Ak4_j6_phi", 35, -3.5, 3.5);
    NN_Ak4_j6_E         = book<TH1F>("NN_Ak4_j6_E", "NN_Ak4_j6_E", 100, 0, 5000);
    NN_Ak4_j6_m         = book<TH1F>("NN_Ak4_j6_m", "NN_Ak4_j6_m", 50, 0, 300);
    NN_Ak4_j6_btag      = book<TH1F>("NN_Ak4_j6_btag", "NN_Ak4_j6_btag", 50, 0, 1);
    NN_M_tt_weighted    = book<TH1F>("NN_M_tt_weighted", "NN_M_tt_weighted", 100, 0, 14000);
    NN_M_tt_notweighted = book<TH1F>("NN_M_tt_notweighted", "NN_M_tt_notweighted", 100, 0, 14000);
    NN_chi2             = book<TH1F>("NN_chi2", "NN_chi2", 100, 0, 100);
    NN_N_HOTVR          = book<TH1F>("NN_N_HOTVR", "NN_N_HOTVR", 20, 0, 20);
    NN_HOTVR_j1_pt      = book<TH1F>("NN_HOTVR_j1_pt", "NN_HOTVR_j1_pt", 150, 0, 3000);
    NN_HOTVR_j1_eta     = book<TH1F>("NN_HOTVR_j1_eta", "NN_HOTVR_j1_eta", 50, -2.5, 2.5);
    NN_HOTVR_j1_phi     = book<TH1F>("NN_HOTVR_j1_phi", "NN_HOTVR_j1_phi", 35, -3.5, 3.5);
    NN_HOTVR_j1_E       = book<TH1F>("NN_HOTVR_j1_E", "NN_HOTVR_j1_E", 100, 0, 5000);
    NN_HOTVR_j1_mSD     = book<TH1F>("NN_HOTVR_j1_mSD", "NN_HOTVR_j1_mSD", 50, 0, 500);
    NN_HOTVR_j1_tau21   = book<TH1F>("NN_HOTVR_j1_tau21", "NN_HOTVR_j1_tau21", 24, 0, 1.2);
    NN_HOTVR_j1_tau32   = book<TH1F>("NN_HOTVR_j1_tau32", "NN_HOTVR_j1_tau32", 24, 0, 1.2);
    NN_HOTVR_j2_pt      = book<TH1F>("NN_HOTVR_j2_pt", "NN_HOTVR_j2_pt", 150, 0, 3000);
    NN_HOTVR_j2_eta     = book<TH1F>("NN_HOTVR_j2_eta", "NN_HOTVR_j2_eta", 50, -2.5, 2.5);
    NN_HOTVR_j2_phi     = book<TH1F>("NN_HOTVR_j2_phi", "NN_HOTVR_j2_phi", 35, -3.5, 3.5);
    NN_HOTVR_j2_E       = book<TH1F>("NN_HOTVR_j2_E", "NN_HOTVR_j2_E", 100, 0, 5000);
    NN_HOTVR_j2_mSD     = book<TH1F>("NN_HOTVR_j2_mSD", "NN_HOTVR_j2_mSD", 50, 0, 500);
    NN_HOTVR_j2_tau21   = book<TH1F>("NN_HOTVR_j2_tau21", "NN_HOTVR_j2_tau21", 24, 0, 1.2);
    NN_HOTVR_j2_tau32   = book<TH1F>("NN_HOTVR_j2_tau32", "NN_HOTVR_j2_tau32", 24, 0, 1.2);
    NN_HOTVR_j3_pt      = book<TH1F>("NN_HOTVR_j3_pt", "NN_HOTVR_j3_pt", 150, 0, 3000);
    NN_HOTVR_j3_eta     = book<TH1F>("NN_HOTVR_j3_eta", "NN_HOTVR_j3_eta", 50, -2.5, 2.5);
    NN_HOTVR_j3_phi     = book<TH1F>("NN_HOTVR_j3_phi", "NN_HOTVR_j3_phi", 35, -3.5, 3.5);
    NN_HOTVR_j3_E       = book<TH1F>("NN_HOTVR_j3_E", "NN_HOTVR_j3_E", 100, 0, 5000);
    NN_HOTVR_j3_mSD     = book<TH1F>("NN_HOTVR_j3_mSD", "NN_HOTVR_j3_mSD", 50, 0, 500);
    NN_HOTVR_j3_tau21   = book<TH1F>("NN_HOTVR_j3_tau21", "NN_HOTVR_j3_tau21", 24, 0, 1.2);
    NN_HOTVR_j3_tau32   = book<TH1F>("NN_HOTVR_j3_tau32", "NN_HOTVR_j3_tau32", 24, 0, 1.2);
    NN_N_Ak8            = book<TH1F>("NN_N_Ak8", "NN_N_Ak8", 20, 0, 20);
    NN_Ak8_j1_pt        = book<TH1F>("NN_Ak8_j1_pt", "NN_Ak8_j1_pt", 150, 0, 3000);
    NN_Ak8_j1_eta       = book<TH1F>("NN_Ak8_j1_eta", "NN_Ak8_j1_eta", 50, -2.5, 2.5);
    NN_Ak8_j1_phi       = book<TH1F>("NN_Ak8_j1_phi", "NN_Ak8_j1_phi", 35, -3.5, 3.5);
    NN_Ak8_j1_E         = book<TH1F>("NN_Ak8_j1_E", "NN_Ak8_j1_E", 100, 0, 5000);
    NN_Ak8_j1_mSD       = book<TH1F>("NN_Ak8_j1_mSD", "NN_Ak8_j1_mSD", 50, 0, 500);
    NN_Ak8_j1_tau21     = book<TH1F>("NN_Ak8_j1_tau21", "NN_Ak8_j1_tau21", 24, 0, 1.2);
    NN_Ak8_j1_tau32     = book<TH1F>("NN_Ak8_j1_tau32", "NN_Ak8_j1_tau32", 24, 0, 1.2);
    NN_Ak8_j1_ttag      = book<TH1F>("NN_Ak8_j1_ttag", "NN_Ak8_j1_ttag", 50, 0, 1);
    NN_Ak8_j2_pt        = book<TH1F>("NN_Ak8_j2_pt", "NN_Ak8_j2_pt", 150, 0, 3000);
    NN_Ak8_j2_eta       = book<TH1F>("NN_Ak8_j2_eta", "NN_Ak8_j2_eta", 50, -2.5, 2.5);
    NN_Ak8_j2_phi       = book<TH1F>("NN_Ak8_j2_phi", "NN_Ak8_j2_phi", 35, -3.5, 3.5);
    NN_Ak8_j2_E         = book<TH1F>("NN_Ak8_j2_E", "NN_Ak8_j2_E", 100, 0, 5000);
    NN_Ak8_j2_mSD       = book<TH1F>("NN_Ak8_j2_mSD", "NN_Ak8_j2_mSD", 50, 0, 500);
    NN_Ak8_j2_tau21     = book<TH1F>("NN_Ak8_j2_tau21", "NN_Ak8_j2_tau21", 24, 0, 1.2);
    NN_Ak8_j2_tau32     = book<TH1F>("NN_Ak8_j2_tau32", "NN_Ak8_j2_tau32", 24, 0, 1.2);
    NN_Ak8_j2_ttag      = book<TH1F>("NN_Ak8_j2_ttag", "NN_Ak8_j2_ttag", 50, 0, 1);
    NN_Ak8_j3_pt        = book<TH1F>("NN_Ak8_j3_pt", "NN_Ak8_j3_pt", 150, 0, 3000);
    NN_Ak8_j3_eta       = book<TH1F>("NN_Ak8_j3_eta", "NN_Ak8_j3_eta", 50, -2.5, 2.5);
    NN_Ak8_j3_phi       = book<TH1F>("NN_Ak8_j3_phi", "NN_Ak8_j3_phi", 35, -3.5, 3.5);
    NN_Ak8_j3_E         = book<TH1F>("NN_Ak8_j3_E", "NN_Ak8_j3_E", 100, 0, 5000);
    NN_Ak8_j3_mSD       = book<TH1F>("NN_Ak8_j3_mSD", "NN_Ak8_j3_mSD", 50, 0, 500);
    NN_Ak8_j3_tau21     = book<TH1F>("NN_Ak8_j3_tau21", "NN_Ak8_j3_tau21", 24, 0, 1.2);
    NN_Ak8_j3_tau32     = book<TH1F>("NN_Ak8_j3_tau32", "NN_Ak8_j3_tau32", 24, 0, 1.2);
    NN_Ak8_j3_ttag      = book<TH1F>("NN_Ak8_j3_ttag", "NN_Ak8_j3_ttag", 50, 0, 1);
    //NN
  }
}


void ZprimeSemiLeptonicHists::fill(const Event & event){

  double weight = event.weight;



  /*
  █      ██ ███████ ████████ ███████
  █      ██ ██         ██    ██
  █      ██ █████      ██    ███████
  █ ██   ██ ██         ██         ██
  █  █████  ███████    ██    ███████
  */


  vector<Jet>* jets = event.jets;
  int Njets = jets->size();
  N_jets->Fill(Njets, weight);

  // for(unsigned int i=0; i<jets->size(); i++){
  //   cout << "Jet pt: " << jets->at(i).pt() << endl
  // }
  // cout << "Jet1 pt: " << jets->at(0).pt() << endl;
  // cout << "Jet2 pt: " << jets->at(0).pt() << endl;
  // cout << "Jet3 pt: " << jets->at(0).pt() << endl;

  for(unsigned int i=0; i<jets->size(); i++){
    pt_jet->Fill(jets->at(i).pt(),weight);
    eta_jet->Fill(jets->at(i).eta(),weight);
    phi_jet->Fill(jets->at(i).phi(),weight);
    m_jet->Fill(jets->at(i).v4().M(),weight);
    deepjetbscore_jet->Fill(jets->at(i).btag_DeepJet(), weight);

    if(i==0){
      pt_jet1->Fill(jets->at(i).pt(),weight);
      eta_jet1->Fill(jets->at(i).eta(),weight);
      phi_jet1->Fill(jets->at(i).phi(),weight);
      m_jet1->Fill(jets->at(i).v4().M(),weight);
      deepjetbscore_jet1->Fill(jets->at(i).btag_DeepJet(), weight);
    }
    else if(i==1){
      pt_jet2->Fill(jets->at(i).pt(),weight);
      eta_jet2->Fill(jets->at(i).eta(),weight);
      phi_jet2->Fill(jets->at(i).phi(),weight);
      m_jet2->Fill(jets->at(i).v4().M(),weight);
      deepjetbscore_jet2->Fill(jets->at(i).btag_DeepJet(), weight);
    }
    else if(i==2){
      pt_jet3->Fill(jets->at(i).pt(),weight);
      eta_jet3->Fill(jets->at(i).eta(),weight);
      phi_jet3->Fill(jets->at(i).phi(),weight);
      m_jet3->Fill(jets->at(i).v4().M(),weight);
      deepjetbscore_jet3->Fill(jets->at(i).btag_DeepJet(), weight);
    }
  }


  int NbjetsDeepJet_loose = 0, NbjetsDeepJet_medium = 0, NbjetsDeepJet_tight = 0;
  DeepJetBTag BtagDeepJet_loose = DeepJetBTag(DeepJetBTag::WP_LOOSE);
  DeepJetBTag BtagDeepJet_medium = DeepJetBTag(DeepJetBTag::WP_MEDIUM);
  DeepJetBTag BtagDeepJet_tight = DeepJetBTag(DeepJetBTag::WP_TIGHT);

  for (unsigned int i =0; i<jets->size(); i++) {
    if(BtagDeepJet_loose(jets->at(i),event))  NbjetsDeepJet_loose++;
    if(BtagDeepJet_medium(jets->at(i),event)) NbjetsDeepJet_medium++;
    if(BtagDeepJet_tight(jets->at(i),event))  NbjetsDeepJet_tight++;
  }

  N_bJetsDeepJet_loose->Fill(NbjetsDeepJet_loose,weight);
  N_bJetsDeepJet_med->Fill(NbjetsDeepJet_medium,weight);
  N_bJetsDeepJet_tight->Fill(NbjetsDeepJet_tight,weight);

  /*
  █ ██   ██    █████  ██████ ██    ██   █████
  █ ██   ██  ██     ██  ██   ██    ██   ██  ██
  █ ███████  ██     ██  ██    ██  ██    █████
  █ ██   ██  ██     ██  ██     ████     ██  ██
  █ ██   ██    █████    ██      ██      ██   ██
  */

  if(ishotvr){
    vector<TopJet>* HOTVRjets = event.topjets;
    unsigned int NHOTVRjets = HOTVRjets->size();
    N_HOTVRjets->Fill(NHOTVRjets, weight);

    for(unsigned int i=0; i<NHOTVRjets; i++){
      double tau21 = HOTVRjets->at(i).tau2_groomed() / HOTVRjets->at(i).tau1_groomed();
      double tau32 = HOTVRjets->at(i).tau3_groomed() / HOTVRjets->at(i).tau2_groomed();

      // Distance to AK8
      double dRmin_Puppi = 99999;
      for(unsigned int j=0; j<event.toppuppijets->size(); j++){
        double dR = deltaR(HOTVRjets->at(i), event.toppuppijets->at(j));
        if(dR < dRmin_Puppi) dRmin_Puppi = dR;
      }

      // Distance to muons
      double dRmin_muon = 99999;
      for(unsigned int j=0; j<event.muons->size(); j++){
        double dR = deltaR(HOTVRjets->at(i), event.muons->at(j));
        if(dR < dRmin_muon) dRmin_muon = dR;
      }

      pt_HOTVRjet->Fill(HOTVRjets->at(i).pt(), weight);
      eta_HOTVRjet->Fill(HOTVRjets->at(i).eta(), weight);
      phi_HOTVRjet->Fill(HOTVRjets->at(i).phi(), weight);
      m_HOTVRjet->Fill(HOTVRjets->at(i).v4().M(), weight);
      dRmin_AK8Puppi_HOTVRjet->Fill(dRmin_Puppi, weight);
      dRmin_mu_HOTVRjet->Fill(dRmin_muon, weight);
      N_subjets_HOTVRjet->Fill(HOTVRjets->at(i).subjets().size(), weight);
      N_daughters_HOTVRjet->Fill(HOTVRjets->at(i).numberOfDaughters(), weight);
      tau1_HOTVRjet->Fill(HOTVRjets->at(i).tau1_groomed(), weight);
      tau2_HOTVRjet->Fill(HOTVRjets->at(i).tau2_groomed(), weight);
      tau3_HOTVRjet->Fill(HOTVRjets->at(i).tau3_groomed(), weight);
      tau21_HOTVRjet->Fill(tau21, weight);
      tau32_HOTVRjet->Fill(tau32, weight);

      if(i==0){
        pt_HOTVRjet1->Fill(HOTVRjets->at(i).pt(), weight);
        eta_HOTVRjet1->Fill(HOTVRjets->at(i).eta(), weight);
        phi_HOTVRjet1->Fill(HOTVRjets->at(i).phi(), weight);
        m_HOTVRjet1->Fill(HOTVRjets->at(i).v4().M(), weight);
        dRmin_AK8Puppi_HOTVRjet1->Fill(dRmin_Puppi, weight);
        dRmin_mu_HOTVRjet1->Fill(dRmin_muon, weight);
        N_subjets_HOTVRjet1->Fill(HOTVRjets->at(i).subjets().size(), weight);
        N_daughters_HOTVRjet1->Fill(HOTVRjets->at(i).numberOfDaughters(), weight);
        tau1_HOTVRjet1->Fill(HOTVRjets->at(i).tau1_groomed(), weight);
        tau2_HOTVRjet1->Fill(HOTVRjets->at(i).tau2_groomed(), weight);
        tau3_HOTVRjet1->Fill(HOTVRjets->at(i).tau3_groomed(), weight);
        tau21_HOTVRjet1->Fill(tau21, weight);
        tau32_HOTVRjet1->Fill(tau32, weight);
      }
      else if(i==1){
        pt_HOTVRjet2->Fill(HOTVRjets->at(i).pt(), weight);
        eta_HOTVRjet2->Fill(HOTVRjets->at(i).eta(), weight);
        phi_HOTVRjet2->Fill(HOTVRjets->at(i).phi(), weight);
        m_HOTVRjet2->Fill(HOTVRjets->at(i).v4().M(), weight);
        dRmin_AK8Puppi_HOTVRjet2->Fill(dRmin_Puppi, weight);
        dRmin_mu_HOTVRjet2->Fill(dRmin_muon, weight);
        N_subjets_HOTVRjet2->Fill(HOTVRjets->at(i).subjets().size(), weight);
        N_daughters_HOTVRjet2->Fill(HOTVRjets->at(i).numberOfDaughters(), weight);
        tau1_HOTVRjet2->Fill(HOTVRjets->at(i).tau1_groomed(), weight);
        tau2_HOTVRjet2->Fill(HOTVRjets->at(i).tau2_groomed(), weight);
        tau3_HOTVRjet2->Fill(HOTVRjets->at(i).tau3_groomed(), weight);
        tau21_HOTVRjet2->Fill(tau21, weight);
        tau32_HOTVRjet2->Fill(tau32, weight);
      }
      else if(i==2){
        pt_HOTVRjet3->Fill(HOTVRjets->at(i).pt(), weight);
        eta_HOTVRjet3->Fill(HOTVRjets->at(i).eta(), weight);
        phi_HOTVRjet3->Fill(HOTVRjets->at(i).phi(), weight);
        m_HOTVRjet3->Fill(HOTVRjets->at(i).v4().M(), weight);
        dRmin_AK8Puppi_HOTVRjet3->Fill(dRmin_Puppi, weight);
        dRmin_mu_HOTVRjet3->Fill(dRmin_muon, weight);
        N_subjets_HOTVRjet3->Fill(HOTVRjets->at(i).subjets().size(), weight);
        N_daughters_HOTVRjet3->Fill(HOTVRjets->at(i).numberOfDaughters(), weight);
        tau1_HOTVRjet3->Fill(HOTVRjets->at(i).tau1_groomed(), weight);
        tau2_HOTVRjet3->Fill(HOTVRjets->at(i).tau2_groomed(), weight);
        tau3_HOTVRjet3->Fill(HOTVRjets->at(i).tau3_groomed(), weight);
        tau21_HOTVRjet3->Fill(tau21, weight);
        tau32_HOTVRjet3->Fill(tau32, weight);
      }
    }

    // Top-tagged HOTVR jets
    vector<TopJet> HOTVRTopTags = event.get(h_AK8TopTags);
    int NHOTVRTaggedjets = 0;
    for(unsigned int i=0; i<HOTVRTopTags.size(); i++){
      NHOTVRTaggedjets++;

      double tau21 = HOTVRTopTags.at(i).tau2_groomed() / HOTVRTopTags.at(i).tau1_groomed();
      double tau32 = HOTVRTopTags.at(i).tau3_groomed() / HOTVRTopTags.at(i).tau2_groomed();

      // Distance to AK8 Puppi
      double dRmin_ak8 = 99999;
      for(unsigned int j=0; j<event.toppuppijets->size(); j++){
        double dR = deltaR(HOTVRTopTags.at(i), event.toppuppijets->at(j));
        if(dR < dRmin_ak8) dRmin_ak8 = dR;
      }

      // Distance to muons
      double dRmin_muon = 99999;
      for(unsigned int j=0; j<event.muons->size(); j++){
        double dR = deltaR(HOTVRTopTags.at(i), event.muons->at(j));
        if(dR < dRmin_muon) dRmin_muon = dR;
      }

      pt_HOTVRTaggedjet->Fill(HOTVRTopTags.at(i).pt(), weight);
      eta_HOTVRTaggedjet->Fill(HOTVRTopTags.at(i).eta(), weight);
      phi_HOTVRTaggedjet->Fill(HOTVRTopTags.at(i).phi(), weight);
      m_HOTVRTaggedjet->Fill(HOTVRTopTags.at(i).v4().M(), weight);
      dRmin_AK8_HOTVRTaggedjet->Fill(dRmin_ak8, weight);
      dRmin_mu_HOTVRTaggedjet->Fill(dRmin_muon, weight);
      N_subjets_HOTVRTaggedjet->Fill(HOTVRTopTags.at(i).subjets().size(), weight);
      N_daughters_HOTVRTaggedjet->Fill(HOTVRTopTags.at(i).numberOfDaughters(), weight);
      tau1_HOTVRTaggedjet->Fill(HOTVRTopTags.at(i).tau1_groomed(), weight);
      tau2_HOTVRTaggedjet->Fill(HOTVRTopTags.at(i).tau2_groomed(), weight);
      tau3_HOTVRTaggedjet->Fill(HOTVRTopTags.at(i).tau3_groomed(), weight);
      tau21_HOTVRTaggedjet->Fill(tau21, weight);
      tau32_HOTVRTaggedjet->Fill(tau32, weight);

      if(i==0){
        pt_HOTVRTaggedjet1->Fill(HOTVRTopTags.at(i).pt(), weight);
        eta_HOTVRTaggedjet1->Fill(HOTVRTopTags.at(i).eta(), weight);
        phi_HOTVRTaggedjet1->Fill(HOTVRTopTags.at(i).phi(), weight);
        m_HOTVRTaggedjet1->Fill(HOTVRTopTags.at(i).v4().M(), weight);
        dRmin_AK8_HOTVRTaggedjet1->Fill(dRmin_ak8, weight);
        dRmin_mu_HOTVRTaggedjet1->Fill(dRmin_muon, weight);
        N_subjets_HOTVRTaggedjet1->Fill(HOTVRTopTags.at(i).subjets().size(), weight);
        N_daughters_HOTVRTaggedjet1->Fill(HOTVRTopTags.at(i).numberOfDaughters(), weight);
        tau1_HOTVRTaggedjet1->Fill(HOTVRTopTags.at(i).tau1_groomed(), weight);
        tau2_HOTVRTaggedjet1->Fill(HOTVRTopTags.at(i).tau2_groomed(), weight);
        tau3_HOTVRTaggedjet1->Fill(HOTVRTopTags.at(i).tau3_groomed(), weight);
        tau21_HOTVRTaggedjet1->Fill(tau21, weight);
        tau32_HOTVRTaggedjet1->Fill(tau32, weight);
      }
      else if(i==1){
        pt_HOTVRTaggedjet2->Fill(HOTVRTopTags.at(i).pt(), weight);
        eta_HOTVRTaggedjet2->Fill(HOTVRTopTags.at(i).eta(), weight);
        phi_HOTVRTaggedjet2->Fill(HOTVRTopTags.at(i).phi(), weight);
        m_HOTVRTaggedjet2->Fill(HOTVRTopTags.at(i).v4().M(), weight);
        dRmin_AK8_HOTVRTaggedjet2->Fill(dRmin_ak8, weight);
        dRmin_mu_HOTVRTaggedjet2->Fill(dRmin_muon, weight);
        N_subjets_HOTVRTaggedjet2->Fill(HOTVRTopTags.at(i).subjets().size(), weight);
        N_daughters_HOTVRTaggedjet2->Fill(HOTVRTopTags.at(i).numberOfDaughters(), weight);
        tau1_HOTVRTaggedjet2->Fill(HOTVRTopTags.at(i).tau1_groomed(), weight);
        tau2_HOTVRTaggedjet2->Fill(HOTVRTopTags.at(i).tau2_groomed(), weight);
        tau3_HOTVRTaggedjet2->Fill(HOTVRTopTags.at(i).tau3_groomed(), weight);
        tau21_HOTVRTaggedjet2->Fill(tau21, weight);
        tau32_HOTVRTaggedjet2->Fill(tau32, weight);
      }
      else if(i==2){
        pt_HOTVRTaggedjet3->Fill(HOTVRTopTags.at(i).pt(), weight);
        eta_HOTVRTaggedjet3->Fill(HOTVRTopTags.at(i).eta(), weight);
        phi_HOTVRTaggedjet3->Fill(HOTVRTopTags.at(i).phi(), weight);
        m_HOTVRTaggedjet3->Fill(HOTVRTopTags.at(i).v4().M(), weight);
        dRmin_AK8_HOTVRTaggedjet3->Fill(dRmin_ak8, weight);
        dRmin_mu_HOTVRTaggedjet3->Fill(dRmin_muon, weight);
        N_subjets_HOTVRTaggedjet3->Fill(HOTVRTopTags.at(i).subjets().size(), weight);
        N_daughters_HOTVRTaggedjet3->Fill(HOTVRTopTags.at(i).numberOfDaughters(), weight);
        tau1_HOTVRTaggedjet3->Fill(HOTVRTopTags.at(i).tau1_groomed(), weight);
        tau2_HOTVRTaggedjet3->Fill(HOTVRTopTags.at(i).tau2_groomed(), weight);
        tau3_HOTVRTaggedjet3->Fill(HOTVRTopTags.at(i).tau3_groomed(), weight);
        tau21_HOTVRTaggedjet3->Fill(tau21, weight);
        tau32_HOTVRTaggedjet3->Fill(tau32, weight);
      }

    }

    N_HOTVRjets->Fill(NHOTVRjets, weight);
    N_HOTVRTaggedjets->Fill(NHOTVRTaggedjets, weight);

  }//end hotvr mode

  /*
  █  █████  ██   ██  █████  ██████  ██    ██ ██████  ██████  ██
  █ ██   ██ ██  ██  ██   ██ ██   ██ ██    ██ ██   ██ ██   ██ ██
  █ ███████ █████    █████  ██████  ██    ██ ██████  ██████  ██
  █ ██   ██ ██  ██  ██   ██ ██      ██    ██ ██      ██      ██
  █ ██   ██ ██   ██  █████  ██       ██████  ██      ██      ██
  */

  if(isdeepAK8){
    vector<TopJet>* AK8Puppijets = event.toppuppijets;
    unsigned int NAK8Puppijets = 0;
    for(unsigned int i=0; i<AK8Puppijets->size(); i++){
      if(AK8Puppijets->at(i).numberOfDaughters()<2) continue;
      NAK8Puppijets++;

      double tau21 = AK8Puppijets->at(i).tau2() / AK8Puppijets->at(i).tau1();
      double tau32 = AK8Puppijets->at(i).tau3() / AK8Puppijets->at(i).tau2();

      deepjet_topscore_jet->Fill(AK8Puppijets->at(i).btag_DeepBoosted_raw_score_top(), weight);
      deepjet_wscore_jet->Fill(AK8Puppijets->at(i).btag_DeepBoosted_raw_score_w(), weight);
      deepjet_zscore_jet->Fill(AK8Puppijets->at(i).btag_DeepBoosted_raw_score_z(), weight);
      deepjet_higgsscore_jet->Fill(AK8Puppijets->at(i).btag_DeepBoosted_raw_score_h(), weight);
      deepjet_qcdscore_jet->Fill(AK8Puppijets->at(i).btag_DeepBoosted_raw_score_qcd(), weight);
      massdecordeepjet_topscore_jet->Fill(AK8Puppijets->at(i).btag_MassDecorrelatedDeepBoosted_raw_score_top(), weight);
      massdecordeepjet_wscore_jet->Fill(AK8Puppijets->at(i).btag_MassDecorrelatedDeepBoosted_raw_score_w(), weight);
      massdecordeepjet_zscore_jet->Fill(AK8Puppijets->at(i).btag_MassDecorrelatedDeepBoosted_raw_score_z(), weight);
      massdecordeepjet_higgsscore_jet->Fill(AK8Puppijets->at(i).btag_MassDecorrelatedDeepBoosted_raw_score_h(), weight);
      massdecordeepjet_qcdscore_jet->Fill(AK8Puppijets->at(i).btag_MassDecorrelatedDeepBoosted_raw_score_qcd(), weight);
      deepjet_TvsQCD_jet->Fill(AK8Puppijets->at(i).btag_DeepBoosted_TvsQCD(), weight);
      massdecordeepjet_TvsQCD_jet->Fill(AK8Puppijets->at(i).btag_MassDecorrelatedDeepBoosted_TvsQCD(), weight);

      // Distance to AK8
      double dRmin_ak8 = 99999;
      for(unsigned int j=0; j<event.topjets->size(); j++){
        double dR = deltaR(AK8Puppijets->at(i), event.topjets->at(j));
        if(dR < dRmin_ak8) dRmin_ak8 = dR;
      }

      // Distance to muons
      double dRmin_muon = 99999;
      for(unsigned int j=0; j<event.muons->size(); j++){
        double dR = deltaR(AK8Puppijets->at(i), event.muons->at(j));
        if(dR < dRmin_muon) dRmin_muon = dR;
      }

      pt_AK8Puppijet->Fill(AK8Puppijets->at(i).pt(), weight);
      eta_AK8Puppijet->Fill(AK8Puppijets->at(i).eta(), weight);
      phi_AK8Puppijet->Fill(AK8Puppijets->at(i).phi(), weight);
      mSD_AK8Puppijet->Fill(AK8Puppijets->at(i).softdropmass(), weight);
      dRmin_AK8_AK8Puppijet->Fill(dRmin_ak8, weight);
      dRmin_mu_AK8Puppijet->Fill(dRmin_muon, weight);
      N_subjets_AK8Puppijet->Fill(AK8Puppijets->at(i).subjets().size(), weight);
      N_daughters_AK8Puppijet->Fill(AK8Puppijets->at(i).numberOfDaughters(), weight);
      tau1_AK8Puppijet->Fill(AK8Puppijets->at(i).tau1(), weight);
      tau2_AK8Puppijet->Fill(AK8Puppijets->at(i).tau2(), weight);
      tau3_AK8Puppijet->Fill(AK8Puppijets->at(i).tau3(), weight);
      tau21_AK8Puppijet->Fill(tau21, weight);
      tau32_AK8Puppijet->Fill(tau32, weight);

      if(i==0){
        pt_AK8Puppijet1->Fill(AK8Puppijets->at(i).pt(), weight);
        eta_AK8Puppijet1->Fill(AK8Puppijets->at(i).eta(), weight);
        phi_AK8Puppijet1->Fill(AK8Puppijets->at(i).phi(), weight);
        mSD_AK8Puppijet1->Fill(AK8Puppijets->at(i).softdropmass(), weight);
        dRmin_AK8_AK8Puppijet1->Fill(dRmin_ak8, weight);
        dRmin_mu_AK8Puppijet1->Fill(dRmin_muon, weight);
        N_subjets_AK8Puppijet1->Fill(AK8Puppijets->at(i).subjets().size(), weight);
        N_daughters_AK8Puppijet1->Fill(AK8Puppijets->at(i).numberOfDaughters(), weight);
        tau1_AK8Puppijet1->Fill(AK8Puppijets->at(i).tau1(), weight);
        tau2_AK8Puppijet1->Fill(AK8Puppijets->at(i).tau2(), weight);
        tau3_AK8Puppijet1->Fill(AK8Puppijets->at(i).tau3(), weight);
        tau21_AK8Puppijet1->Fill(tau21, weight);
        tau32_AK8Puppijet1->Fill(tau32, weight);
        deepjet_topscore_jet1->Fill(AK8Puppijets->at(i).btag_DeepBoosted_raw_score_top(), weight);
        deepjet_wscore_jet1->Fill(AK8Puppijets->at(i).btag_DeepBoosted_raw_score_w(), weight);
        deepjet_zscore_jet1->Fill(AK8Puppijets->at(i).btag_DeepBoosted_raw_score_z(), weight);
        deepjet_higgsscore_jet1->Fill(AK8Puppijets->at(i).btag_DeepBoosted_raw_score_h(), weight);
        deepjet_qcdscore_jet1->Fill(AK8Puppijets->at(i).btag_DeepBoosted_raw_score_qcd(), weight);
        massdecordeepjet_topscore_jet1->Fill(AK8Puppijets->at(i).btag_MassDecorrelatedDeepBoosted_raw_score_top(), weight);
        massdecordeepjet_wscore_jet1->Fill(AK8Puppijets->at(i).btag_MassDecorrelatedDeepBoosted_raw_score_w(), weight);
        massdecordeepjet_zscore_jet1->Fill(AK8Puppijets->at(i).btag_MassDecorrelatedDeepBoosted_raw_score_z(), weight);
        massdecordeepjet_higgsscore_jet1->Fill(AK8Puppijets->at(i).btag_MassDecorrelatedDeepBoosted_raw_score_h(), weight);
        massdecordeepjet_qcdscore_jet1->Fill(AK8Puppijets->at(i).btag_MassDecorrelatedDeepBoosted_raw_score_qcd(), weight);
        deepjet_TvsQCD_jet1->Fill(AK8Puppijets->at(i).btag_DeepBoosted_TvsQCD(), weight);
        massdecordeepjet_TvsQCD_jet1->Fill(AK8Puppijets->at(i).btag_MassDecorrelatedDeepBoosted_TvsQCD(), weight);
      }
      else if(i==1){
        pt_AK8Puppijet2->Fill(AK8Puppijets->at(i).pt(), weight);
        eta_AK8Puppijet2->Fill(AK8Puppijets->at(i).eta(), weight);
        phi_AK8Puppijet2->Fill(AK8Puppijets->at(i).phi(), weight);
        mSD_AK8Puppijet2->Fill(AK8Puppijets->at(i).softdropmass(), weight);
        dRmin_AK8_AK8Puppijet2->Fill(dRmin_ak8, weight);
        dRmin_mu_AK8Puppijet2->Fill(dRmin_muon, weight);
        N_subjets_AK8Puppijet2->Fill(AK8Puppijets->at(i).subjets().size(), weight);
        N_daughters_AK8Puppijet2->Fill(AK8Puppijets->at(i).numberOfDaughters(), weight);
        tau1_AK8Puppijet2->Fill(AK8Puppijets->at(i).tau1(), weight);
        tau2_AK8Puppijet2->Fill(AK8Puppijets->at(i).tau2(), weight);
        tau3_AK8Puppijet2->Fill(AK8Puppijets->at(i).tau3(), weight);
        tau21_AK8Puppijet2->Fill(tau21, weight);
        tau32_AK8Puppijet2->Fill(tau32, weight);
        deepjet_topscore_jet2->Fill(AK8Puppijets->at(i).btag_DeepBoosted_raw_score_top(), weight);
        deepjet_wscore_jet2->Fill(AK8Puppijets->at(i).btag_DeepBoosted_raw_score_w(), weight);
        deepjet_zscore_jet2->Fill(AK8Puppijets->at(i).btag_DeepBoosted_raw_score_z(), weight);
        deepjet_higgsscore_jet2->Fill(AK8Puppijets->at(i).btag_DeepBoosted_raw_score_h(), weight);
        deepjet_qcdscore_jet2->Fill(AK8Puppijets->at(i).btag_DeepBoosted_raw_score_qcd(), weight);
        massdecordeepjet_topscore_jet2->Fill(AK8Puppijets->at(i).btag_MassDecorrelatedDeepBoosted_raw_score_top(), weight);
        massdecordeepjet_wscore_jet2->Fill(AK8Puppijets->at(i).btag_MassDecorrelatedDeepBoosted_raw_score_w(), weight);
        massdecordeepjet_zscore_jet2->Fill(AK8Puppijets->at(i).btag_MassDecorrelatedDeepBoosted_raw_score_z(), weight);
        massdecordeepjet_higgsscore_jet2->Fill(AK8Puppijets->at(i).btag_MassDecorrelatedDeepBoosted_raw_score_h(), weight);
        massdecordeepjet_qcdscore_jet2->Fill(AK8Puppijets->at(i).btag_MassDecorrelatedDeepBoosted_raw_score_qcd(), weight);
        deepjet_TvsQCD_jet2->Fill(AK8Puppijets->at(i).btag_DeepBoosted_TvsQCD(), weight);
        massdecordeepjet_TvsQCD_jet2->Fill(AK8Puppijets->at(i).btag_MassDecorrelatedDeepBoosted_TvsQCD(), weight);
      }
      else if(i==2){
        pt_AK8Puppijet3->Fill(AK8Puppijets->at(i).pt(), weight);
        eta_AK8Puppijet3->Fill(AK8Puppijets->at(i).eta(), weight);
        phi_AK8Puppijet3->Fill(AK8Puppijets->at(i).phi(), weight);
        mSD_AK8Puppijet3->Fill(AK8Puppijets->at(i).softdropmass(), weight);
        dRmin_AK8_AK8Puppijet3->Fill(dRmin_ak8, weight);
        dRmin_mu_AK8Puppijet3->Fill(dRmin_muon, weight);
        N_subjets_AK8Puppijet3->Fill(AK8Puppijets->at(i).subjets().size(), weight);
        N_daughters_AK8Puppijet3->Fill(AK8Puppijets->at(i).numberOfDaughters(), weight);
        tau1_AK8Puppijet3->Fill(AK8Puppijets->at(i).tau1(), weight);
        tau2_AK8Puppijet3->Fill(AK8Puppijets->at(i).tau2(), weight);
        tau3_AK8Puppijet3->Fill(AK8Puppijets->at(i).tau3(), weight);
        tau21_AK8Puppijet3->Fill(tau21, weight);
        tau32_AK8Puppijet3->Fill(tau32, weight);
        massdecordeepjet_topscore_jet3->Fill(AK8Puppijets->at(i).btag_MassDecorrelatedDeepBoosted_raw_score_top(), weight);
        deepjet_topscore_jet3->Fill(AK8Puppijets->at(i).btag_DeepBoosted_raw_score_top(), weight);
        deepjet_TvsQCD_jet3->Fill(AK8Puppijets->at(i).btag_DeepBoosted_TvsQCD(), weight);
        massdecordeepjet_TvsQCD_jet3->Fill(AK8Puppijets->at(i).btag_MassDecorrelatedDeepBoosted_TvsQCD(), weight);
      }
    }

    // Top-tagged jets AK8
    vector<TopJet> AK8PuppiTopTags = event.get(h_AK8TopTags);
    int NAK8PuppiTaggedjets = 0;
    for(unsigned int i=0; i<AK8PuppiTopTags.size(); i++){
      NAK8PuppiTaggedjets++;

      double tau21 = AK8PuppiTopTags.at(i).tau2() / AK8PuppiTopTags.at(i).tau1();
      double tau32 = AK8PuppiTopTags.at(i).tau3() / AK8PuppiTopTags.at(i).tau2();

      // Distance to HOTVR
      double dRmin_HOTVR = 99999;
      for(unsigned int j=0; j<event.topjets->size(); j++){
        double dR = deltaR(AK8PuppiTopTags.at(i), event.topjets->at(j));
        if(dR < dRmin_HOTVR) dRmin_HOTVR = dR;
      }

      // Distance to muons
      double dRmin_muon = 99999;
      for(unsigned int j=0; j<event.muons->size(); j++){
        double dR = deltaR(AK8PuppiTopTags.at(i), event.muons->at(j));
        if(dR < dRmin_muon) dRmin_muon = dR;
      }

      pt_AK8PuppiTaggedjet->Fill(AK8PuppiTopTags.at(i).pt(), weight);
      eta_AK8PuppiTaggedjet->Fill(AK8PuppiTopTags.at(i).eta(), weight);
      phi_AK8PuppiTaggedjet->Fill(AK8PuppiTopTags.at(i).phi(), weight);
      mSD_AK8PuppiTaggedjet->Fill(AK8PuppiTopTags.at(i).softdropmass(), weight);
      dRmin_HOTVR_AK8PuppiTaggedjet->Fill(dRmin_HOTVR, weight);
      dRmin_mu_AK8PuppiTaggedjet->Fill(dRmin_muon, weight);
      N_subjets_AK8PuppiTaggedjet->Fill(AK8PuppiTopTags.at(i).subjets().size(), weight);
      N_daughters_AK8PuppiTaggedjet->Fill(AK8PuppiTopTags.at(i).numberOfDaughters(), weight);
      tau1_AK8PuppiTaggedjet->Fill(AK8PuppiTopTags.at(i).tau1(), weight);
      tau2_AK8PuppiTaggedjet->Fill(AK8PuppiTopTags.at(i).tau2(), weight);
      tau3_AK8PuppiTaggedjet->Fill(AK8PuppiTopTags.at(i).tau3(), weight);
      tau21_AK8PuppiTaggedjet->Fill(tau21, weight);
      tau32_AK8PuppiTaggedjet->Fill(tau32, weight);
      massdecordeepjet_topscore_Taggedjet->Fill(AK8PuppiTopTags.at(i).btag_MassDecorrelatedDeepBoosted_raw_score_top(), weight);
      deepjet_topscore_Taggedjet->Fill(AK8PuppiTopTags.at(i).btag_DeepBoosted_raw_score_top(), weight);
      deepjet_TvsQCD_Taggedjet->Fill(AK8PuppiTopTags.at(i).btag_DeepBoosted_TvsQCD(), weight);
      massdecordeepjet_TvsQCD_Taggedjet->Fill(AK8PuppiTopTags.at(i).btag_MassDecorrelatedDeepBoosted_TvsQCD(), weight);

      if(i==0){
        pt_AK8PuppiTaggedjet1->Fill(AK8PuppiTopTags.at(i).pt(), weight);
        eta_AK8PuppiTaggedjet1->Fill(AK8PuppiTopTags.at(i).eta(), weight);
        phi_AK8PuppiTaggedjet1->Fill(AK8PuppiTopTags.at(i).phi(), weight);
        mSD_AK8PuppiTaggedjet1->Fill(AK8PuppiTopTags.at(i).softdropmass(), weight);
        dRmin_HOTVR_AK8PuppiTaggedjet1->Fill(dRmin_HOTVR, weight);
        dRmin_mu_AK8PuppiTaggedjet1->Fill(dRmin_muon, weight);
        N_subjets_AK8PuppiTaggedjet1->Fill(AK8PuppiTopTags.at(i).subjets().size(), weight);
        N_daughters_AK8PuppiTaggedjet1->Fill(AK8PuppiTopTags.at(i).numberOfDaughters(), weight);
        tau1_AK8PuppiTaggedjet1->Fill(AK8PuppiTopTags.at(i).tau1(), weight);
        tau2_AK8PuppiTaggedjet1->Fill(AK8PuppiTopTags.at(i).tau2(), weight);
        tau3_AK8PuppiTaggedjet1->Fill(AK8PuppiTopTags.at(i).tau3(), weight);
        tau21_AK8PuppiTaggedjet1->Fill(tau21, weight);
        tau32_AK8PuppiTaggedjet1->Fill(tau32, weight);
        massdecordeepjet_topscore_Taggedjet1->Fill(AK8PuppiTopTags.at(i).btag_MassDecorrelatedDeepBoosted_raw_score_top(), weight);
        deepjet_topscore_Taggedjet1->Fill(AK8PuppiTopTags.at(i).btag_DeepBoosted_raw_score_top(), weight);
        deepjet_TvsQCD_Taggedjet1->Fill(AK8PuppiTopTags.at(i).btag_DeepBoosted_TvsQCD(), weight);
        massdecordeepjet_TvsQCD_Taggedjet1->Fill(AK8PuppiTopTags.at(i).btag_MassDecorrelatedDeepBoosted_TvsQCD(), weight);
      }
      else if(i==1){
        pt_AK8PuppiTaggedjet2->Fill(AK8PuppiTopTags.at(i).pt(), weight);
        eta_AK8PuppiTaggedjet2->Fill(AK8PuppiTopTags.at(i).eta(), weight);
        phi_AK8PuppiTaggedjet2->Fill(AK8PuppiTopTags.at(i).phi(), weight);
        mSD_AK8PuppiTaggedjet2->Fill(AK8PuppiTopTags.at(i).softdropmass(), weight);
        dRmin_HOTVR_AK8PuppiTaggedjet2->Fill(dRmin_HOTVR, weight);
        dRmin_mu_AK8PuppiTaggedjet2->Fill(dRmin_muon, weight);
        N_subjets_AK8PuppiTaggedjet2->Fill(AK8PuppiTopTags.at(i).subjets().size(), weight);
        N_daughters_AK8PuppiTaggedjet2->Fill(AK8PuppiTopTags.at(i).numberOfDaughters(), weight);
        tau1_AK8PuppiTaggedjet2->Fill(AK8PuppiTopTags.at(i).tau1(), weight);
        tau2_AK8PuppiTaggedjet2->Fill(AK8PuppiTopTags.at(i).tau2(), weight);
        tau3_AK8PuppiTaggedjet2->Fill(AK8PuppiTopTags.at(i).tau3(), weight);
        tau21_AK8PuppiTaggedjet2->Fill(tau21, weight);
        tau32_AK8PuppiTaggedjet2->Fill(tau32, weight);
        massdecordeepjet_topscore_Taggedjet2->Fill(AK8PuppiTopTags.at(i).btag_MassDecorrelatedDeepBoosted_raw_score_top(), weight);
        deepjet_topscore_Taggedjet2->Fill(AK8PuppiTopTags.at(i).btag_DeepBoosted_raw_score_top(), weight);
        deepjet_TvsQCD_Taggedjet2->Fill(AK8PuppiTopTags.at(i).btag_DeepBoosted_TvsQCD(), weight);
        massdecordeepjet_TvsQCD_Taggedjet2->Fill(AK8PuppiTopTags.at(i).btag_MassDecorrelatedDeepBoosted_TvsQCD(), weight);
      }
      else if(i==2){
        pt_AK8PuppiTaggedjet3->Fill(AK8PuppiTopTags.at(i).pt(), weight);
        eta_AK8PuppiTaggedjet3->Fill(AK8PuppiTopTags.at(i).eta(), weight);
        phi_AK8PuppiTaggedjet3->Fill(AK8PuppiTopTags.at(i).phi(), weight);
        mSD_AK8PuppiTaggedjet3->Fill(AK8PuppiTopTags.at(i).softdropmass(), weight);
        dRmin_HOTVR_AK8PuppiTaggedjet3->Fill(dRmin_HOTVR, weight);
        dRmin_mu_AK8PuppiTaggedjet3->Fill(dRmin_muon, weight);
        N_subjets_AK8PuppiTaggedjet3->Fill(AK8PuppiTopTags.at(i).subjets().size(), weight);
        N_daughters_AK8PuppiTaggedjet3->Fill(AK8PuppiTopTags.at(i).numberOfDaughters(), weight);
        tau1_AK8PuppiTaggedjet3->Fill(AK8PuppiTopTags.at(i).tau1(), weight);
        tau2_AK8PuppiTaggedjet3->Fill(AK8PuppiTopTags.at(i).tau2(), weight);
        tau3_AK8PuppiTaggedjet3->Fill(AK8PuppiTopTags.at(i).tau3(), weight);
        tau21_AK8PuppiTaggedjet3->Fill(tau21, weight);
        tau32_AK8PuppiTaggedjet3->Fill(tau32, weight);
        massdecordeepjet_topscore_Taggedjet3->Fill(AK8PuppiTopTags.at(i).btag_MassDecorrelatedDeepBoosted_raw_score_top(), weight);
        deepjet_topscore_Taggedjet3->Fill(AK8PuppiTopTags.at(i).btag_DeepBoosted_raw_score_top(), weight);
        deepjet_TvsQCD_Taggedjet2->Fill(AK8PuppiTopTags.at(i).btag_DeepBoosted_TvsQCD(), weight);
        massdecordeepjet_TvsQCD_Taggedjet3->Fill(AK8PuppiTopTags.at(i).btag_MassDecorrelatedDeepBoosted_TvsQCD(), weight);
      }

    }

    N_AK8Puppijets->Fill(NAK8Puppijets, weight);
    N_AK8PuppiTaggedjets->Fill(NAK8PuppiTaggedjets, weight);

  }//end deepak8 mode



  /*
  ███    ███ ██    ██  ██████  ███    ██ ███████
  ████  ████ ██    ██ ██    ██ ████   ██ ██
  ██ ████ ██ ██    ██ ██    ██ ██ ██  ██ ███████
  ██  ██  ██ ██    ██ ██    ██ ██  ██ ██      ██
  ██      ██  ██████   ██████  ██   ████ ███████
  */


  vector<Muon>* muons = event.muons;
  int Nmuons = muons->size();
  N_mu->Fill(Nmuons, weight);
  // cout << "N_mu: " << Nmuons << endl;

  for(int i=0; i<Nmuons; i++){

    N_mu_charge->Fill(muons->at(i).charge(), weight);
    

    // if(muons->at(i).charge() > 0){
    //   cout << "muon charge" << muons->at(i).charge() << endl;
    // }
    // cout << "muon charge" << muons->at(i).charge() << endl;

    pt_mu->Fill(muons->at(i).pt(),weight);
    eta_mu->Fill(muons->at(i).eta(),weight);
    phi_mu->Fill(muons->at(i).phi(),weight);
    reliso_mu->Fill(muons->at(i).relIso(),weight);
    reliso_mu_rebin->Fill(muons->at(i).relIso(),weight);
    dRmin_mu_jet->Fill(muons->at(i).get_tag(Muon::twodcut_dRmin), weight);
    dRmin_mu_jet_scaled->Fill(muons->at(i).get_tag(Muon::twodcut_dRmin)*event.jets->at(0).pt(), weight);
    ptrel_mu_jet->Fill(muons->at(i).get_tag(Muon::twodcut_pTrel), weight);
    dRmin_ptrel_mu->Fill(muons->at(i).get_tag(Muon::twodcut_dRmin), muons->at(i).get_tag(Muon::twodcut_pTrel), weight);
    if(i==0){
      pt_mu1->Fill(muons->at(i).pt(),weight);
      eta_mu1->Fill(muons->at(i).eta(),weight);
      phi_mu1->Fill(muons->at(i).phi(),weight);
      reliso_mu1->Fill(muons->at(i).relIso(),weight);
      reliso_mu1_rebin->Fill(muons->at(i).relIso(),weight);
      dRmin_mu1_jet->Fill(muons->at(i).get_tag(Muon::twodcut_dRmin), weight);
      dRmin_mu1_jet_scaled->Fill(muons->at(i).get_tag(Muon::twodcut_dRmin)*event.jets->at(0).pt(), weight);
      ptrel_mu1_jet->Fill(muons->at(i).get_tag(Muon::twodcut_pTrel), weight);
      dRmin_ptrel_mu1->Fill(muons->at(i).get_tag(Muon::twodcut_dRmin), muons->at(i).get_tag(Muon::twodcut_pTrel), weight);
    }
    else if(i==1){
      pt_mu2->Fill(muons->at(i).pt(),weight);
      eta_mu2->Fill(muons->at(i).eta(),weight);
      phi_mu2->Fill(muons->at(i).phi(),weight);
      reliso_mu2->Fill(muons->at(i).relIso(),weight);
      reliso_mu2_rebin->Fill(muons->at(i).relIso(),weight);
    }


    double minDeltaR = 999;
    for(unsigned int k = 0; k < jets->size(); k++){
      if(deltaR(muons->at(i), jets->at(k))<minDeltaR){
        minDeltaR = deltaR(muons->at(i), jets->at(k));
        dR_mu_nearjet->Fill(minDeltaR,weight);
        pt_mu_nearjet->Fill(jets->at(k).pt(),weight);
      }
    }
  }

  for(int i=0; i<Nmuons; i++){
    for(int j=0; j<Nmuons; j++){
      if(j <= i) continue;
      M_mumu->Fill((muons->at(i).v4() + muons->at(j).v4()).M() ,weight);
    }
  }

  /*
  ███████ ██      ███████  ██████ ████████ ██████   ██████  ███    ██ ███████
  ██      ██      ██      ██         ██    ██   ██ ██    ██ ████   ██ ██
  █████   ██      █████   ██         ██    ██████  ██    ██ ██ ██  ██ ███████
  ██      ██      ██      ██         ██    ██   ██ ██    ██ ██  ██ ██      ██
  ███████ ███████ ███████  ██████    ██    ██   ██  ██████  ██   ████ ███████
  */


  vector<Electron>* electrons = event.electrons;
  int Nelectrons = electrons->size();
  N_ele->Fill(Nelectrons, weight);

  for(int i=0; i<Nelectrons; i++){
    N_ele_charge->Fill(electrons->at(i).charge(), weight);
    pt_ele->Fill(electrons->at(i).pt(),weight);
    eta_ele->Fill(electrons->at(i).eta(),weight);
    phi_ele->Fill(electrons->at(i).phi(),weight);
    reliso_ele->Fill(electrons->at(i).relIso(),weight);
    reliso_ele_rebin->Fill(electrons->at(i).relIso(),weight);
    dRmin_ele_jet->Fill(electrons->at(i).get_tag(Electron::twodcut_dRmin), weight);
    dRmin_ele_jet_scaled->Fill(electrons->at(i).get_tag(Electron::twodcut_dRmin)*event.jets->at(0).pt(), weight);
    ptrel_ele_jet->Fill(electrons->at(i).get_tag(Electron::twodcut_pTrel), weight);
    dRmin_ptrel_ele->Fill(electrons->at(i).get_tag(Electron::twodcut_dRmin), electrons->at(i).get_tag(Electron::twodcut_pTrel), weight);
    if(electrons->at(i).pt()<120){
      pt_ele_lowpt->Fill(electrons->at(i).pt(),weight);
      eta_ele_lowpt->Fill(electrons->at(i).eta(),weight);
      phi_ele_lowpt->Fill(electrons->at(i).phi(),weight);
    }
    if(electrons->at(i).pt()>=120 && electrons->at(i).pt()<200){
      pt_ele_midpt->Fill(electrons->at(i).pt(),weight);
      eta_ele_midpt->Fill(electrons->at(i).eta(),weight);
      phi_ele_midpt->Fill(electrons->at(i).phi(),weight);
    }
    if(electrons->at(i).pt()>=200){
      pt_ele_highpt->Fill(electrons->at(i).pt(),weight);
      eta_ele_highpt->Fill(electrons->at(i).eta(),weight);
      phi_ele_highpt->Fill(electrons->at(i).phi(),weight);
    }
    if(i==0){
      pt_ele1->Fill(electrons->at(i).pt(),weight);
      eta_ele1->Fill(electrons->at(i).eta(),weight);
      phi_ele1->Fill(electrons->at(i).phi(),weight);
      reliso_ele1->Fill(electrons->at(i).relIso(),weight);
      reliso_ele1_rebin->Fill(electrons->at(i).relIso(),weight);
      dRmin_ele1_jet->Fill(electrons->at(i).get_tag(Electron::twodcut_dRmin), weight);
      dRmin_ele1_jet_scaled->Fill(electrons->at(i).get_tag(Electron::twodcut_dRmin)*event.jets->at(0).pt(), weight);
      ptrel_ele1_jet->Fill(electrons->at(i).get_tag(Electron::twodcut_pTrel), weight);
      dRmin_ptrel_ele1->Fill(electrons->at(i).get_tag(Electron::twodcut_dRmin), electrons->at(i).get_tag(Electron::twodcut_pTrel), weight);
    }
    else if(i==1){
      pt_ele2->Fill(electrons->at(i).pt(),weight);
      eta_ele2->Fill(electrons->at(i).eta(),weight);
      phi_ele2->Fill(electrons->at(i).phi(),weight);
      reliso_ele2->Fill(electrons->at(i).relIso(),weight);
      reliso_ele2_rebin->Fill(electrons->at(i).relIso(),weight);
    }


    for(unsigned int k = 0; k < jets->size(); k++){
      double minDeltaR = 999;
      if(deltaR(electrons->at(i), jets->at(k))<minDeltaR){
        minDeltaR = deltaR(electrons->at(i), jets->at(k));
        dR_ele_nearjet->Fill(minDeltaR,weight);
        pt_ele_nearjet->Fill(jets->at(k).pt(),weight);
      }
    }

  }

  for(int i=0; i<Nelectrons; i++){
    for(int j=0; j<Nelectrons; j++){
      if(j <= i) continue;
      M_ee->Fill((electrons->at(i).v4() + electrons->at(j).v4()).M() ,weight);
    }
  }

  /*
  ██████  ███████ ███    ██ ███████ ██████   █████  ██
  ██      ██      ████   ██ ██      ██   ██ ██   ██ ██
  ██  ███ █████   ██ ██  ██ █████   ██████  ███████ ██
  ██   ██ ██      ██  ██ ██ ██      ██   ██ ██   ██ ██
  ██████  ███████ ██   ████ ███████ ██   ██ ██   ██ ███████
  */



  int Npvs = event.pvs->size();
  NPV->Fill(Npvs, weight);

  double met = event.met->pt();
  double st = 0., st_jets = 0., st_lep = 0.;
  double ht_lep = 0.;
  for(const auto & jet : *jets)           st_jets += jet.pt();
  for(const auto & electron : *electrons) st_lep += electron.pt();
  for(const auto & muon : *muons)         st_lep += muon.pt();
  st = st_jets + st_lep + met;
  ht_lep = st_lep + met;

  MET->Fill(met, weight);
  MET_rebin->Fill(met, weight);
  MET_rebin2->Fill(met, weight);
  MET_rebin3->Fill(met, weight);
  ST->Fill(st, weight);
  ST_rebin->Fill(st, weight);
  ST_rebin2->Fill(st, weight);
  ST_rebin3->Fill(st, weight);
  STjets->Fill(st_jets, weight);
  STjets_rebin->Fill(st_jets, weight);
  STjets_rebin2->Fill(st_jets, weight);
  STjets_rebin3->Fill(st_jets, weight);
  STlep->Fill(ht_lep, weight);
  STlep_rebin->Fill(ht_lep, weight);
  STlep_rebin2->Fill(ht_lep, weight);
  STlep_rebin3->Fill(ht_lep, weight);

  // Zprime reco
  bool is_zprime_reconstructed_chi2 = event.get(h_is_zprime_reconstructed_chi2);
  bool is_zprime_reconstructed_correctmatch = event.get(h_is_zprime_reconstructed_correctmatch);
  // added "is_mc" to blind data in mttbar hists
  
  if(debug) cout << "Before dY lines:" << endl;


  // ================== DY new =========================================================
  
  if(is_zprime_reconstructed_chi2 && is_tt){
    const auto& genparticles = event.genparticles;
    ZprimeCandidate* BestZprimeCandidate = event.get(h_BestZprimeCandidateChi2);
    float Mreco = BestZprimeCandidate->Zprime_v4().M();
    // float chi2 = BestZprimeCandidate->discriminator("chi2_total");
    ditop_mass->Fill(Mreco, weight);
    M_Zprime->Fill(Mreco, weight);
    M_Zprime_rebin->Fill(Mreco, weight);
    M_Zprime_rebin2->Fill(Mreco, weight);
    M_Zprime_rebin3->Fill(Mreco, weight);

    GenParticle top, antitop;
    for(const GenParticle & gp : *genparticles){
      if(gp.pdgId() == 6){
          top = gp;
      }
      else if(gp.pdgId() == -6){
          antitop = gp;
      }
    }
    if(debug) cout << "after gen 1:" << endl;
    // The Lorentz vectors represent the 4-momenta (energy, and three spatial momentum components) for the leptonic and hadronic tops from the "BestZprimeCandidate" object
    LorentzVector lep_top = BestZprimeCandidate->top_leptonic_v4();
    LorentzVector had_top = BestZprimeCandidate->top_hadronic_v4();

    // vectors to store the deltaR values for the leptonic and hadronic tops with each gen particle
    // this part initializes vectors to store deltaR values with a default of 99.0 and fills in the actual deltaR values by looping over the gen particles (top)
    std::vector<std::pair<double, int>> deltaR_leptonic_values; // ((dR, index), (dR, index), ...)
    std::vector<std::pair<double, int>> deltaR_hadronic_values;

    double deltaR_min_leptonic = 99.0;
    double deltaR_sec_min_leptonic = 99.0;
    int best_gen_for_leptop = -1;
    int sec_best_gen_for_leptop = -1;
    // bool is_leptop_matched = false;

    double deltaR_min_hadronic = 99.0;
    double deltaR_sec_min_hadronic = 99.0;
    int best_gen_for_hadtop = -1;
    int sec_best_gen_for_hadtop = -1;
    // bool is_hadtop_matched = false;

    if(debug) cout << "before gen loop in matching" << endl;
    
    for(unsigned int j=0; j<genparticles->size(); ++j) {
      if(abs(genparticles->at(j).pdgId()) == 6 ){
        if (genparticles->at(j).index() == 2 || genparticles->at(j).index() == 3){
          LorentzVector genparticle_p4(genparticles->at(j).pt(), genparticles->at(j).eta(), genparticles->at(j).phi(), genparticles->at(j).energy());
          deltaR_leptonic_values.push_back(std::make_pair(deltaR(lep_top, genparticle_p4),genparticles->at(j).index() ));
          deltaR_hadronic_values.push_back(std::make_pair(deltaR(had_top, genparticle_p4), genparticles->at(j).index()));
          // cout << "deltaR: " << deltaR(lep_top, genparticle_p4) << j << endl;
        }
      }
    }



    for (const auto& pair_lep : deltaR_leptonic_values) {
      if (pair_lep.first > 0 && pair_lep.first < deltaR_min_leptonic) {
        // deltaR_min_leptonic = pair_lep.first;
        deltaR_sec_min_leptonic = deltaR_min_leptonic;
        deltaR_min_leptonic = pair_lep.first;
        sec_best_gen_for_leptop = best_gen_for_leptop;
        best_gen_for_leptop = pair_lep.second;
        // is_leptop_matched = true;
      }
      else if (pair_lep.first > 0 && pair_lep.first < deltaR_sec_min_leptonic && pair_lep.first != deltaR_min_leptonic && pair_lep.second != best_gen_for_leptop) {
        deltaR_sec_min_leptonic = pair_lep.first;
        sec_best_gen_for_leptop = pair_lep.second;
      }
    }


    for (const auto& pair_had : deltaR_hadronic_values) {
      if (pair_had.first > 0 && pair_had.first < deltaR_min_hadronic) {
        // deltaR_min_hadronic = pair_had.first;
        deltaR_sec_min_hadronic = deltaR_min_hadronic;
        deltaR_min_hadronic = pair_had.first;
        sec_best_gen_for_hadtop = best_gen_for_hadtop;
        best_gen_for_hadtop = pair_had.second;
        // is_hadtop_matched = true;
      }
      else if (pair_had.first > 0 && pair_had.first < deltaR_sec_min_hadronic && pair_had.first != deltaR_min_hadronic && pair_had.second != best_gen_for_hadtop) {
        deltaR_sec_min_hadronic = pair_had.first;
        sec_best_gen_for_hadtop = pair_had.second;
      }
    }

    if(debug) cout << "after dR matching" << endl;



    if (best_gen_for_hadtop == best_gen_for_leptop){
      if(debug) cout << "same index for matched gen" << endl;
      if (deltaR_min_leptonic <= deltaR_min_hadronic){
        if(debug) cout << "deltaR_min_leptonic <= deltaR_min_hadronic" << endl;
        best_gen_for_hadtop = sec_best_gen_for_hadtop;
        deltaR_min_hadronic = deltaR_sec_min_hadronic;
      } else {
        if(debug) cout << "deltaR_min_hadronic <= deltaR_min_leptonic" << endl;
        best_gen_for_leptop = sec_best_gen_for_leptop;
        deltaR_min_leptonic = deltaR_sec_min_leptonic;
      }
      if(debug) cout << "same index for matched gen - done" << endl;
    }
    if(debug) cout << "after the same index statement" << endl;


    GenParticle best_matched_gen_leptop;
    GenParticle best_matched_gen_hadtop;
    
    if(debug) cout <<  genparticles->size() << endl;
    if(debug) cout << static_cast<std::size_t>(best_gen_for_leptop) << endl;
    if(debug) cout << "leptop index: " << best_gen_for_leptop << endl;
    if(debug) cout << "hadtop index :" << best_gen_for_hadtop << endl;
    
    

    float_t DeltaY_gen_best = 99.0;
    float_t DeltaY_reco_best = 99.0;
   

    if(debug) cout << "now will check dR " << endl;
    if(debug) cout << deltaR_min_leptonic << endl;
    if(debug) cout << deltaR_min_hadronic << endl;
    if(debug) cout << best_gen_for_leptop << endl;
    if(debug) cout << best_gen_for_hadtop << endl;

    if (deltaR_min_leptonic < 0.4 && deltaR_min_hadronic < 0.4 && best_gen_for_leptop >= 0 && best_gen_for_hadtop >= 0) {
      if(debug) cout << "in dR loop " << endl;
      if(static_cast<std::size_t>(best_gen_for_leptop) < genparticles->size()) {
        if(debug) cout << "looking for lep match " << endl;
        best_matched_gen_leptop = genparticles->at(best_gen_for_leptop);
      }
      if(debug) cout << "lepton match done" << endl;
      
      if(static_cast<std::size_t>(best_gen_for_hadtop) < genparticles->size()) {
          best_matched_gen_hadtop = genparticles->at(best_gen_for_hadtop);
          if(debug) cout << "hadtop match done" << endl;
      }

      if(debug) cout << "after genparticles->at(best_gen_for_leptop)" << endl;

      N_lep_charge->Fill(BestZprimeCandidate->lepton().charge(),weight);
      // cout <<"Lepton charge is: "<< BestZprimeCandidate->lepton().charge()<<endl;
      // Calculates the delta y (with reco particles) values for the leptonic and hadronic tops depending on the charge of the lepton
      if (BestZprimeCandidate->lepton().charge()>0) {
        DeltaY_reco_best = TMath::Abs(0.5*TMath::Log((lep_top.energy() + lep_top.pt()*TMath::SinH(lep_top.eta()))/(lep_top.energy() - lep_top.pt()*TMath::SinH(lep_top.eta())))) - TMath::Abs(0.5*TMath::Log((had_top.energy() + had_top.pt()*TMath::SinH(had_top.eta()))/(had_top.energy() - had_top.pt()*TMath::SinH(had_top.eta()))));
        DeltaY_gen_best = TMath::Abs(0.5*TMath::Log((best_matched_gen_leptop.energy() + best_matched_gen_leptop.pt()*TMath::SinH(best_matched_gen_leptop.eta()))/(best_matched_gen_leptop.energy() - best_matched_gen_leptop.pt()*TMath::SinH(best_matched_gen_leptop.eta())))) - TMath::Abs(0.5*TMath::Log((best_matched_gen_hadtop.energy() + best_matched_gen_hadtop.pt()*TMath::SinH(best_matched_gen_hadtop.eta()))/(best_matched_gen_hadtop.energy() - best_matched_gen_hadtop.pt()*TMath::SinH(best_matched_gen_hadtop.eta()))));

      } else {
        DeltaY_reco_best = TMath::Abs(0.5*TMath::Log((had_top.energy() + had_top.pt()*TMath::SinH(had_top.eta()))/(had_top.energy() - had_top.pt()*TMath::SinH(had_top.eta())))) - TMath::Abs(0.5*TMath::Log((lep_top.energy() + lep_top.pt()*TMath::SinH(lep_top.eta()))/(lep_top.energy() - lep_top.pt()*TMath::SinH(lep_top.eta()))));
        DeltaY_gen_best = TMath::Abs(0.5*TMath::Log((best_matched_gen_hadtop.energy() + best_matched_gen_hadtop.pt()*TMath::SinH(best_matched_gen_hadtop.eta()))/(best_matched_gen_hadtop.energy() - best_matched_gen_hadtop.pt()*TMath::SinH(best_matched_gen_hadtop.eta())))) - TMath::Abs(0.5*TMath::Log((best_matched_gen_leptop.energy() + best_matched_gen_leptop.pt()*TMath::SinH(best_matched_gen_leptop.eta()))/(best_matched_gen_leptop.energy() - best_matched_gen_leptop.pt()*TMath::SinH(best_matched_gen_leptop.eta()))));

      }
    }
    else {
      if(debug) cout << "not matched" << endl;
      DeltaY_notMatched->Fill(1);
    }

    if(debug) cout << "about to fill response matrix" << endl;

    response_matrix->Fill(DeltaY_reco_best, DeltaY_gen_best, weight);
    DeltaY_reco_best_plot->Fill(DeltaY_reco_best, weight);
    DeltaY_gen_best_plot->Fill(DeltaY_gen_best, weight);
    DeltaY_reco->Fill(DeltaY_reco_best, weight);
    if(debug) cout << "after filling dY hists" << endl;
  
  
  
  
  //begin spin correlation
  // ZprimeCandidate* BestZprimeCandidate = event.get(h_BestZprimeCandidateChi2);
    if (debug) cout << "starting spin corr" << endl;
    bool is_toptag_reconstruction = BestZprimeCandidate->is_toptag_reconstruction(); // Reconstruction process id
    vector <Jet> AK4CHSjets_matched = event.get(h_CHSjets_matched);                  // AK4Puppijets that have been matched to CHSjets
    vector <TopJet> TopTaggedJets = event.get(h_AK8TopTags);                     // AK8Puppi jets TopTagged by DeepAK8TopTagger
    vector <float> jets_hadronic_bscores;                                            // bScores vector for resolved hadronic jets
    float pt_hadTop_thresh = 150;                                                    // Define cut-variable as pt of hadTop for low/high regions
    float btag_WP=0.2783;                                                                   // see https://btv-wiki.docs.cern.ch/ScaleFactors/ for btag WPs
    if (isUL16preVFP) btag_WP = 0.2598;                                              // medium WP for UL16preVFP DeepJet
    if (isUL16postVFP) btag_WP = 0.3657;                                             // medium WP for UL16postVFP DeepJet
    if (isUL17) btag_WP = 0.3040;                                                    // medium WP for UL17 DeepJet
    if (isUL18) btag_WP = 0.2783;                                                    // medium WP for UL18 DeepJet

    // Plot pt of hadronic Top jet
    float pt_hadTop = BestZprimeCandidate->top_hadronic_v4().pt();
    
    // Variable that will hold the highest bscore of the hadronic top's (sub)jets
    float bscore_max = -2;
    if(!is_toptag_reconstruction){
        // Loop over resolved hadronic jets to find their bscore via CHS jets
      for(unsigned int i=0; i<BestZprimeCandidate->jets_hadronic().size(); i++){
        double deltaR_min = 99;
        // Match resolved hadronic jets to CHS jets (which have bscores)
        for(unsigned int j=0; j<AK4CHSjets_matched.size(); j++){
          double deltaR_CHS = deltaR(BestZprimeCandidate->jets_hadronic().at(i), AK4CHSjets_matched.at(j));
          if(deltaR_CHS < deltaR_min) deltaR_min = deltaR_CHS;}
        // Build bScore-vector for resolved hadronic jets whose bscore will correspond by index
        for(unsigned int k=0; k<AK4CHSjets_matched.size(); k++){
          if(deltaR(BestZprimeCandidate->jets_hadronic().at(i), AK4CHSjets_matched.at(k)) == deltaR_min) 
          jets_hadronic_bscores.emplace_back(AK4CHSjets_matched.at(k).btag_DeepJet());} // Using DeepJet btag score
      }
      // Loop over bScores-vector to extract highest bscore
      for(unsigned int i=0; i<jets_hadronic_bscores.size(); i++){
        float bscore = jets_hadronic_bscores.at(i);
        if(bscore > bscore_max) bscore_max = bscore;
      }
    
    }//is not top tag
    if (debug) cout << "done with not top tag" << endl;
    if(is_toptag_reconstruction){
        // Loop over hadronic top's subjets to extract highest bscore
      for(unsigned int i=0; i < BestZprimeCandidate->tophad_topjet_ptr()->subjets().size(); i++){
        float bscore = BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(i).btag_DeepJet(); // Using DeepJet btag score
        if(bscore > bscore_max) bscore_max = bscore;
      }

    }
    if (debug) cout << "done with top tag" << endl;
    if(bscore_max >= btag_WP){
        // event.set(h_bscore_max, bscore_max); // Plot max bscores

        //------------------------------------Define 4vectors of hadronic b-jet and lepton------------------------------------//
        // Hadronic b-jet 4-vector
        TLorentzVector had_top_b(0, 0, 0, 0);

        // Resolved topology
        if(!is_toptag_reconstruction){ // Define hadronic b-jet as hadronic AK4-jet with highest bscore
          for(unsigned int i=0; i< BestZprimeCandidate->jets_hadronic().size(); i++){
            float bscore = jets_hadronic_bscores.at(i);
            if(bscore == bscore_max) had_top_b.SetPtEtaPhiE(BestZprimeCandidate->jets_hadronic().at(i).pt(), 
                                                            BestZprimeCandidate->jets_hadronic().at(i).eta(), 
                                                            BestZprimeCandidate->jets_hadronic().at(i).phi(), 
                                                            BestZprimeCandidate->jets_hadronic().at(i).energy());
          }
        }
        // Merged topology
        if(is_toptag_reconstruction){ // Define hadronic b-jet as hadronic AK8-subjet with highest bscore
          for(unsigned int j=0; j < BestZprimeCandidate->tophad_topjet_ptr()->subjets().size(); j++){
            float bscore = BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).btag_DeepJet();
            if(bscore == bscore_max) had_top_b.SetPtEtaPhiE(BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).pt(), 
                                                            BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).eta(), 
                                                            BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).phi(), 
                                                            BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).energy());
          }
        }

        // Lepton 4-vector
        TLorentzVector lep_top_lep(0, 0, 0, 0);
        LorentzVector lep = BestZprimeCandidate->lepton().v4();
        lep_top_lep.SetPtEtaPhiE(lep.pt(), lep.eta(), lep.phi(), lep.E());
        //------------------------------------Define 4vectors of hadronic b-jet and lepton------------------------------------//


        //-------------------------------- Begin boosting top quarks and their decay products --------------------------------//
        // Define 4vectors of top quarks
        TLorentzVector PosTop(0, 0, 0, 0);
        TLorentzVector NegTop(0, 0, 0, 0);

        // POSITIVE LEPTON CONFIGURATION => Positive charged lepton has Positive Top mother
        if(BestZprimeCandidate->lepton().charge() > 0){
          // Define ttbar system
          PosTop.SetPtEtaPhiE(BestZprimeCandidate->top_leptonic_v4().pt(), 
                              BestZprimeCandidate->top_leptonic_v4().eta(), 
                              BestZprimeCandidate->top_leptonic_v4().phi(), 
                              BestZprimeCandidate->top_leptonic_v4().energy());
          NegTop.SetPtEtaPhiE(BestZprimeCandidate->top_hadronic_v4().pt(), 
                              BestZprimeCandidate->top_hadronic_v4().eta(), 
                              BestZprimeCandidate->top_hadronic_v4().phi(), 
                              BestZprimeCandidate->top_hadronic_v4().energy());
        }
        else if (BestZprimeCandidate->lepton().charge() < 0)
        {
        PosTop.SetPtEtaPhiE(BestZprimeCandidate->top_hadronic_v4().pt(), 
                              BestZprimeCandidate->top_hadronic_v4().eta(), 
                              BestZprimeCandidate->top_hadronic_v4().phi(), 
                              BestZprimeCandidate->top_hadronic_v4().energy());
          NegTop.SetPtEtaPhiE(BestZprimeCandidate->top_leptonic_v4().pt(), 
                              BestZprimeCandidate->top_leptonic_v4().eta(), 
                              BestZprimeCandidate->top_leptonic_v4().phi(), 
                              BestZprimeCandidate->top_leptonic_v4().energy());
        }
        
          TLorentzVector ttbar = PosTop + NegTop;
          // Boost into ttbar CoM-Frame <<<-------//
          lep_top_lep.Boost(-ttbar.BoostVector());
          had_top_b.Boost(-ttbar.BoostVector());
          PosTop.Boost(-ttbar.BoostVector());
          NegTop.Boost(-ttbar.BoostVector());

          // Rotate vectors into Helicity Frame <<<------//
          // Rotate about beamline
          lep_top_lep.RotateZ(-1.*PosTop.Phi());
          had_top_b.RotateZ(-1.*PosTop.Phi());
          PosTop.RotateZ(-1.*PosTop.Phi());
          NegTop.RotateZ(-1.*PosTop.Phi());
          // Rotate about y-axis
          lep_top_lep.RotateY(-1.*PosTop.Theta());
          had_top_b.RotateY(-1.*PosTop.Theta());
          PosTop.RotateY(-1.*PosTop.Theta());
          NegTop.RotateY(-1.*PosTop.Theta());

          // Boost into ttbar Rest-Frame <<<--------//
          lep_top_lep.Boost(-PosTop.BoostVector()); // Positive charged lepton has Positive Top mother
          had_top_b.Boost(-NegTop.BoostVector());   // Positive charged lepton means b-jet has Negative Top mother

        //-------------------------------- End boosting top quarks and their decay products --------------------------------//
    
        // Define angular variables as sum and difference of decay products' phi-coordinates
        // sphi and dphi = PosTopDecayProd_phi +- NegTopDecayProd_phi
        float dphi=0;
        float sphi = lep_top_lep.Phi() + had_top_b.Phi();
        if(BestZprimeCandidate->lepton().charge() > 0){ // lepton is Positive Top's Decay Product
          dphi = lep_top_lep.Phi() - had_top_b.Phi();
        }
        if(BestZprimeCandidate->lepton().charge() < 0)
        {
          dphi = had_top_b.Phi() - lep_top_lep.Phi();
        }
    
        if (debug) cout << "dphi is : " << dphi <<endl;
        if (debug) cout << "sphi is : " << sphi <<endl;
      
          // Map back into original domain if necessary
        if(sphi > TMath::Pi()) sphi = sphi - 2*TMath::Pi();
        if(sphi < -TMath::Pi()) sphi = sphi + 2*TMath::Pi();
        if(dphi > TMath::Pi()) dphi = dphi - 2*TMath::Pi();
        if(dphi < -TMath::Pi()) dphi = dphi + 2*TMath::Pi();
        Sigma_phi->Fill(sphi,weight);
        Delta_phi->Fill(dphi,weight);

          // Plot dphi and sphi for high-pt ranges
        if(pt_hadTop > pt_hadTop_thresh && DeltaY_reco_best >0){
          Sigma_phi_1->Fill(sphi,weight);
          Delta_phi_1->Fill(dphi,weight);
        }
        if(pt_hadTop > pt_hadTop_thresh && DeltaY_reco_best <0){
          Sigma_phi_2->Fill(sphi,weight);
          Delta_phi_2->Fill(dphi,weight);
        }
        if(pt_hadTop < pt_hadTop_thresh && dphi >0){
          DeltaY_reco_d1->Fill(DeltaY_reco_best,weight);
        }
        if(pt_hadTop < pt_hadTop_thresh && dphi <0){
          DeltaY_reco_d2->Fill(DeltaY_reco_best,weight);
        }
        if(pt_hadTop < pt_hadTop_thresh && sphi >0){
          DeltaY_reco_s1->Fill(DeltaY_reco_best,weight);
        }
        if(pt_hadTop < pt_hadTop_thresh && sphi <0){
          DeltaY_reco_s2->Fill(DeltaY_reco_best,weight);
        }
        if(pt_hadTop > pt_hadTop_thresh){
          Sigma_phi_high->Fill(sphi,weight);
          Delta_phi_high->Fill(dphi,weight);
          DeltaY_reco_high->Fill(DeltaY_reco_best,weight);
        }
        // Plot dphi and sphi for low-pt ranges
        if(pt_hadTop < pt_hadTop_thresh){
          Sigma_phi_low->Fill(sphi,weight);
          Delta_phi_low->Fill(dphi,weight);
          DeltaY_reco_low->Fill(DeltaY_reco_best,weight);
        }
      }





  //end spin correlation
  }

  // for backgrounds and DATA
  if(is_zprime_reconstructed_chi2 && !is_tt){
    if (debug) cout << "shouldnt be here if ttbar : " << endl;
   // const auto& genparticles = event.genparticles;
    ZprimeCandidate* BestZprimeCandidate = event.get(h_BestZprimeCandidateChi2);
    float Mreco = BestZprimeCandidate->Zprime_v4().M();
    float chi2 = BestZprimeCandidate->discriminator("chi2_total");
    ditop_mass->Fill(Mreco, weight);
    M_Zprime->Fill(Mreco, weight);
    M_Zprime_rebin->Fill(Mreco, weight);
    M_Zprime_rebin2->Fill(Mreco, weight);
    M_Zprime_rebin3->Fill(Mreco, weight);
    if(Mreco < 6000.) M_Zprime_rebin4->Fill(Mreco, weight);
    else M_Zprime_rebin4->Fill(6050., weight);
    if(Mreco < 5000.) M_Zprime_rebin5->Fill(Mreco, weight);
    else M_Zprime_rebin5->Fill(5050., weight);
    M_Zprime_rebin6->Fill(Mreco, weight);
    M_Zprime_rebin7->Fill(Mreco, weight);
    M_Zprime_rebin8->Fill(Mreco, weight);
    chi2_Zprime->Fill(chi2, weight);
    chi2_Zprime_rebin->Fill(chi2, weight);
    chi2_Zprime_rebin2->Fill(chi2, weight);
    // cout << "the boolean is: "<< isLeptonPositive << endl;
    float_t dyreco = 0.0;
    if (BestZprimeCandidate->lepton().charge()>0) {
      dyreco = TMath::Abs(BestZprimeCandidate->top_leptonic_v4().Rapidity()) - TMath::Abs(BestZprimeCandidate->top_hadronic_v4().Rapidity()); 
    } else {
      dyreco = TMath::Abs(BestZprimeCandidate->top_hadronic_v4().Rapidity()) - TMath::Abs(BestZprimeCandidate->top_leptonic_v4().Rapidity()); 
    }
    N_lep_charge->Fill(BestZprimeCandidate->lepton().charge(),weight);
    // cout <<"Lepton charge is: "<< BestZprimeCandidate->lepton().charge()<<endl;
    DeltaY_reco->Fill(dyreco, weight);
  
  
  //start spin correlation
  // ZprimeCandidate* BestZprimeCandidate = event.get(h_BestZprimeCandidateChi2); 
  bool is_toptag_reconstruction = BestZprimeCandidate->is_toptag_reconstruction(); // Reconstruction process id
  vector <Jet> AK4CHSjets_matched = event.get(h_CHSjets_matched);                  // AK4Puppijets that have been matched to CHSjets
  vector <TopJet> TopTaggedJets = event.get(h_AK8TopTags);                     // AK8Puppi jets TopTagged by DeepAK8TopTagger
  vector <float> jets_hadronic_bscores;                                            // bScores vector for resolved hadronic jets
  float pt_hadTop_thresh = 150;                                                    // Define cut-variable as pt of hadTop for low/high regions                                                   // medium WP for UL18 DeepJet
  float btag_WP=0.2783;  
  if (isUL16preVFP) btag_WP = 0.2598;                                              // medium WP for UL16preVFP DeepJet
  if (isUL16postVFP) btag_WP = 0.3657;                                             // medium WP for UL16postVFP DeepJet
  if (isUL17) btag_WP = 0.3040;                                                    // medium WP for UL17 DeepJet
  if (isUL18) btag_WP = 0.2783;                                                    // medium WP for UL18 DeepJet
   
                                                                // see https://btv-wiki.docs.cern.ch/ScaleFactors/ for btag WPs
                                                // medium WP for UL18 DeepJet

  // Plot pt of hadronic Top jet
  float pt_hadTop = BestZprimeCandidate->top_hadronic_v4().pt();
  
  float bscore_max = -2;
  if(!is_toptag_reconstruction){
      // Loop over resolved hadronic jets to find their bscore via CHS jets
    for(unsigned int i=0; i<BestZprimeCandidate->jets_hadronic().size(); i++){
      double deltaR_min = 99;
      // Match resolved hadronic jets to CHS jets (which have bscores)
      for(unsigned int j=0; j<AK4CHSjets_matched.size(); j++){
        double deltaR_CHS = deltaR(BestZprimeCandidate->jets_hadronic().at(i), AK4CHSjets_matched.at(j));
        if(deltaR_CHS < deltaR_min) deltaR_min = deltaR_CHS;}
      // Build bScore-vector for resolved hadronic jets whose bscore will correspond by index
      for(unsigned int k=0; k<AK4CHSjets_matched.size(); k++){
        if(deltaR(BestZprimeCandidate->jets_hadronic().at(i), AK4CHSjets_matched.at(k)) == deltaR_min) 
        jets_hadronic_bscores.emplace_back(AK4CHSjets_matched.at(k).btag_DeepJet());} // Using DeepJet btag score
    }
    // Loop over bScores-vector to extract highest bscore
    for(unsigned int i=0; i<jets_hadronic_bscores.size(); i++){
      float bscore = jets_hadronic_bscores.at(i);
      if(bscore > bscore_max) bscore_max = bscore;
    }
   
  }//is not top tag
  if(is_toptag_reconstruction){
      // Loop over hadronic top's subjets to extract highest bscore
    for(unsigned int i=0; i < BestZprimeCandidate->tophad_topjet_ptr()->subjets().size(); i++){
      float bscore = BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(i).btag_DeepJet(); // Using DeepJet btag score
      if(bscore > bscore_max) bscore_max = bscore;
    }
   
  }
  if(bscore_max >= btag_WP){
      // event.set(h_bscore_max, bscore_max); // Plot max bscores

      //------------------------------------Define 4vectors of hadronic b-jet and lepton------------------------------------//
      // Hadronic b-jet 4-vector
      TLorentzVector had_top_b(0, 0, 0, 0);

      // Resolved topology
      if(!is_toptag_reconstruction){ // Define hadronic b-jet as hadronic AK4-jet with highest bscore
        for(unsigned int i=0; i< BestZprimeCandidate->jets_hadronic().size(); i++){
          float bscore = jets_hadronic_bscores.at(i);
          if(bscore == bscore_max) had_top_b.SetPtEtaPhiE(BestZprimeCandidate->jets_hadronic().at(i).pt(), 
                                                          BestZprimeCandidate->jets_hadronic().at(i).eta(), 
                                                          BestZprimeCandidate->jets_hadronic().at(i).phi(), 
                                                          BestZprimeCandidate->jets_hadronic().at(i).energy());
        }
      }
      // Merged topology
      if(is_toptag_reconstruction){ // Define hadronic b-jet as hadronic AK8-subjet with highest bscore
        for(unsigned int j=0; j < BestZprimeCandidate->tophad_topjet_ptr()->subjets().size(); j++){
          float bscore = BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).btag_DeepJet();
          if(bscore == bscore_max) had_top_b.SetPtEtaPhiE(BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).pt(), 
                                                          BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).eta(), 
                                                          BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).phi(), 
                                                          BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).energy());
        }
      }

      // Lepton 4-vector
      TLorentzVector lep_top_lep(0, 0, 0, 0);
      LorentzVector lep = BestZprimeCandidate->lepton().v4();
      lep_top_lep.SetPtEtaPhiE(lep.pt(), lep.eta(), lep.phi(), lep.E());
      //------------------------------------Define 4vectors of hadronic b-jet and lepton------------------------------------//


      //-------------------------------- Begin boosting top quarks and their decay products --------------------------------//
      // Define 4vectors of top quarks
      TLorentzVector PosTop(0, 0, 0, 0);
      TLorentzVector NegTop(0, 0, 0, 0);

      // POSITIVE LEPTON CONFIGURATION => Positive charged lepton has Positive Top mother
      if(BestZprimeCandidate->lepton().charge() > 0){
        // Define ttbar system
        PosTop.SetPtEtaPhiE(BestZprimeCandidate->top_leptonic_v4().pt(), 
                            BestZprimeCandidate->top_leptonic_v4().eta(), 
                            BestZprimeCandidate->top_leptonic_v4().phi(), 
                            BestZprimeCandidate->top_leptonic_v4().energy());
        NegTop.SetPtEtaPhiE(BestZprimeCandidate->top_hadronic_v4().pt(), 
                            BestZprimeCandidate->top_hadronic_v4().eta(), 
                            BestZprimeCandidate->top_hadronic_v4().phi(), 
                            BestZprimeCandidate->top_hadronic_v4().energy());
      }
      else if (BestZprimeCandidate->lepton().charge() < 0)
      {
       PosTop.SetPtEtaPhiE(BestZprimeCandidate->top_hadronic_v4().pt(), 
                            BestZprimeCandidate->top_hadronic_v4().eta(), 
                            BestZprimeCandidate->top_hadronic_v4().phi(), 
                            BestZprimeCandidate->top_hadronic_v4().energy());
        NegTop.SetPtEtaPhiE(BestZprimeCandidate->top_leptonic_v4().pt(), 
                            BestZprimeCandidate->top_leptonic_v4().eta(), 
                            BestZprimeCandidate->top_leptonic_v4().phi(), 
                            BestZprimeCandidate->top_leptonic_v4().energy());
      }
      
        TLorentzVector ttbar = PosTop + NegTop;
        // Boost into ttbar CoM-Frame <<<-------//
        lep_top_lep.Boost(-ttbar.BoostVector());
        had_top_b.Boost(-ttbar.BoostVector());
        PosTop.Boost(-ttbar.BoostVector());
        NegTop.Boost(-ttbar.BoostVector());

        // Rotate vectors into Helicity Frame <<<------//
        // Rotate about beamline
        lep_top_lep.RotateZ(-1.*PosTop.Phi());
        had_top_b.RotateZ(-1.*PosTop.Phi());
        PosTop.RotateZ(-1.*PosTop.Phi());
        NegTop.RotateZ(-1.*PosTop.Phi());
        // Rotate about y-axis
        lep_top_lep.RotateY(-1.*PosTop.Theta());
        had_top_b.RotateY(-1.*PosTop.Theta());
        PosTop.RotateY(-1.*PosTop.Theta());
        NegTop.RotateY(-1.*PosTop.Theta());

        // Boost into ttbar Rest-Frame <<<--------//
        lep_top_lep.Boost(-PosTop.BoostVector()); // Positive charged lepton has Positive Top mother
        had_top_b.Boost(-NegTop.BoostVector());   // Positive charged lepton means b-jet has Negative Top mother

       //-------------------------------- End boosting top quarks and their decay products --------------------------------//
  
      // Define angular variables as sum and difference of decay products' phi-coordinates
      // sphi and dphi = PosTopDecayProd_phi +- NegTopDecayProd_phi
      float dphi=0.;
      float sphi = lep_top_lep.Phi() + had_top_b.Phi();
      if(BestZprimeCandidate->lepton().charge() > 0){ // lepton is Positive Top's Decay Product
        dphi = lep_top_lep.Phi() - had_top_b.Phi();
      }
      if(BestZprimeCandidate->lepton().charge() < 0)
      {
        dphi = had_top_b.Phi() - lep_top_lep.Phi();
      }
    
        // Map back into original domain if necessary
      if(sphi > TMath::Pi()) sphi = sphi - 2*TMath::Pi();
      if(sphi < -TMath::Pi()) sphi = sphi + 2*TMath::Pi();
      if(dphi > TMath::Pi()) dphi = dphi - 2*TMath::Pi();
      if(dphi < -TMath::Pi()) dphi = dphi + 2*TMath::Pi();
      Sigma_phi->Fill(sphi,weight);
      Delta_phi->Fill(dphi,weight);

        // Plot dphi and sphi for high-pt ranges
      if(pt_hadTop > pt_hadTop_thresh && dyreco >0){
        Sigma_phi_1->Fill(sphi,weight);
        Delta_phi_1->Fill(dphi,weight);
      }
      if(pt_hadTop > pt_hadTop_thresh && dyreco <0){
        Sigma_phi_2->Fill(sphi,weight);
        Delta_phi_2->Fill(dphi,weight);
      }
      if(pt_hadTop < pt_hadTop_thresh && sphi >0){
        DeltaY_reco_s1->Fill(dyreco,weight);
      }
      if(pt_hadTop < pt_hadTop_thresh && sphi <0){
        DeltaY_reco_s2->Fill(dyreco,weight);
      }
      if(pt_hadTop < pt_hadTop_thresh && dphi >0){
        DeltaY_reco_d1->Fill(dyreco,weight);
      }
      if(pt_hadTop < pt_hadTop_thresh && dphi <0){
        DeltaY_reco_d2->Fill(dyreco,weight);
      }
      if(pt_hadTop > pt_hadTop_thresh){
        Sigma_phi_high->Fill(sphi,weight);
        Delta_phi_high->Fill(dphi,weight);
        DeltaY_reco_high->Fill(dyreco,weight);
      }
      // Plot dphi and sphi for low-pt ranges
      if(pt_hadTop < pt_hadTop_thresh){
        Sigma_phi_low->Fill(sphi,weight);
        Delta_phi_low->Fill(dphi,weight);
        DeltaY_reco_low->Fill(dyreco,weight);
      }
    }










   //end spin correlation
  }

  
  if(is_zprime_reconstructed_chi2 && is_mc){
    ZprimeCandidate* BestZprimeCandidate = event.get(h_BestZprimeCandidateChi2);
    float Mreco = BestZprimeCandidate->Zprime_v4().M();
    float chi2 = BestZprimeCandidate->discriminator("chi2_total");

    LorentzVector toplep = BestZprimeCandidate->top_leptonic_v4();
    LorentzVector tophad = BestZprimeCandidate->top_hadronic_v4();

    toplep_pt->Fill(toplep.Pt(), weight);
    toplep_eta->Fill(toplep.Eta(), weight);
    toplep_phi->Fill(toplep.Phi(), weight);
    toplep_m->Fill(toplep.M(), weight);

    tophad_pt->Fill(tophad.Pt(), weight);
    tophad_eta->Fill(tophad.Eta(), weight);
    tophad_phi->Fill(tophad.Phi(), weight);
    tophad_m->Fill(tophad.M(), weight);

    //ditop_mass->Fill(Mreco, weight);
    ditop_absDeltaPhi->Fill(deltaPhi(toplep, tophad), weight);
    ditop_deltaEta->Fill(toplep.Eta()-tophad.Eta(), weight);
    ditop_absDeltaEta->Fill(abs(toplep.Eta()-tophad.Eta()), weight);
    ditop_deltaR->Fill(deltaR(toplep,tophad), weight);


    TOP_20_001_ditopmass_Fig19->Fill(Mreco, weight);
    TOP_20_001_ditopmass_Fig25->Fill(Mreco, weight);

    LorentzVector had_top = BestZprimeCandidate->top_hadronic_v4();
    LorentzVector lep_top = BestZprimeCandidate->top_leptonic_v4();
    
    TLorentzVector had_top_frame(0,0,0,0);
    had_top_frame.SetPtEtaPhiE(had_top.pt(), had_top.eta(), had_top.phi(), had_top.E());
    TLorentzVector lep_top_frame(0,0,0,0);
    lep_top_frame.SetPtEtaPhiE(lep_top.pt(), lep_top.eta(), lep_top.phi(), lep_top.E());
    TLorentzVector ttbar(0,0,0,0);
    ttbar.SetPtEtaPhiE((had_top+lep_top).pt(), (had_top+lep_top).eta(), (had_top+lep_top).phi(), (had_top+lep_top).E());

    had_top_frame.Boost(-ttbar.BoostVector());
    lep_top_frame.Boost(-ttbar.BoostVector());
    float costhetastar_had = ttbar.Vect().Dot(had_top_frame.Vect())/ttbar.P()/had_top_frame.P();
    float costhetastar_lep = ttbar.Vect().Dot(lep_top_frame.Vect())/ttbar.P()/lep_top_frame.P();

    cos_hadtop_thetastar->Fill(costhetastar_had, weight);
    cos_leptop_thetastar->Fill(costhetastar_lep, weight);
    abs_cos_leptop_thetastar->Fill(fabs(costhetastar_lep), weight);

    mttbar_vs_costhetastar->Fill(costhetastar_lep, Mreco, weight);
    costhetastar_vs_mttbar->Fill(Mreco, costhetastar_lep, weight);

    if(BestZprimeCandidate->is_toptag_reconstruction()){
      M_tophad->Fill(BestZprimeCandidate->tophad_topjet_ptr()->v4().M(), weight);
      M_toplep->Fill(inv_mass(BestZprimeCandidate->top_leptonic_v4()), weight);

      M_Zprime_ttag->Fill(Mreco, weight);
      M_Zprime_ttag_rebin->Fill(Mreco, weight);
      M_Zprime_ttag_rebin2->Fill(Mreco, weight);
      M_Zprime_ttag_rebin3->Fill(Mreco, weight);
      if(Mreco < 6000.) M_Zprime_ttag_rebin4->Fill(Mreco, weight);
      else M_Zprime_ttag_rebin4->Fill(6050., weight);
      if(Mreco < 5000.) M_Zprime_ttag_rebin5->Fill(Mreco, weight);
      else M_Zprime_ttag_rebin5->Fill(5050., weight);
      M_Zprime_ttag_rebin6->Fill(Mreco, weight);
      M_Zprime_ttag_rebin7->Fill(Mreco, weight);
      chi2_Zprime_ttag->Fill(chi2, weight);
      chi2_Zprime_ttag_rebin->Fill(chi2, weight);
      chi2_Zprime_ttag_rebin2->Fill(chi2, weight);

      LorentzVector SumSubjets(0.,0.,0.,0.);
      for(unsigned int k=0; k<BestZprimeCandidate->tophad_topjet_ptr()->subjets().size(); k++) SumSubjets = SumSubjets + BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(k).v4();
      float mhad = inv_mass(SumSubjets);

      M_tophad_ttag->Fill(mhad, weight);
      M_tophad_ttag_v4->Fill(inv_mass(BestZprimeCandidate->top_hadronic_v4()), weight);
      M_tophad_ttag_M->Fill(BestZprimeCandidate->top_hadronic_v4().M(), weight);

      //M_tophad_ttag->Fill(BestZprimeCandidate->tophad_topjet_ptr()->softdropmass(), weight);
      M_toplep_ttag->Fill(inv_mass(BestZprimeCandidate->top_leptonic_v4()), weight);
    }
    else{
      M_tophad->Fill(inv_mass(BestZprimeCandidate->top_hadronic_v4()), weight);
      M_toplep->Fill(inv_mass(BestZprimeCandidate->top_leptonic_v4()), weight);

      M_Zprime_ak4->Fill(Mreco, weight);
      M_Zprime_ak4_rebin->Fill(Mreco, weight);
      M_Zprime_ak4_rebin2->Fill(Mreco, weight);
      M_Zprime_ak4_rebin3->Fill(Mreco, weight);
      if(Mreco < 6000.) M_Zprime_ak4_rebin4->Fill(Mreco, weight);
      else M_Zprime_ak4_rebin4->Fill(6050., weight);
      if(Mreco < 5000.) M_Zprime_ak4_rebin5->Fill(Mreco, weight);
      else M_Zprime_ak4_rebin5->Fill(5050., weight);
      M_Zprime_ak4_rebin6->Fill(Mreco, weight);
      M_Zprime_ak4_rebin7->Fill(Mreco, weight);
      chi2_Zprime_ak4->Fill(chi2, weight);
      chi2_Zprime_ak4_rebin->Fill(chi2, weight);
      chi2_Zprime_ak4_rebin2->Fill(chi2, weight);
      M_tophad_ak4->Fill(inv_mass(BestZprimeCandidate->top_hadronic_v4()), weight);
      M_toplep_ak4->Fill(inv_mass(BestZprimeCandidate->top_leptonic_v4()), weight);
    }
  }
  if(is_zprime_reconstructed_correctmatch){
    // cout << "Correct match is filled" << endl;
    ZprimeCandidate* BestZprimeCandidate = event.get(h_BestZprimeCandidateCorrectMatch);
    float Mreco = BestZprimeCandidate->Zprime_v4().M();
    float dr = BestZprimeCandidate->discriminator("correct_match");
    if(dr < 10.){
      // cout << "dr < 10" << endl;
      if(BestZprimeCandidate->is_toptag_reconstruction()){
        M_tophad_dr_ttag->Fill(BestZprimeCandidate->tophad_topjet_ptr()->softdropmass(), weight);
        M_toplep_dr_ttag->Fill(inv_mass(BestZprimeCandidate->top_leptonic_v4()), weight);
      }
      else{
        M_tophad_dr_ak4->Fill(inv_mass(BestZprimeCandidate->top_hadronic_v4()), weight);
        M_toplep_dr_ak4->Fill(inv_mass(BestZprimeCandidate->top_leptonic_v4()), weight);
      }
      dr_discr_Zprime->Fill(dr, weight);
      M_Zprime_dr->Fill(Mreco, weight);
      M_Zprime_dr_rebin->Fill(Mreco, weight);
      M_Zprime_dr_rebin2->Fill(Mreco, weight);
      M_Zprime_dr_rebin3->Fill(Mreco, weight);
    }
  }

  // Sphericity tensor
  double s11 = -1., s12 = -1., s13 = -1., s22 = -1., s23 = -1., s33 = -1., mag = -1.;
  for(const Jet jet : *event.jets){
    mag += (jet.v4().Px()*jet.v4().Px()+jet.v4().Py()*jet.v4().Py()+jet.v4().Pz()*jet.v4().Pz());
    s11 += jet.v4().Px()*jet.v4().Px();
    s12 += jet.v4().Px()*jet.v4().Py();
    s13 += jet.v4().Px()*jet.v4().Pz();
    s22 += jet.v4().Py()*jet.v4().Py();
    s23 += jet.v4().Py()*jet.v4().Pz();
    s33 += jet.v4().Pz()*jet.v4().Pz();
  }

  s11 = s11 / mag;
  s12 = s12 / mag;
  s13 = s13 / mag;
  s22 = s22 / mag;
  s23 = s23 / mag;
  s33 = s33 / mag;

  S11->Fill(s11, weight);
  S12->Fill(s12, weight);
  S13->Fill(s13, weight);
  S22->Fill(s22, weight);
  S23->Fill(s23, weight);
  S33->Fill(s33, weight);

  sum_event_weights->Fill(1., weight);

  // N_Jets_vs_HT->Fill(Njets, st_jets, weight);

  /*
  ███    ██ ███    ██
  ████   ██ ████   ██
  ██ ██  ██ ██ ██  ██
  ██  ██ ██ ██  ██ ██
  ██   ████ ██   ████
  */

  if(debug) cout << "before NN in hists" << endl;
  if(NN){
    if(debug) cout << "is it going inside NN" << endl;
    for(int i=0; i<Nmuons; i++){
      NN_Mu_pt->Fill(muons->at(i).pt(),weight);
      NN_Mu_eta->Fill(muons->at(i).eta(),weight);
      NN_Mu_phi->Fill(muons->at(i).phi(),weight);
      NN_Mu_E->Fill(muons->at(i).energy(),weight);
    }



    for(int i=0; i<Nelectrons; i++){
      NN_Ele_pt->Fill(electrons->at(i).pt(),weight);
      NN_Ele_eta->Fill(electrons->at(i).eta(),weight);
      NN_Ele_phi->Fill(electrons->at(i).phi(),weight);
      NN_Ele_E->Fill(electrons->at(i).energy(),weight);
    }

    NN_MET_pt->Fill(event.met->pt(),weight);
    NN_MET_phi->Fill(event.met->phi(),weight);

    vector<Jet>* Ak4jets = event.jets;
    int NAk4jets = Ak4jets->size();
    NN_N_Ak4->Fill(NAk4jets,weight);

    for(int i=0; i<NAk4jets; i++){
      if(i==0){
        NN_Ak4_j1_pt->Fill(Ak4jets->at(i).pt(),weight);
        NN_Ak4_j1_eta->Fill(Ak4jets->at(i).eta(),weight);
        NN_Ak4_j1_phi->Fill(Ak4jets->at(i).phi(),weight);
        NN_Ak4_j1_E->Fill(Ak4jets->at(i).energy(),weight);
        NN_Ak4_j1_m->Fill(Ak4jets->at(i).v4().M(),weight);
        NN_Ak4_j1_btag->Fill(Ak4jets->at(i).btag_DeepJet(),weight);
      }
      if(i==1){
        NN_Ak4_j2_pt->Fill(Ak4jets->at(i).pt(),weight);
        NN_Ak4_j2_eta->Fill(Ak4jets->at(i).eta(),weight);
        NN_Ak4_j2_phi->Fill(Ak4jets->at(i).phi(),weight);
        NN_Ak4_j2_E->Fill(Ak4jets->at(i).energy(),weight);
        NN_Ak4_j2_m->Fill(Ak4jets->at(i).v4().M(),weight);
        NN_Ak4_j2_btag->Fill(Ak4jets->at(i).btag_DeepJet(),weight);
      }
      if(i==2){
        NN_Ak4_j3_pt->Fill(Ak4jets->at(i).pt(),weight);
        NN_Ak4_j3_eta->Fill(Ak4jets->at(i).eta(),weight);
        NN_Ak4_j3_phi->Fill(Ak4jets->at(i).phi(),weight);
        NN_Ak4_j3_E->Fill(Ak4jets->at(i).energy(),weight);
        NN_Ak4_j3_m->Fill(Ak4jets->at(i).v4().M(),weight);
        NN_Ak4_j3_btag->Fill(Ak4jets->at(i).btag_DeepJet(),weight);
      }
      if(i==3){
        NN_Ak4_j4_pt->Fill(Ak4jets->at(i).pt(),weight);
        NN_Ak4_j4_eta->Fill(Ak4jets->at(i).eta(),weight);
        NN_Ak4_j4_phi->Fill(Ak4jets->at(i).phi(),weight);
        NN_Ak4_j4_E->Fill(Ak4jets->at(i).energy(),weight);
        NN_Ak4_j4_m->Fill(Ak4jets->at(i).v4().M(),weight);
        NN_Ak4_j4_btag->Fill(Ak4jets->at(i).btag_DeepJet(),weight);
      }
      if(i==4){
        NN_Ak4_j5_pt->Fill(Ak4jets->at(i).pt(),weight);
        NN_Ak4_j5_eta->Fill(Ak4jets->at(i).eta(),weight);
        NN_Ak4_j5_phi->Fill(Ak4jets->at(i).phi(),weight);
        NN_Ak4_j5_E->Fill(Ak4jets->at(i).energy(),weight);
        NN_Ak4_j5_m->Fill(Ak4jets->at(i).v4().M(),weight);
        NN_Ak4_j5_btag->Fill(Ak4jets->at(i).btag_DeepJet(),weight);
      }
      if(i==5){
        NN_Ak4_j6_pt->Fill(Ak4jets->at(i).pt(),weight);
        NN_Ak4_j6_eta->Fill(Ak4jets->at(i).eta(),weight);
        NN_Ak4_j6_phi->Fill(Ak4jets->at(i).phi(),weight);
        NN_Ak4_j6_E->Fill(Ak4jets->at(i).energy(),weight);
        NN_Ak4_j6_m->Fill(Ak4jets->at(i).v4().M(),weight);
        NN_Ak4_j6_btag->Fill(Ak4jets->at(i).btag_DeepJet(),weight);
      }
    }

    if(ishotvr){
      vector<TopJet>* HOTVRjets = event.topjets;
      int N_HOTVRjets = HOTVRjets->size();
      NN_N_HOTVR->Fill(N_HOTVRjets,weight);

      for(int i=0; i<N_HOTVRjets; i++){
        if(i==0){
          NN_HOTVR_j1_pt->Fill(HOTVRjets->at(i).pt(),weight);
          NN_HOTVR_j1_eta->Fill(HOTVRjets->at(i).eta(),weight);
          NN_HOTVR_j1_phi->Fill(HOTVRjets->at(i).phi(),weight);
          NN_HOTVR_j1_E->Fill(HOTVRjets->at(i).energy(),weight);
          NN_HOTVR_j1_mSD->Fill(HOTVRjets->at(i).v4().M(),weight);
          NN_HOTVR_j1_tau21->Fill(HOTVRjets->at(i).tau2_groomed()/HOTVRjets->at(i).tau1_groomed(),weight);
          NN_HOTVR_j1_tau32->Fill(HOTVRjets->at(i).tau3_groomed()/HOTVRjets->at(i).tau2_groomed(),weight);
        }
        if(i==1){
          NN_HOTVR_j2_pt->Fill(HOTVRjets->at(i).pt(),weight);
          NN_HOTVR_j2_eta->Fill(HOTVRjets->at(i).eta(),weight);
          NN_HOTVR_j2_phi->Fill(HOTVRjets->at(i).phi(),weight);
          NN_HOTVR_j2_E->Fill(HOTVRjets->at(i).energy(),weight);
          NN_HOTVR_j2_mSD->Fill(HOTVRjets->at(i).v4().M(),weight);
          NN_HOTVR_j2_tau21->Fill(HOTVRjets->at(i).tau2_groomed()/HOTVRjets->at(i).tau1_groomed(),weight);
          NN_HOTVR_j2_tau32->Fill(HOTVRjets->at(i).tau3_groomed()/HOTVRjets->at(i).tau2_groomed(),weight);
        }
        if(i==2){
          NN_HOTVR_j3_pt->Fill(HOTVRjets->at(i).pt(),weight);
          NN_HOTVR_j3_eta->Fill(HOTVRjets->at(i).eta(),weight);
          NN_HOTVR_j3_phi->Fill(HOTVRjets->at(i).phi(),weight);
          NN_HOTVR_j3_E->Fill(HOTVRjets->at(i).energy(),weight);
          NN_HOTVR_j3_mSD->Fill(HOTVRjets->at(i).v4().M(),weight);
          NN_HOTVR_j3_tau21->Fill(HOTVRjets->at(i).tau2_groomed()/HOTVRjets->at(i).tau1_groomed(),weight);
          NN_HOTVR_j3_tau32->Fill(HOTVRjets->at(i).tau3_groomed()/HOTVRjets->at(i).tau2_groomed(),weight);
        }
      }
    } // end hotvr mode

    if(isdeepAK8){
      vector<TopJet>* Ak8jets = event.toppuppijets;
      int NAk8jets = Ak8jets->size();
      NN_N_Ak8->Fill(NAk8jets,weight);

      for(int i=0; i<NAk8jets; i++){
        if(i==0){
          NN_Ak8_j1_pt->Fill(Ak8jets->at(i).pt(),weight);
          NN_Ak8_j1_eta->Fill(Ak8jets->at(i).eta(),weight);
          NN_Ak8_j1_phi->Fill(Ak8jets->at(i).phi(),weight);
          NN_Ak8_j1_E->Fill(Ak8jets->at(i).energy(),weight);
          NN_Ak8_j1_mSD->Fill(Ak8jets->at(i).softdropmass(),weight);
          NN_Ak8_j1_tau21->Fill(Ak8jets->at(i).tau2()/Ak8jets->at(i).tau1(),weight);
          NN_Ak8_j1_tau32->Fill(Ak8jets->at(i).tau3()/Ak8jets->at(i).tau2(),weight);
          NN_Ak8_j1_ttag->Fill(Ak8jets->at(i).btag_MassDecorrelatedDeepBoosted_TvsQCD(),weight);
        }
        if(i==1){
          NN_Ak8_j2_pt->Fill(Ak8jets->at(i).pt(),weight);
          NN_Ak8_j2_eta->Fill(Ak8jets->at(i).eta(),weight);
          NN_Ak8_j2_phi->Fill(Ak8jets->at(i).phi(),weight);
          NN_Ak8_j2_E->Fill(Ak8jets->at(i).energy(),weight);
          NN_Ak8_j2_mSD->Fill(Ak8jets->at(i).softdropmass(),weight);
          NN_Ak8_j2_tau21->Fill(Ak8jets->at(i).tau2()/Ak8jets->at(i).tau1(),weight);
          NN_Ak8_j2_tau32->Fill(Ak8jets->at(i).tau3()/Ak8jets->at(i).tau2(),weight);
          NN_Ak8_j2_ttag->Fill(Ak8jets->at(i).btag_MassDecorrelatedDeepBoosted_TvsQCD(),weight);
        }
        if(i==2){
          NN_Ak8_j3_pt->Fill(Ak8jets->at(i).pt(),weight);
          NN_Ak8_j3_eta->Fill(Ak8jets->at(i).eta(),weight);
          NN_Ak8_j3_phi->Fill(Ak8jets->at(i).phi(),weight);
          NN_Ak8_j3_E->Fill(Ak8jets->at(i).energy(),weight);
          NN_Ak8_j3_mSD->Fill(Ak8jets->at(i).softdropmass(),weight);
          NN_Ak8_j3_tau21->Fill(Ak8jets->at(i).tau2()/Ak8jets->at(i).tau1(),weight);
          NN_Ak8_j3_tau32->Fill(Ak8jets->at(i).tau3()/Ak8jets->at(i).tau2(),weight);
          NN_Ak8_j3_ttag->Fill(Ak8jets->at(i).btag_MassDecorrelatedDeepBoosted_TvsQCD(),weight);
        }
      }
    } // end deepAK8 mode
  

    if(is_zprime_reconstructed_chi2){
      ZprimeCandidate* BestZprimeCandidate = event.get(h_BestZprimeCandidateChi2);
      float Mass_tt = BestZprimeCandidate->Zprime_v4().M();
      float chi2 = BestZprimeCandidate->discriminator("chi2_total");
      if(is_mc) NN_M_tt_weighted->Fill(Mass_tt,weight);
      if(is_mc) NN_M_tt_notweighted->Fill(Mass_tt);
      NN_chi2->Fill(chi2,weight);
    }
  //NN
  }
  if(debug) cout << "after NN in hists" << endl;
 
  

} //Method



ZprimeSemiLeptonicHists::~ZprimeSemiLeptonicHists(){}