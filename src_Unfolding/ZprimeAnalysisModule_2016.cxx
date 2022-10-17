#include <iostream>
#include <memory>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>
#include "UHH2/common/include/PrintingModules.h"

#include <UHH2/common/include/CleaningModules.h>
#include <UHH2/common/include/NSelections.h>
#include <UHH2/common/include/LumiSelection.h>
#include <UHH2/common/include/TriggerSelection.h>
#include <UHH2/common/include/JetCorrections.h>
#include <UHH2/common/include/ObjectIdUtils.h>
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/ElectronIds.h>
#include <UHH2/common/include/JetIds.h>
#include <UHH2/common/include/TopJetIds.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/AdditionalSelections.h>
#include "UHH2/common/include/LuminosityHists.h"
#include <UHH2/common/include/MCWeight.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/common/include/EventHists.h>

#include <UHH2/ZprimeSemiLeptonic/include/ModuleBASE.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSelections.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicModules.h>
#include <UHH2/ZprimeSemiLeptonic/include/TTbarLJHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicGeneratorHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeCandidate.h>

#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/TTbarReconstruction.h>
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>

using namespace std;
using namespace uhh2;

/*
██████  ███████ ███████ ██ ███    ██ ██ ████████ ██  ██████  ███    ██
██   ██ ██      ██      ██ ████   ██ ██    ██    ██ ██    ██ ████   ██
██   ██ █████   █████   ██ ██ ██  ██ ██    ██    ██ ██    ██ ██ ██  ██
██   ██ ██      ██      ██ ██  ██ ██ ██    ██    ██ ██    ██ ██  ██ ██
██████  ███████ ██      ██ ██   ████ ██    ██    ██  ██████  ██   ████
*/

class ZprimeAnalysisModule_2016 : public ModuleBASE {

public:
  explicit ZprimeAnalysisModule_2016(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  void book_histograms(uhh2::Context&, vector<string>);
  void fill_histograms(uhh2::Event&, string);

protected:

  // Cleaners (to make sure the pre-selection did everything right)
  std::unique_ptr<MuonCleaner> muon_cleaner;
  std::unique_ptr<ElectronCleaner> electron_cleaner;

  // Scale Factors -- Systematics
  unique_ptr<MCMuonScaleFactor> MuonID_module, MuonTrigger_module;

  // AnalysisModules
//  BTagWeight_module.reset(new MCCSVv2ShapeSystematic(ctx, "jets","central","iterativefit","","MCCSVv2ShapeSystematic"));
  unique_ptr<AnalysisModule> LumiWeight_module, PUWeight_module, BTagWeight_module, printer_genparticles;

  // Taggers
  unique_ptr<AK8PuppiTopTagger> TopTaggerPuppi;

  // Mass reconstruction
  unique_ptr<ZprimeCandidateBuilder> CandidateBuilder;

  // Chi2 discriminator
  unique_ptr<ZprimeChi2Discriminator> Chi2DiscriminatorZprime;
  unique_ptr<ZprimeCorrectMatchDiscriminator> CorrectMatchDiscriminatorZprime;

  // Selections
  unique_ptr<Selection> Trigger1_selection, Trigger2_selection, NMuon1_selection, NMuon2_selection, NElectron_selection, TwoDCut_selection, Jet1_selection, Jet2_selection, Met_selection, Chi2_selection, TTbarMatchable_selection, Chi2CandidateMatched_selection, ZprimeTopTag_selection, BlindData_selection;
  std::unique_ptr<uhh2::Selection> met_sel;
  std::unique_ptr<Selection> sel_0btag, sel_1btag, sel_2btag;
  std::unique_ptr<Selection> TopJetBtagSubjet_selection;


  //Handles
  Event::Handle<bool> h_is_zprime_reconstructed_chi2, h_is_zprime_reconstructed_correctmatch;
  Event::Handle<float> h_chi2;  
  Event::Handle<float> h_weight;
  Event::Handle<float> h_MET;
  Event::Handle<int> h_NPV;
  Event::Handle<int> h_ttagN;
  Event::Handle<int> h_btagN;
  Event::Handle<float> h_Mttbar; 
  Event::Handle<float> h_DeltaY;
  Event::Handle<float> h_pTttbar;
  Event::Handle<float> h_Rapidityttbar;
  Event::Handle<float> h_pThad;
  Event::Handle<float> h_pTlep;
  Event::Handle<float> h_Rapidityhad;
  Event::Handle<float> h_Rapiditylep;
  Event::Handle<float> h_DeltaPhi_thad_lepton;
  Event::Handle<float> h_Cos_theta_thad_tlep;
  Event::Handle<float> h_tlead_pT;
  Event::Handle<float> h_tsublead_pT;
  uhh2::Event::Handle<ZprimeCandidate*> h_BestZprimeCandidateChi2;
  uhh2::Event::Handle<std::vector<TopJet>> h_AK8TopTags;


  // Configuration
  bool isMC, ispuppi, islooserselection;
  string Sys_MuonID, Sys_MuonTrigger, Sys_PU, Sys_btag;
  TString sample;
  int runnr_oldtriggers = 299368;

  bool is2016v2, is2016v3, is2017, is2018;
  bool isMuon, isElectron;
};

void ZprimeAnalysisModule_2016::book_histograms(uhh2::Context& ctx, vector<string> tags){
  for(const auto & tag : tags){
    string mytag = tag + "_Skimming";
    // book_HFolder(mytag, new TTbarLJHistsSkimming(ctx,mytag));
    mytag = tag+"_General";
    book_HFolder(mytag, new ZprimeSemiLeptonicHists(ctx,mytag));
    mytag = tag+"_Muons";
    book_HFolder(mytag, new MuonHists(ctx,mytag));
    mytag = tag+"_Electrons";
    book_HFolder(mytag, new ElectronHists(ctx,mytag));
    mytag = tag+"_Jets";
    book_HFolder(mytag, new JetHists(ctx,mytag));
    mytag = tag+"_Event";
    book_HFolder(mytag, new EventHists(ctx,mytag));
    mytag = tag+"_Generator";
    book_HFolder(mytag, new ZprimeSemiLeptonicGeneratorHists(ctx,mytag));
  }
}

void ZprimeAnalysisModule_2016::fill_histograms(uhh2::Event& event, string tag){
  string mytag = tag + "_Skimming";
  // HFolder(mytag)->fill(event);
  mytag = tag+"_General";
  HFolder(mytag)->fill(event);
  mytag = tag+"_Muons";
  HFolder(mytag)->fill(event);
  mytag = tag+"_Electrons";
  HFolder(mytag)->fill(event);
  mytag = tag+"_Jets";
  HFolder(mytag)->fill(event);
  mytag = tag+"_Event";
  HFolder(mytag)->fill(event);
  mytag = tag+"_Generator";
  HFolder(mytag)->fill(event);
}

/*
█  ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████  ██████  ██████
█ ██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█ ██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██    ██    ██ ██████
█ ██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█  ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██     ██████  ██   ██
*/

ZprimeAnalysisModule_2016::ZprimeAnalysisModule_2016(uhh2::Context& ctx){

  for(auto & kv : ctx.get_all()){
    cout << " " << kv.first << " = " << kv.second << endl;
  }

  // Important selection values
  islooserselection = (ctx.get("is_looser_selection") == "true");
  double muon_pt(50.);
  double elec_pt(80.);
  double jet1_pt(150.);
  double jet2_pt(50.);
  //  double stlep_plus_met(150.);
  double chi2_max(30.);
  double mtt_blind(3000.);
  int nmuon_min1, nmuon_max1;
  int nmuon_min2, nmuon_max2;
  int nele_min, nele_max;
  string trigger1,trigger2;// trigger3;
  double MET_cut;
  isMuon = false; isElectron = false;
  if(ctx.get("channel") == "muon") isMuon = true;
  if(ctx.get("channel") == "electron") isElectron = true;

  if(isMuon){//semileptonic muon channel
    trigger1 = "HLT_Mu50_v*";
    trigger2 = "HLT_Mu50_v*";
    //  trigger2 = "HLT_TrkMu50_v*";
  //  trigger3 = "HLT_TkMu50_v*";
    nmuon_min1 = 1, nmuon_max1 = -1;
    nmuon_min2 = 1, nmuon_max2 = 1;
    nele_min = 0; nele_max = 0;
    MET_cut = 50;
    jet1_pt = 150.;
  }
  if(isElectron){//semileptonic electron channel
    nmuon_min1 = 0; nmuon_max1 = 0;
    nmuon_min2 = 0; nmuon_max2 = 0;
    nele_min = 1; nele_max = 1;
    trigger1 = "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v*";
    trigger2 = "HLT_Ele115_CaloIdVT_GsfTrkIdT_v*";
    MET_cut = 120;
    jet1_pt = 185.;
    //    trigger3 = "HLT\_Ele115_v*";
  }

 
  double TwoD_dr = 0.4, TwoD_ptrel = 25.;
  if(islooserselection){
    jet1_pt = 100.;
    TwoD_dr = 0.2;
    TwoD_ptrel = 10.;
    //    stlep_plus_met = 100.;
  }
  const MuonId muonID(PtEtaCut(muon_pt, 2.4));
  const ElectronId electronID(PtEtaSCCut(elec_pt, 2.5));

  // Configuration
  isMC = (ctx.get("dataset_type") == "MC");
  ispuppi = (ctx.get("is_puppi") == "true");
  TString mode = "chs";
  if(ispuppi) mode = "puppi";
  string tmp = ctx.get("dataset_version");
  sample = tmp;
  is2016v2 = (ctx.get("dataset_version").find("2016v2") != std::string::npos);
  is2016v3 = (ctx.get("dataset_version").find("2016v3") != std::string::npos);
  is2017 = (ctx.get("dataset_version").find("2017") != std::string::npos);
  is2018 = (ctx.get("dataset_version").find("2018") != std::string::npos);

  Sys_MuonID = ctx.get("Sys_MuonID");
  Sys_MuonTrigger = ctx.get("Sys_MuonTrigger");
  Sys_PU = ctx.get("Sys_PU");
  Sys_btag = ctx.get("Sys_BTagSF");


  //estos tres funcionan
  BTag::algo btag_algo = BTag::CSVV2;
  BTag::wp btag_wp_tight = BTag::WP_TIGHT;
  JetId id_btag = BTag(btag_algo, btag_wp_tight);

//  BTag::algo btag_algo = BTag::DEEPJET;
//  BTag::wp btag_wp_tight = BTag::WP_TIGHT;
//  JetId id_btag = BTag(btag_algo, btag_wp_tight);


  // Modules
  printer_genparticles.reset(new GenParticlesPrinter(ctx));
  muon_cleaner.reset(new MuonCleaner(muonID));
  electron_cleaner.reset(new ElectronCleaner(electronID));
  LumiWeight_module.reset(new MCLumiWeight(ctx));
  PUWeight_module.reset(new MCPileupReweight(ctx, Sys_PU));
  BTagWeight_module.reset(new MCBTagDiscriminantReweighting(ctx, btag_algo, "jets", Sys_btag));

  if(is2016v3){
     MuonID_module.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/hugobg/UHH2/CMSSW_10_2_10/src/UHH2/common/data/2016/MuonID_EfficienciesAndSF_average_RunBtoH.root", "NUM_HighPtID_DEN_genTracks_eta_pair_newTuneP_probe_pt", 0., "HighPtID", true, Sys_MuonID));
     MuonTrigger_module.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/hugobg/UHH2/CMSSW_10_2_10/src/UHH2/common/data/2017/MuonTrigger_EfficienciesAndSF_RunBtoF_Nov17Nov2017.root", "Mu50_PtEtaBins/pt_abseta_ratio", 0.5, "Trigger", true, Sys_MuonTrigger));
  }
  // if(is2017 || is2018){
  //   MuonID_module.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/karavdia/CMSSW_10_2_11/src/UHH2/common/data/2017/MuonID_94X_RunBCDEF_SF_ID.root", "NUM_HighPtID_DEN_genTracks_pair_newTuneP_probe_pt_abseta", 0., "HighPtID", true, Sys_MuonID));
  //   MuonTrigger_module.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/karavdia/CMSSW_10_2_11/src/UHH2/common/data/2017/MuonTrigger_EfficienciesAndSF_RunBtoF_Nov17Nov2017.root", "Mu50_PtEtaBins/pt_abseta_ratio", 0.5, "Trigger", true, Sys_MuonTrigger));
  // }

  // Selection modules
  Trigger1_selection.reset(new TriggerSelection(trigger1));
  Trigger2_selection.reset(new TriggerSelection(trigger2));
  //  Trigger3_selection.reset(new TriggerSelection(trigger3));
  NMuon1_selection.reset(new NMuonSelection(nmuon_min1, nmuon_max1));
  NMuon2_selection.reset(new NMuonSelection(nmuon_min2, nmuon_max2));
  NElectron_selection.reset(new NElectronSelection(nele_min, nele_max));
  TwoDCut_selection.reset(new TwoDCut(TwoD_dr, TwoD_ptrel));
  Jet1_selection.reset(new NJetSelection(1, -1, JetId(PtEtaCut(jet1_pt, 2.4))));
  Jet2_selection.reset(new NJetSelection(2, -1, JetId(PtEtaCut(jet2_pt, 2.4))));
  //  STlepPlusMet_selection.reset(new STlepPlusMetCut(stlep_plus_met, -1.));
  met_sel.reset(new METCut  (MET_cut   , uhh2::infinity));
  Chi2_selection.reset(new Chi2Cut(ctx, 0., chi2_max));
  TTbarMatchable_selection.reset(new TTbarSemiLepMatchableSelection());
  Chi2CandidateMatched_selection.reset(new Chi2CandidateMatchedSelection(ctx));
  ZprimeTopTag_selection.reset(new ZprimeTopTagSelection(ctx));
  BlindData_selection.reset(new BlindDataSelection(ctx, mtt_blind));

  // Taggers
  TopTaggerPuppi.reset(new AK8PuppiTopTagger(ctx));

  // Zprime candidate builder
  CandidateBuilder.reset(new ZprimeCandidateBuilder(ctx, mode));

  // Zprime discriminators
  Chi2DiscriminatorZprime.reset(new ZprimeChi2Discriminator(ctx));
  h_is_zprime_reconstructed_chi2 = ctx.get_handle<bool>("is_zprime_reconstructed_chi2");
  CorrectMatchDiscriminatorZprime.reset(new ZprimeCorrectMatchDiscriminator(ctx));
  h_is_zprime_reconstructed_correctmatch = ctx.get_handle<bool>("is_zprime_reconstructed_correctmatch");
  h_BestZprimeCandidateChi2 = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestChi2");
  h_AK8TopTags = ctx.get_handle<std::vector<TopJet>>("AK8PuppiTopTags");

  h_chi2 = ctx.declare_event_output<float> ("rec_chi2");
  h_weight = ctx.declare_event_output<float> ("weight"); 
  h_MET = ctx.declare_event_output<float> ("met_pt");
  h_ttagN = ctx.declare_event_output<int>("ttagN");
  h_btagN = ctx.declare_event_output<int>("btagN");
  h_Mttbar = ctx.declare_event_output<float> ("Mttbar");
  h_DeltaY = ctx.declare_event_output<float> ("DeltaY");
  h_pTttbar = ctx.declare_event_output<float> ("pT_ttbar");
  h_Rapidityttbar= ctx.declare_event_output<float> ("Rapidity_ttbar");
  h_pThad= ctx.declare_event_output<float> ("pT_had");
  h_pTlep= ctx.declare_event_output<float> ("pT_lep");
  h_Rapidityhad= ctx.declare_event_output<float> ("Rapidity_had");
  h_Rapiditylep= ctx.declare_event_output<float> ("Rapidity_lep");
  h_DeltaPhi_thad_lepton = ctx.declare_event_output<float> ("DeltaPhi_thad_lepton");
  h_Cos_theta_thad_tlep = ctx.declare_event_output<float> ("Cos_theta_thad_tlep"); 
  h_tlead_pT = ctx.declare_event_output<float> ("tlead_pT");
  h_tsublead_pT = ctx.declare_event_output<float> ("tsublead_pT");
//  h_lep1_pt = ctx.declare_event_output<float> ("lep1_pt");
//  h_lep1_eta = ctx.declare_event_output<float> ("lep1_eta");
//  h_ak4jet1_pt = ctx.declare_event_output<float> ("ak4jet1_pt");
//  h_ak4jet1_eta = ctx.declare_event_output<float> ("ak4jet1_eta");
//  h_ak8jet1_pt = ctx.declare_event_output<float> ("ak8jet1_pt");
//  h_ak8jet1_eta = ctx.declare_event_output<float> ("ak8jet1_eta");



  // btag 
  //CSVBTag::wp btag_wp = CSVBTag::WP_TIGHT; // b-tag workingpoint
  //JetId id_btag = CSVBTag(btag_wp);

  // DeepCSVBTag::wp btag_wp = DeepCSVBTag::WP_TIGHT; // b-tag workingpoint
  // JetId id_btag = DeepCSVBTag(btag_wp);

  //DeepJetBTag::wp btag_wp = DeepJetBTag::WP_TIGHT; // b-tag workingpoint
  //JetId id_btag = DeepJetBTag(btag_wp);

  sel_1btag.reset(new NJetSelection(1, 1, id_btag));
  sel_2btag.reset(new NJetSelection(2, -1, id_btag));


  TopJetBtagSubjet_selection.reset(new ZprimeBTagFatSubJetSelection(ctx));

  // Book histograms
  vector<string> histogram_tags = {"Weights", "Muon1", "Trigger", "Muon2", "Electron1", "TwoDCut", "Jet1", "Jet2", "MET", "MatchableBeforeChi2Cut", "NotMatchableBeforeChi2Cut", "CorrectMatchBeforeChi2Cut", "NotCorrectMatchBeforeChi2Cut", "Chi2", "Matchable", "NotMatchable", "CorrectMatch", "NotCorrectMatch", "TopTagReconstruction", "NotTopTagReconstruction", "Btags2", "Btags1","TopJetBtagSubjet"};
  book_histograms(ctx, histogram_tags);
}

/*
██████  ██████   ██████   ██████ ███████ ███████ ███████
██   ██ ██   ██ ██    ██ ██      ██      ██      ██
██████  ██████  ██    ██ ██      █████   ███████ ███████
██      ██   ██ ██    ██ ██      ██           ██      ██
██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/

bool ZprimeAnalysisModule_2016::process(uhh2::Event& event){

  //  cout << "++++++++++++ NEW EVENT ++++++++++++++" << endl;

  // Initialize reco flags with false
  event.set(h_is_zprime_reconstructed_chi2, false);
  event.set(h_is_zprime_reconstructed_correctmatch, false);
  event.set(h_chi2,-100);
  event.set(h_weight,-100);
  event.set(h_ttagN,-100);
  event.set(h_btagN,-100);  
  event.set(h_MET,-100);
  event.set(h_Mttbar,-100);
  event.set(h_DeltaY,-100);

  event.set(h_pTttbar,-100);
  event.set(h_Rapidityttbar,-100);
  event.set(h_pThad,-100);
  event.set(h_pTlep,-100);
  event.set(h_Rapidityhad,-100);
  event.set(h_Rapiditylep,-100);
  event.set(h_DeltaPhi_thad_lepton,-100);
  event.set(h_Cos_theta_thad_tlep,-100);
  event.set(h_tlead_pT,-100);
  event.set(h_tsublead_pT,-100);

//  event.set(h_lep1_pt,-100);
//  event.set(h_lep1_eta,-100);
//  event.set(h_ak4jet1_pt,-100);
//  event.set(h_ak4jet1_eta,-100);
//  event.set(h_ak8jet1_pt,-100);
//  event.set(h_ak8jet1_eta,-100);
//  event.set(h_NPV,-100);


  // Printing
  // if(!event.isRealData) printer_genparticles->process(event);

  // TODO Apply things that should've been done in the pre-selection already... Fix pre-selection and then remove these steps
  if(isMuon) muon_cleaner->process(event);
  if(isElectron) electron_cleaner->process(event);

  // Weight modules
  LumiWeight_module->process(event);
  PUWeight_module->process(event);
  BTagWeight_module->process(event);
  MuonID_module->process(event);
  
  // Run top-tagging
  TopTaggerPuppi->process(event);
  fill_histograms(event, "Weights");

  if(!NMuon1_selection->passes(event)) return false;
  fill_histograms(event, "Muon1");

  // if(isMC || (!isMC && event.run >= runnr_oldtriggers)){
  //   if(!(Trigger1_selection->passes(event) || Trigger2_selection->passes(event) || Trigger3_selection->passes(event))) return false;
  // }
  // else{
  //   if(!(Trigger1_selection->passes(event))) return false;
  // }
  if(!(Trigger1_selection->passes(event)|| Trigger2_selection->passes(event))) return false;
  MuonTrigger_module->process_onemuon(event, 0);
  fill_histograms(event, "Trigger");

  if(!NMuon2_selection->passes(event)) return false;
  fill_histograms(event, "Muon2");

  if(!NElectron_selection->passes(event)) return false;
  fill_histograms(event, "Electron1");

  if(!TwoDCut_selection->passes(event)) return false;
  fill_histograms(event, "TwoDCut");

  // Here, the Zprime must be reconstructed (we ensured to have >= 2 AK4 jets, >= 1 muon)
  CandidateBuilder->process(event);
  Chi2DiscriminatorZprime->process(event);
  CorrectMatchDiscriminatorZprime->process(event);

  if(sample.Contains("_blinded")){
    if(!BlindData_selection->passes(event)) return false;
  }
  if(!Jet1_selection->passes(event)) return false;
  fill_histograms(event, "Jet1");

  if(!Jet2_selection->passes(event)) return false;
  fill_histograms(event, "Jet2");

  // if(!STlepPlusMet_selection->passes(event)) return false;
  // fill_histograms(event, "STlepPlusMet");

  // MET selection
  if(!met_sel->passes(event)) return false;
  fill_histograms(event, "MET");

  if(TTbarMatchable_selection->passes(event)) fill_histograms(event, "MatchableBeforeChi2Cut");
  else fill_histograms(event, "NotMatchableBeforeChi2Cut");
  if(Chi2CandidateMatched_selection->passes(event)) fill_histograms(event, "CorrectMatchBeforeChi2Cut");
  else fill_histograms(event, "NotCorrectMatchBeforeChi2Cut");

  if(!Chi2_selection->passes(event)) return false;
  fill_histograms(event, "Chi2");

  if(TTbarMatchable_selection->passes(event)) fill_histograms(event, "Matchable");
  else fill_histograms(event, "NotMatchable");

  if(Chi2CandidateMatched_selection->passes(event)) fill_histograms(event, "CorrectMatch");
  else fill_histograms(event, "NotCorrectMatch");

  if(ZprimeTopTag_selection->passes(event)) fill_histograms(event, "TopTagReconstruction");
  else fill_histograms(event, "NotTopTagReconstruction");


  //Test with b-tagging
  if(sel_2btag->passes(event)) fill_histograms(event, "Btags2");
  if(sel_1btag->passes(event)) fill_histograms(event, "Btags1");

  if(TopJetBtagSubjet_selection->passes(event)) fill_histograms(event, "TopJetBtagSubjet");

  // get ttagN
 
  vector<TopJet> TopTags = event.get(h_AK8TopTags);
  event.set(h_ttagN,TopTags.size());

  // get btagN
  
  int btagN = 0;
  int jetbtagN = 0;
  int btag_subjet=0;

  if(TopTags.size() == 0){
     if(sel_2btag->passes(event)) jetbtagN = jetbtagN + 2;
     if(sel_1btag->passes(event)) jetbtagN = jetbtagN + 1; 
  }

  if(TopTags.size() == 1){
     for(auto & topjet : *event.topjets){
        auto subjets = topjet.subjets();
        for (auto & subjet : subjets) {
          if(sel_1btag->passes(event)) btag_subjet++;
        }
     

        for(unsigned int j=0; j<event.jets->size(); j++){
           double dR_jet = deltaR(topjet, event.jets->at(j));
           if(dR_jet > 1.2) continue;
           if(sel_1btag->passes(event)) jetbtagN++;
        }
     }
  }

//  veto on ttagN > 1
//  if(TopTags.size() == 2) return false;
  
  btagN = jetbtagN + btag_subjet;

//  veto on btag > 1 
//  if(btagN > 1) return false;
  event.set(h_btagN,btagN);

  //Fill some vars for monitoring and comparison
  bool is_zprime_reconstructed_chi2 = event.get(h_is_zprime_reconstructed_chi2);

  if(is_zprime_reconstructed_chi2){
    ZprimeCandidate* BestZprimeCandidate = event.get(h_BestZprimeCandidateChi2);
    //    float chi2 = BestZprimeCandidate->discriminator("chi2_total");
    event.set(h_chi2,BestZprimeCandidate->discriminator("chi2_total"));
    event.set(h_Mttbar,BestZprimeCandidate->Zprime_v4().M());

    if (isMuon){
        if(event.muons->at(0).charge() == 1){
                event.set(h_DeltaY,TMath::Abs(BestZprimeCandidate->top_leptonic_v4().Rapidity()) - TMath::Abs(BestZprimeCandidate->top_hadronic_v4().Rapidity()));  
        }
        if(event.muons->at(0).charge() == -1){
       	       	event.set(h_DeltaY,TMath::Abs(BestZprimeCandidate->top_hadronic_v4().Rapidity()) - TMath::Abs(BestZprimeCandidate->top_leptonic_v4().Rapidity()));
        }
          
        event.set(h_DeltaPhi_thad_lepton,deltaPhi(event.muons->at(0),BestZprimeCandidate->top_hadronic_v4()));   
    }
    event.set(h_pTttbar,BestZprimeCandidate->Zprime_v4().Pt());
    event.set(h_Rapidityttbar,BestZprimeCandidate->Zprime_v4().Rapidity()); 
    event.set(h_pThad,BestZprimeCandidate->top_hadronic_v4().Pt());
    event.set(h_pTlep,BestZprimeCandidate->top_leptonic_v4().Pt());
    event.set(h_Rapidityhad,BestZprimeCandidate->top_hadronic_v4().Rapidity());
    event.set(h_Rapiditylep,BestZprimeCandidate->top_leptonic_v4().Rapidity());    

    double mag_lep = TMath::Sqrt(BestZprimeCandidate->top_leptonic_v4().Px()*BestZprimeCandidate->top_leptonic_v4().Px() + BestZprimeCandidate->top_leptonic_v4().Py()*BestZprimeCandidate->top_leptonic_v4().Py() + BestZprimeCandidate->top_leptonic_v4().Pz()*BestZprimeCandidate->top_leptonic_v4().Pz());
    double mag_had = TMath::Sqrt(BestZprimeCandidate->top_hadronic_v4().Px()*BestZprimeCandidate->top_hadronic_v4().Px() + BestZprimeCandidate->top_hadronic_v4().Py()*BestZprimeCandidate->top_hadronic_v4().Py() + BestZprimeCandidate->top_hadronic_v4().Pz()*BestZprimeCandidate->top_hadronic_v4().Pz());
    double dot_prod = BestZprimeCandidate->top_hadronic_v4().Px()*BestZprimeCandidate->top_leptonic_v4().Px() + BestZprimeCandidate->top_hadronic_v4().Py()*BestZprimeCandidate->top_leptonic_v4().Py() + BestZprimeCandidate->top_hadronic_v4().Pz()*BestZprimeCandidate->top_leptonic_v4().Pz();   
    event.set(h_Cos_theta_thad_tlep,dot_prod/(mag_had*mag_lep));

    if( BestZprimeCandidate->top_leptonic_v4().Pt() > BestZprimeCandidate->top_hadronic_v4().Pt()){ 
    	event.set(h_tlead_pT,BestZprimeCandidate->top_leptonic_v4().Pt()); 
    	event.set(h_tsublead_pT,BestZprimeCandidate->top_hadronic_v4().Pt());
    }else{
        event.set(h_tlead_pT,BestZprimeCandidate->top_hadronic_v4().Pt());
        event.set(h_tsublead_pT,BestZprimeCandidate->top_leptonic_v4().Pt());
    }

  }

  event.set(h_weight,event.weight);
  event.set(h_MET,event.met->pt());
  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(ZprimeAnalysisModule_2016)
