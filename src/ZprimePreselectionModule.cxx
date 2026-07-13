#include <iostream>
#include <iomanip>
#include <limits>
#include <memory>
#include <map>
#include <cmath>
#include <algorithm>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>
#include "UHH2/common/include/PrintingModules.h"

#include <UHH2/common/include/CleaningModules.h>
#include <UHH2/common/include/NSelections.h>
#include <UHH2/common/include/LumiSelection.h>
#include <UHH2/common/include/TriggerSelection.h>
#include <UHH2/common/include/JetCorrections.h>
#include <UHH2/common/include/JetCorrectionSets.h>
#include <UHH2/common/include/ObjectIdUtils.h>
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/ElectronIds.h>
#include <UHH2/common/include/JetIds.h>
#include <UHH2/common/include/TopJetIds.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/core/include/GenParticle.h>
#include <set>
#include <UHH2/common/include/AdditionalSelections.h>
#include "UHH2/common/include/LuminosityHists.h"
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/common/include/EventHists.h>
#include <UHH2/common/include/CommonModules.h>
#include <UHH2/common/include/MCWeight.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <stdexcept>
// #include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
// #include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

#include <UHH2/ZprimeSemiLeptonic/include/ModuleBASE.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicModules.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSelections.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicPreselectionHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicGeneratorHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/CHSJetCorrections.h>
#include <UHH2/ZprimeSemiLeptonic/include/TopPuppiJetCorrections.h>
#include "UHH2/HOTVR/include/HOTVRJetCorrectionModule.h"

using namespace std;
using namespace uhh2;

namespace {
  constexpr float kBFragNominalInputRb = 0.855f;

  // True if tau has e or mu in its decay chain (W->tau->e/mu). False if tau->hadronic or unknown.
  bool tau_decays_to_emumu(const std::vector<GenParticle>* gps, int tau_index, std::set<int>& visited) {
    if (!gps || gps->empty() || tau_index < 0) return false;
    if (visited.count(tau_index)) return false;
    visited.insert(tau_index);
    const unsigned short uinv = (unsigned short)(-1);
    for (const auto& gp : *gps) {
      const bool from_m1 = (gp.mother1() != uinv && (int)gp.mother1() == tau_index);
      const bool from_m2 = (gp.mother2() != uinv && (int)gp.mother2() == tau_index);
      if (!from_m1 && !from_m2) continue;
      const int id = std::abs(gp.pdgId());
      if (id == 11 || id == 13) return true;
      if (id == 15 && tau_decays_to_emumu(gps, gp.index(), visited)) return true;
    }
    return false;
  }

  // b-fragmentation : check if the particle is a b-hadron, b-quark, or b-hadron-quark.

  bool is_b_hadron_pdgid(int pdgid) {
    // This removes the sign, so both B hadrons and anti-B hadrons are accepted.
    const int id = std::abs(pdgid);
    // This removes simple particles like quarks and leptons.
    if (id <= 100) return false;
    // Then it checks whether the PDG ID contains a 5 in the hadron quark-content digits.
    return ((id / 100) % 10 == 5) || ((id / 1000) % 10 == 5) || ((id / 10000) % 10 == 5);
  }

  const GenParticle* find_genparticle_by_index(const std::vector<GenParticle>* gps, int index) {
    if (!gps || index < 0) return nullptr;
    for (const auto& gp : *gps) {
      if ((int)gp.index() == index) return &gp;
    }
    return nullptr;
  }

  bool has_ancestor_index(const std::vector<GenParticle>* gps, const GenParticle& gp, int ancestor_index, std::set<int>& visited) {
    if (!gps || ancestor_index < 0) return false;
    if (visited.count(gp.index())) return false;
    visited.insert(gp.index());

    const unsigned short uinv = (unsigned short)(-1);
    for (const auto mother_index : {gp.mother1(), gp.mother2()}) {
      if (mother_index == uinv) continue;
      if ((int)mother_index == ancestor_index) return true;
      const GenParticle* mother = find_genparticle_by_index(gps, mother_index);
      if (mother && has_ancestor_index(gps, *mother, ancestor_index, visited)) return true;
    }
    return false;
  }

  bool has_b_hadron_daughter(const std::vector<GenParticle>* gps, const GenParticle& gp) {
    if (!gps) return false;
    const unsigned short uinv = (unsigned short)(-1);
    for (const auto& candidate : *gps) {
      if (!is_b_hadron_pdgid(candidate.pdgId())) continue;
      const bool from_m1 = (candidate.mother1() != uinv && (int)candidate.mother1() == (int)gp.index());
      const bool from_m2 = (candidate.mother2() != uinv && (int)candidate.mother2() == (int)gp.index());
      if (from_m1 || from_m2) return true;
    }
    return false;
  }

  std::vector<const GenParticle*> find_last_b_hadrons_from_b(const std::vector<GenParticle>* gps, const GenParticle& b_quark) {
    std::vector<const GenParticle*> out;
    if (!gps || b_quark.index() < 0) return out;

    for (const auto& gp : *gps) {
      if (!is_b_hadron_pdgid(gp.pdgId())) continue;
      std::set<int> visited;
      if (!has_ancestor_index(gps, gp, b_quark.index(), visited)) continue;
      if (has_b_hadron_daughter(gps, gp)) continue; // MiniAOD fallback: last-copy B hadron
      out.push_back(&gp);
    }
    return out;
  }

  double lorentz_dot(const LorentzVector& a, const LorentzVector& b) {
    return a.E() * b.E() - a.Px() * b.Px() - a.Py() * b.Py() - a.Pz() * b.Pz();
  }

  // b-fragmentation : compute the xb value for the b-hadron: more detail: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MLReweighting?extralog=-%20caching%20topic#Bfragm
  float compute_bfrag_xb(const GenParticle& b_hadron, const GenParticle& top, const GenParticle& w_boson) {
    const double mt = top.v4().M();
    const double mw = w_boson.v4().M();
    if (!std::isfinite(mt) || mt <= 0.0 || !std::isfinite(mw)) return std::numeric_limits<float>::quiet_NaN();

    const double mt2 = mt * mt;
    const double w = (mw * mw) / mt2;
    const double denom = 1.0 - w;
    if (!std::isfinite(denom) || denom <= 0.0) return std::numeric_limits<float>::quiet_NaN();

    const double xE = 2.0 * lorentz_dot(b_hadron.v4(), top.v4()) / mt2;
    const double xb = xE / denom;
    if (!std::isfinite(xb)) return std::numeric_limits<float>::quiet_NaN();
    return static_cast<float>(xb);
  }

  struct BFragInputs {
    float xb_top = std::numeric_limits<float>::quiet_NaN();
    float xb_antitop = std::numeric_limits<float>::quiet_NaN();
    int valid = 0;
  };

  // b-fragmentation
  BFragInputs compute_bfrag_inputs(const std::vector<GenParticle>* gps, const TTbarGen& ttbargen) {
    BFragInputs out;
    if (!gps || gps->empty() || ttbargen.DecayChannel() == TTbarGen::e_notfound) return out;

    // ttbargen.bTop() bquark from t->bW+ 
    const auto b_hadrons_top = find_last_b_hadrons_from_b(gps, ttbargen.bTop());
    const auto b_hadrons_antitop = find_last_b_hadrons_from_b(gps, ttbargen.bAntitop());

    // check if there is only one b-hadron from the top and one b-hadron from the antitop
    if (b_hadrons_top.size() != 1 || b_hadrons_antitop.size() != 1) return out;
    if (b_hadrons_top.front()->index() == b_hadrons_antitop.front()->index()) return out;

    out.xb_top = compute_bfrag_xb(*b_hadrons_top.front(), ttbargen.Top(), ttbargen.WTop());
    out.xb_antitop = compute_bfrag_xb(*b_hadrons_antitop.front(), ttbargen.Antitop(), ttbargen.WAntitop());

    const bool xb_ok = std::isfinite(out.xb_top) && std::isfinite(out.xb_antitop)
                    && out.xb_top >= 0.f && out.xb_antitop >= 0.f
                    && out.xb_top <= 1.2f && out.xb_antitop <= 1.2f;
    out.valid = xb_ok ? 1 : 0;
    return out;
  }
  // EW uncertainty costheta info at gen level
  TLorentzVector to_tlorentzvector(const LorentzVector& v4) {
    TLorentzVector out;
    out.SetPxPyPzE(v4.Px(), v4.Py(), v4.Pz(), v4.E());
    return out;
  }

  float compute_top_costheta_ttbarframe(const TTbarGen& ttbargen) {
    if (ttbargen.DecayChannel() == TTbarGen::e_notfound) return std::numeric_limits<float>::quiet_NaN();

    const TLorentzVector top = to_tlorentzvector(ttbargen.Top().v4());
    const TLorentzVector antitop = to_tlorentzvector(ttbargen.Antitop().v4());
    const TLorentzVector ttbar = top + antitop;
    if (!std::isfinite(ttbar.E()) || ttbar.E() <= 0.0) return std::numeric_limits<float>::quiet_NaN();

    TLorentzVector top_ttbar_frame = top;
    top_ttbar_frame.Boost(-ttbar.BoostVector());
    if (!std::isfinite(top_ttbar_frame.P()) || top_ttbar_frame.P() <= 0.0) return std::numeric_limits<float>::quiet_NaN();

    double costheta = top_ttbar_frame.Vect().Unit().Dot(TVector3(0.0, 0.0, 1.0));
    if (!std::isfinite(costheta)) return std::numeric_limits<float>::quiet_NaN();
    costheta = std::max(-1.0, std::min(1.0, costheta));
    return static_cast<float>(costheta);
  }

  // b-fragmentation : check if the particle is a b-hadron, b-quark, or b-hadron-quark.
  // True if any tau in the ttbar decay goes tau->hadronic. Exclude such events.
  bool has_any_hadronic_tau_decay(const std::vector<GenParticle>* gps, const TTbarGen& ttbargen) {
    if (!gps || gps->empty()) return true;  // conservative: exclude when uncertain
    const auto dc = ttbargen.DecayChannel();
    if (dc != TTbarGen::e_tauhad && dc != TTbarGen::e_tautau && dc != TTbarGen::e_etau && dc != TTbarGen::e_mutau)
      return false;
    std::set<int> visited;
    if (dc == TTbarGen::e_tauhad) {
      const auto& tau = ttbargen.ChargedLepton();
      if (!tau_decays_to_emumu(gps, tau.index(), visited)) return true;
      return false;
    }
    if (dc == TTbarGen::e_tautau) {
      for (const auto& wd : {ttbargen.Wdecay1(), ttbargen.Wdecay2(), ttbargen.WMinusdecay1(), ttbargen.WMinusdecay2()}) {
        if (std::abs(wd.pdgId()) == 15 && !tau_decays_to_emumu(gps, wd.index(), visited)) return true;
      }
      return false;
    }
    if (dc == TTbarGen::e_etau || dc == TTbarGen::e_mutau) {
      for (const auto& wd : {ttbargen.Wdecay1(), ttbargen.Wdecay2(), ttbargen.WMinusdecay1(), ttbargen.WMinusdecay2()}) {
        if (std::abs(wd.pdgId()) == 15) {
          if (!tau_decays_to_emumu(gps, wd.index(), visited)) return true;
          return false;
        }
      }
    }
    return false;
  }
}

class ZprimePreselectionModule : public ModuleBASE {

public:
  explicit ZprimePreselectionModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  void book_histograms(uhh2::Context&, vector<string>);
  void fill_histograms(uhh2::Event&, string);

protected:
  bool debug;

  // mttbar mass bin lower edges (must match noac_mtt_edges in SystematicsHists for NoAC gen histograms)
  const std::vector<double> mttbar_bin_edges = {0., 500., 750., 1000., 1500.};
  const std::vector<std::string> mttbar_bin_tags = {"0_500", "500_750", "750_1000", "1000_1500", "1500Inf"};
  // Fine bins: [0-500), [500-750), [750-1000), [1000-1500), [1500, Inf)
  // Additionally book/fill mtt_gen_0_750 = merged [0, 750) (overlaps first two fine bins).

  // Helper function to find mttbar bin
  inline int find_mtt_bin(double mtt) {
    for (size_t i = 0; i+1 < mttbar_bin_edges.size(); ++i) {
      if (mtt >= mttbar_bin_edges[i] && mtt < mttbar_bin_edges[i+1]) return static_cast<int>(i);
    }
    if (mtt >= mttbar_bin_edges.back()) return static_cast<int>(mttbar_bin_edges.size() - 1);
    return -1;
  }

  // Standalone lumi weight (applied between raw GEN hists and lumi-scaled GEN hists)
  std::unique_ptr<MCLumiWeight> lumi_weight;

  // Optional NNLO QCD + EWK correction of the GEN xi histograms used to build NoAC weights.
  // Gated by XML flag ApplyNNLOEWKToGenHists (default false): when off, nothing below is
  // constructed, no extra output branches are declared, and the normal preselection ->
  // analysis -> DNN chain is completely unchanged. When on, a parallel "_lumiscaled_nnloewk"
  // folder set is filled with genWeight * lumi * EWK * NNLO-QCD, intended only for producing
  // the NoAC gen-template file (the actual per-event correction stays at the applyNN stage).
  bool apply_nnloewk_to_gen_hists = false;
  std::unique_ptr<AnalysisModule> ttbar_ewk_gen;
  std::unique_ptr<AnalysisModule> ttbar_nnloqcd_gen;
  uhh2::Event::Handle<float> h_w_ewk_nom;
  uhh2::Event::Handle<float> h_w_nnlo_nom;

  // Corrections
  std::unique_ptr<CommonModules> common;
  std::unique_ptr<AnalysisModule> hotvrjetCorr;
  std::unique_ptr<TopPuppiJetCorrections> toppuppijetCorr;
  std::unique_ptr<CHSJetCorrections> CHSjetCorr;

  // Cleaners
  std::unique_ptr<JetCleaner>      jet_IDcleaner, jet_cleaner1, jet_cleaner2;
  std::unique_ptr<AnalysisModule>  hotvrjet_cleaner;
  std::unique_ptr<TopJetCleaner>   topjet_puppi_IDcleaner, topjet_puppi_cleaner;

  // Selections
  std::unique_ptr<uhh2::Selection> genflavor_sel;
  std::unique_ptr<uhh2::Selection> jet1_sel;
  std::unique_ptr<uhh2::Selection> jet2_sel;
  std::unique_ptr<uhh2::Selection> met_sel;
  unique_ptr<Selection> SignSplit;

  bool isMC, isHOTVR;
  string Sys_PU;

  std::unique_ptr<Hists> lumihists;
  TString METcollection;

  bool isUL16preVFP, isUL16postVFP, isUL17, isUL18;

  // additional branch with AK4 CHS jets -> for b-tagging
  Event::Handle<vector<Jet>> h_CHSjets;

  // TTbarGen handle for mttbar calculation
  Event::Handle<TTbarGen> h_ttbargen;
  std::unique_ptr<TTbarGenProducer> ttgenprod;
  uhh2::Event::Handle<float> h_xi_gen;
  uhh2::Event::Handle<float> h_mtt_gen;
  uhh2::Event::Handle<float> h_costheta_gen;
  uhh2::Event::Handle<float> h_DeltaY_gen;
  // b-fragmentation 
  uhh2::Event::Handle<float> h_bfrag_xb_top;
  uhh2::Event::Handle<float> h_bfrag_xb_antitop;
  uhh2::Event::Handle<float> h_bfrag_rb;
  uhh2::Event::Handle<int> h_bfrag_valid;
  uhh2::Event::Handle<float> h_weight_bfrag_nom;
  uhh2::Event::Handle<float> h_weight_bfrag_up;

  // Cutflow histograms: survive hadd across Condor jobs
  // h_cf_raw:         bin content = N_raw (unweighted event count)
  // h_cf_genweight:   bin content = Sum(genWeight) before lumi scaling
  // h_cf_lumiscaled:  bin content = Sum(genWeight * lumiScale) after lumi scaling
  // Bin labels: 0_AllEvents, 1_AfterLumiWeight, 2_GEN_semilep, 3_AfterCommonModules,
  //             4_Lepton1, 5_Jet1, 6_Jet2, 7_MET
  TH1D *h_cf_raw, *h_cf_genweight, *h_cf_lumiscaled;
  static const int NCF = 8; // number of cutflow stages
  double lumi_factor = 1.0; // target_lumi / dataset_lumi, cached for cutflow

  // POWHEG EWK lookup ingredients (sigma_POWHEG vs m_ttbar, costheta*), filled with
  // raw genWeight before any selection and for ALL ttbar decays (no tau veto).
  // Enabled by setting PowhegEWKLookupBinningFile; binning cloned from EWno_0.
  // Harvest offline: sigma_POWHEG_bin = powheg_yield * xsec / powheg_sumw.
  TH2D* h_powheg_lookup_yield = nullptr;
  TH1D* h_powheg_lookup_sumw = nullptr;

};

void ZprimePreselectionModule::book_histograms(uhh2::Context& ctx, vector<string> tags){
  for(const auto & tag : tags){
    string mytag = tag+"_General";
    book_HFolder(mytag, new ZprimeSemiLeptonicPreselectionHists(ctx,mytag));
  }
}

void ZprimePreselectionModule::fill_histograms(uhh2::Event& event, string tag){
  string mytag = tag+"_General";
  HFolder(mytag)->fill(event);
}

ZprimePreselectionModule::ZprimePreselectionModule(uhh2::Context& ctx) {

  debug = false; // true/false

  for(auto & kv : ctx.get_all()){
    cout << " " << kv.first << " = " << kv.second << endl;
  }

  //// CONFIGURATION
  const TString METcollection = ctx.get("METName");
  isMC    = ctx.get("dataset_type") == "MC";
  isHOTVR = ctx.get("is_HOTVR") == "true";
  Sys_PU  = ctx.get("Sys_PU");

  isUL16preVFP  = (ctx.get("dataset_version").find("UL16preVFP")  != std::string::npos);
  isUL16postVFP = (ctx.get("dataset_version").find("UL16postVFP") != std::string::npos);
  isUL17        = (ctx.get("dataset_version").find("UL17")        != std::string::npos);
  isUL18        = (ctx.get("dataset_version").find("UL18")        != std::string::npos);

  // lepton IDs
  ElectronId eleID_veto = ElectronTagID(Electron::mvaEleID_Fall17_noIso_V2_wp90);
  MuonId     muID_veto  = MuonID(Muon::CutBasedIdTight);

  double electron_pt(25.);
  double muon_pt(25.);
  double jet1_pt(30.);
  double jet2_pt(30.);
  double MET(20.);

  // GEN Flavor selection [W+jets flavor-splitting]
  if(ctx.get("dataset_version").find("WJets") != std::string::npos){
    if     (ctx.get("dataset_version").find("_B") != std::string::npos) genflavor_sel.reset(new GenFlavorSelection("b"));
    else if(ctx.get("dataset_version").find("_C") != std::string::npos) genflavor_sel.reset(new GenFlavorSelection("c"));
    else if(ctx.get("dataset_version").find("_L") != std::string::npos) genflavor_sel.reset(new GenFlavorSelection("l"));
    else genflavor_sel.reset(new uhh2::AndSelection(ctx));
  }
  else genflavor_sel.reset(new uhh2::AndSelection(ctx));

  // Cleaning: Mu, Ele, Jets
  const MuonId muonID_veto(AndId<Muon>(PtEtaCut(muon_pt, 2.4), muID_veto));
  const ElectronId electronID_veto(AndId<Electron>(PtEtaSCCut(electron_pt, 2.5), eleID_veto));
  const JetPFID jetID_CHS(JetPFID::WP_TIGHT_CHS);
  const JetPFID jetID_PUPPI(JetPFID::WP_TIGHT_PUPPI);

  jet_IDcleaner.reset(new JetCleaner(ctx, jetID_PUPPI));
  jet_cleaner1.reset(new JetCleaner(ctx, 15., 3.0));
  jet_cleaner2.reset(new JetCleaner(ctx, 20., 2.5));
  hotvrjet_cleaner.reset(new TopJetCleaner(ctx, PtEtaCut(200., 2.5)));
  topjet_puppi_IDcleaner.reset(new TopJetCleaner(ctx, jetID_PUPPI, "toppuppijets"));
  topjet_puppi_cleaner.reset(new TopJetCleaner(ctx, TopJetId(PtEtaCut(200., 2.5)), "toppuppijets"));

  // Split interference signal samples by sign
  if(ctx.get("dataset_version").find("_int") != std::string::npos){
    if     (ctx.get("dataset_version").find("_pos") != std::string::npos) SignSplit.reset(new SignSelection("pos"));
    else if(ctx.get("dataset_version").find("_neg") != std::string::npos) SignSplit.reset(new SignSelection("neg"));
    else SignSplit.reset(new uhh2::AndSelection(ctx));
  }
  else SignSplit.reset(new uhh2::AndSelection(ctx));

  // Standalone lumi weight: applied AFTER raw GEN hists, BEFORE lumi-scaled GEN hists
  if(isMC) lumi_weight.reset(new MCLumiWeight(ctx));

  // Optional: NNLO QCD + EWK correction applied (after lumi) to a parallel GEN folder set,
  // used only to build NoAC gen templates. Off by default so normal jobs are untouched.
  // Reuses the "ttbargen" handle produced above by TTbarGenProducer. Both modules apply to
  // event.weight only for ttbar datasets (gated internally on dataset_version).
  apply_nnloewk_to_gen_hists = isMC && (ctx.get("ApplyNNLOEWKToGenHists", "false") == "true");
  if(apply_nnloewk_to_gen_hists){
    ttbar_ewk_gen.reset(new TTbarEWKCorrection(ctx, "ttbargen"));
    ttbar_nnloqcd_gen.reset(new TTbarNNLOQCDReweighting(ctx, "ttbargen"));
    // Read back the stored nominal weights (set by the modules above) to apply only to the
    // GEN _nnloewk folders, without touching the global event.weight.
    h_w_ewk_nom  = ctx.get_handle<float>("weight_ttbar_ewk_nominal");
    h_w_nnlo_nom = ctx.get_handle<float>("weight_ttbar_nnlo_qcd");
  }

  // Cache lumi factor for cutflow
  if(isMC){
    double dataset_lumi = std::abs(string2double(ctx.get("dataset_lumi")));
    double target_lumi  = string2double(ctx.get("target_lumi"));
    lumi_factor = target_lumi / dataset_lumi;
    cout << "=== LUMI SCALING INFO ===" << endl;
    cout << "  target_lumi  = " << target_lumi  << " pb^-1" << endl;
    cout << "  dataset_lumi = " << dataset_lumi  << " pb^-1" << endl;
    cout << "  lumi_factor  = " << lumi_factor << endl;
    cout << "=========================" << endl;
  }

  // common modules (lumi weight DISABLED - we apply it standalone between the two GEN hist sets)
  common.reset(new CommonModules());
  common->disable_mclumiweight();   // <-- avoid double lumi scaling
  common->switch_jetlepcleaner(true);
  // common->disable_pvfilter();
  common->disable_jetpfidfilter();
  common->switch_jetPtSorter(true);
  common->switch_metcorrection(true);
  common->set_muon_id(muonID_veto);
  common->set_electron_id(electronID_veto);
  common->init(ctx, Sys_PU);

  hotvrjetCorr.reset(new HOTVRJetCorrectionModule(ctx));

  toppuppijetCorr.reset(new TopPuppiJetCorrections());
  toppuppijetCorr->init(ctx);

  CHSjetCorr.reset(new CHSJetCorrections());
  CHSjetCorr->init(ctx);

  // TTbarGen producer. Never throw: this runs on every MC sample, and anything
  // without exactly one top and one antitop (WJets, ST, DY, QCD, toponium EtaT)
  // must yield DecayChannel == e_notfound instead of killing the job.
  if(isMC) ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));

  //// EVENT SELECTION
  jet1_sel.reset(new NJetSelection(1, -1, JetId(PtEtaCut(jet1_pt, 2.5))));
  jet2_sel.reset(new NJetSelection(2, -1, JetId(PtEtaCut(jet2_pt, 2.5))));
  met_sel.reset(new METCut(MET, uhh2::infinity));

  // additional branch with Ak4 CHS jets
  h_CHSjets = ctx.get_handle<vector<Jet>>("jetsAk4CHS");

  // TTbarGen handle for mttbar calculation
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");

  // GEN-level outputs (so they exist in the event and output tree)
  if (isMC) {
    h_xi_gen     = ctx.declare_event_output<float>("xi_gen");
    h_mtt_gen    = ctx.declare_event_output<float>("mtt_gen");
    h_costheta_gen = ctx.declare_event_output<float>("costheta_gen");
    h_DeltaY_gen = ctx.declare_event_output<float>("DeltaY_gen");
    // b-fragmentation 
    h_bfrag_xb_top      = ctx.declare_event_output<float>("bfrag_xb_top");
    h_bfrag_xb_antitop  = ctx.declare_event_output<float>("bfrag_xb_antitop");
    h_bfrag_rb          = ctx.declare_event_output<float>("bfrag_rb");
    h_bfrag_valid       = ctx.declare_event_output<int>("bfrag_valid");
    h_weight_bfrag_nom  = ctx.declare_event_output<float>("weight_bfrag_nom");
    h_weight_bfrag_up   = ctx.declare_event_output<float>("weight_bfrag_up");
  }

  // Book cutflow histograms (survive hadd across Condor jobs)
  {
    const vector<string> cf_labels = {
      "0_AllEvents", "1_GEN_raw", "2_AfterLumiWeight", "3_GEN_lumiscaled",
      "4_AfterCommonModules", "5_Lepton1", "6_Jet1", "7_MET"
    };
    auto book_cf = [&](const char* name, const char* title) -> TH1D* {
      TH1D* h = new TH1D(name, title, NCF, 0, NCF);
      h->Sumw2();
      for(int i = 0; i < NCF; ++i) h->GetXaxis()->SetBinLabel(i+1, cf_labels[i].c_str());
      ctx.put("cutflow", h);  // attach to output file under cutflow/ directory
      return h;
    };
    h_cf_raw        = book_cf("cf_raw",        "Cutflow: unweighted event count");
    h_cf_genweight  = book_cf("cf_genweight",  "Cutflow: sum of genWeight (before lumi)");
    h_cf_lumiscaled = book_cf("cf_lumiscaled", "Cutflow: sum of lumi-scaled weight");
  }

  // Book histograms
  // --- Raw GEN folders (genWeight only, no lumi) ---
  vector<string> histogram_tags = {"Input", "mtt_gen_inclusive"};
    // "CommonModules", "HOTVRCorrections", "PUPPICorrections", "Lepton1", "JetID", "JetCleaner1", "JetCleaner2", "TopjetCleaner", "Jet1", "Jet2", "MET"};
    
  // Add mttbar bin tags
  for (const auto & tag : mttbar_bin_tags) {
    histogram_tags.push_back("mtt_gen_" + tag);
  }
  histogram_tags.push_back("mtt_gen_0_750");

  // --- Lumi-scaled GEN folders (genWeight * lumi, no PU) ---
  histogram_tags.push_back("mtt_gen_inclusive_lumiscaled");
  for (const auto & tag : mttbar_bin_tags) {
    histogram_tags.push_back("mtt_gen_" + tag + "_lumiscaled");
  }
  histogram_tags.push_back("mtt_gen_0_750_lumiscaled");
  // --- NNLO QCD + EWK corrected GEN folders (genWeight * lumi * EWK * NNLOQCD) ---
  // Parallel to the _lumiscaled set above; only booked when ApplyNNLOEWKToGenHists is on.
  if(apply_nnloewk_to_gen_hists){
    histogram_tags.push_back("mtt_gen_inclusive_lumiscaled_nnloewk");
    for (const auto & tag : mttbar_bin_tags) {
      histogram_tags.push_back("mtt_gen_" + tag + "_lumiscaled_nnloewk");
    }
    histogram_tags.push_back("mtt_gen_0_750_lumiscaled_nnloewk");
  }

  // --- Reco-level folders ---
  histogram_tags.push_back("MET");

  book_histograms(ctx, histogram_tags);

  // POWHEG EWK lookup: book yield (HATHOR binning) + sum-of-genWeight histograms.
  // Off unless PowhegEWKLookupBinningFile is set in the XML, so normal runs are unchanged.
  const std::string powheg_lookup_file = ctx.get("PowhegEWKLookupBinningFile", "");
  if (isMC && !powheg_lookup_file.empty()) {
    const std::string powheg_lookup_hist = ctx.get("PowhegEWKLookupBinningHist", "EWno_0");
    TFile lookup_file(powheg_lookup_file.c_str(), "READ");
    if (lookup_file.IsZombie()) {
      throw std::runtime_error("ZprimePreselectionModule: failed to open PowhegEWKLookupBinningFile '" + powheg_lookup_file + "'");
    }
    TH2D* lookup_src = dynamic_cast<TH2D*>(lookup_file.Get(powheg_lookup_hist.c_str()));
    if (!lookup_src) {
      throw std::runtime_error("ZprimePreselectionModule: no TH2 '" + powheg_lookup_hist + "' in '" + powheg_lookup_file + "'");
    }
    h_powheg_lookup_yield = dynamic_cast<TH2D*>(lookup_src->Clone("powheg_yield"));
    h_powheg_lookup_yield->SetDirectory(0);
    h_powheg_lookup_yield->Reset("ICES");
    if (h_powheg_lookup_yield->GetSumw2N() == 0) h_powheg_lookup_yield->Sumw2();
    h_powheg_lookup_yield->SetTitle("POWHEG yield, raw genWeight;M_{t#bar{t}} [GeV];cos#theta*");
    h_powheg_lookup_sumw = new TH1D("powheg_sumw", "sum of genWeight, all events", 1, 0., 1.);
    h_powheg_lookup_sumw->Sumw2();
    ctx.put("powheg_ewk_lookup", h_powheg_lookup_yield);
    ctx.put("powheg_ewk_lookup", h_powheg_lookup_sumw);
    lookup_file.Close();
  }

  lumihists.reset(new LuminosityHists(ctx, "lumi"));
}

bool ZprimePreselectionModule::process(uhh2::Event& event){

  // Process TTbarGen first
  if (isMC && ttgenprod) {
    ttgenprod->process(event);
  }

  // b-fragmentation
  if (isMC) {
    BFragInputs bfrag_inputs;
    if (event.is_valid(h_ttbargen)) {
      bfrag_inputs = compute_bfrag_inputs(event.genparticles, event.get(h_ttbargen));
    }

    event.set(h_bfrag_xb_top, bfrag_inputs.xb_top);
    event.set(h_bfrag_xb_antitop, bfrag_inputs.xb_antitop);
    event.set(h_bfrag_rb, kBFragNominalInputRb);
    event.set(h_bfrag_valid, bfrag_inputs.valid);

    // Real b-fragmentation weights require the trained ONNX models. Keep the
    // branches neutral until model inference or official plugin output is wired.
    event.set(h_weight_bfrag_nom, 1.0f);
    event.set(h_weight_bfrag_up, 1.0f);

    // NNLO QCD + EWK gen-template correction. These modules declare output branches that SFrame
    // requires to be set on EVERY written event, so process() must run for all MC events here
    // (the modules set their weights to 1.0 for non-ttbar). event.weight is restored immediately;
    // the correction itself is applied only to the dedicated _nnloewk GEN folders below, by reading
    // the stored nominal weights from the handles.
    if (apply_nnloewk_to_gen_hists) {
      const double w_keep = event.weight;
      ttbar_ewk_gen->process(event);
      ttbar_nnloqcd_gen->process(event);
      event.weight = w_keep;
    }
  }

  if(debug) cout << "++++++++++++ NEW EVENT ++++++++++++++" << endl;
  if(debug) cout << " run.event: " << event.run << ". " << event.event << endl;
  if(debug) cout << " event.year: " << event.year << ". " << event.event << endl;

  if(!event.isRealData){
    if(!SignSplit->passes(event)) return false;
  }
  if(debug) cout << "beginning: ok" << endl;

  // ==========================================
  // At this point: event.weight = genInfo->weights()[0] (set by framework)
  // NO lumi scaling yet, NO PU, NO reco weights
  // ==========================================

  // --- CUTFLOW: 0_AllEvents ---
  h_cf_raw->Fill(0.5, 1.);
  h_cf_genweight->Fill(0.5, event.weight);
  h_cf_lumiscaled->Fill(0.5, event.weight * lumi_factor);

  // POWHEG EWK lookup: every MC event counts toward the cross-section denominator,
  // independent of decay channel or any later selection (event.weight is still the
  // raw genWeight here).
  if (h_powheg_lookup_sumw && !event.isRealData) {
    h_powheg_lookup_sumw->Fill(0.5, event.weight);
  }

  fill_histograms(event, "Input");
  if(debug) cout << "first plots input: ok" << endl;

  // ==========================================
  // FIRST: fill raw GEN mttbar folders (genWeight only — no lumi)
  // ==========================================
  if (isMC) {
    // Defaults are always set for MC, even if h_ttbargen is not valid
    event.set(h_xi_gen,       std::numeric_limits<float>::quiet_NaN());
    event.set(h_mtt_gen,      std::numeric_limits<float>::quiet_NaN());
    event.set(h_costheta_gen, std::numeric_limits<float>::quiet_NaN());
    event.set(h_DeltaY_gen,   std::numeric_limits<float>::quiet_NaN());

    bool gen_filled = false;
    if (event.is_valid(h_ttbargen)) {
      const auto& ttbargen = event.get(h_ttbargen);
      const auto dc = ttbargen.DecayChannel();

      // POWHEG EWK lookup: fill before any selection, for ALL decay channels
      // (deliberately no hadronic-tau veto — sigma_POWHEG is decay-inclusive).
      // Out-of-range values are clamped into the first/last bin, mirroring how
      // TTbarEWKCorrection::evaluate_histo clamps at evaluation time.
      if (h_powheg_lookup_yield && dc != TTbarGen::e_notfound) {
        const double mtt_lookup = (ttbargen.Top().v4() + ttbargen.Antitop().v4()).M();
        const double ct_lookup = compute_top_costheta_ttbarframe(ttbargen);
        if (std::isfinite(mtt_lookup) && std::isfinite(ct_lookup)) {
          const TAxis* ax = h_powheg_lookup_yield->GetXaxis();
          const TAxis* ay = h_powheg_lookup_yield->GetYaxis();
          h_powheg_lookup_yield->Fill(
            std::min(std::max(mtt_lookup, ax->GetBinCenter(1)), ax->GetBinCenter(ax->GetNbins())),
            std::min(std::max(ct_lookup, ay->GetBinCenter(1)), ay->GetBinCenter(ay->GetNbins())),
            event.weight);
        }
      }
      const bool fill_gen = (dc != TTbarGen::e_notfound && !has_any_hadronic_tau_decay(event.genparticles, ttbargen));
      if (fill_gen) {
        const auto& top  = ttbargen.Top();
        const auto& atop = ttbargen.Antitop();
        double mtt = (top.v4() + atop.v4()).M();
        double dy  = std::abs(top.v4().Rapidity()) - std::abs(atop.v4().Rapidity());
        event.set(h_xi_gen,     std::tanh(dy));
        event.set(h_mtt_gen,    static_cast<float>(mtt));
        event.set(h_costheta_gen, compute_top_costheta_ttbarframe(ttbargen));
        event.set(h_DeltaY_gen, static_cast<float>(dy));

        // --- CUTFLOW: 1_GEN_raw ---
        h_cf_raw->Fill(1.5, 1.);
        h_cf_genweight->Fill(1.5, event.weight);

        // Fill RAW (no lumi) GEN folders
        fill_histograms(event, "mtt_gen_inclusive");
        const int ibin = find_mtt_bin(mtt);
        if (ibin >= 0) {
          fill_histograms(event, "mtt_gen_" + mttbar_bin_tags.at(ibin));
        }
        // Merged [0, 750) GeV: same events as fine bins 0-500 and 500-750
        if (mtt >= mttbar_bin_edges[0] && mtt < mttbar_bin_edges[2]) {
          fill_histograms(event, "mtt_gen_0_750");
        }

        // ==========================================
        // NOW: apply lumi scaling (same events as raw GEN fills above)
        // ==========================================
        lumi_weight->process(event);
        // event.weight = genWeight * (target_lumi / dataset_lumi)

        // --- CUTFLOW: 2_AfterLumiWeight, 3_GEN_lumiscaled ---
        h_cf_raw->Fill(2.5, 1.);
        h_cf_lumiscaled->Fill(2.5, event.weight);
        h_cf_raw->Fill(3.5, 1.);
        h_cf_lumiscaled->Fill(3.5, event.weight);

        // Fill LUMI-SCALED GEN folders (same events, now with lumi weight)
        fill_histograms(event, "mtt_gen_inclusive_lumiscaled");
        if (ibin >= 0) {
          fill_histograms(event, "mtt_gen_" + mttbar_bin_tags.at(ibin) + "_lumiscaled");
        }
        if (mtt >= mttbar_bin_edges[0] && mtt < mttbar_bin_edges[2]) {
          fill_histograms(event, "mtt_gen_0_750_lumiscaled");
        }

        // ==========================================
        // OPTIONAL: NNLO QCD + EWK corrected GEN folders (applied AFTER lumi scaling).
        // Only used to build NoAC gen templates; event.weight is restored afterwards so the
        // downstream reco-level processing (CommonModules onward) sees the lumi-only weight.
        // ==========================================
        if (apply_nnloewk_to_gen_hists) {
          const double w_lumi = event.weight;          // genWeight * lumi
          const float w_ewk  = event.get(h_w_ewk_nom);   // EWK ratio (1.0 if not ttbar)
          const float w_nnlo = event.get(h_w_nnlo_nom);  // sqrt(SF(pt_top)*SF(pt_atop))
          event.weight = w_lumi * w_ewk * w_nnlo;
          fill_histograms(event, "mtt_gen_inclusive_lumiscaled_nnloewk");
          if (ibin >= 0) {
            fill_histograms(event, "mtt_gen_" + mttbar_bin_tags.at(ibin) + "_lumiscaled_nnloewk");
          }
          if (mtt >= mttbar_bin_edges[0] && mtt < mttbar_bin_edges[2]) {
            fill_histograms(event, "mtt_gen_0_750_lumiscaled_nnloewk");
          }
          event.weight = w_lumi;                        // restore lumi-only weight
        }
        gen_filled = true;
      }
    }

    if (!gen_filled) {
      // Non-semileptonic, hadronic-tau, or no ttbargen: still apply lumi weight for reco-level
      lumi_weight->process(event);
    }
  }

  // ==========================================
  // CommonModules: PU reweight + JEC/JER + cleaners
  // (lumi already applied above, disabled in CommonModules)
  // ==========================================
  bool commonResult = common->process(event);
  if (!commonResult) return false;

  // --- CUTFLOW: 4_AfterCommonModules (lumi + PU) ---
  h_cf_raw->Fill(4.5, 1.);
  h_cf_lumiscaled->Fill(4.5, event.weight);

  if(debug) cout << "CommonModules: ok" << endl;
  // fill_histograms(event, "CommonModules");

  sort_by_pt<Muon>(*event.muons);
  sort_by_pt<Electron>(*event.electrons);

  // Correct AK4 CHS jets
  CHSjetCorr->process(event);

  if(isHOTVR){
    hotvrjetCorr->process(event);
    // fill_histograms(event, "HOTVRCorrections");
  }

  toppuppijetCorr->process(event);
  if(debug) cout << "TopPuppiJetCorrections: ok" << endl;
  // fill_histograms(event, "PUPPICorrections");

  // GEN ME quark-flavor selection
  if(!event.isRealData){
    if(!genflavor_sel->passes(event)) return false;
  }
  if(debug) cout << "GenFlavorSelection: ok" << endl;

  const bool pass_lep1 = ((event.muons->size() >= 1) || (event.electrons->size() >= 1));
  if(!pass_lep1) return false;
  h_cf_raw->Fill(5.5, 1.); h_cf_lumiscaled->Fill(5.5, event.weight);
  if(debug) cout << "≥1 leptons: ok" << endl;
  // fill_histograms(event, "Lepton1");

  jet_IDcleaner->process(event);
  // fill_histograms(event, "JetID");
  if(debug) cout << "JetCleaner ID: ok" << endl;

  jet_cleaner1->process(event);
  sort_by_pt<Jet>(*event.jets);
  // fill_histograms(event, "JetCleaner1");
  if(debug) cout << "JetCleaner1: ok" << endl;

  // Lepton-2Dcut variables
  for(auto& muo : *event.muons){
    float    dRmin, pTrel;
    std::tie(dRmin, pTrel) = drmin_pTrel(muo, *event.jets);
    muo.set_tag(Muon::twodcut_dRmin, dRmin);
    muo.set_tag(Muon::twodcut_pTrel, pTrel);
  }

  for(auto& ele : *event.electrons){
    float    dRmin, pTrel;
    std::tie(dRmin, pTrel) = drmin_pTrel(ele, *event.jets);
    ele.set_tag(Electron::twodcut_dRmin, dRmin);
    ele.set_tag(Electron::twodcut_pTrel, pTrel);
  }

  jet_cleaner2->process(event);
  sort_by_pt<Jet>(*event.jets);
  // fill_histograms(event, "JetCleaner2");
  if(debug) cout << "JetCleaner2: ok" << endl;

  hotvrjet_cleaner->process(event);
  sort_by_pt<TopJet>(*event.topjets);

  topjet_puppi_IDcleaner->process(event);
  topjet_puppi_cleaner->process(event);
  sort_by_pt<TopJet>(*event.toppuppijets);

  // fill_histograms(event, "TopjetCleaner");
  if(debug) cout << "TopJetCleaner: ok" << endl;

  // 1st AK4 jet selection
  const bool pass_jet1 = jet1_sel->passes(event);
  if(!pass_jet1) return false;
  h_cf_raw->Fill(6.5, 1.); h_cf_lumiscaled->Fill(6.5, event.weight);

  if(debug) cout << "NJetSelection1: ok" << endl;
  // fill_histograms(event, "Jet1");

  // 2nd AK4 jet selection
  const bool pass_jet2 = jet2_sel->passes(event);
  if(!pass_jet2) return false;
  if(debug) cout << "NJetSelection2: ok" << endl;
  // fill_histograms(event, "Jet2");

  // MET selection
  const bool pass_met = met_sel->passes(event);
  if(!pass_met) return false;
  h_cf_raw->Fill(7.5, 1.); h_cf_lumiscaled->Fill(7.5, event.weight);
  if(debug) cout << "METCut: ok" << endl;
  fill_histograms(event, "MET");

  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(ZprimePreselectionModule)