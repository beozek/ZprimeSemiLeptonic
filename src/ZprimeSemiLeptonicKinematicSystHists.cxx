#include "UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicKinematicSystHists.h"

#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/LorentzVector.h"

#include <cmath>
#include <cassert>
#include <stdexcept>

using namespace std;
using namespace uhh2;

namespace {
  // weight * (var/nom), with division-by-zero / non-finite protection.
  // Same behaviour as safe_syst_ratio() in ZprimeSemiLeptonicSystematicsHists.cxx.
  inline double safe_syst_ratio(double weight, double var, double nom){
    if(!std::isfinite(nom) || nom == 0.0 || !std::isfinite(var)) return weight;
    const double ratio = var / nom;
    if(!std::isfinite(ratio)) return weight;
    return weight * ratio;
  }
}

const std::vector<std::string> & ZprimeSemiLeptonicKinematicSystHists::syst_suffixes(){
  // Order here MUST match the order in which fill() builds the per-event weight vector.
  static const std::vector<std::string> s = {
    // --- simple up/down reweighting (lepton SFs, pileup, prefiring) ---
    "ele_reco_up", "ele_reco_down",
    "ele_id_up", "ele_id_down",
    "ele_trigger_up", "ele_trigger_down",
    "mu_reco_up", "mu_reco_down",
    "mu_iso_stat_up", "mu_iso_stat_down",
    "mu_iso_syst_up", "mu_iso_syst_down",
    "mu_id_stat_up", "mu_id_stat_down",
    "mu_id_syst_up", "mu_id_syst_down",
    "mu_trigger_stat_up", "mu_trigger_stat_down",
    "mu_trigger_syst_up", "mu_trigger_syst_down",
    "pu_up", "pu_down",
    "prefiring_up", "prefiring_down",
    // --- mu_R / mu_F scale (multiplicative weight) ---
    "murmuf_upup", "murmuf_upnone", "murmuf_noneup",
    "murmuf_nonedown", "murmuf_downnone", "murmuf_downdown",
    // --- b-tag shape ---
    "btag_cferr1_up", "btag_cferr1_down",
    "btag_cferr2_up", "btag_cferr2_down",
    "btag_hf_up", "btag_hf_down",
    "btag_hfstats1_up", "btag_hfstats1_down",
    "btag_hfstats2_up", "btag_hfstats2_down",
    "btag_lf_up", "btag_lf_down",
    "btag_lfstats1_up", "btag_lfstats1_down",
    "btag_lfstats2_up", "btag_lfstats2_down",
    // --- top tagging / mistag ---
    "ttag_corr_up", "ttag_corr_down",
    "ttag_uncorr_up", "ttag_uncorr_down",
    "tmistag_up", "tmistag_down",
    // --- ttbar EWK / NNLO-QCD corrections ---
    "ewk_up", "ewk_down",
    "qcd_up", "qcd_down",
    // --- parton shower (ISR/FSR, multiplicative weight) ---
    "isr_up", "isr_down",
    "fsr_up", "fsr_down",
  };
  static const std::vector<std::string> s_full = [](){
    std::vector<std::string> v = s;
    // pdf: 100 NNPDF replica variations (ttbar only), appended after the explicit list.
    for(int i = 1; i <= 100; ++i) v.push_back("pdf_" + std::to_string(i));
    return v;
  }();
  return s_full;
}

ZprimeSemiLeptonicKinematicSystHists::ZprimeSemiLeptonicKinematicSystHists(Context & ctx, const string & dirname):
  Hists(ctx, dirname){

  is_mc      = ctx.get("dataset_type") == "MC";
  isMuon     = ctx.get("channel") == "muon";
  isElectron = ctx.get("channel") == "electron";
  debug      = false;

  // ttbar flag (same definition as ZprimeSemiLeptonicSystematicsHists), used to gate
  // the ttbar-only systematics (mu_R/mu_F scale, ISR/FSR) so non-ttbar gets no deviation.
  const std::string dataset_version = ctx.get("dataset_version");
  is_tt = (dataset_version.find("TTTo") == 0)
       || (dataset_version.find("EFT") != std::string::npos)
       || (dataset_version.find("TTJets") != std::string::npos);

  h_is_zprime_reconstructed_chi2 = ctx.get_handle<bool>("is_zprime_reconstructed_chi2");
  h_BestZprimeCandidateChi2      = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestChi2");

  // top-tagged AK8 PUPPI jets: pick the same collection the nominal Hists uses
  // (DeepAK8 or HOTVR), so the AK8 control plots show the actual hadronic-top proxy.
  const bool ishotvr   = (ctx.get("is_hotvr")   == "true");
  const bool isdeepAK8 = (ctx.get("is_deepAK8") == "true");
  if(isdeepAK8)    h_AK8TopTags = ctx.get_handle<std::vector<TopJet>>("DeepAK8TopTags");
  else if(ishotvr) h_AK8TopTags = ctx.get_handle<std::vector<TopJet>>("HOTVRTopTags");

  // weight handles (verbatim names from ZprimeSemiLeptonicSystematicsHists::init)
  h_ele_reco           = ctx.get_handle<float>("weight_sfelec_reco");
  h_ele_reco_up        = ctx.get_handle<float>("weight_sfelec_reco_up");
  h_ele_reco_down      = ctx.get_handle<float>("weight_sfelec_reco_down");
  h_ele_id             = ctx.get_handle<float>("weight_sfelec_id");
  h_ele_id_up          = ctx.get_handle<float>("weight_sfelec_id_up");
  h_ele_id_down        = ctx.get_handle<float>("weight_sfelec_id_down");
  h_ele_trigger        = ctx.get_handle<float>("weight_sfelec_trigger");
  h_ele_trigger_up     = ctx.get_handle<float>("weight_sfelec_trigger_up");
  h_ele_trigger_down   = ctx.get_handle<float>("weight_sfelec_trigger_down");
  h_mu_reco            = ctx.get_handle<float>("weight_sfmu_reco");
  h_mu_reco_up         = ctx.get_handle<float>("weight_sfmu_reco_up");
  h_mu_reco_down       = ctx.get_handle<float>("weight_sfmu_reco_down");
  h_mu_iso_stat        = ctx.get_handle<float>("weight_sfmu_iso_stat");
  h_mu_iso_stat_up     = ctx.get_handle<float>("weight_sfmu_iso_stat_up");
  h_mu_iso_stat_down   = ctx.get_handle<float>("weight_sfmu_iso_stat_down");
  h_mu_iso_syst        = ctx.get_handle<float>("weight_sfmu_iso_syst");
  h_mu_iso_syst_up     = ctx.get_handle<float>("weight_sfmu_iso_syst_up");
  h_mu_iso_syst_down   = ctx.get_handle<float>("weight_sfmu_iso_syst_down");
  h_mu_id_stat         = ctx.get_handle<float>("weight_sfmu_id_stat");
  h_mu_id_stat_up      = ctx.get_handle<float>("weight_sfmu_id_stat_up");
  h_mu_id_stat_down    = ctx.get_handle<float>("weight_sfmu_id_stat_down");
  h_mu_id_syst         = ctx.get_handle<float>("weight_sfmu_id_syst");
  h_mu_id_syst_up      = ctx.get_handle<float>("weight_sfmu_id_syst_up");
  h_mu_id_syst_down    = ctx.get_handle<float>("weight_sfmu_id_syst_down");
  h_mu_trigger_stat    = ctx.get_handle<float>("weight_sfmu_trigger_stat");
  h_mu_trigger_stat_up = ctx.get_handle<float>("weight_sfmu_trigger_stat_up");
  h_mu_trigger_stat_down = ctx.get_handle<float>("weight_sfmu_trigger_stat_down");
  h_mu_trigger_syst    = ctx.get_handle<float>("weight_sfmu_trigger_syst");
  h_mu_trigger_syst_up = ctx.get_handle<float>("weight_sfmu_trigger_syst_up");
  h_mu_trigger_syst_down = ctx.get_handle<float>("weight_sfmu_trigger_syst_down");
  h_pu                 = ctx.get_handle<float>("weight_pu");
  h_pu_up              = ctx.get_handle<float>("weight_pu_up");
  h_pu_down            = ctx.get_handle<float>("weight_pu_down");
  h_prefiring          = ctx.get_handle<float>("prefiringWeight");
  h_prefiring_up       = ctx.get_handle<float>("prefiringWeightUp");
  h_prefiring_down     = ctx.get_handle<float>("prefiringWeightDown");
  h_murmuf_upup        = ctx.get_handle<float>("weight_murmuf_upup");
  h_murmuf_upnone      = ctx.get_handle<float>("weight_murmuf_upnone");
  h_murmuf_noneup      = ctx.get_handle<float>("weight_murmuf_noneup");
  h_murmuf_nonedown    = ctx.get_handle<float>("weight_murmuf_nonedown");
  h_murmuf_downnone    = ctx.get_handle<float>("weight_murmuf_downnone");
  h_murmuf_downdown    = ctx.get_handle<float>("weight_murmuf_downdown");
  h_isr_up             = ctx.get_handle<float>("weight_isr_2_up");
  h_isr_down           = ctx.get_handle<float>("weight_isr_2_down");
  h_fsr_up             = ctx.get_handle<float>("weight_fsr_2_up");
  h_fsr_down           = ctx.get_handle<float>("weight_fsr_2_down");
  h_btag               = ctx.get_handle<float>("weight_btagdisc_central");
  h_btag_cferr1_up     = ctx.get_handle<float>("weight_btagdisc_cferr1_up");
  h_btag_cferr1_down   = ctx.get_handle<float>("weight_btagdisc_cferr1_down");
  h_btag_cferr2_up     = ctx.get_handle<float>("weight_btagdisc_cferr2_up");
  h_btag_cferr2_down   = ctx.get_handle<float>("weight_btagdisc_cferr2_down");
  h_btag_hf_up         = ctx.get_handle<float>("weight_btagdisc_hf_up");
  h_btag_hf_down       = ctx.get_handle<float>("weight_btagdisc_hf_down");
  h_btag_hfstats1_up   = ctx.get_handle<float>("weight_btagdisc_hfstats1_up");
  h_btag_hfstats1_down = ctx.get_handle<float>("weight_btagdisc_hfstats1_down");
  h_btag_hfstats2_up   = ctx.get_handle<float>("weight_btagdisc_hfstats2_up");
  h_btag_hfstats2_down = ctx.get_handle<float>("weight_btagdisc_hfstats2_down");
  h_btag_lf_up         = ctx.get_handle<float>("weight_btagdisc_lf_up");
  h_btag_lf_down       = ctx.get_handle<float>("weight_btagdisc_lf_down");
  h_btag_lfstats1_up   = ctx.get_handle<float>("weight_btagdisc_lfstats1_up");
  h_btag_lfstats1_down = ctx.get_handle<float>("weight_btagdisc_lfstats1_down");
  h_btag_lfstats2_up   = ctx.get_handle<float>("weight_btagdisc_lfstats2_up");
  h_btag_lfstats2_down = ctx.get_handle<float>("weight_btagdisc_lfstats2_down");
  h_ttag               = ctx.get_handle<float>("weight_toptagsf");
  h_ttag_corr_up       = ctx.get_handle<float>("weight_toptagsf_corr_up");
  h_ttag_corr_down     = ctx.get_handle<float>("weight_toptagsf_corr_down");
  h_ttag_uncorr_up     = ctx.get_handle<float>("weight_toptagsf_uncorr_up");
  h_ttag_uncorr_down   = ctx.get_handle<float>("weight_toptagsf_uncorr_down");
  h_tmistag            = ctx.get_handle<float>("weight_topmistagsf");
  h_tmistag_up         = ctx.get_handle<float>("weight_topmistagsf_up");
  h_tmistag_down       = ctx.get_handle<float>("weight_topmistagsf_down");
  h_ttbar_ewk_nominal  = ctx.get_handle<float>("weight_ttbar_ewk_nominal");
  h_ttbar_ewk_up       = ctx.get_handle<float>("weight_ttbar_ewk_up");
  h_ttbar_ewk_down     = ctx.get_handle<float>("weight_ttbar_ewk_down");
  h_ttbar_nnlo_qcd_nominal = ctx.get_handle<float>("weight_ttbar_nnlo_qcd");
  h_ttbar_nnlo_qcd_up      = ctx.get_handle<float>("weight_ttbar_nnlo_qcd_up");
  h_ttbar_nnlo_qcd_down    = ctx.get_handle<float>("weight_ttbar_nnlo_qcd_down");

  // --- book the kinematic variables. Names and binnings mirror the nominal
  //     ZprimeSemiLeptonicHists so the in-file variations line up bin-for-bin. ---
  // reconstructed tops (best chi2 candidate)
  book_var("toplep_pt",        "p_{T}^{t,lep} [GeV]",       70, 0, 7000);
  book_var("tophad_pt",        "p_{T}^{t,had} [GeV]",       70, 0, 7000);
  book_var("toplep_eta",       "#eta^{t,lep}",              60, -3.0, 3.0);
  book_var("tophad_eta",       "#eta^{t,had}",              60, -3.0, 3.0);
  // MET (name/binning match nominal MET_rebin and NN_MET_phi)
  book_var("MET_rebin",        "missing E_{T} [GeV]",       45, 0, 900);
  book_var("NN_MET_phi",       "#phi^{miss}",               35, -3.5, 3.5);
  // leading DeepAK8 top-tagged AK8 PUPPI jet (hadronic-top proxy)
  book_var("pt_AK8PuppiTaggedjet1",  "p_{T}^{AK8Puppi Tagged jet 1} [GeV]", 45, 0, 900);
  book_var("eta_AK8PuppiTaggedjet1", "#eta^{AK8Puppi Tagged jet 1}",     50, -2.5, 2.5);
  book_var("mSD_AK8PuppiTaggedjet1", "m_{SD}^{AK8Puppi Tagged jet 1}",   50, 0, 500);
  // leading lepton (booked per channel)
  if(isMuon){
    book_var("pt_mu1",         "p_{T}^{#mu 1} [GeV]",       90, 0, 900);
    book_var("eta_mu1",        "#eta^{#mu 1}",              50, -2.5, 2.5);
  }
  if(isElectron){
    book_var("pt_ele1",        "p_{T}^{e 1} [GeV]",         90, 0, 900);
    book_var("eta_ele1",       "#eta^{e 1}",                50, -2.5, 2.5);
  }
  // leading AK4 jet
  book_var("pt_jet1",          "p_{T}^{jet 1} [GeV]",       45, 0, 900);
  book_var("eta_jet1",         "#eta^{jet 1}",              50, -2.5, 2.5);
}

void ZprimeSemiLeptonicKinematicSystHists::book_var(const string & name, const string & title,
                                                    int nbins, double lo, double hi){
  VarHists vh;
  vh.nominal = book<TH1F>(name.c_str(), title.c_str(), nbins, lo, hi);
  vh.nominal->Sumw2();
  const auto & suffixes = syst_suffixes();
  vh.variations.reserve(suffixes.size());
  for(const auto & suf : suffixes){
    const string hname = name + "_" + suf;
    TH1F* h = book<TH1F>(hname.c_str(), title.c_str(), nbins, lo, hi);
    h->Sumw2();
    vh.variations.push_back(h);
  }
  vars_[name] = std::move(vh);
}

void ZprimeSemiLeptonicKinematicSystHists::fill_var(const string & name, double value,
                                                    double w_nom, const vector<double> & w_var){
  auto it = vars_.find(name);
  if(it == vars_.end()) return;  // variable not booked in this channel
  VarHists & vh = it->second;
  vh.nominal->Fill(value, w_nom);
  for(size_t i = 0; i < vh.variations.size(); ++i){
    vh.variations[i]->Fill(value, w_var[i]);
  }
}

void ZprimeSemiLeptonicKinematicSystHists::fill(const Event & event){
  // Systematics only apply to MC; the in-file nominal here is the (MC-only) shift baseline.
  if(!is_mc) return;

  const double weight = event.weight;

  // --- read all systematic weights ---
  const float ele_reco_nom = event.get(h_ele_reco);
  const float ele_id_nom   = event.get(h_ele_id);
  const float ele_trig_nom = event.get(h_ele_trigger);
  const float mu_reco_nom  = event.get(h_mu_reco);
  const float mu_iso_stat_nom = event.get(h_mu_iso_stat);
  const float mu_iso_syst_nom = event.get(h_mu_iso_syst);
  const float mu_id_stat_nom  = event.get(h_mu_id_stat);
  const float mu_id_syst_nom  = event.get(h_mu_id_syst);
  const float mu_trig_stat_nom = event.get(h_mu_trigger_stat);
  const float mu_trig_syst_nom = event.get(h_mu_trigger_syst);
  const float pu_nom       = event.get(h_pu);
  const float prefiring_nom= event.get(h_prefiring);
  const float btag_nom     = event.get(h_btag);
  const float ttag_nom     = event.get(h_ttag);
  const float tmistag_nom  = event.get(h_tmistag);
  const float ewk_nom      = is_tt ? event.get(h_ttbar_ewk_nominal)      : 1.0f;
  const float qcd_nom      = is_tt ? event.get(h_ttbar_nnlo_qcd_nominal) : 1.0f;

  // --- build the per-event weight vector, in the SAME order as syst_suffixes() ---
  vector<double> w;
  w.reserve(syst_suffixes().size());
  // simple up/down reweighting
  w.push_back(safe_syst_ratio(weight, event.get(h_ele_reco_up),   ele_reco_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_ele_reco_down), ele_reco_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_ele_id_up),     ele_id_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_ele_id_down),   ele_id_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_ele_trigger_up),   ele_trig_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_ele_trigger_down), ele_trig_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_mu_reco_up),   mu_reco_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_mu_reco_down), mu_reco_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_mu_iso_stat_up),   mu_iso_stat_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_mu_iso_stat_down), mu_iso_stat_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_mu_iso_syst_up),   mu_iso_syst_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_mu_iso_syst_down), mu_iso_syst_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_mu_id_stat_up),    mu_id_stat_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_mu_id_stat_down),  mu_id_stat_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_mu_id_syst_up),    mu_id_syst_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_mu_id_syst_down),  mu_id_syst_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_mu_trigger_stat_up),   mu_trig_stat_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_mu_trigger_stat_down), mu_trig_stat_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_mu_trigger_syst_up),   mu_trig_syst_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_mu_trigger_syst_down), mu_trig_syst_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_pu_up),   pu_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_pu_down), pu_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_prefiring_up),   prefiring_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_prefiring_down), prefiring_nom));
  // mu_R / mu_F scale (ttbar only): non-ttbar -> nominal weight, no deviation.
  // The stored branch is already the multiplicative factor.
  w.push_back(is_tt ? weight * event.get(h_murmuf_upup)     : weight);
  w.push_back(is_tt ? weight * event.get(h_murmuf_upnone)   : weight);
  w.push_back(is_tt ? weight * event.get(h_murmuf_noneup)   : weight);
  w.push_back(is_tt ? weight * event.get(h_murmuf_nonedown) : weight);
  w.push_back(is_tt ? weight * event.get(h_murmuf_downnone) : weight);
  w.push_back(is_tt ? weight * event.get(h_murmuf_downdown) : weight);
  // b-tag shape
  w.push_back(safe_syst_ratio(weight, event.get(h_btag_cferr1_up),   btag_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_btag_cferr1_down), btag_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_btag_cferr2_up),   btag_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_btag_cferr2_down), btag_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_btag_hf_up),   btag_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_btag_hf_down), btag_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_btag_hfstats1_up),   btag_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_btag_hfstats1_down), btag_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_btag_hfstats2_up),   btag_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_btag_hfstats2_down), btag_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_btag_lf_up),   btag_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_btag_lf_down), btag_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_btag_lfstats1_up),   btag_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_btag_lfstats1_down), btag_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_btag_lfstats2_up),   btag_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_btag_lfstats2_down), btag_nom));
  // top tagging / mistag
  w.push_back(safe_syst_ratio(weight, event.get(h_ttag_corr_up),   ttag_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_ttag_corr_down), ttag_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_ttag_uncorr_up),   ttag_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_ttag_uncorr_down), ttag_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_tmistag_up),   tmistag_nom));
  w.push_back(safe_syst_ratio(weight, event.get(h_tmistag_down), tmistag_nom));
  // ttbar EWK / NNLO-QCD (ttbar only): replace nominal correction by up/down branch.
  // Non-ttbar -> nominal weight (handles never read), so no deviation and crash-safe
  // even if a sample's EWK/QCD weights are not set.
  if(is_tt){
    w.push_back(safe_syst_ratio(weight, event.get(h_ttbar_ewk_up),        ewk_nom));
    w.push_back(safe_syst_ratio(weight, event.get(h_ttbar_ewk_down),      ewk_nom));
    w.push_back(safe_syst_ratio(weight, event.get(h_ttbar_nnlo_qcd_up),   qcd_nom));
    w.push_back(safe_syst_ratio(weight, event.get(h_ttbar_nnlo_qcd_down), qcd_nom));
  }
  else{
    w.push_back(weight);
    w.push_back(weight);
    w.push_back(weight);
    w.push_back(weight);
  }
  // parton shower (ttbar only): non-ttbar -> nominal weight, no deviation.
  w.push_back(is_tt ? weight * event.get(h_isr_up)   : weight);
  w.push_back(is_tt ? weight * event.get(h_isr_down) : weight);
  w.push_back(is_tt ? weight * event.get(h_fsr_up)   : weight);
  w.push_back(is_tt ? weight * event.get(h_fsr_down) : weight);
  // pdf (ttbar only): 100 NNPDF replica reweights, weight * systweight[i+9] / originalXWGTUP.
  // Non-ttbar, or samples lacking the 100 LHE pdf weights -> nominal weight (no deviation).
  // (Same construction as ZprimeSemiLeptonicPDFHists; MY_FIRST_INDEX = 9.)
  const int PDF_FIRST_INDEX = 9;
  const bool have_pdf = is_tt && event.genInfo
      && event.genInfo->systweights().size() >= (size_t)(PDF_FIRST_INDEX + 100);
  const double pdf_orig = have_pdf ? event.genInfo->originalXWGTUP() : 1.0;
  for(int i = 0; i < 100; ++i){
    if(have_pdf && std::isfinite(pdf_orig) && pdf_orig != 0.0)
      w.push_back(weight * (event.genInfo->systweights().at(i + PDF_FIRST_INDEX) / pdf_orig));
    else
      w.push_back(weight);
  }

  assert(w.size() == syst_suffixes().size());
  if(w.size() != syst_suffixes().size()){
    throw std::runtime_error("ZprimeSemiLeptonicKinematicSystHists: weight/suffix size mismatch");
  }

  // --- reconstructed tops (best chi2 candidate) ---
  if(event.get(h_is_zprime_reconstructed_chi2)){
    ZprimeCandidate* cand = event.get(h_BestZprimeCandidateChi2);
    if(cand){
      const LorentzVector toplep = cand->top_leptonic_v4();
      const LorentzVector tophad = cand->top_hadronic_v4();
      fill_var("toplep_pt",  toplep.Pt(),  weight, w);
      fill_var("tophad_pt",  tophad.Pt(),  weight, w);
      fill_var("toplep_eta", toplep.Eta(), weight, w);
      fill_var("tophad_eta", tophad.Eta(), weight, w);
    }
  }

  // --- MET ---
  if(event.met){
    fill_var("MET_rebin",  event.met->pt(),  weight, w);
    fill_var("NN_MET_phi", event.met->phi(), weight, w);
  }

  // --- leading top-tagged AK8 PUPPI jet (match nominal: AK8PuppiTopTags index 0) ---
  if(event.is_valid(h_AK8TopTags)){
    const std::vector<TopJet> & AK8PuppiTopTags = event.get(h_AK8TopTags);
    if(!AK8PuppiTopTags.empty()){
      const TopJet & j0 = AK8PuppiTopTags.at(0);
      fill_var("pt_AK8PuppiTaggedjet1",  j0.pt(),           weight, w);
      fill_var("eta_AK8PuppiTaggedjet1", j0.eta(),          weight, w);
      fill_var("mSD_AK8PuppiTaggedjet1", j0.softdropmass(), weight, w);
    }
  }

  // --- leading lepton ---
  if(isMuon && event.muons && !event.muons->empty()){
    fill_var("pt_mu1",  event.muons->at(0).pt(),  weight, w);
    fill_var("eta_mu1", event.muons->at(0).eta(), weight, w);
  }
  if(isElectron && event.electrons && !event.electrons->empty()){
    fill_var("pt_ele1",  event.electrons->at(0).pt(),  weight, w);
    fill_var("eta_ele1", event.electrons->at(0).eta(), weight, w);
  }

  // --- leading AK4 jet ---
  if(event.jets && !event.jets->empty()){
    fill_var("pt_jet1",  event.jets->at(0).pt(),  weight, w);
    fill_var("eta_jet1", event.jets->at(0).eta(), weight, w);
  }
}

ZprimeSemiLeptonicKinematicSystHists::~ZprimeSemiLeptonicKinematicSystHists(){}
