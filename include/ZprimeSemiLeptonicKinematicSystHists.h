#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicModules.h"

#include <TH1F.h>
#include <string>
#include <vector>
#include <map>

/*
 * ZprimeSemiLeptonicKinematicSystHists
 * ------------------------------------
 * Data/MC kinematic control variables stored together with ALL weight-based
 * systematic variations, in-file. For each variable it books:
 *    <name>               (nominal)
 *    <name>_<suffix>      (one per weight-systematic variation, see syst_suffixes())
 *
 * Object/shape systematics (JEC, JER, hdamp, topmass, ...) are NOT produced here;
 * those come from the dedicated systematic workdirs and are already picked up by the
 * plotter. This class only adds the weight-based systematics (pileup, b-tag shape,
 * lepton SFs, mu_R/mu_F scale, ISR/FSR, top-pt, EWK, NNLO-QCD, ...) so the kinematic
 * data/MC band can include them too.
 *
 * The variable definitions, binnings and the per-systematic weight rules mirror exactly
 * what ZprimeSemiLeptonicHists (nominal) and ZprimeSemiLeptonicSystematicsHists (DeltaY)
 * already do, so the in-file variations line up bin-for-bin with the nominal histograms.
 *
 * Fill this for MC only, at the same selection point as the nominal SR kinematic plots.
 */
class ZprimeSemiLeptonicKinematicSystHists: public uhh2::Hists {
public:
  explicit ZprimeSemiLeptonicKinematicSystHists(uhh2::Context & ctx, const std::string & dirname);
  virtual void fill(const uhh2::Event & event) override;
  virtual ~ZprimeSemiLeptonicKinematicSystHists();

protected:
  // Canonical ordered list of variation suffixes. book_var() and fill() build their
  // histogram lists / weight lists in THIS order, so they stay in lock-step.
  static const std::vector<std::string> & syst_suffixes();

  // One kinematic variable: nominal hist + variation hists (parallel to syst_suffixes()).
  struct VarHists {
    TH1F* nominal = nullptr;
    std::vector<TH1F*> variations;
  };
  void book_var(const std::string & name, const std::string & title, int nbins, double lo, double hi);
  // Fill nominal with w_nom and each variation i with w_var.at(i). w_var must have
  // syst_suffixes().size() entries (asserted at runtime).
  void fill_var(const std::string & name, double value, double w_nom, const std::vector<double> & w_var);

  std::map<std::string, VarHists> vars_;

  bool is_mc;
  bool is_tt;   // true for ttbar datasets; gates the ttbar-only systematics (scale, ISR/FSR)
  bool isMuon, isElectron;
  bool debug;

  // object handles
  uhh2::Event::Handle<bool> h_is_zprime_reconstructed_chi2;
  uhh2::Event::Handle<ZprimeCandidate*> h_BestZprimeCandidateChi2;
  // DeepAK8 / HOTVR top-tagged AK8 PUPPI jets (same handle the nominal Hists uses)
  uhh2::Event::Handle<std::vector<TopJet>> h_AK8TopTags;

  // weight handles (names copied verbatim from ZprimeSemiLeptonicSystematicsHists)
  uhh2::Event::Handle<float> h_ele_reco, h_ele_reco_up, h_ele_reco_down;
  uhh2::Event::Handle<float> h_ele_id, h_ele_id_up, h_ele_id_down;
  uhh2::Event::Handle<float> h_ele_trigger, h_ele_trigger_up, h_ele_trigger_down;
  uhh2::Event::Handle<float> h_mu_reco, h_mu_reco_up, h_mu_reco_down;
  uhh2::Event::Handle<float> h_mu_iso_stat, h_mu_iso_stat_up, h_mu_iso_stat_down;
  uhh2::Event::Handle<float> h_mu_iso_syst, h_mu_iso_syst_up, h_mu_iso_syst_down;
  uhh2::Event::Handle<float> h_mu_id_stat, h_mu_id_stat_up, h_mu_id_stat_down;
  uhh2::Event::Handle<float> h_mu_id_syst, h_mu_id_syst_up, h_mu_id_syst_down;
  uhh2::Event::Handle<float> h_mu_trigger_stat, h_mu_trigger_stat_up, h_mu_trigger_stat_down;
  uhh2::Event::Handle<float> h_mu_trigger_syst, h_mu_trigger_syst_up, h_mu_trigger_syst_down;
  uhh2::Event::Handle<float> h_pu, h_pu_up, h_pu_down;
  uhh2::Event::Handle<float> h_prefiring, h_prefiring_up, h_prefiring_down;
  uhh2::Event::Handle<float> h_murmuf_upup, h_murmuf_upnone, h_murmuf_noneup;
  uhh2::Event::Handle<float> h_murmuf_nonedown, h_murmuf_downnone, h_murmuf_downdown;
  uhh2::Event::Handle<float> h_isr_up, h_isr_down, h_fsr_up, h_fsr_down;
  uhh2::Event::Handle<float> h_btag;
  uhh2::Event::Handle<float> h_btag_cferr1_up, h_btag_cferr1_down, h_btag_cferr2_up, h_btag_cferr2_down;
  uhh2::Event::Handle<float> h_btag_hf_up, h_btag_hf_down, h_btag_hfstats1_up, h_btag_hfstats1_down;
  uhh2::Event::Handle<float> h_btag_hfstats2_up, h_btag_hfstats2_down, h_btag_lf_up, h_btag_lf_down;
  uhh2::Event::Handle<float> h_btag_lfstats1_up, h_btag_lfstats1_down, h_btag_lfstats2_up, h_btag_lfstats2_down;
  uhh2::Event::Handle<float> h_ttag, h_ttag_corr_up, h_ttag_corr_down, h_ttag_uncorr_up, h_ttag_uncorr_down;
  uhh2::Event::Handle<float> h_tmistag, h_tmistag_up, h_tmistag_down;
  uhh2::Event::Handle<float> h_ttbar_ewk_nominal, h_ttbar_ewk_up, h_ttbar_ewk_down;
  uhh2::Event::Handle<float> h_ttbar_nnlo_qcd_nominal, h_ttbar_nnlo_qcd_up, h_ttbar_nnlo_qcd_down;
};
