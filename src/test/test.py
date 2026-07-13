#!/usr/bin/env python2
"""
Make Data/MC plots for the kinematic variables used in the Z' semileptonic
analysis, with one uncertainty band containing MC statistics, luminosity, and
shape systematics from paired external up/down workdirs.

Default input:
  output_DNN/<year>/<channel>/workdir_AnalysisDNN_<year>_<channel>_templatemethod_newmassbins_<suffix>

Default systematics:
  sibling workdirs of the nominal directory whose names contain "_up" or "_down",
  for example workdir_AnalysisDNN_UL17_muon_JEC_up and ..._down.

Examples:
  python tanh_dataMC_withSyst.py --year UL17 --channel muon
  python tanh_dataMC_withSyst.py --year UL17 --channel both --variables all
  python tanh_dataMC_withSyst.py --year UL17 --channel muon --variables met_pt,ak8_pt,mtt
  python tanh_dataMC_withSyst.py --year all  --channel muon     --variables all  # Run 2, muon
  python tanh_dataMC_withSyst.py --year all  --channel combined --variables all  # Run 2, muon+electron merged
  python tanh_dataMC_withSyst.py --year all  --channel all      --variables all  # Run 2, merged PLUS muon-only and electron-only

With --year all (or run2/combined) the four UL years are summed and each systematic is
combined across years per the datacard's correlation scheme (correlated -> linear,
per-year nuisances like JER / stats -> quadrature). Output goes to Run2_<channel>/.
Note: --channel both with --year all makes two plots (Run2_muon, Run2_electron); it does
not merge the two channels into a single histogram.

For systematics stored somewhere else, add one or more parents:
  python tanh_dataMC_withSyst.py --year UL17 --channel muon \
    --syst-parent ../JEC_Analysis

If the systematic files use a different histogram folder, set --syst-folder.
Only use that when the nominal and variation folders correspond to the same
selection step.
"""

from __future__ import print_function

import glob as pyglob
import math
import os
import re
import sys
from array import array
from collections import OrderedDict
from optparse import OptionParser

import ROOT as root
root.gROOT.SetBatch(True)
root.TH1.AddDirectory(False)

from constants import _YEARS
from plotter import NiceStackWithRatio, Process


PREFIX_MC = "uhh2.AnalysisModuleRunner.MC."
PREFIX_DATA = "uhh2.AnalysisModuleRunner.DATA."


SAMPLE_COLORS = OrderedDict([
    ("TTToSemiLeptonic", {"legend": "t#bar{t} semilep.", "color": root.kRed}),
    ("TTToOthers",      {"legend": "t#bar{t} others",   "color": root.kPink + 3}),
    ("Toponium",        {"legend": "Toponium",          "color": root.kViolet + 1}),
    ("WJets",           {"legend": "W+jets",            "color": root.kSpring - 3}),
    ("DY",              {"legend": "Drell-Yan",         "color": root.kOrange}),
    ("ST",              {"legend": "Single t",          "color": root.kBlue}),
    ("QCD",             {"legend": "QCD",               "color": root.kAzure + 10}),
    ("Diboson",         {"legend": "Diboson",           "color": root.kMagenta}),
])

TTBAR_SAMPLES = set(["TTToSemiLeptonic", "TTToOthers"])


class Variable(object):
    def __init__(
        self,
        key,
        title,
        unit,
        hist_names,
        valid_channels=None,
        x_min=None,
        x_max=None,
        projection=None,
        aliases=None,
        rebin=1,
    ):
        self.key = key
        self.title = title
        self.unit = unit
        self.hist_names = hist_names
        self.valid_channels = valid_channels
        self.x_min = x_min
        self.x_max = x_max
        self.projection = projection
        self.aliases = aliases or []
        self.rebin = rebin   # merge this many bins (TH1.Rebin); applied to ALL hists read

    def names_for_channel(self, channel):
        names = self.hist_names
        if isinstance(names, dict):
            return names.get(channel, [])
        return names

    def applies_to_channel(self, channel):
        return self.valid_channels is None or channel in self.valid_channels


VARIABLES = OrderedDict([
    # "Top eta" is plotted as both reconstructed tops in one distribution.
    ("top_eta", Variable(
        "top_eta", "#eta^{top}", None,
        [["toplep_eta"], ["tophad_eta"]],
        x_min=-3.0, x_max=3.0,
        aliases=["topeta", "top_eta_both"],
    )),
    # "Top pt" is plotted as both reconstructed tops in one distribution.
    ("top_pt", Variable(
        "top_pt", "p_{T}^{top}", "GeV",
        [["toplep_pt"], ["tophad_pt"]],
        x_min=0.0, x_max=1500.0,
        aliases=["toppt", "top_pt_both"],
    )),
    ("met_pt", Variable(
        "met_pt", "p_{T}^{miss}", "GeV",
        [["MET_rebin", "MET", "NN_MET_pt"]],
        x_min=0.0, x_max=900.0,
        aliases=["met", "MET_pt"],
    )),
    ("met_phi", Variable(
        "met_phi", "#phi^{miss}", None,
        [["MET_phi", "phi_MET", "NN_MET_phi"]],
        x_min=-3.5, x_max=3.5,
        aliases=["MET_phi"],
    )),
    ("ak8_pt", Variable(
        "ak8_pt", "p_{T}^{AK8 PUPPI top-tagged jet 1}", "GeV",
        [["pt_AK8PuppiTaggedjet1"]],
        x_min=0.0, x_max=1500.0,
        aliases=["AK8jet_pt", "ak8jet_pt", "ak8tagged_pt"],
    )),
    ("ak8_eta", Variable(
        "ak8_eta", "#eta^{AK8 PUPPI top-tagged jet 1}", None,
        [["eta_AK8PuppiTaggedjet1"]],
        x_min=-2.5, x_max=2.5,
        aliases=["AK8jet_eta", "ak8jet_eta", "ak8tagged_eta"],
    )),
    ("ak8_msd", Variable(
        "ak8_msd", "m_{SD}^{AK8 PUPPI top-tagged jet 1}", "GeV",
        [["mSD_AK8PuppiTaggedjet1"]],
        x_min=0.0, x_max=500.0,
        aliases=["ak8_mSD", "AK8jet_mSD", "msd", "ak8tagged_msd"],
    )),
    ("mtt", Variable(
        "mtt", "M_{t#bar{t}}^{reco}", "GeV",
        [["M_Zprime_rebin8", "M_Zprime_rebin7", "M_Zprime_rebin", "M_Zprime", "ditop_mass"]],
        x_min=0.0, x_max=3000.0,
        aliases=["Mtt", "m_tt", "M_Zprime"],
    )),
    # Use M_Zprime_rebin8 (400 bins over 0-10000 -> 25 GeV native bins). It is the same
    # source as the nominal "mtt" and is filled ONCE per event for every sample. Do NOT
    # use M_Zprime / M_Zprime_rebin: those are ALSO filled in a second ttbar-only block in
    # ZprimeSemiLeptonicHists, so ttbar is double-counted and MC overshoots data.
    ("mtt_25gev", Variable(
        "mtt_25gev", "M_{t#bar{t}}^{reco}", "GeV",
        [["M_Zprime_rebin8"]],
        x_min=0.0, x_max=3000.0,
        aliases=["mtt_25", "Mtt_25gev"],
        rebin=1,  # 25 GeV native bins
    )),
    ("mtt_50gev", Variable(
        "mtt_50gev", "M_{t#bar{t}}^{reco}", "GeV",
        [["M_Zprime_rebin8"]],
        x_min=0.0, x_max=3000.0,
        aliases=["mtt_50", "Mtt_50gev"],
        rebin=2,  # 25 GeV native -> 50 GeV bins
    )),
    ("muon_pt", Variable(
        "muon_pt", "p_{T}^{#mu}", "GeV",
        {"muon": [["pt_mu1", "NN_Mu_pt"]]},
        valid_channels=["muon"],
        x_min=0.0, x_max=900.0,
        aliases=["mu_pt", "pt_mu", "pt_mu1"],
        rebin=2,  # pt_mu1 is 10 GeV bins -> 20 GeV to match jet pt
    )),
    ("muon_eta", Variable(
        "muon_eta", "#eta^{#mu}", None,
        {"muon": [["eta_mu1", "NN_Mu_eta"]]},
        valid_channels=["muon"],
        x_min=-2.5, x_max=2.5,
        aliases=["mu_eta", "eta_mu", "eta_mu1"],
    )),
    ("electron_pt", Variable(
        "electron_pt", "p_{T}^{e}", "GeV",
        {"electron": [["pt_ele1", "NN_Ele_pt"]]},
        valid_channels=["electron"],
        x_min=0.0, x_max=900.0,
        aliases=["ele_pt", "pt_ele", "pt_ele1"],
        rebin=2,  # pt_ele1 is 10 GeV bins -> 20 GeV to match jet pt
    )),
    ("electron_eta", Variable(
        "electron_eta", "#eta^{e}", None,
        {"electron": [["eta_ele1", "NN_Ele_eta"]]},
        valid_channels=["electron"],
        x_min=-2.5, x_max=2.5,
        aliases=["ele_eta", "eta_ele", "eta_ele1"],
    )),
    # Flavor-merged lepton (used in the combined-channel plot): muon pt/eta in the muon
    # channel, electron pt/eta in the electron channel, summed into one distribution.
    ("lepton_pt", Variable(
        "lepton_pt", "p_{T}^{lepton}", "GeV",
        {"muon": [["pt_mu1", "NN_Mu_pt"]], "electron": [["pt_ele1", "NN_Ele_pt"]]},
        valid_channels=["muon", "electron"],
        x_min=0.0, x_max=900.0,
        aliases=["lep_pt"],
        rebin=2,  # pt_mu1/pt_ele1 are 10 GeV bins -> 20 GeV bins
    )),
    ("lepton_eta", Variable(
        "lepton_eta", "#eta^{lepton}", None,
        {"muon": [["eta_mu1", "NN_Mu_eta"]], "electron": [["eta_ele1", "NN_Ele_eta"]]},
        valid_channels=["muon", "electron"],
        x_min=-2.5, x_max=2.5,
        aliases=["lep_eta"],
    )),
    ("jet1_pt", Variable(
        "jet1_pt", "p_{T}^{AK4 jet 1}", "GeV",
        [["pt_jet1", "NN_Ak4_j1_pt"]],
        x_min=0.0, x_max=900.0,
        aliases=["jet_pt", "ak4_pt", "pt_jet1"],
    )),
    ("jet1_eta", Variable(
        "jet1_eta", "#eta^{AK4 jet 1}", None,
        [["eta_jet1", "NN_Ak4_j1_eta"]],
        x_min=-2.5, x_max=2.5,
        aliases=["jet_eta", "ak4_eta", "eta_jet1"],
    )),
    ("drmin_pt", Variable(
        "drmin_pt", "#DeltaR_{min}(lepton, jet)", None,
        {
            "muon": [["dRmin_pt_mu1", "dRmin_mu1_jet"]],
            "electron": [["dRmin_pt_ele1", "dRmin_ele1_jet"]],
        },
        valid_channels=["muon", "electron"],
        x_min=0.0, x_max=3.0,
        projection="x",
        aliases=["deltaRmin_pt", "dRmin", "drmin"],
    )),
])


# pdf: 100 NNPDF replica variations (ttbar only). Combined in QUADRATURE to match the
# 100 separate pdf1..pdf100 shape nuisances in the datacard. Defined before
# WEIGHT_SYST_SOURCES because it is appended to that dict just below.
PDF_SOURCE = "pdf"
PDF_MEMBERS = ["pdf_{0}".format(i) for i in range(1, 101)]

# ---------------------------------------------------------------------------
# In-file weight-based systematics (produced by ZprimeSemiLeptonicKinematicSystHists).
# These live in a dedicated folder of the NOMINAL workdir (one histogram per
# variation, named "<base>_<suffix>"), unlike the object/shape systematics which
# come from separate _up/_down workdirs. Each source below groups the variation
# suffixes that belong to one nuisance; the per-bin uncertainty for a source is the
# envelope max_i |variation_i - nominal|. Keep the suffix set in sync with
# ZprimeSemiLeptonicKinematicSystHists::syst_suffixes() in the C++.
# ---------------------------------------------------------------------------
WEIGHT_SYST_SOURCES = OrderedDict([
    ("ele_reco",        ["ele_reco_up", "ele_reco_down"]),
    ("ele_id",          ["ele_id_up", "ele_id_down"]),
    ("ele_trigger",     ["ele_trigger_up", "ele_trigger_down"]),
    ("mu_reco",         ["mu_reco_up", "mu_reco_down"]),
    ("mu_iso_stat",     ["mu_iso_stat_up", "mu_iso_stat_down"]),
    ("mu_iso_syst",     ["mu_iso_syst_up", "mu_iso_syst_down"]),
    ("mu_id_stat",      ["mu_id_stat_up", "mu_id_stat_down"]),
    ("mu_id_syst",      ["mu_id_syst_up", "mu_id_syst_down"]),
    ("mu_trigger_stat", ["mu_trigger_stat_up", "mu_trigger_stat_down"]),
    ("mu_trigger_syst", ["mu_trigger_syst_up", "mu_trigger_syst_down"]),
    ("pu",              ["pu_up", "pu_down"]),
    ("prefiring",       ["prefiring_up", "prefiring_down"]),
    ("murmuf",          ["murmuf_upup", "murmuf_upnone", "murmuf_noneup",
                         "murmuf_nonedown", "murmuf_downnone", "murmuf_downdown"]),
    ("btag_cferr1",     ["btag_cferr1_up", "btag_cferr1_down"]),
    ("btag_cferr2",     ["btag_cferr2_up", "btag_cferr2_down"]),
    ("btag_hf",         ["btag_hf_up", "btag_hf_down"]),
    ("btag_hfstats1",   ["btag_hfstats1_up", "btag_hfstats1_down"]),
    ("btag_hfstats2",   ["btag_hfstats2_up", "btag_hfstats2_down"]),
    ("btag_lf",         ["btag_lf_up", "btag_lf_down"]),
    ("btag_lfstats1",   ["btag_lfstats1_up", "btag_lfstats1_down"]),
    ("btag_lfstats2",   ["btag_lfstats2_up", "btag_lfstats2_down"]),
    ("ttag_corr",       ["ttag_corr_up", "ttag_corr_down"]),
    ("ttag_uncorr",     ["ttag_uncorr_up", "ttag_uncorr_down"]),
    ("tmistag",         ["tmistag_up", "tmistag_down"]),
    ("ewk",             ["ewk_up", "ewk_down"]),
    ("qcd",             ["qcd_up", "qcd_down"]),
    ("isr",             ["isr_up", "isr_down"]),
    ("fsr",             ["fsr_up", "fsr_down"]),
])
# pdf added programmatically (100 members); combined in quadrature, not envelope.
WEIGHT_SYST_SOURCES[PDF_SOURCE] = list(PDF_MEMBERS)

# Map each plotter variable to the in-file base histogram name(s). Two-entry lists
# (the reco tops) are summed, exactly as the central nominal sums them.
INFILE_BASE_NAMES = {
    "top_pt":       ["toplep_pt", "tophad_pt"],
    "top_eta":      ["toplep_eta", "tophad_eta"],
    "met_pt":       ["MET_rebin"],
    "met_phi":      ["NN_MET_phi"],
    "ak8_pt":       ["pt_AK8PuppiTaggedjet1"],
    "ak8_eta":      ["eta_AK8PuppiTaggedjet1"],
    "ak8_msd":      ["mSD_AK8PuppiTaggedjet1"],
    "muon_pt":      ["pt_mu1"],
    "muon_eta":     ["eta_mu1"],
    "electron_pt":  ["pt_ele1"],
    "electron_eta": ["eta_ele1"],
    "jet1_pt":      ["pt_jet1"],
    "jet1_eta":     ["eta_jet1"],
    # Flavor-merged lepton: channel-dependent in-file base names.
    "lepton_pt":    {"muon": ["pt_mu1"],  "electron": ["pt_ele1"]},
    "lepton_eta":   {"muon": ["eta_mu1"], "electron": ["eta_ele1"]},
}

# Variables that only make sense in the combined-channel (flavor-merged) plot;
# excluded from the per-channel modes to avoid duplicating muon_pt/electron_pt.
COMBINED_ONLY_VARS = set(["lepton_pt", "lepton_eta"])


def _infile_base_for(variable, channel):
    """Resolve the in-file base histogram name(s) for a variable in a given channel.
    Supports both a plain list and a {channel: list} mapping (flavor-merged leptons)."""
    base = INFILE_BASE_NAMES.get(variable.key)
    if isinstance(base, dict):
        return base.get(channel, [])
    return base

# Per-process cross-section normalization uncertainties (flat lnN from the datacard),
# keyed by the MC sample names in get_mc_globs(). Value = fractional size (lnN - 1).
#   ttbar 5%, V+jets / ST / QCD / Diboson 20%, Toponium 50%.
CROSS_SECTION_UNC = {
    "TTToSemiLeptonic": 0.05,
    "TTToOthers":       0.05,
    "Toponium":         0.50,
    "WJets":            0.20,
    "DY":               0.20,
    "ST":               0.20,
    "QCD":              0.20,
    "Diboson":          0.20,
}


def get_mc_globs(year):
    y = "_" + year
    return OrderedDict([
        ("TTToSemiLeptonic", [PREFIX_MC + "TTToSemiLeptonic" + y + "_*.root"]),
        ("TTToOthers",      [PREFIX_MC + "TTTo2L2Nu" + y + "_*.root",
                             PREFIX_MC + "TTToHadronic" + y + "_*.root"]),
        ("Toponium",        [PREFIX_MC + "EtaT*" + y + "_*.root"]),
        ("WJets",           [PREFIX_MC + "WJetsToLNu_HT-*" + y + "_*.root"]),
        ("DY",              [PREFIX_MC + "DYJetsToLL_M-50_HT-*" + y + "_*.root"]),
        ("ST",              [PREFIX_MC + "ST_*" + y + "_*.root"]),
        ("QCD",             [PREFIX_MC + "QCD_HT*" + y + "_*.root"]),
        ("Diboson",         [PREFIX_MC + "WW" + y + "_*.root",
                             PREFIX_MC + "WZ" + y + "_*.root",
                             PREFIX_MC + "ZZ" + y + "_*.root"]),
    ])


def get_data_globs(year):
    return [PREFIX_DATA + "*" + "_" + year + "_*.root"]


def default_nominal_dir(opts, channel, year=None):
    year = year or opts.year
    return os.path.join(
        opts.dnn_basedir,
        year,
        channel,
        "workdir_AnalysisDNN_{0}_{1}_templatemethod_newsystematics".format(
            year, channel, opts.dnn_ttbar_suffix
        ),
    )


def expand_files(directory, patterns, max_files=None):
    out = []
    for pat in patterns:
        out.extend(sorted(pyglob.glob(os.path.join(directory, pat))))
    if max_files and max_files > 0:
        return out[:max_files]
    return out


def clone_hist(hist, name):
    out = hist.Clone(name)
    out.SetDirectory(0)
    return out


def project_if_needed(hist, projection, name):
    if hist is None:
        return None
    if hist.InheritsFrom("TH2"):
        if projection == "y":
            out = hist.ProjectionY(name)
        else:
            out = hist.ProjectionX(name)
        out.SetDirectory(0)
        return out
    return clone_hist(hist, name)


def read_first_existing_hist(tfile, folder, hist_names, projection, clone_name):
    for hist_name in hist_names:
        hist = tfile.Get("{0}/{1}".format(folder, hist_name))
        if hist:
            return project_if_needed(hist, projection, clone_name), hist_name
    return None, None


def merge_hist_from_files(directory, patterns, folder, variable, channel, max_files=None):
    hist_names = variable.names_for_channel(channel)
    files = expand_files(directory, patterns, max_files=max_files)
    merged = None
    used_files = 0
    used_hist_names = OrderedDict()

    for ifile, fpath in enumerate(files):
        tfile = root.TFile.Open(fpath, "READ")
        if not tfile or tfile.IsZombie():
            continue

        per_file = None
        for ihist, hist_name_group in enumerate(hist_names):
            candidates = hist_name_group
            if isinstance(candidates, str):
                candidates = [candidates]
            hist, used_name = read_first_existing_hist(
                tfile,
                folder,
                candidates,
                variable.projection,
                "tmp_{0}_{1}_{2}".format(variable.key, ifile, ihist),
            )
            if hist is None:
                continue
            used_hist_names[used_name] = True
            if per_file is None:
                per_file = hist
            else:
                per_file.Add(hist)

        tfile.Close()

        if per_file is None:
            continue
        if merged is None:
            merged = clone_hist(per_file, "merged_" + variable.key)
        else:
            merged.Add(per_file)
        used_files += 1

    if merged is not None and getattr(variable, "rebin", 1) > 1:
        merged.Rebin(variable.rebin)
    return merged, used_files, list(used_hist_names.keys())


def same_binning(a, b):
    if a is None or b is None:
        return False
    if a.GetNbinsX() != b.GetNbinsX():
        return False
    for ibin in range(1, a.GetNbinsX() + 1):
        if abs(a.GetXaxis().GetBinLowEdge(ibin) - b.GetXaxis().GetBinLowEdge(ibin)) > 1e-6:
            return False
        if abs(a.GetXaxis().GetBinUpEdge(ibin) - b.GetXaxis().GetBinUpEdge(ibin)) > 1e-6:
            return False
    return True


def sanitize_name(text):
    return re.sub(r"[^A-Za-z0-9_]+", "_", text).strip("_")


def parse_syst_dir(path, year, channel):
    base = os.path.basename(path.rstrip("/"))
    match = re.search(r"(^|_)(up|down)($|_)", base)
    if not match:
        return None
    sign = match.group(2)
    group = base[:match.start(2)] + base[match.end(2):]
    group = group.replace("__", "_").strip("_")

    prefixes = [
        "workdir_AnalysisDNN_{0}_{1}_".format(year, channel),
        "workdir_Analysis_{0}_{1}_".format(year, channel),
    ]
    label = group
    for prefix in prefixes:
        if label.startswith(prefix):
            label = label[len(prefix):]
            break
    if not label:
        label = group
    return sanitize_name(label), sign


def discover_syst_dirs(nominal_dir, opts, channel, year=None):
    year = year or opts.year
    parents = []
    if not opts.no_auto_systs:
        parents.append(os.path.dirname(nominal_dir))
    if opts.syst_parent:
        parents.extend(opts.syst_parent)

    groups = OrderedDict()
    for parent in parents:
        parent = os.path.abspath(parent)
        if not os.path.isdir(parent):
            print("  WARNING: systematic parent does not exist: {0}".format(parent))
            continue
        pattern = opts.syst_pattern.format(year=year, channel=channel)
        for path in sorted(pyglob.glob(os.path.join(parent, pattern))):
            if not os.path.isdir(path):
                continue
            parsed = parse_syst_dir(path, year, channel)
            if not parsed:
                continue
            label, sign = parsed
            if label not in groups:
                groups[label] = {}
            groups[label][sign] = path

    # Keep one-sided groups too. The missing side will be treated as nominal.
    return groups


def merge_nominal_processes(nominal_dir, folder, variable, year, channel, max_files=None):
    merged = OrderedDict()
    for sample, patterns in get_mc_globs(year).items():
        hist, nfiles, used_names = merge_hist_from_files(
            nominal_dir, patterns, folder, variable, channel, max_files=max_files
        )
        if hist is not None:
            hist.SetName(variable.key + "_" + sample)
            merged[sample] = hist
            print("  {0:<18s} {1:>5d} files  integral={2:>12.2f}  hists={3}".format(
                sample, nfiles, hist.Integral(), ",".join(used_names)
            ))
        else:
            print("  {0:<18s} WARNING: no nominal histogram".format(sample))
    return merged


def merge_data(nominal_dir, folder, variable, year, channel, max_files=None):
    hist, nfiles, used_names = merge_hist_from_files(
        nominal_dir, get_data_globs(year), folder, variable, channel, max_files=max_files
    )
    if hist is not None:
        hist.SetName(variable.key + "_DATA")
        print("  {0:<18s} {1:>5d} files  integral={2:>12.2f}  hists={3}".format(
            "DATA", nfiles, hist.Integral(), ",".join(used_names)
        ))
    else:
        print("  DATA               WARNING: no data histogram")
    return hist


def build_total(hist_map, name):
    total = None
    for sample in SAMPLE_COLORS:
        hist = hist_map.get(sample)
        if hist is None:
            continue
        if total is None:
            total = clone_hist(hist, name)
        else:
            total.Add(hist)
    return total


def merge_syst_total(syst_dir, nominal_hists, folder, variable, year, channel, max_files=None):
    varied = OrderedDict()
    for sample, nominal in nominal_hists.items():
        hist, nfiles, used_names = merge_hist_from_files(
            syst_dir, get_mc_globs(year)[sample], folder, variable, channel, max_files=max_files
        )
        if hist is None:
            varied[sample] = nominal
        elif same_binning(hist, nominal):
            varied[sample] = hist
        else:
            print("    WARNING: {0} has incompatible binning in {1}; using nominal".format(
                sample, os.path.basename(syst_dir)
            ))
            varied[sample] = nominal
    return build_total(varied, "total_" + variable.key + "_syst")


def accumulate_systematics(nominal_total, nominal_hists, syst_groups, folder, variable, year, channel, max_files=None):
    nbins = nominal_total.GetNbinsX()
    err2_up = [0.0] * (nbins + 2)
    err2_down = [0.0] * (nbins + 2)
    used = []

    for label, sides in syst_groups.items():
        up_total = None
        down_total = None
        if "up" in sides:
            up_total = merge_syst_total(
                sides["up"], nominal_hists, folder, variable, year, channel, max_files=max_files
            )
        if "down" in sides:
            down_total = merge_syst_total(
                sides["down"], nominal_hists, folder, variable, year, channel, max_files=max_files
            )
        if up_total is None and down_total is None:
            continue
        if up_total is not None and not same_binning(up_total, nominal_total):
            print("    WARNING: skipping {0} up; incompatible total binning".format(label))
            up_total = None
        if down_total is not None and not same_binning(down_total, nominal_total):
            print("    WARNING: skipping {0} down; incompatible total binning".format(label))
            down_total = None
        if up_total is None and down_total is None:
            continue

        used.append(label)
        for ibin in range(1, nbins + 1):
            nom = nominal_total.GetBinContent(ibin)
            dup = 0.0 if up_total is None else up_total.GetBinContent(ibin) - nom
            ddn = 0.0 if down_total is None else down_total.GetBinContent(ibin) - nom
            up_pos = dup >= 0.0
            dn_pos = ddn >= 0.0

            if up_pos and not dn_pos:
                err2_up[ibin] += dup * dup
                err2_down[ibin] += ddn * ddn
            elif dn_pos and not up_pos:
                err2_up[ibin] += ddn * ddn
                err2_down[ibin] += dup * dup
            elif up_pos and dn_pos:
                err2_up[ibin] += max(dup, ddn) ** 2
            else:
                err2_down[ibin] += min(dup, ddn) ** 2

    return err2_up, err2_down, used


def make_uncertainty_graph(total, syst_up2, syst_down2, lumi_unc, name, ratio=False,
                           fill_color=None, fill_style=3254):
    nbins = total.GetNbinsX()
    x = array("d")
    y = array("d")
    exl = array("d")
    exh = array("d")
    eyl = array("d")
    eyh = array("d")

    for ibin in range(1, nbins + 1):
        nom = total.GetBinContent(ibin)
        stat = total.GetBinError(ibin)
        lumi = lumi_unc * nom
        err_up = math.sqrt(stat * stat + lumi * lumi + syst_up2[ibin])
        err_down = math.sqrt(stat * stat + lumi * lumi + syst_down2[ibin])

        x.append(total.GetBinCenter(ibin))
        exl.append(total.GetBinWidth(ibin) / 2.0)
        exh.append(total.GetBinWidth(ibin) / 2.0)
        if ratio:
            y.append(1.0)
            if nom > 0.0:
                eyl.append(err_down / nom)
                eyh.append(err_up / nom)
            else:
                eyl.append(0.0)
                eyh.append(0.0)
        else:
            y.append(nom)
            eyl.append(err_down)
            eyh.append(err_up)

    graph = root.TGraphAsymmErrors(nbins, x, y, exl, exh, eyl, eyh)
    graph.SetName(name)
    graph.SetLineWidth(0)
    graph.SetMarkerSize(0)
    graph.SetFillColor(root.kGray + 1 if fill_color is None else fill_color)
    graph.SetFillStyle(fill_style)
    return graph


def set_ratio_axis_range(nice, ymin, ymax):
    if not getattr(nice, "show_ratio", False):
        return
    if hasattr(nice, "ratio_null_hist") and nice.ratio_null_hist:
        nice.ratio_null_hist.SetMinimum(ymin)
        nice.ratio_null_hist.SetMaximum(ymax)
    if hasattr(nice, "pad_ratio") and nice.pad_ratio:
        nice.pad_ratio.Modified()
        nice.pad_ratio.Update()


def hide_empty_data_points(nice):
    if nice.blind_data or not hasattr(nice, "data_hist"):
        return
    graphs = []
    if hasattr(nice, "data"):
        graphs.append(nice.data)
    if hasattr(nice, "ratio_data"):
        graphs.append(nice.ratio_data)
    for ibin in range(nice.data_hist.GetNbinsX(), 0, -1):
        if nice.data_hist.GetBinContent(ibin) > 0.0:
            continue
        ipoint = ibin - 1
        for graph in graphs:
            if graph and ipoint < graph.GetN():
                graph.RemovePoint(ipoint)


def write_plot_input(path, variable, nominal_hists, data_hist, total_hist):
    tfile = root.TFile(path, "RECREATE")
    for sample in SAMPLE_COLORS:
        hist = nominal_hists.get(sample)
        if hist is None:
            continue
        out = clone_hist(hist, variable.key + "_" + sample)
        out.Write()
    data_out = clone_hist(data_hist, variable.key + "_DATA")
    data_out.Write()
    total_out = clone_hist(total_hist, "TotalProcs")
    total_out.Write()
    tfile.Close()


def make_processes(variable, nominal_hists):
    out = []
    idx = 0
    for sample, info in SAMPLE_COLORS.items():
        if sample not in nominal_hists:
            continue
        proc = Process(variable.key + "_" + sample, info["legend"], info["color"])
        proc.index = idx
        out.append(proc)
        idx += 1
    return out


def draw_legend(nice, processes, total_unc_graph, unc_label, stat_graph=None, stat_label="MC stat."):
    nice.pad_main.cd()
    n_entries = len(processes) + 2 + (1 if stat_graph is not None else 0)
    # Top-right corner: anchor the top edge at 0.88 and grow downward.
    legend_top = 0.88
    legend = root.TLegend(
        nice.coord.graph_to_pad_x(0.58),
        nice.coord.graph_to_pad_y(legend_top - 0.040 * n_entries),
        nice.coord.graph_to_pad_x(0.95),
        nice.coord.graph_to_pad_y(legend_top),
    )
    legend.SetTextSize(0.023)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    if not nice.blind_data:
        legend.AddEntry(nice.data_hist, "Data", "ep")
    legend.AddEntry(total_unc_graph, unc_label, "f")
    if stat_graph is not None:
        legend.AddEntry(stat_graph, stat_label, "f")
    for proc in processes:
        legend.AddEntry(nice.stack.GetStack().At(proc.index), proc.legend, "f")
    legend.Draw()
    return legend


def read_infile_weight_hists(nominal_dir, infile_folder, base_names, suffixes, year, channel, max_files=None, rebin=1):
    """Open each MC file once and sum the in-file nominal + weight-variation histograms.

    Returns a dict: 'nom' -> summed nominal TH1, '<suffix>' -> summed variation TH1.
    Each entry is the sum over base_names and over all MC samples (no nominal fallback),
    rebinned by `rebin` so it matches the (rebinned) central nominal.
    """
    keys = ["nom"] + list(suffixes)
    out = {}
    for sample, patterns in get_mc_globs(year).items():
        files = expand_files(nominal_dir, patterns, max_files=max_files)
        for fpath in files:
            tfile = root.TFile.Open(fpath, "READ")
            if not tfile or tfile.IsZombie():
                continue
            for key in keys:
                for base in base_names:
                    hname = base if key == "nom" else (base + "_" + key)
                    hist = tfile.Get("{0}/{1}".format(infile_folder, hname))
                    if not hist:
                        continue
                    if out.get(key) is None:
                        clone = hist.Clone("infile_{0}_{1}".format(channel, key))
                        clone.SetDirectory(0)
                        out[key] = clone
                    else:
                        out[key].Add(hist)
    if rebin > 1:
        for key in out:
            if out[key] is not None:
                out[key].Rebin(rebin)
            tfile.Close()
    return out


def accumulate_infile_weight_systematics(nominal_total, infile_hists):
    """Per-bin envelope of the in-file weight systematics, added symmetrically in quadrature.

    For each source the per-bin uncertainty is max_i |variation_i - nominal|; it is added to
    both the up and down sides. Returns (err2_up, err2_down, used_sources).
    """
    nbins = nominal_total.GetNbinsX()
    err2_up = [0.0] * (nbins + 2)
    err2_down = [0.0] * (nbins + 2)
    used = []
    nom_hist = infile_hists.get("nom")
    if nom_hist is None:
        return err2_up, err2_down, used
    if not same_binning(nom_hist, nominal_total):
        print("    WARNING: in-file nominal binning differs from central nominal; "
              "skipping in-file weight systematics")
        return err2_up, err2_down, used
    for source, members in WEIGHT_SYST_SOURCES.items():
        member_hists = []
        for suf in members:
            hist = infile_hists.get(suf)
            if hist is not None and same_binning(hist, nom_hist):
                member_hists.append(hist)
        if not member_hists:
            continue
        used.append(source)
        if source == PDF_SOURCE:
            # pdf: quadrature over the 100 replicas (matches 100 separate datacard nuisances)
            for ibin in range(1, nbins + 1):
                nom = nom_hist.GetBinContent(ibin)
                sumsq = 0.0
                for hist in member_hists:
                    deviation = hist.GetBinContent(ibin) - nom
                    sumsq += deviation * deviation
                err2_up[ibin] += sumsq
                err2_down[ibin] += sumsq
        else:
            # everything else: symmetric per-source envelope
            for ibin in range(1, nbins + 1):
                nom = nom_hist.GetBinContent(ibin)
                env = 0.0
                for hist in member_hists:
                    deviation = abs(hist.GetBinContent(ibin) - nom)
                    if deviation > env:
                        env = deviation
                err2_up[ibin] += env * env
                err2_down[ibin] += env * env
    return err2_up, err2_down, used


def accumulate_crosssection_norm(nominal_total, nominal_hists):
    """Per-bin flat cross-section normalization uncertainty (datacard lnN), processes
    added in quadrature. Each process p contributes unc_p * content_p(bin)."""
    nbins = nominal_total.GetNbinsX()
    err2 = [0.0] * (nbins + 2)
    used = []
    for sample, hist in nominal_hists.items():
        unc = CROSS_SECTION_UNC.get(sample)
        if unc is None or hist is None:
            continue
        used.append(sample)
        for ibin in range(1, nbins + 1):
            contribution = unc * hist.GetBinContent(ibin)
            err2[ibin] += contribution * contribution
    return err2, used


def plot_variable(opts, channel, variable, nominal_dir, syst_groups):
    if not variable.applies_to_channel(channel):
        print("\nSkipping {0}: not defined for {1} channel".format(variable.key, channel))
        return False

    nominal_folder = opts.folder
    syst_folder = opts.syst_folder or opts.folder
    print("\n" + "=" * 80)
    print("Plot {0} ({1}, {2})".format(variable.key, opts.year, channel))
    print("Nominal dir : {0}".format(nominal_dir))
    print("Nominal hist: {0} / {1}".format(nominal_folder, variable.names_for_channel(channel)))
    if syst_groups:
        print("Systematics: {0}".format(", ".join(syst_groups.keys())))
        print("Syst folder: {0}".format(syst_folder))
    else:
        print("Systematics: none found; band will be MC stat + lumi only")
    print("=" * 80)

    nominal_hists = merge_nominal_processes(
        nominal_dir, nominal_folder, variable, opts.year, channel, max_files=opts.max_files
    )
    data_hist = merge_data(
        nominal_dir, nominal_folder, variable, opts.year, channel, max_files=opts.max_files
    )

    if not nominal_hists or data_hist is None:
        print("  ERROR: missing nominal MC or DATA for {0}; skipping".format(variable.key))
        return False

    total = build_total(nominal_hists, "TotalProcs")
    if total is None or total.Integral() <= 0.0:
        print("  ERROR: total MC is empty for {0}; skipping".format(variable.key))
        return False

    syst_up2, syst_down2, used_systs = accumulate_systematics(
        total,
        nominal_hists,
        syst_groups,
        syst_folder,
        variable,
        opts.year,
        channel,
        max_files=opts.max_files,
    )
    print("  Used {0} workdir systematic source(s): {1}".format(
        len(used_systs), ", ".join(used_systs) if used_systs else "none"
    ))

    # In-file weight-based systematics (pileup, b-tag, lepton SFs, scale, PS, EWK/QCD, ...).
    infile_base = _infile_base_for(variable, channel)
    if not opts.no_infile_systs and infile_base:
        all_suffixes = []
        for members in WEIGHT_SYST_SOURCES.values():
            all_suffixes.extend(members)
        infile_hists = read_infile_weight_hists(
            nominal_dir, opts.infile_syst_folder, infile_base,
            all_suffixes, opts.year, channel, max_files=opts.max_files,
            rebin=variable.rebin,
        )
        infile_up2, infile_down2, used_infile = accumulate_infile_weight_systematics(total, infile_hists)
        for ibin in range(len(syst_up2)):
            syst_up2[ibin] += infile_up2[ibin]
            syst_down2[ibin] += infile_down2[ibin]
        if used_infile:
            used_systs = list(used_systs) + ["weight:" + s for s in used_infile]
        print("  Used {0} in-file weight syst source(s): {1}".format(
            len(used_infile), ", ".join(used_infile) if used_infile else "none"
        ))

    # Per-process cross-section normalization (flat lnN from the datacard).
    if not opts.no_xsec_norm:
        xsec_err2, used_xsec = accumulate_crosssection_norm(total, nominal_hists)
        for ibin in range(len(syst_up2)):
            syst_up2[ibin] += xsec_err2[ibin]
            syst_down2[ibin] += xsec_err2[ibin]
        if used_xsec:
            used_systs = list(used_systs) + ["xsec_norm"]
        print("  Added cross-section norm for {0} process(es)".format(len(used_xsec)))

    print("  Total MC integral: {0:.2f}".format(total.Integral()))
    print("  Data integral    : {0:.2f}".format(data_hist.Integral()))
    if total.Integral() > 0.0:
        print("  Data/MC integral : {0:.4f}".format(data_hist.Integral() / total.Integral()))

    lumi_info = _YEARS.get(opts.year, {})
    lumi_text = lumi_info.get("lumi_fb_display", "") + " fb^{#minus1} (13 TeV)"
    lumi_unc = opts.lumi_unc if opts.lumi_unc >= 0.0 else lumi_info.get("lumi_unc", 0.0)
    return _draw_and_save(opts, channel, variable, opts.year, nominal_hists, data_hist, total,
                          syst_up2, syst_down2, used_systs, lumi_text, lumi_unc)


def _draw_and_save(opts, channel, variable, year_label, nominal_hists, data_hist, total,
                   syst_up2, syst_down2, used_systs, lumi_text, lumi_unc, label_text=None):
    """Render the stacked data/MC plot with ratio and the uncertainty band, and save it.
    Shared by the per-year and combined-years paths."""
    channel_outdir = os.path.join(opts.outdir, "{0}_{1}".format(year_label, channel))
    if not os.path.isdir(channel_outdir):
        os.makedirs(channel_outdir)
    tmp_root = os.path.join(channel_outdir, "_tmp_{0}_{1}_{2}.root".format(year_label, channel, variable.key))
    write_plot_input(tmp_root, variable, nominal_hists, data_hist, total)

    if label_text is None:
        label_text = "{0} channel ({1})".format(channel, year_label)

    processes = make_processes(variable, nominal_hists)

    # y-axis floor: --y-min if given (>0), else 1 event for log / 0 for linear.
    if getattr(opts, "y_min", 0.0) and opts.y_min > 0.0:
        y_axis_min = opts.y_min
    elif opts.logy:
        y_axis_min = 1.0
    else:
        y_axis_min = 0.0

    # y-axis ceiling: --y-max if given (>0); else default to 1e9 for log, auto (None) for linear.
    if getattr(opts, "y_max", 0.0) and opts.y_max > 0.0:
        y_axis_max = opts.y_max
    elif opts.logy:
        y_axis_max = 1e9
    else:
        y_axis_max = None

    nice = NiceStackWithRatio(
        infile_path=tmp_root,
        x_axis_title=variable.title,
        x_axis_unit=variable.unit,
        x_axis_min=variable.x_min,
        x_axis_max=variable.x_max,
        prepostfit="prefitRaw",
        processes=processes,
        signals=[],
        syst_names=[],
        lumi_unc=0.0,
        divide_by_bin_width=False,
        data_name=variable.key + "_DATA",
        text_prelim=opts.text_prelim,
        text_simulation=None,
        text_top_left=label_text,
        text_top_right=lumi_text,
        y_axis_min=y_axis_min,
        y_axis_max=y_axis_max,
        nostack=False,
        logy=opts.logy,
        blind_data=False,
        show_ratio=True,
    )
    nice.plot()
    hide_empty_data_points(nice)
    set_ratio_axis_range(nice, opts.ratio_min, opts.ratio_max)

    band = make_uncertainty_graph(total, syst_up2, syst_down2, lumi_unc, "total_unc_" + variable.key, ratio=False)
    ratio_band = make_uncertainty_graph(total, syst_up2, syst_down2, lumi_unc, "ratio_unc_" + variable.key, ratio=True)
    nice._manual_objects = [band, ratio_band]

    # Separate MC-statistics-only band (GetBinError; no syst, no lumi), solid light gray,
    # drawn on top of the total band so the hatched ring shows the systematic part.
    stat_band = None
    ratio_stat_band = None
    if not opts.no_stat_band:
        zeros = [0.0] * (total.GetNbinsX() + 2)
        stat_band = make_uncertainty_graph(total, zeros, zeros, 0.0, "stat_unc_" + variable.key,
                                           ratio=False, fill_color=root.kGray, fill_style=1001)
        ratio_stat_band = make_uncertainty_graph(total, zeros, zeros, 0.0, "ratio_stat_" + variable.key,
                                                 ratio=True, fill_color=root.kGray, fill_style=1001)
        nice._manual_objects += [stat_band, ratio_stat_band]

    nice.pad_main.cd()
    band.Draw("2 same")
    if stat_band is not None:
        stat_band.Draw("2 same")
    if not nice.blind_data:
        nice.data.Draw("pz0 same")
    nice.pad_main.RedrawAxis()

    nice.pad_ratio.cd()
    ratio_band.Draw("2 same")
    if ratio_stat_band is not None:
        ratio_stat_band.Draw("2 same")
    nice.ratio_data.Draw("pz0 same")
    nice.pad_ratio.RedrawAxis()

    unc_label = "MC stat. #oplus syst."
    if not used_systs:
        unc_label = "MC stat. #oplus lumi"
    legend = draw_legend(nice, processes, band, unc_label, stat_graph=stat_band, stat_label="MC stat.")
    nice._manual_objects.append(legend)

    for ext in opts.formats.split(","):
        ext = ext.strip().lstrip(".")
        if not ext:
            continue
        outpath = os.path.join(channel_outdir, variable.key + "_dataMC_withSyst." + ext)
        nice.save_plot(outpath)
        print("  Saved: {0}".format(outpath))

    nice.canvas.Close()
    if (not opts.keep_root) and os.path.isfile(tmp_root):
        os.remove(tmp_root)
    return True


def resolve_variables(text):
    aliases = {}
    for key, variable in VARIABLES.items():
        aliases[key.lower()] = key
        for alias in variable.aliases:
            aliases[alias.lower()] = key

    if text.strip().lower() == "all":
        return list(VARIABLES.keys())

    selected = []
    for raw in text.replace(";", ",").split(","):
        item = raw.strip()
        if not item:
            continue
        key = aliases.get(item.lower())
        if key is None:
            print("ERROR: unknown variable '{0}'. Known variables: {1}".format(
                item, ", ".join(VARIABLES.keys())
            ))
            sys.exit(1)
        if key not in selected:
            selected.append(key)
    return selected


def parse_options():
    parser = OptionParser()
    parser.add_option("-y", "--year", dest="year", default="all", type="str",
                      help="Year, e.g. UL16preVFP, UL16postVFP, UL17, UL18, or all (Run 2, default)")
    parser.add_option("-c", "--channel", dest="channel", default="all", type="str",
                      help="muon, electron, both (separate plots), combined (muon+electron merged "
                           "into one plot), or all (merged plot PLUS separate muon and electron plots)")
    parser.add_option("--variables", dest="variables", default="all", type="str",
                      help="Comma-separated variables or 'all'. Known: {0}".format(", ".join(VARIABLES.keys())))
    parser.add_option("--folder", dest="folder", default="DNN_output0_General", type="str",
                      help="Nominal ROOT folder containing the histograms")
    parser.add_option("--syst-folder", dest="syst_folder", default=None, type="str",
                      help="ROOT folder in systematic workdirs; default is --folder")
    parser.add_option("--infile-syst-folder", dest="infile_syst_folder",
                      default="KinematicSyst_Inclusive_SR", type="str",
                      help="Folder (in the nominal workdir) holding the in-file weight-systematic "
                           "variations from ZprimeSemiLeptonicKinematicSystHists")
    parser.add_option("--no-infile-systs", dest="no_infile_systs", default=False, action="store_true",
                      help="Disable in-file weight systematics; band uses workdir systematics only")
    parser.add_option("--no-xsec-norm", dest="no_xsec_norm", default=False, action="store_true",
                      help="Disable the per-process cross-section normalization uncertainties in the band")
    parser.add_option("--no-stat-band", dest="no_stat_band", default=False, action="store_true",
                      help="Do not draw the separate MC-statistics-only band (only the total stat+syst band)")
    parser.add_option("--dnn-basedir", dest="dnn_basedir",
                      default="/data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/output_DNN",
                      type="str", help="Base directory for output_DNN")
    parser.add_option("--dnn-ttbar-suffix", dest="dnn_ttbar_suffix",
                      default="ttbar1l_combinedgen", type="str",
                      help="Nominal workdir suffix after templatemethod_newmassbins_")
    parser.add_option("--nominal-dir", dest="nominal_dir", default=None, type="str",
                      help="Override nominal workdir. With --channel both this may contain {channel}.")
    parser.add_option("--syst-parent", dest="syst_parent", default=[], action="append", type="str",
                      help="Additional parent directory to scan for paired _up/_down systematic workdirs")
    parser.add_option("--syst-pattern", dest="syst_pattern",
                      default="workdir_*{year}*{channel}*", type="str",
                      help="Glob used under each systematic parent; supports {year} and {channel}")
    parser.add_option("--no-auto-systs", dest="no_auto_systs", default=False, action="store_true",
                      help="Do not scan siblings of the nominal workdir for systematics")
    parser.add_option("--outdir", dest="outdir", default="DataMC_kinematics_withSyst", type="str",
                      help="Output directory")
    parser.add_option("--formats", dest="formats", default="pdf", type="str",
                      help="Comma-separated output formats, e.g. pdf,png")
    parser.add_option("--logy", dest="logy", default=True, action="store_true",
                      help="Log y-axis (default on).")
    parser.add_option("--linear", "--no-logy", dest="logy", action="store_false",
                      help="Linear y-axis (disables the default log scale).")
    parser.add_option("--y-min", dest="y_min", default=0.0, type="float",
                      help="Y-axis minimum. 0 (default) = auto: 1 event for --logy, 0 for linear. "
                           "Set e.g. --y-min 100 for a tighter log axis.")
    parser.add_option("--y-max", dest="y_max", default=0.0, type="float",
                      help="Y-axis maximum. 0 (default) = auto: 1e9 for --logy, auto-scale for linear. "
                           "Set e.g. --y-max 1e7 to override.")
    parser.add_option("--ratio-min", dest="ratio_min", default=0.5, type="float")
    parser.add_option("--ratio-max", dest="ratio_max", default=1.5, type="float")
    parser.add_option("--lumi-unc", dest="lumi_unc", default=-1.0, type="float",
                      help="Override lumi uncertainty. Default uses constants.py.")
    parser.add_option("--text-prelim", dest="text_prelim", default="Private Work", type="str")
    parser.add_option("--max-files", dest="max_files", default=0, type="int",
                      help="Debug option: cap files per sample/glob")
    parser.add_option("--keep-root", dest="keep_root", default=False, action="store_true",
                      help="Keep temporary merged ROOT files")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# All-years (Run 2) combination.
# ---------------------------------------------------------------------------
ALL_YEARS = ["UL16preVFP", "UL16postVFP", "UL17", "UL18"]
RUN2_LUMI_TEXT = "138 fb^{#minus1} (13 TeV)"
RUN2_LUMI_UNC = 0.016  # full Run 2 luminosity uncertainty (override with --lumi-unc)

# Sources that are INDEPENDENT per year -> combine in quadrature across years.
# Everything else is correlated across years -> shifts add linearly. Matches the datacard:
# JER and the *_stats / per-year lepton-SF nuisances are per-year; JEC, pileup, the
# correlated b-tag and lepton-SF *syst* parts, and all ttbar modeling are correlated.
PER_YEAR_SOURCES = set([
    "JER",
    "btag_hfstats1", "btag_hfstats2", "btag_lfstats1", "btag_lfstats2",
    "mu_id_stat", "mu_iso_stat", "mu_trigger_stat",
    "ele_reco", "ele_id", "ele_trigger",
])


def _canonical_source_key(source, year):
    """Per-year sources get a year-specific key (kept independent across years);
    correlated sources keep a year-agnostic key (their shifts add linearly)."""
    if source in PER_YEAR_SOURCES:
        return source + "@" + year
    return source


def workdir_source_sizes(nominal_total, nominal_hists, syst_groups, folder, variable, year, channel, max_files=None):
    """Per-bin symmetric size max(|up-nom|, |down-nom|) for each workdir systematic source."""
    nbins = nominal_total.GetNbinsX()
    sizes = OrderedDict()
    for label, sides in syst_groups.items():
        up_total = merge_syst_total(sides["up"], nominal_hists, folder, variable, year, channel,
                                    max_files=max_files) if "up" in sides else None
        down_total = merge_syst_total(sides["down"], nominal_hists, folder, variable, year, channel,
                                      max_files=max_files) if "down" in sides else None
        if up_total is not None and not same_binning(up_total, nominal_total):
            up_total = None
        if down_total is not None and not same_binning(down_total, nominal_total):
            down_total = None
        if up_total is None and down_total is None:
            continue
        arr = [0.0] * (nbins + 2)
        for ibin in range(1, nbins + 1):
            nom = nominal_total.GetBinContent(ibin)
            dup = abs(up_total.GetBinContent(ibin) - nom) if up_total is not None else 0.0
            ddn = abs(down_total.GetBinContent(ibin) - nom) if down_total is not None else 0.0
            arr[ibin] = max(dup, ddn)
        sizes[label] = arr
    return sizes


def infile_source_sizes(nominal_total, infile_hists):
    """Per-bin size for each in-file weight source (envelope; pdf = quadrature)."""
    nbins = nominal_total.GetNbinsX()
    sizes = OrderedDict()
    nom_hist = infile_hists.get("nom")
    if nom_hist is None or not same_binning(nom_hist, nominal_total):
        return sizes
    for source, members in WEIGHT_SYST_SOURCES.items():
        member_hists = [infile_hists[s] for s in members
                        if infile_hists.get(s) is not None and same_binning(infile_hists[s], nom_hist)]
        if not member_hists:
            continue
        arr = [0.0] * (nbins + 2)
        for ibin in range(1, nbins + 1):
            nom = nom_hist.GetBinContent(ibin)
            if source == PDF_SOURCE:
                sumsq = 0.0
                for hist in member_hists:
                    dev = hist.GetBinContent(ibin) - nom
                    sumsq += dev * dev
                arr[ibin] = math.sqrt(sumsq)
            else:
                env = 0.0
                for hist in member_hists:
                    dev = abs(hist.GetBinContent(ibin) - nom)
                    if dev > env:
                        env = dev
                arr[ibin] = env
        sizes[source] = arr
    return sizes


def xsec_source_sizes(nominal_total, nominal_hists):
    """Per-bin size for each per-process cross-section normalization (one source per process)."""
    nbins = nominal_total.GetNbinsX()
    sizes = OrderedDict()
    for sample, hist in nominal_hists.items():
        unc = CROSS_SECTION_UNC.get(sample)
        if unc is None or hist is None:
            continue
        arr = [0.0] * (nbins + 2)
        for ibin in range(1, nbins + 1):
            arr[ibin] = unc * hist.GetBinContent(ibin)
        sizes["xsec_" + sample] = arr
    return sizes


def plot_variable_combined(opts, variable, units, year_label, channel_label, lumi_text, lumi_unc):
    """Merged data/MC plot summing over the given (year, channel) units, combining each
    systematic per the datacard's correlation scheme: per-year nuisances combine in
    quadrature across years; lepton SFs are channel-specific (kept separate automatically);
    everything else (JEC, pileup, b-tag, ttbar modeling, ...) is correlated across years
    AND channels. units: list of (year, channel, nominal_dir, syst_groups)."""
    nominal_folder = opts.folder
    syst_folder = opts.syst_folder or opts.folder
    all_suffixes = []
    for members in WEIGHT_SYST_SOURCES.values():
        all_suffixes.extend(members)

    combined_nominal_hists = OrderedDict()
    combined_data = None
    combined_sizes = OrderedDict()   # canonical key -> per-bin size (linear sum within a key)
    used_sources = set()
    used_units = []
    nbins = None

    print("\n" + "=" * 80)
    print("Plot {0} [{1}, {2}]  hist: {3}".format(
        variable.key, year_label, channel_label, variable.names_for_channel(units[0][1]) if units else "?"))
    print("=" * 80)

    for (year, channel, nominal_dir, syst_groups) in units:
        if not variable.applies_to_channel(channel):
            continue
        print("  -> {0} / {1}".format(year, channel))
        nh = merge_nominal_processes(nominal_dir, nominal_folder, variable, year, channel, max_files=opts.max_files)
        dh = merge_data(nominal_dir, nominal_folder, variable, year, channel, max_files=opts.max_files)
        if not nh or dh is None:
            print("  WARNING: missing nominal/data for {0} ({1}/{2}); skipping".format(variable.key, year, channel))
            continue
        total_u = build_total(nh, "TotalProcs_{0}_{1}".format(year, channel))
        if total_u is None or total_u.Integral() <= 0.0:
            print("  WARNING: empty MC for {0} ({1}/{2}); skipping".format(variable.key, year, channel))
            continue
        for sample, hist in nh.items():
            print("       {0:<18s} integral={1:>12.2f}".format(sample, hist.Integral()))
        print("       {0:<18s} integral={1:>12.2f}".format("TOTAL MC", total_u.Integral()))
        print("       {0:<18s} integral={1:>12.2f}".format("DATA", dh.Integral()))
        nbins = total_u.GetNbinsX()
        used_units.append("{0}/{1}".format(year, channel))

        sizes = OrderedDict()
        sizes.update(workdir_source_sizes(total_u, nh, syst_groups, syst_folder, variable, year, channel, max_files=opts.max_files))
        infile_base = _infile_base_for(variable, channel)
        if not opts.no_infile_systs and infile_base:
            infile_hists = read_infile_weight_hists(
                nominal_dir, opts.infile_syst_folder, infile_base,
                all_suffixes, year, channel, max_files=opts.max_files,
                rebin=variable.rebin)
            sizes.update(infile_source_sizes(total_u, infile_hists))
        if not opts.no_xsec_norm:
            sizes.update(xsec_source_sizes(total_u, nh))

        for sample, hist in nh.items():
            if combined_nominal_hists.get(sample) is None:
                clone = hist.Clone("comb_{0}_{1}".format(variable.key, sample))
                clone.SetDirectory(0)
                combined_nominal_hists[sample] = clone
            else:
                combined_nominal_hists[sample].Add(hist)
        if combined_data is None:
            combined_data = dh.Clone("comb_{0}_data".format(variable.key))
            combined_data.SetDirectory(0)
        else:
            combined_data.Add(dh)

        for source, arr in sizes.items():
            used_sources.add(source)
            ck = _canonical_source_key(source, year)
            if ck not in combined_sizes:
                combined_sizes[ck] = [0.0] * (nbins + 2)
            for ibin in range(1, nbins + 1):
                combined_sizes[ck][ibin] += arr[ibin]

    if not combined_nominal_hists or combined_data is None:
        print("  Skipping {0}: no usable (year, channel) units".format(variable.key))
        return False

    combined_total = build_total(combined_nominal_hists, "TotalProcs_comb")
    nbins = combined_total.GetNbinsX()

    err2 = [0.0] * (nbins + 2)
    for ck, arr in combined_sizes.items():
        for ibin in range(1, nbins + 1):
            err2[ibin] += arr[ibin] * arr[ibin]

    print("\n" + "=" * 80)
    print("Combined {0} [{1}, {2}]  units: {3}".format(variable.key, year_label, channel_label, ", ".join(used_units)))
    print("  Systematic sources after combination: {0}".format(len(combined_sizes)))
    denom = combined_total.Integral()
    print("  Total MC: {0:.1f}   Data: {1:.1f}   Data/MC: {2:.4f}".format(
        denom, combined_data.Integral(), (combined_data.Integral() / denom) if denom > 0 else 0.0))
    print("=" * 80)

    # Display label uses the actual channels merged, e.g. "muon+electron (Run2)";
    # a single channel keeps the "<channel> channel (<year>)" form.
    channels_used = []
    for unit_label in used_units:
        ch = unit_label.split("/")[1]
        if ch not in channels_used:
            channels_used.append(ch)
    if len(channels_used) > 1:
        label_text = "{0} ({1})".format("+".join(channels_used), year_label)
    elif channels_used:
        label_text = "{0} channel ({1})".format(channels_used[0], year_label)
    else:
        label_text = "{0} ({1})".format(channel_label, year_label)

    return _draw_and_save(opts, channel_label, variable, year_label, combined_nominal_hists, combined_data,
                          combined_total, err2, err2, sorted(used_sources), lumi_text, lumi_unc,
                          label_text=label_text)


def main():
    opts, _ = parse_options()
    opts.channel = opts.channel.strip().lower()
    if opts.channel not in ("muon", "electron", "both", "combined", "all"):
        print("ERROR: --channel must be muon, electron, both, or combined")
        return 1

    selected = resolve_variables(opts.variables)
    combine_years = opts.year.strip().lower() in ("all", "run2", "combined")
    merge_channels = opts.channel in ("combined", "all")

    years = ALL_YEARS if combine_years else [opts.year]
    year_label = "Run2" if combine_years else opts.year

    def build_units(channel_list):
        units = []
        for ch in channel_list:
            for yr in years:
                if opts.nominal_dir:
                    ydir = opts.nominal_dir.format(year=yr, channel=ch)
                else:
                    ydir = default_nominal_dir(opts, ch, year=yr)
                if not os.path.isdir(ydir):
                    print("  WARNING: missing nominal dir ({0}/{1}): {2}".format(yr, ch, ydir))
                    continue
                units.append((yr, ch, ydir, discover_syst_dirs(ydir, opts, ch, year=yr)))
        return units

    def lumi_for():
        if combine_years:
            return RUN2_LUMI_TEXT, (opts.lumi_unc if opts.lumi_unc >= 0.0 else RUN2_LUMI_UNC)
        info = _YEARS.get(opts.year, {})
        return (info.get("lumi_fb_display", "") + " fb^{#minus1} (13 TeV)",
                opts.lumi_unc if opts.lumi_unc >= 0.0 else info.get("lumi_unc", 0.0))

    made_any = False

    def plot_channels_separately(channel_list):
        """One plot per variable per channel (no flavor merging). Combines over years
        when combine_years, else per-year. Used for --channel both and the extra
        per-channel plots produced by --channel all."""
        made = False
        if combine_years:
            lumi_text, lumi_unc = lumi_for()
            for channel in channel_list:
                units = build_units([channel])
                if not units:
                    print("ERROR: no year directories found for combined plot ({0})".format(channel))
                    continue
                print("\nChannel {0}: combining years {1}".format(channel, ", ".join(y for (y, _, _, _) in units)))
                for key in selected:
                    if key in COMBINED_ONLY_VARS:
                        continue
                    if plot_variable_combined(opts, VARIABLES[key], units, year_label, channel, lumi_text, lumi_unc):
                        made = True
        else:
            for channel in channel_list:
                if opts.nominal_dir:
                    nominal_dir = opts.nominal_dir.format(year=opts.year, channel=channel)
                else:
                    nominal_dir = default_nominal_dir(opts, channel)
                if not os.path.isdir(nominal_dir):
                    print("ERROR: nominal directory does not exist: {0}".format(nominal_dir))
                    continue
                syst_groups = discover_syst_dirs(nominal_dir, opts, channel)
                print("\nChannel {0}: found {1} systematic group(s)".format(channel, len(syst_groups)))
                for label, sides in syst_groups.items():
                    print("  {0:<35s} up={1} down={2}".format(
                        label, "yes" if "up" in sides else "no", "yes" if "down" in sides else "no"))
                for key in selected:
                    if key in COMBINED_ONLY_VARS:
                        continue
                    if plot_variable(opts, channel, VARIABLES[key], nominal_dir, syst_groups):
                        made = True
        return made

    if merge_channels:
        # one plot per variable, summing all years x both channels (datacard-style)
        units = build_units(["muon", "electron"])
        if not units:
            print("ERROR: no (year, channel) directories found for combined-channel plot")
            return 1
        lumi_text, lumi_unc = lumi_for()
        print("\nMerging units: {0}".format(", ".join("{0}/{1}".format(y, c) for (y, c, _, _) in units)))
        for key in selected:
            # Skip the flavor-merged lepton (lepton_pt/eta): we keep muon_pt/electron_pt and
            # muon_eta/electron_eta separate instead. Each common variable is merged over all
            # years AND both channels; each lepton variable is restricted to its own channel
            # (via applies_to_channel) and so is combined over years only.
            if key in COMBINED_ONLY_VARS:
                continue
            made = plot_variable_combined(opts, VARIABLES[key], units, year_label, "combined", lumi_text, lumi_unc)
            made_any = made_any or made
        # --channel all: ALSO produce the per-channel muon-only and electron-only plots.
        if opts.channel == "all":
            made_any = plot_channels_separately(["muon", "electron"]) or made_any

    else:
        channels = ["muon", "electron"] if opts.channel == "both" else [opts.channel]
        made_any = plot_channels_separately(channels) or made_any

    if not made_any:
        print("No plots were produced.")
        return 1
    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
