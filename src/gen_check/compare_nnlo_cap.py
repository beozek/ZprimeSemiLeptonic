#!/usr/bin/env python
"""
Gen-level A_C per gen-m_tt bin with the NNLO-QCD SF frozen above pT = 500 GeV (standard,
_lumiscaled_nnloewk folders) and without the freeze (_lumiscaled_nnloewk_nocap folders, from
ZprimePreselectionModule with FillNNLONoCapGenHists=true). Both use the same events, the same
status-22 tops and the same EWK correction; only the NNLO-QCD weight above 500 GeV differs.

Error on the difference: the two weights of an event differ by the factor
r = w_nocap / w_capped in [SF(inf)/SF(500), 1] = [0.9805, 1], so (delta method) the change of
A_C has a statistical uncertainty of at most about (1/r_min - 1) * sigma(A_C) = 0.020 * sigma(A_C);
in practice smaller, since r is close to its bin average.

Usage on the NAF:
  python compare_nnlo_cap.py --workdir-suffix ttbar_NNLOcorrectionatgen
  python compare_nnlo_cap.py --files /path/to/uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_*.root
"""

from __future__ import division, print_function

import argparse
import glob
import math
import os
import sys

DEFAULT_BASE = "/data/dust/group/cms/zprime-uhh"
DEFAULT_YEARS = ["UL16preVFP", "UL16postVFP", "UL17", "UL18"]
DEFAULT_WORKDIR_SUFFIX = "ttbar_NNLOcorrectionatgen"
DEFAULT_SAMPLE_GLOB = "uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_*.root"
TAGS = [("inclusive", "inclusive"), ("0_500", "0-500"), ("500_750", "500-750"),
        ("750_1000", "750-1000"), ("1000_1500", "1000-1500"), ("1500Inf", ">1500")]
FOLDER = "mtt_gen_{0}_lumiscaled_nnloewk{1}_General/{2}"
SF3 = lambda pt: 0.17478303752138802 * math.exp(-0.0045665365705007399 * pt) + 0.89410719420857065
MAX_DELTA = SF3(500.) / SF3(1e9) - 1.   # 1/r_min - 1


def collect_files(args):
    if args.files:
        return sorted(set(f for pattern in args.files for f in glob.glob(pattern)))
    files = []
    for year in args.years:
        workdir = os.path.join(args.base, "Presel_{0}_templatemethod".format(year),
                               "workdir_Preselection_{0}_templatemethod_{1}".format(year, args.workdir_suffix))
        matched = sorted(glob.glob(os.path.join(workdir, args.sample_glob)))
        print("[INFO] {0}: {1} files".format(year, len(matched)) if matched else
              "[WARN] no files in {0}".format(workdir), file=sys.stderr if not matched else sys.stdout)
        files.extend(matched)
    return sorted(set(files))


def use_uproot():
    """uproot >= 4 only; the python2 uproot 3 in CMSSW is fragile, PyROOT is used there instead."""
    try:
        import uproot
        return int(uproot.__version__.split(".")[0]) >= 4
    except Exception:
        return False


def read_hist(path, name):
    """(values, variances) of a TH1 with uproot >= 4, or PyROOT; None if missing; raises if unreadable."""
    if use_uproot():
        import uproot
        f = uproot.open(path)
        try:
            h = f[name]
        except KeyError:
            return None
        return list(h.values()), list(h.variances())
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gErrorIgnoreLevel = ROOT.kError
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie() or f.TestBit(ROOT.TFile.kRecovered):
        raise IOError("cannot open (empty/truncated?)")
    h = f.Get(name)
    if not h:
        f.Close()
        return None
    n = h.GetNbinsX()
    out = [h.GetBinContent(i) for i in range(1, n + 1)], [h.GetBinError(i) ** 2 for i in range(1, n + 1)]
    f.Close()
    return out


def summed(files, name):
    tot_v, tot_e, used = None, None, 0
    for path in files:
        try:
            h = read_hist(path, name)
        except Exception as exc:
            print("[WARN] skipping {0} ({1:.0f} bytes): {2}".format(path, os.path.getsize(path) if os.path.exists(path) else -1, exc), file=sys.stderr)
            continue
        if h is None:
            print("[WARN] {0} has no {1}".format(path, name), file=sys.stderr)
            continue
        v, e = h
        if tot_v is None:
            tot_v, tot_e = list(v), list(e)
        else:
            tot_v = [a + b for a, b in zip(tot_v, v)]
            tot_e = [a + b for a, b in zip(tot_e, e)]
        used += 1
    return tot_v, tot_e, used


def asym(vals, var):
    """A_C from a symmetric binning in xi (first half backward, second half forward), with its error."""
    n = len(vals)
    B, F = sum(vals[:n // 2]), sum(vals[n // 2:])
    vB, vF = sum(var[:n // 2]), sum(var[n // 2:])
    N = F + B
    if N <= 0:
        return float("nan"), float("nan"), N
    return (F - B) / N, 2. * math.sqrt(F * F * vB + B * B * vF) / N ** 2, N


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--files", nargs="*", help="Explicit input files/globs (overrides --base/--years/...).")
    parser.add_argument("--base", default=DEFAULT_BASE)
    parser.add_argument("--years", nargs="*", default=DEFAULT_YEARS)
    parser.add_argument("--workdir-suffix", default=DEFAULT_WORKDIR_SUFFIX)
    parser.add_argument("--sample-glob", default=DEFAULT_SAMPLE_GLOB)
    parser.add_argument("--hist-name", default="DeltaY_xi_gen_300", help="even-binned xi histogram (default %(default)s)")
    args = parser.parse_args()

    files = collect_files(args)
    if not files:
        sys.exit("No input files found.")

    print()
    print("{:<11}{:>22}{:>22}{:>14}{:>12}{:>13}".format(
        "gen m_tt", "A_C capped", "A_C no cap", "nocap-capped", "stat bound", "yield ratio"))
    for tag, label in TAGS:
        cv, ce, n1 = summed(files, FOLDER.format(tag, "", args.hist_name))
        nv, ne, n2 = summed(files, FOLDER.format(tag, "_nocap", args.hist_name))
        if cv is None or nv is None or n1 != n2:
            sys.exit("Missing {0} (run with FillNNLONoCapGenHists=true?) or not in the same files.".format(
                FOLDER.format(tag, "_nocap", args.hist_name)))
        a_c, s_c, y_c = asym(cv, ce)
        a_n, s_n, y_n = asym(nv, ne)
        print("{:<11}{:>12.5f} +- {:.5f}{:>12.5f} +- {:.5f}{:>+14.6f}{:>12.6f}{:>13.4f}".format(
            label, a_c, s_c, a_n, s_n, a_n - a_c, MAX_DELTA * s_c, y_n / y_c if y_c else float("nan")))
    print()
    print("[INFO] summed {0} files".format(len(files)))
    print("stat bound  = approx. upper limit on the statistical uncertainty of (nocap - capped): same events, "
          "weights differ by <= {0:.1%}".format(1. - SF3(1e9) / SF3(500.)))
    print("yield ratio = sum(w_nocap) / sum(w_capped)")


if __name__ == "__main__":
    main()
