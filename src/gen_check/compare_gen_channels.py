#!/usr/bin/env python
"""
Gen-level A_C per ttbar decay channel and gen-m_tt bin, from the gen_channels/ matrices written by
ZprimePreselectionModule with FillGenChannelHists=true (tops taken directly from the GenParticles
collection, all events, no tau veto), for the status-22 tops and for the last-copy (W+b) tops.

Matrix convention: x = decay channel (0 had, 1 ehad, 2 muhad, 3 tauhad, 4 ee, 5 mumu, 6 tautau,
7 emu, 8 etau, 9 mutau, 10 unclassified), y = 2*mttbin + (xi > 0), weights genWeight*lumi
(ac_<def>) or genWeight*lumi*EWK*NNLO (ac_<def>_nnloewk).

Cross-check: the sum of the non-tau channels must reproduce the gen_lastcopy/migration_22_vs_last
matrix (same events, same weights), if that matrix is present in the files.

  python compare_gen_channels.py --files "/data/dust/group/cms/zprime-uhh/Presel_UL16_templatemethod/*/workdir_Preselection_UL*_genchecks2/uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL*_*.root"
  python compare_gen_channels.py --files "..." --weights lumi
"""

from __future__ import division, print_function

import argparse
import glob
import math
import os
import sys

CHANNELS = ["had", "ehad", "muhad", "tauhad", "ee", "mumu", "tautau", "emu", "etau", "mutau", "unclassified"]
GROUPS = [
    ("e/mu + jets (1,2)", [1, 2]),
    ("tau + jets (3)", [3]),
    ("dilepton, no tau (4,5,7)", [4, 5, 7]),
    ("dilepton, with tau (6,8,9)", [6, 8, 9]),
    ("all without tau (0,1,2,4,5,7)", [0, 1, 2, 4, 5, 7]),
    ("all with tau (3,6,8,9)", [3, 6, 8, 9]),
    ("all channels", list(range(11))),
]
BINS = ["inclusive", "0-500", "500-750", "750-1000", "1000-1500", ">1500"]
NBIN = 5


def use_uproot():
    try:
        import uproot
        return int(uproot.__version__.split(".")[0]) >= 4
    except Exception:
        return False


def read_th2(path, name):
    """(values, variances) as nested lists [x][y]; None if missing; raises if unreadable."""
    if use_uproot():
        import uproot
        f = uproot.open(path)
        try:
            h = f[name]
        except KeyError:
            return None
        return [list(r) for r in h.values()], [list(r) for r in h.variances()]
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
    nx, ny = h.GetNbinsX(), h.GetNbinsY()
    out = ([[h.GetBinContent(i + 1, j + 1) for j in range(ny)] for i in range(nx)],
           [[h.GetBinError(i + 1, j + 1) ** 2 for j in range(ny)] for i in range(nx)])
    f.Close()
    return out


def summed(files, name):
    tot_v, tot_e, used = None, None, 0
    for path in files:
        try:
            h = read_th2(path, name)
        except Exception as exc:
            print("[WARN] skipping {0}: {1}".format(path, exc), file=sys.stderr)
            continue
        if h is None:
            continue
        v, e = h
        if tot_v is None:
            tot_v, tot_e = [r[:] for r in v], [r[:] for r in e]
        else:
            for i in range(len(v)):
                for j in range(len(v[i])):
                    tot_v[i][j] += v[i][j]
                    tot_e[i][j] += e[i][j]
        used += 1
    return tot_v, tot_e, used


def ac_from_cells(v, e, channels, b):
    """A_C and its error for a set of channels in m_tt bin b (None = inclusive)."""
    F = B = vF = vB = 0.
    for c in channels:
        for cell in range(len(v[c])):
            if b is not None and cell // 2 != b:
                continue
            if cell % 2 == 1:
                F += v[c][cell]
                vF += e[c][cell]
            else:
                B += v[c][cell]
                vB += e[c][cell]
    N = F + B
    if N <= 0:
        return float("nan"), float("nan"), N
    return (F - B) / N, 2. * math.sqrt(F * F * vB + B * B * vF) / N ** 2, N


def print_table(title, v, e, rows):
    print("\n" + title)
    print("{:<32}{:>10}".format("channel", "yield") + "".join("{:>18}".format(b) for b in BINS))
    for label, chans in rows:
        _, _, n = ac_from_cells(v, e, chans, None)
        if n <= 0:
            continue
        cells = ["{:>+9.5f}+-{:.5f}".format(*ac_from_cells(v, e, chans, b)[:2]) for b in [None] + list(range(NBIN))]
        print("{:<32}{:>10.3g}".format(label, n) + "".join("{:>18}".format(c) for c in cells))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--files", nargs="+", required=True, help="input files/globs")
    parser.add_argument("--weights", choices=["nnloewk", "lumi"], default="nnloewk")
    args = parser.parse_args()
    files = sorted(set(f for pat in args.files for f in glob.glob(pat)))
    if not files:
        sys.exit("No input files found.")
    sfx = "_nnloewk" if args.weights == "nnloewk" else ""
    v22, e22, n22 = summed(files, "gen_channels/ac_22" + sfx)
    vl, el, nl = summed(files, "gen_channels/ac_last" + sfx)
    if v22 is None or vl is None:
        sys.exit("No gen_channels/ac_22{0} / ac_last{0} found (run with FillGenChannelHists=true?).".format(sfx))
    print("[INFO] weights: {0}; summed {1}/{2} files".format(args.weights, n22, len(files)))

    rows = [(CHANNELS[i], [i]) for i in range(11)] + GROUPS
    print_table("A_C, status-22 tops (direct from GenParticles)", v22, e22, rows)
    print_table("A_C, last-copy tops (W+b, direct from GenParticles)", vl, el, rows)

    # cross-check against the TTbarGen-based gen_lastcopy matrix (same events: non-tau channels)
    mig, _, nmig = summed(files, "gen_lastcopy/migration_22_vs_last" + sfx)
    if mig is not None and nmig == n22:
        notau = [0, 1, 2, 4, 5, 7]
        ok = True
        for b in range(NBIN):
            a_ch, _, n_ch = ac_from_cells(v22, e22, notau, b)
            F = sum(mig[2 * b + 1][j] for j in range(len(mig))); B = sum(mig[2 * b][j] for j in range(len(mig)))
            a_lc, n_lc = (F - B) / (F + B), F + B
            if abs(n_ch - n_lc) > 1e-6 * n_lc or abs(a_ch - a_lc) > 1e-9:
                ok = False
                print("[CHECK] bin {0}: channels {1:.6g} / A={2:.6f}  vs gen_lastcopy {3:.6g} / A={4:.6f}".format(BINS[b + 1], n_ch, a_ch, n_lc, a_lc))
        print("[CHECK] non-tau channels vs gen_lastcopy (TTbarGen) matrix: " + ("IDENTICAL" if ok else "DIFFER (see above)"))


if __name__ == "__main__":
    main()
