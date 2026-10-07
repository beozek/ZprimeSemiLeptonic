#!/usr/bin/env python
"""
Gen-level A_C with two top definitions on the SAME events:
  status-22 tops (hard process, what TTbarGen / our NoAC templates use)
  last-copy tops (after ISR recoil / FSR, before decay; t = W + b; TopAnalysis kStable)

Reads the gen_lastcopy/migration_22_vs_last* matrices from Preselection outputs produced with
FillLastCopyGenHists=true (ZprimePreselectionModule), sums them over files and years, and
prints A_C per gen-m_tt bin for both definitions with the paired (same-event) uncertainty
on the difference.

--weights:
  nnloewk  (default) genWeight * lumi * EWK * NNLO-QCD, the same per-event correction as the
           _lumiscaled_nnloewk gen folders (evaluated on the status-22 tops) applied to both
           definitions; needs ApplyNNLOEWKToGenHists=true in the preselection
  nnloewk_nocap  same, but NNLO-QCD without the pT > 500 GeV freeze (needs FillNNLONoCapGenHists=true)
  nnloewk_ownkin        status-22 result with the module weights vs last-copy result with EWK and
                        NNLO-QCD evaluated on the last-copy (W+b) tops, i.e. the dilepton TopAnalysis
                        recipe (needs FillLastCopyOwnKinGenHists=true)
  nnloewk_nocap_ownkin  same with the uncapped NNLO-QCD (needs FillNNLONoCapGenHists=true as well)
  lumi     genWeight * lumi, no higher-order correction

For the *_ownkin modes the two definitions carry different weights: A_22 comes from the
status-22 matrix (w_22), A_last from the *_ownkin matrix (w_last) and the per-cell covariance
sum(w_22 * w_last) from the *_ownkin_cross matrix.

Migration matrix convention (set in ZprimePreselectionModule.cxx):
  x cell = 2*mttbin_22   + (xi_22   > 0)
  y cell = 2*mttbin_last + (xi_last > 0)
Each cell holds sum(w) and sum(w^2) of disjoint event sets, so cells are statistically
independent and any function of them gets its error by linear propagation over cells.

Usage on the NAF (after running the preselection with FillLastCopyGenHists=true):
  python compare_lastcopy_genac.py --workdir-suffix ttbar_NNLOcorrectionatgen                  # EWK*NNLO corrected
  python compare_lastcopy_genac.py --workdir-suffix ttbar_NNLOcorrectionatgen --weights lumi   # genWeight*lumi only
  python compare_lastcopy_genac.py --files /path/to/uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_*.root
  python compare_lastcopy_genac.py --files /data/dust/group/cms/zprime-uhh/Presel_UL16_templatemethod/UL*/workdir_Preselection_*_genchecks2/uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_*.root
--weights nnloewk_ownkin
--weights nnloewk_nocap_ownkin
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
HIST_BASE = "gen_lastcopy/migration_22_vs_last"
# --weights -> (matrix for A_22, matrix for A_last, per-cell sum(w_22*w_last) or None if same weight)
# (suffixes set in ZprimePreselectionModule)
WEIGHT_MATRICES = {
    "nnloewk": (HIST_BASE + "_nnloewk", HIST_BASE + "_nnloewk", None),
    "nnloewk_nocap": (HIST_BASE + "_nnloewk_nocap", HIST_BASE + "_nnloewk_nocap", None),
    "nnloewk_ownkin": (HIST_BASE + "_nnloewk", HIST_BASE + "_nnloewk_ownkin", HIST_BASE + "_nnloewk_ownkin_cross"),
    "nnloewk_nocap_ownkin": (HIST_BASE + "_nnloewk_nocap", HIST_BASE + "_nnloewk_nocap_ownkin", HIST_BASE + "_nnloewk_nocap_ownkin_cross"),
    "lumi":    (HIST_BASE, HIST_BASE, None),
}
MTT_LABELS = ["0-500", "500-750", "750-1000", "1000-1500", ">1500"]


# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------

def build_workdir(base, year, workdir_suffix):
    return os.path.join(
        base,
        "Presel_{0}_templatemethod".format(year),
        "workdir_Preselection_{0}_templatemethod_{1}".format(year, workdir_suffix),
    )


def collect_files(args):
    if args.files:
        files = []
        for pattern in args.files:
            files.extend(glob.glob(pattern))
        return sorted(set(files))
    files = []
    for year in args.years:
        pattern = os.path.join(build_workdir(args.base, year, args.workdir_suffix), args.sample_glob)
        matched = sorted(glob.glob(pattern))
        if not matched:
            print("[WARN] No files matched: {0}".format(pattern), file=sys.stderr)
        else:
            print("[INFO] {0}: {1} files".format(year, len(matched)))
        files.extend(matched)
    return sorted(set(files))


def use_uproot():
    """uproot >= 4 only; the python2 uproot 3 in CMSSW is fragile, PyROOT is used there instead."""
    try:
        import uproot
        return int(uproot.__version__.split(".")[0]) >= 4
    except Exception:
        return False


def read_matrix_root(path, hist_path):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gErrorIgnoreLevel = ROOT.kError
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie() or f.TestBit(ROOT.TFile.kRecovered):
        raise IOError("cannot open (empty/truncated?)")
    h = f.Get(hist_path)
    if not h:
        f.Close()
        return None
    n = h.GetNbinsX()
    content = [[h.GetBinContent(i + 1, j + 1) for j in range(n)] for i in range(n)]
    sumw2 = [[h.GetBinError(i + 1, j + 1) ** 2 for j in range(n)] for i in range(n)]
    f.Close()
    return content, sumw2


def read_matrix_uproot(path, hist_path):
    import uproot
    f = uproot.open(path)
    try:
        h = f[hist_path]
    except KeyError:
        return None
    return [list(row) for row in h.values()], [list(row) for row in h.variances()]


def read_matrix(path, hist_path):
    # uproot >= 4 if available (a broken PyROOT can segfault on import), else PyROOT
    if use_uproot():
        return read_matrix_uproot(path, hist_path)
    return read_matrix_root(path, hist_path)


def sum_matrices(files, hist_path):
    total_c, total_w2, used = None, None, 0
    for path in files:
        try:
            m = read_matrix(path, hist_path)
        except Exception as exc:
            print("[WARN] skipping {0} ({1:.0f} bytes): {2}".format(path, os.path.getsize(path) if os.path.exists(path) else -1, exc), file=sys.stderr)
            continue
        if m is None:
            print("[WARN] {0} has no {1} (run with FillLastCopyGenHists=true; _nnloewk also needs ApplyNNLOEWKToGenHists=true)".format(path, hist_path), file=sys.stderr)
            continue
        c, w2 = m
        if total_c is None:
            total_c = [row[:] for row in c]
            total_w2 = [row[:] for row in w2]
        else:
            for i in range(len(c)):
                for j in range(len(c)):
                    total_c[i][j] += c[i][j]
                    total_w2[i][j] += w2[i][j]
        used += 1
    return total_c, total_w2, used


# ---------------------------------------------------------------------------
# A_C with paired errors
# ---------------------------------------------------------------------------

def asymmetry(fwd_cells, bwd_cells, content):
    """A = (F - B)/(F + B) and dA/dc for every contributing cell."""
    F = sum(content[i][j] for i, j in fwd_cells)
    B = sum(content[i][j] for i, j in bwd_cells)
    N = F + B
    if N == 0:
        return float("nan"), {}
    grad = {}
    for cell in fwd_cells:
        grad[cell] = grad.get(cell, 0.) + 2. * B / N ** 2
    for cell in bwd_cells:
        grad[cell] = grad.get(cell, 0.) - 2. * F / N ** 2
    return (F - B) / N, grad


def propagate(grad, sumw2):
    return math.sqrt(sum(g * g * sumw2[i][j] for (i, j), g in grad.items()))


def paired_error(g22, v22, glast, vlast, cross):
    """Error on A_last - A_22: sum over independent cells of the 2x2 (w_22, w_last) covariance."""
    var = 0.
    for (i, j) in set(g22) | set(glast):
        a, b = g22.get((i, j), 0.), glast.get((i, j), 0.)
        var += a * a * v22[i][j] + b * b * vlast[i][j] - 2. * a * b * cross[i][j]
    return math.sqrt(max(var, 0.))


def analyse(content, sumw2, content_last=None, sumw2_last=None, cross=None):
    """Per gen-m_tt bin (plus inclusive): A_22, A_last, their errors and the paired error on the difference.

    content/sumw2: matrix whose x-marginals give A_22. content_last/sumw2_last: matrix whose
    y-marginals give A_last (default: the same matrix). cross: per-cell sum(w_22 * w_last)
    (default: sumw2, which is exact when both definitions carry the same event weight).
    """
    if content_last is None:
        content_last, sumw2_last, cross = content, sumw2, sumw2
    ncell = len(content)
    nbin = ncell // 2
    cells = [(i, j) for i in range(ncell) for j in range(ncell)]
    rows = []
    for b in list(range(nbin)) + [None]:
        in_22 = (lambda i: True) if b is None else (lambda i, b=b: i // 2 == b)
        a22, g22 = asymmetry([(i, j) for i, j in cells if in_22(i) and i % 2 == 1],
                             [(i, j) for i, j in cells if in_22(i) and i % 2 == 0], content)
        alast, glast = asymmetry([(i, j) for i, j in cells if in_22(j) and j % 2 == 1],
                                 [(i, j) for i, j in cells if in_22(j) and j % 2 == 0], content_last)
        n22 = sum(content[i][j] for i, j in cells if in_22(i))
        if b is None:
            leave, flip = 0., sum(content[i][j] for i, j in cells if i % 2 != j % 2)
        else:
            leave = sum(content[i][j] for i, j in cells if in_22(i) and j // 2 != b)
            flip = sum(content[i][j] for i, j in cells if in_22(i) and j // 2 == b and i % 2 != j % 2)
        rows.append({
            "bin": "inclusive" if b is None else MTT_LABELS[b] if b < len(MTT_LABELS) else str(b),
            "A_22": a22, "err_22": propagate(g22, sumw2),
            "A_last": alast, "err_last": propagate(glast, sumw2_last),
            "diff": alast - a22, "err_diff": paired_error(g22, sumw2, glast, sumw2_last, cross),
            "frac_leave": leave / n22 if n22 else float("nan"),
            "frac_flip": flip / n22 if n22 else float("nan"),
        })
    return rows


def print_table(rows):
    print()
    print("{:<11}{:>22}{:>22}{:>24}{:>9}{:>9}{:>9}".format(
        "gen m_tt", "A_C status-22", "A_C last copy", "last - 22 (paired)", "rel.", "leave", "flip"))
    for r in rows:
        rel = r["diff"] / r["A_22"] if r["A_22"] else float("nan")
        print("{:<11}{:>12.5f} +- {:.5f}{:>12.5f} +- {:.5f}{:>+13.5f} +- {:.5f}{:>+8.1f}%{:>8.1f}%{:>8.1f}%".format(
            r["bin"], r["A_22"], r["err_22"], r["A_last"], r["err_last"], r["diff"], r["err_diff"],
            100. * rel, 100. * r["frac_leave"], 100. * r["frac_flip"]))
    print()
    print("leave = fraction of status-22 events in the bin whose last-copy m_tt falls in another bin")
    print("flip  = fraction that stay in the bin but change the sign of xi")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--files", nargs="*", help="Explicit input files/globs (overrides --base/--years/...).")
    parser.add_argument("--base", default=DEFAULT_BASE)
    parser.add_argument("--years", nargs="*", default=DEFAULT_YEARS)
    parser.add_argument("--workdir-suffix", default=DEFAULT_WORKDIR_SUFFIX)
    parser.add_argument("--sample-glob", default=DEFAULT_SAMPLE_GLOB)
    parser.add_argument("--weights", choices=sorted(WEIGHT_MATRICES), default="nnloewk",
                        help="see the list above. Default: %(default)s.")
    args = parser.parse_args()

    files = collect_files(args)
    if not files:
        sys.exit("No input files found.")
    path_22, path_last, path_cross = WEIGHT_MATRICES[args.weights]
    print("[INFO] weights: {0}".format(args.weights))
    content, sumw2, used = sum_matrices(files, path_22)
    if content is None:
        sys.exit("None of the {0} files contain {1}.".format(len(files), path_22))
    print("[INFO] {0}: summed {1}/{2} files".format(path_22, used, len(files)))
    if path_cross is None:
        print_table(analyse(content, sumw2))
        return
    content_last, sumw2_last, used_last = sum_matrices(files, path_last)
    cross, _, used_cross = sum_matrices(files, path_cross)
    if content_last is None or cross is None or not (used == used_last == used_cross):
        sys.exit("{0} / {1} missing or not in the same files as {2}.".format(path_last, path_cross, path_22))
    print_table(analyse(content, sumw2, content_last, sumw2_last, cross))


if __name__ == "__main__":
    main()
