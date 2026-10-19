#!/usr/bin/env python3
"""
Pooled-null z-score test: standardizes each family/domain's real delta-logL
and its 10 simulated null delta-logLs onto a common per-family scale, then
pools the *null-side* standardized values across every family into one
large empirical reference distribution. Each family's real z-score is then
compared against that pooled reference, rather than against just its own
10 replicates.

Why: the existing permutation p-value (fisher_combined_test.py) ranks the
real delta-logL against only 10 null replicates, so the achievable
resolution is floored at 1/11 ~= 0.091 -- most families (86.8% at last
count) hit that floor and become indistinguishable from each other, even
though their real values differ hugely in magnitude. Standardizing removes
cross-family scale differences (alignment length, taxa count, branch
lengths all push delta-logL up or down regardless of monophyly), and
pooling the null side across ~5000 families turns a 10-point reference
into a ~50,000-point one -- resolution limited by the dataset size, not
by N_REPLICATES.

Per family/domain:
  real_z = (real_delta - mean(null_deltas)) / sd(null_deltas)
    (uses all N null replicates -- the real value isn't part of the null
    set being standardized, so no leave-one-out needed here)

  null z_i, for each of the N null replicates, leave-one-out:
    z_i = (null_i - mean(null without i)) / sd(null without i)
    (excludes replicate i from its own standardization parameters so it
    isn't mechanically shrunk toward 0 by partly standardizing against
    itself, which would make the pooled reference artificially less
    extreme in the tails than it really is)

Pooling (per domain slot, since Domain1/Domain2 may have different
typical null shapes -- parent-gene-1 vs parent-gene-2 alignments, taxa
counts, etc.): every leave-one-out null z from every family, for that
domain, forms the reference. A family's real_z is compared against ITS
OWN domain's pooled reference, never mixed across domains for the
per-family p-value (only the summary/diagnostic numbers mix both slots).

    p = (n_pooled_null_z_as_extreme + 1) / (N_pooled + 1)

Caveat (reported, not enforced): this pooling is only valid if the
*standardized* null shapes are reasonably similar across families --
not that any one family's null is normal, just that different families'
standardized nulls aren't wildly different shapes. Per-domain mean/sd/
skewness of the pooled distribution are printed so this can be eyeballed;
a pooled mean far from 0, sd far from 1, or large |skewness| would be a
reason to distrust the pooled p-values.

Usage:
  python3 src/pooled_null_z_test.py --report data/reports/failure_report_TIMESTAMP.tsv \
      [--out-tsv data/reports/pooled_null_z_test.tsv] [--email you@example.com]
"""
import argparse
import csv
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

from fisher_combined_test import fishers_combined, send_email

# A handful of families have null delta-logL replicates that are all
# ~identical (near-zero sd, e.g. literally [0,0,...,0] or ~1e-12) --
# numerically degenerate, not genuinely low-variance biology. Detected by
# leave-one-out standardizing a family's OWN 10 null replicates against
# each other: if any of THOSE null-vs-null z's exceeds this bound, that
# family's null sd is degenerate and the whole family/domain is excluded
# (from the pool AND from getting a real_z) -- its own real_z is NEVER
# used to decide this, since a huge real_z can also mean a genuine,
# enormous real effect (confirmed: most large real_z's have perfectly
# normal, non-degenerate null sd -- e.g. F1194 Domain1, real_delta=2293,
# null sd=1.7, real_z=1348, legitimate). There is a clean gap in the
# null-vs-null z distribution between well-behaved values (99.9th
# percentile ~20) and degenerate ones (next value up: ~1958) -- this
# bound sits safely in that gap.
OUTLIER_Z_BOUND = 100


def is_degenerate_null(null_deltas):
    """True if this family's null replicates are ~identical (near-zero
    sd), detected via its own leave-one-out null-vs-null z's -- never via
    the real value, which can legitimately be far more extreme than the
    null without the null itself being degenerate."""
    n = len(null_deltas)
    for i in range(n):
        others = null_deltas[:i] + null_deltas[i + 1:]
        sd = statistics.stdev(others)
        if sd == 0:
            return True
        m = statistics.mean(others)
        z = (null_deltas[i] - m) / sd
        if abs(z) > OUTLIER_Z_BOUND:
            return True
    return False


UNRELIABLE_NEGATIVE_FRACTION = 0.5


def is_search_unreliable(null_deltas, threshold=UNRELIABLE_NEGATIVE_FRACTION):
    """True if >=threshold of this domain's null replicates have a
    NEGATIVE delta-logL -- the free (unconstrained) search did worse than
    the constrained tree it was supposed to beat. A search that can't
    reliably outperform a reachable starting point signals too little
    phylogenetic information in this domain's taxa count/alignment
    length/tree length combination to trust any comparison built on its
    null, independent of which way that comparison comes out. Threshold
    chosen from the dataset's own distribution: ~0.9% of domains clear
    0.5, concentrated almost entirely in large-taxa/short-alignment
    combinations -- a real, systematic pattern, not noise."""
    if not null_deltas:
        return False
    n_negative = sum(1 for v in null_deltas if v < 0)
    return (n_negative / len(null_deltas)) >= threshold


def collect_family_data(report_path):
    """Returns by_family: {family: {domain_slot: (real_delta, [null_deltas])}}"""
    by_family = defaultdict(dict)

    with open(report_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        domain_slots = sorted({
            m.group(1) for col in reader.fieldnames
            if (m := re.match(r'^Domain(\d+)_Real_DeltaLogL$', col))
        }, key=int)

        rep_cols = {}
        for d in domain_slots:
            rep_cols[d] = sorted(
                (int(m.group(1)), col) for col in reader.fieldnames
                if (m := re.match(rf'^Domain{d}_DeltaLogL_Rep(\d+)$', col))
            )

        for row in reader:
            fam = row["Family_ID"]
            for d in domain_slots:
                real_raw = row.get(f"Domain{d}_Real_DeltaLogL", "NA")
                if real_raw == "NA" or real_raw == "":
                    continue
                try:
                    real_delta = float(real_raw)
                except ValueError:
                    continue

                null_deltas = []
                for _, col in rep_cols[d]:
                    v = row.get(col, "NA")
                    if v != "NA" and v != "":
                        try:
                            null_deltas.append(float(v))
                        except ValueError:
                            pass

                if len(null_deltas) < 3:
                    continue  # need at least 2 "others" for leave-one-out sd

                by_family[fam][d] = (real_delta, null_deltas)

    return by_family, domain_slots


def leave_one_out_z(values):
    """Standardize each value against the mean/sd of the OTHER values only.
    Caller is expected to have already excluded degenerate families via
    is_degenerate_null -- no bound-checking here."""
    n = len(values)
    zs = []
    for i in range(n):
        others = values[:i] + values[i + 1:]
        sd = statistics.stdev(others)
        if sd == 0:
            continue
        m = statistics.mean(others)
        zs.append((values[i] - m) / sd)
    return zs


def real_z_score(real_delta, null_deltas):
    """No bound-checking here -- a large real_z can be a genuine, enormous
    real effect. Degeneracy is screened separately via is_degenerate_null,
    using only the family's own null replicates, before this is called."""
    sd = statistics.stdev(null_deltas)
    if sd == 0:
        return None
    m = statistics.mean(null_deltas)
    return (real_delta - m) / sd


def empirical_p_upper(z, pooled_null_z):
    n_as_extreme = sum(1 for v in pooled_null_z if v >= z)
    return (n_as_extreme + 1) / (len(pooled_null_z) + 1)


def bh_reject_count(pvalues, alpha=0.05):
    """Benjamini-Hochberg step-up procedure: largest k such that the k-th
    smallest p-value is <= (k/m)*alpha. Returns (k, m, p_cutoff)."""
    m = len(pvalues)
    ps = sorted(pvalues)
    largest_k = 0
    p_cutoff = None
    for k, p in enumerate(ps, start=1):
        if p <= (k / m) * alpha:
            largest_k = k
            p_cutoff = p
    return largest_k, m, p_cutoff


def bh_qvalues(pvalues):
    """Benjamini-Hochberg adjusted p-values (q-values), one per input
    p-value, in the SAME order as the input list (not sorted). Standard
    step-up definition: q_(m) = p_(m); q_(k) = min(q_(k+1), p_(k)*m/k) for
    k = m-1 .. 1, enforced monotone by the running minimum. A q-value is
    the smallest FDR at which that test would be called significant --
    compare directly against alpha (reject if q <= alpha), no further
    correction needed."""
    m = len(pvalues)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: pvalues[i])
    q = [0.0] * m
    running_min = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        k = rank + 1  # 1-indexed rank, smallest p-value = rank 1
        running_min = min(running_min, pvalues[i] * m / k)
        q[i] = running_min
    return q


def pooled_diagnostics(pooled_null_z):
    n = len(pooled_null_z)
    if n < 2:
        return None
    m = statistics.mean(pooled_null_z)
    sd = statistics.stdev(pooled_null_z)
    if sd > 0:
        skew = (sum((v - m) ** 3 for v in pooled_null_z) / n) / (sd ** 3)
    else:
        skew = float("nan")
    return {"n": n, "mean": m, "sd": sd, "skew": skew}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report", required=True, help="generate_failure_report.py output TSV")
    ap.add_argument("--out-tsv", default=None, help="Optional: write per-family, per-domain detail to this file")
    ap.add_argument("--email", default=None, help="If given, email the summary to this address")
    args = ap.parse_args()

    by_family, domain_slots = collect_family_data(args.report)
    if not by_family:
        print("No families with both a real delta-logL and >=3 null replicates found.", file=sys.stderr)
        sys.exit(1)

    # Pass 1: screen for degenerate null variance (own null-vs-null z's
    # only), THEN leave-one-out null z's pooled per domain slot, THEN
    # real_z per family/domain -- a degenerate family contributes neither.
    pooled_null_z = defaultdict(list)
    real_z = defaultdict(dict)  # {domain: {family: z}}
    unreliable = defaultdict(dict)  # {domain: {family: bool}}
    n_degenerate = 0
    n_fam_dom_total = 0

    for fam, domains in by_family.items():
        for d, (real_delta, null_deltas) in domains.items():
            n_fam_dom_total += 1
            if is_degenerate_null(null_deltas):
                n_degenerate += 1
                continue
            pooled_null_z[d].extend(leave_one_out_z(null_deltas))
            rz = real_z_score(real_delta, null_deltas)
            if rz is not None:
                real_z[d][fam] = rz
                unreliable[d][fam] = is_search_unreliable(null_deltas)

    pooled_null_z_all = [v for d in domain_slots for v in pooled_null_z[d]]

    # Pass 2: per-family, per-domain empirical p-value against its own domain's pooled null.
    pooled_p = defaultdict(dict)  # {domain: {family: p}}
    for d in domain_slots:
        ref = pooled_null_z[d]
        for fam, z in real_z[d].items():
            pooled_p[d][fam] = empirical_p_upper(z, ref)

    lines = [
        "Pooled-null z-score test: real delta-logL vs. leave-one-out-standardized,",
        "cross-family-pooled null delta-logL distribution",
        "",
        "Per family/domain: real_z = (real_delta - mean(null)) / sd(null) (own 10 replicates).",
        "Null reference: every family's 10 null replicates, each standardized",
        "leave-one-out against the other 9, pooled across all families for that",
        "domain slot. p = (n_pooled_null_z >= real_z + 1) / (N_pooled + 1).",
        f"Excluded as numerically degenerate null variance (own null-vs-null |z| > {OUTLIER_Z_BOUND}): "
        f"{n_degenerate}/{n_fam_dom_total} family-domain rows.",
        "",
    ]

    for d in sorted(domain_slots, key=int):
        diag = pooled_diagnostics(pooled_null_z[d])
        pvals = list(pooled_p[d].values())
        n_fam = len(pvals)
        lines.append(f"Domain{d}: {n_fam} families, pooled null reference size = {diag['n']}")
        lines.append(f"  Pooled null z diagnostics: mean={diag['mean']:.3f} sd={diag['sd']:.3f} skew={diag['skew']:.3f}"
                     f"  (expect mean~0, sd~1 if standardization behaved; |skew| large -> pooling less trustworthy)")
        min_p = 1.0 / (diag['n'] + 1)
        n_at_floor = sum(1 for p in pvals if p <= min_p)
        n_sig_05 = sum(1 for p in pvals if p <= 0.05)
        n_sig_bonf = sum(1 for p in pvals if p <= 0.05 / n_fam) if n_fam else 0
        lines.append(f"  Achievable minimum p at this pool size: {min_p:.2e}")
        lines.append(f"  Families at that floor: {n_at_floor}/{n_fam} ({100*n_at_floor/n_fam:.1f}%)")
        lines.append(f"  Families at p<=0.05 (uncorrected): {n_sig_05}/{n_fam} ({100*n_sig_05/n_fam:.1f}%)")
        lines.append(f"  Families at p<=0.05/{n_fam} (Bonferroni-corrected): {n_sig_bonf}/{n_fam} ({100*n_sig_bonf/n_fam:.1f}%)")
        bh_k, bh_m, bh_cutoff = bh_reject_count(pvals, 0.05)
        lines.append(f"  Families significant at FDR<=0.05 (Benjamini-Hochberg): {bh_k}/{bh_m} ({100*bh_k/bh_m:.1f}%)"
                     f"  [BH p-value cutoff used: {bh_cutoff:.4g}]" if bh_k else
                     f"  Families significant at FDR<=0.05 (Benjamini-Hochberg): 0/{bh_m} (0.0%)")
        stat, log10_p, k = fishers_combined(pvals)
        lines.append(f"  Fisher's combined (using these finer p-values): stat={stat:.1f} df={2*k}"
                     f"  log10(p)={log10_p:.1f}")

        # Reliability split: >=50% of a family's own null replicates
        # failing to beat the constrained tree means that family's null
        # shouldn't be trusted regardless of which way its p comes out
        # (see is_search_unreliable). Reported separately rather than
        # folded into the headline number above.
        fams = list(pooled_p[d].keys())
        reliable_p = [pooled_p[d][fam] for fam in fams if not unreliable[d][fam]]
        unreliable_p = [pooled_p[d][fam] for fam in fams if unreliable[d][fam]]
        n_bh_rel, m_rel, _ = bh_reject_count(reliable_p, 0.05) if reliable_p else (0, 0, None)
        n_bh_unrel, m_unrel, _ = bh_reject_count(unreliable_p, 0.05) if unreliable_p else (0, 0, None)
        stat_rel, log10_p_rel, k_rel = fishers_combined(reliable_p)
        stat_unrel, log10_p_unrel, k_unrel = fishers_combined(unreliable_p)
        lines.append(f"  RELIABLE   (null search behaved): {m_rel} families; "
                     f"BH<=0.05: {n_bh_rel}/{m_rel} ({100*n_bh_rel/m_rel:.1f}%); "
                     f"Fisher's log10(p)={log10_p_rel:.1f}" if m_rel else
                     "  RELIABLE   (null search behaved): 0 families")
        lines.append(f"  UNRELIABLE (>=50% null reps worse than constrained): {m_unrel} families; "
                     f"BH<=0.05: {n_bh_unrel}/{m_unrel} ({100*n_bh_unrel/m_unrel:.1f}%); "
                     f"Fisher's log10(p)={log10_p_unrel:.1f}" if m_unrel else
                     "  UNRELIABLE (>=50% null reps worse than constrained): 0 families")
        lines.append("")

    diag_all = pooled_diagnostics(pooled_null_z_all)
    all_pvals = [p for d in domain_slots for p in pooled_p[d].values()]
    n_all = len(all_pvals)
    lines.append(f"POOLED across all domain slots: {n_all} family-domain tests, "
                 f"null reference size = {diag_all['n']}")
    lines.append(f"  Pooled null z diagnostics: mean={diag_all['mean']:.3f} sd={diag_all['sd']:.3f} "
                 f"skew={diag_all['skew']:.3f}")
    n_sig_05_all = sum(1 for p in all_pvals if p <= 0.05)
    n_sig_bonf_all = sum(1 for p in all_pvals if p <= 0.05 / n_all) if n_all else 0
    bh_k_all, bh_m_all, bh_cutoff_all = bh_reject_count(all_pvals, 0.05)
    lines.append(f"  Families at p<=0.05 (uncorrected): {n_sig_05_all}/{n_all} ({100*n_sig_05_all/n_all:.1f}%)")
    lines.append(f"  Families at p<=0.05/{n_all} (Bonferroni-corrected): {n_sig_bonf_all}/{n_all} ({100*n_sig_bonf_all/n_all:.1f}%)")
    lines.append(f"  Families significant at FDR<=0.05 (Benjamini-Hochberg): {bh_k_all}/{bh_m_all} ({100*bh_k_all/bh_m_all:.1f}%)")
    stat, log10_p, k = fishers_combined(all_pvals)
    lines.append(f"  Fisher's combined (using these finer p-values): stat={stat:.1f} df={2*k} "
                 f"log10(p)={log10_p:.1f}")

    body = "\n".join(lines)
    print(body)

    if args.out_tsv:
        header = ["Family_ID"]
        for d in domain_slots:
            header += [f"Domain{d}_Real_DeltaLogL", f"Domain{d}_Real_Z", f"Domain{d}_Pooled_Empirical_P"]

        with open(args.out_tsv, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(header)
            for fam in sorted(by_family.keys(), key=lambda s: int(s[1:]) if s[1:].isdigit() else 0):
                row = [fam]
                for d in domain_slots:
                    if fam in real_z[d]:
                        real_delta = by_family[fam][d][0]
                        row += [real_delta, real_z[d][fam], pooled_p[d][fam]]
                    else:
                        row += ["NA", "NA", "NA"]
                writer.writerow(row)
        print(f"\nPer-family, per-domain detail written to: {args.out_tsv}")

    if args.email:
        send_email(args.email, "Pooled-null z-score test: real vs. null delta-logL", body, args.out_tsv)
        print(f"\nEmailed to {args.email}", file=sys.stderr)


if __name__ == "__main__":
    main()
