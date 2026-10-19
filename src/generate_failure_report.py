#!/usr/bin/env python3
"""
Build a running summary of AU-test failure rates across every family
processed so far in a bulk run, and optionally email it.

For each family, scans however many domains currently have a completed
Step 6 domain summary (data/sim_outputs_bulk/<family>/Domain<N>/step6/results.tsv)
-- this is live/partial-safe: a family with only Domain1 done shows Domain1's
numbers and NA for Domain2, rather than waiting for the whole family to
finish. One row per family, one block of columns per domain (a column-name
union across families, since some have more domains than others -- see
src/combine_family_results.py for the same pattern).

"Failed the AU test" is reported two ways per domain, since it's ambiguous
which one a reader wants:
  - AU_Reject_* : raw AU rejections (p < 0.05), regardless of cause
  - Significant_* : AU rejections that also reflect genuine loss of
    monophyly (rejected AND mono=false) -- the definition used everywhere
    else in this pipeline (see docs/METHODS.md), matching how Peter's real
    analysis presumably calls a rejection

Also includes, per domain: every replicate's own delta-logL (not just the
mean/min/max), for downstream comparison against the real delta-logL once
available; and a real-vs-simulated %ID sanity check (skipped by default --
pass --with-pctid to include it; it reads every alignment file and is by
far the slowest part of this script once the dataset reaches thousands of
families).

Usage:
  python3 src/generate_failure_report.py --output-base data/sim_outputs_bulk \
      --manifest data/manifest/manifest_full.json \
      --out-tsv data/reports/failure_report.tsv \
      [--email you@example.com] [--with-pctid]
"""
import argparse
import csv
import glob
import json
import re
import smtplib
import sys
from collections import Counter
from datetime import datetime
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from pathlib import Path

from pooled_null_z_test import (
    is_degenerate_null, is_search_unreliable, leave_one_out_z, real_z_score,
    empirical_p_upper, bh_qvalues,
)
from fisher_combined_test import fishers_combined

DOMAIN_RESULTS_RE = re.compile(r'^(F\d+)/Domain(\d+)/step6/results\.tsv$')


def read_fasta_seqs(path):
    """Ordered list of sequences (ignoring headers) from a FASTA alignment."""
    seqs = []
    cur = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if cur:
                    seqs.append(''.join(cur))
                cur = []
            else:
                cur.append(line)
        if cur:
            seqs.append(''.join(cur))
    return seqs


def pct_identity(seqs):
    """Pooled pairwise identity across an alignment: for every column, count
    matching vs comparable non-gap residue pairs, sum across columns, divide.
    Equivalent to averaging over all taxon pairs but O(taxa x sites) instead
    of O(taxa^2 x sites) -- matters once families have 100+ taxa. Returns
    None if there's nothing comparable (e.g. <2 sequences)."""
    if len(seqs) < 2:
        return None
    total_match = 0
    total_compare = 0
    for col in zip(*seqs):
        counts = Counter(c for c in col if c != '-')
        n = sum(counts.values())
        total_compare += n * (n - 1) // 2
        for cnt in counts.values():
            total_match += cnt * (cnt - 1) // 2
    if total_compare == 0:
        return None
    return 100.0 * total_match / total_compare


def compute_real_pct_id(trim_path):
    if not trim_path or not Path(trim_path).is_file():
        return None
    return pct_identity(read_fasta_seqs(trim_path))


def compute_mean_sim_pct_id(step3_dir):
    """Mean %ID across every rep_*_masked.aln found for this domain (however
    many replicates have made it through Step 3 so far)."""
    rep_files = sorted(glob.glob(str(Path(step3_dir) / "rep_*_masked.aln")))
    if not rep_files:
        return None
    vals = []
    for rf in rep_files:
        pid = pct_identity(read_fasta_seqs(rf))
        if pid is not None:
            vals.append(pid)
    if not vals:
        return None
    return sum(vals) / len(vals)


def read_replicate_deltalogl(step5_dir):
    """Per-replicate delta-logL values (unconstrained - constrained) from
    step5_results.txt, ordered by replicate number (rep_001, rep_002, ...).
    Returns a list, possibly shorter than N_REPLICATES if some replicates
    were skipped -- missing ones are left as None, not silently dropped,
    so position in the list always corresponds to that replicate number."""
    path = Path(step5_dir) / "step5_results.txt"
    if not path.is_file():
        return []
    values = {}
    max_rep = 0
    for line in open(path):
        m = re.match(r'^RESULT rep_(\d+) .*?\bdeltalogl=([^\s]+)', line)
        if not m:
            continue
        rep_num = int(m.group(1))
        max_rep = max(max_rep, rep_num)
        try:
            values[rep_num] = float(m.group(2))
        except ValueError:
            values[rep_num] = None
    return [values.get(i) for i in range(1, max_rep + 1)]


def read_iqtree_treelen(iqtree_path):
    """'Total tree length (sum of branch lengths)' from an .iqtree report.
    Returns None if the file or line doesn't exist."""
    if not Path(iqtree_path).is_file():
        return None
    for line in open(iqtree_path):
        if "Total tree length (sum of branch lengths)" in line:
            m = re.search(r'[-0-9.]+\s*$', line)
            if m:
                try:
                    return float(m.group(0))
                except ValueError:
                    return None
    return None


def read_real_treelens(output_base, fam, domain):
    """Real constrained/unconstrained tree lengths for one domain.
    Constrained: Step 1's brlen_estimate.iqtree (fixed-topology fit on the
    real constrained tree) -- available wherever Step 1 has run, independent
    of this domain having a Peter tree. Unconstrained: the fixed-topology
    evaluation of Peter's real tree (scripts/compute_real_unconstrained_logl.sh)
    where available; for the domains with no matching Peter tree at all,
    falls back to our OWN free-search tree (scripts/build_real_unconstrained_tree.sh,
    real_data/built_unconstrained.iqtree) -- a different, self-built tree, not
    necessarily the one Peter's own AU test used (see review_response.md
    question A1). Returns (constrained_treelen, unconstrained_treelen,
    used_self_built), the first two possibly None."""
    constrained = read_iqtree_treelen(
        output_base / fam / f"Domain{domain}" / "step1" / "brlen_estimate.iqtree")
    real_data_dir = output_base / fam / f"Domain{domain}" / "real_data"
    unconstrained = read_iqtree_treelen(real_data_dir / "real_unconstrained_logl.iqtree")
    used_self_built = False
    if unconstrained is None:
        unconstrained = read_iqtree_treelen(real_data_dir / "built_unconstrained.iqtree")
        used_self_built = unconstrained is not None
    return constrained, unconstrained, used_self_built


def read_iqtree_logl(iqtree_path):
    """'Log-likelihood of the tree: <value> (s.e. ...)' from an .iqtree
    report. Returns None if the file or line doesn't exist."""
    if not Path(iqtree_path).is_file():
        return None
    for line in open(iqtree_path):
        if "Log-likelihood of the tree" in line:
            m = re.search(r'tree:\s*([-0-9.]+)', line)
            if m:
                try:
                    return float(m.group(1))
                except ValueError:
                    return None
    return None


def read_real_logls(output_base, fam, domain):
    """Real constrained/unconstrained log-likelihoods for one domain, and
    their delta (unconstrained - constrained) -- the real-data equivalent
    of the simulated replicates' delta-logL, for direct comparison against
    Domain{d}_Mean_Null_DeltaLogL. Constrained: Step 1's brlen_estimate.iqtree.
    Unconstrained: read straight out of real_unconstrained_logl.txt's
    UNCONSTRAINED_LOGL= line (written by scripts/compute_real_unconstrained_logl.sh)
    where Peter's own tree was available; for domains with no matching Peter
    tree at all, falls back to our OWN free-search tree's
    BUILT_UNCONSTRAINED_LOGL= line (real_data/built_unconstrained.txt, from
    scripts/build_real_unconstrained_tree.sh) -- a different, self-built
    tree, not necessarily the one Peter's own AU test used (see
    review_response.md question A1). Returns (constrained_logl,
    unconstrained_logl, delta_logl, used_self_built), the first three
    possibly None."""
    constrained = read_iqtree_logl(
        output_base / fam / f"Domain{domain}" / "step1" / "brlen_estimate.iqtree")

    real_data_dir = output_base / fam / f"Domain{domain}" / "real_data"
    unconstrained = None
    txt_path = real_data_dir / "real_unconstrained_logl.txt"
    if txt_path.is_file():
        for line in open(txt_path):
            if line.startswith("UNCONSTRAINED_LOGL="):
                try:
                    unconstrained = float(line.strip().split("=", 1)[1])
                except ValueError:
                    pass

    used_self_built = False
    if unconstrained is None:
        built_txt_path = real_data_dir / "built_unconstrained.txt"
        if built_txt_path.is_file():
            for line in open(built_txt_path):
                if line.startswith("BUILT_UNCONSTRAINED_LOGL="):
                    try:
                        unconstrained = float(line.strip().split("=", 1)[1])
                        used_self_built = True
                    except ValueError:
                        pass

    delta = (unconstrained - constrained) if (constrained is not None and unconstrained is not None) else None
    return constrained, unconstrained, delta, used_self_built


def family_sort_key(fam):
    m = re.match(r'^F(\d+)$', fam)
    return int(m.group(1)) if m else fam


def load_domain_results(output_base):
    """Returns {family: {domain_num: {column: value}}}."""
    families = {}
    for path in output_base.glob("*/Domain*/step6/results.tsv"):
        rel = path.relative_to(output_base)
        m = DOMAIN_RESULTS_RE.match(str(rel))
        if not m:
            continue
        fam, dom = m.group(1), m.group(2)
        with open(path, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            row = next(reader, None)
        if row is None:
            continue
        families.setdefault(fam, {})[dom] = row
    return families


def load_manifest_info(manifest_path):
    """Returns ({family: n_domains}, {(family, domain): trim_path})."""
    if not manifest_path or not Path(manifest_path).is_file():
        return {}, {}
    manifest = json.load(open(manifest_path))
    counts = {fam: len(info["domains"]) for fam, info in manifest.items()}
    trim_paths = {
        (fam, dom): dom_info["trim"]
        for fam, info in manifest.items()
        for dom, dom_info in info["domains"].items()
    }
    return counts, trim_paths


def build_table(families, expected_counts, trim_paths, output_base, skip_pctid=True):
    all_domain_nums = sorted({d for doms in families.values() for d in doms.keys()}, key=int)

    # Max replicate count, from the N_Replicates column already loaded in
    # memory (load_domain_results) -- no extra file reads needed. This is
    # a safe upper bound on the actual replicate NUMBER ever seen in a
    # step5_results.txt (N_REPLICATES controls how many get simulated in
    # the first place), so columns never come up short even though
    # N_Replicates itself counts usable (non-skipped) replicates, not the
    # max rep number.
    max_reps = 10
    for doms in families.values():
        for r in doms.values():
            try:
                max_reps = max(max_reps, int(r["N_Replicates"]))
            except (KeyError, ValueError):
                pass

    header = ["Family_ID", "N_Domains_Expected", "N_Domains_Completed"]
    for d in all_domain_nums:
        header += [
            f"Domain{d}_N_Taxa", f"Domain{d}_Alignment_Length", f"Domain{d}_N_Replicates",
            # Real tree branch lengths (rate-of-evolution proxy), near the
            # left alongside the other basic domain descriptors. Constrained:
            # Step 1's fixed-topology fit on the real constrained tree
            # (available wherever Step 1 has run). Unconstrained: the
            # fixed-topology evaluation of Peter's real tree (available only
            # for domains with a matching real tree -- NA otherwise).
            f"Domain{d}_Real_Constrained_TreeLen", f"Domain{d}_Real_Unconstrained_TreeLen",
            # Raw real log-likelihoods, same sources as the tree lengths
            # above -- so Real_DeltaLogL below is verifiable from numbers
            # visible in this same table, not just a black-box subtraction.
            f"Domain{d}_Real_Constrained_LogL", f"Domain{d}_Real_Unconstrained_LogL",
            f"Domain{d}_N_AU_Reject_Total", f"Domain{d}_Pct_AU_Reject_Total",
            f"Domain{d}_N_Significant_Violations", f"Domain{d}_Pct_Significant_Violations",
            f"Domain{d}_Mean_AU_Pvalue",
            f"Domain{d}_Mean_Null_DeltaLogL", f"Domain{d}_Real_DeltaLogL",
            f"Domain{d}_Real_PctID", f"Domain{d}_Mean_Sim_PctID",
        ]
        # Every replicate's own ΔlogL (unconstrained - constrained), not just
        # the mean -- lets downstream analysis compare the real ΔlogL against
        # this null DISTRIBUTION directly (e.g. a permutation p-value or any
        # other summary stat), rather than only the mean/min/max already in
        # Domain{d}_Mean_Null_DeltaLogL. NA for any replicate that was
        # skipped or hasn't run yet.
        header += [f"Domain{d}_DeltaLogL_Rep{i}" for i in range(1, max_reps + 1)]
        # Human-readable, re-derived-from-disk-every-run verification summary
        # for this domain -- what's actually been confirmed complete, and
        # what's known to be missing (not yet run, no real tree, replicates
        # that never produced a usable result, excluded from the pooled-null
        # test). See add_domain_notes() / finalize_notes().
        header.append(f"Domain{d}_Notes")

    rows = []
    for fam in sorted(families.keys(), key=family_sort_key):
        doms = families[fam]
        n_expected = expected_counts.get(fam, "NA")
        row = {
            "Family_ID": fam,
            "N_Domains_Expected": n_expected,
            "N_Domains_Completed": len(doms),
        }
        for d in all_domain_nums:
            r = doms.get(d)
            if r is None:
                for col in ("N_Taxa", "Alignment_Length", "N_Replicates",
                            "Real_Constrained_TreeLen", "Real_Unconstrained_TreeLen",
                            "Real_Constrained_LogL", "Real_Unconstrained_LogL",
                            "N_AU_Reject_Total", "Pct_AU_Reject_Total",
                            "N_Significant_Violations", "Pct_Significant_Violations",
                            "Mean_AU_Pvalue", "Mean_Null_DeltaLogL", "Real_DeltaLogL",
                            "Real_PctID", "Mean_Sim_PctID"):
                    row[f"Domain{d}_{col}"] = "NA"
                for i in range(1, max_reps + 1):
                    row[f"Domain{d}_DeltaLogL_Rep{i}"] = "NA"

                # Cheap (isdir only, no file stats/walks) check of how far
                # this domain got, re-derived from disk every run -- not a
                # live SLURM-state guess. A domain genuinely queued behind a
                # large running job can sit at the same furthest step for
                # days without being "stuck"; this reports the objective
                # fact only, not an active/stalled judgement.
                furthest = "not started"
                for step in ("step1", "step2", "step3", "step4", "step5", "step6"):
                    if (output_base / fam / f"Domain{d}" / step).is_dir():
                        furthest = step
                row[f"Domain{d}_Notes"] = f"Not yet complete (furthest step reached: {furthest})"
            else:
                row[f"Domain{d}_N_Taxa"] = r["N_Taxa"]
                row[f"Domain{d}_Alignment_Length"] = r["Alignment_Length"]
                row[f"Domain{d}_N_Replicates"] = r["N_Replicates"]

                constrained_tl, unconstrained_tl, tl_self_built = read_real_treelens(output_base, fam, d)
                row[f"Domain{d}_Real_Constrained_TreeLen"] = constrained_tl if constrained_tl is not None else "NA"
                row[f"Domain{d}_Real_Unconstrained_TreeLen"] = unconstrained_tl if unconstrained_tl is not None else "NA"

                real_c_logl, real_u_logl, real_delta, logl_self_built = read_real_logls(output_base, fam, d)
                row[f"Domain{d}_Real_Constrained_LogL"] = real_c_logl if real_c_logl is not None else "NA"
                row[f"Domain{d}_Real_Unconstrained_LogL"] = real_u_logl if real_u_logl is not None else "NA"
                used_self_built_tree = tl_self_built or logl_self_built

                row[f"Domain{d}_N_AU_Reject_Total"] = r["N_AU_Rejects_Total"]
                row[f"Domain{d}_Pct_AU_Reject_Total"] = r["Pct_AU_Rejects_Total"]
                row[f"Domain{d}_N_Significant_Violations"] = r["N_Significant_Violations"]
                row[f"Domain{d}_Pct_Significant_Violations"] = r["Pct_Significant_Violations"]
                row[f"Domain{d}_Mean_AU_Pvalue"] = r["Mean_AU_Pvalue"]
                # Mean of ALL replicates' delta-logL (unconstrained - constrained),
                # recorded unconditionally in Step 5 regardless of AU reject/accept
                # -- see docs/METHODS.md. Real_DeltaLogL is the REAL equivalent,
                # computed above from Real_Constrained_LogL/Real_Unconstrained_LogL
                # (the cheap real-tree calculation), not the old expensive
                # real_data_delta_logl.sh pathway (never run in bulk, always NA) --
                # this is the number to compare directly against
                # Mean_Null_DeltaLogL, and the assumption this whole pipeline
                # rests on: that matching alignment size/structure in simulation
                # reproduces a delta-logL comparable to the real one.
                row[f"Domain{d}_Mean_Null_DeltaLogL"] = r["Null_DeltaLogL_Mean"]
                row[f"Domain{d}_Real_DeltaLogL"] = real_delta if real_delta is not None else "NA"

                # Pairwise %ID: real alignment vs. mean across whatever
                # simulated replicates have reached Step 3 so far -- a
                # sanity check that the simulation is actually mimicking the
                # real data's divergence, not just running without error.
                # Reads full alignment files (real + up to N_REPLICATES
                # simulated) per domain -- by far the slowest part of this
                # script once the dataset reaches thousands of families;
                # skipped by default, pass --with-pctid to include it.
                if skip_pctid:
                    row[f"Domain{d}_Real_PctID"] = "NA"
                    row[f"Domain{d}_Mean_Sim_PctID"] = "NA"
                else:
                    trim_path = trim_paths.get((fam, d))
                    real_pid = compute_real_pct_id(trim_path)
                    step3_dir = output_base / fam / f"Domain{d}" / "step3"
                    sim_pid = compute_mean_sim_pct_id(step3_dir)
                    row[f"Domain{d}_Real_PctID"] = round(real_pid, 2) if real_pid is not None else "NA"
                    row[f"Domain{d}_Mean_Sim_PctID"] = round(sim_pid, 2) if sim_pid is not None else "NA"

                step5_dir = output_base / fam / f"Domain{d}" / "step5"
                rep_vals = read_replicate_deltalogl(step5_dir)
                for i in range(1, max_reps + 1):
                    v = rep_vals[i - 1] if i <= len(rep_vals) else None
                    row[f"Domain{d}_DeltaLogL_Rep{i}"] = v if v is not None else "NA"

                # Verification notes for a COMPLETED domain. Stored as a
                # list under a temp key (not in `header`, deleted before
                # writing) -- add_pooled_z_stats() appends one more flag
                # (degenerate-null exclusion) once it's computed later, then
                # finalize_notes() joins everything into the final string.
                notes = []
                if row[f"Domain{d}_Real_Unconstrained_LogL"] == "NA":
                    notes.append("No real (Peter) tree available, and no self-built tree either")
                elif used_self_built_tree:
                    notes.append("Using self-built unconstrained tree (Peter's own tree unavailable) "
                                 "-- not necessarily the same tree Peter's AU test used")

                try:
                    n_replicates = int(r["N_Replicates"])
                    n_skipped = int(r.get("N_Skipped_Replicates", 0))
                except (KeyError, ValueError):
                    n_replicates, n_skipped = None, 0
                total_attempted = (n_replicates + n_skipped) if n_replicates is not None else max_reps

                n_na = sum(
                    1 for i in range(1, total_attempted + 1)
                    if row.get(f"Domain{d}_DeltaLogL_Rep{i}", "NA") == "NA"
                )
                if n_skipped > 0:
                    notes.append(f"{n_skipped} replicate(s) skipped entirely "
                                 f"(missing alignment/trees at Step4/5 input)")
                extra_anomalies = n_na - n_skipped
                if extra_anomalies > 0:
                    notes.append(f"{extra_anomalies} replicate(s) produced a result but with "
                                 f"no usable ΔlogL (see Step5 log)")

                # Null-side search reliability: a NEGATIVE null ΔlogL means
                # the free (unconstrained) search did WORSE than the
                # constrained tree it was supposed to beat -- a search that
                # can't even reliably outperform a reachable starting point
                # signals too little phylogenetic information in this
                # domain's taxa count/alignment length/tree length combo to
                # trust ANY comparison built on it, independent of which way
                # Pooled_Z/P/BH comes out. Threshold (>=50% of valid
                # replicates negative) chosen from the dataset's own
                # distribution: 0.9% of domains clear it, concentrated
                # almost entirely in large-taxa/short-alignment combinations
                # (e.g. 2000+ taxa, <150 sites) -- a real, systematic
                # pattern, not noise. See docs/METHODS.md discussion.
                null_vals = [row[f"Domain{d}_DeltaLogL_Rep{i}"] for i in range(1, total_attempted + 1)
                             if row.get(f"Domain{d}_DeltaLogL_Rep{i}", "NA") != "NA"]
                if null_vals and is_search_unreliable(null_vals):
                    n_null_negative = sum(1 for v in null_vals if v < 0)
                    notes.append(
                        f"{n_null_negative}/{len(null_vals)} null replicates have a negative "
                        f"ΔlogL (free search did worse than the constrained tree) -- "
                        f"simulation search reliability questionable for this domain's "
                        f"taxa count/alignment length/tree length; treat results with caution")

                row[f"Domain{d}_Notes__list"] = notes
        rows.append(row)

    return header, rows, all_domain_nums


def pooled_stats(families, all_domain_nums):
    """Pooled (population-level) significant-violation rate per domain slot,
    across every replicate of every family -- the actual number that matters
    for comparing against the real 41% figure, not a mean-of-percentages."""
    stats = {}
    for d in all_domain_nums:
        n_sig, n_reps, n_fams = 0, 0, 0
        for doms in families.values():
            r = doms.get(d)
            if r is None:
                continue
            n_sig += int(r["N_Significant_Violations"])
            n_reps += int(r["N_Replicates"])
            n_fams += 1
        stats[d] = {
            "n_families": n_fams,
            "n_replicates": n_reps,
            "n_significant": n_sig,
            "pct_significant": round(100.0 * n_sig / n_reps, 2) if n_reps else None,
        }
    return stats


def add_pooled_z_stats(header, rows, all_domain_nums):
    """Adds, per domain, Pooled_Z / Pooled_P / Pooled_P_BH columns -- the
    higher-resolution alternative to the raw per-family rank test, from
    src/pooled_null_z_test.py. Standardizes each family's real delta-logL
    and null replicates using that family's own 10 replicates, pools the
    *null-side* standardized values (leave-one-out) across every family
    in this report into one empirical reference per domain, then compares
    each family's real z-score against that pooled reference. Families
    with numerically degenerate null variance (near-zero sd) are
    excluded, detected from their own null-vs-null z's only -- never from
    the real value, which can legitimately be enormous without the null
    being degenerate. See docs/METHODS.md / pvalue_methods.md for the
    full derivation. Mutates rows in place; returns the updated header
    (new columns inserted right after each domain's Real_DeltaLogL) and a
    per-domain stats dict for the summary email."""
    domain_stats = {}

    for d in all_domain_nums:
        rep_pattern = re.compile(rf'^Domain{d}_DeltaLogL_Rep\d+$')
        candidates = []
        for row in rows:
            real_raw = row.get(f"Domain{d}_Real_DeltaLogL", "NA")
            if real_raw == "NA":
                continue
            null_deltas = [float(v) for k, v in row.items()
                           if rep_pattern.match(k) and v != "NA"]
            if len(null_deltas) < 3:
                continue
            candidates.append((row, float(real_raw), null_deltas))

        pooled_null_z = []
        real_z_rows = []  # (row, z, is_unreliable)
        n_degenerate = 0
        for row, real_delta, null_deltas in candidates:
            if is_degenerate_null(null_deltas):
                n_degenerate += 1
                row.setdefault(f"Domain{d}_Notes__list", []).append(
                    "Excluded from pooled-null test: degenerate null variance "
                    "(own simulated replicates are numerically ~identical)")
                continue
            pooled_null_z.extend(leave_one_out_z(null_deltas))
            rz = real_z_score(real_delta, null_deltas)
            if rz is not None:
                real_z_rows.append((row, rz, is_search_unreliable(null_deltas)))

        p_vals = [empirical_p_upper(z, pooled_null_z) for _, z, _ in real_z_rows]
        q_vals = bh_qvalues(p_vals)

        for (row, z, _), p, q in zip(real_z_rows, p_vals, q_vals):
            row[f"Domain{d}_Pooled_Z"] = round(z, 4)
            row[f"Domain{d}_Pooled_P"] = p
            row[f"Domain{d}_Pooled_P_BH"] = round(q, 6)

        for row in rows:
            row.setdefault(f"Domain{d}_Pooled_Z", "NA")
            row.setdefault(f"Domain{d}_Pooled_P", "NA")
            row.setdefault(f"Domain{d}_Pooled_P_BH", "NA")

        # Reliability split: the Pooled_Z/P/BH columns above are unchanged
        # (same single pool, same global BH ranking, as before) -- this is
        # an additional breakdown layered on top, grouping the ALREADY-
        # computed q-values by whether each family's own null search was
        # reliable (see is_search_unreliable -- >=50% of its null
        # replicates failed to beat the constrained tree). Reliable
        # families are the ones a real-vs-null comparison can actually be
        # trusted for; unreliable ones are reported separately rather than
        # silently mixed into the same headline number.
        reliable_q = [q for (_, _, unreliable), q in zip(real_z_rows, q_vals) if not unreliable]
        unreliable_q = [q for (_, _, unreliable), q in zip(real_z_rows, q_vals) if unreliable]
        reliable_p = [p for (_, _, unreliable), p in zip(real_z_rows, p_vals) if not unreliable]
        unreliable_p = [p for (_, _, unreliable), p in zip(real_z_rows, p_vals) if unreliable]

        fisher_reliable = fishers_combined(reliable_p)
        fisher_unreliable = fishers_combined(unreliable_p)

        n_bh_sig = sum(1 for q in q_vals if q <= 0.05)
        domain_stats[d] = {
            "n_tested": len(real_z_rows),
            "n_degenerate_excluded": n_degenerate,
            "n_pooled_null": len(pooled_null_z),
            "n_bh_sig_05": n_bh_sig,
            "n_reliable": len(reliable_q),
            "n_bh_sig_05_reliable": sum(1 for q in reliable_q if q <= 0.05),
            "fisher_reliable": fisher_reliable,
            "n_unreliable": len(unreliable_q),
            "n_bh_sig_05_unreliable": sum(1 for q in unreliable_q if q <= 0.05),
            "fisher_unreliable": fisher_unreliable,
        }

    new_header = []
    for col in header:
        new_header.append(col)
        m = re.match(r'^Domain(\d+)_Real_DeltaLogL$', col)
        if m:
            d = m.group(1)
            new_header += [f"Domain{d}_Pooled_Z", f"Domain{d}_Pooled_P", f"Domain{d}_Pooled_P_BH"]

    return new_header, domain_stats


def finalize_notes(rows, all_domain_nums):
    """Joins each domain's Domain{d}_Notes__list (built across build_table
    and add_pooled_z_stats) into the final Domain{d}_Notes string, and
    removes the temp list key so it never reaches the CSV writer. A
    completed domain with an empty list -- nothing flagged by any of the
    checks above -- is explicitly confirmed, not just silently blank."""
    for row in rows:
        for d in all_domain_nums:
            key = f"Domain{d}_Notes__list"
            if key in row:
                notes = row.pop(key)
                row[f"Domain{d}_Notes"] = "; ".join(notes) if notes else "Verified complete"


def send_email(to_addr, subject, body, attachment_path):
    msg = MIMEMultipart()
    msg["Subject"] = subject
    msg["From"] = f"au-test-sim@{Path('/etc/hostname').read_text().strip() if Path('/etc/hostname').is_file() else 'localhost'}"
    msg["To"] = to_addr
    msg.attach(MIMEText(body, "plain"))

    with open(attachment_path, "rb") as f:
        part = MIMEApplication(f.read(), Name=Path(attachment_path).name)
    part["Content-Disposition"] = f'attachment; filename="{Path(attachment_path).name}"'
    msg.attach(part)

    with smtplib.SMTP("localhost") as smtp:
        smtp.send_message(msg)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--output-base", required=True)
    ap.add_argument("--manifest", default=None)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--email", default=None, help="If given, email the report to this address")
    ap.add_argument("--with-pctid", action="store_true",
                     help="Include the real-vs-simulated %%ID check (reads every alignment file; "
                          "by far the slowest part once the dataset reaches thousands of families -- "
                          "skipped by default)")
    args = ap.parse_args()

    output_base = Path(args.output_base)
    families = load_domain_results(output_base)
    expected_counts, trim_paths = load_manifest_info(args.manifest)

    if not families:
        print(f"No domain results found yet under {output_base}", file=sys.stderr)
        sys.exit(1)

    header, rows, all_domain_nums = build_table(families, expected_counts, trim_paths, output_base,
                                                 skip_pctid=not args.with_pctid)
    stats = pooled_stats(families, all_domain_nums)
    header, pooled_z_stats = add_pooled_z_stats(header, rows, all_domain_nums)
    finalize_notes(rows, all_domain_nums)

    out_path = Path(args.out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter="\t", restval="NA")
        writer.writeheader()
        writer.writerows(rows)

    n_complete_families = sum(
        1 for fam, doms in families.items()
        if expected_counts.get(fam) is not None and len(doms) >= expected_counts.get(fam, 0)
    )

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    summary_lines = [
        f"AU test simulation -- running failure-rate report, generated {timestamp}",
        "",
        f"Families with at least one domain completed: {len(families)}",
        f"Families with ALL expected domains completed: {n_complete_families}",
        "",
        "Pooled (population-level) significant-violation rate per domain slot",
        "(AU rejected AND composite clade genuinely lost monophyly -- the",
        "headline null false-positive rate, pooled across every replicate of",
        "every family with that domain completed so far, not a mean of",
        "per-family percentages):",
        "",
    ]
    for d in all_domain_nums:
        s = stats[d]
        pct = f"{s['pct_significant']}%" if s['pct_significant'] is not None else "NA"
        summary_lines.append(
            f"  Domain{d}: {s['n_significant']} / {s['n_replicates']} replicates "
            f"({pct}) across {s['n_families']} families"
        )

    summary_lines += [
        "",
        "Pooled-null z-score test (real delta-logL standardized per family,",
        "compared against a pooled, cross-family null reference -- see",
        "pvalue_methods.md for the full method). BH = Benjamini-Hochberg",
        "FDR-adjusted p-value (q-value), the appropriate multiple-testing",
        "correction here (Bonferroni is far too strict given how many",
        "families are tested at once):",
        "",
    ]
    for d in all_domain_nums:
        s = pooled_z_stats[d]
        pct_bh = f"{100*s['n_bh_sig_05']/s['n_tested']:.1f}%" if s['n_tested'] else "NA"
        summary_lines.append(
            f"  Domain{d}: {s['n_tested']} families tested (pooled null size={s['n_pooled_null']}, "
            f"{s['n_degenerate_excluded']} excluded as numerically degenerate); "
            f"significant at BH-FDR<=0.05: {s['n_bh_sig_05']}/{s['n_tested']} ({pct_bh})"
        )

    summary_lines += [
        "",
        "Same test, split by simulation reliability (>=50% of a family's own",
        "null replicates failing to beat the constrained tree means that",
        "family's null -- and so any comparison built on it -- shouldn't be",
        "trusted, regardless of which way it comes out; see each domain's",
        "Notes column). RELIABLE is the number to use for \"is 41% surprising\";",
        "UNRELIABLE is reported separately, not folded into that headline:",
        "",
    ]
    for d in all_domain_nums:
        s = pooled_z_stats[d]
        pct_rel = f"{100*s['n_bh_sig_05_reliable']/s['n_reliable']:.1f}%" if s['n_reliable'] else "NA"
        pct_unrel = f"{100*s['n_bh_sig_05_unreliable']/s['n_unreliable']:.1f}%" if s['n_unreliable'] else "NA"
        _, log10_p_rel, k_rel = s['fisher_reliable']
        _, log10_p_unrel, k_unrel = s['fisher_unreliable']
        fisher_rel_str = f"log10(p)={log10_p_rel:.1f}" if log10_p_rel is not None else "NA"
        fisher_unrel_str = f"log10(p)={log10_p_unrel:.1f}" if log10_p_unrel is not None else "NA"
        summary_lines.append(
            f"  Domain{d} RELIABLE:   {s['n_reliable']} families; "
            f"BH-FDR<=0.05: {s['n_bh_sig_05_reliable']}/{s['n_reliable']} ({pct_rel}); "
            f"Fisher's combined {fisher_rel_str} (k={k_rel})"
        )
        summary_lines.append(
            f"  Domain{d} UNRELIABLE: {s['n_unreliable']} families; "
            f"BH-FDR<=0.05: {s['n_bh_sig_05_unreliable']}/{s['n_unreliable']} ({pct_unrel}); "
            f"Fisher's combined {fisher_unrel_str} (k={k_unrel})"
        )

    summary_lines += [
        "",
        "Mean pairwise %ID, real data vs. simulated replicates (sanity check",
        "that the simulation is mimicking the real data's divergence -- big",
        "gaps here would mean the null distribution isn't trustworthy,",
        "regardless of what the AU-reject rate says):",
        "",
    ]
    for d in all_domain_nums:
        real_vals = [row[f"Domain{d}_Real_PctID"] for row in rows if row[f"Domain{d}_Real_PctID"] != "NA"]
        sim_vals = [row[f"Domain{d}_Mean_Sim_PctID"] for row in rows if row[f"Domain{d}_Mean_Sim_PctID"] != "NA"]
        if real_vals and sim_vals:
            mean_real = sum(real_vals) / len(real_vals)
            mean_sim = sum(sim_vals) / len(sim_vals)
            summary_lines.append(
                f"  Domain{d}: real={mean_real:.2f}%  simulated={mean_sim:.2f}%  "
                f"(diff={mean_real - mean_sim:+.2f} points, n={len(real_vals)} families)"
            )
        else:
            summary_lines.append(f"  Domain{d}: NA (no data yet)")

    summary_lines += [
        "",
        f"Full per-family, per-domain table attached: {out_path.name}",
        f"(also written to: {out_path})",
    ]
    body = "\n".join(summary_lines)

    print(body)

    if args.email:
        subject = f"AU test sim report: {len(families)} families ({n_complete_families} fully complete), {timestamp}"
        send_email(args.email, subject, body, out_path)
        print(f"\nEmailed to {args.email}", file=sys.stderr)


if __name__ == "__main__":
    main()
