#!/usr/bin/env python3
"""
Extract wall-clock timing and basic resource info from a completed (or
in-progress) pipeline run for one family, and append it to a persistent
benchmarking log for HPC capacity planning.

Usage:
  python3 src/extract_timing.py --family-dir data/sim_outputs/F100058 \
      --manifest data/manifest/manifest.json \
      --log data/benchmarks/timing_log.tsv

Run this after (or during -- it only records what has actually finished)
each family's pipeline run. The resulting TSV is the source of truth for
estimating total core-hours and choosing a SLURM array-job batching
strategy: rows are (family, domain, step, replicate)-level wall-clock times,
which can be summed/extrapolated by n_taxa and alignment length across the
full family set without re-running anything.
"""

import argparse
import json
import re
import sys
from pathlib import Path

WALLCLOCK_RE = re.compile(r'Total wall-clock time used:\s*([\d.]+)\s*sec')
THREADS_RE = re.compile(r'Kernel:.*?(\d+)\s*threads', re.IGNORECASE)
# Matches both SeqGen's own banner line and apply_gap_pattern.py's explicit
# "Time taken: N.NNNNNN seconds" line.
TIME_TAKEN_RE = re.compile(r'Time taken:\s*([\d.]+)\s*seconds')


def _last_match(pattern: re.Pattern, text: str):
    """Return the LAST match's first group, not the first. Step logs can
    accumulate output from multiple re-runs (-f, debugging, resumed jobs) via
    shell `tee -a` / Python FileHandler append mode; the most recent run's
    numbers are always the last occurrence in the file."""
    matches = pattern.findall(text)
    return matches[-1] if matches else None


def extract_iqtree_wallclock(iqtree_or_log_path: Path):
    """Look for IQTree's own reported wall-clock time in a .iqtree or .log file."""
    if not iqtree_or_log_path.exists():
        return None
    text = iqtree_or_log_path.read_text(errors='replace')
    val = _last_match(WALLCLOCK_RE, text)
    return float(val) if val is not None else None


def step_log_span_seconds(log_path: Path):
    """Fall back: first-to-last timestamp span in a step*.log file
    (format: [YYYY-MM-DD HH:MM:SS] ...)."""
    if not log_path.exists():
        return None
    ts_re = re.compile(r'^\[(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\]')
    first = last = None
    for line in log_path.read_text(errors='replace').splitlines():
        m = ts_re.match(line)
        if m:
            if first is None:
                first = m.group(1)
            last = m.group(1)
    if first is None or last is None:
        return None
    from datetime import datetime
    fmt = '%Y-%m-%d %H:%M:%S'
    return (datetime.strptime(last, fmt) - datetime.strptime(first, fmt)).total_seconds()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--family-dir', required=True, help='e.g. data/sim_outputs/F100058')
    ap.add_argument('--manifest', default=None, help='manifest.json, used to look up n_taxa/aln_length/model per domain')
    ap.add_argument('--log', required=True, help='Persistent TSV to append rows to (created if missing)')
    ap.add_argument('--threads-per-job', type=int, default=1, help='IQTREE_THREADS_PER_JOB used for this run')
    args = ap.parse_args()

    family_dir = Path(args.family_dir)
    family = family_dir.name

    manifest_info = {}
    if args.manifest and Path(args.manifest).exists():
        with open(args.manifest) as f:
            manifest = json.load(f)
        fam_entry = manifest.get(family, {})
        for dom, info in fam_entry.get('domains', {}).items():
            manifest_info[dom] = info

    rows = []

    for domain_dir in sorted(family_dir.glob('Domain*')):
        domain = domain_dir.name.replace('Domain', '')
        dom_info = manifest_info.get(domain, {})
        n_taxa = dom_info.get('n_taxa', 'NA')
        aln_length = dom_info.get('aln_length', 'NA')
        model = dom_info.get('model', 'NA')

        # Step 2: SeqGen (one shot for all N_REPLICATES, not per-replicate)
        step2_log = domain_dir / 'step2' / 'step2.log'
        if step2_log.exists():
            val = _last_match(TIME_TAKEN_RE, step2_log.read_text(errors='replace'))
            seqgen_secs = float(val) if val is not None else step_log_span_seconds(step2_log)
            rows.append([family, domain, 'step2_seqgen', 'all', n_taxa, aln_length, model,
                         args.threads_per_job, seqgen_secs])

        # Step 3: masking (one shot for all N_REPLICATES)
        step3_log = domain_dir / 'step3' / 'step3.log'
        if step3_log.exists():
            val = _last_match(TIME_TAKEN_RE, step3_log.read_text(errors='replace'))
            mask_secs = float(val) if val is not None else step_log_span_seconds(step3_log)
            rows.append([family, domain, 'step3_masking', 'all', n_taxa, aln_length, model,
                         args.threads_per_job, mask_secs])

        # Step 4: per-replicate unconstrained + constrained tree search
        step4_dir = domain_dir / 'step4'
        for rep_dir in sorted(step4_dir.glob('rep_*')):
            rep = rep_dir.name
            for kind, fname in [('unconstrained', 'masked_unconstrained'), ('constrained', 'masked_constrained')]:
                iqtree_report = rep_dir / f'{fname}.iqtree'
                iqtree_log = rep_dir / f'{fname}.log'
                secs = extract_iqtree_wallclock(iqtree_report) or extract_iqtree_wallclock(iqtree_log)
                if secs is not None:
                    rows.append([family, domain, f'step4_{kind}', rep, n_taxa, aln_length, model,
                                 args.threads_per_job, secs])

        # Step 5: AU test per replicate (includes -zb RELL bootstrap cost)
        step5_dir = domain_dir / 'step5'
        for rep_dir in sorted(step5_dir.glob('rep_*')):
            rep = rep_dir.name
            iqtree_report = rep_dir / 'au_test_result.iqtree'
            iqtree_log = rep_dir / 'au_test_result.log'
            secs = extract_iqtree_wallclock(iqtree_report) or extract_iqtree_wallclock(iqtree_log)
            if secs is not None:
                rows.append([family, domain, 'step5_au_test', rep, n_taxa, aln_length, model,
                             args.threads_per_job, secs])

    if not rows:
        print(f"No completed timing data found under {family_dir}", file=sys.stderr)
        sys.exit(0)

    log_path = Path(args.log)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not log_path.exists()

    from datetime import datetime
    now = datetime.now().isoformat(timespec='seconds')

    with open(log_path, 'a') as f:
        if write_header:
            f.write('\t'.join([
                'family', 'domain', 'step', 'replicate', 'n_taxa', 'aln_length',
                'model', 'threads', 'wall_seconds', 'recorded_at',
            ]) + '\n')
        for row in rows:
            f.write('\t'.join(str(x) for x in row) + f'\t{now}\n')

    print(f"Appended {len(rows)} timing rows to {log_path}", file=sys.stderr)


if __name__ == '__main__':
    main()
