#!/usr/bin/env python3
"""
Scan every family targeted by the bulk <=200-taxa run and report which ones
finished (have step6_family/family_summary.tsv) vs which didn't. Optionally
compile the unfinished families into a fresh set of mop-up chunk files, using
the same one-family-per-line chunk format run_family_chunk.sh expects.

A family counts as "done" purely by output presence -- run_family.sh / every
step script under it already skips completed work, so re-listing an
unfinished family in a mop-up chunk and rerunning it just resumes from
wherever it stopped, it doesn't redo finished work.

Usage:
  check_completion.py --chunks-dir DIR --output-base DIR [--mopup-dir DIR]
                       [--mopup-chunk-size N]

  --chunks-dir        Directory of chunk_NNNN.txt files (the ones already
                       dispatched) -- used to build the full target family
                       list. All of them are scanned, not just one.
  --output-base        Same OUTPUT_BASE passed to run_family.sh /
                       run_family_chunk.sh, e.g. data/sim_outputs_bulk.
  --mopup-dir          If given, write mop-up chunk_NNNN.txt files here for
                       whichever families are not yet done.
  --mopup-chunk-size   Families per mop-up chunk (default 12 -- smaller than
                       the original 36, since the unfinished set is enriched
                       for the slower/larger families that didn't make it
                       through a first-pass chunk in time).
"""
import argparse
import glob
import math
import os
import sys


def load_target_families(chunks_dir):
    families = []
    for chunk_path in sorted(glob.glob(os.path.join(chunks_dir, "chunk_*.txt"))):
        with open(chunk_path) as f:
            for line in f:
                fam = line.strip()
                if fam:
                    families.append(fam)
    return families


def is_done(output_base, family):
    summary_path = os.path.join(output_base, family, "step6_family", "family_summary.tsv")
    if not os.path.isfile(summary_path):
        return False
    # must have a header + at least one data row, not just be touched/empty
    with open(summary_path) as f:
        lines = [l for l in f if l.strip()]
    return len(lines) >= 2


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chunks-dir", required=True)
    ap.add_argument("--output-base", required=True)
    ap.add_argument("--mopup-dir")
    ap.add_argument("--mopup-chunk-size", type=int, default=12)
    args = ap.parse_args()

    families = load_target_families(args.chunks_dir)
    if not families:
        print(f"ERROR: no families found under {args.chunks_dir}/chunk_*.txt", file=sys.stderr)
        sys.exit(1)

    done, missing = [], []
    for fam in families:
        (done if is_done(args.output_base, fam) else missing).append(fam)

    n_total = len(families)
    n_done = len(done)
    n_missing = len(missing)
    pct = 100.0 * n_done / n_total

    print(f"Target families:   {n_total}")
    print(f"Done:               {n_done} ({pct:.1f}%)")
    print(f"Not done:           {n_missing}")

    if not missing:
        print("\nAll target families complete. No mop-up needed.")
        return

    if args.mopup_dir:
        os.makedirs(args.mopup_dir, exist_ok=True)
        # clear any stale mop-up chunks from a previous pass
        for stale in glob.glob(os.path.join(args.mopup_dir, "chunk_*.txt")):
            os.remove(stale)

        chunk_size = args.mopup_chunk_size
        n_chunks = math.ceil(n_missing / chunk_size)
        for i in range(n_chunks):
            chunk_families = missing[i * chunk_size : (i + 1) * chunk_size]
            chunk_path = os.path.join(args.mopup_dir, f"chunk_{i + 1:04d}.txt")
            with open(chunk_path, "w") as f:
                f.write("\n".join(chunk_families) + "\n")

        print(f"\nWrote {n_chunks} mop-up chunk(s) of up to {chunk_size} families "
              f"to {args.mopup_dir}")
        print(f"Resubmit with: sbatch --array=1-{n_chunks} "
              f"hpc_templates/mopup_le200_array.slurm")
    else:
        print("\nNot-done families:")
        for fam in missing:
            print(f"  {fam}")


if __name__ == "__main__":
    main()
