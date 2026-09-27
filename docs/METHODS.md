# Methods

## 1. Scientific purpose

Peter Mulhair's analysis of remodelled/composite gene histories (RHGs —
"reticulate homolog groups", genes formed by fusion of two parent genes)
found that **41% of ~10,829 RHGs appear polyphyletic** for the fused clade,
based on an AU (Approximately Unbiased, Shimodaira 2002) test rejecting the
hypothesis that the fused-gene clade is monophyletic.

The question this pipeline answers: **is 41% surprising?** Simulate
sequence data under the *true* null model — a single fusion event, so the
fused clade really is monophyletic — and measure how often the AU test
*incorrectly* rejects monophyly by chance. That gives the null/expected
false-positive rate to compare 41% against.

## 2. Why domain-split, not single-tree

A gene fusion is a reticulation event: the fused gene has two distinct
ancestries (one segment from Parent gene 1, one from Parent gene 2). A
single bifurcating tree can't represent that without arbitrarily favouring
one parent. The real analysis (and this pipeline) instead builds two
separate trees per fused family — one per parent-gene "domain" — each
using only the alignment region homologous to that parent:

- **Domain *N* side**: that parent gene's own orthologs + the composite
  (fused) taxa, using only the domain-homologous alignment columns.
- Build the unconstrained gene tree, build a monophyly-constrained version
  of the composite clade, run an AU test comparing them — independently per
  domain.

So each fused family yields one AU test *per domain* (2 for most families;
see "Multi-domain families" in [docs/INPUTS.md](INPUTS.md) for the minority
with more), not one. Real data already comes pre-split this way — no
reconstruction from a combined alignment is needed.

## 3. Pipeline architecture

Manifest-driven: `src/build_manifest.py` discovers `(family, domain)` units
from the real data files and computes a deterministic per-family seed, so
all of a family's domains draw from the same underlying simulation
scenario (replicate *i* means the same joint scenario across every domain
of a family). `scripts/run_family.sh` reads the manifest and loops Steps
1–6 over each domain, then the family-level join.

### Step 1 — branch lengths + per-site rates

Fixed-topology IQTree fit (`-te`) on the real alignment under the domain's
real model, plus per-site empirical-Bayes rate categories (`-wsr`), needed
by Step 2 to concentrate simulated rate variation at the same real
positions rather than spreading it uniformly.

### Step 2 — simulate replicates

Groups the real alignment's sites by their Step 1 rate category, then
simulates each category separately (via AliSim's `--branch-scale`) on the
Step 1 branch-length tree and stitches the category blocks back together
at their original column positions. This reproduces the real conservation
*pattern* (which columns are near-invariant), not just its average.

### Step 3 — gap masking

Positional overlay: wherever the real alignment has a gap at column *i*
for a taxon, the simulated sequence gets a gap at column *i* too.

### Step 4 — reconstruct trees

Per replicate: an unconstrained ML search and a constrained fit
(composite clade forced monophyletic). Replicates run concurrently. This
is the dominant compute cost in the pipeline (see §6).

### Step 5 — AU test

Per replicate: compares the two trees' likelihoods via IQTree's AU test,
records the constrained tree's p-value/ΔlogL, and independently checks
composite-clade monophyly in the unconstrained tree (rooting-invariant). A
rejection is classified as reflecting genuine loss of monophyly vs. other
topology changes. ΔlogL and both trees' logL are recorded for every
replicate, unconditionally, enabling a permutation-style comparison against
the real data's own ΔlogL without rerunning anything.

### Step 6 (domain)

Aggregates a domain's replicates: reject rate, mean p-value, mean/SD tree
length, and (optionally) a permutation p-value comparing the real data's
ΔlogL to the null distribution.

### Step 6 (family)

Joins a family's domains' replicate-level results by shared seed/index. A
replicate counts as a significant violation for a domain only if the AU
test rejected **and** that rejection reflects a genuine loss of monophyly.
Reports, per family:

- **OR-rate**: fraction of replicates where *at least one* domain shows a
  significant violation.
- **AND-rate**: fraction of replicates where *all* domains do.

Both rates (and each domain's own rate) are reported side by side rather
than collapsed into one family-level verdict — which one is the right
comparison against the real 41% figure is a downstream analysis decision.
The OR-rate is the more likely candidate, since testing 2+ trees per family
is exactly the kind of multiple-comparisons setup that can inflate a
nominal 5% single-test error rate at the family level.

`scripts/combine_results.sh` combines every family's summary into one
dataset-wide table (a column-name union, since families with 3+ domains
produce extra columns that 2-domain families don't have).

## 4. Resumability

Every step script skips work it finds already done unless a force-rerun
flag is passed. This makes it safe to batch many families into one SLURM
array task with a runtime cap that some batches won't fully finish within —
a killed/resubmitted run resumes rather than redoing completed work.

## 5. Statistical design notes

- Per-replicate p-values need no multiple-testing correction — each
  replicate estimates a single Bernoulli trial (did AU reject on this
  draw), and the summary statistic is a rate, not a significance test in
  itself.
- Per-domain-per-family OR-consolidation *is* a genuine multiple-comparisons
  concern (see §3, Step 6 family).
- `N_REPLICATES` is sized for population-level precision, not per-family
  precision. Pooled across thousands of families, the aggregate null
  FP-rate estimate (the number actually compared against 41%) has a small
  standard error even at a modest per-family replicate count; raising it
  mainly sharpens each individual family's own point estimate, at a
  compute cost that scales close to linearly with replicate count unless
  the parallelism budget (`REPLICATE_JOBS`) is raised to match (see §6).

## 6. HPC operating notes

Step 2 (simulation) is close to replicate-count-independent, since AliSim
batches all replicates for a rate category into one call. Steps 4 and 5 are
not — they run one IQTree invocation per replicate, parallelised via
`REPLICATE_JOBS` (an `xargs -P` pool sized to the job's allocated cores).
Wall-clock time for these steps is approximately
`ceil(N_REPLICATES / REPLICATE_JOBS) × (one replicate's search time)` — so
raising `N_REPLICATES` without raising the parallelism budget to match
increases runtime close to linearly.

`hpc_templates/bulk_le200_array.slurm` batches many families per SLURM
array task (chunked by `data/manifest/chunks_le200/`), sized against a
timing calibration sample so that most chunks finish comfortably within the
partition's time cap. `scripts/check_completion.py` identifies whichever
families didn't finish and compiles them into a smaller mop-up chunk set
for `hpc_templates/mopup_le200_array.slurm` — a family is considered done
once its `step6_family/family_summary.tsv` exists with at least one data
row.

## References

- Shimodaira (2002) An approximately unbiased test of phylogenetic tree
  selection. *Syst Biol* 51:492–508
- Siu-Ting et al. (2019) Inadvertent paralog inclusion drives artifactual
  topologies and timetree estimates in phylogenomics. *Mol Biol Evol*
  36:1344–1356
