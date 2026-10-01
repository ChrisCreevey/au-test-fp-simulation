# AU Test False Positive Rate Simulation Pipeline

A computational pipeline to test whether the AU (Approximately Unbiased)
test produces false positives when detecting convergently evolved
remodelled (fused) gene families.

The pipeline simulates sequence data under the correct null model —
single-origin gene fusion, where the fused clade is truly monophyletic —
using branch lengths, per-site rate variation, and gap patterns calibrated
from real data. It then runs the AU test against the simulated replicates
and measures how often the test incorrectly rejects monophyly, giving an
empirical false-positive rate to compare against the 41% polyphyletic rate
observed in real gene families (Mulhair et al.).

See [docs/METHODS.md](docs/METHODS.md) for the full scientific background
and architecture rationale.

## Dependencies

- **IQ-TREE 2.2+** — tree inference, sequence simulation (AliSim), and the
  AU test. The only external tool the pipeline calls.
- **Python 3.8+** — data processing and summarisation (standard library
  only, no extra packages required)
- **bash**, **xargs** — replicate parallelisation

```bash
conda env create -f environment.yml && conda activate au_test_sim
# or, per-platform installers:
bash install_dependencies_local.sh   # macOS/Homebrew
bash install_dependencies_hpc.sh     # HPC/conda
```

## How input data works

There is no "prepare your inputs" step — the pipeline reads Peter
Mulhair's real RHG data directly. A manifest builder discovers every
(gene family, domain) unit from files already on disk:

```bash
python3 src/build_manifest.py \
    --data-dir /path/to/conv_fusion_data \
    --output-json data/manifest/manifest_full.json \
    --output-tsv  data/manifest/manifest_full.tsv
```

See [docs/INPUTS.md](docs/INPUTS.md) for the exact file format
(`F#####_Domain{N}.trim` / `.constr0.treefile` / `-geneID.tsv` +
`CHG_model.json`) and what the manifest builder validates.

## Running the pipeline

### One family

```bash
bash scripts/run_family.sh \
    -f F100058 \
    -m data/manifest/manifest_full.json \
    -o data/sim_outputs \
    -c config/pipeline.conf
```

Runs Steps 1–6 for every domain the family's manifest entry lists (usually
2), then the family-level OR/AND join. Pass `-r` to also run the (expensive)
real-data permutation comparison, `-x` to force a rerun of every step.

### Many families (SLURM)

```bash
# one family per line
bash scripts/run_family_chunk.sh <chunk_file> <manifest.json> <output_base_dir> config/pipeline.conf
```

`hpc_templates/bulk_le200_array.slurm` is the production array-job template
that drives this at scale (see docs/METHODS.md §6 for the chunk-sizing rationale
and current run status). `hpc_templates/run_single_family.slurm` is a
smaller template for one family. `hpc_templates/mopup_le200_array.slurm` +
`scripts/check_completion.py` handle resubmitting whatever a first pass
didn't finish within its time cap — every step skips work already done, so
a mop-up run resumes rather than redoing anything.

Running on a **different cluster** (different scheduler/partitions, or
someone else's copy of the data)? See
[docs/ADAPTING_TO_YOUR_CLUSTER.md](docs/ADAPTING_TO_YOUR_CLUSTER.md) —
covers both a one-job-per-family approach and the chunked-batch approach
above, plus sizing guidance for larger families.

### Combine results across families

```bash
bash scripts/combine_results.sh -o data/sim_outputs_bulk
```

Writes `combined_family_results.tsv` — every family's `family_summary.tsv`
joined by column-name union (families with more domains produce extra
columns; see `src/combine_family_results.py`).

## Pipeline overview

Per domain (each fused family has one AU test per parent-gene domain, not
one per family — see docs/METHODS.md §2 for why):

1. **Branch lengths + site rates** (`01_estimate_branch_lengths.sh`) — fit
   under the domain's real model on the real alignment
2. **Simulate replicates** (`02_run_alisim.sh`) — AliSim, per-site-rate
   concentrated, `N_REPLICATES` replicates in one batched call
3. **Gap masking** (`03_mask_alignments.sh`) — reapply the real alignment's
   gap pattern positionally
4. **Reconstruct trees** (`04_reconstruct_trees.sh`) — unconstrained ML +
   monophyly-constrained tree per replicate, run concurrently
5. **AU test** (`05_run_au_test.sh`) — compare the two trees, record
   p-value/ΔlogL/monophyly status for every replicate
6. **Summarise** — per-domain (`06_summarise_domain.sh`) then per-family
   OR/AND join (`06_summarise_family.sh`)

## Output structure

```
<output_base>/<family_id>/
├── Domain1/
│   ├── step1/   step2/   step3/   step4/   step5/   step6/
│   │   (step5/step5_results.txt has one RESULT line per replicate;
│   │    step6/results.tsv has the per-domain summary)
├── Domain2/
│   └── ...
├── step6_family/
│   ├── family_summary.tsv            # one row: OR/AND rates, per-domain rates
│   └── family_replicate_detail.tsv   # full per-replicate transparency
└── run_family.log
```

## Key output columns

`step6_family/family_summary.tsv`:
- `Pct_Significant_Domain{N}` — per-domain false-positive rate (AU
  rejected **and** the composite clade actually lost monophyly — not just
  any AU rejection)
- `Pct_OR_Violations` / `Pct_AND_Violations` — family-level rate if *any* /
  *all* domains show a significant violation on the same replicate
- `Real_DeltaLogL_Domain{N}` / `Permutation_P_Upper_Domain{N}` — only
  populated if `-r` was used for that domain

Under a correctly calibrated AU test, per-domain and OR/AND rates should
all be ≈5% (α = 0.05); compare against the real data's 41%.

## Citation

This simulation framework builds on:

- Shimodaira (2002) An approximately unbiased test of phylogenetic tree
  selection. *Syst Biol* 51:492–508
- Siu-Ting et al. (2019) Inadvertent paralog inclusion drives artifactual
  topologies in phylogenomics. *Mol Biol Evol* 36:1344–1356
