# Input Data Specification

The pipeline consumes Peter Mulhair's real RHG (reticulate homolog group)
data directly — there is no separate "prepare your inputs" step. A manifest
builder (`src/build_manifest.py`) discovers every (family, domain) unit from
files already on disk and derives everything else (taxon lists, composite
sets, per-domain model strings, a per-family random seed) from them.

## Directory layout expected by `build_manifest.py`

One flat directory (e.g. `conv_fusion_data/`) containing, per family
`F#####` and per domain `N` (a family has 2 domains in the current data, but
nothing assumes exactly 2 — some families have 3+, see "Multi-domain
families" below):

```
conv_fusion_data/
├── F100058_Domain1.trim               # trimmed FASTA alignment for this domain
├── F100058_Domain1.constr0.treefile   # constrained tree (composite clade forced
│                                       #   monophyletic), branch lengths already
│                                       #   estimated
├── F100058_Domain1-geneID.tsv         # per-taxon label: "composite" or "parent"
├── F100058_Domain2.trim
├── F100058_Domain2.constr0.treefile
├── F100058_Domain2-geneID.tsv
├── ...
└── CHG_model.json                     # {"F#####_Domain{N}": "<IQTree model string>", ...}
                                        # one file, shared across all families
```

### `F#####_Domain{N}.trim`

FASTA alignment for that domain — the region of the fused gene homologous to
one parent gene, spanning both that parent's own orthologs and the
composite/fused taxa's corresponding region. Already trimmed; the pipeline
does not trim.

### `F#####_Domain{N}.constr0.treefile`

Newick tree with the composite clade constrained to be monophyletic.
Branch lengths are expected to already be present (Step 1 re-fits them
under the domain's real model regardless, via a fixed-topology `-te` IQTree
run — see CLAUDE.md — so a tree without branch lengths would also work, just
less accurately).

### `F#####_Domain{N}-geneID.tsv`

Two columns, no header: `taxon<TAB>category`, where `category` is
`composite` (this taxon carries the fused gene) or `parent` (this taxon
carries only this domain's parent gene). This is what Step 5 uses to know
which taxa should form the composite clade when checking monophyly, and
what `build_manifest.py` uses to derive each family's composite taxon set.

### `CHG_model.json`

One JSON object, shared across the whole dataset, mapping
`"F#####_Domain{N}"` → IQTree model string (e.g. `"LG+G4"`, `"JTT+F+R5"`).
`build_manifest.py` looks this up per domain; a domain with no entry here is
skipped (logged as a warning, not a fatal error — see below).

## Building the manifest

```bash
python3 src/build_manifest.py \
    --data-dir /path/to/conv_fusion_data \
    --output-json data/manifest/manifest_full.json \
    --output-tsv  data/manifest/manifest_full.tsv
```

This discovers families/domains from the `.trim` files actually present on
disk (not from `CHG_model.json`'s key set — that file can contain orphaned
entries with no corresponding alignment/tree/geneID triple; those are
skipped with a warning, not an error). For each family it also computes a
deterministic seed (`sha256(family_id)`), so all of a family's domains draw
from the same underlying simulation scenario — replicate *i* means the same
joint scenario across every domain of that family, not independently-drawn
per-domain replicates.

`manifest_full.json` (nested, keyed by family) is what `run_family.sh`
reads at runtime; `manifest_full.tsv` (flat, one row per family+domain) is
for ad hoc shell/`xargs` use, e.g. filtering by `n_taxa` to build a batching
strategy (see `data/manifest/chunks_le200/` for the ≤200-taxa bulk run).

Restrict to one family while testing with `--family F100058`.

## Multi-domain families

Most families have exactly 2 domains (one per parent gene). A minority have
more — `build_manifest.py` places no assumption on domain count, and
`run_family.sh` loops over however many domains a family's manifest entry
lists. Per an earlier project decision, extra-domain families are treated
as ordinary representatives of the family, not a special case: each domain
still gets its own Step 1–6 run, and the family-level OR/AND join (Step 6
family) folds in however many domains are present.

## Validation the manifest builder performs

- **Missing companion file**: a `.trim` with no matching
  `.constr0.treefile` or `-geneID.tsv` — that domain is skipped, family
  keeps whatever other domains are complete.
- **Missing model entry**: a domain with no `CHG_model.json` key — skipped
  the same way.
- **Composite-taxon-set consistency**: a family's domains should agree on
  which taxa are composite (a taxon carrying the fused gene should show
  `composite` in every domain's `-geneID.tsv`, not just some). Inconsistent
  families are still included but flagged
  (`composite_taxa_consistent: false` in the JSON) — worth spot-checking
  before trusting that family's family-level join.

## What downstream steps actually read

`run_family.sh` pulls each domain's `trim`, `tree`, `model`, `geneid` path
and the family's `seed` straight out of `manifest_full.json` — see its
`DOMAIN_INFO` python snippet. No other input file format is supported; if
you have data in a different shape (e.g. separate untrimmed alignments, or
a different taxon-label convention), the manifest builder is the one place
that would need extending, not the six step scripts.
