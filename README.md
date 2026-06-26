# HIV Envelope BF520 Deep Mutational Scanning — Functional Score Analysis

This repository reproduces and extends the functional-score arm of the deep mutational scanning (DMS) analysis of the HIV-1 Envelope protein from strain BF520. The experiment measured how every single amino-acid mutation in Env affects viral entry efficiency by comparing barcode counts before and after a functional selection (VSVG-pseudovirus vs. no-antibody control). The output is a per-mutation **functional score** (log₂ enrichment relative to wildtype) that quantifies whether a mutation is tolerated, neutral, or deleterious for viral entry.

### Original study

The data originate from:

> Radford, C.E., Bloom, J.D., et al. (2023). *Mapping the neutralization landscape of HIV-1 Env with deep mutational scanning.* bioRxiv. https://www.biorxiv.org/content/10.1101/2023.03.23.533993v1

The upstream barcode sequencing and library processing were performed with the [dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline) (Bloom lab). This repository picks up from the processed `variant_counts/` files that pipeline produces and implements a standalone functional-score calculation pipeline.

---

## Linux only

Several dependencies in this pipeline are **not available on Windows**:

- `pysam` — requires POSIX file-descriptor semantics
- `mafft` and `minimap2` — Linux/macOS binaries only
- `entrez-direct` — NCBI command-line tools, Linux/macOS only

Run this pipeline on Linux or macOS (WSL2 is fine).

---

## Environment setup

### Conda (recommended)

```bash
conda env create -f environment.yml
conda activate dms-vep-pipeline1
```

This creates an environment named `dms-vep-pipeline1` with Python 3.11, Snakemake 7, all scientific packages, and the DMS-specific libraries (`dms-variants`, `polyclonal`, `alignparse`, etc.).

It also installs the local package in editable mode (`-e .`) so that `from config import config` works from any script or notebook without path manipulation.

### pip (alternative)

If you prefer pip without conda, first install the non-pip tools (`mafft`, `minimap2`, `entrez-direct`) through your system package manager, then:

```bash
pip install -r requirements.txt
```

> Note: `requirements.txt` lists all pip-installable packages. The system-level tools (`mafft`, `minimap2`, `entrez-direct`) must still be installed separately.

---

## Configuration files

### `config.yaml`

The master configuration for the entire repository. It has two logical sections:

**Original dms-vep-pipeline parameters** (top of file): settings inherited from the Bloom lab's upstream pipeline — PacBio consensus parameters, Illumina barcode parser settings, functional-score thresholds, antibody escape parameters, PDB structures, and all result directory paths used by the original pipeline. These are kept for reference and for the `designed_mutations` path used in the notebook.

**Custom pipeline block** (bottom, under `pipeline:`): paths and parameters consumed by the Snakemake pipeline in this repository:

| Key | Default | Purpose |
|-----|---------|---------|
| `variant_counts_dir` | `data/variant_counts` | Input barcode count CSVs |
| `site_numbering_map` | `data/site_numbering/site_numbering_map.csv` | Sequential → HXB2 position map |
| `functional_selections_clean` | `results/functional_selections_clean.csv` | Output of Step 1 |
| `merged_output_dir` | `results/merged_output` | Output of Step 2 |
| `mapped_output_dir` | `results/mapped_output` | Output of Step 3 |
| `func_scores_dir` | `results/func_scores_output` | Output of Step 4 |
| `pseudocount` | `0.5` | Added to counts before log ratio |
| `min_preselection_counts` | `20` | Minimum pre-selection count per variant |
| `min_preselection_frac` | `0.000001` | Minimum fraction of total pre-selection counts |

To change any path or threshold, edit only the `pipeline:` block.

### `config.py`

A Python module that reads `config.yaml` and exposes the `pipeline:` block as a `SimpleNamespace` object named `config`. This lets notebooks and scripts import paths without hard-coding them:

```python
from config import config

func_scores_dir = config.func_scores_dir   # -> "results/func_scores_output"
site_map        = config.site_numbering_map
```

It also adds `config.designed_mutations` (from the top-level `mutation_design_classification` key) and `config.raw` (the full parsed YAML) for any values not explicitly exposed as attributes.

When imported, `config.py` changes the working directory to the repository root so that all relative paths resolve correctly regardless of where the script is invoked from.

### `pyproject.toml`

Makes the repository itself a pip-installable package (`hiv-dms-analysis`). The only reason this file exists is to support the editable install (`pip install -e .`) so that `from config import config` works as a standard import from anywhere — notebooks, scripts, and tests — without needing to manipulate `sys.path` manually. You do not need to interact with this file directly.

---

## Pipeline workflow

The pipeline is orchestrated by `Snakefile` and has five sequential steps. Each step is a Python script in `src/`.

```
data/variant_counts/          ← input: per-sample barcode count CSVs
        │
        ▼ Step 1 — CleanPairer.py
results/functional_selections_clean.csv
        │
        ▼ Step 2 — variantMerger.py
results/merged_output/        ← one merged CSV per selection pair
        │
        ▼ Step 3 — MutationMapper.py
results/mapped_output/        ← positions remapped to HXB2 numbering
        │
        ▼ Step 4 — FunctionalCalculator.py
results/func_scores_output/   ← functional scores per variant
        │
        ▼ Step 5 — notebooks/functional_score_analysis.ipynb (via papermill)
results/functional_score_analysis.ipynb   ← executed notebook with all plots
```

### `src/CleanPairer.py` — Step 1

Scans `data/variant_counts/` and pairs VSVG-control files (pre-selection) with no-antibody-control files (post-selection) by matching on library, date, virus batch, and replicate extracted from the filename format `{library}_{date}_{batch}_{type}_{replicate}.csv`.

**Output:** `results/functional_selections_clean.csv` — a table where each row is one matched selection pair with columns for `preselection_sample`, `postselection_sample`, `library`, `virus_batch`, `replicate`, and `selection_name`.

### `src/variantMerger.py` — Step 2

Reads the selections CSV from Step 1 and merges each pre- and post-selection count file into a single DataFrame per selection. Validates required columns, strips neutralization-standard barcodes, and writes one merged CSV per selection.

**Output:** `results/merged_output/` — one `{selection_name}_merged.csv` per selection pair containing joined pre- and post-selection counts alongside barcode and variant metadata.

### `src/MutationMapper.py` — Step 3

Reads the merged files from Step 2 and remaps amino-acid substitution positions from sequential (1, 2, 3 …) numbering to the HXB2 HIV reference numbering scheme using `data/site_numbering/site_numbering_map.csv`. Also annotates each variant with `n_aa_substitutions` and `n_codon_substitutions` counts.

**Output:** `results/mapped_output/` — one `{selection_name}_mapped.csv` per file with substitutions written in HXB2 coordinates and a new `aa_substitutions_reference` column alongside the original sequential-coordinate column.

### `src/FunctionalCalculator.py` — Step 4

Computes a **functional score** for every variant in each selection using:

```
func_score = log2( (post_count + pseudocount) / post_wt )
           − log2( (pre_count  + pseudocount) / pre_wt  )
```

Wildtype is identified as rows with `n_codon_substitutions == 0`. Variants that fall below `min_preselection_counts` or `min_preselection_frac` are flagged with a threshold column rather than dropped, so downstream analysis can apply its own filters.

**Output:** `results/func_scores_output/` — one `{selection_name}_merged_mapped_func_score.csv` per selection with columns for pre/post counts, thresholds, and the computed `func_score`.

---

## Analysis notebook

`notebooks/functional_score_analysis.ipynb` is the final pipeline step, executed automatically by Snakemake via [papermill](https://papermill.readthedocs.io/). It performs QC and reproduciblity analysis on the functional scores:

1. **Load pipeline output** — reads all `func_score` CSVs and concatenates them into a single DataFrame.
2. **Pre-selection count QC** — box plots of variant count distributions per sample, verifying that median counts exceed the minimum threshold for all samples.
3. **Fraction above threshold** — bar plots showing what fraction of counts and variants pass the count threshold per sample.
4. **Replicate count correlation** — scatter plots of log₂ barcode counts between replicate 1 and replicate 2 for both VSVG-control and no-antibody-control samples, confirming sequencing reproducibility before score calculation.
5. **Variant class filtering** — filters to variants meeting the pre-count threshold and classifies each as wildtype, synonymous, 1-nonsynonymous, 2+-nonsynonymous, or indel.
6. **Unintended-mutation flagging** — cross-references each variant against the library design to flag substitutions that were not intentionally included.
7. **Functional scores by variant class** — box plot comparing score distributions across variant classes; wildtype/synonymous variants should cluster near zero.
8. **Per-experiment violin plots** — functional score distributions broken out by individual rescue experiment.
9. **Replicate functional score correlation** — scatter plots comparing mutation-level scores between replicates within each library/batch, with Pearson r.
10. **Enrich2 correction** — applies an iterative random-effects shrinkage (the Enrich2 model) to each replicate's scores and repeats the replicate correlation plot to show how the correction improves agreement.
11. **Cross-batch correlation (raw)** — compares mutation scores across different virus rescue batches to assess inter-batch reproducibility.
12. **Cross-batch correlation (Enrich2-corrected)** — same cross-batch comparison after applying the Enrich2 correction.

The executed notebook is saved to `results/functional_score_analysis.html` (rendered) and `results/functional_score_analysis.ipynb` (with outputs).

---

## Running the pipeline with Snakemake

Make sure the conda environment is active and you are in the repository root.

**Full pipeline (all 5 steps):**

```bash
snakemake --cores 1
```

**Dry run (preview what will run without executing):**

```bash
snakemake --cores 1 --dry-run
```

**Run a specific step only:**

```bash
snakemake --cores 1 results/functional_selections_clean.csv   # Step 1
snakemake --cores 1 results/merged_output                     # Steps 1–2
snakemake --cores 1 results/mapped_output                     # Steps 1–3
snakemake --cores 1 results/func_scores_output                # Steps 1–4
```

**Force re-run of all steps (ignore existing outputs):**

```bash
snakemake --cores 1 --forceall
```

**Generate a pipeline graph (requires graphviz):**

```bash
snakemake --filegraph | dot -Tsvg > filegraph.svg
snakemake --rulegraph | dot -Tsvg > rulegraph.svg
```

The final output of the pipeline is `results/functional_score_analysis.ipynb` (the executed analysis notebook). An HTML render is also written alongside it.

---

## Repository layout

```
.
├── config.py                   # Python config loader (imports config.yaml)
├── config.yaml                 # Master configuration
├── pyproject.toml              # Makes repo pip-installable (needed for config import)
├── environment.yml             # Conda environment
├── requirements.txt            # pip requirements
├── Snakefile                   # Pipeline orchestration
├── data/
│   ├── variant_counts/         # Input: per-sample barcode count CSVs
│   └── site_numbering/         # HXB2 position mapping file
├── src/
│   ├── CleanPairer.py          # Step 1: pair pre/post selections
│   ├── variantMerger.py        # Step 2: merge count files
│   ├── MutationMapper.py       # Step 3: remap to HXB2 numbering
│   └── FunctionalCalculator.py # Step 4: compute functional scores
├── notebooks/
│   └── functional_score_analysis.ipynb  # Step 5: QC and analysis
└── results/                    # All pipeline outputs (git-ignored)
```
