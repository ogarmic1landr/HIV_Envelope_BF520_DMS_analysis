configfile: "config.yaml"

# ---------------------------------------------------------------------------
# Paths — all relative to the repo root, defined in config.yaml
# ---------------------------------------------------------------------------
VARIANT_COUNTS_DIR     = config["pipeline"]["variant_counts_dir"]
SITE_NUMBERING_MAP     = config["pipeline"]["site_numbering_map"]
FUNCTIONAL_SELECTIONS  = config["pipeline"]["functional_selections_clean"]
MERGED_OUTPUT_DIR      = config["pipeline"]["merged_output_dir"]
MAPPED_OUTPUT_DIR      = config["pipeline"]["mapped_output_dir"]
FUNC_SCORES_DIR        = config["pipeline"]["func_scores_dir"]

PSEUDOCOUNT            = config["pipeline"].get("pseudocount", 0.5)
MIN_PRE_COUNTS         = config["pipeline"].get("min_preselection_counts", 20)
MIN_PRE_FRAC           = config["pipeline"].get("min_preselection_frac", 1e-6)

ANALYSIS_NOTEBOOK_IN   = "functional_score_analysis.ipynb"
ANALYSIS_NOTEBOOK_OUT  = "results/functional_score_analysis.ipynb"


# ---------------------------------------------------------------------------
# Target rule — running `snakemake` with no arguments builds the final output
# ---------------------------------------------------------------------------
rule all:
    input:
        ANALYSIS_NOTEBOOK_OUT


# ---------------------------------------------------------------------------
# Step 1: Scan variant_counts/ and generate functional selection pairs
# ---------------------------------------------------------------------------
rule clean_pairs:
    input:
        VARIANT_COUNTS_DIR
    output:
        FUNCTIONAL_SELECTIONS
    shell:
        "python CleanPairer.py "
        "--variant-counts-dir {input} "
        "--output-csv {output}"


# ---------------------------------------------------------------------------
# Step 2: Merge pre- and post-selection variant counts for each selection
# ---------------------------------------------------------------------------
rule merge_variants:
    input:
        selections  = FUNCTIONAL_SELECTIONS,
        counts_dir  = VARIANT_COUNTS_DIR
    output:
        directory(MERGED_OUTPUT_DIR)
    shell:
        "python variantMerger.py "
        "--selections {input.selections} "
        "--variant-counts-dir {input.counts_dir} "
        "--output-dir {output}"


# ---------------------------------------------------------------------------
# Step 3: Remap amino acid positions to HXB2 reference numbering
# ---------------------------------------------------------------------------
rule map_mutations:
    input:
        merged_dir = MERGED_OUTPUT_DIR,
        map_file   = SITE_NUMBERING_MAP
    output:
        directory(MAPPED_OUTPUT_DIR)
    shell:
        "python MutationMapper.py "
        "--input-dir {input.merged_dir} "
        "--map-file {input.map_file} "
        "--output-dir {output}"


# ---------------------------------------------------------------------------
# Step 4: Compute functional scores (log2 enrichment relative to WT)
#         Produces one *_merged_mapped_pipeline.csv per input file
# ---------------------------------------------------------------------------
rule compute_func_scores:
    input:
        mapped_dir = MAPPED_OUTPUT_DIR,
        selections  = FUNCTIONAL_SELECTIONS
    output:
        directory(FUNC_SCORES_DIR)
    params:
        pseudocount  = PSEUDOCOUNT,
        min_counts   = MIN_PRE_COUNTS,
        min_frac     = MIN_PRE_FRAC
    shell:
        "python FunctionalCalculator.py "
        "--input-dir {input.mapped_dir} "
        "--output-dir {output} "
        "--selections {input.selections} "
        "--pseudocount {params.pseudocount} "
        "--min-preselection-counts {params.min_counts} "
        "--min-preselection-frac {params.min_frac}"


# ---------------------------------------------------------------------------
# Step 5: Run analysis notebook
# ---------------------------------------------------------------------------
rule run_analysis_notebook:
    input:
        notebook   = ANALYSIS_NOTEBOOK_IN,
        func_scores = FUNC_SCORES_DIR,
        selections  = FUNCTIONAL_SELECTIONS
    output:
        ANALYSIS_NOTEBOOK_OUT
    shell:
        "papermill {input.notebook} {output}"
