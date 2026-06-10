import pandas as pd
import numpy as np
from pathlib import Path

# ==================== CONFIGURATION ====================
VARIANT_COUNT_DIR = Path("results/variant_counts")
MAPPING_FILE = Path("results/func_scores/functional_selections.csv")
PSEUDOCOUNT = 0.5
MIN_PRE_COUNT = 20
NORMALIZE_BY_NEUTRAL = True
OUTPUT_CSV = "functional_scores_all.csv"
# =======================================================

# Read mapping
mapping = pd.read_csv(MAPPING_FILE)

# The columns we need: preselection_library_sample, postselection_library_sample, selection_name
required = ['preselection_library_sample', 'postselection_library_sample', 'selection_name']
if not all(col in mapping.columns for col in required):
    raise ValueError(f"Mapping missing required columns. Found: {mapping.columns.tolist()}")

def compute_pair(pre_file, post_file, sel_name):
    pre_path = VARIANT_COUNT_DIR / f"{pre_file}.csv"
    post_path = VARIANT_COUNT_DIR / f"{post_file}.csv"
    
    if not pre_path.exists():
        print(f"Warning: {pre_path} not found, skipping {sel_name}")
        return None
    if not post_path.exists():
        print(f"Warning: {post_path} not found, skipping {sel_name}")
        return None
    
    pre_df = pd.read_csv(pre_path)
    post_df = pd.read_csv(post_path)
    
    # Merge on barcode
    merged = pre_df.merge(post_df[['barcode', 'count']], on='barcode', suffixes=('_pre', '_post'))
    
    # Filter low pre‑selection counts
    merged = merged[merged['count_pre'] >= MIN_PRE_COUNT].copy()
    
    # Compute raw functional score (log2 enrichment)
    merged['raw_score'] = np.log2((merged['count_post'] + PSEUDOCOUNT) / (merged['count_pre'] + PSEUDOCOUNT))
    
    # Normalize by neutral standards (if requested)
    if NORMALIZE_BY_NEUTRAL:
        neut = merged[merged['aa_substitutions'] == 'neut_standard']
        if len(neut) > 0:
            neut_median = neut['raw_score'].median()
            merged['functional_score'] = merged['raw_score'] - neut_median
        else:
            print(f"Warning: no neut_standard barcodes in {sel_name}, using raw scores")
            merged['functional_score'] = merged['raw_score']
    else:
        merged['functional_score'] = merged['raw_score']
    
    merged['selection_name'] = sel_name
    merged['pre_sample'] = pre_file
    merged['post_sample'] = post_file
    return merged

# Process all pairs
all_results = []
for idx, row in mapping.iterrows():
    pre = row['preselection_library_sample']
    post = row['postselection_library_sample']
    sel = row['selection_name']
    print(f"Processing {pre} vs {post} -> {sel}")
    df = compute_pair(pre, post, sel)
    if df is not None:
        all_results.append(df)

if not all_results:
    raise ValueError("No valid pairs processed. Check file names and mapping.")

# Combine and save
final = pd.concat(all_results, ignore_index=True)

# Reorder columns for clarity
cols = ['barcode', 'selection_name', 'functional_score', 'raw_score',
        'count_pre', 'count_post', 'codon_substitutions', 'aa_substitutions',
        'variant_call_support', 'pre_sample', 'post_sample']
final = final[[c for c in cols if c in final.columns]]
final.to_csv(OUTPUT_CSV, index=False)

print(f"\nSaved {len(final)} rows to {OUTPUT_CSV}")
