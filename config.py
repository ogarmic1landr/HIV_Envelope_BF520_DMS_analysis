import sys
import os
from pathlib import Path
from types import SimpleNamespace
import yaml

_root = Path(__file__).resolve().parent
os.chdir(_root)

_src = str(_root / "src")
if _src not in sys.path:
    sys.path.insert(0, _src)

with open(_root / "config.yaml") as f:
    _raw = yaml.safe_load(f)

_p = _raw["pipeline"]

config = SimpleNamespace(
    variant_counts_dir    = _p["variant_counts_dir"],
    site_numbering_map    = _p["site_numbering_map"],
    functional_selections = _p["functional_selections_clean"],
    merged_output_dir     = _p["merged_output_dir"],
    mapped_output_dir     = _p["mapped_output_dir"],
    func_scores_dir       = _p["func_scores_dir"],
    designed_mutations    = _raw["mutation_design_classification"],
    pipeline              = _p,
    raw                   = _raw,
)
