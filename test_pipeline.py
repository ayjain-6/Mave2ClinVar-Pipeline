import pytest
import numpy as np
import pandas as pd
from pipeline_engine import compute_count_based_oddspath, run_bootstrap_auc, load_and_filter_clinvar

def test_pmid_exact_match():
    # Test that exact PMID boundary matching works and ignores overlapping substrings
    df = pd.DataFrame({"citations": ["PMID:32494009", "1324940099", ""]})
    out, audit = load_and_filter_clinvar(df, ["32494009"])
    assert audit["filtered_submissions"] == 1

def test_oddspath_calculation():
    odds, tier = compute_count_based_oddspath(9, 1, range_type="abnormal", p1=0.10, min_controls=5)
    assert odds > 18.7
    assert tier == "PS3_Strong"

def test_oddspath_insufficient_controls():
    odds, tier = compute_count_based_oddspath(2, 0, range_type="abnormal", p1=0.10, min_controls=11)
    assert odds is None
    assert "Insufficient Controls" in tier

def test_bootstrap_auc_orientation():
    scores = np.array([-2.5, -2.0, 0.5, 1.0])
    labels = np.array([1, 1, 0, 0])
    auc, low, high = run_bootstrap_auc(scores, labels, lower_is_pathogenic=True, n_iterations=100)
    assert auc > 0.85
