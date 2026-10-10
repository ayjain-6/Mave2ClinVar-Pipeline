import os
import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt

def load_and_filter_clinvar(clinvar_df: pd.DataFrame, mave_pmid: str = "PMID:32494009") -> tuple:
    """
    Reduces circular reasoning by filtering out ClinVar submissions 
    that explicitly cite the core MAVE publication.
    Reports total vs. remaining independent counts.
    """
    total_initial = len(clinvar_df)
    if "citations" in clinvar_df.columns:
        independent_df = clinvar_df[~clinvar_df["citations"].str.contains(mave_pmid, na=False)].copy()
    else:
        independent_df = clinvar_df.copy()
    
    retained_count = len(independent_df)
    filtered_out = total_initial - retained_count
    print(f"[Circularity Audit] Initial submissions: {total_initial} | Filtered (cited MAVE): {filtered_out} | Retained Independent: {retained_count}")
    
    return independent_df, retained_count

def separate_variant_types(df: pd.DataFrame) -> tuple:
    """
    Isolates missense-only variants from nonsense/synonymous/intronic controls
    to prevent inflation of AUC metrics.
    """
    if "variant_type" in df.columns:
        missense_df = df[df["variant_type"] == "missense"].copy()
        other_df = df[df["variant_type"] != "missense"].copy()
        return missense_df, other_df
    return df, pd.DataFrame()

def compute_oddspath_ratio(score: float, path_scores: np.ndarray, benign_scores: np.ndarray) -> float:
    """Calculates OddsPath likelihood ratio based on fitted normal distributions."""
    if len(path_scores) < 2 or len(benign_scores) < 2:
        return 1.0
    mu_p, std_p = norm.fit(path_scores)
    mu_b, std_b = norm.fit(benign_scores)
    
    std_p = max(std_p, 1e-5)
    std_b = max(std_b, 1e-5)
    
    prob_path = norm.pdf(score, mu_p, std_p)
    prob_benign = norm.pdf(score, mu_b, std_b)
    
    return float(prob_path / prob_benign) if prob_benign > 0 else 999.0

def map_acmg_strength(odds: float) -> str:
    """
    Maps OddsPath likelihood ratios to ACMG/AMP PS3/BS3 strength tiers 
    following Brnich et al. / Tavtigian framework.
    """
    if odds >= 18.7:
        return "PS3_Strong"
    elif odds >= 4.3:
        return "PS3_Moderate"
    elif odds >= 2.1:
        return "PS3_Supporting"
    elif odds <= 0.05:
        return "BS3_Strong"
    elif odds <= 0.23:
        return "BS3_Moderate"
    elif odds <= 0.48:
        return "BS3_Supporting"
    else:
        return "Indeterminate"

def run_bootstrap_auc(scores: np.ndarray, labels: np.ndarray, n_iterations: int = 1000) -> tuple:
    """Computes non-parametric bootstrap confidence intervals for ROC-AUC (B = 1,000)."""
    np.random.seed(42)
    auc_scores = []
    n = len(scores)
    
    for _ in range(n_iterations):
        indices = np.random.choice(n, size=n, replace=True)
        sample_scores = scores[indices]
        sample_labels = labels[indices]
        
        if len(np.unique(sample_labels)) < 2:
            continue
            
        score_auc = roc_auc_score(sample_labels, sample_scores)
        auc_scores.append(score_auc)
        
    lower_ci = np.percentile(auc_scores, 2.5)
    upper_ci = np.percentile(auc_scores, 97.5)
    mean_auc = np.mean(auc_scores)
    
    return mean_auc, lower_ci, upper_ci
