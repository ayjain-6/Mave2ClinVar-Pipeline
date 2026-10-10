import os
import yaml
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

def load_config(config_path="config.yaml"):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def load_and_filter_clinvar(clinvar_df: pd.DataFrame, gene_pmid: str) -> tuple:
    """
    Reduces circular reasoning by auditing and filtering out ClinVar submissions 
    that explicitly cite the gene-specific core MAVE publication.
    Raises a warning if no citation column is present.
    """
    total_initial = len(clinvar_df)
    if "citations" in clinvar_df.columns:
        independent_df = clinvar_df[~clinvar_df["citations"].str.contains(gene_pmid, na=False)].copy()
        retained_count = len(independent_df)
        filtered_out = total_initial - retained_count
        print(f"[Circularity Audit] PMID {gene_pmid} | Initial: {total_initial} | Filtered: {filtered_out} | Retained: {retained_count}")
    else:
        print("[Warning] 'citations' column missing from ClinVar data. Returning unfiltered set.")
        independent_df = clinvar_df.copy()
        retained_count = total_initial
        
    return independent_df, retained_count

def compute_count_based_oddspath(n_path_in_range: int, n_benign_in_range: int, p1: float = 0.10) -> float:
    """
    Computes true count-based OddsPath likelihood ratio (Brnich et al. / Tavtigian framework).
    OddsPath = [P2 * (1 - P1)] / [(1 - P2) * P1]
    where P2 is the observed fraction of pathogenic variants in the assay bin.
    """
    total_in_bin = n_path_in_range + n_benign_in_range
    if total_in_bin == 0:
        return 1.0
    
    p2 = n_path_in_range / total_in_bin
    if p2 >= 1.0:
        return float("inf")
    if p2 <= 0.0:
        return 0.0
        
    odds_path = (p2 * (1.0 - p1)) / ((1.0 - p2) * p1)
    return float(odds_path)

def map_oddspath_strength(odds: float, is_pathogenic_tier: bool = True) -> str:
    """Maps count-based OddsPath values to ACMG/AMP strength tiers."""
    if is_pathogenic_tier:
        if odds >= 18.7: return "PS3_Strong"
        if odds >= 4.3: return "PS3_Moderate"
        if odds >= 2.1: return "PS3_Supporting"
    else:
        if odds <= 0.053: return "BS3_Strong"
        if odds <= 0.23: return "BS3_Moderate"
        if odds <= 0.48: return "BS3_Supporting"
    return "Indeterminate"

def run_bootstrap_auc(scores: np.ndarray, labels: np.ndarray, n_iterations: int = 1000) -> tuple:
    """
    Computes point estimate ROC-AUC on full sample and non-parametric 
    bootstrap percentile confidence intervals using modern NumPy default_rng.
    """
    rng = np.random.default_rng(42)
    
    # Invert scores if lower functional score indicates pathogenicity
    # to ensure AUC is properly oriented [0.5, 1.0]
    eval_scores = -scores
    
    point_auc = roc_auc_score(labels, eval_scores)
    n = len(eval_scores)
    auc_scores = []
    
    for _ in range(n_iterations):
        indices = rng.choice(n, size=n, replace=True)
        sample_scores = eval_scores[indices]
        sample_labels = labels[indices]
        
        if len(np.unique(sample_labels)) < 2:
            continue
            
        auc_scores.append(roc_auc_score(sample_labels, sample_scores))
        
    if not auc_scores:
        return point_auc, point_auc, point_auc
        
    lower_ci = np.percentile(auc_scores, 2.5)
    upper_ci = np.percentile(auc_scores, 97.5)
    
    return point_auc, lower_ci, upper_ci

if __name__ == "__main__":
    print("Mave2ClinVar Engine Initialized.")
