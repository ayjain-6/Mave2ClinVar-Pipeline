import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.metrics import roc_auc_score

def load_and_filter_clinvar(clinvar_df: pd.DataFrame, mave_pmid: str = "PMID:32494009") -> pd.DataFrame:
    """Isolates independent ClinVar submissions to eliminate circular reasoning."""
    if "citations" in clinvar_df.columns:
        independent_df = clinvar_df[~clinvar_df["citations"].str.contains(mave_pmid, na=False)].copy()
    else:
        independent_df = clinvar_df.copy()
    return independent_df

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

def run_bootstrap_auc(scores: np.ndarray, labels: np.ndarray, n_iterations: int = 1000) -> tuple:
    """Computes non-parametric bootstrap confidence intervals for ROC-AUC."""
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
