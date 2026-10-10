import os
import re
import yaml
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

def load_config(config_path="config.yaml"):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def load_and_filter_clinvar(clinvar_df: pd.DataFrame, gene_pmids: list) -> tuple:
    """
    Rigorously audits and filters out ClinVar submissions that cite 
    gene-specific core MAVE publications using exact boundary matching.
    """
    total_initial = len(clinvar_df)
    if "citations" in clinvar_df.columns:
        pattern = '|'.join([r'\\b' + re.escape(pmid) + r'\\b' for pmid in gene_pmids])
        mask = clinvar_df["citations"].str.contains(pattern, na=False, regex=True)
        independent_df = clinvar_df[~mask].copy()
        retained_count = len(independent_df)
        filtered_out = total_initial - retained_count
    else:
        independent_df = clinvar_df.copy()
        retained_count = total_initial
        filtered_out = 0
        
    audit_info = {
        "initial_submissions": total_initial,
        "filtered_submissions": filtered_out,
        "retained_submissions": retained_count
    }
    return independent_df, audit_info

def separate_variant_types(df: pd.DataFrame) -> tuple:
    """Isolates missense-only variants to prevent metric inflation."""
    if "variant_type" in df.columns:
        missense_df = df[df["variant_type"] == "missense"].copy()
        other_df = df[df["variant_type"] != "missense"].copy()
        return missense_df, other_df
    return df, pd.DataFrame()

def compute_count_based_oddspath(n_path_in_range: int, n_benign_in_range: int, range_type: str = "abnormal", p1: float = 0.10, min_controls: int = 11) -> tuple:
    """
    Computes count-based OddsPath likelihood ratio with Laplace smoothing and small-bin safeguards.
    Evaluates abnormal range against PS3 tiers and normal range against BS3 tiers.
    """
    total_in_bin = n_path_in_range + n_benign_in_range
    if total_in_bin < min_controls:
        return None, "Insufficient Controls (<11)"
    
    # Laplace smoothing (add-half/add-one stabilization)
    p2 = (n_path_in_range + 0.5) / (total_in_bin + 1.0)
    
    if p2 <= 0.0:
        p2 = 1e-4
    if p2 >= 1.0:
        p2 = 1.0 - 1e-4
        
    odds_path = (p2 * (1.0 - p1)) / ((1.0 - p2) * p1)
    
    tier = "Indeterminate"
    if range_type == "abnormal":
        if odds_path >= 18.7: tier = "PS3_Strong"
        elif odds_path >= 4.3: tier = "PS3_Moderate"
        elif odds_path >= 2.1: tier = "PS3_Supporting"
    elif range_type == "normal":
        if odds_path <= 0.053: tier = "BS3_Strong"
        elif odds_path <= 0.23: tier = "BS3_Moderate"
        elif odds_path <= 0.48: tier = "BS3_Supporting"
        
    return float(odds_path), tier

def run_bootstrap_auc(scores, labels, lower_is_pathogenic: bool = True, n_iterations: int = 1000) -> tuple:
    """
    Computes point estimate ROC-AUC and bootstrap percentile confidence intervals 
    using modern NumPy default_rng with safe array casting.
    """
    scores_arr = np.asarray(scores, dtype=float)
    labels_arr = np.asarray(labels, dtype=int)
    
    eval_scores = -scores_arr if lower_is_pathogenic else scores_arr
    
    point_auc = roc_auc_score(labels_arr, eval_scores)
    n = len(eval_scores)
    rng = np.random.default_rng(42)
    auc_scores = []
    
    for _ in range(n_iterations):
        indices = rng.choice(n, size=n, replace=True)
        sample_scores = eval_scores[indices]
        sample_labels = labels_arr[indices]
        
        if len(np.unique(sample_labels)) < 2:
            continue
            
        auc_scores.append(roc_auc_score(sample_labels, sample_scores))
        
    if not auc_scores:
        return point_auc, point_auc, point_auc
        
    lower_ci = float(np.percentile(auc_scores, 2.5))
    upper_ci = float(np.percentile(auc_scores, 97.5))
    
    return float(point_auc), lower_ci, upper_ci

def run_calibration(data_dict: dict, config: dict) -> pd.DataFrame:
    """
    Executes orientation-aware calibration across dual ranges (abnormal & normal), writing summary outputs.
    """
    results = []
    p1 = config.get("calibration", {}).get("prior_pathogenic", 0.10)
    min_ctrls = config.get("calibration", {}).get("min_controls_required", 11)
    
    for gene, df in data_dict.items():
        gene_cfg = config["genes"].get(gene, {})
        pmids = gene_cfg.get("mave_pmids", [gene_cfg.get("mave_pmid", "")])
        lower_path = gene_cfg.get("lower_is_pathogenic", True)
        
        filtered_df, audit = load_and_filter_clinvar(df, pmids)
        missense_df, _ = separate_variant_types(filtered_df)
        
        subsets = {"All Controls": filtered_df, "Missense Only": missense_df}
        for subset_name, sub_df in subsets.items():
            if len(sub_df) < 5:
                continue
                
            scores = sub_df["functional_score"].values
            labels = sub_df["clinvar_class"].values
            
            n_path = int(np.sum(labels == 1))
            n_benign = int(np.sum(labels == 0))
            
            auc, ci_low, ci_high = run_bootstrap_auc(scores, labels, lower_is_pathogenic=lower_path)
            
            ps3_cut = gene_cfg.get("ps3_cutoff", -1.21)
            bs3_cut = gene_cfg.get("bs3_cutoff", -0.18)
            
            # Dynamic binning masks based on orientation (lower_is_pathogenic)
            if lower_path:
                abnormal_mask = (scores <= ps3_cut)
                normal_mask = (scores >= bs3_cut)
            else:
                abnormal_mask = (scores >= ps3_cut)
                normal_mask = (scores <= bs3_cut)
                
            abn_path = int(np.sum(abnormal_mask & (labels == 1)))
            abn_ben = int(np.sum(abnormal_mask & (labels == 0)))
            
            norm_path = int(np.sum(normal_mask & (labels == 1)))
            norm_ben = int(np.sum(normal_mask & (labels == 0)))
            
            odds_abn, tier_abn = compute_count_based_oddspath(abn_path, abn_ben, range_type="abnormal", p1=p1, min_controls=min_ctrls)
            odds_norm, tier_norm = compute_count_based_oddspath(norm_path, norm_ben, range_type="normal", p1=p1, min_controls=min_ctrls)
            
            results.append({
                "Gene": gene,
                "Subset": subset_name,
                "Pathogenic_N": n_path,
                "Benign_N": n_benign,
                "ROC_AUC": round(auc, 3),
                "CI_Lower": round(ci_low, 3),
                "CI_Upper": round(ci_high, 3),
                "PS3_Cutoff": ps3_cut,
                "BS3_Cutoff": bs3_cut,
                "OddsPath_Abnormal": round(odds_abn, 2) if odds_abn is not None else "N/A",
                "ACMG_PS3_Strength": tier_abn,
                "OddsPath_Normal": round(odds_norm, 2) if odds_norm is not None else "N/A",
                "ACMG_BS3_Strength": tier_norm,
                "Filtered_Submissions": audit["filtered_submissions"]
            })
            
    res_df = pd.DataFrame(results)
    os.makedirs("results", exist_ok=True)
    res_df.to_csv("results/summary.csv", index=False)
    return res_df
