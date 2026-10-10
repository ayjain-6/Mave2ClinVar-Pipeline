import numpy as np
import pandas as pd
from scipy.stats import norm

def compute_oddspath_ratio(score: float, path_scores: np.ndarray, benign_scores: np.ndarray) -> float:
    mu_p, std_p = norm.fit(path_scores)
    mu_b, std_b = norm.fit(benign_scores)
    prob_path = norm.pdf(score, mu_p, std_p)
    prob_benign = norm.pdf(score, mu_b, std_b)
    return float(prob_path / prob_benign) if prob_benign > 0 else 999.0

def test_oddspath_arithmetic():
    path_scores = np.array([-2.5, -2.1, -1.9, -2.3, -2.0])
    benign_scores = np.array([-0.1, 0.0, 0.2, -0.2, 0.1])
    
    odds_pathogenic = compute_oddspath_ratio(-2.2, path_scores, benign_scores)
    odds_benign = compute_oddspath_ratio(0.0, path_scores, benign_scores)
    
    assert odds_pathogenic > 10.0
    assert odds_benign < 0.1

def test_hgvs_protein_parsing():
    sample_names = pd.Series([
        "NM_007294.4(BRCA1):c.181T>G (p.Cys61Gly)",
        "NM_000546.6(TP53):c.524G>A (p.Arg175His)"
    ])
    parsed = sample_names.str.extract(r"\((p\.[^\)]+)\)")[0]
    assert parsed[0] == "p.Cys61Gly"
    assert parsed[1] == "p.Arg175His"
