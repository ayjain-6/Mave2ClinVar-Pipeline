# Mave2ClinVar: Functional Evidence Calibration & Benchmarking Framework

[![CI/CD Pipeline](https://github.com/ayjain-6/Mave2ClinVar-Pipeline/actions/workflows/test.yml/badge.svg)](https://github.com/ayjain-6/Mave2ClinVar-Pipeline/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An automated, config-driven framework for calibrating Multiplexed Assays of Variant Effect (MAVE) scores against NCBI ClinVar under the ACMG/AMP PS3/BS3 clinical guidelines using count-based OddsPath likelihood ratios.

---

## Fixed Metadata & Reference Standards
* **Genome Build:** GRCh38 / hg38
* **ClinVar Release Pin:** `2026-09-28`
* **Prior Pathogenicity (P1):** `0.10`

---

## Calibration Controls & OddsPath Performance (`results/summary.csv`)

| Gene | Subset | Pathogenic (N) | Benign (N) | ROC-AUC (95% CI) | PS3 Cutoff | OddsPath (Abnormal) | ACMG Strength |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **BRCA1** | Missense Only | 112 | 198 | **0.979 (0.965–0.991)** | ≤ -1.21 | 52.4 | PS3_Strong |
| **BRCA1** | All Controls | 124 | 215 | **0.984 (0.972–0.994)** | ≤ -1.21 | 61.1 | PS3_Strong |
| **TP53** | Missense Only | 165 | 82 | **0.954 (0.938–0.972)** | ≤ -1.45 | 41.0 | PS3_Strong |
| **PTEN** | Missense Only | 81 | 130 | **0.968 (0.951–0.985)** | ≤ -1.10 | 35.8 | PS3_Strong |

---

## Limitations & Circularity Statement
1. **Circularity Hazard:** Many ClinVar submissions for well-studied genes cite functional assay publications as evidence. To maintain rigor, this pipeline **reduces circular reasoning** by auditing and filtering out independent submissions that cite gene-specific MAVE publications.
2. **Assay Target Boundaries:** Assays evaluate variants strictly within targeted coding/splice regions and cannot evaluate deep intronic or structural rearrangements.

---

## Quick Start
```bash
# Clone Repository
git clone [https://github.com/ayjain-6/Mave2ClinVar-Pipeline.git](https://github.com/ayjain-6/Mave2ClinVar-Pipeline.git)
cd Mave2ClinVar-Pipeline

# Install Dependencies
pip install -r requirements.txt

# Run Unit Tests
pytest

# Launch Streamlit Dashboard
streamlit run app.py
```
