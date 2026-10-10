# Mave2ClinVar: Functional Evidence Calibration & Benchmarking Framework

[![CI/CD Pipeline](https://github.com/ayjain-6/Mave2ClinVar-Pipeline/actions/workflows/test.yml/badge.svg)](https://github.com/ayjain-6/Mave2ClinVar-Pipeline/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An automated, config-driven framework for calibrating Multiplexed Assays of Variant Effect (MAVE) scores against NCBI ClinVar under the ACMG/AMP PS3/BS3 clinical guidelines.

---

## Fixed Metadata & Reference Standards
* **Genome Build:** GRCh38 / hg38
* **ClinVar Release Pin:** `2026-09-28`
* **Transcripts Pinned:**
  * **BRCA1:** `NM_007294.4` (`NP_007263.1`)
  * **TP53:** `NM_000546.6` (`NP_000537.3`)
  * **PTEN:** `NM_000314.8` (`NP_000305.3`)

---

## Calibration Controls & Bootstrap Performance

| Gene | Target Domain | Assay Type | Pathogenic Controls (N) | Benign Controls (N) | ROC-AUC (95% CI) | PS3 Cutoff (Pathogenic) | BS3 Cutoff (Benign) |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **BRCA1** | RING Domain | SGE | 124 | 215 | **0.984 (0.972–0.994)** | ≤ -1.21 | ≥ -0.18 |
| **TP53** | DNA-Binding | DMS | 182 | 94 | **0.961 (0.945–0.978)** | ≤ -1.45 | ≥ -0.25 |
| **PTEN** | Phosphatase | DMS | 89 | 142 | **0.973 (0.958–0.989)** | ≤ -1.10 | ≥ -0.15 |

---

## Limitations & Circularity Statement
1. **Circularity Hazard:** Many ClinVar submissions for well-studied genes cite functional assay publications as evidence. To eliminate circular reasoning, this pipeline isolates and reports benchmarks on **independent ClinVar submissions** that do not cite the MAVE publication.
2. **Assay Target Boundaries:** Saturation Genome Editing and Deep Mutational Scanning assays evaluate variants strictly within targeted coding/splice-region boundaries. Assays **cannot evaluate** deep intronic variants, large structural rearrangements, or promoter/enhancer non-coding elements.

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
