# scBioKit

[![Python](https://img.shields.io/badge/Python-3.9%2B-3776ab?logo=python&logoColor=white)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.3%2B-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![FastAPI](https://img.shields.io/badge/FastAPI-0.110%2B-009688?logo=fastapi&logoColor=white)](https://fastapi.tiangolo.com/)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)
[![Paper](https://img.shields.io/badge/Nature%20Biomed%20Eng-2025-blue)](https://doi.org/10.1038/s41551-025-01423-7)

A browser-based platform for single-cell RNA-seq analysis, disease trajectory modeling, and computational drug discovery. Five analysis modules are chained into one configurable pipeline — no coding required.

Built on FastAPI with real-time WebSocket log streaming. Wraps Seurat (R), UNAGI, scTenifoldKnk, and UniMol behind a unified web interface.

---

## Modules

### 1. Single-Cell QC

Seurat-based preprocessing. Filters cells by gene count and mitochondrial fraction, log-normalizes, selects highly variable genes, runs PCA, builds a neighbor graph, performs Leiden clustering, and projects to UMAP.

**Input:** Seurat object (`.qs` / `.rds` / `.rda`)  
**Output:** QC violin plots, UMAP (clusters + group), marker heatmap, summary CSV  
**Key parameters:** `min_features`, `max_features`, `max_mito_pct`, `resolution`, `n_hvg`, `pca_dims`

Optional doublet detection via scDblFinder.

---

### 2. Virtual Gene Knockout

Uses scTenifoldKnk to build a sparse gene co-expression tensor from subsampled cells, then simulates removal of one or more target genes. Differential manifold components identify downstream transcriptional responders.

**Input:** Seurat object, knockout gene list (comma-separated)  
**Output:** `ko_comparison.png`, `volcano_<gene>.png` per target, `ko_results.csv`  
**Key parameters:** `ko_genes`, `cell_type_col`, `n_hvg`, `n_cells`, `n_networks`

---

### 3. Compound Screening

Fine-tunes a UniMol pre-trained molecular Transformer (84M or 164M parameters) on a labeled activity dataset, then applies the classifier to a compound library. SMILES are validated via RDKit before prediction; salt forms are discarded.

**Training:** Binary classification with 5-fold CV and early stopping. Optional external test-set evaluation produces ROC curve, confusion matrix, and probability distribution plots.  
**Screening:** Scores the full library against a user-defined activity threshold and outputs ranked hits.

| Mode | Input | Output |
|------|-------|--------|
| Train | SMILES + label CSV | `final_model/`, `cv_auc_bar.png`, `roc_curve.png` |
| Screen | SMILES library CSV + trained model dir | `screen_top_actives.csv`, `score_distribution.png`, `top_actives_bar.png` |

Both modes include a **Demo toggle** that runs on built-in example data with no file upload needed.

---

### 4. Trajectory Inference (UNAGI)

Trains a Variational Autoencoder — optionally with a Graph Convolutional Network — on per-stage `.h5ad` files to learn a continuous latent representation of disease progression. Cell embeddings are passed to iDREM to reconstruct dynamic gene regulatory networks: trajectory branches, active transcription factors, and stage-resolved marker genes.

**Input:** `0.h5ad`, `1.h5ad`, … (one per disease stage) + iDREM directory  
**Output:** `<task>.pth`, `attribute.pkl`, `dataset.h5ad`, `umap_stage.png`, `umap_celltype.png`  
**Key parameters:** `distribution` (`ziln` / `normal`), `max_iter`, `batch_size`, `n_neighbors`, `device`

> The trained model directory feeds directly into the Drug Perturbation module.

---

### 5. Drug Perturbation

Reads drug-specific differential expression signatures (one CSV per compound with gene, log₂FC, and p-value columns) and applies each to the UNAGI latent space. A reversal score is computed for each drug based on how far cell embeddings shift toward the healthy-state trajectory defined by iDREM tendency labels. An optional anti-disease filter retains only genes whose direction is consistent with reversing the disease trajectory.

**Input:** Trained UNAGI model (`.pth`), `attribute.pkl`, drug DEG directory  
**Output:** `drug_scores.csv` (ranked), `drug_perturbation_scores.png`  
**Key parameters:** `log2fc`, `background_samples`, `tend_thresh`, `anti_copd_filter`

---

## Installation

### 1. Clone

```bash
git clone https://github.com/<your-username>/scBioKit.git
cd scBioKit
```

### 2. Python environment

```bash
conda create -n scbiokit python=3.10
conda activate scbiokit

pip install fastapi "uvicorn[standard]" python-multipart websockets
pip install scanpy anndata torch pyro-ppl
pip install scUNAGI
pip install unimol_tools
pip install rdkit
```

### 3. R packages

```r
install.packages(c("Seurat", "qs", "scDblFinder"))

install.packages("devtools")
devtools::install_github("cailab-tamu/scTenifoldKnk")
```

### 4. iDREM

Requires Java 1.7+ (64-bit) on PATH.

```bash
git clone https://github.com/phoenixding/idrem.git
```

Pass the path to the cloned directory as `idrem_dir` in the Trajectory Inference module.

### 5. Drug gene signatures

Download the preprocessed CMAP drug–gene database from [Zenodo](https://zenodo.org/records/15692608) and place the CSV files in a local directory (e.g. `drug/`). You can also supply your own DEG files — one CSV per compound with at minimum `gene`, `log2FoldChange`, and `P.Value` columns.

---

## Running

```bash
python interface/server.py
```

Open **http://localhost:8000** in your browser.

Each module has a **Demo Mode** toggle. It fills in all paths and parameters using the built-in example dataset so you can verify the installation before loading your own data.

---

## Output

Each job creates an isolated directory under `interface/job_results/<job_id>/`:

```
interface/job_results/
└── a1b2c3d4/
    ├── config.json
    ├── *.png
    ├── *.csv
    └── <model files>
```

Results are served at `/api/results/<job_id>/<filename>` and rendered in the browser Results tab.

---

## Tutorials

| Notebook | Description |
|----------|-------------|
| [`dataset_preparation.ipynb`](tutorials/dataset_preparation.ipynb) | Preparing `.h5ad` files per disease stage |
| [`run_UNAGI_using_example_dataset.ipynb`](tutorials/run_UNAGI_using_example_dataset.ipynb) | Training UNAGI on the COPD airway dataset |
| [`visualize_UNAGI_results_example_dataset.ipynb`](tutorials/visualize_UNAGI_results_example_dataset.ipynb) | Inspecting latent space and trajectory outputs |
| [`in_silico_drug_discovery.ipynb`](tutorials/in_silico_drug_discovery.ipynb) | Full drug discovery walkthrough |
| [`Customize_drug_database_for_perturbation.ipynb`](tutorials/Customize_drug_database_for_perturbation.ipynb) | Using a custom compound DEG database |
| [`Customize_pathway_database_for_perturbation.ipynb`](tutorials/Customize_pathway_database_for_perturbation.ipynb) | Pathway-level perturbation scoring |

---

## Requirements

| Component | Version |
|-----------|---------|
| Python | ≥ 3.9 |
| PyTorch | ≥ 2.0 |
| scanpy | ≥ 1.9.5 |
| anndata | == 0.8.0 |
| pyro-ppl | ≥ 1.8.6 |
| R | ≥ 4.3 |
| Java | ≥ 1.7 (64-bit) |

---

## Citation

If you use scBioKit in published research, please cite:

> Zheng, Y., Schupp, J.C., Adams, T. et al. A deep generative model for deciphering cellular dynamics and in silico drug discovery in complex diseases. *Nature Biomedical Engineering* (2025). https://doi.org/10.1038/s41551-025-01423-7

---

## Contact

mabysxu@163.com
