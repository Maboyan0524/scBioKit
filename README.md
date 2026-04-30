<div align="center">
  <img src="docs/banner.svg" alt="scBioKit" width="100%"/>
</div>
<div align="center">
![Python](https://img.shields.io/badge/Python-3.9%2B-3776ab?logo=python&logoColor=white)
![R](https://img.shields.io/badge/R-4.3%2B-276DC3?logo=r&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.110%2B-009688?logo=fastapi&logoColor=white)
![License](https://img.shields.io/badge/License-MIT-green)
![Paper](https://img.shields.io/badge/Nature%20Biomed%20Eng-2025-blue)
</div>
---
scBioKit is a browser-based platform that chains five single-cell analysis modules into one configurable pipeline — from raw scRNA-seq QC through disease trajectory modeling to ranked drug candidates — without requiring any coding. It wraps Seurat, UNAGI, scTenifoldKnk, and UniMol behind a FastAPI server with real-time log streaming over WebSocket.
---
Pipeline
```
Raw scRNA-seq data  (.qs / .rds / .rda)
         │
         ▼
 ┌────────────────────┐
 │  1. scRNA-seq QC   │  Filtering · Normalization · UMAP · Clustering      R / Seurat
 └──────────┬─────────┘
            │  filtered Seurat object
            ▼
 ┌────────────────────┐
 │  2. Virtual Gene   │  Co-expression network · KO simulation · DEG         R / scTenifoldKnk
 │     Knockout       │
 └──────────┬─────────┘
            │  validated drug targets
            ▼
 ┌────────────────────┐
 │  3. Compound       │  UniMol Transformer · Activity classifier · Library   Python
 │     Screening      │  screening
 └──────────┬─────────┘
            │  hit compound list
            ▼
 ┌────────────────────┐
 │  4. Trajectory     │  VAE + GCN · iDREM gene regulatory networks           Python / UNAGI
 │     Inference      │
 └──────────┬─────────┘
            │  trained model + iDREM trajectories
            ▼
 ┌────────────────────┐
 │  5. Drug           │  Latent-space perturbation · Reversal scoring          Python
 │     Perturbation   │
 └────────────────────┘
            │
            ▼
   Ranked drug candidates
```
Each module can also be run standalone. Results (figures, CSVs, model files) accumulate in per-job directories and are accessible directly from the browser.
---
Modules
1 · Single-Cell QC
<img src="interface/logo.png" align="right" width="48" style="margin:4px"/>
Seurat-based preprocessing pipeline. Filters cells by gene count and mitochondrial fraction thresholds, log-normalizes, selects highly variable genes, runs PCA, builds a neighbor graph, performs Leiden clustering, and projects into UMAP.
I/O	Files
Input	Seurat object (`.qs` / `.rds` / `.rda`)
Output	QC violin plots, UMAP (clusters + group), marker heatmap, summary CSV
Key parameters: `min_features`, `max_features`, `max_mito_pct`, `resolution`, `n_hvg`, `pca_dims`.  
Doublet detection via scDblFinder is optional.
---
2 · Virtual Gene Knockout
Uses scTenifoldKnk to build a sparse gene co-expression tensor from subsampled cells, then simulates the removal of one or more target genes. Differential manifold components identify downstream responders. Volcano plots and network diagrams are generated automatically.
I/O	Files
Input	Seurat object, knockout gene list (comma-separated)
Output	`ko_comparison.png`, `volcano_<gene>.png` per target, `ko_results.csv`
Key parameters: `ko_genes`, `cell_type_col`, `n_hvg`, `n_cells`, `n_networks`.
---
3 · Compound Screening
Fine-tunes a UniMol pre-trained molecular Transformer (84 M or 164 M parameters) on a user-supplied activity dataset, then applies the classifier to a compound library. Valid SMILES are filtered via RDKit before prediction; salt forms are discarded.
Training — binary classification with 5-fold CV, early stopping, optional external test-set evaluation (ROC, confusion matrix, probability distributions).
Screening — scores the full library, applies a user-defined activity threshold, and outputs ranked hit compounds.
Mode	Input	Output
Train	SMILES + label CSV	`final_model/`, `cv_auc_bar.png`, `roc_curve.png`
Screen	SMILES library CSV + trained model dir	`screen_top_actives.csv`, `score_distribution.png`, `top_actives_bar.png`
Both modes include a demo toggle that runs on built-in example data with no file upload required.
---
4 · Trajectory Inference (UNAGI)
Trains the UNAGI model (VAE optionally augmented with a Graph Convolutional Network) on per-stage `.h5ad` files. The model learns a continuous latent representation of disease progression and passes cell embeddings into iDREM to reconstruct dynamic gene regulatory networks: trajectory branches, active transcription factors, and stage-resolved marker genes.
I/O	Files
Input	`0.h5ad`, `1.h5ad`, … (one per disease stage) + iDREM directory
Output	`<task>.pth`, `attribute.pkl`, `dataset.h5ad`, `umap_stage.png`, `umap_celltype.png`
Key parameters: `distribution` (`ziln` / `normal`), `max_iter`, `batch_size`, `n_neighbors`, `device`.
> The trained model directory is the direct input to the Drug Perturbation module.
---
5 · Drug Perturbation
Reads drug-specific differential expression signatures (one CSV per compound, columns: gene, log₂FC, p-value) and applies each signature to the UNAGI latent space. A reversal score is computed for each drug by comparing how far cell embeddings shift toward healthy-state trajectories as defined by iDREM tendency labels. An optional anti-disease filter retains only genes whose direction is consistent with reversing the disease trajectory.
I/O	Files
Input	Trained UNAGI model (`.pth`), `attribute.pkl`, drug DEG directory
Output	`drug_scores.csv` (ranked), `drug_perturbation_scores.png`
Key parameters: `log2fc`, `background_samples`, `tend_thresh`, `anti_copd_filter`.
---
Installation
Clone
```bash
git clone https://github.com/<your-username>/scBioKit.git
cd scBioKit
```
Python environment
```bash
conda create -n scbiokit python=3.10
conda activate scbiokit

# Core framework
pip install fastapi "uvicorn[standard]" python-multipart websockets

# Single-cell stack
pip install scanpy anndata torch pyro-ppl

# UNAGI
pip install scUNAGI

# UniMol (compound screening)
pip install unimol_tools

# RDKit (SMILES validation)
pip install rdkit
```
R packages
```r
# Single-cell QC
install.packages(c("Seurat", "qs", "scDblFinder"))

# Virtual Knockout
install.packages("devtools")
devtools::install_github("cailab-tamu/scTenifoldKnk")
```
iDREM
iDREM is required for trajectory inference. It runs as a Java subprocess.
```bash
# Java 1.7+ (64-bit) must be on PATH
git clone https://github.com/phoenixding/idrem.git
```
Pass the path to the cloned directory as `idrem_dir` when running the Trajectory Inference module.
Drug gene signatures
Download the preprocessed CMAP drug–gene database from Zenodo and place the CSV files in a directory (e.g. `drug/`). Alternatively, supply your own DEG files — one CSV per compound with at minimum `gene`, `log2FoldChange`, and `P.Value` columns.
---
Running
```bash
python interface/server.py
```
Open http://localhost:8000 in your browser.
Each module has a Demo Mode toggle. Clicking it fills in all paths and parameters using the built-in example dataset so you can verify the installation works before loading your own data.
---
Output structure
Each job creates an isolated directory under `interface/job_results/<job_id>/`:
```
job_results/
└── a1b2c3d4/
    ├── config.json          # parameters used
    ├── *.png                # all generated figures
    ├── *.csv                # tabular results
    └── <model_files>        # .pth / model dirs (training jobs)
```
Results are served at `/api/results/<job_id>/<filename>` and rendered in the browser Results tab.
---
Tutorials
Notebook	Description
`dataset_preparation.ipynb`	Preparing `.h5ad` files per disease stage
`run_UNAGI_using_example_dataset.ipynb`	Training UNAGI on the COPD airway dataset
`visualize_UNAGI_results_example_dataset.ipynb`	Inspecting latent space and trajectory outputs
`in_silico_drug_discovery.ipynb`	Full drug discovery walkthrough
`Customize_drug_database_for_perturbation.ipynb`	Using a custom compound DEG database
`Customize_pathway_database_for_perturbation.ipynb`	Pathway-level perturbation scoring
---
Requirements summary
Component	Version
Python	≥ 3.9
PyTorch	≥ 2.0
scanpy	≥ 1.9.5
anndata	== 0.8.0
pyro-ppl	≥ 1.8.6
R	≥ 4.3
Java	≥ 1.7 (64-bit)
---
Citation
If you use scBioKit in published research, please cite the underlying UNAGI method:
> Zheng, Y., Schupp, J.C., Adams, T. et al. A deep generative model for deciphering cellular dynamics and in silico drug discovery in complex diseases. *Nature Biomedical Engineering* (2025). https://doi.org/10.1038/s41551-025-01423-7
---
Contact
Name	Role	Email
Yumin Zheng	Lead developer	yumin.zheng@mail.mcgill.ca
Naftali Kaminski	PI	naftali.kaminski@yale.edu
Jun Ding	PI	jun.ding@mcgill.ca
