# scBioKit

**scBioKit** is an integrated web-based platform for single-cell RNA-seq analysis, disease trajectory inference, and computational drug discovery. It connects five analysis modules into a unified configurable pipeline, executable through a browser interface without writing code.

---

## Pipeline

```
Raw scRNA-seq data
        │
        ▼
┌───────────────────┐
│  1. Single-cell   │  Cell filtering · Normalization · Clustering     R / Seurat
│     QC            │
└─────────┬─────────┘
          │
          ▼
┌───────────────────┐
│  2. UNAGI Train   │  VAE + GCN · iDREM trajectory inference          Python
└─────────┬─────────┘
          │
          ▼
┌───────────────────┐
│  3. Drug          │  Latent space perturbation · Tendency scoring     Python
│  Perturbation     │
└─────────┬─────────┘
          │
          ▼
┌───────────────────┐
│  4. Virtual       │  Gene regulatory network · KO simulation          R / scTenifoldKnk
│  Knockout         │
└─────────┬─────────┘
          │
          ▼
┌───────────────────┐
│  5. UniMol        │  Molecular representation · Compound ranking      Python
│  Screening        │
└───────────────────┘
          │
          ▼
  Ranked drug candidates
```

---

## Modules

### 1. Single-cell QC
Seurat-based quality control pipeline. Filters low-quality cells by gene count and mitochondrial fraction, normalizes expression, selects highly variable genes, reduces dimensions, and performs graph-based clustering.

- **Input:** Seurat object (`.qs` / `.rds` / `.rda`)
- **Output:** Filtered object, QC violin plots, UMAP clusters

### 2. UNAGI Train
Trajectory inference module. Trains a Variational Autoencoder with Graph Convolutional Network to learn a latent representation of disease progression across stages. Feeds cell embeddings into iDREM to reconstruct dynamic gene regulatory networks and identify trajectory branches, transcription factors, and stage-specific marker genes.

- **Input:** Per-stage `.h5ad` files, iDREM installation
- **Output:** Trained model (`.pth`), UMAP plots, iDREM trajectory, gene scores

### 3. Drug Perturbation
Applies drug-specific differential expression signatures to the UNAGI latent space and scores each drug by how strongly it shifts cells toward healthier states along the inferred trajectory. Supports single drugs, combinations, pathways, and custom compound databases.

- **Input:** Trained UNAGI model, drug DEG directory (one CSV per drug)
- **Output:** Per-drug tendency scores, perturbation UMAPs, ranked drug list

### 4. Virtual Knockout
Uses scTenifoldKnk to construct a gene co-expression network from single-cell data, simulate removal of target genes, and identify downstream transcriptional changes. Useful for validating drug targets identified in the perturbation step.

- **Input:** Seurat object, target gene list
- **Output:** Network visualization, KO effect plots, DEG results

### 5. UniMol Screening
Trains a binary activity classifier using UniMol molecular representations (84M or 164M pre-trained transformer). Screens a compound library and ranks hits by predicted activity probability.

- **Input:** Training CSV (SMILES + label column), compound library CSV
- **Output:** ROC curve, ranked screening results

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/<your-username>/scBioKit.git
cd scBioKit
```

### 2. Python environment

```bash
conda create -n scbiokit python=3.10
conda activate scbiokit
pip install .
pip install fastapi uvicorn websockets
```

### 3. R packages

```r
# For Single-cell QC
install.packages(c("Seurat", "qs"))

# For Virtual Knockout
devtools::install_github("cailab-tamu/scTenifoldKnk")
```

### 4. iDREM (for UNAGI trajectory)

```bash
git clone https://github.com/phoenixding/idrem.git
```

Requires Java 1.7+ (64-bit).

### 5. Drug database (for Drug Perturbation)

Download the preprocessed CMAP database from [Zenodo](https://zenodo.org/records/15692608) and place the files in the project root. Alternatively, provide your own drug–gene DEG files.

---

## Running scBioKit

Start the backend server:

```bash
python interface/server.py
```

Open **http://localhost:8000** in your browser.

---

## Tutorials

- [Dataset preparation](tutorials/dataset_preparation.ipynb)
- [UNAGI training on example dataset](tutorials/run_UNAGI_using_example_dataset.ipynb)
- [Visualizing results](tutorials/visualize_UNAGI_results_example_dataset.ipynb)
- [In-silico drug discovery walkthrough](tutorials/in_silico_drug_discovery.ipynb)
- [Custom drug database](tutorials/Customize_drug_database_for_perturbation.ipynb)
- [Custom pathway database](tutorials/Customize_pathway_database_for_perturbation.ipynb)

---

## Requirements

| Component | Requirement |
|-----------|-------------|
| Python | ≥ 3.9 |
| PyTorch | ≥ 2.0.0 |
| scanpy | ≥ 1.9.5 |
| anndata | == 0.8.0 |
| pyro-ppl | ≥ 1.8.6 |
| R | ≥ 4.3 |
| Java | ≥ 1.7 (64-bit) |

---

## Citation

If you use scBioKit in your research, please cite:

Zheng, Y., Schupp, J.C., Adams, T. et al. A deep generative model for deciphering cellular dynamics and in silico drug discovery in complex diseases. *Nature Biomedical Engineering* (2025). https://doi.org/10.1038/s41551-025-01423-7

---

## Contact

- [Yumin Zheng](mailto:yumin.zheng@mail.mcgill.ca)
- [Naftali Kaminski](mailto:naftali.kaminski@yale.edu)
- [Jun Ding](mailto:jun.ding@mcgill.ca)
