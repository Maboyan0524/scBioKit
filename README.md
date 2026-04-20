# UNAGI + scBioKit

**UNAGI** is a deep generative model for deciphering cellular dynamics and performing unsupervised in-silico drug perturbation in complex diseases.

**scBioKit** is a web-based GUI built on top of UNAGI, providing an integrated interface for the full single-cell analysis pipeline — from QC to trajectory inference, drug perturbation, virtual gene knockout, and compound screening.

> Zheng, Y., Schupp, J.C., Adams, T. et al. *A deep generative model for deciphering cellular dynamics and in silico drug discovery in complex diseases.* **Nature Biomedical Engineering** (2025). https://doi.org/10.1038/s41551-025-01423-7

---

## What UNAGI Does

- Models disease progression across cell states using a **VAE + GCN** architecture
- Reconstructs temporal gene regulatory networks via **iDREM**
- Identifies dynamic marker genes, transcription factors, and trajectory branches
- Performs **unsupervised** in-silico perturbation for drugs, gene sets, and pathways

---

## scBioKit Interface

scBioKit wraps the full pipeline into a browser-based tool with five modules:

| Step | Module | Runtime | Description |
|------|--------|---------|-------------|
| 1 | Single-cell QC | R / Seurat | Filter cells, normalize, cluster |
| 2 | UNAGI Train | Python | Trajectory inference via VAE + GCN + iDREM |
| 3 | Drug Perturbation | Python | Score drugs by tendency to reverse disease trajectory |
| 4 | Virtual Knockout | R / scTenifoldKnk | Gene network knockout simulation |
| 5 | UniMol Screening | Python | Molecular representation-based compound ranking |

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/<your-username>/UNAGI.git
cd UNAGI
```

### 2. Create a Python environment

```bash
conda create -n unagi python=3.10
conda activate unagi
pip install .
```

### 3. Install iDREM

```bash
git clone https://github.com/phoenixding/idrem.git
```

Requires Java 1.7+ (64-bit).

### 4. Install interface dependencies

```bash
pip install fastapi uvicorn websockets
```

For the R modules (QC and Virtual Knockout):

```r
install.packages(c("Seurat", "qs", "scTenifoldKnk"))
```

---

## Running scBioKit

```bash
python interface/server.py
```

Then open **http://localhost:8000** in your browser.

---

## Python API (without GUI)

```python
import UNAGI

# Train
analyst = UNAGI.UNAGI_analyst(
    data_dir='path/to/stages',
    task_name='my_task',
    stage_key='stage',
    label_key='cell_type',
    total_stages=3
)
analyst.train()

# Drug perturbation
analyst.drug_perturbation(drug_dir='path/to/drug_degs')
```

Full documentation: [unagi-docs.readthedocs.io](https://unagi-docs.readthedocs.io/en/latest/)

---

## Tutorials

- [Dataset preparation](tutorials/dataset_preparation.ipynb)
- [Training on example dataset](tutorials/run_UNAGI_using_example_dataset.ipynb)
- [Visualizing results](tutorials/visualize_UNAGI_results_example_dataset.ipynb)
- [In-silico drug discovery walkthrough](tutorials/in_silico_drug_discovery.ipynb)

---

## Requirements

- Python ≥ 3.9
- torch ≥ 2.0.0
- scanpy ≥ 1.9.5
- anndata == 0.8.0
- pyro-ppl ≥ 1.8.6
- matplotlib ≥ 3.7.1

**CMAP drug-gene database** (required for drug perturbation):  
Download from [Zenodo](https://zenodo.org/records/15692608) and place in the project root.

---

## Contact

- [Yumin Zheng](mailto:yumin.zheng@mail.mcgill.ca)
- [Naftali Kaminski](mailto:naftali.kaminski@yale.edu)
- [Jun Ding](mailto:jun.ding@mcgill.ca)
