# Early Motor Imagery Classification with SE-Transformer

This repository implements **early motor imagery (MI) EEG classification** using **functional connectivity graphs** and an **SE-enhanced Transformer**. The goal is to evaluate **early decoding (≤800 ms)** of motor intention using variable-length graph sequences derived from EEG trials.

The workflow consists of **two stages**:

1. Generating connectivity graphs from raw EEG using MATLAB
2. Training and evaluating the SE-Transformer model in Python

---

## Dataset

This project uses the **BCI Competition IV-2a dataset**.

The Python code does **not** operate directly on raw EEG. Instead, it expects **precomputed connectivity graphs** stored as `.mat` files.

---

## Step 1: Generate Connectivity Graphs (MATLAB)

Connectivity graphs must be generated **file by file** from the BCI IV-2a dataset using the provided MATLAB code.

### What the MATLAB code does

For each EEG file (one subject, one session, one task):

1. **Load EEG data**

   * EEG shape: `22 × 750 × n_trials`
   * Segment boundaries (`segpnts`) are provided per trial

2. **Segment the EEG trial-wise**

   * Each trial is split into multiple short segments using SIS segmentation points
   * Segment timing is converted from milliseconds to sample indices

3. **Apply spatial filtering (optional but recommended)**

   * Surface Laplacian / **CSD transformation** using `CSDtoolbox (GetGH)`
   * Uses fixed 22-channel 10–20 electrode locations (θ, φ)

4. **Preprocess each segment**

   * Channel-wise zero-mean
   * Low-rank **SVD reconstruction** (retain ≥95% variance)

5. **Compute functional connectivity**

   * Pearson correlation between channels
   * Produces a **22 × 22 connectivity matrix per segment**

6. **Store metadata**

   * `trial_indices`
   * `segment_indices`
   * `segment_start_ms`, `segment_end_ms`
   * Sample-level timing information

7. **Save output as a `.mat` file**

   ```matlab
   connectivity_graphs_Axx[T/E]_<task>.mat
   ```

Each saved `.mat` file contains:

* `connectivity_graphs` : cell array of 22×22 matrices
* `trial_indices`
* `segment_indices`
* `segment_start_ms`, `segment_end_ms`
* channel `labels`

⚠️ **Important**:
You must run this MATLAB code **separately for each EEG file** (each subject, session, and task).

---

## Step 2: Organize Generated Files

After generating connectivity graphs for **all subjects and sessions**:

1. Upload all `.mat` files to a single directory (e.g., Google Drive).
2. Ensure filenames follow this pattern:

   ```
   connectivity_graphs_A01T_left.mat
   connectivity_graphs_A01E_left.mat
   connectivity_graphs_A02T_right.mat
   ...
   ```
3. Update the `data_dir` path in the Python code to point to this folder.

---

## Step 3: Run the Python SE-Transformer Code

Once all connectivity graphs are available:

1. Open `SE_transformer_best.ipynb`
2. The code:

   * Loads connectivity graphs from **both T and E sessions**
   * Groups segments into **trial-level sequences**
   * Converts each 22×22 graph into a **231-D feature vector** (upper triangle)
   * Trains an **SE-enhanced Transformer** on full trials
   * Evaluates **early decoding** using only segments ending ≤800 ms
   * Performs **8-fold cross-validation**

---

## Key Features

* Variable-length sequence modeling with padding masks
* Squeeze-and-Excitation (SE) blocks for channel-wise recalibration
* Early decoding (≤800 ms) without retraining
* Trial-wise cross-validation
* Reproducible results with fixed random seeds

---

## Requirements

### MATLAB

* MATLAB R2020+
* **CSDtoolbox** (`GetGH.m`, `CSD.m`)

### Python

* Python ≥ 3.8
* PyTorch
* NumPy, SciPy
* scikit-learn

---

## Notes

* Training uses **full trials**, while validation uses **early segments only**
* Early decoding is controlled via `segment_end_ms ≤ 800`
* The code assumes **22 EEG channels** in standard 10–20 layout

---

## Citation

If you use this code or pipeline in your research, please cite the corresponding paper (Sub-Second Motor Imagery Classification via Functional Connectivity Graphs and SE-Transformer) or acknowledge this repository.

---

