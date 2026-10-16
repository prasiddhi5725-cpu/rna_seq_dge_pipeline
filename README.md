# Custom End-to-End RNA-Seq Differential Gene Expression (DGE) Pipeline

> A lightweight, dataset-agnostic Python pipeline for processing raw RNA-Seq count data, running differential expression analysis, and generating publication-grade volcano plots. Available as scripts, an upload-based CLI, and an interactive web GUI.

Built from scratch without black-box packages (DESeq2, edgeR) so normalization, testing, and multiple-testing correction stay fully transparent.

---

## Features

- CPM normalization, log2 fold change, Welch's t-test, Benjamini-Hochberg FDR
- Works with any count matrix in `.csv` or `.tsv` format
- Three ways to run it: fixed scripts, upload-based CLI, interactive GUI
- Publication-ready volcano plots with non-overlapping gene labels
- Downloadable results table and plot

---

## Methodology

### 1. Pre-processing and normalization

**Low-count filtering:** genes with fewer than 10 total reads across all samples are removed.

**Counts Per Million (CPM):** corrects for sequencing depth.

$$\text{CPM}_{i,j} = \left( \frac{c_{i,j}}{\sum_{k} c_{k,j}} \right) \times 10^6$$

**Log2 Fold Change:** a pseudocount ($\epsilon = 1.0$) smooths zero counts.

$$\text{Log}_2\text{FC}_i = \log_2(\bar{x}_{\text{treated},i} + \epsilon) - \log_2(\bar{x}_{\text{control},i} + \epsilon)$$

### 2. Hypothesis testing

- **Welch's t-test** per gene, with no equal-variance assumption.
- **Benjamini-Hochberg FDR** adjustment across all tests:

$$P_{(i)}^{\text{adj}} = \min \left( P_{(i)} \cdot \frac{m}{i}, P_{(i+1)}^{\text{adj}} \right)$$

### 3. Significance criteria

| Status | Condition |
|---|---|
| Upregulated | Log2FC ≥ 1.0 and adjusted p ≤ 0.10 |
| Downregulated | Log2FC ≤ -1.0 and adjusted p ≤ 0.10 |
| No significant change | Fails either threshold |

---

## Three Ways to Run

| Mode | File | Use when |
|---|---|---|
| Fixed pipeline | `gene_analyze.py` + `plot_results.py` | Reproducing the demo on `gene_counts.csv` |
| Upload-based CLI | `dge_for_uploaded_files.py` | Analyzing your own `.csv` / `.tsv` file |
| Interactive GUI | `dge_gui.py` | No-code analysis in the browser |

### 1. Fixed pipeline (demo data)

```bash
python3 generate_data.py     # optional: simulate a count matrix
python3 gene_analyze.py      # writes dge_results.csv
python3 plot_results.py      # writes volcano_plot.png
```

### 2. Upload-based CLI: `dge_for_uploaded_files.py`

Runs DGE on any count matrix you point it to.

```bash
python3 dge_for_uploaded_files.py
# Enter file name in csv or tsv format: my_counts.tsv
```

**Input format**

- First column: gene IDs (used as the index)
- Remaining columns: raw counts, one column per sample
- Delimiter is detected from the extension (`.csv` is comma-separated, anything else is tab-separated)

**Group assignment**

- If column names contain `Control 1` and `Treatment`, columns are grouped by those names (e.g. `Control 1`, `Control 2`, `Treatment 1`, `Treatment 2`).
- Otherwise the columns are split in half: first half is control, second half is treatment.

**Output:** a results table (`dge_results_2.csv`) with Log2FC, p-values, adjusted p-values, and significance status. Plot it with `plot_2.py`.

### 3. Interactive GUI: `dge_gui.py`

A Streamlit app for point-and-click analysis.

```bash
streamlit run dge_gui.py
```

Then open the local URL shown in the terminal (usually `http://localhost:8501`).

**What it does**

- Upload a `.csv` / `.tsv` count matrix
- Runs the full DGE analysis
- Shows an interactive volcano plot with the top genes labeled
- Shows a **Top 20 genes** table sorted by adjusted p-value
- Download buttons for the results (`dge_results.csv`) and the plot (`volcano_plot.jpg`, 300 dpi)

---

## Volcano Plot

![Volcano Plot](volcano_plot.png)

---

## Repository Structure

```text
.
├── generate_data.py            # Simulated count matrix generator
├── gene_analyze.py             # Core DGE pipeline (fixed input)
├── plot_results.py             # Volcano plot generator
├── dge_for_uploaded_files.py   # DGE on any user-supplied csv/tsv
├── dge_gui.py                  # Streamlit GUI
├── gene_counts.csv             # Demo count matrix
├── dge_results.csv             # Demo results
├── volcano_plot.png            # Demo figure
├── requirements.txt
└── README.md
```

---

## Installation

Python 3.8+ required.

```bash
pip install -r requirements.txt
```

`requirements.txt`:

```text
pandas
numpy
scipy
statsmodels
matplotlib
seaborn
adjustText
streamlit
```

---

## Stack

- **Language:** Python 3
- **Computation:** `pandas`, `numpy`, `scipy`, `statsmodels`
- **Visualization:** `matplotlib`, `seaborn`, `adjustText`
- **GUI:** `streamlit`
- **Environment:** Linux / WSL2 / Windows

---

## Author

Prasiddhi, BSc Bioinformatics, Kathmandu University
