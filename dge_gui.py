import io
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests

st.set_page_config(page_title="Differential Gene Expression", page_icon="🧬", layout="wide")

st.markdown("""
<style>
    .banner { background: linear-gradient(135deg, #1E3A8A, #3B82F6); padding: 2rem; border-radius: 14px; margin-bottom: 1.5rem; }
    .banner h1 { color: white; margin: 0; padding: 0; }
    .banner p { color: #DBEAFE; margin: 0.4rem 0 0 0; }
</style>
<div class="banner">
    <h1>🧬 Differential Gene Expression</h1>
    <p>Upload RNA-Seq counts, run the analysis, get a volcano plot and a results table.</p>
</div>
""", unsafe_allow_html=True)

st.sidebar.header("⚙️ Settings")
min_sum = st.sidebar.slider("Minimum total reads per gene", 0, 50, 10)
fc_cut = st.sidebar.slider("Log2 fold change cutoff", 0.0, 4.0, 1.0, 0.1)
padj_cut = st.sidebar.slider("Adjusted p-value cutoff", 0.001, 0.200, 0.05, 0.001)

uploaded_file = st.file_uploader("Upload a count file (CSV/TSV)", type=["csv", "tsv"])

if uploaded_file is None:
    st.info("👆 Upload a file to begin.")
    st.stop()

sep = "," if uploaded_file.name.endswith(".csv") else "\t"
df = pd.read_csv(uploaded_file, sep=sep, index_col=0)
df = df.select_dtypes(include="number")

all_cols = list(df.columns)
if len(all_cols) < 2:
    st.error("Need at least 2 sample columns.")
    st.stop()

st.subheader("Data preview")
st.dataframe(df.head(8), use_container_width=True)

if "Control 1" in all_cols:
    default_ctrl = [c for c in all_cols if "Control" in c]
    default_treat = [c for c in all_cols if "Treatment" in c]
else:
    mid = len(all_cols) // 2
    default_ctrl = all_cols[:mid]
    default_treat = all_cols[mid:]

col1, col2 = st.columns(2)
control_cols = col1.multiselect("Control samples", all_cols, default=default_ctrl)
treat_cols = col2.multiselect("Treatment samples", all_cols, default=default_treat)

if st.button("🚀 Run analysis", type="primary"):

    if len(control_cols) == 0 or len(treat_cols) == 0:
        st.error("Pick at least one control and one treatment sample.")
        st.stop()

    lib_size = df.sum(axis=0)
    cpm = df.div(lib_size, axis=1) * 1_000_000
    cpm = cpm[df.sum(axis=1) >= min_sum]
    log_cpm = np.log2(cpm + 1)

    if len(cpm) == 0:
        st.error("No genes passed the filter. Lower the minimum reads.")
        st.stop()

    res = pd.DataFrame(index=cpm.index)
    res["log2fc"] = log_cpm[treat_cols].mean(axis=1) - log_cpm[control_cols].mean(axis=1)

    if len(control_cols) >= 2 and len(treat_cols) >= 2:
        st.write("Mode: Welch's t-test")
        _, p_vals = stats.ttest_ind(
            log_cpm[control_cols].values,
            log_cpm[treat_cols].values,
            axis=1,
            equal_var=False,
        )
        p_vals = np.nan_to_num(p_vals, nan=1.0)
    else:
        st.write("Mode: single replicate (approximate, BCV = 0.2)")
        mean_c = cpm[control_cols].mean(axis=1).values
        mean_t = cpm[treat_cols].mean(axis=1).values
        se = np.sqrt((0.2 * mean_c) ** 2 / len(control_cols) + (0.2 * mean_t) ** 2 / len(treat_cols))
        se = np.where(se == 0, 1e-6, se)
        p_vals = 2 * stats.norm.sf(np.abs((mean_t - mean_c) / se))

    res["p_value"] = np.clip(p_vals, 0, 1)
    res["padj"] = multipletests(res["p_value"], method="fdr_bh")[1]

    up = (res["log2fc"] >= fc_cut) & (res["padj"] < padj_cut)
    down = (res["log2fc"] <= -fc_cut) & (res["padj"] < padj_cut)
    res["Regulation"] = np.select([up, down], ["Upregulated", "Downregulated"], default="No significant change")
    res["-log10_padj"] = -np.log10(np.clip(res["padj"], 1e-300, 1.0))

    counts = res["Regulation"].value_counts()
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Genes tested", len(res))
    m2.metric("Upregulated", counts.get("Upregulated", 0))
    m3.metric("Downregulated", counts.get("Downregulated", 0))
    m4.metric("Not significant", counts.get("No significant change", 0))

    colors = {"Upregulated": "#e63946", "Downregulated": "#457b9d", "No significant change": "#a8a29e"}

    sns.set_theme(style="ticks")
    fig, ax = plt.subplots(figsize=(9, 6), dpi=150)
    sns.scatterplot(data=res, x="log2fc", y="-log10_padj", hue="Regulation",
                    palette=colors, alpha=0.4, s=12, edgecolor=None, ax=ax)

    ax.axvline(fc_cut, color="black", linestyle="--", linewidth=0.8)
    ax.axvline(-fc_cut, color="black", linestyle="--", linewidth=0.8)
    ax.axhline(-np.log10(padj_cut), color="black", linestyle="--", linewidth=0.8)

    sig = res[res["Regulation"] != "No significant change"].copy()
    sig["abs_fc"] = sig["log2fc"].abs()
    top = sig.sort_values("abs_fc", ascending=False).head(5)
    for gene_id, row in top.iterrows():
        ax.text(row["log2fc"], row["-log10_padj"], str(gene_id), fontsize=8, weight="bold")

    ax.set_title("Volcano Plot", fontsize=14, weight="bold")
    ax.set_xlabel("Log2 Fold Change")
    ax.set_ylabel("-log10(adjusted p-value)")
    ax.legend(title="Status", bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()

    st.pyplot(fig)

    st.subheader("Top 20 genes")
    st.dataframe(res.sort_values("padj")[["log2fc", "p_value", "padj", "Regulation"]].head(20),
                 use_container_width=True)
    
    buf = io.BytesIO()
    fig.savefig(buf, format="jpeg", dpi=300, bbox_inches="tight", facecolor="white")

    d1, d2 = st.columns(2)
    d1.download_button("📥 Results (.csv)", res.to_csv().encode("utf-8"),
                       file_name="dge_results.csv", mime="text/csv", on_click="ignore")
    d2.download_button("📥 Volcano plot (.jpg)", buf.getvalue(),
                       file_name="volcano_plot.jpg", mime="image/jpeg", on_click="ignore")

    