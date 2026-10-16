import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

file_name=input("Enter file name in csv or tsv format: ")

sep="," if file_name.endswith(".csv") else "\t"
df=pd.read_csv(file_name, sep=sep, index_col=0)

df_filtered=df[df.sum(axis=1) >=10].copy()

total_reads=df_filtered.sum(axis=0)
counts_per_million=df_filtered.div(total_reads, axis=1)*1_000_000

all_cols=list(counts_per_million.columns)
n_cols=len(all_cols)

if "Control 1" in all_cols:
    control_cols = [c for c in all_cols if "Control" in c]
    treatment_cols = [c for c in all_cols if "Treatment" in c]
else:
    mid = n_cols // 2
    control_cols = all_cols[:mid]
    treatment_cols = all_cols[mid:]

n_ctrl = len(control_cols)
n_treat = len(treatment_cols)

cpm_mean_control = counts_per_million[control_cols].mean(axis=1)
cpm_mean_treatment = counts_per_million[treatment_cols].mean(axis=1)

log2fc=np.log2((cpm_mean_treatment+1)/(cpm_mean_control+1))

df_filtered["log2fc"]=log2fc

if n_ctrl >= 2 and n_treat >= 2:
    print("Mode: Standard Welch's t-test (Biological Replicates Available)")
    _, p_vals = stats.ttest_ind(
        counts_per_million[control_cols].values,
        counts_per_million[treatment_cols].values,
        axis=1,
        equal_var=False,
    )
    df_filtered["p_value"] = np.nan_to_num(p_vals, nan=1.0)
else:
    print(
        "Mode: Single-replicate dispersion model (BCV = 0.2 estimated)"
    )
    bcv = 0.2  
    ctrl_std = cpm_mean_control * bcv
    treat_std = cpm_mean_treatment * bcv
    se_diff = np.sqrt((ctrl_std**2 / n_ctrl) + (treat_std**2 / n_treat))
    se_diff = np.where(se_diff == 0, 1e-6, se_diff)
    t_stat = (cpm_mean_treatment - cpm_mean_control) / se_diff

    log_p_vals = stats.norm.logsf(np.abs(t_stat)) + np.log(2)
    df_filtered["p_value"] = np.exp(log_p_vals)

_, padj, _, _ = multipletests(df_filtered["p_value"], method="fdr_bh")
df_filtered["padj"] = padj

conditions=[
    (df_filtered["log2fc"] >= 1) & (df_filtered["padj"]<0.10), 
    (df_filtered["log2fc"] <= -1) & (df_filtered["padj"]<0.10)
]

choices=["Upregulated", "Downregulated"]

df_filtered["Regulation"] = np.select(conditions, choices, default="No significant change")

print("---Differential Expression Summary---")
print(df_filtered["Regulation"].value_counts())
print(df_filtered.sort_values("padj")[["log2fc", "p_value", "padj", "Regulation"]].head(10).to_markdown())

df_filtered.to_csv("dge_results_2.csv")
print("Results successfully saved to dge_results_2.csv")