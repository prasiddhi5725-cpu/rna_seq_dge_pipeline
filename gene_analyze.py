import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

df=pd.read_csv("gene_counts.csv", index_col=0)

df_filtered=df[df.sum(axis=1) >=10].copy()

total_reads=df_filtered.sum(axis=0)
counts_per_million=df_filtered.div(total_reads, axis=1)*1_000_000

cpm_mean_control=counts_per_million[["Control 1", "Control 2", "Control 3"]].mean(axis=1)
cpm_mean_treatment=counts_per_million[["Treatment 1", "Treatment 2", "Treatment 3"]].mean(axis=1)

log2fc=np.log2((cpm_mean_treatment+1)/(cpm_mean_control+1))

df_filtered['log2fc']=log2fc

t_stat , p_values=stats.ttest_ind(
    counts_per_million[["Control 1", "Control 2", "Control 3"]].values,
    counts_per_million[["Treatment 1", "Treatment 2", "Treatment 3"]],
    axis=1,
    equal_var=False
)

df_filtered['p_value']=p_values
_, padj, _, _=multipletests(df_filtered['p_value'], method='fdr_bh')
df_filtered['padj']=padj

conditions=[
    (df_filtered['log2fc'] >= 1) & (df_filtered["padj"]<0.10), 
    (df_filtered['log2fc'] <= -1) & (df_filtered["padj"]<0.10)
]

choices=['Upregulated', 'Downregulated']

df_filtered['Regulation'] = np.select(conditions, choices, default='No significant change')

print("---Differential Expression Summary---")
print(df_filtered['Regulation'].value_counts())
print(df_filtered.sort_values('padj')[['log2fc', 'p_value', 'padj', 'Regulation']].head(10).to_markdown())

df_filtered.to_csv("dge_results.csv")
