import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

df = pd.read_csv("dge_results_2.csv", index_col=0)

df["-log10_pvalue"] = -np.log10(np.clip(df["p_value"], 1e-300, 1.0))

fig, ax = plt.subplots(figsize=(9, 6), dpi=300)
sns.set_theme(style="ticks")

colors = {
    "Upregulated": "#e63946",
    "Downregulated": "#457b9d",
    "No significant change": "#a8a29e",
}

sns.scatterplot(
    data=df,
    x="log2fc",
    y="-log10_pvalue",
    hue="Regulation",
    palette=colors,
    alpha=0.35,
    s=12,
    edgecolor=None,
    ax=ax,
)

ax.axvline(x=1.0, color="black", linestyle="--", linewidth=0.8, alpha=0.7)
ax.axvline(x=-1.0, color="black", linestyle="--", linewidth=0.8, alpha=0.7)
ax.axhline(y=-np.log10(0.10), color="black", linestyle="--", linewidth=0.8, alpha=0.7)

sig_genes = df[df["Regulation"] != "No significant change"].copy()
sig_genes["abs_fc"] = sig_genes["log2fc"].abs()
top_genes = sig_genes.sort_values("abs_fc", ascending=False).head(5)

try:
    from adjustText import adjust_text

    texts = [
        ax.text(
            row["log2fc"],
            row["-log10_pvalue"],
            str(gene_id),
            fontsize=8,
            weight="bold",
            color="#1d3557",
        )
        for gene_id, row in top_genes.iterrows()
    ]
    adjust_text(
        texts,
        arrowprops=dict(arrowstyle="->", color="black", lw=0.6),
        expand_points=(2.0, 2.0),
    )
except ImportError:
    for gene_id, row in top_genes.iterrows():
        ax.text(
            row["log2fc"] + 0.1,
            row["-log10_padj"] + 0.1,
            str(gene_id),
            fontsize=8,
            weight="bold",
            color="#1d3557",
        )

ax.set_title(
    "Volcano Plot: Differential Gene Expression",
    fontsize=14,
    weight="bold",
    pad=15,
)
ax.set_xlabel("Log2 Fold Change (Log2FC)", fontsize=11)
ax.set_ylabel("-log10(p-value)", fontsize=11)

ax.legend(
    title="Status",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    borderaxespad=0,
    frameon=True,
)

plt.tight_layout()
plt.savefig("volcano_plot_final.png", bbox_inches="tight")
print("Saved uncapped plot to 'volcano_plot_final.png'!")