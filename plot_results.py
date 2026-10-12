import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

df=pd.read_csv("dge_results.csv",index_col=0)

df["-log10_padj"]=-np.log10(df["padj"]+1e-15)

fig, ax = plt.subplots(figsize=(9, 6), dpi=300)
sns.set_theme(style="ticks")

colors={
    "Upregulated": "#e63946",
    "Downregulated": "#457b9d",
    "No significant change": "#a8a29e"
}

sns.scatterplot(
    data=df,
    x="log2fc",
    y="-log10_padj",
    hue="Regulation",
    palette=colors,
    alpha=0.8,
    s=35,
    edgecolor=None,
    ax=ax
)

ax.axvline(x=1.0, color="black", linestyle="--", linewidth=0.8, alpha=0.7)
ax.axvline(x=-1.0, color="black", linestyle="--", linewidth=0.8, alpha=0.7)
ax.axhline(y=-np.log10(0.10), color="black", linestyle="--", linewidth=0.8, alpha=0.7)
top_genes=df[df["Regulation"]!="No significant change"]

top_genes = df[df["Regulation"] != "No significant change"].copy()

try:
    from adjustText import adjust_text
    texts = []
    for id, row in top_genes.iterrows():
        texts.append(
            ax.text(
                row["log2fc"],
                row["-log10_padj"],
                str(id),
                fontsize=6,
                weight="semibold",
                color="#1d3557"
            )
        )
    adjust_text(texts, arrowprops=dict(arrowstyle="->", color="#333333", lw=0.5), expand_points=(1.8, 1.8))
except ImportError:
    stagger_step = 0
    for id, row in top_genes.iterrows():
        x_offset = 0.04
        y_offset = 0.03 + (stagger_step * 0.04)  
        ax.text(
            row["log2fc"] + x_offset,
            row["-log10_padj"] + y_offset,
            str(id),
            fontsize=9,
            weight="normal",
            color="#1d3557"
        )
        stagger_step += 1

y_max = df["-log10_padj"].max()
ax.set_ylim(-0.05, y_max + 0.25)

ax.set_title("Volcano Plot: Differential Gene Expression", fontsize=14, weight="bold", pad=15)
ax.set_xlabel("Log2 Fold Change (Log2FC)", fontsize=11)
ax.set_ylabel("-log10(Adjusted p-value)", fontsize=11)

ax.legend(
    title="Status",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    borderaxespad=0,
    frameon=True
)

plt.tight_layout()
plt.savefig("volcano_plot.png", bbox_inches="tight")
print("Clean, un-cluttered volcano plot saved to 'volcano_plot.png'!")