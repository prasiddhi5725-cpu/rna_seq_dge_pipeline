import pandas as pd
import numpy as np

generate_data=np.random.negative_binomial(n=10, p=0.01, size=(1000, 6))
row=[f"Gene {i}" for i in range(1, 1001)]
columns=[f"Control {j}" for j in range(1, 4)] + [f"Treatment {k}" for k in range(1,4)]

df=pd.DataFrame(data=generate_data, index=row, columns=columns)

target_rows=row[0:50]
target_columns=columns[3:6]


df.loc[target_rows, target_columns] = df.loc[target_rows, target_columns]*4

df.to_csv("gene_counts.csv")
