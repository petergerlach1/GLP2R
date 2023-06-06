import pandas as pd
from sklearn.preprocessing import MinMaxScaler
import numpy as np

# Read UKBB variant table
UKBB = pd.read_excel("Data/GLP2R_variant_table.xlsx")
print(UKBB.columns)

# Read gnomad variant table
gnomad = pd.read_excel("Data/GLP2R_gnomad_variants.xlsx")
print(gnomad.columns)


df = pd.merge(UKBB, gnomad, on=["Variant", "Segment", "RaSP_pred"], how= "outer")

df = df[["Variant", "Segment", "RaSP_pred", "Allele Frequency", "AF", "Allele Count", "AC"]]

df = df.sort_values(by="AF", ascending=False)
df['number'] = range(1, len(df) + 1)


df["Variant"].to_csv("Data/GEMME_input.txt", index=False, header=False)

# import GEMME input
gemme = pd.read_csv("Data/GEMME_pred.txt", sep= " ")


df = pd.merge(df, gemme, on="Variant", how="left")


d = {'TM1': 'TM', 'TM2': 'TM', 'TM3': 'TM', 'TM4': 'TM', 'TM5': 'TM', 'TM6': 'TM', 'TM7': 'TM',
     'ICL1': 'ICL', 'ICL2': 'ICL', 'ICL3': 'ICL',
     'ECL1': 'ECL', 'ECL2': 'ECL', 'ECL3': 'ECL',}
df = df.replace({"Segment": d})
print(df)


# copy the data
df_sklearn = df.copy()
# apply normalization techniques
column = 'RaSP_pred'
column2 = 'GEMME_pred'
df_sklearn[column] = MinMaxScaler().fit_transform(np.array(df_sklearn[column]).reshape(-1, 1))
df_sklearn[column2] = MinMaxScaler().fit_transform(np.array(df_sklearn[column2]).reshape(-1, 1))


# view normalized data

print(df_sklearn.sort_values(by="RaSP_pred"))

df.to_excel("Data/all_variants.xlsx")



