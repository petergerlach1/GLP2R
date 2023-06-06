import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from results import pheno_search, pval_stars, plot_BT, plot_QT


# Field codings
ukb_coding = pd.read_csv(
    "Data/Data_Dictionary_Showcase.csv",
    error_bad_lines=False,
    warn_bad_lines=False,
    quotechar='"',
    usecols=["FieldID", "Field"],
)

custom_coding = pd.read_csv(
    "Data/Data_Dictionary_Custom.csv",
    error_bad_lines=False,
    warn_bad_lines=False,
    quotechar='"',
)



GENE = "GLP2R"
TRAIT = "BT"
PHENOTYPE = "metabolic"
AAF = None

#from results import pheno_search, pval_stars, plot_BT, plot_QT


list_of_files = ["Data/metabolic.QT.step2_QT_21001.csv",
                 "Data/metabolic.QT.step2_QT_78.csv",
                 "Data/metabolic.QT.step2_QT_48.csv",
                 "Data/metabolic.QT.step2_QT_49.csv",
                 "Data/metabolic.QT.step2_QT_23099.csv",
                 "Data/metabolic.QT.step2_QT_30750.csv",
                 "Data/metabolic.QT.step2_QT_4079.csv",
                 "Data/metabolic.QT.step2_QT_4080.csv",
                 "Data/metabolic.BT.step2_BT_131690.csv",
                 "Data/metabolic.BT.step2_BT_130706.csv",
                 "Data/metabolic.BT.step2_BT_130708.csv",
                 #"Data/metabolic.BT.step2_BT_1687-1.csv",
                 #"Data/metabolic.BT.step2_BT_1687-2.csv",
                 "Data/metabolic.BT.step2_BT_1687-3.csv"]

# %%
# Load raw DF
df_raw = pd.read_csv(list_of_files[0], delimiter=" ", header="infer", comment="#").assign(
    SOURCE=os.path.basename(list_of_files[0])
)
df_raw = pd.concat(
    [df_raw]
    + [
        pd.read_csv(fp, delimiter=" ", comment="#").assign(SOURCE=os.path.basename(fp))
        for fp in list_of_files[1:]
    ],
    axis=0,
)


# only using burden test
df_raw = df_raw.loc[df_raw["TEST"] == "ADD"]

# making the pvalues to real numbers
df_raw.loc[:, "pval"] = np.power(10, -df_raw["LOG10P"])


df = df_raw
df.loc[:, "GENE"] = df.ID.apply(lambda x: x.split(".")[0])
df.loc[:, "MASK"] = df.ALLELE1.apply(lambda x: x.split(".", maxsplit=2)[0])
df.loc[:, "AAF"] = df.ALLELE1.apply(lambda x: x.split(".", maxsplit=1)[-1])
df.loc[:, "TRAIT"] = df.SOURCE.apply(lambda x: x.split(".")[1])
df.loc[:, "PHENO"] = df.SOURCE.apply(lambda x: x.split("_")[-1].split(".")[0])
df = df.drop(["ID", "ALLELE0", "ALLELE1", "EXTRA", "SOURCE", "TEST"], axis=1)


# Filter AAF
if not AAF:
    AAF = df.AAF.max()

df = df.loc[df.AAF == AAF, :]

# Sanity check
df.head()

# %%
# Filters
bt = df.TRAIT == "BT"
qt = df.TRAIT == "QT"

# Fix Binary Traits
df.loc[bt, "OR"] = np.exp(df.loc[bt, "BETA"])
df.loc[bt, "OR_up"] = np.exp(df.loc[bt, "BETA"] + df.loc[bt, "SE"])
df.loc[bt, "OR_low"] = np.exp(df.loc[bt, "BETA"] - df.loc[bt, "SE"])
df.loc[bt, "OR_up_lim"] = df.loc[bt, "OR_up"] - df.loc[bt, "OR"]
df.loc[bt, "OR_low_lim"] = df.loc[bt, "OR"] - df.loc[bt, "OR_low"]

# Fix Quantitative Traits
df.loc[qt, "BETA_up_lim"] = df.loc[qt, "BETA"] + df.loc[qt, "SE"]
df.loc[qt, "BETA_low_lim"] = df.loc[qt, "BETA"] - df.loc[qt, "SE"]

# Final fixes
# Final fixes
#df.loc[:, "Phenotype"] = df.PHENO.apply(
    #lambda x: pheno_search(x, ukb_coding, custom_coding).replace('"', "").strip()
#)




df[["PHENO", "Phenotype_detail"]] = df.PHENO.str.split("-", expand = True)



# Final fixes
df.loc[:, "Phenotype"] = df.PHENO.apply(
    lambda x: pheno_search(x, ukb_coding, custom_coding).replace('"', "").strip()
)


df.loc[:, "pval"] = np.power(10, -df["LOG10P"])
df.loc[:, "pval_stars"] = df["pval"].apply(lambda x: pval_stars(x))
df.loc[:, "N_pos"] = (2 * df["N"] * df["A1FREQ"]).astype(int)

# Singletons
df = df.loc[df.AAF != "singleton", :]
print(len(df))
df.head()

# %%
phenos_to_remove = []
details_to_remove = [1, 2]


plt_df = (
    df.loc[(df.GENE == GENE)]
    .sort_values(by=["Phenotype", "AAF"], ascending=[True, False])  # , "AAF"
    .groupby(["Phenotype", "MASK"])
    .first()
    .reset_index()
)

plt_df = plt_df.loc[~plt_df.Phenotype.astype(str).isin(phenos_to_remove), :]


df.to_excel("test.xlsx")

effect = {"BT": "OR", "QT": "BETA"}[TRAIT]

group_by_mean = (
    pd.DataFrame({"mean": plt_df.groupby(["Phenotype"]).agg("mean")[effect]})
    .sort_values(by="mean", ascending=False)
    .reset_index()
)

sorter = group_by_mean.Phenotype.tolist()

plt_df.loc[:, "Phenotype"] = plt_df.loc[:, "Phenotype"].astype("category")
plt_df.loc[:, "Phenotypes"] = plt_df.loc[:, "Phenotype"].cat.set_categories(sorter)

plt_df = plt_df.sort_values(
    by=["Phenotype", "MASK"], ascending=[True, False]
).reset_index(drop=True)

phenotypes = plt_df.Phenotype.unique()

print(len(plt_df))
plt_df.head()


# %%
plot = plot_BT(
    plt_df, title=f"{GENE} Binary {PHENOTYPE.capitalize()} Traits", xlim=[0, 10]
)

plt.savefig(
    "Figures/BT_test.svg",
    dpi=600,
    bbox_inches="tight",
    format="svg",
)

plot = plot_QT(
    plt_df,
    title=f"{GENE} Quantitative {PHENOTYPE.capitalize()} Traits",
    xlim=[-1, 1],
    height=12,
)


plt.savefig(
    "Figures/QT_test.svg",
    dpi=600,
    bbox_inches="tight",
    format="svg",
)