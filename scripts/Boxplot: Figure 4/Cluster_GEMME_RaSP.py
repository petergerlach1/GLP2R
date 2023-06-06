import pandas as pd

df = pd.read_excel("Data/all_variants.xlsx")
df = df.rename(columns= {"Variant": "mutant"})

df2 = pd.read_csv("Data/regenie_clustering_groups.csv", delimiter=",")
df2 = df2.rename(columns={"GeneID": "mutant"})

# renaming clusters
cluster_mapping = {
    'cluster1': 'cAMP_WT_arrestin_decreased_activity',
    'cluster2': 'cAMP_decreased_activity_arrestin_WT',
    'cluster3': 'WT',
    'cluster4': 'WT',
    'cluster5': 'LoF',
    'cluster6': 'cAMP_increased_activity_arrestin_decreased_activity'
}

df2['Cluster'] = df2['Cluster'].replace(cluster_mapping)



df = pd.merge(df2, df, on = "mutant", how = "left")

df[["mutant", "RaSP_pred", "GEMME_pred", "Cluster"]].to_csv("Data/RaSP_GEMME_cluster.csv", index=False)