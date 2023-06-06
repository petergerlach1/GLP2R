import pandas as pd

# Getting clusters from K-means
df = pd.read_csv("Data/regenie_clustering_groups.csv", delimiter=",")
df = df.rename(columns={"GeneID": "mutant"})

# renaming clusters
cluster_mapping = {
    'cluster1': 'WT/part_LoF',
    'cluster2': 'part_LoF/WT',
    'cluster3': 'WT1',
    'cluster4': 'WT2',
    'cluster5': 'LoF',
    'cluster6': 'part_LoF/part_GoF'
}

df['Cluster'] = df['Cluster'].replace(cluster_mapping)
print(df)

# importing and merging annotation file with cluster file
anno_file = pd.read_csv("/Users/petergerlach/GLP2R_project/Regenie_output/Data/metabolic.txt", delimiter="\t")

GLP2R_annofile = anno_file[anno_file['gene'] == "GLP2R"]
rest_annofile = anno_file[anno_file['gene'] != "GLP2R"]

GLP2R_annofile = pd.merge(GLP2R_annofile, df, how= "left", on= "mutant")

anno_file = pd.concat([GLP2R_annofile, rest_annofile])

anno_file = anno_file.rename(columns={"Cluster": "annotation"})
anno_file["annotation"] = anno_file["annotation"].fillna("None")


anno_file[["ID", "gene", "annotation"]].to_csv("metabolic_final.annotations", sep="\t", index= False, header=False)