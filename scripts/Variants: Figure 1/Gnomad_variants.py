import pandas as pd
import requests

# load data and select relevant columns
df = pd.read_csv("Data/gnomad_variants.csv")

df = df[["Chromosome",
         "Position",
         "rsIDs",
         "Reference",
         "Alternate",
         "Allele Frequency",
         "Allele Count",
         "HGVS Consequence"]]

# making markerID same as for UKBB
df["ID"] = df['Chromosome'].astype(str) +":"+ df["Position"].astype(str) + ":" + \
           df["Reference"].astype(str) + ":" + df["Alternate"].astype(str)


# making the mutant column to be ref/protein_position/alt with 1-letter amino acid codes
df["ref"] = df["HGVS Consequence"].str[2:5]
df["alt"] = df["HGVS Consequence"].str[-3:]
df["Protein_position"] = df["HGVS Consequence"].str[5:-3]
d = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
     'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Asn': 'N',
     'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W',
     'Ala': 'A', 'Val':'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M'}
df = df.replace({"ref": d,"alt": d})
df["mutant"] = df['ref'].astype(str) + df['Protein_position'].astype(str) + df['alt'].astype(str)

# dropping unnecessary columns
df.drop(columns=["Chromosome", "Position", "Reference", "Alternate", "HGVS Consequence", "ref", "alt"], inplace=True)

# drop variant duplicates in the dataset
df = df.drop_duplicates(subset=['mutant'])



# Importing data from VEP
VEP_output = pd.read_csv("Data/VEP_output.txt", sep='\s+')


df["Protein_position"] = df["Protein_position"].astype(int)



### Adding in vitro data to the table
in_vitro_data = pd.read_excel("Data/In_vitro_sumfile.xlsx")
#in_vitro_data = in_vitro_data.dropna(axis = 0)
#in_vitro_data["Protein_position"] = in_vitro_data["Protein_position"].astype(int)



df = pd.merge(df, in_vitro_data,  how='left', on="mutant")



### Adding GPCRdbnumbering, functional annotation and protein segment to the dataframe

# For GPCRdbnumbering and segment
url1 = "https://gpcrdb.org/services/residues/extended/glp2r_human/"
response = requests.get(url1)
protein_data1 = response.json()
SG_df = pd.DataFrame(protein_data1)
SG_df = SG_df.drop('alternative_generic_numbers', axis=1)

# for functional annotation
url2 = "https://gpcrdb.org/services/structure/7D68/interaction/"
response2 = requests.get(url2)
protein_data2 = response2.json()
FA_df = pd.DataFrame(protein_data2)

FA_df = FA_df[["sequence_number", "interaction_type"]]
# Some positions have multiple functional annotations, so I move them into the same data field
FA_df = FA_df.groupby('sequence_number').agg({
                             'interaction_type': ', '.join}).reset_index()

# merging segment, GPCRdbnumering and functional annotation together
SG_FA_df = SG_df.merge(FA_df,how='left', left_on='sequence_number', right_on='sequence_number')
SG_FA_df = SG_FA_df.rename(columns={"sequence_number": "Protein_position", "protein_segment": "Segment",
                                    "display_generic_number": "GPCRdbnumbering", "interaction_type": "Functional Annotation"})


# Merging the Segment, GPCRdb, and functional annotation to the table
df = df.merge(SG_FA_df, on='Protein_position', how='left', indicator=True)






# adding RaSP
# Add predictions based on determined structure


# Add predictions based on alphafold
rasp_alpha = pd.read_csv("Data/RaSP_alphafold_pred.csv")

rasp_alpha = rasp_alpha[["variant", "score_ml_fermi", "score_ml"]]
rasp_alpha.rename(columns={"variant": "mutant", "score_ml": "RaSP_pred"}, inplace = True)

# Merge with RaSP predictor
df = df.merge(rasp_alpha, on='mutant', how='left')

df = df.rename(columns={"mutant": "Variant"})

df = df.set_index("ID")


df = df.sort_values(by="Allele Frequency", ascending=False)

### Making final table
df[["Variant", "Segment", "Allele Frequency", "Allele Count", "Functional Annotation",
    "cAMP_Emax", "cAMP_LogEC50", "arrestin_LogEC50", "arrestin_Emax", "RaSP_pred", "rsIDs"]].\
    to_excel("Data/GLP2R_gnomad_variants.xlsx")



