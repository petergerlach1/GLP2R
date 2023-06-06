import pandas as pd
import requests

### Importing
# Loading file from UK Biobank
ukb_df = pd.read_csv("Data/GLP2R_QC-alt-samples_swiss_army_knife.txt", sep=' ', names=["Chromosome",'Position', 'REF', "ALT", 'Sample',
                                                                 'Genotype', "AC", "AF", "ID"], index_col=None)

# Importing data from VEP
VEP_output = pd.read_csv("Data/VEP_output.txt", sep='\s+')


### Formatting
## UK Biobank file
# Removing variant duplicates from file
ukb_df = ukb_df.drop_duplicates(subset="ID", keep="first")


## VEP output data
# Selecting columns
VEP_output = VEP_output[["#Uploaded_variation", "Amino_acids", "Codons", "Protein_position", "Existing_variation",
         "VEST4_rankscore"]]

# changing data format
VEP_output = VEP_output.rename(columns={"#Uploaded_variation": "ID"})
VEP_output['ID'] = VEP_output['ID'].str.replace('_', ':')
VEP_output['ID'] = VEP_output['ID'].str.replace('/', ':')


### Merging UK Biobank dataframe with VEP output
# merging the two dataframes
df = VEP_output.merge(ukb_df, on='ID', how='left')

## fomatting the dataframe
df["Protein_position"] = df["Protein_position"].astype(int)


df["mutant"] = df['Amino_acids'].str[0] + df['Protein_position'].astype(str) + df['Amino_acids'].str[-1]


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
# Add predictions based on alphafold
rasp_alpha = pd.read_csv("Data/RaSP_alphafold_pred.csv")

rasp_alpha = rasp_alpha[["variant", "score_ml_fermi", "score_ml"]]
rasp_alpha.rename(columns={"variant": "mutant", "score_ml": "RaSP_pred"}, inplace = True)

# Merge with RaSP predictor
df = df.merge(rasp_alpha, on='mutant', how='left')

# sort by allele frequency
df = df.rename(columns={"Existing_variation": "rsID", "mutant": "Variant"})

df = df.set_index("ID")

df = df.sort_values(by="AF", ascending=False)

### Making final table
df[["Variant", "Segment", "AF", "AC", "Functional Annotation",
    "cAMP_Emax", "cAMP_LogEC50", "arrestin_LogEC50", "arrestin_Emax", "RaSP_pred"]].\
    to_excel("Data/GLP2R_variant_table.xlsx")


