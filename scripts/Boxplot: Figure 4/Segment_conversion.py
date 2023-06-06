import pandas as pd


df = pd.read_excel("Data/all_variants.xlsx")


d = {'TM1': 'TM', 'TM2': 'TM', 'TM3': 'TM', 'TM4': 'TM', 'TM5': 'TM', 'TM6': 'TM', 'TM7': 'TM',
     'ICL1': 'ICL', 'ICL2': 'ICL', 'ICL3': 'ICL',
     'ECL1': 'ECL', 'ECL2': 'ECL', 'ECL3': 'ECL'}
df = df.replace({"Segment": d})

df.to_excel("Data/all_variants_segmentconverted.xlsx")