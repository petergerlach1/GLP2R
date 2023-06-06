import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



#%%
camp = pd.read_excel("Data/In_vitro_sumfile.xlsx")
camp = camp.set_index('mutant')


#%% converting "N/A" values to 0
camp = camp.replace('non', np.nan)


#%%
cols0 = list(camp.columns)
cols = cols0[0:4]

camp = camp.dropna()

#for col in cols:
    #col_zscore = col + '_zscore'
    #camp[col_zscore] = (camp[col] - camp[col].mean())/camp[col].std(ddof=0)

camp["cAMP_Emax_zscore"] = (camp["cAMP_Emax"] - 100)/camp["cAMP_Emax"].std(ddof=0)
camp["cAMP_LogEC50_zscore"] = (camp["cAMP_LogEC50"] - (-10.7))/camp["cAMP_LogEC50"].std(ddof=0)


camp["arrestin_Emax_zscore"] = (camp["arrestin_Emax"] - 100)/camp["arrestin_Emax"].std(ddof=0)
camp["arrestin_LogEC50_zscore"] = (camp["arrestin_LogEC50"] - (-8.9))/camp["arrestin_LogEC50"].std(ddof=0)



#%%
camp = camp.drop(["cAMP_Emax","cAMP_LogEC50","arrestin_LogEC50","arrestin_Emax", "n_cAMP","n_arrestin"], axis=1)
print(camp.columns)
camp.to_csv("Data/variant_zscore.csv")
