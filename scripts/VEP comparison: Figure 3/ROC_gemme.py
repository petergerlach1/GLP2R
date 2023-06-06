import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import label_binarize


df = pd.read_excel("Data/RaSP_GEMME.xlsx")
df = df.rename(columns={"Variant":"mutant"})

classification = pd.read_csv("Data/WT_LoF_in_vitro.csv")

df = pd.merge(df, classification, on= "mutant", how= "inner")

df = df[["mutant", "Cluster","GEMME_pred"]]





x = df["GEMME_pred"].values
y = df['Cluster'].values

print(x)

y = label_binarize(y, classes=["LoF", "WT"])
n_classes = y.shape[1]


fpr = dict()
tpr = dict()
roc_auc = dict()
for i in range(n_classes):
    fpr[i], tpr[i], _ = roc_curve(y, x)
    roc_auc[i] = auc(fpr[i], tpr[i])


#plt.figure()
lw = 2
plt.plot(
    fpr[0],
    tpr[0],
    color="darkorange",
    lw=lw,
    label="ROC curve (area = %0.2f)" % roc_auc[0],
)
plt.plot([0, 1], [0, 1], color="navy", lw=lw, linestyle="--")
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Receiver operating characteristic example")
plt.legend(loc="lower right")
plt.show()

