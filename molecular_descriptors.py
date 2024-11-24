
# Install Rdkit
!pip install rdkit

import pandas as pd
import numpy as np
import seaborn as sns
import os
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
from rdkit.Chem import Descriptors
PandasTools.RenderImagesInAllDataFrames(images=True)

# Load data
df = "https://drive.google.com/file/d/1zgjLnw07cK8Xn4wbYzBOaARyRKfWsS73/view?usp=sharing"
# ENLACE AQUÍ, se deben subir a una carpeta de drive *.csv
df='https://drive.google.com/uc?id=' + df.split('/')[-2]
CHEMBL = pd.read_csv(df)
CHEMBL.tail(2)

#df=pd.read_excel("/content/Prueba1.xlsx",  engine='openpyxl')
#print(f"DataFrame shape: {df.shape}.")
#df

#DataFrame info
CHEMBL.info()

# Statistics data
CHEMBL.describe()

""" `RDKit`"""

#mol object
PandasTools.AddMoleculeColumnToFrame(CHEMBL, "canonical smiles")
CHEMBL[0:2]

#Molecular descriptors
CHEMBL["HBD"] = CHEMBL["ROMol"].apply(Descriptors.NumHDonors)
CHEMBL["CSP3"] = CHEMBL["ROMol"].apply(Descriptors.FractionCSP3)
CHEMBL["NumRings"] = CHEMBL["ROMol"].apply(Descriptors.RingCount)
CHEMBL["HetAtoms"] = CHEMBL["ROMol"].apply(Descriptors.NumHeteroatoms)
CHEMBL["RotBonds"] = CHEMBL["ROMol"].apply(Descriptors.NumRotatableBonds)
#Visualizar columnas seleccionadas
CHEMBL[['HBD', 'CSP3','NumRings', 'HetAtoms', 'RotBonds' ]]

#Molecular descriptors
CHEMBL["MW"] = CHEMBL["ROMol"].apply(Descriptors.MolWt)
CHEMBL["HBA"] = CHEMBL["ROMol"].apply(Descriptors.NumHAcceptors)
CHEMBL["HBD"] = CHEMBL["ROMol"].apply(Descriptors.NumHDonors)
CHEMBL["logP"] = CHEMBL["ROMol"].apply(Descriptors.MolLogP)
CHEMBL["TPSA"] = CHEMBL["ROMol"].apply(Descriptors.TPSA)
CHEMBL["CSP3"] = CHEMBL["ROMol"].apply(Descriptors.FractionCSP3)
CHEMBL["NumRings"] = CHEMBL["ROMol"].apply(Descriptors.RingCount)
CHEMBL["HetAtoms"] = CHEMBL["ROMol"].apply(Descriptors.NumHeteroatoms)
CHEMBL["RotBonds"] = CHEMBL["ROMol"].apply(Descriptors.NumRotatableBonds)
#Visualizar columnas seleccionadas
CHEMBL[['MW', 'logP', 'TPSA', 'HBA', 'HBD', 'CSP3','NumRings', 'HetAtoms', 'RotBonds']]

from matplotlib import pyplot as plt
_df_56.plot(kind='scatter', x='TPSA', y='HBA', s=32, alpha=.8)
plt.gca().spines[['top', 'right',]].set_visible(False)

from matplotlib import pyplot as plt
_df_55.plot(kind='scatter', x='logP', y='TPSA', s=32, alpha=.8)
plt.gca().spines[['top', 'right',]].set_visible(False)

from matplotlib import pyplot as plt
_df_54.plot(kind='scatter', x='MW', y='logP', s=32, alpha=.8)
plt.gca().spines[['top', 'right',]].set_visible(False)

from matplotlib import pyplot as plt
_df_53.plot(kind='scatter', x='index', y='MW', s=32, alpha=.8)
plt.gca().spines[['top', 'right',]].set_visible(False)

# Print output
CHEMBL.to_excel('results.xlsx', index=False)
print

#Acumulated graphics
#Graph
hist=sns.histplot(x="HBA", data = CHEMBL)
#Format
hist.axes.set_title("HBA",fontsize=18)
hist.set_xlabel("Value",fontsize=14)
hist.set_ylabel("Frecuency",fontsize=14)
hist.tick_params(labelsize=14)

# Molecular Descriptors
CHEMBL_descriptors = CHEMBL[['MW', 'HBA', 'HBD', 'logP', 'TPSA', 'CSP3', 'NumRings', 'HetAtoms', 'RotBonds']]

# Color graphs
colors = ['maroon', 'blue', 'green', 'purple', 'orange', 'pink', 'brown', 'gray', 'olive']

# Graphs
fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(12, 12))
axes = axes.flatten()

for i, (ax, col) in enumerate(zip(axes, CHEMBL_descriptors.columns)):
    CHEMBL_descriptors[col].hist(ax=ax, color=colors[i])
    ax.set_title(col)
    ax.set_xlabel(f'{col} values')
    ax.set_ylabel('Frequency')

plt.tight_layout()
plt.show()