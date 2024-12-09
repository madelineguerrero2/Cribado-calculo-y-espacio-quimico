
from IPython.utils import io
import tqdm.notebook
import os, sys, random
total = 100
with tqdm.notebook.tqdm(total=total) as pbar:
    with io.capture_output() as captured:
      # Instalar rdkit
      !pip -q install rdkit #.pypi==2021.9.4
      pbar.update(20)
      # Instalar Pillow
      !pip -q install Pillow
      pbar.update(40)
      # Instalar molplotly
      !pip install molplotly
      pbar.update(60)
      # Instalar jupyter-dash
      !pip install jupyter-dash
      pbar.update(80)
      # Instalar el diseño de aplicación dash
      !pip install dash-bootstrap-components
      pbar.update(100)

from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np
from rdkit.Chem import MACCSkeys
from scipy.spatial.distance import pdist, squareform

# compounds
url_data = "/content/Chembl.csv"
data = pd.read_csv(url_data)
data

"""## Morgan2 for compounds

---


"""

fps = pd.DataFrame([[int(y) for y in AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(x), 2, nBits=1024).ToBitString()] for x in data["canonical smiles"]])
fps

# Molecular similarity
SimMat_fps = 1 - pdist(fps[[x for x in range(1024)]], metric="jaccard")#Matriz de similitud_FDA
SimMat_fps
SimMat_fps = squareform(SimMat_fps)
SimMat_fps = 1-SimMat_fps

# Print matrix similarity
SimMat_fps = pd.DataFrame(SimMat_fps)
SimMat_fps

#SimMat_fps=fps

# Data Normalized
from sklearn.preprocessing import StandardScaler
data_std = SimMat_fps.copy()
data_std = StandardScaler().fit_transform(data_std)

from sklearn.decomposition import PCA
pca = PCA(n_components=2)
pca_results = pca.fit_transform(data_std)
pca_results

# Complementary information
label = data[["ID", "canonical smiles"]]
label = label.to_numpy()

# Concat arrays de numpy
arr = np.concatenate((label, pca_results), axis = 1)
# Crreate a new dataframe
pca_dataset = pd.DataFrame(data=arr, columns = ["ID", "SMILES",'component1', 'component2'])
pca_dataset.head(2)

"""# Cumulative variance graph


"""

# Determine explained variance using explained_variance_ration_ attribute
exp_var_pca = pca.explained_variance_ratio_
exp_var_pca

#Graph
import plotly.express as px
import molplotly
fig_pca = px.scatter(pca_dataset,
                            x='component1',
                            y='component2',
                            title='PCA',
                            labels={'component1': f'PC1 {round(exp_var_pca[0] * 100, 2)}%',
                                    'component2': f'PC2 {round(exp_var_pca[1] * 100, 2)}%'},
                            width=700,
                            height=500)

fig_pca.update_traces(marker=dict(color='orange'))

app_marker = molplotly.add_molecules(fig=fig_pca,
                                         df=pca_dataset,
                                         smiles_col='SMILES',
                                         title_col='ID')
#fig_pca.show()
#app_marker.run_server(mode='inline', port=8060, height=1000)
app_marker.run(port=8060)

"""## Morgan3 for compounds

---

"""

fps = pd.DataFrame([[int(y) for y in AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(x), 3, nBits=1024).ToBitString()] for x in data["canonical smiles"]])
fps

# Molecular similarity
SimMat_fps = 1 - pdist(fps[[x for x in range(1024)]], metric="jaccard")#Matriz de similitud_FDA
SimMat_fps
SimMat_fps = squareform(SimMat_fps)
SimMat_fps = 1-SimMat_fps

# Print matrix similarity
SimMat_fps = pd.DataFrame(SimMat_fps)
SimMat_fps

# Data Normalized
from sklearn.preprocessing import StandardScaler
data_std = SimMat_fps.copy()
data_std = StandardScaler().fit_transform(data_std)

from sklearn.decomposition import PCA
pca = PCA(n_components=2)
pca_results = pca.fit_transform(data_std)
pca_results

# Complementary information
label = data[["ID", "canonical smiles"]]
label = label.to_numpy()

# Concat arrays de numpy
arr = np.concatenate((label, pca_results), axis = 1)
# Crreate a new dataframe
pca_dataset = pd.DataFrame(data=arr, columns = ["ID", "SMILES",'component1', 'component2'])
pca_dataset.head(2)

# Determine explained variance using explained_variance_ration_ attribute
exp_var_pca = pca.explained_variance_ratio_
exp_var_pca

#Graph
import plotly.express as px
import molplotly
fig_pca = px.scatter(pca_dataset,
                            x='component1',
                            y='component2',
                            title='PCA',
                            labels={'component1': f'PC1 {round(exp_var_pca[0] * 100, 2)}%',
                                    'component2': f'PC2 {round(exp_var_pca[1] * 100, 2)}%'},
                            width=700,
                            height=500)

fig_pca.update_traces(marker=dict(color='purple'))

app_marker = molplotly.add_molecules(fig=fig_pca,
                                         df=pca_dataset,
                                         smiles_col='SMILES',
                                         title_col='ID')
#fig_pca.show()
#app_marker.run_server(mode='inline', port=8060, height=1000)
app_marker.run(port=8060)

"""## MACCS keys for compounds

---


"""

# Recalculate the MACCS keys
fps = pd.DataFrame([[int(y) for y in MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(x)).ToBitString()] for x in data['canonical smiles']])
# Molecular similarity
SimMat_fps = 1 - pdist(fps[[x for x in range(167)]], metric="jaccard")#Matriz de similitud_FDA
SimMat_fps
SimMat_fps = squareform(SimMat_fps)
SimMat_fps = 1-SimMat_fps

# Print matrix similarity
SimMat_fps = pd.DataFrame(SimMat_fps)
SimMat_fps

# Data Normalized
from sklearn.preprocessing import StandardScaler
data_std = SimMat_fps.copy()
data_std = StandardScaler().fit_transform(data_std)

from sklearn.decomposition import PCA
pca = PCA(n_components=2)
pca_results = pca.fit_transform(data_std)
pca_results

# Complementary information
label = data[["ID", "canonical smiles"]]
label = label.to_numpy()

# Concat arrays de numpy
arr = np.concatenate((label, pca_results), axis = 1)
# Crreate a new dataframe
pca_dataset = pd.DataFrame(data=arr, columns = ["ID", "SMILES",'component1', 'component2'])
pca_dataset.head(2)

# Determine explained variance using explained_variance_ration_ attribute
exp_var_pca = pca.explained_variance_ratio_
exp_var_pca

#Graph
import plotly.express as px
import molplotly
fig_pca = px.scatter(pca_dataset,
                            x='component1',
                            y='component2',
                            title='PCA',
                            labels={'component1': f'PC1 {round(exp_var_pca[0] * 100, 2)}%',
                                    'component2': f'PC2 {round(exp_var_pca[1] * 100, 2)}%'},
                            width=700,
                            height=500)

fig_pca.update_traces(marker=dict(color='red'))

app_marker = molplotly.add_molecules(fig=fig_pca,
                                         df=pca_dataset,
                                         smiles_col='SMILES',
                                         title_col='ID')
#fig_pca.show()
#app_marker.run_server(mode='inline', port=8060, height=1000)
app_marker.run(port=8060)
