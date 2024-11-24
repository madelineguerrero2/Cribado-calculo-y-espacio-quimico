
!pip install -q rdkit

!pip install -q molvs

import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import rdMolDescriptors
from molvs.standardize import Standardizer
from molvs.charge import Uncharger, Reionizer
from molvs.fragment import LargestFragmentChooser
from molvs.tautomer import TautomerCanonicalizer
from rdkit.Chem.rdmolops import GetFormalCharge, RemoveStereochemistry

#DNMT1
url_data = "/content/prueba 4.xlsx"
data = pd.read_excel(url_data, sheet_name='Hoja1')
# data = pd.read_csv(url_data)
data.head(2)  # Muestra las dos primeras filas de información, para corroborar que cargó la información correctamente

print(data.columns)

#Selecionar columna SMILES_chiral
data1 = data[['id', 'smiles']]
data1

set(list(data["id"]))

data1.columns = ['ID','isomeric smiles']

data1

#Seleccionar cinco compuestos
#DNMT1 = DNMT1.iloc[:5,:]

#Información del dataframe
data1.shape

# Definir funciones
STD = Standardizer() # Get the standardized version of a given SMILES string (canonical SMILES).
LFC = LargestFragmentChooser() # Select the largest fragment from a salt (ionic compound).
UC = Uncharger() # Charge corrections are applied to ensure, for example, that free metals are correctly ionized.
RI = Reionizer() # Neutralize molecule by adding/removing hydrogens.
TC = TautomerCanonicalizer()  # Return a tautormer “reasonable” from a chemist’s point, but isn’t guaranteed to be the most energetically favourable.

def pretreatment(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol == None:
            #If rdkit could not parse the smiles, returns Error 1
            return "Error 1"
        else:
            mol = STD(mol)
            mol = LFC(mol)

            allowed_elements = {"H","B","C","N","O","F","Si","P","S","Cl","Se","Br","I"}
            actual_elements = set([atom.GetSymbol() for atom in mol.GetAtoms()])
            if len(actual_elements-allowed_elements) == 0:
                mol = UC(mol)
                mol = RI(mol)
                # paso para remover estereoquímica
                # RemoveStereochemistry(mol)
                mol = TC(mol)
                return Chem.MolToSmiles(mol)
            else:
                # If molecule contains other than the allowed elements, return "Error 2"
                return "Error 2"
    except:
        return "Something else was found"

data1["canonical smiles"] = [pretreatment(x) for x in data1["isomeric smiles"]]
data1

# Borrar smiles que rdkit no puede leer
data1 = data1[data1["canonical smiles"] != "Error 1"]
# Borrar smiles que contienen átomos no permitidos
data1 = data1[data1["canonical smiles"] != "Error 2"]
# Borrar otros errores
data1 = data1[data1["canonical smiles"] != "Something else was found"].reset_index(drop=True)

# Borrar duplicados
data1 = data1.drop_duplicates(subset=["isomeric smiles"], keep="first").reset_index(drop=True)
data1

# Borrar filas con NaN (Not a Number) y valores None del DataFrame.
data1 = data1.dropna()

data1

data1.to_csv("Pruba5.csv", sep=",", index=False)