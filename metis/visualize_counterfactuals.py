# %%
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
import json
import cairosvg
from IPython.display import SVG

# %%


data = pd.read_csv("../data/scaffold_memory.csv")
with open("reinvent_connect/gen_counterfactual.json") as og_file:
    file_contents = json.loads(og_file.read())
mol = Chem.MolFromSmiles(file_contents["parameters"]["scaffolds"][0])
empty_mol = Chem.MolFromSmiles("")
# %%
gen_mols = [Chem.MolFromSmiles(smi) for smi in data.sample(20, replace=False).SMILES]
# %%

img = Draw.MolsToGridImage(
    [empty_mol, empty_mol, mol, empty_mol, empty_mol] + gen_mols,
    molsPerRow=5,
    subImgSize=(300, 300),
    useSVG=True,
)
img
# %%
