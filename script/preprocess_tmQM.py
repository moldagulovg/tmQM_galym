import rdkit
print(rdkit.__version__)
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDetermineBonds
IPythonConsole.drawOptions.addAtomIndices = False
IPythonConsole.molSize = 300,300
from IPython.display import display
from PIL import ImageDraw, ImageFont
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import pandas as pd
import numpy as np

from pandarallel import pandarallel
pandarallel.initialize(nb_workers=40, progress_bar=False)



IPythonConsole.ipython_3d = False
import py3Dmol
import sys
import os

from tqdm import tqdm


import rdmetallics

from IPython.display import display


from util import *

#########################################################
#########################################################
#########################################################



#########################################################
#########################################################
#########################################################






###################################################################################
###################################################################################
###################################################################################


dir_path = '/home/galymzhan/tmQM_galym/tmQM'
tmQM_mols = load_tmQM_dataset(dir_path)
print('Loaded the initial dataset.')


###################################################################################
###################################################################################
###################################################################################


tmQM_mols = tmQM_mols
mnds_list = [int(mol.GetProp('MND')) for mol in tmQM_mols]

temp_df = pd.DataFrame()
temp_df['before_mol'] = tmQM_mols
temp_df['mnd'] = mnds_list

props = [unload_properties(mol) for mol in tmQM_mols]
temp_df['props'] = props


temp_df['after_mol'] = temp_df.parallel_apply(row_mnd, axis=1)

for i in tqdm(range(len(temp_df))):
    tmQM_mols[i] = load_properties(temp_df['after_mol'][i], temp_df['props'][i])



###################################################################################
###################################################################################
###################################################################################

ligands_from_tmQM = []
for mol in tqdm(tmQM_mols):
    if mol is None:
        continue

    fragments = rdmetallics.Get_Ligands_From_Organometallic(mol)
    if len(fragments) >= 2:
        for ligand in fragments[1:]:
            metal_symbol = rdmetallics.get_metal(mol)
            complex_mol = substitute_dummy_atom(ligand, metal_symbol)
            ligands_from_tmQM.append({'sub_complex': complex_mol, 'ligand': remove_dummy_atom(ligand), 'metal_symbol': metal_symbol, })



###################################################################################
###################################################################################
###################################################################################




print('number of all ligand-metal-complex dictionaries: ', len(ligands_from_tmQM))

# Use a dictionary to store unique molecules based on their canonical SMILES
unique_ligands = {}
smiles_from_tmQM = []

for ligand_dict in tqdm(ligands_from_tmQM):
    # Generate canonical SMILES
    ligand_smiles = Chem.MolToSmiles(ligand_dict['ligand'], canonical=True, allHsExplicit=True)
    metal_symbol = ligand_dict['metal_symbol']
    complex_smiles = Chem.MolToSmiles(ligand_dict['sub_complex'], canonical=True, allHsExplicit=True)
    # Add to dictionary (key: SMILES, value: molecule)
    if ligand_smiles not in list(unique_ligands.keys()):
        unique_ligands[ligand_smiles] = [complex_smiles,]
    else:
        unique_ligands[ligand_smiles] += [complex_smiles,]
    smiles_from_tmQM.append({'ligand_smiles': ligand_smiles, 'metal_symbol': metal_symbol, 'subcomplex_smiles': complex_smiles})
        
# Convert the unique molecules back to a list or keep as a set
print('number of unique ligand smiles: ', len(list(unique_ligands.keys())))
max_num_complex_per_ligand = 0
its_smiles = ''
second_smiles = ''
second = 0
for key, value in unique_ligands.items():
    if max_num_complex_per_ligand <= len(value):
        second = max_num_complex_per_ligand
        max_num_complex_per_ligand = len(value)
        second_smiles = its_smiles
        its_smiles = key
print(f'maximum number of complexes per ligand: {max_num_complex_per_ligand}')
print()
print()
print()


###################################################################################
###################################################################################
###################################################################################


import pandas as pd
dataset_batch_df = pd.DataFrame(smiles_from_tmQM)
print(f'len of df before {len(dataset_batch_df)}')
dataset_batch_df = dataset_batch_df.drop_duplicates(subset='subcomplex_smiles', keep='first')
print(f'len of df after {len(dataset_batch_df)}')
print(dataset_batch_df.info())
print()
print()
print()
#get the count of uniques
unique_counts = dataset_batch_df.nunique()
print('count of uniques:\n',unique_counts)



###################################################################################
###################################################################################
###################################################################################

path = '/home/galymzhan/tmQM_galym/tmQM'
file_path = os.path.join(path, 'ligand_metal_complex.csv')
dataset_batch_df.to_csv(file_path)
print(f'saved to {file_path}')
