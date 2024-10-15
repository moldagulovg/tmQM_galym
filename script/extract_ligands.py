<<<<<<< HEAD
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem import rdMolDescriptors
from rdkit import RDLogger

from IPython.display import display
from PIL import ImageDraw, ImageFont


import pandas as pd
import numpy as np
from tqdm import tqdm

import py3Dmol
from collections import defaultdict
import re
from random import sample
import sys
import os
import rdmetallics

from IPython.display import display


from util import *
from global_params import sdf_tmQM_dataset_path, extracted_ligands_path


print(rdkit.__version__)
print()
IPythonConsole.drawOptions.addAtomIndices = False
IPythonConsole.molSize = 300,300
RDLogger.DisableLog('rdApp.*')
IPythonConsole.ipython_3d = False

###################################################################################
###################################################################################
###################################################################################

sdf_tmQM_mols = []

for j in [0,1,2]:
    suppl = Chem.SDMolSupplier(os.path.join(sdf_tmQM_dataset_path, f"tmQM_nonOpt_{j+1}.sdf"), sanitize=False, removeHs=False)
    for mol in tqdm(suppl):
        sdf_tmQM_mols.append(mol)

###################################################################################
###################################################################################
###################################################################################

ligands_from_tmQM = []
for mol in tqdm(sdf_tmQM_mols):
=======
ligands_from_tmQM = []
for mol in tqdm(tmQM_mols):
>>>>>>> 3dc62a990f2991cd2ffa05beef6401a878520e7d
    if mol is None:
        continue

    fragments = rdmetallics.Get_Ligands_From_Organometallic(mol)
    if len(fragments) >= 2:
        for ligand in fragments[1:]:
            metal_symbol = rdmetallics.get_metal(mol)
            complex_mol = substitute_dummy_atom(ligand, metal_symbol)
            ligand = remove_dummy_atom(ligand)
<<<<<<< HEAD
            for at in ligand.GetAtoms():
                at.SetAtomMapNum(at.GetIdx())
            ligands_from_tmQM.append({'sub_complex': complex_mol, 'ligand': ligand, 'metal_symbol': metal_symbol, })

=======
            ligands_from_tmQM.append({'sub_complex': complex_mol, 'ligand': ligand, 'metal_symbol': metal_symbol, })


>>>>>>> 3dc62a990f2991cd2ffa05beef6401a878520e7d
###################################################################################
###################################################################################
###################################################################################

<<<<<<< HEAD
print('number of all ligand-metal-complex dictionaries: ', len(ligands_from_tmQM))

# Use a dictionary to store unique molecules based on their canonical SMILES
unique_ligands = {}
smiles_from_tmQM = []
num_noncanonical_ligand_smiles = []

for ligand_dict in tqdm(ligands_from_tmQM):
    # Generate canonical SMILES
    try:
        ligand_smiles = Chem.MolToSmiles(ligand_dict['ligand'], canonical=True, allHsExplicit=True)
    except:
        ligand_smiles = Chem.MolToSmiles(ligand_dict['ligand'], canonical=False, allHsExplicit=True)
        num_noncanonical_ligand_smiles.append(ligand_smiles)
        print(f'exception #{len(num_noncanonical_ligand_smiles)}')
    metal_symbol = ligand_dict['metal_symbol']
    try:
        complex_smiles = Chem.MolToSmiles(ligand_dict['sub_complex'], canonical=True, allHsExplicit=True)
    except:
        complex_smiles = Chem.MolToSmiles(ligand_dict['sub_complex'], canonical=False, allHsExplicit=True)
    # Add to dictionary (key: SMILES, value: molecule)
    if ligand_smiles not in list(unique_ligands.keys()):
        unique_ligands[ligand_smiles] = [complex_smiles,]
    else:
        unique_ligands[ligand_smiles] += [complex_smiles,]
    smiles_from_tmQM.append({'ligand_smiles': ligand_smiles, 'metal_symbol': metal_symbol, 'subcomplex_smiles': complex_smiles})
=======



print('number of all ligand-metal-complex dictionaries: ', len(ligands_from_tmQM))

# # Use a dictionary to store unique molecules based on their canonical SMILES
# unique_ligands = {}
# smiles_from_tmQM = []

# print('starting ligand labelling')
# print(f'{len(ligands_from_tmQM)} ligands queued')
# progress = []

# for i, ligand_dict in tqdm(enumerate(ligands_from_tmQM)):
#     # Generate canonical SMILES
#     ligand = ligand_dict['ligand']
#     subcomplex = ligand_dict['sub_complex']

#     dent = None

#     ################
#     a,b,c,d = ligand_substructure_search(ligand, subcomplex)
#     if a:
#         dent = c[0]
#     ################
    
#     for atom in ligand.GetAtoms():
#         atom.SetAtomMapNum(atom.GetIdx())
#     ligand_smiles = Chem.MolToSmiles(ligand, canonical=True, allHsExplicit=True, )
#     complex_smiles = Chem.MolToSmiles(subcomplex, canonical=True, allHsExplicit=True)

#     metal_symbol = ligand_dict['metal_symbol']

#     # Add to dictionary (key: SMILES, value: molecule)
#     if ligand_smiles not in list(unique_ligands.keys()):
#         unique_ligands[ligand_smiles] = [complex_smiles,]
#     else:
#         unique_ligands[ligand_smiles] += [complex_smiles,]
#     smiles_from_tmQM.append({'ligand_smiles': ligand_smiles, 'metal_symbol': metal_symbol, 'subcomplex_smiles': complex_smiles, 'dentates': dent})

#     percent = i / len(ligands_from_tmQM) * 100
#     decimal_percents = percent // 10 * 10
#     if decimal_percents not in progress:
#         progress.append(decimal_percents)
#         print(decimal_percents)


>>>>>>> 3dc62a990f2991cd2ffa05beef6401a878520e7d
        
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
<<<<<<< HEAD

print()
print()

# TODO save 'unique_ligands' dictionary as CSV file
# # Check the lengths of the lists
# max_len = max(len(v) for v in unique_ligands.values())

# # Adjust lists to have the same length (filling with None or any other value)
# for key, value in unique_ligands.items():
#     if len(value) < max_len:
#         unique_ligands[key] = value + [None] * (max_len - len(value))

# u_lig_df = pd.DataFrame(unique_ligands)
# file_path = os.path.join(extracted_ligands_path, 'unique_ligands_n_complexes.csv')
# u_lig_df.to_csv(file_path, index=False)
# print(f'saved unique ligands and all their complexes to {file_path}')
=======
print()
print()
print()

>>>>>>> 3dc62a990f2991cd2ffa05beef6401a878520e7d

###################################################################################
###################################################################################
###################################################################################


import pandas as pd
dataset_batch_df = pd.DataFrame(smiles_from_tmQM)
print(f'len of df before {len(dataset_batch_df)}')
dataset_batch_df = dataset_batch_df.drop_duplicates(subset='subcomplex_smiles', keep='first')
print(f'len of df after {len(dataset_batch_df)}')
print(dataset_batch_df.info())
<<<<<<< HEAD

=======
print()
>>>>>>> 3dc62a990f2991cd2ffa05beef6401a878520e7d
print()
print()
#get the count of uniques
unique_counts = dataset_batch_df.nunique()
print('count of uniques:\n',unique_counts)
<<<<<<< HEAD
print()
print()
=======


>>>>>>> 3dc62a990f2991cd2ffa05beef6401a878520e7d

###################################################################################
###################################################################################
###################################################################################

<<<<<<< HEAD
file_path = os.path.join(extracted_ligands_path, 'ligand_metal_complex.csv')
dataset_batch_df.to_csv(file_path)
print(f'saved all ligand-metal-complex to {file_path}')


=======
path = '/home/galymzhan/tmQM_galym/tmQM'
file_path = os.path.join(path, 'ligand_metal_complex.csv')
dataset_batch_df.to_csv(file_path)
print(f'saved to {file_path}')
>>>>>>> 3dc62a990f2991cd2ffa05beef6401a878520e7d
