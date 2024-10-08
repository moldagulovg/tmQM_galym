import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDetermineBonds

from IPython.display import display
from PIL import ImageDraw, ImageFont
from rdkit import RDLogger


import pandas as pd
import numpy as np
import py3Dmol
import sys
import os
from tqdm import tqdm
import rdmetallics
from IPython.display import display
from pandarallel import pandarallel
pandarallel.initialize(nb_workers=40, progress_bar=False)

RDLogger.DisableLog('rdApp.*')
IPythonConsole.drawOptions.addAtomIndices = False
IPythonConsole.molSize = 300,300
IPythonConsole.ipython_3d = False
print(rdkit.__version__)


def parse_mol_properties(prop_string):
    '''
    Extracts assigned properties written in XYZ file.
    q - charge, S - spin, MND - ?
    '''
    
    # Split the string by '|' to get individual properties
    properties = prop_string.split('|')
    
    # Create a dictionary to store the property names and values
    prop_dict = {}
    
    # Iterate over each property and split by '=' to get name and value
    for prop in properties:
        if prop == properties[-1]:
            prop_dict['Source'] = prop.strip()
        else:
            # Strip any leading/trailing spaces and split by '='
            name, value = prop.split('=')
            name, value = name.strip(), value.strip()
            # Strip spaces from name and value, then add to dictionary
            prop_dict[name] = value
    
    return prop_dict

def read_xyz_file(file_path):
    with open(file_path, 'r') as file:
        content = file.read().strip().split("\n")
    
    molecules = []
    i = 0
    while i < len(content):
        num_atoms = int(content[i].strip())
        mol_props = parse_mol_properties(content[i+1])
        xyz_block = "\n".join(content[i:i+num_atoms+2])  # +2 to account for title and atom count
        raw_mol = Chem.MolFromXYZBlock(xyz_block)
        if raw_mol:
            for prop_name, value in mol_props.items():
                raw_mol.SetProp(prop_name, value)
                if prop_name == 'q':
                    charge = value
            mol = Chem.Mol(raw_mol)
            rdDetermineBonds.DetermineConnectivity(mol, charge=int(charge))
            # rdDetermineBonds.DetermineBondOrders(mol, charge=int(charge))
            molecules.append(mol)
        else:
            molecules.append(None)
        i += num_atoms + 3
    
    return molecules

def load_tmQM_dataset(dir_path, subset=None, verbose=True):
    
    # make this tmQM dataset load as a python package through, import something-something
    
    # TODO add checker if directory exists
    # TODO include checker of tmQM_X1-3.xyz file presence.
    tmQM_mols = []
    for i in [1,2,3]:
        file_path = os.path.join(dir_path, f'tmQM_X{i}.xyz')
        tmQM_mols+=read_xyz_file(file_path)
        
    if subset:
        subset_tmQM_mols = []
        for entry_idx in subset:
            subset_tmQM_mols.append(tmQM_mols[entry_idx])
        tmQM_mols = subset_tmQM_mols

    if verbose:
        total_mols, num_valid_mols, num_invalid_mols = len(tmQM_mols), len([mol for mol in tmQM_mols if mol is not None]), len([mol for mol in tmQM_mols if mol is None])
        print(f'total num mols {total_mols}\nvalid mols {num_valid_mols}\nNone mols {num_invalid_mols}\n')
    
    return tmQM_mols

# correct the coordination number to match MND
def correct_coord_num(mol, coord_num):
    if mol is None:
        raise ValueError(f"{mol} does not contain info about any molecule")

    # Get the 3D coordinates of atoms in the first conformer
    conformer = mol.GetConformer(0)
    
    atoms = mol.GetAtoms()
    metals = []
    
    # Identify transition metal atoms
    for atom in atoms:
        if rdmetallics.is_transition_metal(atom):
            metals.append(atom)
    
    editable_mol = Chem.RWMol(mol)
    
    # Process each metal and its neighbors
    for metal in metals:
        nbrs = metal.GetNeighbors()
        metal_idx = metal.GetIdx()
        
        distances = []
        
        # Calculate distances between the metal and its neighbors
        for atom in atoms:
            at_idx = atom.GetIdx()
            if at_idx == metal_idx:
                continue

            pos1 = np.array(conformer.GetAtomPosition(metal_idx))
            pos2 = np.array(conformer.GetAtomPosition(at_idx))
            
            # Calculate Euclidean distance
            distance = np.linalg.norm(pos1 - pos2)
            distances.append((distance, (metal_idx, at_idx)))
        
        distances.sort(key=lambda x: x[0])

        for i, dist_tuple in enumerate(distances):
            metal_idx, at_idx = dist_tuple[1][0], dist_tuple[1][1]
            if (i+1) > coord_num:
                if mol.GetBondBetweenAtoms(metal_idx, at_idx) is not None:
                    editable_mol.RemoveBond(metal_idx, at_idx)
            else:
                if mol.GetBondBetweenAtoms(metal_idx, at_idx) is None:
                    editable_mol.AddBond(metal_idx, at_idx)
                    
    mol = editable_mol.GetMol()  # Convert back to regular Mol
    
    return mol

def row_mnd(row):
    mol = correct_coord_num(row['before_mol'], row['mnd'])
    return mol

def unload_properties(mol):
    # Save assigned properties
    properties = mol.GetPropsAsDict()  # Get properties as a dictionary
    return properties

def load_properties(mol, properties):
    # Transfer the properties
    for key, value in properties.items():
        mol.SetProp(key, str(value))  # Set each property on the copied molecule
    return mol




def ligand_substructure_search(mol_ligand, catalyst_mol):
    has_substructure = False
    
    if catalyst_mol.HasSubstructMatch(mol_ligand):
        has_substructure = True
    else:
        print("No substructure match observed")
        return has_substructure, None, None, None
    
    match_metal_nbrs_idxs = []
    substruct_matches_catalyst = list(catalyst_mol.GetSubstructMatches(mol_ligand))
    
    large_to_small_mapping = {}
    for k, atom_mapping in enumerate(substruct_matches_catalyst):
        # Create a dictionary to map substructure atom indices to the larger molecule's indices
        large_to_small_mapping_ind = {at_idx: i for i, at_idx in enumerate(atom_mapping)}
        large_to_small_mapping[k] = large_to_small_mapping_ind
    
    for match in substruct_matches_catalyst:
        dentate_idxs_match = []
        atoms = [catalyst_mol.GetAtomWithIdx(idx) for idx in match]
        for atom in atoms:
            nbrs = atom.GetNeighbors()
            for nbr in nbrs:
                if rdmetallics.is_transition_metal(nbr):
                    atom_idx = atom.GetIdx()
                    dentate_idxs_match.append(atom_idx)
        match_metal_nbrs_idxs.append(tuple(dentate_idxs_match))
    
    print(large_to_small_mapping)
    
    dentate_in_ligand_idxs = []
    for match_i, dentate_idxs_match in enumerate(match_metal_nbrs_idxs):
        dentate_in_ligand_idxs.append(tuple([large_to_small_mapping[match_i][i] for i in dentate_idxs_match]))
    
    print('large_to_small_atom_mapping: \n', large_to_small_mapping)

    
    mols = [mol_ligand, catalyst_mol, catalyst_mol, catalyst_mol, mol_ligand]
    
    highlighted_atoms_matches = []
    for x in substruct_matches_catalyst:
        highlighted_atoms_matches += list(x)
        
    highlighted_match_metal_nbrs_idxs = []
    for y in match_metal_nbrs_idxs:
        highlighted_match_metal_nbrs_idxs += list(y)
    
    highlights = [
        [],  # Ligand atoms
        [],  # Catalyst atoms
        highlighted_atoms_matches,  # matched atoms
        highlighted_match_metal_nbrs_idxs,
        dentate_in_ligand_idxs[0],
    ]
    # print(highlights)
    
    # draw_options = Draw.MolDrawOptions()
    # draw_options.addAtomIndices = True 
    
    # # Display the molecules with highlighted atoms
    # img = Draw.MolsToGridImage(mols, molsPerRow=5, highlightAtomLists=highlights)
    # display(img)
    
    return has_substructure, substruct_matches_catalyst, dentate_in_ligand_idxs, match_metal_nbrs_idxs


def remove_dummy_atom(mol):
    mol_copy = Chem.Mol(mol)
    
    rwmol = Chem.RWMol(mol_copy)
    atoms = rwmol.GetAtoms()
    for atom in atoms:
        if atom.GetSymbol() == '*':    
            dummy_idx = atom.GetIdx()
    rwmol.RemoveAtom(dummy_idx)
    new_mol = rwmol.GetMol()
    mol_props = new_mol.GetPropsAsDict()
    
    for prop in mol_props.keys():
        if mol_copy.HasProp(prop):
            mol_copy.ClearProp(prop)
            
    return new_mol


def substitute_dummy_atom(mol, atom_symbol):
    mol_copy = Chem.Mol(mol)
    rwmol = Chem.RWMol(mol_copy)
    
    # Create a periodic table to convert symbol to atomic number
    pt = Chem.GetPeriodicTable()
    atomic_num = pt.GetAtomicNumber(atom_symbol)
    
    # Find the index of the dummy atom (usually atomic number 0 for "*")
    for atom in rwmol.GetAtoms():
        if atom.GetAtomicNum() == 0:  # Dummy atom has atomic number 0
            atom_idx = atom.GetIdx()
            break
    else:
        raise ValueError("No dummy atom found in the molecule")
    
    # Replace the dummy atom with the specified atom
    new_atom = Chem.Atom(atomic_num)  # Create the new atom
    rwmol.ReplaceAtom(atom_idx, new_atom)  # Replace atom at the index
    new_mol = rwmol.GetMol()
    return new_mol