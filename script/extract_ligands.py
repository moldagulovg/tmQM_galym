ligands_from_tmQM = []
for mol in tqdm(tmQM_mols):
    if mol is None:
        continue

    fragments = rdmetallics.Get_Ligands_From_Organometallic(mol)
    if len(fragments) >= 2:
        for ligand in fragments[1:]:
            metal_symbol = rdmetallics.get_metal(mol)
            complex_mol = substitute_dummy_atom(ligand, metal_symbol)
            ligand = remove_dummy_atom(ligand)
            ligands_from_tmQM.append({'sub_complex': complex_mol, 'ligand': ligand, 'metal_symbol': metal_symbol, })


###################################################################################
###################################################################################
###################################################################################




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
