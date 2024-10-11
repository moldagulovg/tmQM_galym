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


print(rdkit.__version__)
print()
IPythonConsole.drawOptions.addAtomIndices = False
IPythonConsole.molSize = 300,300
RDLogger.DisableLog('rdApp.*')
IPythonConsole.ipython_3d = False




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

tmQM_properties = {}
for mol in tmQM_mols:
    CSD_code = mol.GetProp('CSD_code')
    prop = unload_properties(mol)
    tmQM_properties[CSD_code] = prop
keys_vals = [(a, b) for a,b in tmQM_properties.items()]
print('Recorded all properties from tmQM mols.')
# print(keys_vals[0])
# print(len(keys_vals))


###################################################################################
###################################################################################
###################################################################################
# Below code reads original CSD derived SDF files (conncetivity info encluded), loads properties and cross checks MND & stoichiometry

CSD_files_path = '/home/galymzhan/Documents/ccdc_structures_sdf'

filenames = os.listdir(CSD_files_path)
print(len(filenames))

# work in batches and save entire dataset into three SDF files
three_cuts = [filenames[:37000], filenames[37000:37000*2], filenames[37000*2:]]
print(len(three_cuts[0])+len(three_cuts[1])+len(three_cuts[2]) == len(filenames))
# sampled_files = sample(filenames, 10000)


############################
verbose = False
ignore_Hs = True
############################


incorrect_mnd_mols = []
incorrect_stoichio_mols = []
failed_upon_check = []
none_mols = []
succesful_mols = []
failed_error_entries = []

for i, sampled_files in enumerate(three_cuts):
    sdf_mols = []

    for sdf in tqdm(sampled_files):
        CSD_code = sdf[:-4]
        try:
            suppl = Chem.SDMolSupplier(os.path.join(CSD_files_path,f'{sdf}'), sanitize=False, removeHs=False)
            for mol in suppl:
                if mol is not None:
                    MND, true_stoichiometry = tmQM_properties[CSD_code]['MND'], tmQM_properties[CSD_code]['Stoichiometry']
                    
                    orgmetal_species = counter_ions_removal(mol)
                    for mol in orgmetal_species:
                        is_mnd_correct, is_stoichio_correct = check_coord_num(mol, MND), check_stoichiometry(mol, true_stoichiometry, ignore_Hs=ignore_Hs)
                        if is_stoichio_correct and is_mnd_correct: break
                        elif is_stoichio_correct: break
                    
                    if is_stoichio_correct:
                        pass
                    else:
                        mol_formula = rdMolDescriptors.CalcMolFormula(mol)
                        incorrect_stoichio_mols.append(tuple([CSD_code, mol_formula, true_stoichiometry]))
                        if verbose: 
                            print(mol_formula, "  >---<  ", true_stoichiometry)
                            print(CSD_code)
                            print()
                    if is_mnd_correct:
                        pass
                    else:
                        incorrect_mnd_mols.append(tuple([CSD_code, MND]))
                        if verbose:
                            print(CSD_code)
                        
                    
                    if is_mnd_correct == True and is_stoichio_correct == True:
                        mol = load_properties(mol, tmQM_properties[CSD_code])
                        sdf_mols.append(mol)
                        succesful_mols.append(CSD_code)
                    else:
                        failed_upon_check.append(CSD_code)
                else:
                    none_mols.append(CSD_code)
                    continue
        except:
            failed_error_entries.append(CSD_code)
            pass

    sdf_writer = Chem.SDWriter(f"/home/galymzhan/tmQM_galym/sdf_tmQM/tmQM_nonOpt_{i+1}.sdf")
    sdf_writer.SetKekulize(False)

    # Write each molecule to the file
    for mol in tqdm(sdf_mols):
        sdf_writer.write(mol)

    # Close the writer
    sdf_writer.close()

total = len(three_cuts[0])+len(three_cuts[1])+len(three_cuts[2])
print(f'total number of entries before {total}\n')

print(f'failed mols that raised errors: {len(failed_error_entries)}')
print(f'none mols: {len(none_mols)}')
print(f'incorrect stoichio: {len(incorrect_stoichio_mols)}')
print(f'incorrect mnd: {len(incorrect_mnd_mols)}')
print(f'incorrect stoichio & mnd : {len(failed_upon_check)}\n')

print(f'successful entries: {len(succesful_mols)}')

total_check = len(failed_error_entries) + len(none_mols) + len(failed_upon_check) + len(succesful_mols)
print(f'Do these numbers add-up to initial number?')
print(total_check == total)



###################################################################################
###################################################################################
###################################################################################
# Below code will convert all tmQM xyz to Mol (adds bonds) and corrects the connectivity at metal by MND


# tmQM_mols = tmQM_mols
# mnds_list = [int(mol.GetProp('MND')) for mol in tmQM_mols]

# temp_df = pd.DataFrame()
# temp_df['before_mol'] = tmQM_mols
# temp_df['mnd'] = mnds_list

# props = [unload_properties(mol) for mol in tmQM_mols]
# temp_df['props'] = props

# print('Start parallel correction of coord nums by MND.')
# temp_df['after_mol'] = temp_df.parallel_apply(row_mnd, axis=1)
# print('Finished correcting all coord numbers by MND.')

# for i in tqdm(range(len(temp_df))):
#     tmQM_mols[i] = load_properties(temp_df['after_mol'][i], temp_df['props'][i])
# print('Reloaded all of the properties.')


###################################################################################
###################################################################################
###################################################################################

