

from rdkit import Chem
from rdkit.Chem import AllChem
import streamlit as st
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
import pandas as pd
import matplotlib.pyplot as plt

def smile_to_mol(smiles):
    mol = Chem.MolFromSmiles(smiles) 
    return mol

def mol_to_image(mol):
    mol2img = Draw.MolToImage(mol)
    return mol2img

def mol_prop(mol):
    mol_formula = Descriptors.CalcMolFormula(mol)
    mol_wt = Descriptors.MolWt(mol)
    return mol_formula,mol_wt
    
def calculate_descriptor(descriptor, mol):
    if descriptor == "Molecular Weight":
        return Descriptors.MolWt(mol)
    elif descriptor == "LogP":
        return Descriptors.MolLogP(mol)
    # Add more descriptors as needed

def smiles_to_mols(smiles_list):
    smiles_list = pretreat_smiles(smiles_list)
    if isinstance(smiles_list,list) and len(smiles_list)> 1:
        mols_list = [smile_to_mol(smile) for smile in smiles_list]
    elif isinstance(smiles_list,list) and len(smiles_list)==1:
        mols_list = smile_to_mol(smiles_list[0])
    return mols_list

def mols_to_imgs(mols):
    if isinstance(mols,list) and len(mols)>1:
        mols2imgs = [mol_to_image(mol) for mol in mols]
    elif isinstance(mols,list) and len(mols)==1:
        mols2imgs = mol_to_image(mols[0])
    return mols2imgs

def pretreat_smiles(input_smiles):
    if isinstance(input_smiles, list):
        if len(input_smiles) > 1:
            return input_smiles
        elif len(input_smiles) == 1:
            return input_smiles[0]
    elif isinstance(input_smiles, str):
        if ',' in input_smiles:
            return input_smiles.split(',')
        else:
            return [input_smiles]
    else:
        return [input_smiles]


def getMolDescriptors(mols, missingVal=None):
    ''' calculate the full list of descriptors for a molecule
    
        missingVal is used if the descriptor cannot be calculated
    '''
    dfs = []
    res = {}
    merged_df = pd.DataFrame()
    
    if isinstance(mols,list) and len(mols)>1:
        for mol in mols:
            for nm,fn in Descriptors._descList:
        # some of the descriptor fucntions can throw errors if they fail, catch those here:
                try:
                    val = fn(mol)
                except:
            # print the error message:
                    import traceback
                    traceback.print_exc()
            # and set the descriptor value to whatever missingVal is
                    val = missingVal
                res[nm] = val
            df = pd.Series(res,name=Chem.MolToSmiles(mol))
            dfs.append(df)
        merged_df = pd.concat(dfs, axis=1)
        return merged_df
    elif isinstance(mols,list) and len(mols)==1:
        for nm,fn in Descriptors._descList:
            try:
                val = fn(mols)
            except:
                import traceback
                traceback.print_exc()
                val = missingVal
            res[nm] = val
            merged_df = pd.Series(res,name=Chem.MolToSmiles(mols[0]))
        return merged_df

def descriptors(mols):
    merged_df = getMolDescriptors(mols)
    return merged_df

import matplotlib.pyplot as plt

def plot_descriptors(df):
    num_columns = len(df.columns)

    # Create a single figure
    fig, ax = plt.subplots()

    for i, column in enumerate(df.columns):
        ax.scatter(range(len(df)), df[column], label=column, s=10)

    # Add labels and title
    ax.set_xlabel('Descriptors')
    ax.set_ylabel('Values')
    ax.set_title('Descriptors of Molecules')

    # Add a legend
    ax.legend()

    # Display the plot
    return fig

    