import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import psi4

from rdkit.Chem import Descriptors
import pandas as pd
import matplotlib.pyplot as plt
from psikit import Psikit

pk1 = Psikit(debug=True)
pk2 = Psikit(debug=True)


smiles = ["c1ccccc1","c1ccccc1N"]
mols = []
smiles_hs = []
for smile in smiles:
    mol = Chem.MolFromSmiles(smile)
    mol = Chem.AddHs(mol)
    smile_h = Chem.MolToSmiles(mol)
    smiles_hs.append(smile_h)
    mols.append(mol)


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


merged_df = descriptors(mols)
merged_df.to_csv("descriptors.csv")

pk1.read_from_smiles(smiles_hs[0])
pk1.optimize()

pk2.read_from_smiles(smiles_hs[1])
pk2.optimize()


merged_df_hs = descriptors([pk1.mol,pk2.mol])
merged_df_hs.to_csv("descriptors_optimized.csv")