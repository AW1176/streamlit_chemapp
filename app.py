
from rdkit import Chem
from rdkit.Chem import Draw
from psikit import Psikit
import streamlit as st
import pandas as pd
import utils
#st.session_state
pk = Psikit(debug=True)

def pretreat_input(text_input):
    input_list = []
    for text in text_input.split(","):
        input_list.append(text)
    return input_list


@st.cache_data
def optmization(smile,basis_sets):
    pk.read_from_smiles(smile)
    opt = pk.optimize(basis_sets=basis_sets)
    return opt

# use this to make a grid of container # not implemented here
def make_grid(cols,rows):
    grid = [0]*cols
    for i in range(cols):
        with st.container():
            grid[i] = st.columns(rows)
    return grid

def opt_callback():
    if 'opt' in st.session_state: st.session_state['opt']= opt

def smiles_callback():
    if 'selected_smiles' in st.session_state: st.session_state['selected_smiles'] = options1

def mols_callback():
    if 'selected_mols' in st.session_state: st.session_state['selected_mols'] = selected_mols
    
def df_merged_callback():
    if 'df_merged' in st.session_state: st.session_state['df_merged'] = df_merged
    st.write("i am still here")
    
col1,col2, col2a= st.columns([3,3,2])
col3,col3a,col4= st.columns([3,2,1])
col5a,col5b = st.columns(2)
col6 = st.columns(1)[0]
col7 = st.columns(1)[0]



if 'opt' not in st.session_state: st.session_state['opt'] = []
if 'images_mol' not in st.session_state: st.session_state['images_mol'] = []
if 'selected_smiles' not in st.session_state: st.session_state['selected_smiles'] = []
if 'selected_mols' not in st.session_state: st.session_state['selected_mols'] = []
if 'df_merged' not in st.session_state: st.session_state['df_merged'] = pd.DataFrame()

with col1:
    st.write("#### Pleas Enter you SMILES here")
    smiles_input = st.text_input(label="",help="Pleas put a common for more than one SMILES")
    pretreated_smiles_input = pretreat_input(smiles_input)
    pretreated_mols_input = [Chem.MolFromSmiles(i) for i in pretreated_smiles_input]
    st.write(smiles_input)
    

with col2:
    st.write("#### Pleas Select SMILES to view")
    selected_smiles = st.selectbox('Pleas select SMILES', options=pretreated_smiles_input)
    
    
with col2a:
    if selected_smiles:
        selected_mols = Chem.MolFromSmiles(selected_smiles)
        mol_image = Draw.MolToImage(selected_mols,size=(100, 100))
        st.image(mol_image)

with col3:
    st.write("#### Optimize Molecule")
    if pretreated_smiles_input:
        options1 = st.multiselect('Select SMILES for optimizations',pretreated_smiles_input,on_change=smiles_callback)
        #smiles_callback()
        selected_mols = [Chem.MolFromSmiles(smile) for smile in options1] 
        mols_callback()

with col3a:
    st.write("""""")
    basis_sets = st.selectbox("#### Select Basis Set", options=["HF/sto-3g","B3LYP/6-31G*", "scf/sto-3g"])
        
with col4:
    opt=[]   
    if st.button("**Optimize Selected molecule**"):
        for smile in options1:
            opt_eng = optmization(smile=smile,basis_sets=basis_sets)
            opt.append(opt_eng)
        opt_callback()

with col5a:
    if st.session_state.opt:
        st.write(f"#### Optimized energy for SMILES {options1} are \n= **{st.session_state.opt}**")
        with col5b:
            images_mol = st.image(Draw.MolsToGridImage(selected_mols))

with col6:
    if st.session_state.opt:
        st.write("# Calculate Descriptors")
        df_merged = utils.descriptors(selected_mols)
        col6a, col6b, col6c = st.columns(3)
        with col6a:
            if st.button("Show descriptors",on_click=df_merged_callback):
                with col7:
                    st.dataframe(df_merged,use_container_width=True)
        with col6b:
            if st.button("# Plot descriptors"):
                with col7:
                    fig = utils.plot_descriptors(df_merged)
                    st.pyplot(fig)
        with col6c:
            if st.button("Export Descriptors"):
                df_merged.to_csv("descriptors.csv")
                st.write("Exported CSV of descriptors")