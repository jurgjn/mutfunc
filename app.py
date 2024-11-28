
import ast, gzip, random, os, tempfile, time, sqlite3, urllib.request
import numpy as np, matplotlib, matplotlib.colors, matplotlib.pyplot as plt, seaborn as sns, pandas as pd, streamlit as st
import py3Dmol, stmol
import streamlit as st
import pandas as pd
import requests
from unipressed import IdMappingClient
from Bio.Data.IUPACData import protein_letters_3to1
import time
#, streamlit_ext as ste, st_aggrid, py3Dmol, stmol

st.cache_resource.clear()

st.set_page_config(
    page_title='mutfunc - prototype',
    page_icon='🔬',
    layout='wide',
)

def query_missense(variants):
    #https://stackoverflow.com/questions/28735213/pandas-read-sql-with-a-list-of-values-for-where-condition
    with sqlite3.connect('resources/missense_23.11.1.sqlite') as db:
        df_ = pd.read_sql_query(sql='SELECT * FROM missense WHERE variant_id in ' + str(tuple(variants)), con=db)
    return df_.replace({
        'pred_ddg': -99,
        'pocketscore': -99,
        'pocketrank': -99,
        'interface': -99,
        'interface_strict': -99,
        'freq': -99,
    }, np.nan)

# Function to query the Ensembl REST API
def fetch_variant_info(variants):
    """
    Fetch variant information from the Ensembl Variant Recoder endpoint.

    Parameters:
    - variants: List of variant identifiers (e.g., rsIDs, HGVS notations).

    Returns:
    - JSON response from the API.
    """
    server = "https://rest.ensembl.org"
    endpoint = "/variant_recoder/homo_sapiens"
    headers = {"Content-Type": "application/json"}
    
    if not variants:
        return []

    # For batch queries, use POST request
    data = {"ids": variants}
    response = requests.post(f"{server}{endpoint}", headers=headers, json=data)

    # Debugging: Print raw response data
    st.subheader("Debugging: Raw API Response")
    st.json(response.json())  # Display raw JSON response in the app

    if response.status_code != 200:
        st.error(f"Error {response.status_code}: {response.json().get('error', 'Unknown error')}")
        return []

    return response.json()


def query_prot_ids(ensemble_ids):
    """
    Query UniProt for protein IDs mapped from Ensembl Protein IDs.

    Parameters:
    - ensemble_ids: List of Ensembl Protein IDs with variants (e.g., "ENSP00000483018.1:p.Gly229Asp").

    Returns:
    - List of UniProt internal IDs with variants.
    """
    # Split the IDs into ENSP and variant
    ensemble_ids = [element.split(":") for element in ensemble_ids]
    variants = [e[1][2:] for e in ensemble_ids] # remove ".p"
    ensemble_ids = [e[0] for e in ensemble_ids]

    # Submit a mapping request to UniProt
    request = IdMappingClient.submit(source="Ensembl_Protein", dest="UniProtKB-Swiss-Prot", ids=set(ensemble_ids))
    result = list(request.each_result())
    
    st.write(result)
    #mapping = {e["from"]: e["to"] for e in result}

    # Get UniProt IDs and convert them to internal IDs with variants
    uniprot_ids = set({e["to"] for e in result})#set(mapping.values())
    st.write(uniprot_ids)
    internal_ids = [convert_to_internal_id(uniprot_id, variants[0]) for uniprot_id in uniprot_ids]

    return internal_ids  # Return results if successful

def convert_to_internal_id(uniprot_id, variant):
    position = variant[3:-3]
    short_variant = protein_letters_3to1[variant[-3:]]+position+protein_letters_3to1[variant[:3]]
    return uniprot_id+"/"+short_variant

# Function to process the API response into a structured DataFrame
def process_variant_info(variant_info):
    # get the proteins 
    rows = []
    for variant in variant_info:
        for nucleotide, details in variant.items():
            input_variant = details.get("input", "N/A")
            hgvsp = details.get("hgvsp", [])
            protein_ids = query_prot_ids(hgvsp)
            """for protein in protein_ids:
                rows.append({
                            "Input Variant": input_variant,
                            "Nucleotide Change": nucleotide,
                            "HGVSp": hgvsp,
                            "Protein ID": protein,
                        })
            """

    return protein_ids #pd.DataFrame(rows)

# Sidebar for input
st.sidebar.header("Input Variants")
uploaded_file = st.sidebar.file_uploader("Upload a file containing genomic variants", type=["txt", "csv"])
manual_input = st.text_area("Enter genomic variants (one per line)")

variants = []

# Handle uploaded file
if uploaded_file is not None:
    try:
        # Assuming the file is a plain text file or CSV with one column of variants
        file_content = pd.read_csv(uploaded_file, header=None)
        variants.extend(file_content[0].tolist())
    except Exception as e:
        st.sidebar.error("Failed to read the uploaded file. Please ensure it contains variants in one column.")

# Handle manual input
if manual_input:
    variants.extend([line.strip() for line in manual_input.split("\n") if line.strip()])

if st.sidebar.button("Analyze Variants"):
    if not variants:
        st.warning("No variants provided. Please upload a file or enter variants manually.")
    else:
        with st.spinner("Fetching variant information from Ensembl..."):
            try:
                # Fetch variant information from the Variant Recoder API
                variant_info = fetch_variant_info(variants)
                if not variant_info:
                    st.error("No information was returned for the provided variants.")
                else:
                    # Process the JSON response into a DataFrame
                    results_df = process_variant_info(variant_info)
                    #st.subheader("Processed Variant Results")
                    #st.dataframe(results_df)
                    #st.download_button("Download Results as CSV", results_df.to_csv(index=False), "variant_results.csv", "text/csv")
            except Exception as e:
                st.error(f"An error occurred: {e}")

    st.write('## Variants')
    #df_ = query_missense(['P09874/S568F', 'P09874/D678H', 'Q96NU1/R28Q', 'P00451/G41C'])
    df_ = query_missense(results_df)
    event = st.dataframe(
        df_,
        column_order = ['variant_id', 'am_pathogenicity', 'am_class', 'pred_ddg', 'pocketscore', 'interface'],
        #column_config=column_configuration,
        use_container_width=True,
        hide_index=True,
        on_select="rerun",
        selection_mode="single-row",
    )
    #st.write(event)
    if len(event['selection']['rows']) > 0:
        sel_variant_index_ = event['selection']['rows']
    else:
        sel_variant_index_ = 0

    r_sel_ = df_.loc[sel_variant_index_].squeeze()
    #st.write(r_sel_.variant_id)

    def parse_variant_id(variant_id):
        # df_var[['uniprot_id', 'aa_pos', 'aa_ref', 'aa_alt']] = df_var.apply(lambda r: parse_varstr(r['protein_variant']), axis=1, result_type='expand')
        aa_pos = int(variant_id[1:-1])
        aa_ref = variant_id[0]
        aa_alt = variant_id[-1]
        #print(uniprot_id, aa_pos, aa_ref, aa_alt)
        return aa_pos, aa_ref, aa_alt

    def parse_varstr(s):
        # df_var[['uniprot_id', 'aa_pos', 'aa_ref', 'aa_alt']] = df_var.apply(lambda r: parse_varstr(r['protein_variant']), axis=1, result_type='expand')
        uniprot_id, variant_id = s.split('/')
        aa_pos = int(variant_id[1:-1])
        aa_ref = variant_id[0]
        aa_alt = variant_id[-1]
        #print(uniprot_id, aa_pos, aa_ref, aa_alt)
        return uniprot_id, aa_pos, aa_ref, aa_alt


    uniprot_id, aa_pos, aa_ref, aa_alt = parse_varstr(r_sel_.variant_id)
    #st.write(uniprot_id)
    #st.write(aa_pos)
    #st.write(aa_ref)
    #st.write(aa_alt)

    col1, col2 = st.columns(2)

    @st.cache_resource
    def read_af2_v4(af2_id):
        url_ = f'https://alphafold.ebi.ac.uk/files/AF-{af2_id}-F1-model_v4.pdb'
        with urllib.request.urlopen(url_) as url:
            return url.read().decode('utf-8')

    with col1:
        st.write('### Structure with variant')
        #st.write(len(read_af2_v4(uniprot_id)))
        pdb_ = read_af2_v4(uniprot_id)
        xyzview = py3Dmol.view()
        #st.write('len2')
        # Add structure
        xyzview.addModel(pdb_, format='pdb')

        xyzview.setStyle({'model': 0}, {    
            'cartoon': {
                #'color': 'grey',
                #'colorscheme': {
                #    'prop': 'resi',
                #    #'map': colors_pocket,
                #},
                'arrows': True,
            }
        })

        xyzview.addStyle({'resi': aa_pos}, {'stick': {'colorscheme': 'purpleCarbon'}})
        xyzview.addLabel( #https://3dmol.org/doc/GLViewer.html#addLabel
            f'{aa_ref}{aa_pos}',
            { #https://3dmol.csb.pitt.edu/doc/LabelSpec.html
                'fontColor': 'red' if r_sel_.am_class == 'pathogenic' else 'blue' if r_sel_.am_class == 'benign' else 'gray',
                'backgroundColor': 'lightgray',
                #'fontOpacity': 0.5,
                'backgroundOpacity': 0.8,
            },
            {'resi': aa_pos},
        )
        xyzview.setBackgroundColor('#eeeeee')
        #xyzview.zoomTo({'resi': aa_pos}) # center exactly onto residue
        xyzview.zoomTo({ # center onto within 30 Angostroms
            'within': {
                'distance': 30,
                'sel': {'resi': aa_pos},
            },
        })
        stmol.showmol(xyzview, height=800, width=800)

    with col2:
        st.write('### Variant details')
        st.dataframe(r_sel_, use_container_width=True, height=800)
