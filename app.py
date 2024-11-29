
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

import util.variant_parser as vp
import util.db_utils as db

st.cache_resource.clear()

st.set_page_config(
    page_title='mutfunc - prototype',
    page_icon='🔬',
    layout='wide',
)

# Streamlit app
st.title("Mutfunc - let's get funky with mutations! 🧬")

# Sidebar for input
st.sidebar.header("Input Variants")
uploaded_file = st.sidebar.file_uploader("Upload a file containing genomic variants", type=["txt", "csv"])
manual_input = st.text_area("Enter genomic variants (one per line)")

prot_variants = []  # List for protein variants
rs_variants = []    # List for rs IDs

# Handle uploaded file
if uploaded_file is not None:
    try:
        # Assuming the file is a plain text file or CSV with one column of variants
        file_content = pd.read_csv(uploaded_file, header=None)
        rs, prot = vp.parse_variants(file_content[0].tolist())
        rs_variants.extend(rs)
        prot_variants.extend(prot)
    except Exception as e:
        st.sidebar.error("Failed to read the uploaded file. Please ensure it contains variants in one column.")

# Handle manual input
if manual_input:
    rs, prot = vp.parse_variants(manual_input.split("\n"))
    rs_variants.extend(rs)
    prot_variants.extend(prot)

if st.sidebar.button("Analyze Variants"):
    if not prot_variants and not rs_variants:
        st.warning("No variants provided. Please upload a file or enter variants manually.")
    else:
        if rs_variants:
            try:
                with st.spinner("Fetching variant information from Ensembl..."):
                    # Fetch variant information from the Variant Recoder API
                    variant_info = vp.fetch_variant_info(rs_variants)
                    if not variant_info:
                        st.error("No information was returned for the provided variants.")
                    else:
                        # Process the JSON response into a DataFrame
                        results_df = vp.process_variant_info(variant_info)
                        prot_variants.extend(results_df)
                        #st.subheader("Processed Variant Results")
                        #st.dataframe(results_df)
                        #st.download_button("Download Results as CSV", results_df.to_csv(index=False), "variant_results.csv", "text/csv")
            except Exception as e:
                st.error(f"An error occurred when transforming the variants: {e}")

    st.write('## Variants')
    #df_ = query_missense(['P09874/S568F', 'P01019/T259M'])#, 'P09874/D678H', 'Q96NU1/R28Q', 'P00451/G41C', 'P01019/T259M'])
    st.write(prot_variants)
    df_ = db.query_missense(prot_variants)
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

    uniprot_id, aa_pos, aa_ref, aa_alt = db.parse_varstr(r_sel_.variant_id)
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
