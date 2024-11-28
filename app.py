
import ast, gzip, random, os, tempfile, time, sqlite3, urllib.request
import numpy as np, matplotlib, matplotlib.colors, matplotlib.pyplot as plt, seaborn as sns, pandas as pd, streamlit as st
import py3Dmol, stmol
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

st.write('## Variants')
df_ = query_missense(['P09874/S568F', 'P09874/D678H', 'Q96NU1/R28Q', 'P00451/G41C'])
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


col1, col2 = st.columns(2)
st.write('## Selected variant')
uniprot_id, aa_pos, aa_ref, aa_alt = parse_varstr(r_sel_.variant_id)
#st.write(uniprot_id)
#st.write(aa_pos)
#st.write(aa_ref)
#st.write(aa_alt)

@st.cache_resource
def read_af2_v4(af2_id):
    url_ = f'https://alphafold.ebi.ac.uk/files/AF-{af2_id}-F1-model_v4.pdb'
    with urllib.request.urlopen(url_) as url:
        return url.read().decode('utf-8')

with col1:
    st.write('### Structure')
    #st.write(len(read_af2_v4(uniprot_id)))
    pdb_ = read_af2_v4(uniprot_id)
    xyzview = py3Dmol.view()
    #st.write('len2')
    # Add structure
    xyzview.addModel(pdb_, format='pdb')

    xyzview.setStyle({'model': 0}, {
        'cartoon': {
            'colorscheme': {
                'prop': 'resi',
                #'map': colors_pocket,
            },
            'arrows': True,
        }
    })

    xyzview.setBackgroundColor('#eeeeee')
    xyzview.zoomTo()
    stmol.showmol(xyzview, height=600, width=600)

with col2:
    st.write('### Variant details')
    st.dataframe(r_sel_, use_container_width=True, height=800)
