import streamlit as st
import py3Dmol, stmol
import util.db_utils as db
import urllib.request
import pandas as pd
from util.footer import show_footer

st.set_page_config(page_title='Mutfunc: Results', page_icon='🧬', layout='wide',)

page_ = '/app/Mutfunc_Home.py'
st.page_link(page_, label='Input Coordinates', icon='⏪')

st.title("Mutfunc - Browse Results 🕵️‍♂️")

#prot_variants = st.session_state.get("prot_variants", {}) # Not needed?
lookup_df = st.session_state.get("lookup_df", pd.DataFrame())

@st.cache_data
def read_af2_v4(af2_id):
    url_ = f'https://alphafold.ebi.ac.uk/files/AF-{af2_id}-F1-model_v4.pdb'
    with urllib.request.urlopen(url_) as url:
        return url.read().decode('utf-8')

if not(len(lookup_df.query('`Inferred type` != "error"')) > 0):
    st.info("No variants available for browsing. Please input variants first.")
else:
    df_ = db.query_missense(lookup_df['Translated variant'])
    df_ = lookup_df.merge(df_, left_on='Translated variant', right_on='variant_id')
    #st.write(df_)

    col1, col2 = st.columns(2)
    with col1:
        st.write('Select a row to display the variant in the structure ☑️')
        event = st.dataframe(
            df_,
            column_order = ['Input variant','variant_id', 'am_pathogenicity', 'am_class', 'pred_ddg', 'pocketscore', 'interface'],
            #column_config=column_configuration,
            use_container_width=True,
            hide_index=True,
            on_select="rerun",
            selection_mode="single-row",
        )
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

    with col2:
        #st.write('### Structure with variant')
        #st.write(len(read_af2_v4(uniprot_id)))
        pdb_ = read_af2_v4(uniprot_id)
        xyzview = py3Dmol.view()
        #st.write('len2')
        # Add structure
        xyzview.addModel(pdb_, format='pdb')

        xyzview.setStyle({'model': 0}, {    
            'cartoon': {
                'color': 'lightgrey',#'white',
                #'colorscheme': {
                #    'prop': 'resi',
                #    #'map': colors_pocket,
                #},
                'arrows': True,
            }
        })

        xyzview.addStyle({'resi': aa_pos}, {'stick': {'colorscheme': 'purpleCarbon'}})
        xyzview.addLabel( #https://3dmol.org/doc/GLViewer.html#addLabel
            f'{aa_ref}{aa_pos}{aa_alt}',
            { #https://3dmol.csb.pitt.edu/doc/LabelSpec.html
                'fontColor': '#ba181b' if r_sel_.am_class == 'pathogenic' else '#0077b6' if r_sel_.am_class == 'benign' else 'gray',
                'backgroundColor': 'lightgray',
                #'fontOpacity': 0.5,
                'backgroundOpacity': 0.7,
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
        stmol.showmol(xyzview, width=800)
        st.write('### Variant details')
        st.dataframe(r_sel_, use_container_width=True, height=800)

# footer 
show_footer()