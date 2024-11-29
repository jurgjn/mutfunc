import streamlit as st
import py3Dmol, stmol
import util.db_utils as db
import urllib.request

st.set_page_config(
    page_title='Browse Results',
    page_icon='🕵️‍♂️',
    layout='wide',
)

st.title("Mutfunc - Browse Results 🧬")

prot_variants = st.session_state.get("prot_variants", {})

#st.write(prot_variants)

@st.cache_resource
def read_af2_v4(af2_id):
    url_ = f'https://alphafold.ebi.ac.uk/files/AF-{af2_id}-F1-model_v4.pdb'
    with urllib.request.urlopen(url_) as url:
        return url.read().decode('utf-8')

if not prot_variants:
    st.info("No variants available for browsing. Please input variants first.")
else:
    internal_ids = list(prot_variants.keys())
    #st.write(prot_variants)
    df_ = db.query_missense(internal_ids)
    df_['input_variant'] = df_['variant_id'].apply(lambda x: prot_variants[x])
    col1, col2 = st.columns(2)
    with col1:
        event = st.dataframe(
            df_,
            column_order = ['input_variant','variant_id', 'am_pathogenicity', 'am_class', 'pred_ddg', 'pocketscore', 'interface'],
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

    with col2:
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
        st.write('### Variant details')
        st.dataframe(r_sel_, use_container_width=True, height=800)