from pprint import pprint
import dash
from dash.exceptions import PreventUpdate
import requests
import util.variant_parser

import pandas as pd
import pyarrow.parquet as pq
from dash import Dash, html, Input, Output, callback, dcc
import dash_ag_grid as dag

import dash_molstar

import dash_bootstrap_components as dbc

from dash_molstar.utils import molstar_helper
from dash_molstar.utils.representations import Representation
import re
import foldcomp

dash.register_page(__name__, path="/variants", title="Browse variants")
#df = pd.read_parquet('data/human-missense.parquet').head(10)

missense_cols_ = [
    'variant_id',
    #'ESM1b_LLR',
    #'ESM1b_is_pathogenic',
    'am_pathogenicity',
    'am_class',
    'pred_ddg',
    #'plddt',
]


def g_convert(query=["CASQ2", "CASQ1", "GSTO1", "DMD", "GSTM2"], target='UNIPROTSWISSPROT_ACC', organism='hsapiens', numeric_namespace='ENTREZGENE_ACC'):
    """
    Query for HGNC gene names using g:convert (https://biit.cs.ut.ee/gprofiler/convert)
    """
    r = requests.post(url='https://biit.cs.ut.ee/gprofiler/api/convert/convert/', json=locals())
    df_ = pd.DataFrame(r.json()['result'])
    return df_

variants_table_ = dag.AgGrid(
    id="variants",
    style={
        "height": '500px',
        "width": "100%",
        "--ag-background-color": "#1a2232",           # --mi-surface
        #"--ag-odd-row-background-color": "#1d2840",   # between surface and surface-2 (stripe effect)
        "--ag-header-background-color": "#1a2232",    # --mi-surface (matches thead)
        "--ag-foreground-color": "#e7ecf5",           # --mi-text
        "--ag-border-color": "#344465",               # --mi-border-strong
        #"--ag-row-hover-color": "#243047",            # --mi-surface-2 (solid proxy for the rgba hover)
        #"--ag-selected-row-background-color": "rgba(25, 93, 230, 0.16)",  # --mi-primary @ 16% (matches .table active)
        "--ag-input-focus-border-color": "#195de6",   # --mi-primary
        "--ag-checkbox-checked-color": "#195de6",     # --mi-primary
        "--ag-header-foreground-color": "#ffffff",    # --mi-text-strong (matches thead th)
        "--ag-secondary-foreground-color": "#93a5c8", # --mi-muted (matches tbody td)
    },
    rowData=[],#df.to_dict('records'),
    columnDefs=[
        #{
        #    "checkboxSelection": True,
        #    "width": 50,
        #},
        {
            "field": "input",
            "headerName": "Input",
            "headerClass": "multiline-header",
            "wrapHeaderText": True,
            "autoHeaderHeight": True,
        },
        {
            "field": "variant_id",
            "headerName": "Missense\nvariant",
            "headerClass": "multiline-header",
            "wrapHeaderText": True,
            "autoHeaderHeight": True,
        },
        {
            "field": "am_pathogenicity",
            "headerName": "Pathogenicity\n(AMS)",
            "headerClass": "multiline-header",
            "wrapHeaderText": True,
            "autoHeaderHeight": True,
            "valueFormatter": {"function": "Number(params.value).toFixed(2)"},
        },
        {
            "field": "pred_ddg",
            "headerName": "Stability\n(ddG)",
            "headerClass": "multiline-header",
            "wrapHeaderText": True,
            "autoHeaderHeight": True,
            "valueFormatter": {"function": "Number(params.value).toFixed(2)"},
        },
        {
            "field": "pocket_score_str",
            "headerName": "Pocket\nscore",
            "headerClass": "multiline-header",
            "wrapHeaderText": True,
            "autoHeaderHeight": True,
        },
        #*[{"field": i, "headerName": i} for i in pq.read_schema('data/missense.parquet').names],
        #*[{"field": i, "headerName": i} for i in ['input'] + missense_cols_ ],
    ],
    columnSize="sizeToFit",
    dashGridOptions={"rowSelection": {"mode": "singleRow"}},
    #dashGridOptions={"rowSelection": "multiple"},
)

structure_viewer_ = dash_molstar.MolstarViewer(
    id='viewer',
    style={
        'width': 'auto',
        'height':'500px',
    },
    layout={
        'backgroundColor': 0x000000,
    }
)

layout = dbc.Container([
    html.Div(id="print-variant-count"),
    html.Div([
        # Left column
        html.Div([
            dcc.Loading(
                id="loading_variants",
                type="default",
                children=dbc.Container([variants_table_], fluid=True, className="dbc dbc-ag-grid"),
            )
        ], style={"flex": "1"}),
    
        # Right column
        html.Div([
            dcc.Loading(
                id="loading_structure",
                type="default",
                children=dbc.Container([structure_viewer_], fluid=True),
            ),
        ], style={"flex": "1"}),

    ], style={"display": "flex", "gap": "10px"}),
], style={"padding": "20px"}, fluid=True)

@callback(
    Output("print-variant-count", "children"),
    #Input("variant-list", "data"),
    Input("variants",  "rowData"),
)
def print_variant_list(data):
    if not data:
        return "variants"
    return f"Showing results for {len(data)} mapped variants"

@callback(
    Output("variants",  "rowData"),
    Output("variants", "selectedRows"),
    Input("variant-list", "data"),
)
def read_variants(variant_list):
    """Propagate variant-list store data into the AgGrid on page load."""

    pattern = r'^.+/[A-Z]\d+[A-Z]$'

    mapping = pd.DataFrame()
    mapping['input'] = variant_list
    mapping = mapping.query('input != ""')
    mapping['vep'] = mapping['input'].map(lambda input: False if re.fullmatch(pattern, input) else True)
    pprint(mapping)

    if len(mapping.query('vep')) > 0:
        query_vep = util.variant_parser.fetch_variants(mapping.query('vep')['input'])['missense']
        mapping_vep = pd.DataFrame.from_records([(r['input'], r['gene_id'], str(r['protein_position']), r['amino_acids'] ) for r in query_vep], columns=['input', 'gene_id', 'protein_position', 'amino_acids']).drop_duplicates(keep='first')
        mapping_vep['aa_ref'] = mapping_vep['amino_acids'].map(lambda s: s.split('/')[0])
        mapping_vep['aa_alt'] = mapping_vep['amino_acids'].map(lambda s: s.split('/')[1])
        mapping_vep = pd.concat([mapping_vep.reset_index(drop=True), g_convert(mapping_vep['gene_id'].tolist(), target='UNIPROTSWISSPROT_ACC').reset_index(drop=True)], axis=1)
        mapping_vep['variant_id'] = mapping_vep['converted'] + '/' + mapping_vep['aa_ref'] + mapping_vep['protein_position'] + mapping_vep['aa_alt']
        mapping_vep = mapping_vep[['input', 'variant_id']]
    else:
        mapping_vep = pd.DataFrame.from_records([], columns=['input', 'variant_id'])

    mapping_merge = mapping.merge(mapping_vep, on='input', how='left')
    mapping_merge['variant_id'] = list(map(lambda input, vep, variant_id: variant_id if vep else input, mapping_merge['input'], mapping_merge['vep'], mapping_merge['variant_id']))

    if not variant_list:
        return [], "No variants loaded."

    df = pd.read_parquet(
        'data/missense.parquet',
        filters=[('variant_id', 'in', mapping_merge['variant_id'])],
    )#[missense_cols_]

    df = mapping_merge[['input', 'variant_id']].merge(df, on='variant_id')
    df[['uniprot_id', 'pos', 'ref', 'alt']] = df.apply(lambda r: parse_varstr(r['variant_id']), axis=1, result_type='expand')

    list_uniprot_id = set(df['variant_id'].map(lambda variant_id: variant_id.split('/')[0]).tolist())
    pockets = pd.read_parquet('data/pockets.parquet', filters=[('uniprot_id', 'in', list_uniprot_id)],)
    print('pockets')
    pprint(pockets)

    pocket_resid = pockets.explode('pocket_resid')[['uniprot_id', 'pocket_resid', 'pocket_id', 'pocket_score_combined_scaled']]\
        .sort_values('pocket_score_combined_scaled', ascending=False)\
        .reset_index(drop=True)\
        .groupby(['uniprot_id', 'pocket_resid'])\
        .head(1)

    df = df.merge(pocket_resid, left_on=['uniprot_id', 'pos'], right_on=['uniprot_id', 'pocket_resid'], how='left')
    df['pocket_score_str'] = df['pocket_score_combined_scaled'].apply(lambda x: f'{x:.2g}')
    df['pocket_score_str'] = df['pocket_score_str'].replace('nan', '')

    return df.to_dict('records'), df.head(1).to_dict('records')
    #return df.to_dict('records'), df.to_dict('records')

def parse_varstr(s):
    # df_var[['uniprot_id', 'aa_pos', 'aa_ref', 'aa_alt']] = df_var.apply(lambda r: parse_varstr(r['protein_variant']), axis=1, result_type='expand')
    uniprot_id, variant_id = s.split('/')
    aa_pos = int(variant_id[1:-1])
    aa_ref = variant_id[0]
    aa_alt = variant_id[-1]
    #print(uniprot_id, aa_pos, aa_ref, aa_alt)
    return uniprot_id, aa_pos, aa_ref, aa_alt

def build_bb_sc(pos_pathogenic, pos_benign):
    bb_rep = Representation(type='cartoon', color='uniform')
    bb_rep.set_type_params({'visuals': ['polymer-trace', 'polymer-gap']})
    bb_rep.set_color_params({'value': 0x888888})

    bb_targets = molstar_helper.get_targets(chain='A')
    bb_comp = molstar_helper.create_component(label='Backbone', targets=bb_targets, representation=[bb_rep])

    ps_rep = Representation(type='ball-and-stick', color='uniform')
    bs_rep = Representation(type='ball-and-stick', color='uniform')
    ps_rep.set_color_params({'value': 0xFF0000})
    bs_rep.set_color_params({'value': 0x0000FF})

    ps_sel  = molstar_helper.get_targets(chain='A', residue=pos_pathogenic)
    bs_sel  = molstar_helper.get_targets(chain='A', residue=pos_benign)

    ps_comp = molstar_helper.create_component(label='p_sidechains', targets=ps_sel, representation=[ps_rep])
    bs_comp = molstar_helper.create_component(label='b_sidechains', targets=bs_sel, representation=[bs_rep])

    component = [ bb_comp ]
    if len(pos_pathogenic) > 0:
        component.append(ps_comp)
    if len(pos_benign) > 0:
        component.append(bs_comp)

    return component
    #return [bb_comp, ps_comp, bs_comp]

@callback(
    Output('viewer', 'data'),
    Input("variants", "selectedRows"),
    #Input("variants", "rowData"),
    prevent_initial_call=True,
)
def show_structure(selected_rows):#, variants_data):
    if not selected_rows:
        raise PreventUpdate

    uniprot_id, pos, ref, alt = parse_varstr(selected_rows[0]['variant_id'])

    df_ = pd.DataFrame(selected_rows)
    df_[['uniprot_id', 'pos', 'ref', 'alt']] = df_.apply(lambda r: parse_varstr(r['variant_id']), axis=1, result_type='expand')
    df_ = df_.set_index('variant_id')

    p_pos = df_.query('am_class == "pathogenic"')['pos'].tolist()
    b_pos = df_.query('am_class == "benign"')['pos'].tolist()

    with foldcomp.open('data/monomers-db', ids=[f'{uniprot_id}-F1']) as db:
        for pdb in db:
            pdb_str = pdb[1]

    structure = molstar_helper.parse_molecule(
        pdb_str,
        fmt='pdb',
        component=build_bb_sc(p_pos, b_pos),
        preset={'kind': 'empty'}
    )

    molstar_data = [ structure ]

    # If selected variant has a pocket, read & show pocket surface
    uniprot_id = selected_rows[0]['uniprot_id']
    pocket_id = selected_rows[0]['pocket_id']
    print('attempt to load', uniprot_id, pocket_id)
    if not(pocket_id is None):
        print('new')
        print(uniprot_id)
        print(pocket_id)

        pockets = pd.read_parquet('data/pockets.parquet', filters=[('uniprot_id', '==', uniprot_id), ('pocket_id', '==', pocket_id)],)
        pocket_str = pockets.head(1)['pocket_cl_str'].squeeze()

        rep = Representation(type='molecular-surface')
        rep.color = 'uniform'
        rep.set_color_params({'value': 0xFFB6C1})
        rep.set_type_params({'alpha': 0.5})
    
        pocket_component = molstar_helper.create_component(
            label="Pocket",
            targets=molstar_helper.get_targets(chain=" "),
            representation=rep,
        )

        pocket = molstar_helper.parse_molecule(
            pocket_str,
            fmt='pdb',
            component=pocket_component,
            preset={'kind': 'empty'},
        )

        molstar_data.append(pocket)

    return molstar_data

    #return molstar_helper.parse_url(
    #    f'https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v6.pdb',
    #    component=build_bb_sc(p_pos, b_pos),
    #    preset={'kind': 'empty'}
    #)

#@app.callback(Output('viewer', 'data'),
#              Input('load_protein', 'n_clicks'),
#              prevent_initial_call=True)
#def display_output(yes):
#    data = molstar_helper.parse_url('https://files.rcsb.org/download/3U7Y.cif')
#    return data

#if __name__ == '__main__':
#    app.run(debug=True)

'''
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
        tab1, tab2 = st.tabs(["Variants", "Pockets"])
        with tab1:
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

        with tab2:
            st.write('Select a row to show/hide pockets ☑️')
            df_pocket = db.query_pockets(uniprot_id)
            sel_pocket = st.dataframe(
                df_pocket,
                column_order = ['uniprot_id', 'pocket_id', 'pocket_score', 'pocket_resid'],
                use_container_width=True,
                hide_index=True,
                on_select='rerun',
                selection_mode='multi-row',
            )

    with col2:
        # Add (monomer) structure
        pdb_ = read_af2_v4(uniprot_id)
        xyzview = py3Dmol.view()
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

        # Add variants to structure
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

        # Add user-selected pockets to structure viewer
        for i, r in df_pocket.loc[sel_pocket['selection']['rows']].iterrows():
            xyzview.addModel(r.pocket_cl, format='pdb')
            xyzview.setStyle({'model': -1}, {})
            xyzview.addSurface(py3Dmol.VDW, {'opacity': 0.7, 'color': 'pink'}, {'model': -1})

        # Set point-of-view/background 
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
'''