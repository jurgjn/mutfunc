import pandas as pd
import sqlite3
import numpy as np
import streamlit as st

@st.cache_resource
def read_clinvar():
    df_ = pd.read_csv('/data/clinvar_snv_vep.missense.tsv', sep='\t', dtype={'CHROM': str})
    return df_

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

def query_missense(variants):
    if len(variants) == 0:
        return pd.DataFrame()
    elif len(variants) == 1:
        variants_tuple = f"('{variants[0]}')"  # Single element tuple as string
    else:
        variants_tuple = str(tuple(variants))  # Normal tuple for multiple elements
    #https://stackoverflow.com/questions/28735213/pandas-read-sql-with-a-list-of-values-for-where-condition
    with sqlite3.connect('/data/missense_23.11.1.sqlite') as db:
        df_ = pd.read_sql_query(sql='SELECT * FROM missense WHERE variant_id in ' + variants_tuple, con=db)
    return df_.replace({
        'pred_ddg': -99,
        'pocketscore': -99,
        'pocketrank': -99,
        'interface': -99,
        'interface_strict': -99,
        'freq': -99,
    }, np.nan)