import ast, gzip, random, os, tempfile, time, sqlite3, urllib.request
import numpy as np, matplotlib, matplotlib.colors, matplotlib.pyplot as plt, seaborn as sns, pandas as pd, streamlit as st
import py3Dmol, stmol
import streamlit as st
import pandas as pd
import requests
from unipressed import IdMappingClient
from Bio.Data.IUPACData import protein_letters_3to1
import time

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
    #st.subheader("Debugging: Raw API Response")
    #st.json(response.json())  # Display raw JSON response in the app

    if response.status_code != 200:
        st.error(f"Error {response.status_code}: {response.json().get('error', 'Unknown error')}")
        return []

    return response.json()


def query_prot_ids(ensemble_ids, input_ids):
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

    ensemble_to_input = dict(zip(ensemble_ids, input_ids))

    # Submit a mapping request to UniProt
    request = IdMappingClient.submit(source="Ensembl_Protein", dest="UniProtKB-Swiss-Prot", ids=set(ensemble_ids))

    max_retries = 10  # Maximum number of retries
    retry_count = 0
    with st.spinner("Looking up Uniprot ids..."):
        while retry_count < max_retries:
            try:
                result = list(request.each_result())
            except Exception as e:
                st.info(f"Failed to fetch UniProt IDs. Retrying... ({retry_count}/{max_retries})")
                retry_count += 1
                time.sleep(2)
                continue
            break

    # TODO: add try except and waiting logic if result takes time

    #mapping = {e["from"]: e["to"] for e in result}

    # Get UniProt IDs and convert them to internal IDs with variants
    uniprot_id_map = list(set({(e["from"], e["to"]) for e in result}))#set(mapping.values())
    input_ids_left = [ensemble_to_input[e[0]] for e in uniprot_id_map]
    uniprot_ids = [e[1] for e in uniprot_id_map]

    internal_ids = [convert_to_internal_id(uniprot_id, variants[0]) for uniprot_id in uniprot_ids]

    return dict(zip(internal_ids,input_ids_left)) # Return results if successful

def convert_to_internal_id(uniprot_id, variant):
    position = variant[3:-3]
    short_variant = protein_letters_3to1[variant[-3:]]+position+protein_letters_3to1[variant[:3]]
    return uniprot_id+"/"+short_variant

# Function to process the API response into a structured DataFrame
def process_variant_info(variant_info):
    # get the proteins 
    rows = []
    all_hgvsp = []
    input_vars = []
    for variant in variant_info:
        for nucleotide, details in variant.items():
            input_variant = details.get("input", "N/A")
            hgvsp = details.get("hgvsp", [])
            all_hgvsp.extend(hgvsp)
            input_vars.extend([input_variant]*len(hgvsp))

    protein_ids = query_prot_ids(all_hgvsp, input_vars)

    return protein_ids #pd.DataFrame(rows)

def parse_variants(variant_lines):
    """
    Parse a list of variant lines into structured lists for rs IDs and protein variants.

    Parameters:
    - variant_lines: List of strings, each representing a variant.

    Returns:
    - A tuple of two lists: rs_variants and prot_variants.
    """
    rs_variants = []
    prot_variants = []
    for line in variant_lines:
        stripped_line = line.strip()
        if stripped_line:
            # Check if it's an rs ID or a protein ID
            if stripped_line.startswith("rs"):
                rs_variants.append(stripped_line)
            else:
                # Split the line into protein ID and mutation
                parts = stripped_line.split()
                if len(parts) == 2:
                    protein_id, mutation = parts
                    #prot_variants.append({"protein": protein_id, "mutation": mutation})
                    prot_variants.append(f"{protein_id}/{mutation}")
                else:
                    st.warning(f"Invalid format in line: {stripped_line}. Expected 'proteinID mutation'.")
    return rs_variants, prot_variants