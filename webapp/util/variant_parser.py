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
st.cache_data(persist='disk')
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

# fetching variant info based on genomic location
st.cache_data(persist='disk')
def fetch_variant_vep(chr, region, mutation, species="human"):
    """
    Fetch variant information from the Ensembl VEP Region endpoint.

    Parameters:
    - species: Species name (e.g., 'homo_sapiens').
    - region: Genomic region in the format 'chr:start-end' (e.g., '1:1000-2000').
    - allele: Allele change in the format 'reference/variant' (e.g., 'A/T').

    Returns:
    - JSON response from the API if successful, otherwise an error message.
    """
    server = "https://rest.ensembl.org"
    endpoint = f"/vep/{species}/region/{chr}:{region}-{region}/{mutation}?content-type=application/json&uniprot=1"
    headers = {"Content-Type": "application/json", "uniprot": "1"}

    url = f"{server}{endpoint}"
    
    try:
        response = requests.get(url, headers=headers)
        
        # Check if the request was successful
        if response.status_code != 200:
            return {"error": f"HTTP {response.status_code}: {response.text}"}

        return response.json()
    except requests.RequestException as e:
        return {"error": f"Request failed: {e}"}
    
def read_vep_result(result):
    # check if result is a list
    try: 
        result = result[0]
    except:
        return
    # check if result is a dict
    if type(result) != dict:
        return
    if result["most_severe_consequence"] != "missense_variant":
        return
    position = result["transcript_consequences"][0]["protein_start"]
    prot_id = result["transcript_consequences"][0]["uniprot_isoform"][0].split("-")[0]
    aa_from = result["transcript_consequences"][0]["amino_acids"][0]
    aa_to = result["transcript_consequences"][0]["amino_acids"][-1]
    # create string
    return f"{prot_id}/{aa_from}{position}{aa_to}"

def translate_genomic_coords(input):
    var_coords = [e.split(" ") for e in input]
    results = [fetch_variant_vep(var[0][3:], var[1], var[2][-1]) for var in var_coords]
    #st.write(results)
    prot_vars = [read_vep_result(result) for result in results]

    # create dict but only if values not none
    return {prot_vars[i]:input[i] for i in range(len(input)) if prot_vars[i] is not None}


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
    #st.write(ensemble_ids)

    ensemble_to_input = dict(zip(ensemble_ids, input_ids))
    ensemble_to_variant = dict(zip(ensemble_ids, variants))

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
    #st.write(result)

    #mapping = {e["from"]: e["to"] for e in result}

    # Get UniProt IDs and convert them to internal IDs with variants
    uniprot_id_map = list(set({(e["from"], e["to"]) for e in result}))#set(mapping.values())
    input_ids_left = [ensemble_to_input[e[0]] for e in uniprot_id_map]
    uniprot_ids = [e[1] for e in uniprot_id_map]
    #st.write(variants)
    filtered_variants = [ensemble_to_variant[e[0]] for e in uniprot_id_map]
    internal_ids = [convert_to_internal_id(uniprot_id, variant) for uniprot_id,variant in zip(uniprot_ids,filtered_variants)]

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
    Parse a list of variant lines into structured lists for rs IDs, genomic positions, and protein variants.

    Parameters:
    - variant_lines: List of strings, each representing a variant.

    Returns:
    - A dictionary with lists of variants:
        {
            "rs_variants": [...],
            "genomic_positions": [...],
            "protein_variants": [...],
            "error_ids": [...],
        }
    """
    rs_variants = []
    genomic_positions = []
    protein_variants = []
    error_ids = []

    for line in variant_lines:
        stripped_line = line.strip()
        if stripped_line:
            if stripped_line.startswith("rs"):
                # rs variant
                rs_variants.append(stripped_line)
            elif stripped_line.startswith("chr"):
                # Genomic position
                parts = stripped_line.split()
                if len(parts) == 3 and "/" in parts[2]:
                    genomic_positions.append(stripped_line)
                else:
                    error_ids.append(stripped_line)
            else:
                # Protein variation
                parts = stripped_line.split(" ")
                if len(parts) == 2 and parts[1][0].isalpha() and parts[1][-1].isalpha():
                    protein_variants.append(f"{parts[0]}/{parts[1]}")
                else:
                    error_ids.append(stripped_line)
        else:
            error_ids.append(stripped_line)

    return {
        "rs_variants": rs_variants,
        "genomic_positions": genomic_positions,
        "protein_variants": protein_variants,
        "error_ids": error_ids,
    }