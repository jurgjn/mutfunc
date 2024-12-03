
import streamlit as st
import pandas as pd
import util.variant_parser as vp
from collections import defaultdict


st.set_page_config(
    page_title='Input variants',
    page_icon='🔬',
    layout='wide',
)

st.title("Mutfunc - Let's get funky with mutations 🧬")

col1, col2 = st.columns(2)

prot_variants = {} #st.session_state.get("prot_variants", {}) # TODO maybe reset this
rs_variants = [] #st.session_state.get("rs_variants", [])
genomic_positions = [] #st.session_state.get("genomic_positions", [])
error_ids = [] #st.session_state.get("error_ids", [])

lookup_df = pd.DataFrame()


with col1:
    manual_input = st.text_area("Enter genomic variants (one per line)", value = "rs699\nrs6265\nP00533 R132C\nP09874 S568F\nP00451 G41C\nchr14 89993420 A/G")#, 'P09874/D678H', 'Q96NU1/R28Q', 'P00451/G41C', 'P01019/T259M'])")

with col2:
    uploaded_file = st.file_uploader("Upload a file containing genomic variants", type=["txt", "csv"])

if uploaded_file is not None:
    try:
        file_content = pd.read_csv(uploaded_file, header=None)
        parsed_variants = vp.parse_variants(file_content[0].tolist())
        rs = parsed_variants.get("rs_variants", [])
        genomic = parsed_variants.get("genomic_positions", [])
        prot = parsed_variants.get("protein_variants", [])
        rs_variants.extend(rs)
        genomic_positions.extend(genomic)
        prot_variants.update({p:p for p in prot})
        error_ids.extend(parsed_variants.get("error_ids", []))
        st.success("File processed successfully!")
    except Exception as e:
        st.error("Failed to read the uploaded file. Please ensure it contains variants in one column.")

if manual_input:
    parsed_variants = vp.parse_variants(manual_input.split("\n"))
    rs = parsed_variants.get("rs_variants", [])
    genomic = parsed_variants.get("genomic_positions", [])
    prot = parsed_variants.get("protein_variants", [])
    rs_variants.extend(rs)
    genomic_positions.extend(genomic)
    prot_variants.update({p:p for p in prot})
    error_ids.extend(parsed_variants.get("error_ids", []))
    st.success("Manual input processed successfully!") # TODO: only print this on input/change 

lookup_df["Input variant"] = rs_variants + genomic_positions + list(prot_variants.keys()) + error_ids
lookup_df["Inferred type"] = ["rsid"]*len(rs_variants) + ["genomic coord"]*len(genomic_positions) + ["protein var"]*len(prot_variants) + ["error"]*len(error_ids)

if st.button("Analyze Variants"):
    if not prot_variants and not rs_variants and not genomic_positions:
        st.warning("No valid variants provided. Please upload a file or enter variants manually.")
    else:
        # translate rs variants
        if rs_variants:
            try:
                with st.spinner("Fetching variant information from Ensembl..."):
                    variant_info = vp.fetch_variant_info(rs_variants)
                    if not variant_info:
                        st.error("No information was returned for the provided rs ids.")
                    else:
                        try:
                            translated_variants = vp.process_variant_info(variant_info)
                            prot_variants.update(translated_variants)
                            st.success("rs ids processed successfully!")
                        except Exception as e:
                            st.error(f"An error occurred when translating to uniprot ids: {e}")
            except Exception as e:
                st.error(f"An error occurred when fetching the rs id information: {e}")
        if genomic_positions:
            try:
                with st.spinner("Fetching coordinate information from Ensembl..."):
                    variant_info = vp.translate_genomic_coords(genomic_positions)
                    prot_variants.update(variant_info)
            except Exception as e:
                st.error(f"An error occurred when fetching the genomic positions: {e}")

        # Step 1: Create a dictionary with lists of all translations
        translations_dict = defaultdict(list)
        for translation, variant in prot_variants.items():
            translations_dict[variant].append(translation)

        lookup_df["Translated variant"] = lookup_df["Input variant"].map(translations_dict)
        lookup_df =  lookup_df.explode("Translated variant", ignore_index=True)
        st.session_state["lookup_df"] = lookup_df

        if len(prot_variants) == 0:
            st.warning("No valid variants provided. Please upload a file or enter variants manually.")
        else:
            st.success("All variants processed!")
            st.session_state["prot_variants"] = prot_variants
            st.write(prot_variants)
            #st.switch_page("pages/2_Browse_Results.py")


#st.session_state["rs_variants"] = rs_variants