
import streamlit as st
import pandas as pd
import util.db_utils as db
import util.variant_parser as vp
import os
from collections import defaultdict
from util.footer import show_footer

st.set_page_config(page_title='Mutfunc: Input', page_icon='🧬', layout='wide',)

st.title("Mutfunc - precomputed mechanistic consequences of mutations")

col1, col2 = st.columns(2)

prot_variants = {} #st.session_state.get("prot_variants", {}) # TODO maybe reset this
rs_variants = [] #st.session_state.get("rs_variants", [])
genomic_positions = [] #st.session_state.get("genomic_positions", [])
error_ids = [] #st.session_state.get("error_ids", [])

lookup_df = pd.DataFrame()

with col1:
    st.write("Enter your variants of interest below or upload a file.  \n  We currently only support human variants. Please use rsids (e.g. rs699); protein variants (e.g. P00533 R132C); or genomic coordinates (e.g. chr14 89993420 A/G).")
    examples = {
        '(User-defined)': '',
        'Mixed genomic/proteomic variants': "rs699\nrs6265\nP00533 R132C\nP09874 S568F\nP00451 G41C\nchr14 89993420 A/G",# 'P09874/D678H', 'Q96NU1/R28Q', 'P00451/G41C', 'P01019/T259M']),
        'SLC25A4/P12235 (Fig. 3e)': 'P12235',
        'MAPK1/P28482 (Fig. 6c)': 'P28482',
    }
    examples_sel = st.selectbox(label='Choose example:', options=examples.keys())
    value_ = examples[examples_sel]
    manual_input = st.text_area(label='Enter variants (one per line):', value=value_)

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

if st.button("Parse Variants"):
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

    # Try to parse input as a single uniprot_id (and fetch all clinvar variants as input)
    if len(manual_input.split("\n")) == 1 and len(prot_variants) == 0:
        if len(manual_input) == 0:
            st.warning("Input looks empty. Please upload a file or enter variants manually.")
            st.stop()

        uniprot_id_ = manual_input.split()[0]
        #st.write(len(db.read_clinvar()))
        #st.write(uniprot_id_)
        #st.write(db.read_clinvar().query('uniprot_id == @uniprot_id_'))
        prot_variants_ = db.read_clinvar().query('uniprot_id == @uniprot_id_')['variant_id'].tolist()
        prot_variants = { var_: var_ for var_ in prot_variants_ }
        error_ids = []

    if not prot_variants and not rs_variants and not genomic_positions:
        st.warning("No valid variants provided. Please upload a file or enter variants manually.")
    else:
        lookup_df["Input variant"] = rs_variants + genomic_positions + list(prot_variants.keys()) + error_ids
        lookup_df["Inferred type"] = ["rsid"]*len(rs_variants) + ["genomic coord"]*len(genomic_positions) + ["protein var"]*len(prot_variants) + ["error"]*len(error_ids)

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
        lookup_df = lookup_df.explode("Translated variant", ignore_index=True)
        st.session_state["lookup_df"] = lookup_df

        in_mutfunc = set(db.query_missense(lookup_df['Translated variant'].dropna())['variant_id'])
        lookup_df['isin_mutfunc'] = lookup_df['Translated variant'].map(lambda var: var in in_mutfunc)

        if len(prot_variants) == 0:
            st.warning("No valid variants provided. Please upload a file or enter variants manually.")
        else:
            st.success("All variants processed!")
            st.session_state["prot_variants"] = prot_variants

            # Show table with coordinate translation outcomes
            st.dataframe(lookup_df, use_container_width=True, hide_index=True)

            # Calculate metrics
            found = lookup_df.query('isin_mutfunc')["Input variant"].nunique()
            total = lookup_df["Input variant"].nunique()
            ratio = found / total

            # Determine emoji based on the ratio
            if ratio >= 1.0:
                emoji = "🎉"  # Great success
            elif 0.5 <= ratio:
                emoji = "🧐"  # Okay success
            else:
                emoji = "😰"  # Poor success

            # Display message with emoji
            st.write(
                "{emoji} **{found} out of {total}** variants were successfully translated and found in the database.".format(
                    emoji=emoji, found=found, total=total
                )
            )

            # Button/link to browsing the results
            page_ = 'pages/2_Browse_Results.py'
            #if st.button('Browse Results'):
            #    st.switch_page(page_)
            st.page_link(page_, label='Browse Results', icon='⏩')

show_footer()