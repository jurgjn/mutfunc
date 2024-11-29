import streamlit as st
import pandas as pd
import util.variant_parser as vp

st.set_page_config(
    page_title='Input variants',
    page_icon='🔬',
    layout='wide',
)

st.title("Mutfunc - Input Variants 🧬")

col1, col2 = st.columns(2)

prot_variants = st.session_state.get("prot_variants", [])
rs_variants = st.session_state.get("rs_variants", [])

with col1:
    manual_input = st.text_area("Enter genomic variants (one per line)", value = "rs699\nP00533 R132C\nP09874 S568F\nP00451 G41C")#, 'P09874/D678H', 'Q96NU1/R28Q', 'P00451/G41C', 'P01019/T259M'])")

with col2:
    uploaded_file = st.file_uploader("Upload a file containing genomic variants", type=["txt", "csv"])

if uploaded_file is not None:
    try:
        file_content = pd.read_csv(uploaded_file, header=None)
        rs, prot = vp.parse_variants(file_content[0].tolist())
        rs_variants.extend(rs)
        prot_variants.extend(prot)
        st.success("File processed successfully!")
    except Exception as e:
        st.error("Failed to read the uploaded file. Please ensure it contains variants in one column.")

if manual_input:
    rs, prot = vp.parse_variants(manual_input.split("\n"))
    rs_variants.extend(rs)
    prot_variants.extend(prot)
    st.success("Manual input processed successfully!")

if st.button("Analyze Variants"):
    if not prot_variants and not rs_variants:
        st.warning("No variants provided. Please upload a file or enter variants manually.")
    else:
        try:
            with st.spinner("Fetching variant information from Ensembl..."):
                variant_info = vp.fetch_variant_info(rs_variants)
                if not variant_info:
                    st.error("No information was returned for the provided variants.")
                else:
                    results_df = vp.process_variant_info(variant_info)
                    prot_variants.extend(results_df)
                    st.session_state["prot_variants"] = prot_variants
                    st.session_state["rs_variants"] = rs_variants
                    st.success("Variants processed successfully!")
                    st.session_state["prot_variants"] = prot_variants
                    st.session_state["rs_variants"] = rs_variants
                    st.switch_page("pages/2_Browse_Results.py")
        except Exception as e:
            st.error(f"An error occurred when processing the variants: {e}")

st.session_state["prot_variants"] = prot_variants
st.session_state["rs_variants"] = rs_variants