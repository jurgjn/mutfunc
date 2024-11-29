
import ast, gzip, random, os, tempfile, time, sqlite3, urllib.request
import numpy as np, matplotlib, matplotlib.colors, matplotlib.pyplot as plt, seaborn as sns, pandas as pd, streamlit as st
import py3Dmol, stmol
import streamlit as st
import pandas as pd
import requests
from unipressed import IdMappingClient
from Bio.Data.IUPACData import protein_letters_3to1
import time

import util.variant_parser as vp
import util.db_utils as db

st.cache_resource.clear()

st.set_page_config(
    page_title='mutfunc - prototype',
    page_icon='🔬',
    layout='wide',
)

st.title("Mutfunc - let's get funky with mutations! 🧬")

st.markdown(
    """
    This is a prototype for the Mutfunc web application. The app is designed to help you analyze the effects of mutations.
    Go to the **Input Variants** page to input your variants, and then proceed to the **Browse Results** page to view the results.
""")

if st.button("Analyze some variants"):
    st.switch_page("pages/1_Input_Variants.py")
