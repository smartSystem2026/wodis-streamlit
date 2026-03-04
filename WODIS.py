# -*- coding: utf-8 -*-
"""
Created on 2026.03.04.
@author: S.M. Kim
"""

import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
from datetime import datetime
import io, sys
import math
from qpcr_app_Claude import run_qpcr_module

# ======================================================
# Session State Initialization
# ======================================================
if "df" not in st.session_state:
    st.session_state.df = None

if "analysis_result" not in st.session_state:
    st.session_state.analysis_result = None

if "p_val" not in st.session_state:
    st.session_state.p_val = None

if "pdf_buffer" not in st.session_state:
    st.session_state.pdf_buffer = None

if "uploaded_filename" not in st.session_state:
    st.session_state.uploaded_filename = None

# ======================================================
# Page Config
# ======================================================
st.set_page_config(page_title="Workflow Optimization System", layout="wide")


# ======================================================
# Sidebar
# ======================================================
with st.sidebar:
    st.title("Workflow Optimization and Data Interpretation System")
    st.divider()

    main_menu = st.selectbox(
        "Analysis Menu",
        ["1. Educational & SOP Support", "2. Data Analysis Functions"]
    )

    if main_menu == "1. Educational & SOP Support":
        sub_menu = st.radio(
            "Sub Menu",
            [
                "1. Regulatory and Institutional Knowledge Support",
                "2. SOP Comprehension and Learning",
                "3. Experimental Design Assistant",
                "4. Protocol Generation and Modification",
                "5. Knowledge Summarization and Literature Synthesis",
                "6. Hypothesis Generation and Experimental Planning"
            ]
        )
    else:
        sub_menu = st.radio(
            "Sub Menu",
            [
                "1. Histological and Imaging Data Interpretation",
                "2. qPCR and Gene Expression Data Analysis",
                "3. Genotyping and Genetic Data Correlation",
                "4. Behavioral Data Interpretation",
                "5. Omics and High-Dimensional Data Summarization",
                "6. Longitudinal and Time-Series Data Interpretation"
            ]
        )


# ======================================================
# Main Content
# ======================================================
st.header(sub_menu)


# ======================================================
# qPCR Analysis Module
# ======================================================
if sub_menu == "2. qPCR and Gene Expression Data Analysis":
	run_qpcr_module()   # ← Call External Module

else:
    st.write("The selected menu feature is currently under development.")
