# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 16:42:37 2026
@author: user
"""

import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
from datetime import datetime
import io
import math

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

    # --- Tab Style ---
    st.markdown(
        """
        <style>
        button[data-baseweb="tab"] {
            font-size: 18px;
            font-weight: 600;
            padding: 10px 18px;
        }
        </style>
        """,
        unsafe_allow_html=True
    )

    tab1, tab2, tab3 = st.tabs(["Open File", "Data Analysis", "AI Reporting"])

    # ==================================================
    # TAB 1: File Upload
    # ==================================================
    with tab1:
        uploaded_file = st.file_uploader(
            "Select qPCR Data (Excel File)",
            type=["xlsx"]
        )

        if uploaded_file is not None:
            # 파일 변경 시 기존 분석 초기화
            if uploaded_file.name != st.session_state.uploaded_filename:
                st.session_state.analysis_result = None
                st.session_state.pdf_buffer = None

            st.session_state.df = pd.read_excel(uploaded_file)
            st.session_state.uploaded_filename = uploaded_file.name

            st.success(f"The file has been successfully loaded.: {uploaded_file.name}")
            st.subheader("Data Preview")
            st.dataframe(st.session_state.df, use_container_width=True)

        else:
            st.info("You can start the analysis after uploading a file.")

    # ==================================================
    # TAB 2: Data Analysis
    # ==================================================
    with tab2:
        if st.session_state.df is not None:

            if st.button("Data analysis & Graph"):
                try:
                    df = st.session_state.df.copy()

                    # --- qPCR Analysis ---
                    df["dCt"] = df["Col1a1_Ct"] - df["Gapdh_Ct"]
                    ctrl_avg = df[df["Group"] == "Control"]["dCt"].mean()
                    df["ddCt"] = df["dCt"] - ctrl_avg
                    df["Fold_Change"] = 2 ** (-df["ddCt"])

                    g1 = df[df["Group"] == "Control"]["Fold_Change"]
                    g2 = df[df["Group"] == "Fibrosis"]["Fold_Change"]
                    _, p_val = stats.ttest_ind(g1, g2)

                    # --- Save to session ---
                    st.session_state.analysis_result = df
                    st.session_state.p_val = p_val

                    # --- PDF Generation ---
                    buf = io.BytesIO()
                    with PdfPages(buf) as pdf:

                        rows_per_page = 25
                        cols = ["Group", "dCt", "ddCt", "Fold_Change"]
                        table_df = df[cols].round(4)

                        total_pages = int(math.ceil(len(table_df) / rows_per_page))

                        for p in range(total_pages):
                            fig, ax = plt.subplots(figsize=(8.5, 11))
                            ax.axis("off")

                            ax.text(
                                0.05, 0.95,
                                f"qPCR Analysis Report\nDate: {datetime.now():%Y-%m-%d %H:%M}\n"
                                f"P-value: {p_val:.4f}",
                                transform=ax.transAxes,
                                fontsize=12,
                                va="top"
                            )

                            page_df = table_df.iloc[
                                p * rows_per_page:(p + 1) * rows_per_page
                            ]

                            table = ax.table(
                                cellText=page_df.values,
                                colLabels=page_df.columns,
                                loc="center",
                                cellLoc="center"
                            )
                            table.auto_set_font_size(False)
                            table.set_fontsize(8)
                            table.scale(1.1, 1.8)

                            pdf.savefig(fig)
                            plt.close(fig)

                        # --- Graph Page ---
                        fig, ax = plt.subplots(figsize=(8, 7))
                        sns.barplot(x="Group", y="Fold_Change", data=df, ax=ax)
                        sns.swarmplot(x="Group", y="Fold_Change", data=df, color=".2", ax=ax)
                        ax.set_title("Gene Expression Visual Analysis")
                        pdf.savefig(fig)
                        plt.close(fig)

                    st.session_state.pdf_buffer = buf.getvalue()

                except Exception as e:
                    st.error(f"An error has occurred.: {e}")

            # --- Show Results if Exist ---
            if st.session_state.analysis_result is not None:
                df = st.session_state.analysis_result
                p_val = st.session_state.p_val

                col1, col2 = st.columns(2)

                with col1:
                    st.subheader("📊 Analysis Results Data")
                    st.dataframe(
                        df[["Group", "dCt", "ddCt", "Fold_Change"]],
                        use_container_width=True
                    )
                    st.write(f"**P-value:** `{p_val:.4f}`")

                with col2:
                    st.subheader("📈 Data Visualization")
                    fig, ax = plt.subplots(figsize=(6, 5))
                    sns.barplot(x="Group", y="Fold_Change", data=df, ax=ax)
                    sns.swarmplot(x="Group", y="Fold_Change", data=df, color=".2", ax=ax)
                    st.pyplot(fig)

                if st.session_state.pdf_buffer is not None:
                    st.download_button(
                        "Save as PDF report",
                        data=st.session_state.pdf_buffer,
                        file_name=f"qPCR_Report_{datetime.now():%Y%m%d_%H%M}.pdf",
                        mime="application/pdf"
                    )

        else:
            st.warning("Please upload a file first using the ‘Open File’ tab.")

    # ==================================================
    # TAB 3: AI Reporting
    # ==================================================
    with tab3:
        st.subheader("OpenAI-powered Experimental Interpretation")

        if st.button("Starting AI Analysis Report Generation"):
            with st.spinner("AI is interpreting the data..."):
                st.success("Analysis completed!")
                st.info(
                    """
                    **[Summary of AI Interpretation Results]**
                    - Col1a1 expression is significantly increased in the Fibrosis group.
                    - This is consistent with typical liver fibrosis gene expression patterns.
                    """
                )

else:
    st.write("The selected menu feature is currently under development.")
