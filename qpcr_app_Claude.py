#!/usr/bin/env python3
"""
qPCR MIQE Statistical Analysis App
Steps 1–8: Data validation → Preprocessing → Statistical Testing →
           FDR Correction → QC → Tables → Figures → Manuscript Text
+ PDF Report Download
"""

import io, warnings, tempfile, base64
from datetime import date
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import shapiro, levene, mannwhitneyu, f as f_dist, kruskal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
warnings.filterwarnings("ignore")

import streamlit as st

# ── ReportLab PDF ────────────────────────────────────────────
from reportlab.lib.pagesizes import letter
from reportlab.lib import colors as rl_colors
from reportlab.lib.units import inch
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                 TableStyle, Image as RLImage,
                                 PageBreak, HRFlowable)
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY

# ═══════════════════════════════════════════════════════════════
# PAGE CONFIG
# ═══════════════════════════════════════════════════════════════
# st.set_page_config(  # Disabled: called in main_app.py
#     page_title="qPCR MIQE Analyzer",
#     page_icon="🧬",
#     layout="wide",
#     initial_sidebar_state="expanded",
# )

# ── Custom CSS ────────────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════
REFERENCE_GENES = {'ACTB','GAPDH','18S','HPRT1','RPLP0','TBP','BACTIN','B-ACTIN','B_ACTIN'}
PALETTE = {
    0: '#4878CF', 1: '#D65F5F', 2: '#6ACC65', 3: '#B47CC7',
    4: '#C4AD66', 5: '#77BEDB', 6: '#F28E2B', 7: '#59A14F',
}

# ═══════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════

def sig_stars(p):
    if p is None: return ''
    if p < 0.001: return '***'
    elif p < 0.01: return '**'
    elif p < 0.05: return '*'
    return 'ns'

def cohen_d(a, b):
    return (np.mean(a) - np.mean(b)) / np.sqrt((np.std(a,ddof=1)**2 + np.std(b,ddof=1)**2)/2)

def rank_biserial(u, n1, n2):
    return 1 - (2*u)/(n1*n2)

def bh_correction(p_vals):
    n = len(p_vals)
    if n == 0: return np.array([])
    order = np.argsort(p_vals)
    q = np.zeros(n)
    for rank, idx in enumerate(order, 1):
        q[idx] = p_vals[idx] * n / rank
    for i in range(n-2, -1, -1):
        q[order[i]] = min(q[order[i]], q[order[i+1]])
    return np.minimum(q, 1.0)

def holm_correction(p_vals):
    n = len(p_vals)
    order = np.argsort(p_vals)
    q = np.zeros(n)
    for rank, idx in enumerate(order):
        q[idx] = min(p_vals[idx] * (n - rank), 1.0)
    for i in range(n-2,-1,-1):
        q[order[i]] = min(q[order[i]], q[order[i+1]])
    return np.minimum(q, 1.0)

def two_way_anova(data, dv, fA, fB):
    grand = data[dv].mean()
    la = sorted(data[fA].unique()); lb = sorted(data[fB].unique())
    a = len(la); b = len(lb)
    cell_m = data.groupby([fA,fB])[dv].mean()
    cell_n = data.groupby([fA,fB])[dv].count()
    marg_a = data.groupby(fA)[dv].mean(); n_a = data.groupby(fA)[dv].count()
    marg_b = data.groupby(fB)[dv].mean(); n_b = data.groupby(fB)[dv].count()
    SS_A  = sum(n_a[ai]*(marg_a[ai]-grand)**2 for ai in la)
    SS_B  = sum(n_b[bi]*(marg_b[bi]-grand)**2 for bi in lb)
    SS_cells = sum(cell_n.get((ai,bi),0)*(cell_m.get((ai,bi),grand)-grand)**2
                   for ai in la for bi in lb)
    SS_AB = SS_cells - SS_A - SS_B
    SS_W  = sum(((data[(data[fA]==ai)&(data[fB]==bi)][dv]-
                  data[(data[fA]==ai)&(data[fB]==bi)][dv].mean())**2).sum()
                for ai in la for bi in lb)
    N = len(data)
    dfA=a-1; dfB=b-1; dfAB=dfA*dfB; dfW=N-a*b
    MSA=SS_A/dfA; MSB=SS_B/dfB; MSAB=SS_AB/dfAB; MSW=SS_W/dfW
    FA=MSA/MSW; FB=MSB/MSW; FAB=MSAB/MSW
    pA=1-f_dist.cdf(FA,dfA,dfW); pB=1-f_dist.cdf(FB,dfB,dfW); pAB=1-f_dist.cdf(FAB,dfAB,dfW)
    return {
        'Source':['Factor 1 (Condition1)','Factor 2 (Condition2)','Interaction (C1×C2)','Within Error'],
        'SS':[SS_A,SS_B,SS_AB,SS_W], 'df':[dfA,dfB,dfAB,dfW],
        'MS':[MSA,MSB,MSAB,MSW], 'F':[FA,FB,FAB,None], 'p':[pA,pB,pAB,None],
        'eta2':[SS_A/(SS_A+SS_W), SS_B/(SS_B+SS_W), SS_AB/(SS_AB+SS_W), None],
        'dfW': dfW,
    }

def fig_to_bytes(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=200, bbox_inches='tight', facecolor='white')
    buf.seek(0)
    return buf.read()

# ═══════════════════════════════════════════════════════════════
# CORE ANALYSIS ENGINE
# ═══════════════════════════════════════════════════════════════

def run_analysis(df):
    """Full MIQE pipeline. Returns dict with all results."""
    res = {'errors': [], 'warnings': [], 'steps': {}}

    # ── STEP 1: Column mapping ────────────────────────────────
    col_map = {}
    for c in df.columns:
        cl = c.lower().strip()
        if 'sample' in cl and 'name' in cl: col_map['Sample_Name'] = c
        elif cl in ('gene','target','target name','target_name'): col_map['Gene'] = c
        elif cl in ('ct','cp','cq','cycle threshold'): col_map['CT'] = c
        elif 'replicate' in cl or cl == 'rep': col_map['Replicate'] = c
        elif 'condition1' in cl or cl == 'condition 1': col_map['Condition1'] = c
        elif 'condition2' in cl or cl == 'condition 2': col_map['Condition2'] = c

    # fallback positional detection
    if 'CT' not in col_map:
        for c in df.columns:
            try:
                pd.to_numeric(df[c], errors='raise')
                if df[c].between(5,45).mean() > 0.5:
                    col_map['CT'] = c; break
            except: pass

    missing = [k for k in ['Sample_Name','Gene','CT','Replicate'] if k not in col_map]
    if missing:
        res['errors'].append(f"Missing required columns: {missing}")
        return res

    dw = df.rename(columns={v:k for k,v in col_map.items()}).copy()
    dw['Sample_Name'] = dw['Sample_Name'].astype(str).str.strip()
    dw['CT'] = pd.to_numeric(dw['CT'], errors='coerce')
    dw['Replicate'] = pd.to_numeric(dw['Replicate'], errors='coerce').astype('Int64')

    has_c1 = 'Condition1' in dw.columns
    has_c2 = 'Condition2' in dw.columns and dw['Condition2'].notna().any()

    # ── STEP 1B: Validation ───────────────────────────────────
    ref_genes = [g for g in dw['Gene'].unique() if str(g).upper().replace('-','_').replace(' ','_') in REFERENCE_GENES]
    tgt_genes = [g for g in dw['Gene'].unique() if g not in ref_genes]

    if not ref_genes:
        res['errors'].append("No reference gene detected. Expected: ACTB, GAPDH, BACTIN, 18S, HPRT1, RPLP0, TBP.")
    if not tgt_genes:
        res['errors'].append("No target gene detected.")
    if dw['CT'].isna().mean() > 0.05:
        res['errors'].append(f"CT missing rate = {dw['CT'].isna().mean()*100:.1f}% (>5% threshold).")
    if has_c1:
        for g in dw['Condition1'].unique():
            n = dw[dw['Condition1']==g]['Sample_Name'].nunique()
            if n < 2: res['errors'].append(f"Group '{g}' has <2 biological replicates.")

    if res['errors']:
        return res

    res['steps']['step1'] = {
        'ref_genes': ref_genes, 'tgt_genes': tgt_genes,
        'has_c1': has_c1, 'has_c2': has_c2,
        'n_samples': dw['Sample_Name'].nunique(),
        'n_genes': len(tgt_genes),
        'conditions': {
            'Condition1': sorted(dw['Condition1'].unique().tolist()) if has_c1 else [],
            'Condition2': sorted(dw['Condition2'].unique().tolist()) if has_c2 else [],
        }
    }

    # ── STEP 2: Preprocessing ─────────────────────────────────
    grp_cols = ['Sample_Name','Gene']
    if has_c1: grp_cols.append('Condition1')
    if has_c2: grp_cols.append('Condition2')

    mc = dw.groupby(grp_cols)['CT'].agg(
        Mean_CT='mean', Min_CT='min', Max_CT='max').reset_index()
    mc['CT_Range'] = mc['Max_CT'] - mc['Min_CT']
    mc['QC_Flag'] = mc['CT_Range'] > 0.5

    ref_ct_df = mc[mc['Gene'].isin(ref_genes)].copy()
    if len(ref_genes) == 1:
        ref_ct_df = ref_ct_df.rename(columns={'Mean_CT':'Ref_CT'})
    else:
        # geometric mean of ref genes
        merge_cols = ['Sample_Name'] + ([c for c in ['Condition1','Condition2'] if c in mc.columns])
        ref_ct_df = ref_ct_df.groupby(merge_cols)['Mean_CT'].apply(
            lambda x: np.exp(np.mean(np.log(x)))).reset_index()
        ref_ct_df.columns = merge_cols + ['Ref_CT']

    ref_merge_on = ['Sample_Name'] + ([c for c in ['Condition1','Condition2'] if c in mc.columns])
    if 'Ref_CT' not in ref_ct_df.columns:
        ref_ct_df = ref_ct_df[ref_merge_on + ['Ref_CT']]

    tgt_df = mc[mc['Gene'].isin(tgt_genes)].copy()
    merged = tgt_df.merge(ref_ct_df[ref_merge_on + ['Ref_CT']], on=ref_merge_on, how='left')
    merged['Delta_CT'] = merged['Mean_CT'] - merged['Ref_CT']

    # Control detection
    ctrl_c1 = None
    if has_c1:
        for candidate in ['Control','WT','Vehicle','control','wt','vehicle']:
            if candidate in merged['Condition1'].values:
                ctrl_c1 = candidate; break
        if ctrl_c1 is None:
            ctrl_c1 = sorted(merged['Condition1'].unique())[0]

    ctrl_c2 = None
    if has_c2:
        for candidate in ['0h','0','Control','control','0min']:
            if candidate in merged['Condition2'].values:
                ctrl_c2 = candidate; break
        if ctrl_c2 is None:
            ctrl_c2 = sorted(merged['Condition2'].unique())[0]

    def get_ctrl_mean(gene):
        if has_c2:
            mask = ((merged['Gene']==gene) &
                    (merged['Condition1']==ctrl_c1) &
                    (merged['Condition2']==ctrl_c2))
        elif has_c1:
            mask = (merged['Gene']==gene) & (merged['Condition1']==ctrl_c1)
        else:
            mask = (merged['Gene']==gene)
        return merged[mask]['Delta_CT'].mean()

    merged['Ctrl_Mean_dCt'] = merged['Gene'].map(lambda g: get_ctrl_mean(g))
    merged['Delta_Delta_CT'] = merged['Delta_CT'] - merged['Ctrl_Mean_dCt']
    merged['Fold_Change'] = 2**(-merged['Delta_Delta_CT'])
    merged['Log2_FC'] = np.log2(merged['Fold_Change'])

    # IQR outliers
    iqr_flags = []
    for gene in tgt_genes:
        for cond in (merged['Condition1'].unique() if has_c1 else [None]):
            msk = (merged['Gene']==gene)
            if cond: msk &= (merged['Condition1']==cond)
            sub = merged[msk]
            if len(sub) >= 4:
                q1,q3 = np.percentile(sub['Delta_CT'],[25,75])
                iqr = q3-q1
                out = sub[(sub['Delta_CT'] < q1-1.5*iqr) | (sub['Delta_CT'] > q3+1.5*iqr)]
                for _, r in out.iterrows():
                    iqr_flags.append({'Sample': r['Sample_Name'], 'Gene': gene,
                                       'dCt': r['Delta_CT'], 'Condition': cond})

    res['steps']['step2'] = {
        'merged': merged, 'mc': mc,
        'flagged': mc[mc['QC_Flag']],
        'ref_genes': ref_genes, 'tgt_genes': tgt_genes,
        'ctrl_c1': ctrl_c1, 'ctrl_c2': ctrl_c2,
        'has_c1': has_c1, 'has_c2': has_c2,
        'iqr_flags': iqr_flags,
    }

    # ── STEP 3 + 4: Statistical Testing + FDR ────────────────
    stat_results = {}

    for gene in tgt_genes:
        gdata = merged[merged['Gene']==gene].copy()
        sr = {'gene': gene, 'mode': None}

        if not has_c1:
            sr['mode'] = 'descriptive'
        elif has_c2:
            # Two-way ANOVA
            sr['mode'] = 'two_way'
            norm_cells = {}
            for c1 in sorted(gdata['Condition1'].unique()):
                for c2 in sorted(gdata['Condition2'].unique()):
                    vals = gdata[(gdata['Condition1']==c1)&(gdata['Condition2']==c2)]['Delta_CT'].values
                    if len(vals) >= 3:
                        w,p = shapiro(vals)
                        norm_cells[f"{c1}×{c2}"] = {'W':w,'p':p,'n':len(vals),'normal':p>0.05}

            groups_all = [gdata[(gdata['Condition1']==c1)&(gdata['Condition2']==c2)]['Delta_CT'].values
                          for c1 in gdata['Condition1'].unique()
                          for c2 in gdata['Condition2'].unique()]
            lev_s, lev_p = levene(*[g for g in groups_all if len(g)>1])

            anova = two_way_anova(gdata, 'Delta_CT', 'Condition1', 'Condition2')
            int_p = anova['p'][2]

            # Post-hoc pairwise
            combos = [(c1,c2) for c1 in sorted(gdata['Condition1'].unique())
                               for c2 in sorted(gdata['Condition2'].unique())]
            pairs = list(combinations(combos,2))
            ph_raw, ph_comps = [], []
            for (a1,b1),(a2,b2) in pairs:
                v1 = gdata[(gdata['Condition1']==a1)&(gdata['Condition2']==b1)]['Delta_CT'].values
                v2 = gdata[(gdata['Condition1']==a2)&(gdata['Condition2']==b2)]['Delta_CT'].values
                if len(v1) >= 2 and len(v2) >= 2:
                    t,p = stats.ttest_ind(v1,v2,equal_var=False)
                    d = cohen_d(v1,v2)
                    ph_comps.append({'g1':f"{a1}/{b1}",'g2':f"{a2}/{b2}",'t':t,'p':p,'d':d})
                    ph_raw.append(p)
            if ph_raw:
                ph_holm = holm_correction(np.array(ph_raw))
                for i,c in enumerate(ph_comps): c['p_holm']=ph_holm[i]; c['sig']=ph_holm[i]<0.05

            sr.update({'anova_2way': anova, 'norm_cells': norm_cells,
                       'levene': (lev_s,lev_p), 'posthoc': ph_comps,
                       'p_main': anova['p'][0], 'p_time': anova['p'][1],
                       'p_int': int_p,
                       'eta2_main': anova['eta2'][0],
                       'F_main': anova['F'][0], 'dfW': anova['dfW']})
        else:
            # One factor
            groups = sorted(gdata['Condition1'].unique())
            n_groups = len(groups)
            group_data = {g: gdata[gdata['Condition1']==g]['Delta_CT'].values for g in groups}
            ns = {g: len(v) for g,v in group_data.items()}

            if n_groups == 2:
                sr['mode'] = 'two_group'
                g1, g2 = groups
                v1, v2 = group_data[g1], group_data[g2]
                n1, n2 = len(v1), len(v2)

                if n1 <= 5 or n2 <= 5:
                    u, p = mannwhitneyu(v1, v2, alternative='two-sided')
                    r = rank_biserial(u, n1, n2)
                    sr.update({'test':'Mann-Whitney U','stat':u,'p':p,
                                'effect_name':'rank-biserial r','effect':r,
                                'normality':None, 'g1':g1,'g2':g2,'n1':n1,'n2':n2})
                else:
                    sw1 = shapiro(v1); sw2 = shapiro(v2)
                    normal = sw1.pvalue > 0.05 and sw2.pvalue > 0.05
                    lev_s, lev_p = levene(v1, v2)
                    sr['normality'] = {g1:{'W':sw1.statistic,'p':sw1.pvalue,'normal':sw1.pvalue>0.05},
                                        g2:{'W':sw2.statistic,'p':sw2.pvalue,'normal':sw2.pvalue>0.05}}
                    sr['levene'] = (lev_s, lev_p)
                    if normal:
                        eq_var = lev_p > 0.05
                        t, p = stats.ttest_ind(v1, v2, equal_var=eq_var)
                        d = cohen_d(v1, v2)
                        test = "Student's t-test" if eq_var else "Welch's t-test"
                        sr.update({'test':test,'stat':t,'p':p,
                                    'effect_name':"Cohen's d",'effect':d,
                                    'g1':g1,'g2':g2,'n1':n1,'n2':n2})
                    else:
                        u, p = mannwhitneyu(v1, v2, alternative='two-sided')
                        r = rank_biserial(u, n1, n2)
                        sr.update({'test':'Mann-Whitney U','stat':u,'p':p,
                                    'effect_name':'rank-biserial r','effect':r,
                                    'g1':g1,'g2':g2,'n1':n1,'n2':n2})
            else:
                sr['mode'] = 'multi_group'
                any_small = any(ns[g] <= 5 for g in groups)
                if any_small:
                    h, p = kruskal(*[group_data[g] for g in groups])
                    n_total = sum(ns.values())
                    eta2_e = (h - len(groups) + 1) / (n_total - len(groups))
                    sr.update({'test':'Kruskal-Wallis','stat':h,'p':p,
                                'effect_name':'epsilon²','effect':max(eta2_e,0),'groups':groups,'group_data':group_data})
                else:
                    all_norm = all(shapiro(group_data[g]).pvalue > 0.05 for g in groups)
                    if all_norm:
                        f_stat, p = stats.f_oneway(*[group_data[g] for g in groups])
                        grand = np.concatenate(list(group_data.values()))
                        ss_between = sum(ns[g]*(np.mean(group_data[g])-np.mean(grand))**2 for g in groups)
                        ss_total = sum((x-np.mean(grand))**2 for x in grand)
                        eta2 = ss_between/ss_total
                        sr.update({'test':'One-Way ANOVA','stat':f_stat,'p':p,
                                    'effect_name':'η²','effect':eta2,'groups':groups,'group_data':group_data})
                    else:
                        h, p = kruskal(*[group_data[g] for g in groups])
                        n_total = sum(ns.values())
                        eta2_e = (h - len(groups) + 1) / (n_total - len(groups))
                        sr.update({'test':'Kruskal-Wallis','stat':h,'p':p,
                                    'effect_name':'epsilon²','effect':max(eta2_e,0),'groups':groups,'group_data':group_data})

        stat_results[gene] = sr

    # BH FDR
    p_vals_for_fdr = []
    for gene in tgt_genes:
        sr = stat_results[gene]
        if sr['mode'] == 'two_way': p_vals_for_fdr.append(sr.get('p_main', 1.0))
        elif sr['mode'] in ('two_group','multi_group'): p_vals_for_fdr.append(sr.get('p', 1.0))
        else: p_vals_for_fdr.append(1.0)

    q_vals = bh_correction(np.array(p_vals_for_fdr))
    for i, gene in enumerate(tgt_genes):
        stat_results[gene]['q'] = q_vals[i]
        stat_results[gene]['sig_fdr'] = q_vals[i] < 0.05

    res['steps']['step3'] = stat_results
    res['steps']['step4'] = {'q_vals': dict(zip(tgt_genes, q_vals))}

    # Descriptive table
    desc_rows = []
    for gene in tgt_genes:
        gdata = merged[merged['Gene']==gene]
        grp_by = (['Condition1'] if has_c1 else []) + (['Condition2'] if has_c2 else [])
        if grp_by:
            for keys, sub in gdata.groupby(grp_by):
                row = {'Gene': gene}
                if not isinstance(keys, tuple): keys = (keys,)
                for k, col in zip(keys, grp_by): row[col] = k
                row.update({'n': len(sub),
                    'Mean_CT': sub['Mean_CT'].mean(), 'Mean_dCt': sub['Delta_CT'].mean(),
                    'SD_dCt': sub['Delta_CT'].std(ddof=1), 'SEM_dCt': sub['Delta_CT'].sem(),
                    'Mean_FC': sub['Fold_Change'].mean(), 'SEM_FC': sub['Fold_Change'].sem(),
                    'Mean_log2FC': sub['Log2_FC'].mean(), 'SEM_log2FC': sub['Log2_FC'].sem()})
                desc_rows.append(row)
        else:
            sub = gdata
            desc_rows.append({'Gene': gene, 'n': len(sub),
                'Mean_CT': sub['Mean_CT'].mean(), 'Mean_dCt': sub['Delta_CT'].mean(),
                'SD_dCt': sub['Delta_CT'].std(ddof=1), 'SEM_dCt': sub['Delta_CT'].sem(),
                'Mean_FC': sub['Fold_Change'].mean(), 'SEM_FC': sub['Fold_Change'].sem(),
                'Mean_log2FC': sub['Log2_FC'].mean(), 'SEM_log2FC': sub['Log2_FC'].sem()})

    res['steps']['step6'] = {'desc_df': pd.DataFrame(desc_rows)}
    return res

# ═══════════════════════════════════════════════════════════════
# FIGURE GENERATORS
# ═══════════════════════════════════════════════════════════════

def make_expression_plot(merged, stat_results, tgt_genes, has_c1, has_c2, ctrl_c1, ctrl_c2):
    n = len(tgt_genes)
    fig, axes = plt.subplots(1, n, figsize=(5.5*n, 5.5), dpi=150)
    if n == 1: axes = [axes]
    fig.patch.set_facecolor('white')

    for ax, gene in zip(axes, tgt_genes):
        gdata = merged[merged['Gene']==gene]
        sr = stat_results[gene]

        if has_c2:
            groups = [(c1,c2) for c1 in sorted(gdata['Condition1'].unique())
                               for c2 in sorted(gdata['Condition2'].unique())]
            labels = [f"{c1}\n{c2}" for c1,c2 in groups]
            means = []; sems = []; pts_list = []
            for c1,c2 in groups:
                sub = gdata[(gdata['Condition1']==c1)&(gdata['Condition2']==c2)]['Log2_FC'].values
                means.append(np.mean(sub)); sems.append(np.std(sub,ddof=1)/np.sqrt(len(sub)) if len(sub)>1 else 0)
                pts_list.append(sub)
            bar_colors = [PALETTE[i % len(PALETTE)] for i in range(len(groups))]
        elif has_c1:
            groups_c1 = sorted(gdata['Condition1'].unique())
            labels = groups_c1
            means = []; sems = []; pts_list = []
            for g in groups_c1:
                sub = gdata[gdata['Condition1']==g]['Log2_FC'].values
                means.append(np.mean(sub)); sems.append(np.std(sub,ddof=1)/np.sqrt(len(sub)) if len(sub)>1 else 0)
                pts_list.append(sub)
            bar_colors = [PALETTE[i % len(PALETTE)] for i in range(len(groups_c1))]
        else:
            sub = gdata['Log2_FC'].values
            labels = [gene]; means = [np.mean(sub)]
            sems = [np.std(sub,ddof=1)/np.sqrt(len(sub)) if len(sub)>1 else 0]
            pts_list = [sub]; bar_colors = [PALETTE[0]]

        x = np.arange(len(labels))
        bars = ax.bar(x, means, width=0.52, color=bar_colors, alpha=0.82,
                      edgecolor='black', linewidth=0.7, zorder=3)
        ax.errorbar(x, means, yerr=sems, fmt='none', color='#333333',
                    capsize=5, capthick=1.2, linewidth=1.2, zorder=4)

        np.random.seed(42)
        for i, pts in enumerate(pts_list):
            jit = np.random.uniform(-0.09,0.09,len(pts))
            ax.scatter(i+jit, pts, s=32, color='#222222', alpha=0.7, zorder=5)

        # Significance annotation
        if sr['mode'] == 'two_group':
            p = sr.get('p',1); q = sr.get('q',1)
            y_top = max(means)+max(sems)+0.3
            ax.plot([0,0,1,1],[y_top,y_top+0.2,y_top+0.2,y_top],'k-',lw=1.1)
            ax.text(0.5, y_top+0.25, f"{sig_stars(p)}  p={p:.3f}\nq={q:.3f}",
                    ha='center',va='bottom',fontsize=8.5,fontweight='bold')
            ax.set_ylim(min(min(p2) for p2 in pts_list)-1.0, y_top+1.4)
        elif sr['mode'] == 'two_way':
            p = sr.get('p_main',1); q = sr.get('q',1)
            ax.text(0.02,0.97,
                f"Treatment: p={p:.4f} {sig_stars(p)}\nTime: p={sr.get('p_time',1):.4f}\nInteraction: p={sr.get('p_int',1):.4f}",
                transform=ax.transAxes,va='top',ha='left',fontsize=8,
                bbox=dict(boxstyle='round,pad=0.3',facecolor='#f8f8f8',edgecolor='#cccccc',alpha=0.9))

        ax.axhline(0, color='#aaaaaa', linestyle='--', linewidth=0.8, alpha=0.6)
        ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9.5)
        ax.set_ylabel('log₂ Fold Change (mean ± SEM)', fontsize=10)
        ax.set_title(gene, fontsize=13, fontweight='bold', pad=8)
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
        ax.tick_params(axis='both', labelsize=9)

    fig.suptitle('Relative Gene Expression  |  Statistics on ΔCt  |  BH FDR corrected',
                 fontsize=10, y=1.01, color='#555555')
    plt.tight_layout()
    return fig

def make_interaction_plot(merged, tgt_genes):
    fig, axes = plt.subplots(1, len(tgt_genes), figsize=(5*len(tgt_genes), 4.5), dpi=150)
    if len(tgt_genes)==1: axes=[axes]
    fig.patch.set_facecolor('white')
    line_c = {'Control':'#4878CF','Drug':'#D65F5F','WT':'#4878CF','Vehicle':'#4878CF'}

    for ax, gene in zip(axes, tgt_genes):
        gdata = merged[merged['Gene']==gene]
        c1_levels = sorted(gdata['Condition1'].unique())
        c2_levels = sorted(gdata['Condition2'].unique())

        for i, c1 in enumerate(c1_levels):
            means_i, sems_i = [],[]
            for c2 in c2_levels:
                sub = gdata[(gdata['Condition1']==c1)&(gdata['Condition2']==c2)]['Delta_CT'].values
                means_i.append(np.mean(sub))
                sems_i.append(np.std(sub,ddof=1)/np.sqrt(len(sub)) if len(sub)>1 else 0)
            col = line_c.get(c1, PALETTE[i % len(PALETTE)])
            mk = ['o','s','^','D'][i % 4]
            ax.plot(range(len(c2_levels)), means_i, color=col, linewidth=2.2,
                    marker=mk, markersize=9, label=c1, zorder=4)
            ax.errorbar(range(len(c2_levels)), means_i, yerr=sems_i, fmt='none',
                        color=col, capsize=5, capthick=1.2, linewidth=1.2)
            np.random.seed(42)
            for j, c2 in enumerate(c2_levels):
                sub = gdata[(gdata['Condition1']==c1)&(gdata['Condition2']==c2)]['Delta_CT'].values
                jit = np.random.uniform(-0.06,0.06,len(sub))
                ax.scatter(j+jit, sub, s=22, color=col, alpha=0.45, zorder=3)

        ax.set_xticks(range(len(c2_levels))); ax.set_xticklabels(c2_levels, fontsize=10)
        ax.set_xlabel('Condition 2 (Time)', fontsize=10); ax.set_ylabel('ΔCt (mean ± SEM)', fontsize=10)
        ax.set_title(f'{gene} — Interaction Plot', fontsize=12, fontweight='bold')
        ax.legend(title='Treatment', fontsize=9, title_fontsize=9, framealpha=0.9)
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

    plt.tight_layout(); return fig

def make_heatmap(merged, tgt_genes, has_c1, has_c2):
    fig, ax = plt.subplots(figsize=(max(4, len(tgt_genes)*2.2), 7), dpi=150)
    fig.patch.set_facecolor('white')

    if has_c2:
        groups = [(c1,c2) for c1 in sorted(merged['Condition1'].unique())
                           for c2 in sorted(merged['Condition2'].unique())]
        group_labels = [f"{c1}/{c2}" for c1,c2 in groups]
    elif has_c1:
        groups = [(c1,) for c1 in sorted(merged['Condition1'].unique())]
        group_labels = [c1 for (c1,) in groups]
    else:
        groups = [()]; group_labels = ['All']

    matrix = []
    row_labels = []
    group_col = []

    for gi, g in enumerate(groups):
        if has_c2:
            sub_all = merged[(merged['Condition1']==g[0])&(merged['Condition2']==g[1])]
        elif has_c1:
            sub_all = merged[merged['Condition1']==g[0]]
        else:
            sub_all = merged

        for _, row in sub_all.drop_duplicates('Sample_Name').iterrows():
            sample_row = []
            for gene in tgt_genes:
                val = merged[(merged['Sample_Name']==row['Sample_Name'])&
                             (merged['Gene']==gene)]['Log2_FC'].values
                sample_row.append(val[0] if len(val)>0 else 0)
            matrix.append(sample_row)
            row_labels.append(f"{row['Sample_Name']} ({group_labels[gi]})")
            group_col.append(gi)

    mat = np.array(matrix)
    vmax = max(abs(mat).max(), 1)
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    im = ax.imshow(mat, aspect='auto', cmap='RdBu_r', norm=norm)

    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            tc = 'white' if abs(mat[i,j]) > vmax*0.65 else 'black'
            ax.text(j, i, f"{mat[i,j]:.2f}", ha='center', va='center',
                    fontsize=8, color=tc, fontweight='bold')

    # Group separators
    prev = None
    for i, gc in enumerate(group_col):
        if gc != prev and i > 0:
            ax.axhline(i-0.5, color='white', linewidth=2.5)
        prev = gc

    ax.set_yticks(range(len(row_labels))); ax.set_yticklabels(row_labels, fontsize=8)
    ax.set_xticks(range(len(tgt_genes))); ax.set_xticklabels(tgt_genes, fontsize=10, fontweight='bold')
    ax.set_title('Heatmap: log₂ Fold Change per Sample\n(blue=down · white=neutral · red=up)',
                 fontsize=11, fontweight='bold', pad=10)
    cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
    cbar.set_label('log₂ Fold Change', fontsize=9)
    plt.tight_layout(); return fig

def make_volcano(merged, stat_results, tgt_genes):
    fig, ax = plt.subplots(figsize=(6.5, 5.5), dpi=150)
    fig.patch.set_facecolor('white')

    xs, ys, labels_v, sigs, qs = [], [], [], [], []
    for gene in tgt_genes:
        sr = stat_results[gene]
        gdata = merged[merged['Gene']==gene]
        if sr['mode'] == 'two_group':
            g1,g2 = sr['g1'], sr['g2']
            v1 = gdata[gdata['Condition1']==g1]['Log2_FC'].values
            v2 = gdata[gdata['Condition1']==g2]['Log2_FC'].values
            fc = np.mean(v2) - np.mean(v1)
            p = sr.get('p',1); q = sr.get('q',1)
        elif sr['mode'] == 'two_way':
            v_ctrl = gdata[(gdata['Condition1']==gdata['Condition1'].mode()[0])]['Log2_FC'].values
            fc = np.mean(gdata['Log2_FC'].values)
            p = sr.get('p_main',1); q = sr.get('q',1)
        elif sr['mode'] == 'multi_group':
            fc = np.mean(gdata['Log2_FC'].values)
            p = sr.get('p',1); q = sr.get('q',1)
        else:
            continue
        xs.append(fc); ys.append(-np.log10(max(p,1e-10)))
        labels_v.append(gene); sigs.append(q < 0.05); qs.append(q)

    col_pts = ['#D65F5F' if s else '#AAAAAA' for s in sigs]
    sz_pts = [90 if s else 55 for s in sigs]
    ax.scatter(xs, ys, c=col_pts, s=sz_pts, alpha=0.9, edgecolors='black', linewidths=0.5, zorder=4)

    for x, y, lbl, sig in zip(xs, ys, labels_v, sigs):
        ax.annotate(lbl, xy=(x,y), xytext=(x+0.08, y+0.12),
                    fontsize=9, fontweight='bold' if sig else 'normal',
                    color='#c0392b' if sig else '#555555',
                    arrowprops=dict(arrowstyle='-', color='#aaaaaa', lw=0.8))

    ax.axhline(-np.log10(0.05), color='#666666', linestyle='--', linewidth=1.2, alpha=0.7, label='p=0.05')
    ax.axvline(0, color='#dddddd', linewidth=0.8)
    ax.axvspan(-1,1, alpha=0.04, color='gray')
    ax.set_xlabel('log₂ Fold Change', fontsize=11)
    ax.set_ylabel('−log₁₀(p-value)', fontsize=11)
    ax.set_title('Volcano Plot (BH FDR q < 0.05 highlighted)', fontsize=12, fontweight='bold')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    legend_elems = [
        mpatches.Patch(facecolor='#D65F5F', edgecolor='black', label='FDR significant (q<0.05)'),
        mpatches.Patch(facecolor='#AAAAAA', edgecolor='black', label='Not significant'),
        Line2D([0],[0], color='#666666', linestyle='--', label='p=0.05'),
    ]
    ax.legend(handles=legend_elems, fontsize=9, framealpha=0.9)
    plt.tight_layout(); return fig

# ═══════════════════════════════════════════════════════════════
# PDF BUILDER
# ═══════════════════════════════════════════════════════════════

def build_pdf(res, fig_bytes_dict):
    buf = io.BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=letter,
        leftMargin=0.85*inch, rightMargin=0.85*inch,
        topMargin=0.9*inch, bottomMargin=0.9*inch)

    S = getSampleStyleSheet()
    TS = ParagraphStyle('TT', parent=S['Title'], fontSize=19, spaceAfter=4,
        textColor=rl_colors.HexColor('#1a3a5c'), alignment=TA_CENTER)
    H1 = ParagraphStyle('H1', parent=S['Heading1'], fontSize=12,
        textColor=rl_colors.HexColor('#1a3a5c'), spaceBefore=14, spaceAfter=6)
    H2 = ParagraphStyle('H2', parent=S['Heading2'], fontSize=10,
        textColor=rl_colors.HexColor('#2c5f8a'), spaceBefore=10, spaceAfter=4)
    NM = ParagraphStyle('NM', parent=S['Normal'], fontSize=9, leading=13.5,
        spaceAfter=4, alignment=TA_JUSTIFY)
    CP = ParagraphStyle('CP', parent=S['Normal'], fontSize=8,
        textColor=rl_colors.HexColor('#555'), alignment=TA_CENTER, spaceAfter=5)
    BN = ParagraphStyle('BN', parent=S['Normal'], fontSize=9, fontName='Helvetica-Bold', leading=13.5)
    WN = ParagraphStyle('WN', parent=S['Normal'], fontSize=8.5,
        textColor=rl_colors.HexColor('#cc5500'), leading=12)

    def mkT(data, cw=None):
        t = Table(data, colWidths=cw, repeatRows=1)
        t.setStyle(TableStyle([
            ('BACKGROUND',(0,0),(-1,0), rl_colors.HexColor('#1a3a5c')),
            ('TEXTCOLOR',(0,0),(-1,0), rl_colors.white),
            ('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),
            ('FONTSIZE',(0,0),(-1,0), 8),
            ('ALIGN',(0,0),(-1,-1),'CENTER'),
            ('VALIGN',(0,0),(-1,-1),'MIDDLE'),
            ('FONTSIZE',(0,1),(-1,-1), 7.5),
            ('ROWBACKGROUNDS',(0,1),(-1,-1),[rl_colors.white, rl_colors.HexColor('#f0f4f8')]),
            ('GRID',(0,0),(-1,-1), 0.35, rl_colors.HexColor('#cccccc')),
            ('TOPPADDING',(0,0),(-1,-1), 3.5),
            ('BOTTOMPADDING',(0,0),(-1,-1), 3.5),
            ('LEFTPADDING',(0,0),(-1,-1), 4),
        ]))
        return t

    def add_fig(story, key, w=5.5*inch, h=4*inch, caption_text=''):
        if key in fig_bytes_dict:
            tmp = tempfile.NamedTemporaryFile(suffix='.png', delete=False)
            tmp.write(fig_bytes_dict[key]); tmp.flush()
            story.append(RLImage(tmp.name, width=w, height=h))
            if caption_text:
                story.append(Paragraph(caption_text, CP))

    s1 = res['steps']['step1']
    s2 = res['steps']['step2']
    s3 = res['steps']['step3']
    s6 = res['steps']['step6']
    merged = s2['merged']
    tgt_genes = s2['tgt_genes']
    ref_genes = s2['ref_genes']
    flagged = s2['flagged']
    has_c1 = s2['has_c1']; has_c2 = s2['has_c2']
    ctrl_c1 = s2['ctrl_c1']; ctrl_c2 = s2['ctrl_c2']

    story = []

    # Title page
    story.append(Spacer(1,0.8*inch))
    story.append(Paragraph('qPCR Statistical Analysis Report', TS))
    story.append(Paragraph('MIQE-Compliant · Steps 1–8 Automated',
        ParagraphStyle('sub',parent=S['Normal'],fontSize=11,
            textColor=rl_colors.HexColor('#2c5f8a'),alignment=TA_CENTER,spaceAfter=4)))
    story.append(HRFlowable(width='100%',thickness=1.5,
        color=rl_colors.HexColor('#1a3a5c'),spaceAfter=12))
    story.append(Spacer(1,0.2*inch))

    mode_str = 'Two-Way ANOVA' if has_c2 else ('Two-Group Comparison' if has_c1 else 'Descriptive')
    gene_mode = f"Multi-Gene ({len(tgt_genes)} targets)" if len(tgt_genes) > 1 else f"Single-Gene ({tgt_genes[0]})"

    sum_rows = [
        ['Field','Value'],
        ['Date', date.today().strftime('%B %d, %Y')],
        ['Analysis Mode', f"{gene_mode} | {mode_str}"],
        ['Reference Gene(s)', ', '.join(ref_genes)],
        ['Target Gene(s)', ', '.join(tgt_genes)],
        ['Normalization', '2^(-ΔΔCt) method'],
        ['Factors', ('Condition1 × Condition2' if has_c2 else 'Condition1' if has_c1 else 'None')],
        ['Control Reference', (f"{ctrl_c1} × {ctrl_c2}" if has_c2 else (ctrl_c1 or 'N/A'))],
        ['Multiple Testing', 'Benjamini-Hochberg FDR'],
    ]
    t = Table([[Paragraph(r[0],BN), Paragraph(str(r[1]),NM)] for r in sum_rows],
              colWidths=[2*inch, 4.4*inch])
    t.setStyle(TableStyle([
        ('BACKGROUND',(0,0),(-1,0),rl_colors.HexColor('#1a3a5c')),
        ('TEXTCOLOR',(0,0),(-1,0),rl_colors.white),
        ('FONTSIZE',(0,0),(-1,-1),9),
        ('ALIGN',(0,0),(-1,-1),'LEFT'),
        ('VALIGN',(0,0),(-1,-1),'MIDDLE'),
        ('ROWBACKGROUNDS',(0,1),(-1,-1),[rl_colors.white,rl_colors.HexColor('#f0f4f8')]),
        ('GRID',(0,0),(-1,-1),0.35,rl_colors.HexColor('#cccccc')),
        ('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5),
        ('LEFTPADDING',(0,0),(-1,-1),7),
    ]))
    story.append(t)
    story.append(PageBreak())

    # QC
    story.append(Paragraph('1. Data & Quality Control', H1))
    story.append(HRFlowable(width='100%',thickness=0.5,color=rl_colors.HexColor('#ccc'),spaceAfter=8))

    story.append(Paragraph('1.1 Technical Replicate QC (Ct range > 0.5)', H2))
    if len(flagged):
        story.append(Paragraph('⚠ Flagged observations (retained, not removed):', WN))
        qd = [['Sample','Gene','Condition1','Condition2','Mean Ct','Range']]
        for _,r in flagged.iterrows():
            qd.append([r['Sample_Name'],r['Gene'],
                        str(r.get('Condition1','')),str(r.get('Condition2','')),
                        f"{r['Mean_CT']:.3f}",f"{r['CT_Range']:.3f}"])
        story.append(mkT(qd))
    else:
        story.append(Paragraph('✓ All technical replicates within acceptable range (≤ 0.5 Ct).', NM))

    story.append(Spacer(1,0.1*inch))
    story.append(Paragraph('1.2 Descriptive Statistics', H2))
    desc_df = s6['desc_df']
    dcols = list(desc_df.columns)
    dd = [dcols]
    for _, r in desc_df.iterrows():
        row_d = []
        for c in dcols:
            v = r[c]
            row_d.append(f"{v:.3f}" if isinstance(v,(float,np.floating)) else str(v))
        dd.append(row_d)
    cw_d = [max(0.5*inch, min(1.2*inch, 6.0*inch/len(dcols))) for _ in dcols]
    story.append(mkT(dd, cw_d))
    story.append(PageBreak())

    # Stats
    story.append(Paragraph('2. Statistical Results', H1))
    story.append(HRFlowable(width='100%',thickness=0.5,color=rl_colors.HexColor('#ccc'),spaceAfter=8))

    for gene in tgt_genes:
        sr = s3[gene]
        story.append(Paragraph(f'Gene: {gene}', H2))

        if sr['mode'] == 'two_way':
            # Normality table
            nc = sr.get('norm_cells',{})
            if nc:
                nd = [['Cell','n','Shapiro-Wilk W','p','Normal']]
                for cell, v in nc.items():
                    nd.append([cell,str(v['n']),f"{v['W']:.4f}",f"{v['p']:.4f}",
                                '✓' if v['normal'] else '✗'])
                story.append(mkT(nd, [1.6*inch,0.4*inch,1.2*inch,1.0*inch,0.8*inch]))
                story.append(Spacer(1,0.06*inch))

            # ANOVA table
            a2 = sr['anova_2way']
            atbl = [['Source','SS','df','MS','F','p','η²p','Sig']]
            for i, src in enumerate(a2['Source']):
                f_s = f"{a2['F'][i]:.3f}" if a2['F'][i] is not None else '—'
                p_s = f"{a2['p'][i]:.4f}" if a2['p'][i] is not None else '—'
                e_s = f"{a2['eta2'][i]:.3f}" if a2['eta2'][i] is not None else '—'
                sig = ''
                if a2['p'][i] is not None:
                    p_v = float(a2['p'][i])
                    sig = '***' if p_v<0.001 else '**' if p_v<0.01 else '*' if p_v<0.05 else 'ns'
                atbl.append([src,f"{a2['SS'][i]:.3f}",str(int(a2['df'][i])),
                               f"{a2['MS'][i]:.3f}",f_s,p_s,e_s,sig])
            at = mkT(atbl,[1.7*inch,0.7*inch,0.35*inch,0.7*inch,0.7*inch,0.65*inch,0.6*inch,0.5*inch])
            story.append(at)
            story.append(Spacer(1,0.06*inch))

            # Post-hoc
            ph = sr.get('posthoc',[])
            if ph:
                story.append(Paragraph('Post-hoc (Holm correction):', BN))
                phd = [['Comparison','t','Cohen\'s d','Raw p','Holm p','Sig']]
                for c in ph:
                    phd.append([f"{c['g1']} vs {c['g2']}",f"{c['t']:.3f}",
                                  f"{abs(c['d']):.3f}",f"{c['p']:.4f}",
                                  f"{c['p_holm']:.4f}",'✓' if c['sig'] else 'ns'])
                story.append(mkT(phd,[2.1*inch,0.6*inch,0.75*inch,0.75*inch,0.75*inch,0.5*inch]))

        elif sr['mode'] == 'two_group':
            if sr.get('normality'):
                nd = [['Group','n','Shapiro-Wilk W','p','Normal']]
                for g, v in sr['normality'].items():
                    nd.append([g,str(sr.get('n1' if g==sr['g1'] else 'n2','')),
                                f"{v['W']:.4f}",f"{v['p']:.4f}",'✓' if v['normal'] else '✗'])
                story.append(mkT(nd,[1.3*inch,0.4*inch,1.2*inch,1.0*inch,0.8*inch]))
                story.append(Spacer(1,0.06*inch))

            itbl = [['Test','Statistic','p-value','Effect Size','Value','BH q','Significant']]
            itbl.append([sr['test'],f"{sr['stat']:.4f}",f"{sr['p']:.4f}",
                          sr['effect_name'],f"{abs(sr['effect']):.4f}",
                          f"{sr['q']:.4f}",'YES ✓' if sr['sig_fdr'] else 'NO'])
            it = mkT(itbl,[1.2*inch,0.8*inch,0.7*inch,1.0*inch,0.7*inch,0.7*inch,0.8*inch])
            story.append(it)

        story.append(Spacer(1,0.1*inch))

    # Global BH table
    story.append(Paragraph('2.X Global FDR Summary', H2))
    gb = [['Gene','Test/Mode','p-value','BH q','Significant']]
    for gene in tgt_genes:
        sr = s3[gene]
        p_v = sr.get('p', sr.get('p_main','—'))
        gb.append([gene, sr.get('test', sr['mode']),
                   f"{p_v:.4f}" if isinstance(p_v,float) else str(p_v),
                   f"{sr['q']:.4f}", 'YES ✓' if sr['sig_fdr'] else 'NO'])
    story.append(mkT(gb,[1.0*inch,1.5*inch,0.9*inch,0.9*inch,0.9*inch]))
    story.append(PageBreak())

    # Figures
    story.append(Paragraph('3. Figures', H1))
    story.append(HRFlowable(width='100%',thickness=0.5,color=rl_colors.HexColor('#ccc'),spaceAfter=8))
    story.append(Paragraph('Figure 1 — Primary Expression Plot', H2))
    add_fig(story,'expr',5.5*inch,4.0*inch,
        'Figure 1. log2 fold change (mean ± SEM) with individual replicate points. Statistics on ΔCt. BH FDR corrected.')
    story.append(Spacer(1,0.12*inch))
    if 'inter' in fig_bytes_dict:
        story.append(Paragraph('Figure 2 — Interaction Plot', H2))
        add_fig(story,'inter',5.0*inch,3.5*inch,
            'Figure 2. Treatment × Time interaction plot (ΔCt ± SEM). Parallel lines indicate no interaction.')
    story.append(PageBreak())
    story.append(Paragraph('3. Figures (continued)', H1))
    story.append(HRFlowable(width='100%',thickness=0.5,color=rl_colors.HexColor('#ccc'),spaceAfter=8))
    story.append(Paragraph('Figure 3 — Sample-Level Heatmap', H2))
    add_fig(story,'heat',5.5*inch,4.2*inch,
        'Figure 3. Heatmap of log2 fold change per sample. Blue=down, white=neutral, red=up.')
    story.append(Spacer(1,0.12*inch))
    story.append(Paragraph('Figure 4 — Volcano Plot', H2))
    add_fig(story,'volc',5.2*inch,3.8*inch,
        'Figure 4. Volcano plot. Red=FDR significant (q<0.05). Dashed line: p=0.05.')
    story.append(PageBreak())

    # Methods + Results
    story.append(Paragraph('4. Manuscript-Ready Text', H1))
    story.append(HRFlowable(width='100%',thickness=0.5,color=rl_colors.HexColor('#ccc'),spaceAfter=8))
    story.append(Paragraph('4.1 Methods', H2))

    method_txt = (
        'Relative mRNA expression was quantified by reverse transcription-quantitative PCR (RT-qPCR) '
        f'following MIQE guidelines. {", ".join(ref_genes)} was used as the reference gene for normalization. '
        'Relative expression was calculated using the 2<super>-ΔΔCt</super> method. '
        '<b>All inferential statistical analyses were performed on ΔCt values.</b> '
    )
    gene_str = tgt_genes[0] if len(tgt_genes)==1 else ', '.join(tgt_genes)
    if has_c2:
        method_txt += (
            f'A 2×2 factorial two-way ANOVA was applied to {gene_str} ΔCt values with '
            'Treatment and Time as independent factors. Normality was assessed per cell using '
            'the Shapiro-Wilk test; homogeneity of variance was verified by Levene\'s test. '
            "Post-hoc pairwise comparisons used Welch's t-tests with Holm correction. "
            'Partial η² was reported as the effect size measure. '
        )
    elif has_c1:
        method_txt += (
            f'For {gene_str}, normality was assessed using the Shapiro-Wilk test and variance equality '
            "by Levene's test. Welch's t-test or Mann-Whitney U was applied depending on assumptions. "
            "Cohen's d or rank-biserial correlation was reported as the effect size. "
        )
    if len(tgt_genes) > 1:
        method_txt += 'Benjamini-Hochberg FDR correction was applied across all target genes. '
    method_txt += 'Technical replicate quality was assessed using a Ct range threshold of 0.5 cycles.'
    story.append(Paragraph(method_txt, NM))

    story.append(Spacer(1,0.12*inch))
    story.append(Paragraph('4.2 Results', H2))
    for gene in tgt_genes:
        sr = s3[gene]
        gdata = merged[merged['Gene']==gene]
        if sr['mode'] == 'two_group':
            fc_treat = gdata[gdata['Condition1']==sr['g2']]['Fold_Change'].mean()
            fc_sem = gdata[gdata['Condition1']==sr['g2']]['Fold_Change'].sem()
            direction = 'upregulated' if fc_treat > 1 else 'downregulated'
            txt = (
                f'<b>{gene}</b> was significantly {direction} in the {sr["g2"]} group '
                f'(mean FC = {fc_treat:.2f} ± {fc_sem:.2f} SEM; '
                f'{sr["test"]}: stat = {sr["stat"]:.3f}, p = {sr["p"]:.4f}, '
                f'{sr["effect_name"]} = {abs(sr["effect"]):.3f}; '
                f'BH q = {sr["q"]:.4f}).'
            )
        elif sr['mode'] == 'two_way':
            txt = (
                f'<b>{gene}</b>: Two-way ANOVA revealed a '
                f'{"significant" if sr["p_main"]<0.05 else "non-significant"} '
                f'main effect of Treatment (F(1,{sr["dfW"]}) = {sr["F_main"]:.3f}, '
                f'p = {sr["p_main"]:.4f}, partial η² = {sr["eta2_main"]:.3f}). '
                f'The main effect of Time was p = {sr["p_time"]:.4f} and '
                f'the interaction was p = {sr["p_int"]:.4f} (BH q = {sr["q"]:.4f}).'
            )
        else:
            txt = f'<b>{gene}</b>: Descriptive analysis only (no experimental conditions provided).'
        story.append(Paragraph(txt, NM))
        story.append(Spacer(1,0.06*inch))

    doc.build(story)
    buf.seek(0)
    return buf.getvalue()

# ═══════════════════════════════════════════════════════════════
# STREAMLIT UI  →  run_qpcr_module() 함수로 래핑
# ═══════════════════════════════════════════════════════════════

def run_qpcr_module():
    """
    qPCR MIQE Analysis Module entry point.
    Called via run_qpcr_module() from the
    sub_menu == '2. qPCR and Gene Expression Data Analysis' block in WODIS.py.
    - st.set_page_config() is handled in WODIS.py and removed here.
    - Sidebar widgets are replaced with st.expander to avoid sidebar conflicts.
    """

    st.markdown("""
    <div class="main-header">
      <h1>🧬 qPCR MIQE Statistical Analyzer</h1>
      <p>MIQE-compliant automated pipeline · Steps 1–8 · Publication-quality figures · PDF report</p>
    </div>
    """, unsafe_allow_html=True)

    # ── Upload & Settings (expander — avoids main_app sidebar conflict) ──
    with st.expander("📁 Data Upload & Settings", expanded=True):
        col_up, col_cfg = st.columns([2, 1])
        with col_up:
            uploaded = st.file_uploader(
                "Upload Excel / CSV file", type=['xlsx','xls','csv'],
                key="qpcr_uploader",
                help="Long-format: Sample_Name, Gene, CT, Replicate, [Condition1], [Condition2]")
            st.caption(
                "**Required columns:** Sample_Name · Gene · CT · Replicate  |  "
                "**Optional columns:** Condition1 · Condition2  |  "
                "**Reference genes auto-detected:** BACTIN, ACTB, GAPDH, 18S, HPRT1, RPLP0, TBP"
            )
        with col_cfg:
            alpha = st.slider("Significance level (α)", 0.01, 0.10, 0.05, 0.01, key="qpcr_alpha")
            ct_range_thresh = st.slider("Ct range flag threshold (cycles)",
                                        0.2, 1.0, 0.5, 0.1, key="qpcr_ct")
            sheet_name = st.text_input("Sheet name (Excel)", value="",
                                       placeholder="Leave blank for first sheet",
                                       key="qpcr_sheet")

    # ── Main Area ─────────────────────────────────────────────────
    if uploaded is None:
        col1, col2, col3 = st.columns(3)
        with col1:
            st.markdown("""
            <div class="step-card">
            <h3>📤 Step 1 · Upload</h3>
            Upload your qPCR data in long format (Excel or CSV). Reference genes are auto-detected.
            </div>""", unsafe_allow_html=True)
        with col2:
            st.markdown("""
            <div class="step-card">
            <h3>⚙️ Steps 2–7 · Auto-Analysis</h3>
            Preprocessing → Statistics → QC → Tables → Figures all run automatically on upload.
            </div>""", unsafe_allow_html=True)
        with col3:
            st.markdown("""
            <div class="step-card">
            <h3>📄 Step 8 · PDF Export</h3>
            Download a publication-ready PDF with methods, results, and all figures.
            </div>""", unsafe_allow_html=True)

        st.info("👈 Upload your Excel/CSV file in the sidebar to begin analysis.")
        st.stop()

    # ── Load Data ─────────────────────────────────────────────────
    @st.cache_data(show_spinner=False)
    def load_data(file_bytes, fname, sheet):
        if fname.endswith('.csv'):
            return pd.read_csv(io.BytesIO(file_bytes))
        else:
            kw = {} if not sheet else {'sheet_name': sheet}
            return pd.read_excel(io.BytesIO(file_bytes), **kw)

    with st.spinner("Loading data…"):
        try:
            raw_df = load_data(uploaded.read(), uploaded.name, sheet_name if sheet_name else None)
            if isinstance(raw_df, dict):
                snames = list(raw_df.keys())
                chosen = st.selectbox("Multiple sheets detected — select one:", snames)
                raw_df = raw_df[chosen]
        except Exception as e:
            st.error(f"Failed to load file: {e}"); st.stop()

    # ── Run Analysis ──────────────────────────────────────────────
    with st.spinner("🔬 Running MIQE pipeline (Steps 1–8)…"):
        res = run_analysis(raw_df.copy())

    if res['errors']:
        st.markdown('<div class="alert-danger"><b>🚨 STOP: Validation Failed</b><br>' +
                    '<br>'.join(res['errors']) + '</div>', unsafe_allow_html=True)
        with st.expander("📋 Raw Data Preview"):
            st.dataframe(raw_df.head(30), use_container_width=True)
        st.stop()

    if res['warnings']:
        for w in res['warnings']:
            st.markdown(f'<div class="alert-warning">⚠️ {w}</div>', unsafe_allow_html=True)

    s1 = res['steps']['step1']
    s2 = res['steps']['step2']
    s3 = res['steps']['step3']
    s6 = res['steps']['step6']
    merged = s2['merged']
    tgt_genes = s2['tgt_genes']
    ref_genes = s2['ref_genes']
    has_c1 = s2['has_c1']; has_c2 = s2['has_c2']

    # ── Summary KPIs ──────────────────────────────────────────────
    n_sig = sum(1 for g in tgt_genes if s3[g].get('sig_fdr', False))
    mode_label = 'Two-Way ANOVA' if has_c2 else ('Group Comparison' if has_c1 else 'Descriptive')
    mode_icon  = '🔬' if has_c2 else ('📊' if has_c1 else '📋')
    sig_color  = '#27ae60' if n_sig > 0 else '#7f8c8d'
    sig_bg     = 'linear-gradient(135deg,#d4edda,#a8d5b5)' if n_sig > 0 else 'linear-gradient(135deg,#f0f0f0,#dcdcdc)'
    sig_border = '#27ae60' if n_sig > 0 else '#aaaaaa'

    st.markdown("""
<div style="display:flex;align-items:center;gap:0.6rem;margin-bottom:1rem;">
  <span style="font-size:1.4rem;">📊</span>
  <span style="font-size:1.2rem;font-weight:700;color:#1a3a5c;letter-spacing:-0.3px;">Analysis Overview</span>
</div>
""", unsafe_allow_html=True)

    k1, k2, k3, k4, k5 = st.columns(5)

    _tgt_label = ', '.join(tgt_genes[:3])
    _factor_label = '2 factors' if has_c2 else ('1 factor' if has_c1 else 'no factor')
    _sig_label = '✓ significant' if n_sig > 0 else 'none significant'

    # Fixed height via flexbox so all 5 cards are identical in size
    _card_h = "130px"
    _card_inner = "display:flex;flex-direction:column;justify-content:center;align-items:center;"

    k1.markdown(f"""
<div style="background:linear-gradient(135deg,#1a3a5c,#2c5f8a);border-radius:14px;
            height:{_card_h};{_card_inner}text-align:center;color:white;
            box-shadow:0 4px 18px rgba(26,58,92,0.28);padding:0 0.8rem;">
  <div style="font-size:0.65rem;font-weight:600;letter-spacing:1.2px;opacity:0.75;
              text-transform:uppercase;margin-bottom:0.35rem;">Total Observations</div>
  <div style="font-size:2.2rem;font-weight:800;line-height:1.1;">{len(raw_df)}</div>
  <div style="font-size:0.65rem;opacity:0.6;margin-top:0.25rem;">rows loaded</div>
</div>""", unsafe_allow_html=True)

    k2.markdown(f"""
<div style="background:linear-gradient(135deg,#1565c0,#1976d2);border-radius:14px;
            height:{_card_h};{_card_inner}text-align:center;color:white;
            box-shadow:0 4px 18px rgba(21,101,192,0.28);padding:0 0.8rem;">
  <div style="font-size:0.65rem;font-weight:600;letter-spacing:1.2px;opacity:0.75;
              text-transform:uppercase;margin-bottom:0.35rem;">Biological Samples</div>
  <div style="font-size:2.2rem;font-weight:800;line-height:1.1;">{s1['n_samples']}</div>
  <div style="font-size:0.65rem;opacity:0.6;margin-top:0.25rem;">unique samples</div>
</div>""", unsafe_allow_html=True)

    k3.markdown(f"""
<div style="background:linear-gradient(135deg,#6a1b9a,#8e24aa);border-radius:14px;
            height:{_card_h};{_card_inner}text-align:center;color:white;
            box-shadow:0 4px 18px rgba(106,27,154,0.28);padding:0 0.8rem;">
  <div style="font-size:0.65rem;font-weight:600;letter-spacing:1.2px;opacity:0.75;
              text-transform:uppercase;margin-bottom:0.35rem;">Target Genes</div>
  <div style="font-size:2.2rem;font-weight:800;line-height:1.1;">{len(tgt_genes)}</div>
  <div style="font-size:0.65rem;opacity:0.6;margin-top:0.25rem;">{_tgt_label}</div>
</div>""", unsafe_allow_html=True)

    k4.markdown(f"""
<div style="background:linear-gradient(135deg,#e65100,#ef6c00);border-radius:14px;
            height:{_card_h};{_card_inner}text-align:center;color:white;
            box-shadow:0 4px 18px rgba(230,81,0,0.28);padding:0 0.8rem;">
  <div style="font-size:0.65rem;font-weight:600;letter-spacing:1.2px;opacity:0.75;
              text-transform:uppercase;margin-bottom:0.35rem;">Analysis Mode</div>
  <div style="font-size:1.1rem;font-weight:800;line-height:1.3;margin:0.2rem 0;">
    {mode_icon} {mode_label}</div>
  <div style="font-size:0.65rem;opacity:0.6;margin-top:0.25rem;">{_factor_label}</div>
</div>""", unsafe_allow_html=True)

    k5.markdown(f"""
<div style="background:{sig_bg};border:2px solid {sig_border};border-radius:14px;
            height:{_card_h};{_card_inner}text-align:center;
            box-shadow:0 4px 18px rgba(0,0,0,0.10);padding:0 0.8rem;">
  <div style="font-size:0.65rem;font-weight:600;letter-spacing:1.2px;color:{sig_color};
              text-transform:uppercase;margin-bottom:0.35rem;">FDR Significant</div>
  <div style="font-size:2.2rem;font-weight:800;line-height:1.1;color:{sig_color};">
    {n_sig}<span style="font-size:1.1rem;font-weight:500;color:#555;">/{len(tgt_genes)}</span></div>
  <div style="font-size:0.65rem;color:{sig_color};opacity:0.8;margin-top:0.25rem;">{_sig_label}</div>
</div>""", unsafe_allow_html=True)

    st.markdown("<div style='margin-top:1.4rem;border-top:2px solid #e8eef5;'></div>",
                unsafe_allow_html=True)

    # ── TABS ──────────────────────────────────────────────────────
    tab_raw, tab_qc, tab_stats, tab_figs, tab_text, tab_dl = st.tabs([
        "📋 Data", "🔍 QC", "📈 Statistics", "📊 Figures", "📝 Manuscript", "⬇️ Download PDF"
    ])

    # ─── TAB 1: Raw Data ─────────────────────────────────────────
    with tab_raw:
        st.markdown('<div class="step-card"><h3>Step 1 · Data Structure Detection</h3></div>',
                    unsafe_allow_html=True)
        c1,c2,c3 = st.columns(3)
        with c1:
            st.markdown(f'<div class="alert-success">✓ Reference genes: <b>{", ".join(ref_genes)}</b></div>',
                        unsafe_allow_html=True)
        with c2:
            st.markdown(f'<div class="alert-success">✓ Target genes: <b>{", ".join(tgt_genes)}</b></div>',
                        unsafe_allow_html=True)
        with c3:
            cond_str = ' × '.join([f"C1({len(s1['conditions']['Condition1'])})",
                                    f"C2({len(s1['conditions']['Condition2'])})"] if has_c2
                                   else ([f"C1: {s1['conditions']['Condition1']}"] if has_c1 else ['No conditions']))
            st.markdown(f'<div class="alert-info">Conditions: {cond_str}</div>',
                        unsafe_allow_html=True)

        st.markdown("**Raw Input Data**")
        st.dataframe(raw_df, use_container_width=True, height=320)
        st.markdown("**Processed Data (ΔCt · ΔΔCt · Fold Change)**")
        show_cols = [c for c in ['Sample_Name','Gene','Condition1','Condition2',
                                   'Mean_CT','Ref_CT','Delta_CT','Delta_Delta_CT','Fold_Change','Log2_FC']
                     if c in merged.columns]
        st.dataframe(merged[show_cols].round(4), use_container_width=True, height=300)

    # ─── TAB 2: QC ───────────────────────────────────────────────
    with tab_qc:
        st.markdown('<div class="step-card"><h3>Steps 2 & 5 · Preprocessing QC</h3></div>',
                    unsafe_allow_html=True)

        flagged = s2['flagged']
        iqr_flags = s2['iqr_flags']

        c1,c2 = st.columns(2)
        with c1:
            st.markdown("#### Technical Replicate QC (Ct range > 0.5)")
            if len(flagged):
                st.markdown(f'<div class="alert-warning">⚠️ <b>{len(flagged)}</b> flagged observation(s)</div>',
                            unsafe_allow_html=True)
                st.dataframe(flagged.round(3), use_container_width=True)
            else:
                st.markdown('<div class="alert-success">✓ All technical replicates pass QC</div>',
                            unsafe_allow_html=True)
        with c2:
            st.markdown("#### ΔCt IQR Outlier Check")
            if iqr_flags:
                st.markdown(f'<div class="alert-warning">⚠️ <b>{len(iqr_flags)}</b> potential outlier(s)</div>',
                            unsafe_allow_html=True)
                st.dataframe(pd.DataFrame(iqr_flags), use_container_width=True)
            else:
                st.markdown('<div class="alert-success">✓ No IQR outliers detected</div>',
                            unsafe_allow_html=True)

        st.markdown("#### Complete Ct QC Table")
        mc_show = s2['mc'].copy()
        mc_show['Status'] = mc_show['QC_Flag'].map({True:'⚠️ FLAGGED', False:'✓ OK'})
        st.dataframe(mc_show.round(3).drop(columns=['QC_Flag']), use_container_width=True, height=300)

    # ─── TAB 3: Statistics ───────────────────────────────────────
    with tab_stats:
        st.markdown('<div class="step-card"><h3>Steps 3 & 4 · Statistical Testing + FDR Correction</h3></div>',
                    unsafe_allow_html=True)

        # Descriptive table
        st.markdown("#### Step 6 · Descriptive Statistics")
        st.dataframe(s6['desc_df'].round(4), use_container_width=True)

        st.markdown("---")
        st.markdown("#### Inferential Statistics")

        for gene in tgt_genes:
            sr = s3[gene]
            with st.expander(f"🧬 Gene: **{gene}** — {sr['mode'].replace('_',' ').title()}", expanded=True):
                if sr['mode'] == 'two_way':
                    # Normality
                    nc = sr.get('norm_cells',{})
                    if nc:
                        st.markdown("**Normality (Shapiro-Wilk per cell):**")
                        nd = pd.DataFrame([
                            {'Cell': cell, 'n': v['n'], 'W': round(v['W'],4),
                             'p': round(v['p'],4), 'Normal': '✓' if v['normal'] else '✗'}
                            for cell, v in nc.items()])
                        st.dataframe(nd, use_container_width=True, hide_index=True)

                    st.markdown(f"**Levene's test:** stat={sr['levene'][0]:.4f}, p={sr['levene'][1]:.4f} → "
                                f"{'Equal variances ✓' if sr['levene'][1]>0.05 else 'Unequal variances ✗'}")

                    # ANOVA table
                    a2 = sr['anova_2way']
                    anova_rows = []
                    for i, src in enumerate(a2['Source']):
                        anova_rows.append({
                            'Source': src,
                            'SS': round(a2['SS'][i],4), 'df': int(a2['df'][i]),
                            'MS': round(a2['MS'][i],4),
                            'F': round(a2['F'][i],4) if a2['F'][i] else '—',
                            'p': round(a2['p'][i],4) if a2['p'][i] is not None else '—',
                            'Partial η²': round(a2['eta2'][i],4) if a2['eta2'][i] is not None else '—',
                            'Significant': sig_stars(a2['p'][i]) if a2['p'][i] is not None else '',
                        })
                    atbl_df = pd.DataFrame(anova_rows)
                    st.markdown("**Two-Way ANOVA Table:**")

                    def color_sig(val):
                        if val in ('*','**','***'): return 'background-color: #d4edda; color: #155724; font-weight:bold'
                        return ''
                    st.dataframe(atbl_df.style.applymap(color_sig, subset=['Significant']),
                                 use_container_width=True, hide_index=True)

                    # Post-hoc
                    ph = sr.get('posthoc',[])
                    if ph:
                        st.markdown("**Post-hoc Pairwise Comparisons (Holm correction):**")
                        ph_df = pd.DataFrame([{
                            'Comparison': f"{c['g1']} vs {c['g2']}",
                            't': round(c['t'],3), "Cohen's d": round(abs(c['d']),3),
                            'Raw p': round(c['p'],4), 'Holm p': round(c['p_holm'],4),
                            'Sig': '✓' if c['sig'] else 'ns',
                        } for c in ph])
                        st.dataframe(ph_df, use_container_width=True, hide_index=True)

                    c1,c2,c3 = st.columns(3)
                    c1.metric("Treatment p", f"{sr['p_main']:.4f}", delta=sig_stars(sr['p_main']))
                    c2.metric("Time p", f"{sr['p_time']:.4f}", delta=sig_stars(sr['p_time']))
                    c3.metric("Interaction p", f"{sr['p_int']:.4f}", delta=sig_stars(sr['p_int']))

                elif sr['mode'] == 'two_group':
                    if sr.get('normality'):
                        nm_df = pd.DataFrame([
                            {'Group': g, 'W': round(v['W'],4), 'p': round(v['p'],4),
                             'Normal': '✓' if v['normal'] else '✗'}
                            for g,v in sr['normality'].items()])
                        st.markdown("**Normality (Shapiro-Wilk):**")
                        st.dataframe(nm_df, use_container_width=True, hide_index=True)

                    c1,c2,c3,c4 = st.columns(4)
                    c1.metric("Test", sr['test'])
                    c2.metric("Statistic", f"{sr['stat']:.4f}")
                    c3.metric("p-value", f"{sr['p']:.4f}", delta=sig_stars(sr['p']))
                    c4.metric(sr['effect_name'], f"{abs(sr['effect']):.4f}")

                    sig_badge = '<span class="sig-badge-yes">✓ Significant</span>' if sr['sig_fdr'] else '<span class="sig-badge-no">ns</span>'
                    st.markdown(f"BH FDR q-value: **{sr['q']:.4f}** → {sig_badge}", unsafe_allow_html=True)

                elif sr['mode'] == 'multi_group':
                    c1,c2,c3 = st.columns(3)
                    c1.metric("Test", sr['test'])
                    c2.metric("Statistic", f"{sr['stat']:.4f}")
                    c3.metric("p-value", f"{sr['p']:.4f}", delta=sig_stars(sr['p']))
                    sig_badge = '<span class="sig-badge-yes">✓ Significant</span>' if sr['sig_fdr'] else '<span class="sig-badge-no">ns</span>'
                    st.markdown(f"BH FDR q-value: **{sr['q']:.4f}** → {sig_badge}", unsafe_allow_html=True)

                else:
                    st.info("Descriptive mode — no inferential test applied (no conditions detected).")

        # Global FDR table
        st.markdown("---")
        st.markdown("#### Step 4 · BH FDR Global Summary")
        fdr_rows = []
        for gene in tgt_genes:
            sr = s3[gene]
            p_v = sr.get('p', sr.get('p_main','—'))
            fdr_rows.append({'Gene': gene,
                'Test': sr.get('test', sr['mode']),
                'p-value': round(p_v,4) if isinstance(p_v,float) else p_v,
                'BH q-value': round(sr['q'],4),
                'Significant': '✓ YES' if sr['sig_fdr'] else 'NO'})
        fdr_df = pd.DataFrame(fdr_rows)
        def fdr_color(val):
            if val == '✓ YES': return 'background-color:#d4edda;color:#155724;font-weight:bold'
            return ''
        st.dataframe(fdr_df.style.applymap(fdr_color, subset=['Significant']),
                     use_container_width=True, hide_index=True)

    # ─── TAB 4: Figures ──────────────────────────────────────────
    with tab_figs:
        st.markdown('<div class="step-card"><h3>Step 7 · Publication-Quality Figures</h3></div>',
                    unsafe_allow_html=True)

        with st.spinner("Generating figures…"):
            fig_expr = make_expression_plot(merged, s3, tgt_genes, has_c1, has_c2,
                                             s2['ctrl_c1'], s2['ctrl_c2'])
            fig_heat = make_heatmap(merged, tgt_genes, has_c1, has_c2)
            fig_volc = make_volcano(merged, s3, tgt_genes)
            fig_inter = make_interaction_plot(merged, tgt_genes) if has_c2 else None

        fig_bytes_dict = {
            'expr': fig_to_bytes(fig_expr),
            'heat': fig_to_bytes(fig_heat),
            'volc': fig_to_bytes(fig_volc),
        }
        if fig_inter:
            fig_bytes_dict['inter'] = fig_to_bytes(fig_inter)

        # Display
        st.markdown("#### Figure 1 · Expression Plot (log₂ FC ± SEM)")
        st.pyplot(fig_expr, use_container_width=True)
        st.caption("*Relative expression calculated using the 2⁻ᐩᐩCt method. Statistics on ΔCt values. BH FDR corrected.*")

        if fig_inter:
            st.markdown("#### Figure 2 · Interaction Plot")
            st.pyplot(fig_inter, use_container_width=True)
            st.caption("*ΔCt mean ± SEM. Parallel lines suggest no Treatment × Time interaction.*")

        col_h, col_v = st.columns(2)
        with col_h:
            st.markdown("#### Figure 3 · Heatmap")
            st.pyplot(fig_heat, use_container_width=True)
            st.caption("*log₂ FC per sample. Blue=down · White=neutral · Red=up*")
        with col_v:
            st.markdown("#### Figure 4 · Volcano Plot")
            st.pyplot(fig_volc, use_container_width=True)
            st.caption("*Red = FDR q < 0.05 · Dashed = p=0.05*")

    # ─── TAB 5: Manuscript Text ──────────────────────────────────
    with tab_text:
        st.markdown('<div class="step-card"><h3>Step 8 · Manuscript-Ready Text</h3></div>',
                    unsafe_allow_html=True)

        # Methods
        method_text = f"""**qPCR Normalization & Statistical Methods**

    Relative mRNA expression was quantified by reverse transcription-quantitative PCR (RT-qPCR) 
    following MIQE guidelines. {", ".join(ref_genes)} served as the reference gene for normalization. 
    Relative expression was calculated using the 2⁻ᐩᐩCt method, where 
    ΔCt = Ct(target) − Ct(reference) and ΔΔCt = ΔCt(sample) − mean ΔCt(control group). 
    **All inferential statistical analyses were performed on ΔCt values.**
    """
        if has_c2:
            method_text += f"""
    A 2×2 factorial two-way ANOVA was applied to ΔCt values with Treatment (Condition1) and 
    Time (Condition2) as independent factors. Within-cell normality was assessed using the 
    Shapiro-Wilk test and homogeneity of variance by Levene's test. Post-hoc pairwise 
    comparisons applied Welch's t-tests with Holm step-down correction. Effect sizes are 
    reported as partial η² (ANOVA) and Cohen's d (pairwise)."""
        elif has_c1:
            method_text += f"""
    Normality was assessed by the Shapiro-Wilk test and variance equality by Levene's test. 
    For normally distributed data with equal variances, Student's t-test was used; 
    unequal variances prompted Welch's t-test; non-normality prompted Mann-Whitney U. 
    Effect sizes are reported as Cohen's d or rank-biserial correlation."""
        if len(tgt_genes) > 1:
            method_text += """
    Benjamini-Hochberg FDR correction was applied across all target genes to control the 
    false discovery rate. Significance was defined as FDR q < 0.05."""
        else:
            method_text += "\nSignificance was defined as p < 0.05."
        method_text += """
    Technical replicates with intra-replicate Ct range > 0.5 cycles were flagged but retained."""

        st.markdown('<div class="manuscript-box">' + method_text.replace('\n','<br>') + '</div>',
                    unsafe_allow_html=True)

        st.markdown("---")
        st.markdown("**Results Summary**")

        results_text = ""
        for gene in tgt_genes:
            sr = s3[gene]
            gdata = merged[merged['Gene']==gene]
            if sr['mode'] == 'two_group':
                g2 = sr['g2']
                fc_m = gdata[gdata['Condition1']==g2]['Fold_Change'].mean()
                fc_s = gdata[gdata['Condition1']==g2]['Fold_Change'].sem()
                direction = 'upregulated' if fc_m > 1 else 'downregulated'
                results_text += (
                    f"**{gene}** mRNA was significantly {direction} in the {g2} group "
                    f"(mean FC = {fc_m:.2f} ± {fc_s:.2f} SEM; {sr['test']}: "
                    f"stat = {sr['stat']:.3f}, p = {sr['p']:.4f}, "
                    f"{sr['effect_name']} = {abs(sr['effect']):.3f}; BH q = {sr['q']:.4f}).\n\n"
                )
            elif sr['mode'] == 'two_way':
                sig_str = 'significant' if sr['p_main'] < 0.05 else 'not significant'
                results_text += (
                    f"**{gene}**: Two-way ANOVA revealed a {sig_str} main effect of Treatment "
                    f"(F(1,{sr['dfW']}) = {sr['F_main']:.3f}, p = {sr['p_main']:.4f}, "
                    f"partial η² = {sr['eta2_main']:.3f}). "
                    f"Time main effect: p = {sr['p_time']:.4f}. "
                    f"Interaction: p = {sr['p_int']:.4f} (BH q = {sr['q']:.4f}).\n\n"
                )
            else:
                desc = s6['desc_df'][s6['desc_df']['Gene']==gene]
                results_text += (
                    f"**{gene}**: Descriptive analysis. "
                    f"Mean FC = {desc['Mean_FC'].mean():.2f}, "
                    f"Mean ΔCt = {desc['Mean_dCt'].mean():.2f}.\n\n"
                )

        st.markdown('<div class="manuscript-box">' +
                    results_text.replace('\n\n','<br><br>').replace('\n','<br>') +
                    '</div>', unsafe_allow_html=True)

    # ─── TAB 6: Download ─────────────────────────────────────────
    with tab_dl:
        st.markdown('<div class="step-card"><h3>📄 Download PDF Report</h3></div>',
                    unsafe_allow_html=True)

        st.markdown("""
        The PDF report includes all 8 steps:
        - **Title page** with analysis summary
        - **QC summary** (technical replicates + IQR outliers)
        - **Descriptive statistics table**
        - **Inferential statistics tables** (normality, ANOVA / t-test, post-hoc)
        - **BH FDR global summary**
        - **All 4 figures** (expression, interaction, heatmap, volcano)
        - **Manuscript-ready Methods & Results text**
        """)

        st.markdown("---")

        # Pre-generate figures if not already cached
        if 'fig_bytes_dict' not in dir():
            with st.spinner("Preparing figures for PDF…"):
                fig_expr = make_expression_plot(merged, s3, tgt_genes, has_c1, has_c2, s2['ctrl_c1'], s2['ctrl_c2'])
                fig_heat = make_heatmap(merged, tgt_genes, has_c1, has_c2)
                fig_volc = make_volcano(merged, s3, tgt_genes)
                fig_inter = make_interaction_plot(merged, tgt_genes) if has_c2 else None
                fig_bytes_dict = {
                    'expr': fig_to_bytes(fig_expr),
                    'heat': fig_to_bytes(fig_heat),
                    'volc': fig_to_bytes(fig_volc),
                }
                if fig_inter:
                    fig_bytes_dict['inter'] = fig_to_bytes(fig_inter)

        with st.spinner("Building PDF…"):
            pdf_bytes = build_pdf(res, fig_bytes_dict)

        fname_gene = tgt_genes[0] if len(tgt_genes)==1 else 'MultiGene'
        pdf_fname = f"qPCR_Report_{fname_gene}_{date.today().strftime('%Y%m%d')}.pdf"

        st.download_button(
            label="⬇️  Download PDF Report",
            data=pdf_bytes,
            file_name=pdf_fname,
            mime="application/pdf",
            use_container_width=True,
            type="primary",
        )

        st.success(f"✅ Report ready: **{pdf_fname}**  ({len(pdf_bytes)//1024} KB)")

        st.markdown("---")
        st.markdown("#### Preview: Report Contents")
        prev_cols = st.columns(2)
        with prev_cols[0]:
            st.markdown("""
            **Pages 1–3: Data & Statistics**
            - Analysis summary table
            - Technical replicate QC
            - Descriptive statistics
            - Normality & variance tests
            - Two-Way ANOVA / t-test / Mann-Whitney
            - Post-hoc comparisons (Holm)
            - BH FDR global table
            """)
        with prev_cols[1]:
            st.markdown("""
            **Pages 4–6: Figures & Text**
            - Figure 1: Expression plot (log₂FC)
            - Figure 2: Interaction plot (if 2-factor)
            - Figure 3: Sample heatmap
            - Figure 4: Volcano plot
            - Methods section (manuscript-ready)
            - Results section (manuscript-ready)
            """)
