"""
app.py — HCV Cross-Neutralization Dashboard (Streamlit, live)
=============================================================
Interactive, filterable rebuild of the V1–V16 heatmaps in one app.

Run locally:
    pip install -r requirements.txt
    streamlit run app.py

Deploy:
    push this folder to GitHub -> share.streamlit.io -> add your service-account
    JSON under Secrets (see .streamlit/secrets.toml.example).

Data sources (pick in the sidebar):
    1. Google Sheet (service account)  -> always-live, needs st.secrets["gcp_service_account"]
    2. Published CSV URL               -> File > Share > Publish to web > CSV
    3. Upload CSV / XLSX               -> quick offline demo
"""

from __future__ import annotations

import io
import pandas as pd
import streamlit as st

import hcv_data as H
import hcv_viz as V

st.set_page_config(page_title="HCV Cross-Neutralization Dashboard",
                   page_icon="🧬", layout="wide")

# Default sheet IDs (from the V16 pipeline) — editable in the sidebar.
DEFAULT_IC50_SHEET = "1pGonKQsnbD4E_-ywm8-XjX_cG6eW0yF5-5hRgHDVWRc"
DEFAULT_CONSTRUCT_SHEET = "1iv9fbnKjvWt_LCKuNPJy58ldJJGfey0SoJSv69cJF3E"


# ============================================================
# DATA LOADING
# ============================================================

@st.cache_data(ttl=600, show_spinner=False)
def load_from_gsheet(sheet_id: str, worksheet: str | None) -> pd.DataFrame:
    import gspread
    if "gcp_service_account" not in st.secrets:
        raise RuntimeError(
            "No service-account credentials found. Add your JSON under "
            "st.secrets['gcp_service_account'] (see .streamlit/secrets.toml.example), "
            "or use the 'Published CSV URL' / 'Upload file' source instead.")
    gc = gspread.service_account_from_dict(dict(st.secrets["gcp_service_account"]))
    sh = gc.open_by_key(sheet_id)
    ws = sh.worksheet(worksheet) if worksheet else sh.sheet1
    values = ws.get_all_values()
    if not values:
        return pd.DataFrame()
    header = [c.strip() for c in values[0]]
    body = [row + [""] * (len(header) - len(row)) for row in values[1:]]
    return pd.DataFrame(body, columns=header)


@st.cache_data(ttl=600, show_spinner=False)
def load_from_csv_url(url: str) -> pd.DataFrame:
    import requests
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return pd.read_csv(io.StringIO(r.text))


@st.cache_data(show_spinner=False)
def load_from_upload(name: str, data: bytes) -> pd.DataFrame:
    if name.lower().endswith((".xlsx", ".xls")):
        return pd.read_excel(io.BytesIO(data))
    return pd.read_csv(io.BytesIO(data))


def sidebar_data_source() -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
    """Returns (df_ic50, df_constructs). df_constructs may be None."""
    st.sidebar.header("1 · Data source")
    source = st.sidebar.radio(
        "Where to read the IC50 sheet from",
        ["Google Sheet (live)", "Published CSV URL", "Upload CSV / XLSX"],
        help="Live reads need a service account; the other two work without one.")

    if st.sidebar.button("🔄 Refresh data (clear cache)"):
        st.cache_data.clear()
        st.rerun()

    df_ic50 = None
    df_const = None

    try:
        if source == "Google Sheet (live)":
            sid = st.sidebar.text_input("IC50 sheet ID", DEFAULT_IC50_SHEET)
            ws  = st.sidebar.text_input("Worksheet name (blank = first)", "")
            cid = st.sidebar.text_input("Constructs sheet ID (optional)",
                                        DEFAULT_CONSTRUCT_SHEET,
                                        help="HCV_Constructs_List sheet — "
                                             "used to resolve full construct names.")
            if sid:
                df_ic50 = load_from_gsheet(sid, ws.strip() or None)
            if cid:
                try:
                    df_const = load_from_gsheet(cid, None)
                except Exception:
                    pass  # constructs sheet is optional

        elif source == "Published CSV URL":
            url = st.sidebar.text_input("Published CSV URL", "")
            if not url:
                st.sidebar.info("Paste a 'Publish to web → CSV' link to load.")
            else:
                df_ic50 = load_from_csv_url(url)

        else:
            up = st.sidebar.file_uploader("Upload the IC50 export",
                                          type=["csv", "xlsx", "xls"])
            if up is not None:
                df_ic50 = load_from_upload(up.name, up.getvalue())

    except Exception as exc:
        st.sidebar.error(f"Could not load data: {exc}")

    return df_ic50, df_const


# ============================================================
# RENDER HELPERS
# ============================================================

def render_heatmap_with_selection(fig, key: str):
    """Render a Plotly heatmap; return a click selection dict if supported."""
    try:
        return st.plotly_chart(fig, use_container_width=True, key=key,
                               on_select="rerun", selection_mode="points")
    except TypeError:
        st.plotly_chart(fig, use_container_width=True, key=key)
        return None


def selection_to_labels(event, value_pivot):
    """Best-effort map a heatmap click back to (construct, psv) raw labels."""
    try:
        pts = event["selection"]["points"]
        if not pts:
            return None, None
        pt = pts[0]
        constructs, psvs = list(value_pivot.index), list(value_pivot.columns)
        # Prefer positional indices (robust to label wrapping); fall back to labels.
        if pt.get("point_indices"):
            r, c = pt["point_indices"]
            disp = list(constructs)[::-1]  # figure reverses y order
            return disp[r], psvs[c]
        return pt.get("y"), pt.get("x")
    except Exception:
        return None, None


def legend_caption(metric: str, mode: str):
    if mode == "threshold":
        st.caption("🟩 ≥ threshold (hit)  ·  ⬜ tested, below threshold  ·  "
                   "✕ not tested  ·  ⬛ No Neutralization")
    elif metric == "log10_ic50":
        st.caption("Color = log₁₀(IC50): red (low) → green (high potency)  ·  "
                   "⬛ No Neutralization  ·  ✕ not tested")
    else:
        st.caption("Color = % neutralization: white → green  ·  ✕ not tested")


def download_view(value_pivot, status_pivot, label: str):
    if value_pivot.empty:
        return
    # Ensure index is named for reset_index to work correctly
    vp = value_pivot.copy()
    sp = status_pivot.copy()
    if vp.index.name is None:
        vp.index.name = "Construct_Description"
    if sp.index.name is None:
        sp.index.name = "Construct_Description"
    
    long = (vp.reset_index()
            .melt(id_vars="Construct_Description", var_name="PSV", value_name="value"))
    stat = (sp.reset_index()
            .melt(id_vars="Construct_Description", var_name="PSV", value_name="status"))
    merged = long.merge(stat, on=["Construct_Description", "PSV"])
    st.download_button("⬇️ Download this view (CSV)",
                       merged.to_csv(index=False).encode(),
                       file_name=f"hcv_view_{label}.csv", mime="text/csv")

# ============================================================
# MAIN
# ============================================================

st.title("🧬 HCV Cross-Neutralization Dashboard")
st.markdown("Constructs × pseudoviruses — every heatmap view, filterable in one place.")

df_raw, df_const_raw = sidebar_data_source()
if df_raw is None or df_raw.empty:
    st.info("⬅️ Choose a data source in the sidebar to begin. "
            "If credentials aren't set up yet, try **Published CSV URL** or **Upload**.")
    st.stop()

construct_lookup = H.build_construct_lookup(df_const_raw) if df_const_raw is not None else {}

ic50_source = st.sidebar.radio(
    "IC50 values",
    ["Background Corrected", "Not Background Corrected"],
    help="Background Corrected = column IC50_corrected (col G)  |  "
         "Not Background Corrected = column IC50 (col F)")
corrected_ic50 = (ic50_source == "Background Corrected")

tidy, info = H.prepare_dataframe(df_raw, construct_lookup=construct_lookup or None,
                                 corrected_ic50=corrected_ic50)
if tidy.empty:
    st.error("No rows survived filtering. Check that the sheet has the expected "
             "columns (Experiment, Group, PSV, Day, PSVX_No, IC50, Dilution, "
             "Avg_Neut_percent_corrected) and HCV experiments.")
    st.stop()

# ---------------- filters ----------------
st.sidebar.header("2 · View")
metric_label = st.sidebar.radio("Metric", ["% neutralization", "log₁₀(IC50)"])
metric = "pct_neut" if metric_label.startswith("%") else "log10_ic50"

dilution = None
if metric == "pct_neut":
    dils = info["all_dilutions"] or [30.0]
    default_idx = next((i for i, d in enumerate(dils) if d == 90.0),
                       next((i for i, d in enumerate(dils) if d == 30.0), 0))
    dilution = st.sidebar.selectbox(
        "Dilution (1:x)", dils, index=default_idx,
        format_func=lambda d: f"1:{int(d)}",
        help="Pick a single dilution point. To compare potency across all dilutions "
             "use log₁₀(IC50) mode instead.")

display_mode = st.sidebar.radio("Cell encoding", ["Gradient", "Threshold (hit map)"])
mode = "threshold" if display_mode.startswith("Threshold") else "gradient"

threshold, ge = (50.0 if metric == "pct_neut" else 3.0), True
if mode == "threshold":
    presets = (["≥50%", "≥75%", "Custom"] if metric == "pct_neut"
               else ["≥3.0", "≥3.5", "Custom"])
    choice = st.sidebar.radio("Threshold", presets, horizontal=True)
    if choice == "Custom":
        threshold = st.sidebar.number_input(
            "Custom threshold", value=float(threshold),
            step=5.0 if metric == "pct_neut" else 0.5)
    else:
        threshold = float(choice.replace("≥", "").replace("%", ""))
    ge = st.sidebar.checkbox("Use ≥ (uncheck for strictly >)", value=True)

view_mode = st.sidebar.radio("Layout", ["Single heatmap", "Small multiples (by subgroup)"])

buckets_present = [b for b in H.BUCKETS if b in set(tidy["Bucket_Type"].unique())]
bucket_opts = buckets_present + (["All (pooled)"] if buckets_present else [])
bucket_choice = st.sidebar.selectbox("Dose window", bucket_opts or ["All (pooled)"])
bucket = None if bucket_choice == "All (pooled)" else bucket_choice

subgroup = "All constructs"
if view_mode == "Single heatmap":
    subgroup = st.sidebar.selectbox("Construct subgroup",
                                    ["All constructs"] + info["subgroups_present"])

with st.sidebar.expander("More filters"):
    exps = st.multiselect("Experiment", info["experiments"], default=info["experiments"])
    grps = st.multiselect("Group", info["groups"], default=info["groups"])
    psvs_sel = st.multiselect("PSV", info["psvs"], default=info["psvs"])
    show_values = st.checkbox("Show values in cells", value=True)
    use_geno = st.checkbox("Group PSV columns by genotype",
                           value=bool(H.PSV_GENOTYPE), disabled=not H.PSV_GENOTYPE)

# ---------------- apply row filters ----------------
f = tidy[tidy["Experiment"].isin(exps) & tidy["PSV"].isin(psvs_sel)]
if grps:
    f = f[f["Group"].isin(grps)]
if subgroup != "All constructs":
    f = f[f["Subgroup"] == subgroup]

view = H.compute_view(f, metric=metric, dilution=dilution)
psv_geno = H.PSV_GENOTYPE if (use_geno and H.PSV_GENOTYPE) else None
thr_pct = threshold if metric == "pct_neut" else 50.0

# ---------------- diagnostics ----------------
with st.expander("📋 Data summary & diagnostics"):
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Raw rows", info["n_raw"])
    c2.metric("After HCV filter", info["n_after_filter"])
    c3.metric("Unmapped bucket", info["n_unknown_bucket"])
    c4.metric("Uncategorized", info["n_uncategorized"])
    if info["n_unknown_bucket"]:
        st.warning(f"{info['n_unknown_bucket']} rows have no Prime/Boost1/Boost2 "
                   "mapping (excluded). Add their (experiment, PSVX, day) keys to "
                   "BUCKET_MAP in hcv_data.py.")
    if info["n_uncategorized"]:
        st.warning(f"{info['n_uncategorized']} rows fell into 'Uncategorized'. "
                   "Extend SUBGROUP_RULES in hcv_data.py to classify them.")
    st.write("Detected columns:", info["columns"])

# ============================================================
# SINGLE HEATMAP
# ============================================================
if view_mode == "Single heatmap":
    value_pivot, status_pivot, counts = H.build_pivots(
        f, view, bucket=bucket, metric=metric, mode=mode, threshold=threshold, ge=ge)

    if metric == "pct_neut":
        _dil_lbl = f"% neut @ 1:{int(dilution)}"
    else:
        _dil_lbl = "log₁₀(IC50)"
    title = (f"{subgroup}  ·  {bucket_choice}  ·  {_dil_lbl}"
             f"{'  ·  ≥' + str(threshold) if mode == 'threshold' else ''}")

    fig = V.build_heatmap_figure(value_pivot, status_pivot, counts, metric, mode,
                                 threshold, title=title, psv_genotype=psv_geno,
                                 show_values=show_values)
    event = render_heatmap_with_selection(fig, key="main_heatmap")
    legend_caption(metric, mode)
    download_view(value_pivot, status_pivot, "single")

    # ---- neutralization-curve detail ----
    st.subheader("🔬 Neutralization curve")
    sel_c, sel_p = selection_to_labels(event, value_pivot) if event else (None, None)
    if not value_pivot.empty:
        constructs = list(value_pivot.index)
        psvs_avail = list(value_pivot.columns)
        ci = constructs.index(sel_c) if sel_c in constructs else 0
        pi = psvs_avail.index(sel_p) if sel_p in psvs_avail else 0
        col_a, col_b = st.columns(2)
        cc = col_a.selectbox("Construct", constructs, index=ci)
        pp = col_b.selectbox("PSV", psvs_avail, index=pi)
        curve = H.get_curve(f, cc, pp, buckets=[bucket] if bucket else None)
        st.plotly_chart(V.build_curve_figure(curve, cc, pp, thr_pct),
                        use_container_width=True)
        if sel_c:
            st.caption("Tip: clicking a heatmap cell pre-selects the curve "
                       "(falls back to the dropdowns if your Streamlit build "
                       "doesn't emit click events).")

# ============================================================
# SMALL MULTIPLES (one heatmap per subgroup, fixed bucket)
# ============================================================
else:
    st.markdown(f"#### Subgroups — **{bucket_choice}** window")
    subs = info["subgroups_present"] or ["All constructs"]
    for sg in subs:
        fsg = f[f["Subgroup"] == sg] if sg != "All constructs" else f
        vsub = H.compute_view(fsg, metric=metric, dilution=dilution)
        vp, sp, cnt = H.build_pivots(vsub, bucket=bucket, metric=metric,
                                     mode=mode, threshold=threshold, ge=ge)
        if vp.empty:
            continue
        with st.expander(f"{sg}  ·  {len(vp)} constructs × {len(vp.columns)} PSVs",
                         expanded=True):
            fig = V.build_heatmap_figure(vp, sp, cnt, metric, mode, threshold,
                                         title="", psv_genotype=psv_geno,
                                         show_values=show_values)
            st.plotly_chart(fig, use_container_width=True, key=f"sm_{sg}")
            legend_caption(metric, mode)
