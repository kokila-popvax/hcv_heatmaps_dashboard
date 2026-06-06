# HCV Cross-Neutralization Dashboard

An interactive, filterable Streamlit rebuild of the V1–V16 heatmaps. One
constructs × pseudoviruses canvas that reproduces every view you've built —
log₁₀(IC50) gradients, % neutralization at any dilution, threshold "hit maps"
with breadth counts, and per-subgroup splits — all driven by live data from
your Google Sheet.

## Files

| File | Purpose |
|------|---------|
| `app.py` | Streamlit UI: data loading, filters, heatmap, curve, exports |
| `hcv_data.py` | Pure data logic (filtering, parsing, bucketing, subgroup classification). **Holds the two blocks you edit.** |
| `hcv_viz.py` | Plotly figure builders (heatmap + breadth bar, neutralization curve) |
| `requirements.txt` | Dependencies |
| `.streamlit/config.toml` | Theme |
| `.streamlit/secrets.toml.example` | Service-account credentials template |

## Run locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

Then pick a **data source** in the sidebar. To try it immediately without
credentials, use **Published CSV URL** (File → Share → Publish to web → CSV in
your sheet) or **Upload CSV / XLSX**.

## Live Google Sheet access (recommended)

1. In Google Cloud, create a **service account**; enable the **Google Sheets**
   and **Google Drive** APIs; download its **JSON key**.
2. Copy `.streamlit/secrets.toml.example` to `.streamlit/secrets.toml` and paste
   the JSON fields under `[gcp_service_account]`.
3. **Share both sheets** (IC50 + construct) with the service-account email
   (`...@...iam.gserviceaccount.com`), Viewer access is enough.
4. Select **Google Sheet (live)** in the sidebar. The default sheet IDs are
   pre-filled; change them if needed. Data is cached for 10 min — hit
   **🔄 Refresh data** to pull the latest.

## Deploy to Streamlit Community Cloud

Push this folder to GitHub → [share.streamlit.io](https://share.streamlit.io) →
point it at `app.py` → paste the service-account JSON into **Settings → Secrets**
(same `[gcp_service_account]` block). The app then always reflects the live sheet.

## ⚠️ Two blocks to make it match V16 exactly

`hcv_data.py` ships with **reconstructed starter versions** so it runs out of the
box, but to match your V16 results exactly, replace these with your authoritative
copies (search the file for the markers):

- **`===== PASTE BUCKET_MAP =====`** — your full `(experiment, PSVX, day) → bucket`
  map. Anything missing is reported in the dashboard's diagnostics panel as
  "Unmapped bucket" and excluded.
- **`===== PASTE SUBGROUP_RULES =====`** — your full ordered regex list (and
  `SUBGROUP_ORDER`, plus `CUSTOM_LABELS` for the cumulative PVXE181 SG4 labels).
  Unmatched constructs land in "Uncategorized" (also reported in diagnostics).

Optionally fill in `PSV_GENOTYPE` (your 15-strain panel: 4×1a, 3×1b, 8 rare) to
enable genotype-grouped PSV columns.

> Or just send me your V16 `.py` file and I'll wire the exact maps in for you.

## Features

- **Metric toggle** — % neutralization ↔ log₁₀(IC50).
- **Dilution selector** — % neutralization at *any* dilution in your array, not
  just 1:30 (read straight from `Avg_Neut_percent_corrected`).
- **Cell encoding** — Gradient (continuous colorscale) or Threshold hit map
  (green ≥ threshold, blank below, ✕ not tested, ⬛ No Neutralization), with
  ≥50% / ≥75% presets and a custom value.
- **Breadth column** — PSVs passing threshold per construct, rows auto-sorted
  most-neutralizing on top.
- **Dose window** — Prime / Boost1 / Boost2 or pooled.
- **Subgroup view** — one subgroup at a time, or **small multiples** showing all
  subgroups at once.
- **Experiment / Group / PSV filters.**
- **Click a cell → neutralization curve** (% neut vs dilution, log-x); dropdown
  fallback always available.
- **CSV export** of the current view; PNG via the Plotly camera button.
- **Diagnostics panel** flags unmapped-bucket and uncategorized rows.

## Data assumptions

The IC50 sheet is expected to have (names are auto-detected, case-insensitive):
`Experiment`, `Group`, `PSVX_No`, `PSV`, `Day`, `Construct_Description`,
`Sample_Type`, `Dilution`, `Avg_Neut_percent_corrected`, and an IC50 column
(`IC50` or `IC50_corrected`). Filtering matches V15/V16: HCV experiments only
(PVXE141/142/167A/167B/174/181), dropping `GNG`/`G1` groups and antibody/plasma
sample types. `NN`, blank, and `0` IC50 cells are treated as No Neutralization.
