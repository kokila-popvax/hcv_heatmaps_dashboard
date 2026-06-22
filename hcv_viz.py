"""
hcv_viz.py
----------
Plotly figure builders for the dashboard. Depends only on plotly / pandas / numpy.

build_heatmap_figure(): constructs x PSVs heatmap + right-hand breadth bar,
    reproducing the V16 encodings:
      gradient mode  -> continuous colorscale (log10IC50: red->green; %neut: white->green)
      threshold mode -> solid green for hits
      NN cells       -> dark gray with '0.00'
      not-tested     -> blank with a small grey 'x'
    Rows are ordered by breadth (most-neutralizing construct on top).

build_curve_figure(): mean % neutralization vs dilution for a clicked cell.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

try:
    from scipy.optimize import curve_fit
    _SCIPY = True
except ImportError:
    _SCIPY = False


# log10(IC50): low potency (red) -> mid (yellow) -> high potency (green)
IC50_COLORSCALE = [
    [0.00, "#C0392B"], [0.30, "#E67E22"], [0.55, "#F4D03F"],
    [0.78, "#82C341"], [1.00, "#1E8449"],
]
# % neutralization: white -> green
NEUT_COLORSCALE = [[0.0, "#FFFFFF"], [0.5, "#A9DFBF"], [1.0, "#1E8449"]]

NN_COLOR = "#363636"
GRID_COLOR = "#dddddd"
HIT_GREEN = "#1E8449"


def _wrap(label: str, width: int = 55) -> str:
    """Soft-wrap long construct labels onto multiple lines for the y-axis."""
    s = str(label)
    if len(s) <= width:
        return s
    words, lines, cur = s.replace("+", " + ").split(), [], ""
    for w in words:
        if len(cur) + len(w) + 1 > width and cur:
            lines.append(cur)
            cur = w
        else:
            cur = f"{cur} {w}".strip()
    if cur:
        lines.append(cur)
    return "<br>".join(lines)


def build_heatmap_figure(value_pivot: pd.DataFrame,
                         status_pivot: pd.DataFrame,
                         counts: pd.Series,
                         metric: str,
                         mode: str,
                         threshold: float,
                         title: str = "",
                         psv_genotype: dict | None = None,
                         show_values: bool = True,
                         row_height: int = 36,
                         label_max_chars: int | None = None) -> go.Figure:
    """Return a clean Plotly heatmap (no breadth bar)."""
    if value_pivot.empty:
        fig = go.Figure()
        fig.add_annotation(text="No data for the current filters",
                           showarrow=False, font=dict(size=16, color="#888"))
        fig.update_layout(height=300, xaxis_visible=False, yaxis_visible=False)
        return fig

    constructs = list(value_pivot.index)[::-1]
    psvs = list(value_pivot.columns)
    value_pivot = value_pivot.reindex(constructs)
    status_pivot = status_pivot.reindex(constructs)

    x_labels = psvs
    if label_max_chars:
        y_labels = [
            (c[:label_max_chars] + "…") if len(c) > label_max_chars else c
            for c in constructs
        ]
    else:
        y_labels = [_wrap(c, width=80) for c in constructs]

    # Build z-matrix
    if mode == "threshold":
        z = np.where(status_pivot.values == "hit", 1.0, np.nan)
        colorscale = [[0, HIT_GREEN], [1, HIT_GREEN]]
        zmin, zmax, showscale = 0, 1, False
        cbar_title = ""
    elif metric == "log10_ic50":
        z = np.where(status_pivot.values == "hit", value_pivot.values, np.nan)
        colorscale, zmin, zmax, showscale = IC50_COLORSCALE, 1.0, 5.0, True
        cbar_title = "log₁₀(IC50)"
    else:
        z = np.where(status_pivot.values == "hit", value_pivot.values, np.nan)
        colorscale, zmin, zmax, showscale = NEUT_COLORSCALE, 0.0, 100.0, True
        cbar_title = "% neutralization"

    customdata = status_pivot.values
    unit = "%" if metric == "pct_neut" else "log₁₀"
    hover = ("<b>%{y}</b><br>PSV: %{x}<br>"
             "value: %{z:.2f} " + unit + "<br>status: %{customdata}<extra></extra>")

    fig = go.Figure()
    fig.add_trace(go.Heatmap(
        z=z, x=x_labels, y=y_labels, customdata=customdata,
        colorscale=colorscale, zmin=zmin, zmax=zmax, showscale=showscale,
        xgap=1.5, ygap=1.5, hovertemplate=hover,
        colorbar=dict(title=cbar_title, thickness=14, len=0.6)
        if showscale else None,
    ))

    # Overlays: NN cells, not-tested 'x', optional value text
    nn_shapes      = []   # filled dark-gray rectangles for NN cells
    nn_annotations = []   # white "0.00" labels — must be annotations, not traces,
                          # so they render above the layer="above" shapes
    nt_x, nt_y = [], []
    txt_x, txt_y, txt = [], [], []  # black value labels on coloured cells

    for yi, c in enumerate(constructs):
        for xi, p in enumerate(psvs):
            st = str(status_pivot.iat[yi, xi])
            if st == "nn":
                nn_shapes.append(dict(
                    type="rect", layer="above",
                    xref="x", yref="y",
                    x0=xi - 0.5, x1=xi + 0.5,
                    y0=yi - 0.5, y1=yi + 0.5,
                    fillcolor=NN_COLOR,
                    line=dict(width=0),
                ))
                nn_annotations.append(dict(
                    x=p, y=y_labels[yi],
                    text="0.00",
                    showarrow=False,
                    font=dict(size=8, color="white"),
                    xref="x", yref="y",
                ))
            elif st == "not_tested":
                nt_x.append(xi); nt_y.append(yi)
            elif st == "hit" and show_values and mode != "threshold":
                v = value_pivot.iat[yi, xi]
                if pd.notna(v):
                    txt_x.append(xi); txt_y.append(yi); txt.append(f"{v:.1f}")
    if nt_x:
        fig.add_trace(go.Scatter(
            x=[x_labels[i] for i in nt_x], y=[y_labels[i] for i in nt_y],
            mode="markers",
            marker=dict(symbol="x-thin", size=7, color="#b0b0b0",
                        line=dict(width=1.2, color="#b0b0b0")),
            hovertemplate="not tested<extra></extra>", showlegend=False))
    if txt_x:
        fig.add_trace(go.Scatter(
            x=[x_labels[i] for i in txt_x], y=[y_labels[i] for i in txt_y],
            mode="text", text=txt,
            textfont=dict(size=8, color="black"),
            hoverinfo="skip", showlegend=False))

    # ── Breadth counts (right-side annotations) ────────────────────────────
    # Numerator:   tested & not NN  → status in {hit, miss}
    # Denominator: all tested PSVs  → status in {hit, miss, nn}
    tested_statuses = {"hit", "miss", "nn"}
    positive_statuses = {"hit", "miss"}
    breadth_annotations = []
    for yi, c in enumerate(constructs):
        row_statuses = [str(status_pivot.iat[yi, xi]) for xi in range(len(psvs))]
        n_tested  = sum(1 for s in row_statuses if s in tested_statuses)
        n_positive = sum(1 for s in row_statuses if s in positive_statuses)
        if n_tested > 0:
            breadth_annotations.append(dict(
                x=1.22, y=y_labels[yi],
                xref="paper", yref="y",
                text=f"<b>{n_positive}</b>/{n_tested}",
                showarrow=False,
                font=dict(size=9, color="#333"),
                xanchor="left", yanchor="middle",
            ))

    sg_annotations = []

    n_rows = len(constructs)
    n_cols = len(psvs)
    # Estimate left margin from the longest (possibly wrapped) y-label
    effective_label_len = label_max_chars if label_max_chars else max(
        (len(lbl.replace("<br>", " ")) for lbl in y_labels), default=40)
    left_margin = max(120, min(effective_label_len * 7, 340))
    # Ensure minimum total width so y-axis labels are never squeezed off-screen
    # Each cell ~60px, plus left margin + right margin (200) + breadth col
    min_width = left_margin + 200 + max(n_cols * 60, 120)
    fig.update_layout(
        title=dict(text=title, font=dict(size=14, color="#222"),
                   y=1.0, yanchor="top", pad=dict(t=8, b=0)),
        height=max(400, 80 + n_rows * row_height),
        width=min_width,
        margin=dict(l=left_margin, r=200, t=200, b=20),
        plot_bgcolor="white", paper_bgcolor="white",
        font=dict(family="Helvetica, Arial, sans-serif", size=10),
        shapes=nn_shapes,
        annotations=nn_annotations + breadth_annotations + sg_annotations,
    )
    fig.update_xaxes(side="top", tickangle=45, showgrid=False,
                     tickfont=dict(size=11, color="black"), automargin=True)
    fig.update_yaxes(autorange="reversed", showgrid=False,
                     tickfont=dict(size=10, color="black"), automargin=True)
    return fig


def _hill4pl(x, top, bottom, ic50, slope):
    """4-parameter logistic (Hill) model.
    x    = dilution (1:x) — higher x means more dilute → less neutralization
    top  = % neut at lowest dilution (most concentrated)
    bottom = % neut at highest dilution (most dilute)
    ic50 = dilution giving 50 % of (top-bottom) drop → neutralization midpoint
    slope = Hill slope (>0; curve falls as x increases)
    """
    return bottom + (top - bottom) / (1.0 + (x / ic50) ** slope)


def build_curve_figure(curve: pd.DataFrame, construct: str, psv: str,
                       threshold_pct: float | None = 50.0) -> go.Figure:
    """4PL-fitted neutralization curve vs dilution (log-x) for one construct × PSV.

    IC50 is the dilution (1:x) at which % neutralization crosses 50 % of the
    dynamic range — i.e. a *higher* IC50 number means the serum is still
    neutralizing at higher dilution → stronger response.
    """
    fig = go.Figure()
    if curve.empty:
        fig.add_annotation(text="No neutralization-curve data for this pair",
                           showarrow=False, font=dict(size=14, color="#888"))
        fig.update_layout(height=300)
        return fig

    x_obs = curve["dilution"].values.astype(float)
    y_obs = curve["pct_neut"].values.astype(float)

    # --- raw data points ---
    fig.add_trace(go.Scatter(
        x=x_obs, y=y_obs,
        mode="markers", marker=dict(size=8, color="#1E8449", opacity=0.8),
        name="observed",
        hovertemplate="1:%{x:.0f}<br>%{y:.1f}%<extra></extra>"))

    fit_label = ""
    if _SCIPY and len(x_obs) >= 4:
        try:
            p0 = [max(y_obs), min(y_obs), float(np.median(x_obs)), 1.0]
            bounds = ([0, -10, x_obs.min() * 0.01, 0.1],
                      [120, 100, x_obs.max() * 100, 10])
            popt, _ = curve_fit(_hill4pl, x_obs, y_obs, p0=p0, bounds=bounds,
                                maxfev=8000)
            top, bottom, ic50_fit, slope = popt

            x_fit = np.logspace(np.log10(x_obs.min() * 0.5),
                                np.log10(x_obs.max() * 2), 300)
            y_fit = _hill4pl(x_fit, *popt)
            fig.add_trace(go.Scatter(
                x=x_fit, y=y_fit,
                mode="lines", line=dict(color="#1E8449", width=2.5),
                name="4PL fit",
                hovertemplate="1:%{x:.0f}<br>%{y:.1f}% (fit)<extra></extra>"))

            # mark IC50 point
            y_ic50 = _hill4pl(ic50_fit, *popt)
            fig.add_trace(go.Scatter(
                x=[ic50_fit], y=[y_ic50],
                mode="markers",
                marker=dict(size=12, color="#E67E22", symbol="diamond"),
                name=f"IC50 = 1:{ic50_fit:.0f}",
                hovertemplate=f"IC50 = 1:{ic50_fit:.0f}<br>{y_ic50:.1f}%<extra></extra>"))
            fit_label = f"  |  IC50 = 1:{ic50_fit:.0f}  (slope={slope:.2f})"
        except Exception:
            # fit failed — just show raw points
            pass
    elif not _SCIPY:
        fig.add_annotation(
            text="Install scipy for curve fitting (pip install scipy)",
            showarrow=False, font=dict(size=11, color="#E67E22"),
            xref="paper", yref="paper", x=0.5, y=0.05)

    if threshold_pct is not None:
        fig.add_hline(y=threshold_pct, line_dash="dash", line_color="#C0392B",
                      annotation_text=f"{threshold_pct:.0f}%",
                      annotation_position="right")

    short_c = (construct[:60] + "…") if len(construct) > 60 else construct
    fig.update_layout(
        title=dict(text=f"{short_c}  ×  {psv}{fit_label}", font=dict(size=12)),
        height=340, margin=dict(l=60, r=20, t=70, b=40),
        plot_bgcolor="white", paper_bgcolor="white",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        xaxis=dict(title="Dilution (1 : x)  — higher = more dilute",
                   type="log", showgrid=True, gridcolor=GRID_COLOR),
        yaxis=dict(title="% neutralization", range=[-5, 110],
                   showgrid=True, gridcolor=GRID_COLOR),
    )
    return fig
