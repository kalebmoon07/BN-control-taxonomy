#!/usr/bin/env python3
if __name__ == "__main__":
    import sys
    from os.path import dirname, abspath
    libdir = dirname(dirname(abspath(__file__)))
    sys.path.insert(0, libdir)

import argparse
import math
import sys
import os

from pathlib import Path
from typing import Dict, Iterable, List, Set

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.ticker import FuncFormatter, MaxNLocator, NullLocator

from bntaxonomy.hierarchy import MultiInputSummary


# ---------------------------------------------------------------------
# Visual constants (inches) — keep bars/pads uniform across inputs
# ---------------------------------------------------------------------
BAR_W_IN = 0.22  # width of each bar (inches)
GAP_IN = 0.12  # gap between neighboring bars (inches)
MARGIN_LR_IN = 0.6  # left/right figure margins (inches)
MARGIN_TB_IN = 0.6  # top/bottom figure margins (inches)
WSPACE_IN = 0.35  # inter-subplot horizontal spacing (inches)
HSPACE_IN = 0.45  # inter-subplot vertical spacing (inches)
FULL_FIG_XLABEL_PAD_IN = 0.55  # reserve bottom space for the full-figure x label

SIGN_COLORS = {1: "tab:blue", 0: "tab:red", -1: "tab:red"}
SIGN_NAME = {1: "Positive", 0: "Negative", -1: "Negative"}
LINEAR_SCORE_TICKS = [-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1]
# Symlog mode: linear region down to 1e-5; ticks cover decades 10^-5 .. 10^0
# on each sign (five decades per half plus 0 at the centre = 11 ticks).
SCORE_LOG_LINTHRESH = 1e-5
SCORE_LOG_EXP_RANGE = (-5, 0)  # exponents of the smallest and largest ticks


def _slot_in() -> float:
    """One x-slot = one bar + one gap, in inches."""
    return BAR_W_IN + GAP_IN


def _bar_width_frac() -> float:
    """Bar width in data-units when 1 x-unit == one slot."""
    return BAR_W_IN / (BAR_W_IN + GAP_IN)


def _compute_figsize_grid(n_bars: int, rows: int, cols: int, panel_h_in: float = 3.6):
    """
    For a grid of subplots where each panel has the same number of bars (n_bars),
    compute a figure size that keeps bar/gap widths constant in inches.
    """
    content_w_in = max(1, n_bars) * _slot_in()
    fig_w = (cols * content_w_in) + (cols - 1) * WSPACE_IN + 2 * MARGIN_LR_IN
    fig_h = (rows * panel_h_in) + (rows - 1) * HSPACE_IN + 2 * MARGIN_TB_IN
    return fig_w, fig_h


# ---------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------
def _format_compact_scientific(value: float) -> str:
    """Single-line compact scientific notation, e.g. ``4e-2`` or ``1.6e-4``.

    Used for y-tick labels where vertical space is tight and two-line
    annotations would double the label height.
    """
    value = abs(float(value))
    if np.isclose(value, 0.0):
        return "0"

    mantissa, exponent = f"{value:.1e}".split("e")
    mantissa = mantissa.rstrip("0").rstrip(".")
    return f"{mantissa}e{int(exponent)}"


def _format_compact_scientific_two_line(value: float) -> str:
    """Two-line compact scientific notation for in-cell bar annotations.

    Mantissa is kept at one decimal (``6.0`` rather than ``6``), followed
    by a newline and the ``e{exponent}`` suffix. So ``0.011`` renders as::

        1.1
        e-2
    """
    value = abs(float(value))
    if np.isclose(value, 0.0):
        return "0"

    mantissa, exponent = f"{value:.1e}".split("e")
    return f"{mantissa}\ne{int(exponent)}"


def format_score_label(value: float, log_scale: bool,
                       *, two_line: bool = False) -> str:
    if not log_scale:
        return f"{abs(value):.2f}"
    return (_format_compact_scientific_two_line(value) if two_line
            else _format_compact_scientific(value))


def annotate_bars(ax, bars, formatter=None, offset_points=3, fontsize=8):
    """
    Put value labels at the end of each bar.
    Offset is applied in screen points so labels remain readable for linear
    and logarithmic axes alike.
    """
    if formatter is None:
        formatter = lambda value: f"{abs(value):.2f}"

    for b in bars:
        h = b.get_height()
        if not h:  # skip zero-height bars
            continue
        x = b.get_x() + b.get_width() / 2.0
        y = b.get_y() + h
        va = "bottom" if h > 0 else "top"
        yo = offset_points if h > 0 else -offset_points
        ax.annotate(
            formatter(h),
            xy=(x, y),
            xytext=(0, yo),
            textcoords="offset points",
            ha="center",
            va=va,
            fontsize=fontsize,
        )


def _gene_sign_wide(df_inst: pd.DataFrame) -> pd.DataFrame:
    """Aggregate scores per (Gene, Sign) and return wide table with pos/neg/total."""
    g: pd.DataFrame = (
        df_inst.reset_index()
        .groupby(["Gene", "Sign"], as_index=False, observed=True)["score"]
        .sum()
    )

    wide = g.pivot(index="Gene", columns="Sign", values="score").fillna(0.0)
    wide = wide.rename(columns={1: "pos", 0: "neg_raw"})

    # Ensure columns exist even if a sign is missing in the data
    if "pos" not in wide.columns:
        wide["pos"] = 0.0
    if "neg_raw" not in wide.columns:
        wide["neg_raw"] = 0.0

    wide["neg"] = wide["neg_raw"].abs()
    wide["total"] = wide["pos"] + wide["neg"]
    return wide


def _sort_genes(df_inst: pd.DataFrame, by: list[str], ascending: list[bool]) -> list[str]:
    wide = _gene_sign_wide(df_inst)
    return (
        wide.sort_values(by=by, ascending=ascending, kind="mergesort")
        .index.to_list()
    )


def sort_by_total_score(df_inst: pd.DataFrame) -> list[str]:
    # primary: total (desc), secondary: pos (desc)
    return _sort_genes(df_inst, by=["total", "pos"], ascending=[False, False])


def sort_by_neg_score(df_inst: pd.DataFrame) -> list[str]:
    # primary: neg (desc), secondary: total (desc)
    return _sort_genes(df_inst, by=["neg", "total"], ascending=[False, False])


def sort_by_pos_score(df_inst: pd.DataFrame) -> list[str]:
    # primary: pos (desc), secondary: total (desc)
    return _sort_genes(df_inst, by=["pos", "total"], ascending=[False, False])



def layout_single_axes(fig_w_in: float, fig_h_in: float):
    """Return (fig, ax) for a single-axes figure with inch-based margins applied."""
    fig, ax = plt.subplots(1, 1, figsize=(fig_w_in, fig_h_in))
    fig.subplots_adjust(
        left=MARGIN_LR_IN / fig_w_in,
        right=1 - MARGIN_LR_IN / fig_w_in,
        top=1 - MARGIN_TB_IN / fig_h_in,
        bottom=MARGIN_TB_IN / fig_h_in,
        wspace=0,
    )
    return fig, ax

def annotate_score_axis(fig: plt.Figure, ax: plt.Axes):
    """Annotate the score axis with sign labels."""
    def _tick_midpoint_in_axes(sign: int) -> float:
        ticks = np.asarray(ax.get_yticks(), dtype=float)
        y0, y1 = ax.get_ylim()
        lo, hi = sorted((y0, y1))
        ticks = ticks[(ticks >= lo) & (ticks <= hi)]

        if sign > 0:
            ticks = ticks[ticks > 0]
            fallback = 0.75
        else:
            ticks = ticks[ticks < 0]
            fallback = 0.25

        if len(ticks) == 0:
            return fallback

        tick_xy = np.column_stack([np.zeros(len(ticks)), ticks])
        tick_disp = ax.transData.transform(tick_xy)[:, 1]
        midpoint_disp = 0.5 * (tick_disp.min() + tick_disp.max())
        return ax.transAxes.inverted().transform((0, midpoint_disp))[1]

    # Put labels just left of the y-axis with fixed physical padding
    trans = ax.transAxes

    # x = 0 is the y-axis in axes coords; move left by 18 points
    pad = mtransforms.ScaledTranslation(-0.55, 0, fig.dpi_scale_trans)

    ax.text(
        0,
        _tick_midpoint_in_axes(1),
        "Activation",
        transform=trans + pad,
        rotation=90,
        va="center",
        ha="center",
        clip_on=False,
    )
    ax.text(
        0,
        _tick_midpoint_in_axes(-1),
        "Inhibition",
        transform=trans + pad,
        rotation=90,
        va="center",
        ha="center",
        clip_on=False,
    )


def annotate_score_axis_fixed(fig: plt.Figure, ax: plt.Axes,
                              pos_y: float = 0.75, neg_y: float = 0.25):
    """Place sign labels at fixed axes-fraction positions.

    Used by the count axis so that 'Activation' / 'Inhibition' sit at
    identical heights regardless of how MaxNLocator chose tick spacing
    for a particular figure.
    """
    trans = ax.transAxes
    pad = mtransforms.ScaledTranslation(-0.55, 0, fig.dpi_scale_trans)
    ax.text(0, pos_y, "Activation", transform=trans + pad,
            rotation=90, va="center", ha="center", clip_on=False)
    ax.text(0, neg_y, "Inhibition", transform=trans + pad,
            rotation=90, va="center", ha="center", clip_on=False)


def _build_symlog_ticks(limit: float, linthresh: float) -> list[float]:
    """Return symmetric powers-of-ten ticks up to the visible y-limit."""
    limit = max(float(limit), float(linthresh))
    positive_ticks = []
    tick = linthresh
    while tick <= limit * 1.0000001:
        positive_ticks.append(tick)
        tick *= 10
    return [-t for t in reversed(positive_ticks)] + [0.0] + positive_ticks


def _fixed_log_ticks() -> list[float]:
    """Symlog ticks at powers of 10 over ``SCORE_LOG_EXP_RANGE``.

    The outer decades are inclusive, so the default (-5, 0) produces
    ticks at 10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 10^0 on each sign, plus
    0 at the centre.
    """
    lo, hi = SCORE_LOG_EXP_RANGE
    positive = [10 ** e for e in range(lo, hi + 1)]
    return [-t for t in reversed(positive)] + [0.0] + positive


def _score_tick_formatter(log_scale: bool) -> FuncFormatter:
    # Suppress tick labels above |1| — the extra headroom (up to ±50) is
    # only there so bar-top annotations don't collide with the axis edge.
    def _suppress_above_one(formatter):
        def wrapped(x, pos):
            if abs(x) > 1.0 + 1e-9:
                return ""
            return formatter(x, pos)
        return wrapped

    if log_scale:
        base = lambda x, pos: format_score_label(x, log_scale=True)
    else:
        base = lambda x, pos: f"{abs(x):.2f}"
    return FuncFormatter(_suppress_above_one(base))


def configure_score_axis(
    ax: plt.Axes,
    fig: plt.Figure,
    values,
    log_scale: bool,
    *,
    show_ticklabels: bool = True,
    show_side_labels: bool = True,
):
    """Apply a consistent score axis style to score plots."""
    max_abs = float(np.max(np.abs(np.asarray(values, dtype=float)))) if len(values) else 0.0
    max_abs = max(max_abs, SCORE_LOG_LINTHRESH)

    if log_scale:
        # Fixed tick ladder: powers of 10 across ``SCORE_LOG_EXP_RANGE``
        # on each sign (10^-5 ... 10^0) plus 0 at the centre. The visible
        # y-limit is intentionally well beyond the top tick so bar-top
        # annotations have breathing room without introducing |value| > 1
        # tick marks.
        ticks = _fixed_log_ticks()
        y_limit = 50.0

        ax.set_yscale("symlog", linthresh=SCORE_LOG_LINTHRESH, linscale=1.0)
        ax.set_ylim(-y_limit, y_limit)
        ax.set_yticks(ticks)
        # symlog's default minor locator keeps emitting decade ticks
        # (10, 100, ...) inside the wider ylim — suppress them.
        ax.yaxis.set_minor_locator(NullLocator())
    else:
        ax.set_ylim(-1.20, 1.20)
        ax.set_yticks(LINEAR_SCORE_TICKS)

    ax.yaxis.set_major_formatter(_score_tick_formatter(log_scale))
    if show_ticklabels:
        ax.tick_params(axis="y", labelleft=True)
    else:
        ax.tick_params(axis="y", labelleft=False)

    if show_side_labels:
        ax.set_ylabel(" ")
        if log_scale:
            annotate_score_axis(fig, ax)
        else:
            annotate_score_axis_fixed(fig, ax)
    else:
        ax.set_ylabel(None)


def configure_count_axis(
    ax: plt.Axes,
    fig: plt.Figure,
    values,
    *,
    show_ticklabels: bool = True,
    show_side_labels: bool = True,
):
    """Apply a symmetric linear axis for signed solution counts."""
    max_abs = int(np.max(np.abs(np.asarray(values, dtype=float)))) if len(values) else 0
    y_limit = max(1, max_abs)
    ax.set_ylim(-1.15 * y_limit, 1.15 * y_limit)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=7, integer=True, symmetric=True))
    ax.yaxis.set_major_formatter(
        FuncFormatter(lambda x, pos: "0" if np.isclose(x, 0.0) else f"{int(round(abs(x)))}")
    )
    if show_ticklabels:
        ax.tick_params(axis="y", labelleft=True)
    else:
        ax.tick_params(axis="y", labelleft=False)

    if show_side_labels:
        ax.set_ylabel(" ")
        annotate_score_axis_fixed(fig, ax)
    else:
        ax.set_ylabel(None)


# ---------------------------------------------------------------------
# MCS pipeline (merged from former evaluate_mcs2)
#
# Per-instance outputs, written under ``{output}/{inst_group}/{inst}/mcs/``:
#   _family_partition.csv
#   _family_jaccard.csv
#   _family_jaccard_heatmap.{fmt}
#   _topk_spearman_{inhibition,activation}.csv
#   _topk_spearman_{inhibition,activation}_heatmap.{fmt}
#
# Family partition is derived by Jaccard = 1.0 union-find over the
# tool-level control-set matrix, ranked in input tool order.
# ---------------------------------------------------------------------
def _jaccard(a: Set, b: Set) -> float:
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def _containment(a: Set, b: Set) -> float:
    """Asymmetric set containment |A ∩ B| / |A|.

    Empty A is treated as vacuously contained (returns 1.0), matching the
    convention that ∅ ⊆ B for any B.
    """
    if not a:
        return 1.0
    return len(a & b) / len(a)


def _geo_mean(values: Iterable[float], eps: float) -> float:
    arr = np.clip(np.asarray(list(values), dtype=float), 0.0, None)
    if arr.size == 0:
        return 0.0
    return float(np.exp(np.mean(np.log(arr + eps))) - eps)


def _derive_families(ctrls: Dict[str, Set]) -> Dict[str, List[str]]:
    tools = list(ctrls)
    pos = {t: i for i, t in enumerate(tools)}
    parent = {t: t for t in tools}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, a in enumerate(tools):
        for b in tools[i + 1:]:
            if _jaccard(ctrls[a], ctrls[b]) == 1.0:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[rb if pos[ra] < pos[rb] else ra] = (
                        ra if pos[ra] < pos[rb] else rb
                    )

    groups: Dict[str, List[str]] = {}
    for t in tools:
        groups.setdefault(find(t), []).append(t)
    fams: Dict[str, List[str]] = {}
    for idx, (rep, members) in enumerate(
        sorted(groups.items(), key=lambda kv: pos[kv[0]]), start=1
    ):
        fams[f"fam{idx}_{rep}"] = members
    return fams


def _ctrls_from_exp(exp, tool_names: List[str]) -> Dict[str, Set]:
    by_name = {r.name: r for r in exp.results}
    return {
        n: {frozenset((g, int(v)) for g, v in c.items())
            for c in by_name[n].d_list}
        for n in tool_names if n in by_name
    }


def _per_family_mcs(score_df: pd.DataFrame, families: Dict[str, List[str]],
                    sign: int, geo: bool, eps: float) -> pd.DataFrame:
    fam_list = list(families)
    tool_to_fam = {t: f for f, ts in families.items() for t in ts}
    df = score_df[score_df["Sign"] == sign].copy()
    df["family"] = df["Algorithm"].map(tool_to_fam)
    df = df.dropna(subset=["family"])

    if geo:
        agg = df.groupby(["Gene", "family"])["score"].apply(
            lambda s: _geo_mean(s, eps=eps)
        )
    else:
        agg = df.groupby(["Gene", "family"])["score"].mean()
    pivot = agg.reset_index().pivot(
        index="Gene", columns="family", values="score"
    ).fillna(0.0)
    return pivot[[f for f in fam_list if f in pivot.columns]]


def _ensemble_series(pivot: pd.DataFrame, geo: bool, eps: float) -> pd.Series:
    if geo:
        s = pivot.apply(lambda r: _geo_mean(r, eps=eps), axis=1)
    else:
        s = pivot.mean(axis=1)
    return s.sort_values(ascending=False)


# Absolute layout constants: fixed fonts, fixed inch margins, figure
# size follows matrix shape. Axes and colorbar placed via add_axes to
# make the margins immune to constrained_layout label re-flow.
_CELL_IN = 0.55
_LEFT_IN = 2.3
_RIGHT_IN = 0.3
_CBAR_W_IN = 0.25
_CBAR_PAD_IN = 0.7
_TOP_IN = 0.9
_BOTTOM_IN = 2.3
_FONT_ANNOT = 7
_FONT_TICK = 9
_FONT_TITLE = 11


def _render_heatmap(mat: pd.DataFrame, outpath: Path, title: str,
                    cmap: str, vmin: float, vmax: float) -> None:
    n_rows, n_cols = len(mat.index), len(mat.columns)
    axes_w = _CELL_IN * n_cols
    axes_h = _CELL_IN * n_rows
    fig_w = _LEFT_IN + axes_w + _RIGHT_IN + _CBAR_W_IN + _CBAR_PAD_IN
    fig_h = _TOP_IN + axes_h + _BOTTOM_IN

    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([
        _LEFT_IN / fig_w, _BOTTOM_IN / fig_h,
        axes_w / fig_w, axes_h / fig_h,
    ])
    im = ax.imshow(mat.values, vmin=vmin, vmax=vmax, cmap=cmap)
    ax.set_xticks(range(n_cols))
    ax.set_yticks(range(n_rows))
    ax.set_xticklabels(mat.columns, rotation=45, ha="right",
                       fontsize=_FONT_TICK)
    ax.set_yticklabels(mat.index, fontsize=_FONT_TICK)
    cmap_obj = plt.get_cmap(cmap)
    span = (vmax - vmin) or 1.0
    for i in range(n_rows):
        for j in range(n_cols):
            v = mat.values[i, j]
            if np.isnan(v):
                ax.text(j, i, "-", ha="center", va="center",
                        color="black", fontsize=_FONT_ANNOT)
                continue
            r, g, b, _ = cmap_obj((v - vmin) / span)
            lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    color="white" if lum < 0.55 else "black",
                    fontsize=_FONT_ANNOT)
    ax.set_title(title, pad=8, fontsize=_FONT_TITLE)

    cax = fig.add_axes([
        (_LEFT_IN + axes_w + _RIGHT_IN) / fig_w, _BOTTOM_IN / fig_h,
        _CBAR_W_IN / fig_w, axes_h / fig_h,
    ])
    cb = fig.colorbar(im, cax=cax)
    cb.ax.tick_params(labelsize=_FONT_TICK)
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _write_family_partition(ctrls: Dict[str, Set],
                            families: Dict[str, List[str]],
                            outdir: Path) -> None:
    rows = []
    for fam, tools in families.items():
        sizes = {len(ctrls[t]) for t in tools if t in ctrls}
        rows.append({
            "family": fam,
            "n_members": len(tools),
            "n_controls": sizes.pop() if len(sizes) == 1 else -1,
            "representative": tools[0] if tools else "",
            "members": ",".join(tools),
        })
    pd.DataFrame(rows).to_csv(outdir / "_family_partition.csv", index=False)


def _write_family_jaccard(ctrls: Dict[str, Set],
                          families: Dict[str, List[str]],
                          outdir: Path, fmt: str, instance: str) -> None:
    fam_sets: Dict[str, Set] = {}
    for fam, tools in families.items():
        members = [ctrls[t] for t in tools if t in ctrls]
        fam_sets[fam] = set.intersection(*members) if members else set()
    names = list(fam_sets)
    M = np.zeros((len(names), len(names)))
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            M[i, j] = _jaccard(fam_sets[a], fam_sets[b])
    df = pd.DataFrame(M, index=names, columns=names)
    df.to_csv(outdir / "_family_jaccard.csv")
    _render_heatmap(
        df, outdir / f"_family_jaccard_heatmap.{fmt}",
        f"Family-level Jaccard ({instance})",
        cmap="YlGn", vmin=0.0, vmax=1.0,
    )


def _write_family_containment(ctrls: Dict[str, Set],
                              families: Dict[str, List[str]],
                              outdir: Path, fmt: str, instance: str) -> None:
    fam_sets: Dict[str, Set] = {}
    for fam, tools in families.items():
        members = [ctrls[t] for t in tools if t in ctrls]
        fam_sets[fam] = set.intersection(*members) if members else set()
    names = list(fam_sets)
    M = np.zeros((len(names), len(names)))
    for i, row in enumerate(names):
        for j, col in enumerate(names):
            # Cell (row=y, col=x) holds |x ∩ y| / |x|.
            M[i, j] = _containment(fam_sets[col], fam_sets[row])
    df = pd.DataFrame(M, index=names, columns=names)
    df.to_csv(outdir / "_family_containment.csv")
    _render_heatmap(
        df, outdir / f"_family_containment_heatmap.{fmt}",
        f"Family-level set containment |x∩y|/|x| ({instance})",
        cmap="YlGn", vmin=0.0, vmax=1.0,
    )


def _write_topk_spearman_vs_ensemble(
    score_df: pd.DataFrame, families: Dict[str, List[str]],
    outdir: Path, sign: int, topk: Iterable[int],
    geo: bool, eps: float, fmt: str, instance: str,
) -> None:
    pivot = _per_family_mcs(score_df, families, sign, geo, eps)
    if pivot.shape[1] < 2:
        return
    ens = _ensemble_series(pivot, geo, eps)
    rows: Dict[int, List[float]] = {}
    # A constant-vector family (e.g. release-empty tools with zero MCS
    # everywhere) gives stddev = 0 inside numpy's corrcoef and raises a
    # "invalid value encountered in divide" RuntimeWarning. The resulting
    # NaN is already handled by .fillna(0.0) below, so silencing the
    # warning is intentional and scoped.
    with np.errstate(invalid="ignore", divide="ignore"):
        for k in topk:
            if k > len(pivot):
                continue
            idx = ens.head(k).index
            sub = pivot.loc[idx]
            ens_sub = ens.loc[idx]
            # Explicit fractional-rank tie handling: Spearman via Pearson
            # on average-method ranks, so tied MCS values share the mean
            # rank.
            ens_rank = ens_sub.rank(method="average")
            rows[k] = [
                sub[f].rank(method="average").corr(ens_rank, method="pearson")
                for f in pivot.columns
            ]
    if not rows:
        return
    df = (pd.DataFrame.from_dict(rows, orient="index", columns=pivot.columns)
            .fillna(0.0).sort_index(ascending=False))
    df.index = [f"top-{k}" for k in df.index]
    label = "inhibition" if sign == 0 else "activation"
    df.to_csv(outdir / f"_topk_spearman_{label}.csv")
    _render_heatmap(
        df, outdir / f"_topk_spearman_{label}_heatmap.{fmt}",
        f"Top-k Spearman rho vs ensemble ({label}) -- {instance}",
        cmap="RdBu_r", vmin=-1.0, vmax=1.0,
    )


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main(argv=None):
    parser = argparse.ArgumentParser(description="Compute/plot control scores.")
    parser.add_argument(
        "-t",
        "--tools",
        nargs="+",
        help="Tools to include. If omitted, include all.",
        default=None,
    )
    parser.add_argument(
        "-ig",
        "--inst_groups",
        nargs="+",
        help="Instance groups to include. Default: all.",
        default=list(),
    )
    parser.add_argument(
        "-i",
        "--instances",
        nargs="+",
        help=("List of specific instances to include."),
        default=None,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output directory for results.",
        default=f"experiments/results",
    )
    parser.add_argument(
        "-g",
        "--genes",
        nargs="+",
        help=(
            "Optional list of genes to evaluate and plot. \n"
            "Only these genes will be included, and they will appear in the exact order given."
        ),
        default=None,
    )
    parser.add_argument(
        "--sort",
        choices=["total", "pos", "neg"],
        help="Sorting method for genes in plots (default: total).",
        default="total",
    )
    parser.add_argument(
        "--format",
        choices=["png", "pdf"],
        help="Output figure format (default: png).",
        default="png",
    )
    parser.add_argument(
        "--geo-mean",
        dest="geo_mean",
        action="store_true",
        help=(
            "Aggregate the per-Gene/per-Sign summary plot with the "
            "epsilon-shifted geometric mean across tools instead of the "
            "arithmetic mean."
        ),
    )
    parser.add_argument(
        "--geo-eps",
        type=float,
        default=1e-6,
        help="Epsilon shift for --geo-mean (default: 1e-6).",
    )
    args = parser.parse_args(argv)

    if args.genes:
        selected_gene_order = list(dict.fromkeys(args.genes))
    else:
        selected_gene_order = None

    # selected_tools = set(args.tools) if args.tools else None
    inst_groups = [p.replace("instances", "results") for p in args.inst_groups]
    instances: list[str] = args.instances

    os.makedirs(args.output, exist_ok=True)
    opath = args.output
    os.makedirs(opath, exist_ok=True)

    # -----------------------------------------------------------------
    # Summaries from experiment groups
    # -----------------------------------------------------------------
    if inst_groups:
        hc = MultiInputSummary.from_inst_groups(inst_groups, "Hierarchy")
    elif args.instances:
        hc = MultiInputSummary.from_instances(instances, "Hierarchy")
    else:
        instances_root = os.path.join("experiments", "instances")
        if not os.path.isdir(instances_root):
            raise FileNotFoundError(f"Instances directory not found: {instances_root}")
        group_dirs = sorted(
            d
            for d in os.listdir(instances_root)
            if os.path.isdir(os.path.join(instances_root, d))
        )
        if not group_dirs:
            raise RuntimeError(f"No instance groups found in {instances_root}")
        inst_groups = [os.path.join("experiments", "results", d) for d in group_dirs]
        hc = MultiInputSummary.from_inst_groups(inst_groups, "Hierarchy")

    if args.tools:
        available = {r.name for e in hc.exp_list for r in e.results}
        missing = [a for a in args.tools if a not in available]
        for m in missing:
            print(
                f"Warning: tool '{m}' not found in {inst_groups}; skipping.",
                file=sys.stderr,
            )
        # keep only those present, in the original order
        selected_tools = [a for a in args.tools if a in available]
    else:
        selected_tools = None
    cat_type = pd.CategoricalDtype(selected_tools, ordered=True)

    # Build records (tight inner loops with local bindings for speed)
    count_list = []  # (Instance, Tool, BN_size, ControlSize, Genes)
    mcs_score_list = []  # (Instance, Tool, Gene, Sign, BN_size, Score)
    solution_gene_list = []  # (Instance, Tool, Gene, Sign)
    for exp in hc.exp_list:
        inst_name = exp.name
        n_nodes = len(exp.bn)
        bn_keys = tuple(exp.bn.keys())  # stable order once

        for tool_result in exp.results:
            tool_name = tool_result.name
            if selected_tools is not None and tool_name not in selected_tools:
                continue

            if not tool_result.d_list:
                count_list.append((inst_name, tool_name, n_nodes, math.inf, dict()))

            gene_set = tool_result.get_controlled_gene_set()
            for gene in gene_set:
                for sign in (0, 1):
                    score, ctrl_list = tool_result.compute_mutation_score(
                        gene, sign, n_nodes
                    )
                    mcs_score_list.append(
                        (inst_name, tool_name, gene, sign, n_nodes, score)
                    )

            base_score, _ = tool_result.compute_mutation_score(
                "_DUMMY_GENE_", 1, n_nodes
            )
            for gene in bn_keys:
                if gene in gene_set:
                    continue
                for sign in (0, 1):
                    mcs_score_list.append(
                        (inst_name, tool_name, gene, sign, n_nodes, base_score)
                    )

            for ctrl_dict in tool_result.d_list:
                csize = len(ctrl_dict)
                count_list.append((inst_name, tool_name, n_nodes, csize, ctrl_dict))
                for gene, sign in ctrl_dict.items():
                    solution_gene_list.append((inst_name, tool_name, gene, sign))

    # -----------------------------------------------------------------
    # Histogram (per Instance)
    # NOTE: Histogram keeps overall control-size counts after filtering control sets
    #       to the selected genes (if provided).
    # -----------------------------------------------------------------
    count_df = pd.DataFrame(
        count_list, columns=["Instance", "Algorithm", "BN_size", "ControlSize", "Genes"]
    )
    mcs_score_df = pd.DataFrame(
        mcs_score_list,
        columns=["Instance", "Algorithm", "Gene", "Sign", "BN_size", "score"],
    )
    solution_gene_df = pd.DataFrame(
        solution_gene_list, columns=["Instance", "Algorithm", "Gene", "Sign"]
    )
    if selected_tools:
        count_df["Algorithm"] = count_df["Algorithm"].astype(cat_type)
        mcs_score_df["Algorithm"] = mcs_score_df["Algorithm"].astype(cat_type)
        if not solution_gene_df.empty:
            solution_gene_df["Algorithm"] = solution_gene_df["Algorithm"].astype(cat_type)
    mcs_score_df = mcs_score_df.sort_values(
        by=["Instance", "Algorithm", "Gene", "Sign"]
    )
    mcs_score_df.to_csv(f"{opath}/score.csv", index=False)

    # -----------------------------------------------------------------
    # MCS outputs per instance: family partition, family Jaccard heatmap,
    # and top-k Spearman vs ensemble heatmaps (one per sign). This block
    # uses ``mcs_score_df`` before the Sign-based sign-flip below -- scores
    # here are still non-negative.
    # -----------------------------------------------------------------
    exp_by_name = {exp.name: exp for exp in hc.exp_list}
    for inst in mcs_score_df["Instance"].unique():
        exp = exp_by_name.get(inst)
        if exp is None:
            continue
        available = {r.name for r in exp.results}
        tool_names = (
            [n for n in selected_tools if n in available]
            if selected_tools is not None
            else [r.name for r in exp.results]
        )
        if len(tool_names) < 2:
            continue
        ctrls = _ctrls_from_exp(exp, tool_names)
        families = _derive_families(ctrls)
        inst_group = hc.get_exp_group_name_from_exp(inst) or "Custom"
        mcs_outdir = Path(opath) / inst_group / inst
        mcs_outdir.mkdir(parents=True, exist_ok=True)

        inst_score_df = mcs_score_df[mcs_score_df["Instance"] == inst]
        # Disabled per user request — outputs not currently needed.
        # _write_family_partition(ctrls, families, mcs_outdir)
        # _write_family_jaccard(ctrls, families, mcs_outdir, args.format, inst)
        # _write_family_containment(ctrls, families, mcs_outdir, args.format, inst)
        # for sign in (0, 1):
        #     _write_topk_spearman_vs_ensemble(
        #         inst_score_df, families, mcs_outdir, sign,
        #         args.topk, args.geo_mean, args.geo_eps, args.format, inst,
        #     )

    for inst, sub_df in count_df.groupby("Instance", sort=False):
        inst_group = hc.get_exp_group_name_from_exp(inst)
        ct1 = pd.crosstab(sub_df["Algorithm"], sub_df["ControlSize"]).sort_index()
        # Ensure columns for both control sizes exist (1, 2)
        full_sizes = [1, 2]
        ct1 = ct1.reindex(columns=full_sizes, fill_value=0)

        tools = ct1.index.to_list()
        sizes = full_sizes

        x = np.arange(len(tools))
        total_width = 0.8
        bar_w = total_width / len(sizes)
        offsets = (np.arange(len(sizes)) - (len(sizes) - 1) / 2.0) * bar_w

        fig, ax = plt.subplots(figsize=(14, 5))
        for i, cs in enumerate(sizes):
            heights = ct1[cs].to_numpy()
            rects = ax.bar(
                x + offsets[i], heights, width=bar_w, label=str(cs), alpha=0.85
            )

            # label every bar, including zeros (slight offset so zeros are visible)
            for r, val in zip(rects, heights):
                ax.annotate(
                    str(int(val)),
                    xy=(r.get_x() + r.get_width() / 2, r.get_height()),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    clip_on=False,
                )

        ax.set_xticks(x, tools, rotation=45)
        ax.grid(True, axis="y", alpha=0.25)
        ax.legend(title="ControlSize", loc="center left", bbox_to_anchor=(1.01, 0.5))
        ax.set_ylim(0, max(1, int(ct1.to_numpy().max())) * 1.15)  # headroom for labels)
        plt.subplots_adjust(bottom=0.22, right=0.85)

        os.makedirs(f"{opath}/{inst_group}/{inst}", exist_ok=True)
        fig.savefig(
            f"{opath}/{inst_group}/{inst}/_histogram_sum_tool.{args.format}",
            dpi=200,
            bbox_inches="tight",
        )
        plt.close(fig)

    # -----------------------------------------------------------------
    # Plotting: per-Instance grids (one panel per gene)
    # -----------------------------------------------------------------
    width_frac = _bar_width_frac()

    mcs_score_df["score"] = mcs_score_df.apply(
        lambda row: row["score"] if row["Sign"] == 1 else -row["score"], axis=1
    )

    for inst, sub_df in mcs_score_df.groupby("Instance", sort=False):
        inst_group = hc.get_exp_group_name_from_exp(inst)
        sub_df: pd.DataFrame = sub_df.groupby(
            ["Algorithm", "Gene", "Sign"], as_index=False, observed=True
        )["score"].sum()

        tools_all = sub_df["Algorithm"].unique()
        # Determine order for this instance
        if selected_gene_order is not None:
            genes_order = [
                g for g in selected_gene_order if g in sub_df["Gene"].unique()
            ]
        else:
            genes_order = sub_df["Gene"].unique()
        sub_df = sub_df[sub_df["Gene"].isin(genes_order)]
        if not len(genes_order):
            print(f"[Instance={inst}] No genes to plot, skipping.", file=sys.stderr)
            continue
        # Apply categorical ordering for stable groupby/sort
        sub_df["Gene"] = pd.Categorical(
            sub_df["Gene"], categories=genes_order, ordered=True
        )
        sub_df = sub_df.sort_values(["Gene", "Sign", "Algorithm"])

        m = len(genes_order)

        total = m
        if m <= 3:
            rows, cols = m, 1
        else:
            cols = max(1, int(math.ceil(math.sqrt(total))))
            rows = int(math.ceil(total / cols))

        full_plot_values = sub_df["score"].to_numpy()

        # Pre-compute per-instance summary aggregation once; the same data
        # feeds both the linear and log-scale summary figures.
        if args.geo_mean:
            eps = args.geo_eps

            def _geo_signed(s: pd.Series) -> float:
                # Scores within a (Gene, Sign) group share a sign after the
                # earlier flip: activation stays positive, inhibition is
                # already negated. Geo-mean on the absolute values, then
                # restore the common sign.
                arr = np.asarray(s.to_numpy(), dtype=float)
                if arr.size == 0:
                    return 0.0
                sign = -1.0 if (arr < 0).any() else 1.0
                mag = np.abs(arr)
                gm = float(np.exp(np.mean(np.log(mag + eps))) - eps)
                return sign * gm

            sum_by_gene = (
                sub_df.groupby(["Gene", "Sign"], observed=True)["score"]
                .apply(_geo_signed)
                .reset_index()
            )
        else:
            sum_by_gene = (
                sub_df.groupby(["Gene", "Sign"], observed=True)["score"]
                .mean()
                .reset_index()
            )

        # Sort genes by the same aggregation that the summary plot
        # displays (geo or arith mean) so the gene order matches what the
        # user sees in _score_summary, instead of the raw per-tool sums.
        if args.sort == "total":
            gene_sorted = sort_by_total_score(sum_by_gene)
        elif args.sort == "pos":
            gene_sorted = sort_by_pos_score(sum_by_gene)
        else:  # args.sort == "neg"
            gene_sorted = sort_by_neg_score(sum_by_gene)
        sum_by_gene["Gene"] = pd.Categorical(
            sum_by_gene["Gene"], categories=gene_sorted, ordered=True
        )
        sum_by_gene = sum_by_gene.sort_values("Gene").set_index("Sign")

        n_gene = len(genes_order)
        content_gene_w = max(1, n_gene) * _slot_in()
        fig_w_s = content_gene_w + 2 * MARGIN_LR_IN
        fig_h_s = 4

        os.makedirs(f"{opath}/{inst_group}/{inst}", exist_ok=True)

        # Emit both linear (default) and log-scale versions of _score_full
        # and _score_summary for every instance.
        for log_scale, suffix in [(False, ""), (True, "_log_scale")]:
            # Figure size derived from desired physical bar/pad sizes
            panel_h_in = 2.9 if log_scale else 3.4
            fig_w, fig_h = _compute_figsize_grid(
                len(tools_all), rows, cols, panel_h_in=panel_h_in
            )
            fig, axes = plt.subplots(
                rows, cols, figsize=(fig_w, fig_h), sharex=True, sharey=True
            )
            axes = np.array(axes).reshape(-1)

            bottom_margin_full_in = MARGIN_TB_IN + FULL_FIG_XLABEL_PAD_IN
            fig.subplots_adjust(
                left=MARGIN_LR_IN / fig_w,
                right=1 - MARGIN_LR_IN / fig_w,
                top=1 - MARGIN_TB_IN / fig_h,
                bottom=bottom_margin_full_in / fig_h,
                wspace=(WSPACE_IN / (len(tools_all) * _slot_in())) if cols > 1 else 0.2,
                hspace=(HSPACE_IN / panel_h_in) if rows > 1 else 0.25,
            )

            for i, (gene, g) in enumerate(
                sub_df.groupby("Gene", sort=True, observed=False)
            ):
                ax: plt.Axes = axes[i]
                g = g.set_index("Sign")
                n_bars = len(tools_all)
                x = np.arange(n_bars) + 0.5
                ax.set_xlim(0, n_bars)
                formatter = lambda value: format_score_label(
                    value, log_scale, two_line=True
                )
                for s in (0, 1):
                    bar = ax.bar(
                        x, g.loc[s, "score"], width_frac,
                        color=SIGN_COLORS[s], alpha=0.85,
                    )
                    annotate_bars(
                        ax, bar, formatter=formatter,
                        fontsize=6 if log_scale else 8,
                    )

                ax.axhline(0, linewidth=1)
                ax.set_title(f"{gene}", fontsize=9)

                row_idx = i // cols
                col_idx = i % cols
                show_bottom = row_idx == rows - 1 or i + cols >= total
                show_left = col_idx == 0
                configure_score_axis(
                    ax, fig, full_plot_values, log_scale,
                    show_ticklabels=show_left,
                    show_side_labels=show_left,
                )

                ax.set_xticks(x)
                if show_bottom:
                    ax.set_xticklabels(
                        tools_all, rotation=35, ha="right", fontsize=8
                    )
                else:
                    ax.tick_params(axis="x", labelbottom=False)
                ax.grid(axis="y", alpha=0.3)

            for j in range(total, len(axes)):
                axes[j].set_visible(False)

            fig.suptitle(f"Instance={inst}", y=0.995, fontsize=12)
            fig.supxlabel("Algorithm", y=0.01)
            fig.savefig(
                f"{opath}/{inst_group}/{inst}/_score_full{suffix}.{args.format}",
                dpi=200, bbox_inches="tight", pad_inches=0.3,
            )
            plt.close(fig)

            # -----------------------------------------------------------------
            # Summary-only (average over algorithms, 1 row figure)
            # -----------------------------------------------------------------
            fig_s, ax_sum = layout_single_axes(fig_w_s, fig_h_s)
            xg = np.arange(n_gene) + 0.5
            ax_sum.set_xlim(0, n_gene)
            for s in (0, 1):
                bars = ax_sum.bar(
                    xg, sum_by_gene.loc[s, "score"], width_frac,
                    color=SIGN_COLORS[s], alpha=0.85,
                )
                annotate_bars(
                    ax_sum, bars,
                    formatter=lambda value: format_score_label(
                        value, log_scale, two_line=True
                    ),
                    fontsize=7 if log_scale else 8,
                )

            configure_score_axis(
                ax_sum, fig_s, sum_by_gene["score"].to_numpy(), log_scale
            )
            ax_sum.axhline(0, linewidth=1)
            ax_sum.set_title(
                "Geometric mean over algorithms" if args.geo_mean
                else "Arithmetic mean over algorithms"
            )
            ax_sum.set_xlabel("Gene")
            ax_sum.set_xticks(xg, gene_sorted, rotation=45, ha="right")
            ax_sum.grid(axis="y", alpha=0.3)

            fig_s.suptitle(f"Instance={inst} — Summary", y=0.98, fontsize=12)
            fig_s.savefig(
                f"{opath}/{inst_group}/{inst}/_score_summary{suffix}.{args.format}",
                dpi=200, bbox_inches="tight", pad_inches=0.3,
            )
            plt.close(fig_s)

        # -----------------------------------------------------------------
        # Per-gene histogram grid (mirrors the _score_full layout): one
        # panel per gene, one bar per tool inside each panel. Positive
        # counts (sign=1) above zero, negative counts (sign=0) below.
        # -----------------------------------------------------------------
        inst_solution_df = solution_gene_df[solution_gene_df["Instance"] == inst]
        tools_list = list(tools_all)
        n_tools = len(tools_list)

        per_tool_counts = (
            inst_solution_df.groupby(
                ["Gene", "Algorithm", "Sign"], observed=True
            )
            .size()
            .unstack("Sign", fill_value=0)
        )
        for sign in (0, 1):
            if sign not in per_tool_counts.columns:
                per_tool_counts[sign] = 0
        per_tool_counts = per_tool_counts[[0, 1]]

        full_idx = pd.MultiIndex.from_product(
            [genes_order, tools_list], names=["Gene", "Algorithm"]
        )
        per_tool_counts = per_tool_counts.reindex(full_idx, fill_value=0)

        m_h = len(genes_order)
        total_h = m_h
        if m_h <= 3:
            rows_h, cols_h = m_h, 1
        else:
            cols_h = max(1, int(math.ceil(math.sqrt(total_h))))
            rows_h = int(math.ceil(total_h / cols_h))

        panel_h_in_h = 3.4
        fig_w_h, fig_h_h = _compute_figsize_grid(
            n_tools, rows_h, cols_h, panel_h_in=panel_h_in_h
        )
        fig_c, axes_h = plt.subplots(
            rows_h, cols_h, figsize=(fig_w_h, fig_h_h),
            sharex=True, sharey=True,
        )
        axes_h = np.array(axes_h).reshape(-1)

        bottom_margin_full_in = MARGIN_TB_IN + FULL_FIG_XLABEL_PAD_IN
        fig_c.subplots_adjust(
            left=MARGIN_LR_IN / fig_w_h,
            right=1 - MARGIN_LR_IN / fig_w_h,
            top=1 - MARGIN_TB_IN / fig_h_h,
            bottom=bottom_margin_full_in / fig_h_h,
            wspace=(WSPACE_IN / (n_tools * _slot_in())) if cols_h > 1 else 0.2,
            hspace=(HSPACE_IN / panel_h_in_h) if rows_h > 1 else 0.25,
        )

        all_pos = per_tool_counts[1].to_numpy(dtype=float)
        all_neg = -per_tool_counts[0].to_numpy(dtype=float)
        global_counts = (
            np.concatenate([all_pos, all_neg]) if m_h else np.array([0])
        )

        x_h = np.arange(n_tools) + 0.5
        for i, gene in enumerate(genes_order):
            ax: plt.Axes = axes_h[i]
            pos_h = per_tool_counts.loc[gene, 1].to_numpy(dtype=float)
            neg_h = -per_tool_counts.loc[gene, 0].to_numpy(dtype=float)
            ax.set_xlim(0, n_tools)
            pos_bars = ax.bar(
                x_h, pos_h, width_frac, color=SIGN_COLORS[1], alpha=0.85
            )
            neg_bars = ax.bar(
                x_h, neg_h, width_frac, color=SIGN_COLORS[0], alpha=0.85
            )
            annotate_bars(
                ax, pos_bars,
                formatter=lambda v: str(int(round(abs(v)))),
            )
            annotate_bars(
                ax, neg_bars,
                formatter=lambda v: str(int(round(abs(v)))),
            )

            ax.axhline(0, linewidth=1)
            ax.set_title(f"{gene}", fontsize=9)

            row_idx = i // cols_h
            col_idx = i % cols_h
            show_bottom = row_idx == rows_h - 1 or i + cols_h >= total_h
            show_left = col_idx == 0
            configure_count_axis(
                ax, fig_c, global_counts,
                show_ticklabels=show_left,
                show_side_labels=show_left,
            )

            ax.set_xticks(x_h)
            if show_bottom:
                ax.set_xticklabels(
                    tools_list, rotation=35, ha="right", fontsize=8
                )
            else:
                ax.tick_params(axis="x", labelbottom=False)
            ax.grid(axis="y", alpha=0.3)

        for j in range(total_h, len(axes_h)):
            axes_h[j].set_visible(False)

        fig_c.suptitle(f"Instance={inst} — Solution histogram",
                       y=0.995, fontsize=12)
        fig_c.supxlabel("Algorithm", y=0.01)
        fig_c.savefig(
            f"{opath}/{inst_group}/{inst}/_histogram_full.{args.format}",
            dpi=200, bbox_inches="tight", pad_inches=0.3,
        )
        plt.close(fig_c)


if __name__ == "__main__":
    main()
