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

SIGN_COLORS = {1: "tab:blue", 0: "tab:red", -1: "tab:red"}
SIGN_NAME = {1: "Positive", 0: "Negative", -1: "Negative"}


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
def annotate_bars(ax, bars, fmt="{:.2f}", dy=0.01):
    """
    Put value labels at the end of each bar.
    - For positive bars, slightly above top.
    - For negative bars, slightly below top (which is a negative number).
    """
    y0, y1 = ax.get_ylim()
    offset = dy * (y1 - y0)

    for b in bars:
        h = b.get_height()
        if not h:  # skip zero-height bars
            continue
        x = b.get_x() + b.get_width() / 2.0
        y = b.get_y() + h
        va = "bottom" if h > 0 else "top"
        yo = offset if h > 0 else -offset
        ax.text(x, y + yo, fmt.format(abs(h)), ha="center", va=va, fontsize=8)


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
    # Put labels just left of the y-axis with fixed physical padding
    trans = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)

    # x = 0 is the y-axis in axes coords; move left by 18 points
    pad = mtransforms.ScaledTranslation(-0.55, 0, fig.dpi_scale_trans)

    ax.text(
        0,
        0.5,
        "Activation",
        transform=trans + pad,
        rotation=90,
        va="center",
        ha="center",
        clip_on=False,
    )
    ax.text(
        0,
        -0.5,
        "Inhibition",
        transform=trans + pad,
        rotation=90,
        va="center",
        ha="center",
        clip_on=False,
    )


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
    for k in topk:
        if k > len(pivot):
            continue
        idx = ens.head(k).index
        sub = pivot.loc[idx]
        ens_sub = ens.loc[idx]
        # Explicit fractional-rank tie handling: Spearman via Pearson on
        # average-method ranks, so tied MCS values share the mean rank.
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
    parser.add_argument(
        "--topk",
        nargs="+",
        type=int,
        default=[3, 5, 10, 15, 20, 30, 60],
        help=(
            "k values for the top-k Spearman vs ensemble heatmap. "
            "Entries exceeding an instance's gene count are skipped."
        ),
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
    if selected_tools:
        count_df["Algorithm"] = count_df["Algorithm"].astype(cat_type)
        mcs_score_df["Algorithm"] = mcs_score_df["Algorithm"].astype(cat_type)
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
        _write_family_partition(ctrls, families, mcs_outdir)
        _write_family_jaccard(ctrls, families, mcs_outdir, args.format, inst)
        for sign in (0, 1):
            _write_topk_spearman_vs_ensemble(
                inst_score_df, families, mcs_outdir, sign,
                args.topk, args.geo_mean, args.geo_eps, args.format, inst,
            )

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
            f"{opath}/{inst_group}/{inst}/_score_histogram.{args.format}",
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

        # Figure size derived from desired physical bar/pad sizes
        fig_w, fig_h = _compute_figsize_grid(len(tools_all), rows, cols, panel_h_in=3.6)
        fig, axes = plt.subplots(
            rows, cols, figsize=(fig_w, fig_h), sharex=False, sharey=False
        )
        axes = np.array(axes).reshape(-1)

        # Convert inch spacings into figure-fractions for precise layout
        fig.subplots_adjust(
            left=MARGIN_LR_IN / fig_w,
            right=1 - MARGIN_LR_IN / fig_w,
            top=1 - MARGIN_TB_IN / fig_h,
            bottom=MARGIN_TB_IN / fig_h,
            wspace=(WSPACE_IN / (len(tools_all) * _slot_in())) if cols > 1 else 0.2,
            hspace=(HSPACE_IN / 3.6) if rows > 1 else 0.25,
        )

        # Per-gene panels
        for i, (gene, g) in enumerate(
            sub_df.groupby("Gene", sort=True, observed=False)
        ):
            ax: plt.Axes = axes[i]
            g = g.set_index("Sign")
            n_bars = len(tools_all)
            x = np.arange(n_bars) + 0.5  # centers at 0.5, 1.5, ...
            ax.set_xlim(0, n_bars)  # 1 data unit == one (bar+gap) slot
            for s in (0, 1):
                bar = ax.bar(
                    x, g.loc[s, "score"], width_frac, color=SIGN_COLORS[s], alpha=0.85
                )
                annotate_bars(ax, bar, fmt="{:.2f}")

            # symmetric y-limit with padding
            ax.set_ylim(-1.15, 1.15)

            ax.axhline(0, linewidth=1)
            ax.set_title(f"{gene}")
            ax.set_ylabel(" ")
            annotate_score_axis(fig, ax)

            ax.set_yticks([-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1])
            ax.set_yticklabels([1, 0.75, 0.5, 0.25, 0, 0.25, 0.5, 0.75, 1])
            ax.yaxis.set_major_formatter(lambda x, pos: f"{abs(x):.2f}")

            ax.set_xticks(x, tools_all, rotation=45, ha="right")
            ax.grid(axis="y", alpha=0.3)

        # Hide any unused axes
        for j in range(total, len(axes)):
            axes[j].set_visible(False)

        fig.suptitle(f"Instance={inst}", y=0.995, fontsize=12)
        plt.tight_layout()
        os.makedirs(f"{opath}/{inst_group}/{inst}", exist_ok=True)
        plt.savefig(
            f"{opath}/{inst_group}/{inst}/_score_full.{args.format}", dpi=200, bbox_inches="tight"
        )
        plt.close(fig)

        # -----------------------------------------------------------------
        # Summary-only (average over algorithms, 1 row figure)
        # -----------------------------------------------------------------

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

        # gene_sorted = sort_by_neg_score(sub_df)
        if args.sort == "total":
            gene_sorted = sort_by_total_score(sub_df)
        elif args.sort == "pos":
            gene_sorted = sort_by_pos_score(sub_df)
        else:  # args.sort == "neg"
            gene_sorted = sort_by_neg_score(sub_df)
        sum_by_gene["Gene"] = pd.Categorical(
            sum_by_gene["Gene"], categories=gene_sorted, ordered=True
        )
        sum_by_gene = sum_by_gene.sort_values("Gene")
        sum_by_gene = sum_by_gene.set_index("Sign")

        n_gene = len(genes_order)
        content_gene_w = max(1, n_gene) * _slot_in()
        fig_w_s = content_gene_w + 2 * MARGIN_LR_IN
        fig_h_s = 4

        fig_s, ax_sum = layout_single_axes(fig_w_s, fig_h_s)

        xg = np.arange(n_gene) + 0.5
        ax_sum.set_xlim(0, n_gene)
        for s in (0, 1):
            bars = ax_sum.bar(
                xg,
                sum_by_gene.loc[s, "score"],
                width_frac,
                color=SIGN_COLORS[s],
                alpha=0.85,
            )
            annotate_bars(ax_sum, bars, fmt="{:.2f}")

        ax_sum.set_ylim(-1.15, 1.15)
        ax_sum.set_yticks([-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1])
        ax_sum.set_yticklabels([1, 0.75, 0.5, 0.25, 0, 0.25, 0.5, 0.75, 1])
        ax_sum.yaxis.set_major_formatter(lambda x, pos: f"{abs(x):.2f}")
        ax_sum.axhline(0, linewidth=1)
        ax_sum.set_title(
            "Geometric mean over algorithms" if args.geo_mean
            else "Arithmetic mean over algorithms"
        )
        ax_sum.set_ylabel(" ")
        ax_sum.set_xlabel("Gene")
        ax_sum.set_xticks(xg, gene_sorted, rotation=45, ha="right")
        ax_sum.grid(axis="y", alpha=0.3)
        
        annotate_score_axis(fig_s, ax_sum)
        fig_s.suptitle(f"Instance={inst} — Summary", y=0.98, fontsize=12)
        os.makedirs(f"{opath}/{inst_group}/{inst}", exist_ok=True)
        suffix = "_geo" if args.geo_mean else ""
        fig_s.savefig(
            f"{opath}/{inst_group}/{inst}/_score_summary{suffix}.{args.format}",
            dpi=200,
            bbox_inches="tight",
        )
        plt.close(fig_s)


if __name__ == "__main__":
    main()
