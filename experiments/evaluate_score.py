#!/usr/bin/env python3
"""
Evaluate and plot control scores for Boolean network control experiments.
Usage:
    python experiments/evaluate_score.py [options]
"""
import argparse
import math
import sys
import os
from collections import OrderedDict

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from bntaxonomy.hierarchy import MultiInputSummary
from automate_test import insts_practical, insts_ce


# ---------------------------------------------------------------------
# Visual constants (inches) — keep bars/pads uniform across inputs
# ---------------------------------------------------------------------
BAR_W_IN = 0.22  # width of each bar (inches)
GAP_IN = 0.12  # gap between neighboring bars (inches)
MARGIN_LR_IN = 0.6  # left/right figure margins (inches)
MARGIN_TB_IN = 0.6  # top/bottom figure margins (inches)
WSPACE_IN = 0.35  # inter-subplot horizontal spacing (inches)
HSPACE_IN = 0.45  # inter-subplot vertical spacing (inches)

SIGN_COLORS = {1: "tab:blue", 0: "tab:red"}
SIGN_NAME = {1: "Positive", 0: "Negative"}


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


def sort_by_lex_score(df_inst: pd.DataFrame):
    # collapse to one row per (Gene, Sign)
    g: pd.DataFrame = (
        df_inst.reset_index()
        .groupby(["Gene", "Sign"], as_index=False, observed=True)["score"]
        .sum()
    )

    # wide form with separate positive/negative columns
    wide = g.pivot(index="Gene", columns="Sign", values="score").fillna(0.0)
    wide = wide.rename(columns={1: "pos", 0: "neg_raw"})
    wide["neg"] = wide["neg_raw"].abs()
    wide["total"] = wide["pos"] + wide["neg"]

    # final sort: primary by total (desc), secondary by positive (desc)
    gene_sorted = wide.sort_values(
        by=["total", "pos"], ascending=[False, False], kind="mergesort"
    ).index.to_list()
    return gene_sorted


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


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main(argv=None):
    cwd = os.path.dirname(os.path.abspath(__file__))

    parser = argparse.ArgumentParser(description="Compute/plot control scores.")
    parser.add_argument(
        "-a",
        "--algorithms",
        "--algs",
        nargs="+",
        help="Algorithms to include. If omitted, include all.",
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
        help=(
            "List of specific instances to include."
        ),
        default=None,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output directory for results.",
        default=f"{cwd}/results/fig_score",
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
    args = parser.parse_args(argv)

    if args.genes:
        selected_gene_order = list(dict.fromkeys(args.genes))
    else:
        selected_gene_order = None

    # selected_algs = set(args.algorithms) if args.algorithms else None
    inst_groups = [p.replace("instances", "results") for p in args.inst_groups]
    instances: list[str] = args.instances

    os.makedirs(args.output, exist_ok=True)
    opath = args.output

    # -----------------------------------------------------------------
    # Summaries from experiment groups
    # -----------------------------------------------------------------
    if inst_groups:
        hc = MultiInputSummary.from_inst_groups(inst_groups, "Hierarchy")
    elif args.instances:
        hc = MultiInputSummary.from_instances(instances, "Hierarchy")
    else:
        raise ValueError("One of --inst_groups or --instances must be provided.")

    if args.algorithms:
        available = {r.name for e in hc.exp_list for r in e.results}
        missing = [a for a in args.algorithms if a not in available]
        for m in missing:
            print(
                f"Warning: algorithm '{m}' not found in {inst_groups}; skipping.",
                file=sys.stderr,
            )
        # keep only those present, in the original order
        selected_algs = [a for a in args.algorithms if a in available]
    else:
        selected_algs = None
    cat_type = pd.CategoricalDtype(selected_algs, ordered=True)

    # Build records (tight inner loops with local bindings for speed)
    count_list = []  # (Instance, Algorithm, BN_size, ControlSize, Genes)
    gene_ctrl_list = (
        []
    )  # (Instance, Algorithm, Gene, Sign, BN_size, ControlSize, Partners)
    for exp in hc.exp_list:
        inst_name = exp.name
        n_nodes = len(exp.bn)
        bn_keys = tuple(exp.bn.keys())  # stable order once

        for alg_result in exp.results:
            alg_name = alg_result.name
            if selected_algs is not None and alg_name not in selected_algs:
                continue

            if not alg_result.d_list:
                count_list.append((inst_name, alg_name, n_nodes, math.inf, dict()))

            ctrl_set = set()
            for ctrl_dict in alg_result.d_list:
                genes = list(ctrl_dict.keys())
                csize = len(ctrl_dict)
                count_list.append((inst_name, alg_name, n_nodes, csize, ctrl_dict))

                # partner tuples are sorted to be hash-stable
                for gene, val in ctrl_dict.items():
                    ctrl_set.add((gene, val))
                    partners = tuple(sorted(g for g in genes if g != gene))
                    gene_ctrl_list.append(
                        (inst_name, alg_name, gene, val, n_nodes, csize, partners)
                    )

            # Fill missing (gene, sign) with ControlSize=-1
            for gene in bn_keys:
                for sign in (0, 1):
                    if (gene, sign) not in ctrl_set:
                        gene_ctrl_list.append(
                            (inst_name, alg_name, gene, sign, n_nodes, -1, ())
                        )

    # -----------------------------------------------------------------
    # Histogram (per Instance)
    # NOTE: Histogram keeps overall control-size counts after filtering control sets
    #       to the selected genes (if provided).
    # -----------------------------------------------------------------
    count_df = pd.DataFrame(
        count_list, columns=["Instance", "Algorithm", "BN_size", "ControlSize", "Genes"]
    )
    if selected_algs:
        count_df["Algorithm"] = count_df["Algorithm"].astype(cat_type)
    count_df.to_csv(f"{opath}/count_control.csv", index=False)

    for inst, sub_df in count_df.groupby("Instance", sort=False):
        ct1 = pd.crosstab(sub_df["Algorithm"], sub_df["ControlSize"]).sort_index()
        # Ensure columns for both control sizes exist (1, 2)
        full_sizes = [1, 2]
        ct1 = ct1.reindex(columns=full_sizes, fill_value=0)

        algs = ct1.index.to_list()
        sizes = full_sizes

        x = np.arange(len(algs))
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

        ax.set_xticks(x, algs, rotation=45)
        ax.grid(True, axis="y", alpha=0.25)
        ax.legend(title="ControlSize", loc="center left", bbox_to_anchor=(1.01, 0.5))
        ax.set_ylim(0, max(1, int(ct1.to_numpy().max())) * 1.15)  # headroom for labels)
        plt.subplots_adjust(bottom=0.22, right=0.85)

        fig.savefig(f"{opath}/{inst}_histogram.png", dpi=200, bbox_inches="tight")
        plt.close(fig)

    # -----------------------------------------------------------------
    # Score calculations
    # -----------------------------------------------------------------
    score_df = pd.DataFrame(
        gene_ctrl_list,
        columns=[
            "Instance",
            "Algorithm",
            "Gene",
            "Sign",
            "BN_size",
            "ControlSize",
            "Partners",
        ],
    )
    score_df.drop_duplicates(
        ["Instance", "Algorithm", "Gene", "Sign", "Partners"], inplace=True
    )
    if selected_algs:
        score_df["Algorithm"] = score_df["Algorithm"].astype(cat_type)

    def compute_score(row):
        if row["ControlSize"] == -1:
            return 0.0
        if row["ControlSize"] == 0:
            return 1.0
        else:
            score = 1 / math.comb(row["BN_size"] - 1, row["ControlSize"] - 1)
            return score if row["Sign"] == 1 else -score

    score_df["score"] = score_df.apply(compute_score, axis=1).fillna(0.0)
    score_df.to_csv(f"{opath}/score_gene_level.csv", index=False)

    # -----------------------------------------------------------------
    # Plotting: per-Instance grids (one panel per gene)
    # -----------------------------------------------------------------
    width_frac = _bar_width_frac()

    for inst, sub_df in score_df.groupby("Instance", sort=False):
        sub_df: pd.DataFrame = sub_df.groupby(
            ["Algorithm", "Gene", "Sign"], as_index=False, observed=True
        )["score"].sum()

        algs_all = sub_df["Algorithm"].unique()
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
            rows, cols = 1, m
        else:
            cols = max(1, int(math.ceil(math.sqrt(total))))
            rows = int(math.ceil(total / cols))

        # Figure size derived from desired physical bar/pad sizes
        fig_w, fig_h = _compute_figsize_grid(len(algs_all), rows, cols, panel_h_in=3.6)
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
            wspace=(WSPACE_IN / (len(algs_all) * _slot_in())) if cols > 1 else 0.2,
            hspace=(HSPACE_IN / 3.6) if rows > 1 else 0.25,
        )

        # Per-gene panels
        for i, (gene, g) in enumerate(
            sub_df.groupby("Gene", sort=True, observed=False)
        ):
            ax: plt.Axes = axes[i]
            g = g.set_index("Sign")
            n_bars = len(algs_all)
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
            ax.set_ylabel("Score")
            ax.set_yticks([-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1])
            ax.set_xticks(x, algs_all, rotation=45, ha="right")
            ax.grid(axis="y", alpha=0.3)

        # Hide any unused axes
        for j in range(total, len(axes)):
            axes[j].set_visible(False)

        fig.suptitle(f"Instance={inst}", y=0.995, fontsize=12)
        plt.tight_layout()
        plt.savefig(f"{opath}/{inst}_full.png", dpi=200, bbox_inches="tight")
        plt.close(fig)

        # -----------------------------------------------------------------
        # Summary-only (average over algorithms, 1 row figure)
        # -----------------------------------------------------------------

        sum_by_gene = (
            sub_df.groupby(["Gene", "Sign"], observed=True)["score"]
            .mean()
            .reset_index()
        )

        gene_sorted = sort_by_lex_score(sub_df)
        sum_by_gene["Gene"] = pd.Categorical(
            sum_by_gene["Gene"], categories=gene_sorted, ordered=True
        )
        sum_by_gene = sum_by_gene.sort_values("Gene")
        sum_by_gene = sum_by_gene.set_index("Sign")

        n_gene = len(genes_order)
        content_gene_w = max(1, n_gene) * _slot_in()
        fig_w_s = content_gene_w + 2 * MARGIN_LR_IN
        fig_h_s = 6

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
        ax_sum.axhline(0, linewidth=1)
        ax_sum.set_title("Average over algorithms")
        ax_sum.set_ylabel("Score")
        ax_sum.set_xlabel("Gene")
        ax_sum.set_xticks(xg, gene_sorted, rotation=45, ha="right")
        ax_sum.grid(axis="y", alpha=0.3)

        fig_s.suptitle(f"Instance={inst} — Summary", y=0.98, fontsize=12)
        plt.savefig(f"{opath}/{inst}_summary.png", dpi=200, bbox_inches="tight")
        plt.close(fig_s)


if __name__ == "__main__":
    main()
