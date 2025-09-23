import math
from bntaxonomy.hierarchy import MultiInputSummary
import os
from automate_test import insts_practical, insts_ce
import pandas as pd

from bntaxonomy.hierarchy import MultiInputSummary
import numpy as np
import matplotlib.pyplot as plt

cwd = os.path.dirname(os.path.abspath(__file__))

inst_groups = [
    result_dir.replace("instances", "results")
    for result_dir in [insts_practical, insts_ce]
]

hc = MultiInputSummary.from_inst_groups(inst_groups, "Hierarchy")

count_list: list[tuple[str, str, str, int, float]] = []
"""inst_name, alg_name, gene, sign, n, score"""

for exp in hc.exp_list:
    # Instance name
    inst_name = exp.name
    n = len(exp.bn)  # number of nodes in the BN

    for alg_result in exp.results:
        # item: CtrlResult for one method
        alg_name = alg_result.name
        for ctrl_dict in alg_result.d_list:
            # ctrl_dict: one control strategy
            for gene, value in ctrl_dict.items():
                count_list.append(
                    (
                        inst_name,
                        alg_result.name,
                        gene,
                        value,
                        n,
                        len(ctrl_dict),
                    )
                )
# print(count_list)
score_df = pd.DataFrame(
    count_list, columns=["Instance", "Algorithm", "Gene", "Sign", "BN_size", "Score"]
)
score_df: pd.DataFrame = (
    score_df.groupby(["Instance", "Algorithm", "Gene", "Sign", "BN_size"])
    .agg(
        List=pd.NamedAgg(column="Score", aggfunc=lambda x: list(x)),
        Count=pd.NamedAgg(column="Score", aggfunc="count"),
        AvgSize=pd.NamedAgg(column="Score", aggfunc=lambda x: sum(x) / (len(x))),
    )
    .reset_index()
)
# count_df["Score"] = count_df[["List","BN_size"]].apply(lambda x: sum(x))
score_df["WeightedSum"] = score_df.apply(
    lambda row: sum(1 / ((row["BN_size"]) ** (x - 1)) for x in row["List"]), axis=1
)
score_df["WeightedSumHalf"] = score_df.apply(
    lambda row: sum(1 / ((row["BN_size"] / 2) ** (x - 1)) for x in row["List"]), axis=1
)
score_df["WeightedSum2Based"] = score_df.apply(
    lambda row: sum(1 / (2 ** (x - 1)) for x in row["List"]), axis=1
)
# count_df["avg_partner"] = count_df["WeightedSum"].apply(lambda x: -log2(x))
score_df.sort_values(by=["Instance", "Gene", "Sign", "Algorithm"], inplace=True)
score_df.to_csv(f"{cwd}/control_score.csv", index=False)


def make_wide_by_gene(df_inst, score_col="WeightedSum"):
    """
    From rows for a single Instance, aggregate <score_col> across BN_size and
    return a wide table indexed by (Gene, Algorithm) with columns:
      - 'pos'  = value for Sign=1  (as-is)
      - 'neg'  = negated value for Sign=0  (so it's <= 0)
    Missing Algorithms are filled with 0.
    """
    agg = df_inst.groupby(["Gene", "Algorithm", "Sign"], as_index=False)[
        score_col
    ].sum()

    wide = agg.pivot_table(
        index=["Gene", "Algorithm"],
        columns="Sign",
        values=score_col,
        aggfunc="sum",
        fill_value=0.0,
    )

    for s in (0, 1):
        if s not in wide.columns:
            wide[s] = 0.0

    wide = wide.rename(columns={1: "pos", 0: "neg_raw"})
    wide["neg"] = -wide["neg_raw"]  # make neg <= 0
    wide = wide.drop(columns=["neg_raw"]).sort_index()

    return wide


def annotate_bars(ax, bars, fmt="{:.2f}", dy=0.01):
    """
    Put value labels at the end of each bar.
    - For positive bars, slightly above top.
    - For negative bars, slightly below top (which is a negative number).
    """
    ylims = ax.get_ylim()
    span = ylims[1] - ylims[0]
    offset = dy * span

    for b in bars:
        h = b.get_height()
        if h == 0:
            continue
        x = b.get_x() + b.get_width() / 2.0
        y = b.get_y() + h
        # choose vertical offset based on sign of bar
        va = "bottom" if h > 0 else "top"
        yo = offset if h > 0 else -offset
        ax.text(x, y + yo, fmt.format(h), ha="center", va=va, fontsize=8, rotation=0)


# ---- plotting: one figure per Instance, subplots per Gene + one summary ----
for inst, sub_df in score_df.groupby("Instance"):
    wide = make_wide_by_gene(sub_df, score_col="WeightedSum")

    genes = list(dict.fromkeys(wide.index.get_level_values("Gene")))  # preserve order
    algs_all = sorted(sub_df["Algorithm"].unique())  # consistent x
    m = len(genes)

    # total panels = per-gene panels + 1 summary panel
    total = m
    cols = math.ceil(math.sqrt(total))
    rows = math.ceil(total / cols)

    fig, axes = plt.subplots(
        rows, cols, figsize=(5 * cols, 5 * rows), sharex=False, sharey=True
    )
    axes = np.array(axes).reshape(-1)

    width = 0.38  # bar width

    # ---- per-gene panels ----
    for i, gene in enumerate(genes):
        ax = axes[i]

        g = (
            wide.loc[gene]
            if gene in wide.index.get_level_values("Gene")
            else wide.iloc[0:0]
        )
        g = g.reindex(algs_all, fill_value=0.0)

        x = np.arange(len(algs_all))
        bars_pos = ax.bar(x, g["pos"].to_numpy(), width, label="Positive")
        bars_neg = ax.bar(x, g["neg"].to_numpy(), width, label="Negative")

        ymax, ymin = 1.0, -1.0
        pad = 0.1 * (ymax - ymin + 1e-9)
        ax.set_ylim(ymin - pad, ymax + pad)

        ax.axhline(0, linewidth=1)
        ax.set_title(f"Gene={gene}")
        ax.set_ylabel("Score (sum 1/n^{|C|-1})")
        ax.set_xticks(x, algs_all, rotation=45, ha="right")
        ax.grid(axis="y", alpha=0.3)

        # annotate values
        annotate_bars(ax, bars_pos, fmt="{:.2f}")
        annotate_bars(ax, bars_neg, fmt="{:.2f}")

        if i == 0:
            ax.legend()

    # hide any unused axes
    for j in range(total, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(f"Instance={inst}", y=0.995, fontsize=12)
    plt.tight_layout()
    plt.savefig(f"{cwd}/fig_score_full_{inst}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # Compute the two summaries
    sum_by_alg = (
        wide.reset_index()
        .groupby("Algorithm", as_index=True)[["pos", "neg"]]
        .sum()
        .reindex(algs_all, fill_value=0.0)
    )
    sum_by_gene = (
        wide.reset_index()
        .groupby("Gene", as_index=True)[["pos", "neg"]]
        .sum()
        .reindex(genes, fill_value=0.0)
    )

    # ---- new figure with summary blocks only ----
    width = 0.38
    fig_s, (ax_sum_alg, ax_sum_gene) = plt.subplots(
        1, 2, figsize=(max(8, 1.5 * len(algs_all) + 4), 6), sharey=False
    )

    # Summary #1: Σ over genes per Algorithm
    x = np.arange(len(algs_all))
    bars_pos_s1 = ax_sum_alg.bar(
        x, sum_by_alg["pos"].to_numpy(), width, label="Positive (Σ genes)"
    )
    bars_neg_s1 = ax_sum_alg.bar(
        x, sum_by_alg["neg"].to_numpy(), width, label="Negative (Σ genes)"
    )

    ymax = float(max(0.0, sum_by_alg["pos"].max())) if len(sum_by_alg) else 0.0
    ymin = float(min(0.0, sum_by_alg["neg"].min())) if len(sum_by_alg) else 0.0
    pad = 0.1 * (ymax - ymin + 1e-9)
    ax_sum_alg.set_ylim(ymin - pad, ymax + pad)

    ax_sum_alg.axhline(0, linewidth=1)
    ax_sum_alg.set_title("Summary: Σ over genes (per Algorithm)")
    ax_sum_alg.set_ylabel("Score (sum 1/n^{|C|-1})")
    ax_sum_alg.set_xticks(x, algs_all, rotation=45, ha="right")
    ax_sum_alg.grid(axis="y", alpha=0.3)
    annotate_bars(ax_sum_alg, bars_pos_s1, fmt="{:.2f}")
    annotate_bars(ax_sum_alg, bars_neg_s1, fmt="{:.2f}")
    ax_sum_alg.legend()

    # Summary #2: Σ over algorithms per Gene
    xg = np.arange(len(genes))
    bars_pos_s2 = ax_sum_gene.bar(
        xg, sum_by_gene["pos"].to_numpy(), width, label="Positive (Σ algs)"
    )
    bars_neg_s2 = ax_sum_gene.bar(
        xg, sum_by_gene["neg"].to_numpy(), width, label="Negative (Σ algs)"
    )

    ymax = float(max(0.0, sum_by_gene["pos"].max())) if len(sum_by_gene) else 0.0
    ymin = float(min(0.0, sum_by_gene["neg"].min())) if len(sum_by_gene) else 0.0
    pad = 0.1 * (ymax - ymin + 1e-9)
    ax_sum_gene.set_ylim(ymin - pad, ymax + pad)

    ax_sum_gene.axhline(0, linewidth=1)
    ax_sum_gene.set_title("Summary: Σ over algorithms (per Gene)")
    ax_sum_gene.set_ylabel("Score (sum 1/n^{|C|-1})")
    ax_sum_gene.set_xticks(xg, genes, rotation=45, ha="right")
    ax_sum_gene.grid(axis="y", alpha=0.3)
    annotate_bars(ax_sum_gene, bars_pos_s2, fmt="{:.2f}")
    annotate_bars(ax_sum_gene, bars_neg_s2, fmt="{:.2f}")
    ax_sum_gene.legend()

    fig_s.suptitle(f"Instance={inst} — Summary", y=0.98, fontsize=12)
    plt.tight_layout()
    plt.savefig(f"{cwd}/fig_score_summary_{inst}.png", dpi=200, bbox_inches="tight")
    plt.close(fig_s)
