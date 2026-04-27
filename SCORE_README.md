# Score outputs — reader's guide

This bundle holds per-instance scores and figures comparing how
multiple Boolean-network control tools agree (or disagree) on which
gene mutations control a given phenotype.

The central quantity is the **Mutation Co-occurrence Score (MCS)**: for
a single tool and a single mutation `(gene = 0/1)`, MCS is the
fraction of that tool's reported control sets in which the mutation
appears. MCS = 1 means every control set the tool produced uses that
mutation; MCS = 0 means none.

## Sign convention (used everywhere)

- **Activation** (gene fixed to 1) — **blue bars above zero**.
- **Inhibition** (gene fixed to 0) — **red bars below zero**.

## Files per instance

```
<group>/<instance>/
├── _histogram_sum_tool.<fmt>     ── how many control sets each tool produced
├── _histogram_full.<fmt>          ── per-gene × per-tool: how often the gene appears
├── _score_full.<fmt>              ── per-gene × per-tool MCS (linear y-axis)
├── _score_full_log_scale.<fmt>    ── same, symmetric-log y-axis
├── _score_summary.<fmt>           ── MCS aggregated across tools (linear)
└── _score_summary_log_scale.<fmt> ── same, symmetric-log
```

### `_histogram_sum_tool`

Bar chart of control set counts by tool, split by control set size
(1 vs 2). A quick view of which tools are prolific vs sparse.

### `_histogram_full`

Per-gene panel grid. Inside each panel, one bar per tool gives the
**number of that tool's control sets containing the gene** with the
given sign. Lets you see whether tools that score a gene high are also
including it in many of their control sets.

### `_score_full` (and `_log_scale`)

Per-gene panel grid of MCS bars, one bar per tool. The linear
version is best for a quick magnitude read; the log-scale version
exposes very small but non-zero scores that vanish on a linear axis.
Both files contain the same data.

### `_score_summary` (and `_log_scale`)

Single-panel summary collapsing the MCS bars across tools. Two
aggregation modes are available; the title indicates which one was
used:

- *Arithmetic mean over algorithms* — straightforward average.
- *Geometric mean over algorithms* — epsilon-shifted geometric mean,
  which collapses to ~0 whenever many tools score 0 for a gene.
  Helpful when scores span several orders of magnitude.

Genes are sorted by the chosen summary value, so when geo and arith
disagree on which gene is "biggest", the gene order can differ between
the two modes.

## `score.csv`

Long-format table, one row per (Instance, Tool, Gene, Sign):

| Column     | Meaning                                                                 |
|------------|-------------------------------------------------------------------------|
| `Instance` | Instance name.                                                          |
| `Tool`     | Tool short name (e.g. `BoNesis[FP]`).                                   |
| `Gene`     | Boolean Network variable.                                               |
| `Sign`     | `1` = activation, `0` = inhibition.                                     |
| `BN_size`  | Number of nodes in the source network.                                  |
| `score`    | MCS in `[0, 1]`.                                                        |

Genes a tool never reports get a constant default score, so every
(Instance, Tool, Gene, Sign) row exists in the table.
