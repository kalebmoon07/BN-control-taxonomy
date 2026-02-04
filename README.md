# Comparing control sets from Boolean network control software tools

This repository provides a framework to run various Boolean network control software tools on a collection of benchmark instances and compare their results. 

The main features are:

- A command line interface (CLI) to run multiple tools on multiple instances in batch mode.
- A **dominance graph** to compare the control sets found by each tool for each instance.
  - If every control found by tool A has another subset control found by tool B, then tool B is said to _dominate_ tool A for that instance.
- The **Mutation Co-occurrence Score** (MCS) for each mutation $(x_i=b)$ predicted by each tool
  - The fraction of the controls containing each mutation are dominated by a control set found by each tool.
  - This metric is consistent with the dominance graph: if tool B dominates tool A, then tool B has higher or equal MCS scores than tool A for all mutations.
  - An ensemble prediction can be made as the average MCS score across a selection of tools.  

## Installation

### Using CoLoMoTo Docker distribution

This is the recommended usage.

You need [`colomoto-docker`](https://github.com/colomoto/colomoto-docker-py) command, see also <https://colomoto.github.io/colomoto-docker>.

### General usage without CoLoMoTo Docker

You need a recent Python. Required dependencies for _result analysis_ are listed in `requirements.txt`, and can be installed using `pip`:

```sh
pip install -r requirements.txt
```

**Important** Each tool must be properly installed, or it will be excluded from the analysis. Most tools can be installed via `conda`, but some may require additional steps.

```sh
conda install --file requirements_tools.txt
```

## Running experiments

General usage:

```sh
python src/bntaxonomy/cli.py {max_size} (--instances PATH [PATH ...] | --inst_groups PATH [PATH ...]) [--tools TOOL [TOOL ...]] [--exclude-targets]
```

**Recommended** Usage within CoLoMoTo Docker distribution:

```sh
colomoto-docker --bind . python src/bntaxonomy/cli.py {max_size} (--instances PATH [PATH ...] | --inst_groups PATH [PATH ...]) [--tools TOOL [TOOL ...]] [--exclude-targets]
```

To analyze the results, use the following commands:

```sh
# run all available tools on B_manually_designed instances with max_size=2
python src/bntaxonomy/cli.py 2 --inst_groups experiments/instances/B_manually_designed

# same, but run only BoNesis[FP] and BoNesis[MTS] tools
python src/bntaxonomy/cli.py 2 --inst_groups experiments/instances/B_manually_designed --tools 'BoNesis[FP]' 'BoNesis[MTS]'

# summarize the results for all instances and make dominance graphs under experiments/results:
python src/bntaxonomy/summarize.py

# compute scores and generating plots for all instances to experiments/results/fig:
python src/bntaxonomy/evaluate_score.py -o experiments/results/fig
```

Outputs are generated under `experiments/results/`, including:

- The dominance graph in `/_summary_tred.dot`.
  - Each arc indicates that the source tool dominates the target tool on all instances.
  - The graph for each instance is located at `/<group>/<instance>/_graph_tred.dot`.
  - Counterexamples for the dominance relations are listed in `counterexamples_first_match.csv` and `counterexamples_full_match.csv`.
- The MCS scores are available as CSV files in `scores.csv`.
- Figures plotting the MCS scores made by `evaluate_score.py` (not synchronized).

## References

If you use this framework in your research, please cite the following paper:

TODO

## Acknowledgements

TODO

## Contributing to the project

See [CONTRIBUTE.md](CONTRIBUTE.md).

## Detailed CLI usage

See [CLI_COMMAND.md](CLI_COMMAND.md).
