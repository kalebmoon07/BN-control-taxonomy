# Comparing control sets from Boolean network control software tools

## Installation

### Using CoLoMoTo Docker distribution

This is the recommended usage.

You need [`colomoto-docker`](https://github.com/colomoto/colomoto-docker-py) command, see also <https://colomoto.github.io/colomoto-docker>.

### General usage

You need a recent Python. Required dependences are listed in `requirements.txt`, and can be installed using `pip`:

```sh
pip install -r requirements.txt
```

**Important** Each tool must be properly installed, or it will be excluded from the analysis.

## Running experiments

General usage:

```sh
python src/bntaxonomy/cli.py {max_size} {instances...} [--tools {tools...}]
```

**Recommended** Usage within CoLoMoTo Docker distribution:

```sh
colomoto-docker --bind . python src/bntaxonomy/cli.py {max_size} {instances...} [--tools {tools...}]
```

## Examples

```sh
# run all available tools on B_manually_designed instances with max_size=2
python src/bntaxonomy/cli.py 2 experiments/instances/B_manually_designed/* --tools 'BoNesis[FP]' 'BoNesis[MTS]'

# same, but run only BoNesis[FP] and BoNesis[MTS] tools
python src/bntaxonomy/cli.py 2 experiments/instances/B_manually_designed/* --tools 'BoNesis[FP]' 'BoNesis[MTS]'

# summarize results and make dominance graphs:
python src/bntaxonomy/summarize.py

# compute scores and generating plots (inside CoLoMoTo Docker):
python src/bntaxonomy/evaluate_score.py -ig experiments/instances/B_manually_designed -o experiments/results/fig
```

## Acknowledgements

TODO

## Detailed CLI usage

The main command is `src/bntaxonomy/cli.py`. The script requires a single positional argument `max_size` (an integer) and lets you select which instances to run on using either `--instances` or `--instance-groups`.

Options (as implemented in `src/bntaxonomy/cli.py`):

- `max_size` (positional): maximum number of perturbations to consider (required).
- `--instances PATH [PATH ...]`: one or more instance folders (each folder must live under an `instances` directory).
- `--instance-groups PATH [PATH ...]`: one or more directories that contain instance subfolders (the CLI will add every subfolder found).
- `--tools TOOL [TOOL ...]`: restrict which analysis tools to run (use the tool short names, e.g. `BoNesis[FP]`).
- `--exclude-targets`: exclude nodes that specify the target phenotype from candidate perturbations.

Behaviour and output:

- For each input instance folder `.../instances/<group>/<instance>` the CLI writes results under the corresponding `.../results/<group>/<instance>` path (the code replaces the `instances` path component with `results`).
- The CLI saves per-instance JSON result files and a `_graph` summary (`.dot`/`.png`) inside each results folder.

Examples (correct usage matching the current code):

```sh
# Run specific instance(s):
python src/bntaxonomy/cli.py 2 --instances experiments/instances/B_manually_designed/ce_long_attr --tools 'BoNesis[FP]' 'BoNesis[MTS]'

# Run all instances inside a group directory:
python src/bntaxonomy/cli.py 2 --instance-groups experiments/instances/B_manually_designed --tools 'BoNesis[FP]' 'BoNesis[MTS]'

# Recommended (inside CoLoMoTo Docker):
colomoto-docker --bind . python src/bntaxonomy/cli.py 2 --instance-groups experiments/instances/A_case_studies/Bladder
```

## Summarize into a dominance graph and detect B_manually_designed (`src/bntaxonomy/summarize.py`)

Generate overall conflict matrices and grouped lists of B_manually_designed.

Basic usage:

```sh
python src/bntaxonomy/summarize.py [options]
```

Important options:

- `-ig`, `--inst_groups` PATH [PATH ...]: instance-group directories under `instances` (they will be mapped to `results`).
- `-i`, `--instances` PATH [PATH ...]: explicit instance folders under `instances`.

Example:

```sh
# Summarize a pair of instance groups and write CSV/graph outputs to experiments/results
python src/bntaxonomy/summarize.py -ig experiments/instances/A_case_studies

# Summarize explicit instances
python src/bntaxonomy/summarize.py -i experiments/instances/B_manually_designed/ce_long_attr experiments/instances/B_manually_designed/ce_yes_P_no_R
```

## Evaluating scores (`src/bntaxonomy/evaluate_score.py`)

This helper computes CSV summaries and generates plots from the per-instance results produced by the CLI.

Basic usage:

```sh
python src/bntaxonomy/evaluate_score.py [options]
```

Important options (see the script for full details):

- `-ig`, `--inst_groups` PATH [PATH ...]: one or more instance-group paths you want to evaluate. These paths should point at the a directory containing multiple `instances` subdirectories (e.g. `experiments/instances/A_case_studies`); the script will automatically translate `instances` -> `results` to find the results folders produced by the CLI.
- `-i`, `--instances` PATH [PATH ...]: explicit instance folders (under `instances`); the script will look for the corresponding `results` folders.
- `-a`, `--algorithms` ALG [ALG ...]: restrict to a subset of algorithms (names must match the algorithm names found in the result files).
- `-g`, `--genes` GENE [GENE ...]: restrict plotting to these genes (also controls order in plots).
- `-o`, `--output` PATH: output directory for CSVs and figures (default `experiments/results`).

What it writes:

- `score.csv` in the chosen output directory.
- Per-instance PNGs: `<instance>/_score_histogram.png`, `<instance>/_score_full.png`, `<instance>/_score_summary.png`.

Example:

```sh
# Evaluate all instances in biological case studies and save plots to the default folder
python src/bntaxonomy/evaluate_score.py -ig experiments/instances/A_case_studies
```

Notes:

- The script expects that you have already run the CLI so that each instance has a corresponding `results` folder (the script maps `instances` -> `results` automatically when you pass `--inst_groups` or `--instances`).
- The generated CSVs are convenient for further analysis or inclusion in figures.
