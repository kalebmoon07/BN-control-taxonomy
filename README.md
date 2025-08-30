## Installation

### Using CoLoMoTo Docker distribution

This is the recommended usage.

You need [`colomoto-docker`](https://github.com/colomoto/colomoto-docker-py) command, see also https://colomoto.github.io/colomoto-docker.

### General usage

You need a recent Python. Required dependences are listed in `requirements.txt`, and can be installed using `pip`:

```
pip install -r requirements.txt
```

**Important** Each tool must be properly installed, or it will be excluded from the analysis.

## Usage

General usage:

```sh
python src/bntaxonomy/cli.py {max_size} {instances...} [--tools {tools...}]
```

**Recommended** Usage within CoLoMoTo Docker distribution:

```sh
colomoto-docker --bind . python src/bntaxonomy/cli.py {max_size} {instances...} [--tools {tools...}]
```


*Examples*

```sh
# run all available tools on counterexamples instances with max_size=2
python src/bntaxonomy/cli.py 2 experiments/instances/counterexamples/* --tools 'BoNesis[FP]' 'BoNesis[MTS]'

# same, but run only BoNesis[FP] and BoNesis[MTS] tools
python src/bntaxonomy/cli.py 2 experiments/instances/counterexamples/* --tools 'BoNesis[FP]' 'BoNesis[MTS]'
```


