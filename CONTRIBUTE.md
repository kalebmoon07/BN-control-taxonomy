# Contributing to the project

## Adding a new input instance

To make a new instance group, create a new subfolder under `experiments/instances/`, and add one or more instance subfolders inside it. Each instance folder must contain:

- `transition_formula.bnet`: A Boolean network model file.
- `setting.json`: A JSON file specifying the control settings. Example format:
  - `inputs`: A dictionary specifying fixed input components and their values (value propagation applied accordingly).
  - `target`: A dictionary specifying the desired target phenotype as a conjunction of components and their values.

```json
// Example settings.json
{
    "inputs": {
        "GrowthInhibitors": 1,
        "EGFR_stimulus": 1,
        "FGFR3_stimulus": 1
    },
    // Apoptosis_b1 AND (NOT RB1)
    "target": {
        "Apoptosis_b1": 1,
        "RB1": 0
    }
}
```

You use an [automated input generator](https://github.com/pauleve/BN-example-generator) for creating inputs with some constraints ([Software Heritage](https://archive.softwareheritage.org/swh:1:rev:2ceb1917c9311e624af5632269fca1d0fc04ef31)).

## Adding a new software tool

To add a new software tool, make an interface to the tool in `src/bntaxonomy/iface`, following the structure of [a template file](src/bntaxonomy/iface/mytool.py.tmpl). Following steps are needed for [cli.py](src/bntaxonomy/cli.py) to grep the tool and automatically run it via `ExperimentHandler.run_tools`.

- Make a class with `@register_tool` decorator for each tool.
  - If more than one configuration is available, create a class for each option.
  - Specify a unique name (e.g., `BoNesis[FP]`), which will be used in the experiment configuration and result files.
  - If a cache is needed, set `uses_cache` to `True` (please refer to `stablemotif.py` for cache usage).
- Define static methods
  - `run`: to execute the tool with given parameters (see the parameters in the template file).
  - `free_experiment`: to clear cached data for a given experiment id.
- Put external files under `src/bntaxonomy/dep` and use relative paths in the `run` method to access them.

A suggested way to interface with the tool is as follows:

- Upload the whole package to anaconda.org or PyPI, so that it can be installed via `conda` or `pip`.
  - Specify the version of your package in `requirements_tools.txt`
- Directly include the tool's source code in this repository.
- (Optional) Integrate it into the [CoLoMoTo Docker](https://github.com/colomoto/colomoto-docker/blob/for-next/CONTRIBUTING.md), and use a Docker image to run experiments.
  - You may refer to `src/bntaxonomy/iface/pyboolnet.py` for an example.
