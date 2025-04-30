"""This module integrates the `optboolnet` library into the project.

`optboolnet` is a The optimization toolbox for Boolean network analysis.
For more details, visit the official repository or PyPI page:

- GitHub: https://github.com/MSOLab/optboolnet
- PyPI: https://pypi.org/project/optboolnet/

To install the library, use the following command:
`pip install optboolnet`
"""

from colomoto.minibn import BooleanNetwork
from optboolnet.algorithm import BendersAttractorControl, BendersFixPointControl
from optboolnet.boolnet import CNFBooleanNetwork
from optboolnet.config import ControlConfig, BendersConfig, LoggingConfig, SolverConfig
from bntaxonomy.utils.control import (
    CtrlResult,
    convert_minibn_perturbation_to_dict,
    refine_pert,
)
from bntaxonomy.utils.log import time_check


@time_check
def make_CNFBooleanNetwork_iface(bn: BooleanNetwork, inputs: dict, target: dict):
    new_bn = BooleanNetwork(bn)
    config = ControlConfig()
    # config.fixed_values = inputs
    config.fixed_values = dict()
    assert len(target) > 0
    config.controllable_vars = list(new_bn.keys())

    if (len(target) == 1) and (target[str(list(target.keys())[0])] == 1):
        config.phenotype = str(list(target.keys())[0])
        config.uncontrollable_vars = list()
    else:
        config.phenotype = "PHENOTYPE"
        cnf_clauses = []
        for var, value in target.items():
            if value == 1:
                cnf_clauses.append(f"{var}")
            else:
                cnf_clauses.append(f"~{var}")
        cnf_formula = " & ".join(cnf_clauses)
        new_bn[config.phenotype] = cnf_formula
        config.uncontrollable_vars = [config.phenotype]
    return CNFBooleanNetwork(new_bn, config, to_cnf=True)


@time_check
def ctrl_optboolnet_sync_attr_iface(
    name: str, cnf_bn: CNFBooleanNetwork, max_size: int, max_length: int = 45
):
    s = BendersAttractorControl(name, cnf_bn).get_control_strategies(
        max_control_size=max_size, max_length=max_length, solve_separation=True
    )
    results = convert_minibn_perturbation_to_dict(s.perturbations())
    return CtrlResult(f"optbn[SA]", results)


@time_check
def ctrl_optboolnet_fp_iface(name, cnf_bn: CNFBooleanNetwork, max_size: int):
    s = BendersFixPointControl(name, cnf_bn).get_control_strategies(
        max_control_size=max_size, max_length=1, solve_separation=False
    )
    results = convert_minibn_perturbation_to_dict(s.perturbations())
    return CtrlResult("optbn[FP]", results)
