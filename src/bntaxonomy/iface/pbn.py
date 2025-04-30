from pyboolnet.file_exchange import bnet2primes

from bntaxonomy.dep.control_strategies import run_control_problem
from bntaxonomy.dep.control_strategies_MC import (
    compute_control_strategies_with_model_checking,
)
from bntaxonomy.utils.control import CtrlResult, suppress_console_output
from bntaxonomy.utils.log import time_check

PRIME_JSON_FILE = "pyboolnet_primes.json"


@time_check
def make_pbn_primes_iface(bnet_fname: str):
    print("Generating primes from pyboolnet")
    return bnet2primes(bnet_fname)


@time_check
def ctrl_pbn_attr_iface(
    pbn_primes: dict,
    inputs: dict[str, int],
    target: dict[str, int],
    update: str,
    max_size: int,
    **kwargs,
):
    with suppress_console_output():
        update_flage = "SA" if update == "synchronous" else "ASA"
        results = compute_control_strategies_with_model_checking(
            primes=pbn_primes,
            avoid_nodes=inputs,
            limit=max_size,
            target=[target],
            update=update,
            **kwargs,
        )
    return CtrlResult(f"PBN[{update_flage}]", results)


@time_check
def ctrl_pbn_heuristics_iface(
    pbn_primes: dict,
    inputs: dict[str, int],
    target: dict[str, int],
    control_type: str,
    max_size: int,
    **kwargs,
):
    with suppress_console_output():
        results = run_control_problem(
            primes=pbn_primes,
            avoid_nodes=inputs,
            limit=max_size,
            target=target,
            control_type=control_type,
            intervention_type="node",
            **kwargs,
        )
        return CtrlResult(f"PBN[{control_type}]", results)
