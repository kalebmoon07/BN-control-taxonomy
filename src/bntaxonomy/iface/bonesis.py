import bonesis
from colomoto.minibn import BooleanNetwork
from bonesis.reprogramming import (
    marker_reprogramming_fixpoints,
    trapspace_reprogramming,
)
from bntaxonomy.utils.control import CtrlResult, suppress_console_output
from bntaxonomy.utils.log import time_check


@time_check
def ctrl_bonesis_mts_iface(
    bn: BooleanNetwork, target: dict[str, int], max_size: int, **kwargs
):
    with suppress_console_output():
        results = list(trapspace_reprogramming(bn, target, max_size))
    return CtrlResult("BoNesis[MTS]", results)


@time_check
def ctrl_bonesis_fp_iface(
    bn: BooleanNetwork, target: dict[str, int], max_size: int, **kwargs
):
    with suppress_console_output():
        results = list(
            marker_reprogramming_fixpoints(bn, target, max_size, at_least_one=False)
        )
    return CtrlResult("BoNesis[FP]", results)
