from colomoto.minibn import BooleanNetwork
import caspo_control

from bntaxonomy.utils.control import CtrlResult, refine_pert, suppress_console_output


def ctrl_caspo_vpts_iface(
    bn: BooleanNetwork, target: dict[str, int], max_size: int, **kwargs
):
    with suppress_console_output():
        model = caspo_control.CaspoControl(bn, {})
        s = model.reprogramming_to_attractor(target, maxsize=max_size, **kwargs)
    return CtrlResult("Caspo", refine_pert(s))
