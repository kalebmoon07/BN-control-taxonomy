from bntaxonomy.iface import register_tool

from colomoto.minibn import BooleanNetwork
import caspo_control

from bntaxonomy.utils.control import refine_pert
from bntaxonomy.utils.log import time_check

@register_tool
class CaspoVPTS:
    bn_type = "colomoto.BooleanNetwork"

    @time_check
    @staticmethod
    def run(bn: BooleanNetwork, max_size: int,
            target: dict[str, int], exclude: list[str]):
        model = caspo_control.CaspoControl(bn, {})
        s = model.reprogramming_to_attractor(target, maxsize=max_size)
        return refine_pert(s)
