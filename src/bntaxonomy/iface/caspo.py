from bntaxonomy.iface import register_tool
from bntaxonomy.utils.control import refine_pert
from bntaxonomy.utils.log import time_check

from colomoto.minibn import BooleanNetwork
import caspo_control

@register_tool
class CaspoVPTS:
    name = "Caspo"
    bn_type = "colomoto.BooleanNetwork"

    @time_check
    @staticmethod
    def run(bn: BooleanNetwork, max_size: int,
            target: dict[str, int], exclude: list[str]):
        model = caspo_control.CaspoControl(bn, {})
        s = model.reprogramming_to_attractor(target, maxsize=max_size)
        return refine_pert(s)
