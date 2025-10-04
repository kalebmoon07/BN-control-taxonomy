from bntaxonomy.iface import register_tool
from bntaxonomy.utils.control import CtrlResult, refine_pert
from bntaxonomy.utils.log import time_check

from actonet import ActoNet


@register_tool
class ActoNet:
    bn_type = "colomoto.BooleanNetwork"

    @time_check
    @staticmethod
    def run(bn, max_size, target, exclude, inputs={}):
        a = ActoNet(bn, inputs)
        r = a.reprogramming_fixpoints(target, maxsize=max_size, ignore=exclude)
        return refine_pert(r)
