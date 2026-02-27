import functools
from bntaxonomy.iface import register_tool
from bntaxonomy.utils.control import refine_pert
from bntaxonomy.utils.log import time_check

from actonet import ActoNet
class myActoNet(ActoNet):
    def __init__(self, bn, inputs=dict()):
        super().__init__(bn, inputs)
        # a fixpoint need not exist; if a control induces no fixpoint, it is a valid control
        self.controls = functools.partial(self.controls, existential=False)


@register_tool
class ActoNetFP:
    name = "ActoNet"
    bn_type = "colomoto.BooleanNetwork"

    @time_check
    @staticmethod
    def run(bn, max_size, target, exclude, inputs={}):
        a = myActoNet(bn, inputs)
        r = a.reprogramming_fixpoints(target, maxsize=max_size)
        return refine_pert(r)
