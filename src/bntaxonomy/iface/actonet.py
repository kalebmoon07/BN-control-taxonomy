import functools
from actonet import ActoNet
from colomoto.minibn import BooleanNetwork

from bntaxonomy.utils.control import CtrlResult, suppress_console_output, refine_pert


class myActoNet(ActoNet):
    def __init__(self, bn: BooleanNetwork, inputs=dict()):
        super().__init__(bn, inputs)
        # a fixpoint need not exist; if a control induces no fixpoint, it is a valid control
        self.controls = functools.partial(self.controls, existential=False)


def ctrl_actonet_fp_iface(
    bn: BooleanNetwork, target: dict[str, int], max_size: int, **kwargs
):
    with suppress_console_output():
        model = myActoNet(bn)
        s = model.reprogramming_fixpoints(target, maxsize=max_size, **kwargs)
    return CtrlResult("ActoNet", refine_pert(s))
