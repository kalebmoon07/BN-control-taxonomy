import bonesis
from colomoto.minibn import BooleanNetwork

from bntaxonomy.utils import CtrlResult, suppress_console_output


def ctrl_bonesis_mts_iface(
    bn: BooleanNetwork, target: dict[str, int], max_size: int, **kwargs
):
    with suppress_console_output():
        bo = bonesis.BoNesis(bn)
        coP = bo.Some(max_size=max_size)
        with bo.mutant(coP):
            x = bo.cfg()
            bo.in_attractor(x)
            x != bo.obs(target)
        results = list(coP.complementary_assignments())
    return CtrlResult("BoNesis", results)
