import os
import tempfile
from colomoto.minibn import BooleanNetwork
import mpbn


def propagate_bn(org_bnet: BooleanNetwork, inputs: dict[str, int]):
    print(f"Propagating with inputs {inputs}")
    fd, bnet_fname = tempfile.mkstemp(suffix=".bnet", prefix="propagated")
    os.close(fd)
    propagated_mbn = mpbn.MPBooleanNetwork(org_bnet)
    propagated_mbn.simplify(in_place=True)
    for k, v in inputs.items():
        propagated_mbn[k] = v
    propagated_mbn.propagate_constants()
    for k, v in propagated_mbn.constants().items():
        propagated_mbn.pop(k)
    propagated_mbn.save(bnet_fname)
    bn = BooleanNetwork(data=bnet_fname)
    return bn, bnet_fname
