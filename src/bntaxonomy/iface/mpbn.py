import os
import tempfile
from colomoto.minibn import BooleanNetwork
import mpbn


def propagate_bn(org_bnet: BooleanNetwork, inputs: dict[str, int]):
    print(f"Propagating with inputs {inputs}")
    f = mpbn.MPBooleanNetwork(org_bnet)
    f.simplify(in_place=True)
    for k, v in inputs.items():
        f[k] = v
    f.propagate_constants()
    for k in f.constants():
        f.pop(k)
    return f
