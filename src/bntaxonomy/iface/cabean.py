import os
import subprocess
import tempfile
from bntaxonomy.utils.control import CtrlResult
from colomoto.minibn import BooleanNetwork
import cabean
from cabean import CabeanInstance
from cabean.iface import CabeanResult

cabean_path = f"{os.path.dirname(os.path.abspath(__file__))}/../dep/cabean_2.0.0"

CABEAN_OUT_MEMORY = False


def make_cabean_iface(bn: BooleanNetwork):
    print("Loading cabean and computing attractors")
    try:
        return cabean.load(bn)
    except Exception as e:
        print(f"Error loading cabean: {e}")
        return CABEAN_OUT_MEMORY


def make_cabean_tempfiles(cabean_obj: CabeanInstance, target: dict[str, int]):
    fd, cabean_bn_fname = tempfile.mkstemp(suffix=".ispl", prefix="cabean_bn")
    os.close(fd)
    fd, cabean_target_fname = tempfile.mkstemp(suffix=".txt", prefix="cabean_target")
    os.close(fd)
    with open(cabean_target_fname, "w") as _f:
        _f.write("node, value\n")
        for k, v in target.items():
            _f.write(f"{k},{v}\n")
    with open(cabean_bn_fname, "w") as _f:
        cabean_obj.iface.write_ispl(_f)
    return cabean_target_fname, cabean_bn_fname


def ctrl_target_control_iface(
    cabean_obj: CabeanInstance,
    target: dict[str, int],
    method: str,
    _debug=False,
    **kwargs,
):

    cabean_target_fname, cabean_bn_fname = make_cabean_tempfiles(cabean_obj, target)
    cmd = [
        cabean_path,
        "-compositional",
        "2",
        "-control",
        method,
        "-tmarker",
        cabean_target_fname,
        cabean_bn_fname,
    ]
    output = subprocess.run(cmd, capture_output=True)
    result = CabeanResult(cabean_obj.iface, output.stdout.decode())
    ctrl_list = []
    line_iter = iter(result.lines)
    line = next(line_iter, "")
    while not "DECOMP" in line:
        line = next(line_iter, "")

    for line in line_iter:
        line: str

        if line.startswith("There is only one attractor."):
            # Trivial solution if the target already satisfies a phenotype
            # Otherwise, the phenotype is not satisfied
            is_phenotype = True
            for state in result.attractors.values():
                if isinstance(state, list):
                    states = state
                elif isinstance(state, dict):
                    states = [state]
                else:
                    raise ValueError(f"Unknown type {type(state)}")
                for state in states:
                    for k, v in target.items():
                        if state[k] != v:
                            is_phenotype = False
            if is_phenotype:
                ctrl_list.append(dict())
            if _debug:
                print(
                    line,
                    cabean_obj.attractors,
                    f"phenotype: {target}, {is_phenotype}",
                )
        if line.startswith("Error:"):
            # Error: could not find any attractor based on the markers of attractors.
            if _debug:
                print(line)
        if line.lower().startswith("control set"):
            p = dict()
            for c in line.strip().split():
                if "=" in c:
                    node, value = c.split("=")
                    p[node] = int(value)
            ctrl_list.append(p)
    return CtrlResult(f"CABEAN-{method}", ctrl_list)
