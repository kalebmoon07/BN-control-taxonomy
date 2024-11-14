from __future__ import annotations
import functools, json, os, re, sys
import tempfile
from actonet import ActoNet
from algorecell_types import *
from typing import Dict, List
import bonesis
from colomoto.minibn import BooleanNetwork
import caspo_control
import biolqm
import pystablemotifs as sm
import subprocess
import networkx as nx
import graphviz
from itertools import combinations, product
from bntaxonomy.dep.control_strategies_MC import (
    compute_control_strategies_with_model_checking,
)
from bntaxonomy.dep.control_strategies import run_control_problem
from pyboolnet.file_exchange import bnet2primes
import mpbn
import cabean
from cabean.iface import CabeanResult
import contextlib
import io


#  TODO: extend ReprogrammingStrategies and Perturbation to do comparison, minimality check, etc.


@contextlib.contextmanager
def suppress_console_output():
    with open(os.devnull, "w") as devnull:
        # suppress stdout and
        orig_stdout_fno = os.dup(sys.stdout.fileno())
        os.dup2(devnull.fileno(), 1)
        orig_stderr_fno = os.dup(sys.stderr.fileno())
        os.dup2(devnull.fileno(), 2)
        try:
            yield
        finally:
            # restore
            os.dup2(orig_stdout_fno, 1)
            os.dup2(orig_stderr_fno, 2)


def check_smaller(p1: dict[str, int], p2: dict[str, int], strict=False):
    is_small = all(p2.get(k, -1) == v for k, v in p1.items())
    if is_small and strict:
        is_small = is_small and not (p1.keys() == p2.keys())
    return is_small


class CtrlResult:
    def __init__(
        self, name: str, d_list: List[Dict[str, int]], sort=True, only_minimal=True
    ) -> None:
        self.name = name
        self.d_list = d_list
        if only_minimal:
            self.drop_nonminimal()
        if sort:
            self.sort_d_list()

    def __repr__(self) -> str:
        return f"CtrlResult({self.name})"

    def __str__(self) -> str:
        return f"{self.d_list}"

    def iter_ctrl_not_included_by(self, other: CtrlResult) -> bool:
        return (
            x
            for x in self.d_list
            if not any(
                # all(x.get(yk, -1) == yv for yk, yv in y.items()) for y in other.d_list
                check_smaller(y, x)
                for y in other.d_list
            )
        )

    def is_stronger_than(self, other: CtrlResult) -> bool:
        for ctrl in other.iter_ctrl_not_included_by(self):
            return False
        return True

    def dump(self, fname):
        with open(fname, "w") as _f:
            json.dump(self.d_list, _f)

    def sort_d_list(self):
        d_list = [dict(sorted(x.items())) for x in self.d_list]
        d_list.sort(key=lambda x: (len(x), sorted(x.items())))
        self.d_list = d_list

    # def drop_duplicates(self):
    #     self.d_list = list({tuple(d.items()) for d in self.d_list})
    #     self.d_list = [dict(x) for x in self.d_list]

    def drop_nonminimal(self):
        d_list = list()
        for ctrl in self.d_list:
            if not any(True for other in d_list if check_smaller(other, ctrl)):
                d_list.append(ctrl)
        self.d_list = d_list


def refine_pert(s: ReprogrammingStrategies):
    text = str(s.perturbations())
    d_list = []
    pattern = re.compile(r"PermanentPerturbation\(([^)]*)\)")
    for match in pattern.findall(text):
        d = dict()
        for item in match.split(", "):
            if "=" not in item:
                continue
            k, v = item.split("=")
            d[k] = int(v)
        d_list.append(d)
    return d_list


class myActoNet(ActoNet):
    def __init__(self, bn, inputs=dict()):
        super().__init__(bn, inputs)
        self.controls = functools.partial(self.controls, existential=False)


class ExperimentRun:
    def __init__(self, results: list[CtrlResult], name: str = "Hierarchy"):
        self.name = name
        self.results = results
        self.G = nx.DiGraph()
        for r1, r2 in combinations(self.results, 2):
            if r1.is_stronger_than(r2):
                self.G.add_edge(r2.name, r1.name)
            if r2.is_stronger_than(r1):
                self.G.add_edge(r1.name, r2.name)

    def save(self, fname: str):
        nx.nx_pydot.write_dot(self.G, f"{fname}.dot")
        tred_cmd = f"tred {fname}.dot | dot -T png > {fname}.png"
        process = subprocess.Popen(tred_cmd, shell=True)
        process.wait()

    @staticmethod
    def from_folder(opath: str, name: str = ""):
        files = [fname for fname in os.listdir(opath) if fname.endswith(".json")]
        sol_list = [
            CtrlResult(fname[:-5], json.load(open(f"{opath}/{fname}")))
            for fname in files
        ]
        if not name:
            name = opath.split("/")[-1]
        return ExperimentRun(sol_list, name)


class ExperimentHandler:
    def __init__(
        self,
        name: str,
        input_path: str,
        output_path: str,
        use_propagated: bool = True,
        precompute_pbn: bool = True,
        precompute_sm: bool = True,
        precompute_cabean: bool = True,
        to_console: bool = True,
        to_file: bool = True,
        only_minimal: bool = True,
        store_results: bool = True,
    ):
        self.name = name
        self.input_path = input_path
        self.output_path = output_path
        self.to_console = to_console
        self.to_file = to_file
        self.only_minimal = only_minimal
        self.store_results = store_results
        self.results = list()
        os.makedirs(output_path, exist_ok=True)

        # Load setting
        with open(f"{input_path}/setting.json") as _f:
            setting = json.load(_f)
        self.inputs: dict[str, int] = setting["inputs"]
        self.target: dict[str, int] = setting["target"]

        # Load the original Boolean network
        self.bnet_fname = f"{self.input_path}/transition_formula.bnet"
        self.org_bnet = BooleanNetwork(data=self.bnet_fname)

        # Load the Boolean network for experiments
        if use_propagated:
            self.bnet_fname = f"{self.input_path}/propagated.bnet"
            if not os.path.exists(self.bnet_fname):
                propagated_mbn = mpbn.MPBooleanNetwork(self.org_bnet)
                for k, v in self.inputs.items():
                    propagated_mbn[k] = v
                propagated_mbn.propagate_constants()
                for k, v in propagated_mbn.constants().items():
                    propagated_mbn.pop(k)
                propagated_mbn.save(self.bnet_fname)
            self.bn = BooleanNetwork(data=self.bnet_fname)
        else:
            self.bn = self.org_bnet

        self.pbn_primes = None
        if precompute_pbn:
            self.make_pbn_primes()

        self.sm_attrs, self.sm_primes = None, None
        if precompute_sm:
            self.make_sm_primes()

        self.cabean = None
        if precompute_cabean:
            self.make_cabean()

    def make_pbn_primes(self):
        if self.pbn_primes is None:
            self.pbn_primes = bnet2primes(self.bnet_fname)

    def make_sm_primes(self):
        if self.sm_primes is None:
            self.sm_primes = sm.format.import_primes(self.bnet_fname)
            self.sm_attrs = sm.AttractorRepertoire.from_primes(self.sm_primes)

    def make_cabean(self):
        if self.cabean is None:
            self.cabean = cabean.load(self.bn)
            self.cabean_bn_fname = f"{self.input_path}/bn.ispl"
            self.cabean_target_fname = f"{self.input_path}/phenotype.txt"
            # fd, tmpfile = tempfile.mkstemp(suffix=".ispl", prefix="cabean")
            # os.close(fd)
            with open(self.cabean_target_fname, "w") as _f:
                _f.write("node, value\n")
                for k, v in self.target.items():
                    _f.write(f"{k},{v}\n")
            with open(self.cabean_bn_fname, "w") as _f:
                self.cabean.iface.write_ispl(_f)

    def postprocess(self, ctrl_result: CtrlResult):
        if self.to_console:
            print(f"{ctrl_result.name:<14}", ctrl_result)
        if self.to_file:
            ctrl_result.dump(f"{self.output_path}/{ctrl_result.name}.json")
        if self.store_results:
            self.results.append(ctrl_result)
        if os.path.exists("program_instance.asp"):
            os.remove("program_instance.asp")
        return ctrl_result

    def ctrl_ActoNet(self, max_size: int, **kwargs):
        with suppress_console_output():
            model = myActoNet(self.bn)
            s = model.reprogramming_fixpoints(self.target, maxsize=max_size, **kwargs)

        return self.postprocess(
            CtrlResult("ActoNet", refine_pert(s), self.only_minimal)
        )

    def ctrl_BoNesis(self, max_size: int, **kwargs):
        # TODO: option for fixed point control
        with suppress_console_output():
            bo = bonesis.BoNesis(self.bn)
            coP = bo.Some(max_size=max_size)
            with bo.mutant(coP):
                x = bo.cfg()
                bo.in_attractor(x)
                x != bo.obs(self.target)
            results = list(coP.complementary_assignments())
        return self.postprocess(CtrlResult("BoNesis", results, self.only_minimal))

    def ctrl_Caspo(self, max_size: int, **kwargs):
        with suppress_console_output():
            model = caspo_control.CaspoControl(self.bn, {})
            s = model.reprogramming_to_attractor(
                self.target, maxsize=max_size, **kwargs
            )
        return self.postprocess(CtrlResult("Caspo", refine_pert(s), self.only_minimal))

    def ctrl_pyboolnet_mc(self, max_size: int, update: str, **kwargs):
        assert update in ["synchronous", "asynchronous"]
        with suppress_console_output():
            result = compute_control_strategies_with_model_checking(
                primes=self.pbn_primes,
                avoid_nodes=self.inputs,
                limit=max_size,
                target=[self.target],
                update=update,
                **kwargs,
            )
        return self.postprocess(
            CtrlResult(f"PBN-mc-{update[:-7]}", result, self.only_minimal)
        )

    def ctrl_pyboolnet_heuristics(self, max_size: int, control_type: str, **kwargs):
        assert control_type in ["percolation", "trap_spaces"]
        with suppress_console_output():
            results = run_control_problem(
                primes=self.pbn_primes,
                avoid_nodes=self.inputs,
                limit=max_size,
                target=self.target,
                control_type=control_type,
                intervention_type="node",
                **kwargs,
            )
        return self.postprocess(
            CtrlResult(f"PBN-{control_type[:4]}", results, self.only_minimal)
        )

    def ctrl_pystablemotif_brute_force(self, max_size: int, **kwargs):
        result = sm.drivers.knock_to_partial_state(
            self.target,
            self.sm_primes,
            min_drivers=0,
            max_drivers=max_size,
            **kwargs,
        )
        return self.postprocess(CtrlResult("SM-bf", result, self.only_minimal))

    def ctrl_pystablemotif_trapspace(
        self, max_size: int, target_method: str, driver_method: str, **kwargs
    ):
        assert target_method in ["merge", "history"]
        assert driver_method in ["minimal", "internal"]

        result = self.sm_attrs.reprogram_to_trap_spaces(
            self.target,
            target_method=target_method,
            driver_method=driver_method,
            max_drivers=max_size,
            **kwargs,
        )
        return self.postprocess(
            CtrlResult(
                f"SM-{target_method}-{driver_method[:3]}", result, self.only_minimal
            )
        )

    def ctrl_cabean_phenotype(self, max_size: int, method: str, **kwargs):
        assert method in ["ITC", "TTC", "PTC"]
        # TODO limit max_size
        cmd = [
            "./cabean_2.0.0",
            "-compositional",
            "2",
            "-control",
            method,
            "-tmarker",
            self.cabean_target_fname,
            self.cabean_bn_fname,
        ]
        # cmd =  f"sh ./cabean_2.0.0 -compositional 2 -control {method} -tmarker {self.cabean_target_fname} {self.cabean_bn_fname}"
        output = subprocess.run(cmd, capture_output=True)
        # return output
        result = CabeanResult(self.cabean.iface, output.stdout.decode())
        ctrl_list = []
        line_iter = iter(result.lines)
        line = next(line_iter, "")
        while not "DECOMP" in line:
            line = next(line_iter, "")

        for line in line_iter:
            print(line)
            line: str
            if line.startswith("There is only one attractor."):
                ctrl_list.append(dict())
            if line.lower().startswith("control set"):
                p = dict()
                for c in line.strip().split():
                    if "=" in c:
                        node, value = c.split("=")
                        p[node] = int(value)
                ctrl_list.append(p)
            #     ctrl_list.append(result.parse_controlset(line))
        ctrl_result = CtrlResult(f"CABEAN-{method}", ctrl_list)
        return self.postprocess(ctrl_result)

        # TODO: parse the output

    def get_run(self, name: str):
        return ExperimentRun(self.results, name)
