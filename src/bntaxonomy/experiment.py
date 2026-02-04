from __future__ import annotations
import json
import os

from colomoto.minibn import BooleanNetwork

from bntaxonomy.utils.control import CtrlResult
from bntaxonomy.utils.log import main_logger
from bntaxonomy.utils.control import suppress_console_output

from bntaxonomy.iface import registered_tools


class ExperimentHandler:
    __next_id = 1

    def __init__(
        self,
        name: str,
        input_path: str,
        output_path: str,
        max_size: int,
        use_propagated: bool = True,
        to_console: bool = True,
        to_file: bool = True,
        only_minimal: bool = True,
        exclude_targets: bool = False,
        dump_full: bool = True,
        load_precompute: bool = False,
    ):
        self.name = name
        self.input_path = input_path
        self.output_path = output_path
        self.max_size = max_size
        self.to_console = to_console
        self.to_file = to_file
        self.only_minimal = only_minimal
        self.dump_full = dump_full
        self.load_precompute = load_precompute
        self.results: list[CtrlResult] = list()
        os.makedirs(output_path, exist_ok=True)

        # Load setting
        with open(f"{input_path}/setting.json") as _f:
            setting = json.load(_f)
        self.inputs: dict[str, int] = setting["inputs"]
        self.target: dict[str, int] = setting["target"]
        self.exclude: list[str] = setting.get("exclude", [])
        if exclude_targets:
            self.exclude.extend(self.target)

        # Load the original Boolean network
        self.bnet_fname = f"{self.input_path}/transition_formula.bnet"
        self.org_bnet = BooleanNetwork(data=self.bnet_fname)

        # Load the Boolean network for experiments
        if use_propagated:
            from bntaxonomy.iface.mpbn import propagate_bn

            self.bn = propagate_bn(self.org_bnet, self.inputs)
            self.inputs = {}
        else:
            self.bn = self.org_bnet
            self.bn |= self.inputs
            self.inputs = {}

        self.cachedir = os.path.join(self.input_path, "cache")
        if not os.path.isdir(self.cachedir):
            os.makedirs(self.cachedir)

        self.bnet_file = os.path.join(self.cachedir, "model.bnet")
        self.bn.save(self.bnet_file)

        self.expid = f"Experiment_{self.__next_id}_{id(self)}"
        self.__class__.__next_id += 1

        self.sm_attrs = None
        self.primes = None
        self.cabean = None

    def postprocess(self, ctrl_result: CtrlResult):
        ctrl_result.sort_d_list()
        if self.dump_full:
            if self.to_console:
                main_logger.info(f"{ctrl_result.name:<14}: {ctrl_result}")
            ctrl_result.dump(f"{self.output_path}/{ctrl_result.name}_full.json")

        # filtering
        ctrl_result.drop_size_limit(self.max_size)
        if self.only_minimal:
            ctrl_result.drop_nonminimal()

        if self.to_console:
            main_logger.info(f"{ctrl_result.name:<14}: {ctrl_result}")
        if self.to_file:
            ctrl_result.dump(f"{self.output_path}/{ctrl_result.name}.json")

        self.results.append(ctrl_result)
        if os.path.exists("program_instance.asp"):
            os.remove("program_instance.asp")
        return ctrl_result

    def run_tools(self, filter_tools):
        for toolcls in registered_tools():
            if filter_tools and toolcls.name not in filter_tools:
                continue

            main_logger.info(f"Running {toolcls.name}")
            args = (self.expid, self.cachedir) if toolcls.uses_cache else ()

            if toolcls.bn_type == "bnet_file":
                bninp = self.bnet_file
            elif toolcls.bn_type == "colomoto.BooleanNetwork":
                bninp = self.bn
            else:
                raise TypeError(
                    f"{toolcls.name}: Unknown BN type input {toolcls.bn_type}"
                )
            with suppress_console_output():
                res = toolcls.run(
                    bninp, self.max_size, self.target, self.exclude, *args
                )
            res = CtrlResult(toolcls.name, res)
            self.postprocess(res)

        for toolcls in registered_tools():
            if toolcls.uses_cache:
                main_logger.info(f"Cleaning cache for {toolcls.name}")
                toolcls.free_experiment(self.expid)
