from __future__ import annotations
import json, os

from colomoto.minibn import BooleanNetwork

from bntaxonomy.utils import CtrlResult


class ExperimentHandler:
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
    ):
        self.name = name
        self.input_path = input_path
        self.output_path = output_path
        self.max_size = max_size
        self.to_console = to_console
        self.to_file = to_file
        self.only_minimal = only_minimal
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
            from bntaxonomy.iface.mpbn import propagate_bn

            self.bn = propagate_bn(self.org_bnet, self.inputs)
        else:
            self.bn = self.org_bnet

        self.pbn_primes = None
        self.sm_attrs, self.sm_primes = None, None
        self.cabean = None

    def postprocess(self, ctrl_result: CtrlResult):
        ctrl_result.drop_size_limit(self.max_size)
        if self.only_minimal:
            ctrl_result.drop_nonminimal()

        if self.to_console:
            print(f"{ctrl_result.name:<14}", ctrl_result)
        if self.to_file:
            ctrl_result.dump(f"{self.output_path}/{ctrl_result.name}.json")

        self.results.append(ctrl_result)
        if os.path.exists("program_instance.asp"):
            os.remove("program_instance.asp")
        return ctrl_result

    ### ActoNet
    def ctrl_actonet_fp(self, **kwargs):
        from bntaxonomy.iface.actonet import ctrl_actonet_fp_iface

        results = ctrl_actonet_fp_iface(self.bn, self.target, self.max_size, **kwargs)
        return self.postprocess(results)

    ### BoNesis
    def ctrl_bonesis_mts(self, **kwargs):
        from bntaxonomy.iface.bonesis import ctrl_bonesis_mts_iface

        results = ctrl_bonesis_mts_iface(self.bn, self.target, self.max_size, **kwargs)
        return self.postprocess(results)

    # TODO: option for fixed point control by BoNesis

    ### Caspo
    def ctrl_caspo_vpts(self, **kwargs):
        from bntaxonomy.iface.caspo import ctrl_caspo_vpts_iface

        return self.postprocess(
            ctrl_caspo_vpts_iface(self.bn, self.target, self.max_size, **kwargs)
        )

    ### PyBoolNet
    def make_pbn_primes(self):
        from bntaxonomy.iface.pbn import make_pbn_primes_iface

        if self.pbn_primes is None:
            self.pbn_primes = make_pbn_primes_iface(self.bnet_fname)

    def ctrl_pyboolnet_model_checking(self, update: str, **kwargs):
        from bntaxonomy.iface.pbn import ctrl_pbn_attr_iface

        assert update in ["synchronous", "asynchronous"]

        self.make_pbn_primes()
        results = ctrl_pbn_attr_iface(
            self.pbn_primes, self.inputs, self.target, update, self.max_size, **kwargs
        )
        return self.postprocess(results)

    def ctrl_pyboolnet_heuristics(self, control_type: str, **kwargs):
        from bntaxonomy.iface.pbn import ctrl_pbn_heuristics_iface

        assert control_type in ["percolation", "trap_spaces"]

        self.make_pbn_primes()
        results = ctrl_pbn_heuristics_iface(
            self.pbn_primes,
            self.inputs,
            self.target,
            control_type,
            self.max_size,
            **kwargs,
        )
        return self.postprocess(results)

    ### pystablemotifs

    def make_sm_primes(self):
        from bntaxonomy.iface.stablemotif import make_sm_primes_iface

        if self.sm_primes is None:
            self.sm_primes = make_sm_primes_iface(self.bnet_fname)

    def make_sm_attrs(self):
        from bntaxonomy.iface.stablemotif import make_sm_attrs_iface

        self.make_sm_primes()
        if self.sm_attrs is None:
            self.sm_attrs = make_sm_attrs_iface(self.sm_primes)

    def ctrl_pystablemotif_brute_force(self, **kwargs):
        from bntaxonomy.iface.stablemotif import ctrl_sm_brute_force_iface

        self.make_sm_primes()

        results = ctrl_sm_brute_force_iface(
            self.sm_primes, self.target, self.max_size, **kwargs
        )
        return self.postprocess(results)

    def ctrl_pystablemotif_trapspace(
        self, target_method: str, driver_method: str, **kwargs
    ):
        from bntaxonomy.iface.stablemotif import ctrl_sm_trapspace_iface

        assert target_method in ["merge", "history"]
        assert driver_method in ["minimal", "internal"]
        self.make_sm_attrs()

        results = ctrl_sm_trapspace_iface(
            self.sm_attrs,
            self.target,
            target_method,
            driver_method,
            self.max_size,
            **kwargs,
        )
        return self.postprocess(results)

    ### CABEAN

    def make_cabean(self):
        from bntaxonomy.iface.cabean import make_cabean_iface

        if self.cabean is None:
            self.cabean = make_cabean_iface(self.bn)

    def ctrl_cabean_phenotype(self, method: str, _debug=False, **kwargs):
        assert method in ["ITC", "TTC", "PTC"]
        self.make_cabean()
        # TODO limit max_size

        from bntaxonomy.iface.cabean import ctrl_target_control_iface

        results = ctrl_target_control_iface(
            self.cabean, self.target, method, _debug, **kwargs
        )
        return self.postprocess(results)
