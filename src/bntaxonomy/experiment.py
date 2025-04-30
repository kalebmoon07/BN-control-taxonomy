from __future__ import annotations
import json, os

from colomoto.minibn import BooleanNetwork

from bntaxonomy.utils.control import CtrlResult
from bntaxonomy.utils.log import main_logger
import bntaxonomy.utils.log as log_utils


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

        # Load the original Boolean network
        self.bnet_fname = f"{self.input_path}/transition_formula.bnet"
        self.org_bnet = BooleanNetwork(data=self.bnet_fname)

        # Load the Boolean network for experiments
        if use_propagated:
            from bntaxonomy.iface.mpbn import propagate_bn

            self.bn, self.bnet_fname = propagate_bn(self.org_bnet, self.inputs)
        else:
            self.bn, self.bnet_fname = self.org_bnet, self.bnet_fname

        self.primes = None
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

    def ctrl_bonesis_fp(self, **kwargs):
        from bntaxonomy.iface.bonesis import ctrl_bonesis_fp_iface

        results = ctrl_bonesis_fp_iface(self.bn, self.target, self.max_size, **kwargs)
        return self.postprocess(results)

    ### optboolnet

    def ctrl_optboolnet_sync_attr(self, **kwargs):
        from bntaxonomy.iface.optboolnet import (
            make_CNFBooleanNetwork_iface,
            ctrl_optboolnet_sync_attr_iface,
        )

        cnf_bn = make_CNFBooleanNetwork_iface(self.bn, self.inputs, self.target)
        results = ctrl_optboolnet_sync_attr_iface(self.name, cnf_bn, self.max_size)

        return self.postprocess(results)

    def ctrl_optboolnet_fp(self, **kwargs):
        from bntaxonomy.iface.optboolnet import (
            make_CNFBooleanNetwork_iface,
            ctrl_optboolnet_fp_iface,
        )

        cnf_bn = make_CNFBooleanNetwork_iface(self.bn, self.inputs, self.target)
        results = ctrl_optboolnet_fp_iface(self.name, cnf_bn, self.max_size)

        return self.postprocess(results)

    ### Caspo
    def ctrl_caspo_vpts(self, **kwargs):
        from bntaxonomy.iface.caspo import ctrl_caspo_vpts_iface

        return self.postprocess(
            ctrl_caspo_vpts_iface(self.bn, self.target, self.max_size, **kwargs)
        )

    ### PyBoolNet
    def make_primes(self):
        from bntaxonomy.iface.pbn import make_pbn_primes_iface, PRIME_JSON_FILE

        if self.primes is None:
            if self.load_precompute:
                try:
                    with open(f"{self.input_path}/{PRIME_JSON_FILE}") as _f:
                        self.primes = json.load(_f)
                    main_logger.info("Loaded precomputed pyboolnet primes successfully")
                except:
                    main_logger.info("Loading precomputed pyboolnet primes fails")
                    self.primes = make_pbn_primes_iface(self.bnet_fname)
            else:
                self.primes = make_pbn_primes_iface(self.bnet_fname)
            with open(f"{self.input_path}/{PRIME_JSON_FILE}", "w") as _f:
                json.dump(self.primes, _f)

    def ctrl_pyboolnet_model_checking(self, update: str, **kwargs):
        from bntaxonomy.iface.pbn import ctrl_pbn_attr_iface

        assert update in ["synchronous", "asynchronous"]

        self.make_primes()
        results = ctrl_pbn_attr_iface(
            self.primes, self.inputs, self.target, update, self.max_size, **kwargs
        )
        return self.postprocess(results)

    def ctrl_pyboolnet_heuristics(self, control_type: str, **kwargs):
        from bntaxonomy.iface.pbn import ctrl_pbn_heuristics_iface

        assert control_type in ["percolation", "trap_spaces"]

        self.make_primes()
        results = ctrl_pbn_heuristics_iface(
            self.primes,
            self.inputs,
            self.target,
            control_type,
            self.max_size,
            **kwargs,
        )
        return self.postprocess(results)

    ### pystablemotifs

    def make_sm_attrs(self):
        from bntaxonomy.iface.stablemotif import make_sm_attrs_iface

        self.make_primes()
        if self.sm_attrs is None:
            self.sm_attrs = make_sm_attrs_iface(self.primes)

    def ctrl_pystablemotif_brute_force(self, **kwargs):
        from bntaxonomy.iface.stablemotif import ctrl_sm_brute_force_iface

        self.make_primes()

        results = ctrl_sm_brute_force_iface(
            self.primes, self.target, self.max_size, **kwargs
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
        from bntaxonomy.iface.cabean import (
            make_cabean_iface,
            CABEAN_OUT_MEMORY,
            ATTR_JSON_FILE,
        )

        # if experienced out of memory, skip retrying
        if self.cabean == CABEAN_OUT_MEMORY:
            return
        if self.cabean is None:
            try:
                self.cabean = make_cabean_iface(self.bn)
                if self.load_precompute:
                    try:
                        with open(f"{self.input_path}/{ATTR_JSON_FILE}") as _f:
                            self.cabean.load_precomputed_attr(json.load(_f))
                        main_logger.info("Loaded precomputed cabean successfully")
                    except:
                        main_logger.info("Loading precomputed cabean fails")
                        self.cabean.compute_attractors()
                else:
                    self.cabean.compute_attractors()
                with open(f"{self.input_path}/{ATTR_JSON_FILE}", "w") as _f:
                    json.dump(self.cabean.attractors, _f)
            except Exception as e:
                main_logger.info(f"Error loading cabean: {e}")
                self.cabean = CABEAN_OUT_MEMORY
                return

    def ctrl_cabean_target_control(self, method: str, _debug=False, **kwargs):
        assert method in ["ITC", "TTC", "PTC"]
        self.make_cabean()

        from bntaxonomy.iface.cabean import ctrl_target_control_iface, CABEAN_OUT_MEMORY

        if self.cabean == CABEAN_OUT_MEMORY:
            # if experienced out of memory, return empty result
            results = CtrlResult(f"CABEAN[{method}]", [])
        else:
            results = ctrl_target_control_iface(
                self.cabean, self.target, method, _debug, **kwargs
            )
        return self.postprocess(results)
