import os, sys
from bntaxonomy.experiment import ExperimentHandler
from bntaxonomy.hierarchy import SingleInputSummary
import bntaxonomy.utils.log as log_utils

cwd = os.path.dirname(os.path.abspath(__file__))

insts_practical = f"{cwd}/instances/practical"
insts_tonello = f"{cwd}/instances/Tonello"
insts_ce = f"{cwd}/instances/counterexamples"


def get_insts(inst_groups: list[str]):
    for inst_selected in inst_groups:
        inst_list = sorted(os.listdir(inst_selected))
        for inst in inst_list:
            yield inst_selected, inst


if __name__ == "__main__":
    log_utils.configure_logging(os.path.basename(__file__))

    max_size = 2
    inst_groups = [insts_ce, insts_tonello, insts_practical]
    # inst_groups = [insts_ce]

    for inst_selected, inst in get_insts(inst_groups):
        log_utils.main_logger.info(f"\tRunning {inst}")
        fpath = f"{inst_selected}/{inst}"
        opath = f"{inst_selected.replace('instances','results')}/{inst}"
        os.makedirs(opath, exist_ok=True)

        exp = ExperimentHandler(
            inst,
            fpath,
            opath,
            max_size,
            to_console=True,
            to_file=True,
            only_minimal=True,
            load_precompute=True,
        )

        exp.ctrl_actonet_fp()
        exp.ctrl_bonesis_mts()
        exp.ctrl_caspo_vpts()
        for update in ["synchronous", "asynchronous"]:
            exp.ctrl_pyboolnet_model_checking(update)
        for control_type in ["percolation", "trap_spaces"]:
            exp.ctrl_pyboolnet_heuristics(control_type)
        exp.ctrl_pystablemotif_brute_force()
        for target_method in ["merge"]:  # ["merge", "history"]:
            for driver_method in ["minimal", "internal"]:
                exp.ctrl_pystablemotif_trapspace(target_method, driver_method)
        for method in ["ITC", "TTC", "PTC"]:
            exp.ctrl_cabean_target_control(method)

        exp_run = SingleInputSummary.from_folder(opath, inst)
        exp_run.save(f"{opath}/_graph")
