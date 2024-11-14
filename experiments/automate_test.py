import os, sys
from bntaxonomy.experiment import ExperimentHandler, ExperimentRun

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

    max_size = 2
    inst_groups = [insts_practical, insts_tonello, insts_ce]
    # inst_groups = [insts_ce]

    for inst_selected, inst in get_insts(inst_groups):
        print(f"\tRunning {inst}")
        fpath = f"{inst_selected}/{inst}"
        opath = f"{inst_selected.replace('instances','results')}/{inst}"
        os.makedirs(opath, exist_ok=True)

        exp = ExperimentHandler(
            inst,
            fpath,
            opath,
            # use_propagated=False,
            precompute_pbn=True,
            precompute_sm=True,
            to_console=False,
            to_file=True,
            only_minimal=True,
            store_results=True,
        )

        exp.ctrl_ActoNet(max_size)
        exp.ctrl_BoNesis(max_size)
        exp.ctrl_Caspo(max_size)
        for update in ["synchronous", "asynchronous"]:
            exp.ctrl_pyboolnet_mc(max_size, update)
        for control_type in ["percolation", "trap_spaces"]:
            exp.ctrl_pyboolnet_heuristics(max_size, control_type)
        exp.ctrl_pystablemotif_brute_force(max_size)
        for target_method in ["merge"]:  # ["merge", "history"]:
            for driver_method in ["minimal", "internal"]:
                exp.ctrl_pystablemotif_trapspace(max_size, target_method, driver_method)
        for method in ["ITC", "TTC", "PTC"]:
            exp.ctrl_cabean_phenotype(max_size, method)

        exp_run = ExperimentRun.from_folder(opath, inst)
        exp_run.save(f"{opath}/_graph")
