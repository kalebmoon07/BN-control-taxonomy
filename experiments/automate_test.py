import os, sys
from bntaxonomy.experiment import ExperimentHandler

cwd = os.path.dirname(os.path.abspath(__file__))

insts_practical = f"{cwd}/instances/practical"
insts_tonello = f"{cwd}/instances/Tonello"
insts_ce = f"{cwd}/instances/counterexamples"

inst_selected = insts_ce # change here

inst_list = sorted(os.listdir(inst_selected))
max_size = 2
for inst in inst_list:
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

    exp_run = exp.get_run(inst)
    exp_run.save(f"{opath}/_graph")
