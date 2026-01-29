if __name__ == "__main__":
    import sys
    from os.path import dirname, abspath
    libdir = dirname(dirname(abspath(__file__)))
    sys.path.insert(0, libdir)

from argparse import ArgumentParser

import os.path

from bntaxonomy.iface import load_tools, tool_names
from bntaxonomy.utils.log import main_logger, configure_logging
from bntaxonomy.experiment import ExperimentHandler
from bntaxonomy.hierarchy import SingleInputSummary
import os

def main():
    configure_logging("cli")
    load_tools()

    ap = ArgumentParser()
    ap.add_argument("--exclude-targets", action="store_true",
            help="Exclude nodes specifying the target phenotype from candidate perturbations")
    ap.add_argument("max_size", type=int,
            help="Maximum number of perturbations")
    ap.add_argument("--instances", type=str, nargs="+", default=list(),
            help="Paths to instances, necessarily in an 'instances' subdirectory")
    ap.add_argument(
        "--instance-groups", type=str, nargs="+", default=list(),
        help="Paths to folders that contain instance subfolders",
)
    
    

    ap.add_argument("--tools", choices=tool_names(), nargs="*")

    args = ap.parse_args()
    for grp in args.instance_groups:
        if not os.path.isdir(grp):
            ap.error(f"Instance group path not a directory: {grp}")
        for name in sorted(os.listdir(grp)):
            path = os.path.join(grp, name)
            if os.path.isdir(path):
                args.instances.append(path)
                
    for inst in args.instances:
        if not os.path.isdir(inst):
            main_logger.warning(f" {inst} is not a directory, ignoring")
            continue
        parts = inst.split(os.path.sep)
        if "instances" not in parts:
            main_logger.warning(f" {inst} is not an 'instances' directory, ignoring")
            continue

        main_logger.info(f"Running {inst}")

        parts[parts.index("instances")] = "results"

        fpath = inst
        opath = os.path.sep.join(parts)
        inst = os.path.basename(inst)

        main_logger.info(f"   output will be in {opath}")
        os.makedirs(opath, exist_ok=True)

        exp = ExperimentHandler(
            inst,
            fpath,
            opath,
            args.max_size,
            to_console=True,
            to_file=True,
            only_minimal=True,
            load_precompute=True,
            exclude_targets=args.exclude_targets,
        )
        exp.run_tools(args.tools)

        """
        exp.ctrl_actonet_fp()
        exp.ctrl_bonesis_fp()
        exp.ctrl_bonesis_mts()
        exp.ctrl_optboolnet_fp()
        exp.ctrl_optboolnet_sync_attr()
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
            exp.ctrl_cabean_target_control(method, _debug=True)
        """

        exp_run = SingleInputSummary.from_folder(opath, inst)
        exp_run.save(f"{opath}/_graph")

if __name__ == "__main__":
    main()
