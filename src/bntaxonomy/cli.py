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
    ap.add_argument("max_size", type=int, help="Maximum number of perturbations")
    ap.add_argument(
        "-ig",
        "--inst_groups",
        nargs="+",
        help="Instance-group directories (under 'instances'). If provided, these will be translated to corresponding 'results' folders.",
        default=list(),
    )
    ap.add_argument(
        "-i",
        "--instances",
        nargs="+",
        help="Explicit instance folders (under 'instances').",
        default=list(),
    )

    ap.add_argument("--tools", choices=tool_names(), nargs="*")

    # Other options
    ap.add_argument(
        "--exclude-targets",
        action="store_true",
        help="Exclude nodes specifying the target phenotype from candidate perturbations",
    )
    ap.add_argument(
        "--print-output",
        action="store_true",
        help="Print console output from tools (except for the final results).",
        )
    ap.add_argument(
        "--clear-cache",
        action="store_true",
        help="Clear any existing cache for each experiment before running tools.",
    )
    
    args = ap.parse_args()
    for grp in args.inst_groups:
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
            print_output=args.print_output,
            clear_cache=args.clear_cache,
        )
        exp.run_tools(args.tools)
        exp_run = SingleInputSummary.from_folder(opath, inst)
        exp_run.save(f"{opath}/_graph")


if __name__ == "__main__":
    main()
